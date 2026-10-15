"""Daily weather forecasts: TimesFM 2.5 when installed, otherwise a statistical baseline.

The weather inputs (temperature, humidity, wind, radiation, rain) are
forecast one variable at a time, and reference ET for the coming days is
computed from the forecast weather with ASCE-PM, so it follows the same
physics as the rest of the model. TimesFM also forecasts the ETo series
directly, which gives quantiles for an uncertainty band.

The baseline is climatology plus today's anomaly decaying with its lag-1
autocorrelation (damped anomaly persistence), the usual reference for
daily weather forecasts.
"""
from datetime import timedelta

import numpy as np

from .et import asce_pm
from .meteo import extraterrestrial_radiation

WEATHER = ("T_max", "T_min", "RH_max", "RH_min", "u2", "Rs", "precip")
HORIZON = 14                    # days
CONTEXT = 1024                  # days of history TimesFM sees
TIMESFM_REPO = "google/timesfm-2.5-200m-pytorch"
# Archived weather-model forecasts (see nwp.py), by backend name
NWP_BACKENDS = {"ecmwf": "ecmwf_ifs025", "gfs": "gfs_seamless"}


def add_eto(records):
    """Copies of the records with their ASCE-PM reference ET under "ETo"."""
    return [dict(d, ETo=asce_pm(d)) for d in records]


def _doy(d):
    return d["timestamp"].timetuple().tm_yday


class Climatology:
    """Smoothed day-of-year mean of each variable and the lag-1 autocorrelation of its anomalies."""

    def __init__(self, records, keys=WEATHER + ("ETo",), half_window=15):
        if len(records) < 365:
            raise ValueError("a climatology needs at least a year of daily records")
        doy = np.array([_doy(d) for d in records])
        self.means, self.phi = {}, {}
        for key in keys:
            values = np.array([d[key] for d in records], dtype=float)
            mean = np.empty(367)
            for day in range(1, 367):
                distance = np.abs(doy - day)
                distance = np.minimum(distance, 366 - distance)
                mean[day] = values[distance <= half_window].mean()
            mean[0] = mean[1]
            anomaly = values - mean[doy]
            self.means[key] = mean
            phi = np.corrcoef(anomaly[:-1], anomaly[1:])[0, 1] if anomaly.std() > 0 else 0.0
            self.phi[key] = float(np.clip(phi, 0.0, 0.99))

    def mean(self, key, doy):
        return self.means[key][doy]


def _future_doys(last_doy, horizon):
    return [(last_doy + k - 1) % 365 + 1 for k in range(1, horizon + 1)]


class WeatherForecaster:
    """Forecasts the weather inputs and reference ET for the next `horizon` days.

    backend: "timesfm" (fail if it cannot load), "baseline", "auto"
    (TimesFM when the package and model load, else the baseline), or
    "ecmwf" / "gfs": archived weather-model forecasts for 1-7 days ahead,
    with the baseline on days that have none.
    """

    def __init__(self, climatology, backend="auto", horizon=HORIZON, context=CONTEXT, batch_size=8):
        self.climatology = climatology
        self.horizon = horizon
        self.context = context
        self.model = None
        self.nwp = None
        self.last_source = None
        self.status = "Baseline: climatology plus damped anomaly persistence."
        if backend in NWP_BACKENDS:
            self.nwp = NWPForecasts(NWP_BACKENDS[backend])
            self.horizon = self.nwp.horizon
            self.status = (f"{self.nwp.name} forecasts archived by Open-Meteo, 1-{self.horizon} days; "
                           "the baseline where none exist.")
        elif backend in ("auto", "timesfm"):
            try:
                self.model = TimesFM(context, horizon, batch_size)
                self.status = f"TimesFM 2.5 (200M), {context}-day context, {horizon}-day horizon."
            except Exception as exc:
                if backend == "timesfm":
                    raise
                self.status = (f"TimesFM not available ({type(exc).__name__}); "
                               "baseline: climatology plus damped anomaly persistence.")

    def forecast_at(self, records, origins):
        """Forecasts issued at the end of each day records[i], i in origins, for days i+1..i+horizon.

        records need the weather variables and "ETo" (see add_eto). Returns a
        dict of arrays shaped (len(origins), horizon): every weather variable,
        "ETo" computed from the forecast weather, "ETo_direct" (the ETo series
        forecast itself), "precip_mean" (expected rain; the "precip" point
        forecast of TimesFM is a median, mostly 0 on a dry climate) and, with
        TimesFM, "ETo_q10" and "ETo_q90".
        """
        keys = WEATHER + ("ETo",)
        series = {k: np.array([d[k] for d in records], dtype=float) for k in keys}
        origins = list(origins)
        if self.nwp is not None:
            out = self._baseline(records, series, origins, keys)
            used = 0
            for j, i in enumerate(origins):
                issued = records[i]["timestamp"].date()
                days = [self.nwp.at(issued, lead) for lead in range(1, self.horizon + 1)]
                if all(days):
                    used += 1
                    for key in WEATHER:
                        out[key][j] = [d[key] for d in days]
            out["precip_mean"] = out["precip"].copy()
            self.last_source = self.nwp.name if used == len(origins) else (
                "baseline (no archived forecast)" if used == 0 else f"{self.nwp.name} for {used} of {len(origins)}")
        elif self.model is None:
            out = self._baseline(records, series, origins, keys)
            out["precip_mean"] = out["precip"].copy()
        else:
            contexts = [series[k][max(0, i + 1 - self.context): i + 1] for i in origins for k in keys]
            points, quantiles = self.model.forecast(contexts, self.horizon)
            n = len(keys)
            out = {k: points[j::n] for j, k in enumerate(keys)}
            out["ETo_q10"] = quantiles[keys.index("ETo")::n, :, 1]
            out["ETo_q90"] = quantiles[keys.index("ETo")::n, :, 9]
            out["precip_mean"] = np.maximum(quantiles[keys.index("precip")::n, :, 0], 0.0)
        out["ETo_direct"] = out.pop("ETo")
        _make_consistent(out, [records[i] for i in origins])
        out["ETo"] = np.array([
            _eto_from_weather(records[i], {k: out[k][j] for k in WEATHER})
            for j, i in enumerate(origins)
        ])
        return out

    def rain_totals(self, records, origins, days):
        """Expected rain (mm) over days i+1..i+days for each origin i; cheaper than forecast_at."""
        rain = np.array([d["precip"] for d in records], dtype=float)
        origins = list(origins)
        if self.model is None:
            out = self._baseline(records, {"precip": rain}, origins, ("precip",))["precip"]
        else:
            contexts = [rain[max(0, i + 1 - self.context): i + 1] for i in origins]
            out = self.model.forecast(contexts, self.horizon)[1][:, :, 0]
        return np.maximum(out[:, :days], 0.0).sum(axis=1)

    def _baseline(self, records, series, origins, keys):
        clim = self.climatology
        out = {k: np.empty((len(origins), self.horizon)) for k in keys}
        lead = np.arange(1, self.horizon + 1)
        for j, i in enumerate(origins):
            doy0 = _doy(records[i])
            doys = _future_doys(doy0, self.horizon)
            for k in keys:
                anomaly = series[k][i] - clim.mean(k, doy0)
                out[k][j] = clim.means[k][doys] + anomaly * clim.phi[k] ** lead
        return out


def _make_consistent(fc, lasts):
    """Repair what independent univariate forecasts can get wrong, in place."""
    lo, hi = np.minimum(fc["T_min"], fc["T_max"]), np.maximum(fc["T_min"], fc["T_max"])
    fc["T_min"], fc["T_max"] = lo, hi
    rh_lo = np.clip(np.minimum(fc["RH_min"], fc["RH_max"]), 1.0, 100.0)
    rh_hi = np.clip(np.maximum(fc["RH_min"], fc["RH_max"]), 1.0, 100.0)
    fc["RH_min"], fc["RH_max"] = rh_lo, rh_hi
    fc["u2"] = np.maximum(fc["u2"], 0.0)
    fc["precip"] = np.maximum(fc["precip"], 0.0)
    for j, last in enumerate(lasts):
        doys = _future_doys(_doy(last), fc["Rs"].shape[1])
        rso = np.array([(0.75 + 2e-5 * last["elev"]) * extraterrestrial_radiation(last["lat"], d) for d in doys])
        fc["Rs"][j] = np.clip(fc["Rs"][j], 0.0, rso)


def _eto_from_weather(last, weather):
    """ASCE-PM reference ET (mm/day) for each forecast day."""
    doys = _future_doys(_doy(last), len(weather["T_max"]))
    eto = []
    for k, doy in enumerate(doys):
        d = {key: float(weather[key][k]) for key in WEATHER}
        d.update(T_mean=(d["T_max"] + d["T_min"]) / 2.0, RH_mean=(d["RH_max"] + d["RH_min"]) / 2.0,
                 P=last["P"], lat=last["lat"], elev=last["elev"], doy=doy)
        eto.append(asce_pm(d))
    return eto


class NWPForecasts:
    """Archived daily forecasts of one weather model for 1-7 days ahead (see nwp.py)."""

    def __init__(self, model, path=None):
        from . import nwp

        self.name = nwp.MODELS[model]
        self.meta, self.table = nwp.load(model, path)
        self.horizon = max(self.table)

    def at(self, issued, lead):
        """The forecast for day issued + lead, made `lead` days before it; None if not archived."""
        return self.table.get(lead, {}).get(issued + timedelta(days=lead))


class TimesFM:
    """TimesFM 2.5 (200M parameters, PyTorch): zero-shot forecasts of many series at once.

    Needs the optional `timesfm` and `torch` packages; the weights (~1 GB) are
    downloaded from Hugging Face on first use.
    """

    def __init__(self, context=CONTEXT, horizon=HORIZON, batch_size=8):
        import timesfm

        self._timesfm = timesfm
        self.model = timesfm.TimesFM_2p5_200M_torch.from_pretrained(TIMESFM_REPO)
        self.compile(context, horizon, batch_size)

    def compile(self, context, horizon, batch_size):
        """(Re)compile for another context length or batch size without reloading the weights."""
        self.model.compile(self._timesfm.ForecastConfig(
            max_context=context,
            max_horizon=horizon,
            normalize_inputs=True,
            use_continuous_quantile_head=True,
            force_flip_invariance=True,
            infer_is_positive=True,
            fix_quantile_crossing=True,
            per_core_batch_size=batch_size,
        ))

    def forecast(self, series, horizon):
        """Median forecasts (n, horizon) and quantiles (n, horizon, 10): mean, then 0.1 ... 0.9."""
        # TimesFM pads the list it is given in place, so pass a fresh one.
        points, quantiles = self.model.forecast(horizon=horizon, inputs=[np.asarray(s, dtype=float) for s in series])
        return np.asarray(points, dtype=float), np.asarray(quantiles, dtype=float)
