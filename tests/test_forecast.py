import math
from datetime import datetime, timedelta

import numpy as np
import pytest

from fao_model import data
from fao_model.et import asce_pm
from fao_model.forecast import (
    TIMESFM_REPO, WEATHER, Climatology, WeatherForecaster, _eto_from_weather, _make_consistent, add_eto,
)


def synthetic_records(years=6, phi=0.7, seed=0):
    """Weather with a seasonal cycle and AR(1) anomalies of known persistence."""
    rng = np.random.default_rng(seed)
    start = datetime(2001, 1, 1)
    anomaly = 0.0
    records = []
    for i in range(365 * years):
        day = start + timedelta(days=i)
        doy = day.timetuple().tm_yday
        season = math.sin(2 * math.pi * (doy - 100) / 365.25)
        anomaly = phi * anomaly + rng.normal(0.0, 2.0)
        t_mean = 15.0 + 10.0 * season + anomaly
        records.append({
            "timestamp": day, "doy": doy, "lat": 37.5, "elev": 1000.0, "P": 90.0,
            "T_max": t_mean + 6.0, "T_min": t_mean - 6.0, "T_mean": t_mean,
            "RH_max": 80.0 - 10.0 * season, "RH_min": 40.0 - 10.0 * season, "RH_mean": 60.0 - 10.0 * season,
            "u2": 2.0, "Rs": 18.0 + 8.0 * season, "precip": max(0.0, rng.normal(0.0, 2.0)),
        })
    return add_eto(records)


def test_climatology_recovers_seasonal_mean_and_persistence():
    records = synthetic_records()
    clim = Climatology(records)
    for doy in (15, 100, 192, 300):
        season = math.sin(2 * math.pi * (doy - 100) / 365.25)
        assert clim.mean("T_max", doy) == pytest.approx(21.0 + 10.0 * season, abs=1.2)
    assert clim.phi["T_max"] == pytest.approx(0.7, abs=0.06)


def test_baseline_decays_from_todays_anomaly_to_climatology():
    records = synthetic_records()
    clim = Climatology(records)
    fc = WeatherForecaster(clim, backend="baseline", horizon=30)
    i = len(records) - 40
    out = fc.forecast_at(records, [i])
    doy0 = records[i]["timestamp"].timetuple().tm_yday
    anomaly = records[i]["T_max"] - clim.mean("T_max", doy0)
    assert out["T_max"][0, 0] == pytest.approx(clim.mean("T_max", doy0 % 365 + 1) + anomaly * clim.phi["T_max"])
    assert out["T_max"][0, -1] == pytest.approx(clim.mean("T_max", (doy0 + 29) % 365 + 1), abs=0.1)
    for key in WEATHER + ("ETo", "ETo_direct"):
        assert out[key].shape == (1, 30), key
    assert "ETo_q10" not in out


def test_forecast_is_made_physically_consistent():
    last = {"timestamp": datetime(2010, 7, 1), "lat": 37.5, "elev": 1000.0}
    fc = {
        "T_max": np.array([[20.0, 30.0]]), "T_min": np.array([[25.0, 10.0]]),
        "RH_max": np.array([[120.0, 50.0]]), "RH_min": np.array([[60.0, 70.0]]),
        "u2": np.array([[-0.5, 2.0]]), "Rs": np.array([[50.0, 20.0]]), "precip": np.array([[-1.0, 3.0]]),
    }
    _make_consistent(fc, [last])
    assert fc["T_min"][0].tolist() == [20.0, 10.0] and fc["T_max"][0].tolist() == [25.0, 30.0]
    assert fc["RH_max"][0].tolist() == [100.0, 70.0] and fc["RH_min"][0].tolist() == [60.0, 50.0]
    assert fc["u2"][0, 0] == 0.0 and fc["precip"][0, 0] == 0.0
    assert fc["Rs"][0, 0] < 35.0      # clipped to clear-sky radiation in July at 37.5° N


def test_eto_from_forecast_weather_is_asce_pm():
    last = {"timestamp": datetime(2010, 7, 1), "lat": 37.5, "elev": 1000.0, "P": 90.0}
    weather = {"T_max": [32.0], "T_min": [16.0], "RH_max": [55.0], "RH_min": [20.0],
               "u2": [2.5], "Rs": [28.0], "precip": [0.0]}
    d = {"T_max": 32.0, "T_min": 16.0, "T_mean": 24.0, "RH_max": 55.0, "RH_min": 20.0, "RH_mean": 37.5,
         "u2": 2.5, "Rs": 28.0, "P": 90.0, "lat": 37.5, "elev": 1000.0, "doy": 183}
    assert _eto_from_weather(last, weather)[0] == pytest.approx(asce_pm(d))


def test_baseline_on_real_data_beats_persistence_for_eto_at_one_week():
    _, records = data.load()
    records = add_eto(records)
    clim = Climatology(records[:7305])
    fc = WeatherForecaster(clim, backend="baseline", horizon=7)
    origins = list(range(7305, len(records) - 7, 10))
    out = fc.forecast_at(records, origins)
    obs = np.array([[records[i + k]["ETo"] for k in range(1, 8)] for i in origins])
    persistence = np.array([[records[i]["ETo"]] * 7 for i in origins])
    assert np.abs(out["ETo_direct"][:, 6] - obs[:, 6]).mean() < np.abs(persistence[:, 6] - obs[:, 6]).mean()


def _timesfm_cached():
    try:
        import timesfm  # noqa: F401
        from huggingface_hub import try_to_load_from_cache
    except ImportError:
        return False
    return isinstance(try_to_load_from_cache(TIMESFM_REPO, "model.safetensors"), str)


@pytest.mark.skipif(not _timesfm_cached(), reason="TimesFM or its downloaded weights are not available")
def test_timesfm_forecasts_weather_and_eto():
    _, records = data.load()
    records = add_eto(records)
    fc = WeatherForecaster(Climatology(records[:7305]), backend="timesfm", horizon=14)
    out = fc.forecast_at(records, [7400, 7500])
    for key in WEATHER + ("ETo", "ETo_direct", "ETo_q10", "ETo_q90"):
        assert out[key].shape == (2, 14) and np.all(np.isfinite(out[key])), key
    assert np.all(out["ETo_q10"] <= out["ETo_direct"] + 1e-6)
    assert np.all(out["ETo_direct"] <= out["ETo_q90"] + 1e-6)
    assert np.all(out["T_min"] <= out["T_max"])


def test_climatology_needs_a_year():
    with pytest.raises(ValueError):
        Climatology(synthetic_records(years=1)[:200])
