import math
import random
from datetime import datetime, timedelta

from .meteo import atmospheric_pressure, extraterrestrial_radiation


class WeatherGenerator:
    """Synthetic daily weather: seasonal cycles plus day-to-day persistence.

    One generate() call is one day, so the daily ET equations get inputs
    with the time scale they were derived for.
    """

    def __init__(self, lat=40.0, elev=100.0, start=None, seed=None):
        self.lat = lat
        self.elev = elev
        self.timestamp = start or datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)
        self.current = None
        self.rng = random.Random(seed)

    def _clip(self, x, lo, hi):
        return max(lo, min(hi, x))

    def _seasonal_temperature(self, doy):
        return 20.0 + 11.0 * math.sin(2 * math.pi * (doy - 100) / 365.25)

    def heat_index(self):
        """Thornthwaite annual heat index of the generator's seasonal climate."""
        from .et import thornthwaite_heat_index
        mid_month_doys = [15 + 30.4 * m for m in range(12)]
        return thornthwaite_heat_index(self._seasonal_temperature(d) for d in mid_month_doys)

    def _new_rain(self, wet_yesterday, doy):
        # Two-state Markov chain, wetter in winter than in summer.
        season = 1.0 + 0.6 * math.cos(2 * math.pi * (doy - 15) / 365.25)
        p_wet = (0.45 if wet_yesterday else 0.15) * season / 1.6
        if self.rng.random() < p_wet:
            return self._clip(self.rng.expovariate(1.0 / 6.0), 0.1, 60.0)
        return 0.0

    def generate(self):
        ts = self.timestamp
        doy = ts.timetuple().tm_yday
        seasonal_T = self._seasonal_temperature(doy)
        if self.current is None:
            T_anom = self.rng.uniform(-1.5, 1.5)
            T_range = self.rng.uniform(8.0, 14.0)
            RH_anom = 0.0
            u2 = self.rng.uniform(1.0, 3.0)
            rain = 0.0
        else:
            prev = self.current
            T_anom = 0.7 * (prev["T_mean"] - prev["T_seasonal"]) + self.rng.gauss(0, 1.8)
            T_range = 11.0 + 0.5 * (prev["T_max"] - prev["T_min"] - 11.0) + self.rng.gauss(0, 1.5)
            RH_anom = 0.6 * prev["RH_anom"] + self.rng.gauss(0, 5.0)
            u2 = 2.0 + 0.5 * (prev["u2"] - 2.0) + self.rng.gauss(0, 0.5)
            rain = self._new_rain(prev["precip"] > 0, doy)
        T_mean = self._clip(seasonal_T + T_anom, -20.0, 40.0)
        T_range = self._clip(T_range - (3.0 if rain > 0 else 0.0), 3.0, 20.0)
        RH_mean = self._clip(
            65.0 - 15.0 * math.sin(2 * math.pi * (doy - 100) / 365.25)
            - 1.0 * T_anom + RH_anom + (12.0 if rain > 0 else 0.0),
            20.0, 98.0,
        )
        RH_max = self._clip(RH_mean + 15.0 + abs(self.rng.gauss(0, 4.0)), RH_mean, 100.0)
        RH_min = self._clip(RH_mean - 15.0 - abs(self.rng.gauss(0, 4.0)), 5.0, RH_mean)
        u2 = self._clip(u2, 0.3, 8.0)
        Rso = (0.75 + 2e-5 * self.elev) * extraterrestrial_radiation(self.lat, doy)
        cloud_factor = self._clip(
            0.85 - (0.35 if rain > 0 else 0.0) - 0.003 * max(0.0, RH_mean - 55.0)
            + self.rng.gauss(0, 0.05),
            0.25, 1.0,
        )
        self.current = {
            "timestamp": ts,
            "doy": doy,
            "lat": self.lat,
            "elev": self.elev,
            "T_mean": T_mean,
            "T_max": T_mean + 0.5 * T_range,
            "T_min": T_mean - 0.5 * T_range,
            "T_seasonal": seasonal_T,
            "RH_mean": RH_mean,
            "RH_max": RH_max,
            "RH_min": RH_min,
            "RH_anom": RH_anom,
            "u2": u2,
            "Rs": Rso * cloud_factor,
            "P": atmospheric_pressure(self.elev),
            "precip": rain,
        }
        self.timestamp += timedelta(days=1)
        return self.current.copy()

    def ensure_required(self, data):
        data["T_mean"] = float(data.get("T_mean", 20.0))
        data["T_max"] = max(data.get("T_max", data["T_mean"] + 5.0), data["T_mean"])
        data["T_min"] = min(data.get("T_min", data["T_mean"] - 5.0), data["T_mean"])
        data["RH_mean"] = self._clip(float(data.get("RH_mean", 65.0)), 0.0, 100.0)
        data["RH_max"] = self._clip(float(data.get("RH_max", data["RH_mean"] + 10)), data["RH_mean"], 100.0)
        data["RH_min"] = self._clip(float(data.get("RH_min", data["RH_mean"] - 10)), 0.0, data["RH_mean"])
        data["u2"] = self._clip(float(data.get("u2", 2.0)), 0.0, 20.0)
        data["Rs"] = max(0.0, float(data.get("Rs", 15.0)))
        data["P"] = max(50.0, float(data.get("P", 101.3)))
        data["precip"] = max(0.0, float(data.get("precip", 0.0)))
        data["doy"] = int(data.get("doy", datetime.now().timetuple().tm_yday))
        data["lat"] = float(data.get("lat", self.lat))
        data["elev"] = float(data.get("elev", self.elev))
        return data
