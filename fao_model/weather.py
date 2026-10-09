import math
import random
from datetime import datetime, timedelta

from .meteo import extraterrestrial_radiation


class WeatherGenerator:
    def __init__(self, lat=40.0, elev=100.0, interval_minutes=10):
        self.lat = lat
        self.elev = elev
        self.interval_minutes = interval_minutes
        self.timestamp = datetime.now()
        self.current = None
        self.rng = random.Random()

    def _clip(self, x, lo, hi):
        return max(lo, min(hi, x))

    def _solar_fraction(self, ts):
        hour = ts.hour + ts.minute / 60.0
        season = math.sin(2 * math.pi * (ts.timetuple().tm_yday - 80) / 365.25)
        sunrise = 6.0 - 1.4 * season
        sunset = 18.0 + 1.4 * season
        if hour <= sunrise or hour >= sunset:
            return 0.0
        x = (hour - sunrise) / max(1e-6, sunset - sunrise)
        return math.sin(math.pi * x)

    def _seasonal_temperature(self, ts):
        return 20.0 + 11.0 * math.sin(
            2 * math.pi * (ts.timetuple().tm_yday - 100) / 365.25
        )

    def _new_rain_state(self, rain):
        if rain > 0.05:
            if self.rng.random() < 0.88:
                return self._clip(rain * self.rng.uniform(0.75, 1.15), 0.0, 35.0)
            return 0.0
        if self.rng.random() < 0.06:
            return self.rng.uniform(0.1, 4.0)
        return 0.0

    def generate(self):
        ts = self.timestamp
        doy = ts.timetuple().tm_yday
        seasonal_T = self._seasonal_temperature(ts)
        solar = self._solar_fraction(ts)
        if self.current is None:
            T_mean = seasonal_T + self.rng.uniform(-1.5, 1.5)
            T_range = self.rng.uniform(6.0, 12.0)
            RH_mean = self.rng.uniform(55.0, 80.0)
            u2 = self.rng.uniform(1.0, 5.0)
            rain = self.rng.uniform(0.0, 1.0)
        else:
            prev = self.current
            T_target = seasonal_T + 2.0 * (solar - 0.45)
            T_mean = 0.94 * prev["T_mean"] + 0.06 * T_target + self.rng.gauss(0, 0.35)
            T_mean = self._clip(T_mean, -20.0, 45.0)
            T_range = 0.90 * (prev["T_max"] - prev["T_min"]) + self.rng.gauss(0, 0.25)
            T_range = self._clip(T_range, 5.0, 16.0)
            humidity_target = 78.0 - 0.55 * (T_mean - seasonal_T) - 18.0 * solar
            RH_mean = 0.92 * prev["RH_mean"] + 0.08 * humidity_target + self.rng.gauss(0, 0.9)
            RH_mean = self._clip(RH_mean, 25.0, 98.0)
            wind_target = 1.0 + 5.0 * solar
            u2 = 0.90 * prev["u2"] + 0.10 * wind_target + self.rng.gauss(0, 0.18)
            u2 = self._clip(u2, 0.1, 12.0)
            rain = self._new_rain_state(prev["precip"])
        T_max = T_mean + 0.5 * T_range
        T_min = T_mean - 0.5 * T_range
        RH_max = self._clip(RH_mean + 6.0 + abs(self.rng.gauss(0, 2.0)), RH_mean, 100.0)
        RH_min = self._clip(RH_mean - 6.0 - abs(self.rng.gauss(0, 2.0)), 5.0, RH_mean)
        Ra = extraterrestrial_radiation(self.lat, doy)
        clear_sky = max(0.0, 0.75 * Ra)
        cloud_factor = self._clip(
            0.92 - 0.004 * max(0.0, RH_mean - 55.0) + self.rng.gauss(0, 0.035),
            0.20, 1.0
        )
        Rs = self._clip(clear_sky * solar * cloud_factor, 0.0, max(0.1, clear_sky))
        P = 101.3 * (1.0 - 2.25577e-5 * self.elev) ** 5.25588
        self.current = {
            "timestamp": ts,
            "doy": doy,
            "lat": self.lat,
            "elev": self.elev,
            "T_mean": T_mean,
            "T_max": max(T_mean, T_max),
            "T_min": min(T_mean, T_min),
            "RH_mean": RH_mean,
            "RH_max": RH_max,
            "RH_min": RH_min,
            "u2": u2,
            "Rs": Rs,
            "P": P,
            "precip": rain,
        }

        self.timestamp += timedelta(minutes=self.interval_minutes)
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
