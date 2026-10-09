from datetime import datetime, timedelta

from fao_model import meteo
from fao_model.weather import WeatherGenerator


def test_generates_consistent_daily_weather():
    gen = WeatherGenerator(lat=37.57, elev=1013.0, start=datetime(2001, 1, 1), seed=7)
    days = [gen.ensure_required(gen.generate()) for _ in range(730)]
    for prev, d in zip(days, days[1:]):
        assert d["timestamp"] - prev["timestamp"] == timedelta(days=1)
    for d in days:
        assert d["T_min"] <= d["T_mean"] <= d["T_max"]
        assert d["RH_min"] <= d["RH_mean"] <= d["RH_max"] <= 100.0
        Rso = (0.75 + 2e-5 * d["elev"]) * meteo.extraterrestrial_radiation(d["lat"], d["doy"])
        assert 0.0 < d["Rs"] <= Rso
        assert d["precip"] >= 0.0
    annual_precip = sum(d["precip"] for d in days) / 2
    assert 100.0 < annual_precip < 1500.0


def test_seed_makes_runs_reproducible():
    a = WeatherGenerator(start=datetime(2001, 1, 1), seed=11)
    b = WeatherGenerator(start=datetime(2001, 1, 1), seed=11)
    assert [a.generate() for _ in range(50)] == [b.generate() for _ in range(50)]


def test_heat_index_of_seasonal_climate():
    # Monthly means between ~9 and ~31 °C give a heat index around 110.
    assert 90.0 < WeatherGenerator().heat_index() < 130.0
