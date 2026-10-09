import math
from datetime import timedelta

import numpy as np
import pytest

from fao_model import data, et, meteo
from reference.pyfao56_refet import ascedaily


@pytest.fixture(scope="module")
def konya():
    return data.load()


def test_default_dataset_is_complete_and_daily(konya):
    meta, records = konya
    assert meta["elevation_m"] == "1013.0"
    assert len(records) == 10958    # 1995-01-01 .. 2024-12-31
    assert records[0]["timestamp"].date().isoformat() == "1995-01-01"
    for prev, d in zip(records, records[1:]):
        assert d["timestamp"] - prev["timestamp"] == timedelta(days=1)
    for d in records:
        assert d["T_min"] <= d["T_mean"] <= d["T_max"]
        assert 0.0 <= d["RH_min"] <= d["RH_mean"] <= d["RH_max"] <= 100.0
        assert d["Rs"] >= 0.0 and d["precip"] >= 0.0 and d["u2"] >= 0.0
        assert 85.0 < d["P"] < 95.0


def test_loader_converts_units(konya):
    _, records = konya
    d = records[0]      # CSV row: 13.8, 6.2 °C; u10 3.31 m/s; 898.6 hPa
    assert d["T_mean"] == pytest.approx(10.0)
    assert d["u2"] == pytest.approx(3.31 * 4.87 / math.log(67.8 * 10 - 5.42))
    assert d["P"] == pytest.approx(89.86)


def test_asce_pm_matches_pyfao56_on_real_days(konya):
    _, records = konya
    for d in records[::7]:
        d = dict(d, P=meteo.atmospheric_pressure(d["elev"]))
        ref = ascedaily("S", d["elev"], d["lat"], d["doy"], d["Rs"], d["T_max"], d["T_min"],
                        rhmax=d["RH_max"], rhmin=d["RH_min"], wndsp=d["u2"], wndht=2.0)
        assert et.asce_pm(d) == pytest.approx(max(0.0, ref), abs=0.005)


def test_asce_pm_tracks_open_meteo_et0(konya):
    # Open-Meteo sums hourly FAO-56 values; a daily computation differs
    # slightly but should follow it closely.
    _, records = konya
    ours = np.array([et.asce_pm(d) for d in records])
    ref = np.array([d["ET0_ref"] for d in records])
    assert np.corrcoef(ours, ref)[0, 1] > 0.99
    assert abs(ours.mean() / ref.mean() - 1.0) < 0.08


def test_csv_round_trip(tmp_path):
    daily = {"time": ["2020-07-01", "2020-07-02"]}
    for var in data.DAILY_VARS:
        daily[var] = [20.0, 25.0]
    daily["temperature_2m_max"] = [30.0, 32.0]
    daily["temperature_2m_min"] = [15.0, 16.0]
    daily["surface_pressure_mean"] = [900.0, 905.0]
    payload = {"latitude": 37.5, "longitude": 32.75, "elevation": 1013.0,
               "timezone": "Europe/Istanbul", "daily": daily}
    path = tmp_path / "site.csv"
    data.write_csv(payload, path, "Test site", 37.57, 32.78)
    meta, records = data.load(path)
    assert meta["site"] == "Test site"
    assert len(records) == 2
    d = records[1]
    assert d["doy"] == 184
    assert d["lat"] == 37.57 and d["elev"] == 1013.0
    assert d["T_mean"] == pytest.approx(24.0)
    assert d["P"] == pytest.approx(90.5)
    assert d["ET0_ref"] == 25.0
