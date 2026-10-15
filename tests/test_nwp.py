from datetime import date, datetime, timedelta

import numpy as np
import pytest

from fao_model import data, nwp
from fao_model.forecast import Climatology, WeatherForecaster, add_eto
from fao_model.meteo import wind_2m
from fao_model.nwp_skill import known_bias


def hourly_two_days():
    """Two local days of hourly forecasts at every lead with known daily statistics."""
    times = [f"2025-07-{d:02d}T{h:02d}:00" for d in (1, 2) for h in range(24)]
    hourly = {"time": times}
    for lead in nwp.LEADS:
        hourly[f"temperature_2m_previous_day{lead}"] = [20.0 + (h % 24) * 0.5 for h in range(48)]
        hourly[f"relative_humidity_2m_previous_day{lead}"] = [80.0 - (h % 24) for h in range(48)]
        hourly[f"wind_speed_10m_previous_day{lead}"] = [3.0] * 48
        hourly[f"shortwave_radiation_previous_day{lead}"] = [500.0 if 6 <= h % 24 < 18 else 0.0 for h in range(48)]
        hourly[f"precipitation_previous_day{lead}"] = [0.5 if h % 24 == 12 else 0.0 for h in range(48)]
        hourly[f"surface_pressure_previous_day{lead}"] = [900.0] * 48
    hourly["temperature_2m_previous_day2"][30] = None      # one missing hour on day 2, lead 2
    return hourly


def test_daily_statistics_by_lead():
    rows = nwp.daily_by_lead(hourly_two_days())
    assert len(rows) == 2 * len(nwp.LEADS) - 1             # day 2 at lead 2 is incomplete
    row = next(r for r in rows if r["date"] == "2025-07-01" and r["lead"] == 3)
    assert row["T_max"] == pytest.approx(31.5) and row["T_min"] == pytest.approx(20.0)
    assert row["RH_max"] == 80.0 and row["RH_min"] == 57.0 and row["RH_mean"] == pytest.approx(68.5)
    assert row["Rs"] == pytest.approx(500.0 * 12 * 3600 / 1e6)          # 21.6 MJ/m2
    assert row["precip"] == pytest.approx(0.5) and row["u10"] == 3.0 and row["P_hPa"] == 900.0


def test_csv_round_trip(tmp_path):
    rows = nwp.daily_by_lead(hourly_two_days())
    path = tmp_path / "model.csv"
    nwp.write_csv(rows, path, "ecmwf_ifs025", {"grid_lat": 37.5, "grid_lon": 32.75, "elevation_m": 1013.0,
                                               "timezone": "Europe/Istanbul"}, 37.57, 32.78)
    meta, table = nwp.load("ecmwf_ifs025", path)
    assert meta["grid_lat"] == "37.5"
    day = table[1][date(2025, 7, 2)]
    assert day["u2"] == pytest.approx(wind_2m(3.0, 10.0))
    assert day["P"] == pytest.approx(90.0) and day["T_mean"] == pytest.approx((31.5 + 20.0) / 2)
    assert date(2025, 7, 2) not in table[2]


def test_archived_forecasts_cover_the_evaluation_period():
    for model in nwp.MODELS:
        _, table = nwp.load(model)
        assert sorted(table) == list(nwp.LEADS)
        assert min(table[1]) == date(2024, 3, 15) and max(table[1]) == date(2026, 9, 30)


def test_nwp_backend_uses_the_archive_and_falls_back_to_the_baseline():
    _, records = data.load_all()
    records = add_eto(records)
    fc = WeatherForecaster(Climatology(records[:10592]), backend="ecmwf")
    assert fc.horizon == 7
    i = next(k for k, d in enumerate(records) if d["timestamp"].date() == date(2025, 7, 14))
    out = fc.forecast_at(records, [i])
    _, table = nwp.load("ecmwf_ifs025")
    assert out["T_max"][0].tolist() == [table[lead][date(2025, 7, 14) + timedelta(days=lead)]["T_max"]
                                        for lead in nwp.LEADS]
    assert fc.last_source == nwp.MODELS["ecmwf_ifs025"]
    fc.forecast_at(records, [400])                   # 1996: nothing archived
    assert fc.last_source.startswith("baseline")


def test_running_bias_learns_a_constant_error_without_looking_ahead():
    start = date(2025, 1, 1)
    days = [start + timedelta(days=k) for k in range(200)]
    observed = {d: {"T_max": 10.0} for d in days}
    table = {1: {d: {"T_max": 12.0} for d in days}}
    table[1][days[150]]["T_max"] = 100.0              # one wild error
    known = known_bias(table, 1, "T_max", observed, weight=0.05)
    assert known[days[149]] == pytest.approx(2.0, abs=0.01)
    assert known[days[150]] > known[days[149]] + 4.0   # it enters only from its own day on
    assert known[days[0]] == pytest.approx(0.05 * 2.0)


def test_report_without_timesfm(tmp_path):
    from fao_model import nwp_skill
    nwp_skill.main(["--no-timesfm", "--first", "2025-06-01", "--last", "2025-07-31", "--out", str(tmp_path)])
    report = (tmp_path / "README.md").read_text(encoding="utf-8")
    for heading in ("## 1.", "## 2.", "## 3.", "## 4.", "## 5."):
        assert heading in report
    assert "ECMWF IFS 0.25°, bias-corrected" in report and "NCEP GFS, raw" in report
    assert (tmp_path / "eto_error.png").stat().st_size > 10_000


def test_load_all_joins_consecutive_files_and_rejects_gaps(tmp_path):
    meta, records = data.load_all()
    assert records[0]["timestamp"] == datetime(1995, 1, 1) and records[-1]["timestamp"] == datetime(2026, 9, 30)
    assert all((b["timestamp"] - a["timestamp"]).days == 1 for a, b in zip(records, records[1:]))

    def write(path, days):
        daily = {"time": days}
        for var in data.DAILY_VARS:
            daily[var] = [20.0] * len(days)
        payload = {"latitude": 37.5, "longitude": 32.75, "elevation": 1013.0, "timezone": "UTC", "daily": daily}
        data.write_csv(payload, path, "site", 37.57, 32.78)

    write(tmp_path / "a.csv", ["2025-01-01", "2025-01-02"])
    write(tmp_path / "b.csv", ["2025-01-04"])
    with pytest.raises(ValueError):
        data.load_all((tmp_path / "a.csv", tmp_path / "b.csv"))
