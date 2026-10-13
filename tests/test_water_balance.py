import pytest

pytest.importorskip("matplotlib")

from fao_model import water_balance   # noqa: E402


def test_report_is_written(tmp_path):
    water_balance.main(["--first-year", "1995", "--last-year", "1996", "--out", str(tmp_path)])
    report = (tmp_path / "README.md").read_text(encoding="utf-8")
    for heading in ("## 1.", "## 2.", "## 3.", "## 4.", "## 5."):
        assert heading in report
    for name in ("example_season", "seasons", "soil_moisture_era5"):
        assert (tmp_path / f"{name}.png").stat().st_size > 10_000


def test_forecast_section_compares_irrigation_with_rain_forecasts(tmp_path):
    water_balance.main(["--first-year", "2014", "--last-year", "2015", "--out", str(tmp_path)])
    report = (tmp_path / "README.md").read_text(encoding="utf-8")
    assert "## 6. Irrigation with rain forecasts (2015-2015)" in report
    assert "Perfect forecast" in report and "Damped anomaly persistence" in report
