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
