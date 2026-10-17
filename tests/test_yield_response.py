from datetime import date

import pytest

from fao_model import data, et
from fao_model import yield_response as yr
from fao_model.crop import Crop
from fao_model.field import FieldModel

CROP = Crop()


def season(deficit_in=None, ratio=0.8, tp=5.0, ea=1.0):
    """Field days of one season with transpiration ratio `ratio` in period `deficit_in` and none elsewhere."""
    days = []
    for t in range(1, CROP.season_length + 1):
        ta = tp * ratio if yr.growth_period(CROP, t) == deficit_in else tp
        days.append({"t": t, "tp": tp, "ta": ta, "ea": ea})
    return days


def test_growth_periods_cover_the_season_in_order():
    periods = [yr.growth_period(CROP, t) for t in range(1, CROP.season_length + 1)]
    assert [periods.count(p) for p in yr.MAIZE_KY] == [70, 15, 35, 30]
    assert periods == sorted(periods, key=list(yr.MAIZE_KY).index)
    stages = {yr.fao56_stage(CROP, t) for t in range(1, CROP.season_length + 1)}
    assert stages == set(yr.MAIZE_KY_BY_FAO56_STAGE)


def test_no_deficit_gives_the_maximum_yield():
    totals = yr.water_use(season() + [{"t": None, "tp": 0.0, "ta": 0.0, "ea": 3.0}], CROP)
    assert yr.stage_yield(totals) == 1.0 and yr.season_yield(totals) == 1.0


def test_a_deficit_at_flowering_costs_the_most():
    # Without soil evaporation, ETa/ETx in the stressed period equals the transpiration ratio.
    losses = {p: 1.0 - yr.stage_yield(yr.water_use(season(p, ea=0.0), CROP)) for p in yr.MAIZE_KY}
    assert losses == pytest.approx({p: 0.2 * ky for p, ky in yr.MAIZE_KY.items()})
    totals = yr.water_use(season("flowering", ea=0.0), CROP)
    assert yr.season_yield(totals) == pytest.approx(1.0 - 1.25 * 0.2 * 15 / 150)


def test_period_yields_multiply_and_stay_between_0_and_1():
    totals = {"vegetative": (50.0, 100.0), "flowering": (80.0, 100.0)}
    assert yr.stage_yield(totals) == pytest.approx((1 - 0.4 * 0.5) * (1 - 1.5 * 0.2))
    assert yr.stage_yield({"flowering": (10.0, 100.0)}) == 0.0
    assert yr.period_yield(0.0, 0.0, 1.5) == 1.0


def test_soil_evaporation_counts_in_both_et_terms():
    day = [{"t": 1, "tp": 4.0, "ta": 2.0, "ea": 2.0}]
    assert yr.water_use(day, CROP) == {"vegetative": (4.0, 6.0)}
    assert yr.water_use(day, CROP, basis="transpiration") == {"vegetative": (2.0, 4.0)}
    with pytest.raises(ValueError):
        yr.water_use(day, CROP, basis="evaporation")


def test_running_estimate_restarts_each_season():
    estimate = yr.SeasonYield(CROP)
    first, second = season("flowering"), season("ripening", ratio=0.5)
    for day in first + [{"t": None, "tp": 0.0, "ta": 0.0, "ea": 1.0}] * 3:
        estimate.add(day)
    assert estimate.last_season == pytest.approx(yr.stage_yield(yr.water_use(first, CROP)))
    for day in second[:80]:
        estimate.add(day)
    assert estimate.relative_yield() == 1.0 and estimate.et_ratio() == 1.0   # no deficit yet
    for day in second[80:]:
        estimate.add(day)
    assert estimate.last_season == pytest.approx(yr.stage_yield(yr.water_use(second, CROP)))


def test_irrigated_maize_yields_far_more_than_rainfed():
    _, records = data.load()
    days = [d for d in records if d["timestamp"].year == 1996]
    relative = {}
    for irrigation in (False, True):
        field = FieldModel(irrigation=irrigation)
        outs = [field.step(d, et.asce_pm(d)) for d in days]
        relative[irrigation] = yr.stage_yield(yr.water_use(outs, field.crop))
    assert relative[False] < 0.2 < 0.7 < relative[True] < 1.0


def test_report_is_written(tmp_path):
    pytest.importorskip("matplotlib")
    from fao_model import yield_report
    yield_report.main(["--first-year", "2009", "--last-year", "2010", "--out", str(tmp_path)])
    report = (tmp_path / "README.md").read_text(encoding="utf-8")
    for heading in ("## 1.", "## 2.", "## 3."):
        assert heading in report
    assert "| 2010 |" in report and "By growth period (FAO-33)" in report
    assert "| Without oxygen stress |" in report
    for name in ("relative_yield", "production_function"):
        assert (tmp_path / f"{name}.png").stat().st_size > 10_000


def test_dashboard_keeps_a_running_yield_estimate():
    from fao_model.simulation import Simulation
    _, records = data.load()
    sim = Simulation(records, forecaster="baseline", start=date(2001, 4, 29))
    for _ in range(5):
        sim.step()
    assert sim.field.last["t"] == 3
    assert set(sim.season_yield.totals) == {"vegetative"}
    assert 0.0 < sim.season_yield.relative_yield() <= 1.0
