from datetime import datetime

import numpy as np
import pytest

from fao_model import data, et
from fao_model.crop import Crop
from fao_model.field import H_FIELD_CAPACITY, FieldModel
from fao_model.richards import RichardsSolver3D


@pytest.fixture(scope="module")
def year_1996():
    _, records = data.load()
    days = [d for d in records if d["timestamp"].year == 1996]
    return days, [et.asce_pm(d) for d in days]


def test_field_water_balance_closes(year_1996):
    days, eto = year_1996
    for irrigation in (False, True):
        field = FieldModel(irrigation=irrigation)
        start = 10.0 * field.soil.storage
        outs = [field.step(d, e) for d, e in zip(days, eto)]
        total = lambda key: sum(o[key] for o in outs)
        change = 10.0 * field.soil.storage - start
        residual = (total("precip") + total("irrigation") - total("ea") - total("ta")
                    - total("drainage") - total("runoff") - change)
        assert abs(residual) < 1e-8
        assert total("ta") <= total("tp") + 1e-9
        assert total("ea") <= total("ep") + 1e-9


def test_irrigation_refills_the_root_zone_in_season_only(year_1996):
    days, eto = year_1996
    field = FieldModel(irrigation=True)
    events = []
    for d, e in zip(days, eto):
        day = d["timestamp"].date()
        zr_cm = 100.0 * field.crop.daily(day, d["u2"], d["RH_min"])["zr"]
        dr, taw = field.root_zone(zr_cm)
        out = field.step(d, e)
        if out["irrigation"] > 0:
            events.append(day)
            assert dr >= field.mad * taw
            assert out["irrigation"] == pytest.approx(dr)
    assert events
    assert all(field.crop.day_of_season(day) is not None for day in events)
    rainfed = FieldModel(irrigation=False)
    assert all(rainfed.step(d, e)["irrigation"] == 0.0 for d, e in zip(days, eto))


def test_root_zone_depletion_at_field_capacity():
    field = FieldModel()
    field.soil.set_head(np.full(field.soil.h.shape, H_FIELD_CAPACITY))
    dr, taw = field.root_zone(100.0)
    assert dr == pytest.approx(0.0, abs=1e-9)
    assert taw == pytest.approx((field.theta_fc - field.theta_wp) * 1000.0)


def test_irrigation_beats_rainfed_transpiration(year_1996):
    days, eto = year_1996
    def season_transpiration(field):
        outs = [field.step(d, e) for d, e in zip(days, eto)]
        return sum(o["ta"] for o in outs if o["t"] is not None)
    wet = season_transpiration(FieldModel(irrigation=True))
    dry = season_transpiration(FieldModel(irrigation=False))
    assert wet > 3.0 * dry


def test_laterally_uniform_grid_matches_one_column(year_1996):
    days, eto = year_1996
    days, eto = days[120:180], eto[120:180]     # May and June
    one = FieldModel(irrigation=True, soil=RichardsSolver3D(nx=1, ny=1, nz=60, dz=2.5))
    grid = FieldModel(irrigation=True, soil=RichardsSolver3D(nx=3, ny=2, nz=60, dz=2.5))
    for d, e in zip(days, eto):
        a, b = one.step(d, e), grid.step(d, e)
        for key in ("ta", "ea", "irrigation", "drainage", "dr"):
            assert a[key] == pytest.approx(b[key], rel=1e-9, abs=1e-12)


def test_bare_soil_outside_the_season():
    field = FieldModel(crop=Crop())
    d = {"timestamp": datetime(2001, 1, 10), "u2": 2.0, "RH_min": 60.0, "precip": 0.0}
    out = field.step(d, 1.0)
    assert out["stage"] == "bare soil" and out["tp"] == 0.0 and out["ta"] == 0.0


def test_rain_forecast_postpones_irrigation(year_1996):
    days, eto = year_1996
    days, eto = days[120:270], eto[120:270]          # the 1996 season
    base = FieldModel(irrigation=True)
    wet_forecast = FieldModel(irrigation=True, rain_forecast=lambda d: 50.0)
    base_out = [base.step(d, e) for d, e in zip(days, eto)]
    wet_out = [wet_forecast.step(d, e) for d, e in zip(days, eto)]
    assert any(o["postponed"] for o in wet_out)
    assert not any(o["postponed"] for o in base_out)
    # Postponing stops at max_depletion: the crop is still irrigated.
    assert sum(o["irrigation"] for o in wet_out) > 0.0
    for o, d in zip(wet_out, days):
        if o["irrigation"] > 0:
            assert o["dr"] >= 0.0
    first_base = next(i for i, o in enumerate(base_out) if o["irrigation"] > 0)
    first_wet = next(i for i, o in enumerate(wet_out) if o["irrigation"] > 0)
    assert first_wet > first_base
