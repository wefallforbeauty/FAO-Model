import numpy as np
import pytest

from fao_model import data, et
from fao_model.crop import Crop, root_weights
from fao_model.field import COLUMN, IRRIGATION_RATE, FieldModel
from fao_model.richards import RichardsSolver3D

CROP = Crop()
TP = 0.6    # cm/day


def column(top_head, bottom_head=-100.0):
    """1.5 m column whose top 30 cm are at top_head and the rest at bottom_head."""
    soil = RichardsSolver3D(**COLUMN)
    h = np.full((1, 1, soil.nz), bottom_head)
    h[:, :, :12] = top_head
    soil.set_head(h)
    return soil


def uptake(soil, omega_c, dt=1e-4):
    weights = root_weights(120.0, soil.dz, soil.nz)
    out = soil.step(0.0, dt=dt, tp=TP, root_weights=weights,
                    alpha=lambda h: CROP.feddes_alpha(h, TP), omega_c=omega_c)
    return out["transpiration"] / dt, weights


def test_without_compensation_uptake_follows_feddes():
    soil = column(-3000.0)
    rate, weights = uptake(soil, omega_c=1.0)
    omega = float(np.sum(CROP.feddes_alpha(column(-3000.0).h[0, 0], TP) * weights))
    assert rate == pytest.approx(TP * omega, rel=1e-3)


def test_compensation_makes_up_for_a_dry_topsoil():
    rates = {omega_c: uptake(column(-3000.0), omega_c)[0] for omega_c in (1.0, 0.8, 0.5)}
    assert rates[1.0] < 0.9 * TP
    assert rates[1.0] < rates[0.8] <= rates[0.5] <= TP * (1 + 1e-9)
    assert rates[0.5] == pytest.approx(TP, rel=1e-3)


def test_below_the_critical_index_uptake_is_omega_over_omega_c():
    soil = column(-7000.0, bottom_head=-3000.0)
    weights = root_weights(120.0, soil.dz, soil.nz)
    omega = float(np.sum(CROP.feddes_alpha(soil.h[0, 0], TP) * weights))
    assert omega < 0.8
    rate, _ = uptake(soil, omega_c=0.8)
    assert rate == pytest.approx(TP * omega / 0.8, rel=1e-2)
    with pytest.raises(ValueError):
        uptake(column(-100.0), omega_c=0.0)


def test_segments_keep_the_water_balance():
    soil = RichardsSolver3D(**COLUMN, h_init=-400.0)
    start = soil.storage
    out = soil.run([(0.1, 9.6), (0.9, -0.4)])
    assert out["infiltration"] == pytest.approx(0.96, abs=1e-6)
    change = soil.storage - start
    assert change == pytest.approx(out["infiltration"] - out["evaporation"] - out["drainage"], abs=1e-9)
    assert abs(soil.mass_balance_error) < 1e-9


def test_irrigation_takes_amount_over_rate_hours():
    _, records = data.load()
    days = [d for d in records if d["timestamp"].year == 1996]
    for rate in (IRRIGATION_RATE, None):
        field = FieldModel(irrigation=True, irrigation_rate=rate)
        start = 10.0 * field.soil.storage
        outs = [field.step(d, et.asce_pm(d)) for d in days]
        events = [o for o in outs if o["irrigation"] > 0]
        assert events
        for o in events:
            expected = 24.0 if rate is None else min(24.0, o["irrigation"] / rate)
            assert o["irrigation_hours"] == pytest.approx(expected)
        total = lambda key: sum(o[key] for o in outs)   # noqa: E731
        residual = (total("precip") + total("irrigation") - total("ea") - total("ta")
                    - total("drainage") - total("runoff") - (10.0 * field.soil.storage - start))
        assert abs(residual) < 1e-8


def test_defaults():
    assert CROP.omega_c == 0.8 and IRRIGATION_RATE == 6.0
    assert FieldModel().irrigation_rate == IRRIGATION_RATE


def test_report_is_written(tmp_path):
    pytest.importorskip("matplotlib")
    from fao_model import field_model_report
    field_model_report.main(["--first-year", "2010", "--last-year", "2010", "--out", str(tmp_path)])
    report = (tmp_path / "README.md").read_text(encoding="utf-8")
    for heading in ("## 1.", "## 2.", "## 3."):
        assert heading in report
    assert "(used now)" in report and "| FAO-56 | 0.55 |" in report
    for name in ("dry_down", "wetting"):
        assert (tmp_path / f"{name}.png").stat().st_size > 10_000
