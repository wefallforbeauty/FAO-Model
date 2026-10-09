from datetime import datetime

import pytest

from fao_model import et
from fao_model.calibration import Calibrator, LITERATURE_COEFFS
from fao_model.weather import WeatherGenerator


def _days(n=200, seed=3):
    gen = WeatherGenerator(lat=37.57, elev=1013.0, start=datetime(2010, 3, 1), seed=seed)
    return [gen.ensure_required(gen.generate()) for _ in range(n)]


def test_starts_from_literature_coefficients():
    assert Calibrator().coeffs == LITERATURE_COEFFS


def test_too_few_rows_leave_coefficients_unchanged():
    cal = Calibrator()
    cal.calibrate(_days(Calibrator.MIN_ROWS - 1))
    assert cal.coeffs == LITERATURE_COEFFS


def test_recovers_known_multiplicative_coefficients():
    days = _days()
    cal = Calibrator()
    cal.calibrate(days, reference=[1.15 * et.hargreaves_base(d) for d in days])
    assert cal.coeffs["hargreaves_a"] == pytest.approx(1.15)
    cal.calibrate(days, reference=[0.9 * et.priestley_taylor_base(d) for d in days])
    assert cal.coeffs["pt_alpha"] == pytest.approx(0.9)


def test_recovers_known_blaney_criddle_line():
    days = _days()
    cal = Calibrator()
    cal.calibrate(days, reference=[0.4 + 0.9 * et.blaney_criddle_f(d) for d in days])
    assert cal.coeffs["bc_a"] == pytest.approx(0.4)
    assert cal.coeffs["bc_b"] == pytest.approx(0.9)


def test_default_reference_is_asce_pm():
    days = _days()
    a, b = Calibrator(), Calibrator()
    a.calibrate(days)
    b.calibrate(days, reference=[et.asce_pm(d) for d in days])
    assert a.coeffs == b.coeffs
