import math
from datetime import datetime

import pytest

from fao_model import et, meteo
from fao_model.calibration import LITERATURE_COEFFS
from fao_model.weather import WeatherGenerator
from fao56 import EXAMPLE_18 as EX18, EXAMPLE_18_RESULTS as R18
from reference.pyfao56_refet import ascedaily

# Brussels-like monthly mean temperatures for the Thornthwaite heat index
BRUSSELS_MONTHLY_T = [3.0, 4.0, 7.0, 9.0, 13.0, 16.0, 18.0, 17.5, 15.0, 11.0, 6.0, 4.0]


def test_asce_pm_fao56_example_18():
    assert et.asce_pm(EX18) == pytest.approx(R18["ETo"], abs=0.05)


def test_asce_pm_matches_pyfao56_reference():
    gen = WeatherGenerator(lat=37.57, elev=1013.0, start=datetime(2001, 1, 1), seed=1)
    for _ in range(730):
        d = gen.ensure_required(gen.generate())
        ref = ascedaily("S", d["elev"], d["lat"], d["doy"], d["Rs"], d["T_max"], d["T_min"],
                        rhmax=d["RH_max"], rhmin=d["RH_min"], wndsp=d["u2"], wndht=2.0)
        assert et.asce_pm(d) == pytest.approx(max(0.0, ref), abs=0.005)


def test_hargreaves_uses_ra_in_mm_per_day():
    d = EX18
    Ra_mm = 0.408 * R18["Ra"]
    expected = 0.0023 * (d["T_mean"] + 17.8) * math.sqrt(d["T_max"] - d["T_min"]) * Ra_mm
    assert et.hargreaves_base(d) == pytest.approx(expected, rel=1e-3)


def test_turc_uses_rs_in_cal_per_cm2():
    d = dict(EX18)
    T = d["T_mean"]
    langley = d["Rs"] / 0.041868        # 1 cal cm-2 = 0.041868 MJ m-2
    expected = 0.013 * T / (T + 15.0) * (langley + 50.0)
    assert et.turc_base(d) == pytest.approx(expected, rel=1e-4)
    d["RH_mean"] = 30.0
    assert et.turc_base(d) == pytest.approx(expected * (1.0 + 20.0 / 70.0), rel=1e-4)


def test_jensen_haise_and_abtew_use_rs_as_evaporation_equivalent():
    d = EX18
    Rs_mm = d["Rs"] / 2.45
    assert et.jensen_haise_base(d) == pytest.approx((0.025 * d["T_mean"] + 0.08) * Rs_mm)
    assert et.abtew_base(d) == pytest.approx(Rs_mm)


def test_priestley_taylor_example_18():
    d = EX18
    ratio = R18["delta"] / (R18["delta"] + R18["gamma"])
    assert et.priestley_taylor_base(d) == pytest.approx(ratio * R18["Rn"] / 2.45, rel=0.005)


def test_blaney_criddle_p_is_percentage_of_annual_daytime_hours():
    d = EX18
    p = 100.0 * R18["N"] / (365 * 12)
    assert et.blaney_criddle_f(d) == pytest.approx(p * (0.46 * d["T_mean"] + 8.13), rel=0.01)


def test_blaney_criddle_fao24_coefficients():
    d = EX18
    n_N = 9.25 / 16.1
    rh, u = d["RH_min"], d["u2"]
    a_expected = 0.0043 * rh - n_N - 1.41
    b_expected = (0.81917 - 0.0040922 * rh + 1.0705 * n_N + 0.065649 * u
                  - 0.0059864 * rh * n_N - 0.0005967 * rh * u)
    a, b = et.blaney_criddle_ab(d)
    assert a == pytest.approx(a_expected, abs=0.005)
    assert b == pytest.approx(b_expected, abs=0.005)


def test_thornthwaite_heat_index_and_monthly_form():
    assert et.thornthwaite_heat_index([20.0] * 12) == pytest.approx(12 * 4.0 ** 1.514)
    assert et.thornthwaite_heat_index([-5.0, 0.0] + [5.0] * 10) == pytest.approx(10.0)
    I = et.thornthwaite_heat_index(BRUSSELS_MONTHLY_T)
    # A standard month: 30 days with 12 h daylight
    assert et.thornthwaite_monthly(16.0, I, 12.0, 30) == pytest.approx(
        et.thornthwaite_unadjusted(16.0, I))
    assert et.thornthwaite_unadjusted(0.0, I) == 0.0
    assert et.thornthwaite_unadjusted(-3.0, I) == 0.0


def test_thornthwaite_willmott_branch_is_continuous():
    for I in (40.0, 60.0, 80.0, 100.0, 120.0):
        below = et.thornthwaite_unadjusted(26.49, I)
        above = et.thornthwaite_unadjusted(26.5, I)
        assert above == pytest.approx(below, rel=0.05)


def test_thornthwaite_daily_rate_is_mm_per_day():
    d = dict(EX18, T_30d=EX18["T_mean"],
             heat_index=et.thornthwaite_heat_index(BRUSSELS_MONTHLY_T))
    monthly = et.thornthwaite_unadjusted(d["T_30d"], d["heat_index"])
    assert et.thornthwaite(d) == pytest.approx(monthly * (R18["N"] / 12.0) / 30.0, rel=0.005)
    # The old version returned ~155 for a day like this (mm/month plotted as mm/day).
    assert 2.0 < et.thornthwaite(d) < 6.0


def test_literature_equations_have_the_magnitude_of_reference_et():
    # Unit-conversion errors show up as factors of 2.45, 23.9 etc.
    d = dict(EX18, T_30d=EX18["T_mean"],
             heat_index=et.thornthwaite_heat_index(BRUSSELS_MONTHLY_T))
    results = et.compute_all(d, LITERATURE_COEFFS)
    eto = results["ET_ASCE_PM"]
    for key, value in results.items():
        if key == "ET_PM_maize":
            continue
        assert 0.6 * eto < value < 1.5 * eto, key
