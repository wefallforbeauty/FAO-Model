import pytest

from fao_model import meteo
from fao56 import EXAMPLE_18 as EX18, EXAMPLE_18_RESULTS as R18


def test_pressure_and_psychrometric_constant_example_2():
    P = meteo.atmospheric_pressure(1800.0)
    assert P == pytest.approx(81.8, abs=0.05)
    assert meteo.psychrometric(P) == pytest.approx(0.054, abs=0.0005)


def test_saturation_vapour_pressure_example_3():
    assert meteo.mean_sat_vp(24.5, 15.0) == pytest.approx(2.39, abs=0.005)


def test_actual_vapour_pressure_example_5():
    d = {"T_max": 25.0, "T_min": 18.0, "RH_max": 82.0, "RH_min": 54.0, "RH_mean": 68.0}
    assert meteo.actual_vp(d) == pytest.approx(1.70, abs=0.005)
    d_mean_only = {"T_max": 25.0, "T_min": 18.0, "RH_mean": 68.0}
    assert meteo.actual_vp(d_mean_only) == pytest.approx(1.78, abs=0.005)


def test_extraterrestrial_radiation_and_day_length_examples_8_9():
    # 3 September at 20°S
    assert meteo.extraterrestrial_radiation(-20.0, 246) == pytest.approx(32.2, abs=0.05)
    assert meteo.day_length(-20.0, 246) == pytest.approx(11.7, abs=0.05)


def test_wind_conversion_example_14():
    assert meteo.wind_2m(3.2, 10.0) == pytest.approx(2.4, abs=0.01)


def test_example_18_intermediate_values():
    d = EX18
    assert meteo.mean_sat_vp(d["T_max"], d["T_min"]) == pytest.approx(R18["es"], abs=0.002)
    ea = meteo.actual_vp(d)
    assert ea == pytest.approx(R18["ea"], abs=0.002)
    assert meteo.slope_vp(d["T_mean"]) == pytest.approx(R18["delta"], abs=0.0005)
    assert meteo.psychrometric(d["P"]) == pytest.approx(R18["gamma"], abs=0.0001)
    Ra = meteo.extraterrestrial_radiation(d["lat"], d["doy"])
    assert Ra == pytest.approx(R18["Ra"], abs=0.02)
    assert meteo.day_length(d["lat"], d["doy"]) == pytest.approx(R18["N"], abs=0.05)
    Rn = meteo.net_radiation(d["Rs"], d["T_max"], d["T_min"], ea, d["lat"], d["doy"], d["elev"])
    assert Rn == pytest.approx(R18["Rn"], abs=0.02)


def test_annual_daylight_hours_close_to_half_a_year():
    # Daylight averages ~12 h over a year at any latitude outside the polar circles.
    for lat in (0.0, 37.57, 50.8, -33.0):
        assert meteo.annual_daylight_hours(lat) == pytest.approx(365 * 12, rel=0.01)


def test_sunshine_fraction_inverts_angstrom():
    Ra = 41.09
    n_N = 9.25 / 16.1
    Rs = (0.25 + 0.5 * n_N) * Ra
    assert meteo.sunshine_fraction(Rs, Ra) == pytest.approx(n_N)
    assert meteo.sunshine_fraction(0.0, Ra) == 0.0
    assert meteo.sunshine_fraction(Ra, Ra) == 1.0
