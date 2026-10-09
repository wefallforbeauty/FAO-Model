"""Daily evapotranspiration equations. Every function returns mm/day.

Inputs are a dict of daily weather (see WeatherGenerator.ensure_required):
temperatures °C, RH %, u2 m s-1 at 2 m, Rs MJ m-2 day-1, P kPa.
The *_base functions return the equation with a unit coefficient; the
coefficient (literature value or calibrated) is applied in compute_all.
"""
import math

from .meteo import (
    LAMBDA, mean_sat_vp, actual_vp, slope_vp, psychrometric,
    extraterrestrial_radiation, day_length, annual_daylight_hours,
    sunshine_fraction, net_radiation,
)

WILLMOTT_T = 26.5   # °C, above this Thornthwaite uses Willmott et al. (1985)

# FAO-24 Blaney-Criddle b regression (Frevert et al., 1983, Table 1)
BC_FREVERT = (0.81917, -0.0040922, 1.0705, 0.065649, -0.0059864, -0.0005967)


def _net_radiation(d):
    return net_radiation(
        d["Rs"], d["T_max"], d["T_min"], actual_vp(d),
        d["lat"], d["doy"], d.get("elev", 100.0), albedo=0.23
    )


def asce_pm(d):
    """FAO-56 Penman-Monteith grass reference ETo (FAO-56 eq. 6), G = 0."""
    T = (d["T_max"] + d["T_min"]) / 2.0
    es = mean_sat_vp(d["T_max"], d["T_min"])
    ea = actual_vp(d)
    delta = slope_vp(T)
    gamma = psychrometric(d["P"])
    Rn = _net_radiation(d)
    eto = (
        0.408 * delta * Rn
        + gamma * (900.0 / (T + 273.0)) * d["u2"] * (es - ea)
    ) / (delta + gamma * (1.0 + 0.34 * d["u2"]))
    return max(0.0, eto)


def thornthwaite_heat_index(monthly_mean_temps):
    """Annual heat index I from the 12 monthly mean temperatures (Thornthwaite, 1948)."""
    return sum((T / 5.0) ** 1.514 for T in monthly_mean_temps if T > 0)


def thornthwaite_unadjusted(T, I):
    """PET of a standard month (30 days of 12 h daylight), mm/month.

    T is a monthly mean temperature. Above 26.5 °C the Willmott et al.
    (1985) form is used, which does not depend on I.
    """
    if T <= 0 or I <= 0:
        return 0.0
    if T >= WILLMOTT_T:
        return -415.85 + 32.24 * T - 0.43 * T ** 2
    a = (
        6.75e-7 * I**3
        - 7.71e-5 * I**2
        + 1.792e-2 * I
        + 0.49239
    )
    return 16.0 * (10.0 * T / I) ** a


def thornthwaite_monthly(T_month, I, mean_day_length, days):
    """Thornthwaite PET of a calendar month, mm/month."""
    return thornthwaite_unadjusted(T_month, I) * (mean_day_length / 12.0) * (days / 30.0)


def thornthwaite(d):
    """Daily rate of the monthly Thornthwaite method.

    Needs d["T_30d"] (trailing 30-day mean temperature, standing in for the
    monthly mean) and d["heat_index"] (annual heat index of the site).
    """
    dl = day_length(d["lat"], d["doy"])
    return max(0.0, thornthwaite_unadjusted(d["T_30d"], d["heat_index"]) * (dl / 12.0) / 30.0)


def hargreaves_base(d):
    """Hargreaves-Samani ETo (FAO-56 eq. 52), Ra converted to mm/day (x 0.408)."""
    Ra = extraterrestrial_radiation(d["lat"], d["doy"])
    return 0.0023 * (d["T_mean"] + 17.8) * math.sqrt(
        max(0.0, d["T_max"] - d["T_min"])
    ) * 0.408 * Ra


def turc_base(d):
    """Turc (1961); Rs converted from MJ m-2 day-1 to cal cm-2 day-1 (x 23.8846)."""
    T = d["T_mean"]
    Rs = 23.8846 * d["Rs"]
    if T <= 0:
        return 0.0
    if d["RH_mean"] >= 50:
        return 0.013 * (T / (T + 15.0)) * (Rs + 50.0)
    return (
        0.013
        * (T / (T + 15.0))
        * (Rs + 50.0)
        * (1.0 + (50.0 - d["RH_mean"]) / 70.0)
    )


def abtew_base(d):
    """Abtew (1996) without its coefficient: Rs / lambda."""
    return d["Rs"] / LAMBDA


def priestley_taylor_base(d):
    """Priestley-Taylor (1972) with alpha = 1, G = 0."""
    T = d["T_mean"]
    delta = slope_vp(T)
    gamma = psychrometric(d["P"])
    Rn = _net_radiation(d)
    return (delta / (delta + gamma)) * max(0.0, Rn) / LAMBDA


def jensen_haise_base(d):
    """Jensen-Haise (1963): (0.025 T + 0.08) Rs, Rs as evaporation equivalent (Rs / lambda)."""
    return (0.025 * d["T_mean"] + 0.08) * d["Rs"] / LAMBDA


def blaney_criddle_f(d):
    """Blaney-Criddle factor p (0.46 T + 8.13), mm/day.

    p is the daily percentage of annual daytime hours (about 0.27 at 12 h).
    """
    p = 100.0 * day_length(d["lat"], d["doy"]) / annual_daylight_hours(d["lat"])
    return p * (0.46 * d["T_mean"] + 8.13)


def blaney_criddle_ab(d):
    """FAO-24 a and b from RH_min, n/N and u2 (Frevert et al., 1983).

    n/N is estimated from Rs/Ra (Angstrom). Returns (a, b).
    """
    rh = d["RH_min"]
    n_N = sunshine_fraction(d["Rs"], extraterrestrial_radiation(d["lat"], d["doy"]))
    u = d["u2"]
    e0, e1, e2, e3, e4, e5 = BC_FREVERT
    a = 0.0043 * rh - n_N - 1.41
    b = e0 + e1 * rh + e2 * n_N + e3 * u + e4 * rh * n_N + e5 * rh * u
    return a, b


def compute_all(d, coeffs):
    """All equations for one day. bc_a/bc_b of None means the FAO-24 regression."""
    eto_asce = asce_pm(d)
    et_pm_maize = coeffs["kc_maize"] * eto_asce
    et_th = thornthwaite(d)
    if coeffs["bc_a"] is None or coeffs["bc_b"] is None:
        bc_a, bc_b = blaney_criddle_ab(d)
    else:
        bc_a, bc_b = coeffs["bc_a"], coeffs["bc_b"]
    et_bc = bc_a + bc_b * blaney_criddle_f(d)
    return {
        "ET_ASCE_PM": eto_asce,
        "ET_PM_maize": et_pm_maize,
        "ET_Thornthwaite": et_th,
        "ET_BlaneyCriddle": max(0.0, et_bc),
        "ET_Turc": max(0.0, coeffs["turc_a"] * turc_base(d)),
        "ET_PriestleyTaylor": max(0.0, coeffs["pt_alpha"] * priestley_taylor_base(d)),
        "ET_Hargreaves": max(0.0, coeffs["hargreaves_a"] * hargreaves_base(d)),
        "ET_JensenHaise": max(0.0, coeffs["jh_a"] * jensen_haise_base(d)),
        "ET_Abtew": max(0.0, coeffs["abtew_a"] * abtew_base(d)),
    }
