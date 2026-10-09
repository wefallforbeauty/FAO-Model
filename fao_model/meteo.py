"""FAO-56 meteorological helpers (Allen et al., 1998).

Units: temperature °C, pressure kPa, radiation MJ m-2 day-1, wind m s-1.
"""
import math
from functools import lru_cache

SIGMA = 4.903e-9    # Stefan-Boltzmann constant, MJ K-4 m-2 day-1
GSC = 0.0820        # solar constant, MJ m-2 min-1
LAMBDA = 2.45       # latent heat of vaporization, MJ kg-1 (1 mm water = 2.45 MJ m-2)


def sat_vp(T):
    """Saturation vapour pressure at T, kPa (FAO-56 eq. 11)."""
    return 0.6108 * math.exp(17.27 * T / (T + 237.3))


def mean_sat_vp(T_max, T_min):
    """Mean saturation vapour pressure of a day, kPa (FAO-56 eq. 12)."""
    return (sat_vp(T_max) + sat_vp(T_min)) / 2.0


def actual_vp(d):
    """Actual vapour pressure, kPa.

    Uses RH_max and RH_min (FAO-56 eq. 17) when both are present,
    otherwise RH_mean (eq. 19).
    """
    if d.get("RH_max") is not None and d.get("RH_min") is not None:
        return (sat_vp(d["T_min"]) * d["RH_max"]
                + sat_vp(d["T_max"]) * d["RH_min"]) / 200.0
    return mean_sat_vp(d["T_max"], d["T_min"]) * d["RH_mean"] / 100.0


def slope_vp(T):
    """Slope of the saturation vapour pressure curve, kPa °C-1 (FAO-56 eq. 13)."""
    es = sat_vp(T)
    return 4098 * es / (T + 237.3) ** 2


def psychrometric(P):
    """Psychrometric constant, kPa °C-1, for pressure P in kPa (FAO-56 eq. 8)."""
    return 0.000665 * P


def atmospheric_pressure(elev):
    """Atmospheric pressure at elevation elev [m], kPa (FAO-56 eq. 7)."""
    return 101.3 * ((293.0 - 0.0065 * elev) / 293.0) ** 5.26


def wind_2m(u, z):
    """Convert wind speed measured at height z [m] to 2 m, m s-1 (FAO-56 eq. 47)."""
    return u * 4.87 / math.log(67.8 * z - 5.42)


def _sunset_hour_angle(phi, decl):
    cos_ws = -math.tan(phi) * math.tan(decl)
    cos_ws = max(-1.0, min(1.0, cos_ws))
    return math.acos(cos_ws)


def extraterrestrial_radiation(lat, doy):
    """Daily extraterrestrial radiation Ra, MJ m-2 day-1 (FAO-56 eq. 21)."""
    phi = math.radians(lat)
    dr = 1 + 0.033 * math.cos(2 * math.pi * doy / 365)
    decl = 0.409 * math.sin(2 * math.pi * doy / 365 - 1.39)
    ws = _sunset_hour_angle(phi, decl)
    Ra = (24 * 60 / math.pi) * GSC * dr * (
        ws * math.sin(phi) * math.sin(decl)
        + math.cos(phi) * math.cos(decl) * math.sin(ws)
    )
    return max(0.0, Ra)


def day_length(lat, doy):
    """Daylight hours N (FAO-56 eq. 34)."""
    phi = math.radians(lat)
    decl = 0.409 * math.sin(2 * math.pi * doy / 365 - 1.39)
    ws = _sunset_hour_angle(phi, decl)
    return 24 * ws / math.pi


@lru_cache(maxsize=64)
def annual_daylight_hours(lat):
    """Sum of daylight hours over a 365-day year, h."""
    return sum(day_length(lat, doy) for doy in range(1, 366))


def sunshine_fraction(Rs, Ra, a_s=0.25, b_s=0.50):
    """Relative sunshine duration n/N estimated from Rs/Ra by inverting the
    Angstrom formula (FAO-56 eq. 35), clipped to [0, 1]."""
    if Ra <= 0:
        return 0.0
    return max(0.0, min(1.0, (Rs / Ra - a_s) / b_s))


def net_radiation(Rs, T_max, T_min, ea, lat, doy, elev=100.0, albedo=0.23):
    """Net radiation Rn = Rns - Rnl, MJ m-2 day-1 (FAO-56 eqs. 37-40).

    Rs/Rso is limited to [0.3, 1.0] as in ASCE-EWRI (2005).
    """
    Ra = extraterrestrial_radiation(lat, doy)
    Rso = (0.75 + 2e-5 * elev) * Ra
    if Rso <= 0:
        Rso = Ra
    ratio = Rs / Rso if Rso > 0 else 1.0
    ratio = max(0.3, min(1.0, ratio))

    Rns = (1 - albedo) * Rs
    Rnl = (
        SIGMA
        * ((T_max + 273.16) ** 4 + (T_min + 273.16) ** 4)
        / 2.0
        * (0.34 - 0.14 * math.sqrt(max(0.0, ea)))
        * (1.35 * ratio - 0.35)
    )
    return Rns - Rnl
