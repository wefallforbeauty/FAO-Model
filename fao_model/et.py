import math

from .meteo import (
    LAMBDA, sat_vp, slope_vp, psychrometric, extraterrestrial_radiation,
    day_length, net_radiation,
)


def asce_pm(d):
    T = d["T_mean"]
    es = sat_vp(T)
    ea = es * d["RH_mean"] / 100.0
    delta = slope_vp(T)
    gamma = psychrometric(d["P"])
    Rn = net_radiation(
        d["Rs"], d["T_max"], d["T_min"], ea,
        d["lat"], d["doy"], d.get("elev", 100.0), albedo=0.23
    )
    eto = (
        0.408 * delta * Rn
        + gamma * (900.0 / (T + 273.0)) * d["u2"] * (es - ea)
    ) / (delta + gamma * (1.0 + 0.34 * d["u2"]))
    return max(0.0, eto)


def thornthwaite(d):
    T = d["T_mean"]
    if T <= 0:
        return 0.0
    I = (T / 5.0) ** 1.514
    a = (
        6.75e-7 * I**3
        - 7.71e-5 * I**2
        + 1.792e-2 * I
        + 0.49239
    )
    dl = day_length(d["lat"], d["doy"])
    pet = 16.0 * (10.0 * T / I) ** a * (dl / 12.0)
    return max(0.0, pet)


def hargreaves_base(d):
    Ra = extraterrestrial_radiation(d["lat"], d["doy"])
    return 0.0023 * (d["T_mean"] + 17.8) * math.sqrt(
        max(0.1, d["T_max"] - d["T_min"])
    ) * Ra


def turc_base(d):
    T = d["T_mean"]
    Rs = d["Rs"]
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


def priestley_taylor_base(d):
    T = d["T_mean"]
    es = sat_vp(T)
    ea = es * d["RH_mean"] / 100.0
    delta = slope_vp(T)
    gamma = psychrometric(d["P"])
    Rn = net_radiation(
        d["Rs"], d["T_max"], d["T_min"], ea,
        d["lat"], d["doy"], d.get("elev", 100.0), albedo=0.23
    )
    return (delta / (delta + gamma)) * max(0.0, Rn) / LAMBDA


def jensen_haise_base(d):
    return (0.025 * d["T_mean"] + 0.08) * d["Rs"] / 28.3


def compute_all(d, coeffs):
    eto_asce = asce_pm(d)
    et_pm_maize = coeffs["kc_maize"] * eto_asce
    et_th = thornthwaite(d)
    dl = day_length(d["lat"], d["doy"])
    p = dl / 24.0
    f = p * (0.46 * d["T_mean"] + 8.13)
    et_bc = coeffs["bc_a"] + coeffs["bc_b"] * f
    return {
        "ET_ASCE_PM": eto_asce,
        "ET_PM_maize": et_pm_maize,
        "ET_Thornthwaite": et_th,
        "ET_BlaneyCriddle": max(0.0, et_bc),
        "ET_Turc": max(0.0, coeffs["turc_a"] * turc_base(d)),
        "ET_PriestleyTaylor": max(0.0, coeffs["pt_alpha"] * priestley_taylor_base(d)),
        "ET_Hargreaves": max(0.0, coeffs["hargreaves_a"] * hargreaves_base(d)),
        "ET_JensenHaise": max(0.0, coeffs["jh_a"] * jensen_haise_base(d)),
        "ET_Abtew": max(0.0, coeffs["abtew_a"] * d["Rs"] / LAMBDA),
    }
