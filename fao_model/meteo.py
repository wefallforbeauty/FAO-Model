import math

SIGMA = 4.903e-9
GSC = 0.0820
LAMBDA = 2.45


def sat_vp(T):
    return 0.6108 * math.exp(17.27 * T / (T + 237.3))


def slope_vp(T):
    es = sat_vp(T)
    return 4098 * es / (T + 237.3) ** 2


def psychrometric(P):
    return 0.000665 * P


def extraterrestrial_radiation(lat, doy):
    phi = math.radians(lat)
    dr = 1 + 0.033 * math.cos(2 * math.pi * doy / 365)
    decl = 0.409 * math.sin(2 * math.pi * doy / 365 - 1.39)
    cos_ws = -math.tan(phi) * math.tan(decl)
    cos_ws = max(-1.0, min(1.0, cos_ws))
    ws = math.acos(cos_ws)
    Ra = (24 * 60 / math.pi) * GSC * dr * (
        ws * math.sin(phi) * math.sin(decl)
        + math.cos(phi) * math.cos(decl) * math.sin(ws)
    )
    return max(0.0, Ra)


def day_length(lat, doy):
    phi = math.radians(lat)
    decl = 0.409 * math.sin(2 * math.pi * doy / 365 - 1.39)
    cos_ws = -math.tan(phi) * math.tan(decl)
    cos_ws = max(-1.0, min(1.0, cos_ws))
    ws = math.acos(cos_ws)
    return 24 * ws / math.pi


def net_radiation(Rs, T_max, T_min, ea, lat, doy, elev=100.0, albedo=0.23):
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
