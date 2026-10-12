"""FAO-56 dual crop coefficients, root growth and Feddes root water uptake.

The defaults describe grain maize: stage lengths from FAO-56 Table 11
(spring planting, Spain/California), crop height from Table 12, basal crop
coefficients from Table 17 and rooting depth and depletion fraction from
Table 22. The water stress thresholds are the corn values of Wesseling
(1991) for the Feddes et al. (1978) uptake function.
"""
from dataclasses import dataclass, replace
from datetime import date

import numpy as np


def _clip(x, lo, hi):
    return max(lo, min(hi, x))


def _climate_term(u2, rh_min, h):
    """FAO-56 eqs. 62, 70 and 72: correction for wind and humidity other than 2 m/s and 45 %."""
    u2 = _clip(u2, 1.0, 6.0)
    rh_min = _clip(rh_min, 20.0, 80.0)
    return (0.04 * (u2 - 2.0) - 0.004 * (rh_min - 45.0)) * (h / 3.0) ** 0.3


@dataclass(frozen=True)
class Crop:
    name: str = "grain maize"
    planting_month: int = 5
    planting_day: int = 1
    l_ini: int = 30                 # stage lengths, days
    l_dev: int = 40
    l_mid: int = 50
    l_late: int = 30
    kcb_ini: float = 0.15
    kcb_mid: float = 1.15
    kcb_end: float = 0.15           # dried in the field; 0.50 for a high-moisture harvest
    h_ini: float = 0.05             # crop height, m
    h_max: float = 2.0
    zr_ini: float = 0.15            # rooting depth, m
    zr_max: float = 1.2             # Table 22 range: 1.0-1.7 m
    p: float = 0.55                 # depletion fraction without stress
    feddes_h1: float = -15.0        # cm; uptake starts below this head (anaerobiosis above)
    feddes_h2: float = -30.0        # optimal uptake from h2 ...
    feddes_h3_high: float = -325.0  # ... down to h3, for potential transpiration >= 0.5 cm/day
    feddes_h3_low: float = -600.0   # h3 for potential transpiration <= 0.1 cm/day
    feddes_h4: float = -8000.0      # no uptake below this head (wilting)

    @property
    def season_length(self):
        return self.l_ini + self.l_dev + self.l_mid + self.l_late

    def planting_date(self, year):
        return date(year, self.planting_month, self.planting_day)

    def day_of_season(self, day):
        """1 on the planting day, ..., season_length on the last day; None outside the season."""
        t = (day - self.planting_date(day.year)).days + 1
        return t if 1 <= t <= self.season_length else None

    def stage(self, t):
        bounds = np.cumsum([self.l_ini, self.l_dev, self.l_mid])
        return ("initial", "development", "mid-season", "late season")[int(np.searchsorted(bounds, t))]

    def basal_kc(self, t):
        """Trapezoidal basal crop coefficient Kcb on day t of the season (FAO-56 Fig. 34)."""
        s1 = self.l_ini
        s2 = s1 + self.l_dev
        s3 = s2 + self.l_mid
        if t <= s1:
            return self.kcb_ini
        if t <= s2:
            return self.kcb_ini + (self.kcb_mid - self.kcb_ini) * (t - s1) / self.l_dev
        if t <= s3:
            return self.kcb_mid
        return self.kcb_mid + (self.kcb_end - self.kcb_mid) * (t - s3) / self.l_late

    def growth_fraction(self, t):
        """0 until the end of the initial stage, rising linearly to 1 at mid-season.

        Equals (Kcb - Kcb_ini) / (Kcb_mid - Kcb_ini) during development, which
        is how pyfao56 grows crop height and roots; neither shrinks afterwards.
        """
        return _clip((t - self.l_ini) / self.l_dev, 0.0, 1.0)

    def height(self, t):
        return self.h_ini + (self.h_max - self.h_ini) * self.growth_fraction(t)

    def root_depth(self, t):
        """Rooting depth Zr, m."""
        return self.zr_ini + (self.zr_max - self.zr_ini) * self.growth_fraction(t)

    def max_kc(self, kcb, h, u2, rh_min):
        """Upper limit of Kc after rain or irrigation (FAO-56 eq. 72)."""
        return max(1.2 + _climate_term(u2, rh_min, h), kcb + 0.05)

    def canopy_cover(self, kcb, kc_max, h):
        """Fraction of the ground covered by the crop (FAO-56 eq. 76), 0-0.99."""
        if kc_max <= self.kcb_ini:
            return 0.0
        ratio = max(0.0, (kcb - self.kcb_ini) / (kc_max - self.kcb_ini))
        return _clip(ratio ** (1.0 + 0.5 * h), 0.0, 0.99)

    def daily(self, day, u2, rh_min):
        """Crop coefficients for one day; outside the season the field is bare soil.

        Returns stage, day of season t, Kcb, crop height h (m), rooting depth
        zr (m), Kc_max, canopy cover fc and exposed soil fraction few
        (FAO-56 eq. 75 with the whole surface wetted, fw = 1).
        """
        t = self.day_of_season(day)
        if t is None:
            kc_max = self.max_kc(0.0, 0.0, u2, rh_min)
            return {"stage": "bare soil", "t": None, "kcb": 0.0, "h": 0.0,
                    "zr": self.zr_ini, "kc_max": kc_max, "fc": 0.0, "few": 1.0}
        kcb = self.basal_kc(t)
        h = self.height(t)
        kc_max = self.max_kc(kcb, h, u2, rh_min)
        fc = self.canopy_cover(kcb, kc_max, h)
        return {"stage": self.stage(t), "t": t, "kcb": kcb, "h": h, "zr": self.root_depth(t),
                "kc_max": kc_max, "fc": fc, "few": _clip(1.0 - fc, 0.01, 1.0)}

    def potential_rates(self, day, eto, u2, rh_min):
        """Potential transpiration Kcb*ETo and soil evaporation Ke*ETo, mm/day.

        Ke = min(Kc_max - Kcb, few * Kc_max) (FAO-56 eq. 71) with Kr = 1:
        the soil model itself reduces evaporation as the surface dries.
        """
        c = self.daily(day, u2, rh_min)
        tp = c["kcb"] * eto
        ep = min(c["kc_max"] - c["kcb"], c["few"] * c["kc_max"]) * eto
        return c, tp, ep

    def feddes_alpha(self, h, tp):
        """Water stress factor 0-1 of root water uptake at pressure head h (cm).

        tp is the potential transpiration in cm/day; it sets h3 between the
        high- and low-demand values.
        """
        h = np.asarray(h, dtype=float)
        r_high, r_low = 0.5, 0.1
        weight = _clip((tp - r_low) / (r_high - r_low), 0.0, 1.0)
        h3 = self.feddes_h3_low + weight * (self.feddes_h3_high - self.feddes_h3_low)
        h1, h2, h4 = self.feddes_h1, self.feddes_h2, self.feddes_h4
        alpha = np.zeros_like(h)
        wet = (h <= h1) & (h > h2)
        alpha[wet] = (h1 - h[wet]) / (h1 - h2)
        alpha[(h <= h2) & (h >= h3)] = 1.0
        dry = (h < h3) & (h > h4)
        alpha[dry] = (h[dry] - h4) / (h3 - h4)
        return alpha


def root_weights(zr, dz, nz):
    """Share of root water uptake taken from each layer of thickness dz (same unit as zr).

    Root density falls linearly from the surface to zero at depth zr
    (Prasad, 1988); the shares sum to 1.
    """
    if zr > nz * dz:
        raise ValueError("rooting depth is below the bottom of the soil column")
    z = np.minimum(np.arange(nz + 1) * dz, zr)
    cumulative = z - z ** 2 / (2.0 * zr)
    return np.diff(cumulative) / (zr / 2.0)


def adjust_for_climate(crop, records):
    """Crop with Kcb_mid (and Kcb_end if >= 0.45) corrected for the site's climate (FAO-56 eq. 70).

    Uses the mean RH_min and u2 of the mid-season (late-season) days over all
    years of the records, so the result is a site constant like the
    Thornthwaite heat index.
    """
    def stage_mean(first, last):
        days = [d for d in records
                if (t := crop.day_of_season(d["timestamp"].date())) is not None and first <= t <= last]
        if not days:
            raise ValueError("the records contain no days of this crop stage")
        return (sum(d["u2"] for d in days) / len(days), sum(d["RH_min"] for d in days) / len(days))

    s2 = crop.l_ini + crop.l_dev
    s3 = s2 + crop.l_mid
    kcb_mid, kcb_end = crop.kcb_mid, crop.kcb_end
    if kcb_mid >= 0.45:
        kcb_mid += _climate_term(*stage_mean(s2 + 1, s3), crop.h_max)
    if kcb_end >= 0.45:
        kcb_end += _climate_term(*stage_mean(s3 + 1, crop.season_length), crop.h_max)
    return replace(crop, kcb_mid=kcb_mid, kcb_end=kcb_end)
