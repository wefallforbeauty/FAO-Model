"""Yield response to water: the FAO-33 water production function.

Doorenbos and Kassam (1979, FAO Irrigation and Drainage Paper 33) relate the
relative yield loss to the relative evapotranspiration deficit,

    1 - Ya/Yx = Ky (1 - ETa/ETx),

where Yx and Ya are the maximum and actual yields, ETx is the crop ET with
the crop's water needs fully met, ETa the actual crop ET, and Ky the yield
response factor. For maize, Ky is 1.25 over the whole season and, by growth
period, 0.40 (vegetative), 1.50 (flowering), 0.50 (yield formation) and 0.20
(ripening) (Steduto et al., 2012, FAO Irrigation and Drainage Paper 66,
Tables 1 and 2).

- **Growth periods:** relative yields of periods with their own Ky are
  multiplied (Rao et al., 1988). The periods are placed on the FAO-56 crop
  stages:
  - vegetative: the initial and development stages;
  - flowering: starts with the mid-season stage (FAO-56: "effective full
    cover for many crops occurs at the initiation of flowering") and is
    assumed to last 15 days;
  - yield formation: the rest of the mid-season stage;
  - ripening: the late-season stage.
- **Which ET:** with dual crop coefficients the deficit is a transpiration
  deficit. Soil evaporation does not depend on the crop's water status, so it
  enters ETx and ETa alike: ETx = Tp + Ea, ETa = Ta + Ea, as in FAO-56, where
  Ke does not depend on Ks. basis="transpiration" uses Ta / Tp instead.

The relation is empirical, and FAO-66 notes that the errors in actual yields
estimated this way are "quite important", so the results are best read as
relative yield losses.
"""
SEASON_KY = 1.25
MAIZE_KY = {"vegetative": 0.40, "flowering": 1.50, "yield formation": 0.50, "ripening": 0.20}
# FAO-66 Table 2 lists the same values against the FAO-56 stages instead.
MAIZE_KY_BY_FAO56_STAGE = {"initial": 0.40, "development": 1.50, "mid-season": 0.50, "late season": 0.20}
FLOWERING_DAYS = 15
BASES = ("et", "transpiration")


def growth_period(crop, t, flowering_days=FLOWERING_DAYS):
    """FAO-33 growth period of day t of the season (see the module docstring)."""
    full_cover = crop.l_ini + crop.l_dev        # last day of the development stage
    if t <= full_cover:
        return "vegetative"
    if t <= full_cover + flowering_days:
        return "flowering"
    if t <= full_cover + crop.l_mid:
        return "yield formation"
    return "ripening"


def fao56_stage(crop, t):
    """FAO-56 crop stage of day t, for the FAO-66 Table 2 placement of Ky."""
    return crop.stage(t)


def _et(day, basis):
    """(ETa, ETx) of one field day, mm."""
    if basis == "et":
        return day["ta"] + day["ea"], day["tp"] + day["ea"]
    if basis == "transpiration":
        return day["ta"], day["tp"]
    raise ValueError(f"basis must be one of {BASES}")


def water_use(days, crop, period=growth_period, basis="et"):
    """{period: (ETa, ETx)} in mm, summed over the in-season days of one season.

    days: FieldModel.step results (with "t", "tp", "ta" and "ea").
    """
    totals = {}
    for day in days:
        if day["t"] is None:
            continue
        name = period(crop, day["t"])
        eta, etx = _et(day, basis)
        a, x = totals.get(name, (0.0, 0.0))
        totals[name] = (a + eta, x + etx)
    return totals


def period_yield(eta, etx, ky):
    """Relative yield Ya/Yx, 0-1, from one period's ET and its Ky."""
    if etx <= 0.0:
        return 1.0
    return min(1.0, max(0.0, 1.0 - ky * (1.0 - eta / etx)))


def stage_yield(totals, ky=MAIZE_KY):
    """Ya/Yx as the product of the periods' relative yields."""
    result = 1.0
    for name, (eta, etx) in totals.items():
        result *= period_yield(eta, etx, ky[name])
    return result


def season_yield(totals, ky=SEASON_KY):
    """Ya/Yx from the season's total ET with the seasonal Ky."""
    eta = sum(a for a, _ in totals.values())
    etx = sum(x for _, x in totals.values())
    return period_yield(eta, etx, ky)


class SeasonYield:
    """Running yield estimate through the season, fed one field day at a time.

    relative_yield() is the estimate from the season so far, as if the rest of
    the season had no water stress; last_season keeps the result of the last
    complete season.
    """

    def __init__(self, crop, ky=MAIZE_KY, period=growth_period, basis="et"):
        self.crop, self.ky, self.period, self.basis = crop, ky, period, basis
        self.totals = {}
        self.last_season = None

    def add(self, day):
        t = day["t"]
        if t is None:
            return
        if t == 1:
            self.totals = {}
        name = self.period(self.crop, t)
        eta, etx = _et(day, self.basis)
        a, x = self.totals.get(name, (0.0, 0.0))
        self.totals[name] = (a + eta, x + etx)
        if t == self.crop.season_length:
            self.last_season = self.relative_yield()

    def relative_yield(self):
        return stage_yield(self.totals, self.ky)

    def et_ratio(self):
        eta = sum(a for a, _ in self.totals.values())
        etx = sum(x for _, x in self.totals.values())
        return eta / etx if etx > 0.0 else 1.0
