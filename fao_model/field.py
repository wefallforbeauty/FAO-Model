"""Daily water balance of a cropped field.

FAO-56 dual crop coefficients split the reference ET into potential
transpiration (Kcb * ETo) and potential soil evaporation (Ke * ETo). The
Richards solver then decides how much of each the soil can actually
supply: roots take water through the Feddes stress function, evaporation
stops when the surface dries. Optional automatic irrigation refills the
root zone when its depletion reaches a threshold.
"""
import numpy as np

from .crop import Crop, root_weights
from .richards import RichardsSolver3D

H_FIELD_CAPACITY = -100.0   # cm (10 kPa); gives FAO-56 loam values for this soil
H_WILTING_POINT = -15000.0  # cm (1500 kPa)
# 1.5 m soil column in 2.5 cm layers. Soil evaporation depends on the top layer
# thickness (5 cm layers gave ~30 % more); the field is laterally uniform, so one
# column gives the same result as a 3D grid.
COLUMN = dict(nx=1, ny=1, nz=60, dz=2.5)


class FieldModel:
    """Soil, crop and irrigation for one field, advanced one day at a time.

    irrigation: True to irrigate automatically during the season whenever
    the root-zone depletion fraction reaches `mad` (management allowed
    depletion, default the crop's p), refilling the root zone to field
    capacity on that day; False for a rainfed field.
    """

    def __init__(self, crop=None, irrigation=False, soil=None, mad=None):
        self.crop = crop or Crop()
        self.soil = soil or RichardsSolver3D(**COLUMN)
        self.irrigation = irrigation
        self.mad = self.crop.p if mad is None else mad
        self.theta_fc = float(self.soil._theta(np.array(H_FIELD_CAPACITY)))
        self.theta_wp = float(self.soil._theta(np.array(H_WILTING_POINT)))
        self.last = None

    def root_zone(self, zr_cm):
        """Depletion below field capacity Dr and total available water TAW of the top zr_cm, mm."""
        dz = self.soil.dz
        tops = np.arange(self.soil.nz) * dz
        overlap = np.clip(zr_cm - tops, 0.0, dz)
        theta = self.soil.theta.mean(axis=(0, 1))
        dr = float(np.sum((self.theta_fc - theta) * overlap)) * 10.0
        taw = float(np.sum((self.theta_fc - self.theta_wp) * overlap)) * 10.0
        return max(dr, 0.0), taw

    def step(self, d, eto):
        """Advance one day of weather d with reference ET eto (mm). Returns the day's terms in mm."""
        day = d["timestamp"].date()
        c, tp, ep = self.crop.potential_rates(day, eto, d["u2"], d["RH_min"])
        zr_cm = 100.0 * c["zr"]
        dr, taw = self.root_zone(zr_cm)
        irrigation = 0.0
        if self.irrigation and c["t"] is not None and dr >= self.mad * taw:
            irrigation = dr
        q_top = (d["precip"] + irrigation - ep) / 10.0
        out = self.soil.step(
            q_top, dt=1.0, tp=tp / 10.0,
            root_weights=root_weights(zr_cm, self.soil.dz, self.soil.nz),
            alpha=lambda h: self.crop.feddes_alpha(h, tp / 10.0),
        )
        dr_end, taw_end = self.root_zone(zr_cm)
        self.last = dict(
            c,
            eto=eto,
            precip=d["precip"],
            irrigation=irrigation,
            tp=tp,
            ep=ep,
            ta=10.0 * out["transpiration"],
            ea=ep - 10.0 * out["evaporation_deficit"],
            drainage=10.0 * out["drainage"],
            runoff=10.0 * out["runoff"],
            dr=dr_end,
            taw=taw_end,
            storage=10.0 * self.soil.storage,
        )
        return self.last
