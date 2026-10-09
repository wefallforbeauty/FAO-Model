"""Van Genuchten-Mualem soil hydraulics and a 3D Richards equation solver.

Units: length cm, time day. Pressure head h is negative in unsaturated soil.
"""
import numpy as np


def vg_theta_arr(h, theta_r, theta_s, alpha, n):
    h_abs = np.abs(h)
    m = 1.0 - 1.0 / n
    se = (1.0 + (alpha * h_abs) ** n) ** (-m)
    se = np.where(h >= 0, 1.0, se)
    return theta_r + (theta_s - theta_r) * se


def vg_head_arr(theta, theta_r, theta_s, alpha, n):
    """Inverse of vg_theta_arr: pressure head for a water content (0 at saturation)."""
    m = 1.0 - 1.0 / n
    se = np.clip((theta - theta_r) / (theta_s - theta_r), 1e-300, 1.0)
    return -np.maximum(se ** (-1.0 / m) - 1.0, 0.0) ** (1.0 / n) / alpha


def vg_K_arr(h, Ks, theta_r, theta_s, alpha, n):
    theta = vg_theta_arr(h, theta_r, theta_s, alpha, n)
    se = (theta - theta_r) / (theta_s - theta_r)
    se = np.clip(se, 1e-8, 1.0)
    m = 1.0 - 1.0 / n
    K = Ks * np.sqrt(se) * (1.0 - (1.0 - se ** (1.0 / m)) ** m) ** 2
    K = np.where(h >= 0, Ks, K)
    return np.maximum(K, 1e-12)


def vg_capacity_arr(h, theta_r, theta_s, alpha, n):
    C = np.zeros_like(h)
    neg = h < 0
    h_abs = np.abs(h[neg])
    m = 1.0 - 1.0 / n
    C[neg] = (
        (theta_s - theta_r)
        * m * n * alpha
        * (alpha * h_abs) ** (n - 1.0)
        * (1.0 + (alpha * h_abs) ** n) ** (-m - 1.0)
    )
    return np.maximum(C, 1e-6)


def _harmonic(a, b):
    return 2.0 * a * b / (a + b + 1e-12)


class RichardsSolver3D:
    """Mass-conservative explicit finite-volume solver for the 3D Richards equation.

    Layer k = 0 is at the surface and depth increases with k. Water content
    is the conserved state and h follows from van Genuchten, so the water
    balance closes to rounding error. Boundaries:
    - surface: flux q_top (cm/day, + into the soil), limited by what a
      near-saturated surface lets in (excess becomes runoff) and by what a
      dry surface at h_min can still evaporate;
    - bottom: free drainage (unit gradient), sides: no flow.
    Sub-steps follow the explicit stability limit, which becomes very small
    in near-saturated soil; that case needs an implicit scheme.
    """

    COURANT = 0.4
    MAX_DTHETA = 0.01           # largest change of water content per sub-step
    MIN_SUBSTEP = 1e-9          # day
    MAX_SUBSTEPS = 200_000      # per call to step()
    SURFACE_SATURATION_GAP = 1e-3   # top layer stays below theta_s minus this

    def __init__(self, nx=5, ny=5, nz=5, dx=10.0, dy=10.0, dz=5.0,
                 h_init=-100.0, h_min=-20000.0):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.dx, self.dy, self.dz = dx, dy, dz
        self.h_min = h_min
        self.params = {
            "theta_r": 0.05,
            "theta_s": 0.43,
            "alpha": 0.036,     # 1/cm
            "n": 1.56,
            "Ks": 21.6,         # cm/day (0.00025 cm/s)
        }
        self.set_head(np.full((nx, ny, nz), float(h_init)))

    def set_head(self, h):
        """Set the pressure head field and restart the water balance."""
        self.theta = self._theta(np.asarray(h, dtype=float))
        self.h = self._head(self.theta)
        self.totals = dict.fromkeys(
            ("infiltration", "evaporation", "drainage", "runoff", "evaporation_deficit"), 0.0)
        self.last_step = dict(self.totals, substeps=0)
        self.initial_storage = self.storage

    def _theta(self, h):
        return vg_theta_arr(h, self.params["theta_r"], self.params["theta_s"],
                            self.params["alpha"], self.params["n"])

    def _head(self, theta):
        return vg_head_arr(theta, self.params["theta_r"], self.params["theta_s"],
                           self.params["alpha"], self.params["n"])

    def _K(self, h):
        return vg_K_arr(h, self.params["Ks"], self.params["theta_r"],
                        self.params["theta_s"], self.params["alpha"], self.params["n"])

    def _capacity(self, h):
        return vg_capacity_arr(h, self.params["theta_r"], self.params["theta_s"],
                               self.params["alpha"], self.params["n"])

    def _internal_rates(self, h, K):
        """dtheta/dt (1/day) from flow between cells and bottom drainage; bottom flux (cm/day)."""
        dx, dy, dz = self.dx, self.dy, self.dz
        rate = np.zeros_like(h)
        qx = _harmonic(K[:-1], K[1:]) * (h[:-1] - h[1:]) / dx
        rate[:-1] -= qx / dx
        rate[1:] += qx / dx
        qy = _harmonic(K[:, :-1], K[:, 1:]) * (h[:, :-1] - h[:, 1:]) / dy
        rate[:, :-1] -= qy / dy
        rate[:, 1:] += qy / dy
        # Downward flux; total head H = h - depth, so gravity adds +1 to the gradient.
        qz = _harmonic(K[:, :, :-1], K[:, :, 1:]) * ((h[:, :, :-1] - h[:, :, 1:]) / dz + 1.0)
        rate[:, :, :-1] -= qz / dz
        rate[:, :, 1:] += qz / dz
        q_bottom = K[:, :, -1]
        rate[:, :, -1] -= q_bottom / dz
        return rate, q_bottom

    def _surface_rate(self, q_top, h, K):
        """Surface flux each column can take (cm/day), before storage limits."""
        half = self.dz / 2.0
        h0, K0 = h[:, :, 0], K[:, :, 0]
        if q_top >= 0:
            # Darcy flux from a saturated surface (h = 0) to the top cell centre
            capacity = 0.5 * (self.params["Ks"] + K0) * (-h0 / half + 1.0)
            return np.minimum(q_top, np.maximum(capacity, 0.0))
        # Upward flux the top cell can deliver to a surface at h_min
        e_max = K0 * ((h0 - self.h_min) / half - 1.0)
        return -np.minimum(-q_top, np.maximum(e_max, 0.0))

    def step(self, q_top, dt=1.0):
        """Advance dt days with surface flux q_top (cm/day, + downward: rain minus ET).

        Returns this step's water balance terms in cm, averaged over the surface.
        """
        p = self.params
        dz = self.dz
        inv_d2 = 1.0 / self.dx ** 2 + 1.0 / self.dy ** 2 + 1.0 / dz ** 2
        theta_cap = p["theta_s"] - self.SURFACE_SATURATION_GAP
        theta_floor = float(self._theta(np.array(self.h_min)))
        step_totals = dict.fromkeys(self.totals, 0.0)
        t, substeps = 0.0, 0
        while t < dt * (1.0 - 1e-12):
            substeps += 1
            if substeps > self.MAX_SUBSTEPS:
                raise RuntimeError("Richards solver: sub-step limit reached "
                                   "(near-saturated soil needs an implicit scheme)")
            h, theta = self.h, self.theta
            K = self._K(h)
            rate_int, q_bottom = self._internal_rates(h, K)
            q_surf_max = self._surface_rate(q_top, h, K)
            diffusivity = K / self._capacity(h)
            sub = min(dt - t, self.COURANT / (2.0 * float(diffusivity.max()) * inv_d2))
            fastest = float(np.abs(rate_int).max() + np.abs(q_surf_max).max() / dz)
            if fastest > 0.0:
                sub = min(sub, self.MAX_DTHETA / fastest)
            while True:
                # Keep the top layer between theta(h_min) and near-saturation;
                # rain that does not fit becomes runoff, unmet demand an ET deficit.
                room_in = ((theta_cap - theta[:, :, 0]) / sub - rate_int[:, :, 0]) * dz
                room_out = ((theta_floor - theta[:, :, 0]) / sub - rate_int[:, :, 0]) * dz
                if q_top >= 0:
                    q_surf = np.minimum(q_surf_max, np.maximum(room_in, 0.0))
                else:
                    q_surf = np.maximum(q_surf_max, np.minimum(room_out, 0.0))
                rate = rate_int.copy()
                rate[:, :, 0] += q_surf / dz
                theta_new = theta + sub * rate
                if np.all(theta_new > p["theta_r"]) and np.all(theta_new < p["theta_s"]):
                    break
                sub *= 0.5
                if sub < self.MIN_SUBSTEP:
                    raise RuntimeError("Richards solver: water content left its valid range")
            self.theta = theta_new
            self.h = self._head(theta_new)
            for key, flux in (
                ("infiltration", np.maximum(q_surf, 0.0)),
                ("evaporation", np.maximum(-q_surf, 0.0)),
                ("drainage", q_bottom),
                ("runoff", np.maximum(q_top, 0.0) - np.maximum(q_surf, 0.0)),
                ("evaporation_deficit", np.maximum(-q_top, 0.0) - np.maximum(-q_surf, 0.0)),
            ):
                step_totals[key] += float(np.mean(flux)) * sub
            t += sub
        for key, value in step_totals.items():
            self.totals[key] += value
        self.last_step = dict(step_totals, substeps=substeps)
        return self.last_step

    @property
    def storage(self):
        """Water stored per unit surface area, cm."""
        return float(self.theta.sum() * self.dz / (self.nx * self.ny))

    @property
    def mass_balance_error(self):
        """Storage change minus net boundary inflow since the start, cm (should be ~0)."""
        t = self.totals
        net_in = t["infiltration"] - t["evaporation"] - t["drainage"]
        return self.storage - self.initial_storage - net_in

    @property
    def theta_mean(self):
        return float(np.mean(self.theta))

    @property
    def pressure_head_mean(self):
        return float(np.mean(self.h))
