import numpy as np

from .et import (
    asce_pm, hargreaves_base, turc_base, abtew_base, priestley_taylor_base,
    jensen_haise_base, blaney_criddle_f,
)

# Coefficients of the published equations, i.e. no calibration.
LITERATURE_COEFFS = {
    "hargreaves_a": 1.0,
    "turc_a": 1.0,
    "abtew_a": 0.53,
    "pt_alpha": 1.26,
    "jh_a": 1.0,
    "bc_a": None,       # None: FAO-24 a and b from RH_min, n/N and u2
    "bc_b": None,
    "kc_maize": 0.35,
}

# Multiplicative coefficient -> equation it scales
SCALED = (
    ("hargreaves_a", hargreaves_base),
    ("turc_a", turc_base),
    ("abtew_a", abtew_base),
    ("pt_alpha", priestley_taylor_base),
    ("jh_a", jensen_haise_base),
)


class Calibrator:
    MIN_ROWS = 30

    def __init__(self):
        self.coeffs = dict(LITERATURE_COEFFS)

    def calibrate(self, data_batch, reference=None):
        """Fit the coefficients so each equation matches a reference ET (mm/day).

        reference defaults to ASCE-PM of each record. Note that this makes
        the calibrated equations agree with the reference by construction;
        judge them on data they were not fitted to (see validate.py).
        """
        if len(data_batch) < self.MIN_ROWS:
            return
        if reference is None:
            reference = [asce_pm(d) for d in data_batch]
        eto = np.asarray(reference, dtype=float)
        for key, base_fn in SCALED:
            base = np.array([base_fn(d) for d in data_batch], dtype=float)
            denom = float(np.dot(base, base))
            if denom > 1e-12:
                # Least-squares coefficient is more stable than sum(y)/sum(x).
                self.coeffs[key] = float(np.dot(base, eto) / denom)
        f = np.array([blaney_criddle_f(d) for d in data_batch], dtype=float)
        A = np.vstack([f, np.ones(len(f))]).T
        sol, _, _, _ = np.linalg.lstsq(A, eto, rcond=None)
        self.coeffs["bc_b"] = float(sol[0])
        self.coeffs["bc_a"] = float(sol[1])
