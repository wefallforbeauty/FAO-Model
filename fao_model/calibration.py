import numpy as np

from .meteo import LAMBDA, day_length
from .et import (
    asce_pm, hargreaves_base, turc_base, priestley_taylor_base, jensen_haise_base,
)


class Calibrator:
    def __init__(self):
        self.coeffs = {
            "hargreaves_a": 1.0,
            "turc_a": 1.0,
            "abtew_a": 0.53,
            "pt_alpha": 1.26,
            "jh_a": 1.0,
            "bc_a": 0.0,
            "bc_b": 1.0,
            "kc_maize": 0.35,
        }

    def calibrate(self, data_batch):
        if len(data_batch) < 3:
            return
        rows = []
        for d in data_batch:
            try:
                eto = asce_pm(d)
                rows.append((
                    eto,
                    hargreaves_base(d),
                    turc_base(d),
                    d["Rs"] / LAMBDA,
                    priestley_taylor_base(d),
                    jensen_haise_base(d),
                    (day_length(d["lat"], d["doy"]) / 24.0)
                    * (0.46 * d["T_mean"] + 8.13),
                ))
            except Exception:
                continue
        if len(rows) < 3:
            return
        a = np.asarray(rows, dtype=float)
        eto = a[:, 0]
        for key, col in (
            ("hargreaves_a", 1),
            ("turc_a", 2),
            ("abtew_a", 3),
            ("pt_alpha", 4),
            ("jh_a", 5),
        ):
            base = a[:, col]
            denom = float(np.dot(base, base))
            if denom > 1e-12:
                # Least-squares coefficient is more stable than sum(y)/sum(x).
                self.coeffs[key] = float(np.dot(base, eto) / denom)
        f = a[:, 6]
        A = np.vstack([f, np.ones(len(f))]).T
        sol, _, _, _ = np.linalg.lstsq(A, eto, rcond=None)
        self.coeffs["bc_b"] = float(sol[0])
        self.coeffs["bc_a"] = float(sol[1])
