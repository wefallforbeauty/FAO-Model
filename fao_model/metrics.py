"""Agreement statistics between an estimate and a reference series."""
import numpy as np


def compare(sim, ref):
    """Bias, relative bias, MAE, RMSE, Pearson r, NSE and the fitted line sim = slope * ref + intercept.

    bias, MAE and RMSE are in the units of the inputs, rel_bias in %.
    """
    sim = np.asarray(sim, dtype=float)
    ref = np.asarray(ref, dtype=float)
    diff = sim - ref
    slope, intercept = np.polyfit(ref, sim, 1)
    return {
        "n": int(len(ref)),
        "mean_ref": float(ref.mean()),
        "mean_sim": float(sim.mean()),
        "bias": float(diff.mean()),
        "rel_bias": float(100.0 * diff.mean() / ref.mean()),
        "mae": float(np.abs(diff).mean()),
        "rmse": float(np.sqrt(np.mean(diff ** 2))),
        "r": float(np.corrcoef(sim, ref)[0, 1]),
        "nse": float(1.0 - np.sum(diff ** 2) / np.sum((ref - ref.mean()) ** 2)),
        "slope": float(slope),
        "intercept": float(intercept),
    }
