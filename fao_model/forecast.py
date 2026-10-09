import math

import numpy as np

try:
    from timesfm import TimesFm, TimesFmHparams
    TIMESFM_AVAILABLE = True
except Exception:
    TIMESFM_AVAILABLE = False


class Forecaster:
    CONTEXT_LEN = 2048
    HORIZON_LEN = 96

    def __init__(self):
        self.model = None
        self.timesfm = False
        self.status = "TimesFM unavailable; using dynamic fallback."
        if TIMESFM_AVAILABLE:
            try:
                hparams = TimesFmHparams(
                    context_len=self.CONTEXT_LEN,
                    horizon_len=self.HORIZON_LEN,
                    input_patch_len=32,
                    output_patch_len=128,
                    num_layers=20,
                    model_dims=1280,
                )
                self.model = TimesFm(hparams)
                self.model.load_from_checkpoint(repo_id="google/timesfm-1.0-200m")
                self.timesfm = True
                self.status = f"TimesFM active: context {self.CONTEXT_LEN}, horizon {self.HORIZON_LEN}"
            except Exception as exc:
                self.model = None
                self.timesfm = False
                self.status = f"TimesFM init failed ({type(exc).__name__}); using fallback."

    def _fallback(self, series, steps):
        if not series:
            return [0.0] * steps
        if len(series) == 1:
            return [float(series[-1])] * steps
        n = min(len(series), self.CONTEXT_LEN)
        y = np.asarray(series[-n:], dtype=float)
        x = np.arange(n, dtype=float)
        w = max(5, min(96, n // 8))
        recent = float(np.mean(y[-w:]))
        older = float(np.mean(y[max(0, n - 2*w):-w])) if n > w else recent
        slope = (recent - older) / max(1.0, w)
        level = float(np.mean(y))
        out = []
        last = float(y[-1])
        for k in range(1, steps + 1):
            trend = last + slope * k
            decay = math.exp(-k / max(12.0, steps / 2.0))
            value = level + (trend - level) * decay
            out.append(float(value))
        return out

    def forecast(self, series, steps=None):
        steps = int(steps or self.HORIZON_LEN)
        steps = max(1, min(steps, self.HORIZON_LEN))
        if self.timesfm and len(series) >= 32:
            try:
                context = np.asarray(series[-self.CONTEXT_LEN:], dtype=np.float32)
                fc = self.model.forecast(context, horizon_len=steps)
                # Different releases return arrays, tuples, or (forecast, ...).
                if isinstance(fc, tuple):
                    fc = fc[0]
                fc = np.asarray(fc, dtype=float).reshape(-1)
                if len(fc) >= steps and np.all(np.isfinite(fc[:steps])):
                    return fc[:steps].tolist()
            except Exception:
                pass
        return self._fallback(series, steps)

    def forecast_all(self, history, steps=None):
        steps = steps or self.HORIZON_LEN
        return {key: self.forecast(series, steps) for key, series in history.items()}
