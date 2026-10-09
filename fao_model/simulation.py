from collections import deque

from .calibration import Calibrator
from .et import compute_all
from .forecast import Forecaster
from .richards import RichardsSolver3D
from .weather import WeatherGenerator


class Simulation:
    """One step is one day: weather -> ET equations -> soil water -> forecasts."""

    CALIBRATION_WINDOW = 365    # days

    def __init__(self, seed=None):
        self.gen = WeatherGenerator(lat=40.0, elev=100.0, seed=seed)
        self.heat_index = self.gen.heat_index()
        self.recent_T = deque(maxlen=30)
        self.calib = Calibrator()
        self.richards = RichardsSolver3D(nx=5, ny=5, nz=5)
        self.forecaster = Forecaster()
        self.data_history = []
        self.data_batch = deque(maxlen=self.CALIBRATION_WINDOW)
        self.history = {}
        self.forecasts = {}
        self.step_number = 0

    def step(self):
        data = self.gen.ensure_required(self.gen.generate())
        self.step_number += 1
        self.recent_T.append(data["T_mean"])
        data["T_30d"] = sum(self.recent_T) / len(self.recent_T)
        data["heat_index"] = self.heat_index
        self.data_history.append(data.copy())
        # Walk-forward: today's values use coefficients fitted on earlier days only.
        self.calib.calibrate(list(self.data_batch))
        results = compute_all(data, self.calib.coeffs)
        self.data_batch.append(data.copy())
        q_top_cm_s = (
            data["precip"] - results["ET_PM_maize"]
        ) * 0.1 / 86400.0
        for _ in range(10):
            self.richards.step(q_top_cm_s, dt=0.1)
        results["Richards_theta_mean"] = self.richards.theta_mean
        results["Richards_head_mean"] = self.richards.pressure_head_mean
        for key, val in results.items():
            self.history.setdefault(key, []).append(float(val))
        self.forecasts = self.forecaster.forecast_all(
            self.history, steps=self.forecaster.HORIZON_LEN
        )
        return data, results, self.forecasts
