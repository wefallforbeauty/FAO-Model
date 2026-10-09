from .calibration import Calibrator
from .et import compute_all
from .forecast import Forecaster
from .richards import RichardsSolver3D
from .weather import WeatherGenerator


class Simulation:
    def __init__(self):
        self.gen = WeatherGenerator(lat=40.0, elev=100.0, interval_minutes=10)
        self.calib = Calibrator()
        self.richards = RichardsSolver3D(nx=5, ny=5, nz=5)
        self.forecaster = Forecaster()
        self.data_history = []
        self.data_batch = []
        self.history = {}
        self.forecasts = {}
        self.step_number = 0

    def step(self):
        data = self.gen.ensure_required(self.gen.generate())
        self.step_number += 1
        self.data_history.append(data.copy())
        self.data_batch.append(data.copy())
        if len(self.data_batch) > 256:
            self.data_batch.pop(0)
        self.calib.calibrate(self.data_batch)
        results = compute_all(data, self.calib.coeffs)
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
