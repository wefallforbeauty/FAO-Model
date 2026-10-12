from collections import deque

from .calibration import Calibrator
from .crop import Crop, adjust_for_climate
from .data import monthly_means
from .et import compute_all, thornthwaite_heat_index
from .field import COLUMN, FieldModel
from .forecast import Forecaster
from .richards import RichardsSolver3D
from .weather import WeatherGenerator


class Simulation:
    """One step is one day: weather -> ET equations -> maize field water balance -> forecasts.

    records: daily weather to replay in order (e.g. fao_model.data.load()[1]);
    None runs the synthetic weather generator instead. step() raises
    StopIteration when the records run out. heat_index (Thornthwaite)
    defaults to the climatology of the records, which must span a year.
    irrigation: True for automatic irrigation of the maize field, False
    for a rainfed field. crop defaults to maize with Kcb_mid adjusted to
    the records' mid-season climate.
    """

    CALIBRATION_WINDOW = 365    # days

    def __init__(self, records=None, seed=None, heat_index=None, irrigation=True, crop=None):
        if records is None:
            self.gen = WeatherGenerator(lat=40.0, elev=100.0, seed=seed)
            self.source = iter(self.gen.generate, None)
            self.heat_index = heat_index or self.gen.heat_index()
            crop = crop or Crop()
        else:
            self.gen = WeatherGenerator(lat=records[0]["lat"], elev=records[0]["elev"])
            self.source = iter(records)
            # Thornthwaite's heat index is a site constant: use the record's climatology.
            self.heat_index = heat_index or thornthwaite_heat_index(monthly_means(records))
            crop = crop or adjust_for_climate(Crop(), records)
        self.recent_T = deque(maxlen=30)
        self.calib = Calibrator()
        soil = RichardsSolver3D(**dict(COLUMN, nx=5, ny=5))
        self.field = FieldModel(crop, irrigation=irrigation, soil=soil)
        self.richards = soil
        self.forecaster = Forecaster()
        self.data_history = []
        self.data_batch = deque(maxlen=self.CALIBRATION_WINDOW)
        self.history = {}
        self.forecasts = {}
        self.step_number = 0

    def step(self):
        data = self.gen.ensure_required(dict(next(self.source)))
        self.step_number += 1
        self.recent_T.append(data["T_mean"])
        data["T_30d"] = sum(self.recent_T) / len(self.recent_T)
        data["heat_index"] = self.heat_index
        self.data_history.append(data.copy())
        # Walk-forward: today's values use coefficients fitted on earlier days only.
        self.calib.calibrate(list(self.data_batch))
        results = compute_all(data, self.calib.coeffs)
        self.data_batch.append(data.copy())
        field = self.field.step(data, results["ET_ASCE_PM"])
        results["ET_crop_potential"] = field["tp"] + field["ep"]
        results["ET_crop_actual"] = field["ta"] + field["ea"]
        results["Richards_theta_mean"] = self.richards.theta_mean
        results["Richards_head_mean"] = self.richards.pressure_head_mean
        for key, val in results.items():
            self.history.setdefault(key, []).append(float(val))
        self.forecasts = self.forecaster.forecast_all(
            self.history, steps=self.forecaster.HORIZON_LEN
        )
        return data, results, self.forecasts
