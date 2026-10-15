from collections import deque
from datetime import datetime

from .calibration import Calibrator
from .crop import Crop, adjust_for_climate
from .data import monthly_means
from .et import compute_all, thornthwaite_heat_index
from .field import COLUMN, FieldModel
from .forecast import Climatology, WeatherForecaster, add_eto
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

    forecaster: "auto" (TimesFM if it loads), "timesfm", "baseline", or
    "ecmwf" / "gfs" (archived weather-model forecasts, 1-7 days). A forecast
    is issued every `forecast_every` days (default 7 with TimesFM, which
    takes about a second, else daily). The baseline's climatology defaults
    to that of the records (a year at least). start: first day to replay;
    the site constants still come from all the records.
    """

    CALIBRATION_WINDOW = 365    # days

    def __init__(self, records=None, seed=None, heat_index=None, irrigation=True, crop=None,
                 forecaster="auto", forecast_every=None, climatology=None, start=None):
        if records is None:
            self.gen = WeatherGenerator(lat=40.0, elev=100.0, seed=seed)
            self.source = iter(self.gen.generate, None)
            self.heat_index = heat_index or self.gen.heat_index()
            crop = crop or Crop()
            sample = WeatherGenerator(lat=40.0, elev=100.0, start=datetime(2001, 1, 1),
                                      seed=None if seed is None else seed + 1)
            climate = [sample.ensure_required(sample.generate()) for _ in range(5 * 365)]
        else:
            self.gen = WeatherGenerator(lat=records[0]["lat"], elev=records[0]["elev"])
            first = 0 if start is None else next(
                i for i, d in enumerate(records) if d["timestamp"].date() >= start)
            self.source = iter(records[first:])
            # Thornthwaite's heat index is a site constant: use the record's climatology.
            self.heat_index = heat_index or thornthwaite_heat_index(monthly_means(records))
            crop = crop or adjust_for_climate(Crop(), records)
            climate = records
        self.recent_T = deque(maxlen=30)
        self.calib = Calibrator()
        soil = RichardsSolver3D(**dict(COLUMN, nx=5, ny=5))
        self.field = FieldModel(crop, irrigation=irrigation, soil=soil)
        self.richards = soil
        # Like the heat index, the climatology behind the baseline forecast is a site constant.
        climatology = climatology or Climatology(add_eto(climate))
        self.forecaster = WeatherForecaster(climatology, backend=forecaster)
        self.forecast_every = forecast_every or (7 if self.forecaster.model is not None else 1)
        self.forecast = None
        self.field_history = []
        self.data_history = []
        self.data_batch = deque(maxlen=self.CALIBRATION_WINDOW)
        self.history = {}
        self.step_number = 0

    def step(self):
        data = self.gen.ensure_required(dict(next(self.source)))
        self.step_number += 1
        self.recent_T.append(data["T_mean"])
        data["T_30d"] = sum(self.recent_T) / len(self.recent_T)
        data["heat_index"] = self.heat_index
        # Walk-forward: today's values use coefficients fitted on earlier days only.
        self.calib.calibrate(list(self.data_batch))
        results = compute_all(data, self.calib.coeffs)
        self.data_batch.append(data.copy())
        data["ETo"] = results["ET_ASCE_PM"]
        self.data_history.append(data.copy())
        field = self.field.step(data, results["ET_ASCE_PM"])
        self.field_history.append(field)
        results["ET_crop_potential"] = field["tp"] + field["ep"]
        results["ET_crop_actual"] = field["ta"] + field["ea"]
        results["Richards_theta_mean"] = self.richards.theta_mean
        results["Richards_head_mean"] = self.richards.pressure_head_mean
        for key, val in results.items():
            self.history.setdefault(key, []).append(float(val))
        if self.forecast is None or self.step_number - self.forecast["step"] >= self.forecast_every:
            issued = self.forecaster.forecast_at(self.data_history, [len(self.data_history) - 1])
            self.forecast = {k: v[0] for k, v in issued.items()}
            self.forecast.update(step=self.step_number, issued=data["timestamp"])
        return data, results, self.forecast
