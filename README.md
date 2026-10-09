# FAO-Model 02

Daily evapotranspiration (ET) equations, a 3D Richards soil-water solver and a
forecasting dashboard, driven by real weather data. What changed in this
version: [UPDATE_LOG.md](UPDATE_LOG.md).

- **ET equations** (mm/day): FAO-56 / ASCE Penman-Monteith reference ET, maize
  ET (Kc × ETo), Thornthwaite, Blaney-Criddle (FAO-24), Turc, Priestley-Taylor,
  Hargreaves, Jensen-Haise and Abtew, checked against FAO-56 worked examples.
- **Calibration** of the empirical equations against a reference ET:
  walk-forward (earlier days only) in the dashboard, training/test split in the
  validation.
- **Soil water**: 3D van Genuchten–Mualem Richards solver (mass-conservative
  explicit finite volumes with adaptive sub-steps) driven by rain minus maize ET.
- **Forecasts** of every series: TimesFM when installed, otherwise a
  damped-trend fallback.
- **Data**: 30 years (1995–2024) of daily ERA5 weather for the Konya Plain,
  Türkiye, from Open-Meteo.

## Quick start

```bash
python3 -m venv .venv
.venv/bin/pip install -e ".[dev]"                 # numpy, matplotlib, pytest
.venv/bin/python -m fao_model                     # dashboard, replays the Konya record
.venv/bin/python -m fao_model --synthetic         # synthetic daily weather instead
.venv/bin/python -m fao_model --interval-ms 200   # faster replay (default 1000 ms per day)
.venv/bin/python -m fao_model.validate            # regenerates reports/validation/
.venv/bin/python -m fao_model.data --help         # download another site or period
.venv/bin/python -m pytest
```

The live graph needs matplotlib with Tk support. If the dashboard reports that
`ImageTk` is missing (system Pillow on Fedora), install `python3-pillow-tk`, or
run `pip install --ignore-installed pillow` inside the venv.

## Layout

| Module | Contents |
|---|---|
| `fao_model/meteo.py` | FAO-56 helpers: vapour pressure, radiation, day length, pressure, wind height |
| `fao_model/et.py` | the ET equations and `compute_all` |
| `fao_model/calibration.py` | literature coefficients and least-squares calibration |
| `fao_model/weather.py` | synthetic daily weather generator |
| `fao_model/data.py` | Open-Meteo download and CSV loader |
| `fao_model/richards.py` | van Genuchten–Mualem functions and the 3D Richards solver |
| `fao_model/forecast.py` | TimesFM forecaster with fallback |
| `fao_model/simulation.py` | the daily model loop, without GUI |
| `fao_model/app.py` | Tkinter dashboard |
| `fao_model/metrics.py`, `fao_model/validate.py` | agreement statistics and the validation report |
| `tests/` | FAO-56 examples, cross-check against pyfao56, solver and data tests |

## Data

`data/konya_cumra_era5_daily_1995_2024.csv`: daily ERA5 weather for the Konya
Plain near Çumra (37.57° N, 32.78° E, 1013 m), 10 958 days without gaps. Weather
data by [Open-Meteo.com](https://open-meteo.com/) (CC BY 4.0), generated using
Copernicus Climate Change Service information. Columns, caveats and citations:
[data/README.md](data/README.md).

## Validation

Calibrated on 1995–2014 and tested on 2015–2024 against this repo's daily
ASCE-PM; full tables and figures in [reports/validation](reports/validation/README.md).

| Method | Relative bias, literature coefficients | Test RMSE (mm/day), literature → calibrated |
|---|---|---|
| Turc | −2.7 % | 0.67 → 0.66 |
| Jensen-Haise | +4.5 % | 0.87 → 0.68 |
| Abtew | +6.0 % | 0.98 → 0.98 |
| Hargreaves | −10.4 % | 0.78 → 0.66 |
| Blaney-Criddle (FAO-24) | +17.0 % | 1.17 → 0.63 |
| Priestley-Taylor | −19.5 % | 0.97 → 0.65 |
| Thornthwaite (calendar months) | −38.6 % | not calibrated |

ASCE-PM matches pyfao56 to within 0.005 mm/day and follows Open-Meteo's
hourly-based ET0 with r = 0.995 (−4.8 %).

## Not done in this update

### From the 01 notes (original wording)

- Then the timesfm forecast should be applicable weather data.
- Richards equation part should be much more complex:
  - A root uptake function should be added to program.
  - The h parameter for Richards eq. should be adjusted to real, plausible soil statistics.
  - Different numerical solutions should be investigated other than FVM like: FDM, FEM, Mixed Finite Element, Explicit Euler, Implicit Euler, Newton-Raphson, Picard.
- Some uncertainty could be added to ET and \theta(x,y,z,t) functions.
- Both the forecast data and the solutions should be an input for another optimization process.
- The soil and plant data could be researched to add more realistic parameters and constants.
- The parameters and constants should be configured in such a way that they would show which equations is most susceptible to which parameter/constant and how.

### Found while making 02

- **TimesFM is not installed** here, so every forecast in 02 comes from the
  damped-trend fallback. Also:
  - The loading code (`TimesFm(hparams)` with `load_from_checkpoint`) has not
    been checked against a current TimesFM release.
  - Forecast errors are swallowed silently.
  - The context length (2048) is larger than timesfm-1.0-200m supports (512).
  - It forecasts the model outputs, not the weather inputs.
- **Crop ET** is drawn only from the top 5 cm layer (no roots), and the maize Kc
  is a constant 0.35 all year. Over 1995–2024 the soil supplied 3577 mm of
  evaporation; another 8722 mm of net demand (ET minus same-day rain) went unmet.
- **Rain and crop ET** are combined into one net daily surface flux. There is no
  irrigation.
- **Near saturation** the explicit solver needs very small sub-steps (heavy rain
  on fine soils, ponding), and it stops with an error if a layer would exceed
  saturation. An implicit scheme removes that limit.
- **Validation** is against computed references (ASCE-PM, Open-Meteo ET0), not
  measured ET, and covers one site. ERA5 is reanalysis, not station data.
- **Blaney-Criddle:** n/N is estimated from Rs/Ra, and the FAO-24 regression is
  used with 24-hour wind instead of daytime wind.
- **Calibration** fits one constant per equation; seasonal or regional
  coefficients are not explored.

## Next steps

1. **Root uptake:** root water uptake (e.g. Feddes) and an FAO-56 maize Kc curve
   with a planting date, so ET comes from the root zone and water stress becomes
   visible. Then irrigation scenarios.
2. **TimesFM:**
   - install the current release and fix the loading code;
   - forecast the weather inputs (temperature, humidity, wind, radiation, rain)
     and compute ET from them;
   - measure forecast skill on 2015–2024.
3. **Implicit Richards solver:** mixed-form modified Picard, then Newton,
   compared with the explicit one. Soil parameters from the Carsel & Parrish
   (1988) texture classes, or SoilGrids with Rosetta for Konya.
4. **Sensitivity analysis:** Morris screening and Sobol indices of each ET
   equation to its inputs and coefficients.
5. **Soil moisture check:** compare the simulated water content with the ERA5
   0–7 cm soil moisture already in the CSV.
6. **Uncertainty:** Monte Carlo on inputs and parameters, and TimesFM quantiles.
7. **Irrigation optimization** driven by the forecasts and the soil water state.
8. **Measured ET:** validation against lysimeter or eddy-covariance data, and
   more sites.
