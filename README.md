# FAO-Model

Daily evapotranspiration (ET) equations, a 3D Richards soil-water solver and a
forecasting dashboard, driven by real weather data.

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

## Validation summary

Konya Plain, calibrated on 1995–2014 and tested on 2015–2024 against this repo's
daily ASCE-PM. Full tables and figures: [reports/validation](reports/validation/README.md).

| Method | Relative bias, literature coefficients | Test RMSE (mm/day), literature → calibrated |
|---|---|---|
| Turc | −2.7 % | 0.67 → 0.66 |
| Jensen-Haise | +4.5 % | 0.87 → 0.68 |
| Abtew | +6.0 % | 0.98 → 0.98 |
| Hargreaves | −10.4 % | 0.78 → 0.66 |
| Blaney-Criddle (FAO-24) | +17.0 % | 1.17 → 0.63 |
| Priestley-Taylor | −19.5 % | 0.97 → 0.65 |
| Thornthwaite (calendar months) | −38.6 % | not calibrated |

- ASCE-PM matches pyfao56 to within 0.005 mm/day and follows Open-Meteo's
  hourly-based ET0 with r = 0.995 (mean difference −4.8 %).
- The earlier convergence of every equation to ASCE-PM came from the method:
  each step refitted all coefficients to ASCE-PM on the same samples. That also
  hid unit errors (calibrated Turc a = 4.25, Jensen-Haise a = 13.9). With correct
  units every equation except Thornthwaite is within ±20 % without calibration.
- Calibration does not remove the bias in the test decade (−7 % to +2 %), and
  the coefficients take on the reference's bias: fitted to Open-Meteo ET0
  instead of ASCE-PM, the multiplicative coefficients shift by 3–4 %
  (Priestley-Taylor α 1.49 vs 1.54).
- Thornthwaite underestimates by about 39 % in this semi-arid climate, in line
  with its known underestimation in arid climates. It no longer sits at the top
  of the graph: it used to return mm/month.

## Data

`data/konya_cumra_era5_daily_1995_2024.csv`: daily ERA5 weather for the Konya
Plain near Çumra (37.57° N, 32.78° E, 1013 m), 10 958 days without gaps. Weather
data by [Open-Meteo.com](https://open-meteo.com/) (CC BY 4.0), generated using
Copernicus Climate Change Service information. Columns, caveats and citations:
[data/README.md](data/README.md).

## Roadmap

Original notes:

- The program should be written for real weather data or real weather data history. Then the timesfm forecast should be applicable weather data.
- Richards equation part should be much more complex:
  - A root uptake function should be added to program.
  - The h parameter for Richards eq. should be adjusted to real, plausible soil statistics.
  - Different numerical solutions should be investigated other than FVM like: FDM, FEM, Mixed Finite Element, Explicit Euler, Implicit Euler, Newton-Raphson, Picard.
  - Adaptive time stamp could be added.
- Some uncertainty could be added to ET and \theta(x,y,z,t) functions.
- Both the forecast data and the solutions should be an input for another optimization process.
- The soil and plant data could be researched to add more realistic parameters and constants.
- Currently the ASCE-PM is taken reference for the calibration for various equations. This seems like it forces them to converge to ASCE-PM, this could be further investigated.
- A validation code should be added.
- Thornthwaite eq. should be adjusted for monthly mean temp values not for daily. This might be the reason why thornthwaite hangs at top of the graph.
- I should also investigate if the convergence of the different equations is because of any bias. Similarly relative bias could be calculated separately which a scatter plot could also help.
- The parameters and constants should be configured in such a way that they would show which equations is most susceptible to which parameter/constant and how.
- All equations, variables and numerical fluxes should be checked for dimensional consistency and unit conversions.
- Add automated water-mass conservation tests for the Richards solver.

Progress (2026-10-09):

| Item | Status |
|---|---|
| Real weather data | Done: the Konya 1995–2024 record drives the dashboard and the validation. TimesFM still forecasts the model outputs, not the weather inputs, and its loading code has not been checked against a current TimesFM release (it falls back silently). |
| Root uptake | Open |
| Plausible soil parameters | Open; the parameters are unchanged (close to the Carsel & Parrish loam) |
| Other numerical schemes | Open; the explicit FVM is now conservative and tested, a baseline to compare against |
| Adaptive time step | Done for the explicit solver (stability-based sub-steps) |
| Uncertainty in ET and θ | Open |
| Forecasts and solutions feeding an optimization | Open |
| Soil and plant data | Open |
| Calibration converging to ASCE-PM | Investigated, see the validation summary |
| Validation code | Done: FAO-56 examples, pyfao56 cross-check, Open-Meteo ET0, training/test report |
| Thornthwaite on monthly means | Done: monthly method with the annual heat index; the old values were mm/month |
| Bias, relative bias, scatter plots | Done: reports/validation, sections 2–3 |
| Parameter sensitivity | Open |
| Dimensional consistency | Done for the ET equations and the Richards fluxes, with unit tests |
| Water-mass conservation tests | Done: `tests/test_richards.py` |
