# FAO-Model 05

Daily evapotranspiration (ET) equations and the water balance of an irrigated
maize field on a 3D Richards soil-water solver, with a forecasting dashboard,
driven by real weather data. What changed in this version:
[UPDATE_LOG.md](UPDATE_LOG.md).

- **ET equations** (mm/day): FAO-56 / ASCE Penman-Monteith reference ET,
  Thornthwaite, Blaney-Criddle (FAO-24), Turc, Priestley-Taylor, Hargreaves,
  Jensen-Haise and Abtew, checked against FAO-56 worked examples.
- **Calibration** of the empirical equations against a reference ET:
  walk-forward (earlier days only) in the dashboard, training/test split in the
  validation.
- **Crop**: grain maize with FAO-56 dual crop coefficients. Potential
  transpiration goes to growing roots through the Feddes stress function,
  potential soil evaporation to the surface.
- **Irrigation**: optional automatic irrigation that refills the root zone
  to field capacity when 55 % of the available water is used; it can wait
  when rain is forecast.
- **Soil water**: 3D van Genuchten–Mualem Richards solver (mass-conservative
  explicit finite volumes with adaptive sub-steps), a 1.5 m column in 2.5 cm
  layers.
- **Forecasts**: forecasts of the weather inputs and reference ET, scored
  against what happened:
  - TimesFM 2.5 (optional) or a statistical baseline, 14 days ahead;
  - archived ECMWF and GFS weather-model forecasts, 1–7 days ahead, with a
    running bias correction.
- **Data**: daily ERA5 weather for the Konya Plain, Türkiye, from Open-Meteo:
  1995–2024 for the reports, extended to September 2026 for the weather-model
  comparison.

## Quick start

Requires Python 3.10 or newer with Tk support. The python.org installers for
Windows and macOS include Tk; some Linux distributions ship it as a separate
package (for example `python3-tk`).

Create a virtual environment and install the package with its development tools:

```sh
python -m venv .venv
# Activate it
#   Windows (PowerShell):  .venv\Scripts\Activate.ps1
#   macOS and Linux:       source .venv/bin/activate
pip install -e ".[dev]"                 # numpy, matplotlib, pytest, pandas, pyfao56
```

TimesFM forecasts are optional. They need PyTorch and download about 1 GB of
weights from Hugging Face on first use; on a CPU a forecast takes about a
second and the model about 2 GB of memory. Without TimesFM everything runs with
the statistical baseline.

```sh
# Linux: the default PyTorch wheel includes GPU libraries (several GB);
# for a CPU-only install, get PyTorch from its CPU index first:
#   pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install -e ".[forecast]"            # timesfm, torch
```

With the environment active:

```sh
python -m fao_model                     # dashboard, replays the Konya record
python -m fao_model --synthetic         # synthetic daily weather instead
python -m fao_model --interval-ms 200   # faster replay (default 1000 ms per day)
python -m fao_model --irrigation none   # rainfed maize instead of automatic irrigation
python -m fao_model --forecaster baseline   # statistical forecasts instead of TimesFM
python -m fao_model --forecaster ecmwf --start 2024-03-15   # archived ECMWF forecasts (or gfs)
python -m fao_model.validate            # regenerates reports/validation/
python -m fao_model.water_balance       # regenerates reports/water_balance/ (about 8 minutes)
python -m fao_model.forecast_skill      # regenerates reports/forecast_skill/ (TimesFM: about 20 minutes)
python -m fao_model.nwp_skill           # regenerates reports/nwp/ (TimesFM: about 20 minutes)
python -m fao_model.nwp --help          # download archived weather-model forecasts
python -m fao_model.data --help         # download another site or period
pytest
```

On systems where `python` points to Python 2 or is missing, use `python3`.
If the dashboard cannot draw the live graph, it shows the import error that
names the missing piece, usually Tk support in Python or in Pillow. If the
TimesFM download fails, setting `HF_HUB_DISABLE_XET=1` makes Hugging Face fall
back to plain HTTP.

## Layout

| Module | Contents |
|---|---|
| `fao_model/meteo.py` | FAO-56 helpers: vapour pressure, radiation, day length, pressure, wind height |
| `fao_model/et.py` | the ET equations and `compute_all` |
| `fao_model/calibration.py` | literature coefficients and least-squares calibration |
| `fao_model/weather.py` | synthetic daily weather generator |
| `fao_model/data.py` | Open-Meteo download and CSV loader |
| `fao_model/richards.py` | van Genuchten–Mualem functions and the 3D Richards solver with root uptake |
| `fao_model/crop.py` | maize parameters, FAO-56 dual crop coefficients, root growth, Feddes stress |
| `fao_model/field.py` | daily field water balance: crop, soil and automatic irrigation |
| `fao_model/forecast.py` | weather and reference ET forecasts: TimesFM 2.5 or climatology with damped anomaly persistence |
| `fao_model/forecast_skill.py` | forecast skill report |
| `fao_model/nwp.py` | download and load archived ECMWF and GFS forecasts by lead day |
| `fao_model/nwp_skill.py` | weather-model forecasts against TimesFM and the baselines |
| `fao_model/simulation.py` | the daily model loop, without GUI |
| `fao_model/app.py` | Tkinter dashboard |
| `fao_model/metrics.py`, `fao_model/validate.py` | agreement statistics and the ET validation report |
| `fao_model/water_balance.py` | maize water balance report: scenarios, pyfao56 and ERA5 comparisons |
| `fao_model/plotting.py` | shared figure style of the reports |
| `tests/` | FAO-56 examples, cross-checks against pyfao56, solver, crop, field and data tests |

## Data

- `data/konya_cumra_era5_daily_1995_2024.csv`: daily ERA5 weather for the
  Konya Plain near Çumra (37.57° N, 32.78° E, 1013 m), 10 958 days without gaps.
- `data/konya_cumra_era5_daily_2025_2026.csv`: the same for 2025-01-01 to
  2026-09-30, kept apart so the 1995–2024 reports stay reproducible.
- `data/nwp/`: daily weather forecast by ECMWF IFS 0.25° and NCEP GFS 1–7 days
  ahead, 2024-03-15 to 2026-09-30, from Open-Meteo's Previous Runs API.

Weather data by [Open-Meteo.com](https://open-meteo.com/) (CC BY 4.0), generated
using Copernicus Climate Change Service information; NCEP GFS (public domain)
and ECMWF IFS open data (CC BY 4.0). Columns, caveats and citations:
[data/README.md](data/README.md).

## Results

### ET equations

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

### Maize water balance

Grain maize planted on 1 May, 30 seasons (1995–2024); full tables and figures
in [reports/water_balance](reports/water_balance/README.md). Season means, mm:

| | Rainfed | Irrigated |
|---|---|---|
| Rain | 51 | 51 |
| Irrigation | 0 | 821 (14.5 times) |
| Actual / potential transpiration | 56 / 722 | 667 / 722 |
| Actual / potential soil evaporation | 42 / 398 | 205 / 398 |
| Drainage below 1.5 m | 1 | 7 |

- **Rainfed:** maize at this site gets only about 8 % of the transpiration it
  needs, so it depends on irrigation.
- **Against pyfao56:** run with the same inputs and start profiles, pyfao56
  gives season means within 5 %: transpiration 696 mm, evaporation 201 mm and
  irrigation 854 mm in the irrigated scenario.
- **Against ERA5:** the simulated topsoil water content follows ERA5's 0–7 cm
  soil moisture (daily r = 0.80, monthly cycle r = 0.97) but is 0.07 m³/m³
  drier; the two soils differ.

### Forecasts

Forecasts issued weekly over 2015–2024 (520), each from data up to its issue
day; full tables and figures in [reports/forecast_skill](reports/forecast_skill/README.md).
Reference ET error (mm/day) by days ahead:

| | 1 day | 3 days | 7 days | 7-day total (mm) |
|---|---|---|---|---|
| Climatology | 0.58 | 0.55 | 0.58 | 2.9 |
| TimesFM 2.5, ETo from forecast weather | 0.40 | 0.52 | 0.58 | 2.5 |
| Damped anomaly persistence | 0.40 | 0.51 | 0.57 | 2.5 |
| Persistence (tomorrow = today) | 0.45 | 0.65 | 0.79 | 3.8 |

- **Zero-shot TimesFM matches the statistical baseline.** It is slightly better
  for minimum temperature and humidity a week ahead. Beyond about three days
  every method is back to climatology.
- **Rain:** no method caught any of the 32 cases with at least 10 mm of rain
  within three days.
- **Uncertainty:** TimesFM's 10–90 % interval is well calibrated (78–84 % of
  observations inside).
- **Why the ceiling:** a single site's history cannot see weather systems
  coming; that takes numerical weather prediction.
- **Irrigation with rain forecasts:** postponing irrigation when rain is
  forecast saves no water at this site, because in-season rain is too rare.
  Over 2015–2024 only 3 irrigations were followed by at least 10 mm of rain
  within three days. A perfect forecast avoids those 3 but uses as much water
  over the season (820 against 818 mm). TimesFM and the baseline never
  postponed. Details:
  [reports/water_balance](reports/water_balance/README.md), section 6.

### Weather-model forecasts

Every day from May 2024 to September 2026 (864 forecasts), each issued the
evening before; full tables and figures in [reports/nwp](reports/nwp/README.md).

| Reference ET error (mm/day) | 1 day | 3 days | 7 days | 7-day total (mm) |
|---|---|---|---|---|
| Climatology | 0.59 | 0.59 | 0.59 | 2.8 |
| Damped anomaly persistence | 0.44 | 0.57 | 0.59 | 2.6 |
| TimesFM 2.5 | 0.43 | 0.57 | 0.60 | 2.7 |
| NCEP GFS, bias-corrected | 0.40 | 0.46 | 0.62 | 2.0 |
| ECMWF IFS 0.25°, bias-corrected | 0.27 | 0.32 | 0.54 | 1.7 |

- **ECMWF has the smallest error at every lead.** Its error is 35–45 % below
  TimesFM's 1–3 days ahead, and it still beats climatology at a week.
- **Rain events:** ECMWF caught 33 of the 48 cases with at least 10 mm of rain
  within three days, with 47 false alarms; GFS caught 28 with 39 false alarms.
  TimesFM and the baseline caught none.
- **Bias correction:** both models have systematic errors at this site (GFS:
  maximum temperature −1.3 °C and wind +0.8 m/s one day ahead). A running bias
  correction removes most of them, cutting GFS's one-day ETo error from 0.53 to
  0.40 mm/day.
- **Caveat:** ERA5, the observations here, comes from a model close to ECMWF's,
  which may favour ECMWF.

## Not done in this update

### From the 01 notes (original wording)

- Richards equation part should be much more complex:
  - The h parameter for Richards eq. should be adjusted to real, plausible soil statistics.
  - Different numerical solutions should be investigated other than FVM like: FDM, FEM, Mixed Finite Element, Explicit Euler, Implicit Euler, Newton-Raphson, Picard.
- Some uncertainty could be added to ET and \theta(x,y,z,t) functions.
- Both the forecast data and the solutions should be an input for another optimization process.
  (Forecasts now feed a simple irrigation rule; there is no optimization yet.)
- The soil and plant data could be researched to add more realistic parameters and constants.
  (Maize parameters now come from FAO-56 and Wesseling (1991); the soil is still
  the loam of 01.)
- The parameters and constants should be configured in such a way that they would show which equations is most susceptible to which parameter/constant and how.

### Known limitations

- **TimesFM is univariate and zero-shot:** each weather variable is forecast
  from its own past at one site, so its skill fades after about three days and
  it misses rain events. TimesFM's covariate support (XReg), which could take
  the weather-model forecasts as input, is not used.
- **Weather-model forecasts:**
  - only 2024-03-15 to 2026-09-30 (the archive), seven days ahead and one run
    per day, with no ensembles;
  - verified against ERA5, which ECMWF makes with a model close to IFS, so the
    comparison may favour ECMWF;
  - the bias correction is a simple running average per variable and lead, and
    rain is not corrected.
- **TimesFM on a CPU** takes about a second per forecast and about 2 GB of
  memory, so with TimesFM the dashboard issues a forecast once a week.
- **Top layer thickness:** soil evaporation depends on it. With the uniform
  grid, 5 cm layers gave about 30 % more evaporation than the 2.5 cm used now. A
  grid that is finer near the surface would remove most of this dependence.
- **Near saturation** the explicit solver needs very small sub-steps (heavy rain
  on fine soils, ponding), and it stops with an error if a layer would exceed
  saturation. An implicit scheme removes that limit.
- **Daily surface flux:** rain, irrigation and soil evaporation are combined into
  one net surface flux per day.
- **Irrigation is idealised:**
  - 100 % efficient and applied within one day;
  - always refills the root zone to field capacity;
  - no limits on amount, interval or available water;
  - forecasts are used only to postpone irrigation before rain.
- **Field capacity** is set at −100 cm. At −330 cm this soil would hold much
  less available water, so irrigation amounts depend on that choice.
- **Crop:**
  - the planting date is fixed at 1 May every year;
  - stage lengths do not respond to temperature;
  - there is no yield model, so water stress does not yet translate into yield.
- **Validation** is against computed references (ASCE-PM, Open-Meteo ET0,
  pyfao56, ERA5 soil moisture), not measurements, and covers one site. ERA5 is
  reanalysis, not station data.
- **Blaney-Criddle:** n/N is estimated from Rs/Ra, and the FAO-24 regression is
  used with 24-hour wind instead of daytime wind.
- **Calibration** fits one constant per equation; seasonal or regional
  coefficients are not explored.

## Next steps

1. **Yield response to water:** relate seasonal water stress to yield with the
   FAO-33 Ky factor, so irrigation strategies can be compared by yield and water
   use.
2. **Irrigation strategies and optimization:**
   - deficit irrigation, application efficiency, limits on amount and interval,
     a seasonal water budget;
   - optimize the trigger and amounts for water productivity, with
     weather-model forecasts of ET for the coming week.
3. **Soil solver:**
   - implicit Richards solver (mixed-form modified Picard, then Newton) on a
     grid refined near the surface, compared with the explicit one;
   - soil parameters from the Carsel & Parrish (1988) texture classes, or
     SoilGrids with Rosetta for the site.
4. **Sensitivity analysis:** Morris screening and Sobol indices of the ET
   equations and of the field water balance to their parameters.
5. **Uncertainty:** Monte Carlo on inputs and parameters; ensemble weather
   forecasts and TimesFM's quantiles in irrigation decisions.
6. **Forecasts:** TimesFM with the weather-model forecasts as covariates
   (XReg), and better bias correction (regression on several predictors).
7. **Measured data:** validation against measured ET (lysimeter or eddy
   covariance) and in-situ soil moisture, and more sites.
