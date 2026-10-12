# FAO-Model 03

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
  to field capacity when 55 % of the available water is used.
- **Soil water**: 3D van Genuchten–Mualem Richards solver (mass-conservative
  explicit finite volumes with adaptive sub-steps), a 1.5 m column in 2.5 cm
  layers.
- **Forecasts** of every series: TimesFM when installed, otherwise a
  damped-trend fallback.
- **Data**: 30 years (1995–2024) of daily ERA5 weather for the Konya Plain,
  Türkiye, from Open-Meteo.

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

With the environment active:

```sh
python -m fao_model                     # dashboard, replays the Konya record
python -m fao_model --synthetic         # synthetic daily weather instead
python -m fao_model --interval-ms 200   # faster replay (default 1000 ms per day)
python -m fao_model --irrigation none   # rainfed maize instead of automatic irrigation
python -m fao_model.validate            # regenerates reports/validation/
python -m fao_model.water_balance       # regenerates reports/water_balance/ (about 3 minutes)
python -m fao_model.data --help         # download another site or period
pytest
```

On systems where `python` points to Python 2 or is missing, use `python3`.
If the dashboard cannot draw the live graph, it shows the import error that
names the missing piece, usually Tk support in Python or in Pillow.

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
| `fao_model/forecast.py` | TimesFM forecaster with fallback |
| `fao_model/simulation.py` | the daily model loop, without GUI |
| `fao_model/app.py` | Tkinter dashboard |
| `fao_model/metrics.py`, `fao_model/validate.py` | agreement statistics and the ET validation report |
| `fao_model/water_balance.py` | maize water balance report: scenarios, pyfao56 and ERA5 comparisons |
| `fao_model/plotting.py` | shared figure style of the reports |
| `tests/` | FAO-56 examples, cross-checks against pyfao56, solver, crop, field and data tests |

## Data

`data/konya_cumra_era5_daily_1995_2024.csv`: daily ERA5 weather for the Konya
Plain near Çumra (37.57° N, 32.78° E, 1013 m), 10 958 days without gaps. Weather
data by [Open-Meteo.com](https://open-meteo.com/) (CC BY 4.0), generated using
Copernicus Climate Change Service information. Columns, caveats and citations:
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

## Not done in this update

### From the 01 notes (original wording)

- Then the timesfm forecast should be applicable weather data.
- Richards equation part should be much more complex:
  - The h parameter for Richards eq. should be adjusted to real, plausible soil statistics.
  - Different numerical solutions should be investigated other than FVM like: FDM, FEM, Mixed Finite Element, Explicit Euler, Implicit Euler, Newton-Raphson, Picard.
- Some uncertainty could be added to ET and \theta(x,y,z,t) functions.
- Both the forecast data and the solutions should be an input for another optimization process.
- The soil and plant data could be researched to add more realistic parameters and constants.
  (Maize parameters now come from FAO-56 and Wesseling (1991); the soil is still
  the loam of 01.)
- The parameters and constants should be configured in such a way that they would show which equations is most susceptible to which parameter/constant and how.

### Known limitations

- **TimesFM is not installed** here, so every forecast comes from the
  damped-trend fallback. Also:
  - The loading code (`TimesFm(hparams)` with `load_from_checkpoint`) has not
    been checked against a current TimesFM release.
  - Forecast errors are swallowed silently.
  - The context length (2048) is larger than timesfm-1.0-200m supports (512).
  - It forecasts the model outputs, not the weather inputs.
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
  - no use of weather forecasts.
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

1. **TimesFM:**
   - install the current release and fix the loading code;
   - forecast the weather inputs (temperature, humidity, wind, radiation, rain)
     and compute ET from them;
   - measure forecast skill on 2015–2024.
2. **Irrigation scheduling with forecasts:** delay or reduce irrigation when rain
   is forecast, add application efficiency and limits on amount and interval,
   and compare strategies, for example full refills against deficit irrigation.
   This leads into the optimization item.
3. **Yield response to water:** relate seasonal water stress to yield with the
   FAO-33 Ky factor, so irrigation strategies can be compared by yield and water
   use.
4. **Soil solver:**
   - implicit Richards solver (mixed-form modified Picard, then Newton) on a grid
     refined near the surface, compared with the explicit one;
   - soil parameters from the Carsel & Parrish (1988) texture classes, or
     SoilGrids with Rosetta for the site.
5. **Sensitivity analysis:** Morris screening and Sobol indices of the ET
   equations and of the field water balance to their parameters.
6. **Uncertainty:** Monte Carlo on inputs and parameters, and TimesFM quantiles.
7. **Measured data:** validation against measured ET (lysimeter or eddy
   covariance) and in-situ soil moisture, and more sites.
