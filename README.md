# FAO Evapotranspiration and Soil-Water Model

A desktop research prototype for exploring the relationship between weather,
reference evapotranspiration (ET), soil-water movement, and time-series
forecasting.

The application generates synthetic weather observations, compares several ET
equations, advances a small three-dimensional Richards-equation model, and
shows the results in a live Tkinter dashboard. When a compatible TimesFM model
is unavailable, the application automatically uses a lightweight trend
forecast instead.

> **Research status:** This project is an exploratory model, not a validated
> irrigation, agronomic, or operational forecasting tool. It currently uses
> simulated weather and simplified soil physics.

## What the application includes

- Synthetic weather generation at configurable time intervals
- ASCE Penman-Monteith reference ET
- Maize crop ET using a crop coefficient
- Thornthwaite, Blaney-Criddle, Turc, Priestley-Taylor, Hargreaves,
  Jensen-Haise, and Abtew estimates
- Online least-squares calibration of several methods against ASCE
  Penman-Monteith
- A `5 x 5 x 5` finite-volume Richards-equation prototype using van Genuchten
  soil relationships
- A 96-step forecast for each result series
- A live text display and, when Matplotlib is installed, a live chart

## Project structure

```text
FAO-Model/
|-- 01.py       # Model equations, simulation, forecasting, and GUI
`-- README.md   # Project documentation
```

The project is currently implemented as one Python script. The main sections
of `01.py` are:

1. Radiation and evapotranspiration functions
2. Synthetic weather generator
3. ET coefficient calibrator
4. van Genuchten functions and the 3D Richards solver
5. TimesFM integration and fallback forecaster
6. Tkinter desktop dashboard

## Requirements

- Python 3.10 or newer is recommended
- Tkinter (normally included with Python on Windows and macOS)
- [NumPy](https://numpy.org/) — required
- [Matplotlib](https://matplotlib.org/) — optional, for the live chart
- [TimesFM](https://github.com/google-research/timesfm) — optional, for model-based
  forecasts

Install the standard Python dependencies:

```bash
python -m pip install numpy matplotlib
```

On Windows, use `py` in place of `python` if that is the launcher available on
your system (for example, `py -m pip install numpy matplotlib` and
`py 01.py`).

On Linux, Tkinter may need to be installed through the operating system's
package manager. TimesFM is intentionally not included in the basic install:
its API and runtime requirements vary by release, and the integration in
`01.py` uses the `TimesFm`/`TimesFmHparams` API and the
`google/timesfm-1.0-200m` checkpoint.

## Run the dashboard

From the project directory:

```bash
python 01.py
```

The application updates once per second. Each update advances the synthetic
weather timestamp by 10 minutes, recalculates all ET methods, advances the soil
model, refreshes forecasts, and updates the display.

If TimesFM cannot be imported or initialized, the status panel reports the
reason and the application continues with its built-in fallback forecaster. If
Matplotlib is missing, the text dashboard still runs without the chart.

## How to read the output

The left panel displays:

- Current synthetic weather inputs
- Coefficients fitted against ASCE Penman-Monteith
- ET estimates in `mm/day`
- Mean Richards-model pressure head in `cm`
- Mean volumetric water content
- Forecast backend status and forecast values

The chart shows recent values for selected ET methods. Dashed lines represent
their forecasts.

## Important limitations

- Weather is synthetic; there is no file, sensor, or weather-API input yet.
- Most ET equations produce daily values, while the simulation clock advances
  in 10-minute steps. Temporal aggregation and unit consistency need to be
  resolved before scientific use.
- Thornthwaite normally uses monthly mean temperatures and an annual heat
  index; the current per-step implementation is experimental.
- The non-ASCE methods are calibrated to ASCE Penman-Monteith, so agreement
  between methods is not independent validation.
- The Richards solver is an explicit, clipped prototype. It has no root-water
  uptake, adaptive time stepping, convergence analysis, or automated
  water-mass conservation test.
- Soil and crop parameters are fixed demonstration values rather than
  site-calibrated measurements.
- Forecasts are made from generated model outputs, not real historical weather.
- Historical data grows for the lifetime of the application, although the
  chart displays only the most recent 1,200 points.

## Development priorities

1. Define and enforce units and time scales throughout the weather, ET, and
   Richards components.
2. Add real historical weather ingestion, cleaning, and validation.
3. Add automated tests for ET equations, boundary conditions, numerical
   stability, and water-mass conservation.
4. Validate each ET method against observations instead of treating ASCE
   Penman-Monteith as ground truth.
5. Improve the Richards model with realistic soil profiles, root-water uptake,
   adaptive time stepping, and alternative implicit/nonlinear solvers.
6. Correct Thornthwaite's monthly heat-index workflow.
7. Quantify uncertainty, bias, and parameter sensitivity.
8. Feed validated forecasts and soil states into an irrigation optimization
   stage.
9. Split the single script into testable modules and add a reproducible
   dependency specification.

## Contributing

Contributions are welcome, especially for validation datasets, dimensional
checks, numerical tests, and clearer separation of the scientific model from
the GUI. When changing an equation, document its source, expected units, time
scale, and valid input range.
