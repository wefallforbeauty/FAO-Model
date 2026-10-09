# Update log

## 02 (2026-10-09)

Pull request #1. The single-file program became a tested Python package that
runs on 30 years of real daily weather. Several unit and sign errors in 01 were
found and fixed.

### Structure

- `01.py` was split into the `fao_model` package: `meteo`, `et`, `calibration`,
  `weather`, `data`, `richards`, `forecast`, `simulation` (model loop without
  GUI), `app` (Tk dashboard), `metrics` and `validate`. The split itself did not
  change behaviour: results, coefficients and forecasts were bit-identical to
  `01.py` over 600 seeded steps.
- Added `pyproject.toml`, `.gitignore` and a pytest suite (51 tests).

### Evapotranspiration equations

All equations now return mm/day and are checked against FAO-56 worked examples
or hand calculations.

- **ASCE-PM:**
  - Saturation vapour pressure now comes from Tmax and Tmin (FAO-56 eq. 12).
  - Actual vapour pressure now comes from RHmax and RHmin (eq. 17).
  - FAO-56 Example 18 now gives 3.88 mm/day (published 3.9; 01 gave 3.75).
- **Thornthwaite:**
  - 01 returned mm/month and built the heat index from the current temperature
    only. That is why it sat at the top of the graph (155 on Example 18).
  - It now uses the annual heat index of the site, the trailing 30-day mean
    temperature, and Willmott et al. (1985) above 26.5 °C, and converts the result
    to mm/day.
  - A calendar-month version is used in the validation.
- **Hargreaves:** Ra is converted to mm/day (× 0.408).
- **Turc:** Rs is converted to cal cm⁻² day⁻¹ (× 23.8846).
- **Jensen-Haise:** Rs / λ instead of the unexplained `/ 28.3`.
- **Blaney-Criddle:**
  - p is now the daily percentage of annual daytime hours (≈ 0.27), not day
    length / 24.
  - Without calibration, a and b come from the FAO-24 regression (Frevert et al.,
    1983).
- **Calibration:**
  - In 01 the calibration hid these errors: Turc a = 4.25, Jensen-Haise a = 13.9,
    Hargreaves a = 0.35, Priestley-Taylor α = 2.17, Blaney-Criddle a = −16.3.
  - It now starts from the literature coefficients, needs 30 days, and only
    uses days before the one being computed (walk-forward).
- **Time scale:** 01 fed 10-minute synthetic weather to daily equations. The
  generator now produces daily weather.

### Real weather data

- `fao_model.data` downloads daily ERA5 weather from the Open-Meteo archive API
  (free, no key, CC BY 4.0) and converts it to the model's units: wind at 2 m,
  pressure in kPa, T_mean = (Tmax + Tmin) / 2.
- Bundled dataset: Konya Plain near Çumra, Türkiye, 1995–2024, 10 958 days
  without gaps. It also contains Open-Meteo's own FAO-56 ET0 and ERA5 0–7 cm
  soil moisture.

### Validation

`python -m fao_model.validate` writes `reports/validation/`: training period
1995–2014, test period 2015–2024, reference = this repo's daily ASCE-PM.

- ASCE-PM matches pyfao56 to within 0.005 mm/day. Against Open-Meteo's
  hourly-based ET0: r = 0.995, bias −0.18 mm/day (−4.8 %).
- Test period results:

  | Method | Relative bias, literature coefficients | RMSE (mm/day), literature → calibrated |
  |---|---|---|
  | Turc | −2.7 % | 0.67 → 0.66 |
  | Jensen-Haise | +4.5 % | 0.87 → 0.68 |
  | Abtew | +6.0 % | 0.98 → 0.98 |
  | Hargreaves | −10.4 % | 0.78 → 0.66 |
  | Blaney-Criddle (FAO-24) | +17.0 % | 1.17 → 0.63 |
  | Priestley-Taylor | −19.5 % | 0.97 → 0.65 |
  | Thornthwaite (calendar months) | −38.6 % | not calibrated |

- In 01 every equation converged to ASCE-PM because all coefficients were refit
  to it at every step on the same samples. That convergence came from the
  method, not from the physics.
- Calibrated coefficients take on the bias of their reference: fitted to
  Open-Meteo ET0 instead of ASCE-PM, they shift by 3–4 %.

### Soil water (Richards solver)

- **Gravity:** it pointed the wrong way. 01 put layer 0 at the bottom, but rain
  entered layer 0. Layer 0 is now the surface.
- **Units:** fluxes between cells and free drainage were not divided by the cell
  size; all terms are now 1/day.
- **Time step:** each weather step advanced the soil by 1 s (10 × 0.1 s). The
  solver now advances one day in cm and day units. Ks is 21.6 cm/day, the same
  value as 0.00025 cm/s. Sub-steps adapt to the explicit stability limit.
- **Mass conservation:**
  - Water content is the conserved state and h follows from an inverse
    van Genuchten function. This replaces the ±2 cm update clip.
  - The water balance closes to about 1e-14 cm over 30 years.
- **Surface boundary:** rain beyond what a near-saturated surface admits
  becomes runoff, and evaporation is limited by a dry surface (h_min = −20000 cm).
- **Conductivity between cells:** now the arithmetic mean. The full 30-year
  replay stopped on 2005-11-04 because the harmonic mean stalled a wetting front
  entering dry soil.
- **Tests:** water balance, drainage direction, steady state under constant
  infiltration (K(h) = q, unit gradient), column symmetry, flux scaling, runoff,
  evaporation limit, and the wetting front of 2005-11-04.

### Dashboard

- `python -m fao_model` replays the Konya record one day per tick.
  - `--synthetic` uses generated weather instead.
  - `--data` loads another CSV.
  - `--interval-ms` sets the speed.
- It stops cleanly at the end of the data.
- The log keeps the last 5000 lines instead of growing without limit. It shows
  Open-Meteo ET0, soil water storage, daily water balance terms and the balance
  error.
- If the live graph cannot load, the dashboard shows the actual import error
  instead of "Install matplotlib".

### Behaviour changes

- One dashboard step is one day (it was 10 minutes).
- `RichardsSolver3D.step(q_top, dt)` takes cm/day and days (it was cm/s and
  seconds). Set the pressure head with `set_head`.
- Literature coefficients are used until 30 days of history exist. Blaney-Criddle
  uses the daily FAO-24 a and b when `bc_a` and `bc_b` are `None`.

### Items addressed from the 01 notes

| 01 note | Done in 02 |
|---|---|
| The program should be written for real weather data or real weather data history. | Konya 1995–2024 ERA5 data drives the dashboard and the validation. |
| Adaptive time stamp could be added. | The Richards solver adapts its sub-steps to the stability limit. |
| Currently the ASCE-PM is taken reference for the calibration … forces them to converge to ASCE-PM … | Investigated: literature against calibrated coefficients, training against test period, and two references in the validation report. |
| A validation code should be added. | FAO-56 examples, pyfao56 cross-check, Open-Meteo ET0, validation report. |
| Thornthwaite eq. should be adjusted for monthly mean temp values … | Monthly method with the annual heat index; 01 returned mm/month. |
| … investigate if the convergence … is because of any bias … relative bias … scatter plot … | Bias, relative bias, RMSE, r and NSE, scatter plots and a relative-bias chart. |
| All equations, variables and numerical fluxes should be checked for dimensional consistency … | Done for the ET equations and the Richards fluxes, with unit tests. |
| Add automated water-mass conservation tests for the Richards solver. | `tests/test_richards.py` |

## 01 (2026-08-19)

- Single-file Tk dashboard `01.py`: synthetic 10-minute weather, nine ET
  equations calibrated to ASCE-PM at every step, an explicit 3D finite-volume
  Richards solver, and TimesFM forecasts with a fallback.
- README: the list of planned improvements.
