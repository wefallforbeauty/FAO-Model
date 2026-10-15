# Data

## `konya_cumra_era5_daily_1995_2024.csv`

Daily weather for the Konya Plain near Çumra, Türkiye (requested point
37.57° N, 32.78° E; ERA5 grid cell 37.50° N, 32.75° E; elevation 1013 m),
1995-01-01 to 2024-12-31 (10 958 days, no gaps), local days (Europe/Istanbul).

Why this site: the Konya Plain is Türkiye's largest irrigated agricultural
area (wheat, sugar beet, maize) in a water-scarce closed basin. Its
semi-arid continental climate has a strong seasonal cycle, including
monthly means below 0 °C, which exercises every branch of the ET
equations. The plain is flat, so one reanalysis grid cell represents it
reasonably well.

| Column | Unit | Description |
|---|---|---|
| `T_max`, `T_min` | °C | daily maximum / minimum air temperature at 2 m |
| `T_mean_24h` | °C | mean of the hourly temperatures (the model uses (T_max + T_min) / 2 as FAO-56 does) |
| `RH_max`, `RH_min`, `RH_mean` | % | relative humidity at 2 m |
| `u10` | m s⁻¹ | mean wind speed at 10 m (the loader converts it to 2 m, FAO-56 eq. 47) |
| `Rs` | MJ m⁻² day⁻¹ | incoming shortwave radiation |
| `precip` | mm day⁻¹ | precipitation |
| `P_hPa` | hPa | mean surface pressure |
| `ET0_openmeteo` | mm day⁻¹ | FAO-56 Penman-Monteith ET0 computed by Open-Meteo from hourly data, an independent reference for this repo's daily ASCE-PM |
| `soil_moisture_0_7cm` | m³ m⁻³ | ERA5 volumetric soil water, 0–7 cm (not used yet; kept for validating the soil model later) |

Regenerate it, or download another site or period:

```sh
python -m fao_model.data                       # this file
python -m fao_model.data --lat 36.9 --lon 39.0 --name "Harran Plain" \
    --start 2000-01-01 --end 2024-12-31 --out data/harran.csv
```

Caveat: these are reanalysis values for a ~25 km grid cell, not station
observations. Temperatures in winter can be milder than at stations in the
basin, and precipitation is spatially smoothed.

### Source and licence

Weather data by [Open-Meteo.com](https://open-meteo.com/) (historical
weather API, `models=era5`), licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Generated using Copernicus Climate Change Service information.

- Zippenfenig, P. (2023). Open-Meteo.com Weather API [Computer software].
  Zenodo. https://doi.org/10.5281/ZENODO.7970649
- Hersbach, H., et al. (2023). ERA5 hourly data on single levels from 1940
  to present [Data set]. ECMWF. https://doi.org/10.24381/cds.adbb2d47

## `konya_cumra_era5_daily_2025_2026.csv`

The same site, variables and source for 2025-01-01 to 2026-09-30 (638 days,
no gaps), kept in a separate file so the 1995-2024 reports stay
reproducible. `fao_model.data.load_all()` joins both files into one series
and checks that they meet without a gap; the dashboard and the weather-model
comparison use it. The most recent months of ERA5 are preliminary (ERA5T)
and can still be revised slightly.

```sh
python -m fao_model.data --start 2025-01-01 --end 2026-09-30 \
    --out data/konya_cumra_era5_daily_2025_2026.csv
```

## `nwp/`: archived weather-model forecasts

`konya_cumra_ecmwf_ifs025_daily_by_lead.csv` (ECMWF IFS 0.25°) and
`konya_cumra_gfs_seamless_daily_by_lead.csv` (NCEP GFS): daily weather that
each model forecast 1 to 7 days ahead, 2024-03-15 to 2026-09-30. ECMWF
misses a few days at some leads in April 2026.

| Column | Unit | Description |
|---|---|---|
| `date` | | the day the forecast is for (local time) |
| `lead` | days | forecast made `lead` × 24 hours before each hour of that day |
| `T_max`, `T_min` | °C | maximum / minimum of the hourly 2 m temperature |
| `RH_max`, `RH_min`, `RH_mean` | % | from the hourly 2 m relative humidity |
| `u10` | m s⁻¹ | mean 10 m wind speed |
| `Rs` | MJ m⁻² day⁻¹ | sum of the hourly mean shortwave radiation |
| `precip` | mm day⁻¹ | precipitation |
| `P_hPa` | hPa | mean surface pressure |

Only days with all 24 hourly values are kept. Regenerate:

```sh
python -m fao_model.nwp                        # both models, 2024-03-15 to 2026-09-30
```

Source: Open-Meteo's
[Previous Runs API](https://open-meteo.com/en/docs/previous-runs-api)
(CC BY 4.0), built from NCEP GFS (public domain, NOAA) and ECMWF IFS open
data (CC BY 4.0).
