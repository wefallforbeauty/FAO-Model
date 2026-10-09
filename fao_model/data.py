"""Daily weather from the Open-Meteo historical weather API (ERA5 reanalysis).

Weather data by Open-Meteo.com (CC BY 4.0), generated from the Copernicus
ERA5 reanalysis. `python -m fao_model.data` (re)downloads the default
dataset into data/; see `--help` for other sites and periods.
"""
import argparse
import csv
import json
import urllib.parse
import urllib.request
from datetime import date, datetime
from pathlib import Path

from .meteo import wind_2m

ARCHIVE_URL = "https://archive-api.open-meteo.com/v1/archive"

# Konya Plain (Cumra): Turkiye's largest irrigated plain, semi-arid continental
# climate, flat terrain that a reanalysis grid cell represents well.
SITE = {"name": "Konya Plain (Cumra), Turkiye", "lat": 37.57, "lon": 32.78}
START, END = "1995-01-01", "2024-12-31"
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DEFAULT_CSV = DATA_DIR / "konya_cumra_era5_daily_1995_2024.csv"

# Open-Meteo daily variable -> CSV column (units as delivered by the API)
DAILY_VARS = {
    "temperature_2m_max": "T_max",                          # °C
    "temperature_2m_min": "T_min",                          # °C
    "temperature_2m_mean": "T_mean_24h",                    # °C, mean of hourly values
    "relative_humidity_2m_max": "RH_max",                   # %
    "relative_humidity_2m_min": "RH_min",                   # %
    "relative_humidity_2m_mean": "RH_mean",                 # %
    "wind_speed_10m_mean": "u10",                           # m s-1 at 10 m
    "shortwave_radiation_sum": "Rs",                        # MJ m-2 day-1
    "precipitation_sum": "precip",                          # mm day-1
    "surface_pressure_mean": "P_hPa",                       # hPa
    "et0_fao_evapotranspiration": "ET0_openmeteo",          # mm day-1, FAO-56 PM from hourly data
    "soil_moisture_0_to_7cm_mean": "soil_moisture_0_7cm",   # m3 m-3
}


def fetch(lat, lon, start, end, model="era5", timeout=900):
    """Download daily data from the Open-Meteo archive API; returns the JSON payload."""
    params = {
        "latitude": lat,
        "longitude": lon,
        "start_date": start,
        "end_date": end,
        "daily": ",".join(DAILY_VARS),
        "models": model,
        "timezone": "auto",
        "wind_speed_unit": "ms",
    }
    url = ARCHIVE_URL + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        payload = json.load(resp)
    if payload.get("error"):
        raise RuntimeError(payload.get("reason", "Open-Meteo request failed"))
    return payload


def write_csv(payload, path, site_name, lat, lon, model="era5"):
    meta = {
        "source": "Open-Meteo historical weather API (https://open-meteo.com/), CC BY 4.0",
        "model": model,
        "site": site_name,
        "requested_lat": lat,
        "requested_lon": lon,
        "grid_lat": payload["latitude"],
        "grid_lon": payload["longitude"],
        "elevation_m": payload["elevation"],
        "timezone": payload["timezone"],
        "downloaded": date.today().isoformat(),
    }
    daily = payload["daily"]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        for key, value in meta.items():
            f.write(f"# {key}: {value}\n")
        writer = csv.writer(f)
        writer.writerow(["date"] + list(DAILY_VARS.values()))
        for i, day in enumerate(daily["time"]):
            values = [daily[var][i] for var in DAILY_VARS]
            writer.writerow([day] + ["" if v is None else v for v in values])


def load(path=DEFAULT_CSV):
    """Read a CSV written by write_csv; returns (metadata, daily records).

    Records use the model's keys and units: T_mean = (T_max + T_min) / 2
    as in FAO-56, wind converted to 2 m, pressure in kPa.
    """
    meta, lines = {}, []
    with open(path, newline="") as f:
        for line in f:
            if line.startswith("#"):
                key, _, value = line[1:].partition(":")
                meta[key.strip()] = value.strip()
            else:
                lines.append(line)
    lat = float(meta["requested_lat"])
    elev = float(meta["elevation_m"])
    records = []
    for row in csv.DictReader(lines):
        day = date.fromisoformat(row["date"])
        T_max, T_min = float(row["T_max"]), float(row["T_min"])
        records.append({
            "timestamp": datetime(day.year, day.month, day.day),
            "doy": day.timetuple().tm_yday,
            "lat": lat,
            "elev": elev,
            "T_mean": (T_max + T_min) / 2.0,
            "T_max": T_max,
            "T_min": T_min,
            "RH_mean": float(row["RH_mean"]),
            "RH_max": float(row["RH_max"]),
            "RH_min": float(row["RH_min"]),
            "u2": wind_2m(float(row["u10"]), 10.0),
            "Rs": float(row["Rs"]),
            "P": float(row["P_hPa"]) / 10.0,
            "precip": float(row["precip"]),
            "ET0_ref": float(row["ET0_openmeteo"]),
        })
    return meta, records


def main():
    parser = argparse.ArgumentParser(description="Download daily ERA5 weather from Open-Meteo.")
    parser.add_argument("--name", default=SITE["name"])
    parser.add_argument("--lat", type=float, default=SITE["lat"])
    parser.add_argument("--lon", type=float, default=SITE["lon"])
    parser.add_argument("--start", default=START)
    parser.add_argument("--end", default=END)
    parser.add_argument("--model", default="era5")
    parser.add_argument("--out", default=str(DEFAULT_CSV))
    args = parser.parse_args()
    payload = fetch(args.lat, args.lon, args.start, args.end, args.model)
    write_csv(payload, args.out, args.name, args.lat, args.lon, args.model)
    print(f"Wrote {len(payload['daily']['time'])} days to {args.out}")


if __name__ == "__main__":
    main()
