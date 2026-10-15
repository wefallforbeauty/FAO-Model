"""Archived forecasts of numerical weather prediction (NWP) models.

Open-Meteo's Previous Runs API keeps, for every hour, the value each model
forecast 1, 2, ... 7 days before (24, 48, ... 168 hours earlier). This
module downloads those hourly series for the site, turns them into daily
weather for each lead day with the same daily statistics as the ERA5
observations, and loads them:

    python -m fao_model.nwp [--model all] [--start 2024-03-15] [--end 2026-09-30]

Forecast data by Open-Meteo.com (CC BY 4.0), from NCEP GFS (public domain)
and ECMWF IFS open data (CC BY 4.0).
"""
import argparse
import csv
import json
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import date, datetime, timedelta

from .data import DATA_DIR, SITE
from .meteo import wind_2m

API_URL = "https://previous-runs-api.open-meteo.com/v1/forecast"
MODELS = {"ecmwf_ifs025": "ECMWF IFS 0.25°", "gfs_seamless": "NCEP GFS"}
LEADS = tuple(range(1, 8))
HOURLY = ("temperature_2m", "relative_humidity_2m", "wind_speed_10m",
          "shortwave_radiation", "precipitation", "surface_pressure")
COLUMNS = ("T_max", "T_min", "RH_max", "RH_min", "RH_mean", "u10", "Rs", "precip", "P_hPa")
# Both models carry every variable at every lead from mid-March 2024 on.
START, END = "2024-03-15", "2026-09-30"
NWP_DIR = DATA_DIR / "nwp"


def csv_path(model):
    return NWP_DIR / f"konya_cumra_{model}_daily_by_lead.csv"


def fetch_hourly(model, lat, lon, start, end, chunk_days=183, timeout=600):
    """Hourly forecasts at every lead, downloaded in chunks: (metadata, {name: values})."""
    names = [f"{v}_previous_day{k}" for v in HOURLY for k in LEADS]
    hourly = defaultdict(list)
    meta = {}
    first, last = date.fromisoformat(start), date.fromisoformat(end)
    while first <= last:
        stop = min(last, first + timedelta(days=chunk_days - 1))
        params = {
            "latitude": lat, "longitude": lon,
            "start_date": first.isoformat(), "end_date": stop.isoformat(),
            "models": model, "timezone": "auto", "wind_speed_unit": "ms",
            "hourly": ",".join(names),
        }
        with urllib.request.urlopen(API_URL + "?" + urllib.parse.urlencode(params), timeout=timeout) as resp:
            payload = json.load(resp)
        if payload.get("error"):
            raise RuntimeError(payload.get("reason", "Open-Meteo request failed"))
        for key, values in payload["hourly"].items():
            hourly[key].extend(values)
        meta = {"grid_lat": payload["latitude"], "grid_lon": payload["longitude"],
                "elevation_m": payload["elevation"], "timezone": payload["timezone"]}
        first = stop + timedelta(days=1)
    return meta, dict(hourly)


def daily_by_lead(hourly):
    """Daily weather for each (local day, lead), only for days with all 24 hourly values.

    Temperature and humidity: maximum and minimum (and mean humidity); wind
    and pressure: mean; radiation: sum of the hourly means, W/m2 -> MJ/m2/day;
    rain: sum.
    """
    by_day = defaultdict(list)
    for i, stamp in enumerate(hourly["time"]):
        by_day[stamp[:10]].append(i)
    rows = []
    for day, idx in sorted(by_day.items()):
        if len(idx) != 24:
            continue
        for lead in LEADS:
            series = {v: [hourly[f"{v}_previous_day{lead}"][i] for i in idx] for v in HOURLY}
            if any(x is None for values in series.values() for x in values):
                continue
            t, rh = series["temperature_2m"], series["relative_humidity_2m"]
            rows.append({
                "date": day, "lead": lead,
                "T_max": max(t), "T_min": min(t),
                "RH_max": max(rh), "RH_min": min(rh), "RH_mean": sum(rh) / 24.0,
                "u10": sum(series["wind_speed_10m"]) / 24.0,
                "Rs": sum(series["shortwave_radiation"]) * 3600.0 / 1e6,
                "precip": sum(series["precipitation"]),
                "P_hPa": sum(series["surface_pressure"]) / 24.0,
            })
    return rows


def write_csv(rows, path, model, meta, lat, lon):
    header = {
        "source": "Open-Meteo Previous Runs API (https://open-meteo.com/), CC BY 4.0",
        "model": f"{model} ({MODELS.get(model, model)})",
        "site": SITE["name"],
        "requested_lat": lat,
        "requested_lon": lon,
        **meta,
        "lead": "days between the forecast and the day it is for (hourly values forecast 24*lead hours earlier)",
        "downloaded": date.today().isoformat(),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        for key, value in header.items():
            f.write(f"# {key}: {value}\n")
        writer = csv.writer(f)
        writer.writerow(("date", "lead") + COLUMNS)
        for row in rows:
            writer.writerow([row["date"], row["lead"]] + [round(row[c], 3) for c in COLUMNS])


def load(model, path=None):
    """(metadata, table) with table[lead][day] = daily forecast weather in the model's units."""
    meta, lines = {}, []
    with open(path or csv_path(model), newline="") as f:
        for line in f:
            if line.startswith("#"):
                key, _, value = line[1:].partition(":")
                meta[key.strip()] = value.strip()
            else:
                lines.append(line)
    table = defaultdict(dict)
    for row in csv.DictReader(lines):
        day = date.fromisoformat(row["date"])
        t_max, t_min = float(row["T_max"]), float(row["T_min"])
        table[int(row["lead"])][day] = {
            "timestamp": datetime(day.year, day.month, day.day),
            "T_max": t_max, "T_min": t_min, "T_mean": (t_max + t_min) / 2.0,
            "RH_max": float(row["RH_max"]), "RH_min": float(row["RH_min"]), "RH_mean": float(row["RH_mean"]),
            "u2": wind_2m(float(row["u10"]), 10.0),
            "Rs": float(row["Rs"]), "precip": float(row["precip"]), "P": float(row["P_hPa"]) / 10.0,
        }
    return meta, dict(table)


def main():
    parser = argparse.ArgumentParser(description="Download archived NWP forecasts from Open-Meteo.")
    parser.add_argument("--model", default="all", choices=("all",) + tuple(MODELS))
    parser.add_argument("--lat", type=float, default=SITE["lat"])
    parser.add_argument("--lon", type=float, default=SITE["lon"])
    parser.add_argument("--start", default=START)
    parser.add_argument("--end", default=END)
    args = parser.parse_args()
    for model in (MODELS if args.model == "all" else (args.model,)):
        meta, hourly = fetch_hourly(model, args.lat, args.lon, args.start, args.end)
        rows = daily_by_lead(hourly)
        write_csv(rows, csv_path(model), model, meta, args.lat, args.lon)
        print(f"{model}: {len(rows)} daily forecasts ({len(rows) // len(LEADS)} days x {len(LEADS)} leads) "
              f"to {csv_path(model)}")


if __name__ == "__main__":
    main()
