"""Crop water balance of grain maize at the data site, rainfed and irrigated.

Runs FieldModel (FAO-56 dual crop coefficients over the Richards solver)
through every season of the record, compares the seasons with pyfao56
(the FAO-56 water balance as a bucket model) when it is installed, and
compares the simulated topsoil moisture with ERA5. Writes a Markdown
report with tables and figures:

    python -m fao_model.water_balance [--data CSV] [--first-year Y] [--last-year Y] [--out DIR]

From 2015 on it also compares irrigation that is postponed when rain is
forecast (perfect forecast, TimesFM, statistical baseline) with the plain rule.
"""
import argparse
import contextlib
import copy
import io
from datetime import date, datetime, timedelta
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from . import data as weather_data
from .crop import Crop, adjust_for_climate
from .et import asce_pm
from .field import COLUMN, H_FIELD_CAPACITY, FieldModel
from .forecast import Climatology, WeatherForecaster, add_eto
from .forecast_skill import SPLIT, cached
from .metrics import compare
from .plotting import INK, INK_2, MUTED, SERIES, SURFACE, style
from .richards import RichardsSolver3D

REPORT_DIR = Path(__file__).resolve().parent.parent / "reports" / "water_balance"
SCENARIOS = {"rainfed": False, "irrigated": True}
TOTALS = ("precip", "irrigation", "eto", "tp", "ta", "ep", "ea", "drainage", "runoff")
TOP_CM = 7.5    # topsoil depth compared with ERA5's 0-7 cm layer


def run_field(records, eto, crop, irrigation, snapshot_at=None):
    """Daily results of one scenario, the soil water profile on each planting day, and a copy
    of the field as it was at the start of `snapshot_at` (None if not reached)."""
    field = FieldModel(crop, irrigation=irrigation, soil=RichardsSolver3D(**COLUMN))
    n_top = int(round(TOP_CM / field.soil.dz))
    days, profiles, snapshot = [], {}, None
    for d, e in zip(records, eto):
        day = d["timestamp"].date()
        if d["timestamp"] == snapshot_at:
            snapshot = copy.deepcopy(field)
        if crop.day_of_season(day) == 1:
            profiles[day.year] = field.soil.theta.mean(axis=(0, 1)).copy()
        out = field.step(d, e)
        theta_top = float(field.soil.theta[:, :, :n_top].mean())
        days.append(dict(out, date=day, theta_top=theta_top))
    return field, days, profiles, snapshot


def seasons(days):
    """Season totals (mm) by year, with the number of irrigations."""
    result = {}
    for d in days:
        if d["t"] is None:
            continue
        s = result.setdefault(d["date"].year, dict.fromkeys(TOTALS + ("irrigations", "postponed"), 0.0))
        for key in TOTALS:
            s[key] += d[key]
        s["irrigations"] += d["irrigation"] > 0
        s["postponed"] += d.get("postponed", False)
    return result


def run_pyfao56(records, eto, crop, field, profiles, irrigation):
    """Season totals from pyfao56 with the same crop, soil limits, weather and start profile.

    Returns None if pyfao56 is not installed. pyfao56 gets a layered soil
    profile (5 cm layers) holding this model's water content on each planting
    day, so both models start the season from the same soil water.
    """
    try:
        import pandas as pd
        import pyfao56 as fao
    except ImportError:
        return None
    wth = fao.Weather()
    wth.wdata = pd.DataFrame(
        {
            "Srad": [d["Rs"] for d in records],
            "Tmax": [d["T_max"] for d in records],
            "Tmin": [d["T_min"] for d in records],
            "Vapr": np.nan,
            "Tdew": np.nan,
            "RHmax": [d["RH_max"] for d in records],
            "RHmin": [d["RH_min"] for d in records],
            "Wndsp": [d["u2"] for d in records],
            "Rain": [d["precip"] for d in records],
            "ETref": eto,
            "MorP": "M",
        },
        index=[d["timestamp"].strftime("%Y-%j") for d in records],
    )
    wth.z, wth.lat, wth.wndht = records[0]["elev"], records[0]["lat"], 2.0
    per_layer = int(round(5.0 / field.soil.dz))
    result = {}
    for year, profile in profiles.items():
        start = crop.planting_date(year)
        end = start + timedelta(days=crop.season_length - 1)
        if end.strftime("%Y-%j") not in wth.wdata.index:
            continue
        sol = fao.SoilProfile()
        theta5 = profile[: len(profile) // per_layer * per_layer].reshape(-1, per_layer).mean(axis=1)
        for i, theta0 in enumerate(theta5):
            sol.sdata.loc[5 * (i + 1)] = [field.theta_fc, field.theta_wp, float(theta0)]
        par = fao.Parameters(
            Kcbini=crop.kcb_ini, Kcbmid=crop.kcb_mid, Kcbend=crop.kcb_end,
            Lini=crop.l_ini, Ldev=crop.l_dev, Lmid=crop.l_mid, Lend=crop.l_late,
            hini=crop.h_ini, hmax=crop.h_max, thetaFC=field.theta_fc, thetaWP=field.theta_wp,
            theta0=float(theta5[0]), Zrini=crop.zr_ini, Zrmax=crop.zr_max, pbase=crop.p,
            Ze=0.10, REW=9.0,   # FAO-56 Table 19, loam
        )
        autoirr = None
        if irrigation:
            autoirr = fao.AutoIrrigate()
            autoirr.addset(start.strftime("%Y-%j"), end.strftime("%Y-%j"), mad=crop.p)
        model = fao.Model(start.strftime("%Y-%j"), end.strftime("%Y-%j"), par, wth,
                          autoirr=autoirr, sol=sol)
        with contextlib.redirect_stdout(io.StringIO()):
            model.run()
        o = model.odata
        result[year] = {"ta": o["T"].sum(), "ea": o["E"].sum(), "irrigation": o["Irrig"].sum(),
                        "irrigations": int((o["Irrig"] > 0).sum()), "drainage": o["DP"].sum()}
    return result


def forecast_scenarios(records, eto, crop, snapshot, start):
    """Irrigated seasons from records[start] on, postponing irrigation when rain is forecast.

    All scenarios continue from the same field state (snapshot). Forecasts of
    the next 3 days' rain are issued the evening before each season day:
    perfect (the rain that fell), TimesFM 2.5 and the statistical baseline.
    """
    in_season = [i for i in range(start, len(records))
                 if crop.day_of_season(records[i]["timestamp"].date()) is not None]
    rain = np.array([d["precip"] for d in records])
    stamps = [records[i]["timestamp"] for i in in_season]
    tables = {"Perfect forecast": dict(zip(stamps, [rain[i:i + 3].sum() for i in in_season]))}
    climatology = Climatology(add_eto(records[:start]))
    for name, backend in (("TimesFM 2.5", "timesfm"), ("Damped anomaly persistence", "baseline")):
        try:
            forecaster = WeatherForecaster(climatology, backend=backend, batch_size=32)
        except Exception as exc:
            print(f"{name} not available ({type(exc).__name__}); skipped.")
            continue
        print(f"{name}: rain forecasts for {len(in_season)} season days")
        key = f"rain3_{backend}_{records[start]['timestamp']:%Y%m%d}_{len(records)}"
        totals = cached(key, lambda: {"total": forecaster.rain_totals(records, [i - 1 for i in in_season], 3)})
        tables[name] = dict(zip(stamps, totals["total"]))
        del forecaster
    results = {}
    for name, table in tables.items():
        field = copy.deepcopy(snapshot)
        field.rain_forecast = lambda d, table=table: table.get(d["timestamp"], 0.0)
        results[name] = [dict(field.step(d, e), date=d["timestamp"].date())
                         for d, e in zip(records[start:], eto[start:])]
    return results


def irrigated_before_rain(days, records_by_date):
    """Irrigations followed by at least 10 mm of rain on that day or the next two."""
    count = 0
    for d in days:
        if d["irrigation"] > 0:
            rain = sum(records_by_date.get(d["date"] + timedelta(days=k), 0.0) for k in range(3))
            count += rain >= 10.0
    return count


# ---------------------------------------------------------------- figures

def _legend(fig_or_ax, handles, **kwargs):
    fig_or_ax.legend(handles=handles, frameon=False, labelcolor=INK_2, fontsize=8, **kwargs)


def example_season_chart(path, year, irrigated_days, rainfed_days):
    fig = Figure(figsize=(10, 6.4), facecolor=SURFACE)
    top, bottom = fig.subplots(2, 1, sharex=True)
    x = [d["date"] for d in irrigated_days]
    top.bar(x, [d["irrigation"] for d in irrigated_days], width=0.9, color=SERIES[1], label="Irrigation")
    top.plot(x, [d["dr"] for d in irrigated_days], color=SERIES[0], linewidth=2, label="Root-zone depletion Dr")
    top.plot(x, [d["taw"] * 0.55 for d in irrigated_days], color=INK_2, linewidth=1,
             label="Irrigation threshold (55 % of TAW)")
    top.plot(x, [d["taw"] for d in irrigated_days], color=MUTED, linewidth=1, label="Total available water TAW")
    top.set_ylabel("mm")
    top.set_title(f"Irrigated maize, {year}: depletion and irrigation", fontsize=10)
    style(top)
    _legend(top, top.get_legend_handles_labels()[0], loc="upper left")
    bottom.plot(x, [d["tp"] for d in irrigated_days], color=INK_2, linewidth=1, label="Potential transpiration")
    bottom.plot(x, [d["ta"] for d in irrigated_days], color=SERIES[0], linewidth=2, label="Actual, irrigated")
    bottom.plot(x, [d["ta"] for d in rainfed_days], color=SERIES[1], linewidth=2, label="Actual, rainfed")
    bottom.set_ylabel("Transpiration (mm/day)")
    bottom.set_title("Transpiration", fontsize=10)
    style(bottom)
    _legend(bottom, bottom.get_legend_handles_labels()[0], loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def seasons_chart(path, by_scenario):
    years = sorted(by_scenario["irrigated"])
    fig = Figure(figsize=(10, 6), facecolor=SURFACE)
    top, bottom = fig.subplots(2, 1, sharex=True)
    top.bar(years, [by_scenario["irrigated"][y]["irrigation"] for y in years], width=0.7, color=SERIES[0])
    top.set_ylabel("Irrigation (mm)")
    top.set_title("Seasonal irrigation, irrigated scenario", fontsize=10)
    style(top)
    for name, color in (("irrigated", SERIES[0]), ("rainfed", SERIES[1])):
        s = by_scenario[name]
        bottom.plot(years, [s[y]["ta"] / s[y]["tp"] for y in years], color=color, linewidth=2,
                    marker="o", markersize=4, label=name.capitalize())
    bottom.set_ylim(0, 1.05)
    bottom.set_ylabel("Actual / potential transpiration")
    bottom.set_title("Water stress by season", fontsize=10)
    style(bottom)
    _legend(bottom, bottom.get_legend_handles_labels()[0], loc="center right")
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def pyfao56_chart(path, ours, theirs):
    fig = Figure(figsize=(10, 4.6), facecolor=SURFACE)
    left, right = fig.subplots(1, 2)
    hi = 0.0
    for name, color in (("irrigated", SERIES[0]), ("rainfed", SERIES[1])):
        years = sorted(theirs[name])
        a = [theirs[name][y]["ta"] + theirs[name][y]["ea"] for y in years]
        b = [ours[name][y]["ta"] + ours[name][y]["ea"] for y in years]
        left.scatter(a, b, s=36, color=color, edgecolors=SURFACE, linewidths=1.5, label=name.capitalize())
        hi = max(hi, max(a), max(b))
    left.plot([0, hi * 1.05], [0, hi * 1.05], color=MUTED, linewidth=1)
    left.set_xlabel("pyfao56 (mm)")
    left.set_ylabel("This model (mm)")
    left.set_title("Seasonal actual ET (T + E)", fontsize=10)
    style(left)
    _legend(left, left.get_legend_handles_labels()[0], loc="upper left")
    years = sorted(theirs["irrigated"])
    a = [theirs["irrigated"][y]["irrigation"] for y in years]
    b = [ours["irrigated"][y]["irrigation"] for y in years]
    right.scatter(a, b, s=36, color=SERIES[0], edgecolors=SURFACE, linewidths=1.5)
    lo, hi = 0.9 * min(a + b), 1.05 * max(a + b)
    right.plot([lo, hi], [lo, hi], color=MUTED, linewidth=1)
    right.set_xlabel("pyfao56 (mm)")
    right.set_ylabel("This model (mm)")
    right.set_title("Seasonal irrigation", fontsize=10)
    style(right)
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def era5_chart(path, days, show_years):
    fig = Figure(figsize=(10, 6.4), facecolor=SURFACE)
    top, bottom = fig.subplots(2, 1)
    shown = [d for d in days if d["date"].year in show_years and d["sm_era5"] is not None]
    x = [d["date"] for d in shown]
    top.plot(x, [d["theta_top"] for d in shown], color=SERIES[0], linewidth=1.5, label=f"This model, 0-{TOP_CM:g} cm")
    top.plot(x, [d["sm_era5"] for d in shown], color=SERIES[1], linewidth=1.5, label="ERA5, 0-7 cm")
    top.set_ylabel("Water content (m³/m³)")
    top.set_title(f"Daily topsoil water content, rainfed field, {min(show_years)}-{max(show_years)}", fontsize=10)
    style(top)
    _legend(top, top.get_legend_handles_labels()[0], loc="upper right")
    months = np.arange(1, 13)
    model_clim, era5_clim = monthly_climatology(days)
    bottom.plot(months, model_clim, color=SERIES[0], linewidth=2, label=f"This model, 0-{TOP_CM:g} cm")
    bottom.plot(months, era5_clim, color=SERIES[1], linewidth=2, label="ERA5, 0-7 cm")
    bottom.set_xticks(months, ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"])
    bottom.set_ylabel("Water content (m³/m³)")
    bottom.set_title("Mean seasonal cycle, all years", fontsize=10)
    style(bottom)
    _legend(bottom, bottom.get_legend_handles_labels()[0], loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def monthly_climatology(days):
    month = np.array([d["date"].month for d in days if d["sm_era5"] is not None])
    model = np.array([d["theta_top"] for d in days if d["sm_era5"] is not None])
    era5 = np.array([d["sm_era5"] for d in days if d["sm_era5"] is not None])
    return (np.array([model[month == m].mean() for m in range(1, 13)]),
            np.array([era5[month == m].mean() for m in range(1, 13)]))


# ----------------------------------------------------------------- report

def _row(label, values, fmt="{:.0f}"):
    return f"| {label} | " + " | ".join(fmt.format(v) for v in values) + " |"


def forecast_section(records, eto, crop, snapshot, start, plain_days):
    scenarios = {"No forecast": [d for d in plain_days if d["date"] >= records[start]["timestamp"].date()]}
    scenarios.update(forecast_scenarios(records, eto, crop, snapshot, start))
    rain_by_date = {d["timestamp"].date(): d["precip"] for d in records}
    summary = {name: seasons(days) for name, days in scenarios.items()}
    years = sorted(summary["No forecast"])
    names = list(scenarios)
    first, last = years[0], years[-1]

    def mean(name, key):
        return float(np.mean([summary[name][y][key] for y in years]))

    lines = [
        "",
        f"## 6. Irrigation with rain forecasts ({first}-{last})",
        "",
        "Same irrigated field and trigger, but a due irrigation is postponed when at least 10 mm of "
        "rain is forecast for the next three days, unless 75 % of the available water is already "
        "used. Forecasts are issued the evening before; \"perfect\" knows the rain that fell and "
        "is the upper limit of what a forecast can save. All scenarios start from the same soil "
        f"state on 1 January {first}. Season means:",
        "",
        "| | " + " | ".join(names) + " |",
        "|---|" + "---|" * len(names),
        "| Irrigation (mm) | " + " | ".join(f"{mean(n, 'irrigation'):.0f}" for n in names) + " |",
        "| Number of irrigations | " + " | ".join(f"{mean(n, 'irrigations'):.1f}" for n in names) + " |",
        "| Postponed irrigation days | " + " | ".join(f"{mean(n, 'postponed'):.1f}" for n in names) + " |",
        "| Actual / potential transpiration | " + " | ".join(
            f"{mean(n, 'ta') / mean(n, 'tp'):.3f}" for n in names) + " |",
        f"| Drainage below {COLUMN['nz'] * COLUMN['dz'] / 100:g} m (mm) | " + " | ".join(
            f"{mean(n, 'drainage'):.1f}" for n in names) + " |",
        f"| Irrigations followed by ≥ 10 mm of rain within 3 days, all {len(years)} seasons | " + " | ".join(
            str(irrigated_before_rain(scenarios[n], rain_by_date)) for n in names) + " |",
        "",
        "Rain of 10 mm or more within three days is rare in this season (see the forecast skill "
        "report), which bounds what forecast-based postponing can save at this site.",
    ]
    return lines


def main(argv=None):
    parser = argparse.ArgumentParser(description="Maize water balance at the data site.")
    parser.add_argument("--data", default=str(weather_data.DEFAULT_CSV))
    parser.add_argument("--first-year", type=int, default=None)
    parser.add_argument("--last-year", type=int, default=None)
    parser.add_argument("--out", default=str(REPORT_DIR))
    args = parser.parse_args(argv)

    meta, records = weather_data.load(args.data)
    crop = adjust_for_climate(Crop(), records)
    if args.first_year:
        records = [d for d in records if d["timestamp"].year >= args.first_year]
    if args.last_year:
        records = [d for d in records if d["timestamp"].year <= args.last_year]
    eto = [asce_pm(d) for d in records]

    days, by_scenario, theirs, snapshots, field = {}, {}, {}, {}, None
    for name, irrigation in SCENARIOS.items():
        field, days[name], profiles, snapshots[name] = run_field(records, eto, crop, irrigation, SPLIT)
        by_scenario[name] = seasons(days[name])
        theirs[name] = run_pyfao56(records, eto, crop, field, profiles, irrigation)
    for d, rec in zip(days["rainfed"], records):
        d["sm_era5"] = rec.get("SM_0_7cm")
    years = sorted(by_scenario["irrigated"])
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    def mean(name, key):
        return float(np.mean([by_scenario[name][y][key] for y in years]))

    first, last = records[0]["timestamp"].date(), records[-1]["timestamp"].date()
    lines = [
        "# Crop water balance: grain maize",
        "",
        f"Generated by `python -m fao_model.water_balance` on {date.today().isoformat()}.",
        "",
        f"- Weather: {meta.get('site')}, {first} to {last} ({len(years)} seasons); reference ET from "
        "this repo's daily ASCE-PM.",
        f"- Crop: {crop.name} planted on {crop.planting_date(2001).strftime('%d %B')}, stages "
        f"{crop.l_ini}/{crop.l_dev}/{crop.l_mid}/{crop.l_late} days, Kcb {crop.kcb_ini}/"
        f"{crop.kcb_mid:.3f}/{crop.kcb_end} (mid-season value adjusted to the site's wind and humidity, "
        f"FAO-56 eq. 70), roots {crop.zr_ini}-{crop.zr_max} m, height up to {crop.h_max} m.",
        f"- Soil: loam (van Genuchten), {COLUMN['nz'] * COLUMN['dz'] / 100:g} m column in "
        f"{COLUMN['dz']:g} cm layers, free drainage at the bottom. Field capacity at "
        f"{H_FIELD_CAPACITY:g} cm (θ = {field.theta_fc:.3f}), wilting point θ = {field.theta_wp:.3f}.",
        "- ET: FAO-56 dual crop coefficients give potential transpiration (Kcb·ETo) and soil "
        "evaporation (Ke·ETo); the Richards solver limits both: roots take water through the "
        "Feddes stress function, evaporation stops as the surface dries.",
        f"- Scenarios: **rainfed**, and **irrigated** whenever the root-zone depletion reaches "
        f"{crop.p:.0%} of the total available water, refilling the root zone to field capacity.",
        "",
        "## 1. Mean seasonal water balance",
        "",
        "Season totals in mm (planting to harvest), averaged over all seasons. The first season "
        "starts from an assumed soil profile, the others from the simulated one.",
        "",
        "| | Rainfed | Irrigated |",
        "|---|---|---|",
    ]
    labels = {
        "precip": "Rain", "irrigation": "Irrigation", "eto": "Reference ET (ETo)",
        "tp": "Potential transpiration", "ta": "Actual transpiration",
        "ep": "Potential soil evaporation", "ea": "Actual soil evaporation",
        "drainage": f"Drainage below {COLUMN['nz'] * COLUMN['dz'] / 100:g} m", "runoff": "Runoff",
    }
    for key, label in labels.items():
        lines.append(_row(label, [mean("rainfed", key), mean("irrigated", key)]))
    lines.append(_row("Number of irrigations", [mean("rainfed", "irrigations"), mean("irrigated", "irrigations")], "{:.1f}"))
    lines.append(_row("Actual / potential transpiration",
                      [mean("rainfed", "ta") / mean("rainfed", "tp"),
                       mean("irrigated", "ta") / mean("irrigated", "tp")], "{:.2f}"))
    lines += [
        "",
        "## 2. Season by season",
        "",
        "| Season | Rain | ETo | Rainfed T/Tp | Irrigation | Irrigations | Irrigated T/Tp | Irrigated T + E |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for y in years:
        r, i = by_scenario["rainfed"][y], by_scenario["irrigated"][y]
        lines.append(f"| {y} | {i['precip']:.0f} | {i['eto']:.0f} | {r['ta'] / r['tp']:.2f} | "
                     f"{i['irrigation']:.0f} | {i['irrigations']:.0f} | {i['ta'] / i['tp']:.2f} | "
                     f"{i['ta'] + i['ea']:.0f} |")

    irrigation_by_year = {y: by_scenario["irrigated"][y]["irrigation"] for y in years}
    median_year = sorted(years, key=lambda y: irrigation_by_year[y])[len(years) // 2]
    in_season = lambda name: [d for d in days[name] if d["t"] is not None and d["date"].year == median_year]
    example_season_chart(out / "example_season.png", median_year, in_season("irrigated"), in_season("rainfed"))
    seasons_chart(out / "seasons.png", by_scenario)
    lines += [
        "",
        f"## 3. An example season ({median_year}, median irrigation)",
        "",
        "![Example season](example_season.png)",
        "",
        "![Seasonal irrigation and water stress](seasons.png)",
    ]

    lines += ["", "## 4. Comparison with pyfao56", ""]
    if theirs["irrigated"] is None:
        lines.append("pyfao56 is not installed (`pip install pyfao56`), so this section was skipped.")
    else:
        common = sorted(theirs["irrigated"])
        lines += [
            "pyfao56 implements the FAO-56 water balance as a bucket model: a two-stage "
            "evaporation layer (Ze = 0.10 m, REW = 9 mm) and a root zone whose transpiration falls "
            "linearly once depletion passes p·TAW. It was given the same weather, ETo, crop "
            "parameters, field capacity and wilting point, the same irrigation trigger, and on "
            "each planting day this model's soil water profile (5 cm layers). Season means in mm:",
            "",
            "| | Rainfed, this model | Rainfed, pyfao56 | Irrigated, this model | Irrigated, pyfao56 |",
            "|---|---|---|---|---|",
        ]
        for key, label in (("ta", "Transpiration"), ("ea", "Soil evaporation"),
                           ("irrigation", "Irrigation"), ("irrigations", "Number of irrigations"),
                           ("drainage", "Drainage / deep percolation")):
            vals = []
            for name in ("rainfed", "irrigated"):
                vals.append(np.mean([by_scenario[name][y][key] for y in common]))
                vals.append(np.mean([theirs[name][y][key] for y in common]))
            lines.append(_row(label, vals, "{:.1f}" if key == "irrigations" else "{:.0f}"))
        refill = crop.p * (field.theta_fc - field.theta_wp) * crop.zr_max * 1000.0
        lines += [
            "",
            "Season by season (this model − pyfao56):",
            "",
            "| | Bias (mm) | Relative bias (%) | r |",
            "|---|---|---|---|",
        ]
        pairs = [(f"Actual ET, {name}",
                  [by_scenario[name][y]["ta"] + by_scenario[name][y]["ea"] for y in common],
                  [theirs[name][y]["ta"] + theirs[name][y]["ea"] for y in common]) for name in SCENARIOS]
        pairs.append(("Irrigation, irrigated",
                      [by_scenario["irrigated"][y]["irrigation"] for y in common],
                      [theirs["irrigated"][y]["irrigation"] for y in common]))
        for label, a, b in pairs:
            m = compare(a, b)
            lines.append(f"| {label} | {m['bias']:+.0f} | {m['rel_bias']:+.1f} | {m['r']:.2f} |")
        lines += [
            "",
            "![Comparison with pyfao56](pyfao56.png)",
            "",
            "The season means agree closely. Single seasons differ more in irrigation because it "
            f"comes as discrete refills of up to about {refill:.0f} mm (55 % of the available water "
            "at full rooting depth): one refill more or less before harvest moves a season's total "
            "by that much. The remaining differences come from the model structure:",
            "",
            "- irrigation here refills the root zone to field capacity, while pyfao56 also adds "
            "the day's estimated ET;",
            "- transpiration stress here follows the Feddes function of the local pressure head, "
            "in pyfao56 a linear function of root-zone depletion;",
            "- here water moves between layers by Richards flow and drains below the roots after "
            "each refill; pyfao56 drains instantly to field capacity.",
            "",
            "Soil evaporation agrees with the 2.5 cm top layer used here; 5 cm layers gave about "
            "30 % more.",
        ]
        pyfao56_chart(out / "pyfao56.png", by_scenario, theirs)

    with_era5 = [d for d in days["rainfed"] if d["sm_era5"] is not None]
    m = compare([d["theta_top"] for d in with_era5], [d["sm_era5"] for d in with_era5])
    model_clim, era5_clim = monthly_climatology(days["rainfed"])
    show = tuple(range(max(first.year, last.year - 2), last.year + 1))
    era5_chart(out / "soil_moisture_era5.png", days["rainfed"], show)
    lines += [
        "",
        "## 5. Topsoil water content against ERA5",
        "",
        f"The rainfed field's water content in the top {TOP_CM:g} cm against ERA5's 0-7 cm soil "
        "moisture (a land-surface model with its own soil and vegetation, no irrigation). The "
        "soils differ, so the timing and the seasonal cycle matter more than the absolute level.",
        "",
        "| Daily values | Bias (m³/m³) | RMSE (m³/m³) | r | r of monthly means |",
        "|---|---|---|---|---|",
        f"| {m['n']} days | {m['bias']:+.3f} | {m['rmse']:.3f} | {m['r']:.2f} | "
        f"{np.corrcoef(model_clim, era5_clim)[0, 1]:.2f} |",
        "",
        "![Topsoil water content against ERA5](soil_moisture_era5.png)",
        "",
    ]
    start = next((i for i, d in enumerate(records) if d["timestamp"] >= SPLIT), 0)
    if snapshots["irrigated"] is not None and start >= 365:
        lines += forecast_section(records, eto, crop, snapshots["irrigated"], start, days["irrigated"])
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out / 'README.md'} and figures")


if __name__ == "__main__":
    main()
