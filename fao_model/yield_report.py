"""Yield response of grain maize to water at the data site (FAO-33 Ky).

Runs the rainfed and irrigated fields of the water balance report through
every season and turns their water use into relative yields with the FAO-33
yield response factors (fao_model.yield_response). Writes a Markdown report
with tables and figures:

    python -m fao_model.yield_report [--data CSV] [--first-year Y] [--last-year Y]
                                     [--max-yield T] [--out DIR]
"""
import argparse
from collections import defaultdict
from dataclasses import replace
from datetime import date
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator

from . import data as weather_data
from . import yield_response as yr
from .crop import Crop, adjust_for_climate
from .et import asce_pm
from .field import FieldModel
from .plotting import CATEGORICAL, INK_2, MUTED, SERIES, SURFACE, style
from .water_balance import SCENARIOS, run_field

REPORT_DIR = Path(__file__).resolve().parent.parent / "reports" / "yield"
MAX_YIELD = 13.0    # t/ha, an assumption; see the report
# Konya province, grain maize in 2024: 2 011 729 t from 1 589 939 da (TÜİK data, in
# Konya Provincial Directorate of Agriculture and Forestry, "Konya Tarım İstatistikleri", 2025).
KONYA_2024 = 12.65  # t/ha
KONYA_SOURCE = "https://konya.tarimorman.gov.tr/Belgeler/kitap/konyatarimistilce2025.pdf"
DEFAULT_METHOD = "By growth period (FAO-33)"


def relative_yields(days, crop):
    """Ya/Yx of one season by each method compared in the report."""
    return {
        DEFAULT_METHOD: yr.stage_yield(yr.water_use(days, crop)),
        f"Whole season, Ky = {yr.SEASON_KY}": yr.season_yield(yr.water_use(days, crop)),
        "Periods on the FAO-56 stages (FAO-66 Table 2)": yr.stage_yield(
            yr.water_use(days, crop, yr.fao56_stage), yr.MAIZE_KY_BY_FAO56_STAGE),
        "By growth period, transpiration only": yr.stage_yield(
            yr.water_use(days, crop, basis="transpiration")),
    }


def by_season(days):
    """In-season field days grouped by year."""
    result = defaultdict(list)
    for d in days:
        if d["t"] is not None:
            result[d["date"].year].append(d)
    return dict(result)


def period_days(crop):
    """{period: (first day, last day)} of the season."""
    spans = {}
    for t in range(1, crop.season_length + 1):
        name = yr.growth_period(crop, t)
        first, _ = spans.get(name, (t, t))
        spans[name] = (first, t)
    return spans


def stress_after_irrigation(seasons):
    """Share of the season days that are irrigation days or the day after, their share of
    the transpiration deficit, and Ta/Tp on those days and on the others."""
    near, deficit, ta, tp = [], [], [], []
    for days in seasons:
        irrigated = [d["irrigation"] > 0 for d in days]
        for i, d in enumerate(days):
            near.append(irrigated[i] or (i > 0 and irrigated[i - 1]))
            deficit.append(d["tp"] - d["ta"])
            ta.append(d["ta"])
            tp.append(d["tp"])
    near, deficit, ta, tp = map(np.asarray, (near, deficit, ta, tp))
    return {
        "day_share": near.mean(),
        "deficit_share": deficit[near].sum() / deficit.sum(),
        "ratio_near": ta[near].sum() / tp[near].sum(),
        "ratio_other": ta[~near].sum() / tp[~near].sum(),
    }


def _legend(ax, **kwargs):
    ax.legend(frameon=False, labelcolor=INK_2, fontsize=8, **kwargs)


def seasons_chart(path, years, by_method):
    fig = Figure(figsize=(10, 4.4), facecolor=SURFACE)
    ax = fig.subplots()
    lowest = 1.0
    for (label, values), color in zip(by_method.items(), CATEGORICAL[:len(by_method)], strict=True):
        ax.plot(years, values, color=color, linewidth=2, marker="o", markersize=5, label=label)
        lowest = min(lowest, min(values))
    ax.set_ylim(max(0.0, lowest - 0.05), 1.0)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.ticklabel_format(axis="x", useOffset=False, style="plain")
    ax.set_ylabel("Relative yield Ya / Yx")
    ax.set_title("Irrigated maize: relative yield by season and method", fontsize=10)
    style(ax)
    _legend(ax, loc="lower left", ncols=2)
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def production_function_chart(path, points):
    """Season ET deficit against yield loss, with the whole-season relation for reference."""
    fig = Figure(figsize=(6.6, 5.2), facecolor=SURFACE)
    ax = fig.subplots()
    x = np.linspace(0.0, 1.0, 101)
    ax.plot(x, np.minimum(yr.SEASON_KY * x, 1.0), color=MUTED, linewidth=1.5,
            label=f"Whole season, Ky = {yr.SEASON_KY}")
    for (label, (deficits, losses)), color in zip(points.items(), SERIES[:len(points)], strict=True):
        ax.scatter(deficits, losses, s=40, color=color, edgecolors=SURFACE, linewidths=1.5, zorder=3, label=label)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel("Season ET deficit, 1 − ETa / ETx")
    ax.set_ylabel("Yield loss, 1 − Ya / Yx")
    ax.set_title("Yield loss against ET deficit, one point per season", fontsize=10)
    style(ax)
    _legend(ax, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Yield response of maize to water at the data site.")
    parser.add_argument("--data", default=str(weather_data.DEFAULT_CSV))
    parser.add_argument("--first-year", type=int, default=None)
    parser.add_argument("--last-year", type=int, default=None)
    parser.add_argument("--max-yield", type=float, default=MAX_YIELD, help="maximum yield Yx, t/ha")
    parser.add_argument("--out", default=str(REPORT_DIR))
    args = parser.parse_args(argv)

    meta, records = weather_data.load(args.data)
    crop = adjust_for_climate(Crop(), records)
    if args.first_year:
        records = [d for d in records if d["timestamp"].year >= args.first_year]
    if args.last_year:
        records = [d for d in records if d["timestamp"].year <= args.last_year]
    eto = [asce_pm(d) for d in records]
    seasons = {}
    for name, irrigation in SCENARIOS.items():
        _, days, _, _ = run_field(records, eto, crop, irrigation)
        seasons[name] = {y: s for y, s in by_season(days).items() if len(s) == crop.season_length}
    # The irrigated field again without the wet end of the Feddes function (oxygen stress).
    no_oxygen = replace(crop, feddes_h1=0.0, feddes_h2=0.0)
    _, days, _, _ = run_field(records, eto, no_oxygen, True)
    seasons["no oxygen stress"] = {y: s for y, s in by_season(days).items() if len(s) == crop.season_length}
    years = sorted(seasons["irrigated"])
    ry = {name: {y: relative_yields(seasons[name][y], crop) for y in years} for name in SCENARIOS}
    totals = {name: {y: yr.water_use(seasons[name][y], crop) for y in years} for name in SCENARIOS}

    def season_sum(name, y, key):
        return sum(d[key] for d in seasons[name][y])

    def et_ratio(name, y):
        return (sum(a for a, _ in totals[name][y].values())
                / sum(x for _, x in totals[name][y].values()))

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    spans = period_days(crop)
    field = FieldModel(crop)
    first, last = records[0]["timestamp"].date(), records[-1]["timestamp"].date()
    lines = [
        "# Yield response to water: grain maize",
        "",
        f"Generated by `python -m fao_model.yield_report` on {date.today().isoformat()}.",
        "",
        f"- Weather: {meta.get('site')}, {first} to {last} ({len(years)} seasons).",
        "- Fields: the rainfed and irrigated maize fields of the "
        "[water balance report](../water_balance/README.md) (FAO-56 dual crop coefficients over the "
        f"Richards solver; irrigation refills the root zone at {crop.p:.0%} depletion). The first "
        "season starts from an assumed soil profile, the others from the simulated one.",
        "- Yield: the FAO-33 water production function (Doorenbos and Kassam, 1979; Steduto et al., "
        "2012, FAO-66), 1 − Ya/Yx = Ky (1 − ETa/ETx). ETx is the crop ET with the crop's water needs "
        "met; with dual crop coefficients the deficit is in transpiration, so soil evaporation counts "
        "the same in both: ETx = Tp + Ea, ETa = Ta + Ea.",
        "- Growth periods and their Ky (FAO-33 values for maize), combined by multiplying the "
        "periods' relative yields (Rao et al., 1988):",
        "",
        "| Period | Days of the season | Ky | Placed on |",
        "|---|---|---|---|",
    ]
    placed = {"vegetative": "initial and development stages",
              "flowering": f"first {yr.FLOWERING_DAYS} days of mid-season (FAO-56: full cover at the "
                           "initiation of flowering); the length is assumed",
              "yield formation": "rest of mid-season", "ripening": "late season"}
    for period, (a, b) in spans.items():
        lines.append(f"| {period.capitalize()} | {a}-{b} | {yr.MAIZE_KY[period]:.2f} | {placed[period]} |")
    lines += [
        "",
        f"The whole-season Ky is {yr.SEASON_KY}. FAO-66 warns that yields estimated this way carry "
        "large errors, so the results are best read as relative yield losses.",
        "",
        "## 1. Relative yield",
        "",
        "Mean Ya/Yx over the seasons by method (lowest and highest season in brackets):",
        "",
        "| Method | Rainfed | Irrigated |",
        "|---|---|---|",
    ]
    for method in ry["irrigated"][years[0]]:
        cells = []
        for name in ("rainfed", "irrigated"):
            values = [ry[name][y][method] for y in years]
            cells.append(f"{np.mean(values):.2f} ({min(values):.2f}-{max(values):.2f})")
        lines.append(f"| {method} | " + " | ".join(cells) + " |")
    rain = np.mean([season_sum("rainfed", y, "precip") for y in years])
    eto_season = np.mean([season_sum("rainfed", y, "eto") for y in years])
    mean_default = {name: np.mean([ry[name][y][DEFAULT_METHOD] for y in years]) for name in SCENARIOS}
    mean_irrigated = {m: np.mean([ry["irrigated"][y][m] for y in years]) for m in ry["irrigated"][years[0]]}
    factors = {p: np.mean([yr.period_yield(*totals["irrigated"][y][p], yr.MAIZE_KY[p]) for y in years])
               for p in spans}
    worst = min(factors, key=factors.get)
    season_method, stage_method = list(mean_irrigated)[1:3]
    lines += [
        "",
        f"- **Rainfed maize fails at this site:** about {rain:.0f} mm of rain falls in the season "
        f"against {eto_season:.0f} mm of reference ET, and the relative yield is "
        f"{mean_default['rainfed']:.2f}.",
        f"- **The irrigated field still loses {1 - mean_default['irrigated']:.0%}** by growth period "
        f"and {1 - mean_irrigated[season_method]:.0%} with the whole-season Ky. The largest loss is in "
        f"the {worst} period (Ky {yr.MAIZE_KY[worst]:.2f}, section 2).",
        f"- **The placement of the periods matters:** with Ky placed on the FAO-56 stages, which puts "
        f"the flowering value on the development stage, the mean is "
        f"{mean_irrigated[stage_method]:.2f} instead of {mean_default['irrigated']:.2f}.",
        "",
        "![Relative yield by season](relative_yield.png)",
        "",
        "![Yield loss against ET deficit](production_function.png)",
        "",
        "Points above the whole-season line lost more yield than their seasonal ET deficit "
        "suggests, because the deficit fell on sensitive periods.",
        "",
        "## 2. Where the irrigated field loses yield",
        "",
        "| Period | Ky | ETa / ETx | Relative yield of the period |",
        "|---|---|---|---|",
    ]
    for period in spans:
        ratios = [totals["irrigated"][y][period][0] / totals["irrigated"][y][period][1] for y in years]
        lines.append(f"| {period.capitalize()} | {yr.MAIZE_KY[period]:.2f} | {np.mean(ratios):.3f} | "
                     f"{factors[period]:.3f} |")
    s = stress_after_irrigation([seasons["irrigated"][y] for y in years])
    h3_depletion = (field.theta_fc - float(field.soil._theta(np.array(crop.feddes_h3_high)))) / (
        field.theta_fc - field.theta_wp)
    lines += [
        "",
        "Why a field irrigated whenever it needs water is still stressed:",
        "",
        f"- **Just after irrigation.** Irrigation days and the day after make up {s['day_share']:.0%} "
        f"of the season but {s['deficit_share']:.0%} of the transpiration deficit: Ta/Tp is "
        f"{s['ratio_near']:.2f} on those days against {s['ratio_other']:.2f} on the others. The "
        "whole refill is applied in one day, which wets the top layers above "
        f"{crop.feddes_h2:g} cm, where the Feddes function reduces uptake for lack of oxygen.",
        f"- **Before irrigation.** At high demand, uptake starts to fall once the head drops below "
        f"{crop.feddes_h3_high:g} cm, which in this soil is {h3_depletion:.0%} depletion, close to the "
        f"{crop.p:.0%} trigger. The "
        "top of the root zone, where most roots are, dries faster than the root-zone average, and "
        "uptake lost there is not made up from wetter, deeper layers.",
        "",
        "The irrigated field run again without oxygen stress (Feddes h1 = h2 = 0 cm) shows how much "
        "of the loss comes from it:",
        "",
        "| Irrigated field | Ta / Tp | Ta / Tp on irrigation days and the day after | Irrigation (mm) | Ya / Yx |",
        "|---|---|---|---|---|",
    ]
    for label, name in (("As modelled", "irrigated"), ("Without oxygen stress", "no oxygen stress")):
        runs = [seasons[name][y] for y in years]
        stress = stress_after_irrigation(runs)
        ratio = sum(d["ta"] for r in runs for d in r) / sum(d["tp"] for r in runs for d in r)
        relative = np.mean([yr.stage_yield(yr.water_use(r, crop)) for r in runs])
        lines.append(f"| {label} | {ratio:.3f} | {stress['ratio_near']:.2f} | "
                     f"{np.mean([season_sum(name, y, 'irrigation') for y in years]):.0f} | {relative:.2f} |")
    lines += [
        "",
        "Both causes are properties of the field model rather than of the yield response: a daily time "
        "step with the whole refill on one day, and root water uptake without compensation.",
        "",
        "## 3. Yield and water productivity",
        "",
        f"Absolute yields need the maximum yield Yx, assumed here to be {args.max_yield:g} t/ha. For "
        f"comparison, grain maize in Konya province averaged {KONYA_2024:g} t/ha in 2024 (TÜİK data in "
        f"[Konya Tarım İstatistikleri]({KONYA_SOURCE}), Konya Provincial Directorate of Agriculture "
        "and Forestry, 2025). All yields and water productivities below scale with Yx.",
        "",
        "| Season | Rain (mm) | Irrigation (mm) | ETa / ETx | Ya / Yx | Yield (t/ha) | "
        "ET water productivity (kg/m³) | Irrigation water productivity (kg/m³) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    wp_et, wp_irr, yields = [], [], []
    for y in years:
        ya = args.max_yield * ry["irrigated"][y][DEFAULT_METHOD]
        ya_rainfed = args.max_yield * ry["rainfed"][y][DEFAULT_METHOD]
        eta = season_sum("irrigated", y, "ta") + season_sum("irrigated", y, "ea")
        irrigation = season_sum("irrigated", y, "irrigation")
        # t/ha over mm: 1 t/ha = 1000 kg/ha and 1 mm = 10 m3/ha.
        wp_et.append(100.0 * ya / eta)
        wp_irr.append(100.0 * (ya - ya_rainfed) / irrigation)
        yields.append(ya)
        lines.append(f"| {y} | {season_sum('irrigated', y, 'precip'):.0f} | {irrigation:.0f} | "
                     f"{et_ratio('irrigated', y):.3f} | {ry['irrigated'][y][DEFAULT_METHOD]:.2f} | "
                     f"{ya:.1f} | {wp_et[-1]:.2f} | {wp_irr[-1]:.2f} |")
    lines += [
        f"| Mean | | | | {mean_default['irrigated']:.2f} | {np.mean(yields):.1f} | {np.mean(wp_et):.2f} | "
        f"{np.mean(wp_irr):.2f} |",
        "",
        "- ET water productivity: yield per cubic metre of actual crop ET (transpiration and soil "
        "evaporation).",
        "- Irrigation water productivity: the yield gained over the rainfed field per cubic metre of "
        "irrigation.",
        "",
    ]
    by_method = {m: [ry["irrigated"][y][m] for y in years] for m in ry["irrigated"][years[0]]}
    seasons_chart(out / "relative_yield.png", years, by_method)
    points = {f"{name.capitalize()} seasons": ([1.0 - et_ratio(name, y) for y in years],
                                               [1.0 - ry[name][y][DEFAULT_METHOD] for y in years])
              for name in ("irrigated", "rainfed")}
    production_function_chart(out / "production_function.png", points)
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out / 'README.md'} and figures")


if __name__ == "__main__":
    main()
