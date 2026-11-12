"""How the field model takes up and applies water: compensated root uptake and timed irrigation.

Three experiments behind the field model's settings, written to a Markdown
report with tables and figures:

1. a root zone at full cover dries from field capacity: where does
   transpiration start to fall, for several critical stress indices omega_c
   of compensated uptake (Jarvis, 1989), and where does FAO-56 put it?
2. one small irrigation on dry bare soil, spread over the whole day or
   applied at the field model's rate: where does the water go?
3. the irrigated maize field of the water balance report, every season,
   with each setting switched on and off.

    python -m fao_model.field_model_report [--data CSV] [--first-year Y] [--last-year Y] [--out DIR]
"""
import argparse
from dataclasses import replace
from datetime import date
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure

from . import data as weather_data
from . import yield_response as yr
from .crop import Crop, adjust_for_climate, root_weights
from .et import asce_pm
from .field import COLUMN, IRRIGATION_RATE, FieldModel
from .plotting import CATEGORICAL, INK_2, MUTED, SURFACE, style
from .richards import RichardsSolver3D
from .water_balance import run_field
from .yield_report import by_season, stress_after_irrigation

REPORT_DIR = Path(__file__).resolve().parent.parent / "reports" / "field_model"
OMEGAS = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5)
DEMANDS = (5.0, 8.0)        # potential transpiration, mm/day
PROFILE = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
SMALL_IRRIGATION = 14.0     # mm, a typical early-season refill
BARE_SOIL_EP = 6.0          # mm/day


def fao56_p(crop, tp):
    """FAO-56 depletion fraction adjusted for the demand (Table 22 footnote), 0.1-0.8."""
    return min(0.8, max(0.1, crop.p + 0.04 * (5.0 - tp)))


def dry_down(crop, omega_c, tp, days=150):
    """Root-zone depletion fraction and Ta/Tp, day by day, from field capacity without rain."""
    field = FieldModel(crop)
    zr = 100.0 * crop.zr_max
    weights = root_weights(zr, field.soil.dz, field.soil.nz)
    depletion, ratio = [], []
    for _ in range(days):
        before, taw = field.root_zone(zr)
        out = field.soil.step(0.0, dt=1.0, tp=tp / 10.0, root_weights=weights,
                              alpha=lambda h: crop.feddes_alpha(h, tp / 10.0), omega_c=omega_c)
        after, _ = field.root_zone(zr)
        depletion.append(0.5 * (before + after) / taw)
        ratio.append(10.0 * out["transpiration"] / tp)
        if depletion[-1] > 0.9:
            break
    return np.array(depletion), np.array(ratio)


def onset(depletion, ratio, level=0.95):
    """Depletion at which Ta/Tp first falls below `level`."""
    below = np.nonzero(ratio < level)[0]
    return float(depletion[below[0]]) if below.size else float("nan")


def one_irrigation(rate, field):
    """A small irrigation on uniformly dry bare soil, then six dry days.

    Returns the water gained by depth after the first day (mm per layer), the
    evaporation of each day (mm) and the water left in the soil at the end.
    """
    theta0 = field.theta_fc - field.crop.p * (field.theta_fc - field.theta_wp)
    soil = RichardsSolver3D(**COLUMN, h_init=float(field.soil._head(np.array(theta0))))
    start = soil.theta.mean(axis=(0, 1)).copy()
    q_day = -BARE_SOIL_EP / 10.0
    hours = 24.0 if rate is None else min(24.0, SMALL_IRRIGATION / rate)
    segments = [(hours / 24.0, q_day + SMALL_IRRIGATION / 10.0 * 24.0 / hours)]
    if hours < 24.0:
        segments.append((1.0 - hours / 24.0, q_day))
    evaporation, gain = [], None
    for day in range(7):
        out = soil.run(segments if day == 0 else [(1.0, q_day)])
        evaporation.append(BARE_SOIL_EP - 10.0 * out["evaporation_deficit"])
        if day == 0:
            gain = (soil.theta.mean(axis=(0, 1)) - start) * soil.dz * 10.0
    left = float((soil.theta.mean(axis=(0, 1)) - start).sum() * soil.dz * 10.0)
    return gain, evaporation, left


def _legend(ax, **kwargs):
    ax.legend(frameon=False, labelcolor=INK_2, fontsize=8, **kwargs)


def dry_down_chart(path, curves, crop, tp):
    fig = Figure(figsize=(7.2, 4.8), facecolor=SURFACE)
    ax = fig.subplots()
    x = np.linspace(0.0, 0.9, 181)
    p = fao56_p(crop, tp)
    ax.plot(x, np.clip((1.0 - x) / (1.0 - p), 0.0, 1.0), color=MUTED, linewidth=1.5,
            label=f"FAO-56 Ks (p = {p:.2f})")
    for (label, (depletion, ratio)), color in zip(curves.items(), CATEGORICAL[:len(curves)], strict=True):
        ax.plot(depletion, ratio, color=color, linewidth=2, label=label)
    ax.set_xlim(0.0, 0.9)
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel("Root-zone depletion, fraction of the available water")
    ax.set_ylabel("Actual / potential transpiration")
    ax.set_title(f"Full-cover maize drying from field capacity, Tp = {tp:g} mm/day", fontsize=10)
    style(ax)
    _legend(ax, loc="lower left")
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def wetting_chart(path, gains, dz):
    fig = Figure(figsize=(6.4, 4.8), facecolor=SURFACE)
    ax = fig.subplots()
    for (label, gain), color in zip(gains.items(), CATEGORICAL[:len(gains)], strict=True):
        depth = dz * (np.arange(len(gain)) + 0.5)
        shown = depth <= 40.0
        ax.plot(gain[shown], depth[shown], color=color, linewidth=2, marker="o", markersize=4, label=label)
    ax.invert_yaxis()
    ax.set_xlabel("Water gained by the end of the day (mm per 2.5 cm layer)")
    ax.set_ylabel("Depth (cm)")
    ax.set_title(f"{SMALL_IRRIGATION:g} mm of irrigation on dry bare soil", fontsize=10)
    style(ax)
    _legend(ax, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Root uptake and irrigation settings of the field model.")
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
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    field = FieldModel(crop)
    first, last = records[0]["timestamp"].date(), records[-1]["timestamp"].date()

    lines = [
        "# Field model: root water uptake and irrigation",
        "",
        f"Generated by `python -m fao_model.field_model_report` on {date.today().isoformat()}.",
        "",
        "Two settings of the field model decide how much a well-irrigated crop is stressed. This "
        "report shows the experiments behind them.",
        "",
        f"- **Compensated root water uptake** (Jarvis, 1989; Šimůnek and Hopmans, 2009), ωc = "
        f"{crop.omega_c:g}. Roots take water where the Feddes function allows it; ω, the "
        "root-weighted mean of the Feddes factor, measures the stress of the whole root zone. While "
        "ω ≥ ωc the roots in favourable layers make up for stressed ones and transpiration stays "
        "potential; below, it is ω/ωc of potential. ωc = 1 is uptake without compensation, as in "
        "versions 03-06.",
        f"- **Irrigation at {IRRIGATION_RATE:g} mm/hour** from the start of the day instead of spread "
        "over all of it.",
        "",
        "## 1. When does a drying root zone start to stress the crop?",
        "",
        "A full-cover root zone (1.2 m, loam) dries from field capacity without rain. The table gives "
        "the depletion at which transpiration falls below 95 % of potential, and Ta/Tp at given "
        "depletions. FAO-56 keeps transpiration potential until the depletion reaches p and then "
        "lowers it linearly to zero at the wilting point (Ks).",
        "",
    ]
    curves_shown, onsets = {}, {}
    for tp in DEMANDS:
        p = fao56_p(crop, tp)
        lines += [
            f"Potential transpiration {tp:g} mm/day (FAO-56 p = {p:.2f}):",
            "",
            "| ωc | Stress from depletion | " + " | ".join(f"Ta/Tp at {x:.1f}" for x in PROFILE) + " |",
            "|---|---|" + "---|" * len(PROFILE),
        ]
        for omega_c in OMEGAS:
            depletion, ratio = dry_down(crop, omega_c, tp)
            at = [np.interp(x, depletion, ratio) if x <= depletion[-1] else float("nan") for x in PROFILE]
            label = f"{omega_c:g}" + (" (no compensation)" if omega_c == 1.0 else "") + (
                " (used here)" if omega_c == crop.omega_c else "")
            onsets[tp, omega_c] = onset(depletion, ratio)
            lines.append(f"| {label} | {onsets[tp, omega_c]:.2f} | " + " | ".join(f"{v:.2f}" for v in at) + " |")
            if tp == DEMANDS[0] and omega_c in (1.0, crop.omega_c, 0.5):
                curves_shown[f"ωc = {omega_c:g}"] = (depletion, ratio)
        lines += ["| FAO-56 | " + f"{p:.2f} | " + " | ".join(
            f"{min(1.0, (1.0 - x) / (1.0 - p)):.2f}" for x in PROFILE) + " |", ""]
    dry_down_chart(out / "dry_down.png", curves_shown, crop, DEMANDS[0])
    lines += [
        "![Dry-down](dry_down.png)",
        "",
        f"- **Without compensation** transpiration falls below 95 % at {onsets[DEMANDS[0], 1.0]:.0%} "
        f"depletion ({DEMANDS[0]:g} mm/day), before the {crop.p:.0%} irrigation trigger: the upper root "
        "zone, where most roots are, dries first.",
        f"- **With ωc = {crop.omega_c:g}** it does so at {onsets[DEMANDS[0], crop.omega_c]:.0%}, against "
        f"a FAO-56 p of {fao56_p(crop, DEMANDS[0]):.2f} for that demand, so the Richards field and the "
        "FAO-56 bucket model agree on when stress begins. No value of ωc "
        "is tabulated for crops; it is chosen here for that agreement. Below the threshold the decline "
        "is steeper than FAO-56's straight line.",
        "",
    ]

    gains, rows = {}, []
    for label, rate in (("Over 24 hours", None), (f"At {IRRIGATION_RATE:g} mm/hour", IRRIGATION_RATE)):
        gain, evaporation, left = one_irrigation(rate, field)
        gains[label] = gain
        depth = field.soil.dz * (np.arange(field.soil.nz) + 1)
        front = depth[np.nonzero(gain > 0.05)[0].max()]
        rows.append(f"| {label} | {gain[:2].sum():.1f} | {front:g} | {evaporation[0]:.1f} | "
                    f"{sum(evaporation):.1f} | {left:+.1f} |")
    wetting_chart(out / "wetting.png", gains, field.soil.dz)
    lines += [
        "## 2. Where does a small irrigation go?",
        "",
        f"{SMALL_IRRIGATION:g} mm of irrigation, a typical early-season refill, on bare loam at "
        f"{crop.p:.0%} depletion, with {BARE_SOIL_EP:g} mm/day of potential evaporation that day and "
        "the six dry days after:",
        "",
        "| Irrigation | Gain in the top 5 cm after the first day (mm) | Wetting front (cm) | "
        "Evaporation, first day (mm) | Evaporation, 7 days (mm) | Left in the soil (mm) |",
        "|---|---|---|---|---|---|",
        *rows,
        "",
        "![Wetting after one irrigation](wetting.png)",
        "",
        "Spread over the whole day, irrigation behaves like a day-long drizzle: the water stays in the "
        "top few centimetres, where evaporation takes it back. Applied in a few hours, it soaks deeper "
        f"and more of it stays. {IRRIGATION_RATE:g} mm/hour is two thirds of this soil's saturated "
        f"conductivity ({field.soil.params['Ks'] * 10.0 / 24.0:.0f} mm/hour), so it infiltrates without "
        "runoff.",
        "",
    ]

    eto = [asce_pm(d) for d in records]
    variants = {
        "As in 06: all day, no compensation": (None, 1.0),
        f"{IRRIGATION_RATE:g} mm/hour, no compensation": (IRRIGATION_RATE, 1.0),
        f"All day, ωc = {crop.omega_c:g}": (None, crop.omega_c),
        f"{IRRIGATION_RATE:g} mm/hour, ωc = {crop.omega_c:g} (used now)": (IRRIGATION_RATE, crop.omega_c),
    }
    lines += [
        "## 3. The irrigated maize field",
        "",
        f"The irrigated field of the water balance report, {first} to {last}, with each setting "
        "switched on and off. Season means:",
        "",
        "| Field model | Ta/Tp | Ta/Tp on irrigation days and the day after | Irrigation (mm) | "
        "Irrigations | Soil evaporation (mm) | Drainage (mm) | Runoff (mm) | Ya/Yx |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for label, (rate, omega_c) in variants.items():
        variant = replace(crop, omega_c=omega_c)
        _, days, _, _ = run_field(records, eto, variant, True, irrigation_rate=rate)
        seasons = [s for s in by_season(days).values() if len(s) == crop.season_length]
        flat = [d for s in seasons for d in s]
        mean = lambda key: sum(d[key] for d in flat) / len(seasons)   # noqa: E731
        near = stress_after_irrigation(seasons)["ratio_near"]
        relative = np.mean([yr.stage_yield(yr.water_use(s, variant)) for s in seasons])
        lines.append(
            f"| {label} | {mean('ta') / mean('tp'):.3f} | {near:.2f} | {mean('irrigation'):.0f} | "
            f"{sum(d['irrigation'] > 0 for d in flat) / len(seasons):.1f} | {mean('ea'):.0f} | "
            f"{mean('drainage'):.0f} | {mean('runoff'):.0f} | {relative:.2f} |")
    lines += [
        "",
        "- **Compensation removes most of the stress between irrigations**, so the field needs a little "
        "more water but uses it for transpiration.",
        "- **Timed irrigation** wets the soil deeper and loses less to evaporation early in the season; "
        "the stress right after irrigation, the Feddes function's oxygen stress in the wet top layers, "
        "remains.",
        "- The yield report and the water balance report use the settings marked \"used now\".",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out / 'README.md'} and figures")


if __name__ == "__main__":
    main()
