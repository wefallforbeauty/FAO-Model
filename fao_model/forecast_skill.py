"""Skill of the daily weather and reference ET forecasts over the test period.

Forecasts are issued every `--every` days of 2015-2024 (rolling origin:
only data up to the issue day is used) for 1-14 days ahead and compared
with what happened. References: climatology (the smoothed day-of-year mean
of 1995-2014) and persistence (tomorrow equals today). TimesFM results are
cached in output/cache/, so the report can be rebuilt quickly.

    python -m fao_model.forecast_skill [--every 7] [--context 1024] [--out DIR]
"""
import argparse
from datetime import date, datetime
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure

from . import data as weather_data
from .forecast import CONTEXT, HORIZON, WEATHER, Climatology, WeatherForecaster, add_eto
from .plotting import BLUE_LIGHT, INK, INK_2, MUTED, SERIES, SURFACE, style

ROOT = Path(__file__).resolve().parent.parent
REPORT_DIR = ROOT / "reports" / "forecast_skill"
CACHE_DIR = ROOT / "output" / "cache"
SPLIT = datetime(2015, 1, 1)
LEADS = (1, 3, 7, 14)
RAIN_DAYS, RAIN_EVENT = 3, 10.0     # rain events: at least 10 mm within 3 days
LABELS = {
    "T_max": "Max temperature (°C)", "T_min": "Min temperature (°C)",
    "RH_max": "Max humidity (%)", "RH_min": "Min humidity (%)", "u2": "Wind (m/s)",
    "Rs": "Radiation (MJ/m²/day)", "precip": "Rain (mm/day)", "ETo": "Reference ET (mm/day)",
}


def cached(name, compute):
    """Load arrays saved under output/cache/<name>.npz, or compute and save them."""
    path = CACHE_DIR / f"{name}.npz"
    if path.exists():
        with np.load(path) as f:
            return {k: f[k] for k in f.files}
    result = compute()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **result)
    return result


def issue_forecasts(forecaster, records, origins, chunk=26):
    """forecast_at in chunks, printing progress (TimesFM takes minutes)."""
    parts = []
    for start in range(0, len(origins), chunk):
        parts.append(forecaster.forecast_at(records, origins[start:start + chunk]))
        print(f"  {min(start + chunk, len(origins))}/{len(origins)} forecasts", flush=True)
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def mae_by_lead(forecast, observed):
    return np.abs(forecast - observed).mean(axis=0)


# ---------------------------------------------------------------- figures

def eto_chart(path, curves, climatology):
    fig = Figure(figsize=(8, 4.4), facecolor=SURFACE)
    ax = fig.subplots()
    leads = np.arange(1, HORIZON + 1)
    ax.plot(leads, climatology, color=MUTED, linewidth=2, label="Climatology")
    for (label, values), color in zip(curves, SERIES + ("#1baf7a", "#eda100")):
        ax.plot(leads, values, color=color, linewidth=2, marker="o", markersize=4, label=label)
    ax.set_xlabel("Days ahead")
    ax.set_ylabel("Mean absolute error (mm/day)")
    ax.set_title("Reference ET forecast error, 2015-2024", fontsize=10)
    ax.set_xticks([1, 3, 5, 7, 9, 11, 13])
    ax.set_ylim(bottom=0)
    style(ax)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_2, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def skill_chart(path, panels):
    """Small multiples: skill score against climatology by lead, two methods per panel."""
    fig = Figure(figsize=(13, 6.4), facecolor=SURFACE)
    axes = fig.subplots(2, 4, sharex=True, sharey=True)
    leads = np.arange(1, HORIZON + 1)
    for ax, (key, curves) in zip(axes.flat, panels):
        ax.axhline(0.0, color=INK_2, linewidth=1)
        for (label, values), color in zip(curves, SERIES):
            ax.plot(leads, values, color=color, linewidth=2, label=label)
        ax.set_title(LABELS[key], fontsize=9)
        ax.set_xticks([1, 7, 14])
        style(ax)
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    axes[0, 0].set_ylim(-0.15, 0.55)
    fig.supxlabel("Days ahead", color=INK_2, fontsize=9)
    fig.supylabel("Skill against climatology (1 − MAE / MAE of climatology)", color=INK_2, fontsize=9)
    fig.suptitle("Where forecasts beat climatology (above 0), 2015-2024", color=INK, fontsize=11)
    fig.tight_layout()
    # The legend goes where the unused last panel would be.
    handles, labels = axes.flat[0].get_legend_handles_labels()
    box = axes.flat[-1].get_position()
    fig.legend(handles, labels, loc="center", bbox_to_anchor=(box.x0 + box.width / 2, box.y0 + box.height / 2),
               frameon=False, fontsize=9, labelcolor=INK_2)
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def examples_chart(path, records, examples):
    """A few forecasts of reference ET with their 10-90 % band against what happened."""
    fig = Figure(figsize=(13, 3.8), facecolor=SURFACE)
    axes = fig.subplots(1, len(examples), sharey=True)
    for ax, (i, fc) in zip(axes, examples):
        past = np.arange(-20, 1)
        ahead = np.arange(1, HORIZON + 1)
        ax.plot(past, [records[i + k]["ETo"] for k in past], color=INK_2, linewidth=1.5, label="Observed")
        ax.plot(ahead, [records[i + k]["ETo"] for k in ahead], color=INK_2, linewidth=1.5, linestyle=(0, (1, 2)),
                label="Observed later")
        if "ETo_q10" in fc:
            ax.fill_between(ahead, fc["ETo_q10"], fc["ETo_q90"], color=BLUE_LIGHT, alpha=0.5, linewidth=0,
                            label="TimesFM 10-90 %")
        ax.plot(ahead, fc["ETo"], color=SERIES[0], linewidth=2, label="Forecast (from forecast weather)")
        ax.axvline(0.5, color=MUTED, linewidth=1)
        ax.set_title(records[i]["timestamp"].strftime("Issued %d %b %Y"), fontsize=9)
        ax.set_xlabel("Days from issue")
        style(ax)
    axes[0].set_ylabel("Reference ET (mm/day)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=8, labelcolor=INK_2)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(path, dpi=110, facecolor=SURFACE)


# ----------------------------------------------------------------- report

def main(argv=None):
    parser = argparse.ArgumentParser(description="Skill of the weather and ETo forecasts.")
    parser.add_argument("--data", default=str(weather_data.DEFAULT_CSV))
    parser.add_argument("--every", type=int, default=7, help="days between forecast issues")
    parser.add_argument("--context", type=int, default=CONTEXT, help="days of history TimesFM sees")
    parser.add_argument("--contexts", default="512,1024,2048",
                        help="context lengths compared for the direct ETo forecast ('' to skip)")
    parser.add_argument("--out", default=str(REPORT_DIR))
    args = parser.parse_args(argv)

    meta, records = weather_data.load(args.data)
    records = add_eto(records)
    split = next(i for i, d in enumerate(records) if d["timestamp"] >= SPLIT)
    clim = Climatology(records[:split])
    origins = list(range(split - 1, len(records) - HORIZON, args.every))
    observed = {k: np.array([[records[i + lead][k] for lead in range(1, HORIZON + 1)] for i in origins])
                for k in WEATHER + ("ETo",)}
    doys = [[(records[i]["timestamp"].timetuple().tm_yday + k - 1) % 365 + 1 for k in range(1, HORIZON + 1)]
            for i in origins]
    climatology = {k: np.array([clim.means[k][row] for row in doys]) for k in observed}
    persistence = {k: np.repeat([[records[i][k]] for i in origins], HORIZON, axis=1) for k in observed}
    baseline = WeatherForecaster(clim, backend="baseline").forecast_at(records, origins)
    tag = f"{Path(args.data).stem}_every{args.every}"

    timesfm, forecaster = None, None
    try:
        forecaster = WeatherForecaster(clim, backend="timesfm", context=args.context, batch_size=26)
    except Exception as exc:
        print(f"TimesFM not available ({type(exc).__name__}); reporting the baselines only.")
    else:
        print("TimesFM forecasts:")
        timesfm = cached(f"timesfm_{tag}_ctx{args.context}",
                         lambda: issue_forecasts(forecaster, records, origins))

    methods = [("Damped anomaly persistence", baseline), ("Persistence", persistence)]
    if timesfm is not None:
        methods.insert(0, ("TimesFM 2.5", timesfm))
    clim_mae = {k: mae_by_lead(climatology[k], observed[k]) for k in observed}

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    first, last = records[origins[0] + 1]["timestamp"].date(), records[origins[-1] + HORIZON]["timestamp"].date()
    lines = [
        "# Forecast skill: weather and reference ET",
        "",
        f"Generated by `python -m fao_model.forecast_skill` on {date.today().isoformat()}.",
        "",
        f"- Data: {meta.get('site')}; forecasts issued every {args.every} days for {first} to {last} "
        f"({len(origins)} forecasts, 1-{HORIZON} days ahead), each using only data up to its issue day.",
        "- **TimesFM 2.5** (200M, zero-shot) forecasts each weather variable from its last "
        f"{args.context} days; reference ET is then computed from the forecast weather with ASCE-PM. "
        "It also forecasts the ETo series directly (\"direct\").",
        "- **Damped anomaly persistence**: climatology plus today's anomaly shrinking by its lag-1 "
        "autocorrelation each day; this is the model's fallback without TimesFM.",
        "- **Persistence**: every coming day equals today. **Climatology**: the smoothed day-of-year "
        "mean of 1995-2014.",
        "- Skill = 1 − MAE / MAE of climatology: above 0 beats climatology, 1 is perfect.",
        "",
        "## 1. Reference ET",
        "",
        "Mean absolute error (mm/day) by days ahead, and of the 7-day total (mm), which matters for "
        "planning irrigation.",
        "",
        "| Method | 1 day | 3 days | 7 days | 14 days | 7-day total |",
        "|---|---|---|---|---|---|",
    ]
    eto_rows = [("Climatology", climatology["ETo"])]
    for name, fc in methods:
        if name == "Persistence":
            eto_rows.append((name, fc["ETo"]))
        else:
            eto_rows.append((f"{name}, from forecast weather", fc["ETo"]))
            eto_rows.append((f"{name}, direct", fc["ETo_direct"]))
    for name, fc in eto_rows:
        mae = mae_by_lead(fc, observed["ETo"])
        week = np.abs(fc[:, :7].sum(axis=1) - observed["ETo"][:, :7].sum(axis=1)).mean()
        lines.append(f"| {name} | " + " | ".join(f"{mae[l - 1]:.2f}" for l in LEADS) + f" | {week:.1f} |")
    shown = [n for n, _ in eto_rows if n.startswith("TimesFM")] + [
        "Damped anomaly persistence, direct", "Persistence"]
    curves = [(n, mae_by_lead(fc, observed["ETo"])) for n, fc in eto_rows if n in shown]
    eto_chart(out / "eto_error.png", curves, mae_by_lead(climatology["ETo"], observed["ETo"]))
    lines += ["", "![Reference ET forecast error](eto_error.png)", ""]

    lines += [
        "## 2. Skill by variable",
        "",
        "Skill against climatology at 1, 3, 7 and 14 days ahead. Rain is left out here: on a dry "
        "climate the absolute error rewards forecasting the median, which is 0 on most days, so it "
        "would make any median forecast look skilful; section 3 scores rain as amounts and events.",
        "",
        "| Variable | Method | 1 day | 3 days | 7 days | 14 days |",
        "|---|---|---|---|---|---|",
    ]
    panels = []
    for key in [k for k in WEATHER if k != "precip"] + ["ETo"]:
        curves = []
        for name, fc in methods[:2]:
            skill = 1.0 - mae_by_lead(fc[key], observed[key]) / clim_mae[key]
            curves.append((name, skill))
            lines.append(f"| {LABELS[key]} | {name} | " + " | ".join(f"{skill[l - 1]:+.2f}" for l in LEADS) + " |")
        panels.append((key, curves))
    skill_chart(out / "skill_by_variable.png", panels)
    lines += ["", "![Skill by variable](skill_by_variable.png)", ""]

    obs_rain = observed["precip"][:, :RAIN_DAYS].sum(axis=1)
    events = obs_rain >= RAIN_EVENT
    lines += [
        "## 3. Rain in the next three days",
        "",
        f"Expected rain over the next {RAIN_DAYS} days against what fell; an event is at least "
        f"{RAIN_EVENT:g} mm (the irrigation rule's threshold). {int(events.sum())} of {len(origins)} "
        "forecast periods had an event. TimesFM and the anomaly baseline use their expected (mean) "
        "rain here, not the median.",
        "",
        "| Method | MAE (mm) | RMSE (mm) | Events caught | False alarms |",
        "|---|---|---|---|---|",
    ]
    rain_methods = [("Climatology", climatology["precip"])] + [
        (name, fc["precip_mean"] if "precip_mean" in fc else fc["precip"]) for name, fc in methods]
    for name, fc in rain_methods:
        total = fc[:, :RAIN_DAYS].sum(axis=1)
        warned = total >= RAIN_EVENT
        hits = int((warned & events).sum())
        lines.append(f"| {name} | {np.abs(total - obs_rain).mean():.1f} | "
                     f"{np.sqrt(((total - obs_rain) ** 2).mean()):.1f} | {hits} of {int(events.sum())} | "
                     f"{int((warned & ~events).sum())} |")

    if timesfm is not None:
        inside = (observed["ETo"] >= timesfm["ETo_q10"]) & (observed["ETo"] <= timesfm["ETo_q90"])
        coverage = inside.mean(axis=0)
        lines += [
            "",
            "## 4. TimesFM uncertainty",
            "",
            "Share of days on which the observed reference ET fell inside TimesFM's 10-90 % interval of "
            "the direct forecast; a calibrated interval holds 80 %.",
            "",
            "| | 1 day | 3 days | 7 days | 14 days |",
            "|---|---|---|---|---|",
            "| Inside the interval | " + " | ".join(f"{coverage[l - 1]:.0%}" for l in LEADS) + " |",
        ]
        picks = [origins.index(min(origins, key=lambda i: abs((records[i]["timestamp"] - target).days)))
                 for target in (datetime(2019, 1, 20), datetime(2019, 5, 10), datetime(2019, 7, 20))]
        examples_chart(out / "examples.png", records,
                       [(origins[j], {k: v[j] for k, v in timesfm.items()}) for j in picks])
        lines += ["", "![Example forecasts](examples.png)", ""]

        contexts = [int(c) for c in args.contexts.split(",") if c.strip()]
        if contexts:
            lines += [
                "## 5. How much history does TimesFM need?",
                "",
                "Direct reference ET forecast with different amounts of history (MAE, mm/day).",
                "",
                "| History (days) | 1 day | 3 days | 7 days | 14 days |",
                "|---|---|---|---|---|",
            ]
            for ctx in contexts:
                print(f"TimesFM direct ETo, context {ctx}:")
                eto_ctx = cached(f"timesfm_eto_{tag}_ctx{ctx}",
                                 lambda: {"ETo": _direct_eto(forecaster.model, ctx, records, origins)})
                mae = mae_by_lead(eto_ctx["ETo"], observed["ETo"])
                lines.append(f"| {ctx} | " + " | ".join(f"{mae[l - 1]:.2f}" for l in LEADS) + " |")
    (out / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {out / 'README.md'} and figures")


def _direct_eto(model, context, records, origins):
    """TimesFM median forecast of the ETo series itself, with `context` days of history."""
    model.compile(context, HORIZON, 32)
    eto = np.array([d["ETo"] for d in records])
    contexts = [eto[max(0, i + 1 - context): i + 1] for i in origins]
    return model.forecast(contexts, HORIZON)[0]


if __name__ == "__main__":
    main()
