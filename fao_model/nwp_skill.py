"""Weather-model forecasts against TimesFM and the statistical baseline, 1-7 days ahead.

Every day of the evaluation period gets forecasts issued the evening before
for the next 7 days, compared with the ERA5 observations:

- NWP: ECMWF IFS 0.25° and NCEP GFS as archived by Open-Meteo, raw and
  bias-corrected (a decaying-average bias per variable and lead, from the
  forecast-observation pairs already known on the issue day);
- TimesFM 2.5 and the damped anomaly persistence baseline;
- climatology (1995-2023) and persistence.

    python -m fao_model.nwp_skill [--no-timesfm] [--first 2024-05-01] [--out DIR]
"""
import argparse
from datetime import date, timedelta
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure

from . import data as weather_data
from . import nwp
from .forecast import (
    CONTEXT, WEATHER, Climatology, WeatherForecaster, _eto_from_weather, _make_consistent, add_eto,
)
from .forecast_skill import LABELS, RAIN_DAYS, RAIN_EVENT, cached, issue_forecasts
from .plotting import BLUE_LIGHT, CATEGORICAL, INK, INK_2, MUTED, SERIES, SURFACE, style

REPORT_DIR = Path(__file__).resolve().parent.parent / "reports" / "nwp"
FIRST = date(2024, 5, 1)        # six weeks after the archive starts: warm-up of the bias correction
BIAS_WEIGHT = 0.05              # memory of about 20 days
CORRECTED = ("T_max", "T_min", "RH_max", "RH_min", "u2", "Rs")
LEADS = nwp.LEADS
H = len(LEADS)
SHOWN_LEADS = (1, 3, 5, 7)
WET_DAY = 1.0                   # mm


def known_bias(table, lead, key, observed, weight=BIAS_WEIGHT):
    """{day: bias of the forecasts of `key` at `lead` known at the end of that day}.

    Decaying average b = (1 - w) b + w (forecast - observed) over the valid
    days observed so far, starting from 0.
    """
    errors = {d: f[key] - observed[d][key] for d, f in table[lead].items() if d in observed}
    day, last = min(errors), max(errors)
    bias, known = 0.0, {}
    while day <= last:
        if day in errors:
            bias = (1.0 - weight) * bias + weight * errors[day]
        known[day] = bias
        day += timedelta(days=1)
    return known


def nwp_forecasts(table, records, origins, observed, correct):
    """Forecast arrays (n, 7) from an archived model, with ETo computed from the forecast weather."""
    issues = [records[i]["timestamp"].date() for i in origins]
    out = {}
    for key in WEATHER:
        values = np.array([[table[lead][o + timedelta(days=lead)][key] for lead in LEADS] for o in issues])
        if correct and key in CORRECTED:
            for j, lead in enumerate(LEADS):
                known = known_bias(table, lead, key, observed)
                values[:, j] -= [known.get(o, 0.0) for o in issues]
        out[key] = values
    _make_consistent(out, [records[i] for i in origins])
    out["precip_mean"] = out["precip"]
    out["ETo"] = np.array([_eto_from_weather(records[i], {k: out[k][j] for k in WEATHER})
                           for j, i in enumerate(origins)])
    return out


def mae(forecast, observed):
    return np.abs(forecast - observed).mean(axis=0)


def detection(forecast_total, observed_total, threshold):
    """(hits, misses, false alarms) of forecasting at least `threshold`."""
    warned, happened = forecast_total >= threshold, observed_total >= threshold
    return int((warned & happened).sum()), int((~warned & happened).sum()), int((warned & ~happened).sum())


# ---------------------------------------------------------------- figures

def eto_chart(path, curves, climatology):
    fig = Figure(figsize=(8, 4.4), facecolor=SURFACE)
    ax = fig.subplots()
    ax.plot(LEADS, climatology, color=MUTED, linewidth=2, label="Climatology")
    for (label, values), color in zip(curves, CATEGORICAL[:len(curves)], strict=True):
        ax.plot(LEADS, values, color=color, linewidth=2, marker="o", markersize=4, label=label)
    ax.set_xlabel("Days ahead")
    ax.set_ylabel("Mean absolute error (mm/day)")
    ax.set_title("Reference ET forecast error", fontsize=10)
    ax.set_xticks(LEADS)
    ax.set_ylim(bottom=0)
    style(ax)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_2, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def skill_chart(path, panels):
    fig = Figure(figsize=(13, 6.4), facecolor=SURFACE)
    axes = fig.subplots(2, 4, sharex=True, sharey=True)
    for ax, (key, curves) in zip(axes.flat, panels):
        ax.axhline(0.0, color=INK_2, linewidth=1)
        for (label, values), color in zip(curves, CATEGORICAL[:len(curves)], strict=True):
            ax.plot(LEADS, values, color=color, linewidth=2, label=label)
        ax.set_title(LABELS[key], fontsize=9)
        ax.set_xticks(LEADS)
        style(ax)
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    axes[0, 0].set_ylim(-0.3, 1.0)
    fig.supxlabel("Days ahead", color=INK_2, fontsize=9)
    fig.supylabel("Skill against climatology (1 − MAE / MAE of climatology)", color=INK_2, fontsize=9)
    fig.suptitle("Where forecasts beat climatology (above 0)", color=INK, fontsize=11)
    fig.tight_layout()
    handles, labels = axes.flat[0].get_legend_handles_labels()
    box = axes.flat[-1].get_position()
    fig.legend(handles, labels, loc="center", bbox_to_anchor=(box.x0 + box.width / 2, box.y0 + box.height / 2),
               frameon=False, fontsize=9, labelcolor=INK_2)
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def rain_chart(path, days, observed, series):
    """Daily rain one day ahead over a wet half-year: observed bars against forecasts."""
    fig = Figure(figsize=(12, 4.2), facecolor=SURFACE)
    ax = fig.subplots()
    ax.bar(days, observed, width=0.9, color=BLUE_LIGHT, label="Observed (ERA5)")
    for (label, values), color in zip(series, CATEGORICAL[:len(series)], strict=True):
        ax.plot(days, values, color=color, linewidth=1.5, label=label)
    ax.set_ylabel("Rain (mm/day)")
    ax.set_title("Rain forecast one day ahead", fontsize=10)
    style(ax)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_2, loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


# ----------------------------------------------------------------- report

def main(argv=None):
    parser = argparse.ArgumentParser(description="Weather-model forecasts against TimesFM and baselines.")
    parser.add_argument("--first", type=date.fromisoformat, default=FIRST, help="first forecast day")
    parser.add_argument("--last", type=date.fromisoformat, default=None, help="last forecast day")
    parser.add_argument("--no-timesfm", action="store_true", help="leave TimesFM out")
    parser.add_argument("--out", default=str(REPORT_DIR))
    args = parser.parse_args(argv)

    meta, records = weather_data.load_all()
    records = add_eto(records)
    index = {d["timestamp"].date(): i for i, d in enumerate(records)}
    observed_by_day = {d["timestamp"].date(): d for d in records}
    last = min(args.last or records[-1]["timestamp"].date(), records[-1]["timestamp"].date())
    models = {name: nwp.load(model) for model, name in nwp.MODELS.items()}
    issues, day = [], args.first - timedelta(days=1)
    while day + timedelta(days=H) <= last:
        if all(t.get(lead, {}).get(day + timedelta(days=lead)) for _, t in models.values() for lead in LEADS):
            issues.append(day)
        day += timedelta(days=1)
    origins = [index[o] for o in issues]
    observed = {k: np.array([[records[i + lead][k] for lead in LEADS] for i in origins])
                for k in WEATHER + ("ETo",)}

    clim = Climatology(records[:index[date(2024, 1, 1)]])
    doys = [[(records[i]["timestamp"].timetuple().tm_yday + lead - 1) % 365 + 1 for lead in LEADS]
            for i in origins]
    methods = {
        "Climatology": {k: np.array([clim.means[k][row] for row in doys]) for k in observed},
        "Persistence": {k: np.repeat([[records[i][k]] for i in origins], H, axis=1) for k in observed},
        "Damped anomaly persistence": WeatherForecaster(clim, backend="baseline", horizon=H)
        .forecast_at(records, origins),
    }
    for name, (_, table) in models.items():
        methods[f"{name}, raw"] = nwp_forecasts(table, records, origins, observed_by_day, correct=False)
        methods[f"{name}, bias-corrected"] = nwp_forecasts(table, records, origins, observed_by_day, correct=True)
    if not args.no_timesfm:
        try:
            forecaster = WeatherForecaster(clim, backend="timesfm", batch_size=26)
        except Exception as exc:
            print(f"TimesFM not available ({type(exc).__name__}); left out.")
        else:
            print("TimesFM forecasts:")
            key = f"timesfm_nwp_{issues[0]:%Y%m%d}_{issues[-1]:%Y%m%d}_{len(issues)}_ctx{CONTEXT}"
            fc = cached(key, lambda: issue_forecasts(forecaster, records, origins))
            methods["TimesFM 2.5"] = {k: v[:, :H] for k, v in fc.items()}
            del forecaster

    ecmwf, gfs = (f"{name}, bias-corrected" for name in nwp.MODELS.values())
    shown = [m for m in (ecmwf, gfs, "TimesFM 2.5", "Damped anomaly persistence") if m in methods]
    clim_mae = {k: mae(methods["Climatology"][k], observed[k]) for k in observed}
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    first_day, last_day = issues[0] + timedelta(days=1), issues[-1] + timedelta(days=H)
    lines = [
        "# Weather-model forecasts against TimesFM",
        "",
        f"Generated by `python -m fao_model.nwp_skill` on {date.today().isoformat()}.",
        "",
        f"- Site: {meta.get('site')}. Forecasts issued every evening from {issues[0]} to {issues[-1]} "
        f"({len(issues)} forecasts) for the next 1-{H} days, compared with ERA5 observations up to {last_day}.",
        "- **Weather models (NWP)**: ECMWF IFS 0.25° and NCEP GFS from Open-Meteo's Previous Runs API: "
        "for each hour, the value forecast 24 × lead hours earlier, aggregated to daily weather like the "
        "observations. Reference ET is computed from the forecast weather with ASCE-PM.",
        f"- **Bias-corrected**: minus a running bias (decaying average with weight {BIAS_WEIGHT}, about "
        f"{1 / BIAS_WEIGHT:.0f} days of memory) of each variable and lead, using only the "
        "forecast-observation pairs known on the issue day; rain is not corrected. The archive starts in "
        f"mid-March 2024, so scoring starts on {args.first} after a warm-up.",
        "- **TimesFM 2.5**, **damped anomaly persistence**, **climatology** (1995-2023) and "
        "**persistence** as in the forecast skill report.",
        "- Caveat: ERA5, the observations here, is made by ECMWF with a model close to IFS, which "
        "may favour ECMWF over GFS.",
        "",
        "## 1. Reference ET",
        "",
        "Mean absolute error (mm/day) by days ahead, and of the 7-day total (mm).",
        "",
        "| Method | " + " | ".join(f"{lead} day{'s' if lead > 1 else ''}" for lead in SHOWN_LEADS) + " | 7-day total |",
        "|---|" + "---|" * (len(SHOWN_LEADS) + 1),
    ]
    order = ["Climatology", "Persistence", "Damped anomaly persistence", "TimesFM 2.5"] + [
        f"{name}, {kind}" for name in nwp.MODELS.values() for kind in ("raw", "bias-corrected")]
    for name in [m for m in order if m in methods]:
        fc = methods[name]
        err = mae(fc["ETo"], observed["ETo"])
        week = np.abs(fc["ETo"].sum(axis=1) - observed["ETo"].sum(axis=1)).mean()
        lines.append(f"| {name} | " + " | ".join(f"{err[lead - 1]:.2f}" for lead in SHOWN_LEADS) + f" | {week:.1f} |")
    eto_chart(out / "eto_error.png", [(m, mae(methods[m]["ETo"], observed["ETo"])) for m in shown],
              clim_mae["ETo"])
    lines += ["", "![Reference ET forecast error](eto_error.png)", ""]

    lines += [
        "## 2. Systematic errors of the raw weather models",
        "",
        "Mean of forecast − observation, which the bias correction removes as it learns it.",
        "",
        "| Variable | " + " | ".join(f"{name}, {lead} day{'s' if lead > 1 else ''}"
                                   for name in nwp.MODELS.values() for lead in (1, 7)) + " |",
        "|---|" + "---|" * (2 * len(nwp.MODELS)),
    ]
    for key in WEATHER:
        cells = []
        for name in nwp.MODELS.values():
            bias = (methods[f"{name}, raw"][key] - observed[key]).mean(axis=0)
            cells += [f"{bias[0]:+.2f}", f"{bias[-1]:+.2f}"]
        lines.append(f"| {LABELS[key]} | " + " | ".join(cells) + " |")

    lines += [
        "",
        "## 3. Skill by variable",
        "",
        "Skill against climatology at 1, 3 and 7 days ahead (rain is scored in section 4).",
        "",
        "| Variable | Method | 1 day | 3 days | 7 days |",
        "|---|---|---|---|---|",
    ]
    panels = []
    for key in [k for k in WEATHER if k != "precip"] + ["ETo"]:
        curves = []
        for name in shown:
            skill = 1.0 - mae(methods[name][key], observed[key]) / clim_mae[key]
            curves.append((name, skill))
            lines.append(f"| {LABELS[key]} | {name} | {skill[0]:+.2f} | {skill[2]:+.2f} | {skill[6]:+.2f} |")
        panels.append((key, curves))
    skill_chart(out / "skill_by_variable.png", panels)
    lines += ["", "![Skill by variable](skill_by_variable.png)", ""]

    obs_wet = observed["precip"][:, 0]
    obs_3day = observed["precip"][:, :RAIN_DAYS].sum(axis=1)
    lines += [
        "## 4. Rain",
        "",
        f"Wet day: at least {WET_DAY:g} mm tomorrow ({int((obs_wet >= WET_DAY).sum())} of {len(issues)} days). "
        f"Rain event: at least {RAIN_EVENT:g} mm within the next {RAIN_DAYS} days "
        f"({int((obs_3day >= RAIN_EVENT).sum())} of {len(issues)} forecasts). TimesFM and the baseline use "
        "their expected (mean) rain.",
        "",
        "| Method | Wet days caught | Wet-day false alarms | Events caught | Event false alarms | 3-day MAE (mm) |",
        "|---|---|---|---|---|---|",
    ]
    for name in [m for m in order if m in methods]:
        rain = methods[name]["precip_mean"] if "precip_mean" in methods[name] else methods[name]["precip"]
        wet = detection(rain[:, 0], obs_wet, WET_DAY)
        total = rain[:, :RAIN_DAYS].sum(axis=1)
        event = detection(total, obs_3day, RAIN_EVENT)
        lines.append(f"| {name} | {wet[0]} of {wet[0] + wet[1]} | {wet[2]} | {event[0]} of {event[0] + event[1]} | "
                     f"{event[2]} | {np.abs(total - obs_3day).mean():.1f} |")
    winter = [j for j, o in enumerate(issues) if date(2024, 10, 1) <= o < date(2025, 5, 1)]
    if winter:
        series = [(n, methods[m]["precip_mean"][winter, 0]) for n, m in
                  ((nwp.MODELS["ecmwf_ifs025"], f"{nwp.MODELS['ecmwf_ifs025']}, raw"), ("TimesFM 2.5", "TimesFM 2.5"),
                   ("Damped anomaly persistence", "Damped anomaly persistence")) if m in methods]
        rain_chart(out / "rain_winter.png", [issues[j] + timedelta(days=1) for j in winter],
                   observed["precip"][winter, 0], series)
        lines += ["", "![Rain forecast one day ahead, winter 2024/25](rain_winter.png)", ""]
    lines += [
        "## 5. What this means for irrigation",
        "",
        "The water balance report (section 6) found that even a perfect 3-day rain forecast saves no "
        "irrigation water at this site in the maize season, because rain is rare then. Better forecasts "
        "cannot beat that bound, so the irrigation rule is not rerun with these models.",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {out / 'README.md'} and figures")


if __name__ == "__main__":
    main()
