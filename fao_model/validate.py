"""Validate the ET equations on real daily weather.

The record is split into a training and a test period. Every equation is
compared with reference ET (this repo's ASCE-PM, and Open-Meteo's
hourly-based ET0) using (1) the published coefficients and (2)
coefficients calibrated on the training period only. Writes a Markdown
report with tables and figures:

    python -m fao_model.validate [--data CSV] [--split 2015-01-01] [--out DIR]
"""
import argparse
import math
from collections import defaultdict, deque
from datetime import date, datetime
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from . import data as weather_data
from .calibration import Calibrator, LITERATURE_COEFFS
from .et import compute_all, thornthwaite_heat_index, thornthwaite_monthly
from .meteo import day_length
from .metrics import compare

METHODS = {
    "ET_Hargreaves": "Hargreaves",
    "ET_Turc": "Turc",
    "ET_PriestleyTaylor": "Priestley-Taylor",
    "ET_JensenHaise": "Jensen-Haise",
    "ET_Abtew": "Abtew",
    "ET_BlaneyCriddle": "Blaney-Criddle",
    "ET_Thornthwaite": "Thornthwaite",
}
UNCALIBRATED = {"ET_Thornthwaite"}
COEFF_LABELS = {
    "hargreaves_a": "Hargreaves a",
    "turc_a": "Turc a",
    "abtew_a": "Abtew k",
    "pt_alpha": "Priestley-Taylor α",
    "jh_a": "Jensen-Haise a",
    "bc_a": "Blaney-Criddle a",
    "bc_b": "Blaney-Criddle b",
}
REPORT_DIR = Path(__file__).resolve().parent.parent / "reports" / "validation"
SPLIT = "2015-01-01"
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# Figure colours: chart surface, ink, hairlines and categorical slots 1-2
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES = ("#2a78d6", "#eb6834")
BLUE_LIGHT, BLUE_DARK = "#86b6ef", "#1c5cab"


def prepare(records, split):
    """Add the Thornthwaite inputs to each record.

    The heat index comes from the training-period monthly climatology,
    T_30d is the trailing 30-day mean temperature (no look-ahead).
    """
    by_month = defaultdict(list)
    for d in records:
        if d["timestamp"] < split:
            by_month[d["timestamp"].month].append(d["T_mean"])
    heat_index = thornthwaite_heat_index(np.mean(t) for t in by_month.values())
    recent = deque(maxlen=30)
    prepared = []
    for d in records:
        recent.append(d["T_mean"])
        prepared.append(dict(d, T_30d=sum(recent) / len(recent), heat_index=heat_index))
    return prepared, heat_index


def run(records, coeffs):
    rows = [compute_all(d, coeffs) for d in records]
    return {key: np.array([r[key] for r in rows]) for key in rows[0]}


def monthly_totals(records, heat_index, asce):
    """Calendar-month Thornthwaite PET and ASCE-PM totals, mm/month."""
    groups = defaultdict(list)
    for i, d in enumerate(records):
        groups[(d["timestamp"].year, d["timestamp"].month)].append(i)
    months, th, ref = [], [], []
    for (year, month), idx in sorted(groups.items()):
        days = [records[i] for i in idx]
        T = np.mean([d["T_mean"] for d in days])
        N = np.mean([day_length(d["lat"], d["doy"]) for d in days])
        months.append(date(year, month, 1))
        th.append(thornthwaite_monthly(T, heat_index, N, len(days)))
        ref.append(asce[idx].sum())
    return months, np.array(th), np.array(ref)


def climatology(records, values):
    """Mean daily value for each calendar month."""
    month = np.array([d["timestamp"].month for d in records])
    return np.array([values[month == m].mean() for m in range(1, 13)])


# ---------------------------------------------------------------- figures

def _style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK_2, labelsize=8)
    ax.xaxis.label.set_color(INK_2)
    ax.yaxis.label.set_color(INK_2)
    ax.title.set_color(INK)


def _hide_unused(axes, n_used):
    """Hide empty panels and show x tick labels on the panels above them."""
    ncols = axes.shape[1]
    for i, ax in enumerate(axes.flat):
        if i >= n_used:
            ax.set_visible(False)
        elif i + ncols >= n_used:
            ax.tick_params(labelbottom=True)


def scatter_grid(path, ref, panels, title, xlabel):
    """One panel per method: estimate against reference, with the 1:1 line."""
    ncols = 4
    nrows = math.ceil(len(panels) / ncols)
    fig = Figure(figsize=(3.1 * ncols, 3.2 * nrows), facecolor=SURFACE)
    axes = fig.subplots(nrows, ncols, sharex=True, sharey=True, squeeze=False)
    lim = 1.05 * max(ref.max(), max(v.max() for _, v, _ in panels))
    for ax, (label, values, m) in zip(axes.flat, panels):
        ax.scatter(ref, values, s=3, color=SERIES[0], alpha=0.25, linewidths=0, rasterized=True)
        ax.plot([0, lim], [0, lim], color=MUTED, linewidth=1)
        ax.set_title(label, fontsize=10)
        ax.text(0.04, 0.96,
                f"bias {m['bias']:+.2f} mm/d ({m['rel_bias']:+.1f}%)\n"
                f"RMSE {m['rmse']:.2f}   r {m['r']:.3f}",
                transform=ax.transAxes, va="top", fontsize=7.5, color=INK_2)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_aspect("equal")
        _style(ax)
    _hide_unused(axes, len(panels))
    fig.supxlabel(xlabel, color=INK_2, fontsize=9)
    fig.supylabel("Estimate (mm/day)", color=INK_2, fontsize=9)
    fig.suptitle(title, color=INK, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def relative_bias_chart(path, rows):
    """Dumbbell: relative bias with literature vs calibrated coefficients."""
    fig = Figure(figsize=(7.5, 0.45 * len(rows) + 1.4), facecolor=SURFACE)
    ax = fig.subplots()
    ys = np.arange(len(rows))[::-1]
    for y, (_, lit, cal) in zip(ys, rows):
        if cal is not None:
            ax.plot([lit, cal], [y, y], color=AXIS, linewidth=2, zorder=1)
            ax.scatter([cal], [y], s=64, color=BLUE_DARK, edgecolors=SURFACE, linewidths=2, zorder=3)
        ax.scatter([lit], [y], s=64, color=BLUE_LIGHT, edgecolors=SURFACE, linewidths=2, zorder=2)
    ax.axvline(0.0, color=INK_2, linewidth=1, zorder=0)
    ax.set_yticks(ys, [label for label, _, _ in rows])
    ax.set_xlabel("Relative bias against ASCE-PM, test period (%)")
    _style(ax)
    ax.grid(False, axis="y")
    handles = [
        Line2D([], [], marker="o", linestyle="", markersize=8, color=BLUE_LIGHT, label="Literature coefficients"),
        Line2D([], [], marker="o", linestyle="", markersize=8, color=BLUE_DARK, label="Calibrated on training period"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=False, labelcolor=INK_2)
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def climatology_chart(path, ref_clim, panels):
    """Small multiples of the mean seasonal cycle, one panel per method."""
    fig = Figure(figsize=(13, 6.2), facecolor=SURFACE)
    axes = fig.subplots(2, 4, sharex=True, sharey=True)
    x = np.arange(1, 13)
    for ax, (label, lit, cal) in zip(axes.flat, panels):
        ax.plot(x, ref_clim, color=INK_2, linewidth=2, solid_capstyle="round")
        ax.plot(x, lit, color=SERIES[0], linewidth=2, solid_capstyle="round")
        if cal is not None:
            ax.plot(x, cal, color=SERIES[1], linewidth=2, solid_capstyle="round")
        ax.set_title(label, fontsize=10)
        ax.set_xticks([1, 4, 7, 10], ["Jan", "Apr", "Jul", "Oct"])
        _style(ax)
    _hide_unused(axes, len(panels))
    fig.supylabel("Mean ET (mm/day)", color=INK_2, fontsize=9)
    fig.suptitle("Mean seasonal cycle, test period", color=INK, fontsize=11)
    handles = [
        Line2D([], [], color=INK_2, linewidth=2, label="ASCE-PM (reference)"),
        Line2D([], [], color=SERIES[0], linewidth=2, label="Literature coefficients"),
        Line2D([], [], color=SERIES[1], linewidth=2, label="Calibrated on training period"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9, frameon=False, labelcolor=INK_2)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(path, dpi=110, facecolor=SURFACE)


def monthly_chart(path, months, ref, th):
    fig = Figure(figsize=(10, 3.8), facecolor=SURFACE)
    ax = fig.subplots()
    ax.plot(months, ref, color=SERIES[0], linewidth=2, label="ASCE-PM, monthly total")
    ax.plot(months, th, color=SERIES[1], linewidth=2, label="Thornthwaite, monthly")
    ax.set_ylabel("ET (mm/month)")
    ax.set_title("Thornthwaite on calendar months, test period", fontsize=11)
    _style(ax)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=8,
              frameon=False, labelcolor=INK_2)
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=SURFACE)


# ----------------------------------------------------------------- report

def _metric_cells(m, unit_digits=2):
    return (f"{m['bias']:+.{unit_digits}f} | {m['rel_bias']:+.1f} | "
            f"{m['rmse']:.{unit_digits}f} | {m['r']:.3f} | {m['nse']:.2f}")


def _coeff(value):
    return "FAO-24 (daily)" if value is None else f"{value:.3f}"


def main(argv=None):
    parser = argparse.ArgumentParser(description="Validate the ET equations on real daily weather.")
    parser.add_argument("--data", default=str(weather_data.DEFAULT_CSV))
    parser.add_argument("--split", default=SPLIT, help="first day of the test period")
    parser.add_argument("--out", default=str(REPORT_DIR))
    args = parser.parse_args(argv)

    meta, records = weather_data.load(args.data)
    split = datetime.fromisoformat(args.split)
    records, heat_index = prepare(records, split)
    train = np.array([d["timestamp"] < split for d in records])
    test = ~train
    train_records = [d for d, t in zip(records, train) if t]

    lit = run(records, LITERATURE_COEFFS)
    asce = lit["ET_ASCE_PM"]
    om = np.array([d["ET0_ref"] for d in records])
    cal_asce = Calibrator()
    cal_asce.calibrate(train_records)
    cal_om = Calibrator()
    cal_om.calibrate(train_records, reference=om[train])
    cal = run(records, cal_asce.coeffs)
    cal_vs_om = run(records, cal_om.coeffs)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    first, last = records[0]["timestamp"].date(), records[-1]["timestamp"].date()
    lines = [
        "# Validation of the ET equations",
        "",
        f"Generated by `python -m fao_model.validate` on {date.today().isoformat()}.",
        "",
        f"- Data: {meta.get('site')} ({meta.get('requested_lat')}° N, {meta.get('requested_lon')}° E, "
        f"{meta.get('elevation_m')} m), {first} to {last}, {len(records)} days; "
        f"{meta.get('source')}, model `{meta.get('model')}`.",
        f"- Training period: {first} to {split.date()} (exclusive), {int(train.sum())} days; "
        f"test period: {split.date()} to {last}, {int(test.sum())} days.",
        "- Reference: daily ASCE-PM grass reference ET (FAO-56) computed by this repo, "
        "unless stated otherwise. Bias = estimate − reference.",
        "- Calibration fits each coefficient by least squares on the training period only; "
        "the test period shows how well that transfers to unseen years.",
        "",
        "## 1. ASCE-PM implementation against Open-Meteo ET0",
        "",
        "Open-Meteo computes FAO-56 ET0 from hourly ERA5 data and sums it per day; "
        "this repo computes it from daily values. Both use the same weather, so this "
        "checks the implementation, not the physics.",
        "",
        "| Period | Bias (mm/day) | Relative bias (%) | RMSE (mm/day) | r | NSE |",
        "|---|---|---|---|---|---|",
    ]
    for name, mask in (("Training", train), ("Test", test), ("All", train | test)):
        lines.append(f"| {name} | {_metric_cells(compare(asce[mask], om[mask]))} |")

    lit_metrics = {k: compare(lit[k][test], asce[test]) for k in METHODS}
    cal_train = {k: compare(cal[k][train], asce[train]) for k in METHODS}
    cal_test = {k: compare(cal[k][test], asce[test]) for k in METHODS}
    lines += [
        "",
        "## 2. Literature coefficients (test period)",
        "",
        "| Method | Bias (mm/day) | Relative bias (%) | RMSE (mm/day) | r | NSE |",
        "|---|---|---|---|---|---|",
    ]
    for k, label in METHODS.items():
        lines.append(f"| {label} | {_metric_cells(lit_metrics[k])} |")
    lines += [
        "",
        "![Literature coefficients against ASCE-PM](scatter_literature.png)",
        "",
        "## 3. Calibrated on the training period",
        "",
        "In-sample RMSE (training period) next to the out-of-sample statistics (test period). "
        "Thornthwaite has no calibration coefficient.",
        "",
        "| Method | Training RMSE | Test bias (mm/day) | Test relative bias (%) | Test RMSE (mm/day) | Test r | Test NSE |",
        "|---|---|---|---|---|---|---|",
    ]
    for k, label in METHODS.items():
        if k in UNCALIBRATED:
            lines.append(f"| {label} | – | – | – | – | – | – |")
        else:
            lines.append(f"| {label} | {cal_train[k]['rmse']:.2f} | {_metric_cells(cal_test[k])} |")
    lines += [
        "",
        "![Calibrated coefficients against ASCE-PM](scatter_calibrated.png)",
        "",
        "![Relative bias before and after calibration](relative_bias.png)",
        "",
        "## 4. Does the calibration reference matter?",
        "",
        "Coefficients fitted on the training period against two references. "
        "Test RMSE is measured against the reference each set was fitted to.",
        "",
        "| Coefficient | Literature | Fitted to ASCE-PM | Fitted to Open-Meteo ET0 |",
        "|---|---|---|---|",
    ]
    for key, label in COEFF_LABELS.items():
        lines.append(f"| {label} | {_coeff(LITERATURE_COEFFS[key])} | "
                     f"{_coeff(cal_asce.coeffs[key])} | {_coeff(cal_om.coeffs[key])} |")
    lines += [
        "",
        "| Method | Test RMSE, fitted to ASCE-PM | Test RMSE, fitted to Open-Meteo ET0 |",
        "|---|---|---|",
    ]
    for k, label in METHODS.items():
        if k in UNCALIBRATED:
            continue
        lines.append(f"| {label} | {cal_test[k]['rmse']:.2f} | "
                     f"{compare(cal_vs_om[k][test], om[test])['rmse']:.2f} |")

    test_records = [d for d, t in zip(records, test) if t]
    months, th_month, ref_month = monthly_totals(test_records, heat_index, asce[test])
    th_m = compare(th_month, ref_month)
    th_daily = lit_metrics["ET_Thornthwaite"]
    lines += [
        "",
        "## 5. Thornthwaite on calendar months",
        "",
        f"Thornthwaite is a monthly method. Heat index I = {heat_index:.1f} from the "
        "training-period monthly climatology. Daily values in the tables above use the "
        "trailing 30-day mean temperature; here the method gets calendar-month means and "
        "is compared with monthly ASCE-PM totals.",
        "",
        "| Time step | Bias | Relative bias (%) | RMSE | r | NSE |",
        "|---|---|---|---|---|---|",
        f"| Calendar month (mm/month) | {_metric_cells(th_m, 1)} |",
        f"| Daily rate (mm/day) | {_metric_cells(th_daily)} |",
        "",
        "![Thornthwaite on calendar months](thornthwaite_monthly.png)",
        "",
        "## 6. Seasonal cycle",
        "",
        "Mean daily ET per calendar month over the test period (mm/day).",
        "",
        "| Series | " + " | ".join(MONTHS) + " |",
        "|---|" + "---|" * 12,
    ]
    test_lit = {k: v[test] for k, v in lit.items()}
    test_cal = {k: v[test] for k, v in cal.items()}
    ref_clim = climatology(test_records, asce[test])
    rows = [("ASCE-PM (reference)", ref_clim),
            ("Open-Meteo ET0", climatology(test_records, om[test]))]
    for k, label in METHODS.items():
        rows.append((f"{label}, literature", climatology(test_records, test_lit[k])))
        if k not in UNCALIBRATED:
            rows.append((f"{label}, calibrated", climatology(test_records, test_cal[k])))
    for label, values in rows:
        lines.append(f"| {label} | " + " | ".join(f"{v:.2f}" for v in values) + " |")
    lines += ["", "![Mean seasonal cycle](monthly_climatology.png)", ""]
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")

    scatter_grid(
        out / "scatter_literature.png", asce[test],
        [(label, lit[k][test], lit_metrics[k]) for k, label in METHODS.items()]
        + [("Open-Meteo ET0", om[test], compare(om[test], asce[test]))],
        "Literature coefficients against ASCE-PM, test period", "ASCE-PM (mm/day)",
    )
    scatter_grid(
        out / "scatter_calibrated.png", asce[test],
        [(label, cal[k][test], cal_test[k]) for k, label in METHODS.items() if k not in UNCALIBRATED],
        "Calibrated on the training period, against ASCE-PM, test period", "ASCE-PM (mm/day)",
    )
    relative_bias_chart(
        out / "relative_bias.png",
        [(label, lit_metrics[k]["rel_bias"], None if k in UNCALIBRATED else cal_test[k]["rel_bias"])
         for k, label in METHODS.items()],
    )
    climatology_chart(
        out / "monthly_climatology.png", ref_clim,
        [(label, climatology(test_records, test_lit[k]),
          None if k in UNCALIBRATED else climatology(test_records, test_cal[k]))
         for k, label in METHODS.items()],
    )
    monthly_chart(out / "thornthwaite_monthly.png", months, ref_month, th_month)
    print(f"Wrote {out / 'README.md'} and figures")


if __name__ == "__main__":
    main()
