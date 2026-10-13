import argparse
import tkinter as tk
from datetime import timedelta
from tkinter import ttk

try:
    import matplotlib.dates as mdates
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from .plotting import BLUE_LIGHT, INK_2, MUTED, SERIES, SURFACE, style
    MATPLOTLIB_AVAILABLE = True
    MATPLOTLIB_ERROR = ""
except Exception as exc:     # e.g. matplotlib missing, or Pillow without ImageTk
    MATPLOTLIB_AVAILABLE = False
    MATPLOTLIB_ERROR = f"{type(exc).__name__}: {exc}"

from . import data as weather_data
from .simulation import Simulation


class App:
    UPDATE_MS = 1000
    PLOT_DAYS = 120
    MAX_TEXT_LINES = 5000

    def __init__(self, root, records=None, source_name="synthetic weather", update_ms=UPDATE_MS,
                 irrigation=True, forecaster="auto", forecast_every=None):
        self.root = root
        field = "irrigated maize" if irrigation else "rainfed maize"
        self.root.title(f"Continuous ET / Richards / TimesFM Dashboard: {source_name}, {field}")
        self.root.geometry("1500x900")
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.update_ms = update_ms
        self.sim = Simulation(records, irrigation=irrigation, forecaster=forecaster,
                              forecast_every=forecast_every)
        self.running = True
        self.after_id = None
        self.build_gui()
        self.update()

    def build_gui(self):
        main = ttk.Panedwindow(self.root, orient=tk.HORIZONTAL)
        main.pack(fill=tk.BOTH, expand=True)
        left = ttk.Frame(main, padding=5)
        right = ttk.Frame(main, padding=5)
        main.add(left, weight=1)
        main.add(right, weight=3)
        self.text = tk.Text(left, height=45, width=75, font=("Consolas", 9))
        scroll = ttk.Scrollbar(left, orient=tk.VERTICAL, command=self.text.yview)
        self.text.configure(yscrollcommand=scroll.set)
        self.text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        if MATPLOTLIB_AVAILABLE:
            self.figure = Figure(figsize=(11, 7), dpi=100, facecolor=SURFACE)
            self.ax, self.ax_rain = self.figure.subplots(
                2, 1, sharex=True, gridspec_kw={"height_ratios": (3, 1)})
            self.canvas = FigureCanvasTkAgg(self.figure, master=right)
            self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            controls = ttk.Frame(right)
            controls.pack(fill=tk.X, pady=4)
            ttk.Label(
                controls,
                text=f"Plot: last {self.PLOT_DAYS} days and the latest "
                     f"{self.sim.forecaster.horizon}-day forecast. {self.sim.forecaster.status}"
            ).pack(side=tk.LEFT)
        else:
            ttk.Label(
                right,
                text="Live graph unavailable (needs matplotlib with Tk support):\n"
                     + MATPLOTLIB_ERROR
            ).pack(expand=True)

    def update(self):
        if not self.running:
            return
        try:
            data, results, forecasts = self.sim.step()
            self.display(data, results, forecasts)
            self.update_plot()
        except StopIteration:
            self.text.insert(tk.END, "\n\nEND OF WEATHER DATA\n")
            self.text.see(tk.END)
            return
        except Exception as exc:
            self.text.insert(
                tk.END,
                "\n\nERROR IN UPDATE: %s: %s\n" % (type(exc).__name__, exc)
            )
        self._trim_text()
        if self.running:
            self.after_id = self.root.after(self.update_ms, self.update)

    def _trim_text(self):
        # Keep the log bounded so a long run does not grow memory without limit.
        excess = int(self.text.index("end-1c").split(".")[0]) - self.MAX_TEXT_LINES
        if excess > 0:
            self.text.delete("1.0", f"{excess + 1}.0")

    def update_plot(self):
        if not MATPLOTLIB_AVAILABLE:
            return
        sim = self.sim
        days = sim.data_history[-self.PLOT_DAYS:]
        fields = sim.field_history[-self.PLOT_DAYS:]
        n = len(days)
        dates = [d["timestamp"] for d in days]
        ax, axr = self.ax, self.ax_rain
        ax.clear()
        axr.clear()
        for i, key in enumerate(("ET_Hargreaves", "ET_PriestleyTaylor", "ET_Thornthwaite")):
            ax.plot(dates, sim.history[key][-n:], color=MUTED, linewidth=0.8,
                    label="Other ET equations" if i == 0 else None)
        ax.plot(dates, sim.history["ET_ASCE_PM"][-n:], color=SERIES[0], linewidth=1.8,
                label="Reference ET (ASCE-PM)")
        ax.plot(dates, sim.history["ET_crop_actual"][-n:], color=SERIES[1], linewidth=1.8,
                label="Maize ET, actual")
        fc = sim.forecast
        ahead = [fc["issued"] + timedelta(days=k) for k in range(1, len(fc["ETo"]) + 1)]
        if "ETo_q10" in fc:
            ax.fill_between(ahead, fc["ETo_q10"], fc["ETo_q90"], color=BLUE_LIGHT, alpha=0.5,
                            linewidth=0, label="TimesFM 10-90 % range")
        ax.plot(ahead, fc["ETo"], color=SERIES[0], linewidth=1.8, linestyle="--",
                label=f"Forecast issued {fc['issued']:%d %b %Y}")
        for a in (ax, axr):
            a.axvline(fc["issued"], color=INK_2, linewidth=1)
            style(a)
        ax.set_ylabel("ET (mm/day)")
        ax.set_title(f"Day {sim.step_number}: evapotranspiration and the latest forecast", fontsize=10)
        axr.bar(dates, [d["precip"] for d in days], width=0.8, color=SERIES[0], label="Rain")
        axr.bar(ahead, fc["precip_mean"], width=0.8, color=BLUE_LIGHT, label="Expected rain (forecast)")
        top = max([1.0] + [d["precip"] for d in days] + list(fc["precip_mean"])) * 1.15
        irrigated = [d for d, f in zip(dates, fields) if f["irrigation"] > 0]
        axr.scatter(irrigated, [top * 0.9] * len(irrigated), marker="v", s=30, color=SERIES[1],
                    label="Irrigation day")
        axr.set_ylim(0, top)
        axr.set_ylabel("Rain (mm/day)")
        # One legend for both panels, where the ET panel has room (its winter values are low).
        handles = ax.get_legend_handles_labels()[0] + axr.get_legend_handles_labels()[0]
        ax.legend(handles=handles, loc="upper left", fontsize=7, frameon=False, labelcolor=INK_2, ncol=2)
        locator = mdates.AutoDateLocator()
        axr.xaxis.set_major_locator(locator)
        axr.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def display(self, data, results, forecasts):
        sim = self.sim
        ts = data["timestamp"].strftime("%Y-%m-%d")
        self.text.insert(
            tk.END,
            "\n" + "=" * 78 + "\n"
            f"DAY {sim.step_number} | date: {ts}\n"
            + "=" * 78 + "\n"
        )
        self.text.insert(
            tk.END,
            "WEATHER INPUTS\n"
            f"T_mean          : {data['T_mean']:.2f} °C\n"
            f"T_max           : {data['T_max']:.2f} °C\n"
            f"T_min           : {data['T_min']:.2f} °C\n"
            f"RH_mean         : {data['RH_mean']:.2f} %\n"
            f"RH_max          : {data['RH_max']:.2f} %\n"
            f"RH_min          : {data['RH_min']:.2f} %\n"
            f"Wind u2         : {data['u2']:.2f} m/s\n"
            f"Solar radiation : {data['Rs']:.2f} MJ/m²/day\n"
            f"Pressure        : {data['P']:.2f} kPa\n"
            f"Precipitation   : {data['precip']:.3f} mm/day\n"
        )
        if "ET0_ref" in data:
            self.text.insert(tk.END, f"Open-Meteo ET0  : {data['ET0_ref']:.2f} mm/day\n")
        self.text.insert(tk.END, "\n")
        coeffs = {k: ("FAO-24" if v is None else f"{v:.5f}") for k, v in sim.calib.coeffs.items()}
        self.text.insert(
            tk.END,
            "CALIBRATED COEFFICIENTS (fitted on earlier days)\n"
            f"Hargreaves a       : {coeffs['hargreaves_a']}\n"
            f"Turc a             : {coeffs['turc_a']}\n"
            f"Abtew a            : {coeffs['abtew_a']}\n"
            f"Priestley-Taylor α : {coeffs['pt_alpha']}\n"
            f"Jensen-Haise a     : {coeffs['jh_a']}\n"
            f"Blaney-Criddle a   : {coeffs['bc_a']}\n"
            f"Blaney-Criddle b   : {coeffs['bc_b']}\n\n"
        )
        self.text.insert(tk.END, "EQUATION RESULTS (mm/day unless noted)\n")
        for key, val in results.items():
            self.text.insert(tk.END, f"{key:25s}: {val:.6f}\n")
        f = sim.field.last
        season = "" if f["t"] is None else f", day {f['t']} of {sim.field.crop.season_length}"
        self.text.insert(
            tk.END,
            f"\nMAIZE FIELD ({'irrigated' if sim.field.irrigation else 'rainfed'})\n"
            f"Stage              : {f['stage']}{season}\n"
            f"Kcb / Kc_max / fc  : {f['kcb']:.3f} / {f['kc_max']:.3f} / {f['fc']:.2f}\n"
            f"Root depth         : {f['zr']:.2f} m\n"
            f"Transpiration      : potential {f['tp']:.2f}, actual {f['ta']:.2f} mm/day\n"
            f"Soil evaporation   : potential {f['ep']:.2f}, actual {f['ea']:.2f} mm/day\n"
            f"Root-zone depletion: {f['dr']:.1f} of {f['taw']:.1f} mm available "
            f"(irrigate at {sim.field.mad:.0%})\n"
            f"Irrigation today   : {f['irrigation']:.1f} mm\n"
        )
        flux = sim.richards.last_step
        self.text.insert(
            tk.END,
            "\nRICHARDS 3D\n"
            f"Grid               : {sim.richards.nx}x{sim.richards.ny}x{sim.richards.nz}\n"
            f"Mean pressure head : {results['Richards_head_mean']:.4f} cm\n"
            f"Mean water content : {results['Richards_theta_mean']:.6f}\n"
            f"Storage            : {10 * sim.richards.storage:.2f} mm\n"
            f"Today (mm)         : infiltration {10 * flux['infiltration']:.2f}, "
            f"evaporation {10 * flux['evaporation']:.2f}, drainage {10 * flux['drainage']:.3f}, "
            f"runoff {10 * flux['runoff']:.2f}, transpiration {10 * flux['transpiration']:.2f}\n"
            f"Water balance error: {10 * sim.richards.mass_balance_error:.2e} mm "
            f"({flux['substeps']} sub-steps)\n\n"
        )
        issued = forecasts["issued"]
        band = "ETo_q10" in forecasts
        self.text.insert(
            tk.END,
            f"FORECAST issued {issued:%Y-%m-%d} (a new one every {sim.forecast_every} days)\n"
            f"{sim.forecaster.status}\n"
            "date         T_max  T_min   rain    ETo" + ("   10-90 %" if band else "") + "\n"
        )
        for k in range(min(7, len(forecasts["ETo"]))):
            day = issued + timedelta(days=k + 1)
            row = (f"{day:%Y-%m-%d}  {forecasts['T_max'][k]:5.1f}  {forecasts['T_min'][k]:5.1f}  "
                   f"{forecasts['precip_mean'][k]:5.1f}  {forecasts['ETo'][k]:5.2f}")
            if band:
                row += f"  {forecasts['ETo_q10'][k]:4.1f}-{forecasts['ETo_q90'][k]:.1f}"
            self.text.insert(tk.END, row + "\n")
        self.text.insert(
            tk.END,
            f"\nTOTAL RETAINED WEATHER SAMPLES : {len(sim.data_history)}\n"
            f"TOTAL RETAINED RESULT POINTS   : {len(next(iter(sim.history.values()), []))}\n"
        )
        self.text.see(tk.END)

    def close(self):
        self.running = False
        if self.after_id is not None:
            try:
                self.root.after_cancel(self.after_id)
            except Exception:
                pass
        self.root.destroy()


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run the ET / soil water dashboard.")
    parser.add_argument("--data", default=str(weather_data.DEFAULT_CSV),
                        help="daily weather CSV written by fao_model.data")
    parser.add_argument("--synthetic", action="store_true",
                        help="use the synthetic weather generator instead of real data")
    parser.add_argument("--interval-ms", type=int, default=App.UPDATE_MS,
                        help="milliseconds between simulated days")
    parser.add_argument("--irrigation", choices=("auto", "none"), default="auto",
                        help="automatic irrigation of the maize field, or a rainfed field")
    parser.add_argument("--forecaster", choices=("auto", "timesfm", "baseline"), default="auto",
                        help="TimesFM 2.5 if installed (auto), or the statistical baseline")
    parser.add_argument("--forecast-every", type=int, default=None,
                        help="days between forecasts (default 7 with TimesFM, else 1)")
    args = parser.parse_args(argv)
    records, source_name = None, "synthetic weather"
    if not args.synthetic:
        try:
            meta, records = weather_data.load(args.data)
            source_name = meta.get("site", args.data)
        except FileNotFoundError:
            print(f"{args.data} not found; using synthetic weather.")
    root = tk.Tk()
    App(root, records, source_name, args.interval_ms, irrigation=args.irrigation == "auto",
        forecaster=args.forecaster, forecast_every=args.forecast_every)
    root.mainloop()
