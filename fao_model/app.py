import argparse
import tkinter as tk
from tkinter import ttk

import numpy as np

try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    MATPLOTLIB_AVAILABLE = True
    MATPLOTLIB_ERROR = ""
except Exception as exc:     # e.g. matplotlib missing, or Pillow without ImageTk
    MATPLOTLIB_AVAILABLE = False
    MATPLOTLIB_ERROR = f"{type(exc).__name__}: {exc}"

from . import data as weather_data
from .simulation import Simulation


class App:
    UPDATE_MS = 1000
    MAX_PLOT_POINTS = 1200
    MAX_TEXT_LINES = 5000

    def __init__(self, root, records=None, source_name="synthetic weather", update_ms=UPDATE_MS):
        self.root = root
        self.root.title(f"Continuous ET / Richards / TimesFM Dashboard: {source_name}")
        self.root.geometry("1500x900")
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.update_ms = update_ms
        self.sim = Simulation(records)
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
            self.figure = Figure(figsize=(11, 7), dpi=100)
            self.ax = self.figure.add_subplot(111)
            self.ax.set_title("Live ET equations and TimesFM forecasts")
            self.ax.set_xlabel("Day")
            self.ax.set_ylabel("ET / forecast value (mm/day)")
            self.ax.grid(True, alpha=0.25)
            self.canvas = FigureCanvasTkAgg(self.figure, master=right)
            self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            controls = ttk.Frame(right)
            controls.pack(fill=tk.X, pady=4)
            ttk.Label(
                controls,
                text="Plot: last %d historical points + %d-step forecast" %
                     (self.MAX_PLOT_POINTS, self.sim.forecaster.HORIZON_LEN)
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
        self.ax.clear()
        self.ax.set_title(
            "Live historical equations + long TimesFM forecast "
            f"(day {self.sim.step_number})"
        )
        self.ax.set_xlabel("Day")
        self.ax.set_ylabel("Value")
        self.ax.grid(True, alpha=0.25)
        plot_keys = [
            "ET_ASCE_PM",
            "ET_PM_maize",
            "ET_Thornthwaite",
            "ET_Hargreaves",
            "ET_PriestleyTaylor",
        ]
        for key in plot_keys:
            series = self.sim.history.get(key, [])
            if not series:
                continue
            start = max(0, len(series) - self.MAX_PLOT_POINTS)
            y = np.asarray(series[start:], dtype=float)
            x = np.arange(start + 1, len(series) + 1)
            self.ax.plot(x, y, linewidth=1.2, label=key)
            fc = self.sim.forecasts.get(key, [])
            if fc:
                fx = np.arange(len(series) + 1, len(series) + len(fc) + 1)
                self.ax.plot(fx, fc, linestyle="--", linewidth=1.1,
                             label=f"{key} forecast")
        self.ax.legend(loc="upper left", fontsize=7, ncol=2)
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
            f"Blaney-Criddle b   : {coeffs['bc_b']}\n"
            f"Maize Kc           : {coeffs['kc_maize']}\n\n"
        )
        self.text.insert(tk.END, "EQUATION RESULTS (mm/day unless noted)\n")
        for key, val in results.items():
            self.text.insert(tk.END, f"{key:25s}: {val:.6f}\n")
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
            f"runoff {10 * flux['runoff']:.2f}, ET deficit {10 * flux['evaporation_deficit']:.2f}\n"
            f"Water balance error: {10 * sim.richards.mass_balance_error:.2e} mm "
            f"({flux['substeps']} sub-steps)\n\n"
        )
        self.text.insert(
            tk.END,
            f"TIMESFM STATUS\n{sim.forecaster.status}\n"
            f"Context length    : {sim.forecaster.CONTEXT_LEN}\n"
            f"Forecast horizon  : {sim.forecaster.HORIZON_LEN} days\n\n"
        )
        self.text.insert(tk.END, "FORECASTS\n")
        for key, fc in forecasts.items():
            if key.startswith("Richards_"):
                continue
            values = ", ".join(f"{v:.3f}" for v in fc)
            self.text.insert(tk.END, f"{key:25s}: {values}\n")
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
    args = parser.parse_args(argv)
    records, source_name = None, "synthetic weather"
    if not args.synthetic:
        try:
            meta, records = weather_data.load(args.data)
            source_name = meta.get("site", args.data)
        except FileNotFoundError:
            print(f"{args.data} not found; using synthetic weather.")
    root = tk.Tk()
    App(root, records, source_name, args.interval_ms)
    root.mainloop()
