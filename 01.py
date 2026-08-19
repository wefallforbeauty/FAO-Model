import math
import random
import tkinter as tk
from tkinter import ttk
from datetime import datetime, timedelta
import numpy as np
try:                                                                                #TIMESFM
    from timesfm import TimesFm, TimesFmHparams
    TIMESFM_AVAILABLE = True
except Exception:
    TIMESFM_AVAILABLE = False
try:                                                                                #MATHPLOTLIB
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    MATPLOTLIB_AVAILABLE = True
except Exception:
    MATPLOTLIB_AVAILABLE = False

SIGMA = 4.903e-9                                                                    #PHYSICALCONSTANTS
GSC = 0.0820
LAMBDA = 2.45
def sat_vp(T):
    return 0.6108 * math.exp(17.27 * T / (T + 237.3))
def slope_vp(T):
    es = sat_vp(T)
    return 4098 * es / (T + 237.3) ** 2
def psychrometric(P):
    return 0.000665 * P
def extraterrestrial_radiation(lat, doy):
    phi = math.radians(lat)
    dr = 1 + 0.033 * math.cos(2 * math.pi * doy / 365)
    decl = 0.409 * math.sin(2 * math.pi * doy / 365 - 1.39)
    cos_ws = -math.tan(phi) * math.tan(decl)
    cos_ws = max(-1.0, min(1.0, cos_ws))
    ws = math.acos(cos_ws)
    Ra = (24 * 60 / math.pi) * GSC * dr * (
        ws * math.sin(phi) * math.sin(decl)
        + math.cos(phi) * math.cos(decl) * math.sin(ws)
    )
    return max(0.0, Ra)
def day_length(lat, doy):
    phi = math.radians(lat)
    decl = 0.409 * math.sin(2 * math.pi * doy / 365 - 1.39)
    cos_ws = -math.tan(phi) * math.tan(decl)
    cos_ws = max(-1.0, min(1.0, cos_ws))
    ws = math.acos(cos_ws)
    return 24 * ws / math.pi
def net_radiation(Rs, T_max, T_min, ea, lat, doy, elev=100.0, albedo=0.23):
    Ra = extraterrestrial_radiation(lat, doy)
    Rso = (0.75 + 2e-5 * elev) * Ra
    if Rso <= 0:
        Rso = Ra
    ratio = Rs / Rso if Rso > 0 else 1.0
    ratio = max(0.3, min(1.0, ratio))

    Rns = (1 - albedo) * Rs
    Rnl = (
        SIGMA
        * ((T_max + 273.16) ** 4 + (T_min + 273.16) ** 4)
        / 2.0
        * (0.34 - 0.14 * math.sqrt(max(0.0, ea)))
        * (1.35 * ratio - 0.35)
    )
    return Rns - Rnl

def asce_pm(d):
    T = d["T_mean"]
    es = sat_vp(T)
    ea = es * d["RH_mean"] / 100.0
    delta = slope_vp(T)
    gamma = psychrometric(d["P"])
    Rn = net_radiation(
        d["Rs"], d["T_max"], d["T_min"], ea,
        d["lat"], d["doy"], d.get("elev", 100.0), albedo=0.23
    )
    eto = (
        0.408 * delta * Rn
        + gamma * (900.0 / (T + 273.0)) * d["u2"] * (es - ea)
    ) / (delta + gamma * (1.0 + 0.34 * d["u2"]))
    return max(0.0, eto)
def thornthwaite(d):                                                                #THORNTHWAITE
    T = d["T_mean"]
    if T <= 0:
        return 0.0
    I = (T / 5.0) ** 1.514
    a = (
        6.75e-7 * I**3
        - 7.71e-5 * I**2
        + 1.792e-2 * I
        + 0.49239
    )
    dl = day_length(d["lat"], d["doy"])
    pet = 16.0 * (10.0 * T / I) ** a * (dl / 12.0)
    return max(0.0, pet)
def hargreaves_base(d):                                                             #HARGREAVES
    Ra = extraterrestrial_radiation(d["lat"], d["doy"])
    return 0.0023 * (d["T_mean"] + 17.8) * math.sqrt(
        max(0.1, d["T_max"] - d["T_min"])
    ) * Ra
def turc_base(d):                                                                   #TURC
    T = d["T_mean"]
    Rs = d["Rs"]
    if T <= 0:
        return 0.0
    if d["RH_mean"] >= 50:
        return 0.013 * (T / (T + 15.0)) * (Rs + 50.0)
    return (
        0.013
        * (T / (T + 15.0))
        * (Rs + 50.0)
        * (1.0 + (50.0 - d["RH_mean"]) / 70.0)
    )
def priestley_taylor_base(d):                                                       #PRIESTLEY-TAYLOR
    T = d["T_mean"]
    es = sat_vp(T)
    ea = es * d["RH_mean"] / 100.0
    delta = slope_vp(T)
    gamma = psychrometric(d["P"])
    Rn = net_radiation(
        d["Rs"], d["T_max"], d["T_min"], ea,
        d["lat"], d["doy"], d.get("elev", 100.0), albedo=0.23
    )
    return (delta / (delta + gamma)) * max(0.0, Rn) / LAMBDA
def jensen_haise_base(d):                                                           #JENSEN-HAISE
    return (0.025 * d["T_mean"] + 0.08) * d["Rs"] / 28.3
class WeatherGenerator:                                                             #WEATHERSIM
    def __init__(self, lat=40.0, elev=100.0, interval_minutes=10):
        self.lat = lat
        self.elev = elev
        self.interval_minutes = interval_minutes
        self.timestamp = datetime.now()
        self.current = None
        self.rng = random.Random()

    def _clip(self, x, lo, hi):
        return max(lo, min(hi, x))
    def _solar_fraction(self, ts):
        hour = ts.hour + ts.minute / 60.0
        season = math.sin(2 * math.pi * (ts.timetuple().tm_yday - 80) / 365.25)
        sunrise = 6.0 - 1.4 * season
        sunset = 18.0 + 1.4 * season
        if hour <= sunrise or hour >= sunset:
            return 0.0
        x = (hour - sunrise) / max(1e-6, sunset - sunrise)
        return math.sin(math.pi * x)
    def _seasonal_temperature(self, ts):
        return 20.0 + 11.0 * math.sin(
            2 * math.pi * (ts.timetuple().tm_yday - 100) / 365.25
        )
    def _new_rain_state(self, rain):
        if rain > 0.05:
            if self.rng.random() < 0.88:
                return self._clip(rain * self.rng.uniform(0.75, 1.15), 0.0, 35.0)
            return 0.0
        if self.rng.random() < 0.06:
            return self.rng.uniform(0.1, 4.0)
        return 0.0
    def generate(self):
        ts = self.timestamp
        doy = ts.timetuple().tm_yday
        seasonal_T = self._seasonal_temperature(ts)
        solar = self._solar_fraction(ts)
        if self.current is None:
            T_mean = seasonal_T + self.rng.uniform(-1.5, 1.5)
            T_range = self.rng.uniform(6.0, 12.0)
            RH_mean = self.rng.uniform(55.0, 80.0)
            u2 = self.rng.uniform(1.0, 5.0)
            rain = self.rng.uniform(0.0, 1.0)
        else:
            prev = self.current
            T_target = seasonal_T + 2.0 * (solar - 0.45)
            T_mean = 0.94 * prev["T_mean"] + 0.06 * T_target + self.rng.gauss(0, 0.35)
            T_mean = self._clip(T_mean, -20.0, 45.0)
            T_range = 0.90 * (prev["T_max"] - prev["T_min"]) + self.rng.gauss(0, 0.25)
            T_range = self._clip(T_range, 5.0, 16.0)
            humidity_target = 78.0 - 0.55 * (T_mean - seasonal_T) - 18.0 * solar
            RH_mean = 0.92 * prev["RH_mean"] + 0.08 * humidity_target + self.rng.gauss(0, 0.9)
            RH_mean = self._clip(RH_mean, 25.0, 98.0)
            wind_target = 1.0 + 5.0 * solar
            u2 = 0.90 * prev["u2"] + 0.10 * wind_target + self.rng.gauss(0, 0.18)
            u2 = self._clip(u2, 0.1, 12.0)
            rain = self._new_rain_state(prev["precip"])
        T_max = T_mean + 0.5 * T_range
        T_min = T_mean - 0.5 * T_range
        RH_max = self._clip(RH_mean + 6.0 + abs(self.rng.gauss(0, 2.0)), RH_mean, 100.0)
        RH_min = self._clip(RH_mean - 6.0 - abs(self.rng.gauss(0, 2.0)), 5.0, RH_mean)
        Ra = extraterrestrial_radiation(self.lat, doy)
        clear_sky = max(0.0, 0.75 * Ra)
        cloud_factor = self._clip(
            0.92 - 0.004 * max(0.0, RH_mean - 55.0) + self.rng.gauss(0, 0.035),
            0.20, 1.0
        )
        Rs = self._clip(clear_sky * solar * cloud_factor, 0.0, max(0.1, clear_sky))
        P = 101.3 * (1.0 - 2.25577e-5 * self.elev) ** 5.25588
        self.current = {
            "timestamp": ts,
            "doy": doy,
            "lat": self.lat,
            "elev": self.elev,
            "T_mean": T_mean,
            "T_max": max(T_mean, T_max),
            "T_min": min(T_mean, T_min),
            "RH_mean": RH_mean,
            "RH_max": RH_max,
            "RH_min": RH_min,
            "u2": u2,
            "Rs": Rs,
            "P": P,
            "precip": rain,
        }

        self.timestamp += timedelta(minutes=self.interval_minutes)
        return self.current.copy()
    def ensure_required(self, data):
        data["T_mean"] = float(data.get("T_mean", 20.0))
        data["T_max"] = max(data.get("T_max", data["T_mean"] + 5.0), data["T_mean"])
        data["T_min"] = min(data.get("T_min", data["T_mean"] - 5.0), data["T_mean"])
        data["RH_mean"] = self._clip(float(data.get("RH_mean", 65.0)), 0.0, 100.0)
        data["RH_max"] = self._clip(float(data.get("RH_max", data["RH_mean"] + 10)), data["RH_mean"], 100.0)
        data["RH_min"] = self._clip(float(data.get("RH_min", data["RH_mean"] - 10)), 0.0, data["RH_mean"])
        data["u2"] = self._clip(float(data.get("u2", 2.0)), 0.0, 20.0)
        data["Rs"] = max(0.0, float(data.get("Rs", 15.0)))
        data["P"] = max(50.0, float(data.get("P", 101.3)))
        data["precip"] = max(0.0, float(data.get("precip", 0.0)))
        data["doy"] = int(data.get("doy", datetime.now().timetuple().tm_yday))
        data["lat"] = float(data.get("lat", self.lat))
        data["elev"] = float(data.get("elev", self.elev))
        return data

class Calibrator:                                                             #CALIBRATIONOFCONSTANTS
    def __init__(self):
        self.coeffs = {
            "hargreaves_a": 1.0,
            "turc_a": 1.0,
            "abtew_a": 0.53,
            "pt_alpha": 1.26,
            "jh_a": 1.0,
            "bc_a": 0.0,
            "bc_b": 1.0,
            "kc_maize": 0.35,
        }
    def calibrate(self, data_batch):
        if len(data_batch) < 3:
            return
        rows = []
        for d in data_batch:
            try:
                eto = asce_pm(d)
                rows.append((
                    eto,
                    hargreaves_base(d),
                    turc_base(d),
                    d["Rs"] / LAMBDA,
                    priestley_taylor_base(d),
                    jensen_haise_base(d),
                    (day_length(d["lat"], d["doy"]) / 24.0)
                    * (0.46 * d["T_mean"] + 8.13),
                ))
            except Exception:
                continue
        if len(rows) < 3:
            return
        a = np.asarray(rows, dtype=float)
        eto = a[:, 0]
        for key, col in (
            ("hargreaves_a", 1),
            ("turc_a", 2),
            ("abtew_a", 3),
            ("pt_alpha", 4),
            ("jh_a", 5),
        ):
            base = a[:, col]
            denom = float(np.dot(base, base))
            if denom > 1e-12:
                # Least-squares coefficient is more stable than sum(y)/sum(x).
                self.coeffs[key] = float(np.dot(base, eto) / denom)
        f = a[:, 6]
        A = np.vstack([f, np.ones(len(f))]).T
        sol, _, _, _ = np.linalg.lstsq(A, eto, rcond=None)
        self.coeffs["bc_b"] = float(sol[0])
        self.coeffs["bc_a"] = float(sol[1])

def vg_theta_arr(h, theta_r, theta_s, alpha, n):                                #VANGENUCHTEN
    h_abs = np.abs(h)
    m = 1.0 - 1.0 / n
    se = (1.0 + (alpha * h_abs) ** n) ** (-m)
    se = np.where(h >= 0, 1.0, se)
    return theta_r + (theta_s - theta_r) * se
def vg_K_arr(h, Ks, theta_r, theta_s, alpha, n):
    theta = vg_theta_arr(h, theta_r, theta_s, alpha, n)
    se = (theta - theta_r) / (theta_s - theta_r)
    se = np.clip(se, 1e-8, 1.0)
    m = 1.0 - 1.0 / n
    K = Ks * np.sqrt(se) * (1.0 - (1.0 - se ** (1.0 / m)) ** m) ** 2
    K = np.where(h >= 0, Ks, K)
    return np.maximum(K, 1e-12)
def vg_capacity_arr(h, theta_r, theta_s, alpha, n):
    C = np.zeros_like(h)
    neg = h < 0
    h_abs = np.abs(h[neg])
    m = 1.0 - 1.0 / n
    C[neg] = (
        (theta_s - theta_r)
        * m * n * alpha
        * (alpha * h_abs) ** (n - 1.0)
        * (1.0 + (alpha * h_abs) ** n) ** (-m - 1.0)
    )
    return np.maximum(C, 1e-6)

class RichardsSolver3D:                                                         #3DFVMRICHARDS
    def __init__(self, nx=5, ny=5, nz=5, dx=10.0, dy=10.0, dz=5.0):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.dx, self.dy, self.dz = dx, dy, dz
        self.h = np.full((nx, ny, nz), -100.0)
        self.params = {
            "theta_r": 0.05,
            "theta_s": 0.43,
            "alpha": 0.036,
            "n": 1.56,
            "Ks": 0.00025,
        }
    def _theta(self, h):
        return vg_theta_arr(h, self.params["theta_r"], self.params["theta_s"],
                            self.params["alpha"], self.params["n"])
    def _K(self, h):
        return vg_K_arr(h, self.params["Ks"], self.params["theta_r"],
                        self.params["theta_s"], self.params["alpha"], self.params["n"])
    def _capacity(self, h):
        return vg_capacity_arr(h, self.params["theta_r"], self.params["theta_s"],
                               self.params["alpha"], self.params["n"])
    def step(self, q_top, dt=0.1):
        h = self.h
        nx, ny, nz = self.nx, self.ny, self.nz
        dx, dy, dz = self.dx, self.dy, self.dz
        C = self._capacity(h)
        K = self._K(h)
        div = np.zeros_like(h)
        for i in range(nx - 1):
            Kf = 2.0 * K[i] * K[i + 1] / (K[i] + K[i + 1] + 1e-12)
            flux = Kf * (h[i] - h[i + 1]) / dx
            div[i] -= flux
            div[i + 1] += flux
        for j in range(ny - 1):
            Kf = 2.0 * K[:, j] * K[:, j + 1] / (K[:, j] + K[:, j + 1] + 1e-12)
            flux = Kf * (h[:, j] - h[:, j + 1]) / dy
            div[:, j] -= flux
            div[:, j + 1] += flux
        for k in range(nz - 1):
            Kf = 2.0 * K[:, :, k] * K[:, :, k + 1] / (K[:, :, k] + K[:, :, k + 1] + 1e-12)
            H_k = h[:, :, k] + k * dz
            H_k1 = h[:, :, k + 1] + (k + 1) * dz
            flux = Kf * (H_k - H_k1) / dz
            div[:, :, k] -= flux
            div[:, :, k + 1] += flux
        div[:, :, 0] += q_top / dz
        div[:, :, -1] -= K[:, :, -1]
        update = dt / C * div
        update = np.clip(update, -2.0, 2.0)
        self.h = np.clip(h + update, -20000.0, 50.0)
    @property
    def theta_mean(self):
        return float(np.mean(self._theta(self.h)))
    @property
    def pressure_head_mean(self):
        return float(np.mean(self.h))

class Forecaster:                                                               #FORECAST
    CONTEXT_LEN = 2048
    HORIZON_LEN = 96
    def __init__(self):
        self.model = None
        self.timesfm = False
        self.status = "TimesFM unavailable; using dynamic fallback."
        if TIMESFM_AVAILABLE:
            try:
                hparams = TimesFmHparams(
                    context_len=self.CONTEXT_LEN,
                    horizon_len=self.HORIZON_LEN,
                    input_patch_len=32,
                    output_patch_len=128,
                    num_layers=20,
                    model_dims=1280,
                )
                self.model = TimesFm(hparams)
                self.model.load_from_checkpoint(repo_id="google/timesfm-1.0-200m")
                self.timesfm = True
                self.status = f"TimesFM active: context {self.CONTEXT_LEN}, horizon {self.HORIZON_LEN}"
            except Exception as exc:
                self.model = None
                self.timesfm = False
                self.status = f"TimesFM init failed ({type(exc).__name__}); using fallback."
    def _fallback(self, series, steps):
        if not series:
            return [0.0] * steps
        if len(series) == 1:
            return [float(series[-1])] * steps
        n = min(len(series), self.CONTEXT_LEN)
        y = np.asarray(series[-n:], dtype=float)
        x = np.arange(n, dtype=float)
        w = max(5, min(96, n // 8))
        recent = float(np.mean(y[-w:]))
        older = float(np.mean(y[max(0, n - 2*w):-w])) if n > w else recent
        slope = (recent - older) / max(1.0, w)
        level = float(np.mean(y))
        out = []
        last = float(y[-1])
        for k in range(1, steps + 1):
            trend = last + slope * k
            decay = math.exp(-k / max(12.0, steps / 2.0))
            value = level + (trend - level) * decay
            out.append(float(value))
        return out
    def forecast(self, series, steps=None):
        steps = int(steps or self.HORIZON_LEN)
        steps = max(1, min(steps, self.HORIZON_LEN))
        if self.timesfm and len(series) >= 32:
            try:
                context = np.asarray(series[-self.CONTEXT_LEN:], dtype=np.float32)
                fc = self.model.forecast(context, horizon_len=steps)
                # Different releases return arrays, tuples, or (forecast, ...).
                if isinstance(fc, tuple):
                    fc = fc[0]
                fc = np.asarray(fc, dtype=float).reshape(-1)
                if len(fc) >= steps and np.all(np.isfinite(fc[:steps])):
                    return fc[:steps].tolist()
            except Exception:
                pass
        return self._fallback(series, steps)
    def forecast_all(self, history, steps=None):
        steps = steps or self.HORIZON_LEN
        return {key: self.forecast(series, steps) for key, series in history.items()}

class App:                                                                      #GUI
    UPDATE_MS = 1000
    MAX_PLOT_POINTS = 1200
    def __init__(self, root):
        self.root = root
        self.root.title("Continuous ET / Richards / TimesFM Dashboard")
        self.root.geometry("1500x900")
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.gen = WeatherGenerator(lat=40.0, elev=100.0, interval_minutes=10)
        self.calib = Calibrator()
        self.richards = RichardsSolver3D(nx=5, ny=5, nz=5)
        self.forecaster = Forecaster()
        self.data_history = []
        self.data_batch = []
        self.history = {}
        self.forecasts = {}
        self.forecast_x = np.arange(1, self.forecaster.HORIZON_LEN + 1)
        self.step_number = 0
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
            self.ax.set_xlabel("Generated step")
            self.ax.set_ylabel("ET / forecast value (mm/day)")
            self.ax.grid(True, alpha=0.25)
            self.canvas = FigureCanvasTkAgg(self.figure, master=right)
            self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            controls = ttk.Frame(right)
            controls.pack(fill=tk.X, pady=4)
            ttk.Label(
                controls,
                text="Plot: last %d historical points + %d-step forecast" %
                     (self.MAX_PLOT_POINTS, self.forecaster.HORIZON_LEN)
            ).pack(side=tk.LEFT)
        else:
            ttk.Label(
                right,
                text="Install matplotlib to enable the live GUI graph."
            ).pack(expand=True)
    def compute_all(self, d):
        eto_asce = asce_pm(d)
        et_pm_maize = self.calib.coeffs["kc_maize"] * eto_asce
        et_th = thornthwaite(d)
        dl = day_length(d["lat"], d["doy"])
        p = dl / 24.0
        f = p * (0.46 * d["T_mean"] + 8.13)
        et_bc = self.calib.coeffs["bc_a"] + self.calib.coeffs["bc_b"] * f
        return {
            "ET_ASCE_PM": eto_asce,
            "ET_PM_maize": et_pm_maize,
            "ET_Thornthwaite": et_th,
            "ET_BlaneyCriddle": max(0.0, et_bc),
            "ET_Turc": max(0.0, self.calib.coeffs["turc_a"] * turc_base(d)),
            "ET_PriestleyTaylor": max(0.0, self.calib.coeffs["pt_alpha"] * priestley_taylor_base(d)),
            "ET_Hargreaves": max(0.0, self.calib.coeffs["hargreaves_a"] * hargreaves_base(d)),
            "ET_JensenHaise": max(0.0, self.calib.coeffs["jh_a"] * jensen_haise_base(d)),
            "ET_Abtew": max(0.0, self.calib.coeffs["abtew_a"] * d["Rs"] / LAMBDA),
        }
    def update(self):
        if not self.running:
            return
        try:
            data = self.gen.ensure_required(self.gen.generate())
            self.step_number += 1
            self.data_history.append(data.copy())
            self.data_batch.append(data.copy())
            if len(self.data_batch) > 256:
                self.data_batch.pop(0)
            self.calib.calibrate(self.data_batch)
            results = self.compute_all(data)
            q_top_cm_s = (
                data["precip"] - results["ET_PM_maize"]
            ) * 0.1 / 86400.0
            for _ in range(10):
                self.richards.step(q_top_cm_s, dt=0.1)
            results["Richards_theta_mean"] = self.richards.theta_mean
            results["Richards_head_mean"] = self.richards.pressure_head_mean
            for key, val in results.items():
                self.history.setdefault(key, []).append(float(val))
            self.forecasts = self.forecaster.forecast_all(
                self.history, steps=self.forecaster.HORIZON_LEN
            )
            self.display(data, results, self.forecasts)
            self.update_plot()
        except Exception as exc:
            self.text.insert(
                tk.END,
                "\n\nERROR IN UPDATE: %s: %s\n" % (type(exc).__name__, exc)
            )
        if self.running:
            self.after_id = self.root.after(self.UPDATE_MS, self.update)
    def update_plot(self):
        if not MATPLOTLIB_AVAILABLE:
            return
        self.ax.clear()
        self.ax.set_title(
            "Live historical equations + long TimesFM forecast "
            f"(step {self.step_number})"
        )
        self.ax.set_xlabel("Generated step")
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
            series = self.history.get(key, [])
            if not series:
                continue
            start = max(0, len(series) - self.MAX_PLOT_POINTS)
            y = np.asarray(series[start:], dtype=float)
            x = np.arange(start + 1, len(series) + 1)
            self.ax.plot(x, y, linewidth=1.2, label=key)
            fc = self.forecasts.get(key, [])
            if fc:
                fx = np.arange(len(series) + 1, len(series) + len(fc) + 1)
                self.ax.plot(fx, fc, linestyle="--", linewidth=1.1,
                             label=f"{key} forecast")
        self.ax.legend(loc="upper left", fontsize=7, ncol=2)
        self.figure.tight_layout()
        self.canvas.draw_idle()
    def display(self, data, results, forecasts):
        ts = data["timestamp"].strftime("%Y-%m-%d %H:%M:%S")
        self.text.insert(
            tk.END,
            "\n" + "=" * 78 + "\n"
            f"STEP {self.step_number} | simulated timestamp: {ts}\n"
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
            f"Precipitation   : {data['precip']:.3f} mm/day\n\n"
        )
        coeffs = self.calib.coeffs
        self.text.insert(
            tk.END,
            "CALIBRATED COEFFICIENTS\n"
            f"Hargreaves a       : {coeffs['hargreaves_a']:.5f}\n"
            f"Turc a             : {coeffs['turc_a']:.5f}\n"
            f"Abtew a            : {coeffs['abtew_a']:.5f}\n"
            f"Priestley-Taylor α : {coeffs['pt_alpha']:.5f}\n"
            f"Jensen-Haise a     : {coeffs['jh_a']:.5f}\n"
            f"Blaney-Criddle a   : {coeffs['bc_a']:.5f}\n"
            f"Blaney-Criddle b   : {coeffs['bc_b']:.5f}\n"
            f"Maize Kc           : {coeffs['kc_maize']:.5f}\n\n"
        )
        self.text.insert(tk.END, "EQUATION RESULTS (mm/day unless noted)\n")
        for key, val in results.items():
            self.text.insert(tk.END, f"{key:25s}: {val:.6f}\n")
        self.text.insert(
            tk.END,
            "\nRICHARDS 3D\n"
            f"Grid               : {self.richards.nx}x{self.richards.ny}x{self.richards.nz}\n"
            f"Mean pressure head : {results['Richards_head_mean']:.4f} cm\n"
            f"Mean water content : {results['Richards_theta_mean']:.6f}\n\n"
        )
        self.text.insert(
            tk.END,
            f"TIMESFM STATUS\n{self.forecaster.status}\n"
            f"Context length    : {self.forecaster.CONTEXT_LEN}\n"
            f"Forecast horizon  : {self.forecaster.HORIZON_LEN} generated steps\n\n"
        )
        self.text.insert(tk.END, "FORECASTS\n")
        for key, fc in forecasts.items():
            if key.startswith("Richards_"):
                continue
            values = ", ".join(f"{v:.3f}" for v in fc)
            self.text.insert(tk.END, f"{key:25s}: {values}\n")
        self.text.insert(
            tk.END,
            f"\nTOTAL RETAINED WEATHER SAMPLES : {len(self.data_history)}\n"
            f"TOTAL RETAINED RESULT POINTS   : {len(next(iter(self.history.values()), []))}\n"
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
if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()
