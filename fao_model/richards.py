import numpy as np


def vg_theta_arr(h, theta_r, theta_s, alpha, n):
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


class RichardsSolver3D:
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
