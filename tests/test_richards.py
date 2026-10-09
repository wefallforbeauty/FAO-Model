import numpy as np
import pytest

from fao_model.richards import RichardsSolver3D, vg_K_arr, vg_head_arr, vg_theta_arr

LOAM = dict(theta_r=0.05, theta_s=0.43, alpha=0.036, n=1.56)
KS = 21.6   # cm/day


def test_head_is_the_inverse_of_water_content():
    h = -np.logspace(-2, 4.3, 50)
    theta = vg_theta_arr(h, **LOAM)
    assert np.allclose(vg_head_arr(theta, **LOAM), h, rtol=1e-6)
    assert vg_head_arr(np.array([LOAM["theta_s"]]), **LOAM)[0] == 0.0


def test_water_balance_closes():
    rng = np.random.default_rng(0)
    s = RichardsSolver3D()
    for _ in range(60):
        s.step(rng.uniform(-0.6, 3.0), dt=1.0)
    t = s.totals
    throughput = t["infiltration"] + t["evaporation"] + t["drainage"]
    assert throughput > 10.0
    assert abs(s.mass_balance_error) < 1e-9 * throughput


def test_gravity_drains_the_profile_downward():
    s = RichardsSolver3D(h_init=-50.0)
    theta_start = s.theta.copy()
    out = s.step(0.0, dt=10.0)
    assert out["drainage"] > 0.0
    assert s.storage == pytest.approx(s.initial_storage - out["drainage"], abs=1e-9)
    assert np.all(s.theta[:, :, 0] < theta_start[:, :, 0])
    # Draining from the bottom with nothing entering at the top: wetter with depth
    assert np.all(np.diff(s.h[2, 2]) > 0.0)


def test_steady_infiltration_gives_unit_gradient():
    # Constant flux q < Ks over free drainage: steady state has K(h) = q in every layer.
    q = 2.0
    s = RichardsSolver3D(nx=1, ny=1)
    for _ in range(60):
        last = s.step(q, dt=1.0)
    assert np.allclose(vg_K_arr(s.h, KS, **LOAM), q, rtol=1e-3)
    assert last["drainage"] == pytest.approx(q, rel=1e-3)
    assert last["runoff"] == 0.0


def test_columns_stay_identical_without_lateral_gradients():
    s = RichardsSolver3D()
    for q in (2.0, -0.4, 0.0, 1.0):
        s.step(q)
    assert np.array_equal(s.h, np.broadcast_to(s.h[:1, :1, :], s.h.shape))


def test_lateral_and_vertical_fluxes_are_divided_by_cell_size():
    s = RichardsSolver3D(nx=2, ny=1, nz=1, dx=10.0, dz=5.0)
    s.set_head(np.array([-50.0, -200.0]).reshape(2, 1, 1))
    K = s._K(s.h)
    rate, q_bottom = s._internal_rates(s.h, K)
    K0, K1 = K[0, 0, 0], K[1, 0, 0]
    qx = 0.5 * (K0 + K1) * (-50.0 + 200.0) / 10.0     # cm/day through the face
    assert rate[0, 0, 0] == pytest.approx(-qx / 10.0 - K0 / 5.0)
    assert rate[1, 0, 0] == pytest.approx(qx / 10.0 - K1 / 5.0)
    assert np.allclose(q_bottom, K[:, :, -1])


def test_dry_surface_limits_evaporation():
    s = RichardsSolver3D(h_init=-15000.0)
    out = s.step(-0.5, dt=5.0)
    assert out["evaporation"] < 0.01
    assert out["evaporation"] + out["evaporation_deficit"] == pytest.approx(2.5)
    assert s.h[:, :, 0].min() >= s.h_min - 1e-6
    assert abs(s.mass_balance_error) < 1e-9


def test_rain_beyond_infiltration_capacity_runs_off():
    s = RichardsSolver3D(nx=1, ny=1, h_init=-5.0)
    out = s.step(100.0, dt=0.2)
    assert out["runoff"] > 10.0
    assert out["infiltration"] + out["runoff"] == pytest.approx(20.0)
    assert np.all(s.theta < s.params["theta_s"])
    assert abs(s.mass_balance_error) < 1e-9


def test_wetting_front_enters_dry_soil():
    # Konya, 2005-11-04: 41.5 mm of rain on a wet top layer over dry soil. With a
    # harmonic mean of K the front stalled and the second layer over-saturated.
    s = RichardsSolver3D(nx=1, ny=1)
    s.set_head(np.array([-38.5, -155.0, -204.1, -202.2, -198.2]).reshape(1, 1, 5))
    theta_start = s.theta.copy()
    s.step(4.12, dt=1.0)
    assert np.all(s.theta < s.params["theta_s"])
    assert s.theta[0, 0, 2] > theta_start[0, 0, 2] + 0.01
    assert abs(s.mass_balance_error) < 1e-9
