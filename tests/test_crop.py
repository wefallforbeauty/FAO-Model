import contextlib
import io
from datetime import date, datetime, timedelta

import numpy as np
import pytest

from fao_model.crop import Crop, adjust_for_climate, root_weights

MAIZE = Crop()


def test_season_calendar():
    assert MAIZE.season_length == 150
    assert MAIZE.day_of_season(date(2020, 5, 1)) == 1
    assert MAIZE.day_of_season(date(2020, 4, 30)) is None
    assert MAIZE.day_of_season(date(2020, 9, 27)) == 150
    assert MAIZE.day_of_season(date(2020, 9, 28)) is None
    assert [MAIZE.stage(t) for t in (1, 30, 31, 70, 71, 120, 121, 150)] == [
        "initial", "initial", "development", "development",
        "mid-season", "mid-season", "late season", "late season"]


def test_basal_kc_trapezoid():
    # FAO-56 Table 17, maize dried in the field: 0.15 / 1.15 / 0.15
    assert MAIZE.basal_kc(1) == MAIZE.basal_kc(30) == 0.15
    assert MAIZE.basal_kc(50) == pytest.approx(0.15 + 1.0 * 20 / 40)
    assert MAIZE.basal_kc(70) == pytest.approx(1.15)
    assert MAIZE.basal_kc(120) == 1.15
    assert MAIZE.basal_kc(135) == pytest.approx(1.15 - 1.0 * 15 / 30)
    assert MAIZE.basal_kc(150) == pytest.approx(0.15)


def test_height_and_roots_grow_during_development_only():
    assert MAIZE.root_depth(1) == MAIZE.root_depth(30) == 0.15
    assert MAIZE.root_depth(50) == pytest.approx(0.15 + 1.05 * 0.5)
    assert MAIZE.root_depth(70) == MAIZE.root_depth(150) == 1.2
    assert MAIZE.height(1) == 0.05 and MAIZE.height(150) == 2.0


def test_max_kc_and_canopy_cover():
    # FAO-56 eq. 72 at the standard climate (u2 = 2 m/s, RHmin = 45 %) is 1.2
    assert MAIZE.max_kc(0.5, 1.0, 2.0, 45.0) == pytest.approx(1.2)
    # windy, dry: 1.2 + (0.04 * 2 - 0.004 * (-25)) * (2/3)^0.3
    assert MAIZE.max_kc(1.15, 2.0, 4.0, 20.0) == pytest.approx(1.2 + 0.18 * (2 / 3) ** 0.3)
    assert MAIZE.max_kc(1.25, 0.5, 2.0, 45.0) == pytest.approx(1.30)    # at least Kcb + 0.05
    # FAO-56 eq. 76
    assert MAIZE.canopy_cover(0.15, 1.2, 0.05) == 0.0
    assert MAIZE.canopy_cover(0.675, 1.2, 1.0) == pytest.approx(0.5 ** 1.5)
    assert MAIZE.canopy_cover(1.2, 1.2, 2.0) == 0.99


def test_potential_rates_split_reference_et():
    c, tp, ep = MAIZE.potential_rates(date(2020, 7, 20), 6.0, 2.0, 45.0)   # mid-season
    assert c["stage"] == "mid-season"
    assert tp == pytest.approx(1.15 * 6.0)
    assert ep == pytest.approx(min(c["kc_max"] - 1.15, c["few"] * c["kc_max"]) * 6.0)
    c, tp, ep = MAIZE.potential_rates(date(2020, 1, 15), 1.0, 2.0, 45.0)    # bare soil
    assert c["stage"] == "bare soil" and tp == 0.0 and ep == pytest.approx(1.2)


def test_feddes_stress_function():
    h = np.array([0.0, -15.0, -22.5, -30.0, -300.0, -325.0, -4162.5, -8000.0, -9000.0])
    alpha = MAIZE.feddes_alpha(h, tp=0.5)   # high demand: h3 = -325 cm
    assert np.allclose(alpha, [0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 0.5, 0.0, 0.0])
    low = MAIZE.feddes_alpha(np.array([-500.0, -600.0]), tp=0.1)   # low demand: h3 = -600 cm
    assert np.allclose(low, 1.0)


def test_root_weights_fall_linearly_with_depth():
    w = root_weights(100.0, 5.0, 30)
    assert w.sum() == pytest.approx(1.0)
    assert np.all(np.diff(w[:20]) < 0) and np.all(w[20:] == 0.0)
    assert w[0] == pytest.approx((5.0 - 25.0 / 200.0) / 50.0)
    with pytest.raises(ValueError):
        root_weights(200.0, 5.0, 30)


def test_climate_adjustment_of_mid_season_kcb():
    days = []
    for year in (2001, 2002):
        day = date(year, 1, 1)
        while day.year == year:
            days.append({"timestamp": datetime(day.year, day.month, day.day), "u2": 3.0, "RH_min": 30.0})
            day += timedelta(days=1)
    adjusted = adjust_for_climate(MAIZE, days)
    # FAO-56 eq. 70: + (0.04 * 1 - 0.004 * (-15)) * (2 / 3)^0.3
    assert adjusted.kcb_mid == pytest.approx(1.15 + 0.10 * (2.0 / 3.0) ** 0.3)
    assert adjusted.kcb_end == 0.15     # below 0.45: not adjusted


def test_crop_curves_match_pyfao56():
    pd = pytest.importorskip("pandas")
    fao = pytest.importorskip("pyfao56")
    start = date(2020, 5, 1)
    n = MAIZE.season_length
    index = [(start + timedelta(days=i)).strftime("%Y-%j") for i in range(n)]
    wth = fao.Weather()
    wth.wdata = pd.DataFrame({"Srad": 25.0, "Tmax": 30.0, "Tmin": 15.0, "Vapr": np.nan, "Tdew": np.nan,
                              "RHmax": 70.0, "RHmin": 25.0, "Wndsp": 3.0, "Rain": 0.0, "ETref": 6.0,
                              "MorP": "M"}, index=index)
    wth.z, wth.lat, wth.wndht = 1000.0, 37.5, 2.0
    par = fao.Parameters(Kcbini=0.15, Kcbmid=1.15, Kcbend=0.15, Lini=30, Ldev=40, Lmid=50, Lend=30,
                         hini=0.05, hmax=2.0, Zrini=0.15, Zrmax=1.2, pbase=0.55)
    model = fao.Model(index[0], index[-1], par, wth)
    with contextlib.redirect_stdout(io.StringIO()):
        model.run()
    o = model.odata
    # pyfao56 counts the planting day as day 0, this model (like FAO-56) as day 1,
    # so pyfao56's row i holds this model's values of the day before.
    for i in range(1, n):
        c = MAIZE.daily(start + timedelta(days=i - 1), 3.0, 25.0)
        assert c["kcb"] == pytest.approx(o["Kcb"].iloc[i], abs=1e-6)
        assert c["h"] == pytest.approx(o["h"].iloc[i], abs=1e-6)
        assert c["zr"] == pytest.approx(o["Zr"].iloc[i], abs=1e-6)
        # pyfao56 converts the 2 m wind to 2 m once more (factor 1.0002), hence 1e-4.
        assert c["kc_max"] == pytest.approx(o["Kcmax"].iloc[i], abs=1e-4)
        assert c["fc"] == pytest.approx(o["fc"].iloc[i], abs=1e-4)
