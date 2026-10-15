from fao_model.calibration import Calibrator, LITERATURE_COEFFS
from fao_model.simulation import Simulation


def test_daily_steps_produce_all_series():
    sim = Simulation(seed=5, forecaster="baseline")
    for _ in range(40):
        data, results, forecasts = sim.step()
    assert sim.step_number == 40
    assert forecasts["issued"] == data["timestamp"]     # the baseline forecasts every day
    for key in ("T_max", "precip_mean", "ETo"):
        assert len(forecasts[key]) == 14, key
    for key, series in sim.history.items():
        assert len(series) == 40, key
    assert 0.0 < results["ET_Thornthwaite"] < 15.0
    assert 0.0 <= results["ET_crop_actual"] <= results["ET_crop_potential"] + 1e-9


def test_calibration_only_uses_earlier_days():
    sim = Simulation(seed=5, forecaster="baseline")
    for _ in range(Calibrator.MIN_ROWS):
        sim.step()
    # Day MIN_ROWS was computed with MIN_ROWS - 1 earlier days: still literature values.
    assert sim.calib.coeffs == LITERATURE_COEFFS
    sim.step()
    assert sim.calib.coeffs != LITERATURE_COEFFS


def test_replays_real_daily_records():
    from fao_model import data
    _, records = data.load()
    sim = Simulation(records[:400], forecaster="baseline")     # 400 days cover every month
    for _ in range(400):
        day, results, _ = sim.step()
    assert day["timestamp"] == records[399]["timestamp"]
    assert 40.0 < sim.heat_index < 90.0
    assert sim.calib.coeffs != LITERATURE_COEFFS
    assert abs(sim.richards.mass_balance_error) < 1e-9
    assert sim.richards is sim.field.soil
    assert sum(sim.history["ET_crop_actual"]) > 0.0
    # Records are copied, not modified in place.
    assert "T_30d" not in records[0]


def test_stops_at_the_end_of_the_records():
    import pytest
    from fao_model import data
    _, records = data.load()
    from fao_model.crop import Crop
    from fao_model.forecast import Climatology, add_eto
    sim = Simulation(records[:3], heat_index=60.0, crop=Crop(), forecaster="baseline",
                     climatology=Climatology(add_eto(records)))
    for _ in range(3):
        sim.step()
    with pytest.raises(StopIteration):
        sim.step()


def test_heat_index_needs_a_full_year_of_records():
    import pytest
    from fao_model import data
    _, records = data.load()
    with pytest.raises(ValueError):
        Simulation(records[:100], forecaster="baseline")


def test_forecasts_are_issued_every_n_days():
    sim = Simulation(seed=5, forecaster="baseline", forecast_every=7)
    issued = {sim.step()[2]["issued"] for _ in range(15)}
    assert len(issued) == 3      # days 1, 8 and 15


def test_dashboard_with_timesfm():
    import pytest
    from test_forecast import _timesfm_cached
    if not _timesfm_cached():
        pytest.skip("TimesFM or its downloaded weights are not available")
    from fao_model import data
    _, records = data.load()
    sim = Simulation(records[:400], forecaster="timesfm")
    assert sim.forecast_every == 7
    for _ in range(10):
        _, _, forecast = sim.step()
    assert "ETo_q10" in forecast and len(forecast["ETo"]) == 14


def test_replay_from_a_start_date_with_weather_model_forecasts():
    from datetime import date
    from fao_model import data
    _, records = data.load_all()
    sim = Simulation(records, forecaster="ecmwf", start=date(2025, 6, 1))
    day, _, forecast = sim.step()
    assert day["timestamp"].date() == date(2025, 6, 1)
    assert sim.forecast_every == 1 and len(forecast["ETo"]) == 7
    assert sim.forecaster.last_source == "ECMWF IFS 0.25°"
