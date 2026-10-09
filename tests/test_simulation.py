from fao_model.calibration import Calibrator, LITERATURE_COEFFS
from fao_model.simulation import Simulation


def test_daily_steps_produce_all_series():
    sim = Simulation(seed=5)
    for _ in range(40):
        data, results, forecasts = sim.step()
    assert sim.step_number == 40
    assert set(forecasts) == set(results)
    for key, series in sim.history.items():
        assert len(series) == 40, key
    assert 0.0 < results["ET_Thornthwaite"] < 15.0


def test_calibration_only_uses_earlier_days():
    sim = Simulation(seed=5)
    for _ in range(Calibrator.MIN_ROWS):
        sim.step()
    # Day MIN_ROWS was computed with MIN_ROWS - 1 earlier days: still literature values.
    assert sim.calib.coeffs == LITERATURE_COEFFS
    sim.step()
    assert sim.calib.coeffs != LITERATURE_COEFFS


def test_replays_real_daily_records():
    from fao_model import data
    _, records = data.load()
    sim = Simulation(records[:400])     # 400 days cover every calendar month
    for _ in range(400):
        day, results, _ = sim.step()
    assert day["timestamp"] == records[399]["timestamp"]
    assert 40.0 < sim.heat_index < 90.0
    assert sim.calib.coeffs != LITERATURE_COEFFS
    assert abs(sim.richards.mass_balance_error) < 1e-9
    # Records are copied, not modified in place.
    assert "T_30d" not in records[0]


def test_stops_at_the_end_of_the_records():
    import pytest
    from fao_model import data
    _, records = data.load()
    sim = Simulation(records[:3], heat_index=60.0)
    for _ in range(3):
        sim.step()
    with pytest.raises(StopIteration):
        sim.step()


def test_heat_index_needs_a_full_year_of_records():
    import pytest
    from fao_model import data
    _, records = data.load()
    with pytest.raises(ValueError):
        Simulation(records[:100])
