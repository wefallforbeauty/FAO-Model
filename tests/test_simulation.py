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
