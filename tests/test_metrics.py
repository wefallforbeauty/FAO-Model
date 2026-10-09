import numpy as np
import pytest

from fao_model.metrics import compare


def test_perfect_agreement():
    ref = np.array([1.0, 2.0, 4.0, 7.0])
    m = compare(ref, ref)
    assert m["bias"] == 0.0 and m["rmse"] == 0.0 and m["mae"] == 0.0
    assert m["r"] == pytest.approx(1.0)
    assert m["nse"] == pytest.approx(1.0)
    assert m["slope"] == pytest.approx(1.0)
    assert m["intercept"] == pytest.approx(0.0, abs=1e-12)


def test_scaled_and_offset_estimates():
    ref = np.array([1.0, 2.0, 3.0, 6.0])     # mean 3
    m = compare(1.1 * ref, ref)
    assert m["rel_bias"] == pytest.approx(10.0)
    assert m["slope"] == pytest.approx(1.1)
    m = compare(ref + 0.5, ref)
    assert m["bias"] == pytest.approx(0.5)
    assert m["rmse"] == pytest.approx(0.5)
    assert m["rel_bias"] == pytest.approx(100.0 * 0.5 / 3.0)
    assert m["r"] == pytest.approx(1.0)
    # NSE = 1 - sum(0.5^2) / sum((ref - 3)^2) = 1 - 1.0 / 14
    assert m["nse"] == pytest.approx(1.0 - 1.0 / 14.0)
