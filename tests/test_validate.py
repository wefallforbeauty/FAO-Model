from datetime import datetime

import numpy as np
import pytest

pytest.importorskip("matplotlib")

from fao_model import data, validate   # noqa: E402


def test_prepare_has_no_look_ahead():
    _, records = data.load()
    split = datetime(2015, 1, 1)
    prepared, heat_index = validate.prepare(records, split)
    assert 40.0 < heat_index < 90.0
    # Trailing 30-day mean: day 30 averages days 1-30 only.
    assert prepared[29]["T_30d"] == pytest.approx(np.mean([d["T_mean"] for d in records[:30]]))
    # Heat index ignores the test period: changing it must not move I.
    shifted = [dict(d, T_mean=d["T_mean"] + 5.0) if d["timestamp"] >= split else d for d in records]
    assert validate.prepare(shifted, split)[1] == pytest.approx(heat_index)


def test_report_is_written(tmp_path):
    validate.main(["--out", str(tmp_path)])
    report = (tmp_path / "README.md").read_text(encoding="utf-8")
    for heading in ("## 1.", "## 2.", "## 3.", "## 4.", "## 5.", "## 6."):
        assert heading in report
    for name in ("scatter_literature", "scatter_calibrated", "relative_bias",
                 "monthly_climatology", "thornthwaite_monthly"):
        assert (tmp_path / f"{name}.png").stat().st_size > 10_000
