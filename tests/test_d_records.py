"""S13 of Implementation.md D.9: the V3 records of D-G14, D-G15 and D-G17 read on the Mac after the V3 `pull.sh`.

The run folder is runs/$OBJDEL_V3_RUN. Without the variable both tests skip (no V3 run pulled); with it a missing file
or key fails. The records are written by D7 step 5 (objdel.train.fast.FastGRPO.log and .train, final_eval).
"""
import json
import os
from pathlib import Path

import pytest

pytestmark = pytest.mark.data
RUNS = Path(__file__).resolve().parents[1] / "runs"


def _run() -> Path:
    name = os.environ.get("OBJDEL_V3_RUN")
    if not name:
        pytest.skip("no V3 run pulled (OBJDEL_V3_RUN unset)")
    return RUNS / name


def test_b04():
    """B04 (D-G14, D-G15): B3 run steps 2-6 make no sync allocation and no pinned growth over the run; per view of
    the B3 val eval, eval_scratch_peak <= max_scratch_bytes and live_device_bytes <= arena_bytes + table_bytes."""
    run = _run()
    log = json.loads((run / "b3-v5" / "rl.json").read_text())["log"]
    steps = {r["step"]: r for r in log if "engine/sync_alloc" in r}
    for key in ("engine/sync_alloc", "engine/sync_free"):
        assert steps[6][key] - steps[1][key] == 0, key
    assert [r["engine/pinned_grows"] for r in log if set(r) == {"engine/pinned_grows", "step"}][-1] == 0
    e = json.loads((run / "eval_b3_val_shifted" / "final_eval.json").read_text())["engine"]
    assert e["eval_scratch_peak"] <= e["max_scratch_bytes"]
    assert e["live_device_bytes"] <= e["arena_bytes"] + e["table_bytes"]


def test_b06():
    """B06 (D-G17 (a), (c)): per view the device process count never rises above its count before scoring (and the
    /proc count is >= 1 before), and peak memory.used <= 30 GB."""
    topo = json.loads((_run() / "d" / "topology.json").read_text())
    assert topo
    for view, r in topo.items():
        assert r["method"] in ("nvidia-smi", "proc"), view
        assert r["before"] >= 1 and r["max_during"] <= r["before"], view
        assert r["peak_memory_used_mib"] <= 30e9 / 2 ** 20, view  # 30 GB in MiB
