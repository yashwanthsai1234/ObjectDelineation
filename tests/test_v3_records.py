"""S12 of Implementation.md D.9 (D7 step 5): the V3 records are written, on the Mac with fakes.

- `FastGRPO.log` on a CUDA build writes `engine/<key>` for every `gc.stats()` key and `engine/max_memory_reserved`
  on steps 1-6 only; `FastGRPO.train` appends `engine/pinned_grows` once at the end of the run (D-G14, D-G15).
- `final_eval` resets the engine counters once before each view's scoring and writes `final_eval.json["engine"]`;
  with `OBJDEL_RUN_DIR` set on a CUDA build its sampler thread writes `d/topology.json` (D-G17 (a), (c)).
"""
import json
import os
import stat
from collections import defaultdict

import pytest

pytest.importorskip("trl")

from transformers import Trainer  # noqa: E402
from trl import GRPOTrainer  # noqa: E402

from experiments.phase2 import final_eval as FE  # noqa: E402
from objdel.geom import engine as _gpucheck  # noqa: E402
from objdel.policy import caches as policy_caches  # noqa: E402
from objdel.rl import trl_adapter as rl_trl_adapter  # noqa: E402

STATS = {"arena_bytes": 11, "sync_alloc": 2, "sync_free": 2, "pool_allocs": 5, "pool_frees": 4, "pinned_grows": 0,
         "table_bytes": 7, "eval_scratch_peak": 3, "live_device_bytes": 18}


def _trainer(step: int):
    t = object.__new__(rl_trl_adapter.FastGRPO)
    t.model = type("M", (), {"training": True})()
    t._metrics = {"train": defaultdict(list), "eval": defaultdict(list)}
    t.state = type("S", (), {"global_step": step, "log_history": []})()
    return t


@pytest.fixture
def cuda_fakes(monkeypatch):
    monkeypatch.setattr(_gpucheck, "build", lambda: "cuda")
    monkeypatch.setattr(_gpucheck, "stats", lambda: dict(STATS))
    monkeypatch.setattr(policy_caches.torch.cuda, "max_memory_reserved", lambda *a, **k: 12345)


@pytest.mark.usefixtures("cuda_fakes")
def test_fast_grpo_logs_engine_counters_on_steps_1_to_6(monkeypatch):
    logged = []
    monkeypatch.setattr(GRPOTrainer, "log", lambda self, logs, start_time=None: logged.append(dict(logs)))
    want = {f"engine/{k}": v for k, v in STATS.items()} | {"engine/max_memory_reserved": 12345}
    for step in range(0, 9):
        _trainer(step).log({"loss": 1.0})
    for step, logs in enumerate(logged):
        engine = {k: v for k, v in logs.items() if k.startswith("engine/")}
        assert engine == (want if 1 <= step <= 6 else {}), step
        assert logs["loss"] == 1.0
    t = _trainer(3)
    t.model.training = False  # an eval log carries no engine keys
    t.log({})
    assert not any(k.startswith("engine/") for k in logged[-1])


@pytest.mark.usefixtures("cuda_fakes")
def test_fast_grpo_logs_pinned_grows_at_the_end_of_the_run(monkeypatch):
    monkeypatch.setattr(Trainer, "train", lambda self, *a, **k: "out")
    t = _trainer(40)
    assert t.train() == "out"
    assert t.state.log_history == [{"engine/pinned_grows": 0, "step": 40}]


def test_fast_grpo_logs_no_engine_keys_on_a_cpu_build(monkeypatch):
    logged = []
    monkeypatch.setattr(GRPOTrainer, "log", lambda self, logs, start_time=None: logged.append(dict(logs)))
    monkeypatch.setattr(Trainer, "train", lambda self, *a, **k: None)
    t = _trainer(1)
    t.log({})
    t.train()
    assert not any(k.startswith("engine/") for k in logged[0]) and t.state.log_history == []


# ------------------------------------------------------------------------------------------------ final_eval

ENGINE_KEYS = {"eval_sets_s", "max_scratch_bytes", "arena_bytes", "table_bytes", "eval_scratch_peak",
               "live_device_bytes", "max_memory_reserved"}


def _fake_smi(tmp_path, apps: str, used: str):
    """An nvidia-smi on PATH that prints fixed CSV rows."""
    d = tmp_path / "bin"
    d.mkdir(exist_ok=True)
    f = d / "nvidia-smi"
    f.write_text(f"#!/bin/sh\ncase \"$1\" in --query-compute-apps*) printf '{apps}';; *) printf '{used}';; esac\n")
    f.chmod(f.stat().st_mode | stat.S_IEXEC)
    return d


def _two_views(ev, tmp_path, monkeypatch):
    import shutil

    from tests.test_final_eval import fake_sampler

    v2 = shutil.copytree(ev, tmp_path / "lidc-val-two")
    monkeypatch.setattr(FE, "load_sampler", fake_sampler(ev))
    return ["--ckpt", "c", "--view", str(ev), "--split", "val", "--view", str(v2), "--split",
            "val", "--out", str(tmp_path / "O")], [ev.name, v2.name]


@pytest.fixture
def ev(tmp_path, clean_store, synthetic_items):
    from tests.test_ceiling import _view

    return _view(tmp_path / "view", synthetic_items)


def test_final_eval_writes_engine_records_and_topology(ev, tmp_path, monkeypatch):
    """Per view: one reset_stats before scoring, the engine record in its final_eval.json, and on a CUDA build with
    OBJDEL_RUN_DIR the sampler's topology record (method nvidia-smi, count before, largest during, peak MiB)."""
    argv, names = _two_views(ev, tmp_path, monkeypatch)
    resets = []
    real = _gpucheck.reset_stats
    monkeypatch.setattr(_gpucheck, "reset_stats", lambda: resets.append(1) or real())
    monkeypatch.setattr(_gpucheck, "build", lambda: "cuda")
    monkeypatch.setattr(FE, "TOPOLOGY_PERIOD_S", 0.01)
    bin_ = _fake_smi(tmp_path, "123, 900 MiB\\n456, 800 MiB\\n", "2048 MiB\\n")
    monkeypatch.setenv("PATH", f"{bin_}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("OBJDEL_RUN_DIR", str(tmp_path / "run"))
    FE.main(argv)
    assert len(resets) == 2
    for name in names:
        engine = json.loads((tmp_path / "O" / name / "final_eval.json").read_text())["engine"]
        assert engine.keys() == ENGINE_KEYS
        assert engine["max_scratch_bytes"] == FE.device.MAX_SCRATCH_BYTES and engine["eval_sets_s"] > 0
        assert 0 < engine["eval_scratch_peak"] <= engine["max_scratch_bytes"]
        assert engine["live_device_bytes"] <= engine["arena_bytes"] + engine["table_bytes"]
    topo = json.loads((tmp_path / "run" / "d" / "topology.json").read_text())
    assert topo.keys() == set(names)
    for rec in topo.values():
        assert rec == {"method": "nvidia-smi", "before": 2, "max_during": 2, "peak_memory_used_mib": 2048}


def test_final_eval_topology_falls_back_to_proc_and_needs_the_run_dir(ev, tmp_path, monkeypatch):
    """An empty compute-apps list: method proc (the container's /proc view); OBJDEL_RUN_DIR unset or a CPU build:
    no sampler thread and no topology record."""
    argv, names = _two_views(ev, tmp_path, monkeypatch)
    monkeypatch.setattr(_gpucheck, "build", lambda: "cuda")
    monkeypatch.setattr(FE, "TOPOLOGY_PERIOD_S", 0.01)
    monkeypatch.setattr(FE, "_proc_count", lambda: 1)
    monkeypatch.setenv("PATH", f"{_fake_smi(tmp_path, '', '10 MiB')}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("OBJDEL_RUN_DIR", str(tmp_path / "run"))
    FE.main(argv)
    topo = json.loads((tmp_path / "run" / "d" / "topology.json").read_text())
    assert {r["method"] for r in topo.values()} == {"proc"} and all(r["before"] == 1 for r in topo.values())
    started = []
    monkeypatch.setattr(FE, "_smi", lambda q: started.append(q))  # the sampler's first act (none: no thread)
    monkeypatch.delenv("OBJDEL_RUN_DIR")
    with FE.topology("x"):
        pass
    monkeypatch.setattr(_gpucheck, "build", lambda: "cpu")
    monkeypatch.setenv("OBJDEL_RUN_DIR", str(tmp_path / "run2"))
    with FE.topology("x"):
        pass
    assert not started and not (tmp_path / "run2").exists()
