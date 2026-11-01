"""V1 readiness of train_sft and the R0 stop rule: deterministic mode without full_determinism (PLAN D24), the run
manifest's timing block, C2's depth read from R2's g_probe.json (branches_from), runner-only config keys, and
project_r1 (§4.7 step 9, §4.9: stop if projected R1 > 2x the high estimate)."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "experiments/phase2/configs"
C2 = sorted(CONFIGS.glob("c2_*.yaml")) + [CONFIGS / "rx8.yaml"]


def _yaml(path):
    return pytest.importorskip("yaml").safe_load(path.read_text())


@pytest.fixture
def TS():
    pytest.importorskip("trl")
    from experiments.phase2 import train_sft

    return train_sft


def _view(tmp_path, epochs):
    (tmp_path / "v").mkdir()
    (tmp_path / "v/manifest.json").write_text(json.dumps({"config": {"epochs": epochs}}))
    return str(tmp_path / "v")


# ---------------------------------------------------------------- 1. deterministic mode (D24)


def test_timing_arms_use_deterministic_not_full_determinism():
    for p in CONFIGS.glob("*.yaml"):
        assert "full_determinism" not in _yaml(p).get("sft", {}), p.name  # it sets CUDA_LAUNCH_BLOCKING=1
    for p in CONFIGS.glob("r0_timing_*.yaml"):
        assert _yaml(p).get("deterministic", False) is p.stem.endswith("_det_on"), p.name


def test_deterministic_is_applied_before_anything_loads(TS, tmp_path, monkeypatch):
    import torch

    monkeypatch.delenv("CUBLAS_WORKSPACE_CONFIG", raising=False)
    a = TS.RunArgs("lfm2.5-vl-450m", _view(tmp_path, 0), str(tmp_path / "o"), vision=False, deterministic=True)
    try:
        with pytest.raises(ValueError, match="view epoch"):
            TS.train(a)
        assert os.environ["CUBLAS_WORKSPACE_CONFIG"] == ":4096:8"
        assert torch.are_deterministic_algorithms_enabled() and torch.is_deterministic_algorithms_warn_only_enabled()
        monkeypatch.setenv("CUBLAS_WORKSPACE_CONFIG", ":16:8")  # an explicit setting wins
        with pytest.raises(ValueError, match="view epoch"):
            TS.train(TS.RunArgs("m", a.view, str(tmp_path / "o"), vision=False, deterministic=True, run_dir=""))
        assert os.environ["CUBLAS_WORKSPACE_CONFIG"] == ":16:8"
    finally:
        torch.use_deterministic_algorithms(False)
    monkeypatch.delenv("CUBLAS_WORKSPACE_CONFIG")
    with pytest.raises(ValueError, match="view epoch"):
        TS.train(TS.RunArgs("m", a.view, str(tmp_path / "o"), vision=False, run_dir=""))
    assert "CUBLAS_WORKSPACE_CONFIG" not in os.environ and not torch.are_deterministic_algorithms_enabled()


# ---------------------------------------------------------------- 2. timing block


def test_timing_excludes_the_warmup_steps(TS, monkeypatch):
    import torch

    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    t = TS.timing([10.0, 12.0, 14.0, 15.0, 16.0, 17.0], 32)  # step-end clocks: steps 4-6 take 1 s each
    assert t == {"steps": 6, "warmup_excluded": 3, "s_per_step": 1.0, "samples_per_s": 32.0, "batch": 32,
                 "max_memory_gib": None}  # fmt: skip
    assert TS.timing([1.0, 2.0, 3.0], 32)["s_per_step"] is None  # nothing after the warm-up
    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setattr(torch.cuda, "max_memory_allocated", lambda: 3 * 2**30)
    assert TS.timing([], 32) == {
        "steps": 0, "warmup_excluded": 3, "s_per_step": None, "samples_per_s": None, "batch": 32, "max_memory_gib": 3.0
    }  # fmt: skip


@pytest.mark.slow
def test_manifest_records_timing(TS, tmp_path, synthetic_items, monkeypatch):
    from test_p2_scripts import _tiny_vlm

    from experiments.phase2 import build_views as BV
    from objdel.data import store

    store.write(synthetic_items(("train",), unanimous=True), tmp_path / "store")
    BV.build(tmp_path / "store", tmp_path / "view", "lidc", "train", epochs=2, seed=3, procs=1)
    sft = {"per_device_train_batch_size": 1, "gradient_accumulation_steps": 1, "use_cpu": True, "bf16": False,
           "dataloader_num_workers": 0, "dataloader_persistent_workers": False,
           "dataloader_prefetch_factor": None}  # fmt: skip
    a = TS.RunArgs("lfm2.5-vl-450m", str(tmp_path / "view"), str(tmp_path / "o"), branches=[], steps_per_epoch=6,
                   model_path=str(_tiny_vlm(tmp_path / "m")), vision=False, rolling_minutes=1e3, sft=sft)  # fmt: skip
    assert TS.train(a) == 0
    t = json.loads((tmp_path / "o/manifest.json").read_text())["timing"]
    assert t["steps"] == 6 and t["s_per_step"] > 0 and t["samples_per_s"] == pytest.approx(1 / t["s_per_step"])
    assert t["batch"] == 1 and t["rows_per_epoch"] > 6 and t["max_memory_gib"] is None


# ---------------------------------------------------------------- 4. C2 depth from R2 (branches_from)


def test_c2_configs_take_their_depth_from_g_probe():
    for p in C2:
        c = _yaml(p)
        assert c["branches_from"] == "r1-lfm/g_probe/g_probe.json" and "branches" not in c, p.name
        assert c["stages"] == ["train", "topology"], p.name


def _probe(tmp_path, **doc):
    (tmp_path / "r1-lfm/g_probe").mkdir(parents=True)
    (tmp_path / "r1-lfm/g_probe/g_probe.json").write_text(json.dumps(doc))
    return "r1-lfm/g_probe/g_probe.json"


def test_branches_from_reads_the_final_verdict_checkpoint(TS, tmp_path):
    ckpt = "/workspace/runs/v1/r1-lfm/branches/e0.5"  # D12(c): the repair branch e0.5 got the final verdict
    rel = _probe(tmp_path, verdict="INCONCLUSIVE", ckpt=ckpt, probed=["e1", "e0.5"])
    a = TS.RunArgs("lfm2.5-vl-450m", _view(tmp_path, 0), "c2", run_dir=str(tmp_path), branches_from=rel, vision=False)
    with pytest.raises(ValueError, match="decays to 0.575 epochs"):  # resolved to [0.5] before the view check
        TS.train(a)
    assert a.branches == [0.5]


@pytest.mark.parametrize(
    "doc", [None, {"verdict": None, "skipped": "gates"}, {"ckpt": "b/e1"}, {"verdict": "PASS", "ckpt": "c"}]
)
def test_branches_from_refuses_without_a_verdict(TS, tmp_path, doc):
    rel = _probe(tmp_path, **doc) if doc else "r1-lfm/g_probe/g_probe.json"
    a = TS.RunArgs("lfm2.5-vl-450m", _view(tmp_path, 3), "c2", run_dir=str(tmp_path), branches_from=rel, vision=False)
    with pytest.raises(SystemExit, match="C2 must not start"):
        TS.train(a)
    assert not (tmp_path / "c2").exists()


@pytest.mark.parametrize("path", C2 + [CONFIGS / "r1_lfm.yaml"], ids=lambda p: p.stem)
def test_runner_only_keys_parse(TS, path, tmp_path):
    from trl import TrlParser

    extra = {"sampler": "hf", "stages": ["train"]}
    (tmp_path / "c.yaml").write_text(json.dumps(_yaml(path) | extra))  # JSON is YAML
    (a,) = TrlParser([TS.RunArgs]).parse_args_and_config(["--config", str(tmp_path / "c.yaml")])
    assert a.sampler == "hf" and a.stages == ["train"]


# ---------------------------------------------------------------- 3. project_r1 (stop rule)


def _manifest(tmp_path, name, s_per_step, model="lfm2.5-vl-450m", workers=6, det=False, batch=32):
    d = tmp_path / name
    d.mkdir()
    cfg = {"model": model, "output_dir": str(d), "deterministic": det, "sft": {"dataloader_num_workers": workers}}
    t = {"s_per_step": s_per_step, "batch": batch, "rows_per_epoch": 22400}
    (d / "manifest.json").write_text(json.dumps({"config": cfg, "timing": t}))
    return str(d / "manifest.json")


def test_project_r1_go_and_stop(tmp_path, capsys):
    pytest.importorskip("trl")
    from experiments.phase2 import project_r1 as P

    off = _manifest(tmp_path, "r0-timing-w6-det-off", 2.0)
    on = _manifest(tmp_path, "r0-timing-w6-det-on", 9.0, det=True)
    w0 = _manifest(tmp_path, "r0-timing-w0-det-off", 9.0, workers=0)
    assert P.main(["--model", "lfm2.5-vl-450m", on, off, w0]) == 0  # basis: det-off workers 6
    out = capsys.readouterr().out.splitlines()
    doc = json.loads(out[1])
    assert out[0].startswith("GO") and doc["basis"] == off and doc["r1_steps"] == 1400 + 52 + 105 + 210
    assert doc["projected_h"] == pytest.approx(1767 * 2.0 / 3600 + 0.3) and doc["stop_h"] == 3.4
    assert doc["arms"][on] == 9.0
    slow = _manifest(tmp_path, "slow", 7.0)  # 1767 x 7 s = 3.44 h + 0.3 > 3.4
    assert P.main(["--model", "lfm2.5-vl-450m", slow]) == 1 and capsys.readouterr().out.startswith("STOP")
    qwen = _manifest(tmp_path, "r0-qwen35-arm", 7.0, model="qwen3.5-0.8b")  # 3.74 h <= 6.8 h
    assert P.main(["--model", "qwen3.5-0.8b", qwen, off]) == 0
    assert json.loads(capsys.readouterr().out.splitlines()[1])["stop_h"] == 6.8


def test_project_r1_refuses_an_unusable_basis(tmp_path):
    pytest.importorskip("trl")
    from experiments.phase2 import project_r1 as P

    for bad in (_manifest(tmp_path, "w0", 2.0, workers=0), _manifest(tmp_path, "b16", 2.0, batch=16),
                _manifest(tmp_path, "none", None)):  # fmt: skip
        with pytest.raises(SystemExit):
            P.main(["--model", "lfm2.5-vl-450m", bad])


def test_project_r1_cli(tmp_path):
    pytest.importorskip("trl")
    m = _manifest(tmp_path, "r0-timing-w6-det-off", 2.0)
    r = subprocess.run([sys.executable, "-I", "-m", "experiments.phase2.project_r1", "--model", "lfm2.5-vl-450m", m],
                       cwd=ROOT, capture_output=True, text=True, timeout=300)  # fmt: skip
    assert r.returncode == 0 and r.stdout.startswith("GO"), r.stderr[-2000:]
