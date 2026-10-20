"""Every Phase 1 script has a command line: `--help` prints usage and exits 0 without touching data (B4-7).

Scripts run as modules in an isolated interpreter, as CLAUDE.md prescribes for anything that reads data/
(`python -I -m experiments.phase1.<script>`; -I drops the current directory from sys.path, so this also checks that
the editable install makes `experiments` importable). One `-I` subprocess runs every script through `runpy`, so the
heavy imports (torch, TRL) are paid once. A script whose optional extra is not installed (the CI `python` job has
only `[dev]`) is skipped, not failed.

References: https://docs.python.org/3/using/cmdline.html#cmdoption-I,
https://docs.python.org/3/library/runpy.html#runpy.run_module, https://docs.python.org/3/library/argparse.html
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = sorted(p.stem for p in (ROOT / "experiments" / "phase1").glob("*.py") if p.stem != "__init__")
OPTIONAL = {"torch", "transformers", "trl", "datasets", "peft", "accelerate", "pydicom", "pylidc"}  # extras

RUNNER = """
import contextlib, io, json, runpy, sys
path, scripts, out = sys.argv[1], sys.argv[2:], {}
for s in scripts:
    buf, sys.argv = io.StringIO(), [s, "--help"]
    try:
        with contextlib.redirect_stdout(buf):
            runpy.run_module("experiments.phase1." + s, run_name="__main__", alter_sys=True)
        out[s] = {"code": "returned without exiting", "stdout": buf.getvalue()}
    except SystemExit as e:
        out[s] = {"code": e.code or 0, "stdout": buf.getvalue()}
    except ModuleNotFoundError as e:
        out[s] = {"missing": e.name}
    except Exception as e:
        out[s] = {"code": repr(e), "stdout": buf.getvalue()}
open(path, "w").write(json.dumps(out))
"""


@pytest.fixture(scope="module")
def helps(tmp_path_factory):
    out = tmp_path_factory.mktemp("help") / "help.json"
    r = subprocess.run([sys.executable, "-I", "-c", RUNNER, str(out), *SCRIPTS], cwd=ROOT, capture_output=True,
                       text=True, timeout=600)
    assert r.returncode == 0 and out.is_file(), r.stderr[-3000:]
    return json.loads(out.read_text())


def test_phase1_scripts_found():
    assert {"census_sizerule", "eval_run", "prepare_rl"} <= set(SCRIPTS)


@pytest.mark.parametrize("script", SCRIPTS)
def test_phase1_script_help(helps, script):
    h = helps[script]
    if "missing" in h:
        if h["missing"].split(".")[0] in OPTIONAL:
            pytest.skip(f"optional dependency {h['missing']} not installed")
        pytest.fail(f"{script}: missing module {h['missing']}")
    assert h["code"] == 0, h
    assert h["stdout"].startswith("usage:")


@pytest.mark.parametrize("cuda, supported, expect", [(False, True, False), (True, False, False), (True, True, True)])
def test_sample_eval_uses_bf16_only_on_cuda_cards_that_support_it(cuda, supported, expect, monkeypatch):
    """A pre-Ampere card (no bf16) must not get bf16=True (TrainingArguments raises); the Mac (no CUDA) samples in
    fp32, as F3 verified the sampler."""
    torch = pytest.importorskip("torch")
    from experiments.phase1 import sample_eval
    monkeypatch.setattr(torch.cuda, "is_available", lambda: cuda)
    monkeypatch.setattr(torch.cuda, "is_bf16_supported", lambda *a, **k: supported)
    assert sample_eval.bf16() is expect


def test_smoke_grpo_group_stats_use_the_gate_definitions():
    """REVIEW #10: smoke.json's GRPO-signal numbers are the Phase 2 gates' (objdel.eval: duplicates = equal rasters,
    sigma of the legality-gated geometric terms), not text duplicates and ungated sigma."""
    pytest.importorskip("trl")
    import numpy as np

    from experiments.phase1 import smoke_grpo
    from objdel.eval import metrics, summarize
    from objdel.reward import ExpertCache

    def square(x0, side=400, dx=0):  # grid units on a 50 px frame (20 per px); dx nudges x within a pixel
        return json.dumps({"nodule": [{"outer": [[x0 + dx, x0], [x0 + dx + side, x0], [x0 + dx + side, x0 + side],
                                                 [x0 + dx, x0 + side]], "holes": []}]}, separators=(",", ":"))

    m = np.zeros((3, 50, 50), bool)
    for i, (a, b) in enumerate([(8, 28), (8, 28), (8, 30)]):
        m[i, a:b, a:b] = True
    caches = {"lidc/a": {"nodule": ExpertCache.build(m)}, "lidc/b": {"nodule": ExpertCache.build(m)}}
    same = [square(160)] * 4 + [square(160, dx=d) for d in (1, 2, 3, 4)]  # 5 texts, one raster
    spread = [square(160 + 20 * i) for i in range(8)]  # 8 rasters
    records = [("lidc/a", t, False, 50, 50, ["nodule"]) for t in same]
    records += [("lidc/b", t, False, 50, 50, ["nodule"]) for t in spread]
    adv = np.r_[np.zeros(8), np.linspace(-1, 1, 8)]
    s = smoke_grpo.group_stats(records, caches.__getitem__, adv)
    assert s["duplicate_rate"] == pytest.approx((7 / 8 + 0) / 2)  # text duplicates would give (3 / 8 + 0) / 2
    items = [{"item_id": k, "dataset": "lidc", **metrics(ts, {"nodule": m}, 50, 50, ["nodule"], [False] * 8, caches[k])}
             for k, ts in (("lidc/a", same), ("lidc/b", spread))]
    gate = summarize(items)["datasets"]["lidc"]
    assert s["duplicate_rate"] == gate["dup_rate"] and s["sigma_above_floor_any_term"] == gate["signal_share"]
    assert s["groups"] == 2 and s["nonzero_adv_groups"] == 1 and s["adv_abs_max"] == 1.0
