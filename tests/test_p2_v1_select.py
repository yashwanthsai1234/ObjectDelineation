"""V1 GPU-session readiness: g_probe stops when R1's §5.1 gates failed, select_ckpt --only topology (C2), the R0
parity driver (experiments/phase2/r0_parity.py, PLAN §4.7 step 7 / G-parity) and the prompt-id echo it needs from
sample_vllm. No GPU, no vLLM: the sampler child and the model scoring are replaced."""

import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("trl")

from experiments.phase2 import g_probe as GP  # noqa: E402
from experiments.phase2 import parity  # noqa: E402
from experiments.phase2 import sample_vllm as V  # noqa: E402
from experiments.phase2 import select_ckpt as SC  # noqa: E402
from objdel.train import serve  # noqa: E402


def _boom(*a, **k):
    raise AssertionError("must not be called")


# ---------------------------------------------------------------- g_probe after failed gates


def _selected(tmp_path, **kw) -> list[str]:
    b = lambda d, e: {"depth": d, "eligible": e}  # noqa: E731
    sel = {"run": str(tmp_path / "run"), "selected": "e2", "gates_branch": "e2", "gates_pass": True,
           "gates_failed": [], "branches": {"e0.5": b(0.5, True), "e1": b(1.0, False), "e2": b(2.0, True)}} | kw
    (tmp_path / "sel").mkdir(exist_ok=True)
    (tmp_path / "sel/selected.json").write_text(json.dumps(sel))
    out = tmp_path / "gp"
    return ["--model", "lfm2.5-vl-450m", "--select", str(tmp_path / "sel"), "--train", str(tmp_path / "prep"),
            "--probe", str(tmp_path / "shifted"), "--out", str(out)]


@pytest.mark.parametrize(("kw", "why"), [
    ({"gates_pass": False, "gates_failed": ["G-format", "G-topology-a"]}, "§5.1 gates failed: G-format, G-topology-a"),
    ({"selected": None, "gates_pass": False, "gates_failed": ["G-vision"]}, "§5.1 gates failed: G-vision"),
    ({"selected": None, "gates_branch": None, "gates_pass": None, "gates_failed": []},
     "§5.1 gates failed: no eligible branch (D12)"),
])
def test_g_probe_skips_when_the_gates_failed(tmp_path, monkeypatch, kw, why):
    monkeypatch.setattr(serve, "preload", _boom)  # skipped before any cache, pool or model
    monkeypatch.setattr(GP, "probe", _boom)
    assert GP.main(_selected(tmp_path, **kw)) == 0  # the run chain ends cleanly; the JSON says why
    doc = json.loads((tmp_path / "gp/g_probe.json").read_text())
    assert doc == {"verdict": None, "skipped": why, "probed": [], "gates_branch": kw.get("gates_branch", "e2")}


def test_g_probe_with_passing_gates_keeps_the_d12c_repair(tmp_path, monkeypatch):
    from objdel.train import dataset

    monkeypatch.setattr(serve, "preload", lambda *a, **k: 0)
    monkeypatch.setattr(dataset, "grpo_rows", lambda *a, **k: [])
    monkeypatch.setattr(GP, "view_rows", lambda *a, **k: [])
    verdicts = iter(["INCONCLUSIVE", "PASS"])
    monkeypatch.setattr(GP, "probe", lambda a, ckpt, *r: {"verdict": next(verdicts), "h0": 0, "gain": 0,
                                                          "ckpt": str(ckpt)})
    assert GP.main(_selected(tmp_path)) == 0
    doc = json.loads((tmp_path / "gp/g_probe.json").read_text())
    assert doc["verdict"] == "PASS" and doc["probed"] == ["e2", "e0.5"] and "skipped" not in doc


def test_select_records_the_failed_gate_names(tmp_path, monkeypatch):
    from tests.test_p2_breakit_runners import _fake_select

    SC_, base, _ = _fake_select(tmp_path, monkeypatch, eligible=False)
    assert SC_.main(base + ["--vision", '{"prepared": "p", "reference": 1}']) == 0
    doc = json.loads((tmp_path / "sel/selected.json").read_text())
    gates = json.loads((tmp_path / "sel/gates.json").read_text())
    failed = [n for n, g in gates["gates"].items() if not g["pass"]]
    assert doc["gates_pass"] is False and failed and doc["gates_failed"] == failed


@pytest.mark.parametrize("ok", [True, False])
def test_vllm_select_gates_on_the_hf_check(tmp_path, monkeypatch, ok):
    """D13 (owner 2026-11-05): with --sampler vllm the judged branch's HF check is a gate; --sampler hf has none."""
    from tests.test_p2_breakit_runners import _fake_select

    SC_, base, _ = _fake_select(tmp_path, monkeypatch, hf_check=ok)
    assert SC_.main(base + ["--vision", '{"prepared": "p", "reference": 1}']) == 0
    doc = json.loads((tmp_path / "sel/selected.json").read_text())
    gates = json.loads((tmp_path / "sel/gates.json").read_text())
    report = json.loads((tmp_path / "sel/report.json").read_text())
    assert report["hf_check"]["branch"] == doc["gates_branch"] == "e1"
    assert gates["gates"]["G-hf-check"]["pass"] is ok and gates["passed"]["G-hf-check"] is ok
    assert ("G-hf-check" in doc["gates_failed"]) is not ok
    SC_, base, _ = _fake_select(tmp_path / "hf", monkeypatch, hf_check=ok)
    assert SC_.main(base + ["--sampler", "hf", "--vision", '{"prepared": "p", "reference": 1}']) == 0
    hf = json.loads((tmp_path / "hf/sel/selected.json").read_text())
    assert "G-hf-check" not in json.loads((tmp_path / "hf/sel/gates.json").read_text())["gates"]
    assert doc["gates_pass"] is (hf["gates_pass"] and ok)  # the other gates are the same in both runs


@pytest.mark.parametrize(("shift", "ok"), [(0.0, True), (0.3, False)])
def test_hf_check_compares_vllm_and_hf_metrics_per_item(tmp_path, monkeypatch, shift, ok):
    """hf_check: HF samples HF_CHECK_ITEMS seeded items of the branch; each metric's per-item vLLM - HF difference
    needs 0 inside its nodule-cluster bootstrap CI."""
    keys = [f"lidc/s{i}/n000/z{i}" for i in range(SC.HF_CHECK_ITEMS)]
    rng = np.random.default_rng(1)
    vals = {i: rng.uniform(0.2, 0.8, 8) for i in keys}

    def samples(i, d):
        return {k: list(vals[i] + d) for k in SC.HF_CHECK_KEYS}

    seen = {}
    monkeypatch.setattr(SC, "view_rows", lambda view, items=None: seen.setdefault("items", items) and
                        [{"key": i} for i in keys])  # fmt: skip
    monkeypatch.setattr(SC, "sample_view", lambda a, b, rows, out, k: seen.update(sampler=a.sampler, out=out) or {})
    monkeypatch.setattr(SC, "score_items", lambda rows, done, root: (
        [{"item_id": r["key"], "samples": samples(r["key"], 0.0)} for r in rows], []))
    (tmp_path / "e1/shifted").mkdir(parents=True)
    (tmp_path / "e1/shifted/items.jsonl").write_text(
        "".join(json.dumps({"item_id": i, "samples": samples(i, shift)}) + "\n" for i in keys))
    a = SimpleNamespace(shifted=tmp_path / "view", run=tmp_path / "run", k=8, sampler="vllm")
    r = SC.hf_check(a, "e1", tmp_path)
    assert seen == {"items": SC.HF_CHECK_ITEMS, "sampler": "hf", "out": tmp_path / "e1/hf-check"}
    assert r["items"] == SC.HF_CHECK_ITEMS and r["pass"] is ok and set(r["vllm_minus_hf"]) == set(SC.HF_CHECK_KEYS)
    assert r["vllm_minus_hf"]["iou_tc"][0] == pytest.approx(shift)


# ---------------------------------------------------------------- select_ckpt --only topology (C2)


def _c2_run(tmp_path, monkeypatch, names=("e1",)):
    from objdel.eval import probes

    run = tmp_path / "c2-a0p_seed2"
    for n in names:
        (run / "branches" / n).mkdir(parents=True)
    calls = []

    def sample_view(a, branch, rows, out, k):
        calls.append({"branch": branch.name, "out": out, "k": k, "rows": len(rows), "sampler": a.sampler})
        return {r["key"]: {"item_id": r["key"], "samples": ["{}"] * k} for r in rows}

    monkeypatch.setattr(SC, "exited", lambda run: {"config": {"model": "lfm2.5-vl-450m"}})
    monkeypatch.setattr(SC, "view_rows", lambda view, items=None, rl_only=True: [{"key": f"{view.name}/{i}"}
                                                                                 for i in range(3)])
    monkeypatch.setattr(SC, "sample_view", sample_view)
    monkeypatch.setattr(SC, "score_view", _boom)
    monkeypatch.setattr(SC, "vision_check", _boom)
    monkeypatch.setattr(serve, "preload", _boom)
    monkeypatch.setattr(probes, "score_precision", _boom)
    topo = {"bins": [[1.0], [0.5], [0.0]], "false_nonsimple": 0.0, "legal": 1.0, "pass": True}
    monkeypatch.setattr(probes, "score_topology", lambda rows, root: topo | {"rows": len(rows), "root": str(root)})
    base = ["--run", str(run), "--probes", str(tmp_path / "eval_views"), "--prepared", str(tmp_path / "prepared"),
            "--only", "topology"]
    return run, base, calls


@pytest.mark.parametrize("extra", [["--sampler", "hf"], []])
def test_only_topology_samples_and_scores_the_topology_probe_alone(tmp_path, monkeypatch, extra):
    run, base, calls = _c2_run(tmp_path, monkeypatch)
    assert SC.main(base + extra) == 0
    assert [(c["branch"], c["out"], c["k"], c["rows"]) for c in calls] == [
        ("e1", run / "select/e1/probe-topology", 8, 3)]
    assert calls[0]["sampler"] == ("hf" if extra else "vllm")
    rep = json.loads((run / "select/report.json").read_text())
    topo = rep["probes"]["topology"]
    assert set(rep["probes"]) == {"topology"} and rep["branch"] == "e1" and rep["only"] == "topology"
    assert topo["rows"] == 3 and topo["root"] == str(tmp_path / "eval_views/probe-topology")
    assert {"bins", "false_nonsimple", "legal"} <= set(topo)  # what c2_analyze reads per seed
    assert not (run / "select/gates.json").exists()


def test_only_topology_needs_a_single_branch_and_probes(tmp_path, monkeypatch):
    _, base, _ = _c2_run(tmp_path, monkeypatch, names=("e0.5", "e1"))
    with pytest.raises(SystemExit, match="one branch"):
        SC.main(base)
    with pytest.raises(SystemExit):
        SC.main([a for a in base if a != "--probes" and not a.endswith("eval_views")])


def test_only_topology_takes_no_vision_flag(tmp_path, monkeypatch):
    """SelectStage: with only set, neither vision nor no_vision (G-vision is not measured under --only)."""
    _, base, _ = _c2_run(tmp_path, monkeypatch)
    with pytest.raises(SystemExit) as e:
        SC.main(base + ["--no-vision"])
    assert e.value.code == 2


def test_full_select_still_needs_shifted(tmp_path, monkeypatch):
    _, base, _ = _c2_run(tmp_path, monkeypatch)
    with pytest.raises(SystemExit):
        SC.main(base[:-2] + ["--no-vision"])


# ---------------------------------------------------------------- the sampler child's extra flags, prompt-id echo


def test_sampler_reuse_is_keyed_on_flags_sampler_and_checkpoint(tmp_path, monkeypatch):
    assert SC.sampler_cmd("vllm", Path("p"), Path("m"), Path("c"), ["--logprobs"])[-1] == "--logprobs"
    prompts, comps, calls = tmp_path / "p.jsonl", tmp_path / "c.jsonl", []
    prompts.write_text(json.dumps({"id": "a", "n": 1}) + "\n")
    monkeypatch.setattr(SC, "run_group", lambda cmd, timeout: calls.append(cmd) or comps.write_text(
        json.dumps({"id": "a", "k": 0}) + "\n"))
    for kind, model, extra in (("vllm", tmp_path, ()), ("vllm", tmp_path, ["--logprobs"]),  # sampled anew
                               ("vllm", tmp_path, ["--logprobs"]),  # reused
                               ("hf", tmp_path, ()),  # the D13 HF fallback must not reuse vLLM's samples
                               ("hf", tmp_path / "other", ())):  # nor another checkpoint's
        SC.run_sampler(kind, prompts, model, comps, 1.0, extra)
    assert [c[-1] == "--logprobs" for c in calls] == [False, True, False, False]


def test_stop_ids_of_base_weights_without_generation_config(tmp_path):
    """Qwen/Qwen3.5-0.8B ships no generation_config.json; its eos is only in text_config (248044)."""
    (tmp_path / "config.json").write_text(json.dumps({"model_type": "qwen3_5", "eos_token_id": None,
                                                      "text_config": {"model_type": "qwen3_5_text",
                                                                      "eos_token_id": 248044}}))
    proc = SimpleNamespace(tokenizer=SimpleNamespace(eos_token_id=248046, pad_token_id=248044))
    assert SC.stop_ids(tmp_path, proc) == {248044, 248046}
    (tmp_path / "generation_config.json").write_text(json.dumps({"eos_token_id": [7, 248044]}))
    assert SC.stop_ids(tmp_path, proc) == {7, 248044, 248046}


def test_sample_vllm_turns_off_the_flashinfer_sampler_before_importing_vllm():
    """R1 2026-11-05: FlashInfer's sampler JIT-compiles with nvcc (absent from /opt/vllm); vLLM reads env at import."""
    import inspect

    src = inspect.getsource(V.main)
    assert V.ENV == {"VLLM_USE_FLASHINFER_SAMPLER": "0"}
    assert -1 < src.index("os.environ.update(ENV)") < src.index("import vllm")


def test_sample_vllm_echoes_prompt_ids_with_logprobs():
    lp = [{9: SimpleNamespace(logprob=-3.0), 5: SimpleNamespace(logprob=-0.5)}, {7: SimpleNamespace(logprob=-0.1)}]
    c = SimpleNamespace(text="t", token_ids=[5, 7], finish_reason="stop", logprobs=lp)
    res = SimpleNamespace(outputs=[c], prompt_token_ids=[1, 2, 3])
    (rec,) = V.outputs([{"id": "a"}], [res])
    assert rec["logprobs"] == [-0.5, -0.1] and rec["prompt_token_ids"] == [1, 2, 3]
    (rec,) = V.outputs([{"id": "a"}], [SimpleNamespace(outputs=[SimpleNamespace(**vars(c) | {"logprobs": None})],
                                                       prompt_token_ids=[1, 2, 3])])
    assert "prompt_token_ids" not in rec and "logprobs" not in rec  # the sampling bridge is unchanged


# ---------------------------------------------------------------- R0 parity driver


class _Proc:
    """A processor whose prompt ids are the prompt's characters (vLLM echoes them unless told otherwise)."""

    def __call__(self, text, images, add_special_tokens, return_tensors):
        assert add_special_tokens is False and len(images[0]) == 1
        return {"input_ids": np.array([[ord(ch) for ch in text[0]]])}


def _r0(tmp_path, monkeypatch, *, bad_ids=(), shift=0.0):
    from PIL import Image

    from experiments.phase2 import r0_parity as R

    Image.new("RGB", (8, 8)).save(tmp_path / "im.png")
    seen = {"teacher": [], "sample": []}

    def sample_view(a, branch, rows, out, k, extra=()):
        seen["sample"].append({"branch": branch, "rows": [r["key"] for r in rows], "k": k, "extra": tuple(extra),
                               "sampler": a.sampler, "model": a.model})
        out.mkdir(parents=True, exist_ok=True)
        ps = [{"id": r["key"], "prompt": f"<{r['key']}>", "image": str(tmp_path / "im.png"),
               "processor_kwargs": {}, "n": k} for r in rows]
        cs = [{"id": p["id"], "k": j, "token_ids": [9, 9, 9], "logprobs": [-1.0 + shift] * 3,
               "prompt_token_ids": [ord(ch) for ch in p["prompt"]][: -1 if p["id"] in bad_ids else None]}
              for p in ps for j in range(k)]
        (out / "prompts.jsonl").write_text("".join(json.dumps(p) + "\n" for p in ps))
        (out / "completions.jsonl").write_text("".join(json.dumps(c) + "\n" for c in cs))
        return {}

    def teacher(model, processor, x, ids):
        seen["teacher"].append(x["prompt"])
        return np.full(len(ids), -1.0)

    monkeypatch.setattr(R, "weights", lambda model: tmp_path / "base")
    monkeypatch.setattr(SC, "view_rows", lambda view, items=None, rl_only=True: [{"key": f"lidc/n{i:02d}/z0"}
                                                                                 for i in range(items)])
    monkeypatch.setattr(SC, "sample_view", sample_view)
    monkeypatch.setattr(parity, "load", lambda model_dir, prompts, dtype: ("model", _Proc()))
    monkeypatch.setattr(parity, "teacher", teacher)
    out = tmp_path / "parity-lfm"
    argv = ["--model", "lfm2.5-vl-450m", "--view", str(tmp_path / "lidc-val-shifted"), "--out", str(out)]
    return R, argv, out, seen


def test_r0_parity_passes_on_matching_ids_and_logprobs(tmp_path, monkeypatch):
    R, argv, out, seen = _r0(tmp_path, monkeypatch)
    assert R.main(argv) == 0
    (s,) = seen["sample"]
    assert s["branch"] == tmp_path / "base" and len(s["rows"]) == 32 and s["k"] == 8
    assert s["extra"] == ("--logprobs",) and s["sampler"] == "vllm" and s["model"] == "lfm2.5-vl-450m"
    assert seen["teacher"] == [f"<{k}>" for k in s["rows"][:8] for _ in range(8)]  # 8 items x n = 8
    doc = json.loads((out / "parity.json").read_text())
    assert doc["pass"] is True and doc["dtype"] == "bfloat16" and doc["sequences"] == 64
    assert doc["prompt_ids"] == {"items": 32, "equal": 32, "mismatched": [], "pass": True}
    assert doc["thresholds"] == {"median_abs": 0.02} and doc["report_only"] == {"seq_p99_abs": 0.5}  # D13, 2026-11-05
    assert doc["model"] == "lfm2.5-vl-450m" and doc["weights"] == str(tmp_path / "base")


def test_r0_parity_fails_on_one_prompt_id_mismatch(tmp_path, monkeypatch):
    R, argv, out, _ = _r0(tmp_path, monkeypatch, bad_ids={"lidc/n20/z0"})  # outside the logprob items
    assert R.main(argv + ["--items", "24"]) == 1
    doc = json.loads((out / "parity.json").read_text())
    assert doc["prompt_ids"] == {"items": 24, "equal": 23, "mismatched": ["lidc/n20/z0"], "pass": False}
    assert doc["pass"] is False and doc["median_abs"] == 0.0


def test_r0_parity_fails_on_a_logprob_shift(tmp_path, monkeypatch):
    R, argv, out, _ = _r0(tmp_path, monkeypatch, shift=0.05)
    assert R.main(argv) == 1
    doc = json.loads((out / "parity.json").read_text())
    assert doc["prompt_ids"]["pass"] is True and doc["median_abs"] == pytest.approx(0.05) and doc["pass"] is False


def test_r0_parity_help_and_no_prepared(tmp_path, monkeypatch):
    R, argv, _, _ = _r0(tmp_path, monkeypatch)
    for extra, code in ((["--help"], 0), (["--prepared", "p"], 2)):  # --prepared was never used: not accepted
        with pytest.raises(SystemExit) as e:
            R.main(argv + extra)
        assert e.value.code == code


def test_sampler_timeout_kills_the_whole_process_group(tmp_path):
    """vLLM runs its engine in a child process (EngineCore) that holds the GPU; subprocess.run's timeout kills only the
    direct child, so a hung sampler would leave the GPU taken for the next job in the chain."""
    import os
    import time

    pid = tmp_path / "grandchild.pid"
    child = "import time; time.sleep(60)"
    code = (f"import subprocess, sys, time; p = subprocess.Popen([sys.executable, '-c', {child!r}]);"
            f" open({str(pid)!r}, 'w').write(str(p.pid)); time.sleep(60)")
    with pytest.raises(subprocess.TimeoutExpired):
        SC.run_group([sys.executable, "-c", code], 2.0)
    gc = int(pid.read_text())
    for _ in range(50):
        try:
            os.kill(gc, 0)
        except ProcessLookupError:
            break
        time.sleep(0.1)
    else:
        os.kill(gc, 9)
        pytest.fail("the sampler's own child survived the timeout")
    with pytest.raises(subprocess.CalledProcessError):
        SC.run_group([sys.executable, "-c", "raise SystemExit(3)"], 10)
