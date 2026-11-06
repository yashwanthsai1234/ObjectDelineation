"""R1 / C2 evaluation on the D13 HF path (V1 R0, RTX 4090: vLLM failed H-A5-2, so `sampler: hf`): hf_sample batches
items with left padding and resumes; sample.sample batches items; select_ckpt kills its sampler child on SIGTERM and
retries without losing finished items; the R0 HF rate step and project_r1's evaluation allowance from it; k3_harness's
exit code; c2_analyze on select_ckpt --only topology's report.json; rx8's sampler; g_probe.json without NaN."""

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "experiments/phase2/configs"


def _yaml(path):
    return pytest.importorskip("yaml").safe_load(path.read_text())


# ---------------------------------------------------------------- hf_sample: batches of items, resume


@pytest.fixture(params=["lfm2.5-vl-450m", "qwen3-vl-2b"])
def tiny_prompts(tmp_path, request):
    """A random tiny model and 4 prompts on images of 3 sizes (different image-token counts, so left padding; for
    Qwen3-VL also different rotary position grids)."""
    pytest.importorskip("trl")
    from PIL import Image
    from test_p2_scripts import _tiny_vlm
    from transformers import AutoProcessor

    model = _tiny_vlm(tmp_path / "tiny", request.param)
    proc = AutoProcessor.from_pretrained(model)
    msg = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": "Outline the nodule."}]}]
    text = proc.apply_chat_template(msg, add_generation_prompt=True, tokenize=False)
    sampling = {"temperature": 1.0, "top_p": 1.0, "top_k": 1, "min_p": 0.0, "repetition_penalty": 1.0,
                "presence_penalty": 0.0, "frequency_penalty": 0.0, "max_tokens": 5,
                "stop_token_ids": [proc.tokenizer.eos_token_id], "seed": 0}  # top_k 1: the argmax, RNG-free
    rng, lines = np.random.default_rng(0), []
    for i, size in enumerate([(64, 64), (200, 120), (300, 300), (64, 64)]):
        Image.fromarray(rng.integers(0, 255, (*size[::-1], 3), dtype=np.uint8)).save(tmp_path / f"im{i}.png")
        lines.append({"id": f"lidc/n{i}/z0", "prompt": text, "image": str(tmp_path / f"im{i}.png"),
                      "processor_kwargs": {}, "sampling": sampling, "n": 2})
    (tmp_path / "p.jsonl").write_text("".join(json.dumps(x) + "\n" for x in lines))
    return model, tmp_path / "p.jsonl"


def test_hf_sample_batched_items_equal_one_item_per_call(tiny_prompts, tmp_path):
    """Left padding with per-item image tokens: a batch of 3 items (+ 1) gives each item what it gets alone."""
    from experiments.phase2 import hf_sample as HS

    model, prompts = tiny_prompts
    outs = {}
    for b in (1, 3):
        outs[b] = tmp_path / f"c{b}.jsonl"
        assert HS.main(["--prompts", str(prompts), "--model", str(model), "--out", str(outs[b]),
                        "--batch_items", str(b), "--timing", str(tmp_path / f"t{b}.json")]) == 0
    one, three = ([json.loads(x) for x in outs[b].read_text().splitlines()] for b in (1, 3))
    assert one == three and [(c["id"], c["k"]) for c in one] == [(f"lidc/n{i}/z0", k) for i in range(4) for k in (0, 1)]
    assert all(set(c) == {"id", "k", "token_ids", "finish_reason", "text"} for c in one)
    t = json.loads((tmp_path / "t3.json").read_text())
    assert t["items"] == 4 and t["samples"] == 8 and t["batch_items"] == 3 and t["generate_s"] > 0
    assert t["s_per_sample"] == pytest.approx(t["generate_s"] / 8) and t["load_s"] > 0 and t["tokens_mean"] > 0


def test_hf_sample_cuts_each_row_at_its_own_first_stop():
    from experiments.phase2 import hf_sample as HS

    ps = [{"id": "a", "n": 2}, {"id": "b", "n": 2}]
    rows = [[5, 7, 0, 0], [5, 5, 5, 5], [7, 9, 9, 0], [9, 2, 7, 5]]  # pad 0 after a stop
    recs = HS.records(rows, ps, {7, 2}, lambda t: "".join(map(str, t)))
    assert [(r["id"], r["k"], r["token_ids"], r["finish_reason"], r["text"]) for r in recs] == [
        ("a", 0, [5, 7], "stop", "57"), ("a", 1, [5, 5, 5, 5], "length", "5555"),
        ("b", 0, [7], "stop", "7"), ("b", 1, [9, 2], "stop", "92")]


def test_hf_sample_resumes_after_the_last_whole_batch(tmp_path, monkeypatch):
    """A killed or timed-out child is restarted on the same --out: whole batches are kept, a half-written batch and
    a cut last line are dropped and resampled with the same seeds (batch boundaries are fixed by --batch_items)."""
    from experiments.phase2 import hf_sample as HS

    prompts = tmp_path / "p.jsonl"
    prompts.write_text("".join(json.dumps({"id": f"i{i}", "n": 2, "sampling": {"seed": 0}, "processor_kwargs": {},
                                           "prompt": "p", "image": "x"}) + "\n" for i in range(5)))
    calls = []

    def generate(model, proc, dev, batch):
        calls.append([i for i, _ in batch])
        return [{"id": p["id"], "k": k, "token_ids": [1], "finish_reason": "stop", "text": "t"}
                for _, p in batch for k in range(p["n"])]

    monkeypatch.setattr(HS, "load", lambda path: ("model", "proc", "cpu"))
    monkeypatch.setattr(HS, "generate", generate)
    out = tmp_path / "c.jsonl"
    argv = ["--prompts", str(prompts), "--model", "m", "--out", str(out), "--batch_items", "2"]
    assert HS.main(argv) == 0 and calls == [[0, 1], [2, 3], [4]]
    full = out.read_text()
    lines = full.splitlines(keepends=True)
    out.write_text("".join(lines[:4] + lines[4:7]) + lines[7][:9])  # batch 0 whole, batch 1 half + a cut line
    calls.clear()
    assert HS.main(argv) == 0 and calls == [[2, 3], [4]]
    assert out.read_text() == full
    calls.clear()
    assert HS.main(argv) == 0 and calls == []  # complete: nothing sampled


# ---------------------------------------------------------------- sample.sample (g_probe): items per call


def test_sample_generates_several_items_per_call(monkeypatch):
    pytest.importorskip("trl")
    from objdel.train import rewards, sample

    calls = []

    def _generate(prompts):
        calls.append(len(prompts))
        return None, [[1, 7]] * len(prompts), None, [f"c-{p}" for p in prompts]

    model = SimpleNamespace(training=True, eval=lambda: None, train=lambda mode: None)
    trainer = SimpleNamespace(model=model, processing_class=None, _generate=_generate,
                              args=SimpleNamespace(per_device_eval_batch_size=16))
    monkeypatch.setattr(rewards, "stop_ids", lambda m, p: {7})
    monkeypatch.setattr(rewards, "_text", lambda c: c)
    rows = [{"prompt": f"p{i}", "key": f"k{i}"} for i in range(5)]
    out = sample.sample(trainer, rows, 8, items=2)
    assert calls == [16, 16, 8]  # 2 items x K = 8 per call, then the last item
    assert [r["item_id"] for r in out] == [f"k{i}" for i in range(5)]
    assert all(r["samples"] == [f"c-p{i}"] * 8 and r["truncated"] == [False] * 8 for i, r in enumerate(out))


# ---------------------------------------------------------------- select_ckpt: SIGTERM, resume across retries


def test_sigterm_to_select_kills_the_sampler_group(tmp_path):
    """The run cap's `timeout` SIGTERMs select_ckpt; its sampler child (own session) must not keep the GPU."""
    pytest.importorskip("trl")
    pid = tmp_path / "child.pid"
    child = f"import os, time; open({str(pid)!r}, 'w').write(str(os.getpid())); time.sleep(60)"
    code = (f"import sys; from experiments.phase2 import select_ckpt as SC; "
            f"SC.run_group([sys.executable, '-c', {child!r}], 60)")
    p = subprocess.Popen([sys.executable, "-I", "-c", code], cwd=ROOT)
    for _ in range(600):
        if pid.exists() and pid.read_text():
            break
        time.sleep(0.1)
    p.send_signal(signal.SIGTERM)
    assert p.wait(30) != 0
    gc = int(pid.read_text())
    for _ in range(50):
        try:
            os.kill(gc, 0)
        except ProcessLookupError:
            return
        time.sleep(0.1)
    os.kill(gc, signal.SIGKILL)
    pytest.fail("the sampler child survived SIGTERM to select_ckpt")


def test_sampler_retry_keeps_finished_items(tmp_path, monkeypatch):
    """The completions file is kept across attempts (the child resumes); an attempt that finished items does not
    count against TRIES; other prompts or another child command start from an empty file."""
    pytest.importorskip("trl")
    from experiments.phase2 import select_ckpt as SC

    prompts, comps, seen = tmp_path / "p.jsonl", tmp_path / "c.jsonl", []
    prompts.write_text("".join(json.dumps({"id": i, "n": 1}) + "\n" for i in "abc"))

    def run(cmd, timeout):  # one more item per attempt, then a timeout
        have = comps.read_text().splitlines() if comps.exists() else []
        seen.append(len(have))
        comps.write_text("".join(x + "\n" for x in have) + json.dumps({"id": "abc"[len(have)], "k": 0}) + "\n")
        raise subprocess.TimeoutExpired(cmd, timeout)

    monkeypatch.setattr(SC, "run_group", run)
    SC.run_sampler("hf", prompts, tmp_path, comps, 1.0)  # 3 timed-out attempts, each finishing an item: done
    assert seen == [0, 1, 2]
    SC.run_sampler("hf", prompts, tmp_path, comps, 1.0)  # reused, no attempt
    assert seen == [0, 1, 2]
    monkeypatch.setattr(SC, "run_group", lambda cmd, timeout: seen.append(-1) or _raise(cmd, timeout))
    prompts.write_text(prompts.read_text() + json.dumps({"id": "d", "n": 1}) + "\n")  # other prompts
    with pytest.raises(RuntimeError, match="failed 2 times"):
        SC.run_sampler("hf", prompts, tmp_path, comps, 1.0)  # two stalled attempts
    assert seen[3:] == [-1, -1] and not comps.exists()  # the old prompts' completions were not resumed


def _raise(cmd, timeout):
    raise subprocess.TimeoutExpired(cmd, timeout)


def test_sampler_timeout_is_sized_per_sample(monkeypatch, tmp_path):
    pytest.importorskip("trl")
    from experiments.phase2 import select_ckpt as SC

    got = []
    monkeypatch.setattr(SC, "run_sampler", lambda kind, p, m, c, timeout, extra=(): got.append(timeout))
    monkeypatch.setattr(SC, "stop_ids", lambda path, proc: {7})
    import transformers

    from objdel.train import sample

    monkeypatch.setattr(transformers.AutoProcessor, "from_pretrained", lambda path: None)
    monkeypatch.setattr(sample, "export_prompts", lambda *a, **k: None)
    monkeypatch.setattr(sample, "import_completions", lambda *a, **k: [])
    a = SimpleNamespace(model="lfm2.5-vl-450m", seed=0, sampler="hf", timeout_s=SC.TIMEOUT_S)
    SC.sample_view(a, tmp_path, [{"key": str(i)} for i in range(100)], tmp_path / "o", 8)
    assert got == [600 + SC.TIMEOUT_S * 800] and SC.TIMEOUT_S == 1.0


# ---------------------------------------------------------------- R0 HF rate step, project_r1


def test_hf_rate_step_samples_base_weights_and_writes_the_rate(tmp_path, monkeypatch):
    pytest.importorskip("trl")
    from experiments.phase2 import hf_rate as HR
    from experiments.phase2 import r0_parity
    from experiments.phase2 import select_ckpt as SC

    seen = []

    def sample_view(a, branch, rows, out, k, extra=()):
        seen.append({"branch": branch, "rows": len(rows), "k": k, "sampler": a.sampler, "extra": extra})
        Path(extra[1]).write_text(json.dumps({"items": len(rows), "samples": len(rows) * k, "generate_s": 64.0,
                                              "s_per_sample": 0.25, "load_s": 9.0, "batch_items": 16,
                                              "tokens_mean": 90.0, "length_share": 0.0}))
        return {}

    monkeypatch.setattr(r0_parity, "weights", lambda model: tmp_path / "base")
    monkeypatch.setattr(SC, "view_rows", lambda view, items=None, rl_only=True: [{"key": str(i)} for i in range(items)])
    monkeypatch.setattr(SC, "sample_view", sample_view)
    out = tmp_path / "hf-rate-lfm"
    assert HR.main(["--model", "lfm2.5-vl-450m", "--view", str(tmp_path / "v"), "--out", str(out)]) == 0
    (s,) = seen
    assert s["branch"] == tmp_path / "base" and s["rows"] == 32 and s["k"] == 8 and s["sampler"] == "hf"
    assert s["extra"][0] == "--timing"
    doc = json.loads((out / "hf_rate.json").read_text())
    assert doc["model"] == "lfm2.5-vl-450m" and doc["s_per_sample"] == 0.25 and doc["load_s"] == 9.0


def test_r0_measures_the_hf_rate_of_the_rung():
    steps = _yaml(CONFIGS / "r0.yaml")["steps"]
    rate = [s for s in steps if "experiments.phase2.hf_rate" in s]
    assert [s[s.index("--model") + 1] for s in rate] == ["qwen3-vl-2b"]
    assert all(s[0] == "~" for s in rate)  # a failure leaves project_r1 without a rate (exit 2), not R0 stopped
    assert all(s[s.index("--view") + 1] == "{data}/eval_views/lidc-val-shifted" for s in rate)


def _manifest(tmp_path, name, s_per_step, model="lfm2.5-vl-450m", workers=6, det=False, batch=32, **cfg):
    d = tmp_path / name
    d.mkdir()
    cfg = {"model": model, "output_dir": str(d), "deterministic": det, "sft": {"dataloader_num_workers": workers}} | cfg
    t = {"s_per_step": s_per_step, "batch": batch, "rows_per_epoch": 22400}
    (d / "manifest.json").write_text(json.dumps({"config": cfg, "timing": t}))
    return str(d / "manifest.json")


def _rate(tmp_path, model="lfm2.5-vl-450m", s=0.05, load=20.0):
    p = tmp_path / f"rate-{model}.json"
    p.write_text(json.dumps({"model": model, "s_per_sample": s, "load_s": load, "samples": 256}))
    return ["--hf_rate", str(p)]


def test_eval_samples_follow_select_ckpt_in_the_chain():
    """policy.chain's R1 select_ckpt with sampler hf (--items 256): 3 branches x 256 x K 8 for selection, then the
    selected branch on the full views: shifted K 8 and K 16, rescaled, rotated, both probes (1,000 + 680)."""
    pytest.importorskip("trl")
    from experiments.phase2 import project_r1 as P

    assert P.eval_samples(3, 256) == 3 * 256 * 8 + 1084 * 8 + 1084 * 16 + 2 * 1084 * 8 + 1680 * 8 == 62944
    assert P.eval_samples(3, None) == 3 * 1084 * 8 + 1084 * 16 + 2 * 1084 * 8 + 1680 * 8  # selection = shifted K 8
    assert P.sampler_calls(3, 256) == 3 + 6 and P.SUBSET == 256


def test_project_r1_go_and_stop(tmp_path, capsys):
    pytest.importorskip("trl")
    from experiments.phase2 import project_r1 as P

    off = _manifest(tmp_path, "r0-timing-w6-det-off", 2.0)
    on = _manifest(tmp_path, "r0-timing-w6-det-on", 9.0, det=True)
    w0 = _manifest(tmp_path, "r0-timing-w0-det-off", 9.0, workers=0)
    assert P.main(["--model", "lfm2.5-vl-450m", *_rate(tmp_path), on, off, w0]) == 0  # basis: det-off workers 6
    out = capsys.readouterr().out.splitlines()
    doc = json.loads(out[1])
    eval_h = (62944 * 0.05 + 9 * 20.0) / 3600  # r1_lfm.yaml: sampler hf
    assert out[0].startswith("GO") and doc["basis"] == off and doc["r1_steps"] == 1400 + 52 + 105 + 210
    assert doc["eval_h"] == pytest.approx(eval_h) and doc["sampler"] == "hf" and doc["eval_samples"] == 62944
    assert doc["projected_h"] == pytest.approx(1767 * 2.0 / 3600 + eval_h) and doc["stop_h"] == 3.4
    assert doc["arms"][on] == 9.0
    slow = _manifest(tmp_path, "slow", 6.0)  # 1767 x 6 s = 2.95 h + 0.92 h > 3.4
    assert P.main(["--model", "lfm2.5-vl-450m", *_rate(tmp_path), slow]) == 1
    assert capsys.readouterr().out.startswith("STOP")
    slow_hf = _rate(tmp_path, s=0.2)  # 62944 x 0.2 s = 3.5 h of HF sampling alone
    assert P.main(["--model", "lfm2.5-vl-450m", *slow_hf, off]) == 1
    capsys.readouterr()
    qwen = _manifest(tmp_path, "r0-qwen3vl-arm", 2.7, model="qwen3-vl-2b")  # LoRA and checkpointing by model (D6, D4)
    assert P.main(["--model", "qwen3-vl-2b", *_rate(tmp_path, "qwen3-vl-2b", s=0.08), qwen, off]) == 0
    doc = json.loads(capsys.readouterr().out.splitlines()[1])
    assert doc["stop_h"] == 6.0 and doc["projected_h"] < 6.0


def test_project_r1_vllm_uses_the_r0_rate_and_adds_the_hf_check(tmp_path, capsys, monkeypatch):
    """D13 (2026-11-05): a vLLM R1 samples every branch on the full shifted view at the R0 vLLM rate, then select_ckpt's
    HF check re-samples HF_CHECK_ITEMS x K with HF."""
    pytest.importorskip("trl")
    from experiments.phase2 import project_r1 as P
    from experiments.phase2 import select_ckpt as SC

    assert P.HF_CHECK_ITEMS == SC.HF_CHECK_ITEMS
    cfg = _yaml(CONFIGS / "r1_lfm.yaml")
    cfg.pop("sampler")
    (tmp_path / "r1.yaml").write_text(json.dumps(cfg))
    monkeypatch.setitem(P.R1, "lfm2.5-vl-450m", (str(tmp_path / "r1.yaml"), 3.4))
    assert P.main(["--model", "lfm2.5-vl-450m", *_rate(tmp_path), _manifest(tmp_path, "a", 2.0)]) == 0
    doc = json.loads(capsys.readouterr().out.splitlines()[1])
    n = P.eval_samples(3, None)
    want = (n * P.VLLM_S_PER_SAMPLE + P.sampler_calls(3, None) * P.VLLM_LOAD_S + 64 * 8 * 0.05 + 20.0) / 3600
    assert doc["sampler"] == "vllm" and doc["eval_samples"] == n and doc["eval_h"] == pytest.approx(want)
    assert doc["hf_check_s"] == pytest.approx(64 * 8 * 0.05 + 20.0)


@pytest.mark.parametrize("why", ["w0", "b16", "none", "no_rate", "rate_model", "freeze", "grad_ckpt", "missing"])
def test_project_r1_exits_2_on_usage_and_data_errors(tmp_path, why, capsys):
    pytest.importorskip("trl")
    from experiments.phase2 import project_r1 as P

    m = {"w0": lambda: _manifest(tmp_path, "w0", 2.0, workers=0),
         "b16": lambda: _manifest(tmp_path, "b16", 2.0, batch=16),
         "none": lambda: _manifest(tmp_path, "none", None),
         "freeze": lambda: _manifest(tmp_path, "f", 2.0, freeze_vision=True),
         "grad_ckpt": lambda: _manifest(tmp_path, "g", 2.0, sft={"dataloader_num_workers": 6,
                                                                  "gradient_checkpointing": True}),
         "missing": lambda: str(tmp_path / "nope.json")}.get(why, lambda: _manifest(tmp_path, "ok", 2.0))()
    rate = {"no_rate": [], "rate_model": _rate(tmp_path, "qwen3-vl-2b")}.get(why, _rate(tmp_path))
    assert P.main(["--model", "lfm2.5-vl-450m", *rate, m]) == 2
    assert "project_r1:" in capsys.readouterr().err


def test_project_r1_cli(tmp_path):
    pytest.importorskip("trl")
    m = _manifest(tmp_path, "r0-timing-w6-det-off", 2.0)
    cmd = [sys.executable, "-I", "-m", "experiments.phase2.project_r1", "--model", "lfm2.5-vl-450m"]
    r = subprocess.run(cmd + [*_rate(tmp_path), m], cwd=ROOT, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0 and r.stdout.startswith("GO"), r.stderr[-2000:]
    r = subprocess.run(cmd + [m], cwd=ROOT, capture_output=True, text=True, timeout=300)
    assert r.returncode == 2 and "hf_rate" in r.stderr


# ---------------------------------------------------------------- k3_harness exit code


@pytest.mark.parametrize(("rollout_ms", "code"), [(1.0, 0), (10.0, 1)])
def test_k3_harness_exits_1_on_a_failed_gate(tmp_path, monkeypatch, rollout_ms, code):
    from experiments.phase2 import k3_harness as K3

    root = str(tmp_path.resolve())
    monkeypatch.setattr(K3.serve, "preload", lambda r, ds: K3.serve.STORE.setdefault(root, {"lidc/a": None}) and 1)
    monkeypatch.setattr(K3.serve, "fork_pool", lambda n: SimpleNamespace(shutdown=lambda: None))
    monkeypatch.setattr(K3, "batch", lambda riga, lidc, share, seed: ["job"])
    monkeypatch.setattr(K3, "warm_wall", lambda jobs, pool=None: (1.0, "h"))
    monkeypatch.setattr(K3, "per_rollout_ms", lambda jobs: [rollout_ms])  # K4: p99 <= 5 ms
    argv = ["--out", str(tmp_path / "k3.json"), "--root", root, "--synth", "0", "--batches", "1", "--threads", "1",
            "--workers", "1"]
    try:
        assert K3.main(argv) == code
    finally:
        K3.serve.STORE.pop(root, None)
        K3.serve.STORE.pop("synthetic", None)  # main's synthetic RIGA caches
    assert json.loads((tmp_path / "k3.json").read_text())["gate"]["K4_pass"] is (code == 0)


# ---------------------------------------------------------------- c2_analyze on select_ckpt's report.json


def _report(tmp_path, name, level, ids, rng):
    topo = {"bins": [np.clip(level + 0.1 * rng.standard_normal(120), 0, 1).tolist() for _ in range(3)],
            "false_nonsimple": 0.01, "legal": 0.95, "ids": ids}
    p = tmp_path / f"{name}.json"
    p.write_text(json.dumps({"only": "topology", "branch": "e1", "sampler": "hf", "k": 8,
                             "probes": {"topology": topo}}))
    return str(p)


def test_c2_analyze_reads_select_ckpts_topology_report(tmp_path, monkeypatch):
    pytest.importorskip("scipy")
    from experiments.phase2 import c2_analyze as C2

    rng, ids = np.random.default_rng(2), [f"lidc/p{i}" for i in range(680)]
    paths = [_report(tmp_path, f"s{i}", v, ids, rng) for i, v in enumerate((0.7, 0.68, 0.72))]
    out = tmp_path / "c2.json"
    monkeypatch.setattr(sys, "argv", ["c2_analyze", "--seeds", *paths, "--out", str(out)])
    C2.main()
    assert json.loads(out.read_text())["outcome"] == "pass"
    bad = _report(tmp_path, "other", 0.7, ids[::-1], rng)  # other items: per-item means could not be paired
    monkeypatch.setattr(sys, "argv", ["c2_analyze", "--seeds", *paths[:2], bad, "--out", str(out)])
    with pytest.raises(SystemExit, match="items"):
        C2.main()


def test_topology_report_carries_the_probe_item_ids(tmp_path, monkeypatch):
    pytest.importorskip("trl")
    from test_p2_v1_select import _c2_run

    from experiments.phase2 import select_ckpt as SC

    run, base, _ = _c2_run(tmp_path, monkeypatch)
    assert SC.main(base + ["--sampler", "hf"]) == 0
    topo = json.loads((run / "select/report.json").read_text())["probes"]["topology"]
    assert topo["ids"] == [f"probe-topology/{i}" for i in range(3)]


# ---------------------------------------------------------------- configs, g_probe.json


def test_every_config_that_samples_on_v1_names_its_sampler():
    """D13 fallback (b37917b): HF generate, except R1 Qwen3-VL-2B: vLLM with select_ckpt's HF check (owner
    2026-11-05)."""
    paths = sorted(CONFIGS.glob("r1_*.yaml")) + sorted(CONFIGS.glob("c2_*.yaml")) + [CONFIGS / "rx8.yaml"]
    want = {p.stem: "hf" for p in paths} | {"r1_qwen3vl2b": "vllm"}
    assert {p.stem: _yaml(p).get("sampler") for p in paths} == want


def test_g_probe_json_writes_null_for_nan(tmp_path):
    pytest.importorskip("scipy")
    from experiments.phase2 import g_probe as GP
    from objdel.eval import signal as S

    gains = np.zeros(40)  # constant gains: Spearman undefined
    rho = S.spearman_ci(np.arange(40.0), gains, [f"lidc/{i}" for i in range(40)], boot=50)
    GP.save(tmp_path, {"verdict": "FAIL", "spearman_gain": {"sd_geom": rho}, "x": np.float32("nan")})

    def strict(c):
        raise ValueError(c)

    doc = json.loads((tmp_path / "g_probe.json").read_text(), parse_constant=strict)
    assert doc["spearman_gain"]["sd_geom"][0] is None and doc["x"] is None and doc["verdict"] == "FAIL"


def test_hf_sample_loads_bf16_on_the_gpu_whatever_the_checkpoint_dtype(monkeypatch):
    """V1 review: from_pretrained's default dtype is "auto", so R0's rate (bf16 base weights) and R1's branches (saved
    fp32 by full FT) would sample at different precisions and the projection would undercount. The GPU always samples
    in bf16, the dtype of D13's vLLM path and of the parity check; CPU smoke runs stay fp32."""
    torch = pytest.importorskip("torch")
    transformers = pytest.importorskip("transformers")

    from experiments.phase2 import hf_sample

    seen = []

    class Fake:
        def __init__(self, dtype):
            self.dtype = dtype

        def to(self, dev):
            return self

        def eval(self):
            return self

    monkeypatch.setattr(transformers.AutoModelForImageTextToText, "from_pretrained",
                        lambda path, **kw: seen.append(kw.get("dtype")) or Fake(kw.get("dtype")))
    proc = type("P", (), {"tokenizer": type("T", (), {"padding_side": "left"})()})()
    monkeypatch.setattr(transformers.AutoProcessor, "from_pretrained", lambda path, **kw: proc)
    for cuda, want in ((True, torch.bfloat16), (False, torch.float32)):
        monkeypatch.setattr(torch.cuda, "is_available", lambda c=cuda: c)
        hf_sample.load("m")
        assert seen[-1] == want
