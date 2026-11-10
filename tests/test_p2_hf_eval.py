"""R1 / C2 evaluation on the D13 HF path (V1 R0, RTX 4090: vLLM failed H-A5-2, so `sampler: hf`): hf_sample batches
items with left padding and resumes; sample.sample batches items; select_ckpt kills its sampler child on SIGTERM and
retries without losing finished items; the R0 HF rate step (its reader project_r1 is archived, E6.6); k3_harness's
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
    from objdel.train import sample as HS

    model, prompts = tiny_prompts
    outs = {}
    for b in (1, 3):
        outs[b] = tmp_path / f"c{b}.jsonl"
        assert HS.main(["--kind", "hf", "--prompts", str(prompts), "--model", str(model), "--out", str(outs[b]),
                        "--batch_items", str(b), "--timing", str(tmp_path / f"t{b}.json")]) is None
    one, three = ([json.loads(x) for x in outs[b].read_text().splitlines()] for b in (1, 3))
    assert one == three and [(c["id"], c["k"]) for c in one] == [(f"lidc/n{i}/z0", k) for i in range(4) for k in (0, 1)]
    assert all(set(c) == {"id", "k", "token_ids", "finish_reason", "text"} for c in one)
    t = json.loads((tmp_path / "t3.json").read_text())
    assert t["items"] == 4 and t["samples"] == 8 and t["batch_items"] == 3 and t["generate_s"] > 0
    assert t["s_per_sample"] == pytest.approx(t["generate_s"] / 8) and t["load_s"] > 0 and t["tokens_mean"] > 0


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
    import transformers

    from objdel.train import sample

    monkeypatch.setattr(sample, "stop_ids", lambda path, proc: {7})
    monkeypatch.setattr(transformers.AutoProcessor, "from_pretrained", lambda path: None)
    monkeypatch.setattr(sample, "export_prompts", lambda *a, **k: None)
    monkeypatch.setattr(sample, "import_completions", lambda *a, **k: [])
    a = SimpleNamespace(model="lfm2.5-vl-450m", seed=0, sampler="hf", timeout_s=None)
    SC.sample_view(a, tmp_path, [{"key": str(i)} for i in range(100)], tmp_path / "o", 8)
    assert got == [600 + 3.0 * 800] and SC.TIMEOUT_S["hf"] == 3.0  # EV-5: 3x HF's ~1 s per answer


# ---------------------------------------------------------------- R0 HF rate step


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
    monkeypatch.setattr(SC, "eval_rows", lambda view, items=None, rl_only=True: [{"key": str(i)} for i in range(items)])
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


# ---------------------------------------------------------------- k3_harness exit code


@pytest.mark.parametrize(("rollout_ms", "code"), [(1.0, 0), (10.0, 1)])
def test_k3_harness_exits_1_on_a_failed_gate(tmp_path, monkeypatch, rollout_ms, code):
    from experiments.phase2 import k3_harness as K3

    root = str(tmp_path.resolve())
    monkeypatch.setattr(K3.serve, "preload", lambda r, ds: K3.serve.STORE.setdefault(root, {"lidc/a": None}) and 1)
    monkeypatch.setattr(K3, "batch", lambda riga, lidc, share, seed: ["job"])
    monkeypatch.setattr(K3, "warm_wall", lambda jobs: (1.0, "h"))
    monkeypatch.setattr(K3, "per_rollout_ms", lambda jobs: [rollout_ms])  # K4: p99 <= 5 ms
    argv = ["--out", str(tmp_path / "k3.json"), "--root", root, "--synth", "0", "--batches", "1"]
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

    from objdel.train import sample as hf_sample

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
        hf_sample.load("m", "cuda" if cuda else "cpu")
        assert seen[-1] == want
