"""Phase 2 runners (M8): train_sft, select_ckpt, g_probe and the run configs.

`--help` for each script; the WSD step plan; the G-vision abort / warn rule (mocked check results); select_ckpt's
"train_sft has exited" and "CUDA untouched" guards; R-B with the eligibility floors on the B2 fixture panels; the
WSD-branch mechanism (stable checkpoint + resumed decay = one uninterrupted WSD run, bit for bit); sample_vllm's
free-memory refusal (mocked); and a slow CPU smoke on a tiny random LFM2-VL (cached processor, no download):
synthetic store -> build_views -> train_sft (2 stable steps, 2 branches) -> select_ckpt with the HF sampler.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from experiments.phase2 import train_sft as TS

ROOT = Path(__file__).resolve().parents[1]
# train_sft configs: not a step list for run.sh (r0.yaml) nor an RL launch config (b1_launch.yaml, test_b1_plans.py)
CONFIGS = sorted(c for c in (ROOT / "experiments/phase2/configs").glob("*.yaml")
                 if not any(k in c.read_text() for k in ("\nsteps:", "\nrows_plan:")))
CPU_SFT = {"per_device_train_batch_size": 1, "gradient_accumulation_steps": 1, "use_cpu": True, "bf16": False,
           "logging_steps": 1, "dataloader_num_workers": 0, "dataloader_persistent_workers": False,
           "dataloader_prefetch_factor": None}  # tiny CPU runs of train_sft


def _module(*args, env=None, timeout=900):
    """python -I -m <args> from the repo root, output captured."""
    return subprocess.run([sys.executable, "-I", "-m", *args], cwd=ROOT, env=env, capture_output=True, text=True,
                          timeout=timeout)


@pytest.mark.parametrize("script", ["train_sft", "select_ckpt", "g_probe", "train_rl", "hf_sample"])
def test_help(script):
    pytest.importorskip("trl")
    r = _module(f"experiments.phase2.{script}", "--help", timeout=300)
    assert r.returncode == 0 and "usage" in r.stdout, r.stderr[-2000:]


@pytest.mark.parametrize("path", CONFIGS, ids=[p.stem for p in CONFIGS])
def test_configs_parse_as_run_args(path):
    yaml = pytest.importorskip("yaml")  # PyYAML comes with the train extra (transformers)
    cfg = yaml.safe_load(path.read_text())
    args = TS.RunArgs(**cfg)
    assert args.model in ("lfm2.5-vl-450m", "qwen3-vl-2b") and args.view.startswith("views/")
    assert args.data_root == "data" and all(b > 0 for b in args.branches)
    assert args.sft.get("gradient_accumulation_steps", 4) * args.sft.get("per_device_train_batch_size", 8) == 32  # D5


def test_config_set_complete():
    yaml = pytest.importorskip("yaml")
    names = {p.stem for p in CONFIGS}
    assert names == {"r1_lfm", "r1_qwen3vl2b", "r0_qwen3vl_arm", "rx8", "c2_a0p_seed2", "c2_a0p_seed3", "c2_a0_seed1"}
    assert yaml.safe_load((ROOT / "experiments/phase2/configs/c2_a0_seed1.yaml").read_text())["a0"] is True


def test_plan_branches_decay_warmup():
    pytest.importorskip("trl")
    from objdel.train import configs

    p = TS.plan(TS.RunArgs("lfm2.5-vl-450m", "v", "o"), 1000)
    assert p["branches"] == {"e0.5": {"at": 500, "decay": 75}, "e1": {"at": 1000, "decay": 150},
                             "e2": {"at": 2000, "decay": 300}}
    assert p["stable"] == 2000 and p["warmup"] == 60  # 3 % of the stable steps
    for b in p["branches"].values():  # the plan is the schedule configs.sft_config builds
        cfg = configs.sft_config("lfm2.5-vl-450m", "o", p["stable"], b["at"])
        assert cfg.max_steps == b["at"] + b["decay"] and cfg.warmup_steps == p["warmup"]


@pytest.mark.parametrize(
    "ratio, tv, action, stop",
    [(0.9, 0.01, "ok", False), (0.3, 0.01, "warn", False), (0.9, 0.001, "warn", False), (0.3, 0.001, "abort", True)],
)
def test_vision_abort_only_when_both_fail(ratio, tv, action, stop):
    """train_sft's reaction to mocked G-vision results (E's decide rule): stop and record a reason only on abort."""
    from objdel.eval import vision

    res = {"ratio": ratio, "tv": tv, "tv_black": tv, "tv_other": 0.5} | vision.decide(ratio, tv, 0.5)
    man = {"vision": []}
    assert res["action"] == action and TS.record_vision(man, 7, res) is stop
    assert man["vision"] == [
        {"step": 7, "action": action, "ratio": ratio, "tv": tv, "tv_black": tv, "tv_other": 0.5, "set_sha256": None}
        | {"device": None, "dtype": None}
    ]
    assert ("abort" in man) is stop and (not stop or "ratio" in man["abort"]["reason"])


def test_select_refuses_until_train_sft_exited(tmp_path):
    from experiments.phase2 import select_ckpt as SC

    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    try:
        for status, pid in (("running", 999999), ("done", child.pid), ("aborted", 999999)):
            (tmp_path / "manifest.json").write_text(json.dumps({"status": status, "pid": pid}))
            with pytest.raises(RuntimeError, match="must exit first"):
                SC.exited(tmp_path)
    finally:
        child.kill()
        child.wait()
    (tmp_path / "manifest.json").write_text(json.dumps({"status": "done", "pid": child.pid}))
    assert SC.exited(tmp_path)["status"] == "done"


def test_select_refuses_sampler_after_cuda_init(monkeypatch):
    torch = pytest.importorskip("torch")
    from experiments.phase2 import select_ckpt as SC

    SC.cuda_untouched()
    monkeypatch.setattr(torch.cuda, "is_initialized", lambda: True)
    with pytest.raises(RuntimeError, match="CUDA is initialized"):
        SC.cuda_untouched()
    cmd = SC.sampler_cmd("vllm", Path("p"), Path("m"), Path("c"))
    assert cmd[:2] == [sys.executable, "-I"] and cmd[2].endswith("sample_vllm.py")  # vllm_native: one venv


def test_sampler_uses_the_second_venv_when_the_image_has_one(monkeypatch, tmp_path):
    from experiments.phase2 import select_ckpt as SC

    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "python").touch()
    monkeypatch.setattr(SC, "VLLM_VENV", tmp_path)
    assert SC.sampler_cmd("vllm", Path("p"), Path("m"), Path("c"))[0] == str(tmp_path / "bin" / "python")


def _summary(fmt=1.0, legal=0.95, dup=0.0, sigma=0.9):
    return {"format": fmt, "legal": legal, "dup_rate": dup, "signal_share": sigma, "group": 8}


OK = {"pass": True}  # a measured, passing G-vision


def _eligible(lidc: dict, vision: dict) -> bool:
    """select_ckpt's D12 eligibility call: gates.eligible on the shifted view's summarize() output."""
    from objdel.eval import gates

    return gates.eligible({"datasets": {"lidc": lidc}}, vision)


def test_eligibility_floors():
    assert _eligible(_summary(), OK)
    for bad in (_summary(fmt=0.97), _summary(legal=0.89), _summary(dup=0.11), _summary(sigma=0.49)):
        assert not _eligible(bad, OK)
    assert not _eligible(_summary(), {"pass": False})


def test_rb_on_fixture_branch_panels():
    """R-B over the B2 fixture panels (LFM shifted val, 36 vs 150 steps) with the eligibility floors applied."""
    from experiments.phase2 import prereg
    from objdel.eval import signal

    fix = json.loads((ROOT / "tests/p2_signal_fixture.json").read_text())["runs"]
    shallow, deep = (signal.panel(signal.rows(fix[r])) for r in ("lfm2.5-vl-450m_s0", "lfm2.5-vl-450m_s0_n150"))
    bs = [
        {"name": "e0.5", "depth": 0.5, "panel": shallow, "eligible": _eligible(_summary(), OK)},
        {"name": "e2", "depth": 2.0, "panel": deep, "eligible": _eligible(_summary(), OK)},
    ]
    name, why = prereg.select_rb(bs)
    assert why["eligible"] == ["e2", "e0.5"] and len(why["pairs"]) == 1
    fires = why["pairs"][0]["sd_geom"]["fires"]
    assert name == ("e0.5" if fires else "e2")  # deepest unless the veto fires
    bs[1]["eligible"] = _eligible(_summary(legal=0.80), OK)  # the deep branch fails G-format
    assert prereg.select_rb(bs)[0] == "e0.5"


def test_wsd_resumed_decay_equals_uninterrupted_run(char_model, tmp_path):
    """The train_sft mechanism: stable run with a full checkpoint at b, then a new trainer resumed from it with
    the decay kwargs = one WSD run decaying at b (weights bit-identical, same LR trace)."""
    torch = pytest.importorskip("torch")
    trl = pytest.importorskip("trl")
    from datasets import Dataset
    from transformers import LlamaForCausalLM, TrainerCallback

    path, tok = char_model
    ds = Dataset.from_list([{"text": "{" + "0123456789"[: 3 + i % 6] + "}"} for i in range(48)])
    warm, b, d = 2, 6, 3

    class SaveAt(TrainerCallback):
        def on_step_end(self, a, state, control, **kw):
            control.should_save = state.global_step == b

    def run(out, steps, stable, decay, resume=None):
        cfg = trl.SFTConfig(
            output_dir=str(out), max_steps=steps, per_device_train_batch_size=4, learning_rate=1e-2,
            lr_scheduler_type="warmup_stable_decay", warmup_steps=warm,
            lr_scheduler_kwargs={"num_stable_steps": stable, "num_decay_steps": decay}, use_cpu=True, bf16=False,
            report_to="none", seed=0, data_seed=1, logging_steps=1, save_strategy="no", max_length=64)
        t = trl.SFTTrainer(model=LlamaForCausalLM.from_pretrained(path), args=cfg, train_dataset=ds,
                           processing_class=tok, callbacks=[SaveAt()])
        t.train(resume_from_checkpoint=resume)
        return t.model.state_dict(), [h["learning_rate"] for h in t.state.log_history if "learning_rate" in h]

    ref, lr_ref = run(tmp_path / "ref", b + d, b - warm, d)
    run(tmp_path / "stable", 12, 12 - warm, 0)
    got, lr_got = run(tmp_path / "decay", b + d, b - warm, d, resume=str(tmp_path / f"stable/checkpoint-{b}"))
    assert lr_got == lr_ref and lr_ref[-1] < lr_ref[b - 1]
    assert all(torch.equal(ref[k], got[k]) for k in ref)


# ---------------------------------------------------------------- FX-select: sampler child, subset, g_probe


def test_sampler_child_is_bounded_checked_and_reused(tmp_path, monkeypatch):
    """BK3-15: a partial output is a failed attempt (retried), a hung child times out (bounded), and completions
    matching the prompts' SHA-256 are reused on a rerun."""
    from experiments.phase2 import select_ckpt as SC

    prompts, comps, calls = tmp_path / "p.jsonl", tmp_path / "c.jsonl", []
    prompts.write_text("".join(json.dumps({"id": i, "n": 2}) + "\n" for i in "ab"))

    def run(cmd, timeout):
        calls.append(timeout)
        rows = [{"id": i, "k": k} for i in "ab" for k in range(2)][: 3 if len(calls) == 1 else 4]
        comps.write_text("".join(json.dumps(r) + "\n" for r in rows))

    monkeypatch.setattr(SC, "run_group", run)
    SC.run_sampler("hf", prompts, tmp_path, comps, 7.0)
    assert calls == [7.0, 7.0]  # the partial first attempt was retried
    SC.run_sampler("hf", prompts, tmp_path, comps, 7.0)
    assert len(calls) == 2  # reused
    prompts.write_text(prompts.read_text().replace('"n": 2', '"n": 1'))
    with pytest.raises(RuntimeError, match="failed 2 times"):
        SC.run_sampler("hf", prompts, tmp_path, comps, 7.0)  # other prompts: resampled, 4 lines != n = 1
    monkeypatch.undo()
    monkeypatch.setattr(SC, "sampler_cmd", lambda *a: [sys.executable, "-c", "import time; time.sleep(30)"])
    with pytest.raises(RuntimeError, match="timed out"):
        SC.run_sampler("hf", prompts, tmp_path, tmp_path / "c2.jsonl", 0.5)


def test_items_subset_is_seeded_not_the_first_n(monkeypatch):
    """S1-12 / BK2-2: --items draws a fixed seeded subset of the scored (rl) items; rl_only=False keeps every record."""
    pytest.importorskip("datasets")
    from experiments.phase2 import select_ckpt as SC
    from objdel.train import dataset as D

    recs = [{"item_id": f"lidc/n{i:03d}/z0", "rl": i % 5 != 0, "disagreement": False, "W": 50, "H": 50}
            for i in range(100)]
    monkeypatch.setattr(D, "records", lambda root, ds, split: recs)
    monkeypatch.setattr(D, "_verified", lambda root, rs: rs)
    rl = [r["item_id"] for r in recs if r["rl"]]
    sub = [r["key"] for r in SC.view_rows(Path("v"), 20)]
    assert len(set(sub)) == 20 and set(sub) <= set(rl) and sub != rl[:20] and sub == sorted(sub)
    assert sub == [r["key"] for r in SC.view_rows(Path("v"), 20)]
    assert len(SC.view_rows(Path("v"))) == 80 and len(SC.view_rows(Path("v"), rl_only=False)) == 100


def test_g_probe_report_form():
    """BK3-16: the training-log form is the mean logged raw composite R (steps 1-5 vs 11-20), and the step-0 legal
    rate is the mean of both step-0 draws."""
    pytest.importorskip("scipy")
    from experiments.phase2 import g_probe as GP

    terms = ("iou_tc", "f5", "legal", "g_area")
    log = [{f"rewards/{t}/mean": 1.0 for t in terms} | {"step": i} for i in range(1, 21)] + [{"train_runtime": 1}]
    for h in log[10:20]:
        h["rewards/iou_tc/mean"] = 0.0
    assert GP.train_log(log) == {"steps_1_5": pytest.approx(1.0), "steps_11_20": pytest.approx(0.7)}
    draw = lambda legal: {"per_item": {"r_mean": np.zeros(40)}, "ids": [f"lidc/{i}/z" for i in range(40)],  # noqa: E731
                          "legal": legal, "sigma_share": 0.9}
    res = GP.verdict(draw(0.96), draw(0.90), draw(0.89))
    assert res["verdict"] == "INCONCLUSIVE"  # legal drop 0.93 - 0.89 = 4 pts; draw a alone would say 7 (FAIL)


def test_g_probe_repair_is_the_next_earlier_eligible_branch():
    """D12(c): a non-PASS selected branch -> the next-earlier eligible branch, probed once."""
    pytest.importorskip("scipy")
    from experiments.phase2 import g_probe as GP

    b = lambda d, e: {"depth": d, "eligible": e}  # noqa: E731
    sel = {"selected": "e2", "branches": {"e0.5": b(0.5, True), "e1": b(1.0, False), "e2": b(2.0, True)}}
    assert GP.candidates(sel) == ["e2", "e0.5"]
    assert GP.candidates(sel | {"selected": "e0.5"}) == ["e0.5"]  # no selection: g_probe skips (test_p2_v1_select)


def test_sample_vllm_refuses_below_90_percent_free():
    """sample_vllm's start check on a mocked torch.cuda.mem_get_info (D13: free >= 90 % before LLM(...))."""
    import runpy

    path = ROOT / "experiments/phase2/sample_vllm.py"
    check = runpy.run_path(str(path), run_name="sample_vllm")["check_memory"]
    with pytest.raises(SystemExit, match="free"):
        check(lambda: (int(0.89 * 24e9), int(24e9)))
    check(lambda: (int(0.95 * 24e9), int(24e9)))  # passes


# ---------------------------------------------------------------- CPU smoke (slow)


TINY = {  # model key -> (cached hub id, config fields of a random 2-layer model)
    "lfm2.5-vl-450m": ("LiquidAI/LFM2.5-VL-450M", dict(
        text=dict(hidden_size=32, intermediate_size=64, block_dim=32, block_ff_dim=64, conv_dim=32, conv_dim_out=32,
                  num_hidden_layers=2, layer_types=["conv", "full_attention"], num_attention_heads=2, num_heads=2,
                  num_key_value_heads=1),
        vision=dict(hidden_size=32, intermediate_size=64, num_hidden_layers=1, num_attention_heads=2))),
    "qwen3-vl-2b": ("Qwen/Qwen3-VL-2B-Instruct", dict(
        text=dict(hidden_size=64, intermediate_size=128, num_hidden_layers=2, num_attention_heads=4,
                  num_key_value_heads=2, head_dim=16,
                  rope_parameters={"rope_type": "default", "rope_theta": 5e6, "mrope_interleaved": True,
                                   "mrope_section": [2, 3, 3]}),
        vision=dict(hidden_size=32, intermediate_size=64, depth=2, num_heads=2, out_hidden_size=64,
                    deepstack_visual_indexes=[0, 1]))),
}


def _tiny_vlm(out: Path, model: str = "lfm2.5-vl-450m") -> Path:
    """A random 2-layer model of `model` with its cached processor (no download)."""
    torch = pytest.importorskip("torch")
    from transformers import AutoConfig, AutoModelForImageTextToText, AutoProcessor

    hub, dims = TINY[model]
    os.environ["HF_HUB_OFFLINE"] = "1"
    try:
        cfg = AutoConfig.from_pretrained(hub)
        proc = AutoProcessor.from_pretrained(hub)
    except OSError:
        pytest.skip(f"{hub} config / processor not in the local HF cache")
    cfg.text_config.update(dims["text"])
    cfg.vision_config.update(dims["vision"])
    if model == "lfm2.5-vl-450m":
        cfg.projector_hidden_size = 32
    torch.manual_seed(0)
    AutoModelForImageTextToText.from_config(cfg).save_pretrained(out)
    proc.save_pretrained(out)
    return out


@pytest.mark.slow
@pytest.mark.parametrize("name, select", [("lfm2.5-vl-450m", True), ("qwen3-vl-2b", False)],
                         ids=["lfm-full-ft", "qwen3-vl-lora"])  # full FT; LoRA with checkpointing (D6)
def test_cpu_smoke_train_then_select(name, select, tmp_path, synthetic_items):
    """train_sft (G-vision at the branch points, WSD branches) then, for the LFM only, select_ckpt and g_probe: a random
    Qwen3-VL never stops, so its 1,536-token samples would take hours on a CPU."""
    pytest.importorskip("trl")
    import yaml

    from experiments.phase2 import build_views as BV
    from objdel.data import store
    from objdel.train import dataset

    model = _tiny_vlm(tmp_path / "tiny", name)
    store.write(synthetic_items(("train", "val"), unanimous=True), tmp_path / "store")
    BV.build(tmp_path / "store", tmp_path / "view", "lidc", "train", epochs=2, seed=3, procs=1)
    dataset.prepare(tmp_path / "store", tmp_path / "evalview", procs=1)  # eval-view layout
    run = tmp_path / "run"
    cfg = {"model": name, "model_path": str(model), "view": str(tmp_path / "view"), "output_dir": str(run),
           "branches": [0.5, 1.0], "steps_per_epoch": 8, "run_dir": str(tmp_path), "sft": CPU_SFT,
           "rolling_minutes": 1000.0,  # saves (and real G-vision checks) at the branch points only
           "vision": {"prepared": str(tmp_path / "evalview"), "reference": f"vision_ref/{name}.json"}}
    (tmp_path / "run.yaml").write_text(yaml.safe_dump(cfg))
    env = os.environ | {"HF_HUB_OFFLINE": "1", "CUDA_VISIBLE_DEVICES": ""}
    ref = tmp_path / f"vision_ref/{name}.json"
    for extra in (["--vision_reference_only", "true"], []):  # R0's G-vision base reference, then the run
        tr = _module("experiments.phase2.train_sft", "--config", str(tmp_path / "run.yaml"), *extra, env=env)
        assert tr.returncode == 0, tr.stderr[-4000:]
    assert json.loads(ref.read_text())["reference"] > 0
    man = json.loads((run / "manifest.json").read_text())
    assert man["status"] == "done" and man["plan"]["branches"]["e1"] == {"at": 8, "decay": 1}
    assert [v["step"] for v in man["vision"]] == [4, 8] and all(v["ratio"] > 0 for v in man["vision"])
    assert {p.name for p in (run / "branches").iterdir()} == {"e0.5", "e1"}
    assert not list((run / "stable").rglob("optimizer.pt"))  # stable checkpoints kept model-only
    for d in ("branches/e0.5", "branches/e1", "stable/e1"):  # LoRA is merged: a plain model dir, no adapter files
        files = {f.name for f in (run / d).iterdir()}
        assert "model.safetensors" in files and not {f for f in files if f.startswith("adapter_")}, (d, files)
    ev = str(tmp_path / "evalview")
    if select:
        vision = json.dumps({"prepared": ev, "reference": str(ref)})  # the run's block
        sel = _module("experiments.phase2.select_ckpt", "--run", str(run), "--shifted", ev, "--sampler", "hf",
                      "--procs", "1", "--prepared", ev, "--vision", vision, env=env)
        assert sel.returncode == 0, sel.stderr[-4000:]
        doc = json.loads((run / "select/selected.json").read_text())
        assert set(doc["branches"]) == {"e0.5", "e1"} and doc["sampler"] == "hf" and doc["cuda_initialized"] is False
        assert doc["checker"] == "cpu"  # the exam scored through objdel._gpucheck (its CPU build on this machine)
        for b in ("e0.5", "e1"):
            lines = (run / "select" / b / "shifted/completions.jsonl").read_text().splitlines()
            assert lines and all(np.isfinite(json.loads(x)["k"]) for x in lines)
        gates = json.loads((run / "select/gates.json").read_text())
        assert "ci_lo/shifted/lidc/nodule" in gates["gates"]["G-prior-b"]["values"]
        assert gates["gates"]["G-vision"]["values"]["ratio"] is not None  # measured on every branch
        k16 = json.loads((run / "select/report.json").read_text())["k16"]["datasets"]["lidc"]
        assert k16["samples"] == 16 * k16["items"] and "hm_iou" in k16["objects"]["nodule"]
    best = run / "branches" / (doc["selected"] or "e1") if select else run / "branches/e1"
    gp = _module("experiments.phase2.g_probe", "--model", name, "--ckpt", str(best), "--train", ev,
                 "--probe", ev, "--out", str(run / "gp"), "--procs", "1", "--grpo", json.dumps(GRPO_SMOKE), env=env)
    assert gp.returncode == 0, gp.stderr[-4000:]
    res = json.loads((run / "gp/g_probe.json").read_text())
    assert res["verdict"] in ("PASS", "FAIL", "INCONCLUSIVE") and res["steps"] == 20
    from objdel.train import configs

    assert res["lr"] == configs.rl_lr(name)


# use_vllm off: vLLM needs a GPU (TRL 1.17 refuses use_vllm without it installed); the box's profile arms run vLLM
GRPO_SMOKE = {"per_device_train_batch_size": 8, "gradient_accumulation_steps": 1, "per_device_eval_batch_size": 8,
              "max_completion_length": 8, "use_cpu": True, "bf16": False, "beta": 0.0, "use_vllm": False}


