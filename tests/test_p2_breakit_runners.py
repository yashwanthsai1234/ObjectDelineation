"""Break-it BK3 (Phase 2 Stage 1): the runner pipeline train_sft -> select_ckpt -> g_probe, its YAML configs and the
configs / vision pieces it wires. Real defects are strict xfails "BK3-<n>" (experiments/phase2/impl/BK3.md);
passing guards stay. Slow tests train a tiny random LFM2-VL on the CPU (test_p2_scripts._tiny_vlm, no download).
"""

import dataclasses
import json
import os
import socket
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from test_p2_scripts import CONFIGS, CPU_SFT, ROOT, stage, stage_of


def _yaml(path):
    yaml = pytest.importorskip("yaml")
    return yaml.safe_load(path.read_text())


def _train_sft():
    pytest.importorskip("trl")
    from experiments.phase2 import train_sft as TS

    return TS


# ---------------------------------------------------------------- configs


@pytest.mark.parametrize("path", CONFIGS, ids=[p.stem for p in CONFIGS])
def test_every_config_builds_its_sft_config(path):
    """train_sft builds SFTConfig as trainer() does (its fixed overrides | the YAML's sft block) before the model."""
    _train_sft()
    from objdel.train import configs

    a = stage_of(path)
    o = {"seed": a.seed, "data_seed": a.seed, "save_strategy": "no", "save_total_limit": 3, "logging_steps": 10}
    for b in (None, 100):
        cfg = configs.sft_config(a.model, "o", 200, b, **(o | a.sft))
        assert cfg.per_device_train_batch_size * cfg.gradient_accumulation_steps == 32


@pytest.mark.parametrize("path", CONFIGS, ids=[p.stem for p in CONFIGS])
def test_configs_carry_the_plan_values(path):
    """PLAN D5 / D6 / D13 / §5.1: LR 5e-5 full FT and 5e-4 LoRA (Qwen3-VL-2B), the model's completion budget, probe
    LRs 5e-6 / 5e-5, WSD."""
    _train_sft()
    from experiments.phase2.prereg import G_PROBE
    from objdel.train import configs

    a = stage_of(path)
    cfg = configs.sft_config(a.model, "o", 100, **({"dataloader_prefetch_factor": None} | a.sft))
    assert cfg.learning_rate == (5e-4 if configs.lora(a.model) else 5e-5)
    assert configs.MAX_COMPLETION[a.model] == configs.sampling(a.model)["max_tokens"] == {
        "lfm2.5-vl-450m": 1024, "qwen3-vl-2b": 1536}[a.model]
    assert G_PROBE["lr_full"] == 5e-6 and G_PROBE["lr_lora"] == 5e-5 and G_PROBE["steps"] == 20
    assert cfg.lr_scheduler_type == "warmup_stable_decay" and not cfg.packing and cfg.max_length is None
    assert set(a.sft) <= {f.name for f in dataclasses.fields(type(cfg))}


def test_unknown_keys_are_rejected():
    from pydantic import ValidationError

    with pytest.raises(ValidationError, match="brnaches"):
        stage("lfm2.5-vl-450m", "v", "o", brnaches=[0.5])


@pytest.mark.parametrize("vision", [{"prepared": "prepared"}, {"prepared": "prepared", "reference": "<R0 value>"}])
def test_vision_config_is_validated_up_front(vision):
    _train_sft()
    with pytest.raises(ValueError, match="reference"):
        stage("lfm2.5-vl-450m", "views/v", "o", vision=vision)


def test_lora_runs_are_merged():
    src = (ROOT / "experiments/phase2/train_sft.py").read_text()
    assert "merge_and_unload" in src  # D6: "merged after SFT"; the vLLM sampler loads a full model dir


def test_train_sft_can_resume_a_crashed_run():
    from objdel.core.spec import SftStage

    assert {"resume", "resume_from_checkpoint"} & set(SftStage.model_fields)


# ---------------------------------------------------------------- schedule math


@pytest.mark.parametrize("spe", [7, 13, 22, 731, 733])
def test_plan_at_odd_epoch_sizes(spe):
    """warmup 3 % of the 2-epoch stable run, decay 15 % of each branch point, every branch point distinct and
    reachable, the deepest decay inside a 3-epoch view and beyond 2 epochs (S1-8)."""
    TS = _train_sft()
    from objdel.train import configs

    p = TS.plan(stage("lfm2.5-vl-450m", "v", "o"), spe)
    at = [b["at"] for b in p["branches"].values()]
    assert len(set(at)) == 3 and p["stable"] == max(at) == 2 * spe
    assert abs(p["warmup"] - 0.03 * p["stable"]) <= 0.5
    for b in p["branches"].values():
        assert abs(b["decay"] - 0.15 * b["at"]) <= 0.5
        cfg = configs.sft_config("lfm2.5-vl-450m", "o", p["stable"], b["at"])
        assert cfg.max_steps == b["at"] + b["decay"] and cfg.warmup_steps == p["warmup"] < b["at"]
    deepest = max(b["at"] + b["decay"] for b in p["branches"].values())
    assert 2 * spe < deepest <= 3 * spe


def test_plan_refuses_colliding_branch_points():
    TS = _train_sft()
    with pytest.raises(ValueError, match="collide"):  # e0.5 and e1 both at step 1
        TS.plan(stage("lfm2.5-vl-450m", "v", "o"), 1)
    with pytest.raises(ValueError, match="decay is 0"):  # e0.5 at step 2, decay round(0.3) = 0
        TS.plan(stage("lfm2.5-vl-450m", "v", "o", branches=[0.5, 1.0]), 4)


def test_timing_runs_are_stable_only():
    """BK3-14: branches [] = one stable epoch, no branch, no decay; the timing / arm configs use it with vision off."""
    TS = _train_sft()
    want = {"steps_per_epoch": 30, "stable": 30, "warmup": 1, "branches": {}}
    assert TS.plan(stage("lfm2.5-vl-450m", "v", "o", branches=[]), 30) == want
    for p in CONFIGS:
        a = stage_of(p)
        timing = p.stem.startswith("r0_")
        assert (a.branches == () and a.vision is False) if timing else (a.branches and a.vision), p.stem


def test_train_sft_refuses_without_g_vision_and_reads_the_r0_reference(tmp_path):
    """BK3-4 / BK3-13: vision {} is refused; every screening config points at the JSON R0 writes, and a reference
    JSON of another model or fixed set is refused."""
    TS = _train_sft()
    (tmp_path / "v").mkdir()
    (tmp_path / "v/manifest.json").write_text(json.dumps({"config": {"epochs": 3}}))
    with pytest.raises(ValueError, match="G-vision"):
        TS.train(stage("lfm2.5-vl-450m", str(tmp_path / "v"), str(tmp_path / "o")))
    for p in CONFIGS:
        v = _yaml(p).get("vision")
        assert v is False or v["reference"] == f"vision_ref/{_yaml(p)['model']}.json", p.stem
    ref = tmp_path / "ref.json"
    ref.write_text(json.dumps({"model": "lfm2.5-vl-450m", "set_sha256": "abc", "reference": 0.47}))
    assert TS.reference(str(ref), "lfm2.5-vl-450m", "abc") == 0.47 and TS.reference(0.3, "x", "y") == 0.3
    for model, sha in (("qwen3-vl-2b", "abc"), ("lfm2.5-vl-450m", "other")):
        with pytest.raises(ValueError, match="not this run"):
            TS.reference(str(ref), model, sha)


def test_train_sft_refuses_a_view_with_too_few_epochs(tmp_path):
    """S1-8: R1's 2-epoch branch decays into epoch 3; a 2-epoch view is refused before anything loads."""
    TS = _train_sft()
    (tmp_path / "v").mkdir()
    (tmp_path / "v/manifest.json").write_text(json.dumps({"config": {"epochs": 2}}))
    a = stage("lfm2.5-vl-450m", str(tmp_path / "v"), str(tmp_path / "o"), vision=False)
    with pytest.raises(ValueError, match="--epochs 3"):
        TS.train(a)
    assert not (tmp_path / "o").exists()


# ---------------------------------------------------------------- param groups / sampling


@pytest.mark.parametrize("name", ["lfm2.5-vl-450m", "qwen3-vl-2b"])
def test_param_groups_cover_every_parameter_once(name):
    """Full-size models on the meta device (cached configs): every trainable parameter in exactly one group, the
    0.1x group = exactly the vision tower + projector parameters."""
    torch = pytest.importorskip("torch")
    pytest.importorskip("trl")
    from transformers import AutoConfig, AutoModelForImageTextToText

    from objdel.train import configs

    os.environ["HF_HUB_OFFLINE"] = "1"
    try:
        cfg = AutoConfig.from_pretrained(configs.MODELS[name])
    except OSError:
        pytest.skip(f"{name} config not cached")
    with torch.device("meta"):
        model = AutoModelForImageTextToText.from_config(cfg)
    lm, vis = configs.param_groups(model, 5e-5)
    ids = [id(p) for g in (lm, vis) for p in g["params"]]
    named = dict(model.named_parameters())
    assert len(ids) == len(set(ids)) == len(named) and vis["lr"] == 5e-6
    want = {n for n in named if any(k in n for k in ("visual", "vision_tower", "multi_modal_projector"))}
    assert {id(named[n]) for n in want} == {id(p) for p in vis["params"]} and want


def test_grpo_and_bridge_get_the_same_sampling_fields():
    """g_probe's GRPOConfig mapping and the exported prompts carry configs.sampling's values (D13)."""
    pytest.importorskip("trl")
    from objdel.train import configs

    for m in ("lfm2.5-vl-450m", "qwen3-vl-2b"):
        s = configs.sampling(m)
        same = {k: s[k] for k in ("temperature", "top_p", "top_k", "min_p", "repetition_penalty")}
        g = configs.grpo_config(m, "o", use_vllm=False, bf16=False, max_completion_length=s["max_tokens"], **same)
        assert (g.temperature, g.top_p, g.top_k, g.min_p, g.repetition_penalty) == (1.0, 1.0, 0, 0.0, 1.0)
        assert g.max_completion_length == configs.MAX_COMPLETION[m]
        assert s["presence_penalty"] == s["frequency_penalty"] == 0


def _hf_sample_kwargs(tmp_path, monkeypatch):
    """objdel.train.sample's HF child on a mocked model / processor -> the kwargs it hands to generate."""
    torch = pytest.importorskip("torch")
    pytest.importorskip("trl")
    import transformers
    from PIL import Image

    from objdel.train import configs
    from objdel.train import sample as HS

    seen = {}

    class Batch(dict):
        def to(self, dev):
            return self

    class Model:
        def to(self, dev):
            return self

        def eval(self):
            return self

        def generate(self, **kw):
            seen.update(kw)
            return torch.zeros((kw["num_return_sequences"], 4), dtype=torch.long)

    class Proc:
        tokenizer = SimpleNamespace(pad_token_id=0, padding_side="left", decode=lambda t, **k: "x")

        def __call__(self, **kw):
            return Batch(input_ids=torch.zeros((1, 2)))

    monkeypatch.setattr(transformers.AutoModelForImageTextToText, "from_pretrained", lambda *a, **k: Model())
    monkeypatch.setattr(transformers.AutoProcessor, "from_pretrained", lambda *a, **k: Proc())
    Image.new("RGB", (8, 8)).save(tmp_path / "i.png")
    s = configs.sampling("lfm2.5-vl-450m") | {"seed": 0, "stop_token_ids": [7]}
    line = {"id": "a", "prompt": "p", "image": str(tmp_path / "i.png"), "processor_kwargs": {}, "sampling": s, "n": 2}
    (tmp_path / "p.jsonl").write_text(json.dumps(line) + "\n")
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    HS.main(["--kind", "hf", "--prompts", str(tmp_path / "p.jsonl"), "--model", "m", "--out",
             str(tmp_path / "c.jsonl")])
    return seen, s


def test_hf_sample_passes_the_sampling_fields(tmp_path, monkeypatch):
    kw, s = _hf_sample_kwargs(tmp_path, monkeypatch)
    assert (kw["temperature"], kw["top_p"], kw["top_k"], kw["repetition_penalty"]) == (1.0, 1.0, 0, 1.0)
    assert kw["max_new_tokens"] == s["max_tokens"] == 1024 and kw["do_sample"] is True


def test_hf_sample_min_p_is_explicit(tmp_path, monkeypatch):
    kw, _ = _hf_sample_kwargs(tmp_path, monkeypatch)
    assert kw["min_p"] == 0.0


# ---------------------------------------------------------------- select_ckpt guards


def _sc():
    pytest.importorskip("trl")
    from experiments.phase2 import select_ckpt as SC

    return SC


def test_exited_treats_permission_error_as_alive(tmp_path, monkeypatch):
    SC = _sc()
    (tmp_path / "manifest.json").write_text(json.dumps({"status": "done", "pid": 12345}))

    def kill(pid, sig):
        raise PermissionError

    monkeypatch.setattr(SC.os, "kill", kill)
    with pytest.raises(RuntimeError, match="must exit first"):
        SC.exited(tmp_path)


def _fake_select(tmp_path, monkeypatch, *, names=("e0.5", "e1"), eligible=True, hf_check=True):
    """select_ckpt.main with the sampler / scorer replaced: records every score_view call; the HF check (D13, run
    with the default --sampler vllm) returns `hf_check` as its pass."""
    SC = _sc()
    from experiments.phase2 import prereg
    from objdel.train import serve

    run = tmp_path / "run"
    for n in names:
        (run / "branches" / n).mkdir(parents=True)
    calls = []

    def score_view(a, branch, view, out, k):
        calls.append({"branch": branch.name, "view": view.name, "items": a.items, "k": k})
        out.mkdir(parents=True, exist_ok=True)  # as sample_view does
        lidc = {"format": 1.0, "legal": 1.0 if eligible else 0.5, "dup_rate": 0.0, "signal_share": 0.9, "group": 8}
        return {"summary": {"datasets": {"lidc": lidc}}, "panel": {"summary": {}}, "iou": {}}

    monkeypatch.setattr(SC, "exited", lambda run: {"config": {"model": "lfm2.5-vl-450m"}})
    monkeypatch.setattr(serve, "preload", lambda *a, **k: 0)
    monkeypatch.setattr(SC, "score_view", score_view)
    monkeypatch.setattr(SC, "vision_check", lambda b, cfg: {"pass": True} if cfg else {})
    monkeypatch.setattr(SC, "hf_check", lambda a, name, out: {"branch": name, "vllm_minus_hf": {},
                                                                    "pass": hf_check})
    monkeypatch.setattr(
        prereg, "select_rb", lambda bs: (next((b["name"] for b in bs[::-1] if b["eligible"]), None), {"eligible": []})
    )
    base = ["--run", str(run), "--shifted", str(tmp_path / "shifted"), "--out", str(tmp_path / "sel")]
    return SC, base, calls


def test_select_with_every_branch_ineligible_gates_the_deepest(tmp_path, monkeypatch):
    SC, base, calls = _fake_select(tmp_path, monkeypatch, names=("e2", "e0.5", "e1"), eligible=False)
    assert SC.main(base + ["--vision", '{"prepared": "p", "reference": 1}']) == 0
    doc = json.loads((tmp_path / "sel/selected.json").read_text())
    assert doc["selected"] is None and doc["gates_branch"] == "e2"
    assert [c["branch"] for c in calls if c["k"] == 8] == ["e0.5", "e1", "e2"]  # then K = 16 on e2


def test_select_with_an_empty_branch_list(tmp_path, monkeypatch):
    SC, base, calls = _fake_select(tmp_path, monkeypatch, names=())
    (tmp_path / "run" / "branches").mkdir(parents=True)
    assert SC.main(base + ["--vision", '{"prepared": "p", "reference": 1}']) == 0
    doc = json.loads((tmp_path / "sel/selected.json").read_text())
    assert doc["selected"] is None and doc["gates_branch"] is None and not calls


def test_select_refuses_to_judge_eligibility_without_g_vision(tmp_path, monkeypatch, capsys):
    SC, base, _ = _fake_select(tmp_path, monkeypatch)
    with pytest.raises(SystemExit) as e:
        SC.main(base)
    assert e.value.code == 2 and "vision" in capsys.readouterr().err


def test_items_subset_is_for_selection_only(tmp_path, monkeypatch):
    SC, base, calls = _fake_select(tmp_path, monkeypatch)
    rescaled = ["--rescaled", str(tmp_path / "rescaled"), "--items", "256", "--sampler", "hf"]
    SC.main(base + rescaled + ["--vision", '{"prepared": "p", "reference": 1}'])
    sel = [c for c in calls if c["branch"] == "e1"]
    assert {c["view"] for c in sel if c["items"] is None} == {"shifted", "rescaled"}


def test_import_completions_refuses_a_partial_group(tmp_path):
    pytest.importorskip("trl")
    from objdel.train import sample

    rows = [{"id": "a", "k": k, "text": "{}", "token_ids": [1, 7], "finish_reason": "stop"} for k in range(4)]
    (tmp_path / "c.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
    with pytest.raises(ValueError):
        sample.import_completions(tmp_path / "c.jsonl", {7})  # prompts asked for n = 8


def test_sampler_child_failure_stops_select(tmp_path, monkeypatch):
    """A crashing sampler child is not swallowed (subprocess.run(check=True))."""
    SC = _sc()
    cmd = SC.sampler_cmd("hf", tmp_path / "p", tmp_path / "m", tmp_path / "c")
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(cmd[:2] + ["-c", "raise SystemExit(1)"], check=True)
    assert cmd[2:6] == ["-m", "objdel.train.sample", "--kind", "hf"]


# ---------------------------------------------------------------- g_probe null


def test_g_probe_zero_variance_null_is_finite():
    pytest.importorskip("scipy")
    from experiments.phase2.prereg import g_probe

    ids = [f"lidc/{i // 2}/{i}" for i in range(40)]
    r0 = np.linspace(0.1, 0.5, 40)
    for r20, want in ((r0, "INCONCLUSIVE"), (r0 - 0.1, "FAIL"), (r0 + 0.1, "PASS")):
        res = g_probe(r20, r0, r0.copy(), ids, 0.95, 0.95, 0.9)
        assert res["h0"] == 0.0 and res["verdict"] == want


# ---------------------------------------------------------------- slow: in-process train_sft on a tiny model


@pytest.fixture
def tiny_run(tmp_path, synthetic_items, monkeypatch):
    pytest.importorskip("trl")
    from test_p2_scripts import _tiny_vlm

    from experiments.phase2 import build_views as BV
    from objdel.data import store

    model = _tiny_vlm(tmp_path / "tiny")
    store.write(synthetic_items(("train",), unanimous=True), tmp_path / "store")
    BV.build(tmp_path / "store", tmp_path / "view", "lidc", "train", epochs=2, seed=3, procs=1)
    TS = _train_sft()

    def args(out, **kw):
        base = dict(model="lfm2.5-vl-450m", model_path=str(model), view=str(tmp_path / "view"),
                    output_dir=str(tmp_path / out), branches=[0.5, 1.0], rolling_minutes=0.0, vision=False, sft=CPU_SFT,
                    steps_per_epoch=8)  # e0.5 at 4, e1 at 8, each with a 1-step decay (decay 0 is refused, BK3-10)
        return stage(**(base | kw))

    from objdel.eval import vision

    monkeypatch.setattr(vision, "fixed_set", lambda prepared: ([], "set-sha"))  # vision_check is mocked where used
    return TS, args, tmp_path


def _weights(d):
    from safetensors.numpy import load_file

    return {k: v for f in sorted(Path(d).glob("*.safetensors")) for k, v in load_file(f).items()}


def _same_weights(a, b) -> bool:
    wa, wb = _weights(a), _weights(b)
    return wa.keys() == wb.keys() and all(np.array_equal(wa[k], wb[k]) for k in wa)


@pytest.mark.slow
def test_two_runs_are_bit_identical_and_rolling_keeps_three(tiny_run, monkeypatch):
    TS, args, tmp = tiny_run
    rolling = []

    def check(model, proc, cfg):
        rolling.append(len(list((tmp / cur[0] / "rolling").glob("checkpoint-*"))))
        return {"action": "ok", "ratio": 1.0, "tv": 0.1, "tv_black": 0.1, "tv_other": 0.1}

    monkeypatch.setattr(TS, "vision_check", check)
    vis = {"prepared": "unused", "reference": 1.0}
    cur = ["a"]
    assert TS.train(args("a", vision=vis)) == 0
    cur[0] = "b"
    assert TS.train(args("b", vision=vis)) == 0
    assert rolling and max(rolling) <= 3  # rotation (save_total_limit 3) while the stable run saves every step
    assert all(_same_weights(tmp / "a/branches" / b, tmp / "b/branches" / b) for b in ("e0.5", "e1"))
    man = json.loads((tmp / "a/manifest.json").read_text())
    assert man["status"] == "done" and not (tmp / "a/rolling").exists() and len(man["vision"]) >= 4
    assert {(v["device"], v["dtype"]) for v in man["vision"]} == {("cpu", "torch.float32")}  # BK3-12


@pytest.mark.slow
def test_deepest_branch_decay_moves_the_weights(tiny_run, monkeypatch):
    """e1 is the deepest branch (= the stable run's end, plan at 8, decay 1): its one decay step must use the WSD
    factor 1.0 as an uninterrupted run does, so branches/e1 != stable/e1. The stable run's scheduler has stepped past
    its end (factor min_lr_ratio = 0) when the checkpoint is saved; the resumed optimizer keeps that lr for one step."""
    TS, args, tmp = tiny_run
    assert TS.train(args("d")) == 0
    man = json.loads((tmp / "d/manifest.json").read_text())
    assert man["plan"]["branches"]["e1"] == {"at": 8, "decay": 1}
    wa, ws = _weights(tmp / "d/branches/e1"), _weights(tmp / "d/stable/e1")
    assert any(not np.array_equal(wa[k], ws[k]) for k in wa)


@pytest.mark.slow
def test_g_vision_abort_stops_before_the_branches(tiny_run, monkeypatch):
    TS, args, tmp = tiny_run
    calls = []

    def check(model, proc, cfg):
        calls.append(model.training)
        bad = len(calls) >= 4  # after the save at e0.5's branch point (step 4)
        return {"action": "abort" if bad else "warn", "ratio": 0.1, "tv": 0.001 if bad else 0.1}

    monkeypatch.setattr(TS, "vision_check", check)
    assert TS.train(args("ab", vision={"prepared": "unused", "reference": 1.0})) == TS.ABORTED
    man = json.loads((tmp / "ab/manifest.json").read_text())
    assert man["status"] == "aborted" and man["abort"]["step"] == 4 and "freeze_vision" in man["abort"]["reason"]
    assert [v["action"] for v in man["vision"]] == ["warn"] * 3 + ["abort"]  # warn alone did not stop training
    assert not (tmp / "ab/branches").exists() and (tmp / "ab/stable/e0.5").is_dir()
    assert not list((tmp / "ab/stable").rglob("optimizer.pt"))


@pytest.mark.slow
def test_resume_after_a_crash_equals_an_uninterrupted_run(tiny_run, monkeypatch):
    """BK3-8: a crash in the stable phase (here: at the step-6 G-vision check), then --resume: the newest rolling
    checkpoint, then the branches; the branches equal an uninterrupted run's bit for bit. A finished run, or a
    second start without --resume, is refused."""
    TS, args, tmp = tiny_run
    ok = {"action": "ok", "ratio": 1.0, "tv": 0.1, "tv_black": 0.1, "tv_other": 0.1}
    calls = []

    def crash(model, proc, cfg):
        calls.append(1)
        if len(calls) == 6:
            raise RuntimeError("host lost")
        return ok

    vis = {"prepared": "unused", "reference": 1.0}
    monkeypatch.setattr(TS, "vision_check", lambda m, p, c: ok)
    assert TS.train(args("ref", vision=vis)) == 0
    monkeypatch.setattr(TS, "vision_check", crash)
    with pytest.raises(RuntimeError, match="host lost"):
        TS.train(args("cr", vision=vis))
    assert json.loads((tmp / "cr/manifest.json").read_text())["status"] == "running"
    with pytest.raises(ValueError, match="--resume"):
        TS.train(args("cr", vision=vis))
    assert TS.train(args("cr", vision=vis, resume=True)) == 0
    man = json.loads((tmp / "cr/manifest.json").read_text())
    assert man["status"] == "done" and len(man["resumed"]) == 1 and man["host"] == socket.gethostname()
    assert all(_same_weights(tmp / "ref/branches" / b, tmp / "cr/branches" / b) for b in ("e0.5", "e1"))
    with pytest.raises(ValueError, match="--resume"):
        TS.train(args("cr", vision=vis, resume=True))


@pytest.mark.slow
def test_lora_branches_are_saved_merged(tiny_run, monkeypatch):
    """BK3-3 / D6: a LoRA model trains LM-only adapters; branches and stable checkpoints are full merged models.
    The tiny LFM stands in for Qwen3-VL-2B (the tests/test_p2_scripts.py CPU smoke runs the real architecture)."""
    from transformers import AutoModelForImageTextToText

    from objdel.train import configs

    monkeypatch.setattr(configs, "_LORA", {"lfm2.5-vl-450m"})
    TS, args, tmp = tiny_run
    assert TS.train(args("lo")) == 0
    for d in ("branches/e0.5", "branches/e1", "stable/e0.5", "stable/e1"):
        files = {f.name for f in (tmp / "lo" / d).iterdir()}
        assert "model.safetensors" in files and not {f for f in files if f.startswith("adapter_")}, (d, files)
        AutoModelForImageTextToText.from_pretrained(tmp / "lo" / d)
    from safetensors.torch import load_file  # the tiny base is bf16 (numpy cannot read it)

    base, merged = load_file(tmp / "tiny/model.safetensors"), load_file(tmp / "lo/branches/e1/model.safetensors")
    assert base.keys() == merged.keys() and any(not base[k].float().equal(merged[k].float()) for k in base)


@pytest.mark.slow
def test_vision_reference_only_writes_the_r0_reference(tiny_run, monkeypatch):
    """BK3-13: --vision-reference-only computes the base projector distance (CPU float32) into vision.reference,
    under run_dir; a later run reads it."""
    from objdel.eval import vision

    TS, args, tmp = tiny_run
    from PIL import Image

    imgs = [Image.new("RGB", (64, 64), (40 * i, 80, 120)) for i in range(4)]
    monkeypatch.setattr(vision, "fixed_set", lambda prepared: ([{"image": im} for im in imgs], "set-sha"))
    monkeypatch.setattr(vision, "projector_distance", lambda m, p, images: 0.25)
    a = args("unused", run_dir=str(tmp / "box"), vision={"prepared": "p", "reference": "vision_ref/lfm.json"})
    assert TS.train(a.model_copy(update={"vision_reference_only": True})) == 0
    doc = json.loads((tmp / "box/vision_ref/lfm.json").read_text())
    assert doc["reference"] == 0.25 and doc["dtype"] == "float32" and doc["model"] == "lfm2.5-vl-450m"
    assert TS.reference(str(tmp / "box/vision_ref/lfm.json"), "lfm2.5-vl-450m", "set-sha") == 0.25
