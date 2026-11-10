"""E4: one sampler interface (objdel.train.sample; policy/sampler.py after the E5 move): params_for, build_prompt,
the vLLM request builder, HfSampler on a tiny offline model, the child-process entry and a light import (no
trl.trainer, no vllm)."""
import json
import os
import subprocess
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from objdel.train import configs, rewards
from objdel.train import sample as S

ROOT = Path(__file__).resolve().parents[2]
VIEW = ROOT / "data/eval_views/lidc-val-shifted"


class SamplingParams:
    def __init__(self, **kw):
        self.kw = kw


def _lines(tmp_path, kwargs=({}, {"max_pixels": 9}, {})):
    return [{"id": f"lidc/1.2.1/n00{i}/z000", "prompt": f"p{i}", "image": str(tmp_path / f"im{i}.png"),
             "processor_kwargs": kw, "sampling": {"temperature": 1.0, "seed": 3}, "n": 2}
            for i, kw in enumerate(kwargs)]


def test_stop_ids_of_base_weights_without_generation_config(tmp_path):
    """Qwen/Qwen3.5-0.8B ships no generation_config.json; its eos is only in text_config (248044)."""
    (tmp_path / "config.json").write_text(json.dumps({"model_type": "qwen3_5", "eos_token_id": None,
                                                      "text_config": {"model_type": "qwen3_5_text",
                                                                      "eos_token_id": 248044}}))
    proc = SimpleNamespace(tokenizer=SimpleNamespace(eos_token_id=248046, pad_token_id=248044))
    assert S.stop_ids(tmp_path, proc) == {248044, 248046}
    (tmp_path / "generation_config.json").write_text(json.dumps({"eos_token_id": [7, 248044]}))
    assert S.stop_ids(tmp_path, proc) == {7, 248044, 248046}


def test_params_for_is_the_one_sampling_source_key_for_key():
    s = (7, 2)
    got = S.params_for("qwen3-vl-2b", n=8, seed=0, stop_token_ids=s)
    want = configs.sampling("qwen3-vl-2b", "vllm") | {"n": 8, "seed": 0, "stop_token_ids": [7, 2]}
    assert got == want and list(got) == list(want)
    assert {"min_p", "repetition_penalty", "presence_penalty", "frequency_penalty"} <= set(got)


def test_requests_carry_processor_kwargs_and_a_seed_per_item(tmp_path):
    lines = _lines(tmp_path)
    reqs, params = S.requests(lines, SamplingParams, lambda p: f"image:{p}")
    assert ["mm_processor_kwargs" in r for r in reqs] == [False, True, False]
    assert reqs[1]["mm_processor_kwargs"] == {"max_pixels": 9}
    assert [p.kw["seed"] for p in params] == [S.request_seed(3, x["id"]) for x in lines]
    assert len({p.kw["seed"] for p in params}) == 3 and all(0 <= p.kw["seed"] < 2**31 for p in params)
    assert all(p.kw["n"] == 2 and "logprobs" not in p.kw for p in params)
    assert S.requests(lines, SamplingParams, str, logprobs=True)[1][0].kw["logprobs"] == 0


def test_outputs_echo_prompt_ids_only_with_logprobs(tmp_path):
    lines = _lines(tmp_path)[:1]
    lp = [{5: types.SimpleNamespace(logprob=-0.5)}, {1: types.SimpleNamespace(logprob=-0.1)}]
    res = [types.SimpleNamespace(prompt_token_ids=[9, 9], outputs=[
        types.SimpleNamespace(text="a", token_ids=[5, 1], finish_reason="stop", logprobs=lp),
        types.SimpleNamespace(text="b", token_ids=[5], finish_reason="length", logprobs=None)])]
    a, b = S.outputs(lines, res)
    assert a["logprobs"] == [-0.5, -0.1] and a["prompt_token_ids"] == [9, 9] and a["k"] == 0
    assert set(b) == {"id", "k", "text", "token_ids", "finish_reason"}


def test_refuses_to_start_on_a_busy_gpu():
    with pytest.raises(SystemExit, match="89 %"):
        S.check_memory(lambda: (89, 100))
    S.check_memory(lambda: (90, 100))


def test_samples_follow_the_line_order_and_the_rewards_truncation(tmp_path):
    lines = _lines(tmp_path)[:2]
    recs = [{"id": x["id"], "k": k, "text": f"{x['id']}{k}", "token_ids": ids}
            for x in reversed(lines) for k, ids in enumerate([[4, 1], [4, 4]])]
    got = S.samples(lines, recs, {1})
    assert [g.item_id for g in got] == [x["id"] for x in lines]
    assert got[0].texts == (f"{lines[0]['id']}0", f"{lines[0]['id']}1") and got[0].token_ids == ((4, 1), (4, 4))
    assert got[0].truncated == (rewards.truncated([4, 1], {1}), rewards.truncated([4, 4], {1})) == (False, True)


def _fake_vllm(monkeypatch, seen):
    import builtins

    class LLM:
        def __init__(self, **kw):
            seen["llm"] = kw

    def fake_import(name, *args, **kw):
        if name == "vllm" or name.startswith("vllm."):
            seen["env"] = os.environ.get("VLLM_USE_FLASHINFER_SAMPLER")
            return types.SimpleNamespace(LLM=LLM, SamplingParams=SamplingParams, lora=types.SimpleNamespace(
                request=types.SimpleNamespace(LoRARequest=lambda *a: ("lora", *a))))
        return real(name, *args, **kw)

    real = builtins.__import__
    torch = pytest.importorskip("torch")
    monkeypatch.setattr(torch.cuda, "mem_get_info", lambda: (95, 100))
    monkeypatch.setenv("VLLM_USE_FLASHINFER_SAMPLER", "1")
    monkeypatch.setattr(builtins, "__import__", fake_import)


def test_vllm_sampler_sets_the_sampler_env_before_importing_vllm(monkeypatch):
    seen = {}
    _fake_vllm(monkeypatch, seen)
    s = S.VllmSampler(Path("ckpt"), "qwen3-vl-2b", 0.5, 1024, logprobs=True)
    assert seen["env"] == "0" and s.lora is None
    assert seen["llm"] == {"model": "ckpt", "dtype": "bfloat16", "gpu_memory_utilization": 0.5, "max_model_len": 1024,
                           "enable_prefix_caching": True, "seed": 0,
                           "limit_mm_per_prompt": {"image": {"count": 1, "width": 512, "height": 512}},
                           "logprobs_mode": "processed_logprobs"}


def test_vllm_sampler_serves_an_adapter_unmerged(monkeypatch, tmp_path):
    """Option F (final_eval --adapter): the engine serves LoRA (the trainer's vllm_engine_kwargs) and every request
    names the adapter folder."""
    seen = {}
    _fake_vllm(monkeypatch, seen)
    s = S.VllmSampler(Path("ckpt"), "qwen3-vl-2b", adapter=tmp_path)
    assert seen["llm"]["enable_lora"] is True and seen["llm"]["max_lora_rank"] == 16
    assert s.lora == ("lora", "adapter", 1, str(tmp_path))


def test_import_loads_neither_trl_trainer_nor_vllm():
    code = ("import sys, objdel.train.sample, objdel.train.evaluate\n"
            "print(sorted(m for m in sys.modules if m.startswith(('trl.trainer', 'vllm'))))")
    out = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "[]"


def test_help_runs_without_vllm():
    out = subprocess.run([sys.executable, "-I", "-m", "objdel.train.sample", "--help"], capture_output=True,
                         text=True, check=True)
    assert "--kind" in out.stdout and "--logprobs" in out.stdout and "--batch_items" in out.stdout


@pytest.mark.data
def test_build_prompt_ids_equal_trls_own_on_64_validation_rows():
    """TRL's GRPOTrainer._tokenize_prompts: apply_chat_template(tokenize=True) on prepare_multimodal_messages; the
    sampler sends build_prompt's text and the image to vLLM, which tokenizes text + image like processor(...)."""
    pytest.importorskip("trl")
    if not (VIEW / "manifest.json").is_file():
        pytest.skip(f"{VIEW} absent")
    from PIL import Image
    from transformers import AutoProcessor
    from trl.data_utils import prepare_multimodal_messages

    from objdel.train import dataset

    os.environ["HF_HUB_OFFLINE"] = "1"
    try:
        proc = AutoProcessor.from_pretrained(configs.MODELS["qwen3-vl-2b"],
                                             **configs.processor_kwargs("qwen3-vl-2b", "lidc"))
    except OSError:
        pytest.skip("Qwen3-VL-2B processor not in the local HF cache")
    rows = dataset.prompts(VIEW, "lidc", "val")[:64]
    assert len(rows) == 64
    for r in rows:
        image = Image.open(r["image"]).convert("RGB")
        trl_ids = proc.apply_chat_template([prepare_multimodal_messages(r["prompt"], images=[image])],
                                           add_generation_prompt=True, tokenize=True, return_dict=True)["input_ids"][0]
        text, paths = S.build_prompt(r, proc)
        assert paths == [r["image"]]
        ours = proc(text=[text], images=[[Image.open(paths[0]).convert("RGB")]])["input_ids"][0]
        assert list(ours) == list(trl_ids), r["key"]


@pytest.fixture
def tiny(tmp_path):
    pytest.importorskip("trl")
    from PIL import Image
    from test_p2_scripts import _tiny_vlm

    from objdel.train.dataset import _user

    model = _tiny_vlm(tmp_path / "tiny", "qwen3-vl-2b")
    rng, rows = np.random.default_rng(0), []
    for i, size in enumerate([(64, 64), (200, 120), (64, 64)]):
        Image.fromarray(rng.integers(0, 255, (*size[::-1], 3), dtype=np.uint8)).save(tmp_path / f"im{i}.png")
        rows.append({"prompt": [_user("lidc")], "image": str(tmp_path / f"im{i}.png"),
                     "key": f"lidc/1.2.1/n00{i}/z000"})
    return model, rows


def test_hf_sampler_returns_samples_of_the_right_shape(tiny):
    model, rows = tiny
    s = S.HfSampler(model, "qwen3-vl-2b", batch_items=2)
    stop = sorted(rewards.stop_ids(s.net, s.processor))
    params = S.params_for("qwen3-vl-2b", n=2, seed=0, stop_token_ids=stop) | {"max_tokens": 4}
    got = s.sample(rows, params)
    assert [g.item_id for g in got] == [r["key"] for r in rows]
    for g in got:
        assert len(g.texts) == len(g.token_ids) == len(g.truncated) == 2
        assert all(len(t) <= 4 for t in g.token_ids)
        assert g.truncated == tuple(rewards.truncated(list(t), set(stop)) for t in g.token_ids)


def test_hf_child_resumes_after_the_last_whole_batch(tmp_path, monkeypatch):
    """A killed or timed-out child is restarted on the same --out: whole batches are kept, a half-written batch and
    a cut last line are dropped and resampled with the same seeds (batch boundaries are fixed by --batch_items)."""
    prompts = tmp_path / "p.jsonl"
    prompts.write_text("".join(json.dumps({"id": f"i{i}", "n": 2, "sampling": {"seed": 0}, "processor_kwargs": {},
                                           "prompt": "p", "image": "x"}) + "\n" for i in range(5)))
    calls = []

    def generate(model, proc, dev, batch):
        calls.append([i for i, _ in batch])
        return [{"id": p["id"], "k": k, "token_ids": [1], "finish_reason": "stop", "text": "t"}
                for _, p in batch for k in range(p["n"])]

    monkeypatch.setattr(S, "load", lambda path, dev: ("model", "proc"))
    monkeypatch.setattr(S, "generate", generate)
    out = tmp_path / "c.jsonl"
    argv = ["--kind", "hf", "--prompts", str(prompts), "--model", "m", "--out", str(out), "--batch_items", "2"]
    S.main(argv)
    assert calls == [[0, 1], [2, 3], [4]]
    full = out.read_text()
    lines = full.splitlines(keepends=True)
    out.write_text("".join(lines[:4] + lines[4:7]) + lines[7][:9])  # batch 0 whole, batch 1 half + a cut line
    calls.clear()
    S.main(argv)
    assert calls == [[2, 3], [4]] and out.read_text() == full


def test_hf_records_cut_each_row_at_its_own_first_stop():
    ps = [{"id": "a", "n": 2}, {"id": "b", "n": 2}]
    rows = [[5, 7, 0, 0], [5, 5, 5, 5], [7, 9, 9, 0], [9, 2, 7, 5]]  # pad 0 after a stop
    recs = S.records(rows, ps, {7, 2}, lambda t: "".join(map(str, t)))
    assert [(r["id"], r["k"], r["token_ids"], r["finish_reason"], r["text"]) for r in recs] == [
        ("a", 0, [5, 7], "stop", "57"), ("a", 1, [5, 5, 5, 5], "length", "5555"),
        ("b", 0, [7], "stop", "7"), ("b", 1, [9, 2], "stop", "92")]
