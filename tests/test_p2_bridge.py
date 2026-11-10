"""Phase 2 M5: the JSONL sampling bridge (C-BRIDGE) between the training env (objdel.policy.sampler) and the vLLM env
(objdel.policy.sampler's vLLM request builder and output records)."""
import json
import types
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


class Processor:
    """apply_chat_template as a processor renders text (tokenize=False); records what it was given."""

    def __init__(self):
        self.seen = []

    def apply_chat_template(self, conversation, add_generation_prompt, tokenize):
        assert add_generation_prompt and not tokenize
        self.seen.append(conversation)
        parts = conversation[0][0]["content"]  # a batch of one conversation, as B5 called it
        return [f"<user>{''.join('<image>' if p['type'] == 'image' else p['text'] for p in parts)}<assistant>"]


@pytest.fixture
def rows(tmp_path):
    pytest.importorskip("datasets")  # BK4-1: the CI [dev] job has no datasets
    from objdel.policy.prompt import _user
    return [{"prompt": [_user("lidc")], "image": str(tmp_path / f"im{i}.png"), "key": f"lidc/1.2.1/n000/z00{i}"}
            for i in (1, 3)]


@pytest.fixture
def V():
    """objdel.policy.sampler, which imports TRL's data utils: the CI [dev]-only job skips the tests that use it."""
    return pytest.importorskip("objdel.policy.sampler")


def export(rows, tmp_path, model="qwen3-vl-2b", **kw):
    pytest.importorskip("trl")
    from objdel.policy import sampler as sample
    path = tmp_path / "prompts.jsonl"
    sample.export_prompts(rows, path, model, "lidc", 8, processor=Processor(), **kw)
    return path, [json.loads(line) for line in path.read_text().splitlines()]


def test_export_schema_kwargs_and_sampling(rows, tmp_path):
    from objdel.policy import models as configs
    from objdel.policy.prompt import PROMPTS
    _, lines = export(rows, tmp_path, stop_ids={2, 1}, seed=3)
    assert [set(x) for x in lines] == [{"id", "prompt", "image", "processor_kwargs", "sampling", "n"}] * 2
    x = lines[0]
    assert (x["id"], x["image"], x["n"]) == (rows[0]["key"], rows[0]["image"], 8)
    assert x["prompt"] == f"<user><image>{PROMPTS['lidc']}<assistant>"
    assert x["processor_kwargs"] == {}  # configs.processor_kwargs minus padding_side
    assert configs.processor_kwargs("qwen3-vl-2b", "lidc")["padding_side"] == "left"
    assert x["sampling"] == configs.sampling("qwen3-vl-2b") | {"seed": 3, "stop_token_ids": [1, 2]}


def test_export_caps_the_completion_only_when_asked(rows, tmp_path):
    """CLEAN8: select_ckpt's max_tokens (the CPU smoke's cap on a random Qwen3-VL) replaces the model's budget."""
    from objdel.policy import models as configs
    assert export(rows, tmp_path, max_tokens=32)[1][0]["sampling"]["max_tokens"] == 32
    assert export(rows, tmp_path)[1][0]["sampling"]["max_tokens"] == configs.MAX_COMPLETION["qwen3-vl-2b"]


def test_export_takes_image_paths_from_a_dataset(rows, tmp_path):
    datasets = pytest.importorskip("datasets")
    from PIL import Image
    for r in rows:
        Image.new("RGB", (4, 4)).save(r["image"])
    ds = datasets.Dataset.from_list(rows).cast_column("image", datasets.Image())  # as grpo_rows
    _, lines = export(ds, tmp_path)
    assert [x["image"] for x in lines] == [r["image"] for r in rows]
    with pytest.raises(ValueError, match="twice"):
        export([rows[0], rows[0]], tmp_path)


class SamplingParams:
    def __init__(self, **kw):
        self.kw = kw

    def __repr__(self):
        return f"SamplingParams({self.kw})"


def completion(text, ids, finish):
    return types.SimpleNamespace(text=text, token_ids=ids, finish_reason=finish, logprobs=None)


def test_requests_echo_the_exported_sampling(V, rows, tmp_path):
    _, lines = export(rows, tmp_path, stop_ids={1})
    opened = []
    reqs, params = V.requests(lines, SamplingParams, lambda p: opened.append(p) or f"image:{p}")
    assert opened == [r["image"] for r in rows]
    assert reqs[0] == {"prompt": lines[0]["prompt"], "multi_modal_data": {"image": f"image:{rows[0]['image']}"}}
    seeds = [V.request_seed(x["sampling"]["seed"], x["id"]) for x in lines]
    assert [p.kw for p in params] == [{"n": 8, **x["sampling"], "seed": s} for x, s in zip(lines, seeds, strict=True)]
    assert seeds[0] != seeds[1] and seeds == [V.request_seed(0, x["id"]) for x in lines]  # per item, reproducible
    assert V.request_seed(1, lines[0]["id"]) != seeds[0] and all(0 <= s < 2**31 for s in seeds)
    assert V.requests(lines, SamplingParams, str, logprobs=True)[1][0].kw["logprobs"] == 0


def test_round_trip_and_truncation_is_the_rewards_rule(V, rows, tmp_path):
    from objdel.policy import models as rewards
    from objdel.policy import sampler as sample
    _, lines = export(rows, tmp_path)
    stop = {1, 2}
    outs = [types.SimpleNamespace(outputs=[completion(f"t{i}{k}", ids, f) for k, (ids, f) in enumerate(
        [([5, 1], "stop"), ([5, 6], "length"), ([], "length"), ([7, 2], "stop")] * 2)]) for i in range(2)]
    out = tmp_path / "completions.jsonl"
    V.write(V.outputs(lines, outs), out)
    recs = [json.loads(line) for line in out.read_text().splitlines()]
    assert set(recs[0]) == {"id", "k", "text", "token_ids", "finish_reason"} and len(recs) == 16
    got = sample.import_completions(out, stop)
    assert [g["item_id"] for g in got] == [r["key"] for r in rows]
    assert got[0]["samples"] == [f"t0{k}" for k in range(8)]
    assert got[1]["truncated"] == [rewards.truncated(c.token_ids, stop) for c in outs[1].outputs]
    assert got[0]["truncated"][:4] == [False, True, True, False]
    assert got[0]["tokens"][:4] == [2, 2, 0, 2] and set(got[0]) == {"item_id", "samples", "truncated", "tokens"}
    shuffled = tmp_path / "shuffled.jsonl"
    shuffled.write_text("".join(reversed(out.read_text().splitlines(keepends=True))))
    key = lambda r: r["item_id"]  # noqa: E731
    assert sorted(sample.import_completions(shuffled, stop), key=key) == sorted(got, key=key)  # samples by k
    shuffled.write_text("".join(out.read_text().splitlines(keepends=True)[1:]))  # k = 0 of the first id lost
    with pytest.raises(ValueError, match="samples"):
        sample.import_completions(shuffled, stop)


def test_import_refuses_short_duplicate_or_missing_items(tmp_path):
    """BK2-3 / BK3-7: every exported id comes back with exactly n samples (default G = 8), each (id, k) once."""
    pytest.importorskip("trl")
    from objdel.policy import sampler as sample

    def run(recs, **kw):
        path = tmp_path / "c.jsonl"
        path.write_text("".join(json.dumps({"id": i, "k": k, "text": "", "token_ids": [1],
                                            "finish_reason": "stop"}) + "\n" for i, k in recs))
        return sample.import_completions(path, {1}, **kw)

    full = [(i, k) for i in "ab" for k in range(8)]
    assert [r["item_id"] for r in run(full, ids=["b", "a"])] == ["a", "b"]
    assert len(run(full[:3] + full[8:11], n=3)[1]["samples"]) == 3
    for recs, kw, msg in ((full[:-1], {}, "b: 7 samples"), (full + [("a", 0)], {}, "a: sample k = 0 twice"),
                          (full, {"ids": ["a", "b", "c"]}, r"missing: \['c'\]"),
                          (full, {"ids": ["a"]}, r"not exported: \['b'\]"), (full, {"n": 16}, "8 samples")):
        with pytest.raises(ValueError, match=msg):
            run(recs, **kw)


def test_refuses_to_start_on_a_busy_gpu(V):
    V.check_memory(lambda: (95, 100))
    with pytest.raises(SystemExit, match="89"):
        V.check_memory(lambda: (89, 100))


def test_parity_help():
    from experiments.phase2 import parity
    with pytest.raises(SystemExit) as e:
        parity.main(["--help"])
    assert e.value.code == 0


@pytest.mark.slow
def test_logprob_parity_on_a_tiny_model(tmp_path):
    """S1-16: experiments/phase2/parity.py teacher-forces sample_vllm --logprobs completions. The "vLLM" logprobs
    here come from the lab's path (apply_chat_template(tokenize=True) on the conversation, B5 hf_teacher), so a pass
    shows the exported prompt text gives the same inputs; a 0.1-nat shift per token fails."""
    torch = pytest.importorskip("torch")
    pytest.importorskip("trl")
    from PIL import Image
    from transformers import AutoModelForImageTextToText, AutoProcessor
    from trl.data_utils import prepare_multimodal_messages

    from experiments.phase2 import parity
    from objdel.policy import sampler as sample
    from objdel.policy.prompt import _user
    from tests.test_p2_scripts import _tiny_vlm

    model_dir = _tiny_vlm(tmp_path / "tiny")
    proc = AutoProcessor.from_pretrained(model_dir)
    model = AutoModelForImageTextToText.from_pretrained(model_dir, dtype=torch.float32)  # "auto" would be bf16
    rng = np.random.default_rng(0)
    rows = []
    for i in range(2):
        Image.fromarray(rng.integers(0, 256, (64, 64, 3), np.uint8)).save(tmp_path / f"im{i}.png")
        rows.append({"prompt": [_user("lidc")], "image": str(tmp_path / f"im{i}.png"), "key": f"k{i}"})
    sample.export_prompts(rows, tmp_path / "p.jsonl", "qwen3-vl-2b", "lidc", 2, processor=proc)
    comps = []
    for r in rows:
        conv = prepare_multimodal_messages(r["prompt"], images=[Image.open(r["image"]).convert("RGB")])
        enc = proc.apply_chat_template([conv], add_generation_prompt=True, tokenize=True, return_dict=True,
                                       return_tensors="pt")
        n = enc["input_ids"].shape[1]
        for k in range(2):
            ids = rng.integers(10, 1000, 6).tolist()
            kw = dict(enc, input_ids=torch.cat([enc["input_ids"], torch.tensor([ids])], 1),
                      attention_mask=torch.ones(1, n + len(ids), dtype=torch.long))
            with torch.no_grad():
                lp = torch.log_softmax(model(**kw).logits[0, n - 1:-1].float(), -1)[range(len(ids)), ids]
            comps.append({"id": r["key"], "k": k, "text": "", "token_ids": ids, "finish_reason": "length",
                          "logprobs": lp.tolist(), "prompt_token_ids": enc["input_ids"][0].tolist()})
    base = ["--prompts", str(tmp_path / "p.jsonl"), "--model", str(model_dir), "--out", str(tmp_path / "r.json")]
    for shift, ok in ((0.0, True), (0.1, False)):
        (tmp_path / "c.jsonl").write_text("".join(json.dumps(dict(c, logprobs=[v + shift for v in c["logprobs"]]))
                                                  + "\n" for c in comps))
        r = parity.main(base + ["--completions", str(tmp_path / "c.jsonl")])
        assert (r["pass"], r["tokens"], r["sequences"]) == (ok, 24, 4) and r["max_abs"] < shift + 1e-5
        assert r["prompt_ids"] == {"items": 2, "equal": 2, "mismatched": [], "pass": True}
    assert json.loads((tmp_path / "r.json").read_text())["median_abs"] == pytest.approx(0.1, abs=1e-4)
