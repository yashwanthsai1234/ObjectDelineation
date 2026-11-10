"""policy.caches.PromptMemo (C2 step 5): TRL pass 2's processor call with the text part once per (text, image sizes).

The output must be the processor's own, bit for bit, so the RL step does not change: U-PM1..PM5, P-PM1, I-PM1 on
the real Qwen3-VL processor (skipped when it is not in the local HF cache), U-IM8 (the ImageMemo under it), and a
two-step tiny FastGRPO run with and without the memo that must give identical logs and weights."""

import copy
import hashlib
import pickle
import threading

import numpy as np
import PIL.Image
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from transformers import AutoProcessor  # noqa: E402

from objdel.policy import caches  # noqa: E402

PAD = "<|vision_start|><|image_pad|><|vision_end|>"


def load():
    try:
        return AutoProcessor.from_pretrained("Qwen/Qwen3-VL-2B-Instruct", local_files_only=True)
    except OSError:
        pytest.skip("the Qwen3-VL processor is not in the local HF cache")


def image(seed, h=64, w=64, mode="RGB"):
    img = PIL.Image.fromarray(np.random.default_rng(seed).integers(0, 256, (h, w, 3), dtype=np.uint8))
    return img.convert(mode)


def installed():
    """(processor with ImageMemo and PromptMemo, a never-installed copy, the memo, a list of original-call sizes)."""
    proc, plain = load(), load()
    seen, real = [], type(proc).__call__

    class Counting(type(proc)):
        def __call__(self, *a, **k):
            seen.append(len(k.get("text") or []))
            return real(self, *a, **k)

    Counting.__name__ = Counting.__qualname__ = type(proc).__name__
    proc.__class__ = Counting
    memo = caches.PromptMemo.install(proc, caches.ImageMemo.install(proc))
    return proc, plain, memo, seen


def call(p, texts, imgs):
    return p(images=[[x] for x in imgs], text=texts, padding=True, return_tensors="pt")


def assert_same(got, want):
    assert list(got.keys()) == list(want.keys())
    for k in want:
        assert got[k].dtype == want[k].dtype and torch.equal(got[k], want[k]), k


def test_u_pm1_other_families_and_text_only_classes_get_nothing():
    class Other:
        image_processor = object()

    class Text:
        pass

    for p in (Other(), Text()):
        image_memo = caches.ImageMemo.install(p)  # the trainer installs it first; PromptMemo must add nothing
        before, state = type(p), copy.deepcopy(p.__dict__)
        assert caches.PromptMemo.install(p, image_memo) is None
        assert type(p) is before and p.__dict__.keys() == state.keys()
    proc = load()
    assert caches.PromptMemo.install(proc, None) is None and type(proc).__name__ == "Qwen3VLProcessor"


def test_u_pm2_isinstance_attributes_saved_files_copies_and_one_wrap(tmp_path):
    proc, plain = load(), load()
    mro, tok, tmpl = type(proc).__mro__, proc.tokenizer, proc.chat_template
    memo = caches.PromptMemo.install(proc, caches.ImageMemo.install(proc))
    assert memo is not None and all(isinstance(proc, c) for c in mro)
    assert proc.tokenizer is tok and proc.chat_template == tmpl
    assert caches.PromptMemo.install(proc, proc.image_processor) is memo and type(proc).__mro__[1] is mro[0]
    plain.save_pretrained(tmp_path / "a")
    proc.save_pretrained(tmp_path / "b")
    digest = lambda d: {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(d.iterdir())}  # noqa: E731
    assert digest(tmp_path / "a") == digest(tmp_path / "b")
    for back in (copy.deepcopy(proc), pickle.loads(pickle.dumps(proc))):
        assert type(back) is mro[0]
        assert_same(call(back, [PAD + "x"], [image(0)]), call(plain, [PAD + "x"], [image(0)]))


def test_u_pm3_other_signatures_run_the_processor_and_store_nothing():
    proc, plain, memo, seen = installed()
    img, msg = image(0), [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": "x"}]}]
    cases = [dict(text=[PAD + "x"]), dict(images=[img], return_tensors="pt"),
             dict(images=[[img]], text=[PAD + "x"], padding="max_length", max_length=40, return_tensors="pt"),
             dict(images=[[img]], text=[PAD + "x"], padding=True, return_tensors=None),
             dict(images=[[img]], text=[PAD + "x"], padding=True, return_tensors="pt", add_special_tokens=False),
             dict(images=[img], text=PAD + "x", padding=True, return_tensors="pt"),
             dict(images=img, text=[PAD + "x"], padding=True, return_tensors="pt")]  # fmt: skip
    for kw in cases:
        n = len(seen)
        got, want = proc(**kw), plain(**kw)
        assert len(seen) == n + 1 and memo.rows == {}
        for k in want:
            assert np.array_equal(np.asarray(got[k]), np.asarray(want[k])), (kw, k)
    assert proc.apply_chat_template(msg, add_generation_prompt=True) == plain.apply_chat_template(
        msg, add_generation_prompt=True)


def test_u_pm4_the_original_runs_once_per_key_and_a_failure_is_not_cached(monkeypatch):
    proc, _, memo, seen = installed()
    imgs = [image(i) for i in range(32)]
    call(proc, [PAD + "x"] * 256, [imgs[i % 32] for i in range(256)])
    proc.image_processor.clear()
    call(proc, [PAD + "x"] * 256, [imgs[i % 32] for i in range(256)])
    assert seen.count(1) == 1  # one row call (key miss); the rest are the batch's image steps (text None: 0)
    seen.clear()
    big = [image(i, 96, 64) for i in range(2)]
    call(proc, [PAD + "x", PAD + "y", PAD + "x", PAD + "y"], [big[0], big[1], imgs[0], imgs[1]])
    assert seen.count(1) == 3  # (x, 96 x 64), (y, 96 x 64), (y, 64 x 64) are new; (x, 64 x 64) is not
    real, n = type(proc).__mro__[1].__call__, [0]

    def flaky(self, *a, **k):
        n[0] += 1
        if n[0] == 1:
            raise RuntimeError("tokenizer down")
        return real(self, *a, **k)

    monkeypatch.setattr(memo, "call", flaky)
    with pytest.raises(RuntimeError, match="tokenizer down"):
        call(proc, [PAD + "z"], [imgs[0]])
    assert (PAD + "z", ((64, 64),)) not in memo.rows
    call(proc, [PAD + "z"], [imgs[0]])
    assert n[0] == 3  # the retry runs the row again (2) and the image step (3)


@pytest.mark.parametrize("side", ["left", "right"])
def test_u_pm5_rows_are_padded_as_the_tokenizer_pads_and_results_do_not_alias(side):
    proc, plain, memo, _ = installed()
    proc.tokenizer.padding_side = plain.tokenizer.padding_side = side
    texts, imgs = [PAD + "a", PAD + "a much longer prompt text"] * 2, [image(0), image(1), image(2, 96, 64), image(0)]
    got = call(proc, texts, imgs)
    assert_same(got, call(plain, texts, imgs))
    got["input_ids"][0, -1] = 7
    assert_same(call(proc, texts, imgs), call(plain, texts, imgs))
    proc.image_processor.per_image = True
    assert len(call(proc, texts, imgs)["image_grid_thw"]) == 4


def test_p_pm1_mixed_batches_equal_the_processor():
    """150 cases: 4 texts x 3 image sizes, random rows per case; every output key equals the plain processor's."""
    proc, plain, memo, seen = installed()
    rng = np.random.default_rng(0)
    texts = [PAD + t for t in ("outline", "outline the optic disc", "cup", "outline the nodule in this view")]
    sizes = [(64, 64), (96, 64), (64, 128)]
    pool = [image(i, *sizes[i % 3], mode="L" if i % 5 == 0 else "RGB") for i in range(12)]
    for _ in range(150):
        n = int(rng.integers(1, 9))
        t = [texts[i] for i in rng.integers(0, 4, n)]
        im = [pool[i] for i in rng.integers(0, 12, n)]
        assert_same(call(proc, t, im), call(plain, t, im))
    assert seen.count(1) == len(memo.rows) <= 12  # one row call per (text, size) key


def test_i_pm1_real_prompts_and_the_digit_image_equal_the_processor():
    proc, plain, _, _ = installed()
    msg = lambda q: [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": q}]}]  # noqa: E731
    texts = [plain.apply_chat_template(msg(q), add_generation_prompt=True) for q in
             ("Outline the nodule.", "Outline the optic disc.", "Outline the cup.", "Give the outline as points.",
              "Outline the lesion in this CT view.", "Draw the boundary.", "Trace the region.", "Outline it.")]
    digit = PIL.Image.fromarray((np.arange(64 * 64 * 3) & 255).astype(np.uint8).reshape(64, 64, 3))
    imgs = [image(i, 512, 512) for i in range(8)]
    assert_same(call(proc, texts + texts[:1], imgs + [digit]), call(plain, texts + texts[:1], imgs + [digit]))


def test_u_im8_no_double_count_inside_prompt_memo_and_other_threads_raise():
    proc, _, _, _ = installed()
    memo, imgs = proc.image_processor, [image(i) for i in range(4)]
    h, m = memo.hits, memo.misses
    call(proc, [PAD + "x"] * 8, [imgs[i % 4] for i in range(8)])
    assert (memo.hits - h, memo.misses - m) == (4, 4)  # 8 rows counted once: the row calls skip the memo
    err = []
    t = threading.Thread(target=lambda: err.append(pytest.raises(RuntimeError, memo, imgs, return_tensors="pt")))
    t.start()
    t.join()
    assert str(memo.owner) in str(err[0].value)
    t = threading.Thread(target=lambda: err.append(pytest.raises(RuntimeError, memo.prefill, imgs)))
    t.start()
    t.join()
    assert len(err) == 2


def test_u_im8_the_reward_prefetch_never_calls_the_image_memo(monkeypatch):
    from objdel.rl import trl_adapter

    monkeypatch.setattr(caches.ImageMemo, "__call__", lambda *a, **k: pytest.fail("memo called"))
    monkeypatch.setattr(caches.ImageMemo, "prefill", lambda *a, **k: pytest.fail("memo prefilled"))
    src = __import__("inspect").getsource(trl_adapter.FastGRPO.prefetch_rewards)
    assert "images" not in src and "processing_class" not in src


# The RL step, with and without the memo ----------------------------------------------------------------------------


def tiny_vlm_run(tmp_path, mode: str, steps: int = 2):
    """Two FastGRPO steps (CPU, fp32, TRL's own generate) of a tiny random Qwen3-VL with the real processor: 4 rows
    of two image sizes x 4 answers a step. mode: "off" (flag off), "flag" (the --prompt-memo flag as at 13f5cb7:
    pass-1 override and prefill, no PromptMemo) or "memo" (use_prompt_memo: the flag and PromptMemo). Returns
    (train logs, parameters, the PromptMemo's row-call count); steps 0: the untrained parameters."""
    from datasets import Dataset, Features, Image, Value
    from transformers import Qwen3VLConfig, Qwen3VLForConditionalGeneration
    from trl import GRPOConfig

    from objdel.rl import trl_adapter

    proc = load()
    cfg = proc.tokenizer
    cfg.padding_side = "left"  # as the run loads it (policy.models: padding_side left; TRL pads prompt ids left)
    torch.manual_seed(0)
    model = Qwen3VLForConditionalGeneration(Qwen3VLConfig(
        text_config=dict(vocab_size=len(cfg), hidden_size=32, intermediate_size=64, num_hidden_layers=2,
                         num_attention_heads=2, num_key_value_heads=1, head_dim=16,
                         rope_scaling={"rope_type": "default", "mrope_section": [2, 3, 3], "mrope_interleaved": True}),
        vision_config=dict(depth=1, hidden_size=32, intermediate_size=64, num_heads=2, out_hidden_size=32,
                           patch_size=16, spatial_merge_size=2, temporal_patch_size=2, deepstack_visual_indexes=[0]),
        image_token_id=proc.image_token_id, video_token_id=proc.video_token_id,
        vision_start_token_id=cfg.convert_tokens_to_ids("<|vision_start|>"),
        vision_end_token_id=cfg.convert_tokens_to_ids("<|vision_end|>")))  # fmt: skip
    model.generation_config.eos_token_id, model.generation_config.pad_token_id = cfg.eos_token_id, cfg.pad_token_id
    msg = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": "Outline the nodule."}]}]
    rows = Dataset.from_dict({"prompt": [msg] * 4, "image": [image(i, *[(64, 64), (96, 64)][i % 2]) for i in range(4)]},
                             features=Features({"prompt": [{"role": Value("string"), "content": [{
                                 "type": Value("string"), "text": Value("string")}]}], "image": Image()}))  # fmt: skip
    args = GRPOConfig(output_dir=str(tmp_path / mode), max_steps=max(steps, 1), per_device_train_batch_size=16,
                      num_generations=4, max_completion_length=4, use_cpu=True, bf16=False, report_to="none",
                      logging_steps=1, save_strategy="no", seed=0, remove_unused_columns=False, learning_rate=1e-2,
                      beta=0.0, loss_type="dr_grpo", use_vllm=False)  # fmt: skip
    t = trl_adapter.FastGRPO(model=model, reward_funcs=lambda completions, **kw: [len(str(c)) % 3 * 1.0 for c in
                             completions], args=args, train_dataset=rows, processing_class=proc)  # fmt: skip
    if mode == "flag":
        t.prompt_memo = True
    elif mode == "memo":
        assert t.use_prompt_memo()
    if steps:
        t.train()
    logs = [{k: v for k, v in h.items() if "time" not in k and "per_second" not in k}  # wall clock: time/*, step_time
            for h in t.state.log_history]  # fmt: skip
    memo = type(proc).__dict__.get("prompt_memo")
    return logs, {n: p.detach().clone() for n, p in t.model.named_parameters()}, memo and len(memo.rows)


def test_rl_step_is_bit_identical_with_and_without_prompt_memo(tmp_path):
    """PromptMemo changes no bit of the RL step: two steps with it give the logs (loss, rewards, grad norm, lengths)
    and every weight of the flag-off run and of the flag as it was (pass-1 override only); it served pass 2 from 2
    row calls (one prompt text, two image sizes)."""
    off_logs, off_w, _ = tiny_vlm_run(tmp_path, "off")
    flag_logs, flag_w, _ = tiny_vlm_run(tmp_path, "flag")
    memo_logs, memo_w, rows = tiny_vlm_run(tmp_path, "memo")
    assert rows == 2 and len(off_logs) >= 2 and "loss" in off_logs[0]
    assert memo_logs == off_logs == flag_logs
    assert all(torch.equal(memo_w[n], off_w[n]) and torch.equal(flag_w[n], off_w[n]) for n in off_w)
    assert any(not torch.equal(off_w[n], w) for n, w in tiny_vlm_run(tmp_path, "off", steps=0)[1].items())
