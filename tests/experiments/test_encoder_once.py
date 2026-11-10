"""C3 (--encoder-once): the generation batch's distinct images run through the frozen vision tower once; every forward
gets each row's features as mm_encoder_outputs, never pixel values (plan C3, tests U-EN1..U-EN7, I-EN1..I-EN4, P-EN1).

The reference is the pixel_values path with the same grouping: VisionCache holds the distinct images from one tower
call, so the equality checks the plumbing, not the GEMM shape (C3 step 0 measures that on the box)."""

import dataclasses

import numpy as np
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from transformers.models.qwen3_vl.modeling_qwen3_vl import Qwen3VLModel  # noqa: E402

from objdel.policy import caches as policy_caches  # noqa: E402
from objdel.rl import trl_adapter as rl_trl_adapter  # noqa: E402
from tests.experiments.test_trl_adapter_equivalence import (  # noqa: E402
    PATCH,
    Fast,
    GRPOTrainer,
    group_batch,
    meta_for,
    tiny_model,
    trainer,
)
from tests.policy.test_caches import counted_tower  # noqa: E402

G = 4
BUCKETS = (4, 8, 12, 16, 24)


def four_groups():
    """4 prompts x 4 answers (answers of 0 and 1 tokens included), one distinct random image per prompt."""
    return group_batch(groups=4, G=G, lengths=[0, 1, 5, 3, 7, 2, 4, 6, 10, 1, 3, 8, 2, 9, 0, 5], pads=[2, 0, 1, 3])


def distinct(kw):
    """The batch's distinct images (one per group of G rows): pixel rows and grid."""
    n = len(kw["num_images"])
    return kw["pixel_values"].view(n, -1, PATCH)[::G].reshape(-1, PATCH), kw["image_grid_thw"][::G]


def base_of(m):
    return next(x for x in m.modules() if isinstance(x, Qwen3VLModel))


def pixel_trainer(m, kw, **attrs):
    """The pixel_values path with the encoder-once grouping: one tower call on the distinct images fills VisionCache."""
    vc = policy_caches.VisionCache.install(m)
    vc.forward(*distinct(kw))
    return trainer(Fast, buckets=BUCKETS, vision=vc, **attrs)


def frozen(m):
    """m with its vision tower frozen, as configs.lora leaves it (a bare tiny model trains every weight)."""
    base_of(m).visual.requires_grad_(False)
    return m


class MetaFast(Fast):
    """Fast with FastGRPO's rule: a forward given no host_meta (reference and old passes) takes the generation
    batch's RowMeta `meta`."""

    meta = None

    def _get_per_token_logps_and_entropies(self, *args, **kw):
        if self.meta is not None:
            kw.setdefault("host_meta", self.meta)
        return super()._get_per_token_logps_and_entropies(*args, **kw)


def encoder_trainer(m, kw, meta, **attrs):
    """The encoder-once path; `meta` is the generation batch's RowMeta (its image_id field gives each row's image)."""
    enc = policy_caches.EncoderCache(policy_caches.VisionCache.install(frozen(m)), base_of(m))
    enc.encode(*distinct(kw))
    return trainer(MetaFast, buckets=BUCKETS, encoder=enc, meta=meta, **attrs)


def shared_attrs(m, on):
    return {"shared_prefix": True, "rope": policy_caches.RopeCache.install(m)} if on else {}


def calls(seen):
    return len(seen)


@pytest.mark.parametrize("shared_prefix", [False, True])
@pytest.mark.parametrize("attn", ["sdpa", "eager"])
@pytest.mark.parametrize("chunk", [8, 4])
def test_encoder_once_logps_and_entropies_equal_the_pixel_path(shared_prefix, attn, chunk):
    """I-EN1: torch.equal on every position; the tower runs once per generation batch on the new path (reference and
    update passes), and more than once on the old path's chunked passes (the counter cannot always read 1)."""
    ids, mask, T, kw = four_groups()
    out, towers = {}, {}
    for path in ("pixel", "encoder"):
        m = tiny_model()
        m.set_attn_implementation(attn)
        seen = counted_tower(m)
        attrs = shared_attrs(m, shared_prefix)
        t = (pixel_trainer(m, kw, **attrs) if path == "pixel"
             else encoder_trainer(m, kw, meta_for(mask, T, G), **attrs))  # fmt: skip
        for _ in range(2):  # reference pass and an update pass of one generation batch
            out[path] = t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=chunk, compute_entropy=True,
                                                             **kw)[:2]  # fmt: skip
        towers[path] = calls(seen)
    assert torch.equal(out["encoder"][0], out["pixel"][0]) and torch.equal(out["encoder"][1], out["pixel"][1])
    assert towers["encoder"] == 1
    m = tiny_model()
    seen = counted_tower(m)
    policy_caches.VisionCache.install(m)
    trainer(Fast, buckets=BUCKETS)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=chunk, **kw)
    assert calls(seen) > 1  # control: today's path encodes per chunk


def test_encoder_once_lora_gradients_equal_the_pixel_path():
    """I-EN2: one micro-batch, every LoRA A/B grad torch.equal; the tower has no grad and the features none."""
    from peft import get_peft_model

    from objdel.policy import models as configs

    ids, mask, T, kw = four_groups()
    real = mask[:, -T:].float()
    grads = []
    for path in ("pixel", "encoder"):
        m = get_peft_model(tiny_model(), configs.lora("qwen3-vl-2b"))
        g = torch.Generator().manual_seed(1)
        for n, p in m.named_parameters():
            if "lora_" in n:
                p.data = torch.randn(p.shape, generator=g) * 0.05
        t = (pixel_trainer(m, kw) if path == "pixel"
             else encoder_trainer(m, kw, meta_for(mask, T, G)))  # fmt: skip
        lp, _, _ = t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=16, **kw)
        (lp * real).sum().backward()
        grads.append({n: p.grad for n, p in m.named_parameters() if p.grad is not None})
        assert not [n for n in grads[-1] if ".visual." in n]
    assert grads[0].keys() == grads[1].keys() and all("lora_" in n for n in grads[0])
    assert all(torch.equal(grads[0][n], grads[1][n]) for n in grads[0])
    assert not any(x.requires_grad for x in t.encoder.outputs([0, 1]).pooler_output)


def test_every_forward_gets_encoder_outputs_and_never_pixel_values():
    """U-EN2 and I-EN3: a spy on the base's forward sees mm_encoder_outputs {"image"} and no pixel_values for a
    reference pass (2 chunks) and 4 update micro-batches; the tower ran once, on the 4 distinct images (4 x 16 rows,
    not 16 x 16); a new generation batch encodes again (2). With the path off the spy sees pixel_values only."""
    ids, mask, T, kw = four_groups()
    for on in (True, False):
        m = tiny_model()
        seen, rows = counted_tower(m), []
        orig = m.model.visual.forward
        m.model.visual.forward = lambda pv, grid_thw=None, _o=orig, _r=rows, **k: _r.append(len(pv)) or _o(pv, grid_thw,
                                                                                                        **k)
        t = (encoder_trainer(m, kw, meta_for(mask, T, G)) if on
             else trainer(Fast, buckets=BUCKETS, vision=policy_caches.VisionCache.install(m)))  # fmt: skip
        got = []
        m.model.register_forward_pre_hook(lambda mod, a, k, _g=got: _g.append((k.get("mm_encoder_outputs"),
                                                                               k.get("pixel_values"))),
                                          with_kwargs=True)  # fmt: skip
        t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=8, **kw)
        for s in range(0, 16, 4):  # update micro-batches: shuffled subsets carrying their host_meta
            perm = np.random.default_rng(s).permutation(16)[:4]
            meta = meta_for(mask, T, G)
            sub = {"pixel_values": kw["pixel_values"].view(16, -1, PATCH)[perm].reshape(-1, PATCH),
                   "image_grid_thw": kw["image_grid_thw"][perm], "num_images": [1] * 4,
                   "mm_token_type_ids": kw["mm_token_type_ids"][perm]}  # fmt: skip
            pick = policy_caches.RowMeta(*(a[perm] for a in dataclasses.astuple(meta)))
            t._get_per_token_logps_and_entropies(m, ids[perm], mask[perm], T, host_meta=pick, **sub)
        assert len(got) == 6
        if on:
            assert all(e is not None and set(e) == {"image"} and p is None for e, p in got)
            assert calls(seen) == 1 and rows[0] == 4 * 16
            t.encoder.encode(*distinct(kw))
            assert calls(seen) == 2
        else:
            assert all(e is None and p is not None for e, p in got)


def test_image_ids_follow_shuffled_rows_and_a_wrong_permutation_changes_the_logps():
    """P-EN1: an update micro-batch of shuffled rows takes each row's own features through host_meta (equal to the
    pixel path on the same rows); permuting the id field but not the rows changes the logprobs (the oracle has
    teeth: the images differ a lot)."""
    ids, mask, T, kw = four_groups()
    perm = np.random.default_rng(3).permutation(16)
    sub = {"pixel_values": kw["pixel_values"].view(16, -1, PATCH)[perm].reshape(-1, PATCH),
           "image_grid_thw": kw["image_grid_thw"][perm], "num_images": [1] * 16,
           "mm_token_type_ids": kw["mm_token_type_ids"][perm]}  # fmt: skip
    meta = meta_for(mask, T, G)
    pick = policy_caches.RowMeta(*(a[perm] for a in dataclasses.astuple(meta)))
    m0 = tiny_model()
    want, _, _ = pixel_trainer(m0, kw)._get_per_token_logps_and_entropies(
        m0, ids[perm], mask[perm], T, batch_size=4, host_meta=pick, **sub)
    m = tiny_model()
    t = encoder_trainer(m, kw, None)
    got, _, _ = t._get_per_token_logps_and_entropies(m, ids[perm], mask[perm], T, batch_size=4, host_meta=pick, **sub)
    assert torch.equal(got, want)
    wrong = dataclasses.replace(pick, image_id=pick.image_id[np.random.default_rng(4).permutation(16)])
    bad, _, _ = t._get_per_token_logps_and_entropies(m, ids[perm], mask[perm], T, batch_size=4, host_meta=wrong,
                                                     **sub)  # fmt: skip
    assert not torch.equal(bad, want)


@pytest.mark.parametrize("image_ids, error", [(None, RuntimeError), (np.arange(15) // G, ValueError),
                                              (np.r_[np.arange(15) // G, 99], KeyError),
                                              (np.r_[-1, np.arange(15) // G], KeyError)])  # fmt: skip
def test_missing_or_bad_image_ids_raise(image_ids, error):
    """U-EN3: ids missing, of length B - 1, 99 of 4, -1: the forward raises naming "image" and never falls back to
    pixel values (no tower call afterwards). With the path off and no ids the forward works."""
    ids, mask, T, kw = four_groups()
    m = tiny_model()
    seen = counted_tower(m)
    meta = None if image_ids is None else dataclasses.replace(meta_for(mask, T, G), image_id=image_ids.astype(np.int32))
    t = encoder_trainer(m, kw, meta)
    with pytest.raises(error, match="image"):
        t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=8, **kw)
    assert calls(seen) == 1 and t.encoder.pending is None
    m = tiny_model()
    trainer(Fast, buckets=BUCKETS, vision=policy_caches.VisionCache.install(m))._get_per_token_logps_and_entropies(
        m, ids, mask, T, batch_size=8, **kw)


def test_text_only_rows_skip_the_encoder_and_a_token_mismatch_is_an_error():
    """U-EN4: an all-text batch never calls the encoder and equals TRL; a grid of 4 merged tokens against 3 image
    pads in the prompt raises the model's own error, no truncation."""
    ids, mask, T, kw = four_groups()
    m = tiny_model()
    seen = counted_tower(m)
    t = encoder_trainer(m, kw, meta_for(mask, T, G))
    text = torch.where(ids == 100, 7, ids)
    want, _, _ = trainer(GRPOTrainer)._get_per_token_logps_and_entropies(m, text, mask, T, batch_size=8)
    got, _, _ = t._get_per_token_logps_and_entropies(m, text, mask, T, batch_size=8)
    torch.testing.assert_close(got, want, atol=1e-5, rtol=1e-5)
    assert calls(seen) == 1  # the encode call only
    short = ids.clone()
    first = (short == 100).int().argmax(1)
    short[torch.arange(16), first] = 7  # 3 image pads per row, the grid says 4
    with pytest.raises(ValueError, match="Image features and image tokens do not match"):
        t._get_per_token_logps_and_entropies(m, short, mask, T, batch_size=8,
                                             **(kw | {"mm_token_type_ids": (short == 100).long()}))


def test_per_image_follows_the_flag_and_is_cleared_before_each_generation(monkeypatch):
    """U-EN5: the image memo returns one block per distinct image only under --encoder-once; set and cleared before
    every generation, per trainer (two trainers sharing one processor each get their own setting)."""
    memo = policy_caches.ImageMemo(lambda images, **k: None)
    memo.memo["old"] = None
    out = ("pids", [[1]], None, ["a"], None, None, {}, None, [])
    monkeypatch.setattr(GRPOTrainer, "_generate", lambda self, prompts: out)
    on = trainer(rl_trl_adapter.FastGRPO, images=memo, encoder=object(), use_vllm=False, reward_funcs=[])
    off = trainer(rl_trl_adapter.FastGRPO, images=memo, use_vllm=False, reward_funcs=[])
    on._generate(["p"])
    assert memo.per_image is True and memo.memo == {}
    off._generate(["p"])
    assert memo.per_image is False


def test_a_trainable_tower_is_refused():
    """U-EN6: features computed once are valid only for a frozen tower."""
    m = tiny_model()
    for p in m.model.visual.parameters():
        p.requires_grad_(False)
    next(m.model.visual.parameters()).requires_grad_(True)
    with pytest.raises(ValueError, match="frozen"):
        policy_caches.EncoderCache(policy_caches.VisionCache.install(m), m.model)


def test_shape_guard_raises_on_pixels_cut_by_the_row_grid():
    """U-EN7: under --encoder-once a tower call with 4 images of pixel rows and a 16-row grid raises EncoderMismatch;
    equal shapes pass; with the flag off the tower wrapper is VisionCache's."""
    ids, mask, T, kw = four_groups()
    pv, grid = distinct(kw)
    m = tiny_model()
    vc = policy_caches.VisionCache.install(frozen(m))
    assert m.model.visual.forward == vc.forward
    policy_caches.EncoderCache(vc, m.model)
    with pytest.raises(policy_caches.EncoderMismatch):
        m.model.visual(pv, grid_thw=kw["image_grid_thw"])
    assert m.model.visual(pv, grid_thw=grid).pooler_output.shape == (16, 32)


def tiny_vl_trainer(monkeypatch, tmp_path, proc, encoder_once):
    """FastGRPO on a tiny Qwen3-VL with the real Qwen3-VL processor's token ids, TRL's vLLM path with FakeVLLM, a
    PEFT adapter and beta > 0 (TRL's reference pass runs inside _generate_and_score_completions)."""
    from datasets import Dataset
    from trl import GRPOConfig

    from objdel.policy import models as configs
    from tests.experiments.test_trl_adapter_equivalence import FakeVLLM

    monkeypatch.setattr("trl.trainer.grpo_trainer.VLLMGeneration", FakeVLLM)
    monkeypatch.setattr(rl_trl_adapter.FastGRPO, "sync_without_merge", lambda self: False)  # FakeVLLM has no engine
    tok = proc.tokenizer
    ids = [tok.convert_tokens_to_ids(t) for t in ("<|image_pad|>", "<|video_pad|>", "<|vision_start|>",
                                                  "<|vision_end|>")]  # fmt: skip
    model = tiny_model(len(tok), *ids)
    args = GRPOConfig(output_dir=str(tmp_path), max_steps=1, per_device_train_batch_size=4, num_generations=2,
                      max_completion_length=4, use_cpu=True, bf16=False, report_to="none", save_strategy="no",
                      seed=0, remove_unused_columns=False, beta=0.04, loss_type="dr_grpo", use_vllm=True,
                      vllm_mode="colocate", vllm_importance_sampling_correction=False)  # fmt: skip
    rows = Dataset.from_list([{"prompt": [{"role": "user", "content": "x"}]}])
    t = rl_trl_adapter.FastGRPO(model=model, reward_funcs=lambda completions, **kw: [0.0] * len(completions), args=args,
                      train_dataset=rows, processing_class=proc, peft_config=configs.lora("qwen3-vl-2b"))
    t.host_meta = True
    if encoder_once:
        assert t.use_encoder_once()
    return t


def test_trls_full_generate_and_score_runs_under_the_guard(monkeypatch, tmp_path):
    """U-EN7 second half and U-EN1: TRL's own _generate_and_score_completions on the tiny model with the guard on does
    not raise; the batch holds no pixel_values, a row-aligned host_meta of length B and one grid row per row; the
    tower ran once (the encode). The reference logprobs equal the pixel path's within float noise."""
    import PIL.Image
    from transformers import AutoProcessor

    try:
        proc = AutoProcessor.from_pretrained("Qwen/Qwen3-VL-2B-Instruct", local_files_only=True)
    except OSError:
        pytest.skip("the Qwen3-VL processor is not in the local HF cache")
    imgs = [PIL.Image.fromarray(np.random.default_rng(i).integers(0, 255, (64, 64, 3), dtype=np.uint8))
            for i in range(2)]  # fmt: skip
    rows = [{"prompt": [{"role": "user", "content": "outline"}], "image": imgs[i // 2]} for i in range(4)]
    out = {}
    for on in (False, True):
        t = tiny_vl_trainer(monkeypatch, tmp_path, proc, on)
        seen = []
        if on:  # the tower behind the guard
            t.encoder.tower = lambda *a, _t=t.encoder.tower, _s=seen, **k: _s.append(1) or _t(*a, **k)
        with torch.no_grad():
            out[on] = t._generate_and_score_completions(rows)
        assert calls(seen) == on
    b = out[True]
    assert "pixel_values" not in b and "pixel_values" in out[False]
    assert len(b["host_meta"]) == 4 and [r.image_id for r in b["host_meta"]] == [0, 0, 1, 1]
    assert b["image_grid_thw"].shape == (4, 3)
    torch.testing.assert_close(b["ref_per_token_logps"], out[False]["ref_per_token_logps"], atol=1e-5, rtol=1e-5)


def test_train_rl_and_g_probe_carry_encoder_once(capsys):
    from pathlib import Path

    from experiments.phase2 import train_rl

    with pytest.raises(SystemExit):
        train_rl.main(["--help"])
    assert "--encoder-once" in capsys.readouterr().out
    src = (Path(__file__).resolve().parents[2] / "src/objdel/rl/build.py").read_text()
    assert "flags.encoder_once" in src and "trainer.use_encoder_once()" in src


def test_log_reads_the_peak_memory_on_cuda(monkeypatch):
    """G-MEM (C3 step 4): perf/max_memory_gb after the update, with and without --encoder-once (CUDA only)."""
    monkeypatch.setattr(GRPOTrainer, "log", lambda self, logs, start_time=None: None)
    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setattr(torch.cuda, "max_memory_allocated", lambda: 3 * 2**30)
    t = trainer(rl_trl_adapter.FastGRPO)
    t.log({})
    assert t._metrics["train"]["perf/max_memory_gb"] == [3.0]
