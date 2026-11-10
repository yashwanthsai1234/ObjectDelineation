"""objdel.train.fast: per-chunk trimming of the GRPO forwards and the per-step vision cache (perf C1, C6, RL-6).

The reference for every equality test is TRL 1.17's own `_get_per_token_logps_and_entropies` on the untrimmed batch,
run on a tiny random Qwen3-VL (CPU, fp32) with TRL's fused LM head.
"""

import collections
import contextlib
import functools
import json
import math
import time

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from datasets import Dataset  # noqa: E402
from transformers import Qwen3VLConfig, Qwen3VLForConditionalGeneration  # noqa: E402
from trl import GRPOTrainer  # noqa: E402
from trl.trainer.utils import add_fused_lm_head  # noqa: E402

from objdel.train import fast  # noqa: E402

IMG, START, END = 100, 102, 103
PATCH = 3 * 2 * 16 * 16  # channels x temporal patch x patch x patch


def tiny_model(vocab_size=128, img=IMG, video=101, start=START, end=END):
    torch.manual_seed(0)
    cfg = Qwen3VLConfig(
        text_config=dict(vocab_size=vocab_size, hidden_size=32, intermediate_size=64, num_hidden_layers=3,
                         num_attention_heads=2, num_key_value_heads=1, head_dim=16,
                         rope_scaling={"rope_type": "default", "mrope_section": [2, 3, 3], "mrope_interleaved": True}),
        vision_config=dict(depth=2, hidden_size=32, intermediate_size=64, num_heads=2, out_hidden_size=32,
                           patch_size=16, spatial_merge_size=2, temporal_patch_size=2, deepstack_visual_indexes=[0]),
        image_token_id=img, video_token_id=video, vision_start_token_id=start, vision_end_token_id=end)  # fmt: skip
    m = Qwen3VLForConditionalGeneration(cfg).eval()
    add_fused_lm_head(m, outputs=("log_probs", "entropy"))
    return m


class _Accel:
    autocast = staticmethod(contextlib.nullcontext)
    is_main_process = True


def trainer(cls, **attrs):
    """A GRPOTrainer subclass instance with only what the forward and log need (no model load, no vLLM); its state
    is a Trainer's at step 1 (Trainer.__init__ sets one), so log on a CUDA build takes the D7 engine/* keys."""
    from collections import defaultdict
    from types import SimpleNamespace

    t = object.__new__(cls)
    t.accelerator, t._metrics = _Accel(), {"train": defaultdict(list), "eval": defaultdict(list)}
    t.state = SimpleNamespace(global_step=1, log_history=[])
    t.model = type("M", (), {"training": True})()
    for k, v in attrs.items():
        setattr(t, k, v)
    return t


class Plain(GRPOTrainer):
    pass


class Fast(fast.FastForward, GRPOTrainer):
    pass


def batch(pads, lengths, P=12, T=10, images=None, seed=1):
    """Rows with `pads` left-pad tokens, one 4x4-patch image (4 merged tokens) per row, completions of `lengths`."""
    g = torch.Generator().manual_seed(seed)
    n = len(pads)
    ids = torch.randint(0, 99, (n, P + T), generator=g)
    mask = torch.ones_like(ids)
    for r, (p, L) in enumerate(zip(pads, lengths, strict=True)):
        mask[r, :p] = 0
        ids[r, p], ids[r, p + 1:p + 5], ids[r, p + 5] = START, IMG, END
        mask[r, P + L:] = 0
    img = torch.randn(n if images is None else max(images) + 1, 16, PATCH, generator=g)
    pv = img[list(range(n)) if images is None else images].reshape(-1, PATCH)
    kw = {"pixel_values": pv, "image_grid_thw": torch.tensor([[1, 4, 4]] * n), "num_images": [1] * n,
          "mm_token_type_ids": (ids == IMG).long()}  # fmt: skip
    return ids, mask, T, kw


def test_window_cuts_to_the_longest_real_answer_rounded_to_a_bucket():
    # prompt width 12 (left pads 3 and 5 -> 9 real columns), completion width 10, longest real answer 4
    mask = torch.ones(2, 22, dtype=torch.long)
    mask[0, :3], mask[1, :5] = 0, 0
    mask[0, 12 + 4:], mask[1, 12 + 2:] = 0, 0
    assert fast.window(mask, 10, buckets=(8, 16, 20)) == (3, 7)  # need 9 + 4 = 13 -> bucket 16: keep 16 - 9 = 7
    assert fast.window(mask, 10, buckets=(13,)) == (3, 4)  # an exact fit keeps the 3 left pads, no more
    assert fast.window(mask, 10, buckets=(64,)) == (0, 10)  # a bucket wider than the batch: nothing to cut
    mask[:, 12:] = 0  # no real answer token: keep 1 column (TRL's labels[:, :-keep] needs keep >= 1)
    assert fast.window(mask, 10, buckets=(9, 10)) == (3, 1)


@pytest.mark.parametrize("batch_size", [None, 1, 2])
def test_trimmed_logps_equal_trls_untrimmed(batch_size):
    m = tiny_model()
    ids, mask, T, kw = batch(pads=[2, 0, 5], lengths=[3, 7, 1])
    ref_lp, ref_ent, _ = trainer(Plain)._get_per_token_logps_and_entropies(
        m, ids, mask, T, batch_size=batch_size, compute_entropy=True, **kw)
    lp, ent, _ = trainer(Fast, buckets=(4, 8, 12, 16))._get_per_token_logps_and_entropies(
        m, ids, mask, T, batch_size=batch_size, compute_entropy=True, **kw)
    real = mask[:, -T:].bool()
    assert lp.shape == ref_lp.shape == (3, T)
    torch.testing.assert_close(lp[real], ref_lp[real], atol=1e-5, rtol=1e-5)
    torch.testing.assert_close(ent[real], ref_ent[real], atol=1e-5, rtol=1e-5)
    assert (lp[~real] == 0).all()  # cut positions come back as 0, as TRL's own masked positions


def counted_tower(m):
    """Counts the images the real vision tower encodes (grid rows per call)."""
    seen, orig = [], m.model.visual.forward

    def forward(pixel_values, grid_thw=None, **kw):
        seen.append(len(grid_thw))
        return orig(pixel_values, grid_thw=grid_thw, **kw)

    m.model.visual.forward = forward
    return seen


def test_vision_cache_gives_trls_logps_and_encodes_each_image_once():
    m = tiny_model()
    ids, mask, T, kw = batch(pads=[2, 0, 5, 1, 0], lengths=[3, 7, 1, 2, 5], images=[0, 0, 1, 1, 0])  # 2 images
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
    seen = counted_tower(m)
    cache = fast.VisionCache.install(m)
    t = trainer(Fast, buckets=(4, 8, 12, 16))
    for _ in range(3):  # old-logprob, reference and update passes of one generation batch
        lp, _, _ = t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
        torch.testing.assert_close(lp, ref, atol=1e-5, rtol=1e-5)
    assert sum(seen) == 2  # 5 rows x 3 passes = 15 image uses, 2 distinct images encoded
    cache.clear()  # a new generation batch: the next pass encodes again
    t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
    assert sum(seen) == 4


def test_rope_cache_gives_trls_logps_and_computes_each_row_once():
    """C5: Qwen3-VL's 3-D rope positions (get_rope_index: per-row host reads) are computed once per row per generation
    batch at the full width and sliced per chunk, not once per chunk in each of the old, reference and update passes.
    The update's micro-batches are other row subsets (TRL shuffles): rows are found by content."""
    m = tiny_model()
    ids, mask, T, kw = batch(pads=[2, 0, 5, 1, 0], lengths=[3, 7, 1, 2, 5])
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
    rows, orig = [], m.model.get_rope_index
    m.model.get_rope_index = lambda input_ids, **k: rows.append(len(input_ids)) or orig(input_ids, **k)
    t = trainer(Fast, buckets=(4, 8, 12, 16), rope=fast.RopeCache.install(m))
    for _ in range(2):  # old-logprob and reference passes over the whole batch
        lp, _, _ = t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
        torch.testing.assert_close(lp, ref, atol=1e-5, rtol=1e-5)
    perm = [3, 0, 4]  # an update micro-batch: a shuffled subset
    sub = {"pixel_values": kw["pixel_values"].reshape(5, 16, PATCH)[perm].reshape(-1, PATCH), "num_images": [1] * 3,
           "image_grid_thw": kw["image_grid_thw"][perm], "mm_token_type_ids": kw["mm_token_type_ids"][perm]}
    lp, _, _ = t._get_per_token_logps_and_entropies(m, ids[perm], mask[perm], T, **sub)
    torch.testing.assert_close(lp, ref[perm], atol=1e-5, rtol=1e-5)
    assert sum(rows) == 5  # 5 rows x 3 passes, each row's positions computed once
    t.rope.clear()  # a new generation batch
    t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
    assert sum(rows) == 10 and t.rope.pending is None


def test_rope_cache_keys_rows_by_their_exact_tokens_and_image_grids():
    """The prompt text is fixed per dataset, so a row's ids carry only its pad count t*h/2*w/2: grids (1, 4, 8) and
    (1, 8, 4) give the same ids but other image positions. Rows are keyed exactly, ids, mask and grids."""

    class Base(torch.nn.Module):
        calls = []

        def get_rope_index(self, input_ids, image_grid_thw, **k):
            self.calls.append(image_grid_thw.tolist())
            pos = image_grid_thw[:, 1:2].float().expand(len(input_ids), input_ids.size(1))  # h as the position
            return pos.unsqueeze(0).expand(3, -1, -1), None

    cache = fast.RopeCache(Base())
    ids, mask, mm = torch.ones(3, 4, dtype=torch.long), torch.ones(3, 4, dtype=torch.long), torch.zeros(3, 4)
    grid = torch.tensor([[1, 4, 8], [1, 8, 4], [1, 4, 8]])
    pos = cache.positions(ids, mask, mm, grid, [1, 1, 1])
    assert Base.calls == [[[1, 4, 8], [1, 8, 4]]]  # row 2 repeats row 0 exactly: computed once
    assert pos[0, :, 0].tolist() == [4.0, 8.0, 4.0]



def test_rope_positions_are_computed_on_host_copies_and_returned_on_the_rows_device():
    """get_rope_index reads the host several times per row (~9 syncs on CUDA tensors); on the host copies keys()
    already makes, the same function makes none, and one copy returns the positions to the rows' device."""
    seen = []

    class Base(torch.nn.Module):
        def get_rope_index(self, input_ids, image_grid_thw, attention_mask, mm_token_type_ids):
            seen.append({t.device.type for t in (input_ids, image_grid_thw, attention_mask, mm_token_type_ids)})
            return torch.zeros(3, *input_ids.shape), None

    ids = torch.ones(2, 4, dtype=torch.long)
    pos = fast.RopeCache(Base()).positions(ids, torch.ones_like(ids), torch.zeros_like(ids),
                                           torch.tensor([[1, 4, 8], [1, 8, 4]]), [1, 1])  # fmt: skip
    assert seen == [{"cpu"}] and pos.device == ids.device and pos.shape == (3, 2, 4)


def test_lora_leaves_the_vision_tower_frozen():
    """The cache is valid only because no trainable weight sits in the tower or the projector (configs.lora, D6)."""
    from peft import get_peft_model

    from objdel.train import configs

    m = get_peft_model(tiny_model(), configs.lora("qwen3-vl-2b"))
    trainable = [n for n, p in m.named_parameters() if p.requires_grad]
    assert trainable and not [n for n in trainable if ".visual." in n]


def test_fast_grpo_clears_the_vision_cache_per_generation_batch(monkeypatch):
    m = tiny_model()
    t = trainer(fast.FastGRPO, vision=fast.VisionCache.install(m))
    t.vision.store["old image"] = None
    monkeypatch.setattr(GRPOTrainer, "_generate_and_score_completions", lambda self, inputs: {"done": inputs})
    assert t._generate_and_score_completions(["prompts"]) == {"done": ["prompts"]}
    assert t.vision.store == {}


def test_fast_grpo_trains_on_open_once_rows_from_the_main_process(monkeypatch):
    from objdel.train.dataset import OpenOnce

    seen = []
    monkeypatch.setattr(GRPOTrainer, "get_train_dataloader", lambda self: seen.append(self.train_dataset) or "dl")
    rows = Dataset.from_list([{"prompt": "q", "image": None, "key": "k"}])
    t = trainer(fast.FastGRPO, train_dataset=rows, args=type("A", (), {"dataloader_num_workers": 0})())
    assert fast.FastGRPO.get_train_dataloader(t) == "dl"
    assert isinstance(seen[0], OpenOnce) and t.train_dataset is rows  # wrapped only for the loader
    t.args.dataloader_num_workers = 2  # a worker would keep its own last stack: the 16 fetches spread over workers
    with pytest.raises(ValueError, match="dataloader_num_workers"):
        fast.FastGRPO.get_train_dataloader(t)


def test_fast_grpo_skips_flop_counting():
    assert trainer(fast.FastGRPO).floating_point_ops({"input_ids": torch.ones(2, 3)}) == 0


def test_compiled_text_layers_give_eager_logps():
    """Regional compilation (Stage F arm): the text decoder layers compile in place; names and outputs unchanged."""
    m = tiny_model()
    ids, mask, T, kw = batch(pads=[2, 0], lengths=[3, 7])
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, **kw)
    names = [n for n, _ in m.named_parameters()]
    assert fast.compile_text_layers(m, "default") == 3
    lp, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, **kw)
    torch.testing.assert_close(lp, ref, atol=1e-4, rtol=1e-4)
    assert [n for n, _ in m.named_parameters()] == names  # the vLLM weight sync reads these names
    import torch._dynamo.config as dc  # 8 widths x grad on/off x adapter on/off (ref pass): past the default 8, eager
    assert dc.recompile_limit >= 4 * len(fast.BUCKETS) and dc.accumulated_recompile_limit >= 1024


def test_vllm_prompts_carry_one_image_pad_per_image():
    """TRL sends processor-expanded prompt ids plus the raw image (vllm_generation.py 723-727) and trains on its own
    ids, not vLLM's (grpo_trainer.py 1844). One pad per image is the input vLLM's processor expands itself."""
    rows = [{"prompt_token_ids": [1, START, IMG, IMG, IMG, END, 7, IMG, IMG, 8], "multi_modal_data": {"image": "a"}},
            {"prompt_token_ids": [1, 2, 3]}]  # fmt: skip
    assert fast.one_pad_per_image(rows, IMG) == [
        {"prompt_token_ids": [1, START, IMG, END, 7, IMG, 8], "multi_modal_data": {"image": "a"}},
        {"prompt_token_ids": [1, 2, 3]}]  # fmt: skip


def vllm_trainer(expand):
    """A FastGRPO whose vLLM engine records its prompts and expands each pad to `expand` pads, as vLLM's processor."""
    sent = []

    class LLM:
        def generate(self, prompts, **kw):
            sent.append(prompts)
            ids = [[x for t in r["prompt_token_ids"] for x in ([t] * expand if t == IMG else [t])] for r in prompts]
            return [type("O", (), {"prompt_token_ids": i})() for i in ids]

    t = trainer(fast.FastGRPO, vision=None, vllm_generation=type("G", (), {"llm": LLM()})(), use_vllm=True,
                model=type("M", (), {"config": type("C", (), {"image_token_id": IMG})()})())  # fmt: skip
    return t, sent


def test_fast_grpo_sends_vllm_one_pad_per_image():
    t, sent = vllm_trainer(expand=2)
    fast.FastGRPO.patch_vllm_prompts(t)
    fast.FastGRPO.patch_vllm_prompts(t)  # idempotent
    rows = [{"prompt_token_ids": [START, IMG, IMG, END], "multi_modal_data": {"image": "a"}}]
    out = t.vllm_generation.llm.generate(rows, sampling_params=None)
    assert [o.prompt_token_ids for o in out] == [[START, IMG, IMG, END]]
    assert sent == [[{"prompt_token_ids": [START, IMG, END], "multi_modal_data": {"image": "a"}}]]


def test_fast_grpo_fails_loudly_when_vllm_expands_the_image_differently():
    """vLLM then samples from other image features than the trainer scores: off-policy without an error."""
    t, _ = vllm_trainer(expand=3)
    fast.FastGRPO.patch_vllm_prompts(t)
    rows = [{"prompt_token_ids": [START, IMG, IMG, END], "multi_modal_data": {"image": "a"}}]
    with pytest.raises(RuntimeError, match="vLLM expanded the image"):
        t.vllm_generation.llm.generate(rows, sampling_params=None)


def test_vllm_generate_checked_is_the_one_prompt_path():
    """B2.2 F step 4b.3: one pad per image sent, the keywords passed through, a different expansion raises."""
    sent = []

    def generate(rows, expand=2, **kw):
        sent.append((rows, kw))
        ids = [[x for t in r["prompt_token_ids"] for x in ([t] * expand if t == IMG else [t])] for r in rows]
        return [type("O", (), {"prompt_token_ids": i})() for i in ids]

    rows = [{"prompt_token_ids": [START, IMG, IMG, END], "multi_modal_data": {"image": "a"}}]
    fast.vllm_generate_checked(generate, rows, IMG, lora_request=None)
    assert sent == [([{"prompt_token_ids": [START, IMG, END], "multi_modal_data": {"image": "a"}}],
                     {"lora_request": None})]
    with pytest.raises(RuntimeError, match="vLLM expanded the image to 3 pads"):
        fast.vllm_generate_checked(generate, rows, IMG, expand=3)


def test_fast_grpo_prefetch_waits_for_a_generation_batch(monkeypatch):
    """g_probe's step-0 sampling calls _generate outside a generation batch: no rows, nothing to prefetch."""
    f = type("F", (), {"__call__": lambda self, **kw: None})()
    f.trl_reward = type("R", (), {"prefetch": lambda self, *a: pytest.fail("prefetch outside a batch")})()
    t = object.__new__(fast.FastGRPO)  # rows as __init__ leaves them
    t.accelerator, t.model, t.use_vllm, t.reward_funcs = _Accel(), type("M", (), {"training": True})(), False, [f]
    t._metrics, t.rows = {"train": __import__("collections").defaultdict(list)}, fast.FastGRPO.rows
    out = ("pids", [[1]], None, ["a"], None, None, {}, None, [])
    monkeypatch.setattr(GRPOTrainer, "_generate", lambda self, prompts: out)
    assert t._generate(["p"]) == out


def test_fast_grpo_logs_the_phase_times(monkeypatch):
    """time/generate_s (vLLM), time/rewards_s (checker + scoring), time/batch_s (generate + old/ref logprobs +
    rewards); TRL's step_time adds the update. They reach log_history through TRL's metric means."""
    t = trainer(fast.FastGRPO, vision=None, use_vllm=False, reward_funcs=[])
    out = ("pids", [[1]], None, ["a"], None, None, {}, None, [])
    monkeypatch.setattr(GRPOTrainer, "_generate", lambda self, prompts: out)
    monkeypatch.setattr(GRPOTrainer, "_calculate_rewards", lambda self, *a: "rewards")

    def batch(self, inputs):
        assert self._generate(["p"]) == out and self._calculate_rewards(1, 2, 3, 4) == "rewards"
        return {}

    monkeypatch.setattr(GRPOTrainer, "_generate_and_score_completions", batch)
    t._generate_and_score_completions([{"prompt": []}])
    assert sorted(t._metrics["train"]) == ["time/batch_s", "time/generate_s", "time/rewards_s"]  # fmt: skip
    assert all(len(v) == 1 and v[0] >= 0 for v in t._metrics["train"].values())


def test_fast_grpo_starts_scoring_when_vllm_returns(monkeypatch):
    """The answers go to TRLReward.prefetch as soon as _generate returns, with the dataset columns of the same rows,
    so scoring runs while TRL computes the old and reference logprobs on the GPU."""
    got = []
    reward = type("R", (), {"prefetch": lambda self, *a: got.append(a)})()
    f = type("F", (), {"__call__": lambda self, **kw: None})()
    f.trl_reward = reward
    t = trainer(fast.FastGRPO, vision=None, use_vllm=False, reward_funcs=[f])
    rows = [{"prompt": [], "key": k, "W": 50, "H": 40, "names": ["nodule"]} for k in ("a", "a", "b", "b")]
    answers = [[{"role": "assistant", "content": c}] for c in "wxyz"]
    ids = [[1, 7], [2, 7], [3], [4, 7]]
    out = ("pids", ids, None, answers, None, None, {}, None, [])
    monkeypatch.setattr(GRPOTrainer, "_generate", lambda self, prompts: out)
    monkeypatch.setattr(GRPOTrainer, "_generate_and_score_completions",
                        lambda self, inputs: self._generate([r["prompt"] for r in inputs]) and {})
    t._generate_and_score_completions(rows)
    assert got == [(answers, ids, ["a", "a", "b", "b"], [50] * 4, [40] * 4, [["nodule"]] * 4)]
    t.reward_funcs = [lambda **kw: None]  # no TRLReward among the reward functions: nothing to prefetch
    t._generate_and_score_completions(rows)
    assert len(got) == 1


def test_is_ratio_follows_trls_modes_on_a_worked_example():
    """TRL 1.17's vLLM importance-sampling ratio (grpo_trainer.py 2776-2824): token gaps 0, 0.5 and a masked one;
    the sequence gap is 0.5; exp(0.5) = 1.6487."""
    lp, smp, m = torch.tensor([[0.0, -1.0, -2.0]]), torch.tensor([[0.0, -1.5, -1.0]]), torch.tensor([[1.0, 1.0, 0.0]])
    e = 1.6487212707
    cases = {("token_truncate", None, 1.5): [[1.0, 1.5, 1.0]], ("token_mask", None, 1.5): [[1.0, 0.0, 1.0]],
             ("sequence_mask", None, 3.0): [[e]], ("sequence_mask", None, 1.5): [[0.0]],
             ("sequence_truncate", None, 1.5): [[1.5]], ("sequence_mask", 2.0, None): [[0.0]]}  # fmt: skip
    for (mode, lo, hi), want in cases.items():
        torch.testing.assert_close(fast.is_ratio(lp, smp, m, mode, lo, hi), torch.tensor(want), msg=mode)
    nan = torch.tensor([[0.0, float("nan"), -1.0]])  # a token vLLM could not score: no correction
    torch.testing.assert_close(fast.is_ratio(lp, nan, m, "token_mask", None, 3.0), torch.tensor([[1.0, 1.0, 1.0]]))


def c10_trainer(steps_per_generation=16, num_iterations=1):
    args = type("A", (), {"gradient_accumulation_steps": 16, "steps_per_generation": steps_per_generation})()
    return trainer(fast.FastGRPO, vision=None, use_vllm=True, vllm_importance_sampling_correction=True, args=args,
                   num_iterations=num_iterations, loss_type="dr_grpo", vllm_importance_sampling_mode="sequence_mask",
                   vllm_importance_sampling_clip_min=None, vllm_importance_sampling_clip_max=3.0,
                   off_policy_mask_threshold=None)  # fmt: skip


def test_fast_grpo_skips_the_old_logprob_pass_when_generation_and_update_align(monkeypatch):
    """C10: one generation per optimizer step and one iteration: the weights do not change between the old-logprob
    pass and the update, so TRL's old pass (run only for the vLLM IS ratio) is skipped."""
    seen = []
    monkeypatch.setattr(GRPOTrainer, "_generate_and_score_completions",
                        lambda self, inputs: seen.append(self.vllm_importance_sampling_correction) or {})
    t = c10_trainer()
    t._generate_and_score_completions([{"prompt": []}])
    assert seen == [False] and t.vllm_importance_sampling_correction is True  # restored for the loss
    t = c10_trainer(num_iterations=2)  # the weights move between iterations: TRL's old pass is needed
    t._generate_and_score_completions([{"prompt": []}])
    assert seen == [False, True]


def test_update_forward_gives_trl_its_importance_ratio(monkeypatch, tmp_path):
    """The update micro-batch's own logprobs (detached) stand in for the old ones: the ratio is set before TRL's loss
    multiplies by inputs["importance_sampling_ratio"]. TRL's sampling metrics cover the whole generation batch (min,
    max and token-weighted means over every micro-batch), kept on the device and read once, at log time."""
    lp = torch.tensor([[0.0, -1.0, -2.0], [-0.5, -0.5, -0.5]], requires_grad=True)
    monkeypatch.setattr(fast.FastForward, "_get_per_token_logps_and_entropies", lambda self, *a, **k: (lp, None, None))
    logged = []
    monkeypatch.setattr(GRPOTrainer, "log", lambda self, logs, start=None: logged.append({**self._metrics["train"]}))
    t = c10_trainer()
    t.args.output_dir = str(tmp_path)  # log appends the step's seq_gaps rows there
    a = {"sampling_per_token_logps": torch.tensor([[0.0, -1.5, -1.0], [-0.5, -0.5, -0.5]]),
         "completion_mask": torch.tensor([[1, 1, 0], [1, 1, 1]]),
         "advantages": torch.tensor([1.0, -1.0])}  # token gaps 0, 0.5 | 0, 0, 0: ratios e^0.5, 1
    b = {"sampling_per_token_logps": torch.tensor([[0.0, -1.0, -2.0], [-0.5, -0.5, float("nan")]]),
         "completion_mask": torch.tensor([[1, 0, 0], [1, 1, 1]]),
         "advantages": torch.tensor([1.0, -1.0])}  # gaps 0 | 0, 0, NaN (not counted): ratios 1, 1
    for inputs in (a, b):  # two micro-batches of one generation batch
        t.loss_inputs = inputs
        assert t._get_per_token_logps_and_entropies(None, None, None, 3)[0] is lp
    torch.testing.assert_close(a["importance_sampling_ratio"], torch.tensor([[1.6487212707], [1.0]]))
    t.loss_inputs = None  # outside the loss (old and reference passes): nothing set
    t._get_per_token_logps_and_entropies(None, None, None, 3)
    assert not t._metrics["train"]  # nothing read from the device before the log
    t.log({})
    m = logged[0]
    assert m["sampling/sampling_logp_difference/max"] == [pytest.approx(0.5)]
    assert m["sampling/sampling_logp_difference/mean"] == [pytest.approx(0.5 / 8)]  # 8 counted tokens
    assert m["sampling/importance_sampling_ratio/min"] == [1.0]
    assert m["sampling/importance_sampling_ratio/max"] == [pytest.approx(1.6487212707)]
    assert m["sampling/importance_sampling_ratio/mean"] == [pytest.approx((1.6487212707 + 3) / 4)]  # 4 sequences
    t.log({})  # the stats were reset: nothing more is appended (TRL clears _metrics after its own log)
    assert logged[1]["sampling/importance_sampling_ratio/min"] == [1.0]


def test_by_length_orders_rows_longest_first_and_moves_every_field_with_them():
    """TRL shuffles the generation batch before splitting it into update micro-batches, so each chunk is as wide as
    its longest row; sorted, rows of like length share a chunk and FastForward trims each to its own width."""
    d = {"prompt_mask": torch.tensor([[0, 1], [1, 1], [1, 1], [0, 1]]),
         "completion_mask": torch.tensor([[1, 0, 0], [1, 1, 1], [1, 0, 0], [1, 1, 0]]),
         "pixel_values": ["p0", "p1", "p2", "p3"], "num_items_in_batch": torch.tensor(9), "tool_mask": None}
    out = fast.by_length(d)  # real lengths 2, 5, 3, 3: ties keep their order
    assert out["pixel_values"] == ["p1", "p2", "p3", "p0"]
    assert out["completion_mask"].tolist() == [[1, 1, 1], [1, 0, 0], [1, 1, 0], [1, 0, 0]]
    assert out["num_items_in_batch"] == 9 and out["tool_mask"] is None


def test_update_micro_batches_are_length_sorted_only_when_that_leaves_the_gradient_unchanged(monkeypatch):
    """dr_grpo divides each micro-batch's sum by a constant (micro-batch rows x max_completion_length, then the
    accumulation steps: grpo_trainer.py 3285-3289), so when every micro-batch of a generation batch lands in one
    optimizer step (gradient_accumulation_steps a multiple of steps_per_generation), the partition does not change
    the update. Otherwise TRL's shuffle stays."""
    import trl.trainer.grpo_trainer as G

    seen = []
    monkeypatch.setattr(GRPOTrainer, "_prepare_inputs", lambda self, batch: seen.append(G.shuffle_sequence_dict))
    shuffle = G.shuffle_sequence_dict
    t = c10_trainer()
    t.model = type("M", (), {"training": True})()
    t._prepare_inputs({})
    assert seen == [fast.by_length] and G.shuffle_sequence_dict is shuffle  # restored after the call
    t.args.steps_per_generation = 32  # micro-batches of one generation batch span two optimizer steps
    t._prepare_inputs({})
    t.args.steps_per_generation, t.loss_type = 16, "bnpo"  # a per-micro-batch token-count normaliser
    t._prepare_inputs({})
    assert seen[1:] == [shuffle, shuffle]


def test_lora_input_cast_is_dropped_only_under_bf16_autocast_where_it_is_bit_identical():
    """PEFT keeps the adapters in fp32 on a bf16 base and casts every LoRA input bf16 -> fp32; under bf16 autocast
    (TRL's bf16 training: Accelerate wraps the forward) the matmul casts it straight back, so the copy changes nothing
    and FastGRPO turns it off. Without autocast the cast is what makes the dtypes match, so it stays."""
    from peft import get_peft_model

    from objdel.train import configs

    m = get_peft_model(tiny_model().to(torch.bfloat16), configs.lora("qwen3-vl-2b"), autocast_adapter_dtype=True)
    assert next(p for n, p in m.named_parameters() if "lora_A" in n).dtype == torch.float32
    ids = torch.randint(0, 100, (2, 9))
    with torch.autocast("cpu", torch.bfloat16):
        ref = m(input_ids=ids).logits
        n = fast.no_lora_input_cast(m)
        assert n > 0 and all(not mod.cast_input_dtype_enabled for _, mod in m.named_modules()
                             if hasattr(mod, "cast_input_dtype_enabled"))  # fmt: skip
        assert torch.equal(m(input_ids=ids).logits, ref)
    t = trainer(fast.FastGRPO, vision=None, rope=None, args=type("A", (), {"bf16": False})(), model=m)
    assert fast.FastGRPO.drop_lora_input_cast(t) is False  # no autocast: nothing changed (and the model left as is)


def _lora_bf16():
    """tiny_model with configs.lora's LoRA (all-linear, vision excluded) in bf16, LoRA B set non-zero so the merge
    moves the weights (a fresh adapter has B = 0)."""
    from peft import get_peft_model

    from objdel.train import configs

    m = get_peft_model(tiny_model(), configs.lora("qwen3-vl-2b")).to(torch.bfloat16)
    g = torch.Generator().manual_seed(1)
    for n, p in m.named_parameters():
        if "lora_B" in n:
            p.data = (torch.randn(p.shape, generator=g) * 0.05).to(p.dtype)
    return m


def _trl_iter(model):
    """TRL 1.17's VLLMGeneration._iter_named_params on one process (in-place merge_adapter / unmerge_adapter)."""
    from types import MethodType, SimpleNamespace

    from trl.generation.vllm_generation import VLLMGeneration

    dist = SimpleNamespace(is_fsdp=False, gather_params=lambda ps: contextlib.nullcontext())
    host = SimpleNamespace(model=model, _dist=dist)
    host._fix_param_name_to_vllm = MethodType(VLLMGeneration._fix_param_name_to_vllm, host)
    return VLLMGeneration._iter_named_params(host)


def test_vllm_sync_pushes_the_same_weights_as_trls_merge():
    m = _lora_bf16()
    ours = {n: t.clone() for n, t in fast.lora_merged_params(m)}
    trl = {n: t.clone() for n, t in _trl_iter(m)}
    assert ours.keys() == trl.keys() and not any("lora_" in n for n in ours)
    assert all(torch.equal(ours[n], trl[n]) for n in ours)


def test_vllm_sync_never_writes_the_base_weights():
    """TRL merges LoRA into the bf16 base weights in place and unmerges: W + d - d != W in bf16 (trl #6688), so the
    frozen base drifts every step. The sync must leave every parameter bit-identical."""
    m = _lora_bf16()
    before = {n: p.detach().clone() for n, p in m.named_parameters()}
    for _ in range(20):
        for _ in fast.lora_merged_params(m):
            pass
    assert all(torch.equal(p, before[n]) for n, p in m.named_parameters())
    for _ in range(20):  # the bug this replaces, so the test above is not vacuous
        for _ in _trl_iter(m):
            pass
    assert any(not torch.equal(p, before[n]) for n, p in m.named_parameters())


def test_fast_grpo_routes_trls_vllm_sync_through_the_no_merge_path():
    from types import SimpleNamespace

    m = _lora_bf16()
    t = object.__new__(fast.FastGRPO)
    t.model, t.accelerator = m, SimpleNamespace(num_processes=1)
    t.vllm_generation = SimpleNamespace(model=m, _dist=SimpleNamespace(is_fsdp=False))
    assert t.sync_without_merge()
    before = {n: p.detach().clone() for n, p in m.named_parameters()}
    pushed = dict(t.vllm_generation._iter_named_params())
    assert pushed.keys() == dict(fast.lora_merged_params(m)).keys()
    assert all(torch.equal(p, before[n]) for n, p in m.named_parameters())
    t.vllm_generation = None  # no vLLM (use_vllm False): nothing to route
    assert not t.sync_without_merge()


def test_uncompiled_chunks_round_widths_to_64_not_to_the_compile_buckets():
    """The 8 BUCKETS bound torch.compile's recompiles; uncompiled they only pad (430 tokens ran 512 wide)."""
    mask = torch.zeros(2, 1000, dtype=torch.long)  # prompt width P = 400, 600 answer columns
    mask[:, 70:400] = 1  # 330 prompt tokens
    mask[0, 400:500] = 1  # the longest answer: 100 tokens
    start, keep = fast.window(mask, 600, fast.FastForward.buckets)
    assert fast.FastForward.buckets == tuple(range(64, 2049, 64))
    assert (400 + keep - start, keep) == (448, 118)  # 330 + 100 = 430 -> 448
    start, keep = fast.window(mask, 600, fast.BUCKETS)
    assert 400 + keep - start == 512


def test_compiled_layers_go_back_to_the_compile_buckets():
    from pathlib import Path

    src = (Path(__file__).resolve().parents[1] / "experiments/phase2/g_probe.py").read_text()
    assert "fast.compile_text_layers(trainer.model, a.compile)" in src and "trainer.buckets = fast.BUCKETS" in src


def test_vllm_sync_wake_and_sleep_are_timed():
    from collections import defaultdict
    from types import SimpleNamespace

    calls = []
    llm = SimpleNamespace(wake_up=lambda tags=None: calls.append(("wake", tags)),
                          sleep=lambda level=1: calls.append(("sleep", level)))
    vg = SimpleNamespace(llm=llm, sync_weights=lambda: calls.append(("sync",)))
    t = object.__new__(fast.FastGRPO)
    t.model, t.vllm_generation = SimpleNamespace(training=True), vg
    t._metrics = {"train": defaultdict(list), "eval": defaultdict(list)}
    t.time_vllm_phases()
    t.time_vllm_phases()  # idempotent: one timer layer
    vg.sync_weights(), llm.wake_up(tags=["kv_cache"]), llm.sleep(level=2)
    assert calls == [("sync",), ("wake", ["kv_cache"]), ("sleep", 2)]
    timed = {k: len(v) for k, v in t._metrics["train"].items()}
    assert timed == {"time/sync_s": 1, "time/wake_kv_s": 1, "time/sleep_s": 1}


def test_wake_time_is_split_by_tags_and_no_wake_time_is_lost(monkeypatch):
    """U-MS3: ["weights"] charges time/wake_weights_s only, ["kv_cache"] time/wake_kv_s only; any other call (no tags,
    both tags) charges time/wake_s, so the sum of the wake keys is the total."""
    from collections import defaultdict
    from types import SimpleNamespace

    clock = iter(range(100))
    monkeypatch.setattr(fast.time, "perf_counter", lambda: next(clock) * 1.0)
    llm = SimpleNamespace(wake_up=lambda tags=None: None, sleep=lambda level=1: None)
    t = object.__new__(fast.FastGRPO)
    t.model, t.vllm_generation = SimpleNamespace(training=True), SimpleNamespace(llm=llm, sync_weights=lambda: None)
    t._metrics = {"train": defaultdict(list), "eval": defaultdict(list)}
    t.time_vllm_phases()
    llm.wake_up(tags=["weights"]), llm.wake_up(["kv_cache"]), llm.wake_up(), llm.wake_up(tags=["weights", "kv_cache"])
    assert {k: len(v) for k, v in t._metrics["train"].items()} == {
        "time/wake_weights_s": 1, "time/wake_kv_s": 1, "time/wake_s": 2}  # fmt: skip


def group_batch(groups, G, lengths, pads=None, P=12, T=10, seed=3):
    """GRPO layout: `groups` prompts x G answers, rows of a group consecutive with the same prompt tokens, pads and
    image; answers of `lengths` (one per row) differ."""
    g = torch.Generator().manual_seed(seed)
    n = groups * G
    ids = torch.randint(0, 99, (n, P + T), generator=g)
    mask = torch.ones_like(ids)
    pads = pads or [0] * groups
    for q in range(groups):
        p, rows = pads[q], slice(q * G, q * G + G)
        prompt = torch.randint(0, 99, (P,), generator=g)
        prompt[p], prompt[p + 1:p + 5], prompt[p + 5] = START, IMG, END
        ids[rows, :P] = prompt
        mask[rows, :p] = 0
    for r, L in enumerate(lengths):
        mask[r, P + L:] = 0
    img = torch.randn(groups, 16, PATCH, generator=g)
    pv = img[[r // G for r in range(n)]].reshape(-1, PATCH)
    kw = {"pixel_values": pv, "image_grid_thw": torch.tensor([[1, 4, 4]] * n), "num_images": [1] * n,
          "mm_token_type_ids": (ids == IMG).long()}  # fmt: skip
    return ids, mask, T, kw


def shared(m, **attrs):
    return trainer(Fast, shared_prefix=True, rope=fast.RopeCache.install(m), **attrs)


@pytest.mark.parametrize("batch_size", [8, 4])
def test_shared_prefix_logps_equal_trls(batch_size):
    """Each group's prompt runs once: [prompt minus its last token] + G x [last prompt token, answer], answers masked
    from each other. Same log-probs and entropies as TRL's one-row-per-answer forward."""
    m = tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[2, 0])
    ref_lp, ref_ent, _ = trainer(Plain)._get_per_token_logps_and_entropies(
        m, ids, mask, T, batch_size=batch_size, compute_entropy=True, **kw)
    lp, ent, _ = shared(m)._get_per_token_logps_and_entropies(
        m, ids, mask, T, batch_size=batch_size, compute_entropy=True, **kw)
    real = mask[:, -T:].bool()
    assert lp.shape == ref_lp.shape == (8, T)
    torch.testing.assert_close(lp[real], ref_lp[real], atol=1e-5, rtol=1e-5)
    torch.testing.assert_close(ent[real], ref_ent[real], atol=1e-5, rtol=1e-5)
    assert (lp[~real] == 0).all()


def test_shared_prefix_gradients_equal_trls():
    m = tiny_model().train()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[1, 3])
    real = mask[:, -T:].float()
    grads = []
    for t in (trainer(Plain), shared(m)):
        m.zero_grad()
        lp, _, _ = t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=8, **kw)
        (lp * real).sum().backward()
        grads.append({n: p.grad.clone() for n, p in m.named_parameters() if p.grad is not None})
    assert grads[0].keys() == grads[1].keys()
    for n in grads[0]:
        torch.testing.assert_close(grads[1][n], grads[0][n], atol=1e-4, rtol=1e-4)


def test_shared_prefix_runs_far_fewer_tokens():
    m = tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10])
    seen, orig = [], m.model.language_model.forward

    def forward(*a, **k):
        seen.append(k["inputs_embeds"].shape[:2])
        return orig(*a, **k)

    m.model.language_model.forward = forward
    t = shared(m)
    t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=8, **kw)
    (rows, width), = seen
    assert rows == 2  # one packed row per group
    assert width == 11 + 4 * 1 + (2 + 6 + 4 + 10)  # prompt minus its last token + 4 x (last prompt token + answer)
    assert rows * width < 8 * (12 + 10) / 2  # 74 tokens, not TRL's 8 rows x 22
    real = 2 * 4 * 12 + (3 + 7 + 1 + 5 + 2 + 6 + 4 + 10)  # every row's prompt + its answer: 134
    assert t._metrics["train"]["shared/token_ratio"] == [(2 * 11 + 8 + 38) / real]


def test_shared_prefix_falls_back_when_rows_do_not_share_prompts():
    m = tiny_model()
    ids, mask, T, kw = batch(pads=[2, 0, 5], lengths=[3, 7, 1])  # three different prompts
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, **kw)
    lp, _, _ = shared(m)._get_per_token_logps_and_entropies(m, ids, mask, T, **kw)
    torch.testing.assert_close(lp, ref, atol=1e-5, rtol=1e-5)


def test_by_group_length_keeps_each_group_together_longest_group_first():
    n = torch.tensor([3, 9, 2, 4, 5, 5, 1, 1, 7, 2, 2, 2])  # 3 groups of 4: longest members 9, 5, 7
    batch_ = {"prompt_mask": torch.ones(12, 1, dtype=torch.long),
              "completion_mask": torch.zeros(12, 9, dtype=torch.long),
              "tag": list(range(12)), "advantages": torch.arange(12.0)}
    for r, L in enumerate(n.tolist()):
        batch_["completion_mask"][r, :L] = 1
    out = fast.by_group_length(batch_, 4)
    assert out["tag"] == [0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7]
    assert out["advantages"].tolist() == [float(i) for i in out["tag"]]


def test_shared_prefix_update_micro_batches_keep_whole_groups(monkeypatch):
    """With shared_prefix the update's order keeps each prompt's G answers together (by_group_length), so a
    micro-batch of whole groups packs; micro-batches that cannot hold whole groups keep by_length."""
    import trl.trainer.grpo_trainer as G

    seen = []
    monkeypatch.setattr(GRPOTrainer, "_prepare_inputs", lambda self, batch: seen.append(G.shuffle_sequence_dict))
    t = c10_trainer()
    t.model, t.shared_prefix, t.num_generations = type("M", (), {"training": True})(), True, 8
    t.args.per_device_train_batch_size = 16
    t._prepare_inputs({})
    assert seen[0].func is fast.by_group_length and seen[0].keywords == {"G": 8}
    t.args.per_device_train_batch_size = 12  # not whole groups of 8
    t._prepare_inputs({})
    assert seen[1] is fast.by_length


def test_shared_prefix_normalizes_over_the_sampled_support_as_trl():
    m = tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[2, 0])
    g = torch.Generator().manual_seed(5)
    sup = torch.randint(0, 99, (8, T, 6), generator=g)
    sup[..., 0] = ids[:, -T:]  # the sampled token is in its support
    sup[..., 5] = -1  # padded support slot
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, sampling_mask=sup, **kw)
    lp, _, _ = shared(m)._get_per_token_logps_and_entropies(m, ids, mask, T, sampling_mask=sup, **kw)
    real = mask[:, -T:].bool()
    torch.testing.assert_close(lp[real], ref[real], atol=1e-5, rtol=1e-5)


def test_train_rl_and_g_probe_turn_on_shared_prefix_by_flag(capsys):
    from pathlib import Path

    from experiments.phase2 import train_rl

    with pytest.raises(SystemExit):
        train_rl.main(["--help"])
    assert "--shared-prefix" in capsys.readouterr().out
    src = (Path(__file__).resolve().parents[1] / "experiments/phase2/g_probe.py").read_text()
    assert 'trainer.shared_prefix = bool(getattr(a, "shared_prefix", False))' in src


@pytest.mark.parametrize("impl", ["eager", "sdpa"])
def test_shared_prefix_mask_works_under_both_attention_kinds(impl):
    m = tiny_model()
    m.set_attn_implementation(impl)
    assert m.config._attn_implementation == impl
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[2, 0])
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, **kw)
    lp, _, _ = shared(m)._get_per_token_logps_and_entropies(m, ids, mask, T, **kw)
    real = mask[:, -T:].bool()
    torch.testing.assert_close(lp[real], ref[real], atol=1e-5, rtol=1e-5)


def test_packing_is_chosen_only_when_it_costs_less():
    """Packed attention is dense over the whole packed row (prompt + 8 answers): with long answers its square
    outgrows the 8 prompts it saves. Real sizes: 325-token prompts, 16 rows, linear cost ~12,000x an attention pair."""
    short = [127] * 16
    assert fast.packing_wins(325, short, groups=2)
    assert not fast.packing_wins(325, [1536] * 16, groups=2)  # answers at the cap: per-row is cheaper
    assert not fast.packing_wins(325, short, groups=16)  # nothing shared


def test_long_answers_fall_back_to_the_per_row_forward(monkeypatch):
    m = tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10])
    monkeypatch.setattr(fast, "packing_wins", lambda *a, **k: False)
    t = shared(m)
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, **kw)
    lp, _, _ = t._get_per_token_logps_and_entropies(m, ids, mask, T, **kw)
    torch.testing.assert_close(lp, ref, atol=1e-5, rtol=1e-5)
    assert not t._metrics["train"]["shared/token_ratio"]


def test_a_model_without_the_rope_hook_takes_the_per_row_forward():
    """A separate reference model (full fine-tuning, beta > 0) has no RopeCache hook: it must not get a packed row."""
    m, other = tiny_model(), tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10])
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(other, ids, mask, T, **kw)
    t = shared(m)
    lp, _, _ = t._get_per_token_logps_and_entropies(other, ids, mask, T, **kw)
    torch.testing.assert_close(lp, ref, atol=1e-5, rtol=1e-5)
    assert not t._metrics["train"]["shared/token_ratio"]


def test_update_order_keeps_trls_shuffle_under_a_per_micro_batch_entropy_cut(monkeypatch):
    import trl.trainer.grpo_trainer as G

    seen = []
    monkeypatch.setattr(GRPOTrainer, "_prepare_inputs", lambda self, batch: seen.append(G.shuffle_sequence_dict))
    shuffle = G.shuffle_sequence_dict
    t = c10_trainer()
    t.model, t.top_entropy_quantile = type("M", (), {"training": True})(), 0.2
    t._prepare_inputs({})
    assert seen == [shuffle]


def test_shared_prefix_is_off_when_the_text_layers_are_compiled():
    from pathlib import Path

    src = (Path(__file__).resolve().parents[1] / "experiments/phase2/g_probe.py").read_text()
    assert 'trainer.shared_prefix = bool(getattr(a, "shared_prefix", False)) and not getattr(a, "compile", None)' in src


# B2.1: what the IS correction loses, and what every other mode would do ------------------------------------------


class FakeVLLM:
    """TRL's VLLMGeneration without an engine (fixture tiny_loop_vllm_fake, B.9.0). Each row's answer and planted
    per-sequence gap c are drawn from a generator seeded by the row's prompt and its place in the group, so an
    eval in between does not change later train rows. Its sampling logprobs are the model's own plus c on every
    token: seq_kl (vLLM minus trainer, mean over the answer) of row i is c_i up to float error. planted: (c, length)
    per generated row."""

    VOCAB = 43

    def __init__(self, model, max_completion_length, **kw):
        self.model, self.max_len, self.planted = model, max_completion_length, []
        noop = lambda *a, **k: None  # noqa: E731
        self.llm = type("LLM", (), {"generate": staticmethod(noop), "wake_up": staticmethod(noop),
                                    "sleep": staticmethod(noop)})()  # fmt: skip

    def sync_weights(self):
        pass

    @torch.no_grad()
    def generate(self, prompts, images, num_generations, profiler=None):
        ids, lps = [], []
        for j, p in enumerate(prompts):
            g = torch.Generator().manual_seed(hash((tuple(p), j % num_generations)) % 2**31)
            n = int(torch.randint(1, self.max_len + 1, (1,), generator=g))
            c = float(torch.rand(1, generator=g) * 0.4 - 0.2)
            a = torch.randint(3, self.VOCAB, (n,), generator=g)
            x = torch.cat([torch.tensor(p), a])[None]
            lp = self.model(input_ids=x).logits[0, len(p) - 1:-1].float().log_softmax(-1)[torch.arange(n), a]
            ids.append(a.tolist())
            lps.append([[v + c] for v in lp.tolist()])
            self.planted.append((c, n))
        return prompts, ids, lps, None, None


def tiny_loop_vllm_fake(monkeypatch, out, steps=3, eval_steps=None):
    """FastGRPO on a 1-layer Llama (CPU) for `steps` steps, logging every step, with TRL's vLLM path on (IS
    correction, sequence_mask cap 3) and FakeVLLM in place of the engine; 8 rows per step (rows_per_step); TRL's
    own train loader (OpenOnce needs an image column).
    eval_steps: a held-out eval every eval_steps steps. Returns the trained trainer."""
    import trl.trainer.grpo_trainer as G
    from tokenizers import Tokenizer, models, pre_tokenizers
    from transformers import LlamaConfig, LlamaForCausalLM, PreTrainedTokenizerFast
    from trl import GRPOConfig

    monkeypatch.setattr(G, "VLLMGeneration", FakeVLLM)
    monkeypatch.setattr(fast.FastGRPO, "get_train_dataloader", GRPOTrainer.get_train_dataloader)  # rows: no image
    monkeypatch.setattr(fast.FastGRPO, "patch_vllm_prompts", lambda self: None)  # text only: no image pads
    vocab = {w: i for i, w in enumerate(["<pad>", "<eos>", "prompt", *map(str, range(FakeVLLM.VOCAB - 3))])}
    tk = Tokenizer(models.WordLevel(vocab, unk_token="<pad>"))
    tk.pre_tokenizer = pre_tokenizers.WhitespaceSplit()
    tok = PreTrainedTokenizerFast(tokenizer_object=tk, pad_token="<pad>", eos_token="<eos>", padding_side="left")
    torch.manual_seed(0)
    model = LlamaForCausalLM(LlamaConfig(vocab_size=len(vocab), hidden_size=16, intermediate_size=32,
                                         num_hidden_layers=1, num_attention_heads=2, num_key_value_heads=2,
                                         pad_token_id=0, eos_token_id=1, bos_token_id=1))  # fmt: skip
    ev = {"eval_strategy": "steps", "eval_steps": eval_steps, "per_device_eval_batch_size": 8} if eval_steps else {}
    args = GRPOConfig(output_dir=str(out), max_steps=steps, per_device_train_batch_size=8, num_generations=8,
                      max_completion_length=4, use_cpu=True, bf16=False, report_to="none", logging_steps=1,
                      save_strategy="no", seed=0, remove_unused_columns=False, learning_rate=1e-2, beta=0.0,
                      loss_type="dr_grpo", use_vllm=True, vllm_mode="colocate",
                      vllm_importance_sampling_correction=True, vllm_importance_sampling_mode="sequence_mask",
                      **ev)  # fmt: skip
    rows = lambda *n: Dataset.from_list([{"prompt": f"prompt {i}"} for i in n])  # noqa: E731
    t = fast.FastGRPO(model=model, reward_funcs=lambda completions, **kw: [len(c) % 3 * 1.0 for c in completions],
                      args=args, train_dataset=rows(1, 2, 3, 4), eval_dataset=rows(7) if eval_steps else None,
                      processing_class=tok)  # fmt: skip
    t.train()
    return t


OLD5 = ("sampling_logp_difference/mean", "sampling_logp_difference/max", "importance_sampling_ratio/min",
        "importance_sampling_ratio/mean", "importance_sampling_ratio/max")


def train_logs(t) -> list[dict]:
    return [h for h in t.state.log_history if "loss" in h]


R01_GOLDEN = {  # captured on this fixture before B2.1 (commit f5e0a05 + A6), CPU fp32
    "sampling_logp_difference/mean": [0.07079750015622094, 0.11064745903015137, 0.09729175269603729],
    "sampling_logp_difference/max": [0.17806434631347656, 0.19784259796142578, 0.19302105903625488],
    "importance_sampling_ratio/min": [0.8013280034065247, 0.632602870464325, 0.4812847077846527],
    "importance_sampling_ratio/mean": [1.092457890510559, 1.2020281553268433, 1.005218505859375],
    "importance_sampling_ratio/max": [1.706070899963379, 2.050137758255005, 1.2129082679748535]}


def test_r01_old_sampling_metrics_unchanged(monkeypatch, tmp_path):
    logs = train_logs(tiny_loop_vllm_fake(monkeypatch, tmp_path))
    assert {k: [h[f"sampling/{k}"] for h in logs] for k in OLD5} == R01_GOLDEN


B21_KEYS = ("is_masked_share", "is_clamped_share", "seq_logdiff_mean", "seq_logdiff_max", "tok_absdiff_gt_0.5_share",
            "masked_len_mean", "kept_len_mean", "modified_token_share", "ess", "cf/seqmask3_share",
            "cf/tok_clamp2_share", "cf/geo_0.99_1.01_share", "cf/opsm_share@0.05", "cf/opsm_share@0.1",
            "cf/opsm_share@0.5", "cf/opsm0.1_bigtok_unmasked_share", "cf/seq_ratio_lt_0.1_share", "cf/ess_C",
            "cf/ess_B", "cf/ess_E@0.05", "cf/ess_E@0.1", "cf/ess_E@0.5", "seq_kl/q50", "seq_kl/q99",
            "seq_kl/q99.5")  # fmt: skip


def rows(*gaps, adv=None, kl_nan=()):
    """One micro-batch from per-row trainer-minus-vLLM token gaps (rows right-padded): logps 0, sampling = -gap;
    kl_nan: (row, token) places where vLLM could not score (sampling NaN)."""
    T = max(map(len, gaps))
    g = torch.zeros(len(gaps), T)
    m = torch.zeros(len(gaps), T, dtype=torch.long)
    for i, x in enumerate(gaps):
        g[i, : len(x)], m[i, : len(x)] = torch.tensor(x, dtype=torch.float32), 1
    s = -g
    for i, j in kl_nan:
        s[i, j] = float("nan")
    a = torch.tensor(adv if adv is not None else [1.0] * len(gaps))
    return torch.zeros_like(g), {"sampling_per_token_logps": s, "completion_mask": m, "advantages": a}


def is_trainer(mode, lo, hi, out, opsm=None):
    from types import SimpleNamespace

    return trainer(fast.FastGRPO, use_vllm=True, vllm_importance_sampling_correction=True, args=SimpleNamespace(
        output_dir=str(out)), vllm_importance_sampling_mode=mode, vllm_importance_sampling_clip_min=lo,
        vllm_importance_sampling_clip_max=hi, off_policy_mask_threshold=opsm)  # fmt: skip


def logged(monkeypatch, out, mode, lo, hi, *micro, opsm=None) -> dict:
    """set_is_ratio on each micro-batch of one generation batch, then log: the train metrics (one value each)."""
    got = []
    monkeypatch.setattr(GRPOTrainer, "log", lambda self, logs, start=None: got.append(dict(self._metrics["train"])))
    t = is_trainer(mode, lo, hi, out, opsm)
    for lp, inputs in micro:
        t.set_is_ratio(lp, inputs)
    t.log({})
    return {k.removeprefix("sampling/"): v[0] for k, v in got[0].items() if k.startswith("sampling/")}


def test_u37_active_mode_stats_on_a_hand_case(monkeypatch, tmp_path):
    """sequence_mask, cap 3: row 1 sums 1.8 (e^1.8 = 6.05 > 3: masked, 3 tokens); the others are kept (0.2 over 2
    tokens, 0 over 1, -2.0 over 4). Tokens with |gap| > 0.5: three 0.6 and one -0.6; -0.5 is not (strict)."""
    m = logged(monkeypatch, tmp_path, "sequence_mask", None, 3.0,
               rows([0.1, 0.1], [0.6, 0.6, 0.6], [0.0], [-0.5, -0.5, -0.6, -0.4]))  # fmt: skip
    assert m["is_masked_share"] == 0.25 and m["is_clamped_share"] == 0.0
    assert m["seq_logdiff_mean"] == pytest.approx(0.0, abs=1e-6) and m["seq_logdiff_max"] == pytest.approx(1.8)
    assert m["tok_absdiff_gt_0.5_share"] == pytest.approx(4 / 10)
    assert m["masked_len_mean"] == 3.0 and m["kept_len_mean"] == pytest.approx(7 / 3)
    m = logged(monkeypatch, tmp_path, "sequence_mask", None, 3.0, rows([0.1], [0.2]))
    assert m["masked_len_mean"] == 0.0 and m["is_masked_share"] == 0.0  # none masked: 0, not NaN


def test_u37_clamped_counts_only_outside_the_bounds(monkeypatch, tmp_path):
    """token_truncate: e^0.8 = 2.23 > 2 is clamped; a ratio exactly equal to hi is not; with lo None nothing counts
    below (e^-5)."""
    hi = float(torch.exp(torch.tensor(0.3)))  # the float32 ratio of a 0.3 gap, exactly
    m = logged(monkeypatch, tmp_path, "token_truncate", None, hi, rows([0.8, 0.3, -5.0, 0.0]))
    assert m["is_clamped_share"] == 0.25 and m["is_masked_share"] == 0.0
    m = logged(monkeypatch, tmp_path, "token_truncate", 0.1, hi, rows([0.8, 0.3, -5.0, 0.0]))
    assert m["is_clamped_share"] == 0.5  # e^-5 < 0.1 now counts


def test_u37_summary_is_float32_on_the_input_device():
    lp, x = rows([0.1, 0.2], [0.3])
    s = fast.is_summary(lp, x["sampling_per_token_logps"], x["completion_mask"].bool(),
                        fast.is_ratio(lp, x["sampling_per_token_logps"], x["completion_mask"], "token_mask", None, 3.0),
                        x["advantages"], "token_mask", None, 3.0)  # fmt: skip
    assert s.dtype == torch.float32 and s.device == lp.device and s.shape == (len(fast.SUMMARY),)


U54 = rows([0.005] * 800,  # 0: sum 4.0 (e^4 > 3), mean 0.005 inside the geo band, seq_kl -0.005
           [-0.2] * 5,  # 1: seq_kl +0.2, negative advantage: OPSM drops it at 0.05 and 0.1, not at 0.5
           [0.2] * 5,  # 2: mirror, seq_kl -0.2: OPSM keeps it
           [-0.8] * 2,  # 3: seq_kl 0.8 but positive advantage: never OPSM; 2 big tokens
           [-0.15, -0.15],  # 4: first token unscored (NaN): seq_kl (0 + 0.15) / 2 = 0.075, as TRL
           [math.log(0.05)],  # 5: weight 0.05 under sequence_mask (no lower bound): kept, down-weighted; 1 big token
           adv=[-1.0, -1.0, -1.0, 1.0, -1.0, 1.0], kl_nan=[(4, 0)])  # fmt: skip


def test_u54_counterfactual_shares_on_planted_rows(monkeypatch, tmp_path):
    m = logged(monkeypatch, tmp_path, "sequence_mask", None, 3.0, U54)
    assert m["cf/seqmask3_share"] == pytest.approx(1 / 6)  # row 0 only
    assert m["cf/geo_0.99_1.01_share"] == pytest.approx(5 / 6)  # all but row 0 (row 4: e^-0.075 < 0.99)
    assert m["cf/opsm_share@0.05"] == pytest.approx(2 / 6)  # rows 1 and 4
    assert m["cf/opsm_share@0.1"] == pytest.approx(1 / 6)  # row 1; the mirror row 2 is kept
    assert m["cf/opsm_share@0.5"] == 0.0
    assert m["cf/seq_ratio_lt_0.1_share"] == pytest.approx(1 / 6)  # row 5
    assert m["cf/tok_clamp2_share"] == 0.0  # no token gap above ln 2
    assert m["cf/opsm0.1_bigtok_unmasked_share"] == pytest.approx(3 / 815)  # rows 3 and 5 kept; row 1 is not big
    assert m["is_masked_share"] == pytest.approx(1 / 6)  # active sequence_mask drops row 0


def test_u54_seq_kl_quantiles_equal_numpy(monkeypatch, tmp_path):
    m = logged(monkeypatch, tmp_path, "sequence_mask", None, 3.0, U54)
    kl = np.array([-0.005, 0.2, -0.2, 0.8, 0.075, -math.log(0.05)], dtype=np.float32)
    for q, k in ((0.5, "q50"), (0.99, "q99"), (0.995, "q99.5")):
        assert m[f"seq_kl/{k}"] == pytest.approx(float(np.quantile(kl, q)), abs=1e-6)


def test_u54_ess_active_and_counterfactual_at_token_level(monkeypatch, tmp_path):
    """Row 0: gaps (1, 1), seq_kl -1; row 1: gaps (0, -1), seq_kl 0.5; both negative advantage. Token weights:
    active sequence_mask cap 3 and C: (0, 0, a, a) with a = e^-1 (e^2 > 3 dropped); B, token TIS C = 2: (2, 2, 1, a);
    E@0.1: OPSM drops row 1, (2, 2, 0, 0); E@0.5 keeps it (0.5 <= 0.5). ESS = (sum w)^2 / (N sum w^2), N = 4."""
    a = math.exp(-1)
    m = logged(monkeypatch, tmp_path, "sequence_mask", None, 3.0, rows([1.0, 1.0], [0.0, -1.0], adv=[-1.0, -1.0]))
    assert m["ess"] == pytest.approx(0.5) and m["cf/ess_C"] == pytest.approx(0.5)
    assert m["cf/ess_B"] == pytest.approx((5 + a) ** 2 / (4 * (9 + a * a)))
    assert m["cf/ess_E@0.1"] == pytest.approx(16 / (4 * 8))
    assert m["cf/ess_E@0.5"] == pytest.approx(m["cf/ess_B"]) and m["cf/ess_E@0.05"] == pytest.approx(0.5)
    assert logged(monkeypatch, tmp_path, "token_mask", None, 3.0, rows([0.0] * 4))["ess"] == 1.0  # equal weights
    one = logged(monkeypatch, tmp_path, "token_mask", None, 3.0, rows([5.0, 5.0, 0.0, 5.0]))  # one non-zero weight
    assert one["ess"] == pytest.approx(1 / 4)


def test_u54_modified_token_share_counts_a_mask_and_a_2pct_clamp_not_half_a_pct(monkeypatch, tmp_path):
    m = logged(monkeypatch, tmp_path, "token_mask", None, 3.0, rows([math.log(4), 0.0, 0.0, 0.0]))
    assert m["modified_token_share"] == 0.25
    m = logged(monkeypatch, tmp_path, "token_truncate", 1.0, 3.0,
               rows([math.log(0.98), math.log(0.995), 0.0, 0.0]))  # fmt: skip
    assert m["modified_token_share"] == 0.25


def test_u62_counterfactual_opsm_equals_trls_off_policy_mask():
    """Oracle: TRL's own GRPOTrainer.get_off_policy_mask (imported here only), drops = rows it returns 0 for."""
    g = torch.Generator().manual_seed(0)
    names = list(fast.SUMMARY)
    for _ in range(200):
        lp = -torch.rand(6, 9, generator=g) * 3
        sp = lp + torch.randn(6, 9, generator=g) * torch.rand(6, 1, generator=g)
        sp[0, 2] = float("nan")
        mask = (torch.rand(6, 9, generator=g) < 0.8).long()
        mask[1] = 0
        adv = torch.randn(6, generator=g)
        r = fast.is_ratio(lp, sp, mask, "sequence_mask", None, 3.0)
        s = fast.is_summary(lp, sp, mask.bool(), r, adv, "sequence_mask", None, 3.0).tolist()
        for d in fast.OPSM_DELTAS:
            keep = GRPOTrainer.get_off_policy_mask(adv[:, None], lp, sp, mask.float(), d)
            assert s[names.index(f"cf_opsm_{d}")] == float((keep == 0).sum())


def test_u67_sixteen_micro_batches_give_the_quantiles_of_one_batch(monkeypatch, tmp_path):
    g = torch.Generator().manual_seed(1)
    gaps = torch.randn(256, 20, generator=g) * 0.1
    batch = lambda sl: (torch.zeros(len(gaps[sl]), 20), {  # noqa: E731
        "sampling_per_token_logps": -gaps[sl], "completion_mask": torch.ones(len(gaps[sl]), 20),
        "advantages": torch.ones(len(gaps[sl]))})  # fmt: skip
    one = logged(monkeypatch, tmp_path / "a", "sequence_mask", None, 3.0, batch(slice(None)))
    many = logged(monkeypatch, tmp_path / "b", "sequence_mask", None, 3.0,
                  *(batch(slice(16 * i, 16 * i + 16)) for i in range(16)))  # fmt: skip
    for k in ("q50", "q99", "q99.5"):
        assert many[f"seq_kl/{k}"] == pytest.approx(one[f"seq_kl/{k}"], abs=1e-6)


def test_u67_the_buffer_stays_on_the_device_until_log_and_is_cleared(monkeypatch, tmp_path):
    t = is_trainer("sequence_mask", None, 3.0, tmp_path)
    monkeypatch.setattr(GRPOTrainer, "log", lambda self, logs, start=None: None)
    reads = []
    for name in ("item", "cpu", "tolist", "numpy", "__float__", "__bool__"):
        orig = getattr(torch.Tensor, name)
        monkeypatch.setattr(torch.Tensor, name, lambda self, *a, _o=orig, _n=name: reads.append(_n) or _o(self, *a))
    for _ in range(3):
        t.set_is_ratio(*U54)
    assert reads == [] and isinstance(t.is_stats, torch.Tensor) and len(t.seq_gaps) == 3
    t.log({})
    assert reads == ["tolist"] and t.seq_gaps == [] and t.is_stats is None


def test_u75_active_opsm_logs_its_drops(monkeypatch, tmp_path):
    """OPSM at 0.1: rows 1 (seq_kl 0.2) and 3 (0.3) have negative advantages and are dropped; row 2 (0.3, positive)
    and row 0 (0.05) are kept."""
    x = rows([-0.05], [-0.2], [-0.3], [-0.3], adv=[-1.0, -0.5, 2.0, -0.25])
    m = logged(monkeypatch, tmp_path, "sequence_mask", None, 3.0, x, opsm=0.1)
    assert m["opsm_masked_share"] == 0.5 and m["opsm_masked_adv_sum"] == -0.75
    assert "opsm_masked_share" not in logged(monkeypatch, tmp_path, "sequence_mask", None, 3.0, x)


def test_u36_is_ratio_and_summary_edges():
    lp, x = rows([0.0, 0.0], [1e4], [-1e4], [0.3])
    x["completion_mask"][0] = 0  # an all-masked row
    sp, mask = x["sampling_per_token_logps"], x["completion_mask"]
    for mode in ("sequence_mask", "token_mask", "token_truncate", "sequence_truncate"):
        r = fast.is_ratio(lp, sp, mask, mode, None, 3.0)
        s = fast.is_summary(lp, sp, mask.bool(), r, x["advantages"], mode, None, 3.0)
        assert torch.isfinite(r).all() and torch.isfinite(s[: len(fast.SUMS)]).all(), mode
    with pytest.raises(ValueError, match="clip_max"):
        fast.is_ratio(lp, sp, mask, "token_truncate", 0.5, None)


def test_p13_the_summary_never_changes_the_ratio(monkeypatch, tmp_path):
    """The ratio TRL multiplies by is is_ratio's (A7 pins is_ratio to TRL's block), whatever the summary computes."""
    t = is_trainer("sequence_mask", None, 3.0, tmp_path)
    lp, x = U54
    t.set_is_ratio(lp, x)
    assert torch.equal(x["importance_sampling_ratio"],
                       fast.is_ratio(lp, x["sampling_per_token_logps"], x["completion_mask"], "sequence_mask", None,
                                     3.0))  # fmt: skip


def test_eval_rows_stay_out_of_the_seq_kl_buffer(monkeypatch, tmp_path):
    """R01c unit part: in eval mode the stats still go to the eval metrics, but no seq_kl row is kept."""
    t = is_trainer("sequence_mask", None, 3.0, tmp_path)
    t.model.training = False
    t.set_is_ratio(*U54)
    assert t.seq_gaps == [] and t.is_stats is not None


def test_r01_every_b21_key_once_per_logging_step(monkeypatch, tmp_path):
    counts, log = [], GRPOTrainer.log

    def counted(self, logs, start_time=None):
        if self.model.training and "loss" in logs:  # train steps, not the end-of-train summary
            counts.append({k: len(v) for k, v in self._metrics["train"].items() if k.startswith("sampling/")})
        return log(self, logs, start_time)

    monkeypatch.setattr(GRPOTrainer, "log", counted)
    logs = train_logs(tiny_loop_vllm_fake(monkeypatch, tmp_path))
    assert len(counts) == len(logs) == 3
    for c, h in zip(counts, logs, strict=True):
        assert all(c[f"sampling/{k}"] == 1 for k in B21_KEYS), c  # one value per step, not one per micro-batch
        assert not any(k.startswith("sampling/opsm_") for k in h)  # OPSM off


def read_seq_gaps(out):
    meta = json.loads((out / "seq_gaps.json").read_text())
    a = np.fromfile(out / "seq_gaps.f32", dtype="<f4")
    return meta, a.reshape(meta["steps"], meta["rows_per_step"], 3)


def test_r01b_seq_gaps_file_holds_each_steps_rows_in_order(monkeypatch, tmp_path):
    t = tiny_loop_vllm_fake(monkeypatch, tmp_path)
    rows = t.args.per_device_train_batch_size * t.args.gradient_accumulation_steps  # rows_per_step
    meta, a = read_seq_gaps(tmp_path)
    assert meta == {"rows_per_step": rows, "cols": ["seq_kl", "advantage", "valid_tokens"], "steps": 3}
    assert (tmp_path / "seq_gaps.f32").stat().st_size == 3 * rows * 3 * 4
    planted = t.vllm_generation.planted  # (c, n) per generated row, in generation order
    for k in range(3):  # rows are length-sorted inside a step (by_length): compare as sets of (seq_kl, length)
        got = sorted(zip(a[k, :, 0].tolist(), a[k, :, 2].tolist(), strict=True))
        want = sorted(planted[k * rows:(k + 1) * rows])
        assert [n for _, n in got] == [n for _, n in want]
        np.testing.assert_allclose([c for c, _ in got], [c for c, _ in want], atol=1e-5)
        assert abs(a[k, :, 1].sum()) < 1e-4  # group advantages sum to 0


def test_r01c_an_eval_between_steps_adds_no_rows(monkeypatch, tmp_path):
    q = lambda logs: [[h[f"sampling/seq_kl/{k}"] for k in ("q50", "q99", "q99.5")] for h in logs]  # noqa: E731
    plain = q(train_logs(tiny_loop_vllm_fake(monkeypatch, tmp_path / "plain", steps=2)))
    t = tiny_loop_vllm_fake(monkeypatch, tmp_path / "eval", steps=2, eval_steps=1)
    evals = [h for h in t.state.log_history if "eval_loss" in h]
    assert len(evals) == 2 and all("eval_sampling/importance_sampling_ratio/mean" in h for h in evals)
    assert not any(k.startswith("eval_sampling/seq_kl") for h in evals for k in h)
    meta, a = read_seq_gaps(tmp_path / "eval")
    assert meta["steps"] == 2 and a.shape == (2, 8, 3)
    assert q(train_logs(t)) == plain
    np.testing.assert_array_equal(a, read_seq_gaps(tmp_path / "plain")[1])


CUDA = pytest.mark.skipif(not torch.cuda.is_available(), reason="no CUDA device: runs on the box (b1_rows.yaml step 2)")


def g02_run(device="cuda"):
    """16 update micro-batches of 16 x 1536 through set_is_ratio under sync-debug "error", then one log."""
    import tempfile

    g = torch.Generator(device=device).manual_seed(0)
    t = is_trainer("sequence_mask", None, 3.0, tempfile.mkdtemp())
    torch.cuda.set_sync_debug_mode("error")
    try:
        for _ in range(16):
            lp = -torch.rand(16, 1536, device=device, generator=g) * 3
            x = {"sampling_per_token_logps": lp + torch.randn(16, 1536, device=device, generator=g) * 1e-3,
                 "completion_mask": (torch.rand(16, 1536, device=device, generator=g) < 0.7).long(),
                 "advantages": torch.randn(16, device=device, generator=g)}  # fmt: skip
            t.set_is_ratio(lp, x)
        assert isinstance(t.is_stats, torch.Tensor) and t.is_stats.is_cuda
    finally:
        torch.cuda.set_sync_debug_mode(0)
    return t


@pytest.mark.gpu
@CUDA
def test_g02_is_summary_no_sync(monkeypatch):
    """G02: is_ratio + is_summary on the device raise nothing under sync-debug "error"; the log reads once."""
    t = g02_run()
    monkeypatch.setattr(GRPOTrainer, "log", lambda self, logs, start=None: None)
    t.log({})
    assert t.is_stats is None and t.seq_gaps == []


def test_g02_log_on_a_cuda_build_takes_the_engine_keys(monkeypatch, tmp_path):
    """The G02 trainer (trainer()) logged on a CUDA build at step 1 (faked on the Mac): the D7 engine/* keys come in;
    the box found the helper had no Trainer state there."""
    from objdel import _gpucheck

    got = []
    monkeypatch.setattr(_gpucheck, "build", lambda: "cuda")
    monkeypatch.setattr(_gpucheck, "stats", lambda: {"sync_alloc": 2})
    monkeypatch.setattr(fast.torch.cuda, "max_memory_reserved", lambda *a, **k: 7)
    monkeypatch.setattr(GRPOTrainer, "log", lambda self, logs, start=None: got.append(logs))
    t = is_trainer("sequence_mask", None, 3.0, tmp_path)
    t.set_is_ratio(*rows([0.1, 0.1]))
    t.log({})
    assert got == [{"engine/sync_alloc": 2, "engine/max_memory_reserved": 7}]


@pytest.mark.gpu
@CUDA
def test_g08a_smoke(monkeypatch):
    """G08a: the CUDA checker build (fails, never skips, on the box) and the G02 no-sync run, counterfactuals on."""
    from objdel import _gpucheck

    assert _gpucheck.build() == "cuda"
    test_g02_is_summary_no_sync(monkeypatch)


@pytest.mark.gpu
@CUDA
def test_f04_is_summary_costs_under_2pct_of_a_step():
    """F04: the summary of a step's 16 micro-batches (16 x 1536 each) costs < 2 % of the R5 step (23.6 s, measured
    with --shared-prefix on an RTX 5090)."""
    g = torch.Generator(device="cuda").manual_seed(0)
    lp = -torch.rand(16, 1536, device="cuda", generator=g) * 3
    sp, adv = lp + torch.randn(16, 1536, device="cuda", generator=g) * 1e-3, torch.randn(16, device="cuda")
    mask = torch.rand(16, 1536, device="cuda", generator=g) < 0.7
    r = fast.is_ratio(lp, sp, mask, "sequence_mask", None, 3.0)
    run = lambda: [fast.is_summary(lp, sp, mask, r, adv, "sequence_mask", None, 3.0) for _ in range(16)]  # noqa: E731
    run()
    torch.cuda.synchronize()
    t0 = time.perf_counter()  # host launches and device time together
    run()
    torch.cuda.synchronize()
    assert time.perf_counter() - t0 < 0.02 * 23.6


def rope_rows(prompt, answers, pads, W):
    """Rows laid out as TRL builds them: left pads, the shared prompt, an answer of text tokens, right pads."""
    ids, mask = torch.zeros(len(answers), W, dtype=torch.long), torch.zeros(len(answers), W, dtype=torch.long)
    for r, (a, p) in enumerate(zip(answers, pads, strict=True)):
        row = prompt + a
        ids[r, p:p + len(row)], mask[r, p:p + len(row)] = torch.tensor(row), 1
    return ids, mask, (ids == IMG).long() * mask


def test_rope_positions_equal_get_rope_index_and_flag_rows_off_the_template():
    """C2: one get_rope_index call on the bare prompt gives a template; every row's positions follow from its real
    token count by one cumsum, gather and where. Rows that leave the template set `bad` (no host read inside)."""
    m, g = tiny_model(), torch.Generator().manual_seed(3)
    prompt = [5, 6, START, *[IMG] * 4, END, 7, 8, 9]
    grid = torch.tensor([[1, 4, 4]])
    ans = lambda n: torch.randint(0, 99, (n,), generator=g).tolist()  # noqa: E731
    answers, pads = [[], ans(1), ans(13), ans(20), ans(9)], [3, 0, 2, 0, 1]
    ids, mask, mm = rope_rows(prompt, answers, pads, W=34)
    mask[4, 1 + len(prompt) + 4] = 0  # a hole in an answer
    mask[3, len(prompt):] = 0  # a truncated row whose answer is masked
    mm = mm * mask
    tmpl = fast.rope_template(m.model, torch.tensor(prompt), (torch.tensor(prompt) == IMG).long(), grid)
    pos, bad = fast.rope_positions(tmpl, ids, mask, mm)
    want, _ = m.model.get_rope_index(ids, mm_token_type_ids=mm, image_grid_thw=grid.repeat(5, 1), attention_mask=mask)
    assert torch.equal(pos, want) and pos.dtype == torch.int64 and not bad.any()
    odd = ids.clone()
    odd[0, 3 + 1] = 9  # a prompt token differs
    mm_odd = mm.clone()
    mm_odd[2, 2 + len(prompt) + 3] = 1  # an image token in an answer
    short = mask.clone()
    short[1, 4:] = 0  # fewer real tokens than the prompt
    assert fast.rope_positions(tmpl, odd, mask, mm)[1].tolist() == [True, False, False, False, False]
    assert fast.rope_positions(tmpl, ids, mask, mm_odd)[1].tolist() == [False, False, True, False, False]
    assert fast.rope_positions(tmpl, ids, short, mm * short)[1].tolist() == [False, True, False, False, False]


@functools.cache
def rope_model():
    return tiny_model()


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(rows=st.lists(st.tuples(st.integers(0, 40), st.integers(0, 6), st.booleans()), min_size=1, max_size=6),
       seed=st.integers(0, 2**16))  # fmt: skip
def test_rope_positions_equal_get_rope_index_on_random_rows(rows, seed):
    """C2 property: any answer lengths, left pads, a mask hole or a masked (truncated) answer: rope_positions equals
    get_rope_index on every column and no row is bad."""
    m, g = rope_model(), torch.Generator().manual_seed(seed)
    prompt = [5, 6, START, *[IMG] * 4, END, 7, 8, 9]
    answers = [torch.randint(0, 99, (n,), generator=g).tolist() for n, _, _ in rows]
    ids, mask, mm = rope_rows(prompt, answers, [p for _, p, _ in rows], W=6 + len(prompt) + 40)
    for r, (n, p, odd) in enumerate(rows):
        if odd and n > 1:
            mask[r, p + len(prompt) + int(torch.randint(0, n, (1,), generator=g))] = 0  # a hole in the answer
        elif odd:
            mask[r, p + len(prompt):] = 0  # a truncated answer, masked
    mm = mm * mask
    grid = torch.tensor([[1, 4, 4]])
    tmpl = fast.rope_template(m.model, torch.tensor(prompt), (torch.tensor(prompt) == IMG).long(), grid)
    pos, bad = fast.rope_positions(tmpl, ids, mask, mm)
    want, _ = m.model.get_rope_index(ids, mm_token_type_ids=mm, image_grid_thw=grid.repeat(len(rows), 1),
                                     attention_mask=mask)  # fmt: skip
    assert torch.equal(pos, want) and not bad.any()


def test_rope_positions_equal_get_rope_index_on_the_real_lidc_prompt():
    """C2: the real processor's LIDC prompt (325 tokens, grid (1, 32, 32)) on a model with the real vision token ids:
    answers of 0, 1, 131 and 1,536 tokens, left pads, a mask hole and a truncated masked row equal get_rope_index;
    an image pad in an answer sets bad."""
    import PIL.Image
    from transformers import AutoConfig, AutoProcessor

    from objdel.train import dataset

    try:
        name = "Qwen/Qwen3-VL-2B-Instruct"
        proc, cfg = AutoProcessor.from_pretrained(name, local_files_only=True), AutoConfig.from_pretrained(
            name, local_files_only=True)
    except OSError:
        pytest.skip("the Qwen3-VL processor is not in the local HF cache")
    text = proc.apply_chat_template([dataset._user("lidc")], add_generation_prompt=True, tokenize=False)
    one = proc(text=[text], images=[PIL.Image.new("RGB", (512, 512))], return_tensors="pt")
    prompt, grid = one["input_ids"][0].tolist(), one["image_grid_thw"]
    assert len(prompt) == 325 and grid.tolist() == [[1, 32, 32]]
    m = tiny_model(vocab_size=cfg.text_config.vocab_size, img=cfg.image_token_id, video=cfg.video_token_id,
                   start=cfg.vision_start_token_id, end=cfg.vision_end_token_id)  # fmt: skip
    g = torch.Generator().manual_seed(5)
    answers = [torch.randint(0, 151000, (n,), generator=g).tolist() for n in (0, 1, 131, 1536, 131)]
    ids, mask, _ = rope_rows(prompt, answers, [3, 0, 2, 0, 1], W=325 + 1536 + 3)
    mm = (ids == cfg.image_token_id).long() * mask
    mask[2, 2 + 325 + 7] = 0  # a hole
    mask[4, 1 + 325:] = 0  # truncated and masked
    mm = mm * mask
    tmpl = fast.rope_template(m.model, torch.tensor(prompt), one["mm_token_type_ids"][0], grid)
    pos, bad = fast.rope_positions(tmpl, ids, mask, mm)
    want, _ = m.model.get_rope_index(ids, mm_token_type_ids=mm, image_grid_thw=grid.repeat(5, 1), attention_mask=mask)
    assert torch.equal(pos, want) and not bad.any()
    mm[3, 325 + 10] = 1  # an image token in an answer
    assert fast.rope_positions(tmpl, ids, mask, mm)[1].tolist() == [False, False, False, True, False]


def test_rope_bad_rows_routes_on_host_data_only():
    """U-RP7: a prompt off its template, an answer holding the image pad id or the vision start id; a clean row."""
    tmpl, vision = (1, 2, 3), (IMG, 101, START, END)
    bad = fast.rope_bad_rows([(1, 2, 3), (1, 2, 4), (1, 2, 3), (1, 2, 3)], [[7, 8], [7], [7, IMG], [START]],
                             [tmpl] * 4, vision)  # fmt: skip
    assert bad.dtype == bool and bad.tolist() == [False, True, True, True]


def test_vllm_kwargs_wraps_trls_llm_only_inside_the_block(monkeypatch):
    import trl.generation.vllm_generation as vg

    seen = []
    monkeypatch.setattr(vg, "LLM", lambda **kw: seen.append(kw) or "engine", raising=False)
    outer = vg.LLM
    with fast.vllm_kwargs(dtype="float16"):
        assert vg.LLM(model="m") == "engine"
    assert seen == [{"dtype": "float16", "model": "m"}] and vg.LLM is outer
    with fast.vllm_kwargs():  # nothing to add: TRL's name is left alone
        assert vg.LLM is outer


def _fp16_trainer(engine_dtype, fp16=True):
    from types import SimpleNamespace

    t = object.__new__(fast.FastGRPO)
    t.args = SimpleNamespace(fp16=fp16)
    t.vllm_generation = SimpleNamespace(llm=SimpleNamespace(
        llm_engine=SimpleNamespace(model_config=SimpleNamespace(dtype=engine_dtype))))  # fmt: skip
    return t


def test_fp16_run_refuses_an_engine_that_is_not_float16():
    _fp16_trainer(torch.float16).check_engine_dtype()
    with pytest.raises(RuntimeError, match="float16"):
        _fp16_trainer(torch.bfloat16).check_engine_dtype()
    _fp16_trainer(torch.bfloat16, fp16=False).check_engine_dtype()  # other presets are not checked
    t = _fp16_trainer(torch.bfloat16)
    t.vllm_generation = None  # use_vllm False
    t.check_engine_dtype()


def test_lora_input_cast_is_dropped_under_fp16_autocast_too():
    m = _lora_bf16()
    t = trainer(fast.FastGRPO, vision=None, rope=None, args=type("A", (), {"bf16": False, "fp16": True})(), model=m)
    assert fast.FastGRPO.drop_lora_input_cast(t) is True
    assert all(not mod.cast_input_dtype_enabled for _, mod in m.named_modules()
               if hasattr(mod, "cast_input_dtype_enabled"))  # fmt: skip


@pytest.mark.gpu
@CUDA
def test_g08b_fp16_dtype():
    """G08b fp16 part: vllm_kwargs makes TRL's LLM build a float16 engine on the start checkpoint (a missing
    folder fails, never skips)."""
    import os

    import trl.generation.vllm_generation as vg

    ckpt = os.environ["OBJDEL_RUN_DIR"] + "/ckpt"
    with fast.vllm_kwargs(dtype="float16"):
        llm = vg.LLM(model=ckpt, gpu_memory_utilization=0.3, max_model_len=2048, enforce_eager=True)
    assert llm.llm_engine.model_config.dtype == torch.float16


@pytest.mark.gpu
@CUDA
def test_g08b_smoke():
    """G08b smoke (b3_v5.yaml step 1): the G01 checks on the first 1,000 crops of checker_eval/masks.npz, so the
    reward v5 sums are right on the box before they train anything: sum_x and sum_y equal np.nonzero's, and every
    other field equals the `windows_1000` entry of tests/goldens/b_g01.json."""
    from experiments.phase2.perf.make_b_goldens import windows
    from tests.test_scorer_sums import GOLDEN, _masks_npz, _old_fields_ok, _sums_ok

    gold = json.loads(GOLDEN.read_text())
    sha, ms, r = windows(_masks_npz())
    assert len(ms) == 1000 and sha == gold["windows_1000"]["input_sha256"]
    _sums_ok(r, ms)
    _old_fields_ok(r, gold, "windows_1000")


def test_log_by_prompt_adds_the_per_prompt_rewards_to_train_logs_only(monkeypatch):
    from collections import defaultdict
    from types import SimpleNamespace

    seen = []
    monkeypatch.setattr(GRPOTrainer, "log", lambda self, logs, start_time=None: seen.append(dict(logs)))
    reward = SimpleNamespace(by_prompt=[0.5, 0.25])
    t = object.__new__(fast.FastGRPO)
    t.model, t.is_stats, t.cuda_times, t.last_log = SimpleNamespace(training=True), None, None, None
    t._metrics = {"train": defaultdict(list), "eval": defaultdict(list)}
    t.reward_funcs = [SimpleNamespace(trl_reward=reward)]
    t.log({})  # the default keeps the B1 log
    t.log_by_prompt = True
    t.log({})
    t.model.training = False
    t.log({})
    assert seen == [{}, {"rewards/iou_tc/by_prompt": [0.5, 0.25]}, {}]


def capped_trainer(caps):
    """vllm_trainer with the runaway cap on (train_rl --runaway-cap), its engine recording the sampling params."""
    t, _ = vllm_trainer(expand=1)
    t.runaway_cap = True
    t.rows = [{"cap": c} for c in caps] if caps is not None else [{}] * 2
    seen, llm = [], t.vllm_generation.llm
    generate = llm.generate
    llm.generate = lambda prompts, sampling_params=None, **kw: seen.append(sampling_params) or generate(prompts, **kw)
    fast.FastGRPO.patch_vllm_prompts(t)
    return t, seen


def test_runaway_cap_gives_each_prompt_its_own_max_tokens():
    class Params:
        def __init__(self, max_tokens):
            self.max_tokens = max_tokens

        def clone(self):
            return Params(self.max_tokens)

    t, seen = capped_trainer([5, 9])
    prompts = [{"prompt_token_ids": [START, END]}] * 2
    base = Params(1536)
    t.vllm_generation.llm.generate(prompts, sampling_params=base)
    assert [p.max_tokens for p in seen[0]] == [5, 9] and base.max_tokens == 1536
    t.rows = [{"cap": 1}]
    with pytest.raises(ValueError, match="1 caps for 2 prompts"):
        t.vllm_generation.llm.generate(prompts, sampling_params=base)
    t.rows = None  # g_probe samples outside a generation batch: the global limit stays
    t.vllm_generation.llm.generate(prompts, sampling_params=base)
    assert seen[1] is base


def test_capped_share_counts_answers_that_reached_their_cap():
    t = object.__new__(fast.FastGRPO)
    t.runaway_cap, t._metrics = True, {"train": collections.defaultdict(list)}
    t.model = type("M", (), {"training": True})()
    t.rows = [{"cap": 3}, {"cap": 3}, {"cap": 4}, {"cap": 2}]
    t.log_capped([[1, 2, 3], [1], [1, 2, 3, 4], [1, 2]])
    assert t._metrics["train"]["completions/capped_share"] == [0.75]


def test_beta_zero_runs_only_the_update_forward(monkeypatch, tmp_path):
    """C8: beta 0 loads no reference adapter and runs no reference pass; with one generation per step the old-logprob
    pass is skipped too (C10): one forward per update micro-batch (8 rows, micro-batch 8: one) and nothing else."""
    calls, real = [], fast.FastGRPO._get_per_token_logps_and_entropies
    monkeypatch.setattr(fast.FastGRPO, "_get_per_token_logps_and_entropies",
                        lambda self, *a, **k: calls.append(self.loss_inputs is not None) or real(self, *a, **k))
    t = tiny_loop_vllm_fake(monkeypatch, tmp_path, steps=2)
    assert t.args.beta == 0.0 and calls == [True, True]  # U-CO1: every call is an update forward (loss_inputs set)


def test_dr_grpo_loss_scales_by_the_completion_limit_on_a_one_row_example():
    """U-CO2: dr_grpo divides the summed token loss by rows x max_completion_length (TRL 1.17, grpo_trainer.py 3285),
    so lowering the limit from 1,536 to 1,200 scales the loss by 1,536 / 1,200 = 1.28. One row, 3 answer tokens,
    ratio 1 (the old logprobs are the new ones), advantage 2: the token loss is -2 each, the sum -6."""
    loss = {}
    for limit in (1536, 1200):
        t = trainer(Plain, top_entropy_quantile=1.0, aux_loss_enabled=False, off_policy_mask_threshold=None,
                    importance_sampling_level="token", beta=0.0, loss_type="dr_grpo", epsilon_low=0.2,
                    epsilon_high=0.28, args=type("A", (), {"delta": None})(), use_vllm=False,
                    max_completion_length=limit, current_gradient_accumulation_steps=1,
                    _entropy_bonus_enabled=False)  # fmt: skip
        t._get_per_token_logps_and_entropies = lambda *a, **k: (torch.zeros(1, 3), torch.zeros(1, 3), None)
        t.accelerator.gather, t.accelerator.reduce = (lambda x: x), (lambda x, **_: x)
        loss[limit] = t._compute_loss(None, {"prompt_ids": torch.ones(1, 2, dtype=torch.long),
                                             "prompt_mask": torch.ones(1, 2, dtype=torch.long),
                                             "completion_ids": torch.ones(1, 3, dtype=torch.long),
                                             "completion_mask": torch.ones(1, 3, dtype=torch.long),
                                             "advantages": torch.tensor([2.0])}).item()  # fmt: skip
    assert loss[1536] == -0.00390625 and loss[1200] == pytest.approx(-0.005)  # -6 / 1536 = -1 / 256; -6 / 1200


def test_tokenize_prompts_with_prompt_memo_gives_trls_prompt_ids_tokenizing_each_layout_once():
    """C2 step 1 (--prompt-memo): rows are grouped by (rendered text, image sizes); TRL tokenizes one row per group.
    Prompt ids and images equal TRL's own; the vLLM path drops the multimodal fields."""
    import PIL.Image
    from transformers import AutoProcessor

    try:
        proc = AutoProcessor.from_pretrained("Qwen/Qwen3-VL-2B-Instruct", local_files_only=True)
    except OSError:
        pytest.skip("the Qwen3-VL processor is not in the local HF cache")
    imgs = [PIL.Image.new("RGB", s, (i, 0, 0)) for i, s in enumerate([(64, 64)] * 3 + [(96, 64)])]
    ask = lambda img, q: [{"role": "user", "content": [{"type": "image", "image": img},  # noqa: E731
                                                       {"type": "text", "text": q}]}]  # fmt: skip
    prompts = [ask(imgs[i % 4], "outline the disc" if i % 3 else "outline the cup") for i in range(24)]
    attrs = dict(processing_class=proc, chat_template=None, chat_template_kwargs={}, _is_vlm=True, tools=None,
                 environment_factories=None, use_vllm=True)  # fmt: skip
    want_ids, want_images, _ = trainer(Plain, **attrs)._tokenize_prompts(prompts)
    calls, real = [], GRPOTrainer._tokenize_prompts
    t = trainer(fast.FastGRPO, prompt_memo=True, **attrs)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(GRPOTrainer, "_tokenize_prompts", lambda self, p: calls.append(len(p)) or real(self, p))
        ids, images, fields = t._tokenize_prompts(prompts)
    assert [list(x) for x in ids] == [list(x) for x in want_ids] and fields == {}
    assert all(a[0] is b[0] for a, b in zip(images, want_images, strict=True))
    assert calls == [1] * 4  # 2 texts x 2 image sizes


def meta_for(mask, T, G, bad=()):
    """The RowMeta host lists of a test batch: real prompt and answer tokens; image, group and prompt layout per G
    rows; rows `bad` marked for the HF rope path; 16 pixel rows (one 4 x 4-patch image) per row."""
    n, P = mask.size(0), mask.size(1) - T
    grp = np.arange(n, dtype=np.int32) // G
    return fast.RowMeta(mask[:, :P].sum(1).numpy().astype(np.int32), mask[:, P:].sum(1).numpy().astype(np.int32),
                        grp, grp, grp, np.isin(np.arange(n), bad).astype(np.int32),
                        np.full(n, 16, dtype=np.int32))  # fmt: skip


def no_device_reads(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("device read under host_meta")

    monkeypatch.setattr(fast, "window", boom)
    monkeypatch.setattr(fast, "prompt_groups", boom)


@pytest.mark.parametrize("batch_size", [8, 4, 3])
def test_host_meta_forward_takes_windows_runs_and_lengths_from_the_host(monkeypatch, batch_size):
    """C4 step 3: with host_meta the chunk windows, packed runs and lengths come from RowMeta, same logprobs."""
    m = tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[2, 0])
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=batch_size, **kw)
    no_device_reads(monkeypatch)
    for t in (trainer(Fast, buckets=(4, 8, 12, 16)), shared(m)):
        lp, _, _ = t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=batch_size,
                                                        host_meta=meta_for(mask, T, 4), **kw)  # fmt: skip
        real = mask[:, -T:].bool()
        torch.testing.assert_close(lp[real], ref[real], atol=1e-5, rtol=1e-5)
        assert (lp[~real] == 0).all()


@pytest.mark.parametrize("shared_prefix", [False, True])
def test_host_meta_forward_never_reads_the_grid(monkeypatch, shared_prefix):
    """C4 step 3: under host_meta each row's pixel rows come from RowMeta (image pads x merge^2), so the forward makes
    no Tensor.tolist read of its own; without host_meta it reads the grid once, in grid_rows (control)."""
    import sys

    m = tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[2, 0])
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=4, **kw)
    reads, real = [], torch.Tensor.tolist

    def spy(x):
        f = sys._getframe(1).f_code
        reads.append(f.co_filename.endswith("fast.py") and f.co_name in ("_get_per_token_logps_and_entropies",
                                                                         "grid_rows"))  # fmt: skip
        return real(x)

    monkeypatch.setattr(torch.Tensor, "tolist", spy)
    t = shared(m) if shared_prefix else trainer(Fast, buckets=(4, 8, 12, 16))
    lp, _, _ = t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=4, host_meta=meta_for(mask, T, 4), **kw)
    torch.testing.assert_close(lp, ref, atol=1e-5, rtol=1e-5)
    assert not any(reads)
    t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=4, **kw)
    assert any(reads)


def test_compute_loss_hands_the_micro_batch_host_meta_to_the_forward(monkeypatch):
    """C4: _compute_loss reads the micro-batch's host_meta rows (moved by TRL's split and shuffle) for its forward;
    without them under --host-meta it raises; the old and reference passes get the generation batch's RowMeta (C2:
    their rope positions, windows and pixel rows from the host too), none outside a generation batch."""
    seen = []
    monkeypatch.setattr(fast.FastForward, "_get_per_token_logps_and_entropies",
                        lambda self, *a, **k: seen.append(k.get("host_meta")) or (None, None, None))
    monkeypatch.setattr(GRPOTrainer, "_compute_loss",
                        lambda self, model, inputs: self._get_per_token_logps_and_entropies(model, None, None, 1))
    t = trainer(fast.FastGRPO, host_meta=True, use_vllm=False)
    rows = [fast.MetaRow(5, 3, 0, 0, 0, 0, 16), fast.MetaRow(4, 0, 1, 1, 1, 0, 16)]
    t._compute_loss(None, {"host_meta": rows})
    t._get_per_token_logps_and_entropies(None, None, None, 1)
    t.meta = fast.meta_of({"host_meta": rows[::-1]})
    t._get_per_token_logps_and_entropies(None, None, None, 1)
    assert fast.host_rows(seen[0]) == rows and seen[1] is None and seen[2] is t.meta
    t.meta = None
    with pytest.raises(KeyError, match="host_meta"):
        t._compute_loss(None, {})
    t.host_meta = False
    t._compute_loss(None, {"host_meta": rows})
    assert seen[-1] is None


def test_packed_chunks_stage_their_indices_in_one_fresh_buffer_each(monkeypatch):
    """U-RM6 (C4 step 4): each packed chunk fills its own torch.empty staging buffer (pinned on CUDA) from numpy and
    copies it once; two chunks never share one, and the trainer keeps none."""
    m = tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[2, 0])
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=4, **kw)
    bufs, real_empty = [], torch.empty

    def spy(*a, pin_memory=False, **k):
        out = real_empty(*a, **k)
        if a and isinstance(a[0], int) and k.get("dtype") is torch.int64:
            bufs.append((out, pin_memory))
        return out

    monkeypatch.setattr(fast.torch, "empty", spy)
    t = shared(m)
    lp, _, _ = t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=4, **kw)
    torch.testing.assert_close(lp, ref, atol=1e-5, rtol=1e-5)
    assert len(bufs) == 2 and bufs[0][0].data_ptr() != bufs[1][0].data_ptr()
    assert [p for _, p in bufs] == [False, False]  # CPU here; pin_memory=True when the rows are on CUDA
    assert not any(isinstance(v, torch.Tensor) for v in vars(t).values())


def test_packed_chunks_gather_their_rows_with_index_select(monkeypatch):
    """C4 step 4: packed_chunk gathers tok, mmt and posp from the row tensors with index_select (host indices, one
    staging copy) instead of per-row torch.full / torch.cat loops; same logprobs. Per chunk: 3 gathers, the labels'
    one torch.full, and the 2 image cats (pixel_values, image_grid_thw)."""
    import sys

    m = tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[2, 0])
    ref, ent_ref, _ = shared(m)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=4,
                                                                   compute_entropy=True, **kw)  # fmt: skip
    calls, real = collections.Counter(), {n: getattr(torch, n) for n in ("full", "cat", "index_select")}
    for name in real:
        def spy(*a, _name=name, **k):
            calls[_name] += sys._getframe(1).f_code.co_name == "packed_chunk"
            return real[_name](*a, **k)

        monkeypatch.setattr(torch, name, spy)
    lp, ent, _ = shared(m)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=4, compute_entropy=True,
                                                              **kw)  # fmt: skip
    assert torch.equal(lp, ref) and torch.equal(ent, ent_ref)
    assert calls == {"index_select": 6, "full": 2, "cat": 4}  # 2 packed chunks


def counted_rope(m):
    calls, orig = [], m.model.get_rope_index

    def get_rope_index(input_ids, **k):
        calls.append(input_ids.size(0))
        return orig(input_ids, **k)

    m.model.get_rope_index = get_rope_index
    return calls


def test_host_meta_rope_positions_come_from_one_template_per_layout_bad_rows_from_hf():
    """C2 step 2: under host_meta each layout's template is built by one HF call, good rows take rope_positions, the
    host's bad rows take get_rope_index in the same chunk (rope/fallback_rows); positions equal HF's on every row."""
    m = tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[2, 0])
    rope = fast.RopeCache.install(m)
    ref = rope.positions(ids, mask, kw["mm_token_type_ids"], kw["image_grid_thw"], kw["num_images"])
    rope.clear()
    calls = counted_rope(m)
    pos, nbad = rope.template_positions(ids, mask, kw["mm_token_type_ids"], kw["image_grid_thw"], kw["num_images"],
                                        ids.size(1) - T, meta_for(mask, T, 4, bad=[5]), check=False)  # fmt: skip
    real = mask.bool()
    assert torch.equal(pos[:, real], ref[:, real]) and nbad == 1
    assert calls == [1, 1, 1]  # two layout templates, one bad row
    rope.clear()
    with pytest.MonkeyPatch.context() as mp:  # C2: the keyed cache is the no-host_meta path only
        mp.setattr(fast.RopeCache, "keys", None)
        assert torch.equal(rope.template_positions(ids, mask, kw["mm_token_type_ids"], kw["image_grid_thw"],
                                                   kw["num_images"], ids.size(1) - T, meta_for(mask, T, 4, bad=[5]),
                                                   check=True)[0], pos)  # fmt: skip
    plain = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=4, **kw)[0]
    t = trainer(Fast, shared_prefix=True, rope=rope)
    rope.clear()
    lp = t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=4, host_meta=meta_for(mask, T, 4, bad=[5]),
                                              **kw)[0]  # fmt: skip
    torch.testing.assert_close(lp, plain, atol=1e-5, rtol=1e-5)
    assert t._metrics["train"]["rope/fallback_rows"] == [1]


def test_check_rope_compares_every_row_with_hf_on_check_steps(monkeypatch):
    """C2 step 3 (--check-rope-every): on a check step a template position that differs from HF raises; a row bad on
    the device but good on the host raises; other steps never compare."""
    m = tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[2, 0])
    args = (ids, mask, kw["mm_token_type_ids"], kw["image_grid_thw"], kw["num_images"], ids.size(1) - T,
            meta_for(mask, T, 4))  # fmt: skip
    rope = fast.RopeCache.install(m)
    rope.template_positions(*args, check=True)
    real = fast.rope_positions
    monkeypatch.setattr(fast, "rope_positions", lambda *a: (real(*a)[0] + 1, real(*a)[1]))
    rope.clear()
    with pytest.raises(AssertionError, match="rope positions differ from HF"):
        rope.template_positions(*args, check=True)
    rope.clear()
    rope.template_positions(*args, check=False)
    monkeypatch.setattr(fast, "rope_positions", lambda *a: (real(*a)[0], torch.ones_like(real(*a)[1])))
    rope.clear()
    with pytest.raises(AssertionError, match="bad on the device"):
        rope.template_positions(*args, check=True)


def test_vllm_engine_kwargs_follow_the_config_and_default_to_none():
    from types import SimpleNamespace as NS

    assert fast.vllm_engine_kwargs(NS(fp16=False)) == {}  # the B1 control: TRL's engine untouched
    assert fast.vllm_engine_kwargs(NS(fp16=True)) == {"dtype": "float16"}
    assert fast.vllm_engine_kwargs(NS(fp16=False, vllm_lora=True)) == {"enable_lora": True, "max_lora_rank": 16}


def test_generate_carries_the_lora_request_only_under_option_f(monkeypatch, tmp_path):
    import sys
    from types import SimpleNamespace as NS

    monkeypatch.setitem(sys.modules, "vllm", NS())
    monkeypatch.setitem(sys.modules, "vllm.lora", NS())
    monkeypatch.setitem(sys.modules, "vllm.lora.request", NS(LoRARequest=lambda *a: ("lora", *a)))
    for lora_dir, want in ((None, {}), (tmp_path, {"lora_request": ("lora", "policy", 1, str(tmp_path))})):
        t, _ = vllm_trainer(expand=1)
        t.lora_dir, seen = lora_dir, []
        gen = t.vllm_generation.llm.generate
        t.vllm_generation.llm.generate = lambda p, seen=seen, gen=gen, **kw: seen.append(kw) or gen(p, **kw)
        fast.FastGRPO.patch_vllm_prompts(t)
        t.vllm_generation.llm.generate([{"prompt_token_ids": [START, END]}], sampling_params=None)
        assert seen == [{"sampling_params": None, **want}]
