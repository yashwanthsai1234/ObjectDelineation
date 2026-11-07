"""objdel.train.fast: per-chunk trimming of the GRPO forwards and the per-step vision cache (perf C1, C6, RL-6).

The reference for every equality test is TRL 1.17's own `_get_per_token_logps_and_entropies` on the untrimmed batch,
run on a tiny random Qwen3-VL (CPU, fp32) with TRL's fused LM head.
"""

import contextlib

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from transformers import Qwen3VLConfig, Qwen3VLForConditionalGeneration  # noqa: E402
from trl import GRPOTrainer  # noqa: E402
from trl.trainer.utils import add_fused_lm_head  # noqa: E402

from objdel.train import fast  # noqa: E402

IMG, START, END = 100, 102, 103
PATCH = 3 * 2 * 16 * 16  # channels x temporal patch x patch x patch


def tiny_model():
    torch.manual_seed(0)
    cfg = Qwen3VLConfig(
        text_config=dict(vocab_size=128, hidden_size=32, intermediate_size=64, num_hidden_layers=3,
                         num_attention_heads=2, num_key_value_heads=1, head_dim=16,
                         rope_scaling={"rope_type": "default", "mrope_section": [2, 3, 3], "mrope_interleaved": True}),
        vision_config=dict(depth=2, hidden_size=32, intermediate_size=64, num_heads=2, out_hidden_size=32,
                           patch_size=16, spatial_merge_size=2, temporal_patch_size=2, deepstack_visual_indexes=[0]),
        image_token_id=IMG, video_token_id=101, vision_start_token_id=START, vision_end_token_id=END)  # fmt: skip
    m = Qwen3VLForConditionalGeneration(cfg).eval()
    add_fused_lm_head(m, outputs=("log_probs", "entropy"))
    return m


class _Accel:
    autocast = staticmethod(contextlib.nullcontext)
    is_main_process = True


def trainer(cls, **attrs):
    """A GRPOTrainer subclass instance with only what the forward needs (no model load, no vLLM)."""
    from collections import defaultdict

    t = object.__new__(cls)
    t.accelerator, t._metrics = _Accel(), {"train": defaultdict(list), "eval": defaultdict(list)}
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
                   vllm_importance_sampling_clip_min=None, vllm_importance_sampling_clip_max=3.0)  # fmt: skip


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


def test_update_forward_gives_trl_its_importance_ratio(monkeypatch):
    """The update micro-batch's own logprobs (detached) stand in for the old ones: the ratio and TRL's sampling
    metrics are set from them before TRL's loss multiplies by inputs["importance_sampling_ratio"]."""
    lp = torch.tensor([[0.0, -1.0, -2.0], [-0.5, -0.5, -0.5]], requires_grad=True)
    monkeypatch.setattr(fast.FastForward, "_get_per_token_logps_and_entropies", lambda self, *a, **k: (lp, None, None))
    t = c10_trainer()
    inputs = {"sampling_per_token_logps": torch.tensor([[0.0, -1.5, -1.0], [-0.5, -0.5, -0.5]]),
              "completion_mask": torch.tensor([[1, 1, 0], [1, 1, 1]])}  # fmt: skip
    t.loss_inputs = inputs
    assert t._get_per_token_logps_and_entropies(None, None, None, 3)[0] is lp
    torch.testing.assert_close(inputs["importance_sampling_ratio"], torch.tensor([[1.6487212707], [1.0]]))
    m = t._metrics["train"]
    assert m["sampling/sampling_logp_difference/max"] == [pytest.approx(0.5)]
    assert m["sampling/importance_sampling_ratio/min"] == [1.0]
    t.loss_inputs = None  # outside the loss (old and reference passes): nothing set
    t._get_per_token_logps_and_entropies(None, None, None, 3)
    assert len(m["sampling/importance_sampling_ratio/min"]) == 1
