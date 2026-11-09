"""fast.RowMeta (plan C4): the generation batch's row lengths, image ids and groups from the host lists vLLM returns,
and the window / order derivations from them, equal to today's device derivations (G-EQ3 on synthetic batches)."""

import numpy as np
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from objdel.train import fast  # noqa: E402

EOS, PAD = 7, 0


def rows(seed, n=16, G=4, P=12, T=10):
    """Host lists (prompt ids, completion ids) and the TRL-padded device masks of the same rows."""
    rng = np.random.default_rng(seed)
    plens = rng.integers(3, P + 1, n)
    prompt = [list(rng.integers(10, 99, k)) for k in plens]
    comp = []
    for _ in range(n):
        L = int(rng.integers(1, T + 1))
        c = list(rng.integers(10, 99, L))
        if rng.random() < 0.7:
            c[-1] = EOS  # ended; else truncated at its last token
        comp.append(c)
    pm = torch.tensor([[0] * (P - len(p)) + [1] * len(p) for p in prompt])
    return prompt, comp, pm, P, T, G


def device_masks(prompt, comp, pm, T, mask_truncated):
    cm = torch.tensor([[1] * len(c) + [0] * (T - len(c)) for c in comp])
    if mask_truncated:
        cm *= torch.tensor([c[-1] in (EOS, PAD) for c in comp]).long()[:, None]  # TRL: grpo_trainer.py 2604-2608
    return cm


@pytest.mark.parametrize("mask_truncated", [True, False])
@pytest.mark.parametrize("seed", range(5))
def test_row_meta_lengths_equal_the_device_mask_sums(seed, mask_truncated):
    prompt, comp, pm, P, T, G = rows(seed)
    meta = fast.row_meta(prompt, comp, list(range(len(prompt))), G, mask_truncated, [EOS, PAD])
    cm = device_masks(prompt, comp, pm, T, mask_truncated)
    assert meta.prompt_len.tolist() == pm.sum(1).tolist() and meta.completion_len.tolist() == cm.sum(1).tolist()
    assert meta.group.tolist() == [i // G for i in range(len(prompt))]


def test_row_meta_refuses_an_empty_answer():
    with pytest.raises(ValueError):
        fast.row_meta([[1, 2]], [[]], [0], 1, True, [EOS, PAD])


@pytest.mark.parametrize("seed", range(5))
def test_window_order_and_group_order_from_meta_equal_the_device_derivations(seed):
    prompt, comp, pm, P, T, G = rows(seed)
    meta = fast.row_meta(prompt, comp, list(range(len(prompt))), G, True, [EOS, PAD])
    cm = device_masks(prompt, comp, pm, T, True)
    mask = torch.cat([pm, cm], 1)
    for sel in (list(range(8)), list(range(8, 16)), [3, 1, 12]):
        for buckets in (fast.FINE, (8, 16, 20), (13,)):
            want = fast.window(mask[sel], T, buckets)
            assert fast.window_from_meta(meta, np.array(sel), P, T, buckets) == want
    batch = {"prompt_mask": pm, "completion_mask": cm, "i": list(range(len(prompt)))}
    assert fast.length_order(meta).tolist() == fast.by_length(batch)["i"]
    assert fast.group_order(meta, G).tolist() == fast.by_group_length(batch, G)["i"]


def test_image_ids_are_dense_by_object_in_first_appearance_order():
    a, b, c = object(), object(), object()
    assert fast.image_ids([{"image": x} for x in (b, b, a, b, c, a)]) == [0, 0, 1, 0, 2, 1]


def test_host_meta_orders_without_the_device_masks_and_moves_with_the_rows():
    """U-RM4 / C4 step 3: by_length and by_group_length take the order from host_meta when the batch has it."""
    prompt, comp, pm, P, T, G = rows(3)
    meta = fast.row_meta(prompt, comp, list(range(len(prompt))), G, True, [EOS, PAD])
    cm = device_masks(prompt, comp, pm, T, True)
    host = fast.host_rows(meta)
    assert [r.prompt_len for r in host] == pm.sum(1).tolist()
    for order in (fast.by_length, lambda b: fast.by_group_length(b, G)):
        want = order({"prompt_mask": pm, "completion_mask": cm, "i": list(range(len(prompt)))})
        got = order({"prompt_mask": None, "completion_mask": None, "host_meta": host, "i": list(range(len(prompt)))})
        assert got["i"] == want["i"] and [r.prompt_len for r in got["host_meta"]] == pm.sum(1)[want["i"]].tolist()


def generation_trainer(monkeypatch, every, step=0):
    """A FastGRPO with TRL's batch builder replaced: _generate returns host lists, the batch has their masks."""
    from collections import defaultdict

    prompt, comp, pm, P, T, G = rows(5)
    cm = device_masks(prompt, comp, pm, T, True)
    imgs = [object() for _ in range(len(prompt) // G)]
    out = (prompt, comp, None, ["x"] * len(comp), None, None, {}, None, None)
    monkeypatch.setattr(fast.GRPOTrainer, "_generate", lambda self, prompts: out)
    t = object.__new__(fast.FastGRPO)
    t.host_meta, t.check_meta_every, t.use_vllm, t.images = True, every, False, None
    t.num_generations, t.mask_truncated_completions, t.eos_token_ids = G, True, [EOS]
    t._tokenizer = type("Tok", (), {"pad_token_id": PAD})()
    t.reward_funcs, t._metrics, t.model = [], {"train": defaultdict(list)}, torch.nn.Linear(1, 1).train()
    t.state = type("S", (), {"global_step": step})()
    t.vllm_importance_sampling_correction = False
    t.vision = t.rope = None

    def build(self, inputs):
        self._generate(prompt)
        return {"prompt_mask": pm, "completion_mask": cm.clone()}

    monkeypatch.setattr(fast.GRPOTrainer, "_generate_and_score_completions", build)
    inputs = [{"image": imgs[i // G]} for i in range(len(prompt))]
    return t, inputs, G


def test_generate_puts_row_aligned_host_meta_in_the_batch(monkeypatch):
    """C4 steps 1-2: --host-meta builds RowMeta in _generate and the batch carries one MetaRow per row."""
    t, inputs, G = generation_trainer(monkeypatch, every=1)
    batch = t._generate_and_score_completions(inputs)
    host = batch["host_meta"]
    assert [r.prompt_len for r in host] == batch["prompt_mask"].sum(1).tolist()
    assert [r.completion_len for r in host] == batch["completion_mask"].sum(1).tolist()
    assert [r.image_id for r in host] == [i // G for i in range(len(host))] == [r.group for r in host]
    assert t.meta is None


def test_generate_without_host_meta_keeps_the_batch(monkeypatch):
    t, inputs, _ = generation_trainer(monkeypatch, every=1)
    t.host_meta = False
    assert "host_meta" not in t._generate_and_score_completions(inputs)


@pytest.mark.parametrize(("every", "step", "checks"), [(1, 7, True), (50, 0, True), (50, 50, True), (50, 1, False),
                                                       (0, 0, False)])  # fmt: skip
def test_check_meta_every_compares_the_device_derivation(monkeypatch, every, step, checks):
    """U-RM5: a host length one off raises naming field and row on a checked step; other steps never compare."""
    t, inputs, _ = generation_trainer(monkeypatch, every, step)
    real = fast.host_rows
    monkeypatch.setattr(fast, "host_rows", lambda m: [r._replace(completion_len=r.completion_len + (i == 3))
                                                      for i, r in enumerate(real(m))])  # fmt: skip
    if checks:
        with pytest.raises(AssertionError, match="row 3: completion_len"):
            t._generate_and_score_completions(inputs)
    else:
        t._generate_and_score_completions(inputs)
