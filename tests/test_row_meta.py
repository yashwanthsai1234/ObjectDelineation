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
