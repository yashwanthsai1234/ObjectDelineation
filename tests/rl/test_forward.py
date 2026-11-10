"""objdel.rl.forward: the pure parts of the GRPO forwards (the window, the update orders, packing, the importance-ratio
modes and their summary); from the old test_train_fast.py (E5)."""

import pytest

torch = pytest.importorskip("torch")
from objdel.rl import forward as rl_forward  # noqa: E402


def test_window_cuts_to_the_longest_real_answer_rounded_to_a_bucket():
    # prompt width 12 (left pads 3 and 5 -> 9 real columns), completion width 10, longest real answer 4
    mask = torch.ones(2, 22, dtype=torch.long)
    mask[0, :3], mask[1, :5] = 0, 0
    mask[0, 12 + 4:], mask[1, 12 + 2:] = 0, 0
    assert rl_forward.window(mask, 10, buckets=(8, 16, 20)) == (3, 7)  # need 9 + 4 = 13 -> bucket 16: keep 16 - 9 = 7
    assert rl_forward.window(mask, 10, buckets=(13,)) == (3, 4)  # an exact fit keeps the 3 left pads, no more
    assert rl_forward.window(mask, 10, buckets=(64,)) == (0, 10)  # a bucket wider than the batch: nothing to cut
    mask[:, 12:] = 0  # no real answer token: keep 1 column (TRL's labels[:, :-keep] needs keep >= 1)
    assert rl_forward.window(mask, 10, buckets=(9, 10)) == (3, 1)


def test_is_ratio_follows_trls_modes_on_a_worked_example():
    """TRL 1.17's vLLM importance-sampling ratio (grpo_trainer.py 2776-2824): token gaps 0, 0.5 and a masked one;
    the sequence gap is 0.5; exp(0.5) = 1.6487."""
    lp, smp, m = torch.tensor([[0.0, -1.0, -2.0]]), torch.tensor([[0.0, -1.5, -1.0]]), torch.tensor([[1.0, 1.0, 0.0]])
    e = 1.6487212707
    cases = {("token_truncate", None, 1.5): [[1.0, 1.5, 1.0]], ("token_mask", None, 1.5): [[1.0, 0.0, 1.0]],
             ("sequence_mask", None, 3.0): [[e]], ("sequence_mask", None, 1.5): [[0.0]],
             ("sequence_truncate", None, 1.5): [[1.5]], ("sequence_mask", 2.0, None): [[0.0]]}  # fmt: skip
    for (mode, lo, hi), want in cases.items():
        torch.testing.assert_close(rl_forward.is_ratio(lp, smp, m, mode, lo, hi), torch.tensor(want), msg=mode)
    nan = torch.tensor([[0.0, float("nan"), -1.0]])  # a token vLLM could not score: no correction
    torch.testing.assert_close(rl_forward.is_ratio(lp, nan, m, "token_mask", None, 3.0), torch.tensor([[1.0, 1.0,
        1.0]]))


def test_by_length_orders_rows_longest_first_and_moves_every_field_with_them():
    """TRL shuffles the generation batch before splitting it into update micro-batches, so each chunk is as wide as
    its longest row; sorted, rows of like length share a chunk and FastForward trims each to its own width."""
    d = {"prompt_mask": torch.tensor([[0, 1], [1, 1], [1, 1], [0, 1]]),
         "completion_mask": torch.tensor([[1, 0, 0], [1, 1, 1], [1, 0, 0], [1, 1, 0]]),
         "pixel_values": ["p0", "p1", "p2", "p3"], "num_items_in_batch": torch.tensor(9), "tool_mask": None}
    out = rl_forward.by_length(d)  # real lengths 2, 5, 3, 3: ties keep their order
    assert out["pixel_values"] == ["p1", "p2", "p3", "p0"]
    assert out["completion_mask"].tolist() == [[1, 1, 1], [1, 0, 0], [1, 1, 0], [1, 0, 0]]
    assert out["num_items_in_batch"] == 9 and out["tool_mask"] is None


def test_by_group_length_keeps_each_group_together_longest_group_first():
    n = torch.tensor([3, 9, 2, 4, 5, 5, 1, 1, 7, 2, 2, 2])  # 3 groups of 4: longest members 9, 5, 7
    batch_ = {"prompt_mask": torch.ones(12, 1, dtype=torch.long),
              "completion_mask": torch.zeros(12, 9, dtype=torch.long),
              "tag": list(range(12)), "advantages": torch.arange(12.0)}
    for r, L in enumerate(n.tolist()):
        batch_["completion_mask"][r, :L] = 1
    out = rl_forward.by_group_length(batch_, 4)
    assert out["tag"] == [0, 1, 2, 3, 8, 9, 10, 11, 4, 5, 6, 7]
    assert out["advantages"].tolist() == [float(i) for i in out["tag"]]


def test_packing_is_chosen_only_when_it_costs_less():
    """Packed attention is dense over the whole packed row (prompt + 8 answers): with long answers its square
    outgrows the 8 prompts it saves. Real sizes: 325-token prompts, 16 rows, linear cost ~12,000x an attention pair."""
    short = [127] * 16
    assert rl_forward.packing_wins(325, short, groups=2)
    assert not rl_forward.packing_wins(325, [1536] * 16, groups=2)  # answers at the cap: per-row is cheaper
    assert not rl_forward.packing_wins(325, short, groups=16)  # nothing shared


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


def test_u37_summary_is_float32_on_the_input_device():
    lp, x = rows([0.1, 0.2], [0.3])
    s = rl_forward.is_summary(lp, x["sampling_per_token_logps"], x["completion_mask"].bool(),
                        rl_forward.is_ratio(lp, x["sampling_per_token_logps"], x["completion_mask"], "token_mask",
                            None, 3.0),
                        x["advantages"], "token_mask", None, 3.0)  # fmt: skip
    assert s.dtype == torch.float32 and s.device == lp.device and s.shape == (len(rl_forward.SUMMARY),)


def test_u62_counterfactual_opsm_equals_trls_off_policy_mask():
    """Oracle: TRL's own GRPOTrainer.get_off_policy_mask (imported here only), drops = rows it returns 0 for."""
    GRPOTrainer = pytest.importorskip("trl").GRPOTrainer
    g = torch.Generator().manual_seed(0)
    names = list(rl_forward.SUMMARY)
    for _ in range(200):
        lp = -torch.rand(6, 9, generator=g) * 3
        sp = lp + torch.randn(6, 9, generator=g) * torch.rand(6, 1, generator=g)
        sp[0, 2] = float("nan")
        mask = (torch.rand(6, 9, generator=g) < 0.8).long()
        mask[1] = 0
        adv = torch.randn(6, generator=g)
        r = rl_forward.is_ratio(lp, sp, mask, "sequence_mask", None, 3.0)
        s = rl_forward.is_summary(lp, sp, mask.bool(), r, adv, "sequence_mask", None, 3.0).tolist()
        for d in rl_forward.OPSM_DELTAS:
            keep = GRPOTrainer.get_off_policy_mask(adv[:, None], lp, sp, mask.float(), d)
            assert s[names.index(f"cf_opsm_{d}")] == float((keep == 0).sum())


def test_u36_is_ratio_and_summary_edges():
    lp, x = rows([0.0, 0.0], [1e4], [-1e4], [0.3])
    x["completion_mask"][0] = 0  # an all-masked row
    sp, mask = x["sampling_per_token_logps"], x["completion_mask"]
    for mode in ("sequence_mask", "token_mask", "token_truncate", "sequence_truncate"):
        r = rl_forward.is_ratio(lp, sp, mask, mode, None, 3.0)
        s = rl_forward.is_summary(lp, sp, mask.bool(), r, x["advantages"], mode, None, 3.0)
        assert torch.isfinite(r).all() and torch.isfinite(s[: len(rl_forward.SUMS)]).all(), mode
    with pytest.raises(ValueError, match="clip_max"):
        rl_forward.is_ratio(lp, sp, mask, "token_truncate", 0.5, None)
