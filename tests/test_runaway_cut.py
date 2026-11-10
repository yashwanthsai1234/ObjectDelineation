"""C7 outline limit and C11 per-row answer cap: the cap rule, per-row sampling params, the run guards and the sizing
script (U-RC1, U-RC2, U-RC4, U-RC5, U-CO4)."""

import json

import pytest

from experiments.phase2 import runaway_cut, train_rl
from objdel.train import configs, runaway


def test_caps_hand_values_min_with_limit_and_end_token_counted():
    """U-RC1: cap = min(L, ceil(factor x T) + margin); T already counts the end token."""
    assert runaway.runaway_caps([100, 1000, 1182], 1536, 1.5, 64) == [214, 1536, 1536]
    assert runaway.runaway_caps([101], 1536, 1.5, 64) == [216]  # ceil(151.5) = 152
    assert runaway.runaway_caps([1000], 1200, 1.5, 64) == [1200]
    assert runaway.runaway_caps([], 1536, 1.5, 64) == []


def test_outline_tokens_counts_the_end_token_and_skips_none():
    assert runaway.outline_tokens(len, ["abcd", None, "ab"]) == 5


def test_per_row_params_order_length_and_no_mutation():
    """U-RC2."""
    class Params:  # vLLM's SamplingParams (not installed on the Mac): clone() returns an independent copy
        def __init__(self, **kw):
            self.__dict__.update(kw)

        def clone(self):
            return Params(**self.__dict__)

    base = Params(max_tokens=1536, temperature=1.0)
    rows = runaway.per_row_params(base, [10, 20, 30])
    assert [r.max_tokens for r in rows] == [10, 20, 30]
    assert all(r.temperature == 1.0 for r in rows) and base.max_tokens == 1536
    assert len({id(r) for r in rows}) == 3


def args(**kw):
    return train_rl.parse(["--model", "qwen3-vl-2b", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "1",
                           *[x for k, v in kw.items() for x in (f"--{k.replace('_', '-')}", v)]])


def test_runaway_cap_is_off_by_default():
    assert args().runaway_cap == 0
    train_rl.check_flags(args())  # no raise


def test_runaway_cap_needs_unmasked_truncated_rows():
    """U-RC5: a capped answer is cut mid-text; masking it would bias the update toward more capping."""
    with pytest.raises(ValueError, match="mask_truncated_completions"):
        train_rl.check_flags(args(runaway_cap="1.5"))
    with pytest.raises(ValueError, match="mask_truncated_completions"):
        train_rl.check_flags(args(runaway_cap="1.5", grpo='{"mask_truncated_completions": true}'))


def test_outline_limit_is_never_raised_above_the_model_budget():
    """C7 step 1: L is set from the count, never raised."""
    with pytest.raises(ValueError, match="1536"):
        train_rl.check_flags(args(grpo='{"max_completion_length": 1800}'))
    train_rl.check_flags(args(grpo='{"max_completion_length": 1200}'))


def test_soft_overlong_formula_values():
    """U-CO4: TRL's DAPO penalty at L = 1200, c = 200 is 0 up to 1,000 tokens, then falls linearly to -1 at 1,200."""
    from trl.rewards import get_soft_overlong_punishment

    penalty = get_soft_overlong_punishment(max_completion_len=1200, soft_punish_cache=200)
    assert penalty([[1] * n for n in (999, 1000, 1100, 1200, 1201)]) == [0.0, 0.0, -0.5, -1.0, -1.0]


def test_soft_overlong_rule():
    """U-CO4: L - c >= longest doctor outline + 32 (1,214 today)."""
    configs.check_soft_overlong(1536, 322)
    with pytest.raises(ValueError, match="1214"):
        configs.check_soft_overlong(1536, 323)
    with pytest.raises(ValueError):
        configs.check_soft_overlong(1200, 200)  # the penalty would start at 1,000 tokens
    assert configs.LONGEST_OUTLINE == 1182


def test_beta_zero_builds_and_b1_keeps_1e3(tmp_path):
    """U-CO5 (config half): beta 0 is a plain override; b1_rows.yaml never sets it."""
    from pathlib import Path

    assert configs.grpo_config("qwen3-vl-2b", str(tmp_path), use_cpu=True, use_vllm=False, beta=0.0).beta == 0.0
    plan = (Path(__file__).resolve().parents[1] / "experiments/phase2/configs/b1_rows.yaml").read_text()
    assert "beta" not in plan.replace("beta 1e-3", "")  # the comment may name it; no key overrides it
    assert configs.grpo_config("qwen3-vl-2b", str(tmp_path), use_cpu=True, use_vllm=False).beta == 1e-3


def fixture(tmp_path, tokens, ious, name="eval_val"):
    run = tmp_path / name
    run.mkdir()
    (run / "answers.jsonl").write_text(json.dumps({"item_id": "A", "tokens": tokens}) + "\n")
    (run / "metrics.jsonl").write_text(json.dumps(
        {"item_id": "A", "samples": {"legal": [1] * len(tokens), "iou_tc": ious}}) + "\n")
    return run


def test_sizing_hand_values():
    """U-RC4: T = 100 (cap 214 at 1.5, 264 at 2.0). Found answers are IoU_tc >= 0.1."""
    recs = [{"item_id": "A", "targets": ["x" * 99, None]}]
    ans = [{"item_id": "A", "tokens": [50, 230, 300, 900]}]
    met = [{"item_id": "A", "samples": {"legal": [1, 1, 1, 1], "iou_tc": [0.8, 0.5, 0.05, 0.9]}}]
    r = runaway_cut.size(recs, ans, met, len, 1536, [1.5, 2.0], block=32)
    assert r[1.5]["cap_min"] == 214 and r[1.5]["found"] == 3
    assert r[1.5]["found_above_cap"] == 2 and r[1.5]["share"] == pytest.approx(2 / 3)  # 230 and 900
    assert r[2.0]["found_above_cap"] == 1  # only 900
    assert r[1.5]["seconds_removed_per_step"] == pytest.approx(0.0043 * (900 - 214))


def test_sizing_asserts_no_outline_above_its_cap():
    recs = [{"item_id": "A", "targets": ["x" * 2000]}]  # 2001 tokens: the cap is L = 1536 < outline
    ans = [{"item_id": "A", "tokens": [10]}]
    met = [{"item_id": "A", "samples": {"legal": [1], "iou_tc": [0.5]}}]
    with pytest.raises(AssertionError):
        runaway_cut.size(recs, ans, met, len, 1536, [1.5], block=32)


def test_sizing_refuses_a_test_split(tmp_path):
    run = fixture(tmp_path, [10], [0.5], name="eval_test")
    with pytest.raises(ValueError, match="test"):
        runaway_cut.refuse_test(run / "answers.jsonl", "val")
    runaway_cut.refuse_test(fixture(tmp_path, [10], [0.5]) / "answers.jsonl", "val")
    with pytest.raises(ValueError, match="test"):
        runaway_cut.refuse_test(tmp_path / "eval_val" / "answers.jsonl", "test")


def test_row_caps_follow_each_rows_own_outline():
    """C11 start-up: one cap per RL row, from the record of the row's key (a view has its own record and targets)."""
    records = [{"item_id": "a", "targets": ["x" * 10, None]}, {"item_id": "b", "targets": ["y" * 1000, "y"]}]
    caps = runaway.row_caps(["b", "a", "a"], records, len, 1536, 1.5, 64)
    assert caps == [1536, 81, 81]  # a: ceil(1.5 x 11) = 17, + 64; b: capped by L


def test_runaway_cap_and_group_requests_go_together():
    """U-RC3: one n = G request per block carries the block's one cap (vllm_requests.group_params)."""
    a = args(runaway_cap="1.5", grpo='{"mask_truncated_completions": false}')
    a.group_requests = True
    train_rl.check_flags(a)  # no raise
