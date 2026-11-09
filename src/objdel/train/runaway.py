"""C11 runaway cut: each answer is capped by its own doctor outlines instead of one global limit.

cap_i = min(L, ceil(factor x T_i) + margin), with T_i the longest doctor outline of the row's view in tokens, end
token counted (the sft_rows rule). A capped answer is cut mid-text: it is unparseable, gets R = 0 and a negative
advantage, so `--runaway-cap` needs mask_truncated_completions=False (train_rl.check_flags)."""
import math


def outline_tokens(tokens_of, texts) -> int:
    """T_i: the longest target in tokens, end token counted; tokens_of(text) counts none (None targets are skipped)."""
    return max(tokens_of(t) + 1 for t in texts if t is not None)


def runaway_caps(target_tokens, limit: int, factor: float, margin: int) -> list[int]:
    """Per-row max_tokens cap: factor x the row's target token count plus margin, never above limit."""
    return [min(limit, math.ceil(factor * t) + margin) for t in target_tokens]


def per_row_params(sampling_params, caps) -> list:
    """One clone of sampling_params per row with max_tokens = that row's cap; the input is not changed."""
    rows = []
    for cap in caps:
        p = sampling_params.clone()
        p.max_tokens = cap
        rows.append(p)
    return rows
