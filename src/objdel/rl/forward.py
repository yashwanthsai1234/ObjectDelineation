"""The pure parts of the fast GRPO forwards for Qwen3-VL (perf C1, C4, RL-6; from the old train/fast.py, E5).

TRL 1.17 pads the 256 answers of a generation batch to the longest one and runs every forward at that width: the
old-logprob and reference passes (chunks of per_device_train_batch_size rows) and each update micro-batch. All three
go through `GRPOTrainer._get_per_token_logps_and_entropies`, so `rl.trl_adapter.FastForward` overrides that one method
with the functions here:

    per chunk: [pad pad P P P P | A A A . . . . . .]      cut to  [pad P P P P | A A A .]  (a BUCKETS width)
               left pads of the whole chunk ─┘   └─ cut columns come back as 0 (masked, as TRL's own padding)

Widths are rounded up to BUCKETS (8 shapes: the torch.compile recompile limit) so a compiled or graphed layer sees
few shapes. The window is computed from the attention mask with one host read per chunk. The image caches of the
same forwards are in `policy.caches`.
"""

import dataclasses
import functools
import math
import sys
import typing

import numpy as np
import torch

from objdel.policy.caches import RowMeta, rope_bad_rows

BUCKETS = (256, 384, 512, 768, 1024, 1280, 1536, 2048)  # total prompt + completion widths, compiled layers


def window(attention_mask: torch.Tensor, logits_to_keep: int, buckets=BUCKETS) -> tuple[int, int]:
    """(start, keep): the chunk's columns [start, P + keep) hold every real token, where P is the prompt width.
    keep >= the longest real answer; the width P + keep - start is the smallest bucket that fits (or the full width)."""
    W = attention_mask.size(1)
    P = W - logits_to_keep
    first, longest = torch.stack([attention_mask[:, :P].any(0).int().argmax(),
                                  attention_mask[:, P:].sum(1).max()]).tolist()  # fmt: skip
    prompt = P - first
    longest = max(longest, 1)  # keep >= 1: TRL scores labels[:, -keep:]
    width = min(next((b for b in buckets if b >= prompt + longest), W), W)
    keep = min(logits_to_keep, width - prompt)
    return P + keep - width, keep


def row_meta(prompt_ids_list, completion_ids_list, image_ids, G: int, mask_truncated: bool,
             eos_and_pad_ids, sizes=None, vision_ids=(), image_pad=None) -> RowMeta:  # fmt: skip
    """RowMeta from the host lists TRL already holds. A row is truncated when its last id is not in eos_and_pad_ids
    (TRL's own rule, grpo_trainer.py 2604-2606, with TRL's list). image_pad: (image pad id, spatial merge size), or
    None for a text model (pixel rows 0)."""
    if any(not c for c in completion_ids_list):
        raise ValueError("an empty answer: TRL's truncation rule reads its last id")
    ends = set(eos_and_pad_ids)
    clen = [0 if mask_truncated and c[-1] not in ends else len(c) for c in completion_ids_list]
    i32 = functools.partial(np.asarray, dtype=np.int32)
    seen: dict[tuple, int] = {}
    layout = [seen.setdefault((tuple(p), z), len(seen))
              for p, z in zip(prompt_ids_list, sizes or [None] * len(clen), strict=True)]  # fmt: skip
    bad = rope_bad_rows(prompt_ids_list, completion_ids_list, prompt_ids_list, vision_ids)
    pad, merge = image_pad or (None, 0)
    return RowMeta(i32([len(p) for p in prompt_ids_list]), i32(clen), i32(image_ids),
                   i32([i // G for i in range(len(clen))]), i32(layout), i32(bad),
                   i32([list(p).count(pad) * merge**2 for p in prompt_ids_list]))  # fmt: skip


class MetaRow(typing.NamedTuple):
    """One row of a RowMeta: the batch's host_meta list field, which TRL's split and shuffle move with the row."""

    prompt_len: int
    completion_len: int
    image_id: int
    group: int
    layout: int
    rope_bad: int
    pixel_rows: int


def host_rows(meta: RowMeta) -> list[MetaRow]:
    """meta as one MetaRow per row (Python ints)."""
    return [MetaRow(*r) for r in np.stack(dataclasses.astuple(meta), 1).tolist()]


def meta_of(batch: dict) -> RowMeta | None:
    """The RowMeta of a batch's host_meta rows (in their current order), or None without host_meta."""
    rows = batch.get("host_meta")
    return None if rows is None else RowMeta(*(np.asarray(c, dtype=np.int32) for c in zip(*rows, strict=True)))


def image_ids(rows) -> list[int]:
    """Each dataset row's image as a dense id in first-appearance order, by object (OpenOnce shares one object per
    distinct image)."""
    seen: dict[int, int] = {}
    return [seen.setdefault(id(r.get("image")), len(seen)) for r in rows]


def check_meta(batch: dict) -> None:
    """C4 step 5 (--check-meta-every): host_meta's lengths equal the device mask sums and its pixel rows the grid's;
    raises naming field and row."""
    dev = {f: batch[m].sum(1).tolist() for f, m in (("prompt_len", "prompt_mask"),
                                                     ("completion_len", "completion_mask"))}  # fmt: skip
    if batch.get("image_grid_thw") is not None:
        dev["pixel_rows"] = grid_rows(batch["image_grid_thw"], batch["num_images"])
    for field, vals in dev.items():
        for i, (r, d) in enumerate(zip(batch["host_meta"], vals, strict=True)):
            if getattr(r, field) != d:
                raise AssertionError(f"host_meta row {i}: {field} {getattr(r, field)} != device {d}")


def window_from_meta(meta: RowMeta, rows, P: int, logits_to_keep: int, buckets=BUCKETS) -> tuple[int, int]:
    """window() of the chunk `rows` (indices into meta) without a device read; P is the padded prompt width."""
    prompt, longest = int(meta.prompt_len[rows].max()), max(int(meta.completion_len[rows].max()), 1)
    W = P + logits_to_keep
    width = min(next((b for b in buckets if b >= prompt + longest), W), W)
    keep = min(logits_to_keep, width - prompt)
    return P + keep - width, keep


def meta_runs(meta: RowMeta, rows: slice) -> list[tuple[int, int]] | None:
    """prompt_groups() of the chunk `rows` from the host: runs of consecutive rows of one group (one prompt and
    image), or None when no run has two rows."""
    g = meta.group[rows]
    cut = [0, *(np.flatnonzero(g[1:] != g[:-1]) + 1).tolist(), len(g)]
    runs = list(zip(cut[:-1], cut[1:], strict=True))
    return runs if any(b - a > 1 for a, b in runs) else None


def length_order(meta: RowMeta) -> np.ndarray:
    """by_length's order: real prompt + answer tokens, longest first, stable."""
    return np.argsort(-(meta.prompt_len.astype(np.int64) + meta.completion_len), kind="stable")


def group_order(meta: RowMeta, G: int) -> np.ndarray:
    """by_group_length's order: groups of G by their longest member, longest first, stable; rows kept in place."""
    n = (meta.prompt_len.astype(np.int64) + meta.completion_len).reshape(-1, G).max(1)
    return (np.argsort(-n, kind="stable")[:, None] * G + np.arange(G)).ravel()


def grid_rows(grid, nimg) -> list[int]:
    """Each row's pixel rows (its images' t*h*w) from the image grid: one host read."""
    per_image, cut = grid.prod(-1).tolist(), np.cumsum([0, *nimg])
    return [sum(per_image[a:b]) for a, b in zip(cut[:-1], cut[1:], strict=True)]


def prompt_groups(prompt_ids, prompt_mask, pix, image_ids=None) -> list[tuple[int, int]] | None:
    """Runs [a, b) of consecutive rows with the same prompt tokens, padding and pixels, or image ids under
    --encoder-once (exact compares, one host read), or None when no run has two rows (nothing to share)."""
    same = (prompt_ids[1:] == prompt_ids[:-1]).all(1) & (prompt_mask[1:] == prompt_mask[:-1]).all(1)
    if image_ids is not None:
        same &= torch.as_tensor(image_ids[1:] == image_ids[:-1], device=same.device)
    if pix is not None:
        if len({p.shape for p in pix}) == 1:
            px = torch.stack(pix).flatten(1)
            same &= (px[1:] == px[:-1]).all(1)
        else:
            eq = [p.shape == q.shape and torch.equal(p, q) for p, q in zip(pix[1:], pix[:-1], strict=True)]
            same &= torch.tensor(eq, dtype=torch.bool, device=same.device)
    runs, a = [], 0
    for r, eq in enumerate(same.tolist() + [False], 1):
        if not eq:
            runs.append((a, r))
            a = r
    return runs if any(b - a > 1 for a, b in runs) else None


PAIR = 12_000  # linear-layer cost of one token / attention cost of one query-key pair (Qwen3-VL-2B: ~2.8 GFLOP /


def packing_wins(prompt_len, lens, groups) -> bool:
    """Whether `groups` packed rows (prompt once + each answer and its repeated last prompt token) cost less than
    TRL's rows (prompt + the longest answer each), counting n + n^2 / PAIR per row of n tokens."""
    def cost(n):
        return n + n * n / PAIR

    size = len(lens) // groups
    packed = sum(cost(prompt_len + sum(L + 1 for L in lens[g * size:(g + 1) * size])) for g in range(groups))
    return packed < len(lens) * cost(prompt_len + max(lens))


def raw_ratio(logps, sampling, mask, mode: str) -> torch.Tensor:
    """exp of the trainer - vLLM logprob gap per token (`token_*`) or summed per sequence (`sequence_*`, shape
    (B, 1)), before TRL's bounds. Tokens vLLM could not score (NaN) get no correction."""
    diff = torch.nan_to_num((logps - sampling) * mask, nan=0.0)
    if mode.startswith("sequence"):
        diff = diff.sum(-1, keepdim=True)
    return torch.exp(diff)


def bounded(r, mode: str, lo: float | None, hi: float | None) -> torch.Tensor:
    """TRL's bounds on a raw ratio: `*_truncate` clamps to [lo, hi], `*_mask` zeroes ratios outside it. A truncate
    without clip_max would pass an overflowed exp (inf) to the loss, so it is refused."""
    if mode.endswith("truncate"):
        if hi is None:
            raise ValueError(f"{mode} needs vllm_importance_sampling_clip_max (an exp overflow would reach the loss)")
        return torch.clamp(r, min=lo, max=hi)
    lo, hi = -math.inf if lo is None else lo, math.inf if hi is None else hi
    return r.masked_fill((r < lo) | (r > hi), 0.0)


def is_ratio(logps, sampling, mask, mode: str, lo: float | None, hi: float | None) -> torch.Tensor:
    """TRL 1.17's vLLM importance-sampling ratio (grpo_trainer.py 2776-2824): bounded(raw_ratio(...))."""
    return bounded(raw_ratio(logps, sampling, mask, mode), mode, lo, hi)


def clip_counts(raw, r, keep, lo: float | None, hi: float | None) -> torch.Tensor:
    """(zero, hi, lo) counts over the counted entries `keep`, on the device: the bounded ratio r is 0; the raw ratio
    is above hi; the raw ratio is below lo or 0 (an fp32 underflow of exp: trainer << vLLM), each entry once."""
    lo, hi = -math.inf if lo is None else lo, math.inf if hi is None else hi
    return torch.stack([((r == 0) & keep).sum(), ((raw > hi) & keep).sum(), (((raw < lo) | (raw == 0)) & keep).sum()])


OPSM_DELTAS = (0.05, 0.1, 0.5)  # counterfactual OPSM thresholds (B2.1; B2.4 rule 4)

# is_summary's fields: SUMS merge over micro-batches by addition, MAXES by maximum, MINS by minimum. The first 7 sums,
# d_max, r_max and r_min are TRL's sampling metrics and the A6 zero / hi / lo counts.
SUMS = ("d_sum", "d_n", "r_sum", "r_n", "zero", "hi", "lo", "masked", "clamped", "seqs", "seq_gap_sum", "tok",
        "tok_big", "masked_seqs", "masked_len", "kept_len", "modified", "w_sum", "w_sq", "cf_seqmask3",
        "cf_tok_clamp2", "cf_geo", *(f"cf_opsm_{d}" for d in OPSM_DELTAS), "cf_bigtok_unmasked", "cf_lt_0.1",
        "C_sum", "C_sq", "B_sum", "B_sq", *(f"E_{d}_{x}" for d in OPSM_DELTAS for x in ("sum", "sq")), "opsm_masked",
        "opsm_adv_sum")  # fmt: skip

MAXES, MINS = ("d_max", "r_max", "seq_gap_max"), ("r_min",)

SUMMARY = SUMS + MAXES + MINS


def seq_kl(logps, sampling, mask) -> torch.Tensor:
    """TRL's OPSM gap per sequence (grpo_trainer.py get_off_policy_mask): the vLLM-minus-trainer logprob mean over
    the valid tokens (an unscored, NaN token counts as 0 and stays in the count), shape (B,)."""
    m = mask.to(logps.dtype)
    return (torch.nan_to_num(sampling - logps, nan=0.0) * m).sum(-1) / m.sum(-1).clamp(min=1.0)


def is_summary(logps, sampling, mask, ratio, advantages, mode: str, lo: float | None, hi: float | None,
               opsm: float | None = None) -> torch.Tensor:
    """What the IS correction of `mode` does to one micro-batch, and what other modes would do to the same tensors,
    as a float32 vector of SUMMARY fields on the input device (no host read). mask: valid answer tokens (bool, B x
    T); ratio: is_ratio's output; advantages: (B,); opsm: the active off_policy_mask_threshold. Token-level weights
    broadcast a sequence weight to its tokens; every token count is over valid tokens.
    - masked / clamped: ratio units (sequences or tokens) outside [lo, hi] (strict) that a `*_mask` mode zeroes / a
      `*_truncate` mode clamps; masked_seqs, masked_len, kept_len: sequences with a masked unit, and lengths;
    - seq_gap_*: the signed summed trainer-minus-vLLM gap per sequence; tok_big: tokens with |gap| > 0.5;
    - modified: tokens whose weight differs from the raw ratio by more than 1 %; w_*: ESS sums of the weights;
    - cf_*: sequence_mask cap 3 (seqmask3, C), token TIS C = 2 (tok_clamp2, B), the geometric mask [0.99, 1.01],
      TRL's OPSM at each OPSM_DELTAS (drop: advantage < 0 and seq_kl > delta), the big tokens OPSM at 0.1 keeps,
      kept sequences with weight in (0, 0.1); E = token TIS times OPSM;
    - opsm_*: the active OPSM's dropped sequences and the sum of their advantages (0 when opsm is None)."""
    f = lambda x: x.sum().float()  # noqa: E731
    seqmode, lo_, hi_ = mode.startswith("sequence"), -math.inf if lo is None else lo, math.inf if hi is None else hi
    ok = mask & ~torch.isnan(sampling)
    gap = torch.nan_to_num((logps - sampling) * mask, nan=0.0)  # raw_ratio's gap
    raw = raw_ratio(logps, sampling, mask, mode)
    keep = torch.ones_like(ratio, dtype=torch.bool) if seqmode else mask
    out = ((raw > hi_) | (raw < lo_)) & keep
    none = torch.zeros_like(out)
    masked, clamped = (out, none) if mode.endswith("mask") else (none, out)
    seq_masked = masked.any(-1)
    lens = mask.sum(-1).float()
    seq_gap = gap.sum(-1)
    w = ratio.expand_as(gap) * mask
    kl = seq_kl(logps, sampling, mask)
    neg = advantages.reshape(-1) < 0
    drop = {d: neg & (kl > d) for d in OPSM_DELTAS}
    big = (gap.abs() > 0.5) & mask
    tis = torch.exp(gap).clamp(max=2.0) * mask
    seq_raw = torch.exp(seq_gap)
    w_c = torch.where(seq_raw > 3.0, 0.0, seq_raw)[:, None] * mask
    geo = torch.exp(seq_gap / lens.clamp(min=1))  # the geometric mean ratio
    sw = w.sum(-1) / lens.clamp(min=1)  # a sequence's mean weight
    act = neg & (kl > opsm) if opsm is not None else torch.zeros_like(neg)
    absd = (logps - sampling).abs()
    ess = lambda x: (f(x), f(x * x))  # noqa: E731
    sums = [f(absd.masked_fill(~ok, 0)), f(ok), f(ratio.masked_fill(~keep, 0)), f(keep),
            *clip_counts(raw, ratio, keep, lo, hi).float(), f(masked), f(clamped),
            f(torch.ones_like(seq_gap)), f(seq_gap), f(mask), f(big), f(seq_masked),
            f(lens * seq_masked), f(lens * ~seq_masked), f(((w / raw.expand_as(w) - 1).abs() > 0.01) & mask),
            *ess(w), f(seq_raw > 3.0), f((torch.exp(gap) > 2.0) & mask), f((geo < 0.99) | (geo > 1.01)),
            *(f(drop[x]) for x in OPSM_DELTAS), f(big & ~drop[0.1][:, None]), f((sw > 0) & (sw < 0.1)), *ess(w_c),
            *ess(tis), *(v for x in OPSM_DELTAS for v in ess(tis * ~drop[x][:, None])), f(act),
            f(advantages.reshape(-1) * act)]  # fmt: skip
    return torch.cat([torch.stack(sums), torch.stack([
        absd.masked_fill(~ok, -math.inf).amax(), ratio.masked_fill(~keep, -math.inf).amax(), seq_gap.amax(),
        ratio.masked_fill(~keep, math.inf).amin()]).float()])  # fmt: skip


def merge_summary(old, new) -> torch.Tensor:
    """Two is_summary vectors of one generation batch as one (on the device)."""
    if old is None:
        return new
    a, b = len(SUMS), len(SUMS) + len(MAXES)
    return torch.cat([old[:a] + new[:a], torch.maximum(old[a:b], new[a:b]), torch.minimum(old[b:], new[b:])])


SEQ_KL_Q = (0.5, 0.99, 0.995)  # logged seq_kl quantiles: sampling/seq_kl/{q50,q99,q99.5}


def summary_metrics(s: dict, opsm: float | None) -> dict:
    """The logged `sampling/*` values from a merged is_summary (field name -> float); a ratio of an empty count is
    0, never NaN. The opsm_* keys only when OPSM is active."""
    div = lambda a, b: a / b if b else 0.0  # noqa: E731
    ess = lambda k: div(s[f"{k}_sum"] ** 2, s["tok"] * s[f"{k}_sq"])  # noqa: E731
    out = {"sampling_logp_difference/mean": div(s["d_sum"], s["d_n"]),
           "sampling_logp_difference/max": s["d_max"] if s["d_n"] else 0.0,
           "importance_sampling_ratio/min": s["r_min"] if s["r_n"] else 0.0,
           "importance_sampling_ratio/mean": div(s["r_sum"], s["r_n"]),
           "importance_sampling_ratio/max": s["r_max"] if s["r_n"] else 0.0,
           **{f"importance_sampling_ratio/{k}_share": div(s[k], s["r_n"]) for k in ("zero", "hi", "lo")},
           "is_masked_share": div(s["masked"], s["r_n"]), "is_clamped_share": div(s["clamped"], s["r_n"]),
           "seq_logdiff_mean": div(s["seq_gap_sum"], s["seqs"]),
           "seq_logdiff_max": s["seq_gap_max"] if s["seqs"] else 0.0,
           "tok_absdiff_gt_0.5_share": div(s["tok_big"], s["tok"]),
           "masked_len_mean": div(s["masked_len"], s["masked_seqs"]),
           "kept_len_mean": div(s["kept_len"], s["seqs"] - s["masked_seqs"]),
           "modified_token_share": div(s["modified"], s["tok"]), "ess": ess("w"),
           "cf/seqmask3_share": div(s["cf_seqmask3"], s["seqs"]),
           "cf/tok_clamp2_share": div(s["cf_tok_clamp2"], s["tok"]),
           "cf/geo_0.99_1.01_share": div(s["cf_geo"], s["seqs"]),
           **{f"cf/opsm_share@{d}": div(s[f"cf_opsm_{d}"], s["seqs"]) for d in OPSM_DELTAS},
           "cf/opsm0.1_bigtok_unmasked_share": div(s["cf_bigtok_unmasked"], s["tok"]),
           "cf/seq_ratio_lt_0.1_share": div(s["cf_lt_0.1"], s["seqs"]), "cf/ess_C": ess("C"), "cf/ess_B": ess("B"),
           **{f"cf/ess_E@{d}": ess(f"E_{d}") for d in OPSM_DELTAS}}  # fmt: skip
    if opsm is not None:
        out |= {"opsm_masked_share": div(s["opsm_masked"], s["seqs"]), "opsm_masked_adv_sum": s["opsm_adv_sum"]}
    return out


def write_seq_gaps(out_dir: str, values: list[float], rows: int) -> None:
    """Appends one logged step's (rows x 3) per-sequence values to <out_dir>/seq_gaps.f32 (raw little-endian float32,
    steps in order) and rewrites seq_gaps.json; B2.4 rules 3 and 4 read the file (B1 runs once)."""
    import json
    from array import array
    from pathlib import Path

    d = Path(out_dir)
    d.mkdir(parents=True, exist_ok=True)
    a = array("f", values)
    if sys.byteorder != "little":
        a.byteswap()
    with open(d / "seq_gaps.f32", "ab") as f:
        a.tofile(f)
    steps = (d / "seq_gaps.f32").stat().st_size // (rows * 3 * 4)
    (d / "seq_gaps.json").write_text(json.dumps({"rows_per_step": rows, "cols": ["seq_kl", "advantage",
                                                 "valid_tokens"], "steps": steps}) + "\n")  # fmt: skip


def by_length(batch: dict) -> dict:
    """The generation batch's rows longest first (prompt + answer real tokens; stable), every per-row field moved
    with them, in place of TRL's shuffle_sequence_dict (same field rules: None and 0-d tensors stay, tensors are
    indexed, lists reordered). One host read for the order (none with host_meta, C4)."""
    meta = meta_of(batch)
    if meta is not None:
        order = torch.from_numpy(length_order(meta))
    else:
        n = batch["prompt_mask"].sum(1) + batch["completion_mask"].sum(1)
        order = torch.sort(n, descending=True, stable=True).indices
    idx = order.tolist()

    def move(v):
        if v is None or isinstance(v, torch.Tensor) and v.ndim == 0:
            return v
        return v[order.to(v.device)] if isinstance(v, torch.Tensor) else [v[i] for i in idx]

    return {k: move(v) for k, v in batch.items()}


def by_group_length(batch: dict, G: int) -> dict:
    """by_length with each group of G consecutive rows (one prompt's answers) kept together: groups longest member
    first (stable), so update micro-batches of whole groups can share their prompt (FastForward.shared_prefix)."""
    meta = meta_of(batch)
    if meta is not None:  # C4: no device read
        order = torch.from_numpy(group_order(meta, G))
    else:
        n = batch["prompt_mask"].sum(1) + batch["completion_mask"].sum(1)
        groups = torch.sort(n.view(-1, G).max(1).values, descending=True, stable=True).indices
        order = (groups[:, None] * G + torch.arange(G, device=groups.device)).flatten()
    idx = order.tolist()

    def move(v):
        if v is None or isinstance(v, torch.Tensor) and v.ndim == 0:
            return v
        return v[order.to(v.device)] if isinstance(v, torch.Tensor) else [v[i] for i in idx]

    return {k: move(v) for k, v in batch.items()}
