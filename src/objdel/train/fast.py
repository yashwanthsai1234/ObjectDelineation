"""Fast GRPO forwards for Qwen3-VL (perf C1, RL-6, C6; experiments/phase2/perf/implementation.md).

TRL 1.17 pads the 256 answers of a generation batch to the longest one and runs every forward at that width: the
old-logprob and reference passes (chunks of per_device_train_batch_size rows) and each update micro-batch. All three
go through `GRPOTrainer._get_per_token_logps_and_entropies`, so `FastForward` overrides that one method:

    per chunk: [pad pad P P P P | A A A . . . . . .]      cut to  [pad P P P P | A A A .]  (a BUCKETS width)
               left pads of the whole chunk ─┘   └─ cut columns come back as 0 (masked, as TRL's own padding)

Widths are rounded up to BUCKETS (8 shapes: the torch.compile recompile limit) so a compiled or graphed layer sees
few shapes. The window is computed from the attention mask with one host read per chunk.

`VisionCache` (C6): the G = 8 answers of a prompt share its image, so one generation batch holds 32 images 256 times,
and the frozen vision tower ran on every copy in every pass (old, reference, update). The cache encodes each distinct
image once per generation batch. Valid because LoRA leaves the tower and projector out (configs.lora, D6).
`EncoderCache` (C3, --encoder-once) takes over its tower: the batch's distinct images in one call, no pixel values in
any forward.
"""

import contextlib
import dataclasses
import functools
import math
import sys
import time
import typing

import numpy as np
import torch
import torch.nn.functional as F
import trl.trainer.grpo_trainer as trl_grpo
from trl import GRPOTrainer

from objdel.train import nvtx, runaway, trl_adapter, vllm_requests
from objdel.train.dataset import OpenOnce

BUCKETS = (256, 384, 512, 768, 1024, 1280, 1536, 2048)  # total prompt + completion widths, compiled layers
FINE = tuple(range(64, 2049, 64))  # uncompiled: round to 64 only (BUCKETS would pad a 430-wide chunk to 512)
ROW_KEYS = ("mm_token_type_ids", "token_type_ids")  # per-token, prompt + completion wide: cut like input_ids
OTHER_VISION = ("pixel_attention_mask", "spatial_shapes", "num_tiles", "image_sizes", "image_position_ids")


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



@dataclasses.dataclass(frozen=True)
class RowMeta:
    """Host facts of a generation batch's rows (C4), int32 arrays of shape (B,): real prompt tokens, real answer
    tokens (0 for a truncated row under mask_truncated_completions, as TRL zeroes its completion_mask), the row's
    distinct image (first-appearance index), its group of G answers, its prompt layout (dense id of prompt ids and
    image size: one rope template each, C2) and the host's rope routing (1: the HF path, rope_bad_rows)."""

    prompt_len: np.ndarray
    completion_len: np.ndarray
    image_id: np.ndarray
    group: np.ndarray
    layout: np.ndarray
    rope_bad: np.ndarray


def row_meta(prompt_ids_list, completion_ids_list, image_ids, G: int, mask_truncated: bool,
             eos_and_pad_ids, sizes=None, vision_ids=()) -> RowMeta:  # fmt: skip
    """RowMeta from the host lists TRL already holds. A row is truncated when its last id is not in eos_and_pad_ids
    (TRL's own rule, grpo_trainer.py 2604-2606, with TRL's list)."""
    if any(not c for c in completion_ids_list):
        raise ValueError("an empty answer: TRL's truncation rule reads its last id")
    ends = set(eos_and_pad_ids)
    clen = [0 if mask_truncated and c[-1] not in ends else len(c) for c in completion_ids_list]
    i32 = functools.partial(np.asarray, dtype=np.int32)
    seen: dict[tuple, int] = {}
    layout = [seen.setdefault((tuple(p), z), len(seen))
              for p, z in zip(prompt_ids_list, sizes or [None] * len(clen), strict=True)]  # fmt: skip
    bad = rope_bad_rows(prompt_ids_list, completion_ids_list, prompt_ids_list, vision_ids)
    return RowMeta(i32([len(p) for p in prompt_ids_list]), i32(clen), i32(image_ids),
                   i32([i // G for i in range(len(clen))]), i32(layout), i32(bad))  # fmt: skip


class MetaRow(typing.NamedTuple):
    """One row of a RowMeta: the batch's host_meta list field, which TRL's split and shuffle move with the row."""

    prompt_len: int
    completion_len: int
    image_id: int
    group: int
    layout: int
    rope_bad: int


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
    """C4 step 5 (--check-meta-every): host_meta's lengths equal the device mask sums; raises naming field and row."""
    for field, mask in (("prompt_len", "prompt_mask"), ("completion_len", "completion_mask")):
        for i, (r, d) in enumerate(zip(batch["host_meta"], batch[mask].sum(1).tolist(), strict=True)):
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

def per_row(pv, rows, nimg) -> list | None:
    """Each row's pixel rows (its images' slices of the chunk's pixel_values), or None without images."""
    if pv is None:
        return None
    out, i, at = [], 0, 0
    for k in nimg:
        r = sum(rows[i:i + k])
        out.append(pv[at:at + r])
        i, at = i + k, at + r
    return out


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
#                ~0.23 MFLOP); packed attention is dense over the packed row (SDPA's masked kernel skips no tiles)


def packing_wins(prompt_len, lens, groups) -> bool:
    """Whether `groups` packed rows (prompt once + each answer and its repeated last prompt token) cost less than
    TRL's rows (prompt + the longest answer each), counting n + n^2 / PAIR per row of n tokens."""
    def cost(n):
        return n + n * n / PAIR

    size = len(lens) // groups
    packed = sum(cost(prompt_len + sum(L + 1 for L in lens[g * size:(g + 1) * size])) for g in range(groups))
    return packed < len(lens) * cost(prompt_len + max(lens))


class FastForward:
    """GRPOTrainer mixin: each forward chunk runs on its own window (see the module docstring)."""

    buckets = FINE  # make_trainer sets BUCKETS when it compiles the text layers
    shared_prefix = False  # pack each group's shared prompt once per chunk (packed_chunk); needs the RopeCache

    def _get_per_token_logps_and_entropies(self, model, input_ids, attention_mask, logits_to_keep, batch_size=None,
                                           compute_entropy=False, compute_aux_loss=False, **kw):  # fmt: skip
        meta = kw.pop("host_meta", None)  # C4: the micro-batch's RowMeta (FastGRPO._compute_loss), None elsewhere
        if any(kw.get(k) is not None for k in OTHER_VISION) or (
                kw.get("pixel_values") is not None and kw.get("image_grid_thw") is None):
            return super()._get_per_token_logps_and_entropies(  # not a grid (Qwen-VL style) model: TRL's path
                model, input_ids, attention_mask, logits_to_keep, batch_size, compute_entropy, compute_aux_loss, **kw)
        n = input_ids.size(0)
        size = batch_size or n
        pv, grid, nimg = kw.pop("pixel_values", None), kw.pop("image_grid_thw", None), kw.pop("num_images", None)
        enc = getattr(self, "encoder", None) if grid is not None else None
        if enc is not None:  # C3: each row's features from the generation batch's one tower call, never pixels
            pv, ids = None, enc.check_rows(meta.image_id if meta is not None else self.row_image_ids, nimg, model)
        rows = grid.prod(-1).tolist() if pv is not None else None  # pixel rows per image (one host read)
        P = input_ids.size(1) - logits_to_keep
        rope = getattr(self, "rope", None)
        pos = None
        if rope is not None and grid is not None and kw.get("mm_token_type_ids") is not None:
            if meta is None:
                pos = rope.positions(input_ids, attention_mask, kw["mm_token_type_ids"], grid, nimg)
            else:  # C2: templates per prompt layout, the host's bad rows on HF's path
                pos, nbad = rope.template_positions(input_ids, attention_mask, kw["mm_token_type_ids"], grid, nimg, P,
                                                    meta, getattr(self, "rope_check", False))  # fmt: skip
                self._metrics["train" if self.model.training else "eval"]["rope/fallback_rows"].append(nbad)
        lps, ents, auxs = [], [], []
        img = row = 0
        for s in range(0, n, size):
            e = min(s + size, n)
            if self.shared_prefix and pos is not None and not compute_aux_loss and rope.hooks(model):
                k = sum(nimg[s:e])
                r = 0 if rows is None else sum(rows[img:img + k])
                pix = None if enc is not None else per_row(pv[row:row + r], rows[img:img + k], nimg[s:e])
                chunk_ids = None if enc is None else ids[s:e]
                if meta is None:
                    runs = prompt_groups(input_ids[s:e, :P], attention_mask[s:e, :P], pix, chunk_ids)
                    lengths = torch.cat([attention_mask[s:e, :P].sum(1), attention_mask[s:e, P:].sum(1)]).tolist()
                    plen, lens = lengths[:e - s], lengths[e - s:]
                else:
                    runs = meta_runs(meta, slice(s, e))
                    plen, lens = meta.prompt_len[s:e].tolist(), meta.completion_len[s:e].tolist()
                if runs is not None and packing_wins(sum(plen) / len(plen), lens, len(runs)):
                    smask = kw["sampling_mask"][s:e] if kw.get("sampling_mask") is not None else None
                    lp, ent = self.packed_chunk(model, input_ids[s:e], P, plen, lens, pos[:, s:e],
                                                kw["mm_token_type_ids"][s:e], runs, pix, grid[img:img + k],
                                                nimg[s:e], smask, compute_entropy, chunk_ids)  # fmt: skip
                    img, row = img + k, row + r
                    lps.append(lp)
                    if compute_entropy:
                        ents.append(ent)
                    continue
            start, keep = (window(attention_mask[s:e], logits_to_keep, self.buckets) if meta is None
                           else window_from_meta(meta, slice(s, e), P, logits_to_keep, self.buckets))
            cols = slice(start, P + keep)
            sub = {k: kw[k][s:e, cols] for k in ROW_KEYS if kw.get(k) is not None}
            if kw.get("sampling_mask") is not None:
                sub["sampling_mask"] = kw["sampling_mask"][s:e, :keep]
            if grid is not None:
                k = sum(nimg[s:e])
                if enc is not None:  # taken by the base's forward (EncoderCache.inject): TRL passes no such kwarg
                    enc.pending = (enc.outputs(ids[s:e]), grid[img:img + k])
                else:
                    r = sum(rows[img:img + k])
                    sub |= {"pixel_values": pv[row:row + r], "image_grid_thw": grid[img:img + k],
                            "num_images": nimg[s:e]}  # fmt: skip
                    row += r
                img += k
            if pos is not None:
                rope.pending = pos[:, s:e, cols]  # the chunk's slice, taken by the model's forward (RopeCache.inject)
            try:
                lp, ent, aux = super()._get_per_token_logps_and_entropies(
                    model, input_ids[s:e, cols], attention_mask[s:e, cols], keep, None, compute_entropy,
                    compute_aux_loss, **sub)  # fmt: skip
            finally:
                if pos is not None:
                    rope.pending = None
                if enc is not None:
                    enc.pending = None
            lps.append(F.pad(lp, (0, logits_to_keep - keep)))
            if compute_entropy:
                ents.append(F.pad(ent, (0, logits_to_keep - keep)))
            if compute_aux_loss:
                auxs.append(aux)
        return (torch.cat(lps), torch.cat(ents) if compute_entropy else None,
                torch.stack(auxs).mean() if compute_aux_loss else None)  # fmt: skip

    def packed_chunk(self, model, ids, P, plen, lens, pos, mm, runs, pix, grid, nimg, sampling_mask, compute_entropy,
                     image_ids=None):  # fmt: skip
        """One packed row per group of rows sharing a prompt and image (prompt_groups' runs):
            [prompt minus its last token] + for each answer: [last prompt token, answer tokens]
        The last prompt token is repeated at the head of each answer block with its own position, so its hidden
        state (same token, position and visible keys) predicts each answer's first token, and TRL's fused head
        scores every block with its usual one-token shift. A 4-D mask from segment numbers (0 = prompt, k = answer
        block k, -1 = padding) lets the prompt attend causally and each block attend to the prompt and causally to
        itself, never to another answer. Positions are the rows' own (RopeCache), sliced. `plen` and `lens` are the
        rows' real prompt and answer lengths (host lists). Returns (rows, answer width) log-probs and entropies laid
        out as TRL's, 0 off the answers. `image_ids` (--encoder-once): the rows' image ids, each group's features from
        EncoderCache instead of `pix`."""
        dev, Wr, T = ids.device, ids.size(1), ids.size(1) - P
        srcs, segs, gi, pi, ri, ci = [], [], [], [], [], []
        for g, (a, b) in enumerate(runs):
            f, lp0 = P - plen[a], plen[a] - 1
            src, seg = [*range(a * Wr + f, a * Wr + P - 1)], [0] * lp0  # flat (row, column) of each packed token
            off = lp0
            for k, r in enumerate(range(a, b), 1):
                L = lens[r]
                src += range(r * Wr + P - 1, r * Wr + P + L)
                seg += [k] * (L + 1)
                gi += [g] * L  # the fused head's outputs are shifted: index q - 1 scores token q
                pi += range(off, off + L)
                ri += [r] * L
                ci += range(L)
                off += L + 1
            srcs.append(src), segs.append(seg)
        W, c = max(map(len, segs)), len(runs)
        host = np.concatenate([np.asarray([x + [-1] * (W - len(x)) for x in segs]).ravel(),
                               np.asarray([x + [0] * (W - len(x)) for x in srcs]).ravel(),
                               np.asarray([gi, pi, ri, ci]).ravel()]).astype(np.int64)  # fmt: skip
        stage = torch.empty(host.size, dtype=torch.int64, pin_memory=dev.type == "cuda")  # fresh per chunk (C4 step 4)
        stage.numpy()[:] = host
        flat = stage.to(dev, non_blocking=True)  # one copy; the caching host allocator keeps `stage` until it ends
        seg, src, idx = flat[:c * W].view(c, W), flat[c * W:2 * c * W], flat[2 * c * W:].view(4, -1)
        pad = seg < 0
        tok = torch.index_select(ids.reshape(-1), 0, src).view(c, W).masked_fill(pad, 0)
        mmt = torch.index_select(mm.reshape(-1), 0, src).view(c, W).masked_fill(pad, 0)
        posp = torch.index_select(pos.reshape(pos.size(0), -1), 1, src).view(-1, c, W).masked_fill(pad, 0)
        lab = torch.full((c, W), -100, dtype=ids.dtype, device=dev)
        lab[idx[0], idx[1] + 1] = tok[idx[0], idx[1] + 1]  # each answer token labels itself; prompt and heads -100
        q, k_ = seg[:, :, None], seg[:, None, :]
        causal = torch.ones(W, W, dtype=torch.bool, device=dev).tril()
        allow = causal & (k_ >= 0) & (q >= 0) & ((k_ == 0) | (k_ == q))
        allow |= torch.eye(W, dtype=torch.bool, device=dev)  # padding rows see only themselves
        impl = getattr(model.config, "_attn_implementation", "sdpa")
        if impl == "eager":
            m4 = torch.zeros(c, 1, W, W, dtype=model.dtype, device=dev)
            m4 = m4.masked_fill(~allow[:, None], torch.finfo(model.dtype).min)
        elif impl == "sdpa":
            m4 = allow[:, None]
        else:
            raise NotImplementedError(f"shared_prefix needs sdpa or eager attention, not {impl}")
        inputs = {"input_ids": tok, "attention_mask": m4, "mm_token_type_ids": mmt, "labels": lab}
        if pix is not None or image_ids is not None:
            heads = [a for a, _ in runs]
            img0 = [sum(nimg[:a]) for a in heads]
            inputs["image_grid_thw"] = torch.cat([grid[i:i + nimg[a]] for i, a in zip(img0, heads, strict=True)])
            if image_ids is None:
                inputs["pixel_values"] = torch.cat([pix[a] for a in heads])
            else:
                inputs["mm_encoder_outputs"] = {"image": self.encoder.outputs(image_ids[heads])}
        if sampling_mask is not None:
            sup = sampling_mask.new_full((c, W, sampling_mask.size(-1)), -1)
            sup[idx[0], idx[1] + 1] = sampling_mask[idx[2], idx[3]]  # aligned with the tokens, shifted like labels
            inputs["support_ids"] = sup
        mode = "train" if self.model.training else "eval"  # packed tokens / the rows' own (prompt + answer) tokens
        self._metrics[mode]["shared/token_ratio"].append(sum(map(len, segs)) / (sum(plen) + sum(lens)))
        self.rope.pending = posp
        try:
            with self.accelerator.autocast():
                out = model(**inputs, fused_lm_head=True)
        finally:
            self.rope.pending = None

        def scatter(x):
            return x.new_zeros(ids.size(0), T).index_put((idx[2], idx[3]), x[idx[0], idx[1]])

        return scatter(out.log_probs), (scatter(out.entropy) if compute_entropy else None)


INT = {4: torch.int32, 2: torch.int16}  # float bits as integers: an exact, order-free hash


class VisionCache:
    """Wraps a Qwen-VL vision tower's forward (in place: parameter names, and so the vLLM weight sync, unchanged).
    Images are keyed by an exact integer hash of their pixel rows and their grid; `clear()` per generation batch."""

    def __init__(self, visual):
        self.visual, self.encode, self.store = visual, visual.forward, {}
        visual.forward = self.forward

    @classmethod
    def install(cls, model):
        """The cache on the model's vision tower (a PEFT wrapper is searched through), or None without one."""
        visual = next((m for n, m in model.named_modules() if n.endswith("visual") and hasattr(m, "merger")), None)
        return cls(visual) if visual is not None else None

    def clear(self):
        """Empty the cache."""
        self.store.clear()

    def keys(self, pixel_values, grid_thw) -> list[tuple]:
        """Exact integer hash and (t, h, w) grid of each image: one tuple per image, with one host read per call."""
        rows = grid_thw.prod(-1)
        image = torch.repeat_interleave(torch.arange(len(rows), device=rows.device), rows)
        first = torch.cumsum(rows, 0) - rows
        g = torch.Generator().manual_seed(0)
        w = torch.randint(1, 2**31, (pixel_values.size(1),), generator=g).to(pixel_values.device)
        bits = pixel_values.contiguous().view(INT[pixel_values.element_size()])
        per_row = torch.cat([(c.long() * w).sum(1) for c in bits.split(8192)])
        place = torch.arange(len(per_row), device=rows.device) - first[image] + 1
        h = torch.zeros(len(rows), dtype=torch.long, device=rows.device).index_add_(0, image, per_row * place)
        return list(zip(h.tolist(), map(tuple, grid_thw.tolist()), strict=True))  # one host read per call

    def forward(self, pixel_values, grid_thw=None, **kw):
        """Vision tower forward that encodes only images missing from the cache and returns the pooled and deepstack
        features of all, in order.
        """
        with nvtx.range("vision_cache"):
            keys = self.keys(pixel_values, grid_thw)
            rows = grid_thw.prod(-1).tolist()
            starts = [sum(rows[:i]) for i in range(len(rows))]
            new = list({k: i for i, k in reversed(list(enumerate(keys))) if k not in self.store}.values())
            if new:
                px = torch.cat([pixel_values[starts[i]:starts[i] + rows[i]] for i in new])
                with torch.no_grad():
                    out = self.encode(px, grid_thw=grid_thw[new], **kw)
                self.kind = type(out)
                merged = [rows[i] // self.visual.spatial_merge_size**2 for i in new]
                embeds = out.pooler_output.split(merged)
                levels = [f.split(merged) for f in out.deepstack_features]
                for j, i in enumerate(new):
                    self.store[keys[i]] = (embeds[j], [lv[j] for lv in levels])
            hits = [self.store[k] for k in keys]
            return self.kind(pooler_output=torch.cat([e for e, _ in hits]),
                             deepstack_features=[torch.cat([d[j] for _, d in hits]) for j in range(len(hits[0][1]))])


class RopeCache:
    """Qwen3-VL's 3-D rope positions per row (C5). `get_rope_index` reads the host several times per row; TRL runs
    it in every forward chunk of the old, reference and update passes. Here each row's positions are computed once per
    generation batch at the full width (rows found by their exact ids, mask and image grids: the update's micro-batches
    are shuffled subsets) and each chunk's slice goes to the model through a forward pre-hook. Positions depend only
    on a row's real tokens, so a column slice equals the trimmed chunk's own positions on every real token."""

    def __init__(self, base):
        self.base, self.store, self.pending = base, {}, None
        base.register_forward_pre_hook(self.inject, with_kwargs=True)

    @classmethod
    def install(cls, model):
        """The cache on the model's Qwen-VL base (a PEFT wrapper is searched through), or None without one."""
        base = next((m for _, m in model.named_modules()
                     if hasattr(m, "get_rope_index") and hasattr(m, "compute_3d_position_ids")), None)  # fmt: skip
        return cls(base) if base is not None else None

    def clear(self):
        """Empty the cache."""
        self.store.clear()

    def inject(self, module, args, kwargs):
        """Forward pre-hook: give the pending rope positions as `position_ids` when the call has none."""
        if self.pending is not None and kwargs.get("position_ids") is None:
            kwargs["position_ids"] = self.pending
        return args, kwargs

    @staticmethod
    def keys(input_ids, attention_mask, grid, nimg) -> list[tuple]:
        """Exact row keys: the row's ids and mask bytes and its images' (t, h, w). The prompt text is fixed per
        dataset, so ids alone carry only the pad count t*h*w/4, not h and w (one host read per tensor)."""
        ids, mask, g = input_ids.cpu().numpy(), attention_mask.cpu().numpy(), grid.tolist()
        first = [sum(nimg[:i]) for i in range(len(nimg))]
        return [(ids[i].tobytes(), mask[i].tobytes(), tuple(map(tuple, g[first[i]:first[i] + nimg[i]])))
                for i in range(len(nimg))]  # fmt: skip

    def hooks(self, model) -> bool:
        """Whether `model` contains the hooked base (a separate reference model does not: positions not injected)."""
        seen = self.__dict__.setdefault("hooked", {})
        if id(model) not in seen:
            seen[id(model)] = any(m is self.base for m in model.modules())
        return seen[id(model)]

    def positions(self, input_ids, attention_mask, mm_token_type_ids, grid, nimg) -> torch.Tensor:
        """(3, rows, width) positions of these rows, computing only rows not seen in this generation batch."""
        with nvtx.range("rope_cache"):
            keys = self.keys(input_ids, attention_mask, grid, nimg)
            first = [sum(nimg[:i]) for i in range(len(nimg))]
            miss = list({k: i for i, k in enumerate(keys) if k not in self.store}.values())
            if miss:  # on host copies: get_rope_index's per-row host reads are then no GPU syncs; one copy back
                g = torch.cat([grid[first[i]:first[i] + nimg[i]] for i in miss]).cpu()
                pos, _ = self.base.get_rope_index(
                    input_ids[miss].cpu(), mm_token_type_ids=mm_token_type_ids[miss].cpu(), image_grid_thw=g,
                    attention_mask=attention_mask[miss].cpu())  # fmt: skip
                pos = pos.to(input_ids.device)
                for j, i in enumerate(miss):
                    self.store[keys[i]] = pos[:, j]
            return torch.stack([self.store[k] for k in keys], dim=1)


    def template_positions(self, input_ids, attention_mask, mm_token_type_ids, grid, nimg, P: int, meta: RowMeta,
                           check: bool) -> tuple[torch.Tensor, int]:  # fmt: skip
        """C2: (positions (3, rows, width), bad row count). One rope_template per new prompt layout (built from
        its first row's real prompt), rope_positions for the rows meta marks good, HF's path (positions) for the
        rest. With `check` (--check-rope-every) every row is compared with HF and the device `bad` with the host."""
        with nvtx.range("rope_template"):
            first = np.cumsum([0, *nimg])
            bad = meta.rope_bad.astype(bool)
            out = torch.zeros(3, *input_ids.shape, dtype=torch.long, device=input_ids.device)
            for lay in np.unique(meta.layout[~bad]):
                rows = np.flatnonzero((meta.layout == lay) & ~bad)
                key = ("layout", int(lay))
                if key not in self.store:
                    i, f = rows[0], P - int(meta.prompt_len[rows[0]])
                    self.store[key] = rope_template(self.base, input_ids[i, f:P], mm_token_type_ids[i, f:P],
                                                    grid[first[i]:first[i + 1]])  # fmt: skip
                r = torch.as_tensor(rows, device=input_ids.device)
                out[:, r], dev_bad = rope_positions(self.store[key], input_ids[r], attention_mask[r],
                                                    mm_token_type_ids[r])  # fmt: skip
                if check and dev_bad.any():
                    raise AssertionError(f"rope rows {rows[dev_bad.cpu().numpy()].tolist()} bad on the device, "
                                         "good on the host")  # fmt: skip
            rows = np.arange(len(nimg)) if check else np.flatnonzero(bad)
            if rows.size:
                r = torch.as_tensor(rows, device=input_ids.device)
                g = torch.cat([grid[first[i]:first[i + 1]] for i in rows])
                hf = self.positions(input_ids[r], attention_mask[r], mm_token_type_ids[r], g, [nimg[i] for i in rows])
                if check:
                    real, good = attention_mask.bool(), torch.as_tensor(~bad, device=input_ids.device)
                    if not torch.equal(torch.where(real, out, 0)[:, good], torch.where(real, hf, 0)[:, good]):
                        raise AssertionError("rope positions differ from HF")
                    hf = hf[:, torch.as_tensor(bad, device=input_ids.device)]
                out[:, torch.as_tensor(np.flatnonzero(bad), device=input_ids.device)] = hf
            return out, int(bad.sum())


class EncoderMismatch(RuntimeError):
    """A vision tower call whose pixel rows are not its grid's (C3: e.g. 32 images of pixels cut by a 256-row grid)."""


class EncoderCache:
    """C3 (--encoder-once): `encode` runs the generation batch's distinct images through the frozen tower in one call
    (FastGRPO._generate, after vLLM returns); `outputs(image_ids)` gives each row's features as Qwen3VLModel.forward
    takes them in `mm_encoder_outputs`. TRL's forward passes no such kwarg, so `pending` reaches the base through a
    pre-hook (like RopeCache's positions). Takes over VisionCache's tower wrapper: any other tower call passes the
    shape guard and runs uncached."""

    hooks = RopeCache.hooks  # whether a model holds the hooked base

    def __init__(self, vision: VisionCache, base):
        if any(p.requires_grad for p in vision.visual.parameters()):
            raise ValueError("--encoder-once needs a frozen vision tower: its features are computed once a batch")
        self.visual, self.tower, self.base = vision.visual, vision.encode, base
        self.feats = self.pending = None
        vision.visual.forward = self.guard
        base.register_forward_pre_hook(self.inject, with_kwargs=True)

    def clear(self):
        """Drop the previous generation batch's features."""
        self.feats = None

    def guard(self, pixel_values, grid_thw=None, **kw):
        """The tower, refusing pixel rows that are not the grid's (a TRL change that cuts the per-image pixels by
        the per-row grid)."""
        if pixel_values.size(0) != int(grid_thw.prod(-1).sum()):
            raise EncoderMismatch(f"{pixel_values.size(0)} pixel rows for a grid of {len(grid_thw)} images")
        return self.tower(pixel_values, grid_thw=grid_thw, **kw)

    def encode(self, patches, grid) -> None:
        """The distinct images' pixel rows and (t, h, w) grid rows (image id order) in one no-grad tower call."""
        with nvtx.range("encoder_once"), torch.no_grad():
            out = self.tower(patches.to(self.visual.device, self.visual.dtype), grid_thw=grid.to(self.visual.device),
                             return_dict=True)  # fmt: skip
        merged = (grid.prod(-1) // self.visual.spatial_merge_size**2).tolist()
        self.kind = type(out)
        self.feats = out.pooler_output.split(merged), [f.split(merged) for f in out.deepstack_features]

    def check_rows(self, image_ids, nimg, model):
        """The forward's per-row image ids, refusing a missing or misaligned field, rows of other than one image and
        a model without the hooked base (a separate reference model would get no image features)."""
        if image_ids is None:
            raise RuntimeError("--encoder-once: no image ids for these rows")
        if len(image_ids) != len(nimg) or any(k != 1 for k in nimg):
            raise ValueError(f"--encoder-once: {len(image_ids)} image ids for {len(nimg)} rows of one image each")
        if not self.hooks(model):
            raise RuntimeError("--encoder-once: the model does not hold the hooked Qwen-VL base")
        return image_ids

    def outputs(self, image_ids):
        """The rows' features in row order: a tuple per row (pooler_output) and per DeepStack level."""
        if self.feats is None:
            raise RuntimeError("--encoder-once: no image encoded in this generation batch")
        embeds, levels = self.feats
        bad = [int(i) for i in image_ids if not 0 <= i < len(embeds)]
        if bad:
            raise KeyError(f"image ids {bad} outside the {len(embeds)} encoded images")
        return self.kind(pooler_output=tuple(embeds[i] for i in image_ids),
                         deepstack_features=[tuple(lv[i] for i in image_ids) for lv in levels])  # fmt: skip

    def inject(self, module, args, kwargs):
        """Forward pre-hook: the pending chunk's features and grid in place of pixel values."""
        if self.pending is not None:
            kwargs["mm_encoder_outputs"], kwargs["image_grid_thw"] = {"image": self.pending[0]}, self.pending[1]
        return args, kwargs


class RopeTemplate(typing.NamedTuple):
    """One prompt layout's rope positions (C2): `pos` (3, Lp), the next position `end`, its ids and token types."""

    pos: torch.Tensor
    end: int
    ids: torch.Tensor
    types: torch.Tensor


def rope_template(base, input_ids, mm_token_type_ids, image_grid_thw) -> RopeTemplate:
    """One HF get_rope_index call on one bare prompt row (Lp,), on host copies; the template lives on the row's
    device. Text after the prompt continues from `end` = max + 1 on all three axes."""
    pos, _ = base.get_rope_index(input_ids[None].cpu(), mm_token_type_ids=mm_token_type_ids[None].cpu(),
                                 image_grid_thw=image_grid_thw.cpu())  # fmt: skip
    pos = pos[:, 0]
    return RopeTemplate(pos.to(input_ids.device), int(pos.max()) + 1, input_ids, mm_token_type_ids.long())


def rope_positions(tmpl: RopeTemplate, input_ids, attention_mask, mm_token_type_ids):
    """(positions (3, B, W), bad (B,)) for rows that start with the template's prompt; no host read. Real token r
    (cumsum of the mask - 1) takes tmpl.pos[:, r] inside the prompt, tmpl.end + r - Lp after it; masked columns 0.
    `bad` marks rows whose first Lp real tokens are not the template's, with a non-text token after the prompt, or
    with fewer than Lp real tokens."""
    n = tmpl.ids.numel()
    real = attention_mask != 0
    r = real.long().cumsum(1) - 1
    inside, idx = r < n, r.clamp(0, n - 1)
    pos = torch.where(inside, tmpl.pos[:, idx], tmpl.end + r - n)
    pos = torch.where(real, pos, 0)
    head = real & inside & ((input_ids != tmpl.ids[idx]) | (mm_token_type_ids != tmpl.types[idx]))
    tail = real & ~inside & (mm_token_type_ids != 0)
    return pos, (head | tail).any(1) | (real.sum(1) < n)


def rope_bad_rows(prompt_ids, completion_ids, tmpl_ids, vision_ids) -> np.ndarray:
    """The host's routing decision per row: its prompt ids differ from its template's, or its answer holds an image,
    video or vision start/end id (TRL types every answer token as text)."""
    vision = set(vision_ids)
    return np.array([tuple(p) != tuple(t) or not vision.isdisjoint(c)
                     for p, c, t in zip(prompt_ids, completion_ids, tmpl_ids, strict=True)], dtype=bool)  # fmt: skip


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




def no_lora_input_cast(model) -> int:
    """Turns off PEFT's cast of each LoRA layer's input to its fp32 adapter dtype; returns the layers changed. Only
    under bf16 autocast: there the adapter matmul casts the input straight back to bf16, so the fp32 copy (one per
    LoRA layer per forward and recompute) changes nothing; without autocast the cast is what matches the dtypes."""
    layers = [m for _, m in model.named_modules() if hasattr(m, "cast_input_dtype_enabled")]
    for m in layers:
        m.cast_input_dtype_enabled = False
    return len(layers)


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


def one_pad_per_image(prompts: list[dict], pad: int) -> list[dict]:
    """vLLM prompts with each run of image pads cut to one pad (rows without multi_modal_data unchanged)."""
    out = []
    for row in prompts:
        if "multi_modal_data" in row:
            ids = row["prompt_token_ids"]
            row = row | {"prompt_token_ids": [t for i, t in enumerate(ids) if t != pad or not i or ids[i - 1] != pad]}
        out.append(row)
    return out


def prompt_images(prompts) -> list[list]:
    """Each conversational prompt's image objects, in order, as TRL extracts them (grpo_trainer.py 1760-1771)."""
    return [[c["image"] for m in p if isinstance(m["content"], list)
             for c in m["content"] if c["type"] == "image"] for p in prompts]  # fmt: skip


class ImageMemo:
    """The processor's image step once per distinct image per generation batch (perf: images prepared 2-3 times).

    TRL runs the image processor over every row twice a step: for the prompt ids sent to vLLM (`_tokenize_prompts`,
    whose pixel values it then drops) and for the update's pixel values. A step has 256 rows but 32 distinct images,
    and dataset.OpenOnce hands both passes the same object for the 8 rows of an image, so this memo (keyed by the object
    and the call's kwargs, holding the object so its id stays unique) prepares each image once and assembles each
    pass's output from the per-image slices, in row order. FastGRPO clears it before each generation. Every other
    attribute is the wrapped processor's."""

    KEYS = {"pixel_values", "image_grid_thw"}  # Qwen-VL image processors: patch rows and one (t, h, w) row per image

    hits = misses = fallback_calls = 0  # C1 counters: rows served from the memo, images prepared, CPU fallbacks

    def __init__(self, inner, gpu: tuple[str, str] | None = None, per_image: bool = False):
        """gpu: (gpu_pixels backend, device) to prepare missed images with gpu_pixels.patchify (C1), or None for the
        processor itself (the B1 control). per_image (--encoder-once, C3): pixel_values holds each distinct image of
        the call once (first-appearance order) instead of once per row; image_grid_thw stays one row per row."""
        self.inner, self.gpu, self.per_image = inner, gpu, per_image
        self.memo: dict[tuple, tuple] = {}  # (id(image), kwargs) -> (image, its patch rows, its grid row)

    @classmethod
    def install(cls, processing_class):
        """Wrap processing_class.image_processor once; None for a processing class without one."""
        ip = getattr(processing_class, "image_processor", None)
        if ip is None or isinstance(ip, cls):
            return ip
        processing_class.image_processor = memo = cls(ip)
        return memo

    def __getattr__(self, name):
        if "inner" not in self.__dict__:  # copy / pickle probe a bare instance before __init__ state exists
            raise AttributeError(name)
        return getattr(self.__dict__["inner"], name)

    def __deepcopy__(self, memo):
        """The plain image processor: TRL saves the processor with every checkpoint, and ProcessorMixin.to_dict
        deep-copies its attributes, so the saved files are the processor's own, without the cached pixels."""
        import copy

        return copy.deepcopy(self.inner, memo)

    def clear(self) -> None:
        """Drop the cached images and zero the counters (one generation batch: 32 misses, 256 hits expected)."""
        self.memo.clear()
        self.hits = self.misses = self.fallback_calls = 0

    def prefill(self, images, **kwargs) -> None:
        """Prepare the distinct images of a generation batch before TRL's processor calls (C1: on the GPU); counts
        one miss per distinct image and no hit."""
        from transformers.image_utils import make_flat_list_of_images

        self(list({id(x): x for x in make_flat_list_of_images(images)}.values()), **kwargs)

    def prepare(self, miss: list, kwargs: dict):
        """(pixel_values, image_grid_thw) of the missed images: gpu_pixels.patchify when self.gpu is set and every
        call kwarg is the processor's own value (else the processor: a counted fallback), or the processor."""
        if self.gpu is not None:
            from objdel.train import gpu_pixels

            own = all(v is None or getattr(self.inner, k, object()) == v for k, v in kwargs.items()
                      if k != "return_tensors")  # fmt: skip
            try:
                if own:
                    return gpu_pixels.patchify(miss, self.inner, device=self.gpu[1], backend=self.gpu[0])
            except gpu_pixels.UnsupportedImage:
                pass
            self.fallback_calls += 1
        out = self.inner(miss, **kwargs)
        return None if set(out.keys()) != self.KEYS else (out["pixel_values"], out["image_grid_thw"])

    def __call__(self, images, **kwargs):
        """Image processor call that processes only images not yet in the memo and returns pixel values and grids of
        all, in order.
        """
        with nvtx.range("image_memo"):
            from transformers.feature_extraction_utils import BatchFeature
            from transformers.image_utils import make_flat_list_of_images

            flat, kw = make_flat_list_of_images(images), repr(sorted(kwargs.items()))
            miss = list({id(x): x for x in flat if (id(x), kw) not in self.memo}.values())
            self.misses, self.hits = self.misses + len(miss), self.hits + len(flat) - len(miss)
            if miss:
                out = self.prepare(miss, kwargs)
                if out is None:  # another processor family: no memo
                    return self.inner(images, **kwargs)
                pv, grid = out
                rows = torch.split(pv, grid.prod(-1).tolist())
                for x, r, g in zip(miss, rows, grid, strict=True):
                    self.memo[(id(x), kw)] = (x, r, g)
            parts = [self.memo[(id(x), kw)] for x in flat]
            once = {id(x): p for x, p in zip(flat, parts, strict=True)}.values() if self.per_image else parts
            data = {"pixel_values": torch.cat([r for _, r, _ in once]),
                    "image_grid_thw": torch.stack([g for *_, g in parts])}
            return BatchFeature(data, tensor_type=kwargs.get("return_tensors"))


@contextlib.contextmanager
def vllm_kwargs(**kw):
    """Inside the block TRL's module-level `LLM` (imported at trl/generation/vllm_generation.py top) gets these
    keyword arguments; no `kw` leaves it alone."""
    if not kw:
        yield
        return
    import trl.generation.vllm_generation as vg

    orig = vg.LLM
    vg.LLM = functools.partial(orig, **kw)
    try:
        yield
    finally:
        vg.LLM = orig


class FastGRPO(FastForward, GRPOTrainer):
    """GRPOTrainer with the per-chunk windows, the vision cache (cleared per generation batch), rewards scored while
    the old and reference logprobs run (C9), phase timers and no FLOP count."""

    rows = None  # the generation batch's dataset rows; None outside one (g_probe's step-0 sampling calls _generate)
    loss_inputs = None  # the micro-batch inside _compute_loss (C10: its forward sets the IS ratio)
    loss_meta = None  # the micro-batch's RowMeta inside _compute_loss under host_meta (C4: its forward reads it)
    is_stats = None  # C10, B2.1: the generation batch's merged is_summary, on the device until `log`
    seq_gaps = []  # B2.1: the step's per-sequence (seq_kl, advantage, valid tokens) rows (train mode), on the device
    vision = rope = None  # VisionCache, RopeCache (set in __init__ when the model has a Qwen-VL tower and base)
    encoder = None  # C3, train_rl --encoder-once: EncoderCache in place of VisionCache (use_encoder_once)
    row_image_ids = None  # C3: the generation batch's per-row image ids, for its reference and old passes
    images = None  # ImageMemo on the processing class's image processor (set in __init__ when it has one)
    host_meta = False  # C4, train_rl --host-meta: RowMeta from _generate's host lists, row-aligned in the batch
    check_meta_every = 50  # C4, train_rl --check-meta-every N: compare host_meta with the device masks (0: never)
    check_rope_every = 50  # C2, train_rl --check-rope-every N: compare template rope positions with HF (0: never)
    rope_check = False  # this generation batch's forwards run the rope check (FastForward reads it)
    meta = None  # the generation batch's RowMeta, from _generate to the batch
    prompt_memo = False  # C2, train_rl --prompt-memo: _tokenize_prompts tokenizes one row per prompt layout
    count_syncs = False  # C0, train_rl --count-syncs: perf/host_syncs/<phase> per timed phase
    cuda_times = None  # OBJDEL_NVTX_SYNC=0: (mode, phase, start, end) CUDA event pairs of `timed`, read by `log`
    last_log = None  # perf_counter at the previous `log` (time/step_s)
    group_requests = False  # C6, train_rl --group-requests: vLLM gets one n = G request per prompt (vllm_requests)
    opsm_watch = None  # B2.4 rule 4: callbacks.OpsmDelta, fed the train-mode per-sequence (seq_kl, advantage) rows
    log_by_prompt = False  # B2.3: train logs also hold rewards/iou_tc/by_prompt (one mean per prompt of the step)
    runaway_cap = False  # C11, train_rl --runaway-cap: max_tokens per prompt from the rows' `cap` (runaway)

    def __init__(self, *args, **kwargs):
        trl_adapter.check()  # the private TRL names this class relies on, before TRL builds anything
        fp16 = getattr(kwargs.get("args"), "fp16", False)
        with vllm_kwargs(**({"dtype": "float16"} if fp16 else {})):  # TRL builds its vLLM engine without a dtype
            super().__init__(*args, **kwargs)
        self.check_engine_dtype()
        self.vision, self.rope = VisionCache.install(self.model), RopeCache.install(self.model)
        self.images = ImageMemo.install(self.processing_class)
        self.drop_lora_input_cast()
        self.sync_without_merge()

    def use_gpu_pixels(self, images, allow_cpu: bool = False) -> str:
        """C1 start-up self-check (train_rl --gpu-pixels auto): the first gpu_pixels backend equal to the CPU
        processor on `images` (training views) plus an in-memory digit image prepares the memo's misses on CUDA.
        None passing raises, unless allow_cpu (then the CPU processor stays: returns "cpu")."""
        from objdel.train import gpu_pixels

        if self.images is None or not torch.cuda.is_available():
            raise RuntimeError("--gpu-pixels auto needs CUDA and a Qwen-VL image processor")
        digit = (np.arange(64 * 64 * 3) & 255).astype(np.uint8).reshape(64, 64, 3)
        backend = gpu_pixels.pick_backend([*images, digit], self.images.inner, "cuda")
        if backend is None and not allow_cpu:
            raise RuntimeError("no GPU pixel backend equals the CPU processor (pass --allow-cpu-pixels to run on it)")
        self.images.gpu = None if backend is None else (backend, "cuda")
        return backend or "cpu"

    def use_encoder_once(self) -> bool:
        """C3, train_rl --encoder-once (needs --host-meta: the update's rows carry their image ids): EncoderCache takes
        over VisionCache's tower. False without a Qwen-VL tower, base and image processor."""
        if not self.host_meta:
            raise ValueError("--encoder-once needs --host-meta")
        if any(x is None for x in (self.vision, self.rope, self.images)):
            return False
        self.encoder, self.vision = EncoderCache(self.vision, self.rope.base), None
        return True

    def sync_without_merge(self) -> bool:
        """vLLM weight sync through lora_merged_params (the base weights are never written), on one process; FSDP and
        ZeRO-3 keep TRL's gathered path."""
        vg = getattr(self, "vllm_generation", None)
        if vg is None or not hasattr(self.model, "peft_config"):
            return False
        if vg._dist.is_fsdp or self.accelerator.num_processes > 1:
            return False
        vg._iter_named_params = lambda: lora_merged_params(vg.model)
        return True

    def check_engine_dtype(self) -> None:
        """fp16 runs only: the vLLM engine must hold float16 (the vllm_kwargs patch took effect)."""
        vg = getattr(self, "vllm_generation", None)
        if vg is None or not getattr(self.args, "fp16", False):
            return
        got = vg.llm.llm_engine.model_config.dtype
        if got != torch.float16:
            raise RuntimeError(f"fp16 run, but the vLLM engine dtype is {got}, not torch.float16")

    def drop_lora_input_cast(self) -> bool:
        """no_lora_input_cast when TRL trains in bf16 or fp16 (Accelerate then runs every forward under autocast)."""
        if not (getattr(self.args, "bf16", False) or getattr(self.args, "fp16", False)):
            return False
        no_lora_input_cast(self.model)
        return True

    def get_train_dataloader(self):
        """TRL's loader over dataset.OpenOnce rows: one stack of 256 rows and 32 image opens a step, not 16 x 256."""
        if self.args.dataloader_num_workers:
            raise ValueError("OpenOnce keeps the step's stack in the main process: dataloader_num_workers must be 0, "
                             f"not {self.args.dataloader_num_workers}")
        rows, self.train_dataset = self.train_dataset, OpenOnce(self.train_dataset)
        try:
            return super().get_train_dataloader()
        finally:
            self.train_dataset = rows

    def _generate_and_score_completions(self, inputs):
        for cache in (self.vision, self.rope, self.encoder):
            if cache is not None:
                cache.clear()
        self.rows = inputs
        skip = self.skip_old_logps()
        if skip:  # TRL then skips the old-logprob pass and its IS ratio; the update's forward sets it (C10)
            self.vllm_importance_sampling_correction = False
        try:
            with self.timed("batch"):
                batch = super()._generate_and_score_completions(inputs)
            if self.encoder is not None:  # C3: the update takes the encoder's features; TRL would cut these by grid
                batch.pop("pixel_values", None)
            if self.meta is not None:
                batch["host_meta"] = host_rows(self.meta)
                n = self.check_meta_every
                if n and self.state.global_step % n == 0:  # steps 1, N + 1, 2N + 1, ...
                    check_meta(batch)
                n = self.check_rope_every
                self.rope_check = bool(n) and self.state.global_step % n == 0
            return batch
        finally:
            self.rows, self.meta, self.row_image_ids = None, None, None
            if skip:
                self.vllm_importance_sampling_correction = True

    def skip_old_logps(self) -> bool:
        """C10: with one generation per optimizer step and one iteration (TRL's own aligned case), the weights do not
        change between the old-logprob pass and the update, so the update's logprobs (detached) are the old ones.
        TRL still runs the old pass under vLLM, only for the IS ratio: the update's forward computes it instead."""
        a = getattr(self, "args", None)
        return bool(a is not None and self.use_vllm and self.vllm_importance_sampling_correction
                    and self.loss_type != "vespo"
                    and a.gradient_accumulation_steps % (a.steps_per_generation * self.num_iterations) == 0)

    def _prepare_inputs(self, generation_batch):
        """Update micro-batches cut from the length-sorted generation batch (by_length) when the partition leaves the
        update unchanged: dr_grpo's constant normaliser and every micro-batch in one optimizer step. Each chunk is
        then trimmed to rows of like length instead of the batch's longest answer."""
        a = self.args
        if not (self.model.training and self.loss_type == "dr_grpo" and getattr(self, "top_entropy_quantile", 1.0) >= 1
                and a.gradient_accumulation_steps % a.steps_per_generation == 0):  # fmt: skip (entropy cut: per batch)
            return super()._prepare_inputs(generation_batch)
        order = by_length
        G = getattr(self, "num_generations", None)
        if self.shared_prefix and G and a.per_device_train_batch_size % G == 0:
            order = functools.partial(by_group_length, G=G)  # whole groups per micro-batch: packed_chunk shares them
        shuffle, trl_grpo.shuffle_sequence_dict = trl_grpo.shuffle_sequence_dict, order
        try:
            return super()._prepare_inputs(generation_batch)
        finally:
            trl_grpo.shuffle_sequence_dict = shuffle

    def training_step(self, *args, **kwargs):
        """Run the micro-batch step inside the `micro` NVTX range."""
        with nvtx.range("micro"):  # the micro-batch's forward and backward (the first one also generates: `batch`)
            return super().training_step(*args, **kwargs)

    def get_batch_samples(self, *args, **kwargs):
        """Fetch batch samples inside the `loader` NVTX range."""
        with nvtx.range("loader"):  # transformers fetches the step's 16 micro-batches before on_step_begin
            return super().get_batch_samples(*args, **kwargs)

    def _compute_loss(self, model, inputs):
        if "importance_sampling_ratio" not in inputs and self.use_vllm and self.vllm_importance_sampling_correction:
            self.loss_inputs = inputs
        if self.host_meta:
            self.loss_meta = meta_of({"host_meta": inputs["host_meta"]})  # KeyError: rows lost their host_meta
        try:
            with nvtx.range("forward"):
                return super()._compute_loss(model, inputs)
        finally:
            self.loss_inputs = self.loss_meta = None

    def _get_per_token_logps_and_entropies(self, *args, **kwargs):
        if self.loss_meta is not None:
            kwargs["host_meta"] = self.loss_meta
        with nvtx.range("logps"):
            out = super()._get_per_token_logps_and_entropies(*args, **kwargs)
        if self.loss_inputs is not None:
            self.set_is_ratio(out[0].detach(), self.loss_inputs)
        return out

    def set_is_ratio(self, logps, inputs):
        """inputs["importance_sampling_ratio"] from the update's own logprobs. The micro-batch's is_summary is merged
        into the generation batch's on the device, and in train mode its per-sequence (seq_kl, advantage, valid
        tokens) rows are kept for `log` (no host read here)."""
        sampling = inputs.get("sampling_per_token_logps")
        if sampling is None:
            raise RuntimeError("C10 needs vLLM's sampling logprobs (vllm_generation logprobs=0) for the IS ratio")
        mask = (inputs["completion_mask"] * inputs.get("tool_mask", 1)).bool()
        mode, lo, hi = (self.vllm_importance_sampling_mode, self.vllm_importance_sampling_clip_min,
                        self.vllm_importance_sampling_clip_max)  # fmt: skip
        r = inputs["importance_sampling_ratio"] = is_ratio(logps, sampling, mask, mode, lo, hi)
        adv = inputs["advantages"].reshape(-1)
        new = is_summary(logps, sampling, mask, r, adv, mode, lo, hi, self.off_policy_mask_threshold)
        self.is_stats = merge_summary(self.is_stats, new)
        if self.model.training:  # eval rows inform the eval metrics only; they never enter seq_gaps (B2.1, R01c)
            kl = seq_kl(logps, sampling, mask)
            rows = torch.stack([kl, adv.to(logps.dtype), mask.sum(-1).to(logps.dtype)], 1)
            if self.opsm_watch is not None:
                self.opsm_watch.add(kl, adv)
            self.seq_gaps = [*self.seq_gaps, rows.float()]

    def log(self, logs, start_time=None):
        """TRL's sampling metrics of the generation batch, the A6 zero / hi / lo shares and the B2.1 keys, from the
        device in one read (C10, A6, B2.1), time/step_s (the host-clock gap since the previous log: one log per step
        at logging_steps 1) and, with OBJDEL_NVTX_SYNC=0, the phase times of the CUDA event pairs; then TRL's log.
        In train mode the same read takes the step's seq_kl quantiles and its per-sequence rows, appended to
        <output_dir>/seq_gaps.f32 (write_seq_gaps)."""
        with nvtx.range("log"):
            train = self.model.training
            m = self._metrics["train" if train else "eval"]
            now = time.perf_counter()
            if self.last_log is not None:
                m["time/step_s"].append(now - self.last_log)
            self.last_log = now
            if torch.cuda.is_available():  # G-MEM (C3 step 4): the peak so far, read after the update
                m["perf/max_memory_gb"].append(torch.cuda.max_memory_allocated() / 2**30)
            if self.is_stats is not None:
                rows = torch.cat(self.seq_gaps) if train and self.seq_gaps else None
                parts = [self.is_stats]
                if rows is not None:  # scalar q: no host-to-device copy
                    parts += [torch.stack([torch.quantile(rows[:, 0], q) for q in SEQ_KL_Q]), rows.flatten()]
                vals, n = torch.cat(parts).tolist(), len(SUMMARY)
                self.is_stats, self.seq_gaps = None, []
                s = dict(zip(SUMMARY, vals[:n], strict=True))
                for k, v in summary_metrics(s, self.off_policy_mask_threshold).items():
                    m[f"sampling/{k}"].append(v)
                if rows is not None:
                    for k, v in zip(("q50", "q99", "q99.5"), vals[n:n + 3], strict=True):
                        m[f"sampling/seq_kl/{k}"].append(v)
                    if self.accelerator.is_main_process:
                        write_seq_gaps(self.args.output_dir, vals[n + 3:], len(rows))
            for mode, phase, start, end in self.cuda_times or ():
                end.synchronize()  # long done: TRL's loss read and the stats read above waited for the stream
                self._metrics[mode][f"time/{phase}_s"].append(start.elapsed_time(end) / 1e3)
            self.cuda_times = None
            if self.log_by_prompt and train:  # B2.3; reward_funcs read only when asked (fake trainers have none)
                reward = next((f.trl_reward for f in self.reward_funcs if hasattr(f, "trl_reward")), None)
                if reward is not None:
                    logs = {**logs, "rewards/iou_tc/by_prompt": reward.by_prompt}
            return super().log(logs, start_time)

    def _tokenize_prompts(self, prompts):
        """With prompt_memo on the vLLM path: TRL's own tokenization once per (rendered text, image sizes), its
        prompt ids reused by every row of that layout, each row's own image objects. vLLM drops the multimodal
        fields (C1 note S1), so none are built."""
        if not (self.prompt_memo and self.use_vllm and self._is_vlm):
            return super()._tokenize_prompts(prompts)
        rows = [trl_grpo.prepare_multimodal_messages(p) for p in prompts]
        images = prompt_images(rows)
        memo, ids = {}, []
        for prompt, row, imgs in zip(prompts, rows, images, strict=True):
            text = self.processing_class.apply_chat_template(
                row, chat_template=self.chat_template, tokenize=False, add_generation_prompt=True,
                **self.chat_template_kwargs)  # fmt: skip
            key = (text, tuple(i.size for i in imgs))
            if key not in memo:
                memo[key] = super()._tokenize_prompts([prompt])[0][0]
            ids.append(memo[key])
        return ids, ([i or None for i in images] if any(images) else None), {}

    def _generate(self, prompts):
        if self.images is not None:  # the previous batch's images
            self.images.per_image = self.encoder is not None  # C3: one pixel block per distinct image
            self.images.clear()
            flat = [x for row in prompt_images(prompts) for x in row] if self.prompt_memo else []
            if flat:  # pass 1 prepares no image with --prompt-memo: fill the memo with pass 2's kwargs (C1)
                self.images.prefill(flat, return_tensors="pt")
        if self.use_vllm:
            self.patch_vllm_prompts()
            self.time_vllm_phases()
        with self.timed("generate"):
            out = super()._generate(prompts)
        if self.host_meta and self.rows is not None and len(self.rows) == len(out[1]):
            eos_and_pad = [*self.eos_token_ids, self._tokenizer.pad_token_id]  # TRL's truncation list
            cfg = getattr(self.model, "config", None)
            vision = [getattr(cfg, k, None) for k in ("image_token_id", "video_token_id", "vision_start_token_id",
                                                      "vision_end_token_id")]  # fmt: skip
            self.meta = row_meta(out[0], out[1], image_ids(self.rows), self.num_generations,
                                 self.mask_truncated_completions, eos_and_pad,
                                 [getattr(r.get("image"), "size", None) for r in self.rows],
                                 [v for v in vision if v is not None])  # fmt: skip
        if self.encoder is not None and self.rows is not None:
            self.encode_images()
        self.prefetch_rewards(out[3], out[1])
        self.log_capped(out[1])
        return out

    def encode_images(self) -> None:
        """C3: the generation batch's distinct images (the memo's pass-2 entries, prepared here when missing) through
        the tower in one call; row_image_ids then gives the reference and old passes each row's image. Nothing for
        rows without an image."""
        imgs = list({id(r.get("image")): r.get("image") for r in self.rows}.values())
        if any(x is None for x in imgs):
            return
        self.row_image_ids = np.asarray(image_ids(self.rows), dtype=np.int32)
        px = self.images(imgs, return_tensors="pt")
        self.encoder.encode(px["pixel_values"], px["image_grid_thw"])

    def log_capped(self, completion_ids):
        """C11: completions/capped_share, the share of this batch's answers that reached their row's cap."""
        if self.runaway_cap and self.rows is not None:
            hit = [len(ids) >= r["cap"] for ids, r in zip(completion_ids, self.rows, strict=True)]
            self._metrics["train" if self.model.training else "eval"]["completions/capped_share"].append(
                sum(hit) / len(hit))

    def prefetch_rewards(self, completions, completion_ids):
        """Start scoring (objdel.train.rewards.TRLReward.prefetch) on the answers vLLM just returned, with the
        dataset columns of the same rows, so the checker and the CPU scorers overlap TRL's old and reference
        logprob forwards. The reward call later checks it got the same answers."""
        reward = next((f.trl_reward for f in self.reward_funcs if hasattr(f, "trl_reward")), None)
        if reward is not None and self.rows is not None and len(self.rows) == len(completions):
            cols = [[r[k] for r in self.rows] for k in ("key", "W", "H", "names")]
            reward.prefetch(completions, completion_ids, *cols)

    def _calculate_rewards(self, inputs, prompts, completions, completion_ids_list):
        with self.timed("rewards"):
            return super()._calculate_rewards(inputs, prompts, completions, completion_ids_list)

    @contextlib.contextmanager
    def timed(self, phase: str):
        """phase_clock(phase); with count_syncs (train_rl --count-syncs, C0) also perf/host_syncs/<phase>, torch's
        synchronizing operations inside the phase and outside its nested phases (nvtx.count_syncs). The vLLM range
        (`generate`) logs perf/host_syncs_vllm instead, for information only."""
        if not self.count_syncs:
            with self.phase_clock(phase):
                yield
            return
        n = None
        try:
            with nvtx.count_syncs() as n, self.phase_clock(phase):
                yield
        finally:
            if n is not None:
                key = "perf/host_syncs_vllm" if phase == "generate" else f"perf/host_syncs/{phase}"
                self._metrics["train" if self.model.training else "eval"][key].append(n.syncs)

    @contextlib.contextmanager
    def phase_clock(self, phase: str):
        """time/<phase>_s in TRL's metrics (logged as the mean per logging step), inside the NVTX range `phase`, also
        when the phase raises. By default one device sync at the end, so queued GPU work is charged to the phase that
        queued it (a few per step: negligible). With OBJDEL_NVTX_SYNC=0 and CUDA: no sync, a CUDA event pair on the
        current stream instead, read by `log` (a trace then shows the real overlap between phases)."""
        mode = "train" if self.model.training else "eval"
        with nvtx.range(phase):
            if not nvtx.SYNC and torch.cuda.is_available():
                start, end = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
                start.record()
                try:
                    yield
                finally:
                    end.record()
                    self.cuda_times = [*(self.cuda_times or ()), (mode, phase, start, end)]
                return
            t = time.perf_counter()
            try:
                yield
            finally:
                if nvtx.SYNC and torch.cuda.is_available():
                    torch.cuda.synchronize()
                self._metrics[mode][f"time/{phase}_s"].append(time.perf_counter() - t)

    def time_vllm_phases(self):
        """time/sync_s (weight push, with its weights wake-up), time/wake_s and time/sleep_s: the colocate
        handovers where the GPU idles between generation and training. Idempotent; the engine exists only after TRL's
        vLLM init."""
        vg = self.vllm_generation
        hooks = ((vg, "sync_weights", "sync"), (vg.llm, "wake_up", "wake"), (vg.llm, "sleep", "sleep"))
        for host, name, phase in hooks:
            fn = getattr(host, name)
            if getattr(fn, "timed_phase", None):
                continue

            def wrapped(*args, _fn=fn, _phase=phase, **kwargs):
                with self.timed(_phase):
                    return _fn(*args, **kwargs)

            wrapped.timed_phase = phase
            setattr(host, name, wrapped)

    def patch_vllm_prompts(self):
        """vLLM gets one image pad per image (one_pad_per_image) and expands it itself. TRL trains on its own
        processor-expanded ids (grpo_trainer.py 1832, `_generate_single_turn`), so vLLM's expansion must give exactly
        those ids, or vLLM samples from other image features than the trainer scores (off-policy, silently): checked
        per prompt.
        Idempotent; the engine exists only after TRL's vLLM init."""
        llm = self.vllm_generation.llm
        if getattr(llm.generate, "one_pad", False):
            return
        generate, pad = llm.generate, self.model.config.image_token_id

        def one_pad(prompts, *args, **kwargs):
            outputs = generate(one_pad_per_image(prompts, pad), *args, **kwargs)
            for row, o in zip(prompts, outputs, strict=True):
                if "multi_modal_data" in row and list(o.prompt_token_ids) != list(row["prompt_token_ids"]):
                    n, m = list(o.prompt_token_ids).count(pad), list(row["prompt_token_ids"]).count(pad)
                    raise RuntimeError(f"vLLM expanded the image to {n} pads, the trainer's processor to {m}: "
                                       "vLLM would sample from other image features than the trainer scores")
            return outputs

        one_pad.one_pad = True
        llm.generate = one_pad
        if self.runaway_cap:  # C11: the generation batch's rows hold the caps; outside one (g_probe) the limit stays

            def capped(prompts, *args, sampling_params=None, **kwargs):
                if self.rows is not None:
                    if len(self.rows) != len(prompts):
                        raise ValueError(f"{len(self.rows)} caps for {len(prompts)} prompts")
                    sampling_params = runaway.per_row_params(sampling_params, [r["cap"] for r in self.rows])
                return one_pad(prompts, *args, sampling_params=sampling_params, **kwargs)

            capped.one_pad = True
            llm.generate = capped
        if self.group_requests:  # C6: outside one_pad, so its check compares the 32 sent prompts

            def by_group(prompts, *args, **kwargs):
                G = self.num_generations if self.model.training else self.num_generations_eval
                return vllm_requests.grouped(one_pad, G)(prompts, *args, **kwargs)

            by_group.one_pad = True
            llm.generate = by_group

    def floating_point_ops(self, inputs):
        """FLOP count: always 0, which skips the Trainer's per-micro-batch module walk."""
        return 0  # perf F3: the Trainer default walks ~2.5k PEFT modules per micro-batch


def lora_merged_params(model, with_delta: bool = True):
    """(vLLM name, tensor) for every weight of a PEFT model, as TRL 1.17's VLLMGeneration._iter_named_params names
    them, with each LoRA layer's delta added into a new tensor. TRL instead runs merge_adapter / unmerge_adapter on the
    bf16 base weights in place, and W + d - d != W in bf16, so the frozen base drifts at every sync (trl #6688,
    #7423). Same values pushed (W + delta in W's dtype, as PEFT's merge computes it); no parameter is written.
    with_delta=False: the same names with the base W alone (vllm_sync's option-F base push)."""
    from peft.tuners.lora import LoraLayer

    deltas = {}
    for name, m in model.named_modules():
        if isinstance(m, LoraLayer) and not m.merged:
            ads = [a for a in m.active_adapters if a in m.lora_A]
            if ads and (getattr(m, "use_dora", {}).get(ads[0]) or getattr(m, "lora_bias", {}).get(ads[0])):
                raise NotImplementedError(f"{name}: DoRA / LoRA bias are not merged here")
            if ads and with_delta:
                deltas[f"{name}.base_layer.weight"] = sum(m.get_delta_weight(a) for a in ads)
    for name, p in model.named_parameters():
        out = name.removeprefix("base_model.model.").replace(".base_layer", "")
        if model.prefix in out or "original_module" in out:
            continue
        out = out.replace("_checkpoint_wrapped_module.", "").replace("modules_to_save.default.", "")
        d = deltas.get(name)
        yield out, (p.data if d is None else p.data + d.to(p.dtype))


def compile_text_layers(model, mode: str) -> int:
    """Stage F arm: torch.compile each text decoder layer in place (regional compilation: one code object, BUCKETS
    shapes; parameter names unchanged). mode: "default" or "reduce-overhead" (CUDA Graphs). Returns the layer count."""
    import torch._dynamo.config as dc

    # one code object serves every layer; its variants are BUCKETS widths x grad on/off (update vs old/ref) x adapter
    # on/off (the ref pass): past recompile_limit (default 8) Dynamo runs eager with only a warning
    dc.recompile_limit = max(dc.recompile_limit, 64)
    dc.accumulated_recompile_limit = max(dc.accumulated_recompile_limit, 1024)
    layers = next(m for n, m in model.named_modules() if n.endswith("language_model.layers"))
    for layer in layers:
        layer.compile(mode=mode, dynamic=False)
    return len(layers)
