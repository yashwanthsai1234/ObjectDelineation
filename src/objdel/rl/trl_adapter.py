"""The one TRL adapter (E5, ADR 0019): FastGRPO and FastForward (the speed patches of TRL's GRPOTrainer), sample (the
trainer's own generation), and the start-up check of the private TRL names and the silent couplings they rely on
(TRL pinned to 1.17.1). No other module uses a private TRL name (tests/core/test_layers.py).

`check()` runs first in `FastGRPO.__init__` and in the preflight ("CODE TRL private names"). It reads signatures and
source text only, so it runs on the Mac without vLLM or a GPU. Under lora-seqmask, `check_lora(llm)` also checks the
option-F vLLM names (LORA_POINTS) on the live engine, in `FastGRPO.__init__` once TRL built it (the preflight's
configs are train_sft runs with no numerics preset, and the Mac has no vLLM engine). A TRL change that breaks one
item fails in seconds, naming the item and the place in this module that overrides or reads it. The version is checked
last: with OBJDEL_TRL_UNPINNED=1 a different TRL version only warns, so an upgrade can be tried on the contract alone.
"""

import contextlib
import functools
import inspect
import os
import time
import warnings
from collections.abc import Iterable

import numpy as np
import torch
import torch.nn.functional as F
import trl.trainer.grpo_trainer as trl_grpo
from trl import GRPOTrainer
from trl.data_utils import prepare_multimodal_messages

from objdel.geom import engine as _gpucheck
from objdel.policy import marks as nvtx
from objdel.policy import models
from objdel.policy.caches import EncoderCache, ImageMemo, RopeCache, VisionCache
from objdel.policy.sampler import G
from objdel.rl import record, runaway, vllm_requests
from objdel.rl.forward import (
    SEQ_KL_Q,
    SUMMARY,
    by_group_length,
    by_length,
    check_meta,
    grid_rows,
    host_rows,
    image_ids,
    is_ratio,
    is_summary,
    merge_summary,
    meta_of,
    meta_runs,
    packing_wins,
    prompt_groups,
    row_meta,
    seq_kl,
    summary_metrics,
    window,
    window_from_meta,
    write_seq_gaps,
)
from objdel.rl.rows import OpenOnce

PIN = "1.17.1"

# name -> (parameter names in order, the place that overrides or calls it)
METHODS = {  # GRPOTrainer
    "_generate_and_score_completions": (("self", "inputs"), "FastGRPO._generate_and_score_completions"),
    "_prepare_inputs": (("self", "generation_batch"), "FastGRPO._prepare_inputs"),
    "_compute_loss": (("self", "model", "inputs"), "FastGRPO._compute_loss"),
    "_get_per_token_logps_and_entropies": (
        ("self", "model", "input_ids", "attention_mask", "logits_to_keep", "batch_size", "compute_entropy",
         "compute_aux_loss", "pixel_values", "image_grid_thw", "num_images", "pixel_attention_mask", "spatial_shapes",
         "num_tiles", "image_sizes", "token_type_ids", "mm_token_type_ids", "image_position_ids", "sampling_mask"),
        "FastForward._get_per_token_logps_and_entropies"),  # fmt: skip
    "_generate": (("self", "prompts"), "FastGRPO._generate"),
    "_calculate_rewards": (("self", "inputs", "prompts", "completions", "completion_ids_list"),
                           "FastGRPO._calculate_rewards"),  # fmt: skip
}
GLOBALS = {"shuffle_sequence_dict": (("seq_dict",), "FastGRPO._prepare_inputs")}  # trl.trainer.grpo_trainer
VLLM_METHODS = {"_iter_named_params": (("self",), "FastGRPO.sync_without_merge")}  # VLLMGeneration

# name -> (sources that must contain the needle, needle, the place that relies on it)
PRIVATE_ATTRS = {
    "_metrics": (("GRPOTrainer.__init__",), "self._metrics", "FastGRPO.timed, FastGRPO.log"),
    "_dist": (("VLLMGeneration",), "self._dist", "FastGRPO.sync_without_merge"),
}
ATTRS = {
    **{n: (("GRPOTrainer",), f"self.{n} =", p) for n, p in (
        ("use_vllm", "FastGRPO.skip_old_logps"), ("loss_type", "FastGRPO._prepare_inputs"),
        ("num_iterations", "FastGRPO.skip_old_logps"), ("num_generations", "FastGRPO._prepare_inputs"),
        ("top_entropy_quantile", "FastGRPO._prepare_inputs"),
        ("vllm_importance_sampling_correction", "FastGRPO._generate_and_score_completions"),
        ("vllm_importance_sampling_mode", "FastGRPO.set_is_ratio"),
        ("vllm_importance_sampling_clip_min", "FastGRPO.set_is_ratio"),
        ("vllm_importance_sampling_clip_max", "FastGRPO.set_is_ratio"),
        ("vllm_generation", "FastGRPO.sync_without_merge"), ("off_policy_mask_threshold", "FastGRPO.set_is_ratio"))},
    "llm": (("VLLMGeneration",), "self.llm =", "FastGRPO.time_vllm_phases, FastGRPO.patch_vllm_prompts"),
    "model": (("VLLMGeneration",), "self.model =", "FastGRPO.sync_without_merge"),
    "sync_weights": (("VLLMGeneration",), "def sync_weights", "FastGRPO.time_vllm_phases"),
    "_llm_weights_sleeping": (("VLLMGeneration",), "self._llm_weights_sleeping =", "vllm_sync.push_weights, sync_lora"),
}  # fmt: skip
_GASC = "GRPOTrainer._generate_and_score_completions"
KEYS = {
    "sampling_per_token_logps": ((_GASC,), '"sampling_per_token_logps"', "FastGRPO.set_is_ratio"),
    "completion_mask": ((_GASC,), '"completion_mask"', "FastGRPO.set_is_ratio"),
    "tool_mask": ((_GASC,), '"tool_mask"', "FastGRPO.set_is_ratio"),
    "advantages": ((_GASC,), '"advantages"', "FastGRPO.set_is_ratio"),
    "importance_sampling_ratio": ((_GASC, "GRPOTrainer._compute_loss"), '"importance_sampling_ratio"',
                                  "FastGRPO.set_is_ratio, FastGRPO._compute_loss"),
}  # fmt: skip
BEHAVIOURS = {  # C10: TRL skips the old-logprob pass when the IS correction is off; forward.is_ratio copies this block
    "old-logprob pass gated on the IS correction": (
        (_GASC,), "self.use_vllm and self.vllm_importance_sampling_correction", "FastGRPO.skip_old_logps"),
    "IS ratio block": ((_GASC,), "per_token_logps_diff = (old_per_token_logps - sampling_per_token_logps) * mask",
                       "forward.is_ratio"),
    "log_metric to reward functions": (("GRPOTrainer._calculate_rewards",), 'reward_kwargs["log_metric"]',
                                       "TRLReward.funcs"),
    "log_extra to reward functions": (("GRPOTrainer._calculate_rewards",), 'reward_kwargs["log_extra"]',
                                      "TRLReward.funcs"),
    "completion_ids to reward functions": (("GRPOTrainer._calculate_rewards",), "completion_ids=completion_ids_list",
                                           "TRLReward.funcs"),
}  # fmt: skip
# Option F (lora-seqmask, B2.2 F step 3): vLLM names on the live engine (vg.llm), checked by check_lora at FastGRPO
# start-up under lora-seqmask only (vLLM is not on the Mac, so check() cannot read them).
# name -> (attribute path from vg.llm, what must hold, the place that relies on it)
LORA_POINTS = {
    "add_lora": ("llm_engine.add_lora", "callable", "vllm_sync.sync_lora, train_rl --f-check (load count)"),
    "remove_lora": ("llm_engine.remove_lora", "callable", "f_check.run_comparator (restore)"),
    "lora_index_to_id": ("llm_engine.model_executor.driver_worker.model_runner.lora_manager._adapter_manager"
                         ".lora_index_to_id", "a sequence", "f_check.lora_slot_index"),
    "lora_a_stacked / lora_b_stacked": ("llm_engine.model_executor.driver_worker.model_runner.model",
                                        "a module holding both", "f_check.lora_equal"),
}  # fmt: skip
_SIGNED = {"method": METHODS, "global": GLOBALS, "VLLMGeneration method": VLLM_METHODS}
_TEXT = {"private attribute": PRIVATE_ATTRS, "attribute": ATTRS, "key": KEYS, "behaviour": BEHAVIOURS}
_LIVE = {"option-F vLLM name": LORA_POINTS}


class TRLContractError(RuntimeError):
    """A private TRL name or coupling that this adapter relies on changed."""


def _fail(kind: str, name: str, expected, found, place: str):
    raise TRLContractError(f"TRL {kind} {name}: expected {expected!r}, found {found!r}; relied on by {place}. "
                           f"Pin is {PIN}.")  # fmt: skip


def check() -> None:
    """Raise TRLContractError on the first broken item (methods, globals, private attributes, attributes, keys,
    behaviours), then on a version other than PIN unless OBJDEL_TRL_UNPINNED=1 (then one warning)."""
    import trl
    import trl.trainer.grpo_trainer as grpo
    from trl.generation.vllm_generation import VLLMGeneration  # imports no vLLM module

    owners = {"method": grpo.GRPOTrainer, "global": grpo, "VLLMGeneration method": VLLMGeneration}
    for kind, items in _SIGNED.items():
        for name, (params, place) in items.items():
            fn = getattr(owners[kind], name, None)
            found = "missing" if fn is None else tuple(inspect.signature(fn).parameters)
            if found != params:
                _fail(kind, name, params, found, place)
    sources = {}
    for kind, items in _TEXT.items():
        for name, (where, needle, place) in items.items():
            for w in where:
                if w not in sources:
                    obj = {"GRPOTrainer": grpo.GRPOTrainer, "VLLMGeneration": VLLMGeneration}[w.split(".")[0]]
                    for attr in w.split(".")[1:]:
                        obj = getattr(obj, attr)
                    sources[w] = inspect.getsource(obj)
                if needle not in sources[w]:
                    _fail(kind, name, f"{needle} in the source of {w}", "absent", place)
    if trl.__version__ != PIN:
        msg = f"TRL {trl.__version__} is installed; every contract item passed. Pin is {PIN}."
        if os.environ.get("OBJDEL_TRL_UNPINNED") != "1":
            raise TRLContractError(msg + " Set OBJDEL_TRL_UNPINNED=1 to try this version.")
        warnings.warn(msg, UserWarning, stacklevel=2)


def check_lora(llm) -> None:
    """Under lora-seqmask (B2.2 F step 3): raise TRLContractError on the first option-F vLLM name missing from the
    live engine `llm` (vg.llm): add_lora and remove_lora callable, the adapter manager's lora_index_to_id a
    sequence, and at least one engine module holding lora_a_stacked and lora_b_stacked."""
    for name, (path, want, place) in LORA_POINTS.items():
        obj = llm
        for attr in path.split("."):
            obj = getattr(obj, attr, None)
        if name == "lora_a_stacked / lora_b_stacked":
            ok = obj is not None and any(hasattr(m, "lora_a_stacked") and hasattr(m, "lora_b_stacked")
                                         for _, m in obj.named_modules())  # fmt: skip
        elif want == "callable":
            ok = callable(obj)
        else:
            ok = obj is not None and hasattr(obj, "__len__") and hasattr(obj, "__getitem__")
        if not ok:
            _fail("option-F vLLM name", name, f"{path}: {want}", "absent" if obj is None else type(obj).__name__,
                  place)


FINE = tuple(range(64, 2049, 64))  # uncompiled: round to 64 only (BUCKETS would pad a 430-wide chunk to 512)

WAKE_PHASE = {("weights",): "wake_weights", ("kv_cache",): "wake_kv"}  # llm.wake_up(tags) -> its time/<phase>_s

ROW_KEYS = ("mm_token_type_ids", "token_type_ids")  # per-token, prompt + completion wide: cut like input_ids

OTHER_VISION = ("pixel_attention_mask", "spatial_shapes", "num_tiles", "image_sizes", "image_position_ids")


class FastForward:
    """GRPOTrainer mixin: each forward chunk runs on its own window (see the module docstring)."""

    buckets = FINE  # build_trainer sets BUCKETS when it compiles the text layers
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
            pv, ids = None, enc.check_rows(None if meta is None else meta.image_id, nimg, model)
        rows = None if pv is None else meta.pixel_rows.tolist() if meta is not None else grid_rows(grid, nimg)
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
                r = 0 if rows is None else sum(rows[s:e])
                pix = None if pv is None else list(pv[row:row + r].split(rows[s:e]))
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
                    r = sum(rows[s:e])
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


def no_lora_input_cast(model) -> int:
    """Turns off PEFT's cast of each LoRA layer's input to its fp32 adapter dtype; returns the layers changed. Only
    under bf16 autocast: there the adapter matmul casts the input straight back to bf16, so the fp32 copy (one per
    LoRA layer per forward and recompute) changes nothing; without autocast the cast is what matches the dtypes."""
    layers = [m for _, m in model.named_modules() if hasattr(m, "cast_input_dtype_enabled")]
    for m in layers:
        m.cast_input_dtype_enabled = False
    return len(layers)


def one_pad_per_image(prompts: list[dict], pad: int) -> list[dict]:
    """vLLM prompts with each run of image pads cut to one pad (rows without multi_modal_data unchanged)."""
    out = []
    for row in prompts:
        if "multi_modal_data" in row:
            ids = row["prompt_token_ids"]
            row = row | {"prompt_token_ids": [t for i, t in enumerate(ids) if t != pad or not i or ids[i - 1] != pad]}
        out.append(row)
    return out


def vllm_generate_checked(generate, prompts: list[dict], pad: int, *args, **kw):
    """generate(one_pad_per_image(prompts, pad), ...) with vLLM's own `generate`, then each multimodal row's
    expanded ids checked against the trainer's processor-expanded ids (B2.2 F step 4b.3: the patch_vllm_prompts
    wrapper and the F comparator both send prompts through here, never raw)."""
    outputs = generate(one_pad_per_image(prompts, pad), *args, **kw)
    for row, o in zip(prompts, outputs, strict=True):
        if "multi_modal_data" in row and list(o.prompt_token_ids) != list(row["prompt_token_ids"]):
            n, m = list(o.prompt_token_ids).count(pad), list(row["prompt_token_ids"]).count(pad)
            raise RuntimeError(f"vLLM expanded the image to {n} pads, the trainer's processor to {m}: "
                               "vLLM would sample from other image features than the trainer scores")
    return outputs


def prompt_images(prompts) -> list[list]:
    """Each conversational prompt's image objects, in order, as TRL extracts them (grpo_trainer.py 1760-1771)."""
    return [[c["image"] for m in p if isinstance(m["content"], list)
             for c in m["content"] if c["type"] == "image"] for p in prompts]  # fmt: skip


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


def vllm_engine_kwargs(args) -> dict:
    """The LLM(...) keywords TRL leaves out: float16 under fp16, LoRA serving under option F (args.vllm_lora)."""
    return ({"dtype": "float16"} if getattr(args, "fp16", False) else {}) | (
        {"enable_lora": True, "max_lora_rank": 16} if getattr(args, "vllm_lora", False) else {})


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
    images = None  # ImageMemo on the processing class's image processor (set in __init__ when it has one)
    host_meta = False  # C4, train_rl --host-meta: RowMeta from _generate's host lists, row-aligned in the batch
    check_meta_every = 50  # C4, train_rl --check-meta-every N: compare host_meta with the device masks (0: never)
    check_rope_every = 50  # C2, train_rl --check-rope-every N: compare template rope positions with HF (0: never)
    checks_from_step = 1  # C0, train_rl --checks-from-step S: --check-*-every count from step S; recording starts there
    record_dir = None  # C0, train_rl --record-batches DIR (record.open_dir): each generation batch from step S on
    rope_check = False  # this generation batch's forwards run the rope check (FastForward reads it)
    meta = None  # the generation batch's RowMeta, from _generate to the batch
    prompt_memo = False  # C2, train_rl --prompt-memo: _tokenize_prompts tokenizes one row per prompt layout
    count_syncs = False  # C0, train_rl --count-syncs: perf/host_syncs/<phase> per timed phase
    cuda_times = None  # OBJDEL_NVTX_SYNC=0: (mode, phase, start, end) CUDA event pairs of `timed`, read by `log`
    last_log = None  # perf_counter at the previous `log` (time/step_s)
    sync_calls = None  # weight syncs since the previous `log`, counted once time_vllm_phases wraps sync_weights
    group_requests = False  # C6, train_rl --group-requests: vLLM gets one n = G request per prompt (vllm_requests)
    opsm_watch = None  # B2.4 rule 4: callbacks.OpsmDelta, fed the train-mode per-sequence (seq_kl, advantage) rows
    log_by_prompt = False  # B2.3: train logs also hold rewards/iou_tc/by_prompt (one mean per prompt of the step)
    runaway_cap = False  # C11, train_rl --runaway-cap: max_tokens per prompt from the rows' `cap` (runaway)
    lora_dir = None  # option F: the tmpfs folder of the adapter vLLM serves; every generate carries its LoRARequest

    def __init__(self, *args, **kwargs):
        check()  # the private TRL names this class relies on, before TRL builds anything
        with vllm_kwargs(**vllm_engine_kwargs(kwargs.get("args"))):  # TRL builds its vLLM engine without these
            super().__init__(*args, **kwargs)
        self.check_engine_dtype()
        if getattr(self.args, "vllm_lora", False) and getattr(self, "vllm_generation", None) is not None:
            check_lora(self.vllm_generation.llm)  # option F's vLLM names (B2.2 F step 3)
        self.vision, self.rope = VisionCache.install(self.model), RopeCache.install(self.model)
        self.images = ImageMemo.install(self.processing_class)
        self.drop_lora_input_cast()
        self.sync_without_merge()

    def use_gpu_pixels(self, images, allow_cpu: bool = False) -> str:
        """C1 start-up self-check (train_rl --gpu-pixels auto): the first gpu_pixels backend equal to the CPU
        processor on `images` (training views) plus an in-memory digit image prepares the memo's misses on CUDA.
        None passing raises, unless allow_cpu (then the CPU processor stays: returns "cpu")."""
        from objdel.policy import gpu_pixels

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
        if self.host_meta:  # C2: check steps (check_due) check every pass of this generation batch
            self.rope_check = self.check_due(self.check_rope_every)
        try:
            with self.timed("batch"):
                batch = super()._generate_and_score_completions(inputs)
            if self.encoder is not None:  # C3: the update takes the encoder's features; TRL would cut these by grid
                batch.pop("pixel_values", None)
            if self.meta is not None:
                if self.host_meta:
                    batch["host_meta"] = host_rows(self.meta)
                    if self.check_due(self.check_meta_every):
                        check_meta(batch)
                if self.record_dir and self.model.training and self.state.global_step + 1 >= self.checks_from_step:
                    torch.save(record.batch_record(batch, self.meta, self.rope and self.rope.base),
                               self.record_dir / record.file_name(self.state.global_step + 1))
            return batch
        finally:
            self.rows, self.meta = None, None
            if skip:
                self.vllm_importance_sampling_correction = True

    def check_due(self, every: int) -> bool:
        """Whether this generation batch is a check step of a --check-*-every flag: steps S, S + every, S + 2 every,
        ... where S is --checks-from-step (default 1; every 0 never)."""
        late = self.state.global_step + 1 - self.checks_from_step  # steps since S; negative before it
        return bool(every) and late >= 0 and late % every == 0

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
        meta = self.loss_meta if self.loss_meta is not None else self.meta  # C2: the reference and old passes too
        if meta is not None:
            kwargs["host_meta"] = meta
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
        <output_dir>/seq_gaps.f32 (write_seq_gaps). On a CUDA build, train steps 1-6 add `engine/<key>` for every
        gc.stats() key and `engine/max_memory_reserved` (host reads of counters, no CUDA sync; D7 step 5)."""
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
            if self.sync_calls is not None:
                m["time/sync_calls"].append(self.sync_calls)
                self.sync_calls = 0
            if self.log_by_prompt and train:  # B2.3; reward_funcs read only when asked (fake trainers have none)
                reward = next((f.trl_reward for f in self.reward_funcs if hasattr(f, "trl_reward")), None)
                if reward is not None:
                    logs = {**logs, "rewards/iou_tc/by_prompt": reward.by_prompt}
            if train and _gpucheck.build() == "cuda" and 1 <= self.state.global_step <= 6:  # D-G14, D-G15 (D7)
                logs = logs | {f"engine/{k}": v for k, v in _gpucheck.stats().items()} | {
                    "engine/max_memory_reserved": torch.cuda.max_memory_reserved()}
            return super().log(logs, start_time)

    def train(self, *args, **kwargs):
        """Trainer.train; on a CUDA build, then one more log_history row with the run's `engine/pinned_grows`
        (D-G15 reads it over the whole run)."""
        out = super().train(*args, **kwargs)
        if _gpucheck.build() == "cuda":
            self.state.log_history.append({"engine/pinned_grows": _gpucheck.stats()["pinned_grows"],
                                           "step": self.state.global_step})
        return out

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

    def decode_tail(self):
        """C11: the NVTX `generate/tail` range of the vLLM decode loop, opened once fewer than num_generations answers
        run. Nothing without NVTX, without vLLM, or under --group-requests (one request then holds G answers)."""
        if not (nvtx.enabled() and self.use_vllm) or self.group_requests:
            return contextlib.nullcontext()
        G = self.num_generations if self.model.training else self.num_generations_eval
        return nvtx.tail(self.vllm_generation.llm.llm_engine, G)

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
        with self.timed("generate"), self.decode_tail():
            out = super()._generate(prompts)
        if (self.host_meta or self.record_dir) and self.rows is not None and len(self.rows) == len(out[1]):
            eos_and_pad = [*self.eos_token_ids, self._tokenizer.pad_token_id]  # TRL's truncation list
            cfg = getattr(self.model, "config", None)
            vision = [getattr(cfg, k, None) for k in ("image_token_id", "video_token_id", "vision_start_token_id",
                                                      "vision_end_token_id")]  # fmt: skip
            merge = getattr(getattr(cfg, "vision_config", None), "spatial_merge_size", None)
            self.meta = row_meta(out[0], out[1], image_ids(self.rows), self.num_generations,
                                 self.mask_truncated_completions, eos_and_pad,
                                 [getattr(r.get("image"), "size", None) for r in self.rows],
                                 [v for v in vision if v is not None],
                                 (vision[0], merge) if vision[0] is not None and merge else None)  # fmt: skip
        if self.encoder is not None and self.rows is not None:
            self.encode_images()
        self.prefetch_rewards(out[3], out[1])
        self.log_capped(out[1])
        return out

    def encode_images(self) -> None:
        """C3: the generation batch's distinct images (the memo's pass-2 entries, prepared here when missing) through
        the tower in one call, in image id order (RowMeta.image_id). Nothing for rows without an image."""
        imgs = list({id(r.get("image")): r.get("image") for r in self.rows}.values())
        if any(x is None for x in imgs):
            return
        px = self.images(imgs, return_tensors="pt")
        self.encoder.encode(px["pixel_values"], px["image_grid_thw"])

    def log_capped(self, completion_ids):
        """C11: completions/capped_share, the share of this batch's answers that reached their row's cap."""
        if self.runaway_cap and self.rows is not None:
            hit = [len(ids) >= r["cap"] for ids, r in zip(completion_ids, self.rows, strict=True)]
            self._metrics["train" if self.model.training else "eval"]["completions/capped_share"].append(
                sum(hit) / len(hit))

    def prefetch_rewards(self, completions, completion_ids):
        """Start scoring (objdel.rl.reward.TRLReward.prefetch) on the answers vLLM just returned, with the
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
        """time/sync_s (weight push, with its weights wake-up), time/sync_calls (pushes per step, logged by `log`),
        time/wake_weights_s and time/wake_kv_s (by the wake_up `tags`; any other call charges time/wake_s, so no wake
        time is lost) and time/sleep_s: the colocate handovers where the GPU idles between generation and training.
        Idempotent; the engine exists only after TRL's vLLM init."""
        vg = self.vllm_generation
        self.sync_calls = self.sync_calls or 0
        hooks = ((vg, "sync_weights", "sync"), (vg.llm, "wake_up", "wake"), (vg.llm, "sleep", "sleep"))
        for host, name, phase in hooks:
            fn = getattr(host, name)
            if getattr(fn, "timed_phase", None):
                continue

            def wrapped(*args, _fn=fn, _phase=phase, **kwargs):
                if _phase == "sync":
                    self.sync_calls += 1
                if _phase == "wake":
                    _phase = WAKE_PHASE.get(tuple((args[0] if args else kwargs.get("tags")) or ()), _phase)
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
            if self.lora_dir is not None:
                from vllm.lora.request import LoRARequest

                from objdel.rl.vllm_sync import LORA_NAME

                kwargs["lora_request"] = LoRARequest(LORA_NAME, 1, str(self.lora_dir))
            return vllm_generate_checked(generate, prompts, pad, *args, **kwargs)

        one_pad.one_pad = True
        send = one_pad
        if self.group_requests:  # C6: outside one_pad, so its check compares the 32 sent prompts
            plain = send

            def send(prompts, *args, **kwargs):
                G = self.num_generations if self.model.training else self.num_generations_eval
                return vllm_requests.grouped(plain, G)(prompts, *args, **kwargs)

        if self.runaway_cap:  # C11: the generation batch's rows hold the caps; outside one (g_probe) the limit stays
            sent = send

            def send(prompts, *args, sampling_params=None, **kwargs):
                if self.rows is not None:
                    if len(self.rows) != len(prompts):
                        raise ValueError(f"{len(self.rows)} caps for {len(prompts)} prompts")
                    sampling_params = runaway.per_row_params(sampling_params, [r["cap"] for r in self.rows])
                return sent(prompts, *args, sampling_params=sampling_params, **kwargs)

        send.one_pad = True
        llm.generate = send

    def comparator_prompts(self, rows: list[dict], G: int) -> list[dict]:
        """B2.2 F step 4, the 32 x 8 batch: each dataset row G times, through TRL's own prompt path (the image into
        the conversational prompt, then `_tokenize_prompts`, ImageMemo on), in TRL's vLLM form: the trainer's
        processor-expanded ids and multi_modal_data (send them only through vllm_generate_checked)."""
        rows = [r for r in rows for _ in range(G)]
        prompts = [trl_grpo.prepare_multimodal_messages(r["prompt"], images=[r["image"]]) for r in rows]
        ids, images, _ = self._tokenize_prompts(prompts)
        return [{"prompt_token_ids": list(i)} | ({"multi_modal_data": {"image": im[0] if len(im) == 1 else im}}
                if im else {}) for i, im in zip(ids, images or [None] * len(ids), strict=True)]  # fmt: skip

    @torch.no_grad()
    def way_b_gap(self, prompts: list[dict], outputs, adapter: bool) -> float:
        """B2.3b way (b): mean |vLLM logp - trainer logp| over every answer token of `outputs` (vLLM RequestOutputs of
        `prompts`), the trainer scoring with TRL's plain forward (no FastGRPO patch), the active adapter on
        (adapter True) or disabled. One forward per run of rows sharing a prompt."""
        def same_prompt(a, b):  # same ids and the same image object (comparator_prompts repeats each row G times)
            img = ("multi_modal_data" in a and a["multi_modal_data"]["image"]) is (
                "multi_modal_data" in b and b["multi_modal_data"]["image"])
            return img and a["prompt_token_ids"] == b["prompt_token_ids"]

        pad = self.model.config.image_token_id
        total = n = 0.0
        i = 0
        while i < len(prompts):
            j = i
            while j < len(prompts) and same_prompt(prompts[j], prompts[i]):
                j += 1
            ids_p = torch.tensor(prompts[i]["prompt_token_ids"])
            outs = [o.outputs[0] for o in outputs[i:j]]
            T = max(len(o.token_ids) for o in outs)
            comp = torch.full((j - i, T), self._tokenizer.pad_token_id)
            cmask = torch.zeros(j - i, T, dtype=torch.long)
            vlp = torch.zeros(j - i, T, dtype=torch.float64)
            for r, o in enumerate(outs):
                k = len(o.token_ids)
                comp[r, :k], cmask[r, :k] = torch.tensor(list(o.token_ids)), 1
                vlp[r, :k] = torch.tensor([lp[t].logprob for t, lp in zip(o.token_ids, o.logprobs, strict=True)])
            ids = torch.cat([ids_p.expand(j - i, -1), comp], 1)
            mask = torch.cat([torch.ones(j - i, len(ids_p), dtype=torch.long), cmask], 1)
            kw = {"mm_token_type_ids": (ids == pad).long()}
            image = prompts[i].get("multi_modal_data", {}).get("image")
            if image is not None:
                px = self.processing_class.image_processor(images=[image] * (j - i), return_tensors="pt")
                kw |= {"pixel_values": px["pixel_values"], "image_grid_thw": px["image_grid_thw"],
                       "num_images": [1] * (j - i)}  # fmt: skip
            dev = self.accelerator.device
            kw = {k: v.to(dev) if torch.is_tensor(v) else v for k, v in kw.items()}
            off = self.model.disable_adapter() if not adapter else contextlib.nullcontext()
            with off:
                lp = GRPOTrainer._get_per_token_logps_and_entropies(self, self.model, ids.to(dev), mask.to(dev), T,
                                                                    batch_size=j - i, **kw)[0]
            m = cmask.bool()
            total += (lp.double().cpu() - vlp).abs()[m].sum().item()
            n += int(m.sum())
            i = j
        return total / max(n, 1)

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


def sample(trainer, rows: Iterable[dict], k: int, items: int = 16) -> list[dict]:
    """k completions (a multiple of G) per row of objdel.rl.rows.grpo_rows (prompt, image, key), `items` rows
    per generation call (a 450M model on a 4090 is overhead-bound at one row per call). Each row: item_id, samples,
    truncated and tokens (each answer's completion id count, as TRL's completion length)."""
    if k <= 0 or k % G:
        raise ValueError(f"k = {k}: sample whole groups of {G}")
    stop = models.stop_ids(trainer.model, trainer.processing_class)  # the reward's rule
    rows = list(rows)
    was_training = trainer.model.training
    trainer.model.eval()  # TRL then generates as in held-out evaluation (num_generations_eval)
    out = []
    try:
        for start in range(0, len(rows), items):
            chunk = rows[start:start + items]
            prompts = [r["prompt"] if r.get("image") is None else
                       prepare_multimodal_messages(r["prompt"], images=[r["image"]]) for r in chunk]
            _, ids, _, completions, *_ = trainer._generate([p for p in prompts for _ in range(k)])
            for j, r in enumerate(chunk):
                part = slice(j * k, (j + 1) * k)
                out.append({"item_id": r["key"], "samples": [models._text(c) for c in completions[part]],
                            "truncated": [models.truncated(x, stop) for x in ids[part]],
                            "tokens": [len(x) for x in ids[part]]})
    finally:
        trainer.model.train(was_training)
    return out
