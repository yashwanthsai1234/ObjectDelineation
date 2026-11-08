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
"""

import contextlib
import math
import time

import torch
import torch.nn.functional as F
import trl.trainer.grpo_trainer as trl_grpo
from trl import GRPOTrainer

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


class FastForward:
    """GRPOTrainer mixin: each forward chunk runs on its own window (see the module docstring)."""

    buckets = FINE  # make_trainer sets BUCKETS when it compiles the text layers

    def _get_per_token_logps_and_entropies(self, model, input_ids, attention_mask, logits_to_keep, batch_size=None,
                                           compute_entropy=False, compute_aux_loss=False, **kw):  # fmt: skip
        if any(kw.get(k) is not None for k in OTHER_VISION) or (
                kw.get("pixel_values") is not None and kw.get("image_grid_thw") is None):
            return super()._get_per_token_logps_and_entropies(  # not a grid (Qwen-VL style) model: TRL's path
                model, input_ids, attention_mask, logits_to_keep, batch_size, compute_entropy, compute_aux_loss, **kw)
        n = input_ids.size(0)
        size = batch_size or n
        pv, grid, nimg = kw.pop("pixel_values", None), kw.pop("image_grid_thw", None), kw.pop("num_images", None)
        rows = grid.prod(-1).tolist() if grid is not None else None  # pixel rows per image (one host read)
        P = input_ids.size(1) - logits_to_keep
        rope = getattr(self, "rope", None)
        pos = (rope.positions(input_ids, attention_mask, kw["mm_token_type_ids"], grid, nimg)
               if rope is not None and grid is not None and kw.get("mm_token_type_ids") is not None else None)
        lps, ents, auxs = [], [], []
        img = row = 0
        for s in range(0, n, size):
            e = min(s + size, n)
            start, keep = window(attention_mask[s:e], logits_to_keep, self.buckets)
            cols = slice(start, P + keep)
            sub = {k: kw[k][s:e, cols] for k in ROW_KEYS if kw.get(k) is not None}
            if kw.get("sampling_mask") is not None:
                sub["sampling_mask"] = kw["sampling_mask"][s:e, :keep]
            if grid is not None:
                k = sum(nimg[s:e])
                r = sum(rows[img:img + k])
                sub |= {"pixel_values": pv[row:row + r], "image_grid_thw": grid[img:img + k], "num_images": nimg[s:e]}
                img, row = img + k, row + r
            if pos is not None:
                rope.pending = pos[:, s:e, cols]  # the chunk's slice, taken by the model's forward (RopeCache.inject)
            try:
                lp, ent, aux = super()._get_per_token_logps_and_entropies(
                    model, input_ids[s:e, cols], attention_mask[s:e, cols], keep, None, compute_entropy,
                    compute_aux_loss, **sub)  # fmt: skip
            finally:
                if pos is not None:
                    rope.pending = None
            lps.append(F.pad(lp, (0, logits_to_keep - keep)))
            if compute_entropy:
                ents.append(F.pad(ent, (0, logits_to_keep - keep)))
            if compute_aux_loss:
                auxs.append(aux)
        return (torch.cat(lps), torch.cat(ents) if compute_entropy else None,
                torch.stack(auxs).mean() if compute_aux_loss else None)  # fmt: skip


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
        self.store.clear()

    def keys(self, pixel_values, grid_thw) -> list[tuple]:
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
        self.store.clear()

    def inject(self, module, args, kwargs):
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

    def positions(self, input_ids, attention_mask, mm_token_type_ids, grid, nimg) -> torch.Tensor:
        """(3, rows, width) positions of these rows, computing only rows not seen in this generation batch."""
        keys = self.keys(input_ids, attention_mask, grid, nimg)
        first = [sum(nimg[:i]) for i in range(len(nimg))]
        miss = list({k: i for i, k in enumerate(keys) if k not in self.store}.values())
        if miss:  # on host copies: get_rope_index's per-row host reads are then no GPU syncs; one copy back
            g = torch.cat([grid[first[i]:first[i] + nimg[i]] for i in miss]).cpu()
            pos, _ = self.base.get_rope_index(input_ids[miss].cpu(), mm_token_type_ids=mm_token_type_ids[miss].cpu(),
                                              image_grid_thw=g, attention_mask=attention_mask[miss].cpu())  # fmt: skip
            pos = pos.to(input_ids.device)
            for j, i in enumerate(miss):
                self.store[keys[i]] = pos[:, j]
        return torch.stack([self.store[k] for k in keys], dim=1)


def is_ratio(logps, sampling, mask, mode: str, lo: float | None, hi: float | None) -> torch.Tensor:
    """TRL 1.17's vLLM importance-sampling ratio (grpo_trainer.py 2776-2824): exp of the trainer - vLLM logprob gap per
    token (`token_*`) or summed per sequence (`sequence_*`, shape (B, 1)); `*_truncate` clamps to [lo, hi], `*_mask`
    zeroes ratios outside it. Tokens vLLM could not score (NaN) get no correction."""
    diff = torch.nan_to_num((logps - sampling) * mask, nan=0.0)
    if mode.startswith("sequence"):
        diff = diff.sum(-1, keepdim=True)
    r = torch.exp(diff)
    if mode.endswith("truncate"):
        return torch.clamp(r, min=lo, max=hi)
    lo, hi = -math.inf if lo is None else lo, math.inf if hi is None else hi
    return r.masked_fill((r < lo) | (r > hi), 0.0)


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
    indexed, lists reordered). One host read for the order."""
    n = batch["prompt_mask"].sum(1) + batch["completion_mask"].sum(1)
    order = torch.sort(n, descending=True, stable=True).indices
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


class ImageMemo:
    """The processor's image step once per distinct image per generation batch (perf: images prepared 2-3 times).

    TRL runs the image processor over every row twice a step: for the prompt ids sent to vLLM (`_tokenize_prompts`,
    whose pixel values it then drops) and for the update's pixel values. A step has 256 rows but 32 distinct images,
    and dataset.OpenOnce hands both passes the same object for the 8 rows of an image, so this memo (keyed by the object
    and the call's kwargs, holding the object so its id stays unique) prepares each image once and assembles each
    pass's output from the per-image slices, in row order. FastGRPO clears it before each generation. Every other
    attribute is the wrapped processor's."""

    KEYS = {"pixel_values", "image_grid_thw"}  # Qwen-VL image processors: patch rows and one (t, h, w) row per image

    def __init__(self, inner):
        self.inner = inner
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
        self.memo.clear()

    def __call__(self, images, **kwargs):
        from transformers.feature_extraction_utils import BatchFeature
        from transformers.image_utils import make_flat_list_of_images

        flat, kw = make_flat_list_of_images(images), repr(sorted(kwargs.items()))
        miss = list({id(x): x for x in flat if (id(x), kw) not in self.memo}.values())
        if miss:
            out = self.inner(miss, **kwargs)
            if set(out.keys()) != self.KEYS:  # another processor family: no memo
                return self.inner(images, **kwargs)
            grid = out["image_grid_thw"]
            rows = torch.split(out["pixel_values"], grid.prod(-1).tolist())
            for x, r, g in zip(miss, rows, grid, strict=True):
                self.memo[(id(x), kw)] = (x, r, g)
        parts = [self.memo[(id(x), kw)] for x in flat]
        data = {"pixel_values": torch.cat([r for _, r, _ in parts]),
                "image_grid_thw": torch.stack([g for *_, g in parts])}
        return BatchFeature(data, tensor_type=kwargs.get("return_tensors"))


class FastGRPO(FastForward, GRPOTrainer):
    """GRPOTrainer with the per-chunk windows, the vision cache (cleared per generation batch), rewards scored while
    the old and reference logprobs run (C9), phase timers and no FLOP count."""

    rows = None  # the generation batch's dataset rows; None outside one (g_probe's step-0 sampling calls _generate)
    loss_inputs = None  # the micro-batch inside _compute_loss (C10: its forward sets the IS ratio)
    is_stats = None  # C10: the generation batch's IS-metric sums, counts, min and max, on the device until `log`
    vision = rope = None  # VisionCache, RopeCache (set in __init__ when the model has a Qwen-VL tower and base)
    images = None  # ImageMemo on the processing class's image processor (set in __init__ when it has one)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.vision, self.rope = VisionCache.install(self.model), RopeCache.install(self.model)
        self.images = ImageMemo.install(self.processing_class)
        self.drop_lora_input_cast()
        self.sync_without_merge()

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

    def drop_lora_input_cast(self) -> bool:
        """no_lora_input_cast when TRL trains in bf16 (Accelerate then runs every forward under bf16 autocast)."""
        if not getattr(self.args, "bf16", False):
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
        for cache in (self.vision, self.rope):
            if cache is not None:
                cache.clear()
        self.rows = inputs
        skip = self.skip_old_logps()
        if skip:  # TRL then skips the old-logprob pass and its IS ratio; the update's forward sets it (C10)
            self.vllm_importance_sampling_correction = False
        try:
            with self.timed("batch"):
                return super()._generate_and_score_completions(inputs)
        finally:
            self.rows = None
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
        if not (self.model.training and self.loss_type == "dr_grpo"
                and a.gradient_accumulation_steps % a.steps_per_generation == 0):  # fmt: skip
            return super()._prepare_inputs(generation_batch)
        shuffle, trl_grpo.shuffle_sequence_dict = trl_grpo.shuffle_sequence_dict, by_length
        try:
            return super()._prepare_inputs(generation_batch)
        finally:
            trl_grpo.shuffle_sequence_dict = shuffle

    def _compute_loss(self, model, inputs):
        if "importance_sampling_ratio" not in inputs and self.use_vllm and self.vllm_importance_sampling_correction:
            self.loss_inputs = inputs
        try:
            return super()._compute_loss(model, inputs)
        finally:
            self.loss_inputs = None

    def _get_per_token_logps_and_entropies(self, *args, **kwargs):
        out = super()._get_per_token_logps_and_entropies(*args, **kwargs)
        if self.loss_inputs is not None:
            self.set_is_ratio(out[0].detach(), self.loss_inputs)
        return out

    def set_is_ratio(self, logps, inputs):
        """inputs["importance_sampling_ratio"] from the update's own logprobs. TRL's sampling metrics are kept on the
        device over the generation batch's micro-batches (no host read here) and logged once by `log`."""
        sampling = inputs.get("sampling_per_token_logps")
        if sampling is None:
            raise RuntimeError("C10 needs vLLM's sampling logprobs (vllm_generation logprobs=0) for the IS ratio")
        mask = (inputs["completion_mask"] * inputs.get("tool_mask", 1)).bool()
        mode = self.vllm_importance_sampling_mode
        r = is_ratio(logps, sampling, mask, mode, self.vllm_importance_sampling_clip_min,
                     self.vllm_importance_sampling_clip_max)  # fmt: skip
        inputs["importance_sampling_ratio"] = r
        ok = mask & ~torch.isnan(sampling)
        d = (logps - sampling).abs()
        rk = torch.ones_like(r, dtype=torch.bool) if mode.startswith("sequence") else mask
        new = torch.stack([d.masked_fill(~ok, 0).sum(), ok.sum(), d.masked_fill(~ok, -math.inf).amax(),
                           r.masked_fill(~rk, 0).sum(), rk.sum(), r.masked_fill(~rk, math.inf).amin(),
                           r.masked_fill(~rk, -math.inf).amax()]).float()  # fmt: skip
        old = self.is_stats
        if old is not None:  # sums and counts add; min and max combine
            new = torch.stack([new[0] + old[0], new[1] + old[1], torch.maximum(new[2], old[2]), new[3] + old[3],
                               new[4] + old[4], torch.minimum(new[5], old[5]), torch.maximum(new[6], old[6])])
        self.is_stats = new

    def log(self, logs, start_time=None):
        """TRL's sampling metrics of the generation batch, from the device in one read (C10), then TRL's log."""
        if self.is_stats is not None:
            d_sum, d_n, d_max, r_sum, r_n, r_min, r_max = self.is_stats.tolist()
            self.is_stats = None
            m = self._metrics["train" if self.model.training else "eval"]
            for k, v in (("sampling_logp_difference/mean", d_sum / d_n if d_n else 0.0),
                         ("sampling_logp_difference/max", d_max if d_n else 0.0),
                         ("importance_sampling_ratio/min", r_min if r_n else 0.0),
                         ("importance_sampling_ratio/mean", r_sum / r_n if r_n else 0.0),
                         ("importance_sampling_ratio/max", r_max if r_n else 0.0)):  # fmt: skip
                m[f"sampling/{k}"].append(v)
        return super().log(logs, start_time)

    def _generate(self, prompts):
        if self.images is not None:  # the previous batch's images
            self.images.clear()
        if self.use_vllm:
            self.patch_vllm_prompts()
            self.time_vllm_phases()
        with self.timed("generate"):
            out = super()._generate(prompts)
        self.prefetch_rewards(out[3], out[1])
        return out

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
        """time/<phase>_s in TRL's metrics (logged as the mean per logging step). One device sync at the end, so
        queued GPU work is charged to the phase that queued it (a few per step: negligible)."""
        t = time.perf_counter()
        yield
        if torch.cuda.is_available():
            torch.cuda.synchronize()
        self._metrics["train" if self.model.training else "eval"][f"time/{phase}_s"].append(time.perf_counter() - t)

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
        processor-expanded ids (grpo_trainer.py 1844), so vLLM's expansion must give exactly those ids, or vLLM
        samples from other image features than the trainer scores (off-policy, silently): checked per prompt.
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

    def floating_point_ops(self, inputs):
        return 0  # perf F3: the Trainer default walks ~2.5k PEFT modules per micro-batch


def lora_merged_params(model):
    """(vLLM name, tensor) for every weight of a PEFT model, as TRL 1.17's VLLMGeneration._iter_named_params names
    them, with each LoRA layer's delta added into a new tensor. TRL instead runs merge_adapter / unmerge_adapter on the
    bf16 base weights in place, and W + d - d != W in bf16, so the frozen base drifts at every sync (trl #6688,
    #7423). Same values pushed (W + delta in W's dtype, as PEFT's merge computes it); no parameter is written."""
    from peft.tuners.lora import LoraLayer

    deltas = {}
    for name, m in model.named_modules():
        if isinstance(m, LoraLayer) and not m.merged:
            ads = [a for a in m.active_adapters if a in m.lora_A]
            if ads and (getattr(m, "use_dora", {}).get(ads[0]) or getattr(m, "lora_bias", {}).get(ads[0])):
                raise NotImplementedError(f"{name}: DoRA / LoRA bias are not merged here")
            if ads:
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
