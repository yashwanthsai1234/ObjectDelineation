"""Per-batch caches of the GRPO forwards (perf C2, C3, C6; from the old train/fast.py, E5): VisionCache,
EncoderCache, RopeCache, RopeTemplate, ImageMemo, and RowMeta (the host metadata of a batch's rows)."""

import dataclasses
import typing

import numpy as np
import torch

from objdel.policy import marks as nvtx


@dataclasses.dataclass(frozen=True)
class RowMeta:
    """Host facts of a generation batch's rows (C4), int32 arrays of shape (B,): real prompt tokens, real answer
    tokens (0 for a truncated row under mask_truncated_completions, as TRL zeroes its completion_mask), the row's
    distinct image (first-appearance index), its group of G answers, its prompt layout (dense id of prompt ids and
    image size: one rope template each, C2), the host's rope routing (1: the HF path, rope_bad_rows) and its pixel
    rows (its images' t*h*w: image pads x merge^2, C4 step 3)."""

    prompt_len: np.ndarray
    completion_len: np.ndarray
    image_id: np.ndarray
    group: np.ndarray
    layout: np.ndarray
    rope_bad: np.ndarray
    pixel_rows: np.ndarray


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
        its first row's real prompt), rope_positions for the rows meta marks good, HF get_rope_index for the
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
                g = torch.cat([grid[first[i]:first[i + 1]] for i in rows]).cpu()
                hf = self.base.get_rope_index(input_ids[r].cpu(), mm_token_type_ids=mm_token_type_ids[r].cpu(),
                                              image_grid_thw=g, attention_mask=attention_mask[r].cpu())[0]
                hf = hf.to(input_ids.device)  # on host copies, like rope_template: no GPU sync per row
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
        """Drop the cached images. The counters run on (a generation batch adds 32 misses and 256 hits: read deltas)."""
        self.memo.clear()

    def prefill(self, images, **kwargs) -> None:
        """Prepare the distinct images of a generation batch before TRL's processor calls (C1: on the GPU); counts
        one miss per image not yet in the memo and no hit (an empty list, or images already there: nothing)."""
        from transformers.image_utils import make_flat_list_of_images

        kw = repr(sorted(kwargs.items()))
        miss = self.missing(make_flat_list_of_images(images) if images else [], kw)
        self.misses += len(miss)
        if miss:
            self.store(miss, kw, kwargs)

    def missing(self, flat: list, kw: str) -> list:
        """The distinct images of `flat` that the memo does not hold for these kwargs (kw: their repr)."""
        return list({id(x): x for x in flat if (id(x), kw) not in self.memo}.values())

    def store(self, miss: list, kw: str, kwargs: dict) -> bool:
        """Prepare the missed images into the memo; False (nothing stored) for another processor family."""
        out = self.prepare(miss, kwargs)
        if out is None:
            return False
        pv, grid = out
        for x, r, g in zip(miss, torch.split(pv, grid.prod(-1).tolist()), grid, strict=True):
            self.memo[(id(x), kw)] = (x, r, g)
        return True

    def prepare(self, miss: list, kwargs: dict):
        """(pixel_values, image_grid_thw) of the missed images: gpu_pixels.patchify when self.gpu is set and every
        call kwarg is the processor's own value (else the processor: a counted fallback), or the processor."""
        if self.gpu is not None:
            from objdel.policy import gpu_pixels

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
            miss = self.missing(flat, kw)
            self.misses, self.hits = self.misses + len(miss), self.hits + len(flat) - len(miss)
            if miss and not self.store(miss, kw, kwargs):  # another processor family: no memo
                return self.inner(images, **kwargs)
            parts = [self.memo[(id(x), kw)] for x in flat]
            once = {id(x): p for x, p in zip(flat, parts, strict=True)}.values() if self.per_image else parts
            data = {"pixel_values": torch.cat([r for _, r, _ in once]),
                    "image_grid_thw": torch.stack([g for *_, g in parts])}
            return BatchFeature(data, tensor_type=kwargs.get("return_tensors"))
