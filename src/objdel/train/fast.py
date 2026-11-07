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
import time

import torch
import torch.nn.functional as F
from trl import GRPOTrainer

BUCKETS = (256, 384, 512, 768, 1024, 1280, 1536, 2048)  # total prompt + completion widths
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

    buckets = BUCKETS

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
            lp, ent, aux = super()._get_per_token_logps_and_entropies(
                model, input_ids[s:e, cols], attention_mask[s:e, cols], keep, None, compute_entropy, compute_aux_loss,
                **sub)  # fmt: skip
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


def one_pad_per_image(prompts: list[dict], pad: int) -> list[dict]:
    """vLLM prompts with each run of image pads cut to one pad (rows without multi_modal_data unchanged)."""
    out = []
    for row in prompts:
        if "multi_modal_data" in row:
            ids = row["prompt_token_ids"]
            row = row | {"prompt_token_ids": [t for i, t in enumerate(ids) if t != pad or not i or ids[i - 1] != pad]}
        out.append(row)
    return out


class FastGRPO(FastForward, GRPOTrainer):
    """GRPOTrainer with the per-chunk windows, the vision cache (cleared per generation batch), rewards scored while
    the old and reference logprobs run (C9), phase timers and no FLOP count."""

    rows = None  # the generation batch's dataset rows; None outside one (g_probe's step-0 sampling calls _generate)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.vision = VisionCache.install(self.model)

    def _generate_and_score_completions(self, inputs):
        if self.vision is not None:
            self.vision.clear()
        self.rows = inputs
        try:
            with self.timed("batch"):
                return super()._generate_and_score_completions(inputs)
        finally:
            self.rows = None

    def _generate(self, prompts):
        if self.use_vllm:
            self.patch_vllm_prompts()
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
