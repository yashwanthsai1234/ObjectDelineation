"""B1 shared pieces: the three arms, fixed eval items, centred and shifted views, the image-blind prior, K-sample
generation through TRL's own GRPOTrainer path, and item-bootstrap statistics.

Arms (local config, because objdel.train.configs knows only SmolVLM-256M, LFM2.5-VL-450M and Qwen3-VL-2B): Qwen3.5-0.8B
reuses Qwen3-VL-2B's settings (processor defaults plus left padding, the configured LoRA adapter); the other two use
their own configs entries. Every arm trains the same adapter (configs.lora: r16, alpha 32, all linear layers).

Views (coordinator, after A4 findings 1-3): "centred" is targets.frame (the crop centred on the raters' union);
"shifted" moves the same S x S window by (dx, dy), drawn uniformly so that the union stays inside the window (A4's
null_prior_jitter procedure), one fixed draw per item (seed, item index). Image and masks are cut from the store item
exactly as targets.frame does (lung window, channels-first crop, Pillow bicubic upscale to image_side for LIDC); at
(0, 0) the view is asserted equal to targets.frame and to the prepared PNG's pixel hash.

Prior: A4's image-blind centred disc (diameter = train median for the frame side S, per structure;
experiments/phase2/lab/A4/null_prior.json), scored with objdel.reward.score_mask on each view.

References:
- A4 null prior and shift: experiments/phase2/lab/A4/null_prior.py, null_prior_jitter.py; A4_eval_data.md §0.1
- TRL 1.16 GRPOTrainer._generate (HF generate path, GenerationConfig from GRPOConfig: temperature, top_p 1, top_k 0):
  https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/grpo_trainer.py
- Bootstrap percentile CIs: B. Efron, R. Tibshirani, An Introduction to the Bootstrap (1993), ch. 13
"""
import json
import re
import time
from pathlib import Path

import numpy as np
import torch
from PIL import Image

from objdel.checker import as_mask
from objdel.data import store
from objdel.reward import ExpertCache, score_mask
from objdel.targets import frame
from objdel.targets.frame import FRAMES, LUNG_WINDOW, _crop
from objdel.train import configs, rewards
from objdel.train import dataset as D

ROOT = Path(__file__).resolve().parents[3]
SCRATCH = ROOT / "runs/scratch/phase2/B1"  # was a temporary folder of the session that ran it
PREPARED, STORE = ROOT / "data/prepared", ROOT / "data/store"
PRIOR = json.loads((ROOT / "experiments/phase2/lab/A4/null_prior.json").read_text())
ARMS = {
    "smolvlm-256m": {"hf": "HuggingFaceTB/SmolVLM-256M-Instruct", "cfg": "smolvlm-256m"},
    "lfm2.5-vl-450m": {"hf": "LiquidAI/LFM2.5-VL-450M", "cfg": "lfm2.5-vl-450m"},
    "qwen3.5-0.8b": {"hf": "Qwen/Qwen3.5-0.8B", "cfg": "qwen3-vl-2b"},  # local config: Qwen3-VL-2B's settings
}
LORA = configs.lora("qwen3-vl-2b")  # the configured adapter, for every arm
VISION = r".*(visual|vision_model|vision_tower|connector|multi_modal_projector).*"  # vision tower + projector


def lora(scope: str = "all"):
    """The configured adapter ("all": all linear layers, vision included), or the same on the language model only
    ("lm": vision tower and projector frozen, via PEFT exclude_modules)."""
    if scope == "all":
        return LORA
    from peft import LoraConfig
    return LoraConfig(r=LORA.r, lora_alpha=LORA.lora_alpha, target_modules="all-linear", exclude_modules=VISION,
                      task_type="CAUSAL_LM")
MPS_FRACTION = 0.55  # of recommended_max_memory (12.7 GB on this 16 GB Mac): 7.0 GB


def cap_mps(fraction=MPS_FRACTION):
    if torch.backends.mps.is_available():
        torch.mps.set_per_process_memory_fraction(fraction)


def mps_gb() -> float:
    return torch.mps.driver_allocated_memory() / 1e9 if torch.backends.mps.is_available() else 0.0


def processor(arm: str, dataset: str):
    from transformers import AutoProcessor
    a = ARMS[arm]
    return AutoProcessor.from_pretrained(a["hf"], **configs.processor_kwargs(a["cfg"], dataset))


def eval_records(dataset: str, n: int, seed: int = 0) -> list[dict]:
    """n RL-eligible val records (>= 2 raters), a seeded random subset, in item_id order."""
    recs = [r for r in D.records(PREPARED, dataset, "val") if r["rl"]]
    idx = np.random.default_rng(seed).permutation(len(recs))[:n]
    return sorted((recs[i] for i in idx), key=lambda r: r["item_id"])


def _view(it, dx: int, dy: int) -> tuple[np.ndarray, dict, int]:
    """targets.frame's image and masks with the window moved by (dx, dy) (scale-1 frames only)."""
    fr = frame(it)
    if fr.scale != 1:
        raise ValueError(f"{it.item_id}: scale {fr.scale}")
    S, cfg = fr.W, FRAMES[it.dataset]
    x0, y0 = fr.box[0] + dx, fr.box[1] + dy
    masks = {n: _crop(as_mask(s, 3), x0, y0, S) for n, s in it.masks.items()}
    image = it.image
    if image.dtype != np.uint8:
        lo, hi = LUNG_WINDOW
        image = ((np.clip(image.astype(np.int64), lo, hi) - lo) * 255 + (hi - lo) // 2) // (hi - lo)
    image = _crop(np.atleast_3d(image)[..., :3].transpose(2, 0, 1), x0, y0, S)
    image = np.ascontiguousarray(np.repeat(image, 3 // len(image), 0).transpose(1, 2, 0), np.uint8)
    if cfg.image_side and cfg.image_side != S:
        image = np.array(Image.fromarray(image).resize((cfg.image_side,) * 2, Image.Resampling.BICUBIC))
    if dx == dy == 0:
        assert np.array_equal(image, fr.image) and all(np.array_equal(masks[n], fr.masks[n]) for n in masks)
    return image, masks, S


def _shift(it, rng) -> tuple[int, int]:
    """A4's draw: uniform window offset that keeps every rater's union inside the S x S window."""
    fr = frame(it)
    S, (x0, y0, _, _) = fr.W, fr.box
    u = np.logical_or.reduce([as_mask(s, 3).any(0) for s in it.masks.values()])
    ys, xs = np.nonzero(u)
    dx = int(rng.integers(-(x0 + S - xs.max() - 1), xs.min() - x0 + 1))
    dy = int(rng.integers(-(y0 + S - ys.max() - 1), ys.min() - y0 + 1))
    return dx, dy


def views(dataset: str, recs: list[dict], seed: int = 0) -> dict[str, list[dict]]:
    """{"centred": [...], "shifted": [...]}: per record, item_id, PIL image, W, H, names, caches, shift."""
    want = {r["item_id"]: (i, r) for i, r in enumerate(recs)}
    out = {"centred": [None] * len(recs), "shifted": [None] * len(recs)}
    for it in store.load(dataset, "val", STORE):
        if it.item_id not in want:
            continue
        i, rec = want[it.item_id]
        for v, (dx, dy) in (("centred", (0, 0)), ("shifted", _shift(it, np.random.default_rng([seed, i])))):
            image, masks, S = _view(it, dx, dy)
            if v == "centred" and D._pixels_sha(image) != rec["image_sha256"]:
                raise ValueError(f"{it.item_id}: centred view differs from the prepared image")
            out[v][i] = {"item_id": it.item_id, "image": Image.fromarray(image), "W": S, "H": S,
                         "names": list(masks), "shift": [dx, dy],
                         "caches": {n: ExpertCache.build(m) for n, m in masks.items()}}
    missing = [r["item_id"] for r, x in zip(recs, out["centred"], strict=True) if x is None]
    if missing:
        raise ValueError(f"not in the store: {missing[:3]}")
    return out


def prior_mask(dataset: str, S: int, name: str) -> np.ndarray:
    """A4's image-blind disc at the frame centre, diameter = train median for frame side S."""
    by = {int(k): v for k, v in PRIOR[dataset]["prior_diam_by_S"].items()}
    d = (by.get(S) or by[min(by, key=lambda s: abs(s - S))])[name]
    yy, xx = np.mgrid[:S, :S] + 0.5
    return (xx - S / 2) ** 2 + (yy - S / 2) ** 2 <= (d / 2) ** 2


def prior_iou(dataset: str, view: dict) -> dict[str, float]:
    return {n: float(score_mask(prior_mask(dataset, view["W"], n), view["caches"][n]).iou_tc) for n in view["names"]}


def box_of(mask: np.ndarray, S: int) -> list[float] | None:
    """[x1, y1, x2, y2] of a mask, normalised to [0, 1] of the frame (pixel edges)."""
    ys, xs = np.nonzero(mask)
    return None if not ys.size else [xs.min() / S, ys.min() / S, (xs.max() + 1) / S, (ys.max() + 1) / S]


def box_iou(a, b) -> float:
    ix = max(0.0, min(a[2], b[2]) - max(a[0], b[0]))
    iy = max(0.0, min(a[3], b[3]) - max(a[1], b[1]))
    inter = ix * iy
    union = (a[2] - a[0]) * (a[3] - a[1]) + (b[2] - b[0]) * (b[3] - b[1]) - inter
    return inter / union if union > 0 else 0.0


ZERO_LED = re.compile(r"^0\d+$")


def zero_led_ids(tokenizer) -> set[int]:
    """Token ids whose text is a multi-digit number starting with 0 (LFM's '042'; none on digit tokenizers)."""
    return {i for t, i in tokenizer.get_vocab().items() if ZERO_LED.match(t)}


def sample_k(trainer, rows: list[dict], k: int, zero_ids: set[int]) -> list[dict]:
    """objdel.train.sample.sample's loop (TRL GRPOTrainer._generate, eval mode, whole groups of G), keeping the
    completion ids for token counts and zero-led number tokens."""
    from trl.data_utils import prepare_multimodal_messages
    stop = rewards.stop_ids(trainer.model, trainer.processing_class)
    trainer.model.eval()
    out = []
    for r in rows:
        prompt = prepare_multimodal_messages(r["prompt"], images=[r["image"]])
        t = time.time()
        _, ids, _, completions, *_ = trainer._generate([prompt] * k)
        out.append({"item_id": r["key"], "samples": [rewards._text(c) for c in completions],
                    "truncated": [rewards.truncated(x, stop) for x in ids],
                    "tokens": [len(x) for x in ids], "zero_led": [sum(i in zero_ids for i in x) for x in ids],
                    "s": round(time.time() - t, 2)})
    return out


def eval_trainer(model, proc, arm: str, out_dir: str, max_completion: int):
    """A GRPOTrainer used only for its generation path (as experiments/phase1/sample_eval.py), HF generate."""
    from datasets import Dataset
    from trl import GRPOTrainer
    dummy = Dataset.from_list([{"prompt": [{"role": "user", "content": "x"}]}] * 8)
    args = configs.grpo_config(ARMS[arm]["cfg"], out_dir, reward_weights=None, log_completions=False, bf16=False,
                               use_vllm=False, max_completion_length=max_completion, per_device_eval_batch_size=8,
                               per_device_train_batch_size=8, gradient_accumulation_steps=1, report_to="none",
                               beta=0.0)  # no reference model: generation only
    return GRPOTrainer(model=model, reward_funcs=lambda completions, **_: [0.0] * len(completions), args=args,
                       train_dataset=dummy, processing_class=proc)


def boot(stat, n: int, B: int = 2000, seed: int = 0) -> np.ndarray:
    """B bootstrap replicates of stat(idx) over n items (idx resampled with replacement)."""
    rng = np.random.default_rng(seed)
    return np.array([stat(rng.integers(0, n, n)) for _ in range(B)])


def ci(stat, n: int, B: int = 2000) -> list[float]:
    """[point, lo, hi]: the statistic on all items and its 95 % percentile bootstrap interval."""
    reps = boot(stat, n, B)
    return [round(float(stat(np.arange(n))), 4), round(float(np.quantile(reps, 0.025)), 4),
            round(float(np.quantile(reps, 0.975)), 4)]
