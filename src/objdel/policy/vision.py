"""G-vision: does the model still look at the image? Two measures, both on the CPU in float32:

(a) the mean pairwise cosine distance between the projector outputs (get_image_features, what the language model
    sees; all image tokens of one image flattened into one vector) of a fixed image set, the 4 centred LIDC val crops
    plus a black image, as a ratio to the base model's value (`projector_distance(base, ...)`, measured once);
(b) the teacher-forced next-token total variation over the first T = 24 tokens of the first rater target, on 8 fixed
    LIDC val items: real crop vs a black image and real crop vs the next item's crop; per item the mean over
    positions of 0.5 * sum |p - q|, then the mean over items.

Thresholds: ratio >= 0.5 and both TVs >= 0.002 (collapsed checkpoints measured 0.001-0.064 and <= 0.0007, intact ones
0.68-1.0 and >= 0.004; the thresholds hold only under this exact protocol). `decide` is the rule of select_ckpt
(eligibility needs both) and of train_sft at each rolling checkpoint: both below threshold is "abort", one is "warn",
none is "ok". The fixed sets (`fixed_set`) are a seed-0 permutation of the dataset's RL-eligible val records, its first
POOL in item-id order, the first n of those; every image is checked against its prepared image_sha256 and the set's
SHA-256 (ids, image hashes, targets) is returned for the run manifest. torch is imported lazily.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from objdel.core.io import read_jsonl
from objdel.policy.prompt import user_turn

RATIO, TV = 0.5, 0.002
N_IMAGES, N_ITEMS, T = 4, 8, 24
POOL, SEED = 64, 0  # a seed-0 draw of 64 RL-eligible val items, then [:4] images and [:8] items


def decide(ratio: float, tv_black: float, tv_other: float) -> dict:
    """Thresholds and the in-training action (module docstring)."""
    ok_ratio, ok_tv = ratio >= RATIO, min(tv_black, tv_other) >= TV
    action = "ok" if ok_ratio and ok_tv else ("warn" if ok_ratio or ok_tv else "abort")
    return {"pass_ratio": bool(ok_ratio), "pass_tv": bool(ok_tv), "pass": action == "ok", "action": action}


def fixed_set(prepared_root, dataset: str, n: int = N_ITEMS) -> tuple[list[dict], str]:
    """([{"item_id", "image": PIL RGB, "target"}], SHA-256 of the set) of the fixed val items of `dataset` (module
    docstring); ValueError when an image is not the one recorded."""
    from PIL import Image
    root = Path(prepared_root)
    recs = [r for r in read_jsonl(root / dataset / "val.jsonl") if r["rl"]]
    idx = np.random.default_rng(SEED).permutation(len(recs))[:POOL]
    out = []
    for r in sorted((recs[i] for i in idx), key=lambda r: r["item_id"])[:n]:
        a = np.asarray(Image.open(root / "images" / f"{r['item_id']}.png").convert("RGB"))
        if hashlib.sha256(repr(a.shape).encode() + a.tobytes()).hexdigest() != r["image_sha256"]:
            raise ValueError(f"{r['item_id']}: image differs from its image_sha256")
        out.append({"item_id": r["item_id"], "image": Image.fromarray(a),
                    "target": next(t for t in r["targets"] if t), "image_sha256": r["image_sha256"]})
    key = [[it["item_id"], it["image_sha256"], it["target"]] for it in out]
    return out, hashlib.sha256(json.dumps(key).encode()).hexdigest()


def black(image):
    """An all-black RGB image of the same size."""
    from PIL import Image
    return Image.new("RGB", image.size)


def _batch(model, processor, image, prompt):
    msgs = [user_turn(prompt)]
    chat = processor.apply_chat_template(msgs, tokenize=False, add_generation_prompt=True)
    b = processor(text=[chat], images=[[image]], return_tensors="pt")
    return {k: v.to(model.device, model.dtype) if v.is_floating_point() else v.to(model.device) for k, v in b.items()}


def _features(model, processor, image, prompt):
    import torch
    m = model.get_base_model() if hasattr(model, "get_base_model") else model  # PEFT
    get = getattr(m, "get_image_features", None) or m.model.get_image_features
    kw = {k: v for k, v in _batch(model, processor, image, prompt).items()
          if k not in ("input_ids", "attention_mask", "mm_token_type_ids")}
    with torch.inference_mode():
        f = get(**kw)
    f = getattr(f, "pooler_output", f)
    f = torch.cat([x.reshape(-1) for x in f]) if isinstance(f, (list, tuple)) else f.reshape(-1)
    return f.double()


def projector_distance(model, processor, images, prompt: str) -> float:
    """Mean pairwise cosine distance of the projector outputs of images + a black image."""
    import torch
    F = [_features(model, processor, im, prompt) for im in [*images, black(images[0])]]
    d = [1.0 - float(torch.nn.functional.cosine_similarity(F[i], F[j], 0))
         for i in range(len(F)) for j in range(i + 1, len(F))]
    return float(np.mean(d))


def extend(batch: dict, ids) -> dict:
    """A one-prompt processor batch with answer ids (a (1, T) tensor) appended for a teacher-forced forward: every
    token-aligned (1, n) tensor is extended over the answer (input_ids by the ids, attention_mask by ones, the others,
    e.g. mm_token_type_ids, by zeros), each in its own dtype and device; other entries are unchanged."""
    import torch
    n, out = batch["input_ids"].shape[1], {}
    for k, v in batch.items():
        if isinstance(v, torch.Tensor) and v.shape == (1, n):
            t = ids.to(v.device)
            tail = t if k == "input_ids" else torch.ones_like(t) if k == "attention_mask" else torch.zeros_like(t)
            v = torch.cat([v, tail.to(v.dtype)], 1)
        out[k] = v
    return out


def next_token_probs(model, processor, image, target: str, prompt: str, T: int = T) -> np.ndarray:
    """(positions, vocab) float64 next-token distributions over the first T teacher-forced target tokens."""
    import torch
    b = _batch(model, processor, image, prompt)
    n = b["input_ids"].shape[1]
    tgt = processor.tokenizer(target, add_special_tokens=False, return_tensors="pt")["input_ids"][:, :T]
    b = extend(b, tgt)
    with torch.inference_mode():
        logits = model(**b).logits[0, n - 1:n - 1 + tgt.shape[1]]
    return torch.softmax(logits.double(), -1).cpu().numpy()


def image_use(model, processor, items, prompt: str, T: int = T) -> dict:
    """{"black", "other"}: mean TV between the real crop and a black image / the next item's crop."""
    tv = {"black": [], "other": []}
    for i, it in enumerate(items):
        p = next_token_probs(model, processor, it["image"], it["target"], prompt, T)
        for name, im in (("black", black(it["image"])), ("other", items[(i + 1) % len(items)]["image"])):
            q = next_token_probs(model, processor, im, it["target"], prompt, T)
            tv[name].append(float(0.5 * np.abs(p - q).sum(-1).mean()))
    return {k: float(np.mean(v)) for k, v in tv.items()}


def check(model, processor, images, items, *, reference: float, prompt: str, T: int = T) -> dict:
    """G-vision of one checkpoint: images = the fixed crops (a black image is added), items = [{"image", "target"}],
    reference = projector_distance of the base model, prompt = the dataset's prompt. -> {"ratio", "distance",
    "reference", "tv" (the smaller TV), "tv_black", "tv_other", "pass_ratio", "pass_tv", "pass", "action"}. Runs in eval
    mode, then gives every module back the mode it had (S1-9: train_sft calls it mid-training)."""
    modes = {m: m.training for m in model.modules()}
    model.eval()
    try:
        distance = projector_distance(model, processor, images, prompt)
        tv = image_use(model, processor, items, prompt, T)
    finally:
        for m, mode in modes.items():
            m.training = mode
    ratio = distance / reference
    return {"ratio": ratio, "distance": distance, "reference": reference, "tv": min(tv.values()),
            "tv_black": tv["black"], "tv_other": tv["other"], **decide(ratio, tv["black"], tv["other"])}
