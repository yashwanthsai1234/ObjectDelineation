"""Phase 0b baseline B: zero-shot small VLM asked for disc and cup polygons on RIGA.

The image is cropped around the experts' disc (grown by half its size per side),
the same crop MedSAM gets, but with no box prompt. Coordinates are requested in
0-1000 relative to the crop, Qwen-VL's grounding convention.

    .venv/bin/python -m phase0.baseline_vlm polygon data/raw/riga data/riga_true phase0/results/vlm_polygon 30
    .venv/bin/python -m phase0.baseline_vlm box data/raw/riga data/riga_true phase0/results/vlm_box 30

"box" asks for Qwen-VL's native bbox_2d grounding; each box is scored as its inscribed ellipse.
"""
import json
import random
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
from PIL import Image
from transformers import AutoModelForImageTextToText, AutoProcessor

from checker import AnnotatorSet
from checker.polygon_io import as_points, extract_json, rasterize, ring_problems
from phase0.baseline_medsam import TRUE, score_row, write_rows
from phase0.medsam import box_of
from phase0.riga_outlines import jobs as riga_jobs

MODEL = "Qwen/Qwen3-VL-2B-Instruct"
SIDE = 448
PROMPTS = {"polygon": (
    "This is a crop of a retinal fundus photograph centred on the optic nerve head. "
    "Outline the optic disc and the optic cup. Reply with JSON only, in this form: "
    '{"disc": [[x, y], ...], "cup": [[x, y], ...]} '
    "with 12 to 24 points per outline, in order around the boundary. "
    "Coordinates are integers from 0 to 1000, relative to the image width and height."
), "box": (
    "This is a crop of a retinal fundus photograph centred on the optic nerve head. "
    "Locate the optic disc and the optic cup. Reply with JSON only, in this form: "
    '{"disc": {"bbox_2d": [x1, y1, x2, y2]}, "cup": {"bbox_2d": [x1, y1, x2, y2]}} '
    "with coordinates from 0 to 1000 relative to the image width and height."
)}


def load():
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    model = AutoModelForImageTextToText.from_pretrained(MODEL, dtype=torch.bfloat16).to(device).eval()
    return model, AutoProcessor.from_pretrained(MODEL), device


@torch.no_grad()
def ask(model, processor, device, image, prompt):
    msgs = [{"role": "user", "content": [{"type": "image", "image": image}, {"type": "text", "text": prompt}]}]
    inputs = processor.apply_chat_template(msgs, tokenize=True, add_generation_prompt=True,
                                           return_dict=True, return_tensors="pt").to(device)
    out = model.generate(**inputs, max_new_tokens=500, do_sample=False)
    return processor.batch_decode(out[:, inputs["input_ids"].shape[1]:], skip_special_tokens=True)[0]


def ellipse_points(box, n=32):
    x1, y1, x2, y2 = box
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.stack([(x1 + x2) / 2 + (x2 - x1) / 2 * np.cos(t), (y1 + y2) / 2 + (y2 - y1) / 2 * np.sin(t)], axis=1)


def outline(obj, name, mode):
    """Model output for one structure -> (points in 0-1000 or None, problems)."""
    v = obj.get(name) if isinstance(obj, dict) else None
    if mode == "box":
        b = as_points(v.get("bbox_2d")) if isinstance(v, dict) else None
        if b is None or b.size != 4:
            return None, ["unparseable"]
        x1, y1, x2, y2 = b.ravel()
        probs = ["out_of_range"] if b.min() < 0 or b.max() > 1000 else []
        if x2 <= x1 or y2 <= y1:
            return None, probs + ["empty_box"]
        return ellipse_points((x1, y1, x2, y2)), probs
    pts = as_points(v)
    return pts, ring_problems(pts, bounds=(0, 1000))


def main(mode, raw_root, true_root, out_dir, n_images):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    primes = {(sub.replace("/", "_"), idx): p for sub, idx, p, _, _ in riga_jobs(raw_root, None)}
    experts = defaultdict(dict)
    for f in Path(true_root).glob("*.png"):
        if m := TRUE.match(f.name):
            experts[(m.group(1), int(m.group(2)))][int(m.group(3))] = f
    keys = sorted(k for k in experts if k in primes and len(experts[k]) >= 2)
    keys = sorted(random.Random(0).sample(keys, min(int(n_images), len(keys))))
    model, processor, device = load()
    rows, replies = [], []
    for n, key in enumerate(keys, 1):
        ex = [np.asarray(Image.open(f)) for _, f in sorted(experts[key].items())]
        disc, cup = [e >= 100 for e in ex], [e >= 200 for e in ex]
        x0, y0, x1, y1 = box_of(np.any(disc, 0), 0)
        w, h = x1 - x0, y1 - y0
        img = np.asarray(Image.open(primes[key]).convert("RGB"))
        cy0, cx0 = max(0, y0 - h // 2), max(0, x0 - w // 2)
        cy1, cx1 = min(img.shape[0], y1 + h // 2 + 1), min(img.shape[1], x1 + w // 2 + 1)
        crop = Image.fromarray(img[cy0:cy1, cx0:cx1])
        H, W = cy1 - cy0, cx1 - cx0
        scale = SIDE / max(W, H)
        t = time.time()
        reply = ask(model, processor, device, crop.resize((round(W * scale), round(H * scale))), PROMPTS[mode])
        replies.append({"subset": key[0], "image": key[1], "reply": reply, "seconds": round(time.time() - t, 1)})
        obj = extract_json(reply) or {}
        preds, valid = {}, {}
        for name in ("disc", "cup"):
            pts, probs = outline(obj, name, mode)
            valid[name] = ";".join(probs) or "ok"
            ok = pts is not None and len(pts) >= 3
            preds[name] = rasterize(pts * np.array([W / 1000, H / 1000]), W, H) if ok else np.zeros((H, W), bool)
        cd = [d[cy0:cy1, cx0:cx1] for d in disc]
        cc = [c[cy0:cy1, cx0:cx1] for c in cup]
        for name, pred, gt in (("disc", preds["disc"], cd), ("cup", preds["cup"], cc),
                               ("rim", preds["disc"] & ~preds["cup"], [d & ~c for d, c in zip(cd, cc)])):
            rows.append({"subset": key[0], "image": key[1], "structure": name, "experts": len(gt),
                         "parsed": bool(obj), "valid": valid.get(name, f"{valid['disc']}|{valid['cup']}"),
                         **score_row(AnnotatorSet(gt, min_size=100), pred)})
        print(f"vlm {n}/{len(keys)} {replies[-1]['seconds']}s valid={valid}", flush=True)
        write_rows(out / "scores.csv", rows)
        (out / "replies.jsonl").write_text("\n".join(json.dumps(r) for r in replies))


if __name__ == "__main__":
    main(*sys.argv[1:6])
