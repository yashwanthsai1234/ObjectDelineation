"""H-A1-1: zero-shot box (and polygon) prior of each arm on centred and shifted val crops, greedy, no training.

    HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B1.zeroshot --arm qwen3.5-0.8b [--lidc 64]

Box prompts: "native" is the model card's grounding convention (Qwen3.5: bbox_2d, 0-1000; LFM2.5-VL: bbox, [0, 1]);
"generic" asks every model for {"bbox": [x1, y1, x2, y2]} in integers 0-1000 (SmolVLM has no convention of its own).
Parse (strict): a "bbox_2d" / "bbox" key followed by 4 numbers in the asked range with x1 < x2 and y1 < y2. Box IoU
against the box of the raters' union for the structure (disc for RIGA, nodule for LIDC), in frame-normalised
coordinates; an unparsed answer scores 0. Reference row: the box of A4's image-blind centred disc.
Polygon: the D10 prompt (objdel.train.dataset.PROMPTS), greedy, max 384 new tokens, on the first --poly items of
each dataset (centred), scored with objdel.eval.metrics (two copies of the greedy answer).
Reads untrusted data: run with python -I.

References:
- Qwen3.5-0.8B card (non-thinking default; 0-1000 grounding): https://huggingface.co/Qwen/Qwen3.5-0.8B ;
  Qwen3-VL grounding cookbook (bbox_2d JSON, relative 0-1000):
  https://github.com/QwenLM/Qwen3-VL/blob/main/cookbooks/2d_grounding.ipynb
- LFM2.5-VL-450M card (grounding prompt, coordinates normalised to [0, 1]): https://huggingface.co/LiquidAI/LFM2.5-VL-450M
"""
import argparse
import json
import re
import time

import numpy as np
import torch

from objdel.eval import metrics, summarize
from objdel.train import dataset as D

from . import common as C

QUERY = {"riga": ("optic disc", "fundus", "disc"), "lidc": ("lung nodule", "CT", "nodule")}
NUM = r"-?\d+(?:\.\d+)?"


def prompts(arm: str, ds: str) -> dict[str, tuple[str, str, float]]:
    """{name: (prompt, key, scale)}; coordinates / scale are frame-normalised."""
    q, modality, _ = QUERY[ds]
    out = {"generic": (f'Find the {q} in this {modality} image. Answer with JSON only: {{"bbox":[x1,y1,x2,y2]}}, '
                       'integers 0-1000 relative to the image.', "bbox", 1000.0)}
    if arm == "qwen3.5-0.8b":
        out["native"] = (f"Locate the {q} in this {modality} image. Report its bbox coordinates in JSON format.",
                         "bbox_2d", 1000.0)
    if arm == "lfm2.5-vl-450m":
        out["native"] = (f'Detect all instances of: {q}. Response must be a JSON array: [{{"label": ..., "bbox": '
                         '[x1, y1, x2, y2]}, ...]. Coordinates are normalized to [0,1].', "bbox", 1.0)
    return out


def parse_box(text: str, key: str, scale: float):
    m = re.search(rf'"{key}"\s*:\s*\[\s*({NUM})\s*,\s*({NUM})\s*,\s*({NUM})\s*,\s*({NUM})', text)
    if not m:
        return None
    b = [float(x) / scale for x in m.groups()]
    ok = all(0 <= v <= 1 for v in b) and b[0] < b[2] and b[1] < b[3]
    return b if ok else None


def generate(model, proc, images, prompt_texts, max_new: int, batch: int = 8) -> list[str]:
    out = []
    for i in range(0, len(images), batch):
        convs = [[{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": p}]}]
                 for p in prompt_texts[i:i + batch]]
        texts = [proc.apply_chat_template(c, add_generation_prompt=True, tokenize=False) for c in convs]
        inputs = proc(text=texts, images=[[im] for im in images[i:i + batch]], return_tensors="pt", padding=True)
        inputs = {k: v.to(model.device) if hasattr(v, "to") else v for k, v in inputs.items()}
        with torch.no_grad():
            ids = model.generate(**inputs, max_new_tokens=max_new, do_sample=False)
        out += proc.batch_decode(ids[:, inputs["input_ids"].shape[1]:], skip_special_tokens=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=list(C.ARMS))
    ap.add_argument("--lidc", type=int, default=64)
    ap.add_argument("--riga", type=int, default=64)
    ap.add_argument("--poly", type=int, default=16)
    a = ap.parse_args()
    C.cap_mps()
    from transformers import AutoModelForImageTextToText
    out_dir = C.SCRATCH / "zeroshot"
    out_dir.mkdir(parents=True, exist_ok=True)
    model = AutoModelForImageTextToText.from_pretrained(C.ARMS[a.arm]["hf"], dtype=torch.float32).to("mps").eval()
    res, raw = {"arm": a.arm}, {}
    for ds, n in (("riga", a.riga), ("lidc", a.lidc)):
        proc = C.processor(a.arm, ds)
        recs = C.eval_records(ds, n)
        V = C.views(ds, recs)
        struct = QUERY[ds][2]
        for view, items in V.items():
            gt = [C.box_of(x["caches"][struct].union, x["W"]) for x in items]
            pm = [C.box_of(C.prior_mask(ds, x["W"], struct), x["W"]) for x in items]
            prior = np.array([C.box_iou(p, g) for p, g in zip(pm, gt, strict=True)])
            res[f"{ds}/{view}/prior_box_iou_median"] = round(float(np.median(prior)), 4)
            for name, (p, key, scale) in prompts(a.arm, ds).items():
                t = time.time()
                texts = generate(model, proc, [x["image"] for x in items], [p] * len(items), 96)
                boxes = [parse_box(s, key, scale) for s in texts]
                iou = np.array([C.box_iou(b, g) if b else 0.0 for b, g in zip(boxes, gt, strict=True)])
                parsed = np.array([b is not None for b in boxes], float)
                diff = iou - prior
                res[f"{ds}/{view}/{name}"] = {
                    "items": len(items), "parse": C.ci(lambda i: parsed[i].mean(), len(items)),
                    "box_iou_median": C.ci(lambda i: np.median(iou[i]), len(items)),
                    "box_iou_mean": C.ci(lambda i: iou[i].mean(), len(items)),
                    "minus_prior_mean": C.ci(lambda i: diff[i].mean(), len(items)),
                    "win_share_vs_prior": round(float((diff > 0).mean()), 4), "s": round(time.time() - t, 1)}
                raw[f"{ds}/{view}/{name}"] = texts
                print(ds, view, name, json.dumps(res[f"{ds}/{view}/{name}"]), flush=True)
        # polygon prior (D10 prompt), centred, greedy
        items = V["centred"][:a.poly]
        t = time.time()
        texts = generate(model, proc, [x["image"] for x in items], [D.PROMPTS[ds]] * len(items), 384)
        its = [{"item_id": x["item_id"], "dataset": ds,
                **metrics([s, s], {k: c.raters for k, c in x["caches"].items()}, x["W"], x["H"], x["names"],
                          None, x["caches"])} for x, s in zip(items, texts, strict=True)]
        s = summarize(its)["datasets"][ds]
        res[f"{ds}/polygon"] = {"items": len(items), "format": s["format"], "legal": s["legal"], "errors": s["errors"],
                                "iou_tc_median": {k: v["iou_tc_median"] for k, v in s["objects"].items()},
                                "s": round(time.time() - t, 1)}
        raw[f"{ds}/polygon"] = texts
        print(ds, "polygon", json.dumps(res[f"{ds}/polygon"]), flush=True)
    res["peak_mps_gb"] = round(C.mps_gb(), 2)
    (out_dir / f"{a.arm}.json").write_text(json.dumps(res, indent=1))
    (out_dir / f"{a.arm}_raw.json").write_text(json.dumps(raw, indent=1))


if __name__ == "__main__":
    main()
