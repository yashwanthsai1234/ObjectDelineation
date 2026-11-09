"""Test 2: zero-shot optic-disc bbox + polygon on jittered RIGA crops, scored against expert masks.

    HF_HOME=./hf HF_HUB_OFFLINE=1 .venv/bin/python -I t2_zeroshot.py <model> [n_images=5] [device=mps] [nested|simple]
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import torch  # noqa: E402

from common import (HERE, PROMPT_POLY, Timer, box_iou, mask_bbox, model_path, mps_mem_gb,  # noqa: E402
                    peak_rss_gb, riga_samples)
from objdel.checker.polygon_io import as_points, extract_json, rasterize, ring_problems  # noqa: E402
from objdel.reward import AnnotatorSet  # noqa: E402

BOX_PROMPT = {
    "qwen": ('Locate the optic disc in this retinal fundus image. Reply with JSON only: '
             '{"disc": {"bbox_2d": [x1, y1, x2, y2]}} with coordinates from 0 to 1000 relative to the image width and height.'),
    "lfm": ('Detect all instances of: optic disc. Response must be a JSON array: '
            '[{"label": ..., "bbox": [x1, y1, x2, y2]}, ...]. Coordinates are normalized to [0,1].'),
}


def load(name, device):
    from transformers import AutoModelForImageTextToText, AutoProcessor
    path = model_path(name)
    dtype = torch.bfloat16 if name in ("qwen3-vl-2b", "qwen3.5-0.8b") else torch.float32
    kw = {"size": {"longest_edge": 512}} if name.startswith("smolvlm") else {}
    proc = AutoProcessor.from_pretrained(path, **kw)
    model = AutoModelForImageTextToText.from_pretrained(path, dtype=dtype).to(device).eval()
    return model, proc


@torch.no_grad()
def ask(model, proc, device, image, prompt, max_new_tokens):
    msgs = [{"role": "user", "content": [{"type": "image", "image": image}, {"type": "text", "text": prompt}]}]
    kw = {"enable_thinking": False}
    inputs = proc.apply_chat_template(msgs, tokenize=True, add_generation_prompt=True, return_dict=True,
                                      return_tensors="pt", **kw).to(device)
    out = model.generate(**inputs, max_new_tokens=max_new_tokens, do_sample=False, repetition_penalty=1.05)
    return proc.batch_decode(out[:, inputs["input_ids"].shape[1]:], skip_special_tokens=True)[0]


def _walk(o):
    if isinstance(o, dict):
        for k, v in o.items():
            yield k, v
            yield from _walk(v)
    elif isinstance(o, list):
        for v in o:
            yield from _walk(v)


def to_1000(a):
    a = np.asarray(a, dtype=float)
    return a * 1000 if a.size and np.nanmax(np.abs(a)) <= 1.5 else a


def parse_box(reply):
    obj = extract_json(reply)
    cands = []
    if obj is not None:
        cands = [v for k, v in _walk(obj) if k in ("bbox_2d", "bbox", "box")]
    if not cands:  # bare list of 4 numbers
        import re
        m = re.search(r"\[\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)\s*\]", reply)
        if m:
            cands = [[float(g) for g in m.groups()]]
    for c in cands:
        p = as_points(c)
        if p is not None and p.size == 4:
            return to_1000(p.ravel())
    return None


def parse_polygon(reply):
    obj = extract_json(reply)
    if obj is None:
        return None
    for k, v in _walk(obj):
        if k in ("outer", "polygon", "disc", "points"):
            p = as_points(v)
            if p is not None and len(p) >= 3:
                return to_1000(p)
            if isinstance(v, list) and v and isinstance(v[0], dict) and "outer" in v[0]:
                p = as_points(v[0]["outer"])
                if p is not None and len(p) >= 3:
                    return to_1000(p)
    return None


PROMPT_POLY_SIMPLE = ("Outline the optic disc in this retinal fundus image. Reply with JSON only, in this form: "
                      '{"disc": [[x, y], ...]} with 12 to 24 points in order around the boundary. '
                      "Coordinates are integers from 0 to 1000, relative to the image width and height.")


def main(name, n=5, device="mps", poly_prompt="nested"):
    n = int(n)
    prompt_poly = PROMPT_POLY if poly_prompt == "nested" else PROMPT_POLY_SIMPLE
    samples = riga_samples(n, seed=0)
    with Timer() as tl:
        model, proc = load(name, device)
    print(f"{name}: loaded in {tl.s:.1f}s rss={peak_rss_gb():.2f}GB mps={mps_mem_gb():.2f}GB", flush=True)
    fam = "lfm" if name.startswith("lfm") else "qwen"
    rows = []
    for s in samples:
        side = s["image"].size[0]
        aset = AnnotatorSet(s["disc_masks"])
        vote = aset.vote
        gt_box = np.array(mask_bbox(vote)) / (side - 1) * 1000
        ys, xs = np.nonzero(vote)
        cxy = (xs.mean() / (side - 1) * 1000, ys.mean() / (side - 1) * 1000)
        row = {"key": s["key"], "gt_box": gt_box.round().tolist()}
        t = time.time()
        rb = ask(model, proc, device, s["image"], BOX_PROMPT[fam], 96)
        row["box_s"] = round(time.time() - t, 1)
        b = parse_box(rb)
        row["box_reply"] = rb[:300]
        row["box_parsed"] = b is not None
        if b is not None:
            x1, y1, x2, y2 = b
            row["box"] = b.round().tolist()
            row["box_iou"] = round(box_iou([x1, y1, x2, y2], gt_box), 3)
            row["box_has_centroid"] = bool(x1 <= cxy[0] <= x2 and y1 <= cxy[1] <= y2)
        t = time.time()
        rp = ask(model, proc, device, s["image"], prompt_poly, 400)
        row["poly_s"] = round(time.time() - t, 1)
        p = parse_polygon(rp)
        row["poly_reply"] = rp[:400]
        row["poly_parsed"] = p is not None
        if p is not None:
            row["poly_n"] = len(p)
            row["poly_problems"] = ";".join(ring_problems(p, bounds=(0, 1000))) or "ok"
            m = rasterize(p * (side - 1) / 1000, side, side)
            row["poly_iou"] = round(aset.score(m).iou, 3)
            row["poly_has_centroid"] = bool(m[int(round(cxy[1] * (side - 1) / 1000)), int(round(cxy[0] * (side - 1) / 1000))])
        rows.append(row)
        print(json.dumps({k: v for k, v in row.items() if not k.endswith("reply")}), flush=True)
        print("   box reply:", repr(rb[:160]), "\n   poly reply:", repr(rp[:200]), flush=True)
    # trivial baseline: box predicted at the expected centre of a jittered crop
    centre = [333, 333, 667, 667]
    base = np.mean([box_iou(centre, r["gt_box"]) for r in rows])
    agg = {
        "model": name, "poly_prompt": poly_prompt, "n": n, "device": device, "load_s": round(tl.s, 1),
        "box_parse_rate": np.mean([r["box_parsed"] for r in rows]),
        "box_iou_mean": float(np.mean([r.get("box_iou", 0) for r in rows])),
        "box_centroid_rate": float(np.mean([r.get("box_has_centroid", False) for r in rows])),
        "poly_parse_rate": np.mean([r["poly_parsed"] for r in rows]),
        "poly_valid_rate": float(np.mean([r.get("poly_problems") == "ok" for r in rows])),
        "poly_iou_mean": float(np.mean([r.get("poly_iou", 0) for r in rows])),
        "poly_centroid_rate": float(np.mean([r.get("poly_has_centroid", False) for r in rows])),
        "centre_box_baseline_iou": round(float(base), 3),
        "sec_per_box": float(np.mean([r["box_s"] for r in rows])),
        "sec_per_poly": float(np.mean([r["poly_s"] for r in rows])),
        "peak_rss_gb": round(peak_rss_gb(), 2), "mps_driver_gb": round(mps_mem_gb(), 2),
    }
    agg = {k: (round(float(v), 3) if isinstance(v, (float, np.floating)) else v) for k, v in agg.items()}
    print("SUMMARY", json.dumps(agg), flush=True)
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    (out / f"t2_{name}_{poly_prompt}.json").write_text(json.dumps({"summary": agg, "rows": rows}, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:])
