"""A4 sandbox: is Phase 0's MedSAM-on-LIDC IoU (0.33) a pipeline artefact? Box-prompted MedSAM on OUR frames.

Same frame, same split (val), same scoring path as the VLM: frame(item).image (RIGA crop; LIDC 512-px upscale of the
lung-windowed crop), box = union of the raters' masks (+ margin) mapped to the image, mask resampled to the frame at
pixel centres, scored with objdel.reward.score_mask (IoU_tc) on ExpertCache.build(raters), plus IoU vs vote / D_max.
Null baseline with the same information: the ellipse inscribed in the same box.
Variants: 'hf' = the transformers processor as in Phase 0 (rescale 1/255, no min-max); 'minmax' = the official
MedSAM inference (min-max to [0, 1] after resize).
Run: .venv/bin/python -I medsam_check.py --lidc-every 9 --riga-every 3 --out medsam_check.json
"""
import argparse
import json
import time

import numpy as np
import torch
from transformers import SamModel, SamProcessor

from objdel.data import store
from objdel.reward import ExpertCache, score_mask
from objdel.targets import frame

MODEL = "flaviagiammarino/medsam-vit-base"


def iou(a, b):
    u = np.count_nonzero(a | b)
    return 1.0 if u == 0 else np.count_nonzero(a & b) / u


class MedSAM:
    def __init__(self):
        self.dev = "cpu"  # MPS peak footprint 4.5 GB > the 3 GB sandbox budget
        self.model = SamModel.from_pretrained(MODEL).to(self.dev).eval()
        self.proc = SamProcessor.from_pretrained(MODEL)

    @torch.no_grad()
    def embed(self, img, minmax):
        inp = self.proc(img, return_tensors="pt")
        pv = inp["pixel_values"]
        if minmax:  # official MedSAM: min-max of the resized image (padding region excluded by using the content)
            h, w = inp["reshaped_input_sizes"][0].tolist()
            c = pv[..., :h, :w]
            pv = pv.clone()
            pv[..., :h, :w] = (c - c.min()) / max(float(c.max() - c.min()), 1e-8)
        emb = self.model.get_image_embeddings(pv.to(self.dev))
        return emb, inp

    @torch.no_grad()
    def segment(self, emb, inp, box):
        # scale the box like the processor does (longest edge -> 1024)
        H, W = inp["original_sizes"][0].tolist()
        s = 1024 / max(H, W)
        bx = torch.tensor([[[box[0] * s, box[1] * s, box[2] * s, box[3] * s]]], dtype=torch.float32, device=self.dev)
        out = self.model(image_embeddings=emb, input_boxes=bx, multimask_output=False)
        probs = self.proc.image_processor.post_process_masks(out.pred_masks.sigmoid().cpu(), inp["original_sizes"],
                                                             inp["reshaped_input_sizes"], binarize=False)[0]
        return probs[0, 0].numpy() > 0.5


def to_frame(mask_img, S):
    """Image-resolution mask -> S x S frame mask, sampled at frame pixel centres."""
    n = mask_img.shape[0]
    idx = ((np.arange(S) + 0.5) * n / S).astype(int)
    return mask_img[np.ix_(idx, idx)]


def ellipse(box, S):
    yy, xx = np.mgrid[:S, :S] + 0.5
    x0, y0, x1, y1 = box
    return (((xx - (x0 + x1) / 2) / ((x1 - x0) / 2)) ** 2 + ((yy - (y0 + y1) / 2) / ((y1 - y0) / 2)) ** 2) <= 1


def run(ds, every, model, margin_px):
    rows = []
    for i, it in enumerate(store.load(ds, "val", "data/store")):
        if i % every or it.image is None:
            continue
        fr = frame(it)
        S, n = fr.W, fr.image.shape[0]
        k = n / S
        embs = {mm: model.embed(fr.image, mm) for mm in (False, True)}
        for name, st in fr.masks.items():
            cache = ExpertCache.build(st)
            u = st.any(0)
            if not u.any():
                continue
            ys, xs = np.nonzero(u)
            fbox = (xs.min() - margin_px, ys.min() - margin_px, xs.max() + 1 + margin_px, ys.max() + 1 + margin_px)
            ibox = [c * k for c in fbox]
            preds = {f"medsam_{'minmax' if mm else 'hf'}": to_frame(model.segment(*embs[mm], ibox), S)
                     for mm in (False, True)}
            preds["box_ellipse"] = ellipse(fbox, S)
            vote = cache.vote
            for p, m in preds.items():
                t = score_mask(m, cache)
                rows.append({"ds": ds, "obj": name, "pred": p, "iou_tc": float(t.iou_tc), "iou_vote": iou(m, vote),
                             "d_max": max(iou(m, r) for r in st), "diam": float(2 * np.sqrt(u.sum() / np.pi / 1.0)),
                             "S": S})
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--lidc-every", type=int, default=9)
    ap.add_argument("--riga-every", type=int, default=3)
    ap.add_argument("--out", default="medsam_check.json")
    a = ap.parse_args()
    torch.set_num_threads(2)
    m = MedSAM()
    t0 = time.time()
    rows = run("lidc", a.lidc_every, m, 1) + run("riga", a.riga_every, m, 5)
    summ = {}
    for key in sorted({(r["ds"], r["obj"], r["pred"]) for r in rows}):
        R = [r for r in rows if (r["ds"], r["obj"], r["pred"]) == key]
        summ["/".join(key)] = {"n": len(R), **{f"{x}_median": round(float(np.median([r[x] for r in R])), 3)
                                               for x in ("iou_tc", "iou_vote", "d_max")},
                               "iou_tc_ge_0.5": round(float(np.mean([r["iou_tc"] >= 0.5 for r in R])), 3)}
    out = {"summary": summ, "rows": rows, "wall_s": time.time() - t0, "device": m.dev}
    print(json.dumps(summ, indent=1))
    json.dump(out, open(a.out, "w"))
