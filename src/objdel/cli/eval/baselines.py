"""D17 baselines on the eval views: image-blind prior, MedSAM (oracle box, plus a jittered box) and SAM 2.1-tiny.

    .venv/bin/python -I -m objdel.cli.eval.baselines --out data/eval/baselines.json \
        [--views centred=data/prepared shifted=data/eval_views/lidc-val-shifted
        rescaled=data/eval_views/lidc-val-rescaled] \
        [--datasets lidc riga] [--models prior medsam sam2] [--every 1] [--threads 2]

Every row is labelled with its view and scored by objdel.verify.reward.score.score_mask (IoU_tc) on the view's prepared
caches, the gate path (§5.1 G-prior, G-geometry). Promoted from archive/phase2_lab/A4 null_prior.py,
null_prior_jitter.py (the prior on shifted frames is now simply the prior on the shifted view) and medsam_check.py
(official min-max variant only).
- prior: a disc at the frame centre whose diameter is the median train-target diameter for that frame side S and
  structure (fit on the centred train records' certified targets, >= 5 px); never looks at the image.
- medsam: flaviagiammarino/medsam-vit-base, box = the view's rater-union box + 1 px (LIDC) / + 5 px (RIGA) mapped to
  the model image; "medsam_jitter" moves each box side by U(0, 20) px at 1024 scale in a random direction.
- sam2: facebook/sam2.1-hiera-tiny (Apache-2.0), the same oracle box.
For G-prior (a) the report gives, per dataset and structure, the per-item centred win share over the prior (K = 1)
of MedSAM and SAM 2.1-tiny, with a two-sided sign test. Torch / transformers are imported only by the SAM models
(CPU; ~3.5 h for all frames, A4 §3.3). Reads untrusted data: run with -I.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

from objdel.core import datasets
from objdel.core.io import read_jsonl
from objdel.geom.shape import to_mask
from objdel.geom.text import parse
from objdel.verify.reward.cache import ExpertCache
from objdel.verify.reward.score import score_mask

MEDSAM, SAM2 = "flaviagiammarino/medsam-vit-base", "facebook/sam2.1-hiera-tiny"
MARGIN = {"lidc": 1, "riga": 5}


def _records(root, dataset, split):
    return read_jsonl(Path(root, dataset, f"{split}.jsonl"))


def _caches(root, rec):
    return {n: ExpertCache.load(Path(root, "caches", rec["item_id"], f"{n}.npz"), d) for n, d in rec["caches"].items()}


def fit_prior(root, dataset, split="train"):
    """{S: {structure: median equivalent diameter}} of the certified targets (>= 5 px) per frame side."""
    by = {}
    for r in _records(root, dataset, split):
        for t in r["targets"]:
            p = parse(t, datasets.get(dataset).names) if t else None
            for n, parts in p.objs.items() if p and p.ok else ():
                m = to_mask(parts, r["W"], r["H"]).mask
                if m.sum() >= 5:
                    by.setdefault(r["W"], {}).setdefault(n, []).append(2 * np.sqrt(m.sum() / np.pi))
    return {S: {n: float(np.median(v)) for n, v in d.items()} for S, d in by.items()}


def _disc(S, d):
    yy, xx = np.mgrid[:S, :S] + 0.5
    return (xx - S / 2) ** 2 + (yy - S / 2) ** 2 <= (d / 2) ** 2


def _score_rows(root, dataset, split, predict, every=1):
    """[{"item_id", "iou_tc": {structure: v}}] of predict(rec, caches) -> {structure: mask} on RL items."""
    rows = []
    for i, rec in enumerate(r for r in _records(root, dataset, split) if r["rl"]):
        if i % every:
            continue
        caches = _caches(root, rec)
        masks = predict(rec, caches)
        rows.append(
            {"item_id": rec["item_id"], "iou_tc": {n: float(score_mask(masks[n], c).iou_tc) for n, c in caches.items()}}
        )
    return rows


def score_prior(root, dataset, split, prior, every=1):
    def predict(rec, caches):
        S = rec["W"]
        pS = prior.get(S) or prior[min(prior, key=lambda s: abs(s - S))]
        return {n: _disc(S, pS[n]) for n in caches}

    return _score_rows(root, dataset, split, predict, every)


def win_share(a_rows, b_rows, name):
    """Share of the common items where a's IoU_tc beats b's (ties are not wins); sign test on the non-ties."""
    b = {r["item_id"]: r["iou_tc"][name] for r in b_rows}
    d = np.array([r["iou_tc"][name] - b[r["item_id"]] for r in a_rows if r["item_id"] in b])
    wins, losses = int((d > 0).sum()), int((d < 0).sum())
    p = binomtest(wins, wins + losses).pvalue if wins + losses else 1.0
    return {"n": len(d), "share": float(wins / len(d)) if len(d) else None, "sign_p": float(p)}


def summary(rows):
    names = rows[0]["iou_tc"] if rows else {}
    return {
        n: {
            "items": len(rows),
            "iou_tc_median": float(np.median([r["iou_tc"][n] for r in rows])),
            "iou_tc_ge_0.5": float(np.mean([r["iou_tc"][n] >= 0.5 for r in rows])),
        }
        for n in names
    }


class _Sam:
    """Box-prompted SAM-family model on the CPU (lazy torch / transformers)."""

    def __init__(self, kind, threads=2):
        import torch

        torch.set_num_threads(threads)
        self.torch, self.kind = torch, kind
        if kind == "sam2":
            from transformers import Sam2Model, Sam2Processor

            self.model, self.proc = Sam2Model.from_pretrained(SAM2).eval(), Sam2Processor.from_pretrained(SAM2)
        else:
            from transformers import SamModel, SamProcessor

            self.model, self.proc = SamModel.from_pretrained(MEDSAM).eval(), SamProcessor.from_pretrained(MEDSAM)

    def __call__(self, img, box):
        """Boolean mask at image resolution for an image-pixel box (x0, y0, x1, y1)."""
        with self.torch.no_grad():
            inp = self.proc(images=img, input_boxes=[[list(map(float, box))]], return_tensors="pt")
            if self.kind == "medsam":  # official MedSAM: min-max of the resized content
                h, w = inp["reshaped_input_sizes"][0].tolist()
                c = inp["pixel_values"][..., :h, :w]
                inp["pixel_values"][..., :h, :w] = (c - c.min()) / max(float(c.max() - c.min()), 1e-8)
            out = self.model(**inp, multimask_output=False)
            pp = self.proc.post_process_masks if self.kind == "sam2" else self.proc.image_processor.post_process_masks
            args = (
                (inp["original_sizes"],)
                if self.kind == "sam2"
                else (inp["original_sizes"], inp["reshaped_input_sizes"])
            )
            return np.asarray(pp(out.pred_masks.cpu(), *args)[0][0, 0], bool)


def _to_frame(mask_img, S):
    idx = ((np.arange(S) + 0.5) * mask_img.shape[0] / S).astype(int)  # sample at frame pixel centres
    return mask_img[np.ix_(idx, idx)]


def sam_predict(model, root, dataset, jitter_rng=None):
    from objdel.policy.rows import image_path, rgb

    def predict(rec, caches):
        img = np.asarray(rgb(image_path(root, rec["item_id"])))
        S, k = rec["W"], img.shape[0] / rec["W"]
        out = {}
        for n, c in caches.items():
            y0, y1, x0, x1 = c.box
            m = MARGIN[dataset]
            box = np.array([x0 - m, y0 - m, x1 + m, y1 + m], float) * k
            if jitter_rng is not None:
                box += jitter_rng.uniform(0, 20, 4) * jitter_rng.choice([-1, 1], 4) * img.shape[0] / 1024
            out[n] = _to_frame(model(img, np.clip(box, 0, img.shape[0])), S)
        return out

    return predict


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--views",
        nargs="+",
        default=[
            "centred=data/prepared",
            "shifted=data/eval_views/lidc-val-shifted",
            "rescaled=data/eval_views/lidc-val-rescaled",
        ],
    )
    ap.add_argument("--prepared", default="data/prepared", help="centred train records the prior is fit on")
    ap.add_argument("--datasets", nargs="+", default=["lidc"])
    ap.add_argument("--models", nargs="+", default=["prior", "medsam", "sam2"])
    ap.add_argument("--split", default="val")
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    t0, res = time.time(), {}
    sams = {k: _Sam(k, a.threads) for k in ("medsam", "sam2") if k in a.models}
    for ds in a.datasets:
        prior = fit_prior(a.prepared, ds)
        res[ds] = {"prior_diam_by_S": {str(k): v for k, v in sorted(prior.items())}, "views": {}}
        for view, root in (v.split("=", 1) for v in a.views):
            if not Path(root, ds, f"{a.split}.jsonl").is_file():
                continue
            rows = {}
            if "prior" in a.models:
                rows["prior"] = score_prior(root, ds, a.split, prior, a.every)
            for k, m in sams.items():
                rows[k] = _score_rows(root, ds, a.split, sam_predict(m, root, ds), a.every)
                if k == "medsam":
                    rows["medsam_jitter"] = _score_rows(
                        root, ds, a.split, sam_predict(m, root, ds, np.random.default_rng(0)), a.every
                    )
            res[ds]["views"][view] = {"summary": {k: summary(v) for k, v in rows.items()}, "rows": rows}
            if view == "centred" and "prior" in rows:
                res[ds]["centred_win_share_vs_prior"] = {
                    k: {n: win_share(rows[k], rows["prior"], n) for n in datasets.get(ds).names} for k in sams if k in rows
                }
            print(ds, view, json.dumps(res[ds]["views"][view]["summary"]), f"{time.time() - t0:.0f} s", flush=True)
    res["wall_s"] = time.time() - t0
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(res) + "\n")


if __name__ == "__main__":
    main()
