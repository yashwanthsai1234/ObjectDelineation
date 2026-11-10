"""Build the SFT targets of every train/val item and rater from the data store, and measure K5-K7 and R1 (target part).

    .venv/bin/python -I -m experiments.phase1.build_targets [--store data/store] [--out data/targets] [--procs 3]

Per item: objdel.targets.frame (crop policy), then one target per rater with the reward's size filter
(min_area of the raters' union per structure, and of the rim union for RIGA). Every target is re-checked from its
text alone (parse -> polycore -> signature, IoU), and its answer length is counted with the tokenizers of the three
screened models. Writes <out>/<dataset>/<split>.jsonl (item_id, rater, W, H, box, scale, text), <out>/manifest.json
(target-set SHA-256 per file and overall, pinned versions, config) and <out>/summary.json.
Reads untrusted data: run with python -I.

References:
- Plan §4.4 (SFT targets, D11) and §7 (K5-K7, R1): experiments/phase1/PLAN.md; implementation note
  experiments/phase1/impl/M5_targets.md. The target algorithms and their sources: objdel.targets (build, simplify).
- Tokenizers (answer length without special tokens): transformers AutoTokenizer.from_pretrained and
  PreTrainedTokenizer.__call__(add_special_tokens=False):
  https://huggingface.co/docs/transformers/main/en/model_doc/auto ,
  https://huggingface.co/docs/transformers/main/en/main_classes/tokenizer
- multiprocessing.Pool.imap (ordered results, chunksize): https://docs.python.org/3/library/multiprocessing.html
"""
import argparse
import hashlib
import json
import platform
import time
from collections import Counter
from importlib.metadata import version
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from objdel.data import store
from objdel.geom.shape import to_mask
from objdel.geom.text import SCHEMAS, parse, serialize
from objdel.geom.topology import signature
from objdel.targets import FRAMES, BuildConfig, build_target, frame, min_areas, target_set_hash
from objdel.train.configs import MODELS  # the three screened models; their tokenizers measure answer lengths

PINS = ("numpy", "numba", "opencv-python-headless", "objdel", "transformers")
_TOK = {}


def _tokens(text):
    if not _TOK:
        from transformers import AutoTokenizer
        _TOK.update({k: AutoTokenizer.from_pretrained(v) for k, v in MODELS.items()})
    return {k: len(t(text, add_special_tokens=False)["input_ids"]) for k, t in _TOK.items()}


def _job(job):
    """All raters of one framed item -> rows (one per rater) with the text and its independent re-check."""
    item_id, ds, raters, masks = job
    names = SCHEMAS[ds]
    W = next(iter(masks.values())).shape[2]
    derived = {"rim": masks["disc"] & ~masks["cup"]} if ds == "riga" else {}
    ma = min_areas(masks)
    rows = []
    for r, rater in enumerate(raters):
        src = {n: masks[n][r] for n in names}
        row = {"item_id": item_id, "rater": rater, "W": W}
        try:
            t = build_target(src, names, BuildConfig(), ma)
        except ValueError as e:
            rows.append({**row, "error": str(e)})
            continue
        p = parse(t.text, names)
        row.update(text=t.text, parse=p.ok, time=t.info["time"], tokens=_tokens(t.text), obj={})
        out = {}
        for n in names:
            g = to_mask(p.objs[n], W, W)
            out[n] = g.mask
            u = np.count_nonzero(g.mask | src[n])
            row["obj"][n] = {
                "legal": g.legal, "iou": 1.0 if u == 0 else np.count_nonzero(g.mask & src[n]) / u,
                "k5": signature(g.mask, ma[n]) == signature(src[n], ma[n]),
                "raw": signature(g.mask) == signature(src[n]), "tau": t.info[n]["tau"],
                "vertices": t.info[n]["vertices"], "tokens": _tokens(serialize({n: p.objs[n]}, (n,)))}
        if derived:
            rim, rim_src = out["disc"] & ~out["cup"], derived["rim"][r]
            row["obj"]["rim"] = {"k5": signature(rim, ma["rim"]) == signature(rim_src, ma["rim"]),
                                 "raw": signature(rim) == signature(rim_src)}
        rows.append(row)
    return rows


def _jobs(root, ds, split, frames):
    for it in store.load(ds, split, root):
        f = frame(it)
        frames.append({"item_id": it.item_id, "W": f.W, "H": f.H, "box": [int(v) for v in f.box], "scale": f.scale})
        yield it.item_id, ds, it.raters, f.masks


def _pct(a, q):
    return float(np.percentile(a, q)) if len(a) else None


def summarize(rows, frames):
    ok = [r for r in rows if "error" not in r]
    s = {"items": len(frames), "targets": len(rows), "errors": [(r["item_id"], r["rater"], r["error"])
                                                                for r in rows if "error" in r],
         "parse": float(np.mean([r["parse"] for r in ok])),
         "frame_side": dict(sorted(Counter(f["W"] for f in frames).items())),
         "scale": dict(sorted(Counter(f["scale"] for f in frames).items())),
         "build_s": {"median": _pct([r["time"] for r in ok], 50), "p95": _pct([r["time"] for r in ok], 95),
                     "max": max(r["time"] for r in ok), "total": sum(r["time"] for r in ok)},
         "answer_tokens": {k: {q: _pct([r["tokens"][k] for r in ok], q) for q in (50, 95, 99, 100)}
                           for k in MODELS}, "objects": {}}
    for n in sorted({n for r in ok for n in r["obj"]}):
        o = [r["obj"][n] for r in ok]
        e = {"n": len(o), "K5": float(np.mean([x["k5"] for x in o])),
            "raw_topology": float(np.mean([x["raw"] for x in o]))}
        if "iou" in o[0]:
            iou = np.array([x["iou"] for x in o])
            e.update(legal=float(np.mean([x["legal"] for x in o])), iou_min=float(iou.min()),
                     iou_median=float(np.median(iou)),
                     vertices={q: _pct([x["vertices"] for x in o], q) for q in (50, 95, 100)},
                     tau=dict(sorted(Counter(str(x["tau"]) for x in o).items())),
                     tokens={k: {q: _pct([x["tokens"][k] for x in o], q) for q in (50, 95, 100)} for k in MODELS})
        s["objects"][n] = e
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default="data/store")
    ap.add_argument("--out", default="data/targets")
    ap.add_argument("--procs", type=int, default=3)
    a = ap.parse_args()
    out = Path(a.out)
    manifest = {"config": {"build": BuildConfig().__dict__, "frames": {k: v.__dict__ for k, v in FRAMES.items()}},
                "versions": {"python": platform.python_version(), **{p: version(p) for p in PINS}},
                "store_manifest_sha256": hashlib.sha256((Path(a.store) / "manifest.json").read_bytes()).hexdigest(),
                "files": {}}
    summary, all_rows = {}, []
    t0 = time.time()
    with Pool(a.procs) as pool:
        for ds in ("riga", "lidc"):
            for split in ("train", "val"):
                frames, rows = [], []
                for rs in pool.imap(_job, _jobs(a.store, ds, split, frames), chunksize=4):
                    rows += rs
                rel = f"{ds}/{split}.jsonl"
                (out / ds).mkdir(parents=True, exist_ok=True)
                fr = {f["item_id"]: f for f in frames}
                with open(out / rel, "w") as fh:
                    for r in rows:
                        if "text" in r:
                            f = fr[r["item_id"]]
                            fh.write(json.dumps({"item_id": r["item_id"], "rater": r["rater"], "W": f["W"],
                                                 "H": f["H"], "box": f["box"], "scale": f["scale"],
                                                 "text": r["text"]}) + "\n")
                triples = [(r["item_id"], r["rater"], r["text"]) for r in rows if "text" in r]
                manifest["files"][rel] = {"targets": len(triples), "sha256": target_set_hash(triples)}
                all_rows += triples
                summary[f"{ds}/{split}"] = summarize(rows, frames)
                print(rel, json.dumps(summary[f"{ds}/{split}"]["objects"]), f"{time.time() - t0:.0f}s", flush=True)
    manifest["target_set_sha256"] = target_set_hash(all_rows)
    summary["wall_s"] = round(time.time() - t0, 1)
    summary["target_set_sha256"] = manifest["target_set_sha256"]
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
