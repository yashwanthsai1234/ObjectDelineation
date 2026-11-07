"""Frozen checker eval set (LIDC only): real masks and polygons with the CPU checker's answers as references.

Run: .venv/bin/python -I -m experiments.cuda_checker.prepare_data  (writes data/checker_eval/)

- masks.npz: every LIDC Phase 0 mask (data/fixtures/lidc_full.npz, lidc_experts.npz), packed as `bits` (flat uint8),
  `hw` (n, 2) and `start` (n + 1); references from objdel.checker.certify (b0, b1, area, perimeter, corners, chi8)
  and `ring_sha`, the SHA-256 of objdel.checker.polygonize's rings (vertices int64 + signed area, in output order).
- polygons.npz: LIDC polygons on their crop (W, H): training targets, eval-view targets, probe targets and the R1
  model's outputs (illegal ones included: they are the hard cases). Packed like polycore.pack: `xy` (P, 2) int32,
  `ring_start` (R + 1), `shape_start` (n + 1, into rings), `roles` (R,) as geometry.shape.to_mask builds them;
  references from polycore (G = 1000): rasterize areas, px_multi, rings_mixed_sign, validate totals (`defects`,
  DEFECTS order) and `evenodd_sha` / `nonzero_sha`, the SHA-256 of np.packbits of each (H, W) mask.
- manifest.json: counts per source and the SHA-256 of every input and output file.

RIGA is never read (its licence forbids rented hosts; this set goes to a GPU box).
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from objdel.checker import certify, polycore, polygonize
from objdel.geometry import parse

OUT = Path("data/checker_eval")
MASK_FILES = ["data/fixtures/lidc_full.npz", "data/fixtures/lidc_experts.npz"]
BUNDLE = Path("data/views/lidc-bundle/rows.jsonl")
VIEWS = Path("data/eval_views")
R1 = Path("runs/v1/r1-lfm/select/e1")
R1_VIEW = {"shifted": "lidc-val-shifted", "shifted-full": "lidc-val-shifted", "shifted-k16": "lidc-val-shifted",
           "rescaled": "lidc-val-rescaled", "rotated": "lidc-val-rotated", "probe-precision": "probe-precision",
           "probe-topology": "probe-topology"}
SRC = {"train_target": 0, "eval_target": 1, "probe_target": 2, "model_output": 3}


def sha(b: bytes) -> bytes:
    return hashlib.sha256(b).hexdigest().encode()


def file_sha(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def jsonl(p):
    with open(p) as f:
        return [json.loads(line) for line in f if line.strip()]


def masks():
    bits, hw, src, ref = [], [], [], {k: [] for k in ("b0", "b1", "area", "perimeter", "corners", "chi8")}
    ring_sha = []
    for i, f in enumerate(MASK_FILES):
        z = np.load(f, allow_pickle=False)
        for key in sorted(z.files):
            for m in z[key].astype(bool):
                r = certify(m)
                assert r.certified, (f, key, r.failures)  # the CPU checker agrees with itself on every reference
                for k in ref:
                    ref[k].append(getattr(r, k))
                ring_sha.append(sha(b"".join(np.ascontiguousarray(g.pts, np.int64).tobytes()
                                             + int(g.area).to_bytes(8, "little", signed=True) for g in polygonize(m))))
                bits.append(m.astype(np.uint8).ravel())
                hw.append(m.shape)
                src.append(i)
    start = np.concatenate([[0], np.cumsum([b.size for b in bits])]).astype(np.int64)
    return dict(bits=np.concatenate(bits), hw=np.array(hw, np.int32), start=start, src=np.array(src, np.uint8),
                ring_sha=np.array(ring_sha, "S64"), **{k: np.array(v, np.int64) for k, v in ref.items()})


def polygon_cases():
    """(text, W, H, src) of every distinct LIDC polygon text on its crop size."""
    seen = {}
    for r in jsonl(BUNDLE):
        seen.setdefault((r["target"], r["side"], r["side"]), SRC["train_target"])
    meta = {}
    for v in sorted(p.name for p in VIEWS.iterdir()):
        kind = SRC["probe_target"] if v.startswith("probe") else SRC["eval_target"]
        for r in jsonl(VIEWS / v / "lidc/val.jsonl"):
            meta[(v, r["item_id"])] = (r["W"], r["H"])
            for t in r["targets"]:
                seen.setdefault((t, r["W"], r["H"]), kind)
    for sub, v in R1_VIEW.items():
        for c in jsonl(R1 / sub / "completions.jsonl"):
            seen.setdefault((c["text"], *meta[(v, c["id"])]), SRC["model_output"])
    return seen


def polygons():
    xy, ring_start, shape_start, roles, W, H, src, unparsed = [], [0], [0], [], [], [], [], 0
    ref = {k: [] for k in ("area_evenodd", "area_nonzero", "px_multi", "rings_mixed_sign", "defects",
                           "evenodd_sha", "nonzero_sha")}
    for (text, w, h), s in polygon_cases().items():
        p = parse(text, ["nodule"])
        if not p.ok:
            unparsed += 1
            continue
        rings, rl = [], []
        for outer, holes in p.objs["nodule"]:
            rl += [-1] + [len(rings)] * len(holes)
            rings += [outer, *holes]
        if not rings:
            continue  # "no object": nothing for the checker to do
        shape = polycore.pack(rings)
        ras = polycore.rasterize(shape, w, h, G=1000)
        for k in ("area_evenodd", "area_nonzero", "px_multi", "rings_mixed_sign"):
            ref[k].append(int(ras[k]))
        ref["defects"].append(polycore.validate(shape, roles=np.array(rl, np.int32), G=1000)["total"])
        ref["evenodd_sha"].append(sha(np.packbits(np.asarray(ras["evenodd"], bool)).tobytes()))
        ref["nonzero_sha"].append(sha(np.packbits(np.asarray(ras["nonzero"], bool)).tobytes()))
        xy += rings
        ring_start += list(ring_start[-1] + np.cumsum([len(r) for r in rings]))
        shape_start.append(shape_start[-1] + len(rings))
        roles += rl
        W.append(w), H.append(h), src.append(s)
    out = dict(xy=np.concatenate(xy).astype(np.int32), ring_start=np.array(ring_start, np.int64),
               shape_start=np.array(shape_start, np.int64), roles=np.array(roles, np.int32),
               W=np.array(W, np.int32), H=np.array(H, np.int32), src=np.array(src, np.uint8),
               defects=np.array(ref.pop("defects"), np.int64),
               evenodd_sha=np.array(ref.pop("evenodd_sha"), "S64"), nonzero_sha=np.array(ref.pop("nonzero_sha"), "S64"),
               **{k: np.array(v, np.int64) for k, v in ref.items()})
    return out, unparsed


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    m = masks()
    np.savez_compressed(OUT / "masks.npz", **m)
    p, unparsed = polygons()
    np.savez_compressed(OUT / "polygons.npz", **p)
    inputs = [*MASK_FILES, BUNDLE, *sorted(VIEWS.glob("*/lidc/val.jsonl")), *sorted(R1.glob("*/completions.jsonl"))]
    illegal = (p["defects"].sum(1) > 0)
    manifest = {
        "licence": "LIDC-IDRI (CC BY 3.0) only; no RIGA",
        "masks": {"n": int(m["hw"].shape[0]), "per_file": {f: int((m["src"] == i).sum()) for i, f in enumerate(MASK_FILES)},
                  "holes": int((m["b1"] > 0).sum()), "multi_part": int((m["b0"] > 1).sum())},
        "polygons": {"n": int(p["W"].size), "unparsed_skipped": unparsed,
                     "per_src": {k: int((p["src"] == v).sum()) for k, v in SRC.items()},
                     "illegal": int(illegal.sum()), "illegal_model_outputs": int((illegal & (p["src"] == 3)).sum()),
                     "with_holes": int(sum((p["roles"][p["shape_start"][i]:p["shape_start"][i + 1]] >= 0).any()
                                           for i in range(p["W"].size)))},
        "files": {f.name: file_sha(f) for f in (OUT / "masks.npz", OUT / "polygons.npz")},
        "inputs": {str(f): file_sha(f) for f in inputs},
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(json.dumps({k: manifest[k] for k in ("masks", "polygons", "files")}, indent=1))


if __name__ == "__main__":
    main()
