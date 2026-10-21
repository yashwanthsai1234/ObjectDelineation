"""Build the fixed-seed shifted, rescaled and rotated val views (plan §5.0 "Views", §6.1; C-VIEW eval layout).

    .venv/bin/python -I -m experiments.phase2.build_eval_views [--store data/store] [--out data/eval_views]
        [--datasets lidc] [--split val] [--seed 0] [--procs 2] [--limit N]

Per dataset the directories <out>/<dataset>-<split>-shifted, -rescaled and (LIDC) -rotated, each in the data/prepared layout
(images/, caches/, <dataset>/<split>.jsonl, manifest.json with each jsonl's SHA-256 and the target-set hash), so
objdel.train.dataset.records / cache_loader / grpo_rows, objdel.train.serve.preload and eval_run read them
unchanged. Each record is objdel.train.dataset._prepare's (frame image as PNG, one certified target per rater,
ExpertCache.build of the view's rater stack) plus "view": {"offset": [dx, dy], "side": S, "kept_side": bool}.
- shifted: the policy side S, one uniform offset per item (augment.draw with side=S: union inside with 2 px
  clearance, no other outlined nodule);
- rescaled: S moved to a neighbouring divisor (augment.neighbour_side; LIDC 50 <-> 100 <-> 125, RIGA 200 <-> 250 <->
  500) plus a shift; items with no admissible neighbour keep S (kept_side, counted in the manifest; plan §5.0: if
  > 30 %, the rescaled check is weaker than stated);
- rotated (LIDC only, D10: fundus orientation is physiological): the shifted view's windows (same draws) with the D4
  op r90 applied at frame resolution to image and masks before the Pillow upscale; targets and caches are rebuilt
  from the rotated masks. select_ckpt --rotated compares it with the shifted view (report-only attribution check,
  plan §6.1 R1).
Draws: augment.row_rng(seed, view, item_id), independent of the item order. An item with no admissible shifted
window keeps the policy frame and is counted (fallback). objdel.train.dataset (and with it `datasets`) is imported
only when building, so --help and the CI [dev] job run without the train extra (BK4-2). Reads untrusted data: run
with python -I.

References: experiments/phase2/PLAN.md §5.0, §5.1; A4 §2.5 proposal 1; objdel.train.dataset.prepare (layout).
"""
import argparse
import hashlib
import json
import platform
from dataclasses import replace
from importlib.metadata import version
from itertools import batched
from multiprocessing import Pool
from pathlib import Path

from objdel.data import store
from objdel.targets import FRAMES, augment, frame, target_set_hash
from objdel.targets.frame import rgb

KINDS = ("shifted", "rescaled", "rotated")
ROTATION = {"lidc": "r90"}  # the rotated view's op per dataset; no rotated RIGA view


def _rotated(it, offset, side, op):
    """frame(it, offset, side) with `op` applied to masks and image at frame resolution, then the upscale."""
    f = frame(it, replace(FRAMES[it.dataset], image_side=None), offset=offset, side=side)
    image = rgb(augment.apply(f.image.transpose(2, 0, 1), op), FRAMES[it.dataset].image_side)
    return f._replace(image=image, masks={n: augment.apply(m, op) for n, m in f.masks.items()})


def _record(job):
    from objdel.train import dataset as D  # lazy: needs the train extra (BK4-2)
    out, it, kind, seed = job
    bare = replace(it, image=None)
    S0 = frame(bare).W
    rng = augment.row_rng(seed, "shifted" if kind == "rotated" else kind, it.item_id)  # rotated: shifted's windows
    got = augment.neighbour_side(bare, S0, rng) if kind == "rescaled" else None
    kept = got is None
    got = got or augment.draw(bare, rng, side=S0)
    (dx, dy), S = got or ((0, 0), S0)
    f = _rotated(it, (dx, dy), S, ROTATION[it.dataset]) if kind == "rotated" else frame(it, offset=(dx, dy), side=S)
    rec = D._prepare((out, it.item_id, it.dataset, it.group_id, it.raters, f))  # S1-2: needs a public helper (FX-T)
    rec["view"] = {"offset": [dx, dy], "side": S, "kept_side": kept and kind == "rescaled", "fallback": got is None,
                   "op": ROTATION[it.dataset] if kind == "rotated" else "id"}
    return rec


def build(store_root, out_root, datasets=("lidc",), split="val", seed=0, procs=2, limit=None) -> dict:
    from objdel.reward import fingerprint
    from objdel.train import dataset as D  # lazy: needs the train extra (BK4-2)
    manifests = {}
    with Pool(procs) as pool:
        for ds in datasets:
            for kind in (k for k in KINDS if k != "rotated" or ds in ROTATION):
                name = f"{ds}-{split}-{kind}"
                out = Path(out_root, name)
                (out / ds).mkdir(parents=True, exist_ok=True)
                its = [it for i, it in enumerate(store.load(ds, split, store_root))
                       if (limit is None or i < limit) and it.image is not None]
                recs = []
                for chunk in batched(((str(out), it, kind, seed) for it in its), 12):  # bounds parent memory
                    recs += pool.map(_record, chunk)
                text = "".join(json.dumps(r, sort_keys=True) + "\n" for r in recs)
                rel = f"{ds}/{split}.jsonl"
                (out / rel).write_text(text)
                triples = [(r["item_id"], k, t) for r in recs for k, t in zip(r["raters"], r["targets"], strict=True)
                           if t]
                m = {"view": name, "seed": seed,
                     "store_manifest_sha256": hashlib.sha256(Path(store_root, "manifest.json").read_bytes()).hexdigest(),
                     "reward_fingerprint": fingerprint(), "prompts": D.PROMPTS,
                     "versions": {"python": platform.python_version(), **{p: version(p) for p in D._PINS}},
                     "files": {rel: {"sha256": hashlib.sha256(text.encode()).hexdigest()}},
                     "sets": {f"{ds}/{split}": {
                         "items": len(recs), "skipped_no_image": 0, "targets": len(triples),
                         "target_set_sha256": target_set_hash(triples),
                         "target_errors": [e for r in recs for e in r["errors"]],
                         "rl_items": sum(r["rl"] for r in recs),
                         "rl_targets": sum(t is not None for r in recs if r["rl"] for t in r["targets"]),
                         "disagreement_items": sum(r["disagreement"] for r in recs),
                         "kept_side": sum(r["view"]["kept_side"] for r in recs),
                         "fallback": sum(r["view"]["fallback"] for r in recs),
                         "side": {str(s): sum(r["view"]["side"] == s for r in recs)
                                  for s in sorted({r["view"]["side"] for r in recs})}}}}
                (out / "manifest.json").write_text(json.dumps(m, indent=1, sort_keys=True) + "\n")
                manifests[name] = m
    return manifests


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--store", default="data/store")
    ap.add_argument("--out", default="data/eval_views")
    ap.add_argument("--datasets", default="lidc")
    ap.add_argument("--split", default="val")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--limit", type=int, default=None, help="first n items per dataset (smoke test)")
    a = ap.parse_args(argv)
    ms = build(a.store, a.out, a.datasets.split(","), a.split, a.seed, a.procs, a.limit)
    print(json.dumps({k: m["sets"] for k, m in ms.items()}, indent=1))


if __name__ == "__main__":
    main()
