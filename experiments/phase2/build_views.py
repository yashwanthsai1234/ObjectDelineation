"""Build a Phase 2 training view (plan §6.0 bundle, D9-D11; C-VIEW): per-epoch row specs with certified targets.

    .venv/bin/python -I -m experiments.phase2.build_views --out data/views/lidc-bundle-v2 [--store data/store]
        [--dataset lidc] [--split train] [--epochs 3] [--seed 0] [--nonsimple-repeat 4] [--procs 2] [--limit N]

Per item with pixels (objdel.data.store): the policy frame gives the reward's size rule (min_areas: shift-, scale-
and D4-invariant while the union stays inside) and the voting raters (objdel.reward.voters: non-empty size-filtered
signature; I27, B2-4). Per epoch, voting rater and repeat (non-simple targets, i.e. any object other than one part
without holes, are repeated `nonsimple_repeat` times, each repeat with its own draw), one row spec drawn from
augment.row_rng(seed, epoch, item_id, rater_index, rep): the window (augment.draw: S uniform over the policy side
and its ladder neighbours, the sides of the shifted and rescaled eval views, uniform shift with the union inside and
2 px clearance, no other outlined nodule), the op
(D4 for LIDC, id / fh for RIGA), and for LIDC the CT window (L -300 +- 60, W 1400 +- 80) and noise (sigma <= 10 HU).
The target is build_target on the op applied to the rater's mask in that window (certified, and pixel-exact for
LIDC, or ValueError); its raster is compared pixel for pixel with the mask for the per-op exact share (RIGA < 1).
The answer-token budget (I5) is applied when rows are read (objdel.train.dataset.view_rows), per tokenizer.

Writes <out>/rows.jsonl (sorted by epoch, item_id, rater_index, rep), <out>/context/<item_id>.npz (augment.
write_context: int16 HU for LIDC, uint8 RGB for RIGA) and <out>/manifest.json (config, source manifests' SHA-256,
counts, per-op pixel-exact counts, view_hash = SHA-256 of rows.jsonl, contexts = {item_id: SHA-256 of its npz},
checked by dataset.view_rows); two builds are byte-identical. Timing and disk go to <out>/../build_reports/<view>.json,
outside the view (the bundle stays deterministic, BK1-7). Default 3 epochs: R1's 2-epoch branch decays into epoch 3
(S1-8). Reads untrusted data: run with python -I.

References:
- plan experiments/phase2/PLAN.md D9, D10, D11, §6.0; promoted from archive/phase2_lab/B3/dihedral.py (P3) and archive/phase2_lab/A4, archive/phase2_lab/B1 shifts
- multiprocessing.Pool.imap (ordered results): https://docs.python.org/3/library/multiprocessing.html
"""
import argparse
import hashlib
import json
import time
from collections import Counter
from dataclasses import replace
from itertools import islice
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from objdel.checker import signature
from objdel.data import store
from objdel.geometry import SCHEMAS, to_mask
from objdel.reward import voters
from objdel.targets import BuildConfig, augment, build_target, frame, min_areas

NO_JITTER = {"window": None, "noise_sigma": 0.0, "noise_seed": 0}


def _item(job):
    it, out, epochs, seed, repeat = job
    if it.image is None:
        return None  # counted as skipped
    ds, names = it.dataset, SCHEMAS[it.dataset]
    policy = frame(replace(it, image=None))
    ma = min_areas(policy.masks)
    sigs = {n: [signature(m, ma[n]) for m in policy.masks[n]] for n in names}
    vote = [r for r in range(len(it.raters)) if all(r in voters(sigs[n]) for n in names)] or range(len(it.raters))
    ctx = Path(out, "context", it.item_id + ".npz")
    augment.write_context(ctx, it, policy.box)
    bare, rows, stats = replace(it, image=None), [], Counter()
    for e in range(epochs):
        for r in vote:
            simple = all(sigs[n][r] == (0,) for n in names)
            for rep in range(1 if simple else repeat):
                rng = augment.row_rng(seed, e, it.item_id, r, rep)
                got = augment.draw(bare, rng)
                if got is None:  # no admissible moved window: the policy frame (always valid)
                    stats["fallback"] += 1
                (dx, dy), S = got or ((0, 0), policy.W)
                ops = augment.DATASET_OPS[ds]
                op = ops[int(rng.integers(len(ops)))]
                jit = augment.jitter(rng) if ds == "lidc" else NO_JITTER
                f = frame(bare, offset=(dx, dy), side=S)
                src = {n: augment.apply(f.masks[n][r], op) for n in names}
                row = {"row_id": f"{e}/{it.item_id}/{r}/{rep}", "epoch": e, "item_id": it.item_id, "dataset": ds,
                       "rater": it.raters[r], "rater_index": r, "op": op, "offset": [dx, dy], "side": S,
                       "box": list(f.box), **jit, "simple": simple, "rep": rep}
                try:
                    t = build_target(src, names, BuildConfig(), ma)
                except ValueError as err:
                    stats[f"error/{op}"] += 1
                    stats["errors"] += 1
                    print(f"{row['row_id']}: {err}", flush=True)
                    continue
                row["target"] = t.text
                stats[f"rows/{op}"] += 1
                stats[f"exact/{op}"] += all(np.array_equal(to_mask(t.objs[n], S, S).mask, src[n]) for n in names)
                stats[f"side/{S}"] += 1
                rows.append(row)
    return it.item_id, rows, stats, hashlib.sha256(ctx.read_bytes()).hexdigest()


def build(store_root, out, dataset="lidc", split="train", epochs=3, seed=0, nonsimple_repeat=4, procs=2,
          limit=None, prepared=None) -> dict:
    out, t0, cpu0 = Path(out), time.time(), time.process_time()
    out.mkdir(parents=True, exist_ok=True)
    rows, stats, contexts, skipped = [], Counter(), {}, 0
    its = islice(store.load(dataset, split, store_root), limit)
    jobs = ((it, str(out), epochs, seed, nonsimple_repeat) for it in its)
    with Pool(procs) as pool:
        for got in pool.imap(_item, jobs, chunksize=8):
            if got is None:
                skipped += 1
                continue
            item_id, rs, st, contexts[item_id] = got
            rows += rs
            stats.update(st)
    rows.sort(key=lambda r: (r["epoch"], r["item_id"], r["rater_index"], r["rep"]))
    text = "".join(json.dumps(r, sort_keys=True) + "\n" for r in rows)
    (out / "rows.jsonl").write_text(text)
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest() if p and Path(p).is_file() else None  # noqa: E731
    per_op = {op: {k: stats[f"{k}/{op}"] for k in ("rows", "exact", "error")}
              for op in augment.DATASET_OPS[dataset]}
    manifest = {
        "view": out.name, "dataset": dataset, "split": split,
        "config": {
            "epochs": epochs, "seed": seed, "nonsimple_repeat": nonsimple_repeat,
            "ops": augment.DATASET_OPS[dataset],
            "scale": "policy side or a LADDER neighbour, uniform", "clearance": augment.CLEARANCE,
            "window_jitter": augment.WINDOW_JITTER if dataset == "lidc" else None,
            "noise_max": augment.NOISE_MAX if dataset == "lidc" else None, "build": BuildConfig().__dict__},
        "sources": {"store_manifest_sha256": sha(Path(store_root, "manifest.json")),
                    "prepared_manifest_sha256": sha(Path(prepared, "manifest.json")) if prepared else None},
        "counts": {"items": len({r["item_id"] for r in rows}), "skipped_no_image": skipped, "rows": len(rows),
                   "nonsimple_rows": sum(not r["simple"] for r in rows), "errors": stats["errors"],
                   "fallback_windows": stats["fallback"], "per_op": per_op,
                   "side": {k.split("/")[1]: v for k, v in sorted(stats.items()) if k.startswith("side/")}},
        "view_hash": hashlib.sha256(text.encode()).hexdigest(), "contexts": contexts}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    ctx = [Path(out, "context", i + ".npz") for i in contexts]
    report = {"view": out.name, "wall_s": round(time.time() - t0, 1),
              "parent_cpu_s": round(time.process_time() - cpu0, 1), "procs": procs, "context_files": len(ctx),
              "context_mb": round(sum(p.stat().st_size for p in ctx) / 2**20, 1), "rows_mb": round(len(text) / 2**20, 1),
              "exact_share": sum(v["exact"] for v in per_op.values()) / max(len(rows), 1)}
    report_path(out).parent.mkdir(exist_ok=True)
    report_path(out).write_text(json.dumps(report, indent=1) + "\n")
    return manifest


def report_path(out) -> Path:
    """Timing and disk of a build, outside the view dir (not bundled: wall clock is not deterministic, BK1-7)."""
    return Path(out).parent / "build_reports" / f"{Path(out).name}.json"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--store", default="data/store")
    ap.add_argument("--prepared", default="data/prepared")
    ap.add_argument("--out", required=True, help="e.g. data/views/lidc-bundle-v2")
    ap.add_argument("--dataset", default="lidc", choices=sorted(SCHEMAS))
    ap.add_argument("--split", default="train")
    ap.add_argument("--epochs", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--nonsimple-repeat", type=int, default=4)
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--limit", type=int, default=None, help="first n items (smoke test)")
    a = ap.parse_args(argv)
    m = build(a.store, a.out, a.dataset, a.split, a.epochs, a.seed, a.nonsimple_repeat, a.procs, a.limit, a.prepared)
    print(json.dumps(m["counts"], indent=1))
    print(report_path(a.out).read_text())


if __name__ == "__main__":
    main()
