"""Hash gate of the polygonizer seam (B1.4a; D3 step 3 re-runs it): rebuild every recorded target set through the
normal builders and compare each target_set_sha256 byte for byte.

    .venv/bin/python -I -m experiments.phase2.verify_target_hashes [manifest.json ...] [--store data/store] [--procs 6]

Without manifests it checks the recorded sets: data/prepared, data/targets, data/eval_views/lidc-val-{shifted,
rescaled,rotated} and data/views/lidc-train-{shifted,rescaled}. Two layouts:
- prepared (objdel.train.dataset.prepare, build_eval_views): sets["<ds>/<split>"].target_set_sha256 over the records
  of <dir>/<ds>/<split>.jsonl; each record's store item is framed as it was built (the policy frame, or the record's
  "view": offset, side and op; a rotated window through build_eval_views._rotated);
- targets (experiments.phase1.build_targets): files["<ds>/<split>.jsonl"].sha256 over every store item of the split,
  and target_set_sha256 over all four.
Per rater: objdel.targets.build_target with the reward's size rule (min_areas of the frame masks), as both builders
call it; a rater whose build raises ValueError has no target, as there. A store whose manifest SHA-256 differs from
the recorded one is a mismatch (its targets differ for another reason). Prints one line per set; exit 1 on any
mismatch. Writes nothing. Reads untrusted data: run with python -I.
"""
import argparse
import hashlib
import json
from itertools import batched
from multiprocessing import Pool
from pathlib import Path

from experiments.phase2.build_eval_views import _rotated
from objdel.data import store
from objdel.geometry import SCHEMAS
from objdel.targets import build_target, frame, min_areas, target_set_hash

MANIFESTS = ("data/prepared/manifest.json", "data/targets/manifest.json",
             *(f"data/eval_views/lidc-val-{k}/manifest.json" for k in ("shifted", "rescaled", "rotated")),
             *(f"data/views/lidc-train-{k}/manifest.json" for k in ("shifted", "rescaled")))


def _triples(job) -> list[tuple[str, str, str]]:
    it, view = job
    if view is None:
        f = frame(it)
    elif view["op"] != "id":
        f = _rotated(it, tuple(view["offset"]), view["side"], view["op"])
    else:
        f = frame(it, offset=tuple(view["offset"]), side=view["side"])
    names, ma, out = SCHEMAS[it.dataset], min_areas(f.masks), []
    for r, rater in enumerate(it.raters):
        try:
            out.append((it.item_id, rater, build_target({n: f.masks[n][r] for n in names}, names, None, ma).text))
        except ValueError:
            pass
    return out


def _rebuild(pool, store_root, ds, split, views=None) -> list:
    """Target triples of one split: the items of `views` ({item_id: view or None}), or every store item. Items go
    to the pool in batches, which bounds memory."""
    jobs = ((it, None if views is None else views[it.item_id]) for it in store.load(ds, split, store_root)
            if views is None or it.item_id in views)
    out, seen = [], set()
    for batch in batched(jobs, 12):  # prepare._BATCH: RIGA items hold 6 x 1584 x 2376 masks
        seen.update(it.item_id for it, _ in batch)
        out += [t for ts in pool.map(_triples, batch) for t in ts]
    if views is not None and seen != set(views):
        raise ValueError(f"{ds}/{split}: {len(set(views) - seen)} recorded items are not in the store")
    return out


def check(path, store_root, pool) -> list[tuple[str, str, str]]:
    """(name, recorded, rebuilt) per target set of one manifest."""
    path = Path(path)
    m = json.loads(path.read_text())
    have = hashlib.sha256(Path(store_root, "manifest.json").read_bytes()).hexdigest()
    if m["store_manifest_sha256"] != have:
        return [(f"{path} store", m["store_manifest_sha256"], have)]
    if "sets" in m:
        rows = []
        for key, s in m["sets"].items():
            with open(path.parent / f"{key}.jsonl") as f:
                views = {r["item_id"]: r.get("view") for r in map(json.loads, f)}
            rows.append((f"{path} {key}", s["target_set_sha256"],
                         target_set_hash(_rebuild(pool, store_root, *key.split("/"), views))))
        return rows
    rows, every = [], []
    for rel, e in m["files"].items():
        t = _rebuild(pool, store_root, *rel.removesuffix(".jsonl").split("/"))
        every += t
        rows.append((f"{path} {rel}", e["sha256"], target_set_hash(t)))
    return rows + [(f"{path} all", m["target_set_sha256"], target_set_hash(every))]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("manifests", nargs="*", default=MANIFESTS)
    ap.add_argument("--store", default="data/store")
    ap.add_argument("--procs", type=int, default=6)
    a = ap.parse_args(argv)
    bad = 0
    with Pool(a.procs) as pool:
        for p in a.manifests:
            for name, recorded, rebuilt in check(p, a.store, pool):
                bad += recorded != rebuilt
                print(f"{'ok' if recorded == rebuilt else 'MISMATCH':8} {recorded[:16]} {rebuilt[:16]} {name}",
                      flush=True)
    print(f"{bad} mismatches" if bad else "every target set is byte-identical")
    return int(bad > 0)


if __name__ == "__main__":
    raise SystemExit(main())
