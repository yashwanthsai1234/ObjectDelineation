"""Checker benchmark (Phase D, D1): wall time, per-stage device time and engine counters of objdel._gpucheck.

Run as a module with -I (it reads data, which is untrusted):
    python -I -m experiments.phase2.perf.bench_checker --out runs/d1/bench.json [--engine current|old] [--stages]
        [--workload NAME ...] [--warmup 2] [--repeats 10]
    python -I -m experiments.phase2.perf.bench_checker --dump-replay OUT.npz [--synthetic N]
    python -I -m experiments.phase2.perf.bench_checker --replay OUT.npz --out replay.json [--engine current|old]
    python -I -m experiments.phase2.perf.bench_checker --rungs runs/d1/rungs.json

Workloads (check_masks): lesion512_b1, lesion512_b256, lesion512_b8192 (512 x 512 frames with two lesion-like blobs,
generated, seed 20261108) and crops (the 37,093 crops of checker_eval/masks.npz). Each runs `--warmup` untimed calls,
then `--repeats` timed ones: median and IQR of the wall time, the counters of one call (gc.stats() after a reset),
and with --stages the median per-stage device time (timing=True).

Inputs resolve under ROOT = $OBJDEL_DATA_ROOT, else `data` (the Mac); on the box the code runs in /workspace/src and
the LIDC bundle sits in /workspace/data/lidc, so the box sets the variable. Every input is checked before anything
loads (exit 1 naming the missing path), and its resolved path and SHA-256 go into the output JSON.

The replay (D-G8): --dump-replay, in the main venv, parses the first 256 rl rows of prepared/lidc/val.jsonl (sorted by
item_id; completion = the row's first target text), loads their caches from prepared/caches/lidc/ and writes the
score_polygons arguments (pack_objects arrays, GRID), obj_cache (object -> cache) and per cache the 13 add_tables
arguments to one .npz, with its SHA-256 beside it (--synthetic N: N two-rater disc items, no data). --replay needs
only numpy and _gpucheck: it checks the SHA-256, uploads the tables and times score_polygons; the per-object results
go into the JSON so two engines can be compared.

--engine old reruns this script by path with the pre-D2 engine's Python (--old-python, default /opt/old/bin/python,
built by run_box.sh --gpu-check from old_engine.tar; Implementation.md D-G7) and this environment, so both engines
read the same root. That venv has no objdel.geom, objdel.verify.reward, cv2 or scipy: the module-level imports below are
the standard library, numpy and objdel._gpucheck only; --dump-replay and --rungs import the rest inside their functions.

--rungs (D1 step 3, input to the D6 grid decision): the vertex count of every certified rung (the target texts of
the first 500 rows of data/targets/{lidc,riga}/train.jsonl under ROOT) and its validity time on this build.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from objdel import _gpucheck as gc

SEED = 20261108
REPLAY_ROWS = 256
LESION = {"lesion512_b1": 1, "lesion512_b256": 256, "lesion512_b8192": 8192}
WORKLOADS = (*LESION, "crops")
MASKS, VAL, CACHES = "checker_eval/masks.npz", "prepared/lidc/val.jsonl", "prepared/caches/lidc"
# the 13 add_tables arguments, in order (ExpertCache.table_args)
TABLE_ARGS = ("raters", "uni", "core", "vote", "dist", "rater_area", "voter", "sig_start", "sig", "n_vote",
              "perimeter", "t5", "min_area")


def data_root() -> Path:
    return Path(os.environ.get("OBJDEL_DATA_ROOT") or "data").resolve()


def file_sha256(path: Path) -> str:
    """SHA-256 hex digest of the file. Not objdel.core.hashing: --engine old runs this module in the pre-D2 venv, whose
    objdel has no hashing module (only the standard library, numpy and objdel._gpucheck; tests/experiments/test_bench_checker
    S10)."""
    with open(path, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def need(root: Path, rels) -> dict:
    """{rel: {path, sha256}} of every input; exit 1 naming the first missing one, before anything loads."""
    for rel in rels:
        if not (root / rel).exists():
            sys.exit(f"bench_checker: missing input {root / rel} (OBJDEL_DATA_ROOT={root})")
    return {rel: {"path": str(root / rel), "sha256": file_sha256(root / rel)} for rel in rels if (root / rel).is_file()}


def device_name() -> str:
    if gc.build() != "cuda":
        return "cpu"
    q = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader", "-i", "0"], capture_output=True,
                       text=True)
    return q.stdout.strip() or "unknown"


def pack(masks) -> tuple:
    hw = np.array([m.shape for m in masks], np.int32).reshape(-1, 2)
    start = np.cumsum([0, *(m.size for m in masks)], dtype=np.int64)
    return np.concatenate([m.ravel() for m in masks]).astype(np.uint8), hw, start


def lesion_frames(n: int, side: int = 512) -> tuple:
    """n frames, each with two filled ellipses (radii 6-60 px), one in three with an elliptic hole; seed SEED."""
    rng = np.random.default_rng(SEED)
    bits = np.zeros((n, side, side), np.uint8)
    for f in range(n):
        for _ in range(2):
            cy, cx = rng.uniform(64, side - 64, 2)
            ry, rx = rng.uniform(6, 60, 2)
            y0, y1, x0, x1 = int(cy - ry), int(cy + ry) + 1, int(cx - rx), int(cx + rx) + 1
            y, x = np.ogrid[y0:y1, x0:x1]
            e = ((y + 0.5 - cy) / ry) ** 2 + ((x + 0.5 - cx) / rx) ** 2
            inside = (e <= 1) & ~((e <= 0.2) & (rng.random() < 1 / 3))
            bits[f, y0:y1, x0:x1] |= inside.astype(np.uint8)
    return bits.reshape(-1), np.full((n, 2), side, np.int32), np.arange(n + 1, dtype=np.int64) * side * side


def timed(call, warmup: int, repeats: int, stages: bool) -> dict:
    """Median and IQR of `call`'s wall time, one call's counters, and with `stages` the median stage times."""
    for _ in range(warmup):
        call(False)
    wall, times = [], []
    for _ in range(repeats):
        gc.reset_stats()
        t0 = time.perf_counter()
        out = call(False)
        wall.append(1e3 * (time.perf_counter() - t0))
        stats = gc.stats()
        if stages:
            times.append(call(True)["timing"])
    q1, med, q3 = np.percentile(wall, [25, 50, 75])
    r = {"repeats": repeats, "median_ms": float(med), "iqr_ms": float(q3 - q1), "wall_ms": wall, "stats": stats}
    if stages:
        r["stages_ms"] = {k: float(np.median([t[k] for t in times])) for k in times[0]}
    return r, out


def bench(a, root: Path) -> dict:
    inputs = need(root, [MASKS] if "crops" in a.workload else [])
    out = {}
    for w in a.workload:
        if w == "crops":
            with np.load(root / MASKS) as z:
                args = (z["bits"], z["hw"], z["start"])
        else:
            args = lesion_frames(LESION[w])
        r, res = timed(lambda t, args=args: gc.check_masks(*args, timing=t), a.warmup, a.repeats, a.stages)
        out[w] = r | {"masks": len(args[1]), "rings": len(res["ring_area"]), "corners": len(res["pts"]),
                      "results_sha256": digest(res)}
    return {"inputs": inputs, "workloads": out}


def digest(res: dict) -> str:
    h = hashlib.sha256()
    for k in sorted(res):
        if k != "timing":
            h.update(k.encode() + np.ascontiguousarray(res[k]).tobytes())
    return h.hexdigest()


# ------------------------------------------------------------------ replay (D-G8)

def synthetic_replay(n: int) -> list:
    """n (completion, cache) pairs: two-rater disc caches on 24-48 px frames and square answers; fixed seed."""
    from objdel.verify.reward.cache import ExpertCache

    rng = np.random.default_rng(SEED)
    out = []
    for _ in range(n):
        side = int(rng.integers(24, 49))
        y, x = np.mgrid[:side, :side]
        r = rng.uniform(4, side / 3)
        c = ExpertCache.build(np.stack([np.hypot(y - side / 2, x - side / 2) <= r + d for d in (0, 1.5)]))
        a, b = sorted(int(v) for v in rng.integers(100, 900, 2))
        text = json.dumps({"nodule": [{"outer": [[a, a], [b + 1, a], [b + 1, b + 1], [a, b + 1]], "holes": []}]},
                          separators=(",", ":"))
        out.append((text, c))
    return out


def dump_replay(path: Path, synthetic: int | None) -> None:
    from objdel.geom.shape import pack_objects
    from objdel.geom.text import GRID, parse
    from objdel.verify.reward.cache import ExpertCache

    root, meta = data_root(), {}
    if synthetic:
        pairs = [((t, ("nodule",), {"nodule": c}), f"synthetic/{i}")
                 for i, (t, c) in enumerate(synthetic_replay(synthetic))]
    else:
        inputs = need(root, [VAL, CACHES])
        with open(root / VAL) as f:
            rows = sorted((r for r in map(json.loads, f) if r["rl"]), key=lambda r: r["item_id"])[:REPLAY_ROWS]
        pairs = []
        for r in rows:
            key = r["item_id"].removeprefix("lidc/")
            caches = {}
            for n, d in r["caches"].items():
                p = root / CACHES / key / f"{n}.npz"
                inputs[f"{CACHES}/{key}/{n}.npz"] = {"path": str(p), "sha256": file_sha256(p)}
                caches[n] = ExpertCache.load(p, d)
            pairs.append(((r["targets"][0], tuple(r["caches"]), caches), r["item_id"]))
        meta = {"root": str(root), "inputs": inputs}
    objs, obj_cache, arrays, caches = [], [], {}, []
    for (text, names, cs), _ in pairs:
        p = parse(text, names)
        if not p.ok:
            continue
        for n in names:
            c = cs[n]
            objs.append((p.objs[n], c.union.shape[1], c.union.shape[0]))
            obj_cache.append(len(caches))
            caches.append(c)
    for i, c in enumerate(caches):
        arrays |= {f"t{i}_{k}": np.asarray(v) for k, v in zip(TABLE_ARGS, c.table_args(), strict=True)}
    names = ("xy", "ring_start", "shape_start", "roles", "W", "H")
    arrays |= dict(zip(names, pack_objects(objs), strict=True))
    meta |= {"item_ids": [i for _, i in pairs], "caches": len(caches)}
    arrays |= {"GRID": np.int64(GRID), "obj_cache": np.array(obj_cache, np.int64), "meta": np.array(json.dumps(meta))}
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        np.savez(f, **arrays)
    with open(path.with_name(path.name + ".sha256"), "w") as f:
        f.write(f"{file_sha256(path)}  {path.name}\n")


def replay(a) -> dict:
    path = Path(a.replay)
    with open(path.with_name(path.name + ".sha256")) as f:
        want = f.read().split()[0]
    if file_sha256(path) != want:
        sys.exit(f"bench_checker: {path}: SHA-256 differs from {path.name}.sha256")
    with np.load(path) as z:
        d = {k: z[k] for k in z.files}
    meta = json.loads(str(d["meta"]))
    handles = []
    try:
        for i in range(meta["caches"]):
            t = [d[f"t{i}_{k}"] for k in TABLE_ARGS]
            handles.append(gc.add_tables(*t[:9], int(t[9]), int(t[10]), float(t[11]), int(t[12])))
        args = (d["xy"], d["ring_start"], d["shape_start"], d["roles"], d["W"], d["H"], int(d["GRID"]),
                np.array(handles, np.int64)[d["obj_cache"]])
        r, res = timed(lambda t: gc.score_polygons(*args, timing=t), a.warmup, a.repeats, a.stages)
    finally:
        for h in handles:
            gc.drop_tables(h)
    res = {k: v.tolist() for k, v in res.items() if k != "timing"}
    inputs = {path.name: {"path": str(path.resolve()), "sha256": want}} | meta.get("inputs", {})
    return {"inputs": inputs, "replay": r | {"objects": len(d["W"]), "items": meta["item_ids"], "results": res}}


# ------------------------------------------------------------------ rungs (D1 step 3, D6)

def rungs(path: Path) -> None:
    from objdel.geom.shape import pack_objects
    from objdel.geom.text import GRID, parse

    root, out = data_root(), []
    rels = [f"targets/{ds}/train.jsonl" for ds in ("lidc", "riga")]
    inputs = need(root, rels)
    for rel in rels:
        with open(root / rel) as f:
            rows = [json.loads(line) for _, line in zip(range(500), f)]
        for row in rows:
            doc = json.loads(row["text"])
            p = parse(row["text"], tuple(doc))
            for n in doc:
                parts = p.objs[n]
                ring_sizes = [len(r) for o, hs in parts for r in (o, *hs)]
                args = (*pack_objects([(parts, row["W"], row["H"])]), GRID)
                ms = np.median([gc.check_polygons(*args, timing=True)["timing"]["validity"] for _ in range(3)])
                out.append({"dataset": rel.split("/")[1], "item_id": row["item_id"], "rater": row["rater"],
                            "object": n, "vertices": sum(ring_sizes), "max_ring": max(ring_sizes, default=0),
                            "validity_ms": float(ms)})
    big = [r for r in out if r["vertices"] >= 512]
    total = sum(r["validity_ms"] for r in out)
    result = {"build": gc.build(), "root": str(root), "inputs": inputs, "rungs": out,
              "share_ge512": len(big) / len(out), "time_share_ge512": sum(r["validity_ms"] for r in big) / total,
              "largest": sorted(out, key=lambda r: -r["vertices"])[:10]}
    write(path, result)


def write(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(result, f, indent=1)


def run_old(a, argv) -> int:
    """This script by path under the old engine's Python, this environment, --engine current."""
    if not Path(a.old_python).is_file():
        sys.exit(f"bench_checker: old engine not found: {a.old_python} (run_box.sh --gpu-check builds it from "
                 "old_engine.tar; Implementation.md D-G7)")
    rest, skip = [], False
    for x in argv:
        if skip:
            skip = False
        elif x in ("--engine", "--old-python"):
            skip = True
        elif not x.startswith(("--engine=", "--old-python=")):
            rest.append(x)
    return subprocess.run([a.old_python, "-I", str(Path(__file__).resolve()), *rest, "--engine", "current"],
                          env=os.environ).returncode


def main(argv=None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--engine", choices=("current", "old"), default="current")
    ap.add_argument("--old-python", default="/opt/old/bin/python")
    ap.add_argument("--workload", action="append", choices=WORKLOADS, help="repeatable; default all")
    ap.add_argument("--stages", action="store_true", help="per-stage device times (timing=True)")
    ap.add_argument("--warmup", type=int, default=2)
    ap.add_argument("--repeats", type=int, default=10)
    ap.add_argument("--dump-replay", type=Path, metavar="OUT")
    ap.add_argument("--synthetic", type=int, metavar="N")
    ap.add_argument("--replay", type=Path, metavar="FILE")
    ap.add_argument("--rungs", type=Path, metavar="OUT")
    ap.add_argument("--out", type=Path, required=not ({"--dump-replay", "--rungs"} & {x.split("=")[0] for x in argv}))
    a = ap.parse_args(argv)
    if a.engine == "old":
        sys.exit(run_old(a, argv))
    if a.dump_replay:
        return dump_replay(a.dump_replay, a.synthetic)
    if a.rungs:
        return rungs(a.rungs)
    a.workload = a.workload or list(WORKLOADS)
    head = {"build": gc.build(), "device": device_name(), "engine": a.engine, "root": str(data_root())}
    write(a.out, head | (replay(a) if a.replay else bench(a, data_root())))


if __name__ == "__main__":
    main()
