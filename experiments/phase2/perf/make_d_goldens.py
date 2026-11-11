"""Goldens of Phase D's CPU-vs-GPU checks, written on the Mac CPU build (one venv holds one build), LIDC only:

- tests/goldens/d_gpu02.json (GPU02, D-G3/D-G10 box half; tests/experiments/test_d_goldens.py compares on the box): per
  entry the SHA-256 of its inputs, then {key: dtype, shape, SHA-256 of tobytes()} of every output array. Entries:
  check_masks and topology on the 6,024 fuzz masks (fuzz_masks(2000) of seeds 0-2) and the 37,093 crops of
  checker_eval/masks.npz; check_polygons and certify_polygons on the fuzz polygons (fuzz_polygons(150, seed 0), 2,853
  cases, minus the 150 of the K1 "holes" family, whose generator rounds math.cos/sin with a platform-dependent last
  bit) and the 144,133 LIDC polygons of checker_eval/polygons.npz (certify in chunks of 20,000; each object against
  its own even-odd raster from check_polygons); score_polygons and eval_sets on the D1 replay set
  (bench_checker.dump_replay: prepared/lidc/val.jsonl and its caches); add_tables_from_raters and edt on the first
  1,000 caches of prepared/caches/lidc by sorted path. Section "gpu03": the sanitizer corpus of
  tests/geom/test_gpu_sanitize.py (GPU03 compares its GPU outputs with these).
- tests/goldens/d_r01_lidc.json (D-R01, always-run copy; tests/experiments/test_d_goldens.py): the first 50 items
  of the LIDC val store, framed as objdel.cli.data.targets frames them: per rater the mask of each object as
  run lengths and the target text stored in data/targets/lidc/val.jsonl (built before D).
- tests/goldens/d_r04_scorer.json (D-R04; tests/experiments/test_d_goldens.py): score_polygons and score_masks of 200
  synthetic items (d_r04_items) as the hex of tobytes() of every output.

Inputs resolve under OBJDEL_DATA_ROOT, else `data` (bench_checker.data_root). Reads untrusted data: run with -I.

    .venv/bin/python -I -m experiments.phase2.perf.make_d_goldens [--store data/store]
"""
import argparse
import functools
import hashlib
import json
import tempfile
from pathlib import Path

import numpy as np

from experiments.phase2.perf import bench_checker as B
from objdel.core.io import write_json
from objdel.geom import engine as gc
from objdel.geom.shape import pack_objects, to_mask, unpack_rows
from objdel.geom.text import GRID, parse
from objdel.verify import device as D
from objdel.verify.reward.cache import ExpertCache
from tests.oracles import cuda_parity as C

GOLDENS = Path(__file__).resolve().parents[3] / "tests/goldens"
POLYGONS, CHUNK, N_CACHES = "checker_eval/polygons.npz", 20000, 1000
POLY_KEYS = ("xy", "ring_start", "shape_start", "roles", "W", "H")
MIN_AREAS = np.array([0, 1, 3, 10], np.int64)  # topology's size rules, cycled over the masks


def record(a) -> dict:
    a = np.asarray(a)
    return {"dtype": str(a.dtype), "shape": list(a.shape),
            "sha256": hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()}


def outputs(r: dict) -> dict:
    """{key: record} of an engine result; table handles are process state, not results."""
    return {k: record(v) for k, v in sorted(r.items()) if k != "handles"}


def arrays_sha(*arrays) -> str:
    h = hashlib.sha256()
    for a in arrays:
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


# ------------------------------------------------------------------ GPU02 corpora (cached: several entries share one)

@functools.cache
def fuzz_masks() -> dict:
    return C.pack_masks([m for s in (0, 1, 2) for _, m in C.fuzz_masks(2000, seed=s)])


@functools.cache
def crops(root: Path) -> dict:
    with np.load(root / B.MASKS) as z:
        return {k: z[k] for k in ("bits", "hw", "start")}


@functools.cache
def fuzz_polygons() -> dict:
    """{G: packed polygons} of the fuzz cases without the "holes" family."""
    cases = [c for c in C.fuzz_polygons(150, seed=0) if c[0] != "holes"]
    return {g: C.pack_polygons([(r, rl, W, H) for _, r, rl, W, H, G in cases if G == g])
            for g in sorted({c[5] for c in cases})}


@functools.cache
def lidc_polygons(root: Path) -> dict:
    with np.load(root / POLYGONS) as z:
        return {k: z[k] for k in POLY_KEYS}


def chunk(p: dict, a: int, b: int) -> dict:
    """Shapes a..b-1 of a packed polygon set, re-based."""
    ss, rs = p["shape_start"], p["ring_start"]
    r0, r1 = int(ss[a]), int(ss[b])
    return {"xy": p["xy"][rs[r0]:rs[r1]], "ring_start": rs[r0:r1 + 1] - rs[r0], "shape_start": ss[a:b + 1] - r0,
            "roles": p["roles"][r0:r1], "W": p["W"][a:b], "H": p["H"][a:b]}


def certify_self(p: dict, G: int) -> tuple[str, dict]:
    """(SHA-256 of the references, certify_polygons of each shape against its own even-odd raster)."""
    ev = gc.check_polygons(*(p[k] for k in POLY_KEYS), G)
    ws = ev["word_start"]
    refs = C.pack_masks([unpack_rows(ev["evenodd"][ws[i]:ws[i + 1]], int(p["H"][i]), int(p["W"][i]))
                         for i in range(len(p["W"]))])
    n = len(p["W"])
    refs = (refs["bits"], refs["hw"], refs["start"])
    return arrays_sha(*refs), gc.certify_polygons(*(p[k] for k in POLY_KEYS), G, *refs, np.arange(n, dtype=np.int32),
                                                  MIN_AREAS[np.arange(n) % 4])


@functools.cache
def replay(root: Path) -> tuple[dict, dict]:
    """(inputs, arrays) of the D1 replay set, dumped by bench_checker.dump_replay under OBJDEL_DATA_ROOT = root."""
    with tempfile.TemporaryDirectory() as t:
        path = Path(t) / "replay_inputs.npz"
        B.dump_replay(path, None)
        with np.load(path) as z:
            d = {k: z[k] for k in z.files}
    files = json.loads(str(d["meta"]))["inputs"]
    caches = sorted((k, v["sha256"]) for k, v in files.items() if k.startswith(B.CACHES + "/"))
    return {B.VAL: files[B.VAL]["sha256"], B.CACHES: arrays_sha(np.array([f"{k} {s}" for k, s in caches]))}, d


def replay_tables(d: dict) -> list[int]:
    out = []
    for i in range(json.loads(str(d["meta"]))["caches"]):
        t = [d[f"t{i}_{k}"] for k in B.TABLE_ARGS]
        out.append(gc.add_tables(*t[:9], int(t[9]), int(t[10]), float(t[11]), int(t[12])))
    return out


def with_replay(root: Path, call) -> tuple[dict, dict]:
    inputs, d = replay(root)
    handles = replay_tables(d)
    try:
        return inputs, call(d, np.array(handles, np.int64)[d["obj_cache"]])
    finally:
        for h in handles:
            gc.drop_tables(h)


def score_replay(d: dict, tables) -> dict:
    return gc.score_polygons(*(d[k] for k in POLY_KEYS), int(d["GRID"]), tables)


def eval_replay(d: dict, tables) -> dict:
    """Each replay object is a set of K = 1 sample (legal), cuts at 0.1/0.5/0.9, NSD tolerance 2 px, band 2 px. One call
    per rater count (the engine wants one count per call); outputs are named "<key>@k<raters>"."""
    ks = np.array([len(d[f"t{c}_raters"]) for c in d["obj_cache"]])
    out = {}
    for k in sorted(set(ks.tolist())):
        sel = np.flatnonzero(ks == k)
        p = chunk_of(d, sel)
        n = len(sel)
        cuts = np.array([[[D.soft_cut(1, t), D.soft_cut(k, t)] for t in (0.1, 0.5, 0.9)]] * n, np.int32)
        r = gc.eval_sets(*(p[x] for x in POLY_KEYS), int(d["GRID"]), tables[sel], 1, np.ones(n, np.uint8), cuts,
                         np.full(n, D.nsd_cut(2.0), np.int64), np.full(n, 2, np.int32))
        out |= {f"{key}@k{k}": v for key, v in r.items()}
    return out


def chunk_of(p: dict, idx) -> dict:
    """The shapes idx (ascending) of a packed polygon set, re-based."""
    ss, rs = p["shape_start"], p["ring_start"]
    rings = np.concatenate([np.arange(ss[i], ss[i + 1]) for i in idx]) if len(idx) else np.zeros(0, np.int64)
    pts = [np.arange(rs[r], rs[r + 1]) for r in rings]
    counts = np.array([len(x) for x in pts], np.int64)
    per = np.array([ss[i + 1] - ss[i] for i in idx], np.int64)
    return {"xy": p["xy"][np.concatenate(pts)] if pts else p["xy"][:0],
            "ring_start": np.r_[0, np.cumsum(counts)].astype(rs.dtype),
            "shape_start": np.r_[0, np.cumsum(per)].astype(ss.dtype),
            "roles": p["roles"][rings], "W": p["W"][idx], "H": p["H"][idx]}


def cache_paths(root: Path) -> list[Path]:
    return sorted((root / B.CACHES).rglob("*.npz"))[:N_CACHES]


@functools.cache
def caches(root: Path) -> tuple[dict, list]:
    paths = cache_paths(root)
    sums = [f"{p.relative_to(root)} {B.file_sha256(p)}" for p in paths]
    return {B.CACHES + "[:1000]": arrays_sha(np.array(sums))}, [ExpertCache.load(p) for p in paths]


def tables_of(root: Path) -> tuple[dict, dict]:
    inputs, cs = caches(root)
    r = gc.add_tables_from_raters([c.raters.view(np.uint8) for c in cs], np.array([c.t5 for c in cs], np.float64),
                                  np.array([c.min_area for c in cs], np.int64))
    for h in r["handles"]:
        gc.drop_tables(int(h))
    return inputs, r


def edt_maps(cs) -> list[np.ndarray]:
    """The two distance maps of every cache: outside the union, inside the core."""
    return [m for c in cs for m in (~c.union, c.core)]


def edt_of(root: Path) -> tuple[dict, dict]:
    inputs, cs = caches(root)
    return inputs, {"edt": np.concatenate([gc.edt(m.view(np.uint8)).ravel() for m in edt_maps(cs)])}


def masks_entry(corpus, call):
    def run(root):
        p = corpus(root) if corpus is crops else corpus()
        return {"packed": arrays_sha(p["bits"], p["hw"], p["start"])}, call(p)
    return run


def check(p):
    return gc.check_masks(p["bits"], p["hw"], p["start"])


def topo(p):
    return gc.topology(p["bits"], p["hw"], p["start"], MIN_AREAS[np.arange(len(p["hw"])) % 4])


def polygon_entry(p_of, G: int, certify: bool):
    def run(root):
        p = p_of(root)
        inputs = {"packed": arrays_sha(*(p[k] for k in POLY_KEYS))}
        if not certify:
            return inputs, gc.check_polygons(*(p[k] for k in POLY_KEYS), G)
        refs, r = certify_self(p, G)
        return inputs | {"refs": refs}, r
    return run


def entries() -> dict:
    """{name: run(root) -> (inputs {name: SHA-256}, engine result)} of every GPU02 entry."""
    out = {"check_masks/fuzz": masks_entry(fuzz_masks, check), "check_masks/crops": masks_entry(crops, check),
           "topology/fuzz": masks_entry(fuzz_masks, topo), "topology/crops": masks_entry(crops, topo)}
    for g in fuzz_polygons():
        for kind in ("check_polygons", "certify_polygons"):
            out[f"{kind}/fuzz_G{g}"] = polygon_entry(lambda _, g=g: fuzz_polygons()[g], g, kind == "certify_polygons")
    out["check_polygons/lidc"] = polygon_entry(lidc_polygons, GRID, False)
    n = 144133  # the payload's polygon count; a different file fails on its input SHA-256 first
    for a in range(0, n, CHUNK):
        out[f"certify_polygons/lidc_{a}"] = polygon_entry(
            lambda root, a=a: chunk(lidc_polygons(root), a, min(a + CHUNK, n)), GRID, True)
    out |= {"score_polygons/replay": lambda root: with_replay(root, score_replay),
            "eval_sets/replay": lambda root: with_replay(root, eval_replay),
            "add_tables_from_raters/caches": tables_of, "edt/caches": edt_of}
    return out


# ------------------------------------------------------------------ D-R01: 50 LIDC items, masks and pre-D texts

def rle(m: np.ndarray) -> list[int]:
    """Run lengths of the row-major mask, starting with a run of zeros (possibly empty)."""
    f = m.ravel()
    runs = np.diff(np.r_[0, np.flatnonzero(f[1:] != f[:-1]) + 1, f.size]).tolist()
    return [0] + runs if f.size and f[0] else runs


def unrle(runs, H: int, W: int) -> np.ndarray:
    return np.repeat(np.arange(len(runs)) % 2 == 1, runs).reshape(H, W)


def d_r01(store_root: Path, root: Path) -> dict:
    from objdel.data import store
    from objdel.data.targets.build import build_target, min_areas
    from objdel.data.targets.frame import frame

    with open(root / "targets/lidc/val.jsonl") as f:
        stored = {(r["item_id"], r["rater"]): r["text"] for r in map(json.loads, f) if "text" in r}
    items = []
    for it in store.load("lidc", "val", store_root):
        if not it.item_id.startswith("lidc/") or "riga" in it.item_id.lower():
            raise ValueError(f"D-R01 golden holds LIDC items only, got {it.item_id}")
        f = frame(it)
        ma = min_areas(f.masks)
        raters = []
        for r, rater in enumerate(it.raters):
            text = stored.get((it.item_id, rater))
            if text is None:
                continue
            masks = {"nodule": f.masks["nodule"][r]}
            if build_target(masks, ("nodule",), None, ma).text != text:
                raise ValueError(f"{it.item_id} rater {rater}: rebuilt text differs from the stored one")
            raters.append({"rater": rater, "nodule": rle(masks["nodule"]), "text": text})
        if raters:
            H, W = f.masks["nodule"].shape[1:]
            items.append({"item_id": it.item_id, "H": int(H), "W": int(W), "min_area": ma, "raters": raters})
        if len(items) == 50:
            return {"items": items}
    raise ValueError("fewer than 50 LIDC val items with stored texts")


# ------------------------------------------------------------------ D-R04: 200 scorer items

KINDS = ("nothing", "small", "hole", "everything")  # what the answer's size filter removes


def _sq(x0, y0, x1, y1, hole=False) -> list:
    r = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
    return r[::-1] if hole else r


def d_r04_items() -> list[tuple[str, list, ExpertCache]]:
    """200 (kind, parts, cache) on 40 x 40 frames (25 grid units per pixel): a cache of 3 disc raters (tau and T5
    fixed, so no host float enters) and an answer whose size filter (min_area 5 here) removes nothing (one square),
    1-4 px specks beside the square (in pixels 0-7, the square starts at 10+), fills a 1-4 px hole of the square, or removes everything (specks only)."""
    rng = np.random.default_rng(20261110)
    y, x = np.mgrid[:40, :40]
    out = []
    for i in range(200):
        kind = KINDS[i % 4]
        cy, cx, r = rng.uniform(14, 26, 3)
        c = ExpertCache(np.stack([np.hypot(y - cy, x - cx) <= r * s for s in (0.6, 0.65, 0.7)]), 2.0, 2.0)
        a = int(rng.integers(10, 19))
        a, b = a * 25, int(rng.integers(a + 6, 35)) * 25  # a square of 6+ px
        cells = rng.choice(9, 3, replace=False)  # specks of 1-4 px on a 3-px lattice in px 0-7: never touching
        specks = [(int(v) % 3 * 75, int(v) // 3 * 75, int(rng.integers(1, 3)) * 25) for v in cells]
        main = {"outer": _sq(a, a, b, b), "holes": []}
        if kind == "hole":
            h = (a + b) // 50 * 25
            main["holes"] = [_sq(h, h, h + int(rng.integers(1, 3)) * 25, h + 25, hole=True)]
        parts = [] if kind == "everything" else [main]
        if kind in ("small", "everything"):
            parts += [{"outer": _sq(sx, sy, sx + s, sy + s), "holes": []} for sx, sy, s in specks]
        text = json.dumps({"nodule": parts}, separators=(",", ":"))
        out.append((kind, parse(text, ("nodule",)).objs["nodule"], c))
    return out


def d_r04_scores(items) -> dict:
    """{"score_polygons": result, "score_masks": result} of the items, each against its own cache."""
    tables = np.array([c.device() for _, _, c in items], np.int64)
    objs = [(parts, 40, 40) for _, parts, _ in items]
    masks = [to_mask(*o).mask for o in objs]
    p = C.pack_masks(masks)
    return {"score_polygons": gc.score_polygons(*pack_objects(objs), GRID, tables),
            "score_masks": gc.score_masks(p["bits"], p["hw"], p["start"], tables)}


def floats(r: dict) -> dict:
    return {k: {"dtype": str(v.dtype), "shape": list(v.shape), "hex": np.ascontiguousarray(v).tobytes().hex()}
            for k, v in sorted(r.items())}


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--store", type=Path, default=Path("data/store"))
    a = ap.parse_args(argv)
    if gc.build() != "cpu":
        raise SystemExit("make_d_goldens: the goldens come from the CPU build")
    from tests.geom.test_gpu_sanitize import CORPUS

    root = B.data_root()
    gpu02 = {}
    for name, run in entries().items():
        inputs, r = run(root)
        gpu02[name] = {"inputs": inputs, "outputs": outputs(r)}
        print(name, flush=True)
    write_json(GOLDENS / "d_gpu02.json", {"entries": gpu02, "gpu03": {n: outputs(f()) for n, f in CORPUS.items()}})
    write_json(GOLDENS / "d_r01_lidc.json", d_r01(a.store, root))
    write_json(GOLDENS / "d_r04_scorer.json", {k: floats(v) for k, v in d_r04_scores(d_r04_items()).items()})


if __name__ == "__main__":
    main()
