"""Goldens of the scorer sums gate (B3.1, B.9.5 G01): the scorer's results on fixed inputs, written on the Mac CPU build
and compared on the box (tests/test_scorer_sums.py), since one venv holds one build.

Sets (`b_sets`): `random` = 200 random masks plus all-zero, all-one and one-pixel 64 x 64 masks; `polygons_<G>` = the
first 10,000 fuzz polygons of experiments/cuda_checker/compare.py but the K1 "holes" family (`polygon_cases`), one
score_polygons call per grid G; `windows_1000` =
the first 1,000 crops of checker_eval/masks.npz. Tables rule: every object m is scored against
`b_tables(m)` = ExpertCache.build of three rolls of m, three raters in m's own frame (m itself is not a rater).
The golden holds per set the SHA-256 of the packed call inputs and {field: [dtype, shape, sha256]} of every result field
except `sum_x`, `sum_y`, which the test checks against np.nonzero sums instead.

    .venv/bin/python -I -m experiments.phase2.perf.make_b_goldens [--masks PATH] [--out tests/goldens/b_g01.json]
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path

import numpy as np

from experiments.cuda_checker import compare as C
from objdel import _gpucheck as G
from objdel.geometry.shape import unpack_rows
from objdel.reward import ExpertCache

ROOT = Path(__file__).resolve().parents[3]
SUMS = ("sum_x", "sum_y")


def b_tables(m) -> ExpertCache:
    """G01 tables rule: three raters in m's own frame, m not one of them."""
    return ExpertCache.build(np.stack([np.roll(m, 1, axis=1), np.roll(m, 1, axis=0), np.roll(m, -1, axis=1)]))


def random_masks(rng):
    yield np.zeros((64, 64), bool)
    yield np.ones((64, 64), bool)
    one = np.zeros((64, 64), bool)
    one[63, 63] = True
    yield one
    for _ in range(200):
        h, w = rng.integers(8, 513, 2)
        yield rng.random((h, w)) < rng.random()


def _sha(*arrays) -> str:
    h = hashlib.sha256()
    for a in arrays:
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def score_masks_of(masks, packed=None) -> tuple[str, dict]:
    """(input SHA-256, score_masks result) of the masks, each against its b_tables, in the layout of
    objdel.reward.score.score_mask (one byte per pixel, hw int32, start int64)."""
    if packed is None:
        packed = (np.concatenate([np.ascontiguousarray(m, bool).view(np.uint8).reshape(-1) for m in masks]),
                  np.array([m.shape for m in masks], np.int32),
                  np.concatenate([[0], np.cumsum([m.size for m in masks])]).astype(np.int64))
    caches = [b_tables(m) for m in masks]  # alive until the call returns (tables freed on GC)
    return _sha(*packed), G.score_masks(*packed, np.array([c.device() for c in caches], np.int64))


def polygon_cases() -> list:
    """The fuzz cases of the `polygons_<G>` sets: the first 10,000 without the K1 "holes" family, whose generator rounds
    math.cos/sin (last bit platform-dependent: 7 of its 530 cases differ on the box's glibc, compare.polygon_reference)."""
    return [c for c in C.fuzz_polygons(530, seed=0)[:10000] if c[0] != "holes"]


def b_sets(masks_npz, n_masks=203, n_polys=None, n_crops=1000) -> dict:
    """{name: (input SHA-256, objects m, result)}; the defaults are G01's sets, G07 passes 32, 64, 32."""
    out = {}
    ms = list(itertools.islice(random_masks(np.random.default_rng(1)), n_masks))
    sha, r = score_masks_of(ms)
    out["random"] = (sha, ms, r)
    cases = polygon_cases()
    for g in sorted({c[5] for c in cases}):  # G = 7, 1000, 2^30: score_polygons takes one G
        sub = [c for c in cases if c[5] == g][:n_polys]
        p = C.pack_polygons([(r, rl, W, H) for _, r, rl, W, H, _ in sub])
        args = (p["xy"], p["ring_start"], p["shape_start"], p["roles"], p["W"], p["H"])
        ch = G.check_polygons(*args, g)
        ms = [unpack_rows(ch["evenodd"][ch["word_start"][i]:ch["word_start"][i + 1]], H, W)
              for i, (_, _, _, W, H, _) in enumerate(sub)]
        caches = [b_tables(m) for m in ms]
        out[f"polygons_{g}"] = (_sha(*args, np.int64(g)), ms,
                                G.score_polygons(*args, g, np.array([c.device() for c in caches], np.int64)))
    out["windows_1000"] = windows(masks_npz, n_crops)
    return out


def windows(masks_npz, n_crops=1000) -> tuple:
    """(input SHA-256, the first n_crops crops of checker_eval/masks.npz, score_masks result): G01's `windows_1000`."""
    with np.load(masks_npz) as z:  # bits: one byte (0/1) per pixel; hw int32; start int64
        st = z["start"][:n_crops + 1]
        packed = (z["bits"][:st[n_crops]], z["hw"][:n_crops], st)
    ms = [packed[0][st[k]:st[k + 1]].view(bool).reshape(packed[1][k]) for k in range(n_crops)]
    sha, r = score_masks_of(ms, packed)
    return sha, ms, r


def fields(r: dict) -> dict:
    """{field: [dtype, shape, sha256]} of every result field except the sums."""
    out = {}
    for k in sorted(set(r) - set(SUMS)):
        a = np.ascontiguousarray(r[k])
        out[k] = [str(a.dtype), list(a.shape), hashlib.sha256(a.tobytes()).hexdigest()]
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    data = Path(os.environ.get("OBJDEL_DATA_ROOT", ROOT / "data"))
    ap.add_argument("--masks", type=Path, default=data / "checker_eval" / "masks.npz")
    ap.add_argument("--out", type=Path, default=ROOT / "tests" / "goldens" / "b_g01.json")
    a = ap.parse_args(argv)
    if G.build() != "cpu":
        raise SystemExit("make_b_goldens writes the CPU-build golden: run it on the Mac build")
    sets = b_sets(a.masks)
    gold = {name: {"input_sha256": sha, "fields": fields(r)} for name, (sha, _, r) in sets.items()}
    a.out.write_text(json.dumps(gold, indent=1, sort_keys=True) + "\n")
    print(f"wrote {a.out}: " + ", ".join(f"{name} ({len(ms)} objects)" for name, (_, ms, _) in sets.items()))


if __name__ == "__main__":
    main()
