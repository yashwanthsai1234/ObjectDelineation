"""Cross-check K1 (via the k1_dump binary) against objdel.checker.polygonize and objdel.checker.bitquads.

Run with the repo's venv and -I (the masks file is data, untrusted):
    REPO/.venv/bin/python -I cross_numpy.py <k1_dump binary> <workdir> [--masks REPO/data/checker_eval/masks.npz]
        [--limit N] [--random N] [--seed-rows K]
Writes <workdir>/in.{bits,hw,start}, runs k1_dump, and compares every mask: ring count, every ring's corners and
signed area, and the five bit-quad counts. Exit code 0 only if all are identical.
"""
import argparse
import dataclasses
import subprocess
import sys
from pathlib import Path

import numpy as np
from objdel.checker import bitquads, polygonize


def _quads(b):
    t = dataclasses.astuple(b) if dataclasses.is_dataclass(b) else tuple(b)
    return np.asarray(t, np.int64)


def random_masks(n, rng):
    out = []
    for k in range(n):
        h, w = int(rng.integers(1, 90)), int(rng.integers(1, 90))
        if k % 3 == 0:  # blobs: discs on a grid, minus smaller discs
            yy, xx = np.mgrid[0:h, 0:w]
            m = np.zeros((h, w), bool)
            for t in range(int(rng.integers(1, 8))):
                cx, cy, r = rng.uniform(0, w), rng.uniform(0, h), rng.uniform(1, 15)
                disc = (xx + 0.5 - cx) ** 2 + (yy + 0.5 - cy) ** 2 <= r * r
                m = m & ~disc if t % 3 == 2 else m | disc
        else:
            m = rng.random((h, w)) < rng.uniform(0.05, 0.95)
        out.append(m.astype(np.uint8))
    return out


def pack(masks):
    hw = np.array([m.shape for m in masks], np.int32).reshape(-1, 2)
    start = np.zeros(len(masks) + 1, np.int64)
    start[1:] = np.cumsum([m.size for m in masks])
    bits = np.concatenate([m.ravel() for m in masks]) if masks else np.zeros(0, np.uint8)
    return bits.astype(np.uint8), hw, start


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("workdir")
    ap.add_argument("--masks", default=None)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--random", type=int, default=0)
    ap.add_argument("--seed-rows", type=int, default=16)
    a = ap.parse_args()
    wd = Path(a.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    if a.masks:
        z = np.load(a.masks, allow_pickle=False)
        bits, hw, start = z["bits"], z["hw"].astype(np.int32), z["start"].astype(np.int64)
        if a.limit:
            n = min(a.limit, len(hw))
            hw, start = hw[:n], start[: n + 1]
            bits = bits[: start[-1]]
        bits = (bits != 0).astype(np.uint8)
    else:
        bits, hw, start = pack(random_masks(a.random or 2000, np.random.default_rng(7)))
    n = len(hw)
    bits.tofile(wd / "in.bits")
    hw.tofile(wd / "in.hw")
    start.tofile(wd / "in.start")
    subprocess.run([a.dump, str(wd / "in"), str(wd / "out.bin"), str(a.seed_rows)], check=True)

    raw = (wd / "out.bin").read_bytes()
    n2, R, P = np.frombuffer(raw[:24], np.int64)
    off = 24

    def take(dtype, count):
        nonlocal off
        v = np.frombuffer(raw[off: off + np.dtype(dtype).itemsize * count], dtype)
        off += np.dtype(dtype).itemsize * count
        return v

    quads = take(np.int64, 5 * n).reshape(n, 5)
    ring_start = take(np.int64, n + 1)
    pts_start = take(np.int64, R + 1)
    ring_area = take(np.int64, R)
    pts = take(np.int32, 2 * P).reshape(P, 2)
    assert n2 == n and off == len(raw)

    bad = 0
    rings = 0
    for m in range(n):
        H, W = int(hw[m, 0]), int(hw[m, 1])
        mask = bits[start[m]: start[m + 1]].reshape(H, W).astype(bool)
        want = polygonize(mask) if mask.size else []
        q = _quads(bitquads(mask)) if mask.size else np.zeros(5, np.int64)
        got_r = range(int(ring_start[m]), int(ring_start[m + 1]))
        ok = len(want) == len(got_r) and np.array_equal(q, quads[m])
        if ok:
            for ring, r in zip(want, got_r):
                g = pts[pts_start[r]: pts_start[r + 1]]
                if ring.area != ring_area[r] or not np.array_equal(ring.pts, g):
                    ok = False
                    break
        rings += len(want)
        if not ok:
            bad += 1
            if bad <= 5:
                print(f"mismatch at mask {m} ({H}x{W}): want {len(want)} rings, got {len(got_r)}", file=sys.stderr)
    print(f"cross_numpy: {n} masks, {rings} rings, {bad} mismatches -> {'PASS' if bad == 0 else 'FAIL'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
