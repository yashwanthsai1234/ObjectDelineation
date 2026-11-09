"""K1 runner: polycore.rasterize vs the oracle on deterministic generated cases.

Usable from pytest (test_k1_raster.py) and as a script for the exhaustive
pure-oracle run:  python tests/k1_runner.py --n 10000 --pure --workers 8
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(os.environ.get("POLYCORE_ORACLE_DIR", HERE.parents[1] / "common"))))
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True

import polycore as pc  # noqa: E402
from helpers import k1_case, on_boundary_centres, oracle_winding  # noqa: E402


def run_chunk(args):
    lo, hi, pure, seed = args
    st = collections.Counter()
    fails = []
    for i in range(lo, hi):
        fam, rings, W, H = k1_case(i, seed)
        res = pc.rasterize(rings, W, H, winding=True)
        w, c = oracle_winding(rings, W, H, pure=pure)
        bad_w = int((res["winding"] != w).sum())
        bad_eo = int((res["evenodd"] != (c % 2 == 1)).sum())
        bad_nz = int((res["nonzero"] != (w != 0)).sum())
        bad_ring = 0
        if len(rings) > 1:
            for k, ring in enumerate(rings):
                wr, _ = oracle_winding([ring], W, H, pure=pure)
                bad_ring += int(min(0, wr.min()) != res["ring_wmin"][k]) + int(max(0, wr.max()) != res["ring_wmax"][k])
        else:
            bad_ring += int(min(0, w.min()) != res["ring_wmin"][0]) + int(max(0, w.max()) != res["ring_wmax"][0])
        st["cases"] += 1
        st[f"fam:{fam}"] += 1
        st["pixels"] += W * H
        st["bad_winding_px"] += bad_w
        st["bad_evenodd_px"] += bad_eo
        st["bad_nonzero_px"] += bad_nz
        st["bad_ring_flags"] += bad_ring
        st["mixed_sign_rings"] += int(((res["ring_wmin"] < 0) & (res["ring_wmax"] > 0)).sum())
        st["eo_ne_nz_cases"] += int(((c % 2 == 1) != (w != 0)).any())
        if i % 10 == 0:  # tie coverage on a 10 % sample (exact on-edge centre test)
            st["sampled_on_edge_centres"] += int(on_boundary_centres(rings, W, H).sum())
        if bad_w or bad_eo or bad_nz or bad_ring:
            fails.append(i)
    return dict(st), fails


def run(n, pure=False, workers=None, seed=12345, start=0, chunk=50):
    workers = workers or os.cpu_count()
    jobs = [(a, min(a + chunk, start + n), pure, seed) for a in range(start, start + n, chunk)]
    total = collections.Counter()
    fails = []
    with ProcessPoolExecutor(workers) as ex:
        for st, f in ex.map(run_chunk, jobs):
            total.update(st)
            fails += f
    return dict(total), fails


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--pure", action="store_true", help="use pure-Python oracle.py (slow, the spec itself)")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    t = time.time()
    st, fails = run(a.n, a.pure, a.workers, a.seed)
    st["seconds"] = round(time.time() - t, 1)
    st["oracle"] = "oracle.py (pure)" if a.pure else "numpy vectorisation of oracle.py"
    st["failing_cases"] = fails[:50]
    print(json.dumps(st, indent=1, sort_keys=True))
    if a.out:
        Path(a.out).write_text(json.dumps(st, indent=1, sort_keys=True))
    sys.exit(1 if fails else 0)
