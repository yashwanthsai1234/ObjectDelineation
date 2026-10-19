"""Batch API correctness and K10 determinism (runs, processes, thread counts, pools)."""
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

import polycore as pc
from rollouts import batch_args, make_rollouts, stats_hash

ROLLS = make_rollouts(256, seed=7)


def test_batch_equals_single_calls():
    xys, offs, Ws, Hs, roles, gts = batch_args(ROLLS)
    b = pc.score_batch(xys, offs, Ws, Hs, roles, gts, threads=4)
    for i, r in enumerate(ROLLS):
        ras = pc.rasterize((r["xy"], r["off"]), r["S"], r["S"])
        val = pc.validate((r["xy"], r["off"]), r["roles"])
        gt = r["gt"].astype(bool)
        assert b["area_evenodd"][i] == ras["evenodd"].sum() and b["area_nonzero"][i] == ras["nonzero"].sum()
        assert b["inter_evenodd"][i] == (ras["evenodd"] & gt).sum()
        assert b["inter_nonzero"][i] == (ras["nonzero"] & gt).sum()
        assert b["gt_area"][i] == gt.sum()
        assert (b["defects"][i] == val["total"]).all()
        s = pc.score(r["xy"], r["off"], r["S"], r["S"], r["roles"], r["gt"])
        assert s["inter_evenodd"] == b["inter_evenodd"][i] and (s["defects"] == b["defects"][i]).all()


def _batch_hash(threads):
    return stats_hash(pc.score_batch(*batch_args(ROLLS), threads=threads))


def test_k10_thread_counts_and_repeats():
    ref = _batch_hash(1)
    for t in (0, 1, 2, 3, 4, 8, 16):
        for _ in range(2):
            assert _batch_hash(t) == ref


def test_k10_python_thread_pool():
    xys, offs, Ws, Hs, roles, gts = batch_args(ROLLS)
    ref = pc.score_batch(xys, offs, Ws, Hs, roles, gts, threads=1)
    for nt in (2, 8):
        with ThreadPoolExecutor(nt) as ex:
            res = list(ex.map(lambda i: pc.score(xys[i], offs[i], Ws[i], Hs[i], roles[i], gts[i]), range(len(ROLLS))))
        assert [r["inter_evenodd"] for r in res] == ref["inter_evenodd"].tolist()
        assert np.array_equal(np.stack([r["defects"] for r in res]), ref["defects"])
        # full masks too
        with ThreadPoolExecutor(nt) as ex:
            masks = list(ex.map(lambda r: pc.rasterize((r["xy"], r["off"]), r["S"], r["S"])["nonzero"], ROLLS))
        for r, m in zip(ROLLS, masks):
            assert np.array_equal(m, pc.rasterize((r["xy"], r["off"]), r["S"], r["S"])["nonzero"])


def test_k10_across_processes():
    here = Path(__file__).resolve().parent
    code = ("import sys; sys.path.insert(0, %r); sys.dont_write_bytecode=True\n"
            "import polycore as pc\nfrom rollouts import batch_args, make_rollouts, stats_hash\n"
            "R = make_rollouts(256, seed=7)\nprint(stats_hash(pc.score_batch(*batch_args(R), threads=%d)))")
    hashes = {subprocess.run([sys.executable, "-c", code % (str(here), t)], capture_output=True, text=True,
                             check=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}).stdout.strip()
              for t in (1, 8)}
    assert hashes == {_batch_hash(1)}
