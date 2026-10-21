"""A6 (d'): ExpertCache.load cost per item (the GRPO reward's cache_loader on an LRU miss) by stage, its memory,
and whether scipy's EDT releases the GIL (so a prefetch thread could overlap it with generation).

    .venv/bin/python -I -m experiments.phase2.lab.A6.prof_cache --out <dir>
"""
import argparse
import threading
import time
from pathlib import Path

import numpy as np
from scipy import ndimage

from objdel.reward import ExpertCache
from objdel.reward.cache import _band, _digest, edt
from objdel.reward.spec import fingerprint

from .common import PREPARED, Timer, dump, env, profile, records, stats_ms


def staged_load(path, digest, t: Timer):
    a = time.perf_counter()
    with np.load(Path(path), allow_pickle=False) as z:
        raters, tau, t5 = z["raters"].astype(bool), float(z["tau"]), float(z["t5"])
        fp, dg = str(z["fingerprint"]), str(z["sha256"])
    b = time.perf_counter()
    assert _digest(raters, tau, t5, fp) == dg == digest
    c = time.perf_counter()
    fingerprint()
    d = time.perf_counter()
    _band(raters)
    e = time.perf_counter()
    ExpertCache(raters, tau, t5)  # the whole constructor (band again, vote, signatures, box, digest)
    f = time.perf_counter()
    t.add("npz read+inflate", b - a)
    t.add("sha256 digest check", c - b)
    t.add("fingerprint()", d - c)
    t.add("_band (2 EDT + boundary)", e - d)
    t.add("ExpertCache.__init__ total", f - e)
    return raters.shape


def nbytes(c: ExpertCache) -> int:
    return sum(x.nbytes for x in (c.raters, c.union, c.core, c.vote, c.dist_out_union, c.dist_in_core))


def gil_test(shape=(500, 500), n=40):
    """Spinner thread's Python-loop rate alone vs while the main thread runs scipy EDTs."""
    stop = threading.Event()
    count = [0]

    def spin():
        while not stop.is_set():
            count[0] += 1

    m = np.zeros(shape, np.uint8)
    m[100:400, 120:380] = 1
    th = threading.Thread(target=spin)
    th.start()
    time.sleep(0.3)
    c0 = count[0]
    time.sleep(0.5)
    solo = (count[0] - c0) / 0.5
    c1 = count[0]
    t0 = time.perf_counter()
    for _ in range(n):
        ndimage.distance_transform_edt(m)
    dt = time.perf_counter() - t0
    during = (count[0] - c1) / dt
    stop.set()
    th.join()
    return {"spinner_solo_per_s": solo, "spinner_during_edt_per_s": during, "ratio": during / solo,
            "edt_ms_500": 1e3 * dt / n}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = {"env": env()}
    for ds, n in (("riga", 24), ("lidc", 200)):
        recs = [r for r in records(ds, "train") if r["rl"]][:n]
        t = Timer()
        full, mem, shapes = [], [], []
        for r in recs:
            for name, dg in r["caches"].items():
                p = PREPARED / "caches" / r["item_id"] / f"{name}.npz"
                t0 = time.perf_counter()
                c = ExpertCache.load(p, dg)
                full.append(1e3 * (time.perf_counter() - t0))
                mem.append(nbytes(c))
                shapes.append(c.raters.shape)
                staged_load(p, dg, t)
        k = len(full)
        res[ds] = {"objects": k, "load_ms_per_object": stats_ms(full),
                   "stages_ms_per_object": {x: v / k for x, v in t.t.items()},
                   "mem_MB_per_object": float(np.mean(mem)) / 2**20,
                   "shapes": sorted({str(s) for s in shapes})[:6]}
        rec = recs[0]
        _, rows, total = profile(lambda: [ExpertCache.load(PREPARED / "caches" / rec["item_id"] / f"{nm}.npz", d)
                                          for nm, d in rec["caches"].items()], top=12)
        res[ds]["cprofile_one_item"] = [(round(tt, 3), round(ct, 3), nc, nm) for tt, ct, nc, nm in rows]
        print(ds, res[ds]["load_ms_per_object"]["median"], res[ds]["stages_ms_per_object"], flush=True)
    # EDT alternatives on a RIGA-size mask (cropped union box vs full frame)
    m = np.zeros((500, 500), bool)
    yy, xx = np.mgrid[:500, :500]
    m[(yy - 250) ** 2 / 120**2 + (xx - 250) ** 2 / 110**2 < 1] = True
    res["edt_500_ms"] = {}
    for label, fn in (("scipy edt", lambda: edt(~m)), ("scipy edt core", lambda: edt(m))):
        ts = []
        for _ in range(20):
            t0 = time.perf_counter()
            fn()
            ts.append(1e3 * (time.perf_counter() - t0))
        res["edt_500_ms"][label] = float(np.median(ts))
    try:
        import cv2
        cv2.setNumThreads(0)
        ts = []
        for _ in range(20):
            t0 = time.perf_counter()
            cv2.distanceTransform((~m).view(np.uint8), cv2.DIST_L2, cv2.DIST_MASK_PRECISE)
            ts.append(1e3 * (time.perf_counter() - t0))
        res["edt_500_ms"]["opencv precise (not used: x86 IPP nondeterminism)"] = float(np.median(ts))
    except Exception as e:  # pragma: no cover
        res["edt_500_ms"]["opencv"] = str(e)
    res["gil"] = gil_test()
    print(res["edt_500_ms"], res["gil"], flush=True)
    res["env_end"] = env()
    dump(res, f"{a.out}/cache.json")


if __name__ == "__main__":
    main()
