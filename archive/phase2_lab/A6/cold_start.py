"""A6: cold-start costs on a fresh box: numba JIT of objdel.targets.simplify with an empty cache (run with
NUMBA_CACHE_DIR=<empty dir>), then the warm call; and scipy EDT in 1 vs 2 Python threads (does it release the GIL?).

    NUMBA_CACHE_DIR=$(mktemp -d) .venv/bin/python -I -m experiments.phase2.lab.A6.cold_start --out <dir>
"""
import argparse
import os
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from .common import dump, env


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    res = {"env": env(), "numba_cache_dir": os.environ.get("NUMBA_CACHE_DIR")}
    t0 = time.perf_counter()
    from objdel.targets import build_target, min_areas  # noqa: PLC0415 - timed import
    res["import_targets_s"] = time.perf_counter() - t0
    yy, xx = np.mgrid[:500, :500]
    disc = (yy - 250) ** 2 / 120**2 + (xx - 250) ** 2 / 110**2 < 1
    cup = (yy - 250) ** 2 / 50**2 + (xx - 255) ** 2 / 45**2 < 1
    for label in ("first_call_s", "second_call_s"):
        t0 = time.perf_counter()
        build_target({"disc": disc, "cup": cup}, ("disc", "cup"), None, min_areas({"disc": disc[None], "cup": cup[None]}))
        res[label] = time.perf_counter() - t0
    from scipy import ndimage  # noqa: PLC0415
    m = ~disc
    jobs = [m] * 32
    t0 = time.perf_counter()
    for j in jobs:
        ndimage.distance_transform_edt(j)
    serial = time.perf_counter() - t0
    with ThreadPoolExecutor(2) as ex:
        t0 = time.perf_counter()
        list(ex.map(ndimage.distance_transform_edt, jobs))
        two = time.perf_counter() - t0
    res["edt_32x500_serial_s"], res["edt_32x500_2threads_s"] = serial, two
    res["edt_thread_speedup"] = serial / two
    res["env_end"] = env()
    print(res, flush=True)
    dump(res, f"{a.out}/cold_start{'_warm' if 'warm' in (os.environ.get('NUMBA_CACHE_DIR') or '') else ''}.json")


if __name__ == "__main__":
    main()
