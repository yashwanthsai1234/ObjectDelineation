"""A6: per-image decode cost of the prepared PNG frames: PIL (what datasets.Image uses), OpenCV, torchvision
decode_png, and a RAM-resident uint8 array (copy only); file sizes. Interleaved per image.

    .venv/bin/python -I -m experiments.phase2.lab.A6.decode --out <dir>
"""
import argparse
import io
import time
from pathlib import Path

import cv2
import numpy as np
from PIL import Image

from .common import PREPARED, dump, env, records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=60)
    a = ap.parse_args()
    cv2.setNumThreads(0)
    try:
        import torch  # noqa: PLC0415
        import torchvision.io as tvio  # noqa: PLC0415
        torch.set_num_threads(1)
    except ImportError:  # pragma: no cover
        tvio = None
    res = {"env": env()}
    for ds in ("riga", "lidc"):
        recs = records(ds, "train")[:: max(1, len(records(ds, "train")) // a.n)][: a.n]
        paths = [PREPARED / "images" / (r["item_id"] + ".png") for r in recs]
        T = {"pil": [], "opencv": [], "torchvision": [], "ram_copy": []}
        sizes = [p.stat().st_size for p in paths]
        for p in paths:
            data = Path(p).read_bytes()
            t0 = time.perf_counter()
            im = np.asarray(Image.open(io.BytesIO(data)).convert("RGB"))
            t1 = time.perf_counter()
            cv2.imdecode(np.frombuffer(data, np.uint8), cv2.IMREAD_COLOR)
            t2 = time.perf_counter()
            if tvio is not None:
                tvio.decode_png(torch.frombuffer(bytearray(data), dtype=torch.uint8))
            t3 = time.perf_counter()
            im.copy()
            t4 = time.perf_counter()
            T["pil"].append(1e3 * (t1 - t0))
            T["opencv"].append(1e3 * (t2 - t1))
            T["torchvision"].append(1e3 * (t3 - t2))
            T["ram_copy"].append(1e3 * (t4 - t3))
        res[ds] = {"images": len(paths), "shape": str(im.shape), "png_kB_mean": float(np.mean(sizes)) / 1e3,
                   "decode_ms_median": {k: float(np.median(v)) for k, v in T.items()}}
        print(ds, res[ds], flush=True)
    res["env_end"] = env()
    dump(res, f"{a.out}/decode.json")


if __name__ == "__main__":
    main()
