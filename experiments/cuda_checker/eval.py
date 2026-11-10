"""Eval of the batched checker objdel._gpucheck against its references: the CPU checker it replaced (polycore, archived
in archive/cpu_checker; its results recorded) for polygons, the numpy polygonizer (tests/oracles/polygonize_oracle.py)
and bitquad (run here) for masks.

Run (Mac: CPU build; CUDA box: after `make cuda`, the GPU build):
    .venv/bin/python -I -m experiments.cuda_checker.eval [--data data/checker_eval] [--quick] [--no-speed]

Correctness, every output the reward and certify use: legality verdict and the 13 defect counts, the even-odd mask,
bit-quad counts, rings (Betti numbers b0/b1 and the ring hash), on
  (a) seeded fuzz cases from the repo's generators (compare.fuzz_polygons / fuzz_masks: K1 raster families, K2
      rings, multi-ring shapes with roles, rollouts, int32-extreme coordinates, shapes at the caps; ten mask families
      and edge cases), checked against polycore's results recorded in tests/goldens/polycore.npz and against
      polygonize run here;
  (b) every LIDC polygon of data/checker_eval/polygons.npz (training, eval and probe targets, R1 model outputs);
  (c) every LIDC Phase 0 mask of data/checker_eval/masks.npz;
  (b) and (c) are checked against the references stored by prepare_data.py, after verifying the files' SHA-256 against
  manifest.json. Every mismatch is listed with its category and first differing field, and saved to
  data/checker_eval/repro/<build>/ (git-ignored: may hold LIDC geometry); a mask is first shrunk greedily to the
  smallest one that still differs (a polygon cannot be: its sub-cases have no reference since polycore is archived).
Speed: shapes/s of this build, and masks/s against the Python certify path and OpenCV, by batch size and shape size.

Writes experiments/cuda_checker/results/<build>.json and the build's section of experiments/cuda_checker/RESULTS.md.
Licence rule: RIGA must never reach a rented host, so any path naming RIGA, or a manifest that is not LIDC-only, is
refused.
"""
import argparse
import collections
import hashlib
import json
import os
import platform
import random
import re
import shutil
import subprocess
import time
from pathlib import Path

import cv2
import numpy as np

from objdel.geom import engine as _gpucheck
from objdel.geom.topology import betti
from tests.oracles import cuda_parity as C
from tests.oracles.bitquad_oracle import bitquads
from tests.oracles.polygonize_oracle import polygonize

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PIXELS_PER_CALL = 1 << 26  # polygon batches are split at 2^26 pixels per call (8 MB of packed mask words)
SRC = ("train_target", "eval_target", "probe_target", "model_output")  # polygons.npz src codes


def refuse_riga(*paths):
    for p in paths:
        if "riga" in str(Path(p).resolve()).lower():
            raise SystemExit(f"refused: {p} names RIGA, whose licence forbids rented hosts")


def check_data(data: Path) -> dict:
    """manifest.json: LIDC only, and the SHA-256 of masks.npz / polygons.npz as recorded."""
    man = json.loads((data / "manifest.json").read_text())
    refuse_riga(*man.get("inputs", {}), *man.get("masks", {}).get("per_file", {}))
    if "no RIGA" not in man.get("licence", ""):
        raise SystemExit(f"refused: {data}/manifest.json does not declare the set LIDC-only")
    for name, want in man["files"].items():
        h = hashlib.sha256((data / name).read_bytes()).hexdigest()
        if h != want:
            raise SystemExit(f"{data / name}: SHA-256 {h} != manifest {want}")
    return man


# ------------------------------------------------------------------ shrinking a mask mismatch
def mask_diff_one(m):
    return C.mask_diff(C.run_masks(C.pack_masks([m])), 0, C.mask_reference(m))


def shrink_mask(m):
    """Greedy: drop border rows and columns while the GPU still differs from polygonize/bitquad."""
    changed = True
    while changed and m.size:
        changed = False
        for cut in (np.s_[1:, :], np.s_[:-1, :], np.s_[:, 1:], np.s_[:, :-1]):
            if mask_diff_one(m[cut]):
                m, changed = m[cut], True
                break
    return m


def save_repro(out_dir: Path, name: str, case: dict):
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{name}.json").write_text(json.dumps(case, default=lambda o: np.asarray(o).tolist()) + "\n")


# ------------------------------------------------------------------ correctness
class Tally:
    """Counts per category, mismatches per (category, first differing field), and the first repro of each."""

    def __init__(self, name, repro_dir, max_repro=10):
        self.name, self.repro_dir, self.max_repro = name, repro_dir, max_repro
        self.n, self.bad = collections.Counter(), collections.Counter()
        self.examples, self.saved = [], collections.Counter()

    def add(self, cat, field, index, repro):
        self.n[cat] += 1
        if field is None:
            return
        self.bad[(cat, field)] += 1
        if self.saved[(cat, field)] < self.max_repro:
            self.saved[(cat, field)] += 1
            case = repro()
            save_repro(self.repro_dir, f"{self.name}-{cat}-{field}-{index}", case)
            self.examples.append({"index": index, "category": cat, "field": field, "repro_size": case["size"]})

    def summary(self):
        cats = sorted(self.n)
        rows = {c: {"n": self.n[c], "mismatches": sum(v for (k, _), v in self.bad.items() if k == c)} for c in cats}
        total = sum(self.n.values())
        bad = sum(self.bad.values())
        return {"n": total, "mismatches": bad, "exact_rate": (total - bad) / total if total else None,
                "per_category": rows, "first_field": {f"{c}/{f}": v for (c, f), v in sorted(self.bad.items())},
                "examples": self.examples}


def polygon_batches(n, W, H):
    """Shape index ranges whose masks fit PIXELS_PER_CALL."""
    px = np.asarray(W, np.int64) * np.asarray(H, np.int64)
    a = 0
    while a < n:
        b, tot = a, 0
        while b < n and (b == a or tot + px[b] <= PIXELS_PER_CALL):
            tot += px[b]
            b += 1
        yield a, b
        a = b


def subset(p, a, b):
    r0, r1 = p["shape_start"][a], p["shape_start"][b]
    v0, v1 = p["ring_start"][r0], p["ring_start"][r1]
    return dict(xy=p["xy"][v0:v1], ring_start=p["ring_start"][r0:r1 + 1] - v0,
                shape_start=p["shape_start"][a:b + 1] - r0, roles=p["roles"][r0:r1], W=p["W"][a:b], H=p["H"][a:b])


def eval_fuzz_polygons(tally, n_per_family, seeds):
    for seed in seeds:
        cases = C.fuzz_polygons(n_per_family, seed)
        for G in sorted({c[5] for c in cases}):
            sub = [c for c in cases if c[5] == G]
            p = C.pack_polygons([(r, rl, W, H) for _, r, rl, W, H, _ in sub])
            for a, b in polygon_batches(len(sub), p["W"], p["H"]):
                q = subset(p, a, b)
                out = C.run_polygons(q, G)
                for s in range(b - a):
                    fam, rings, roles, W, H, _ = sub[a + s]
                    roles = [-1] * len(rings) if roles is None else list(roles)
                    ref = C.polygon_reference(rings, roles, W, H, G)
                    d = C.polygon_diff(C.gpu_polygon(out, q, s), ref)
                    tally.add(fam, d, f"s{seed}-{a + s}", lambda rings=rings, roles=roles, W=W, H=H, G=G, ref=ref:
                              repro_polygon(rings, roles, W, H, G, ref))


def repro_polygon(rings, roles, W, H, G, ref):
    p = C.pack_polygons([(rings, roles, W, H)])
    gpu = C.gpu_polygon(C.run_polygons(p, G), p, 0)
    return {"rings": [np.asarray(r).tolist() for r in rings], "roles": list(roles), "W": W, "H": H, "G": G,
            "size": sum(map(len, rings)), "gpu_defects": gpu["defects"], "ref_defects": ref["defects"],
            "evenodd_equal": C.mask_sha(gpu["evenodd"]) == ref["evenodd"]}


def eval_polygon_set(tally, p, per_class):
    """polygons.npz against its stored references; per_class counts cases and mismatches per DEFECTS class."""
    n = len(p["W"])
    for a, b in polygon_batches(n, p["W"], p["H"]):
        q = subset(p, a, b)
        out = C.run_polygons(q)
        for s in range(b - a):
            i = a + s
            ref = {"defects": p["defects"][i], "evenodd": p["evenodd_sha"][i]}
            d = C.polygon_diff(C.gpu_polygon(out, q, s), ref)
            classes = [c for c, v in zip(_gpucheck.DEFECTS, ref["defects"], strict=True) if v] or ["legal"]
            for c in classes:
                per_class[c][0] += 1
                per_class[c][1] += d is not None
            rings, roles, W, H = C.shape_of(p, i)
            tally.add(SRC[p["src"][i]], d, i, lambda rings=rings, roles=roles, W=W, H=H, ref=ref:
                      repro_polygon(rings, roles, W, H, 1000, ref))


def repro_mask(m):
    m = shrink_mask(m)
    return {"mask": ["".join(".#"[v] for v in row) for row in m.astype(int)], "size": int(m.size)}


def eval_fuzz_masks(tally, n, seeds):
    for seed in seeds:
        cases = C.fuzz_masks(n, seed)
        p = C.pack_masks([m for _, m in cases])
        out = C.run_masks(p)
        for i, (fam, m) in enumerate(cases):
            m = np.asarray(m, bool)
            tally.add(fam, C.mask_diff(out, i, C.mask_reference(m)), f"s{seed}-{i}", lambda m=m: repro_mask(m))


def eval_mask_set(tally, z):
    out = C.run_masks(z)
    names = ("lidc_full", "lidc_experts")
    for i in range(len(z["hw"])):
        ref = {k: z[k][i] for k in ("b0", "b1", "area", "perimeter", "corners", "chi8", "ring_sha")}
        d = C.mask_diff(out, i, ref)
        tally.add(names[z["src"][i]], d, i, lambda i=i: repro_mask(C.mask_of(z, i)))


# ------------------------------------------------------------------ speed
def timed(f, min_time=0.3, max_reps=50):
    """Median seconds per call over repetitions filling ~min_time (one warm-up call first)."""
    f()
    ts, t0 = [], time.perf_counter()
    while len(ts) < 3 or (time.perf_counter() - t0 < min_time and len(ts) < max_reps):
        t = time.perf_counter()
        f()
        ts.append(time.perf_counter() - t)
    return float(np.median(ts))


def speed_polygons(quick):
    rng = np.random.default_rng(0)
    rows = []
    sizes = [(16, 64), (128, 256), (1000, 512)]
    for nv, S in sizes:
        for B in ([1, 64] if quick else [1, 64, 1024]):
            p = C.pack_polygons([([C.blob(rng, nv, 500, 500, 300)], None, S, S) for _ in range(B)])
            rows.append({"vertices": nv, "crop": S, "batch": B, "gpucheck": B / timed(lambda: C.run_polygons(p))})
    return rows


def speed_masks(quick, lidc):
    rng = np.random.default_rng(0)
    rows = []
    sets = [("lidc", lidc[:256]), ("lidc", lidc[:8192])] if lidc else []
    sets += [("noise512", [rng.random((512, 512)) < .5 for _ in range(4)])]
    if quick:
        sets = sets[:1] + sets[-1:]
    for name, masks in sets:
        p = C.pack_masks(masks)
        def python_certify():
            for m in masks:
                bitquads(m), betti(m), polygonize(m)
        def opencv():
            for m in masks:
                u = m.astype(np.uint8)
                cv2.findContours(u, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE), cv2.connectedComponents(u, connectivity=8)
        B = len(masks)
        rows.append({"masks": name, "shape": list(masks[0].shape), "batch": B,
                     "gpucheck": B / timed(lambda: C.run_masks(p)),
                     "python_certify": B / timed(python_certify, max_reps=3),
                     "opencv_contours_ccl": B / timed(opencv)})
    return rows


# ------------------------------------------------------------------ report
def host_info():
    gpu = ""
    if shutil.which("nvidia-smi"):
        gpu = subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"],
                             capture_output=True, text=True).stdout.strip()
    commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"], capture_output=True,
                            text=True).stdout.strip()
    return {"build": _gpucheck.build(), "devices": _gpucheck.device_count(), "gpu": gpu,
            "machine": f"{platform.system()} {platform.machine()}", "cpu_count": os.cpu_count(),
            "python": platform.python_version(), "numpy": np.__version__, "opencv": cv2.__version__, "commit": commit,
            "date": time.strftime("%Y-%m-%d %H:%M")}


def markdown(res) -> str:
    h = res["host"]
    lines = [f"Run {h['date']}, commit {h['commit']}, build **{h['build']}**"
             + (f" on {h['gpu']}" if h["gpu"] else "") + f", {h['machine']}, {h['cpu_count']} CPUs; "
             f"wall time {res['wall_s']:.0f} s.", "",
             "| Case set | Cases | Mismatches | Exact match | First differing field (count) |", "|---|---|---|---|---|"]
    for name, s in res["correctness"].items():
        if s is None:
            lines.append(f"| {name} | absent | | | |")
            continue
        ff = ", ".join(f"{k} ({v})" for k, v in s["first_field"].items()) or "none"
        rate = f"{100 * s['exact_rate']:.3f} %" if s["exact_rate"] is not None else ""
        lines.append(f"| {name} | {s['n']:,} | {s['mismatches']} | {rate} | {ff} |")
    if res.get("per_class"):
        lines += ["", "LIDC polygons by reference DEFECTS class (a case counts once per class it has):", "",
                  "| Class | Cases | Mismatches |", "|---|---|---|"]
        lines += [f"| {c} | {n:,} | {b} |" for c, (n, b) in res["per_class"].items()]
    if res.get("speed_polygons"):
        lines += ["", "Polygons per second (defects + even-odd mask):", "",
                  "| Vertices | Crop | Batch | gpucheck |", "|---|---|---|---|"]
        lines += [f"| {r['vertices']} | {r['crop']} | {r['batch']} | {r['gpucheck']:,.0f} |"
                  for r in res["speed_polygons"]]
    if res.get("speed_masks"):
        lines += ["", "Masks per second (bit-quads + rings; python_certify = bitquad + betti + polygonize; OpenCV = "
                  "findContours + connectedComponents, a reference point, not the same output):", "",
                  "| Masks | Shape | Batch | gpucheck | Python certify path | OpenCV |", "|---|---|---|---|---|---|"]
        lines += [f"| {r['masks']} | {r['shape'][0]}x{r['shape'][1]} | {r['batch']} | {r['gpucheck']:,.0f} | "
                  f"{r['python_certify']:,.0f} | {r['opencv_contours_ccl']:,.0f} |" for r in res["speed_masks"]]
    return "\n".join(lines)


def write_section(md: str, build: str):
    doc = HERE / "RESULTS.md"
    text = doc.read_text()
    begin, end = f"<!-- BEGIN eval:{build} -->", f"<!-- END eval:{build} -->"
    block = f"{begin}\n{md}\n{end}"
    if begin in text:
        text = re.sub(re.escape(begin) + ".*?" + re.escape(end), lambda _: block, text, flags=re.S)
    else:
        text = text.rstrip() + "\n\n" + block + "\n"
    doc.write_text(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=str(ROOT / "data" / "checker_eval"))
    ap.add_argument("--quick", action="store_true", help="small fuzz set and speed grid (smoke run)")
    ap.add_argument("--no-speed", action="store_true")
    a = ap.parse_args()
    refuse_riga(a.data)
    t0 = time.time()
    data = Path(a.data)
    build = _gpucheck.build()
    if build == "cuda" and _gpucheck.device_count() == 0:
        raise SystemExit("CUDA build but no CUDA device")
    repro = data / "repro" / build
    res = {"host": host_info(), "correctness": {}}
    tp, tm = Tally("fuzz_polygons", repro), Tally("fuzz_masks", repro)
    eval_fuzz_polygons(tp, 20 if a.quick else 150, [0])  # the sets whose references were recorded
    eval_fuzz_masks(tm, 300 if a.quick else 2000, [0] if a.quick else [0, 1, 2])
    res["correctness"]["(a) fuzz polygons"] = tp.summary()
    res["correctness"]["(a) fuzz masks"] = tm.summary()
    lidc_masks = []
    if (data / "manifest.json").exists():
        man = check_data(data)
        res["manifest"] = {k: man[k] for k in ("masks", "polygons", "files")}
        p = dict(np.load(data / "polygons.npz", allow_pickle=False))
        t, per_class = Tally("lidc_polygons", repro), collections.defaultdict(lambda: [0, 0])
        eval_polygon_set(t, p, per_class)
        res["correctness"]["(b) LIDC polygons"] = t.summary()
        res["per_class"] = {c: per_class[c] for c in ("legal", *_gpucheck.DEFECTS) if c in per_class}
        z = dict(np.load(data / "masks.npz", allow_pickle=False))
        t = Tally("lidc_masks", repro)
        eval_mask_set(t, z)
        res["correctness"]["(c) LIDC Phase 0 masks"] = t.summary()
        sample = random.Random(0).sample(range(len(z["hw"])), min(8192, len(z["hw"])))
        lidc_masks = [C.mask_of(z, i) for i in sample]
    else:
        res["correctness"]["(b) LIDC polygons"] = res["correctness"]["(c) LIDC Phase 0 masks"] = None
    if not a.no_speed:
        res["speed_polygons"] = speed_polygons(a.quick)
        res["speed_masks"] = speed_masks(a.quick, lidc_masks)
    res["wall_s"] = time.time() - t0
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    (out / f"{build}{'-quick' if a.quick else ''}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    md = markdown(res)
    if not a.quick:
        write_section(md, build)
    print(md)
    bad = sum(s["mismatches"] for s in res["correctness"].values() if s)
    if bad:
        print(f"\n{bad} mismatches; smallest repro cases in {repro}")
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
