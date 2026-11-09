"""B4: aggregate shape statistics of the RIGA K3 prompt set (the first 64 RL train items), used only to calibrate
the synthetic RIGA-like generator (synth.py). Prints and writes quantiles only; no mask, image or item id leaves the
Mac. With --synthetic, the same statistics of synth.corpus (no RIGA data read).

    .venv/bin/python -I -m experiments.phase2.lab.B4.riga_stats --out <scratch>/phase2/B4 [--synthetic]
"""
import argparse

import cv2
import numpy as np

from experiments.phase1.reward_gates import make_jobs
from experiments.phase2.lab.A6.common import caches, dump, records
from objdel.checker.topology import instance_topology
from objdel.reward import ExpertCache

from .synth import corpus


def q(xs) -> list[float]:
    return [round(float(v), 3) for v in np.percentile(np.asarray(xs, np.float64), (5, 25, 50, 75, 95))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=64)
    ap.add_argument("--synthetic", action="store_true")
    a = ap.parse_args()
    cv2.setNumThreads(0)
    if a.synthetic:
        pairs = [{n: ExpertCache.build(it[n]) for n in ("disc", "cup")} for it in corpus(a.n)]
    else:
        pairs = [caches(r) for r in [r for r in records("riga", "train") if r["rl"]][:a.n]]
    out = shape_stats(pairs)
    print(out)
    dump(out, f"{a.out}/riga_stats{'_synthetic' if a.synthetic else ''}.json")


def shape_stats(pairs) -> dict:
    st = {k: [] for k in ("k", "frame_h", "frame_w", "disc_r", "cup_over_disc", "disc_cy", "disc_cx", "cup_dy",
                          "cup_dx", "band_disc", "band_cup", "disc_comps", "cup_comps", "disc_holes", "cup_holes",
                          "rater_r_sd_disc", "rater_r_sd_cup", "aspect_disc", "aspect_cup")}
    prompts = []
    for cs in pairs:
        prompts.append((("disc", "cup"), cs))
        d, c = cs["disc"], cs["cup"]
        H, W = d.union.shape
        st["k"].append(d.k)
        st["frame_h"].append(H)
        st["frame_w"].append(W)
        rd = [np.sqrt(np.count_nonzero(m) / np.pi) for m in d.raters]
        rc = [np.sqrt(np.count_nonzero(m) / np.pi) for m in c.raters]
        st["disc_r"].append(np.mean(rd))
        st["cup_over_disc"].append(np.mean(rc) / max(1e-9, np.mean(rd)))
        st["rater_r_sd_disc"].append(np.std(rd) / np.mean(rd))
        st["rater_r_sd_cup"].append(np.std(rc) / max(1e-9, np.mean(rc)))
        yd, xd = np.nonzero(d.vote)
        yc, xc = np.nonzero(c.vote) if c.n_vote else (yd, xd)
        st["disc_cy"].append(yd.mean() / H)
        st["disc_cx"].append(xd.mean() / W)
        st["cup_dy"].append((yc.mean() - yd.mean()) / np.mean(rd))
        st["cup_dx"].append((xc.mean() - xd.mean()) / np.mean(rd))
        st["band_disc"].append((np.count_nonzero(d.union) - np.count_nonzero(d.core)) / d.perimeter)
        st["band_cup"].append((np.count_nonzero(c.union) - np.count_nonzero(c.core)) / c.perimeter)
        for n, cache in (("disc", d), ("cup", c)):
            for m in cache.raters:
                t = instance_topology(m)
                st[f"{n}_comps"].append(len(t.areas))
                st[f"{n}_holes"].append(int(np.sum(t.holes)) if len(t.areas) else 0)
            _, _, w, h = cv2.boundingRect(cache.vote.view(np.uint8))
            st[f"aspect_{n}"].append(h / max(1, w))
    jobs = [j for b in range(4) for j in make_jobs(prompts[8 * b:8 * b + 32] if len(prompts) >= 32 else prompts, b)]
    out = {k: q(v) for k, v in st.items()}
    out["max_comps"] = {"disc": int(max(st["disc_comps"])), "cup": int(max(st["cup_comps"]))}
    out["max_holes"] = {"disc": int(max(st["disc_holes"])), "cup": int(max(st["cup_holes"]))}
    out["chars_per_completion"] = q([len(j[0]) for j in jobs])
    out["n_items"] = len(pairs)
    return out


if __name__ == "__main__":
    main()
