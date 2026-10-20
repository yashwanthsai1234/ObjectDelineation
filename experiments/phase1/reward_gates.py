"""Reward v4 gates (plan §7: K3, K4, K8, K8b, K9, K10, R3, R4, R5) and the M6 studies (per-rollout local
gradient, leave-one-rater-out g_area = AB4, hole-matching spike = AB3, Phase 0 / v2 / v4 = AB2) on the fixtures
and the full-LIDC train disagreement items (239; 251 before sub-filter raters abstained, B2-4). Writes the result
tables into experiments/phase1/impl/M6_reward.md (between the `gates` markers) and prints them.

    .venv/bin/python -I -m experiments.phase1.reward_gates [--quick] [--workers 3]

Data: data/fixtures (untrusted, hence -I). The harness (honest rollouts, adversaries, K9 pairs) and the A/B
variants live in tests/oracles/reward_harness.py and reward_variants.py.

References:
- Gate definitions: experiments/phase1/PLAN.md §7 (K3, K4, K8, K8b, K9, K10, R3, R4, R5); reward v4:
  objdel.reward (its module references); lab protocol: experiments/phase1/lab/t_reward/results/ITERATION4.md.
- Spearman rho (K8): https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html
- Cluster (item-level) percentile bootstrap CIs (K8, K9, R3, R4): C. A. Field, A. H. Welsh, "Bootstrapping clustered
  data", JRSS B 69(3), 2007, https://doi.org/10.1111/j.1467-9868.2007.00593.x
- Timing (K3, K4): https://docs.python.org/3/library/time.html#time.perf_counter_ns
- Determinism in a fresh interpreter (K10): https://docs.python.org/3/using/cmdline.html#cmdoption-I
"""
import argparse
import hashlib
import os
import pickle
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from multiprocessing import Pool
from pathlib import Path

import cv2
import numpy as np
from scipy import stats

from objdel.geometry import canonicalize, serialize
from objdel.reward import NUMERIC, ExpertCache, advantages_v4, score_batch, score_completion
from tests.oracles import reward_harness as H
from tests.oracles import reward_variants as V

ROOT = Path(__file__).resolve().parents[2]
NOTE = ROOT / "experiments" / "phase1" / "impl" / "M6_reward.md"
BEGIN, END = "<!-- gates:begin -->", "<!-- gates:end -->"
OUT: list[str] = []


def emit(line: str = "") -> None:
    print(line, flush=True)
    OUT.append(line)


def pct(x) -> str:
    return f"{100 * x:.1f} %"


def ci_str(lo, hi) -> str:
    return f"[{pct(lo)}, {pct(hi)}]"


def rho_ci(v, a, items, rng, n_boot=1000):
    """Spearman rho and a cluster-bootstrap (items) 95 % CI."""
    items = np.array(items)
    rows = [np.flatnonzero(items == k) for k in sorted(set(items))]
    boots = []
    for _ in range(n_boot):
        sel = np.concatenate([rows[i] for i in rng.integers(0, len(rows), len(rows))])
        boots.append(stats.spearmanr(v[sel], a[sel]).statistic)
    return stats.spearmanr(v, a).statistic, *np.percentile(boots, [2.5, 97.5])


# ------------------------------------------------------------------ band: K8, K8b, R5

def band_gates(pool, items, rng):
    emit("### K8, K8b, R5 (band)\n")
    eligible = [(ds, k, st) for ds, k, st in items if len(st) >= 3]
    rows = [r for rs in pool.starmap(H.k8_rows, eligible, chunksize=8) for r in rs]
    emit("| stratum (leave one rater out) | n | Spearman ρ(F5, area) [95 % CI] | median F5 | K8 (abs ρ ≤ 0.10) |")
    emit("|---|---|---|---|---|")
    for name, sub, k in (("riga/disc", "riga/disc", None), ("riga/cup", "riga/cup", None),
                         ("lidc k=3", "lidc/nodule", 3), ("lidc k=4", "lidc/nodule", 4)):
        g = [r for r in rows if r["subset"] == sub and k in (None, r["k"])]
        if not g:
            continue
        v, a = np.array([r["f5"] for r in g]), np.array([r["area"] for r in g])
        rho, lo, hi = rho_ci(v, a, [r["key"] for r in g], rng)
        emit(f"| {name} | {len(g)} | {rho:+.3f} [{lo:+.3f}, {hi:+.3f}] | {np.median(v):.3f} | "
             f"{'PASS' if abs(rho) <= 0.10 else 'FAIL'} |")
    v, a = np.array([r["f5"] for r in rows]), np.array([r["area"] for r in rows])
    emit(f"\nPooled ρ (reported, not gated): {stats.spearmanr(v, a).statistic:+.3f}. LIDC k = 2 items have no "
         "leave-one-out case (one remaining rater is not RL-eligible; production has no single-rater prior).\n")
    k8b = dict(zip([k for _, k, _ in items], pool.map(H.k8b_item, [st for _, _, st in items], chunksize=8), strict=True))
    held = {}
    for r in rows:
        held.setdefault(r["key"], []).append(r["f5"])
    mono = np.array([d["dil_mono"] and d["ero_mono"] for d in k8b.values()])
    keys = sorted(held)
    med = np.array([np.median(held[k]) for k in keys])
    dil = np.array([k8b[k]["dil"] for k in keys])
    ero = np.array([k8b[k]["ero"] for k in keys])
    m50, p10 = np.median(v), np.percentile(v, 10)
    emit("| K8b criterion | result | target | verdict |")
    emit("|---|---|---|---|")
    emit(f"| (i) F5 monotone in excursion (0.25–2 τ, dilate and erode), all {len(mono)} items | {pct(mono.mean())} | ≥ 98 % | "
         f"{'PASS' if mono.mean() >= 0.98 else 'FAIL'} |")
    emit(f"| (ii) 2τ-dilated vote below the median held-out rater ({len(keys)} items, k ≥ 3) | {pct(np.mean(dil < med))} "
         f"(eroded: {pct(np.mean(ero < med))}) | ≥ 95 % | {'PASS' if np.mean(dil < med) >= 0.95 else 'FAIL'} |")
    emit(f"| (iii) held-out F5 median / p10 ({len(v)} cases) | {m50:.3f} / {p10:.3f} | ≥ 0.90 / ≥ 0.70 | "
         f"{'PASS' if m50 >= 0.9 and p10 >= 0.7 else 'FAIL'} |")
    riga = [(st, j) for j, (ds, _, st) in enumerate(items) if ds == "riga"]
    sub = pool.starmap(H.r5_item, riga, chunksize=8)
    emit("\n| R5: m reference raters (RIGA) | median F5 of a held-out expert | p10 | verdict ([0.92, 0.98]) |")
    emit("|---|---|---|---|")
    for m in (2, 3, 4, 5):
        x = np.concatenate([s[m] for s in sub])
        emit(f"| {m} | {np.median(x):.3f} | {np.percentile(x, 10):.3f} | {'PASS' if 0.92 <= np.median(x) <= 0.98 else 'FAIL'} |")
    emit("\nSingle-rater items: `ExpertCache.build` refuses k = 1 (`rl_eligible`); tested.\n")


# ------------------------------------------------------------------ K9

def _k9_item(args):
    key, stack, two_part, seed = args
    return [(g, hi, key) for g, hi in H.k9_groups(stack, np.random.default_rng(seed), two_part, pack_fn=V.pack_ab)]


def k9_gate(pool, item_args, label, rng):
    groups = [g for gs in pool.map(_k9_item, item_args, chunksize=4) for g in gs]
    T = np.array([g for g, _, _ in groups])
    hi = np.array([h for _, h, _ in groups])
    keys = [k for _, _, k in groups]
    order = rng.permutation(len(T))
    res = []
    for name in V.SCHEMES:
        a = np.zeros((len(T), H.G))
        a[order] = H.in_steps(V.SCHEMES[name], T[order])
        ok = a[:, 0] > a[:, 1]
        lo, up = H.cluster_ci(ok, keys, rng)
        hlo, hup = H.cluster_ci(ok[hi], [k for k, h in zip(keys, hi, strict=True) if h], rng)
        res.append(f"{pct(ok.mean())} {ci_str(lo, up)}; IoU ≥ 0.5: {pct(ok[hi].mean())}")
    emit(f"| {label} | {len(T)} pairs, {len(set(keys))} items | " + " | ".join(res) + " |")


# ------------------------------------------------------------------ R3 adversaries, R4 neutrality

INFLATION = ("dilated_union", "band_speck", "hull", "union", "core")


def adversary_table(data_fix, data_all, rng, mode):
    emit(f"\nHonest rollouts: **{mode}**. Cell = share of groups where the adversary's advantage is below the median "
         "honest advantage [95 % cluster-bootstrap CI]; hole rows: win rate against the median honest rollout of the "
         "same topology class (re-registered KPI ≤ 60 %, CI upper ≤ 65 %), then the original below-median share.\n")
    emit("| adversary | items | " + " | ".join(V.SCHEMES) + " |")
    emit("|---|---|" + "---|" * len(V.SCHEMES))
    for name in INFLATION + H.HOLES:
        data = data_fix if name in INFLATION else data_all
        T, keys, same = H.adversary_groups(data, name, mode, rng)
        if not len(T):
            continue
        cells = []
        for adv in V.SCHEMES.values():
            a = H.in_steps(adv, T)
            below = a[:, 0] < np.median(a[:, 1:], axis=1)
            lo, hi = H.cluster_ci(below, keys, rng)
            if name in H.HOLES:
                win = np.array([a[i, 0] > np.median(a[i, idx]) for i, idx in enumerate(same)])
                wlo, whi = H.cluster_ci(win, keys, rng)
                cells.append(f"win {pct(win.mean())} {ci_str(wlo, whi)}; below {pct(below.mean())}")
            else:
                cells.append(f"{pct(below.mean())} {ci_str(lo, hi)}")
        emit(f"| {name} | {len(set(keys))} | " + " | ".join(cells) + " |")


def neutrality_table(data, rng):
    emit(f"\n| R4 neutrality (minority − mode advantage), {len(data)} items, 300 batches | " + " | ".join(V.SCHEMES)
         + " |")
    emit("|---|" + "---|" * len(V.SCHEMES))
    for mode in ("exact", "jit1px", "real"):
        cells = []
        for adv in V.SCHEMES.values():
            m, lo, hi = H.mean_ci(H.neutrality(data, mode, rng, adv), rng)
            ok = lo <= 0 <= hi and (hi - lo) / 2 <= 0.10
            cells.append(f"{m:+.3f} [{lo:+.3f}, {hi:+.3f}] {'ACCEPT' if ok else 'reject'}")
        emit(f"| {mode} | " + " | ".join(cells) + " |")


# ------------------------------------------------------------------ local gradient, LOO g_area

def local_gradient(pool, args, label, rng):
    res = pool.map(H.local_gradient_item, args, chunksize=4)
    emit(f"\n| local gradient, {label} | moved rollouts | ΔA < 0 (per rollout) | mean ΔA [95 % CI over items] |")
    emit("|---|---|---|---|")
    for name, _, _ in H.MOVES:
        d = np.concatenate([r["deltas"][name] for r in res if r["deltas"][name]])
        per_item = [np.mean(r["deltas"][name]) for r in res if r["deltas"][name]]
        m, lo, hi = H.mean_ci(per_item, rng)
        emit(f"| {name} | {len(d)} | {pct((d < 0).mean())} | {m:+.3f} [{lo:+.3f}, {hi:+.3f}] |")
    allv = np.concatenate([np.concatenate([r["deltas"][n] for n, _, _ in H.MOVES]) for r in res])
    away = np.concatenate([np.concatenate([np.array(r["away"][n], bool) for n, _, _ in H.MOVES]) for r in res])
    emit(f"| all moves | {len(allv)} | **{pct((allv < 0).mean())}** (target ≥ 95 %) | |")
    emit(f"| moves away from the median rater area (dilating a rollout ≥ it, eroding one ≤ it) | {away.sum()} | "
         f"{pct((allv[away] < 0).mean())} | |")
    emit(f"| moves toward the median rater area | {(~away).sum()} | {pct((allv[~away] < 0).mean())} | |")


def loo_g_area(pool, stacks, label):
    pairs = np.array([p for ps in pool.map(H.g_area_loo, stacks, chunksize=8) for p in ps])
    z, t = pairs[:, 0], pairs[:, 1]
    emit(f"| {label} | {len(pairs)} | {np.median(z):.3f} / {np.percentile(z, 10):.3f} / {np.mean(z < 1):.1%} | "
         f"{np.median(t):.3f} / {np.percentile(t, 10):.3f} / {np.mean(t < 1):.1%} |")
    return np.median(z)


# ------------------------------------------------------------------ K3 / K4 / K10 (completion text)

def to_parts(m, W, Hh, eps=1.0):
    """Mask -> polygon parts on the 0..1000 grid (findContours + approxPolyDP 1 px): realistic completions."""
    cs, hier = cv2.findContours(m.astype(np.uint8), cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
    if hier is None:
        return []

    def ring(c):
        a = cv2.approxPolyDP(c, eps, True).reshape(-1, 2).astype(np.float64)
        g = np.stack([np.round((a[:, 0] + 0.5) * 1000 / W), np.round((a[:, 1] + 0.5) * 1000 / Hh)], 1)
        return np.clip(g, 0, 1000).astype(np.int64) if len(g) >= 3 else None

    parts, idx = [], {}
    for i, (c, h) in enumerate(zip(cs, hier[0], strict=True)):
        if h[3] < 0 and (r := ring(c)) is not None:
            idx[i] = len(parts)
            parts.append((r, []))
    for c, h in zip(cs, hier[0], strict=True):
        if h[3] in idx and (r := ring(c)) is not None:
            parts[idx[h[3]]][1].append(r)
    return parts


def make_jobs(prompts, seed):
    """32 prompts x 8 completions; RIGA prompts carry disc + cup of one rater, perturbed."""
    rng = np.random.default_rng(seed)
    jobs = []
    for names, caches in prompts:
        Hh, W = caches[names[0]].union.shape
        for _ in range(H.G):
            r = int(rng.integers(caches[names[0]].k))
            kind = int(rng.integers(0, 4))
            objs = {}
            for n in names:
                m = caches[n].raters[r]
                m = (H.jitter(m, float(rng.uniform(0.5, 3)), rng, 6.0) if kind == 0 else
                     H.shift(m, *map(int, rng.integers(-3, 4, 2))) if kind == 1 else
                     H.grow(m, float(rng.uniform(-2, 2))) if kind == 2 else m)
                objs[n] = to_parts(m, W, Hh)
            jobs.append((serialize(canonicalize(objs), names), False, names, caches, W, Hh))
    return jobs


def throughput(items, rng, n_batches):  # items: all fixture items (RIGA disc + cup pairs needed)
    cv2.setNumThreads(0)  # OpenCV runs sequentially (on macOS's GCD backend only values <= 0 take effect)
    riga = {}
    for ds, k, st in items:
        if ds == "riga":
            riga.setdefault(k.rsplit("__", 1)[0], {})[k.rsplit("__", 1)[1]] = st
    prompts_r = [(("disc", "cup"), {n: ExpertCache.build(v[n]) for n in ("disc", "cup")}) for v in list(riga.values())[:64]]
    prompts_l = [(("nodule",), {"nodule": ExpertCache.build(st)}) for ds, _, st in items if ds == "lidc"][:200]
    emit("\n| K3 / K4 workload (256 completions = 32 prompts × 8) | executor | batch wall median / max (ms) | "
         "per-completion p50 / p99 (ms) |")
    emit("|---|---|---|---|")
    for label, share_riga in (("mixed (40 % RIGA disc + cup, 60 % LIDC)", 0.4), ("RIGA only", 1.0)):
        batches = []
        for b in range(n_batches):
            n_r = int(round(32 * share_riga))
            pr = [prompts_r[i] for i in rng.choice(len(prompts_r), n_r, replace=False)]
            pl = [prompts_l[i] for i in rng.choice(len(prompts_l), 32 - n_r, replace=False)]
            batches.append(make_jobs(pr + pl, 1000 + b))
        lat = []
        for jobs in batches:
            for j in jobs:
                t0 = time.perf_counter_ns()
                score_completion(*j)
                lat.append((time.perf_counter_ns() - t0) / 1e6)
        for ex_name, ex in (("serial", None), ("8 threads", ThreadPoolExecutor(8))):
            walls = []
            for jobs in batches:
                for rep in range(2):
                    t0 = time.perf_counter()
                    advantages_v4(score_batch(jobs, ex))
                    if rep:
                        walls.append(1e3 * (time.perf_counter() - t0))
            if ex is not None:
                ex.shutdown()
            emit(f"| {label} | {ex_name} | {np.median(walls):.1f} / {np.max(walls):.1f} | "
                 f"{np.median(lat):.2f} / {np.percentile(lat, 99):.2f} |")
    emit(f"\nLoad average during the run: {os.getloadavg()[0]:.1f} (8-core Apple M3, shared with other agents).")
    return batches[0]


def digest(t) -> str:
    """SHA-256 of a batch's numeric terms and its advantages (K10)."""
    h = hashlib.sha256()
    for k in NUMERIC:
        h.update(np.ascontiguousarray(t[k], np.float64).tobytes())
    h.update(np.ascontiguousarray(advantages_v4(t), np.float64).tobytes())
    return h.hexdigest()


K10_CODE = """import pickle, sys
from experiments.phase1.reward_gates import digest
from objdel.reward import score_batch
print(digest(score_batch(pickle.load(open(sys.argv[1], "rb")))))
"""


def k10(jobs, tmp: Path):
    out = {"serial run 1": digest(score_batch(jobs)), "serial run 2": digest(score_batch(jobs))}
    with ThreadPoolExecutor(8) as ex:
        out["8 threads"] = digest(score_batch(jobs, ex))
    with ProcessPoolExecutor(3) as ex:
        out["3 processes"] = digest(score_batch(jobs, ex))
    cv2.setNumThreads(8)
    out["OpenCV 8 threads"] = digest(score_batch(jobs))
    cv2.setNumThreads(0)
    tmp.write_bytes(pickle.dumps(jobs))
    env = dict(os.environ, PYTHONHASHSEED="123")
    out["fresh interpreter"] = subprocess.run([sys.executable, "-I", "-c", K10_CODE, str(tmp)], capture_output=True,
                                              text=True, check=True, env=env).stdout.strip()
    same = len(set(out.values())) == 1
    emit(f"\nK10: SHA-256 of 256 completions' terms + advantages identical across {', '.join(out)}: "
         f"**{'PASS' if same else 'FAIL'}** ({next(iter(out.values()))[:16]}…).")


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixtures", type=Path, default=ROOT / "data" / "fixtures")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--quick", action="store_true", help="every 3rd fixture item (smoke run)")
    ap.add_argument("--no-note", action="store_true", help="print only; do not rewrite the impl note")
    args = ap.parse_args()
    rng = np.random.default_rng(0)
    t_start = time.time()
    items = H.load_fixture_items(args.fixtures)
    lidc = H.load_lidc_disagreement(args.fixtures)
    if args.quick:
        items = items[::3]
    emit(f"Fixture items with ≥ 2 raters: {len(items)}; full-LIDC train disagreement items: {len(lidc)}.\n")
    with Pool(args.workers) as pool:
        band_gates(pool, items, rng)

        emit("### K9 (topology ranking; v4 gated ≥ 95 %)\n")
        emit("| item set | pairs | " + " | ".join(V.SCHEMES) + " |")
        emit("|---|---|" + "---|" * len(V.SCHEMES))
        fix_args = [(k, st, False, j) for j, (_, k, st) in enumerate(items)]
        for sub in ("lidc/nodule", "riga/cup"):
            grp = [(k, st) for ds, k, st in items if H.subset(ds, k) == sub]
            fix_args += [(f"{ka}+{kb}", H.composite(a, b, 3 + j % 6), True, 10_000 + j)
                         for j, ((ka, a), (kb, b)) in enumerate(zip(grp[0::2], grp[1::2], strict=False))]
        k9_gate(pool, fix_args, "fixtures + 2-part composites", rng)
        k9_gate(pool, [(k, st, False, j) for j, (k, st) in enumerate(lidc.items())], "full-LIDC train disagreement",
                rng)

        emit("\n### R3 adversaries (dilated union and in-band speck gated ≥ 95 %; others reported) and AB2/AB3/AB4\n")
        fix_data = pool.map(H.item_pools, [(k, st, j, 12, V.pack_ab, True) for j, (_, k, st) in enumerate(items)],
                            chunksize=4)
        lidc_data = pool.map(H.item_pools, [(k, st, 50_000 + j, 12, V.pack_ab, True) for j, (k, st) in
                                            enumerate(lidc.items())], chunksize=4)
        all_data = fix_data + lidc_data
        for mode in ("real", "jit1px"):
            adversary_table(fix_data, all_data, rng, mode)
        neutrality_table(lidc_data, rng)

        emit("\n### Local gradient (dilating or eroding one honest rollout, others fixed)")
        local_gradient(pool, [(k, st, 70_000 + j, 20) for j, (k, st) in enumerate(lidc.items())],
                       f"{len(lidc)} full-LIDC items × 20 groups", rng)
        local_gradient(pool, [(k, st, 80_000 + j, 6) for j, (_, k, st) in enumerate(items)],
                       "fixture items × 6 groups", rng)

        emit("\n### AB4: leave-one-rater-out g_area (held-out rater vs the other raters' area range)\n")
        emit("| item set | held-out cases | zero tolerance: median / p10 / share < 1 | τ·P tolerance: median / p10 / share < 1 |")
        emit("|---|---|---|---|")
        meds = [loo_g_area(pool, [st for ds, _, st in items if ds == "riga"], "RIGA fixtures"),
                loo_g_area(pool, [st for ds, _, st in items if ds == "lidc"], "LIDC fixtures"),
                loo_g_area(pool, list(lidc.values()), "full-LIDC disagreement")]
        emit(f"\nSwitch rule (median zero-tolerance g_area < 0.9 → τ·P): medians {', '.join(f'{m:.3f}' for m in meds)} → "
             f"**{'switch to τ·P' if min(meds) < 0.9 else 'keep zero tolerance'}**.")

    emit("\n### K3, K4, K10 (completion text → reward → advantages)")
    jobs = throughput(H.load_fixture_items(args.fixtures), rng, 4 if args.quick else 8)
    with tempfile.TemporaryDirectory() as tmp:
        k10(jobs, Path(tmp) / "jobs.pkl")
    emit(f"\nTotal run time {time.time() - t_start:.0f} s, {args.workers} worker processes.")
    if not args.no_note:
        text = NOTE.read_text() if NOTE.exists() else f"{BEGIN}\n{END}\n"
        a, b = text.index(BEGIN) + len(BEGIN), text.index(END)
        NOTE.write_text(text[:a] + "\n" + "\n".join(OUT) + "\n" + text[b:])


if __name__ == "__main__":
    main()
