"""Iteration 3 (addendum): does clipping the DM-loo bracket at eps = 0.1 preserve the anti-collapse
behaviour? Uses ../t_dynamics/grpo_toy.py read-only (exact GRPO on tabular topology policies).

r = 1[valid] * max(eps, 1 + lam (log f(a) - [b(f)] - log pi_loo(a))), pi_loo = max(n-1, 0.5)/(G-1), G = 8;
b(f) = exact binomial bias of the LOO estimator at pi = f (the "_debias" variant, frozen in iteration 3).
Suite = t_dynamics SUITE + "hole 1/6" (6 raters, f = 1/6: the clip can bind at lam >= 0.5).
Settings: topology-only and IoU-gap 0.1 (sum, sigma 0.05), both mean-centred (Dr. GRPO, lr 1.0),
SFT calibrated or flattened (sharpen 0.5: over-produces minority answers). 20 seeds, 300 visits.
Reported (final, mean over seeds and prompts, 95 % CI over seeds): TV to the rater histogram,
invalid mass, min-mode retention, and the share of sampled valid rollouts whose bracket was clipped.

    ../.venv/bin/python scripts/dynamics_clip.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "t_dynamics"))

import numpy as np  # noqa: E402

import grpo_toy as gt  # noqa: E402

EPS = 0.1
CLIPPED = {"n": 0, "valid": 0}
_orig = gt.topo_reward


def _bias(fa, G):
    from scipy import stats as st

    n = np.arange(G)
    out = np.zeros_like(fa)
    for f in np.unique(fa[fa > 0]):
        pmf = st.binom.pmf(n, G - 1, f)
        out[fa == f] = np.log(f) - np.sum(pmf * np.log(np.maximum(n, 0.5) / (G - 1)))
    return out


def topo_reward(name, a, freq, logp_old, lam=0.5, floor=1e-3):
    if not name.startswith("set+dm_loo_clip"):
        return _orig(name, a, freq, logp_old, lam, floor)
    S, P, G = a.shape
    f = np.broadcast_to(freq, (S, P, gt.K))
    fa = np.take_along_axis(f, a, -1)
    valid = (fa > 0).astype(float)
    n, _ = gt.counts_in_group(a)
    bias = _bias(fa, G) if name.endswith("_debias") else 0.0
    raw = 1 + lam * (np.log(np.maximum(fa, floor)) - bias - np.log(np.maximum(n - 1, 0.5) / (G - 1)))
    CLIPPED["n"] += int(((raw < EPS) & (valid > 0)).sum())
    CLIPPED["valid"] += int(valid.sum())
    return valid * np.maximum(EPS, raw)


gt.topo_reward = topo_reward
SUITE = dict(gt.SUITE, **{"hole 1/6": {(1, 0): 5, (1, 1): 1}})


def iou_table(gap, base=0.85, inv_drop=0.02):
    mu = np.full((len(SUITE), gt.K), base - inv_drop)
    for p, n in enumerate(SUITE):
        f = gt.hist(SUITE[n])
        for k in np.where(f > 0)[0]:
            mu[p, k] = base - (gap if gt.ANSWERS[k] != (1, 0) else 0)
    return mu


def final(o):
    out = {}
    for key in ("tv", "invalid", "minret"):
        v = o[key][-1]  # (S, P)
        v = np.where(np.isfinite(v), v, np.nan)
        per_seed = np.nanmean(v, 1)
        m = per_seed.mean()
        h = 1.96 * per_seed.std(ddof=1) / np.sqrt(len(per_seed))
        out[key] = f"{m:.3f}±{h:.3f}"
    return out


def main():
    print(f"{'setting':44s} {'variant':18s} {'TV':>12s} {'invalid':>12s} {'min-mode ret':>13s} {'clipped share':>13s}")
    for setting, kw, sharpen in (("topology only, calibrated SFT", dict(combine="none"), 1.0),
                                 ("topology only, flattened SFT (sharpen 0.5)", dict(combine="none"), 0.5),
                                 ("IoU gap 0.1, calibrated SFT", dict(combine="sum", iou_mu=iou_table(0.1), iou_sd=0.05), 1.0)):
        for lam in (0.5, 1.0):
            for variant in ("set+dm_loo", "set+dm_loo_clip", "set+dm_loo_clip_debias"):
                CLIPPED.update(n=0, valid=0)
                cfg = gt.Cfg(reward=variant, lam=lam, adv="mean", lr=1.0, **kw)
                o = gt.run(cfg, suite=SUITE, S=20, T=300, seed=3, sharpen=sharpen)
                r = final(o)
                share = CLIPPED["n"] / max(1, CLIPPED["valid"]) if "clip" in variant else float("nan")
                print(f"{setting:44s} {variant + f' λ={lam}':30s} {r['tv']:>12s} {r['invalid']:>12s} {r['minret']:>13s} {share:13.4f}")


if __name__ == "__main__":
    main()
