"""Honest GRPO rollouts for the adversary and neutrality tests (iteration 4).

A pool entry perturbs EVERY rater of the item with the SAME amplitude, so all topology classes share
one perturbation-quality distribution. Iteration 4 found that per-rater rejection sampling (retry
until the size-filtered topology survives, else fall back to the exact mask) fell back 16.5 % of
the time for minority-topology raters (holes and small parts are fragile) vs 2.8 % for the mode
class, which gave minority answers systematically better honest copies.

modes
  jit1px : amplitude 1 px, correlation length 3 px (iteration 3 construction)
  real   : with probability 0.25 every rater's exact mask; otherwise amplitude ~ U(0, tau) and
           correlation length max(3, tau) px, tau = the item's P95 rater-pair contour distance
For an amplitude, each rater gets up to 20 noise draws to keep its topology; if any rater fails, a
new amplitude is drawn for all raters (up to 30 times, then all exact).
"""
from __future__ import annotations

from rewardcore import perturb
from rewardcore.topology import filtered_signature


def honest_pools(stack, sigs, cache, rng, mode: str, n: int, pack, score) -> list[list]:
    """pools[j] = list of n packed Terms for rater j (pool entry t is matched across raters)."""
    k = len(stack)
    exact = [pack(score(m, cache)) for m in stack]
    pools = [[] for _ in range(k)]
    for _ in range(n):
        if mode == "real" and rng.random() < 0.25:
            for j in range(k):
                pools[j].append(exact[j])
            continue
        entry = None
        for _amp_try in range(30):
            amp = 1.0 if mode == "jit1px" else float(rng.uniform(0, cache.tau))
            corr = 3.0 if mode == "jit1px" else max(3.0, cache.tau)
            out = []
            for j, m in enumerate(stack):
                for _t in range(20):
                    x = perturb.jitter(m, amp, rng, corr=corr)
                    if filtered_signature(x, cache.min_area) == sigs[j]:
                        out.append(x)
                        break
                else:
                    break
            if len(out) == k:
                entry = [pack(score(x, cache)) for x in out]
                break
        for j in range(k):
            pools[j].append(entry[j] if entry is not None else exact[j])
    return pools
