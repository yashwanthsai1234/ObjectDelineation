"""Phase 2 pre-registration: every pre-registered decision rule as data plus pure functions (PLAN §5.0-§5.2, §6.3,
D12; M8). Authoritative where it and the plan's prose differ (§5.0). Fixed before V1; a rule is never changed after
a run sees data (a later re-specification gets a new name and an AMEND note here and in PLAN §9).

    .venv/bin/python -I -m experiments.phase2.prereg   # prints the rules as JSON

Rules:
- C1 stop rule (§5.2): rungs smallest first; C1 = the first rung passing every C1 gate. Report-only: G-prior (a)
  (unless the owner sets G_PRIOR_A_THRESHOLD before V1, never after seeing C1), G-topology (b); G-probe never climbs
  the ladder (not PASS after the D12(c) repair: the rung passes with a "G-probe warning").
- D12 rule R-B: `objdel.eval.signal.select_rb` (veto statistic sd_geom, post-hoc; r_std reported, disagreement
  flagged), constants in RB.
- G-probe (§5.1): PASS / FAIL / INCONCLUSIVE against the step-0 null (two step-0 draws, different sampler seeds).
- C2 (§5.2, §5.0): clauses (1)-(5), the 3-seed rule (>= 2 of 3 seeds >= 0.50, none < 0.30) and the fourth-seed
  clause (>= 3 of 4), outcome pass / partial / fail with at most one conditional run (owner approval required).
- §6.3 fallback: the component removed by the single bundle-level diagnostic.
"""

import argparse
import json

import numpy as np

from objdel.eval import gates, signal

# ---------------------------------------------------------------- C1 (§5.1, §5.2)
G_PRIOR_A_THRESHOLD = None  # owner-set before V1 only (<= SAM 2.1-tiny's centred share - 5 pts); None = report-only
C1_GATES = ("G-format", "G-geometry", "G-prior-b", "G-precision", "G-topology-a", "G-signal", "G-vision")
REPORT_ONLY = ("G-prior-a", "G-topology-b", "G-probe")


def c1_gates(prior_a_threshold: float | None = G_PRIOR_A_THRESHOLD) -> tuple[str, ...]:
    return C1_GATES + (("G-prior-a",) if prior_a_threshold is not None else ())


def entry(report: dict, datasets=None) -> dict:
    """gates.phase2_entry under the pre-registered G-prior (a) threshold (read at call time): a rung's C1 decision,
    whose "passed" is c1_stop's input (BK2-12: a set threshold then has its "G-prior-a" gate)."""
    return gates.phase2_entry(report, datasets, G_PRIOR_A_THRESHOLD)


def c1_stop(rungs, prior_a_threshold: float | None = G_PRIOR_A_THRESHOLD) -> tuple[str | None, list[str]]:
    """rungs = [(name, {gate: passed})] smallest first -> (C1 rung or None, its warnings). Stops at the first pass."""
    for name, passed in rungs:
        if all(passed.get(g, False) for g in c1_gates(prior_a_threshold)):
            return name, (["G-probe warning"] if passed.get("G-probe") is False else [])
    return None, []


# ---------------------------------------------------------------- D12 R-B (signal.select_rb)
RB = {
    "veto_stat": signal.VETO_STAT,
    "reported_stat": signal.ALT_STAT,
    "veto_rel": signal.VETO_REL,
    "boot": signal.BOOT,
    "level": signal.LEVEL,
    "eligible": {"parse": 0.98, "legal": 0.90, "dup_max": 0.10, "sigma_share": 0.50, "g_vision": True},
}
select_rb = signal.select_rb

# ---------------------------------------------------------------- G-probe (§5.1)
G_PROBE = {
    "steps": 20,
    "lr_full": 5e-6,
    "lr_lora": 5e-5,
    "items": 128,
    "K": 8,
    "legal_drop_max": 0.05,
    "sigma_share_min": 0.50,
}


def g_probe_verdict(gain_ci, null_ci, legal_drop: float, sigma_share: float) -> str:
    """gain_ci / null_ci = (lo, hi) of the paired per-item gain and of the step-0 null; legal_drop in fractions."""
    if gain_ci[1] < null_ci[0] or legal_drop > G_PROBE["legal_drop_max"] or sigma_share < G_PROBE["sigma_share_min"]:
        return "FAIL"
    return "PASS" if gain_ci[0] > 0 else "INCONCLUSIVE"


def g_probe(r20, r0a, r0b, ids, legal0: float, legal20: float, sigma20: float, boot=signal.BOOT, seed=signal.SEED):
    """Per-item raw composite R at step 20 and the two step-0 draws -> null CI, h0, gain CI and the verdict."""
    r20, r0a, r0b = (np.asarray(x, float) for x in (r20, r0a, r0b))
    D = signal.draws(ids, boot, seed)
    null, gain = signal.ci(r0a - r0b, D), signal.ci(r20 - (r0a + r0b) / 2, D)
    verdict = g_probe_verdict(gain[1:], null[1:], legal0 - legal20, sigma20)
    return {
        "null": null,
        "h0": (null[2] - null[1]) / 2,
        "gain": gain,
        "legal_drop": legal0 - legal20,
        "sigma_share": sigma20,
        "verdict": verdict,
    }


# ---------------------------------------------------------------- C2 (§5.2, D20)
C2 = {
    "pass_line": 0.50,
    "floor": 0.30,
    "false_nonsimple_max": 0.05,
    "legal_min": 0.90,
    "seeds": 3,
    "seeds_needed": {3: 2, 4: 3},
    "bin_weights": (0.4386, 0.3614, 0.2000),
    "per_bin": 120,
}


def c2_seed_clause(seed_points) -> bool:
    """Clause (2): >= 2 of 3 (or, with the fourth seed, >= 3 of 4) seed estimates >= 0.50 and none < 0.30."""
    x = np.asarray(seed_points, float)
    if len(x) not in C2["seeds_needed"]:
        raise ValueError("C2 is judged on 3 seeds, or 4 with the fourth-seed clause")
    return bool((x >= C2["pass_line"]).sum() >= C2["seeds_needed"][len(x)] and (x >= C2["floor"]).all())


def c2_clauses(pooled, lo, seed_points, bin_points, false_nonsimple, legal) -> dict[str, bool]:
    return {
        "1": bool(pooled >= C2["pass_line"] and lo > C2["floor"]),
        "2": c2_seed_clause(seed_points),
        "3": bool(min(bin_points) >= C2["floor"]),
        "4": bool(false_nonsimple <= C2["false_nonsimple_max"]),
        "5": bool(legal >= C2["legal_min"]),
    }


def c2_outcome(pooled: float, clauses: dict[str, bool]) -> tuple[str, str | None]:
    """(outcome, conditional run): pass -> none; fail (pooled < 0.30) -> R-ctrl; partial -> one x8 A0' seed (R-x8),
    not run when clause (4) failed."""
    if pooled < C2["floor"]:
        return "fail", "R-ctrl"
    if all(clauses.values()):
        return "pass", None
    return "partial", ("R-x8" if clauses["4"] else None)


def fourth_seed_eligible(seed_points, pooled: float, other: dict[str, bool]) -> bool:
    """The 3-seed rule is decided by a single seed's side of 0.50: one seed crossing 0.50 flips clause (2) (1 or 2 of
    3 at >= 0.50, none < 0.30) and the flip changes the outcome. Adding the seed also needs owner approval."""
    x = np.asarray(seed_points, float)
    if len(x) != 3 or (x < C2["floor"]).any() or (x >= C2["pass_line"]).sum() not in (1, 2):
        return False
    return c2_outcome(pooled, other | {"2": True}) != c2_outcome(pooled, other | {"2": False})


# ---------------------------------------------------------------- §6.3 bundle-level fallback
FALLBACK_TRIGGERS = ("G-format", "G-geometry", "G-prior-b", "G-precision", "G-topology-a")
HOLE_DEFECT_SHARE = 0.50


def fallback_component(failed, hole_defect_share: float, legal: float, used: bool = False) -> str | None:
    """Component removed by the single §6.3 diagnostic, or None (no trigger, or the diagnostic was already used).
    failed = the selected checkpoint's failed gates (G-precision after its capacity run); hole_defect_share = share
    of illegal samples with a hole-related defect; legal = legal share at T = 1."""
    failed = set(failed) & set(FALLBACK_TRIGGERS)
    if used or not failed:
        return None
    if hole_defect_share >= HOLE_DEFECT_SHARE or "G-topology-a" in failed:
        return "hole_oversampling"
    if failed & {"G-geometry", "G-prior-b", "G-precision"} and legal >= C2["legal_min"]:
        return "dihedral"
    return "window_jitter"


def rules() -> dict:
    return {
        "C1": {"gates": c1_gates(), "report_only": REPORT_ONLY, "G_PRIOR_A_THRESHOLD": G_PRIOR_A_THRESHOLD},
        "RB": RB,
        "G_PROBE": G_PROBE,
        "C2": C2,
        "fallback": {
            "triggers": FALLBACK_TRIGGERS,
            "hole_defect_share": HOLE_DEFECT_SHARE,
            "order": ("hole_oversampling", "dihedral", "window_jitter"),
        },
    }


def main():
    argparse.ArgumentParser(description="Print Phase 2's pre-registered rules as JSON.").parse_args()
    print(json.dumps(rules(), indent=1, default=str))


if __name__ == "__main__":
    main()
