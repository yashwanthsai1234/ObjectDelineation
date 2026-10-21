"""Phase 2 entry gates (experiments/phase2/PLAN.md §5.1, D18) decided from one report of a selected checkpoint.

Gated, per dataset (B4-5): G-format on the shifted val (parse >= 98 %, legal >= 90 %); G-geometry on the shifted
**and** the rescaled val (median IoU_tc per structure, illegal = 0, and the share of samples legal with IoU_tc >= 0.5
on every object); G-prior (b) on both views (paired mean (model - prior) IoU_tc, 95 % CI lower bound > 0, every
structure); G-signal on the shifted val at G = 8 (duplicates <= 10 %, sigma-share >= 50 %); on LIDC the precision
probe (zero-shot pass, else the capacity run's pass, with a warning) and the topology probe (a) (simple controls; the
probe's own pass rule, objdel.eval.probes). Gated once per model: G-vision. Reported, never gating: G-prior (a)
(centred win share, revision 7; gated per dataset as "G-prior-a", every structure's win share >= the threshold, only
when `prior_a_threshold` = prereg.G_PRIOR_A_THRESHOLD is set before V1), topology-probe hole recall (decided by C2),
truncation share and mean / p95 tokens beside G-format, the signal panel, and G-probe (not PASS -> a warning handed
to Phase 3, §5.1). The unanimous-item topology gate (topo_any_unanimous) is removed (D18: vacuous, 0 unanimous
non-simple val items). G-x86, G-cost, G-parity and housekeeping are measured outside the eval harness.

A model passes on a dataset when every gate value of that dataset passes, and overall when it passes on every judged
dataset and G-vision passes. Judged datasets: `datasets` (fixed in advance; a model the precision probe dropped from
LIDC is judged on RIGA alone, datasets=("riga",), and no LIDC value is read) or those in the shifted view's summary.
Anything missing (a view, a dataset, a structure, a probe, the vision check) fails its gate; with no judged dataset
every dataset gate fails (BK2-10). `eligible(summary, vision)` is D12 R-B's branch eligibility (prereg.RB["eligible"]):
G-format, G-signal and G-vision of phase2_entry on one shifted-view summary, so selection and the gate table agree.

The report (C-EVAL; agents G, P, R fill it; every part optional, missing parts fail their gates):

    {"views": {"centred" | "shifted" | "rescaled":           # eval_run's report.json on that view
                   {"summary": objdel.eval.summarize(items), ...}},
     "prior": {"shifted" | "rescaled":                        # G-prior (b), per item IoU_tc = K-sample mean
                   {dataset: {structure: {"diff": mean, "ci": [lo, hi], ...}}},   # paired, nodule-cluster bootstrap
               "centred": {dataset: {structure: {"win_share": share, ...}}}},       # G-prior (a), report-only
     "probes": {"precision": probes.score_precision(rows)     # {"pass": bool, ..., "capacity": {"pass": bool, ...}}
                "topology": probes.score_topology(rows)},     # {"pass": bool (rule (a)), "hole_recall": x, ...}
     "vision": vision.check(...),                             # {"ratio", "tv", "pass", ...}
     "signal": signal.panel(groups),                          # report-only (C-SIGNAL)
     "g_probe": {"verdict": "PASS" | "FAIL" | "INCONCLUSIVE", ...}}

-> {"pass", "datasets": {ds: pass}, "gates": {name: {"requirement", "values": {key: value}, "pass"}},
    "passed": {name: bool} (prereg.c1_stop's input: the C1 gates plus "G-probe" = verdict is PASS),
    "reported": {name: value or {key: value}}, "warnings": [str]}. Gate names are prereg.py's C1_GATES ("G-format",
"G-geometry", "G-prior-b", "G-precision", "G-topology-a", "G-signal", "G-vision"); report-only names its REPORT_ONLY
("G-prior-a", "G-topology-b", "G-probe"). Value keys are "<metric>/<view>/<dataset>[/<structure>]".

The geometry gates read objdel.eval.summarize, where an illegal rollout scores IoU_tc = 0 (objdel.reward.legal_gate;
REVIEW #2).

References:
- Gate table and thresholds: experiments/phase2/PLAN.md §5.1, D17-D20; per-dataset decision: B4-5.
- Duplicate and sigma gates on groups of G = 8 at temperature 1: a group with equal rewards has zero advantage (DAPO
  dynamic sampling, https://arxiv.org/abs/2503.14476); sigma floors as in objdel.reward.advantages_v4.
"""
from collections.abc import Sequence

from objdel.reward import REWARD_V4

FORMAT, LEGAL = 0.98, 0.90
MEDIAN_IOU_TC = {("riga", "disc"): 0.85, ("riga", "cup"): 0.75, ("lidc", "nodule"): 0.60}
IOU_TC_PASS, DUPLICATES, SIGNAL = 0.70, 0.10, 0.50
GEOMETRY_VIEWS = ("shifted", "rescaled")  # G-geometry and G-prior (b)
PROBE_DATASET = "lidc"  # the precision and topology probes are LIDC patches
G = REWARD_V4.group
REQUIREMENTS = {  # keys = experiments/phase2/prereg.py C1_GATES
    "G-format": f"shifted val: parse >= {FORMAT:.0%}, legal >= {LEGAL:.0%}",
    "G-geometry": f"shifted and rescaled val: median IoU_tc (illegal = 0) >= 0.85 disc, 0.75 cup, 0.60 nodule; "
                  f">= {IOU_TC_PASS:.0%} of rollouts legal with IoU_tc >= 0.5 on every object",
    "G-prior-b": "shifted and rescaled val: paired (model - image-blind prior) IoU_tc, 95 % CI lower bound > 0",
    "G-precision": "LIDC precision probe passes (zero-shot, else the capacity run)",
    "G-topology-a": "LIDC topology probe (a): simple controls correct >= 95 %",
    "G-signal": f"shifted val, groups of {G}: duplicate outputs <= {DUPLICATES:.0%}, >= {SIGNAL:.0%} of groups with a "
                "geometric sigma above its floor",
    "G-vision": "projector cosine ratio >= 0.5 and image-use TV >= 0.002 (vs black, vs another item)",
}
PRIOR_A = "G-prior-a: centred val, every structure's win share over the image-blind prior >= {:.0%}"
ELIGIBILITY = ("G-format", "G-signal", "G-vision")  # D12 R-B eligible branch (prereg.RB["eligible"])


def _dig(d, *path):
    for p in path:
        if not isinstance(d, dict):
            return None
        d = d.get(p)
    return d


def _ge(v, t) -> bool:
    return v is not None and v >= t


def _checks(report: dict, ds: str, prior_a: float | None = None) -> list[tuple[str, str, object, bool]]:
    """(gate, value name, value, passed) of one dataset; G-prior-a rows only when its threshold is set."""
    summ = lambda view, *p: _dig(report, "views", view, "summary", "datasets", ds, *p)  # noqa: E731
    names = [n for x, n in MEDIAN_IOU_TC if x == ds]
    rows = [("G-format", f"format/shifted/{ds}", summ("shifted", "format"), FORMAT),
            ("G-format", f"legal/shifted/{ds}", summ("shifted", "legal"), LEGAL)]
    out = [(g, k, v, _ge(v, t)) for g, k, v, t in rows]
    for view in GEOMETRY_VIEWS:
        out += [("G-geometry", f"median/{view}/{ds}/{n}", v := summ(view, "objects", n, "iou_tc_median"),
                 _ge(v, MEDIAN_IOU_TC[ds, n])) for n in names]
        out.append(("G-geometry", f"pass_share/{view}/{ds}", v := summ(view, "iou_tc_pass"), _ge(v, IOU_TC_PASS)))
        for n in names:
            lo = (_dig(report, "prior", view, ds, n, "ci") or [None])[0]
            out.append(("G-prior-b", f"ci_lo/{view}/{ds}/{n}", lo, lo is not None and lo > 0))
    if prior_a is not None:
        out += [("G-prior-a", f"win_share/centred/{ds}/{n}", v := _dig(report, "prior", "centred", ds, n, "win_share"),
                 _ge(v, prior_a)) for n in names]
    g8 = summ("shifted", "group") == G  # the GRPO gates need groups of G
    sig, dup = (summ("shifted", "signal_share"), summ("shifted", "dup_rate")) if g8 else (None, None)
    out += [("G-signal", f"sigma_share/shifted/{ds}", sig, _ge(sig, SIGNAL)),
            ("G-signal", f"duplicates/shifted/{ds}", dup, dup is not None and dup <= DUPLICATES)]
    if ds == PROBE_DATASET:
        p, t = _dig(report, "probes", "precision") or {}, _dig(report, "probes", "topology") or {}
        ok = p.get("pass") is True or _dig(p, "capacity", "pass") is True
        out += [("G-precision", f"zero_shot/{ds}", p.get("pass"), ok),
                ("G-topology-a", f"pass/{ds}", t.get("pass"), t.get("pass") is True)]
    return out


def _reported(report: dict, judged: Sequence[str]) -> dict:
    out = {"G-prior-a": {}, "truncated": {}, "tokens_mean": {}, "tokens_p95": {}}
    for ds in judged:
        for n, v in (_dig(report, "prior", "centred", ds) or {}).items():
            out["G-prior-a"][f"{ds}/{n}"] = _dig(v, "win_share")
        for view in ("centred", *GEOMETRY_VIEWS):
            d = _dig(report, "views", view, "summary", "datasets", ds)
            if d:
                out["truncated"][f"{view}/{ds}"] = d.get("truncated")
                out["tokens_mean"][f"{view}/{ds}"] = _dig(d, "tokens", "mean")
                out["tokens_p95"][f"{view}/{ds}"] = _dig(d, "tokens", "p95")
    lidc = PROBE_DATASET in judged
    out["G-topology-b"] = _dig(report, "probes", "topology", "hole_recall") if lidc else None
    out["G-probe"] = _dig(report, "g_probe", "verdict")
    return out


def phase2_entry(report: dict, datasets: Sequence[str] | None = None, prior_a_threshold: float | None = None) -> dict:
    """The §5.1 entry decision of one report (module docstring: schema and rules); prior_a_threshold =
    prereg.G_PRIOR_A_THRESHOLD (None: G-prior (a) is report-only and no "G-prior-a" gate is written)."""
    shifted = _dig(report, "views", "shifted", "summary", "datasets") or {}
    judged = tuple(datasets) if datasets is not None else tuple(sorted(shifted))
    rows = {ds: _checks(report, ds, prior_a_threshold) for ds in judged}
    reqs = REQUIREMENTS | ({} if prior_a_threshold is None else {"G-prior-a": PRIOR_A.format(prior_a_threshold)})
    gates = {name: {"requirement": req, "values": {}, "pass": bool(judged)} for name, req in reqs.items()}
    for ds_rows in rows.values():
        for gate, key, value, ok in ds_rows:
            gates[gate]["values"][key] = value
            gates[gate]["pass"] &= ok
    vision = _dig(report, "vision") or {}
    gates["G-vision"]["values"] = {k: vision.get(k) for k in ("ratio", "tv")}
    gates["G-vision"]["pass"] = vision.get("pass") is True
    per = {ds: all(ok for *_, ok in r) for ds, r in rows.items()}
    reported = _reported(report, judged)
    warnings = []
    if reported["G-probe"] != "PASS":
        warnings.append(f"G-probe {reported['G-probe'] or 'not measured'}: warning handed to Phase 3 (§5.1)")
    if PROBE_DATASET in judged and _dig(report, "probes", "precision", "pass") is False and \
            _dig(report, "probes", "precision", "capacity", "pass") is True:
        warnings.append("precision probe: zero-shot failed, capacity run passed (LIDC SFT is the bottleneck)")
    passed = {name: g["pass"] for name, g in gates.items()} | {"G-probe": reported["G-probe"] == "PASS"}
    return {"pass": bool(per) and all(per.values()) and gates["G-vision"]["pass"], "datasets": per, "gates": gates,
            "passed": passed, "reported": reported, "warnings": warnings}


def eligible(summary: dict, vision: dict | None, dataset: str = PROBE_DATASET) -> bool:
    """D12 R-B eligibility of one branch: G-format, G-signal (groups of G) and G-vision of phase2_entry on its
    shifted-view `summary` (objdel.eval.summarize output) and `vision` (vision.check output; None or {} = unmeasured,
    not eligible)."""
    passed = phase2_entry({"views": {"shifted": {"summary": summary}}, "vision": vision}, [dataset])["passed"]
    return all(passed[g] for g in ELIGIBILITY)
