"""Phase 2 entry gates (plan §7 table) that an evaluation report can decide, checked on the validation split after
SFT: format, geometry (medians per structure and the IoU_tc >= 0.5 share), topology on unanimous items and the GRPO
signal (duplicates and sigma above the floor, on groups of G = 8 sampled at temperature 1, which the caller ensures).

Per dataset (B4-5): a model passes on a dataset when every gate value of that dataset passes; it passes overall when
it passes on every dataset it is judged on. By default these are the datasets in the report; the caller names them
with `datasets` when the set is fixed in advance (a requested dataset missing from the report fails). Plan §7 drops
a model that fails the precision probe from LIDC; such a model is judged on RIGA alone (`datasets=("riga",)`), and
the run records the drop. A structure missing from a judged dataset fails its gate. The remaining Phase 2 gates
(reward speed on x86, cost, precision probe, housekeeping) are measured outside the eval harness.

The geometry gates read objdel.eval.summarize, where an illegal rollout scores IoU_tc = 0 (the reward's legality gate,
objdel.reward.legal_gate; REVIEW #2): the medians count it as 0 and the pass share needs legal and IoU_tc >= 0.5 for
every object.

References:
- Gate table: experiments/phase1/PLAN.md §7 "Phase 2 entry gates"; per-dataset decision and the LIDC drop: §7
  precision probe, EXECUTION.md §3 (M8, B4-5).
- Duplicate and sigma gates on groups of G = 8 at temperature 1: a group with equal rewards has zero advantage (DAPO
  dynamic sampling, https://arxiv.org/abs/2503.14476); sigma floors as in objdel.reward.advantages_v4.
"""
from collections.abc import Sequence

from objdel.reward import REWARD_V4

FORMAT, LEGAL = 0.98, 0.90
MEDIAN_IOU_TC = {("riga", "disc"): 0.85, ("riga", "cup"): 0.75, ("lidc", "nodule"): 0.60}
IOU_TC_PASS, TOPOLOGY, DUPLICATES, SIGNAL = 0.70, 0.90, 0.10, 0.50
G = REWARD_V4.group
REQUIREMENTS = {
    "format": f"parse >= {FORMAT:.0%}, legal >= {LEGAL:.0%}",
    "geometry_median": "median IoU_tc (illegal = 0) >= 0.85 disc, 0.75 cup, 0.60 nodule",
    "geometry_pass": f">= {IOU_TC_PASS:.0%} of rollouts legal with IoU_tc >= 0.5 (every object)",
    "topology": f"topology-any >= {TOPOLOGY:.0%} on unanimous items",
    "grpo_duplicates": f"duplicate outputs <= {DUPLICATES:.0%} within groups of {G}",
    "grpo_sigma": f">= {SIGNAL:.0%} of groups of {G} with a geometric sigma above its floor",
}


def _checks(d: dict | None, ds: str) -> list[tuple[str, str, object, bool]]:
    """(gate, value name, value, passed) of one dataset's summary d (None when the dataset is missing)."""
    get = lambda *path: _dig(d, path)  # noqa: E731
    g8 = get("group") == G  # the GRPO gates need groups of G
    rows = [("format", f"{ds}/format", get("format"), FORMAT), ("format", f"{ds}/legal", get("legal"), LEGAL)]
    rows += [("geometry_median", f"{ds}/{n}", get("objects", n, "iou_tc_median"), t)
             for (x, n), t in MEDIAN_IOU_TC.items() if x == ds]
    rows += [("geometry_pass", ds, get("iou_tc_pass"), IOU_TC_PASS),
             ("topology", ds, get("topo_any_unanimous"), TOPOLOGY),
             ("grpo_sigma", ds, get("signal_share") if g8 else None, SIGNAL)]
    out = [(gate, k, v, v is not None and v >= t) for gate, k, v, t in rows]
    dup = get("dup_rate") if g8 else None
    return out + [("grpo_duplicates", ds, dup, dup is not None and dup <= DUPLICATES)]


def phase2_entry(report: dict, datasets: Sequence[str] | None = None) -> dict:
    """report["summary"] = objdel.eval.summarize(items); datasets: the datasets the model is judged on (default: those
    in the report). -> {"pass", "datasets": {ds: pass}, "gates": {name: {pass, values, requirement}}}."""
    d = report["summary"]["datasets"]
    judged = tuple(datasets) if datasets is not None else tuple(sorted(d))
    rows = {ds: _checks(d.get(ds), ds) for ds in judged}
    gates = {name: {"requirement": req, "values": {}, "pass": True} for name, req in REQUIREMENTS.items()}
    for ds_rows in rows.values():
        for gate, key, value, ok in ds_rows:
            gates[gate]["values"][key] = value
            gates[gate]["pass"] &= ok
    per = {ds: all(ok for *_, ok in r) for ds, r in rows.items()}
    return {"pass": bool(per) and all(per.values()), "datasets": per, "gates": gates}


def _dig(d, path):
    for p in path:
        if d is None:
            return None
        d = d.get(p)
    return d
