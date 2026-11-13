"""Gate B1 and Gate B3 (Implementation.md B.6.2) on the runs' per-item values, and the after-V1 decision table
(front section 1.1, B.4.1) as a pure function of the Gate B1 rule results. The CLI is src/objdel/cli/eval/gates.py.

Inputs per run: final_eval's `items.jsonl` (the K-mean IoU_tc and legal share per item; any other column, a medoid
for one, is ignored) and `final_eval.json` (`view`, `view_content_sha256`); per view the frozen view's records
(patient `group_id`, `view.offset_frac`, `view.padded`, `view.fallback`) and its frozen_views.json entry (content hash,
`known_label_issues`). Known label issues leave both arms of every rule; the report prints the IoU line with them.

Decision interval of a mean rule: the 95 % percentile patient bootstrap (objdel.verify.stats.draws with the patient
key, B = 2,000, seed 0; one set of draws per view, so the difference in differences uses the same draws) when its
subset holds >= 50 patients, else the CR1 cluster_t interval (A1 step 5). Decision values are rounded to
stats.DIGITS before a threshold.

Gate B1 (whole val; hard rules, all must hold): `iou` CI lower end of (B1 - base) > 0; `far_miss` B1 miss share at
offset > 0.3 <= 20.5 % and the CI upper end of the paired miss difference < 0; `centred` IoU at offset < 0.15, CI lower
end > -0.02; `did` CI lower end of (far gain - near gain) > 0; `edge` edge view, non-padded and non-fallback items, CI
lower end > -0.02; `legal` B1 legal share >= 0.99; `view_hash` every B1 and A3 `view_content_sha256` equals its frozen
entry. Baseline (front 6.2): A3's item 5 when it ran; else R5's items.jsonl when the hash chain of B1.0 holds
("hash only, correlation not run"), else the gate waits for V2 item 1b. The edge rule pairs with item 4; without it
the gate is "open: edge baseline missing", never a pass. The view-identity correlation is not a hard rule: when it
fails, only the R5 report line is dropped. A near-subset headroom upper bound below 0.013 marks the
difference-in-differences rule uninformative and names the result "B1 vs R5", not "data effect"; pass/fail is
unchanged. Report-only lines never change the result.

Gate B3 (val-select patients; B3 vs the B1 re-eval of the same rental): `miss_share` falls >= 25 % relative (point)
and the CI upper end of (B3 - B1) miss difference < 0; `both_found_iou` IoU on items found by both runs, CI lower end
> -0.01; `edge` as Gate B1; `legal` B3 legal share >= 0.99.
"""

from pathlib import Path
from typing import NamedTuple

import numpy as np

from objdel.core.io import read_json, read_jsonl, read_split
from objdel.verify import stats as signal

MISS, DIGITS = signal.MISS, signal.DIGITS
FAR, NEAR = 0.3, 0.15
FAR_MISS_MAX, NEAR_FLOOR, EDGE_FLOOR, LEGAL_MIN = 0.205, -0.02, -0.02, 0.99
B3_FALL, B3_BOTH_FLOOR = 0.75, -0.01
HEADROOM_MIN = 0.013  # A2's near-subset headroom upper bound below which the DiD rule is uninformative
MIN_PATIENTS = 50
PLANNED_STEPS = 260
TARGET = 0.66  # val IoU target: reported, not gated
HELD_OUT = {"250": 0.634, "770": 0.659}  # rl-full-3ep held-out 128-item points at steps 250 and 770
B1_RULES = ("iou", "far_miss", "centred", "did", "edge", "legal", "view_hash")
STOP_RULES = ("centred", "edge", "legal", "view_hash")  # front 1.1 row 4
BEST = "best of B1 and B3 on val-select by the Gate B3 rules (a tie goes to B1)"
R5_FINAL = "R5, unless the 4B arm wins on val"


class Run(NamedTuple):
    """Per-item scores of one run."""
    iou: dict[str, float]  # item_id -> K-mean IoU_tc
    legal: dict[str, float]
    view: str  # the view folder's name
    content_sha256: str | None  # final_eval.json view_content_sha256 (R5's run has none)


class GateB1(NamedTuple):
    """The Gate B1 verdict and its per-rule results."""
    ok: bool
    rules: dict[str, bool | None]  # None: the edge rule has no baseline
    r5_line_dropped: bool
    reported: dict


def load_run(out: Path) -> Run | None:
    """A final_eval output folder; None when the step gave no output (no items.jsonl)."""
    out = Path(out)
    if not (out / "items.jsonl").is_file():
        return None
    rows, fe = read_jsonl(out / "items.jsonl"), read_json(out / "final_eval.json")
    return Run({r["item_id"]: r["iou_tc"] for r in rows}, {r["item_id"]: r["legal"] for r in rows},
               Path(fe["view"]).name, fe.get("view_content_sha256"))


def view_meta(view: Path, split: str = "val") -> dict[str, dict]:
    """item_id -> {group, offset, padded, fallback} from a view's LIDC records."""
    rows = read_split(view, "lidc", split)
    return {r["item_id"]: {"group": r["group_id"], "offset": r["view"]["offset_frac"], "padded": r["view"]["padded"],
                           "fallback": r["view"]["fallback"]} for r in rows}


def hash_chain_ok(pre: dict, entry: dict) -> bool:
    """B1.0's chain for the R5 fallback: the R5 stamp's data hash equals the recomputed one and the pre-rebuild view
    hashes equal the frozen entry (a rebuild that changed the content leaves item 5 the only baseline)."""
    return (pre.get("r5_stamp_check", {}).get("equal") is True
            and all(pre.get(k) == entry[k] for k in ("content_sha256", "target_set_sha256")))


def weights(ids: list[str], groups: np.ndarray) -> np.ndarray:
    """(B, N) item multiplicities of the patient draws."""
    D = signal.draws(ids, key=dict(zip(ids, groups.tolist(), strict=True)).__getitem__)
    return np.stack([np.bincount(d, minlength=len(ids)) for d in D]).astype(float)


def _reps(x: np.ndarray, mask: np.ndarray, W: np.ndarray) -> np.ndarray:
    Wm = W[:, mask]
    with np.errstate(invalid="ignore", divide="ignore"):
        return Wm @ x[mask] / Wm.sum(1)


def _round(v) -> list[float]:
    return [round(float(x), DIGITS) for x in v]


def interval(x: np.ndarray, groups: np.ndarray, mask: np.ndarray, W: np.ndarray) -> list[float]:
    """[mean, lo, hi] of x over mask: percentile patient interval from W at >= 50 patients, CR1 cluster_t below."""
    if len(np.unique(groups[mask])) < MIN_PATIENTS:
        return _round(signal.cluster_t(x[mask], groups[mask]))
    q = (1 - signal.LEVEL) / 2
    return _round([x[mask].mean(), *np.nanquantile(_reps(x, mask, W), [q, 1 - q])])


class _Pair(NamedTuple):
    new: np.ndarray
    old: np.ndarray
    groups: np.ndarray
    W: np.ndarray
    ids: list[str]


def _pair(new: Run, old: Run, meta: dict, drop=frozenset(), keep=None) -> _Pair:
    """Both runs' values on new's items minus `drop` (and within `keep`); ValueError unless old holds them all."""
    ids = [i for i in new.iou if i not in drop and (keep is None or keep(i))]
    if lost := [i for i in ids if i not in old.iou]:
        raise ValueError(f"{len(lost)} items of {new.view} are not in the baseline, first {lost[0]}")
    groups = np.array([meta[i]["group"] for i in ids])
    return _Pair(np.array([new.iou[i] for i in ids]), np.array([old.iou[i] for i in ids]), groups,
                 weights(ids, groups), ids)


def _issues(frozen: dict, view: str) -> frozenset:
    return frozenset(frozen.get(view, {}).get("known_label_issues", ()))


def _edge(new: Run, old: Run | None, meta: dict, frozen: dict, keep=None) -> tuple[list | None, dict]:
    """The edge rule's interval on non-padded, non-fallback items, and the padded / fallback report line."""
    if old is None:
        return None, {"items": None, "mean": None}
    p = _pair(new, old, meta, _issues(frozen, new.view), keep)
    side = np.array([meta[i]["padded"] or meta[i]["fallback"] for i in p.ids])
    d = p.new - p.old
    return interval(d, p.groups, ~side, p.W), {"items": int(side.sum()),
                                                "mean": float(d[side].mean()) if side.any() else None}


def _legal(run: Run, ids: list[str]) -> float:
    return round(float(np.mean([run.legal[i] for i in ids])), DIGITS)


def gate_b1(*, iou_ci_lo, far_miss, far_miss_ci_hi, centred_ci_lo, did_ci_lo, edge_ci_lo, legal, view_hash_ok,
            identity_ok=None, **reported) -> GateB1:
    """The Gate B1 hard rules on their values; keyword values beyond them are report lines and never decide."""
    rules = {"iou": iou_ci_lo > 0, "far_miss": far_miss <= FAR_MISS_MAX and far_miss_ci_hi < 0,
             "centred": centred_ci_lo > NEAR_FLOOR, "did": did_ci_lo > 0,
             "edge": None if edge_ci_lo is None else edge_ci_lo > EDGE_FLOOR, "legal": legal >= LEGAL_MIN,
             "view_hash": bool(view_hash_ok)}
    return GateB1(all(v is True for v in rules.values()), rules, identity_ok is False, reported)


def decide(rules: dict, *, partial_steps: int | None = None, chain_broken: bool = False,
           raise_planned: bool = True) -> dict:
    """Front section 1.1's row for Gate B1 rule results ({rule: True | False | None}; edge None = no baseline)."""
    def row(n, step, final, reason=None):
        b18 = "B1.8 waits for the owner review" if n == 4 else "B1.8 runs" if n else None
        return {"row": n, "next": step, "b18": b18, "final": final, "reason": reason}

    if chain_broken:
        return row(None, "V2 item 1b", None, "R5 hash chain broken: A3-val comes from V2 item 1b")
    if stop := [r for r in STOP_RULES if rules[r] is False]:
        return row(4, "stop", "R5", "failed: " + ", ".join(stop))
    if rules["edge"] is None:
        return row(None, "V2 item 1b", None, "edge baseline missing") if raise_planned else row(
            4, "stop", "R5", "edge baseline missing")
    if not (rules["iou"] and rules["did"]):
        return row(3, "no B3 spend", R5_FINAL)
    if not rules["far_miss"]:
        return row(2, "B3 on B1 rows", BEST)
    partial = partial_steps is not None and partial_steps < PLANNED_STEPS
    return row(1, "B-S2 after the owner reviews the partial run" if partial else "B-S2", BEST)


def evaluate_b1(b1: Run, b1_edge: Run, *, a3: Run | None, a3_edge: Run | None, r5: Run, meta: dict, edge_meta: dict,
                frozen: dict, chain_ok: bool, identity_ok: bool | None, near_headroom_hi: float | None = None,
                raise_planned: bool = True, partial_steps: int | None = None, r5_selection=(),
                report_pairs: dict | None = None, noise_mean_sd: float | None = None) -> dict:
    """Gate B1 on the runs (module docstring); identity_ok: whether A3 and R5 gave identical answers on the shared
    views (None: not checked)."""
    base = a3 or r5
    issues = _issues(frozen, b1.view)
    p = _pair(b1, base, meta, issues)
    d, ones = p.new - p.old, np.ones(len(p.ids), bool)
    off = np.array([meta[i]["offset"] for i in p.ids])
    far, near = off > FAR, off < NEAR
    miss = (p.new < MISS).astype(float) - (p.old < MISS)
    did_reps = _reps(d, far, p.W) - _reps(d, near, p.W)
    q = (1 - signal.LEVEL) / 2
    edge, side = _edge(b1_edge, a3_edge, edge_meta, frozen)
    values = {"items": len(p.ids), "patients": len(np.unique(p.groups)), "iou": interval(d, p.groups, ones, p.W),
              "far_miss": round(float((p.new[far] < MISS).mean()), DIGITS),
              "far_miss_diff": interval(miss, p.groups, far, p.W), "centred": interval(d, p.groups, near, p.W),
              "did": _round([d[far].mean() - d[near].mean(), *np.nanquantile(did_reps, [q, 1 - q])]),
              "edge": edge, "legal": _legal(b1, p.ids)}
    hashed = [r for r in (b1, b1_edge, a3, a3_edge) if r is not None]
    view_hash = all(r.content_sha256 == frozen.get(r.view, {}).get("content_sha256") for r in hashed)
    gate = gate_b1(iou_ci_lo=values["iou"][1], far_miss=values["far_miss"], far_miss_ci_hi=values["far_miss_diff"][2],
                   centred_ci_lo=values["centred"][1], did_ci_lo=values["did"][1],
                   edge_ci_lo=None if edge is None else edge[1], legal=values["legal"], view_hash_ok=view_hash,
                   identity_ok=identity_ok if a3 else None)
    chain_broken = a3 is None and not chain_ok
    decision = decide(gate.rules, partial_steps=partial_steps, chain_broken=chain_broken, raise_planned=raise_planned)
    result = ("wait for V2 item 1b" if chain_broken else "fail" if False in gate.rules.values()
              else "open: edge baseline missing" if gate.rules["edge"] is None else "pass")
    uninformative = near_headroom_hi is not None and near_headroom_hi < HEADROOM_MIN

    def iou_on(**kw):
        q = _pair(b1, base, meta, **kw)
        return {"items": len(q.ids), "iou": interval(q.new - q.old, q.groups, np.ones(len(q.ids), bool), q.W)}

    def r5_line():
        q = _pair(b1, r5, meta, issues)
        return interval(q.new - q.old, q.groups, np.ones(len(q.ids), bool), q.W)

    pairs = {}
    for name, (new, old) in (report_pairs or {}).items():
        q = _pair(new, old, meta, issues)
        pairs[name] = interval(q.new - q.old, q.groups, np.ones(len(q.ids), bool), q.W)
    report = {"effect": "B1 vs R5" if uninformative else "data effect",
              "did_note": (f"uninformative: A2 near headroom upper bound {near_headroom_hi} < {HEADROOM_MIN}"
                           if uninformative else None),
              "r5_line": r5_line() if a3 and identity_ok else None,
              "known_label_issues": iou_on() | {"items": len(issues)},
              "without_r5_selection": iou_on(drop=issues | set(r5_selection)),
              "gains": {"far": float(d[far].mean()), "near": float(d[near].mean())},
              "edge_padded_or_fallback": side, "pairs": pairs, "held_out": HELD_OUT, "noise_mean_sd": noise_mean_sd,
              "target": {"iou": TARGET, "b1": float(p.new.mean())}}
    if not a3:
        identity = "hash only, correlation not run"
    else:
        identity = "view identity holds" if identity_ok else "view identity fails: R5 line dropped"
    return {"result": result, "label": f"partial run, {partial_steps} steps" if decision["row"] and partial_steps
            is not None and partial_steps < PLANNED_STEPS else None,
            "baseline": "A3 (item 5)" if a3 else "R5 items.jsonl", "identity": identity,
            "edge_baseline": "A3 (item 4)" if a3_edge else "missing", "rules": gate.rules, "decision": decision,
            "values": values, "report": report}


def evaluate_b3(b3: Run, b3_edge: Run, *, b1: Run, b1_edge: Run, meta: dict, edge_meta: dict, frozen: dict,
                select: set[str]) -> dict:
    """Gate B3 on val-select (patients in `select`), B3 vs the B1 re-eval of the same rental."""
    keep = lambda i: meta[i]["group"] in select  # noqa: E731
    p = _pair(b3, b1, meta, _issues(frozen, b3.view), keep)
    m3, m1 = p.new < MISS, p.old < MISS
    ones = np.ones(len(p.ids), bool)
    share = [round(float(m3.mean()), DIGITS), round(float(m1.mean()), DIGITS)]
    miss = interval(m3.astype(float) - m1, p.groups, ones, p.W)
    both = interval(p.new - p.old, p.groups, ~m3 & ~m1, p.W)
    edge, side = _edge(b3_edge, b1_edge, edge_meta, frozen, lambda i: edge_meta[i]["group"] in select)
    legal = _legal(b3, p.ids)
    rules = {"miss_share": share[0] <= round(B3_FALL * share[1], DIGITS) and miss[2] < 0,
             "both_found_iou": both[1] > B3_BOTH_FLOOR, "edge": edge[1] > EDGE_FLOOR, "legal": legal >= LEGAL_MIN}
    patients = len(np.unique(p.groups))
    return {"result": "pass" if all(rules.values()) else "fail", "rules": rules,
            "values": {"items": len(p.ids), "patients": patients,
                       "interval": "percentile" if patients >= MIN_PATIENTS else "cluster_t", "miss_share": share,
                       "miss_diff": miss, "both_found_iou": both, "edge": edge, "legal": legal},
            "report": {"edge_padded_or_fallback": side}}
