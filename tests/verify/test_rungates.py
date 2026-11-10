"""objdel.verify.rungates: the Gate B1 and Gate B3 evaluators on planted paired item values (U49) and the after-V1
decision table of front section 1.1 as a pure function of the Gate B1 rule results (I43)."""
import json

import numpy as np
import pytest

from objdel.verify import rungates as G
from objdel.verify import stats as signal

FAR_MISSES = (8, 9, 10)  # per patient: items 0-3 near (offset 0.05), 4-7 mid (0.25), 8-13 far (0.4)
FROZEN = {"lidc-val-shifted": {"content_sha256": "a" * 64, "known_label_issues": []},
          "lidc-val-edge": {"content_sha256": "b" * 64, "known_label_issues": []}}


def world(patients=60, seed=0):
    """(meta, baseline K-means): far items 8-10 of every patient are misses (0.02), the rest found."""
    rng = np.random.default_rng(seed)
    meta, base = {}, {}
    for p in range(patients):
        for j, off in enumerate([0.05] * 4 + [0.25] * 4 + [0.4] * 6):
            i = f"lidc/p{p:02d}/n000/z{j:03d}"
            meta[i] = {"group": f"P{p:02d}", "offset": off, "padded": j == 1, "fallback": False}
            base[i] = 0.02 if j in FAR_MISSES else float(rng.uniform(0.5, 0.9))
    return meta, base


def converted(base, n=2):
    """Base with the first n far misses of each patient found (0.7)."""
    return {i: 0.7 if int(i[-3:]) in FAR_MISSES[:n] else v for i, v in base.items()}


def run(iou, view="lidc-val-shifted", sha=None, legal=1.0):
    return G.Run(iou=dict(iou), legal=dict.fromkeys(iou, legal), view=view,
                 content_sha256=FROZEN[view]["content_sha256"] if sha is None else sha)


def b1_case(b1=None, frozen=FROZEN, **kw):
    meta, base = world()
    args = dict(b1=run(b1 or converted(base)), b1_edge=run(base, "lidc-val-edge"), a3=run(base),
                a3_edge=run(base, "lidc-val-edge"), r5=run(base, sha="r5"), meta=meta, edge_meta=meta, frozen=frozen,
                chain_ok=True, identity_ok=True)
    return G.evaluate_b1(**(args | kw))


GOOD = dict(iou_ci_lo=0.01, far_miss=0.20, far_miss_ci_hi=-0.01, centred_ci_lo=-0.01, did_ci_lo=0.001,
            edge_ci_lo=-0.01, legal=0.995, view_hash_ok=True)


def test_gate_b1_needs_all_hard_rules():
    assert G.gate_b1(**GOOD, rescaled_ci_lo=-0.05).ok  # a report-only line never decides
    for k, bad in [("iou_ci_lo", -0.001), ("far_miss", 0.206), ("far_miss_ci_hi", 0.0), ("centred_ci_lo", -0.03),
                   ("did_ci_lo", 0.0), ("edge_ci_lo", -0.021), ("legal", 0.98), ("view_hash_ok", False)]:
        assert not G.gate_b1(**{**GOOD, k: bad}).ok, k
    assert G.gate_b1(**{**GOOD, "far_miss": 0.205, "legal": 0.99}).ok  # both boundaries pass
    low = G.gate_b1(**GOOD, identity_ok=False)  # the correlation is not a hard rule (B.6.2)
    assert low.ok and low.r5_line_dropped
    assert G.gate_b1(**{**GOOD, "edge_ci_lo": None}).rules["edge"] is None  # no edge baseline: never a pass
    assert not G.gate_b1(**{**GOOD, "edge_ci_lo": None}).ok


def test_planted_effect_passes_and_names_the_a3_baseline():
    d = b1_case()
    assert d["result"] == "pass" and all(d["rules"].values())
    assert d["decision"]["row"] == 1 and d["decision"]["next"] == "B-S2"
    assert d["baseline"] == "A3 (item 5)" and d["edge_baseline"] == "A3 (item 4)"
    assert d["values"]["far_miss"] == pytest.approx(1 / 6)
    assert d["report"]["r5_line"] is not None and d["report"]["effect"] == "data effect"


def test_view_hash_mismatch_stops():
    frozen = FROZEN | {"lidc-val-edge": FROZEN["lidc-val-edge"] | {"content_sha256": "c" * 64}}
    d = b1_case(frozen=frozen)
    assert d["rules"]["view_hash"] is False and d["result"] == "fail"
    assert d["decision"]["row"] == 4 and d["decision"]["next"] == "stop"


def test_equal_gain_on_near_and_far_fails_the_difference_in_differences():
    meta, base = world()
    d = b1_case(b1={i: v + 0.05 for i, v in base.items()})
    assert d["rules"]["iou"] and d["rules"]["did"] is False
    assert d["decision"]["row"] == 3 and d["decision"]["final"] == "R5, unless the 4B arm wins on val"


def test_low_near_headroom_marks_the_did_rule_uninformative():
    d, ref = b1_case(near_headroom_hi=0.0129), b1_case(near_headroom_hi=0.013)
    assert d["report"]["effect"] == "B1 vs R5" and "uninformative" in d["report"]["did_note"]
    assert ref["report"]["effect"] == "data effect"
    assert d["result"] == ref["result"] == "pass"  # the pass/fail is unchanged


def test_item_5_missing_falls_back_to_r5_with_the_hash_chain():
    d = b1_case(a3=None)
    assert d["result"] == "pass" and d["baseline"] == "R5 items.jsonl"
    assert d["identity"] == "hash only, correlation not run" and d["edge_baseline"] == "A3 (item 4)"
    broken = b1_case(a3=None, chain_ok=False)
    assert broken["result"] == "wait for V2 item 1b" and broken["decision"]["row"] is None


def test_items_4_and_5_missing_leave_the_gate_open():
    d = b1_case(a3=None, a3_edge=None)
    assert d["result"] == "open: edge baseline missing" and d["rules"]["edge"] is None
    assert all(d["rules"][r] for r in ("iou", "far_miss", "centred", "did", "legal", "view_hash"))  # still computed
    assert d["decision"]["next"] == "V2 item 1b"
    no_raise = b1_case(a3=None, a3_edge=None, raise_planned=False)
    assert no_raise["decision"]["row"] == 4 and no_raise["decision"]["reason"] == "edge baseline missing"


def test_failing_correlation_keeps_a3_and_drops_the_r5_line():
    d = b1_case(identity_ok=False)
    assert d["result"] == "pass" and d["baseline"] == "A3 (item 5)" and d["report"]["r5_line"] is None


def test_known_label_issues_are_excluded_from_both_arms():
    meta, base = world()
    bad = [i for i in base if int(i[-3:]) < 4 and int(i.split("/")[1][1:]) < 10]  # 40 near items of 10 patients
    b1 = converted(base) | dict.fromkeys(bad, 0.0)
    frozen = FROZEN | {"lidc-val-shifted": FROZEN["lidc-val-shifted"] | {"known_label_issues": bad}}
    assert b1_case(b1=b1)["rules"]["centred"] is False  # the planted items would fail the guard
    d = b1_case(b1=b1, frozen=frozen)
    assert d["result"] == "pass" and d["values"]["items"] == len(base) - len(bad)
    assert d["report"]["known_label_issues"]["items"] == len(bad)
    assert d["report"]["known_label_issues"]["iou"][0] < d["values"]["iou"][0]  # the line with them


def test_report_lines_never_change_the_result():
    meta, base = world()
    sel = [i for i in base if i.startswith("lidc/p0")]
    pairs = {"lidc-val-rescaled": (run({i: 0.0 for i in base}, sha="x"), run(base, sha="y"))}
    d = b1_case(r5_selection=sel, report_pairs=pairs, noise_mean_sd=0.004)
    assert d["result"] == b1_case()["result"] == "pass"
    r = d["report"]
    assert r["pairs"]["lidc-val-rescaled"][0] == pytest.approx(-np.mean(list(base.values())))
    assert r["without_r5_selection"]["items"] == len(base) - len(sel)
    assert r["held_out"] == {"250": 0.634, "770": 0.659} and r["noise_mean_sd"] == 0.004
    assert set(r["gains"]) == {"far", "near"} and r["edge_padded_or_fallback"]["items"] == 60


def test_partial_run_is_labelled_and_needs_review():
    d = b1_case(partial_steps=210)
    assert d["result"] == "pass" and d["label"] == "partial run, 210 steps"
    assert d["decision"]["next"] == "B-S2 after the owner reviews the partial run"
    assert b1_case(partial_steps=260)["label"] is None


def test_interval_is_percentile_at_50_patients_and_cluster_t_below():
    meta, base = world(patients=60)
    ids = list(base)
    groups = np.array([meta[i]["group"] for i in ids])
    x = np.array([base[i] for i in ids])
    W = G.weights(ids, groups)
    allm = np.ones(len(ids), bool)
    reps = W @ x / W.sum(1)
    assert G.interval(x, groups, allm, W)[1:] == pytest.approx(np.quantile(reps, [0.025, 0.975]).tolist())
    few = groups < "P49"  # 49 patients
    assert G.interval(x, groups, few, W) == pytest.approx(signal.cluster_t(x[few], groups[few]))


def test_load_run_reads_the_k_mean_and_the_view_identity(tmp_path):
    (tmp_path / "items.jsonl").write_text("".join(json.dumps(r) + "\n" for r in [
        {"item_id": "lidc/a/n000/z001", "iou_tc": 0.5, "legal": 1.0, "sel_medoid": 0.9}]))
    (tmp_path / "final_eval.json").write_text(json.dumps({"view": "/x/eval_views/lidc-val-shifted",
                                                          "view_content_sha256": "d" * 64}))
    r = G.load_run(tmp_path)
    assert r == G.Run({"lidc/a/n000/z001": 0.5}, {"lidc/a/n000/z001": 1.0}, "lidc-val-shifted", "d" * 64)
    assert G.load_run(tmp_path / "missing") is None  # a step that gave no output


def b3_case(b3=None, b1=None, legal=1.0, **kw):
    meta, base = world(patients=76)
    select = {f"P{p:02d}" for p in range(0, 76, 2)}  # 38 patients: cluster_t
    args = dict(b3=run(b3 or converted(base, 2), legal=legal), b3_edge=run(base, "lidc-val-edge"), b1=run(b1 or base),
                b1_edge=run(base, "lidc-val-edge"), meta=meta, edge_meta=meta, frozen=FROZEN, select=select)
    return G.evaluate_b3(**(args | kw))


def test_gate_b3_rules():
    d = b3_case()
    assert d["result"] == "pass" and d["values"]["patients"] == 38 and d["values"]["interval"] == "cluster_t"
    assert d["values"]["miss_share"] == pytest.approx([1 / 14, 3 / 14])
    meta, base = world(patients=76)
    assert b3_case(b3=converted(base, 0))["rules"]["miss_share"] is False  # no fall
    assert b3_case(b3=converted(base, 1))["result"] == "pass"  # 2/14 = 0.143 <= 0.75 x 3/14 = 0.161
    worse = {i: (v - 0.05 if v > 0.3 else v) for i, v in converted(base, 2).items()}
    assert b3_case(b3=worse)["rules"]["both_found_iou"] is False
    lower_edge = run({i: v - 0.05 for i, v in base.items()}, "lidc-val-edge")
    assert b3_case(b3_edge=lower_edge)["rules"]["edge"] is False
    assert b3_case(legal=0.99)["result"] == "pass"
    low = b3_case(legal=0.9899)
    assert low["rules"]["legal"] is False and low["result"] == "fail"  # one failing rule fails the gate


DECISIONS = [
    (dict(), 1, "B-S2", "B1.8 runs", "best of B1 and B3 on val-select by the Gate B3 rules (a tie goes to B1)"),
    (dict(far_miss=False), 2, "B3 on B1 rows", "B1.8 runs",
     "best of B1 and B3 on val-select by the Gate B3 rules (a tie goes to B1)"),
    (dict(iou=False), 3, "no B3 spend", "B1.8 runs", "R5, unless the 4B arm wins on val"),
    (dict(did=False, far_miss=False), 3, "no B3 spend", "B1.8 runs", "R5, unless the 4B arm wins on val"),
    (dict(edge=False), 4, "stop", "B1.8 waits for the owner review", "R5"),
    (dict(legal=False, iou=False), 4, "stop", "B1.8 waits for the owner review", "R5"),
    (dict(centred=False), 4, "stop", "B1.8 waits for the owner review", "R5"),
    (dict(view_hash=False), 4, "stop", "B1.8 waits for the owner review", "R5"),
]


@pytest.mark.parametrize("fail,row,step,b18,final", DECISIONS)
def test_decision_table(fail, row, step, b18, final):
    rules = dict.fromkeys(G.B1_RULES, True) | fail
    d = G.decide(rules)
    assert (d["row"], d["next"], d["b18"], d["final"]) == (row, step, b18, final)


def test_decision_table_partial_and_edge_missing():
    ok = dict.fromkeys(G.B1_RULES, True)
    assert G.decide(ok, partial_steps=210)["next"] == "B-S2 after the owner reviews the partial run"
    assert G.decide(ok | {"iou": False}, partial_steps=210)["row"] == 3  # the row the rules give applies
    assert G.decide(ok | {"edge": None})["next"] == "V2 item 1b"
    stop = G.decide(ok | {"edge": None}, raise_planned=False)
    assert stop["row"] == 4 and stop["reason"] == "edge baseline missing"
    assert G.decide(ok | {"edge": None, "legal": False})["row"] == 4  # a failed stop rule stops at once
