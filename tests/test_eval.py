"""Evaluation harness (plan §4.7, gate R10): metrics, hack audit and Phase 2 entry gates on synthetic items."""
import json
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from oracles import reward_harness as H

from objdel.eval import gates, hack_audit, metrics, select, summarize
from objdel.eval.metrics import _pairwise_iou, ged, min_mode_retention, pass_at_k, soft_dice, tv_js
from objdel.reward import ExpertCache
from objdel.targets import BuildConfig, build_target

S = 50  # frame side (divides 1000)


def disk(cy, cx, r, side=S):
    yy, xx = np.ogrid[:side, :side]
    return (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r


def nodule_raters():
    """4 raters, two of them with a hole: rater topology histogram {(0,): 1/2, (1,): 1/2}."""
    return {"nodule": np.stack([disk(25, 25, 12), disk(25, 27, 14), disk(24, 25, 11) & ~disk(24, 25, 3),
                                disk(26, 24, 13) & ~disk(25, 25, 3)])}


def riga_raters():
    """3 raters of a disc with a cup inside, one topology answer."""
    side = 100
    return {"disc": np.stack([disk(50, 50, r, side) for r in (30, 32, 34)]),
            "cup": np.stack([disk(50, 51, r, side) for r in (12, 14, 16)])}


def targets(raters, names):
    ma = {n: ExpertCache.build(raters[n]).min_area for n in names}
    k = len(raters[names[0]])
    return [build_target({n: raters[n][r] for n in names}, names, min_area=ma).text for r in range(k)]


def text_of(masks: dict, names, raters):
    ma = {n: ExpertCache.build(raters[n]).min_area for n in names}
    return build_target(masks, names, min_area=ma).text


@pytest.fixture(scope="module")
def nod():
    raters = nodule_raters()
    return raters, targets(raters, ("nodule",))


@pytest.fixture(scope="module")
def riga():
    raters = riga_raters()
    return raters, targets(raters, ("disc", "cup"))


# ---------------------------------------------------------------- building blocks

@settings(max_examples=300, deadline=None)
@given(st.integers(1, 12).flatmap(lambda n: st.tuples(st.just(n), st.integers(0, n), st.integers(1, n))))
def test_pass_at_k_equals_brute_force(nck):
    n, c, k = nck
    ok = [True] * c + [False] * (n - c)
    brute = sum(any(ok[i] for i in s) for s in combinations(range(n), k)) / comb(n, k)
    assert pass_at_k(n, c, k) == pytest.approx(brute, abs=1e-12)


def test_tv_js_hand_cases():
    f = {"a": 0.75, "b": 0.25}
    assert tv_js(f, f) == (0.0, 0.0)
    tv, js = tv_js({"a": 1.0}, f)
    want_js = 0.5 * np.log2(1 / 0.875) + 0.5 * (0.75 * np.log2(0.75 / 0.875) + 0.25 * np.log2(0.25 / 0.125))
    assert tv == pytest.approx(0.25) and js == pytest.approx(want_js)
    assert tv_js({"c": 1.0}, f) == pytest.approx((1.0, 1.0))  # disjoint supports: both metrics at their maximum


def test_min_mode_retention():
    f = {"a": 0.75, "b": 0.25}
    assert min_mode_retention(f, f) == 1.0
    assert min_mode_retention({"a": 1.0}, f) == 0.0  # collapse to the majority loses the minority mode
    assert min_mode_retention({"a": 0.5, "b": 0.5}, f) == pytest.approx(2 / 3)  # the majority keeps 2/3 of its share
    assert min_mode_retention({"a": 0.375, "b": 0.125, -1: 0.5}, f) == pytest.approx(1.0)  # invalid mass is excluded
    assert min_mode_retention({-1: 1.0}, f) == 0.0


def test_soft_dice_mrnet_protocol_hand_case():
    """Raters {0, 1} and {0}: mean 1 / 0.5. One sample {0}: Dice 2/3 at t = 0.1, 0.3, 0.5 and 1 at 0.7, 0.9."""
    raters = np.array([[[1, 1]], [[1, 0]]], bool)
    assert soft_dice(np.array([[[1, 0]]], bool), raters) == pytest.approx(0.8)
    assert soft_dice(raters, raters) == 1.0
    assert soft_dice(np.zeros((1, 1, 2), bool), np.zeros((2, 1, 2), bool)) == 1.0  # empty vs empty agrees


def test_ged_hand_cases():
    a, b = np.eye(4, dtype=bool)[None], np.ones((1, 4, 4), bool)
    Y = np.concatenate([a, b])
    assert ged(Y, Y) == pytest.approx(0.0)
    assert ged(np.concatenate([a, a]), Y) == pytest.approx(0.75 / 2)  # d(a, b) = 1 - 4/16; GED = d / 2
    z = np.zeros((1, 3, 3), bool)
    assert ged(z, np.concatenate([z, z])) == 0.0  # two empty masks are at distance 0 (Kohl et al., App. B)


@settings(max_examples=100, deadline=None)
@given(st.integers(0, 2**31 - 1))
def test_pairwise_iou_equals_brute_force(seed):
    """Random masks with duplicates and empty masks; the dedup + crop + matrix product path is exact."""
    rng = np.random.default_rng(seed)
    base = rng.random((rng.integers(1, 5), 7, 9)) < rng.uniform(0, 0.6)
    masks = base[rng.integers(0, len(base), rng.integers(1, 9))]
    want = np.array([[np.sum(a & b) / np.sum(a | b) if np.any(a | b) else 1.0 for b in masks] for a in masks])
    assert np.array_equal(_pairwise_iou(masks), want)


# ---------------------------------------------------------------- R10: perfect values on the raters' own targets

def assert_perfect(m, names):
    s = m["samples"]
    assert min(s["format_ok"]) == 1 and min(s["legal"]) == 1 and min(s["topo_any"]) == 1 and min(s["passed"]) == 1
    assert min(s["iou_tc"]) >= 0.99
    for n in names:
        o = m["objects"][n]
        assert min(o["iou_tc"]) >= 0.99 and min(o["d_max"]) >= 0.99 and min(o["f5"]) >= 0.95
        assert o["soft_dice"] >= 0.98 and o["ged"] <= 0.02
    assert m["tv"] == 0.0 and m["js"] == 0.0 and m["min_mode"] == 1.0 and m["coverage"] == 1.0
    assert all(v == 1.0 for v in m["pass_at_k"].values())


def test_raters_own_targets_score_perfectly(nod, riga):
    raters, texts = nod
    m = metrics(texts, raters, S, S)
    assert_perfect(m, ("nodule",))
    assert m["modes"] == 2 and m["K"] == 4 and set(m["pass_at_k"]) == {"1", "2", "4"}
    raters, texts = riga
    m = metrics(texts, raters, 100, 100)
    assert_perfect(m, ("disc", "cup"))
    assert m["modes"] == 1


def test_raters_own_masks_are_exact():
    """Lossless targets (tau = 0 rung) reproduce the masks: every metric is exactly at its perfect value."""
    raters = nodule_raters()
    ma = {"nodule": ExpertCache.build(raters["nodule"]).min_area}
    exact = [build_target({"nodule": r}, ("nodule",), BuildConfig(tau=0.0), ma).text for r in raters["nodule"]]
    o = metrics(exact, raters, S, S)["objects"]["nodule"]
    assert min(o["iou_tc"]) == 1.0 and min(o["d_max"]) == 1.0 and min(o["f5"]) == 1.0
    assert o["soft_dice"] == 1.0 and o["ged"] == pytest.approx(0.0, abs=1e-12)


def test_tv_zero_when_samples_reproduce_the_histogram(nod):
    raters, texts = nod
    m = metrics(texts * 2, raters, S, S)  # every rater twice: the K-sample histogram equals the rater histogram
    assert m["tv"] == 0.0 and m["js"] == 0.0 and m["min_mode"] == 1.0


def test_known_tv_for_a_hand_built_case(nod):
    raters, texts = nod
    m = metrics([texts[0]] * 3 + ["junk"], raters, S, S)  # 3/4 the no-hole answer, 1/4 a format failure
    assert m["tv"] == pytest.approx(0.5)  # |3/4 - 1/2| + |0 - 1/2| + |1/4 - 0|, halved
    assert m["min_mode"] == 0.0 and m["coverage"] == 0.5
    assert m["samples"]["format_ok"] == [1.0, 1.0, 1.0, 0.0] and m["samples"]["error"][3] == "syntax"


def test_failed_samples_score_zero(nod):
    raters, texts = nod
    m = metrics([texts[0], texts[0][:-3], texts[1]], raters, S, S, truncated=[False, False, True])
    s, o = m["samples"], m["objects"]["nodule"]
    assert s["format_ok"] == [1.0, 0.0, 0.0] and s["legal"] == [1.0, 0.0, 0.0] and s["passed"] == [1.0, 0.0, 0.0]
    assert s["error"][1:] == ["incomplete", "truncated"]
    for key in ("iou_tc", "f5", "g_area", "iou_vote", "dice_vote", "d_max"):
        assert o[key][1:] == [0.0, 0.0]
    assert m["pass_at_k"]["1"] == pytest.approx(1 / 3) and m["pass_at_k"]["2"] == pytest.approx(2 / 3)


def test_illegal_samples_give_no_topology_answer(nod):
    raters, texts = nod
    spike = '{"nodule":[{"outer":[[300,300],[700,300],[700,700],[300,700],[300,500],[200,500],[300,500]],"holes":[]}]}'
    m = metrics([texts[0], texts[2], spike, spike], raters, S, S)  # illegal square whose raster is a rater answer
    assert m["samples"]["format_ok"] == [1.0] * 4 and m["samples"]["legal"] == [1.0, 1.0, 0.0, 0.0]
    assert m["samples"]["passed"] == [1.0, 1.0, 0.0, 0.0]
    assert m["tv"] == pytest.approx(0.5) and m["min_mode"] == 1.0 and m["coverage"] == 1.0  # half the mass on "fail"


def test_illegal_answers_score_zero_on_geometry():
    """REVIEW #2 repro: 8 copies of a square with a collinear spike (parsed, illegal) against two near-identical
    raters. Illegal scores 0 on every geometric score, as the reward's legality gate zeroes it, so the Phase 2
    geometry gates read 0; the audit still measures the drawn area."""
    sq = np.zeros((S, S), bool)
    sq[15:35, 15:35] = True
    wide = sq.copy()
    wide[15:35, 35] = True
    raters = {"nodule": np.stack([sq, wide])}
    spike = '{"nodule":[{"outer":[[300,300],[700,300],[700,700],[300,700],[300,500],[200,500],[300,500]],"holes":[]}]}'
    it = {"item_id": "lidc/x", "dataset": "lidc", **metrics([spike] * 8, raters, S, S)}
    s, o = it["samples"], it["objects"]["nodule"]
    assert s["format_ok"] == [1.0] * 8 and s["legal"] == [0.0] * 8 and s["passed"] == [0.0] * 8
    for key in ("iou_tc", "f5", "g_area"):
        assert s[key] == [0.0] * 8, key
    for key in ("iou_tc", "f5", "g_area", "iou_vote", "dice_vote", "d_max"):
        assert o[key] == [0.0] * 8, key
    assert o["soft_dice"] == 0.0 and o["ged"] > 0.9  # the probability map holds no illegal raster
    assert o["a_amax"] == [400 / 420] * 8 and o["components"] == [1.0] * 8  # audit: the parsed raster's area
    d = summarize([it])["datasets"]["lidc"]
    assert d["iou_tc_pass"] == 0.0 and d["pass_rate"] == 0.0
    assert d["objects"]["nodule"]["iou_tc_median"] == 0.0 and d["objects"]["nodule"]["f5_median"] == 0.0
    g = gates.phase2_entry({"views": {v: {"summary": summarize([it])} for v in ("shifted", "rescaled")}})
    vals = g["gates"]["G-geometry"]["values"]
    assert g["gates"]["G-geometry"]["pass"] is False
    assert vals["pass_share/shifted/lidc"] == 0.0 and vals["median/rescaled/lidc/nodule"] == 0.0
    assert hack_audit([it])["datasets"]["lidc"]["rhr"] == 1.0  # every parsed output is a failure


def test_pass_needs_every_object(riga):
    """A RIGA answer with a rater's disc and a misplaced cup (same topology) does not pass: IoU_tc(cup) < 0.5."""
    raters, texts = riga
    bad = text_of({"disc": raters["disc"][0], "cup": disk(30, 50, 6, 100)}, ("disc", "cup"), raters)
    m = metrics([texts[0], bad], raters, 100, 100)
    assert m["samples"]["legal"] == [1.0, 1.0] and m["samples"]["topo_any"] == [1.0, 1.0]
    assert m["objects"]["cup"]["iou_tc"][1] < 0.5 <= m["samples"]["iou_tc"][1]  # the completion mean would pass
    assert m["samples"]["passed"] == [1.0, 0.0]


def test_iou_and_dice_vs_majority_and_dmax(nod):
    raters, _ = nod
    cache = ExpertCache.build(raters["nodule"])
    vote = cache.vote
    m = metrics([text_of({"nodule": vote}, ("nodule",), raters)] * 2, raters, S, S)
    o = m["objects"]["nodule"]
    assert o["iou_vote"][0] >= 0.99 and o["dice_vote"][0] >= 0.99
    ious = [np.sum(vote & r) / np.sum(vote | r) for r in raters["nodule"]]
    assert o["d_max"][0] == pytest.approx(max(ious), abs=0.01)


def test_groups_duplicates_and_signal(nod):
    raters, texts = nod
    m = metrics([texts[0]] * 8, raters, S, S)
    assert m["group"] == 8 and m["groups"]["dup"] == [7 / 8] and m["groups"]["signal"] == [False]
    m = metrics(texts * 2, raters, S, S)  # K = 8: one group of 8 with 4 distinct texts
    assert m["groups"]["dup"] == [0.5]
    m = metrics(texts * 4, raters, S, S)  # K = 16: two groups of 8, as GRPO samples them
    assert m["group"] == 8 and m["groups"]["dup"] == [0.5, 0.5]
    m = metrics(texts[:3], raters, S, S)  # K not a multiple of 8: one group of K
    assert m["group"] == 3 and len(m["groups"]["dup"]) == 1


def test_metrics_is_deterministic_and_json_native(nod):
    import json
    raters, texts = nod
    a, b = metrics(texts, raters, S, S), metrics(texts, raters, S, S)
    assert json.dumps(a, allow_nan=False) == json.dumps(b, allow_nan=False)


def test_metrics_refuses_single_rater_and_bad_inputs(nod):
    raters, texts = nod
    with pytest.raises(ValueError):
        metrics(texts[:1], {"nodule": raters["nodule"][:1]}, S, S)
    with pytest.raises(ValueError):
        metrics(texts, raters, S, S, truncated=[False])


# ---------------------------------------------------------------- hack audit

def item(texts, raters, side, dataset="lidc"):
    return {"item_id": f"{dataset}/x", "dataset": dataset, **metrics(texts, raters, side, side)}


@pytest.fixture(scope="module")
def adversary_texts(nod):
    raters, _ = nod
    cache = ExpertCache.build(raters["nodule"])
    adv = H.adversaries(raters["nodule"], cache, np.random.default_rng(0))
    return {a: text_of({"nodule": m}, ("nodule",), raters) for a, m in adv.items()}


def test_audit_does_not_flag_the_raters_own_targets(nod, riga):
    items = [item(nod[1], nod[0], S), item(riga[1], riga[0], 100, "riga")]
    a = hack_audit(items)
    assert a["fired"] == []
    assert a["alarms"]["components"]["fired"] is None  # no baseline given
    assert a["flags"]["topology_collapse"]["fired"] is None  # K < 64: not evaluated
    assert a["datasets"]["lidc"]["rhr"] == 0.0


EXPECTED = {  # adversary -> the per-output audit flag that must catch it (duplicates aside)
    "dilated_union": {"area_median", "area_tail", "inflation"},
    "hull": {"inflation"},
    "union": {"inflation"},
    "core": {"deflation"},
    "band_speck": {"fragmentation"},
    "holes_stamped": {"topology_collapse"},
    "minority_hole": {"topology_collapse"},
}


def test_audit_flags_every_r3_adversary(nod, adversary_texts):
    raters, _ = nod
    assert {"dilated_union", "hull", "union", "core", "band_speck", "holes_stamped"} <= set(adversary_texts)
    for name, text in adversary_texts.items():
        a = hack_audit([item([text] * 64, raters, S)])
        assert EXPECTED[name] <= set(a["fired"]), (name, a["fired"])
        assert "duplicates" in a["fired"]


def test_component_alarm_against_a_baseline(nod):
    raters, texts = nod
    cache = ExpertCache.build(raters["nodule"])
    speck = H.band_speck(cache, np.random.default_rng(1))
    assert speck is not None
    base = hack_audit([item(texts, raters, S)])
    a = hack_audit([item([text_of({"nodule": speck}, ("nodule",), raters)] * 4, raters, S)], baseline=base)
    assert a["alarms"]["components"]["fired"] is True and "components" in a["fired"]
    assert hack_audit([item(texts, raters, S)], baseline=base)["alarms"]["components"]["fired"] is False


def test_envelope_flags_rise_over_a_baseline(nod):
    """One union per group of five is below the absolute 50 % line but 20 points above a baseline without hedges."""
    raters, texts = nod
    union = text_of({"nodule": ExpertCache.build(raters["nodule"]).union}, ("nodule",), raters)
    mixed = [item([union] + texts, raters, S)]
    assert hack_audit(mixed)["flags"]["inflation"]["fired"] is False
    a = hack_audit(mixed, baseline=hack_audit([item(texts, raters, S)]))
    assert a["flags"]["inflation"]["fired"] is True and a["flags"]["inflation"]["rise"]["lidc/nodule"] == 0.2
    assert a["flags"]["deflation"]["fired"] is False


def test_newly_rewarded_failures(nod):
    """I18: a failure fails the pass criterion or lies outside [A_min / 1.3, 1.3 A_max]. A 1.34x dilation of a rater
    passes (IoU_tc 0.75) but is a failure; the raters' own masks and the union (inside 1.3 A_max) are not. When the
    reward favours the dilation, NRFR > RHR."""
    raters, texts = nod
    big = text_of({"nodule": H.grow(raters["nodule"][1], 2.5)}, ("nodule",), raters)
    union = text_of({"nodule": ExpertCache.build(raters["nodule"]).union}, ("nodule",), raters)
    it = item([big, union] + texts, raters, S)
    assert it["samples"]["passed"] == [1.0] * 6 and 1.3 < it["objects"]["nodule"]["a_amax"][0] < 1.4
    it["samples"]["advantage"] = [1.0, -1.0, 0.5, -0.5, 0.5, -0.5]  # the reward favours the dilation
    a = hack_audit([it])
    d = a["datasets"]["lidc"]
    assert (d["rhr"], d["nrf"], d["nrfr"]) == (1 / 6, 1, 1 / 3) and "newly_rewarded_failures" in a["fired"]
    it["samples"]["advantage"] = [-1.0, 1.0, 0.5, -0.5, 0.5, -0.5]  # the reward favours rater-like outputs
    assert "newly_rewarded_failures" not in hack_audit([it])["fired"]


def test_area_tail_is_relative_to_the_baseline(nod):
    """I17: with a baseline, the tail alarm fires only above the baseline's share + 10 points (honest output reaches
    15 % above 1.3 A_max on small nodules); without one, at the absolute 10 % line."""
    raters, texts = nod
    big = text_of({"nodule": H.grow(raters["nodule"][1], 2.5)}, ("nodule",), raters)  # 1.34 A_max
    base_items = [item([big] + texts * 2, raters, S)]  # 1 of 9 outputs above 1.3 A_max: 11 %
    base = hack_audit(base_items)
    assert base["alarms"]["area_tail"]["fired"] is True  # absolute rule
    assert hack_audit(base_items, baseline=base)["alarms"]["area_tail"]["fired"] is False
    hedged = hack_audit([item([big] * 3 + texts * 2, raters, S)], baseline=base)  # 3 / 11 = 27 %: +16 points
    assert hedged["alarms"]["area_tail"]["fired"] is True


def test_alarms_for_training_evaluations(nod):
    """B4-4: alarms_for audits what a training evaluation holds (per completion: key, text, truncated, W, H, names;
    a cache loader) exactly as hack_audit audits eval.metrics items; the baseline is a previous alarms_for result."""
    from objdel.eval import alarms_for
    raters, texts = nod
    caches = {"nodule": ExpertCache.build(raters["nodule"])}
    union = text_of({"nodule": caches["nodule"].union}, ("nodule",), raters)
    base = alarms_for([("lidc/a", t, False, S, S, ["nodule"]) for t in texts * 2], {"lidc/a": caches}.__getitem__)
    assert base == hack_audit([{"item_id": "lidc/a", "dataset": "lidc", **metrics(texts * 2, raters, S, S)}])
    assert base["fired"] == []
    a = alarms_for([("lidc/a", union, False, S, S, ["nodule"])] * 8, {"lidc/a": caches}.__getitem__, baseline=base)
    assert {"duplicates", "inflation"} <= set(a["fired"])


# ---------------------------------------------------------------- summary and Phase 2 entry gates

def good_report():
    """A one-view (shifted) report as eval_run writes it; the full §5.1 schema is tested in test_p2_eval.py."""
    def ds(objs, iou):
        return {"items": 10, "samples": 640, "format": 0.99, "legal": 0.95, "topo_any_unanimous": 0.95,
                "iou_tc_pass": 0.8, "group": 8, "dup_rate": 0.05, "signal_share": 0.6,
                "objects": {n: {"iou_tc_median": iou[n]} for n in objs}}
    return {"views": {"shifted": {"summary": {"datasets": {
        "riga": ds(("disc", "cup"), {"disc": 0.9, "cup": 0.8}), "lidc": ds(("nodule",), {"nodule": 0.7})}}}}}


def test_phase2_gates_on_one_view_pass_what_it_measures_and_fail_the_rest():
    g = gates.phase2_entry(good_report())
    assert set(g["gates"]) == {"G-format", "G-geometry", "G-prior-b", "G-precision", "G-topology-a", "G-signal",
                               "G-vision"}
    assert all(g["gates"][k]["pass"] for k in ("G-format", "G-signal"))
    assert g["pass"] is False  # the rescaled view, the prior, the probes and G-vision are not in this report


@pytest.mark.parametrize("path,value,gate", [
    (("riga", "format"), 0.97, "G-format"), (("lidc", "legal"), 0.89, "G-format"),
    (("lidc", "dup_rate"), 0.11, "G-signal"), (("riga", "signal_share"), 0.49, "G-signal"),
    (("lidc", "group"), 4, "G-signal"),
])
def test_phase2_gates_fail_cases(path, value, gate):
    r = good_report()
    d = r["views"]["shifted"]["summary"]["datasets"]
    for p in path[:-1]:
        d = d[p]
    d[path[-1]] = value
    g = gates.phase2_entry(r)
    assert g["gates"][gate]["pass"] is False
    assert all(g["gates"][k]["pass"] for k in ("G-format", "G-signal") if k != gate)


def test_phase2_gates_per_dataset():
    """B4-5: gates are decided per dataset; a model dropped from LIDC (plan §5.1 precision probe) is judged on RIGA."""
    r = good_report()
    r["views"]["shifted"]["summary"]["datasets"]["lidc"]["format"] = 0.5
    g = gates.phase2_entry(r, datasets=("riga",))
    assert g["datasets"] == {"riga": False} and g["gates"]["G-format"]["pass"] is True
    assert "format/shifted/lidc" not in g["gates"]["G-format"]["values"]
    assert gates.phase2_entry({"views": {"shifted": {"summary": {"datasets": {}}}}})["pass"] is False


def test_summary_of_the_raters_own_targets(nod, riga):
    s = summarize([item(nod[1], nod[0], S), item(riga[1], riga[0], 100, "riga")])
    lidc, rg = s["datasets"]["lidc"], s["datasets"]["riga"]
    assert lidc["format"] == 1.0 and lidc["legal"] == 1.0 and lidc["pass_rate"] == 1.0 and lidc["tv"] == 0.0
    assert lidc["disagreement"]["items"] == 1 and lidc["disagreement"]["min_mode"] == 1.0
    assert rg["topo_any_unanimous"] == 1.0 and rg["objects"]["cup"]["iou_tc_median"] >= 0.99
    assert rg["disagreement"]["items"] == 0 and rg["disagreement"]["tv"] is None
    assert lidc["pass_at_k"] == {"1": 1.0, "2": 1.0, "4": 1.0}


# ---------------------------------------------------------------- eval_demo

def test_eval_demo_builds_targets_with_the_prepare_size_rule(riga, monkeypatch):
    """REVIEW #11: the demo's rater, honest and adversary targets use targets.min_areas (with the RIGA rim), the rule
    build_targets and train.dataset.prepare use, so the "raters" model is the rater's own SFT target."""
    from experiments.phase1 import eval_demo
    from objdel.targets import min_areas
    raters, texts = riga
    seen = []

    def spy(masks, names, cfg=None, min_area=None):
        seen.append(min_area)
        return build_target(masks, names, cfg, min_area)

    monkeypatch.setattr(eval_demo, "build_target", spy)
    _, out = eval_demo._texts(("riga/x", {"dataset": "riga", "masks": raters}))
    assert seen and all(m == min_areas(raters) for m in seen) and "rim" in seen[0]
    assert out["raters"] == texts


def test_precision_and_recall_vs_the_majority_outline(nod):
    """Precision = answer pixels inside the doctors' majority outline / answer pixels; recall = majority pixels the
    answer covers / majority pixels. An answer around the majority has recall 1, one inside it precision 1."""
    raters, _ = nod
    vote = ExpertCache.build(raters["nodule"]).vote
    big, small = disk(25, 25, 18), disk(25, 32, 3)  # small: inside every rater, clear of the majority's hole
    m = metrics([text_of({"nodule": big}, ("nodule",), raters), text_of({"nodule": small}, ("nodule",), raters),
                 "not json"], raters, S, S)  # fmt: skip
    o = m["objects"]["nodule"]
    assert o["recall_vote"][0] >= 0.99 and o["precision_vote"][0] == pytest.approx(vote.sum() / big.sum(), rel=0.05)
    assert o["precision_vote"][1] >= 0.99 and o["recall_vote"][1] == pytest.approx(small.sum() / vote.sum(), rel=0.15)
    assert o["precision_vote"][2] == o["recall_vote"][2] == 0.0  # a failed answer scores 0, as the other terms


# ---------------------------------------------------------------- A4 / A5: answer IoU matrix and boundary IoU

F = 128
BOWTIE = '{"nodule":[{"outer":[[300,300],[700,700],[700,300],[300,700]],"holes":[]}]}'  # self-intersecting
GOLDEN_I02 = Path(__file__).parent / "goldens" / "a_i02_metrics_before.json"


def i02_item(k):
    """I02's item on a 128 px frame: K = 4 answers (rater 0's copy, it shifted 3 px, a self-intersecting ring,
    unparseable) against k = 2 or 3 raters."""
    r = {"nodule": np.stack([disk(64, 64, 20, F), disk(64, 67, 22, F), disk(62, 64, 18, F)][:k])}
    texts = [text_of({"nodule": m}, ("nodule",), r) for m in (r["nodule"][0], np.roll(r["nodule"][0], 3, 1))]
    return r, [*texts, BOWTIE, "not json"]


@pytest.mark.parametrize("k", [2, 3])
def test_metrics_additions_keep_the_old_keys(k):
    """I02: sample_iou is K x K, symmetric, in [0, 1], 1 on the diagonal; the keys and values before A4 / A5 are
    unchanged (golden dumped from commit df4c9f7 by a one-off script; hm_iou left out: A5 moves it to the LCM
    protocol)."""
    raters, texts = i02_item(k)
    m = metrics(texts, raters, F, F)
    o = m["objects"]["nodule"]
    s = np.array(o["sample_iou"])
    assert s.shape == (4, 4) and np.array_equal(s, s.T) and ((0 <= s) & (s <= 1)).all()
    assert np.diag(s).tolist() == [1.0] * 4
    old = json.loads(GOLDEN_I02.read_text())[str(k)]
    del old["objects"]["nodule"]["hm_iou"], o["hm_iou"]
    for key, v in old.items():
        if key == "objects":
            assert v["nodule"].items() <= o.items()
        else:
            assert m[key] == v, key


def test_illegal_failed_and_empty_answers_in_sample_iou(nod):
    """I03 (A4 part): two failed or illegal answers have IoU 1 with each other; no NaN."""
    raters, texts = nod
    o = metrics([texts[0], "not json", BOWTIE, ""], raters, S, S)
    s = np.array(o["objects"]["nodule"]["sample_iou"])
    assert o["samples"]["legal"][1:] == [0.0] * 3
    assert s[1, 2] == 1.0 and s[0, 1] == 0.0 and np.isfinite(s).all()


def test_sample_iou_is_averaged_over_the_objects(riga):
    """U16: the item's answer matrix is (disc + cup) / 2, exactly."""
    raters, texts = riga
    o = metrics(texts * 2, raters, 100, 100)["objects"]
    A, B = np.array(o["disc"]["sample_iou"]), np.array(o["cup"]["sample_iou"])
    assert np.array_equal(select.sample_iou(o), (A + B) / 2)

