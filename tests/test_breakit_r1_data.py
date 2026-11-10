"""Break-it round 1 (data and eval): edge-case tests; each states the correct behaviour."""
import json

import numpy as np
import pytest
from test_data_windows import lidc_item

from objdel.data import coverage as Cov
from objdel.data import splits as Sp
from objdel.data import windows as Wn
from objdel.verify import ceiling as C
from objdel.verify import gates as G
from objdel.verify import select
from objdel.verify import stats as S


def _gate_items(n_patients):
    items = []
    for i in range(4):
        m = (0.5, 0.6, 0.4, 0.7)[i]
        d = {k: 0.5 for k in FE().METRICS}
        d.update(item_id=f"lidc/p{i % n_patients}/n00{i}/z1", iou_tc=m, sel_mean_of_k=m,
                 sel_medoid=m + 0.05 + 0.03 * (i % 2))
        items.append(d)
    return items, {it["item_id"]: it["item_id"].split("/")[1] for it in items}


def FE():
    from objdel.policy import final_eval
    return final_eval


# --- bootstrap with one cluster -------------------------------------------------------------------------------------
def test_bootstrap_one_cluster_is_refused_or_not_degenerate():
    ids = ["lidc/a/n1/z1", "lidc/a/n1/z2", "lidc/a/n1/z3"]
    x = np.array([0.1, 0.5, 0.9])
    try:
        mean, lo, hi = S.ci(x, S.draws(ids, 50, 0))
    except ValueError:
        return  # refusing is right (cluster_t refuses one patient too)
    assert lo < mean < hi


def test_gate_a_needs_more_than_one_patient():
    items, groups = _gate_items(1)
    try:
        s = FE().summarize(items, groups, 200, 0)
    except ValueError:
        return
    assert FE().gate_a(s, None)["change_default"] is False


# --- medoid ties ----------------------------------------------------------------------------------------------------
def test_medoid_identical_answers_pick_lowest_index():
    rng = np.random.default_rng(0)
    for _ in range(20000):
        a, b, c, d = (round(float(v), 3) for v in rng.integers(1, 999, 4) / 1000)
        K = 6
        M = np.eye(K)
        others = {0: a, 2: b, 4: c, 5: d}
        M[1, 3] = M[3, 1] = 1.0
        for j, v in others.items():
            M[1, j] = M[j, 1] = M[3, j] = M[j, 3] = v
        for i in (0, 2, 4, 5):
            for j in (0, 2, 4, 5):
                if i < j:
                    M[i, j] = M[j, i] = 0.001
        if select.medoid(M, [True] * K) == 3:
            M2 = M.copy()
            assert select.medoid(M2, [True] * K) == 1  # answers 1 and 3 are the same mask: lowest index
            return
    pytest.skip("no float tie found")


# --- pool and schedule ----------------------------------------------------------------------------------------------
def test_pool_duplicate_items_give_unique_window_ids():
    it = lidc_item("lidc/1.2.3/n000/z001", [(100, 100, 8), (100, 101, 8)])
    p = Wn.pool([it, it], draws=2, seed=0)
    ids = [w.window_id for w in p.windows]
    assert len(ids) == len(set(ids))


def test_cell_of_nan_raises():
    with pytest.raises(ValueError):
        Wn.cell(float("nan"), 0.0)


# --- ceiling --------------------------------------------------------------------------------------------------------
def test_ci_with_nan_item_is_not_half_finite():
    x = np.array([0.2, np.nan, 0.6, 0.4])
    ids = [f"lidc/p{i}/n000/z1" for i in range(4)]
    mean, lo, hi = S.ci(x, S.draws(ids, 200, 0))
    assert np.isnan(mean) == np.isnan(lo) == np.isnan(hi)


def test_ceiling_found_miss_partition_all_rows():
    def row(i, m):
        return {"item_id": f"lidc/s{i}/n000/z001", "k": 3, "n_voters": 3, "single": False, "doctor": .5,
                "legacy_dmax": [.5, .5, .5], "legacy_loo_dmax_item": .5, "model": m, "model_all_k": m, "diff": m - .5}
    rows = [row(0, .6), row(1, .05), row(2, float("nan")), row(3, .7), row(4, .2)]
    groups = {r["item_id"]: f"p{i % 2}" for i, r in enumerate(rows)}
    try:
        s = C.summarize(rows, groups, boot=20)
    except ValueError:
        return
    assert s["found"]["n"] + s["miss"]["n"] == s["all"]["n"]


def test_val_half_table_of_nothing_raises():
    with pytest.raises(ValueError):
        Sp.val_half_table([])


# --- coverage -------------------------------------------------------------------------------------------------------
def test_props_of_empty_object_is_none():
    assert Cov.props(json.dumps({}), 100, "op") is None


# --- gates ----------------------------------------------------------------------------------------------------------
def test_unknown_dataset_does_not_pass_geometry_without_structures():
    summ = {"format": 1.0, "legal": 1.0, "iou_tc_pass": 1.0, "group": 8, "signal_share": 1.0, "dup_rate": 0.0,
            "objects": {}}
    rep = {"views": {v: {"summary": {"datasets": {"foo": summ}}} for v in ("shifted", "rescaled")},
           "vision": {"pass": True, "ratio": 1, "tv": 1}}
    out = G.phase2_entry(rep, ("foo",))
    assert out["gates"]["G-geometry"]["pass"] is False
