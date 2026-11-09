"""Break-it round 2 (data and eval): edge-case tests; each states the correct behaviour."""
import numpy as np
import pytest

from objdel.data import coverage as Cov
from objdel.data import quality as Q
from objdel.data import windows as Wn
from objdel.eval import select


def _fe():
    from experiments.phase2 import final_eval
    return final_eval


def _report(prefixes):
    return {"leakage": {"patients": 0, "series": 0, "items": 0, "slice_sha": 0, "crop_sha": 0, "near": []},
            "rows": {"not_train": [], "val_test_patient": [], "view_image": [], "excluded": [],
                     "exclusions_stale": False},
            "windows": {"windows": 10, "inside": 10, "clear": 10, "fallback": 0, "kept_side": 0, "mismatch": [],
                        "padded": 0},
            "schedule": {"prefixes": prefixes, "mix": {"tv": 0.0}, "s0_share": 0.9, "count_mismatch": [],
                         "unknown": 0, "disagreement_share": 0.1, "val": {"val": {"coverage": {"tv": 0.0}}}},
            "tokens": {"max": 100, "limit": 2000, "errors": 0}, "agreement": {"share_lt_0.3": 0.0}, "frozen": {},
            "z": []}


@pytest.mark.xfail(strict=True, reason="BR2-data-1: final_eval.summarize raises ValueError for the whole summary when "
                   "the found items (K-mean IoU_tc >= MISS) all come from one patient or nodule; signal.draws "
                   "refuses one cluster (BR1 fix) and iou_found does not catch it; the CI should be None")
def test_summarize_one_found_patient_does_not_crash():
    fe = _fe()
    items = []
    for i, (p, iou) in enumerate([("a", 0.05), ("b", 0.05), ("c", 0.9), ("c", 0.8)]):
        d = {k: 0.5 for k in fe.METRICS}
        d.update(item_id=f"lidc/{p}/n00{i}/z1", iou_tc=iou, sel_mean_of_k=iou, sel_medoid=iou)
        items.append(d)
    groups = {it["item_id"]: it["item_id"].split("/")[1] for it in items}
    out = fe.summarize(items, groups, 50, 0)
    assert out["iou_found"]["ci95_patient"] is None


@pytest.mark.xfail(strict=True, reason="BR2-data-2: quality.gates takes min/max of the prefix coverage values; a NaN "
                   "that is not first is dropped (min([0.9, nan]) = 0.9), so DG7, DG8 and DG9 pass on NaN prefixes")
def test_gates_dg7_dg9_fail_on_nan_prefix():
    nan = float("nan")
    good = {"worst_ratio": 0.9, "tv": 0.05, "ks": 0.02}
    gs = {g.name: g for g in Q.gates(_report([good, {"worst_ratio": nan, "tv": nan, "ks": nan}]))}
    assert not gs["DG7"].ok and not gs["DG8"].ok and not gs["DG9"].ok


@pytest.mark.xfail(strict=True, reason="BR2-data-3: mix_coverage with no target offsets returns ks = NaN (scipy "
                   "SmallSampleWarning); the docstring says a failure is never NaN")
def test_mix_coverage_empty_target_offsets_is_not_nan():
    w = Wn.Window("lidc/a/n1/z1", "a", "lidc/a/n1", "shifted", 0, 0, 0, 64, 64, 0.1, 20, False)
    out = Cov.mix_coverage([w], np.full((4, 4), 1 / 16), [])
    assert not np.isnan(out["ks"])


@pytest.mark.xfail(strict=True, reason="BR2-data-4: coverage.gaps against an empty exam distribution returns [] (no "
                   "gap): distribution([]) gives empty shares instead of raising, as cell_shares does for no windows")
def test_empty_exam_is_not_full_coverage():
    row = {k: 1 for k in Cov.KEYS}
    train = Cov.distribution([row])
    try:
        exam = Cov.distribution([])
    except ValueError:
        return
    assert Cov.gaps(train, exam) != []


@pytest.mark.xfail(strict=True, reason="BR2-data-5: medoid with a NaN in the IoU matrix returns the first NaN-scored "
                   "answer (np.argmax picks NaN) with no error; a non-finite matrix should be refused")
def test_medoid_nan_matrix_is_refused():
    M = np.array([[1, np.nan, 0.9, 0.9], [np.nan, 1, 0.9, 0.9], [0.9, 0.9, 1, 0.95], [0.9, 0.9, 0.95, 1]])
    with pytest.raises(ValueError):
        select.medoid(M, [True] * 4)
