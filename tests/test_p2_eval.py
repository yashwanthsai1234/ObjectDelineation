"""Phase 2 M4a + M9 (eval): HM-IoU and NSD in eval.metrics, truncation / token stats in summarize, G-vision
(eval.vision) on a tiny random VLM, the revised phase2_entry (plan §5.1) and eval_run's single-rater and prepared-cache
paths (I45)."""
import json
import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from oracles.eval_oracle import nsd
from scipy import ndimage
from scipy.optimize import linear_sum_assignment, linprog

from objdel.data.targets.build import build_target, min_areas
from objdel.verify import gates
from objdel.verify.metrics import NSD_TAU, hm_iou, metrics, min_mode_retention, summarize, tv_js
from objdel.verify.reward.cache import ExpertCache

S = 50


def disk(cy, cx, r, side=S):
    yy, xx = np.ogrid[:side, :side]
    return (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r


def raters_and_targets():
    raters = {"nodule": np.stack([disk(25, 25, 12), disk(25, 27, 14), disk(24, 25, 11) & ~disk(24, 25, 3)])}
    ma = {"nodule": ExpertCache.build(raters["nodule"]).min_area}
    texts = [build_target({"nodule": m}, ("nodule",), min_area=ma).text for m in raters["nodule"]]
    return raters, texts


# ---------------------------------------------------------------- HM-IoU

def transport_hm_iou(iou: np.ndarray) -> float:
    """The LCM protocol as a transport LP (math-stats §4): sample i sends L / K units, rater j takes L / k units,
    maximise the matched IoU; the transport polytope has integer vertices, so the LP optimum / L is the duplicated
    sets' Hungarian optimum."""
    K, k = iou.shape
    L = math.lcm(K, k)
    rows = np.kron(np.eye(K), np.ones(k))  # sum_j x_ij = L / K
    cols = np.kron(np.ones(K), np.eye(k))  # sum_i x_ij = L / k
    r = linprog(-iou.ravel(), A_eq=np.vstack([rows, cols]),
                b_eq=np.r_[np.full(K, L // K), np.full(k, L // k)], bounds=(0, None), method="highs")
    return -r.fun / L


@settings(max_examples=150, deadline=None)
@given(st.integers(1, 5).flatmap(lambda k: st.tuples(st.just(k), st.integers(max(k, 2), 8))), st.integers(0, 2**31))
def test_hm_iou_is_the_lcm_protocol(kK, seed):
    """HPU-Net App. B.1: "duplicate both sets so that their number of elements matches their least common
    multiple" (math-stats §11 item 3)."""
    k, K = kK
    rng = np.random.default_rng(seed)
    iou = rng.random((K, k))
    if seed % 3 == 0:
        iou = np.round(iou, 1)  # ties
    assert hm_iou(iou) == pytest.approx(transport_hm_iou(iou), abs=1e-9)


def test_hm_iou_k3_is_not_the_floor_ceil_tiling():
    """K = 8, k = 3: each sample may split over the raters in thirds. Two samples near rater 0, six near rater 1:
    rater 0 takes both (0.9) and 2/3 of one other (0.6), raters 1 and 2 take 8/3 each (0.8, 0.5): 17/24. The old
    floor / ceil tiling (each rater 2 or 3 whole samples) read 0.725 here."""
    iou = np.array([[.9, .5, .1]] * 2 + [[.6, .8, .5]] * 6)
    assert hm_iou(iou) == pytest.approx((2 * .9 + 2 / 3 * .6 + 8 / 3 * .8 + 8 / 3 * .5) / 8, abs=1e-12)
    assert hm_iou(iou) == pytest.approx(17 / 24, abs=1e-12)


@settings(max_examples=50, deadline=None)
@given(st.sampled_from([(1, 4), (2, 4), (4, 8), (2, 8), (4, 16)]), st.integers(0, 2**31))
def test_hm_iou_is_the_standard_tiled_square_matching_when_k_divides_K(kK, seed):
    """HPU-Net / CAR: raters duplicated K / k times, square Hungarian matching, mean matched IoU."""
    k, K = kK
    iou = np.random.default_rng(seed).random((K, k))
    sq = np.tile(iou, K // k)
    r, c = linear_sum_assignment(sq, maximize=True)
    assert hm_iou(iou) == pytest.approx(sq[r, c].mean(), abs=1e-12)


def test_hm_iou_in_metrics_on_the_raters_own_targets():
    raters, texts = raters_and_targets()
    m = metrics(texts * 2, raters, S, S)  # K = 6, k = 3: each rater matched to its own two copies
    assert m["objects"]["nodule"]["hm_iou"] >= 0.99
    m = metrics([texts[0]] * 6, raters, S, S)  # one mode only: matched pairs include the other raters
    assert m["objects"]["nodule"]["hm_iou"] < 0.9


# ---------------------------------------------------------------- NSD: DeepMind surface-distance reference

def _dm_table():
    """surface_distance.lookup_tables.create_table_neighbour_code_to_contour_length((1, 1)), verbatim."""
    first, second = 1.0, 1.0
    diag = 0.5 * math.sqrt(first ** 2 + second ** 2)
    t = np.zeros([16])
    t[int("00000001", 2)] = diag
    t[int("00000010", 2)] = diag
    t[int("00000011", 2)] = second
    t[int("00000100", 2)] = diag
    t[int("00000101", 2)] = first
    t[int("00000110", 2)] = 2 * diag
    t[int("00000111", 2)] = diag
    t[int("00001000", 2)] = diag
    t[int("00001001", 2)] = 2 * diag
    t[int("00001010", 2)] = first
    t[int("00001011", 2)] = diag
    t[int("00001100", 2)] = second
    t[int("00001101", 2)] = diag
    t[int("00001110", 2)] = diag
    return t


def dm_surface_dice(mask_gt, mask_pred, tol):
    """google-deepmind/surface-distance metrics.compute_surface_distances + compute_surface_dice_at_tolerance (2D,
    spacing (1, 1)), following the library code; nan when both masks are empty."""
    table, kernel = _dm_table(), np.array([[8, 4], [2, 1]])
    both = mask_gt | mask_pred
    if not both.any():
        return float("nan")
    rows, cols = np.flatnonzero(both.any(1)), np.flatnonzero(both.any(0))
    bmin, bmax = np.array([rows[0], cols[0]]), np.array([rows[-1], cols[-1]])

    def crop(m):
        c = np.zeros((bmax - bmin) + 2, np.uint8)
        c[0:-1, 0:-1] = m[bmin[0]:bmax[0] + 1, bmin[1]:bmax[1] + 1]
        return c

    code_gt = ndimage.correlate(crop(mask_gt), kernel, mode="constant", cval=0)
    code_pred = ndimage.correlate(crop(mask_pred), kernel, mode="constant", cval=0)
    b_gt, b_pred = (code_gt != 0) & (code_gt != 15), (code_pred != 0) & (code_pred != 15)
    d_gt = ndimage.distance_transform_edt(~b_gt, sampling=(1, 1)) if b_gt.any() else np.inf * np.ones(b_gt.shape)
    d_pred = ndimage.distance_transform_edt(~b_pred, sampling=(1, 1)) if b_pred.any() else np.inf * np.ones(b_gt.shape)
    dist_gt_to_pred, dist_pred_to_gt = d_pred[b_gt], d_gt[b_pred]
    area_gt, area_pred = table[code_gt][b_gt], table[code_pred][b_pred]
    overlap_gt = np.sum(area_gt[dist_gt_to_pred <= tol])
    overlap_pred = np.sum(area_pred[dist_pred_to_gt <= tol])
    return float((overlap_gt + overlap_pred) / (np.sum(area_gt) + np.sum(area_pred)))


masks = st.integers(0, 2**31).map(lambda s: np.random.default_rng(s))


def blob(rng, side):
    m = ndimage.binary_opening(rng.random((side, side)) < rng.uniform(0.2, 0.7), iterations=int(rng.integers(0, 3)))
    return m if rng.random() > 0.1 else np.zeros_like(m)


@settings(max_examples=200, deadline=None)
@given(masks, st.integers(3, 40), st.sampled_from([0.0, 0.5, 1.0, 2.0, 3.5, 9.0]))
def test_nsd_equals_the_deepmind_reference(rng, side, tau):
    a, b = blob(rng, side), blob(rng, side)
    ref = dm_surface_dice(a, b, tau)
    got = nsd(a, b, tau)
    assert got == (1.0 if math.isnan(ref) else pytest.approx(ref, abs=1e-6))
    assert 0.0 <= got <= 1.0 and nsd(a, a, tau) == 1.0


def test_values_rounding_pushes_out_of_range_are_clamped():
    """S1-15 (Hypothesis falsifying example: NSD -2.2e-16 at tau = 0) and the same rounding class in TV / JS /
    minimum-mode retention (1 + 2.2e-16 on these exact histograms)."""
    rng = np.random.default_rng(703608005)
    assert nsd(blob(rng, 39), blob(rng, 39), 0.0) == 0.0
    f = [0.4035760681074393, 0.06121156431979067, 0.402805720374044, 0.13240664719872622]
    assert tv_js(dict(enumerate(f)), {i + 100: v for i, v in enumerate(f[::-1])}) == (1.0, 1.0)
    g = dict(enumerate([0.1915658702295475, 0.37455625100036816, 0.1851729057940717, 0.2487049729760125]))
    assert min_mode_retention(g, g) == 1.0


def test_nsd_hand_cases():
    a = disk(25, 25, 10)
    assert nsd(a, np.zeros_like(a), 2.0) == 0.0 and nsd(np.zeros_like(a), np.zeros_like(a), 2.0) == 1.0
    shifted = np.roll(a, 3, axis=1)
    assert nsd(a, shifted, 3.0) == 1.0 and nsd(a, shifted, 0.0) < 0.5  # every surfel moves by 3 px or less


def test_nsd_in_metrics_is_one_on_the_raters_own_targets():
    raters, texts = raters_and_targets()
    m = metrics(texts + ["{", texts[0]], raters, S, S)
    v = m["objects"]["nodule"]["nsd"]
    assert v[:3] == [1.0] * 3 and v[3] == 0.0 and v[4] == 1.0  # a failed sample is an empty mask
    riga = {"disc": np.stack([disk(50, 50, r, 100) for r in (30, 32)]),
            "cup": np.stack([disk(50, 51, r, 100) for r in (12, 14)])}
    names = ("disc", "cup")
    ma = min_areas(riga)  # with the rim's rule
    t = [build_target({n: riga[n][r] for n in names}, names, min_area=ma).text for r in range(2)]
    m = metrics(t, riga, 100, 100, names)
    assert all(m["objects"][n]["nsd"] == [1.0, 1.0] for n in names)
    assert NSD_TAU == {"nodule": 2.0, "disc": 9.0, "cup": 16.0}


# ---------------------------------------------------------------- summarize: truncation and tokens, new metrics

def test_summarize_reports_truncation_tokens_hm_iou_and_nsd():
    raters, texts = raters_and_targets()
    a = {"item_id": "lidc/a", "dataset": "lidc",
         **metrics(texts + ["{"], raters, S, S, truncated=[0, 0, 0, 1], tokens=[10, 20, 30, 100])}
    b = {"item_id": "lidc/b", "dataset": "lidc", **metrics(texts * 2, raters, S, S, tokens=[10] * 6)}
    d = summarize([a, b])["datasets"]["lidc"]
    assert d["truncated"] == pytest.approx(0.1)
    assert d["tokens"]["mean"] == pytest.approx(220 / 10) and d["tokens"]["p95"] == pytest.approx(
        np.percentile([10, 20, 30, 100] + [10] * 6, 95))
    o = d["objects"]["nodule"]
    assert o["hm_iou"] == pytest.approx(np.mean([a["objects"]["nodule"]["hm_iou"], b["objects"]["nodule"]["hm_iou"]]))
    assert o["nsd_mean"] == pytest.approx(0.9) and o["nsd_median"] == 1.0
    c = {"item_id": "lidc/c", "dataset": "lidc", **metrics(texts, raters, S, S)}  # no token counts
    assert summarize([a, c])["datasets"]["lidc"]["tokens"] is None
    with pytest.raises(ValueError, match="token"):
        metrics(texts, raters, S, S, tokens=[1])


# ---------------------------------------------------------------- revised phase2_entry (plan §5.1)

def _ds(objs, iou):
    return {"items": 10, "samples": 640, "format": 0.99, "legal": 0.95, "iou_tc_pass": 0.8, "group": 8,
            "dup_rate": 0.05, "signal_share": 0.6, "truncated": 0.01, "tokens": {"mean": 80.0, "p95": 150.0},
            "objects": {n: {"iou_tc_median": iou[n]} for n in objs}}


def _summary():
    return {"datasets": {"riga": _ds(("disc", "cup"), {"disc": 0.9, "cup": 0.8}),
                         "lidc": _ds(("nodule",), {"nodule": 0.7})}}


def full_report():
    prior = lambda lo: {"riga": {"disc": {"diff": 0.1, "ci": [lo, 0.2]}, "cup": {"diff": 0.1, "ci": [lo, 0.2]}},  # noqa: E731
                        "lidc": {"nodule": {"diff": 0.1, "ci": [lo, 0.2]}}}
    centred = {"riga": {"disc": {"win_share": 0.7}, "cup": {"win_share": 0.65}},
               "lidc": {"nodule": {"win_share": 0.55}}}
    return {"views": {v: {"summary": _summary()} for v in ("centred", "shifted", "rescaled")},
            "prior": {"shifted": prior(0.01), "rescaled": prior(0.02), "centred": centred},
            "probes": {"precision": {"pass": True, "legal": 0.97}, "topology": {"pass": True, "hole_recall": 0.2}},
            "vision": {"ratio": 0.9, "tv": 0.01, "pass": True}, "g_probe": {"verdict": "PASS"},
            "signal": {"report": "only"}}


GATES = {"G-format", "G-geometry", "G-prior-b", "G-precision", "G-topology-a", "G-signal", "G-vision"}


def test_phase2_entry_passes_a_full_report_and_reports_the_report_only_values():
    g = gates.phase2_entry(full_report())
    assert g["pass"] is True and g["datasets"] == {"lidc": True, "riga": True} and set(g["gates"]) == GATES
    assert g["warnings"] == []
    r = g["reported"]
    assert r["G-prior-a"]["lidc/nodule"] == 0.55 and r["G-topology-b"] == 0.2 and r["G-probe"] == "PASS"
    assert r["truncated"]["shifted/lidc"] == 0.01 and r["tokens_p95"]["shifted/riga"] == 150.0
    assert {"median/shifted/lidc/nodule", "median/rescaled/riga/cup", "pass_share/rescaled/lidc"} <= set(
        g["gates"]["G-geometry"]["values"])
    assert g["passed"] == dict.fromkeys(GATES | {"G-probe"}, True)


def _set(r, path, value):
    d = r
    for p in path[:-1]:
        d = d[p]
    d[path[-1]] = value


@pytest.mark.parametrize("path,value,gate", [
    (("views", "shifted", "summary", "datasets", "riga", "format"), 0.97, "G-format"),
    (("views", "shifted", "summary", "datasets", "lidc", "legal"), 0.89, "G-format"),
    (("views", "shifted", "summary", "datasets", "riga", "objects", "cup", "iou_tc_median"), 0.74, "G-geometry"),
    (("views", "rescaled", "summary", "datasets", "lidc", "objects", "nodule", "iou_tc_median"), 0.59,
     "G-geometry"),
    (("views", "rescaled", "summary", "datasets", "riga", "iou_tc_pass"), 0.69, "G-geometry"),
    (("prior", "shifted", "lidc", "nodule", "ci"), [0.0, 0.1], "G-prior-b"),
    (("prior", "rescaled", "riga", "cup", "ci"), [-0.01, 0.1], "G-prior-b"),
    (("probes", "precision", "pass"), False, "G-precision"),
    (("probes", "topology", "pass"), False, "G-topology-a"),
    (("views", "shifted", "summary", "datasets", "lidc", "dup_rate"), 0.11, "G-signal"),
    (("views", "shifted", "summary", "datasets", "riga", "signal_share"), 0.49, "G-signal"),
    (("vision", "pass"), False, "G-vision"),
])
def test_phase2_entry_fail_cases(path, value, gate):
    r = full_report()
    _set(r, path, value)
    g = gates.phase2_entry(r)
    assert g["pass"] is False and g["gates"][gate]["pass"] is False
    assert all(v["pass"] for k, v in g["gates"].items() if k != gate)


@pytest.mark.parametrize("path,value", [
    (("prior", "centred", "lidc", "nodule", "win_share"), 0.1),  # G-prior (a): report-only (revision 7)
    (("probes", "topology", "hole_recall"), 0.0),  # G-topology (b): report-only, decided by C2
    (("views", "centred", "summary", "datasets", "lidc", "format"), 0.1),  # centred view gates nothing
    (("views", "shifted", "summary", "datasets", "lidc", "truncated"), 0.9),  # reported beside G-format
])
def test_report_only_values_never_fail(path, value):
    r = full_report()
    _set(r, path, value)
    assert gates.phase2_entry(r)["pass"] is True


def test_topo_any_unanimous_is_not_read():
    g = gates.phase2_entry(full_report())  # summarize no longer gives it (X2-5); no gate asks for it
    assert g["pass"] is True and "G-topology-a" in g["gates"]
    assert not any("unanimous" in k for gate in g["gates"].values() for k in gate["values"])


def test_g_probe_not_pass_is_a_warning_not_a_failure():
    for verdict in ("FAIL", "INCONCLUSIVE", None):
        r = full_report()
        r["g_probe"] = {"verdict": verdict} if verdict else None
        g = gates.phase2_entry(r)
        assert g["pass"] is True and any("G-probe" in w for w in g["warnings"]) and g["passed"]["G-probe"] is False


def test_precision_capacity_run_rescues_a_zero_shot_failure_with_a_warning():
    r = full_report()
    r["probes"]["precision"] = {"pass": False, "capacity": {"pass": True}}
    g = gates.phase2_entry(r)
    assert g["pass"] is True and any("capacity" in w for w in g["warnings"])
    r["probes"]["precision"]["capacity"]["pass"] = False
    assert gates.phase2_entry(r)["gates"]["G-precision"]["pass"] is False


def test_a_model_dropped_from_lidc_is_judged_on_riga_alone():
    """Plan §5.1 G-precision: a capacity failure drops the rung from LIDC; RIGA gates still decide."""
    r = full_report()
    r["probes"]["precision"] = {"pass": False, "capacity": {"pass": False}}
    r["probes"]["topology"] = {"pass": False}
    r["views"]["shifted"]["summary"]["datasets"]["lidc"]["legal"] = 0.5
    r["prior"]["shifted"]["lidc"]["nodule"]["ci"] = [-0.2, -0.1]
    g = gates.phase2_entry(r)
    assert g["pass"] is False and g["datasets"] == {"lidc": False, "riga": True}
    g = gates.phase2_entry(r, datasets=("riga",))
    assert g["pass"] is True and g["datasets"] == {"riga": True}
    assert not any("lidc" in k for gate in g["gates"].values() for k in gate["values"])
    assert g["gates"]["G-precision"]["values"] == {} and g["gates"]["G-precision"]["pass"] is True


@pytest.mark.parametrize("drop", [("views", "rescaled"), ("prior", "rescaled"), ("prior", "shifted"),
                                  ("probes", "precision"), ("probes", "topology"), ("vision",), ("views", "shifted")])
def test_a_missing_part_fails_its_gate(drop):
    r = full_report()
    d = r
    for p in drop[:-1]:
        d = d[p]
    del d[drop[-1]]
    g = gates.phase2_entry(r)
    assert g["pass"] is False
    assert gates.phase2_entry({})["pass"] is False  # nothing measured


def test_a_requested_dataset_or_structure_missing_fails():
    r = full_report()
    del r["views"]["rescaled"]["summary"]["datasets"]["riga"]
    assert gates.phase2_entry(r, datasets=("riga", "lidc"))["datasets"]["riga"] is False
    r = full_report()
    del r["prior"]["shifted"]["riga"]["cup"]
    assert gates.phase2_entry(r)["gates"]["G-prior-b"]["pass"] is False


def test_prior_a_gates_only_when_a_threshold_is_set():
    r = full_report()  # centred win shares: disc 0.7, cup 0.65, nodule 0.55
    assert "G-prior-a" not in gates.phase2_entry(r)["passed"]
    rep = gates.phase2_entry(r, prior_a_threshold=0.55)
    assert rep["pass"] and rep["passed"]["G-prior-a"]
    assert rep["gates"]["G-prior-a"]["values"]["win_share/centred/lidc/nodule"] == 0.55
    rep = gates.phase2_entry(r, prior_a_threshold=0.6)
    assert not rep["pass"] and not rep["passed"]["G-prior-a"] and rep["datasets"] == {"riga": True, "lidc": False}
    del r["prior"]["centred"]["riga"]["cup"]  # missing: fails
    assert not gates.phase2_entry(r, prior_a_threshold=0.1)["passed"]["G-prior-a"]


def test_eligible_is_g_format_g_signal_g_vision_of_phase2_entry():
    """select_ckpt's R-B eligibility (prereg.RB["eligible"]) from the gate rules, on one shifted-view summary."""
    s, ok = _summary(), {"pass": True}
    assert gates.eligible(s, ok) and gates.eligible(s, ok, "riga")
    assert not gates.eligible(s, {}) and not gates.eligible(s, None)  # unmeasured G-vision (BK2-5)
    assert not gates.eligible(s, {"pass": False})
    for key, bad in (("group", 12), ("format", 0.97), ("legal", 0.89), ("dup_rate", 0.11), ("signal_share", 0.49)):
        t = _summary()
        t["datasets"]["lidc"][key] = bad
        assert not gates.eligible(t, ok), key  # BK2-6: K != 8 is not eligible
    t = _summary()
    t["datasets"]["lidc"]["objects"]["nodule"]["iou_tc_median"] = 0.0  # geometry is not an eligibility gate
    assert gates.eligible(t, ok)
    assert not gates.eligible({"datasets": {}}, ok)


def test_the_report_from_summarize_fits_the_schema():
    """A real summarize() output in every view slot is read without KeyError."""
    raters, texts = raters_and_targets()
    s = summarize([{"item_id": "lidc/a", "dataset": "lidc", **metrics(texts * 8, raters, S, S)}])
    r = full_report()
    r["views"] = {v: {"summary": s} for v in ("shifted", "rescaled")}
    g = gates.phase2_entry(r, datasets=("lidc",))
    assert g["gates"]["G-format"]["pass"] is True and g["gates"]["G-geometry"]["pass"] is True
    json.dumps(g, allow_nan=False)


# ---------------------------------------------------------------- M9: eval_run (I45)

@pytest.fixture(scope="module")
def stores(tmp_path_factory, synthetic_items):
    from objdel.data import store
    root = tmp_path_factory.mktemp("evalrun")
    store.write(synthetic_items(("val",)), root / "store")
    return root


def _rows(frames):
    out = []
    for item_id, fr in sorted(frames.items()):
        names = tuple(fr["masks"])
        k = len(fr["masks"][names[0]])
        ma = min_areas(fr["masks"])
        texts = [build_target({n: fr["masks"][n][r] for n in names}, names, min_area=ma).text for r in range(k)]
        out.append({"item_id": item_id, "samples": (texts * 2)[:max(2, k)]})
    return out


def test_eval_run_skips_and_counts_single_rater_items(stores):
    """I45: a single-rater item (synthetic z003) no longer stops eval_run; it is skipped and counted."""
    from experiments.phase1 import eval_run
    frames = eval_run.load_frames(stores / "store", "val")
    single = [i for i, f in frames.items() if len(next(iter(f["masks"].values()))) < 2]
    assert single == ["lidc/1.2.1/n000/z003"]
    rows = _rows(frames)
    results = eval_run.evaluate(rows, frames)
    assert {r["item_id"] for r in results} == {r["item_id"] for r in rows} - set(single)
    report = eval_run.build_report(results, skipped=eval_run.skipped(rows, frames))
    assert report["skipped"] == {"single_rater": 1, "ids": single}
    assert "skipped" in eval_run.markdown(report, "t")


def test_eval_run_scores_with_the_prepared_caches(stores):
    """I45 / D25: with a prepared root, metrics get the prepared (digest-checked) caches and equal the store path."""
    pytest.importorskip("datasets")
    pytest.importorskip("cv2")
    from experiments.phase1 import eval_run
    from objdel.data import prepare as data_prepare
    from objdel.train import dataset as D
    from objdel.verify import store as verify_store
    data_prepare.prepare(stores / "store", stores / "prep", splits=("val",), procs=1, prompts=D.PROMPTS)
    frames = eval_run.load_prepared(stores / "prep", "val")
    assert frames and all(f["prepared"] == str((stores / "prep").resolve()) for f in frames.values())
    ref = eval_run.load_frames(stores / "store", "val", list(frames))
    rows = _rows(ref)
    got, want = eval_run.evaluate(rows, frames), eval_run.evaluate(rows, ref)
    assert json.dumps(got) == json.dumps(want)
    seen = []

    def loader(key):  # any key -> {object: ExpertCache} loader, e.g. serve.loader(root) (coordinator request)
        seen.append(key)
        return verify_store.cache_loader(stores / "prep")(key)

    assert json.dumps(eval_run.evaluate(rows, frames, caches=loader)) == json.dumps(want)
    assert sorted(seen) == sorted(r["item_id"] for r in rows if frames[r["item_id"]]["k"] >= 2)
    with pytest.raises(KeyError, match="not in the val"):
        eval_run.load_prepared(stores / "prep", "val", ["lidc/9.9.9/n000/z001"])


# ---------------------------------------------------------------- G-vision on a tiny random VLM

@pytest.fixture(scope="module")
def vlm():
    """A tiny random Llava (CLIP vision tower, 2-layer projector, 1-layer Llama) and its processor, built offline."""
    torch = pytest.importorskip("torch")
    pytest.importorskip("transformers")
    from tokenizers import Regex, Tokenizer, decoders, models, pre_tokenizers
    from transformers import (
        CLIPImageProcessor,
        CLIPVisionConfig,
        LlamaConfig,
        LlavaConfig,
        LlavaForConditionalGeneration,
        LlavaProcessor,
        PreTrainedTokenizerFast,
    )
    chars = sorted(set('{}[]":,0123456789abcdefghijklmnopqrstuvwxyz .'))
    vocab = {w: i for i, w in enumerate(["<pad>", "<eos>", "<image>", *chars])}
    tk = Tokenizer(models.WordLevel(vocab, unk_token="<pad>"))
    tk.pre_tokenizer = pre_tokenizers.Split(Regex("<image>|."), behavior="isolated")
    tk.decoder = decoders.Fuse()
    tok = PreTrainedTokenizerFast(tokenizer_object=tk, pad_token="<pad>", eos_token="<eos>")
    tok.add_special_tokens({"additional_special_tokens": ["<image>"]})
    template = ("{% for m in messages %}{% for c in m['content'] %}{% if c['type'] == 'image' %}<image>{% else %}"
                "{{ c['text'] }}{% endif %}{% endfor %}{% endfor %}{% if add_generation_prompt %}:{% endif %}")
    proc = LlavaProcessor(image_processor=CLIPImageProcessor(size={"shortest_edge": 16},
                                                             crop_size={"height": 16, "width": 16}),
                          tokenizer=tok, patch_size=8, vision_feature_select_strategy="full", chat_template=template,
                          num_additional_image_tokens=1)

    def build(seed=0):
        torch.manual_seed(seed)
        cfg = LlavaConfig(
            vision_config=CLIPVisionConfig(hidden_size=16, intermediate_size=32, num_hidden_layers=1,
                                           num_attention_heads=2, image_size=16, patch_size=8),
            text_config=LlamaConfig(vocab_size=len(tok), hidden_size=16, intermediate_size=32, num_hidden_layers=1,
                                    num_attention_heads=2, num_key_value_heads=2, pad_token_id=0, eos_token_id=1,
                                    bos_token_id=1),
            image_token_index=vocab["<image>"], vision_feature_select_strategy="full", vision_feature_layer=-1)
        return LlavaForConditionalGeneration(cfg).eval()

    return build, proc


def _images(n, seed=0):
    from PIL import Image
    rng = np.random.default_rng(seed)
    return [Image.fromarray(rng.integers(0, 256, (50, 50, 3), dtype=np.uint8)) for _ in range(n)]


PROMPT = "outline the nodule"


def _items(n=4):
    return [{"image": im, "target": '{"nodule":[{"outer":[[1,2],[3,4]]}]}'} for im in _images(n, 1)]


def test_vision_check_of_a_model_against_itself_has_ratio_one(vlm):
    from objdel.eval import vision
    build, proc = vlm
    model = build()
    ref = vision.projector_distance(model, proc, _images(4), PROMPT)
    out = vision.check(model, proc, _images(4), _items(), reference=ref, prompt=PROMPT)
    assert out["ratio"] == pytest.approx(1.0) and 0 < out["distance"] <= 2
    assert 0.0 <= out["tv_black"] <= 1.0 and 0.0 <= out["tv_other"] <= 1.0 and out["tv"] == min(
        out["tv_black"], out["tv_other"])
    json.dumps(out, allow_nan=False)


def test_vision_check_restores_the_training_mode(vlm):
    """S1-9: check runs in eval mode and leaves the model in the mode it found (train mode mid-SFT)."""
    from objdel.eval import vision
    build, proc = vlm
    model = build()
    for mode in (True, False):
        model.train(mode)
        vision.check(model, proc, _images(2), _items(2), reference=1.0, prompt=PROMPT)
        assert model.training is mode and all(m.training is mode for m in model.modules())


def test_vision_check_fails_when_the_projector_is_constant(vlm):
    import torch

    from objdel.eval import vision
    build, proc = vlm
    base = build()
    ref = vision.projector_distance(base, proc, _images(4), PROMPT)
    model = build()
    lin = model.model.multi_modal_projector.linear_2
    with torch.no_grad():
        lin.weight.zero_()
        lin.bias.fill_(0.5)
    out = vision.check(model, proc, _images(4), _items(), reference=ref, prompt=PROMPT)
    assert out["ratio"] == pytest.approx(0.0, abs=1e-6) and out["tv_black"] == out["tv_other"] == 0.0
    assert out["pass"] is False and out["action"] == "abort"


def test_vision_tv_equals_a_numpy_reference(vlm):
    """TV = mean over items and the first T teacher-forced positions of 0.5 * sum |p - q|, p and q softmaxed in
    float64 numpy from logits of a forward pass built here."""
    import torch

    from objdel.eval import vision
    build, proc = vlm
    model = build(3)
    items, T = _items(3), 5

    def probs(image, target):
        msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": PROMPT}]}]
        text = proc.apply_chat_template(msgs, tokenize=False, add_generation_prompt=True)
        b = proc(text=[text], images=[[image]], return_tensors="pt")
        n = b["input_ids"].shape[1]
        tgt = proc.tokenizer(target, add_special_tokens=False, return_tensors="pt")["input_ids"][:, :T]
        b["input_ids"] = torch.cat([b["input_ids"], tgt], 1)
        b["attention_mask"] = torch.ones_like(b["input_ids"])
        with torch.no_grad():
            z = model(**b).logits[0, n - 1:n - 1 + tgt.shape[1]].double().numpy()
        e = np.exp(z - z.max(-1, keepdims=True))
        return e / e.sum(-1, keepdims=True)

    black = [vision.black(it["image"]) for it in items]
    want_black = np.mean([0.5 * np.abs(probs(it["image"], it["target"]) - probs(bl, it["target"])).sum(-1).mean()
                          for it, bl in zip(items, black, strict=True)])
    want_other = np.mean([0.5 * np.abs(probs(it["image"], it["target"]) -
                                       probs(items[(i + 1) % 3]["image"], it["target"])).sum(-1).mean()
                          for i, it in enumerate(items)])
    got = vision.image_use(model, proc, items, PROMPT, T)
    assert got["black"] == pytest.approx(want_black, abs=1e-6) and got["other"] == pytest.approx(want_other, abs=1e-6)
    assert 0.0 < got["black"] <= 1.0


def test_next_token_probs_keeps_every_input_on_the_model_device(monkeypatch):
    """R1 2026-11-05: on CUDA, mm_token_type_ids was extended with CPU zeros (the target ids stayed on the CPU) and
    G-vision crashed at the first rolling checkpoint. "meta" stands in for the GPU on the Mac."""
    pytest.importorskip("torch")
    from types import SimpleNamespace

    import torch

    from objdel.eval import vision

    dev = torch.device("meta")
    monkeypatch.setattr(vision, "_batch", lambda *a: {
        "input_ids": torch.zeros(1, 4, dtype=torch.long, device=dev),
        "attention_mask": torch.ones(1, 4, dtype=torch.long, device=dev),
        "mm_token_type_ids": torch.zeros(1, 4, dtype=torch.long, device=dev)})
    proc = SimpleNamespace(tokenizer=lambda *a, **k: {"input_ids": torch.tensor([[5, 6, 7]])})  # on the CPU

    def model(**b):
        assert {k: v.device for k, v in b.items()} == {k: dev for k in b} and b["mm_token_type_ids"].shape == (1, 7)
        return SimpleNamespace(logits=torch.zeros(1, 7, 11))

    p = vision.next_token_probs(model, proc, None, "x", T=3)
    assert p.shape == (3, 11) and np.allclose(p.sum(-1), 1.0)


@pytest.mark.parametrize("ratio,tv_black,tv_other,action", [
    (0.9, 0.01, 0.01, "ok"), (0.4, 0.01, 0.01, "warn"), (0.9, 0.001, 0.01, "warn"), (0.9, 0.01, 0.0019, "warn"),
    (0.49, 0.0, 0.0, "abort"), (0.5, 0.002, 0.002, "ok"),
])
def test_vision_decision_rule(ratio, tv_black, tv_other, action):
    """Plan §5.1 / M8: thresholds ratio >= 0.5 and both TVs >= 0.002; in training, both failing aborts, one warns."""
    from objdel.eval import vision
    d = vision.decide(ratio, tv_black, tv_other)
    assert d["action"] == action and d["pass"] is (action == "ok")


def test_vision_module_imports_without_torch():
    import subprocess
    import sys
    code = ("import sys\n"
            "class Block:\n"
            "    def find_spec(self, name, path=None, target=None):\n"
            "        if name.split('.')[0] in ('torch', 'transformers'): raise ImportError(name)\n"
            "sys.meta_path.insert(0, Block())\n"
            "from objdel.eval import vision; print(vision.decide(1, 1, 1)['action'])")
    out = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "ok"


def test_vision_fixed_sets_are_pinned_by_hash(tmp_path):
    """fixed_set: B1's selection (seeded permutation of the RL LIDC val records, item-id order), each image checked
    against its image_sha256; the set's SHA-256 changes with any id, image or target."""
    import hashlib

    from PIL import Image

    from objdel.eval import vision
    recs = []
    (tmp_path / "lidc").mkdir()
    for i, im in enumerate(_images(12)):
        item_id = f"lidc/x/n000/z{i:03d}"
        p = tmp_path / "images" / f"{item_id}.png"
        p.parent.mkdir(parents=True, exist_ok=True)
        im.save(p)
        a = np.asarray(im)
        sha = hashlib.sha256(repr(a.shape).encode() + a.tobytes()).hexdigest()
        recs.append({"item_id": item_id, "rl": i != 5, "image_sha256": sha, "targets": [None, f'{{"t":{i}}}']})
    (tmp_path / "lidc" / "val.jsonl").write_text("".join(json.dumps(r) + "\n" for r in recs))
    items, h = vision.fixed_set(tmp_path, 4)
    again, h2 = vision.fixed_set(tmp_path, 4)
    assert h == h2 and [it["item_id"] for it in items] == [it["item_id"] for it in again]
    rl = [r for r in recs if r["rl"]]
    want = sorted((rl[j] for j in np.random.default_rng(0).permutation(len(rl))[:vision.POOL]),
                  key=lambda r: r["item_id"])[:4]
    assert [it["item_id"] for it in items] == [r["item_id"] for r in want]
    assert all(it["target"] == r["targets"][1] for it, r in zip(items, want, strict=True))
    Image.fromarray(np.zeros((50, 50, 3), np.uint8)).save(tmp_path / "images" / f"{items[0]['item_id']}.png")
    with pytest.raises(ValueError, match="image_sha256"):
        vision.fixed_set(tmp_path, 4)


PREPARED = __import__("pathlib").Path(__file__).resolve().parents[1] / "data" / "prepared"
B1_SET_SHA256 = "b776dd12a085b69affc2cd1a07e4b723a9333f350c7409a76c0c562a5ea4f30e"  # 8 items, data/prepared (M4.md)


@pytest.mark.data
def test_vision_fixed_set_on_the_prepared_val_is_pinned():
    """The G-vision sets of B1 §2.5 (8 LIDC val items, the first 4 also the image set), fixed by hash."""
    from objdel.eval import vision
    if not (PREPARED / "lidc" / "val.jsonl").is_file():
        pytest.skip("data/prepared not found")
    items, h = vision.fixed_set(PREPARED, vision.N_ITEMS)
    assert h == B1_SET_SHA256 and len(items) == 8 and all(it["image"].size == items[0]["image"].size for it in items)
