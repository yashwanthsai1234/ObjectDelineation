"""Fair ceiling (Implementation.md A2; tests U19-U24, P10-P12, I16, I19): each doctor and the model are scored
against the same k - 1 doctors. CPU build of the checker, synthetic masks."""

import json
from types import SimpleNamespace

import numpy as np
import pytest

from objdel.verify import ceiling as C
from objdel.verify.reward import score
from objdel.verify.reward.cache import ExpertCache

H = W = 64


def disc(cy, cx, r, hole=0):
    y, x = np.ogrid[:H, :W]
    d2 = (y - cy) ** 2 + (x - cx) ** 2
    return (d2 <= r * r) & (d2 > hole * hole)


def rect(y0, y1, x0, x1):
    m = np.zeros((H, W), bool)
    m[y0:y1, x0:x1] = True
    return m


def plain_iou(a, b):
    u = np.count_nonzero(a | b)
    return np.count_nonzero(a & b) / u if u else 1.0


FOUR = np.stack([disc(30, 30, 10), disc(31, 32, 11), disc(29, 30, 9), disc(32, 29, 12)])


# U19 -----------------------------------------------------------------------------------------------------------------

def test_k4_leave_one_out_drops_one_rater_and_keeps_order():
    loos = C.loo_caches(ExpertCache.build(FOUR))
    assert [lo.index for lo in loos] == [0, 1, 2, 3]
    for lo in loos:
        assert isinstance(lo, tuple) and lo._fields == ("index", "cache", "single")
        assert not lo.single and lo.cache.k == 3
        assert lo.cache.digest == ExpertCache.build(np.delete(FOUR, lo.index, 0)).digest


def test_k3_leave_one_out_keeps_two_raters():
    loos = C.loo_caches(ExpertCache.build(FOUR[:3]))
    assert len(loos) == 3 and all(not lo.single and lo.cache.k == 2 for lo in loos)


def test_k2_leave_one_out_is_the_other_rater_twice():
    loos = C.loo_caches(ExpertCache.build(FOUR[:2]))
    assert [lo.index for lo in loos] == [0, 1]
    for lo, other in zip(loos, (FOUR[1], FOUR[0]), strict=True):
        assert lo.single and lo.cache.k == 2
        assert np.array_equal(lo.cache.raters[0], other) and np.array_equal(lo.cache.raters[1], other)


def test_a_non_voter_gets_no_entry_but_stays_in_every_stack():
    speck = np.zeros((H, W), bool)
    speck[5, 5] = True  # below the size filter of the union: no topology vote
    stack = np.concatenate([FOUR[:3], speck[None]])
    cache = ExpertCache.build(stack)
    assert cache.voters == (0, 1, 2)
    loos = C.loo_caches(cache)
    assert [lo.index for lo in loos] == [0, 1, 2]
    assert all(lo.cache.k == 3 and lo.cache.raters[2].sum() == 1 and lo.cache.raters[2][5, 5] for lo in loos)
    row = C.item(cache, None, None, W, H)
    assert row["k"] == 4 and row["n_voters"] == 3 and not row["single"]


def test_k1_raises_a_clear_value_error():
    one = SimpleNamespace(raters=FOUR[:1], k=1, voters=(0,))
    with pytest.raises(ValueError, match="k >= 2"):
        C.loo_caches(one)


# P10 -----------------------------------------------------------------------------------------------------------------

def test_k2_leave_one_out_scores_equal_plain_iou():
    A, B = disc(30, 30, 10), rect(22, 40, 25, 44)
    loos = C.loo_caches(ExpertCache.build(np.stack([A, B])))
    rng = np.random.default_rng(0)
    masks = [disc(*rng.integers(15, 50, 2), rng.integers(1, 15)) for _ in range(200)]
    masks += [np.zeros((H, W), bool), A, B, rect(0, 2, 0, 64), rect(50, 60, 50, 60), rect(10, 50, 31, 32)]
    for m in masks:
        assert score.score_mask(m, loos[0].cache).iou_tc == pytest.approx(plain_iou(m, B), abs=1e-12)
        assert score.score_mask(m, loos[1].cache).iou_tc == pytest.approx(plain_iou(m, A), abs=1e-12)
    assert score.score_mask(B, loos[0].cache).iou_tc == 1.0
    assert score.score_mask(rect(50, 60, 50, 60), loos[0].cache).iou_tc == 0.0


def test_k2_both_empty_scores_one_not_nan():
    empty = np.zeros((H, W), bool)
    (lo, _) = C.loo_caches(ExpertCache.build(np.stack([empty, empty])))
    assert score.score_mask(empty, lo.cache).iou_tc == 1.0


# item(): doctor side and legacy values (U20) -------------------------------------------------------------------------

def text_of(mask):
    """The repo's target text of a mask (exact for LIDC nodules)."""
    from objdel.data.targets.build import build_target, min_areas

    return build_target({"nodule": mask}, ("nodule",), None, min_areas({"nodule": mask[None]})).text


def test_identical_doctors_score_one():
    row = C.item(ExpertCache.build(np.stack([FOUR[0]] * 3)), None, None, W, H)
    assert row["doctor"] == 1.0 and row["legacy_dmax"] == [1.0, 1.0, 1.0] and row["legacy_loo_dmax_item"] == 1.0


def test_k2_doctor_is_the_plain_iou_of_the_pair():
    A, B = FOUR[0], rect(22, 40, 25, 44)
    row = C.item(ExpertCache.build(np.stack([A, B])), None, None, W, H)
    x = plain_iou(A, B)
    assert row["single"] and row["doctor"] == pytest.approx(x, abs=1e-12) and row["legacy_dmax"] == [x, x]


def test_one_signature_doctor_equals_the_legacy_max_iou():
    row = C.item(ExpertCache.build(FOUR), None, None, W, H)
    assert row["doctor"] == pytest.approx(row["legacy_loo_dmax_item"], abs=1e-12)


def test_legacy_dmax_is_the_max_plain_iou_with_any_other_doctor():
    # intervals [0,20), [0,40), [0,50) on rows 10-20: IoU A-B .5, A-C .4, B-C .8 (the plan's .5/.2/.8 breaks the
    # Jaccard triangle inequality, so no masks have it)
    stack = np.stack([rect(10, 20, 0, 20), rect(10, 20, 0, 40), rect(10, 20, 0, 50)])
    row = C.item(ExpertCache.build(stack), None, None, W, H)
    assert row["legacy_dmax"] == [0.5, 0.8, 0.8] and row["legacy_loo_dmax_item"] == pytest.approx(0.7, abs=1e-12)


# item(): answers (U21-U23, P12) --------------------------------------------------------------------------------------

BOWTIE = '{"nodule":[{"outer":[[0,0],[1000,1000],[1000,0],[0,1000]],"holes":[]}]}'  # self-intersecting: illegal


def test_without_answers_the_model_keys_are_null():
    row = C.item(ExpertCache.build(FOUR), None, None, W, H)
    assert {k: row[k] for k in ("model", "model_all_k", "diff")} == dict.fromkeys(("model", "model_all_k", "diff"))
    assert {"doctor", "legacy_dmax", "legacy_loo_dmax_item", "k", "single", "n_voters"} <= set(row)


def test_illegal_failed_and_truncated_answers_score_zero():
    cache = ExpertCache.build(FOUR)
    good = text_of(FOUR[0])
    row = C.item(cache, [good, BOWTIE, "not json"], [False] * 3, W, H)
    ref = [score.score_mask(FOUR[0], lo.cache).iou_tc for lo in C.loo_caches(cache)]
    assert row["model"] == pytest.approx(np.mean(ref) / 3, abs=1e-12)
    assert row["model_all_k"] == pytest.approx(score.score_mask(FOUR[0], cache).iou_tc / 3, abs=1e-12)
    assert C.item(cache, [good], [True], W, H)["model"] == 0.0  # parse refuses a truncated answer
    with pytest.raises(ValueError, match="truncated"):
        C.item(cache, [good], [False, False], W, H)


def test_an_answer_equal_to_the_doctors_scores_like_them():
    row = C.item(ExpertCache.build(np.stack([FOUR[1]] * 3)), [text_of(FOUR[1])], [False], W, H)
    assert row["model"] == row["doctor"] == 1.0 and row["diff"] == 0.0


def test_the_entry_r_term_of_an_answer_equal_to_doctor_r_is_the_doctor_term():
    cache = ExpertCache.build(FOUR)
    for lo in C.loo_caches(cache):
        (t, _), = score.score_many([(text_of(FOUR[lo.index]), False, ("nodule",), {"nodule": lo.cache}, W, H)])
        assert t["legal"] == 1.0 and t["iou_tc"] == score.score_mask(FOUR[lo.index], lo.cache).iou_tc


def test_more_references_can_raise_model_all_k_above_model():
    row = C.item(ExpertCache.build(FOUR[:3]), [text_of(FOUR[0])], [False], W, H)
    assert row["model_all_k"] == 1.0 > row["model"]


def test_doctor_and_model_are_scored_on_the_same_k_minus_1_doctors(monkeypatch):
    cache = ExpertCache.build(FOUR)
    seen = {"doctor": [], "model": []}
    real_mask, real_many = score.score_mask, score.score_many

    def spy_mask(mask, c):
        seen["doctor"].append(c.digest)
        return real_mask(mask, c)

    def spy_many(jobs):
        seen["model"] += [job[3]["nodule"].digest for job in jobs]
        out = real_many(jobs)
        assert all(t["iou_tc"] <= 1.0 for t, _ in out)
        return out

    monkeypatch.setattr(score, "score_mask", spy_mask)
    monkeypatch.setattr(score, "score_many", spy_many)
    texts = [text_of(m) for m in FOUR] * 2
    C.item(cache, texts, [False] * 8, W, H)
    loo = [lo.cache.digest for lo in C.loo_caches(cache)]
    assert seen["doctor"] == loo  # one doctor call per left-out r, on its cache
    assert seen["model"] == [d for d in loo for _ in texts] + [cache.digest] * 8  # every LOO cache, then model_all_k


# P11: the doctor side is backed by a check ---------------------------------------------------------------------------

def _named_shapes(S=125):
    y, x = np.ogrid[:S, :S]
    d2 = (y - 60) ** 2 + (x - 60) ** 2
    ring = (d2 <= 25 ** 2) & (d2 > 10 ** 2)
    return {"disc": d2 <= 20 ** 2, "ellipse": ((y - 60) / 30) ** 2 + ((x - 60) / 12) ** 2 <= 1, "ring (hole)": ring,
            "C with a hole": ring & ~((y > 55) & (y < 65) & (x > 60)), "thin diagonal": (y == x) & (x < 90),
            "touches all four borders": ~(d2 <= 10 ** 2), "1-px sliver": (y == 40) & (x > 10) & (x < 90),
            "one pixel": (y == 3) & (x == 4)}


def test_target_text_rasterizes_to_its_rater_mask():
    """A doctor's target text (the repo's writer), parsed and rasterized, is the doctor's mask pixel for pixel: the
    doctor gains nothing from raw masks over the model's grid text. Holes round-trip (no known loss)."""
    from objdel.geom.shape import to_mask, to_masks
    from objdel.geom.text import parse

    rng = np.random.default_rng(0)
    masks = list(_named_shapes().values())
    while len(masks) < 208:  # 200 random simple (star-shaped) polygons on the LIDC frame sides
        S, k = int(rng.choice([50, 100, 125])), int(rng.integers(3, 9))
        ang, rad = np.sort(rng.uniform(0, 2 * np.pi, k)), rng.uniform(50, 300, k)
        c = rng.uniform(200, 800, 2)
        pts = np.clip(np.round(c + np.c_[np.cos(ang) * rad, np.sin(ang) * rad]), 0, 1000).astype(int)
        g = to_mask([(pts, [])], S, S)
        if g.legal and g.mask.any():
            masks.append(g.mask)
    for m in masks:
        p = parse(text_of(m), ("nodule",))
        g, = to_masks([(p.objs["nodule"], m.shape[1], m.shape[0])])
        assert np.array_equal(g.mask, m)


# summarize (U24, I24) ------------------------------------------------------------------------------------------------

def _row(i, k, doctor, dmax, model=None, n_voters=None):
    m = None if model is None else {"model": model, "model_all_k": model + .05, "diff": model - doctor}
    return {"item_id": f"lidc/s{i}/n000/z{i:03d}", "k": k, "n_voters": n_voters or k, "single": k == 2,
            "doctor": doctor, "legacy_dmax": dmax, "legacy_loo_dmax_item": float(np.mean(dmax)),
            **(m or dict.fromkeys(("model", "model_all_k", "diff")))}


KS = (3, 3, 3, 4, 4, 4, 2, 2, 2, 2)


def _rows(model=True):
    rng = np.random.default_rng(1)
    return [_row(i, k, float(rng.uniform(.5, .9)), list(rng.uniform(.5, 1, k)),
                 float(rng.uniform(0, .8)) if model else None, n_voters=k - (i == 4)) for i, k in enumerate(KS)]


def _groups(rows):
    return {r["item_id"]: f"LIDC-IDRI-{i % 4:04d}" for i, r in enumerate(rows)}


def _walk(x):
    if isinstance(x, dict):
        return all(isinstance(k, str) and _walk(v) for k, v in x.items())
    if isinstance(x, list):
        return all(map(_walk, x))
    return type(x) in (float, int, bool, str, type(None))


def test_summarize_strata_add_up_and_use_patient_then_nodule_draws():
    from objdel.verify import stats as signal

    rows, groups = _rows(), _groups(_rows())
    s = C.summarize(rows, groups, boot=300)
    assert s["all"]["n"] == s["k>=3"]["n"] + s["k=2"]["n"] == 10 and s["k=2"]["n"] == 4
    assert s["found"]["n"] + s["miss"]["n"] == 10
    assert s["found"]["n"] == sum(r["model_all_k"] >= signal.MISS for r in rows) > 0 and s["miss"]["n"] > 0
    assert s["all"]["n_patients"] == 4 and s["all"]["n_voters_lt_raters"] == 1 and s["k=2"]["n_voters_lt_raters"] == 0
    assert s["k=2"]["single_reference"] is True and "single_reference" not in s["k>=3"]
    ids = [r["item_id"] for r in rows]
    for m in C.METRICS:
        x = np.array([r[m] for r in rows])
        mean, lo, hi = s["all"][m]
        assert mean == pytest.approx(x.mean(), abs=1e-12) and lo <= mean <= hi
        assert s["all"][m] == signal.ci(x, signal.draws(ids, 300, 0, key=groups.__getitem__))
        assert s["all"]["nodule"][m] == signal.ci(x, signal.draws(ids, 300, 0))
    d = np.array([r["model"] - r["doctor"] for r in rows])  # paired: the CI of the per-item differences
    assert s["all"]["diff"] == signal.ci(d, signal.draws(ids, 300, 0, key=groups.__getitem__))
    assert _walk(s) and json.dumps(s, allow_nan=False)


def test_summarize_pools_the_legacy_values_over_doctors():
    rows = [_row(0, 2, .5, [.5, .5]), _row(1, 4, 1., [1., 1., 1., 1.])]
    s = C.summarize(rows, boot=50)
    assert s["all"]["legacy_loo_dmax_pooled"] == pytest.approx(5 / 6, abs=1e-12)
    assert s["all"]["legacy_loo_dmax_item"][0] == pytest.approx(.75, abs=1e-12)
    same_k = [r for r in _rows(model=False) if r["k"] == 3]
    t = C.summarize(same_k, boot=50)["all"]
    assert t["legacy_loo_dmax_pooled"] == pytest.approx(t["legacy_loo_dmax_item"][0], abs=1e-12)


def test_summarize_without_answers_or_groups():
    s = C.summarize([r for r in _rows(model=False) if r["k"] == 2], boot=50)
    assert set(s) == {"all", "k>=3", "k=2"}  # no found / miss without answers
    assert s["k>=3"]["n"] == 0 and s["k>=3"]["doctor"] is None and s["k>=3"]["legacy_loo_dmax_pooled"] is None
    assert s["all"]["model"] is None and s["all"]["nodule"]["diff"] is None
    assert s["all"]["n_patients"] is None and s["all"]["doctor"][1:] == [None, None]  # no patient line
    assert json.dumps(s, allow_nan=False)


def test_summarize_refuses_answers_without_groups_and_mixed_answers():
    with pytest.raises(ValueError, match="groups"):
        C.summarize(_rows(), boot=50)
    rows = _rows()
    rows[3] = _row(3, 4, .7, [.8] * 4)
    with pytest.raises(ValueError, match=rows[3]["item_id"]):
        C.summarize(rows, _groups(rows), boot=50)


# CLI (I16, I19) ------------------------------------------------------------------------------------------------------

def _view(tmp_path, synthetic_items):
    """A prepared val root: 6 LIDC items of 3 patients, k = 2, 3, 4, one with a non-voting doctor (and conftest's
    RIGA item: prepare builds every dataset)."""
    from objdel.data import prepare as data_prepare
    from objdel.data import store
    from objdel.data.store import Item
    from objdel.train import dataset as D

    S, rng = 100, np.random.default_rng(0)
    y, x = np.ogrid[:S, :S]

    def d(cy, cx, r):
        return (y - cy) ** 2 + (x - cx) ** 2 <= r * r

    speck = (y == 5) & (x == 5)
    stacks = [[d(50, 50, 12), d(51, 50, 11)], [d(40, 45, 9), d(41, 44, 10), d(40, 46, 8)],
              [d(60, 50, 14), d(59, 51, 13), d(60, 49, 15), d(61, 50, 12)], [d(50, 50, 12), d(52, 49, 12), speck],
              [d(45, 55, 7), d(46, 55, 8)], [d(50, 40, 10), d(50, 41, 11), d(49, 40, 10)]]
    items = [Item(f"lidc/1.2.{i // 2}/n000/z00{i}", "lidc", f"LIDC-IDRI-000{i // 2}", "val",
                  rng.integers(-1000, 400, (S, S), dtype=np.int16), {"nodule": np.stack(st)},
                  [str(j + 1) for j in range(len(st))], {}) for i, st in enumerate(stacks)]
    riga = [it for it in synthetic_items(("val",)) if it.dataset == "riga"]
    store.write(riga + items, tmp_path / "store")
    data_prepare.prepare(tmp_path / "store", tmp_path / "view", splits=("val",), procs=1, prompts=D.PROMPTS)
    return tmp_path / "view"


@pytest.fixture(scope="module")
def view_root(tmp_path_factory, synthetic_items):
    return _view(tmp_path_factory.mktemp("ceiling"), synthetic_items)


@pytest.fixture
def view(view_root, clean_store):
    return view_root


def _answers(view, path, drop=None):
    from objdel.train import dataset as D

    recs = [r for r in D.records(view, "lidc", "val") if r["rl"] and r["item_id"] != drop]
    lines = [{"item_id": r["item_id"], "group_id": r["group_id"], "W": r["W"], "H": r["H"],
              "samples": ([t for t in r["targets"] if t] * 8)[:7] + ["garbage"], "truncated": [False] * 8,
              "tokens": [10] * 8} for r in recs]
    path.write_text("".join(json.dumps(x) + "\n" for x in lines))
    return path


def _run(view, out, *extra):
    from experiments.phase2 import ceiling as CLI

    CLI.main(["--view", str(view), "--split", "val", "--out", str(out), *extra])
    return json.loads((out / "ceiling.json").read_text()), (out / "ceiling_items.jsonl").read_text()


def test_cli_without_answers_writes_the_doctor_side(view, tmp_path):
    doc, items = _run(view, tmp_path / "out")
    rows = [json.loads(x) for x in items.splitlines()]
    assert len(rows) == 6 and all(r["model"] is r["model_all_k"] is r["diff"] is None for r in rows)
    assert [r["k"] for r in rows] == [2, 3, 4, 3, 2, 3] and rows[3]["n_voters"] == 2
    assert doc["all"]["n"] == 6 and doc["k=2"]["n"] == 2 and doc["all"]["n_patients"] == 3
    assert doc["all"]["n_voters_lt_raters"] == 1 and doc["seconds"] > 0 and doc["answers"] is None
    assert "found" not in doc and json.dumps(doc, allow_nan=False)


def test_cli_joins_the_answers_by_item_id(view, tmp_path):
    doc, items = _run(view, tmp_path / "out", "--answers", str(_answers(view, tmp_path / "answers.jsonl")))
    rows = [json.loads(x) for x in items.splitlines()]
    assert all(0 < r["model"] <= 1 and r["diff"] == r["model"] - r["doctor"] for r in rows)
    assert doc["found"]["n"] + doc["miss"]["n"] == 6


def test_cli_refuses_answers_that_miss_an_item(view, tmp_path):
    lost = "lidc/1.2.1/n000/z003"
    with pytest.raises(ValueError, match=lost):
        _run(view, tmp_path / "out", "--answers", str(_answers(view, tmp_path / "a.jsonl", drop=lost)))
    assert not (tmp_path / "out").exists()


def test_cli_scores_in_one_process(view, tmp_path):
    """I19 after D7: two runs of the CLI, which scores in its own process, give the same rows."""
    answers = str(_answers(view, tmp_path / "answers.jsonl"))
    got = [_run(view, tmp_path / f"out{n}", "--answers", answers)[1] for n in (1, 2)]
    assert got[0] == got[1]
