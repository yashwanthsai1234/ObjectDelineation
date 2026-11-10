"""Break-it B2 (reward, TRL wiring, eval): reproducers of the confirmed breakdowns, now fixed (ids B2-<n>, see
experiments/phase1/impl/breakit/B2_reward.md) and the robustness checks that held (plain tests). No data needed."""
import hashlib
import json

import numpy as np
import pytest

from objdel.data.audit import hack_audit
from objdel.verify import gates
from objdel.verify.metrics import metrics, summarize
from objdel.verify.reward import cache as cache_mod
from objdel.verify.reward.advantages import advantages_v4, topology_reward
from objdel.verify.reward.cache import ExpertCache
from objdel.verify.reward.score import score_batch, score_completion, score_mask, sig_id

W = 50
NOD = ("nodule",)


def square(x0, y0, x1, y1):
    m = np.zeros((W, W), bool)
    m[y0:y1, x0:x1] = True
    return m


def ring_pts(x0, y0, x1, y1, ccw=False):
    r = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
    return r[::-1] if ccw else r


def doc(**objs):
    return json.dumps({n: [{"outer": o, "holes": h} for o, h in parts] for n, parts in objs.items()},
                      separators=(",", ":"))


@pytest.fixture(scope="module")
def nod():
    """Three raters, squares of ~20 px on a 50 x 50 crop (20 grid units per pixel)."""
    st = np.stack([square(10, 10, 30, 30), square(11, 11, 31, 31), square(10, 10, 32, 32)])
    return st, {"nodule": ExpertCache.build(st)}


def terms_of(masks, cache):
    """Single-object completion terms of legal masks (what score_completion returns for them)."""
    rows = []
    for m in masks:
        t = score_mask(m, cache)
        rows.append(dict(format_ok=1.0, legal=1.0, iou_tc=t.iou_tc, f5=t.f5, g_area=t.g_area,
                         topo_any=float(t.topo_freq > 0), topo_freq=t.topo_freq, sig_id=sig_id((t.sig,))))
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


# ------------------------------------------------------------------ reward semantics

def test_b2_4_sub_filter_rater_does_not_make_a_speck_a_valid_answer():
    st = np.stack([square(20, 20, 22, 22), square(14, 14, 21, 20)])  # 4 px and 42 px (LIDC: 18 such items)
    c = ExpertCache.build(st)
    assert c.rater_sigs == ((), (0,))
    honest = [st[1]] * 7
    T = terms_of(honest + [st[0]], c)  # the eighth rollout: a 4 px speck on the tiny rater
    a = advantages_v4(T)
    assert topology_reward(T)[7] == 0.0  # not a valid topology answer
    assert a[7] <= np.median(a[:7])  # and it does not beat the honest median (today: +2.77 vs +0.29 on real data)


def test_b2_7_missed_cup_is_not_a_valid_topology_answer():
    disc = np.stack([square(10, 10, 40, 40)] * 3)
    cup = np.stack([square(20, 20, 30, 30), square(20, 20, 31, 31), square(19, 19, 30, 30)])
    caches = {"disc": ExpertCache.build(disc), "cup": ExpertCache.build(cup)}
    good = doc(disc=[(ring_pts(200, 200, 800, 800), [])], cup=[(ring_pts(400, 400, 600, 600), [])])
    missed = doc(disc=[(ring_pts(200, 200, 800, 800), [])], cup=[(ring_pts(0, 0, 200, 200), [])])
    T = score_batch([(good if i < 7 else missed, False, ("disc", "cup"), caches, W, W) for i in range(8)])
    m = metrics([good, missed], {"disc": disc, "cup": cup}, W, W, caches=caches)
    assert m["samples"]["passed"][1] == 0.0  # the eval harness fails it (per object) ...
    assert topology_reward(T)[7] == 0.0  # ... the reward should not grant it the topology bonus


def test_b2_9_caches_with_different_rater_counts_are_refused():
    disc3 = ExpertCache.build(np.stack([square(10, 10, 40, 40)] * 2 + [square(11, 10, 40, 40)]))
    cup2 = ExpertCache.build(np.stack([square(18, 18, 30, 30), square(19, 19, 30, 30)]))
    t = doc(disc=[(ring_pts(200, 200, 800, 800), [])], cup=[(ring_pts(360, 360, 600, 600), [])])
    for caches in ({"disc": cup2, "cup": disc3}, {"disc": disc3, "cup": cup2}):
        with pytest.raises(ValueError):
            score_completion(t, False, ("disc", "cup"), caches, W, W)


def test_b2_10_cache_normalises_bool_bytes(tmp_path):
    y, x = np.mgrid[:40, :40]
    r = np.hypot(y - 19.5, x - 19.5)
    ring = (r <= 12) & (r > 4)
    st = np.stack([ring, ring, r <= 12])
    bad = (st.view(np.uint8) * 255).view(bool)  # e.g. np.asarray(PIL mode "1" image, bool), issue I1
    good, wrong = ExpertCache.build(st), ExpertCache.build(bad)
    assert good.rater_sigs == ((1,), (1,), (0,))
    assert wrong.digest == good.digest  # the integrity hash cannot tell them apart ...
    wrong.save(tmp_path / "c.npz")
    assert ExpertCache.load(tmp_path / "c.npz").rater_sigs == good.rater_sigs  # ... so this must hold


# ------------------------------------------------------------------ TRL wiring and cache loading

def test_b2_5_group_size_other_than_eight_is_refused(tmp_path, nod):
    pytest.importorskip("trl")
    from objdel.rl import config as configs
    from objdel.rl import reward as R
    _, caches = nod
    texts = [doc(nodule=[(ring_pts(200 + 5 * i, 200 + 5 * i, 600 + 7 * i, 600 + 7 * i), [])]) for i in range(16)]
    funcs, _ = R.TRLReward({"a": caches}.__getitem__, {7}).funcs()
    adv = np.array(funcs[-1](completions=texts, completion_ids=[[1, 7]] * 16, key=["a"] * 16, W=[W] * 16,
                             H=[W] * 16, names=[list(NOD)] * 16))
    T = score_batch([(t, False, NOD, caches, W, W) for t in texts])
    assert not np.allclose(adv, advantages_v4(T, 16))  # what TRL's group of 16 would need: not what it gets
    with pytest.raises(ValueError):
        configs.grpo_config("qwen3-vl-2b", str(tmp_path), num_generations=16, bf16=False, use_cpu=True)


def prepared(tmp_path, n=3):
    recs = []
    for i in range(n):
        m = np.zeros((2, 30, 30), bool)
        m[0, 5:20 + i, 5:20] = True
        m[1, 6:21, 5:20 + i] = True
        d = tmp_path / "caches" / f"it{i}"
        d.mkdir(parents=True)
        digest = ExpertCache.build(m).save(d / "nodule.npz")
        recs.append({"item_id": f"it{i}", "rl": True, "caches": {"nodule": digest}})
    (tmp_path / "lidc").mkdir()
    text = "".join(json.dumps(r) + "\n" for r in recs)
    (tmp_path / "lidc" / "train.jsonl").write_text(text)
    files = {"lidc/train.jsonl": {"sha256": hashlib.sha256(text.encode()).hexdigest()}}  # S1-14: cache_loader (I45)
    (tmp_path / "manifest.json").write_text(json.dumps({"files": files}))
    return tmp_path


def test_b2_8_cache_loader_after_a_fingerprint_change(tmp_path, monkeypatch):
    pytest.importorskip("datasets")
    from objdel.verify import store as dataset
    root = prepared(tmp_path)
    monkeypatch.setattr(cache_mod, "fingerprint", lambda: "numpy or opencv pin changed")
    try:
        caches = dataset.cache_loader(root)("it0")
    except ValueError as e:
        assert "fingerprint" in str(e)
    else:
        assert caches["nodule"].k == 2


def test_cache_loader_lru_eviction_reloads_identically(tmp_path):
    pytest.importorskip("datasets")
    from objdel.verify import store as dataset
    load = dataset.cache_loader(prepared(tmp_path), maxsize=2)
    first = load("it0")["nodule"].digest
    for k in ("it1", "it2", "it1"):
        load(k)
    assert load.cache_info().currsize == 2 and load("it0")["nodule"].digest == first


# ------------------------------------------------------------------ eval / hack audit

def test_b2_6_duplicate_alarm_sees_identical_masks(nod):
    st, caches = nod
    nudged = [doc(nodule=[([[200, 200], [600 + d, 200], [600, 600], [200, 600]], [])]) for d in range(8)]
    spaced = [nudged[0] + " " * i for i in range(8)]
    for texts in (nudged, spaced):
        m = metrics(texts, {"nodule": st}, W, W, caches=caches)
        assert set(m["samples"]["advantage"]) == {0.0}  # the reward sees one answer: the group is dropped
        assert m["groups"]["dup"][0] >= 0.5
        assert hack_audit([{"item_id": "x", "dataset": "lidc", **m}])["alarms"]["duplicates"]["fired"]


# ------------------------------------------------------------------ attacks that held (kept as regression checks)

@pytest.mark.parametrize("parts, legal", [
    ([(ring_pts(200, 200, 600, 600), [])] * 2, False),  # duplicated part
    ([(ring_pts(200, 200, 600, 600), []), (ring_pts(300, 300, 400, 400), [])], False),  # nested outers
    ([(ring_pts(200, 200, 400, 600), []), (ring_pts(400, 200, 600, 600), [])], False),  # parts sharing an edge
    ([(ring_pts(200, 200, 600, 600), [ring_pts(200, 200, 600, 600, ccw=True)])], False),  # hole = outer
    ([(ring_pts(200, 200, 600, 600), [ring_pts(300, 300, 400, 400)])], False),  # hole wound like an outer
    ([(ring_pts(200, 200, 600, 600, ccw=True), [])], False),  # counter-clockwise outer
    ([([[200, 200], [600, 200], [400, 200]], [])], False),  # zero-area spike
    ([(ring_pts(0, 0, 1000, 1000), [])], True),  # the whole crop
])
def test_text_hacks_are_illegal_or_score_low(nod, parts, legal):
    _, caches = nod
    honest = score_completion(doc(nodule=[(ring_pts(200, 200, 600, 600), [])]), False, NOD, caches, W, W)
    r = score_completion(doc(nodule=parts), False, NOD, caches, W, W)
    assert r["format_ok"] == 1.0 and r["legal"] == float(legal)
    T = {k: np.array([honest[k]] * 7 + [r[k]]) for k in honest if k not in ("error", "defects")}
    a = advantages_v4(T)
    assert a[7] < a[:7].min()


def test_ring_cap_and_truncation_are_format_failures(nod):
    _, caches = nod
    eleven = doc(nodule=[(ring_pts(200 + 30 * i, 200, 210 + 30 * i, 210), []) for i in range(11)])
    assert score_completion(eleven, False, NOD, caches, W, W)["error"] == "cap_rings"
    ok = doc(nodule=[(ring_pts(200, 200, 600, 600), [])])
    assert score_completion(ok, True, NOD, caches, W, W)["format_ok"] == 0.0


def _g(n=8, **kw):
    base = dict(format_ok=np.ones(n), legal=np.ones(n), iou_tc=np.full(n, .9), f5=np.full(n, .9), g_area=np.full(n, .9),
                topo_any=np.ones(n), topo_freq=np.ones(n), sig_id=np.full(n, sig_id(((0,),))))
    return base | {k: np.asarray(v, float) for k, v in kw.items()}


@pytest.mark.parametrize("terms", [
    _g(legal=np.zeros(8), iou_tc=np.linspace(0, 1, 8)),  # every rollout illegal
    _g(),  # identical
    _g(format_ok=np.zeros(8), legal=np.zeros(8), iou_tc=np.zeros(8), f5=np.zeros(8), g_area=np.zeros(8),
       topo_any=np.zeros(8), topo_freq=np.zeros(8), sig_id=-np.ones(8)),  # every rollout a format failure
    _g(legal=np.r_[1, np.zeros(7)]),  # one legal among seven illegal
    _g(iou_tc=np.r_[np.full(4, .85), np.full(4, .95)]),  # sigma exactly at the 0.05 floor
    _g(sig_id=np.r_[sig_id(((1,),)), np.full(7, sig_id(((0,),)))], topo_freq=np.r_[1 / 6, np.full(7, 5 / 6)]),
])
def test_edge_groups_give_finite_centred_advantages(terms):
    a = advantages_v4(terms)
    assert np.isfinite(a).all() and abs(a.sum()) < 1e-9 and np.isfinite(topology_reward(terms)).all()


@pytest.mark.parametrize("group", [2, 3, 4, 16])
def test_other_group_sizes_are_finite(group):
    t = _g(group, iou_tc=np.linspace(.5, .9, group), topo_freq=np.full(group, .5),
           sig_id=np.r_[np.full(group - 1, sig_id(((0,),))), sig_id(((1,),))])
    assert np.isfinite(advantages_v4(t, group)).all()


def test_eval_on_all_truncated_and_mixed_k_items(nod):
    st, caches = nod
    texts = [doc(nodule=[(ring_pts(200, 200, 600, 600), [])])] * 8
    allt = metrics(texts, {"nodule": st}, W, W, truncated=[True] * 8, caches=caches)
    assert allt["tv"] == 1.0 and allt["min_mode"] == 0.0 and set(allt["samples"]["advantage"]) == {0.0}
    items = [{"item_id": "a", "dataset": "lidc", **allt},
             {"item_id": "b", "dataset": "lidc", **metrics(texts[:2], {"nodule": st}, W, W, caches=caches)}]
    s = summarize(items)["datasets"]["lidc"]
    assert s["group"] is None and set(s["pass_at_k"]) == {"1", "2"}
    audit = hack_audit(items)
    assert audit["alarms"]["components"]["fired"] is None  # no baseline: not measured
    assert gates.phase2_entry({"summary": summarize(items)})["pass"] is False  # 8 of 10 samples truncated: format fails
