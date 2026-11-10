"""objdel.train.rewards.TRLReward without a model: terms, advantages, truncation, group checks, typed logs."""
from types import SimpleNamespace

import numpy as np
import pytest

from objdel.geom.engine import DEFECTS
from objdel.train import rewards as R
from objdel.verify.reward import advantages as reward_advantages
from objdel.verify.reward.advantages import advantages_v4, legal_gate, topology_reward
from objdel.verify.reward.cache import ExpertCache
from objdel.verify.reward.score import score_batch

G, W = 8, 50
NAMES = ("nodule",)
EOS, PAD = 7, 0


def square(x0, y0, x1, y1):
    return f'{{"nodule":[{{"outer":[[{x0},{y0}],[{x1},{y0}],[{x1},{y1}],[{x0},{y1}]],"holes":[]}}]}}'


def rater_stack():
    m = np.zeros((3, W, W), bool)
    for i, (a, b) in enumerate([(10, 30), (11, 31), (10, 32)]):
        m[i, a:b, a:b] = True
    return m


CACHES = {"a": {"nodule": ExpertCache.build(rater_stack())}, "b": {"nodule": ExpertCache.build(rater_stack()[:, ::-1])}}
BOWTIE = '{"nodule":[{"outer":[[200,200],[600,600],[600,200],[200,600]],"holes":[]}]}'


def group_texts(rng):
    out = [square(200, 200, 600, 600), square(220, 220, 620, 620), square(200, 200, 640, 640), BOWTIE, "nonsense",
           square(100, 100, 900, 900), square(200, 200, 610, 600), square(200, 200, 600, 600)]
    return [out[i] for i in rng.permutation(G)]


def batch(seed=0, truncate=()):
    rng = np.random.default_rng(seed)
    texts = group_texts(rng) + group_texts(rng)
    keys = ["a"] * G + ["b"] * G
    ids = [[5, 6, EOS] for _ in texts]
    for i in truncate:
        ids[i] = [5, 6, 6]
    return dict(prompts=[[{"role": "user", "content": "p"}]] * 2 * G,
                completions=[[{"role": "assistant", "content": t}] for t in texts], completion_ids=ids,
                key=keys, W=[W] * 2 * G, H=[W] * 2 * G, names=[list(NAMES)] * 2 * G), texts, keys


def reward():
    return R.TRLReward(CACHES.__getitem__, stop_ids={EOS, PAD})


def call_all(rw, kw, **extra):
    funcs, _ = rw.funcs()
    return {f.__name__: f(**kw, **extra) for f in funcs}


def test_funcs_are_weight_zero_loggers_plus_one_aggregator():
    funcs, weights = reward().funcs()
    assert [f.__name__ for f in funcs] == [*reward_advantages.LOGGED, R.AGGREGATOR]
    assert weights == [0.0] * len(reward_advantages.LOGGED) + [1.0]


def test_values_equal_score_batch_and_advantages_v4():
    kw, texts, keys = batch(truncate=(1,))
    got = call_all(reward(), kw)
    jobs = [(t, i == 1, NAMES, CACHES[k], W, W) for i, (t, k) in enumerate(zip(texts, keys, strict=True))]
    ref = score_batch(jobs)
    legal = legal_gate(ref)
    for name in ("format_ok", "legal", "topo_any"):
        assert got[name] == ref[name].tolist()
    for name in reward_advantages.GATED:  # the geometric terms as advantages_v4 sees them: 0 for an illegal rollout
        assert got[name] == (ref[name] * legal).tolist()
    assert got["r_topo"] == topology_reward(ref).tolist()
    assert got[R.AGGREGATOR] == advantages_v4(ref).tolist()
    assert np.ptp(got[R.AGGREGATOR]) > 0  # a non-degenerate batch
    assert got["format_ok"][1] == 0.0  # truncated: no EOS at the end


def test_one_computation_per_batch(monkeypatch):
    calls = []
    monkeypatch.setattr(R, "score_batch", lambda jobs: calls.append(1) or score_batch(jobs))
    rw = reward()
    kw, _, _ = batch()
    call_all(rw, kw)
    assert len(calls) == 1
    kw2, _, _ = batch(seed=1)
    call_all(rw, kw2)
    assert len(calls) == 2


def test_truncation_comes_from_completion_ids():
    kw, texts, _ = batch()
    good = texts.index(square(200, 200, 600, 600))
    for ids, expect in (([5, EOS], 1.0), ([5, PAD], 1.0), ([5, 6], 0.0)):
        kw["completion_ids"] = [list(i) for i in kw["completion_ids"]]
        kw["completion_ids"][good] = ids
        assert call_all(reward(), kw)["format_ok"][good] == expect


@pytest.mark.parametrize("cut", [slice(0, 12), slice(4, 12)])
def test_partial_groups_are_refused(cut):
    kw, _, _ = batch()
    kw = {k: v[cut] for k, v in kw.items()}
    with pytest.raises(ValueError, match="whole groups"):
        call_all(reward(), kw)


def test_interleaved_groups_are_refused():
    kw, _, _ = batch()
    order = [i for pair in zip(range(G), range(G, 2 * G), strict=True) for i in pair]
    kw = {k: [v[i] for i in order] for k, v in kw.items()}
    with pytest.raises(ValueError, match="whole groups"):
        call_all(reward(), kw)


def test_typed_error_logs():
    metrics, extra = {}, {}
    kw, texts, _ = batch()
    call_all(reward(), kw, log_metric=metrics.__setitem__, log_extra=extra.__setitem__)
    assert extra["error"].count("") == sum(t != "nonsense" for t in texts)
    assert metrics["error/syntax"] == pytest.approx(2 / 16)
    assert {f"defect/{d}" for d in DEFECTS} <= set(metrics)
    assert sum(metrics[f"defect/{d}"] for d in DEFECTS) > 0  # the bow-tie crosses itself
    assert all(0 <= v <= 1 for v in metrics.values())


def test_stop_ids_follow_trl():
    model = SimpleNamespace(generation_config=SimpleNamespace(eos_token_id=[3, 4]))
    proc = SimpleNamespace(tokenizer=SimpleNamespace(eos_token_id=9, pad_token_id=0))
    assert R.stop_ids(model, proc) == {3, 4, 9, 0}
    model.generation_config.eos_token_id = 9
    assert R.stop_ids(model, proc.tokenizer) == {9, 0}


def test_one_truncation_rule():
    """TRL's rule (`ids[-1] not in eos_and_pad`; an empty completion has no end either), shared by the reward and
    objdel.train.sample (REVIEW #7)."""
    stop = {EOS, PAD}
    assert [R.truncated(ids, stop) for ids in ([], [5], [5, EOS], [5, PAD], [EOS, 5])] == [True, True, False, False,
                                                                                            True]


SPIKE = '{"nodule":[{"outer":[[200,200],[400,200],[400,100],[400,200],[600,200],[600,600],[200,600]],"holes":[]}]}'


def test_an_illegal_high_iou_answer_does_not_raise_the_logged_geometry():
    """R2 follow-up: the best-checkpoint metric (eval_rewards/iou_tc/mean) is the iou_tc logger's mean, so an illegal
    answer (a square with a collinear spike: IoU_tc 1, legal 0) must log 0, as advantages_v4 scores it; the
    aggregator's values are unchanged (reward v4 is frozen)."""
    kw, texts, keys = batch()
    honest = call_all(reward(), kw)
    swap = [i for i, (t, k) in enumerate(zip(texts, keys, strict=True))
            if k == "a" and t == square(220, 220, 620, 620)]  # legal, IoU_tc < 1
    assert swap
    kw["completions"] = [[{"role": "assistant", "content": SPIKE if i in swap else t}] for i, t in enumerate(texts)]
    hacked = call_all(reward(), kw)
    jobs = [(SPIKE if i in swap else t, False, NAMES, CACHES[k], W, W)
            for i, (t, k) in enumerate(zip(texts, keys, strict=True))]
    ref = score_batch(jobs)
    assert all(ref["legal"][i] == 0 and ref["iou_tc"][i] == 1.0 for i in swap)  # the precondition
    for name in reward_advantages.GATED:
        assert all(hacked[name][i] == 0.0 for i in swap)
        if name != "loc":  # the honest answers are hits: their location term is 0 too
            assert np.mean(hacked[name]) < np.mean(honest[name])
    assert hacked[R.AGGREGATOR] == advantages_v4(ref).tolist()


HIT, MISS = square(200, 200, 600, 600), square(700, 700, 900, 900)  # IoU_tc about 1 and 0 against cache "a"


def view_logs(groups, kinds=None, cells=None) -> dict:
    """log_metric output of one reward call over whole groups of 8 texts (all on cache "a")."""
    texts = [t for g in groups for t in g]
    n = len(texts)
    kw = dict(prompts=[[{"role": "user", "content": "p"}]] * n,
              completions=[[{"role": "assistant", "content": t}] for t in texts], completion_ids=[[5, EOS]] * n,
              key=["a"] * n, W=[W] * n, H=[W] * n, names=[list(NAMES)] * n)
    if kinds:
        kw |= {"kind": [k for k in kinds for _ in range(G)], "cell": [c for c in cells for _ in range(G)]}
    metrics = {}
    call_all(reward(), kw, log_metric=metrics.__setitem__)
    return metrics


def test_kind_and_cell_logs():
    """U44 (B1.5, reward v4): per-kind and per-cell IoU_tc and the share of all-miss groups per cell."""
    m = view_logs([[MISS] * G, [HIT] * G], kinds=["shifted", "centred"], cells=["21", "00"])
    hit = score_batch([(HIT, False, NAMES, CACHES["a"], W, W)])
    assert hit["iou_tc"][0] > 0.9 and hit["legal"][0] == 1
    assert m["kind/shifted/iou_tc"] == 0.0 and m["cell/21/iou_tc"] == 0.0
    assert m["kind/centred/iou_tc"] == pytest.approx(hit["iou_tc"][0]) == m["cell/00/iou_tc"]
    assert m["cell/21/allmiss_groups"] == 1.0 and m["cell/00/allmiss_groups"] == 0.0
    assert m["allmiss_share"] == 0.5


def test_allmiss_counts_whole_groups():
    m = view_logs([[MISS] * 7 + [HIT]], kinds=["shifted"], cells=["21"])
    assert m["cell/21/allmiss_groups"] == 0.0  # 7 misses + 1 hit is not an all-miss group
    m = view_logs([[MISS] * G, [MISS] * 7 + [HIT], [HIT] * G], kinds=["shifted"] * 3, cells=["21"] * 3)
    assert m["cell/21/allmiss_groups"] == pytest.approx(1 / 3) and m["allmiss_share"] == pytest.approx(1 / 3)


def test_no_kind_or_cell_logs_without_the_columns():
    m = view_logs([[MISS] * G])
    assert m["allmiss_share"] == 1.0 and not [k for k in m if k.startswith(("kind/", "cell/"))]


# CLEAN6 (CT-09): the reward without the Controller's record list; pre-CLEAN literals captured at 08e6455
ADVANTAGE_PRE_CLEAN = [1.6219269363403064, -2.92319357363354, -2.92319357363354, 1.6219269363403064,
                       -2.2632475344344516, 1.6219269363403064, 1.6219269363403064, 1.6219269363403064,
                       0.6228655995292476, 0.8391300385069551, 0.6228655995292476, -1.5441340302455195,
                       -0.46553721920509716, 0.8460784426014387, 0.6228655995292476, -1.5441340302455195]
IOU_PRE_CLEAN = [1.0, 0.0, 0.0, 1.0, 0.3025, 1.0, 1.0, 1.0, 0.37267080745341613, 0.4666666666666667,
                 0.37267080745341613, 0.0, 0.3025, 0.42857142857142855, 0.37267080745341613, 0.0]


def test_values_equal_the_pre_clean_literal_and_nothing_is_recorded():
    kw, _, _ = batch(truncate=(1,))
    rw = reward()
    size = len(vars(rw))
    got = call_all(rw, kw)
    assert got[R.AGGREGATOR] == ADVANTAGE_PRE_CLEAN and got["iou_tc"] == IOU_PRE_CLEAN  # bit-equal
    assert not hasattr(rw, "recorded")
    for _ in range(1000):  # a new completions list each time: a full reward computation per call
        call_all(rw, kw | {"completions": list(kw["completions"])})
    assert not hasattr(rw, "recorded") and len(vars(rw)) == size


def test_every_function_points_back_to_its_own_reward():
    import copy

    r1, r2 = reward(), reward()
    assert all(f.trl_reward is r1 for f in r1.funcs()[0]) and all(f.trl_reward is r2 for f in r2.funcs()[0])
    assert all(copy.copy(f).trl_reward is r1 for f in r1.funcs()[0])


def test_fastgrpo_prefetch_finds_the_reward_through_the_back_reference():
    """FastGRPO.prefetch_rewards reaches TRLReward.prefetch through f.trl_reward; without it prefetch would be off
    with no error."""
    pytest.importorskip("trl")
    from objdel.train.fast import FastGRPO

    kw, _, _ = batch()
    rw = reward()
    rows = [dict(zip(("key", "W", "H", "names"), c, strict=True))
            for c in zip(kw["key"], kw["W"], kw["H"], kw["names"], strict=True)]
    FastGRPO.prefetch_rewards(SimpleNamespace(reward_funcs=rw.funcs()[0], rows=rows), kw["completions"],
                              kw["completion_ids"])
    assert rw._ahead is not None
    assert call_all(rw, kw)[R.AGGREGATOR] == call_all(reward(), kw)[R.AGGREGATOR]


def test_by_prompt_is_the_iou_tc_mean_of_each_group():
    r = reward()
    texts = [HIT] * G + [MISS] * G
    n = len(texts)
    kw = dict(prompts=[[{"role": "user", "content": "p"}]] * n,
              completions=[[{"role": "assistant", "content": t}] for t in texts], completion_ids=[[5, EOS]] * n,
              key=["a"] * n, W=[W] * n, H=[W] * n, names=[list(NAMES)] * n)
    call_all(r, kw, log_metric=lambda *_: None)
    assert len(r.by_prompt) == 2 and r.by_prompt[0] > 0.9 and r.by_prompt[1] == 0.0
