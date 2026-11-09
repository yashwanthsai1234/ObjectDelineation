"""objdel.train.rewards.TRLReward without a model: terms, advantages, truncation, group checks, typed logs."""
from types import SimpleNamespace

import numpy as np
import pytest

from objdel._gpucheck import DEFECTS
from objdel.reward import ExpertCache, advantages_v4, legal_gate, score_batch, topology_reward
from objdel.train import rewards as R

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
    assert [f.__name__ for f in funcs] == [*R.LOGGED, R.AGGREGATOR]
    assert weights == [0.0] * len(R.LOGGED) + [1.0]


def test_values_equal_score_batch_and_advantages_v4():
    kw, texts, keys = batch(truncate=(1,))
    got = call_all(reward(), kw)
    jobs = [(t, i == 1, NAMES, CACHES[k], W, W) for i, (t, k) in enumerate(zip(texts, keys, strict=True))]
    ref = score_batch(jobs)
    legal = legal_gate(ref)
    for name in ("format_ok", "legal", "topo_any"):
        assert got[name] == ref[name].tolist()
    for name in R.GATED:  # the geometric terms as advantages_v4 sees them: 0 for an illegal rollout
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
    for name in R.GATED:
        assert all(hacked[name][i] == 0.0 for i in swap)
        assert np.mean(hacked[name]) < np.mean(honest[name])
    assert hacked[R.AGGREGATOR] == advantages_v4(ref).tolist()
