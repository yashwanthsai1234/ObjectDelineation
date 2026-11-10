"""objdel.verify.select (Implementation.md A4): the medoid answer and the pixel vote, chosen with no doctor masks."""

import inspect

import numpy as np
import pytest
from test_eval import disk, text_of

from objdel.verify.metrics import metrics
from objdel.verify.reward.cache import ExpertCache
from objdel.verify.select import medoid, vote_mask

M3 = np.array([[1, .5, .5], [.5, 1, .9], [.5, .9, 1]])


def test_medoid_hand_cases():
    """U13."""
    assert medoid(M3, [1, 1, 1]) == 1  # scores .5, .7, .7: a tie goes to the lowest index
    assert medoid(M3, [0, 0, 0]) == 0
    assert (medoid(M3, [0, 1, 0]), medoid(M3, [0, 0, 1]), medoid(M3, [1, 0, 1])) == (1, 2, 0)
    for diag in (0.0, 1.0, 7.0):  # the diagonal never counts (j != i)
        m = M3.copy()
        np.fill_diagonal(m, diag)
        assert medoid(m, [1, 1, 1]) == 1
    # three legal empty answers (IoU 1 among them) win over two distinct non-empty ones: the first empty one
    m5 = np.array([[1, 0, 0, 0, 0], [0, 1, .3, 0, 0], [0, .3, 1, 0, 0], [0, 0, 0, 1, 1], [0, 0, 0, 1, 1]], float)
    m5[0, 3] = m5[3, 0] = m5[0, 4] = m5[4, 0] = 1.0
    assert medoid(m5, [1] * 5) == 0
    assert medoid(np.array([[1, .2], [.2, 1]]), [1, 1]) == 0
    # the mean over eligible j != i only: over all j the answer 0 wins (illegal answer 3 lifts it), over the legal 2
    m4 = np.array([[1, .3, .4, 1], [.3, 1, .5, 0], [.4, .5, 1, 0], [1, 0, 0, 1]])
    assert np.argmax((m4.sum(1) - 1) / 3) == 0
    assert medoid(m4, [1, 1, 1, 0]) == 2


def test_medoid_contract():
    """U14: no rater input; a bad shape raises; K = 0 raises (pin)."""
    assert list(inspect.signature(medoid).parameters) == ["sample_iou", "legal"]
    with pytest.raises(ValueError):
        medoid(np.ones((2, 3)), [1, 1])
    with pytest.raises(ValueError):
        medoid(M3, [1, 1])
    with pytest.raises(ValueError):
        medoid(np.ones((0, 0)), [])


def test_vote_mask_hand_cases():
    """U15: strict majority of the legal masks."""
    rng = np.random.default_rng(0)
    m = rng.random((4, 6, 5)) < 0.5
    assert np.array_equal(vote_mask(m[:3], [1, 1, 1]), m[:3].sum(0) >= 2)
    assert np.array_equal(vote_mask(m, [1, 1, 1, 1]), m.sum(0) >= 3)  # a 2-2 tie is False
    same = np.stack([m[0], m[0], m[2], m[3]])
    assert np.array_equal(vote_mask(same, [1, 1, 0, 0]), m[0])
    assert np.array_equal(vote_mask(m, [1, 1, 0, 0]), m[0] & m[1])
    none = vote_mask(m, [0, 0, 0, 0])
    assert none.shape == (6, 5) and none.dtype == bool and not none.any()
    assert np.array_equal(vote_mask(m, [0, 0, 1, 0]), m[2])
    full = np.concatenate([m[:3], np.ones((5, 6, 5), bool)])
    assert np.array_equal(vote_mask(full, [1, 1, 1, 0, 0, 0, 0, 0]), m[:3].sum(0) >= 2)
    assert not np.array_equal(vote_mask(full, [1] * 8), m[:3].sum(0) >= 2)


def oracle(M, legal):
    eligible = [i for i in range(len(M)) if legal[i]]
    if not eligible:
        return 0
    if len(eligible) == 1:
        return eligible[0]
    score = {i: np.mean([M[i, j] for j in eligible if j != i]) for i in eligible}
    return max(eligible, key=lambda i: (score[i], -i))


def test_medoid_equals_a_brute_force_oracle():
    """P06: 2,000 cases on a 0.25 grid (ties), and permutation equivariance on float matrices."""
    rng = np.random.default_rng(6)
    for _ in range(2000):
        K = int(rng.integers(1, 9))
        a = rng.integers(0, 5, (K, K)) / 4
        M = np.triu(a) + np.triu(a, 1).T
        legal = rng.random(K) < 0.7
        assert medoid(M, legal) == oracle(M, legal), (M, legal)
    for _ in range(300):
        K = int(rng.integers(2, 9))
        a = rng.random((K, K))
        M = np.triu(a) + np.triu(a, 1).T
        legal = rng.random(K) < 0.7
        if legal.sum() < 2:  # none legal returns 0, not a permuted index
            continue
        p = rng.permutation(K)
        i, j = medoid(M, legal), p[medoid(M[np.ix_(p, p)], legal[p])]
        e = np.flatnonzero(legal)
        score = {a: np.mean([M[a, b] for b in e if b != a]) for a in e}
        assert i == j if len(set(score.values())) == len(score) else score.get(i) == score.get(j)


def test_vote_mask_invariants():
    """P07."""
    rng = np.random.default_rng(7)
    m = rng.random((5, 7)) < 0.4
    assert np.array_equal(vote_mask(np.stack([m, m, m]), [1, 1, 1]), m)
    assert not vote_mask(np.stack([m, ~m]), [1, 1]).any()
    for _ in range(300):
        K, H, W = (int(x) for x in rng.integers(1, 10, 3))
        K = min(K, 8)
        ms = rng.random((K, H, W)) < rng.random()
        legal = rng.random(K) < 0.7
        before = ms.copy()
        v = vote_mask(ms, legal)
        assert np.array_equal(ms, before)
        assert np.array_equal(v, ms[legal].sum(0) * 2 > legal.sum())
        if legal.any():  # between the intersection and the union of the legal masks
            assert (ms[legal].all(0) <= v).all() and (v <= ms[legal].any(0)).all()
        p = rng.permutation(K)
        assert np.array_equal(vote_mask(ms[p], legal[p]), v)
    raters = rng.random((3, 9, 9)) < 0.5
    raters[:, 4, 4] = True
    assert np.array_equal(vote_mask(raters, [1, 1, 1]), ExpertCache.build(raters).vote)


def test_rater_data_never_reaches_the_medoid():
    """P14: permuted or replaced rater masks leave sample_iou, the medoid and the vote answer's IoU_tc as they are
    (the vote answer is scored against the doctors, so only the permutation keeps it)."""
    S = 50
    raters = np.stack([disk(25, 25, 10), disk(25, 27, 12), disk(24, 25, 9)])
    texts = [text_of({"nodule": m}, ("nodule",), {"nodule": raters}) for m in (*raters, disk(20, 20, 6))] + ["x"]
    texts = texts * 2
    rng = np.random.default_rng(14)
    m = metrics(texts, {"nodule": raters}, S, S)
    a, legal = m["objects"]["nodule"], m["samples"]["legal"]
    b = metrics(texts, {"nodule": raters[[2, 0, 1]]}, S, S)["objects"]["nodule"]
    noise = np.stack([disk(*rng.integers(15, 35, 2), int(rng.integers(5, 12))) for _ in range(3)])
    c = metrics(texts, {"nodule": noise}, S, S)["objects"]["nodule"]
    assert a["sample_iou"] == b["sample_iou"] == c["sample_iou"]
    assert medoid(np.array(a["sample_iou"]), legal) == medoid(np.array(c["sample_iou"]), legal)
    assert a["pixel_vote_iou_tc"] == b["pixel_vote_iou_tc"]
