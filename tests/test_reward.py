import numpy as np

from objdel.reward import AnnotatorSet, group_advantages


def make(W, H, ops):
    g = np.zeros((H, W), bool)
    for v, x0, y0, x1, y1 in ops:
        g[max(y0, 0):y1, max(x0, 0):x1] = bool(v)
    return g


def test_reproduces_artifact_grpo_demo():
    W, H = 22, 14
    annot = [make(W, H, [[1, 3, 2, 19, 12], [0, 8, 5, 13, 9]]),
             make(W, H, [[1, 3, 2, 18, 12], [0, 8, 5, 14, 10]]),
             make(W, H, [[1, 4, 3, 19, 12], [0, 9, 5, 13, 9]])]
    samples = {
        "faithful": [[1, 4, 2, 19, 12], [0, 8, 5, 13, 9]],
        "hole_filled": [[1, 3, 2, 19, 12]],
        "extra_speck": [[1, 4, 2, 19, 12], [0, 8, 5, 13, 9], [1, 0, 0, 2, 1]],
        "over_dilated": [[1, 2, 1, 20, 13], [0, 9, 6, 12, 8]],
        "broken_ring": [[1, 4, 2, 19, 12], [0, 8, 5, 13, 9], [0, 15, 2, 16, 12]],
        "slightly_thin": [[1, 4, 3, 18, 11], [0, 8, 5, 13, 9]],
    }
    aset = AnnotatorSet(annot)
    scores = [aset.score(make(W, H, ops)) for ops in samples.values()]
    adv_iou = group_advantages(scores, 1, 0, 0)
    adv = group_advantages(scores)
    expect_iou = [1.01, 0.52, 0.89, -1.38, 0.35, -1.38]
    expect = [1.94, -0.45, 0.16, -1.36, -0.36, 0.08]
    assert np.allclose(adv_iou, expect_iou, atol=0.01), adv_iou
    assert np.allclose(adv, expect, atol=0.01), adv
    assert scores[1].errors["missing_holes"] == 1
    assert scores[2].errors["extra_objects"] == 1
