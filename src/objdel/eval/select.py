"""Choose one answer of K with no doctor masks (Implementation.md A4): the medoid under IoU, and the pixel vote.

- medoid: the legal answer with the highest mean IoU to the other legal answers (sample-based minimum-Bayes-risk
  decoding with utility IoU, Bertsch et al. 2023, https://arxiv.org/abs/2310.01387; the sample 1-median under the
  Jaccard distance). It reads only the answers' IoU matrix and legality, so it cannot read labels.
- vote_mask: the pixel strict majority of the legal masks, the exact minimiser of the summed Hamming distance (the
  pixel median; not the Jaccard medoid). With an even count, a pixel in exactly half of the masks is dropped. A mask
  is not an answer: deploying the vote needs the CUDA polygonizer (Phase D, kernels/K1).
"""

from collections.abc import Mapping, Sequence

import numpy as np


def medoid(sample_iou: np.ndarray, legal: Sequence[bool]) -> int:
    """Index of the legal answer with the highest mean IoU to the other legal answers; ties go to the lowest index;
    0 when no answer is legal. Legal empty answers have IoU 1 with each other, so a majority of empty answers wins:
    intended (the samples agree that there is nothing to outline). ValueError for a non-square matrix, K = 0, or one
    legal flag per answer missing."""
    M, legal = np.asarray(sample_iou, float), np.asarray(legal, bool)
    if M.ndim != 2 or M.shape[0] != M.shape[1] or len(M) == 0 or legal.shape != (len(M),):
        raise ValueError(f"medoid: need a K x K matrix (K >= 1) and K legal flags, got {M.shape} and {legal.shape}")
    eligible = np.flatnonzero(legal)
    if len(eligible) < 2:
        return int(eligible[0]) if len(eligible) else 0
    scores = [np.mean([M[i, j] for j in eligible if j != i]) for i in eligible]
    return int(eligible[np.argmax(scores)])


def vote_mask(masks: np.ndarray, legal: Sequence[bool]) -> np.ndarray:
    """Pixel strict majority (2 * count > n_legal) of the legal masks (K, H, W); all False when none is legal."""
    m = np.asarray(masks, bool)[np.asarray(legal, bool)]
    return 2 * m.sum(0) > len(m)


def sample_iou(objects: Mapping[str, dict]) -> np.ndarray:
    """The item's K x K answer IoU matrix: metrics' per-object `sample_iou`, averaged over the objects (RIGA)."""
    return np.mean([o["sample_iou"] for o in objects.values()], 0)
