"""REWARD_V4: the frozen constants of reward v4 (plan §4.3; lab/t_reward/results/ITERATION4.md).

Changing any value is a change of the reward definition, not of its implementation: it needs a new lab
round, and it invalidates every stored ExpertCache (their fingerprint includes these values).

References:
- Values: experiments/phase1/PLAN.md §4.3 and experiments/phase1/lab/t_reward/results/ITERATION4.md (weights, floors,
  c_topo, lam, eps, iou_gate), ITERATION3.md (debias, clip, no batch z).
- Frozen dataclass, astuple for the fingerprint: https://docs.python.org/3/library/dataclasses.html ; SHA-256:
  https://docs.python.org/3/library/hashlib.html
- Library versions in the fingerprint (plan §6, caches rebuilt when a pin changes): numpy.__version__,
  scipy.__version__, cv2.__version__ (https://numpy.org/doc/stable/reference/routines.version.html).
"""
import hashlib
from dataclasses import astuple, dataclass

import cv2
import numpy as np
import scipy

from objdel.checker.topology import min_area

TERMS = ("iou_tc", "f5", "legal", "g_area")  # the z-scored terms, in weight / floor order


@dataclass(frozen=True)
class RewardV4:
    weights: tuple[float, ...] = (0.3, 0.3, 0.1, 0.3)  # IoU_tc, F5, legal, g_area
    floors: tuple[float, ...] = (0.05, 0.05, 0.25, 0.05)  # sigma floors of the same terms
    c_topo: float = 2.0  # weight of the mean-centred topology reward (advantage units)
    lam: float = 0.5  # lambda of the topology bracket
    eps: float = 0.1  # clip floor of the topology bracket: a valid rollout always beats an invalid one
    group: int = 8  # G, rollouts per prompt
    drop_spread: float = 1e-2  # groups whose terms all span less than this (and one topology answer) get 0
    iou_gate: float = 0.5  # a rater topology counts only with IoU_tc >= this
    expert_score: float = 0.95  # s*: expected F5 of a new expert's mask
    tau_floor: float = 1.0  # px, floor of tau and of the F5 tolerance
    tau_quantile: float = 95.0  # tau = this percentile of pooled rater-pair contour distances


REWARD_V4 = RewardV4()


def fingerprint() -> str:
    """SHA-256 of everything a stored ExpertCache value depends on besides the rater masks: REWARD_V4, the size rule
    (probed at fixed union areas) and the versions of the libraries that compute tau and T5: numpy (percentile,
    sums), OpenCV (contour erosion) and scipy (distance transforms; quad and the normal law in S_m beyond the frozen
    m <= 8). Plan §6: caches are rebuilt when a pin changes. Signatures, the vote and the band are rederived on load
    and need no fingerprint."""
    parts = (astuple(REWARD_V4), [min_area(a) for a in (0, 999, 1001, 10**4, 10**6)], np.__version__, cv2.__version__,
             scipy.__version__)
    return hashlib.sha256(repr(parts).encode()).hexdigest()
