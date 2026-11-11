"""REWARD_V4: the frozen constants of reward v4;
CACHE_SPEC: the constants stored ExpertCache values (tau, T5) depend on.

Changing a REWARD_V4 value is a change of the reward definition, not of its implementation: it needs a new round of
tuning, but it leaves stored caches valid. Changing a CACHE_SPEC value invalidates every stored ExpertCache (their
fingerprint covers CACHE_SPEC, not the reward weights).

References:
- Frozen dataclasses, astuple for the fingerprint: https://docs.python.org/3/library/dataclasses.html ; SHA-256:
  https://docs.python.org/3/library/hashlib.html
- Library versions in the fingerprint (caches rebuilt when a pin changes): numpy.__version__,
  scipy.__version__, cv2.__version__ (https://numpy.org/doc/stable/reference/routines.version.html).
"""
import hashlib
from dataclasses import astuple, dataclass

import cv2
import numpy as np
import scipy

from objdel.geom.topology import min_area

TERMS = ("iou_tc", "f5", "legal", "g_area")  # the z-scored terms, in weight / floor order


@dataclass(frozen=True)
class CacheSpec:
    """Parameters of the expert cache: expected expert F5, tau floor (pixels) and tau percentile."""
    expert_score: float = 0.95  # s*: expected F5 of a new expert's mask
    tau_floor: float = 1.0  # px, floor of tau and of the F5 tolerance
    tau_quantile: float = 95.0  # tau = this percentile of pooled rater-pair contour distances


CACHE_SPEC = CacheSpec()


@dataclass(frozen=True)
class RewardV4:
    """Parameters of reward v4: term weights, sigma floors, topology bracket, group size and gates."""
    weights: tuple[float, ...] = (0.3, 0.3, 0.1, 0.3)  # IoU_tc, F5, legal, g_area
    floors: tuple[float, ...] = (0.05, 0.05, 0.25, 0.05)  # sigma floors of the same terms
    c_topo: float = 2.0  # weight of the mean-centred topology reward (advantage units)
    lam: float = 0.5  # lambda of the topology bracket
    eps: float = 0.1  # clip floor of the topology bracket: a valid rollout always beats an invalid one
    group: int = 8  # G, rollouts per prompt
    drop_spread: float = 1e-2  # groups whose terms all span less than this (and one topology answer) get 0
    iou_gate: float = 0.5  # a rater topology counts only with IoU_tc >= this
    name: str = "reward_v4"


REWARD_V4 = RewardV4()


@dataclass(frozen=True)
class RewardV5:
    """Reward v5 (Implementation.md B3.3): R = legal (IoU_tc + alpha_loc loc + w_f5 F5 + w_area g_area) + c_topo
    r_topo, RLOO advantages (advantages_v5). alpha_loc < miss_iou keeps a hit strictly above a miss with full
    location for equal F5 and g_area (math-stats notes s11 item 2: at alpha_loc = miss_iou they tie)."""
    c_topo: float  # weight of r_topo (reward units); B3.7 rule 8 sets it
    name: str = "reward_v5"
    alpha_loc: float = 0.09  # weight of the gated location term; < miss_iou
    miss_iou: float = 0.1  # a missed answer: IoU_tc below this; the location ramp ends here
    w_f5: float = 0.1
    w_area: float = 0.1
    loc_width: float = 1.0  # Gaussian kernel width in units of the union's radius (B3.7 rule 5)
    area_gate: bool = False  # True scales the area term by min(1, IoU_tc / miss_iou) (B3.7 rule 4)
    lam: float = 0.5
    eps: float = 0.1
    group: int = 8
    drop_spread: float = 1e-2
    iou_gate: float = 0.5

    def __post_init__(self):
        if self.alpha_loc >= self.miss_iou:
            raise ValueError(f"alpha_loc {self.alpha_loc} must be < miss_iou {self.miss_iou}: a miss with full "
                             f"location would tie or beat a hit at IoU_tc = miss_iou")


REWARD_V5 = RewardV5(c_topo=1.0)  # provisional c_topo (c * eps >= a 0.1 margin) until the B3.7 lab round
REWARDS = {"v4": REWARD_V4, "v5": REWARD_V5}  # train_rl --reward


def reward_fingerprint(spec) -> str:
    """SHA-256 of a reward spec's constants (rl.json provenance)."""
    return hashlib.sha256(repr(astuple(spec)).encode()).hexdigest()


def fingerprint() -> str:
    """SHA-256 of everything a stored ExpertCache value depends on besides the rater masks: CACHE_SPEC, the size rule
    (probed at fixed union areas) and the versions of the libraries that compute tau and T5: numpy (percentile,
    sums), OpenCV (contour erosion) and scipy (quad and the normal law in S_m beyond the frozen m <= 8; until the
    engine took over, also the distance transforms: its bit-identical EDT replaced them without a fingerprint
    change).
    Caches are rebuilt when a pin changes. Signatures, the vote and the band are rederived on load and need
    no fingerprint; reward weights are not stored, so they are not in it."""
    parts = (astuple(CACHE_SPEC), [min_area(a) for a in (0, 999, 1001, 10**4, 10**6)], np.__version__, cv2.__version__,
             scipy.__version__)
    return hashlib.sha256(repr(parts).encode()).hexdigest()
