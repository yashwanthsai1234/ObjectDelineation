"""Reward-v4 group advantages, for a batch ordered prompt by prompt.

    A = sum_k w_k z_k  +  c_topo (r_topo - mean_group r_topo)

- Legality is a hard gate: every term of an illegal (or unparseable) rollout is set to 0 before anything else.
- z_k is the per-group z-score of term k in (IoU_tc, F5, legal, g_area) with sigma_k <- max(sigma_k, floor_k).
  The floors keep a near-constant term from amplifying numeric noise (the GDPO "dark room"); there is no
  batch-level z-score, so one rollout's geometric advantage is bounded by sum_k w_k sqrt(G - 1) and
  c_topo = 2 dominates it.
- r_topo is rater-frequency matching of the topology answer with a group leave-one-out estimate of the
  policy's answer probability, an instance of MaD-RL's gated reverse-KL with the rater histogram as target:
      r_topo = valid * max(eps, 1 + lam (log f(a) - b(f(a), G) - log pi_loo(a))),
      pi_loo(a) = max(n(a) - 1, 0.5) / (G - 1),  n(a) = rollouts of the group with answer a,
      b(f, G) = log f - E_{N ~ Bin(G-1, f)} log(max(N, 0.5) / (G - 1)),
  valid = legal, a rater drew the answer and IoU_tc >= 0.5 for every object (iou_min; the eval pass
  criterion; terms without iou_min are single-object, where it equals iou_tc). b removes the estimator's bias,
  so at pi = f no rater answer is favoured; the clip eps keeps every valid rollout strictly above every invalid one. The
  term is only mean-centred (Dr. GRPO), never divided by its group std.
- A group whose terms all span < 1e-2 and whose rollouts give one topology answer gets advantage 0 (DAPO's
  dynamic sampling, by a spread threshold instead of exact equality).

Reward v5 (advantages_v5; Implementation.md B3.4) sums the terms into one reward and takes RLOO advantages:
    R = legal (IoU_tc + alpha_loc loc + w_f5 F5 + w_area g_area) + c_topo r_topo,
    A_i = R_i - mean_{j != i} R_j   (= G/(G-1) (R_i - mean_g R): group mean exactly 0, linear in R, no std).
loc is the scorer's gated location term (objdel.verify.reward.score.gated_loc). The drop rule tests each unweighted term
(legal, IoU_tc, F5, g_area, loc, r_topo) and one topology answer: a test on the weighted R would drop the far-miss
groups the location term is for.

References:
- GDPO, per-reward group normalization then a batch z-score (here: floors, no batch z):
  https://arxiv.org/abs/2601.05242 ; https://github.com/NVlabs/GDPO ; TRL multi_objective_aggregation:
  https://huggingface.co/docs/trl/main/en/grpo_trainer
- Dr. GRPO, no std / length normalization: https://arxiv.org/abs/2503.20783
- DAPO, dynamic sampling of zero-variance groups: https://arxiv.org/abs/2503.14476
- RLOO, leave-one-out baselines: https://arxiv.org/abs/2402.14740 ; Kool et al. 2019,
  https://openreview.net/forum?id=r1lgTGL5DE
- FlowRL, reward-distribution matching: https://arxiv.org/abs/2509.15207
- MaD-RL (Kulkarni et al. 2026): reverse-KL reward log P*(z) - log P_pi(z) (§2.3, Eq. 10) with finite-group
  frequency estimates, gated through correct rollouts R_correct (1 + lambda R_div) (§2.5, Eq. 17; checked in the
  full text, https://arxiv.org/html/2609.31644): https://arxiv.org/abs/2609.31644
- Bias of plug-in log-frequency estimators at finite group size (Jensen): Mohri, Schneider, Zhong 2026,
  https://arxiv.org/abs/2605.02435 ; b(f, G) here is the exact correction at pi = f for the LOO estimator.
- "The Dark Room in the Reward Channel" (z-scored dense rewards amplify noise): https://arxiv.org/abs/2607.21273
"""
from functools import cache
from math import comb, log

import numpy as np

from objdel.verify.reward.spec import REWARD_V4, REWARD_V5, TERMS

_INPUTS = (*TERMS, "topo_any", "topo_freq", "sig_id")


def _groups(terms: dict, group: int, extra: tuple[str, ...] = ()) -> dict[str, np.ndarray]:
    out = {k: np.asarray(terms[k], np.float64).reshape(-1) for k in (*_INPUTS, *extra)}
    out["iou_min"] = np.asarray(terms.get("iou_min", terms["iou_tc"]), np.float64).reshape(-1)
    n = {v.size for v in out.values()}
    if len(n) != 1 or n.pop() % group:
        raise ValueError(f"every term needs the same length, a multiple of group={group}")
    return {k: v.reshape(-1, group) for k, v in out.items()}


@cache
def loo_bias(f: float, group: int) -> float:
    """b(f, G) = log f - E[log(max(N, 0.5) / (G - 1))], N ~ Binomial(G - 1, f) (exact sum)."""
    n1 = group - 1
    expect = sum(comb(n1, n) * f**n * (1 - f) ** (n1 - n) * log(max(n, 0.5) / n1) for n in range(group))
    return log(f) - expect


def legal_gate(T: dict) -> np.ndarray:
    """Elementwise hard legality gate (module docstring): True for a legal rollout. The one rule by which
    advantages_v4, r_topo's validity and eval (objdel.verify.metrics: an illegal answer scores 0) treat illegal
    output."""
    return np.asarray(T["legal"]) >= 0.5


def topology_valid(T: dict, spec=REWARD_V4) -> np.ndarray:
    """Elementwise validity gate of r_topo (module docstring), also eval's pass criterion: legal, a rater drew the
    answer, IoU_tc >= iou_gate for every object (T["iou_min"])."""
    return legal_gate(T) & (T["topo_any"] > 0) & (T["iou_min"] >= spec.iou_gate) & (T["topo_freq"] > 0)


def _answers(T: dict) -> np.ndarray:
    """Topology answers with every illegal rollout mapped to one shared FAIL answer (as objdel.verify.metrics), so an
    illegal rollout's signature never moves the legal rollouts' r_topo or the flat-group check."""
    return np.where(legal_gate(T), T["sig_id"], -1.0)  # sig_ids are CRC32 values >= 0


def _topology(T: dict, group: int, spec=REWARD_V4) -> np.ndarray:
    s = _answers(T)
    pi = np.maximum((s[:, :, None] == s[:, None, :]).sum(-1) - 1, 0.5) / (group - 1)
    valid = topology_valid(T, spec)
    f = np.where(valid, T["topo_freq"], 1.0)
    u, inv = np.unique(f, return_inverse=True)
    b = np.array([loo_bias(float(x), group) for x in u])[inv].reshape(f.shape)
    bracket = np.maximum(spec.eps, 1 + spec.lam * (np.log(f) - b - np.log(pi)))
    return np.where(valid, bracket, 0.0)


def topology_reward(terms: dict, group: int = REWARD_V4.group, spec=REWARD_V4) -> np.ndarray:
    """r_topo (B*G,) of score_completion terms (module docstring); 0 for invalid or illegal rollouts."""
    return _topology(_groups(terms, group), group, spec).reshape(-1)


def advantages_v4(terms: dict, group: int = REWARD_V4.group) -> np.ndarray:
    """Reward-v4 advantages (B*G,) from per-rollout terms (keys of score_completion; flat, prompt by prompt)."""
    T = _groups(terms, group)
    R = np.stack([T[k] for k in TERMS], -1)  # (B, G, 4)
    X = np.where(legal_gate(T)[..., None], R, 0.0)  # where, not *: NaN * 0 is NaN
    z = (X - X.mean(1, keepdims=True)) / np.maximum(X.std(1, keepdims=True), REWARD_V4.floors)
    r = _topology(T, group)
    a = (z * REWARD_V4.weights).sum(-1) + REWARD_V4.c_topo * (r - r.mean(1, keepdims=True))
    flat = (np.ptp(R, axis=1).max(-1) < REWARD_V4.drop_spread) & (_answers(T) == _answers(T)[:, :1]).all(1)
    a[flat] = 0.0
    return a.reshape(-1)


def reward_parts_v5(T: dict, r_topo: np.ndarray, spec=REWARD_V5) -> dict[str, np.ndarray]:
    """The weighted summands of reward-v5 R (module docstring), by term; their sum is R."""
    area = spec.w_area * T["g_area"]
    if spec.area_gate:
        area = area * np.minimum(1.0, T["iou_tc"] / spec.miss_iou)
    legal = legal_gate(T)
    return {"iou": legal * T["iou_tc"], "loc": legal * spec.alpha_loc * T["loc"], "f5": legal * spec.w_f5 * T["f5"],
            "area": legal * area, "topo": spec.c_topo * r_topo}


def rewards_v5(T: dict, r_topo: np.ndarray, spec=REWARD_V5) -> np.ndarray:
    """Reward-v5 R of grouped terms T (any shape) and r_topo (module docstring)."""
    return sum(reward_parts_v5(T, r_topo, spec).values())


def advantages_v5(terms: dict, group: int = REWARD_V5.group, spec=REWARD_V5) -> np.ndarray:
    """Reward-v5 RLOO advantages (B*G,) from per-rollout terms (keys of score_completion; flat, prompt by prompt)."""
    T = _groups(terms, group, ("loc",))
    r = _topology(T, group, spec)
    R = rewards_v5(T, r, spec)
    a = R - (R.sum(1, keepdims=True) - R) / (group - 1)
    raw = np.stack([T["legal"], T["iou_tc"], T["f5"], T["g_area"], T["loc"], r], -1)
    flat = (np.ptp(raw, axis=1).max(-1) < spec.drop_spread) & (T["sig_id"] == T["sig_id"][:, :1]).all(1)
    a[flat] = 0.0
    return a.reshape(-1)


LOGGED = ("format_ok", "legal", "iou_tc", "f5", "g_area", "topo_any", "loc", "r_topo")

GATED = ("iou_tc", "f5", "g_area", "loc")  # times legal_gate: 0 for an illegal rollout, as the advantages see them


def logged_terms(T: dict) -> dict[str, np.ndarray]:
    """The LOGGED terms of score_batch's T but r_topo (which needs whole groups): the GATED ones times legal_gate."""
    legal = legal_gate(T)
    return {n: T[n] * legal if n in GATED else T[n] for n in LOGGED[:-1]}
