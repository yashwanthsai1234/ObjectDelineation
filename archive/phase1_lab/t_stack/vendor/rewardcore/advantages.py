"""Group-advantage schemes for multi-term rewards (GRPO / GDPO variants).

Input: a batch of groups, terms R with shape (B, G, K) (B prompts, G rollouts, K reward terms
in [0, 1]), weights w (K,). Every scheme returns (A, C): advantages (B, G) and the per-term
contributions C (B, G, K) with A = sum_k C[..., k] (before any final batch rescaling, which
is applied to both).

Schemes
- grpo:        A = z_group(sum_k w_k r_k)                      (sum, then group-normalize)
- drgrpo:      A = sum_k w_k r_k - group mean                  (Dr. GRPO: no std division)
- gdpo:        A = z_group(sum_k w_k z_group(r_k))             (= objdel group_advantages, no gate)
- gdpo_batch:  A = z_batch(sum_k w_k z_group(r_k))             (GDPO paper: batch-level final z)
- gdpo_floor:  z_group with sigma_k <- max(sigma_k, floor_k), then batch-level final z
- gdpo_mean:   mean-only centring (scaled by a fixed floor_k) for the low-variance terms,
               z_group for the others, then batch-level final z
- gated:       r = legal * (w_iou IoU + w_band band) * (1 if topo else gamma), then
               (r - mean) / max(std, floor_total) per group (gates instead of additive terms)
- gdpo_floor_gated (recommended): legality is a hard gate (illegal -> every term 0), the topology
               term only counts when IoU >= iou_gate, then gdpo_floor
"""
from __future__ import annotations

import numpy as np

EPS = 1e-8


def _z(x: np.ndarray, axis: int = -1, floor: float | np.ndarray = 0.0) -> np.ndarray:
    mu = x.mean(axis=axis, keepdims=True)
    sd = x.std(axis=axis, keepdims=True)
    sd = np.maximum(sd, floor)
    return np.where(sd > EPS, (x - mu) / np.where(sd > EPS, sd, 1.0), 0.0)


def _batch_rescale(a: np.ndarray, c: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mu, sd = a.mean(), a.std()
    s = sd if sd > EPS else 1.0
    # centring is distributed over terms in proportion to their weight share (sum stays exact)
    return (a - mu) / s, (c - mu / c.shape[-1]) / s


def grpo(R, w, **_):
    s = R @ w  # (B, G)
    mu, sd = s.mean(1, keepdims=True), s.std(1, keepdims=True)
    scale = np.where(sd > EPS, sd, 1.0)
    a = (s - mu) / scale
    c = (R - R.mean(1, keepdims=True)) * w / scale[..., None]
    return a, c


def drgrpo(R, w, **_):
    c = (R - R.mean(1, keepdims=True)) * w
    return c.sum(-1), c


def gdpo(R, w, **_):
    c = _z(R, axis=1) * w
    s = c.sum(-1)
    mu, sd = s.mean(1, keepdims=True), s.std(1, keepdims=True)
    scale = np.where(sd > EPS, sd, 1.0)
    return (s - mu) / scale, (c - mu[..., None] / c.shape[-1]) / scale[..., None]


def gdpo_batch(R, w, **_):
    c = _z(R, axis=1) * w
    return _batch_rescale(c.sum(-1), c)


def gdpo_floor(R, w, floors, **_):
    c = _z(R, axis=1, floor=np.asarray(floors)) * w
    return _batch_rescale(c.sum(-1), c)


def gdpo_mean(R, w, floors, low_var, **_):
    z = _z(R, axis=1)
    centred = (R - R.mean(1, keepdims=True)) / np.asarray(floors)
    c = np.where(np.asarray(low_var), centred, z) * w
    return _batch_rescale(c.sum(-1), c)


def gated(R, w, names, floor_total=0.05, gamma=0.5, iou_gate=0.5, **_):
    i = {n: j for j, n in enumerate(names)}
    iou, band, topo, legal = (R[..., i[n]] for n in ("iou", "band", "topo", "legal"))
    base = (w[i["iou"]] * iou + w[i["band"]] * band) / (w[i["iou"]] + w[i["band"]])
    topo_ok = (topo >= 0.5) & (iou >= iou_gate)
    r = (legal >= 0.5) * base * np.where(topo_ok, 1.0, gamma)
    mu, sd = r.mean(1, keepdims=True), r.std(1, keepdims=True)
    a = (r - mu) / np.maximum(sd, floor_total)
    c = np.zeros(R.shape)
    c[..., i["iou"]] = a  # contributions are not separable under gates; attribute all to the base
    return a, c


def gdpo_floor_gated(R, w, floors, names, iou_gate=0.5, **_):
    i = {n: j for j, n in enumerate(names)}
    G = R * (R[..., i["legal"]] >= 0.5)[..., None]
    G[..., i["topo"]] *= G[..., i["iou"]] >= iou_gate
    return gdpo_floor(G, w, floors)


SCHEMES = {"grpo": grpo, "drgrpo": drgrpo, "gdpo": gdpo, "gdpo_batch": gdpo_batch,
           "gdpo_floor": gdpo_floor, "gdpo_mean": gdpo_mean, "gated": gated,
           "gdpo_floor_gated": gdpo_floor_gated}


# ---- iteration 2: anti-collapse topology term (t_dynamics) + topology-neutral geometry ----

def dm_loo_topology(sig_ids: np.ndarray, freq: np.ndarray, valid: np.ndarray, lam: float = 0.5,
                    floor: float = 1e-3) -> np.ndarray:
    """r = 1[valid] * (1 + lam * (log f(a) - log pi_loo(a))) per rollout, shape (B, G).

    sig_ids: integer id of each rollout's topology answer (equal ids = same answer within a group);
    freq: rater frequency f(a) of that answer; valid: answer matches a rater (and passes gates).
    pi_loo(a) = max(n(a) - 1, 0.5) / (G - 1), n(a) = rollouts in the group with the same answer.
    """
    sig_ids = np.asarray(sig_ids)
    G = sig_ids.shape[1]
    n = (sig_ids[:, :, None] == sig_ids[:, None, :]).sum(-1)
    pi = np.maximum(n - 1, 0.5) / (G - 1)
    return np.asarray(valid, float) * (1 + lam * (np.log(np.maximum(freq, floor)) - np.log(pi)))


def reward_v2(R, w, floors, names, r_topo, c_topo=1.0, **_):
    """Iteration-2 advantages.

    Geometric terms R (columns `names`, must include 'legal'): legality gate, then per-term group z
    with sigma floors, weighted sum, batch-level z (= gdpo_floor on the gated terms). The DM-loo
    topology reward r_topo (B, G) is legality-gated and only MEAN-CENTRED per group, then added
    with weight c_topo (advantage units): it is never divided by its own group std.
    """
    i = {n: j for j, n in enumerate(names)}
    gate = R[..., i["legal"]] >= 0.5
    a_geo, c_geo = gdpo_floor(R * gate[..., None], w, floors)
    t = np.asarray(r_topo, float) * gate
    a_topo = c_topo * (t - t.mean(1, keepdims=True))
    return a_geo + a_topo, np.concatenate([c_geo, a_topo[..., None]], -1)
