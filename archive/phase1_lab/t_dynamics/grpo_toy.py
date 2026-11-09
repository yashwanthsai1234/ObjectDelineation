"""Toy GRPO on a per-prompt categorical policy over discrete topology answers (beta0, beta1).

Everything is vectorised over (S seeds, P prompts, G group members). The policy for prompt p
in seed s is softmax(theta[s, p] + bias[s]) over K answers; `bias` is an optional shared
parameter (generalisation across prompts). Updates are exact GRPO on that tabular policy:
group-normalised advantages, PPO-clipped ratio with mu inner iterations against pi_old,
optional k3 KL estimator to a reference policy (DeepSeekMath form), plain SGD on logits.
"""
from dataclasses import dataclass, field

import numpy as np

# answer space: (b0, b1) with b0 in {0,1,2}, b1 in {0,1,2}
ANSWERS = [(b0, b1) for b0 in range(3) for b1 in range(3)]
K = len(ANSWERS)
IDX = {a: i for i, a in enumerate(ANSWERS)}


def hist(d):
    """rater histogram {answer: count} -> frequency vector over ANSWERS."""
    f = np.zeros(K)
    for a, c in d.items():
        f[IDX[a]] = c
    return f / f.sum()


# 4 raters per prompt (LIDC-like). Mixed prompts are the ones that matter.
SUITE = {
    "unanimous": {(1, 0): 4},
    "hole 2/4": {(1, 0): 2, (1, 1): 2},
    "hole 1/4": {(1, 0): 3, (1, 1): 1},
    "hole 3/4": {(1, 0): 1, (1, 1): 3},
    "3-way 2/1/1": {(1, 0): 2, (1, 1): 1, (2, 0): 1},
    "split 1/4": {(1, 0): 3, (2, 0): 1},
}


def sft_policy(freq, invalid_mass=0.15, sharpen=1.0):
    """SFT init: rater frequencies (optionally sharpened, i.e. under-producing minority modes)
    mixed with uniform mass over all K answers (so some invalid mass exists)."""
    f = np.asarray(freq, float) ** sharpen
    f = f / f.sum(-1, keepdims=True)
    return (1 - invalid_mass) * f + invalid_mass / K


def softmax(x):
    z = x - x.max(-1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(-1, keepdims=True)


# ---------------------------------------------------------------- advantages

def standardize(v, axis=-1, mode="std", eps=0.0):
    """Group normalisation along `axis`.

    mode='std'  : (v-mean)/std, exactly like the repo's annotators.standardize when eps=0
                  (std==0 -> divide by 1, i.e. all zeros).
    mode='mean' : Dr. GRPO, mean-centring only.
    """
    v = np.asarray(v, float)
    c = v - v.mean(axis, keepdims=True)
    if mode == "mean":
        return c
    sd = v.std(axis, keepdims=True)
    if eps > 0:
        return c / (sd + eps)
    return c / np.where(sd > 0, sd, 1.0)


def gdpo_advantages(iou, band, topo, w_iou=0.3, w_band=0.3, w_topo=0.4, gate=0.5):
    """Vectorised replica of objdel.reward.annotators.group_advantages (last axis = group)."""
    t = np.where(iou >= gate, topo, 0.0)
    mix = w_iou * standardize(iou) + w_band * standardize(band) + w_topo * standardize(t)
    return standardize(mix)


# ---------------------------------------------------------------- topology rewards

def counts_in_group(a):
    """n[s,p,g] = number of group members with the same answer as member g."""
    oh = np.eye(K)[a]  # S,P,G,K
    n = oh.sum(-2)  # S,P,K
    return np.take_along_axis(n, a, -1), n


def topo_reward(name, a, freq, logp_old, lam=0.5, floor=1e-3):
    """Topology reward per rollout. a: (S,P,G) answer indices; freq: (P,K) rater freq."""
    S, P, G = a.shape
    f = np.broadcast_to(freq, (S, P, K))
    fa = np.take_along_axis(f, a, -1)
    valid = (fa > 0).astype(float)
    if name == "set":  # any rater counts
        return valid
    if name == "freq":  # share of raters (annotators.topo_freq)
        return fa
    if name == "majority":  # only the modal answer(s)
        return (fa >= f.max(-1, keepdims=True) - 1e-12).astype(float)
    if name == "near":  # repo default annotators.topo_near: max_r 1/(1+|db0|+|db1|)
        ans = np.array(ANSWERS)
        d = np.abs(ans[:, None, :] - ans[None, :, :]).sum(-1)  # K,K
        near = 1.0 / (1 + d)  # near[a, r]
        tab = np.where(np.broadcast_to(freq, (P, K))[:, None, :] > 0, near[None], 0).max(-1)  # P,K
        return np.take_along_axis(np.broadcast_to(tab, (S, P, K)), a, -1)
    if name == "logfreq":  # plain log rater frequency (no entropy term)
        return np.log(np.maximum(fa, floor))
    if name == "cov":  # set + split coverage credit: lam/n(a) per valid member (sum = lam*#distinct)
        n, _ = counts_in_group(a)
        return valid * (1 + lam / n)
    if name == "uniq":  # set + leave-one-out marginal coverage: lam if the only one with this answer
        n, _ = counts_in_group(a)
        return valid * (1 + lam * (n == 1))
    if name == "share":  # set + fitness sharing toward rater freq: lam*(f(a)*G/n(a) - 1), capped
        n, _ = counts_in_group(a)
        return valid * (1 + lam * (np.minimum(fa * G / n, 4.0) - 1))
    if name == "dm":  # FlowRL-style reverse-KL distribution matching: log f~ - log pi_old
        return np.log(np.maximum(fa, floor)) - logp_old
    if name == "set+dm":  # set reward + lam*(log f - log pi_old) on valid answers only
        return valid * (1 + lam * (np.log(np.maximum(fa, floor)) - logp_old))
    if name == "set+dm_emp":  # deployable: pi(topology) estimated by the group's own counts n(a)/G
        n, _ = counts_in_group(a)
        return valid * (1 + lam * (np.log(np.maximum(fa, floor)) - np.log(n / G)))
    if name == "set+dm_loo":  # leave-one-out count estimate (n-1)/(G-1), floored at 0.5/(G-1)
        n, _ = counts_in_group(a)
        return valid * (1 + lam * (np.log(np.maximum(fa, floor)) - np.log(np.maximum(n - 1, 0.5) / (G - 1))))
    raise ValueError(name)


# ---------------------------------------------------------------- GRPO step

@dataclass
class Cfg:
    reward: str = "set"
    lam: float = 0.5
    beta: float = 0.0  # KL coefficient
    kl: str = "k3"  # 'k3' (sampled, DeepSeekMath) or 'exact'
    G: int = 8
    lr: float = 0.5
    mu: int = 2  # inner PPO iterations per batch
    eps_low: float = 0.2
    eps_high: float = 0.28
    adv: str = "std"  # 'std' (GRPO) or 'mean' (Dr. GRPO)
    adv_eps: float = 0.0
    # continuous term
    combine: str = "none"  # 'none' | 'sum' | 'gdpo'
    w_iou: float = 0.5  # for 'sum': r = w_iou*iou + (1-w_iou)*topo
    iou_mu: np.ndarray = field(default_factory=lambda: None)  # (P,K) mean IoU per answer
    iou_sd: float = 0.05
    gate: float = 0.5
    shared_lr: float = 0.0  # lr of shared bias (0 = no shared parameter)


def grad_surrogate(theta_tot, a, A, logp_old, ref, cfg):
    """d/dlogits of (1/G) sum_g [min(rho A, clip(rho) A) - beta*k3]. Returns (S,P,K)."""
    pi = softmax(theta_tot)
    logpi = np.log(pi)
    lp = np.take_along_axis(logpi, a, -1)
    rho = np.exp(lp - logp_old)
    flows = ~(((A > 0) & (rho > 1 + cfg.eps_high)) | ((A < 0) & (rho < 1 - cfg.eps_low)))
    coef = np.where(flows, A * rho, 0.0)  # d(rho*A)/dlogpi(a) = rho*A
    if cfg.beta > 0 and cfg.kl == "k3":
        lref = np.broadcast_to(np.log(ref), theta_tot.shape)
        r = np.exp(np.take_along_axis(lref, a, -1) - lp)
        coef = coef - cfg.beta * (1 - r)
    oh = np.eye(K)[a]  # S,P,G,K
    G = a.shape[-1]
    # sum_g coef_g (e_{a_g} - pi) / G
    g = (coef[..., None] * oh).sum(-2) / G - coef.sum(-1, keepdims=True) / G * pi
    if cfg.beta > 0 and cfg.kl == "exact":
        lr_ = logpi - np.log(ref)
        kl = (pi * lr_).sum(-1, keepdims=True)
        g = g - cfg.beta * pi * (lr_ - kl)
    return g


def advantages(cfg, a, freq, logp_old, rng):
    topo = topo_reward(cfg.reward, a, freq, logp_old, lam=cfg.lam)
    if cfg.combine == "none":
        return standardize(topo, mode=cfg.adv, eps=cfg.adv_eps), topo, None
    mu = np.take_along_axis(np.broadcast_to(cfg.iou_mu, a.shape[:2] + (K,)), a, -1)
    iou = np.clip(mu + cfg.iou_sd * rng.standard_normal(a.shape), 0, 1)
    if cfg.combine == "sum":
        r = cfg.w_iou * iou + (1 - cfg.w_iou) * topo
        return standardize(r, mode=cfg.adv, eps=cfg.adv_eps), topo, iou
    if cfg.combine == "gdpo":
        band = np.clip(iou + 0.05 + 0.02 * rng.standard_normal(a.shape), 0, 1)
        return gdpo_advantages(iou, band, topo, gate=cfg.gate), topo, iou
    raise ValueError(cfg.combine)


def sample(pi, G, rng):
    c = np.cumsum(pi, -1)
    u = rng.random(pi.shape[:-1] + (G,))
    return np.minimum((u[..., None] > c[..., None, :]).sum(-1), K - 1)


def step(theta, bias, freq, ref, cfg, rng):
    """One GRPO batch per (seed, prompt). Mutates theta/bias in place."""
    pi_old = softmax(theta + bias[:, None, :])
    a = sample(pi_old, cfg.G, rng)
    logp_old = np.take_along_axis(np.log(pi_old), a, -1)
    A, topo, iou = advantages(cfg, a, freq, logp_old, rng)
    for _ in range(cfg.mu):
        g = grad_surrogate(theta + bias[:, None, :], a, A, logp_old, ref, cfg)
        theta += cfg.lr * g
        if cfg.shared_lr > 0:
            bias += cfg.shared_lr * g.mean(1)
    return a, A, topo


# ---------------------------------------------------------------- metrics

def metrics(pi, freq):
    f = np.broadcast_to(freq, pi.shape)
    invalid = (pi * (f == 0)).sum(-1)
    tv = 0.5 * np.abs(pi - f).sum(-1)
    m = 0.5 * (pi + f)

    def kl(p, q):
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(p > 0, p * np.log2(p / q), 0.0)
        return t.sum(-1)

    js = 0.5 * kl(pi, m) + 0.5 * kl(f, m)
    ent = -(pi * np.log(np.maximum(pi, 1e-300))).sum(-1)
    # minority retention: min over valid answers of pi(a)/f(a), restricted to the valid simplex
    pv = pi * (f > 0)
    pv = pv / np.maximum(pv.sum(-1, keepdims=True), 1e-300)
    ratio = np.where(f > 0, pv / np.where(f > 0, f, 1), np.inf).min(-1)
    return dict(invalid=invalid, tv=tv, js=js, ent=ent, minret=ratio)


def run(cfg, suite=None, S=20, T=300, seed=0, invalid_mass=0.15, sharpen=1.0, ref_kind="sft",
        log_every=5):
    suite = SUITE if suite is None else suite
    names = list(suite)
    freq = np.stack([hist(suite[n]) for n in names])  # P,K
    sft = sft_policy(freq, invalid_mass, sharpen)
    ref = sft if ref_kind == "sft" else sft_policy(freq, 0.01, 1.0)
    rng = np.random.default_rng(seed)
    theta = np.broadcast_to(np.log(sft), (S,) + sft.shape).copy()
    bias = np.zeros((S, K))
    logs = {k: [] for k in ["invalid", "tv", "js", "ent", "minret"]}
    steps = []
    for t in range(T + 1):
        if t % log_every == 0 or t == T:
            m = metrics(softmax(theta + bias[:, None, :]), freq)
            for k in logs:
                logs[k].append(m[k])
            steps.append(t)
        if t < T:
            step(theta, bias, freq, ref, cfg, rng)
    out = {k: np.stack(v, 0) for k, v in logs.items()}  # (n_log, S, P)
    out["steps"] = np.array(steps)
    out["names"] = names
    out["freq"] = freq
    out["pi"] = softmax(theta + bias[:, None, :])
    return out
