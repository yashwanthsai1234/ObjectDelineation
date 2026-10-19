"""Question B: synthetic RLVR curriculum simulator with typed checker errors.

Latent model (per run):
  * K error types: extra_hole, missing_hole, extra_part, missing_part, boundary.
  * Each item i has, per type k, a relevance flag (prevalence q_k) and a difficulty d_ik; for
    irrelevant types d_ik = -8 (error essentially impossible).
  * Policy state: a global skill s_k per type (transfers to held-out), a per-item memorisation m_i
    (train items only, does NOT transfer: overfitting), and an optional hack variable H.
  * p_ik = floor_k*rel_ik + (1 - floor_k*rel_ik) * sigmoid(d_ik - s_k - m_i)   (EP also gets the hack term)
  * A GRPO group of G rollouts draws typed error indicators e_gk ~ Bernoulli(p_ik); the scalar
    train reward is r_g = 1 - sum_k w_k e_gk; advantages are group z-scores.
  * The learning signal for type k from item i is c_ik = max(0, -mean_g A_g e_gk): GRPO only
    learns to fix an error type when, inside the group, rollouts with that error got lower advantage.
  * s_k += alpha_k * mean_i(c_ik) * exp(-s_k / s_cap) + transfer; s_k *= (1 - forget) every step.
  * m_i += alpha_m * sum_k c_ik on every visit.
Samplers choose which B items of the train set form each batch.
"""
from dataclasses import dataclass, field, replace

import numpy as np

TYPES = ["extra_hole", "missing_hole", "extra_part", "missing_part", "boundary"]
K = len(TYPES)
EP, BD = 2, 4


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


@dataclass
class World:
    n_train: int = 512
    n_val: int = 128
    n_test: int = 256
    prevalence: tuple = (0.12, 0.12, 0.08, 0.12, 0.70)
    d_mu: tuple = (0.5, 1.0, 0.5, 0.5, 0.0)
    d_sd: float = 1.0
    alpha: tuple = (0.4, 0.12, 0.4, 0.4, 0.4)  # learnability per type
    s_cap: float = 4.0  # diminishing returns scale
    floor: tuple = (0, 0, 0, 0, 0)  # irreducible typed error on relevant items (rater ambiguity)
    forget: float = 0.001
    alpha_m: float = 0.03  # memorisation per visit
    transfer: float = 0.0  # share of other types' skill gain that transfers
    hack: bool = False
    alpha_h: float = 1.5
    kappa: float = 1.0
    w: tuple = (0.2, 0.2, 0.2, 0.2, 0.2)  # train reward weights
    w_ep_hack: float = 0.02  # EP weight in train reward when the hack scenario is on
    G: int = 8
    B: int = 32
    T: int = 600
    eval_every: int = 10
    val_rollouts: int = 8


def make_items(n, world, rng):
    rel = rng.random((n, K)) < np.array(world.prevalence)
    d = rng.normal(np.array(world.d_mu), world.d_sd, (n, K))
    return rel, np.where(rel, d, -8.0)


class Policy:
    def __init__(self, world, n_train):
        self.w = world
        self.s = np.zeros(K)
        self.m = np.zeros(n_train)
        self.H = 0.0
        self.guard = False

    def probs(self, rel, d, m=None):
        w = self.w
        logit = d - self.s[None]
        if m is not None:
            logit = logit - m[:, None]
        p = sigmoid(logit)
        fl = np.array(w.floor)[None] * rel
        p = fl + (1 - fl) * p
        if w.hack:
            ph = sigmoid(-4 + w.kappa * self.H)
            p[:, EP] = 1 - (1 - p[:, EP]) * (1 - ph)
            p[:, BD] = p[:, BD] * sigmoid(-self.H + 4) / sigmoid(4)  # hack suppresses boundary misses
        return p

    def reward_weights(self):
        w = np.array(self.w.w, float)
        if self.w.hack and not self.guard:
            w[EP] = self.w.w_ep_hack
        return w


def rollouts(p, G, rng):
    return (rng.random((p.shape[0], G, K)) < p[:, None, :]).astype(float)  # B,G,K


def group_signal(e, w):
    r = 1 - (e * w).sum(-1)  # B,G
    sd = r.std(-1, keepdims=True)
    A = (r - r.mean(-1, keepdims=True)) / np.where(sd > 0, sd, 1.0)
    c = np.maximum(0.0, -(A[..., None] * e).mean(1))  # B,K
    return c, r


def fail_rate(p):
    return 1 - np.prod(1 - p, -1)


# ------------------------------------------------------------------ samplers

class Sampler:
    name = "base"

    def __init__(self, n, init_e, rng, world):
        self.n = n
        self.e_hat = init_e.copy()  # n,K last observed typed error rates
        self.rng = rng
        self.world = world
        self.val_hist = []  # list of K-vectors of val typed error

    def weights(self):
        return np.ones(self.n)

    def draw(self, B):
        w = np.maximum(self.weights(), 0) + 1e-12
        return self.rng.choice(self.n, B, replace=False, p=w / w.sum())

    def observe(self, idx, e):
        self.e_hat[idx] = e.mean(1)

    def observe_val(self, val_typed):
        self.val_hist.append(val_typed)


class Uniform(Sampler):
    name = "uniform"


class DiffFilter(Sampler):
    """Online difficulty filtering on stale pass rates: keep items with pass in [lo, hi]."""
    name = "difficulty-filter"
    lo, hi = 0.2, 0.8

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.passrate = 1 - fail_rate(self.e_hat)  # approx from typed rates (indep.)

    def observe(self, idx, e):
        super().observe(idx, e)
        self.passrate[idx] = 1 - e.max(-1).mean(1)

    def weights(self):
        keep = (self.passrate >= self.lo) & (self.passrate <= self.hi)
        if keep.sum() < 2 * self.world.B:  # fall back: nearest to 0.5
            order = np.argsort(np.abs(self.passrate - 0.5))[: 2 * self.world.B]
            keep = np.zeros(self.n, bool)
            keep[order] = True
        return keep.astype(float)


class Learnability(Sampler):
    """p(1-p) on the scalar pass rate (ODF / PODS-style)."""
    name = "learnability p(1-p)"

    def __init__(self, *a, **k):
        super().__init__(*a, **k)
        self.passrate = 1 - fail_rate(self.e_hat)

    def observe(self, idx, e):
        super().observe(idx, e)
        self.passrate[idx] = 1 - e.max(-1).mean(1)

    def weights(self):
        return self.passrate * (1 - self.passrate) + 0.01


class Typed(Sampler):
    """Item weight = sum_k u_k * e_hat_ik (1 - e_hat_ik) + eps: typed learnability with type weights u_k.

    mode='balance': u_k = 1 / (total learnable mass of type k)  (equal budget per still-learnable type)
    mode='stall'  : u_k multiplied up when the type's val error stopped falling (as specified in the brief)
    mode='lp'     : u_k proportional to recent learning progress of the type (Matiisen et al. TSCL)
    mode='flat'   : u_k = 1 (typed learnability without type reweighting)
    """
    W = 3  # val evals per window

    def __init__(self, *a, mode="balance", **k):
        super().__init__(*a, **k)
        self.mode = mode
        self.name = f"typed-{mode}"
        self.u = np.ones(K)

    def observe_val(self, val_typed):
        super().observe_val(val_typed)
        h = np.array(self.val_hist)
        if len(h) <= self.W:
            return
        prev, cur = h[-self.W - 1], h[-1]
        rel_drop = (prev - cur) / np.maximum(prev, 1e-3)
        if self.mode == "stall":
            stalled = (rel_drop < 0.05) & (cur > 0.02)
            self.u = np.clip(np.where(stalled, self.u * 1.5, self.u * 0.9 + 0.1), 0.25, 4.0)
        elif self.mode == "lp":
            lp = np.maximum(prev - cur, 0)
            self.u = 0.2 + lp / max(lp.max(), 1e-6)

    def weights(self):
        L = self.e_hat * (1 - self.e_hat)  # n,K
        u = self.u
        if self.mode == "balance":
            mass = L.sum(0)
            u = np.where(mass > 1e-6, 1.0 / np.maximum(mass, 1e-6), 0.0)
            u = u / max(u.max(), 1e-12)
        return (L * u).sum(1) + 1e-3 * (L.sum(1).max() + 1e-9)


SAMPLERS = {
    "uniform": lambda n, e, r, w: Uniform(n, e, r, w),
    "difficulty-filter": lambda n, e, r, w: DiffFilter(n, e, r, w),
    "learnability": lambda n, e, r, w: Learnability(n, e, r, w),
    "typed-flat": lambda n, e, r, w: Typed(n, e, r, w, mode="flat"),
    "typed-balance": lambda n, e, r, w: Typed(n, e, r, w, mode="balance"),
    "typed-stall": lambda n, e, r, w: Typed(n, e, r, w, mode="stall"),
    "typed-lp": lambda n, e, r, w: Typed(n, e, r, w, mode="lp"),
}


# ------------------------------------------------------------------ alarm

class Alarm:
    """Overfitting/hacking alarm: train reward rising while some val typed error rises above its best.

    For each type k: cur = mean of the last `win` val evals, best = lowest earlier `win`-eval rolling mean.
    Type k is 'rising' when cur - best > max(min_rise, z * se), se = sqrt(2 cur(1-cur) / (n_k * win)),
    n_k = number of val rollouts on items relevant to k. Fires when some type is rising AND the train
    reward EMA is higher than `win` evals ago, for `persist` consecutive evals."""

    def __init__(self, win=5, z=4.0, min_rise=0.03, persist=3):
        self.win, self.z, self.min_rise, self.persist = win, z, min_rise, persist
        self.train_r, self.val, self.roll = [], [], []
        self.streak = 0
        self.fired_at = None

    def update(self, step, train_reward_ema, val_typed, n_k=100):
        self.train_r.append(train_reward_ema)
        self.val.append(np.asarray(val_typed, float))
        if len(self.val) >= self.win:
            self.roll.append(np.mean(self.val[-self.win:], 0))
        if self.fired_at is not None or len(self.roll) <= self.win:
            return False
        cur = self.roll[-1]
        best = np.min(self.roll[:-self.win], 0)
        se = np.sqrt(2 * np.clip(cur, 0.01, 0.99) * (1 - np.clip(cur, 0.01, 0.99)) / (np.asarray(n_k) * self.win))
        rising = (cur - best) > np.maximum(self.min_rise, self.z * se)
        up = self.train_r[-1] > self.train_r[-1 - self.win]
        self.streak = self.streak + 1 if (rising.any() and up) else 0
        if self.streak >= self.persist:
            self.fired_at = step
            return True
        return False


# ------------------------------------------------------------------ run

def run(world, sampler_name, seed, alarm=False, track_best=True, alarm_kw=None):
    rng = np.random.default_rng(seed)
    wrng = np.random.default_rng(10_000 + seed)  # world identical across samplers for a seed
    tr_rel, tr_d = make_items(world.n_train, world, wrng)
    va_rel, va_d = make_items(world.n_val, world, wrng)
    te_rel, te_d = make_items(world.n_test, world, wrng)
    pol = Policy(world, world.n_train)
    # SFT-time estimation pass over the train set (same for every sampler)
    init_e = rollouts(pol.probs(tr_rel, tr_d, pol.m), world.G, wrng).mean(1)
    smp = SAMPLERS[sampler_name](world.n_train, init_e, rng, world)
    al = Alarm(**(alarm_kw or {})) if alarm else None
    alpha = np.array(world.alpha)
    log = dict(step=[], test_fail=[], test_typed=[], test_worst=[], train_fail=[], train_reward=[], H=[])
    r_ema = None
    best = (np.inf, None)
    visits = np.zeros(world.n_train)
    for t in range(world.T + 1):
        if t % world.eval_every == 0:
            pt = pol.probs(te_rel, te_d)
            typed = np.array([pt[te_rel[:, k], k].mean() if te_rel[:, k].any() else 0 for k in range(K)])
            log["step"].append(t)
            log["test_fail"].append(fail_rate(pt).mean())
            log["test_typed"].append(typed)
            log["test_worst"].append(typed.max())
            log["train_fail"].append(fail_rate(pol.probs(tr_rel, tr_d, pol.m)).mean())
            log["train_reward"].append(r_ema if r_ema is not None else np.nan)
            log["H"].append(pol.H)
            # noisy val evaluation (what a controller could actually see)
            ev = rollouts(pol.probs(va_rel, va_d), world.val_rollouts, rng).mean(1)
            val_typed = np.array([ev[va_rel[:, k], k].mean() if va_rel[:, k].any() else 0 for k in range(K)])
            smp.observe_val(val_typed)
            val_fail = fail_rate(ev).mean()
            if track_best and val_fail < best[0]:
                best = (val_fail, (log["test_fail"][-1], t))
            n_k = np.maximum(va_rel.sum(0), 1) * world.val_rollouts
            if al is not None and r_ema is not None and al.update(t, r_ema, val_typed, n_k):
                pol.guard = True
        if t == world.T:
            break
        idx = smp.draw(world.B)
        visits[idx] += 1
        p = pol.probs(tr_rel[idx], tr_d[idx], pol.m[idx])
        e = rollouts(p, world.G, rng)
        c, r = group_signal(e, pol.reward_weights())
        smp.observe(idx, e)
        gain = alpha * c.mean(0) * np.exp(-pol.s / world.s_cap)
        if world.transfer > 0:
            gain = gain + world.transfer * (gain.sum() - gain) / (K - 1)
        pol.s = (pol.s + gain) * (1 - world.forget)
        pol.m[idx] += world.alpha_m * c.sum(1)
        if world.hack:
            if pol.guard:
                pol.H *= 0.97
            else:
                pol.H += world.alpha_h * c[:, BD].mean() * 0.1
        rm = r.mean()
        r_ema = rm if r_ema is None else 0.9 * r_ema + 0.1 * rm
    out = {k: np.array(v) for k, v in log.items()}
    out["alarm_at"] = al.fired_at if al is not None else None
    out["best_val_ckpt"] = best[1]
    out["visits_gini"] = gini(visits)
    return out


def gini(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    if x.sum() == 0:
        return 0.0
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


SCENARIOS = {
    "S0 homogeneous types": dict(prevalence=(0.3,) * 5, d_mu=(0.5,) * 5, alpha=(0.4,) * 5),
    "S1 default (heterogeneous)": dict(),
    "S2 + unlearnable EP (floor 0.3)": dict(floor=(0, 0, 0.3, 0, 0)),
    "S3 + strong forgetting": dict(forget=0.004),
    "S4 + strong memorisation": dict(alpha_m=0.15),
    "S5 + reward hack": dict(hack=True),
    "S6 + cross-type transfer 0.5": dict(transfer=0.5),
    "S7 small train set (128)": dict(n_train=128, alpha_m=0.08),
    "S8 rare type also slow (MH α=0.5, q=0.05)": dict(prevalence=(0.12, 0.05, 0.08, 0.12, 0.7),
                                                      alpha=(0.4, 0.06, 0.4, 0.4, 0.4)),
    "S9 no diminishing returns (s_cap=∞)": dict(s_cap=1e9),
    "S10 ambiguous type (EP q=0.35, floor 0.45)": dict(prevalence=(0.12, 0.12, 0.35, 0.12, 0.7),
                                                       floor=(0, 0, 0.45, 0, 0)),
    "S11 slow reward hack (α_h=0.6)": dict(hack=True, alpha_h=0.6),
}


def world_for(name, **over):
    return replace(World(), **SCENARIOS[name], **over)
