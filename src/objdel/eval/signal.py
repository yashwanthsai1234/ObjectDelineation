"""GRPO-signal panel of K = 8 sampled groups per item, D12 checkpoint rule R-B and the per-item Spearman (PLAN §5.1
G-signal / G-probe, D12, M11; promoted from archive/phase2_lab/B2 score.py, analyze.py, posthoc.py `parts`).

Per item (one GRPO group of G = 8 samples at T = 1), from `objdel.eval.metrics` output (`group`):
- R = 0.3 IoU_tc + 0.3 F5 + 0.1 legal + 0.3 g_area, the G-probe raw composite on the legality-gated terms
  (illegal or unparseable samples: R = 0);
- r_std (population std of R), sigma (some gated geometric term above its REWARD_V4 floor; σ-share = its mean),
  mixed (>= 1 legal and >= 1 illegal), headroom (max R - mean R), distinct (distinct texts / K), raster_distinct
  (1 - duplicate share), legal_wrong (share legal with IoU_tc < 0.5), active (v4 advantages not all 0), abs_adv
  (mean |A_v4|), legal, r_mean, iou_mean;
- the exact split (law of total variance; illegal R = 0): Var(R) = p(1 - p) mu_L^2 + p Var_L(R), with p the legal
  share, mu_L / Var_L the legal samples' mean / population variance; sd_geom = sqrt(p Var_L) (D12's veto statistic,
  post-hoc in B2 and labelled so), geom_share = var_geom / var.

Statistics (PLAN §5.0): 95 % percentile bootstrap, B = 2,000, resampling LIDC nodules (item id up to the slice; any
other id is its own cluster); paired differences resample the same nodules for both arms (items aligned by id).
The draws reproduce B2's analyze.draws exactly (same seed, same order), so B2's CIs are reproduced to the digit.

R-B (D12): among eligible branches take the deepest; walk downward pairwise and move to the next shallower eligible
branch only if its sd_geom is higher with the paired CI excluding 0 AND by >= 15 % of the deeper branch's value AND
its median IoU_tc is not significantly lower (paired CI upper bound < 0 = significant loss); stop at the first pair
that does not fire (ties stay deep). The same walk under r_std (B2's pre-registered statistic) is reported, with
`disagree` when the two statistics would select differently or fire differently on an evaluated pair.
"""

from collections.abc import Callable, Sequence

import numpy as np
from scipy import stats

WEIGHTS = {"iou_tc": 0.3, "f5": 0.3, "legal": 0.1, "g_area": 0.3}  # REWARD_V4.weights on the raw gated terms
IOU_OK = 0.5  # legal_wrong: legal and IoU_tc < this (REWARD_V4.iou_gate)
BOOT, LEVEL, SEED = 2000, 0.95, 0
VETO_REL = 0.15
MISS = 0.1  # an item is a miss when its K-mean IoU_tc is below this (final_eval, ceiling)
DIGITS = 12  # decision values are rounded to 12 decimals before a >= / <= line (BK2-4, BK2-7): float sums of exact
#              shares land 1e-16 off (0.345 / 0.30 - 1 = 0.1499999...), real differences are far above 1e-12
VETO_STAT, ALT_STAT = "sd_geom", "r_std"


def group(item_id: str, texts: Sequence[str], m: dict) -> dict:
    """One item's group record from `objdel.eval.metrics(texts, ...)` output (B2 score.py's scored form)."""
    s = m["samples"]
    R = sum(w * np.asarray(s[t]) for t, w in WEIGHTS.items())
    return {
        "item_id": item_id,
        "R": [float(v) for v in R],
        "legal": s["legal"],
        "iou_tc": s["iou_tc"],
        "advantage": s["advantage"],
        "sigma": float(np.mean(m["groups"]["signal"])),
        "dup": float(np.mean(m["groups"]["dup"])),
        "distinct": len(set(texts)) / len(texts),
    }


def rows(columns: dict) -> list[dict]:
    """B2's column form ({key: [per item]}) to a list of group records."""
    return [dict(zip(columns, v, strict=True)) for v in zip(*columns.values(), strict=True)]


def split(R: np.ndarray, legal: np.ndarray) -> dict[str, np.ndarray]:
    """Var(R) = p(1-p) mu_L^2 + p Var_L(R) per row of R (n, K); illegal entries of R must be 0."""
    R, L = np.asarray(R, float), np.asarray(legal, float) > 0
    p, n = L.mean(1), np.maximum(L.sum(1), 1)
    mu = (R * L).sum(1) / n
    out = {
        "var": R.var(1),
        "var_legality": p * (1 - p) * mu**2,
        "var_geom": p * ((R - mu[:, None]) ** 2 * L).sum(1) / n,
    }
    out["sd_geom"] = np.sqrt(out["var_geom"])
    v = out["var"]
    out["geom_share"] = np.clip(np.where(v > 0, out["var_geom"] / np.where(v > 0, v, 1), 0.0), 0.0, 1.0)
    return out


def per_item(groups: Sequence[dict]) -> dict[str, np.ndarray]:
    """Per-item panel metrics (module docstring) of group records."""
    R, legal, iou = (np.array([g[k] for g in groups], float) for k in ("R", "legal", "iou_tc"))
    adv = np.abs(np.array([g["advantage"] for g in groups], float))
    L = legal > 0
    return {
        "r_std": R.std(1),
        "sigma": np.array([g["sigma"] for g in groups], float),
        "mixed": (L.any(1) & ~L.all(1)).astype(float),
        "headroom": np.maximum(R.max(1) - R.mean(1), 0.0),
        "distinct": np.array([g["distinct"] for g in groups], float),
        "raster_distinct": 1 - np.array([g["dup"] for g in groups], float),
        "legal_wrong": (L & (iou < IOU_OK)).mean(1),
        "active": (adv.max(1) > 1e-12).astype(float),
        "abs_adv": adv.mean(1),
        "legal": legal.mean(1),
        "r_mean": R.mean(1),
        "iou_mean": iou.mean(1),
    } | split(R, L)


def cluster(item_id: str) -> str:
    return item_id.rsplit("/", 1)[0] if item_id.startswith("lidc/") else item_id


def draws(ids: Sequence[str], boot: int = BOOT, seed: int = SEED, key: Callable[[str], str] = cluster
          ) -> list[np.ndarray]:
    """B cluster bootstrap index draws: each draw picks C of the C clusters (key(id); first-appearance order) with
    replacement and keeps all their items. The default key (nodule) reproduces B2's analyze.draws bit for bit;
    key = {item_id: group_id}.__getitem__ draws patients (A1's gate CI)."""
    if len(ids) == 0:
        raise ValueError("draws: no ids to resample")
    nod = [key(i) for i in ids]
    groups = [np.array([j for j, x in enumerate(nod) if x == u]) for u in dict.fromkeys(nod)]
    rng = np.random.default_rng(seed)
    return [np.concatenate([groups[g] for g in rng.integers(0, len(groups), len(groups))]) for _ in range(boot)]


def _ci(f, n: int, D: list[np.ndarray]) -> list[float]:
    reps, point = np.array([f(d) for d in D]), float(f(np.arange(n)))
    if np.isnan(point):  # a NaN item: no finite interval from the draws that missed it
        return [point, np.nan, np.nan]
    q = (1 - LEVEL) / 2
    return [point, float(np.nanquantile(reps, q)), float(np.nanquantile(reps, 1 - q))]


def ci(x: np.ndarray, D: list[np.ndarray]) -> list[float]:
    """[mean, lo, hi] of x under the draws D."""
    x = np.asarray(x, float)
    return _ci(lambda d: x[d].mean(), len(x), D)


def cluster_t(x, groups: Sequence[str], level: float = LEVEL) -> list[float]:
    """[mean, lo, hi]: the CR1 cluster-robust t interval of the item-weighted mean of x, with groups[i] the cluster
    (patient) of item i. se^2 = P/(P-1) sum_p (S_p - m N_p)^2 / N^2, t with P - 1 degrees of freedom (Cameron and
    Miller 2015, eq. 11-12 at one regressor). The decision interval of a gate on fewer than 50 patients (A1 step 5)."""
    x = np.asarray(x, float)
    _, inv = np.unique(np.asarray(groups), return_inverse=True)
    P = len(np.unique(inv))
    if P < 2:
        raise ValueError(f"cluster_t: needs at least 2 patients, got {P}")
    S, N = np.bincount(inv, x), np.bincount(inv).astype(float)
    m = S.sum() / N.sum()
    h = stats.t.ppf((1 + level) / 2, P - 1) * np.sqrt(P / (P - 1) * ((S - m * N) ** 2).sum() / N.sum() ** 2)
    return [float(m), float(m - h), float(m + h)]


def panel(groups: Sequence[dict], boot: int = BOOT, seed: int = SEED) -> dict:
    """Signal panel of one branch / view: per-item metrics, their means with nodule-cluster CIs, and the per-sample
    IoU_tc kept for paired median comparisons. JSON-native."""
    ids = [g["item_id"] for g in groups]
    m = per_item(groups)
    D = draws(ids, boot, seed)
    iou = [list(map(float, g["iou_tc"])) for g in groups]
    return {
        "items": len(ids),
        "nodules": len(set(map(cluster, ids))),
        "ids": ids,
        "per_item": {k: v.tolist() for k, v in m.items()},
        "iou_tc": iou,
        "summary": {k: ci(v, D) for k, v in m.items()} | {"iou_median": [float(np.median(iou))]},
    }


def paired(a: dict, b: dict, key: str, boot: int = BOOT, seed: int = SEED) -> list[float]:
    """[a - b, lo, hi] on the items both panels have; key = a per-item metric (difference of means) or "iou_median"
    (difference of the median IoU_tc over the resampled items' samples)."""
    ib = {i: j for j, i in enumerate(b["ids"])}
    ia = [j for j, i in enumerate(a["ids"]) if i in ib]
    ids = [a["ids"][j] for j in ia]
    jb = [ib[i] for i in ids]
    D = draws(ids, boot, seed)
    if key == "iou_median":
        A, B = np.asarray(a["iou_tc"], float)[ia], np.asarray(b["iou_tc"], float)[jb]
        return _ci(lambda d: np.median(A[d]) - np.median(B[d]), len(ids), D)
    x = np.asarray(a["per_item"][key], float)[ia] - np.asarray(b["per_item"][key], float)[jb]
    return ci(x, D)


def _pair(shallow: dict, deep: dict) -> dict:
    """Veto evidence for moving from `deep` to the next shallower `shallow`, under both statistics."""
    s, d = shallow["panel"], deep["panel"]
    iou = paired(s, d, "iou_median")
    out = {"deep": deep["name"], "shallow": shallow["name"], "iou_median": iou}
    ix = {i: j for j, i in enumerate(d["ids"])}
    for stat in (VETO_STAT, ALT_STAT):
        diff = paired(s, d, stat)
        base = float(np.mean([d["per_item"][stat][ix[i]] for i in s["ids"] if i in ix]))  # deeper branch's value
        rel = round(diff[0] / base, DIGITS) if base > 0 else (np.inf if diff[0] > 0 else 0.0)
        out[stat] = {"ci": diff, "rel": float(rel), "fires": bool(diff[1] > 0 and rel >= VETO_REL and iou[2] >= 0)}
    return out


def select_rb(branches: Sequence[dict]) -> tuple[str | None, dict]:
    """D12 rule R-B over branches {"name", "depth", "eligible", "panel"} (panel: `panel` output on the shifted val).
    Returns (selected name, or None when no branch is eligible; reasons). ValueError on a repeated name or depth
    (the walk order would be the input order, BK2-8)."""
    for key in ("name", "depth"):
        values = [b[key] for b in branches]
        if len(set(values)) != len(values):
            raise ValueError(f"R-B: repeated branch {key} in {values}")
    order = sorted((b for b in branches if b["eligible"]), key=lambda b: -b["depth"])
    why = {"eligible": [b["name"] for b in order], "pairs": [], "selected_r_std": None, "disagree": False}
    if not order:
        return None, why
    memo = {}

    def walk(stat):
        cur, seen = 0, []
        for k in range(1, len(order)):
            if (cur, k) not in memo:
                memo[cur, k] = _pair(order[k], order[cur])
            p = memo[cur, k]
            seen.append(p)
            if not p[stat]["fires"]:
                break
            cur = k
        return order[cur]["name"], seen

    name, pairs = walk(VETO_STAT)
    alt, _ = walk(ALT_STAT)
    why |= {"pairs": pairs, "selected_r_std": alt}
    why["disagree"] = alt != name or any(p[VETO_STAT]["fires"] != p[ALT_STAT]["fires"] for p in pairs)
    return name, why


def spearman_ci(x, y, ids: Sequence[str] | None = None, boot: int = BOOT, seed: int = SEED) -> list[float]:
    """[rho, lo, hi]: Spearman's rank correlation (scipy) with a nodule-cluster (ids given) or item bootstrap CI."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    D = draws([str(i) for i in range(len(x))] if ids is None else ids, boot, seed)
    return _ci(lambda d: stats.spearmanr(x[d], y[d]).statistic, len(x), D)
