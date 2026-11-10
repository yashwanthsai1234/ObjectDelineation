"""Gate power check before any spend (Implementation.md B1.9; $0, Mac).

Input: R5's per-item K-means of IoU_tc on val (`runs/r5/eval/eval_val/items.jsonl`), each item's offset from the
frozen view (`view.offset_frac`), its size (the doctors' union-box extent, from the store) and patient.

Simulated B1 effect: a share of the far misses (offset > 0.3, K-mean IoU_tc < 0.1), taken patient by patient in a
random patient order, become found images with an IoU drawn from R5's found images of the same B1.1 cell (size x
offset). Noise: each arm gets independent N(0, c mu_i (1 - mu_i) / 8) per item, mu_i the arm's true K-mean, not
clipped. c = 0.077 is the one measured pair (two K = 8 samplings of rl-full-3ep, retest variance 0.00324); c = 1 is
the Bernoulli bound on both arms (the stress line); noise 0 is the third line. With --item-sd (A3's split-half
`noise_floor.item_sd`, the SD of a retest difference) one arm gets item_sd^2 / 2 (math-stats §11 item 1).

Gate B1 rules (B.6.2, paired, patient percentile interval with signal.draws' B and seed; every subset has >= 50
patients, so no rule needs cluster_t): (1) val IoU CI lower end > 0; (2) far-miss share <= 20.5 % and the CI upper
end of the paired miss difference < 0; (3) IoU at offset < 0.15 CI lower end > -0.02; (4) the CI lower end of
(far gain - near gain) > 0 on the same draws. Joint P(pass) = the share of runs that pass all four, nothing else.
Not simulated: the edge-view rule (no R5 edge values before V1), the legal share and the frozen-view hash.

Rule (2) margin (math-stats §11 item 5): at the 25 % effect the B1 far-miss point is (111 - 28) / 406 = 20.44 %
against 20.5 %, a quarter of an item: a run fails rule (2) as soon as noise moves one net far item below 0.1, so
the joint value is capped by P(rule 2), and a 20 % effect fails it by construction. The smallest effect with joint
P(pass) >= 0.8 (searched in 5 % steps) is the useful number.

Gate B3 power (report only, never in the joint value): on val-select, the B3 arm converts a further share of the
B1 arm's remaining misses (same rule); both arms get the B1 arm's c; rules with signal.cluster_t: the miss share
falls >= 25 % relative with the CI upper end of (B3 - B1) < 0; IoU on items found by both, CI lower end > -0.01.

.venv/bin/python -I -m experiments.phase2.gate_power [--item-sd X] --out runs/b1/gate_power
"""

import argparse
import json
from pathlib import Path
from typing import NamedTuple

import numpy as np

from objdel.data.windows import OFFSET_EDGES, SIZE_EDGES  # B1.1 cells
from objdel.io import read_json, read_jsonl
from objdel.verify import stats as signal

K = 8
MISS = signal.MISS
FAR, NEAR = 0.3, 0.15
FAR_MISS_MAX, NEAR_FLOOR = 0.205, -0.02
B3_FALL, B3_BOTH_FLOOR = 0.75, -0.01
C_MEASURED = 0.077
EFFECT = 0.25
B3_EFFECTS = (0.25, 0.35, 0.5)
GRID = tuple(round(0.05 * i, 2) for i in range(1, 21))  # effects searched for joint P(pass) >= 0.8
POWER = 0.8
LINES = {"c=0.077": (C_MEASURED, C_MEASURED), "c=1": (1.0, 1.0), "noise 0": (0.0, 0.0)}  # (baseline, B1) arm
B1_RULES = ("iou", "far_miss", "near_iou", "did")
B3_RULES = ("miss_share", "both_found_iou")
NOT_SIMULATED = ("edge view", "legal share", "frozen-view hash")
ARMS = ("base", "b1", "b1_re", "b3")
MIN_PATIENTS = 50  # below this a gate set uses cluster_t (A1 step 5)


class Items(NamedTuple):
    ids: list[str]
    mu: np.ndarray  # R5 K-mean IoU_tc
    groups: np.ndarray  # patient code 0..P-1
    cell: np.ndarray  # 4 x size bin + offset bin
    far: np.ndarray
    near: np.ndarray
    select: np.ndarray  # val-select half (Gate B3)


def cell_of(offset, extent) -> np.ndarray:
    return 4 * np.searchsorted(SIZE_EDGES, extent, side="right") + np.searchsorted(OFFSET_EDGES, offset, side="right")


def weights(ids: list[str], groups: np.ndarray) -> np.ndarray:
    """(B, N) item multiplicities of the gate's patient draws (signal.draws, B and seed as the gate)."""
    D = signal.draws(ids, key=dict(zip(ids, groups.tolist(), strict=True)).__getitem__)
    return np.stack([np.bincount(d, minlength=len(ids)) for d in D]).astype(float)


def streams(seed: int, runs: int, n: int):
    """The effect generator and one independent standard normal (runs, n) per arm."""
    ss = np.random.SeedSequence(seed).spawn(1 + len(ARMS))
    return np.random.default_rng(ss[0]), {a: np.random.default_rng(s).standard_normal((runs, n))
                                          for a, s in zip(ARMS, ss[1:], strict=True)}


def noisy(true: np.ndarray, c: float, z: np.ndarray) -> np.ndarray:
    return true + np.sqrt(c * true * (1 - true) / K) * z


def convert(true: np.ndarray, eligible: np.ndarray, items: Items, share: float, rng) -> np.ndarray:
    """Copy of true (runs, N): per run, round(share x eligible count) eligible items, taken patient by patient in a
    random patient order, get an IoU drawn from R5's found items of their cell."""
    out = true.copy()
    pools = {c: items.mu[(items.cell == c) & (items.mu >= MISS)] for c in np.unique(items.cell)}
    for r in range(len(out)):
        idx = np.flatnonzero(eligible[r])
        n = round(share * len(idx))
        if not n:
            continue
        rank = rng.permutation(items.groups.max() + 1)
        pick = idx[np.lexsort((rng.random(len(idx)), rank[items.groups[idx]]))][:n]
        for c in np.unique(items.cell[pick]):
            at = pick[items.cell[pick] == c]
            if not len(pools[c]):
                raise ValueError(f"cell {c} has no found item to draw from")
            out[r, at] = rng.choice(pools[c], len(at))
    return out


def _arms(items: Items, share: float, c_base: float, c_b1: float, rng, z):
    mu = np.broadcast_to(items.mu, z["b1"].shape)
    true = convert(mu, (mu < MISS) & items.far, items, share, rng)
    return true, noisy(true, c_b1, z["b1"]), noisy(mu, c_base, z["base"])


def arms(items: Items, share: float, c_base: float, c_b1: float, runs: int, seed: int):
    """(B1 true K-means, B1 arm, baseline arm), each (runs, N)."""
    return _arms(items, share, c_base, c_b1, *streams(seed, runs, len(items.mu)))


def _reps(x: np.ndarray, f: np.ndarray, W: np.ndarray) -> np.ndarray:
    """(runs, B) item mean of x[:, f] under each draw; NaN for a draw with no item of f (signal._ci drops it)."""
    Wf = W[:, f]
    with np.errstate(invalid="ignore", divide="ignore"):
        return x[:, f] @ Wf.T / Wf.sum(1)


def _ends(reps: np.ndarray) -> np.ndarray:
    q = (1 - signal.LEVEL) / 2
    return np.nanquantile(reps, [q, 1 - q], axis=1)


def b1_rules(b1: np.ndarray, base: np.ndarray, far: np.ndarray, near: np.ndarray, W: np.ndarray) -> dict:
    """Gate B1 rules (1)-(4) per run (module docstring) on the paired arms (runs, N)."""
    d = b1 - base
    miss = (b1 < MISS).astype(float) - (base < MISS)
    lo = _ends(_reps(d, np.ones(d.shape[1], bool), W))[0]
    far_point = np.round((b1[:, far] < MISS).mean(1), signal.DIGITS)
    far_hi = _ends(_reps(miss, far, W))[1]
    near_lo = _ends(_reps(d, near, W))[0]
    did = _reps(d, far, W) - _reps(d, near, W)
    did_lo = _ends(did)[0]
    nan = lambda r: np.isnan(r).sum(1)  # noqa: E731  draws with no subset item, dropped by _ends (math-stats s1)
    return {"iou": lo > 0, "far_miss": (far_point <= FAR_MISS_MAX) & (far_hi < 0), "near_iou": near_lo > NEAR_FLOOR,
            "did": did_lo > 0, "far_nan_draws": nan(_reps(miss, far, W)), "near_nan_draws": nan(_reps(d, near, W)),
            "did_nan_draws": nan(did)}


def b3_rules(b3: np.ndarray, b1: np.ndarray, groups: np.ndarray) -> dict:
    """The two Gate B3 rules on item values per run, with signal.cluster_t (val-select has < 50 patients)."""
    out = {r: np.zeros(len(b3), bool) for r in B3_RULES}
    for r, (x3, x1) in enumerate(zip(b3, b1, strict=True)):
        m3, m1 = x3 < MISS, x1 < MISS
        point = round(m3.mean(), signal.DIGITS) <= round(B3_FALL * m1.mean(), signal.DIGITS)
        out["miss_share"][r] = point and signal.cluster_t(m3.astype(float) - m1, groups)[2] < 0
        both = ~m3 & ~m1  # under 2 patients there: no interval, the rule fails
        out["both_found_iou"][r] = (len(np.unique(groups[both])) >= 2
                                    and signal.cluster_t((x3 - x1)[both], groups[both])[1] > B3_BOTH_FLOOR)
    return out


def simulate(items: Items, W: np.ndarray, share: float, c_base: float, c_b1: float, runs: int, seed: int,
             b3_shares=()) -> dict:
    """Per-run Gate B1 rule flags at one effect, and the Gate B3 flags per B3 share (both B3 arms at c_b1)."""
    rng, z = streams(seed, runs, len(items.mu))
    true, b1, base = _arms(items, share, c_base, c_b1, rng, z)
    s, b3 = items.select, {}
    for share3 in b3_shares:
        t3 = convert(true, (true < MISS) & s, items, share3, rng)
        b3[share3] = b3_rules(noisy(t3, c_b1, z["b3"])[:, s], noisy(true, c_b1, z["b1_re"])[:, s], items.groups[s])
    return {"rules": b1_rules(b1, base, items.far, items.near, W), "b3": b3}


def summary(flags: dict) -> dict:
    """P(pass) per rule and jointly; the max empty-subset draw count per B1 subset rule (math-stats s1)."""
    nan = {k: int(v.max()) for k, v in flags.items() if k.endswith("_nan_draws")}
    flags = {k: v for k, v in flags.items() if k not in nan}
    return ({"nan_draws_max": nan} if nan else {}) | {
        "p_pass": {k: float(v.mean()) for k, v in flags.items()}, "joint": float(np.all(list(flags.values()), axis=0).mean())}


def c_from_item_sd(item_sd: float, mu: np.ndarray) -> float:
    """c whose mean per-arm variance equals item_sd^2 / 2 (item_sd: SD of a retest difference of two K-means)."""
    return item_sd**2 / 2 / float(np.mean(mu * (1 - mu) / K))


def report(items: Items, W: np.ndarray, runs: int, seed: int, lines: dict, b3_shares=B3_EFFECTS) -> dict:
    fm = int((items.far & (items.mu < MISS)).sum())
    doc = {"runs": runs, "seed": seed, "effect": EFFECT, "items": len(items.ids),
           "patients": {k: int(len(np.unique(items.groups[m]))) for k, m in
                        (("all", np.ones(len(items.ids), bool)), ("far", items.far), ("near", items.near),
                         ("select", items.select))},
           "far_items": int(items.far.sum()), "far_misses": fm,
           "far_miss_point_at_effect": (fm - round(EFFECT * fm)) / int(items.far.sum()),
           "not_simulated": list(NOT_SIMULATED), "lines": {}}
    for name, (cb, c1) in lines.items():
        sim = simulate(items, W, EFFECT, cb, c1, runs, seed, b3_shares)
        joint = {str(s): summary(simulate(items, W, s, cb, c1, runs, seed)["rules"])["joint"] for s in GRID}
        doc["lines"][name] = summary(sim["rules"]) | {
            "c": [cb, c1], "b3": {str(s): summary(f) for s, f in sim["b3"].items()},
            "joint_by_effect": joint, "smallest_effect": next((float(s) for s in joint if joint[s] >= POWER), None),
            "flags": {k: v.astype(int).tolist() for k, v in sim["rules"].items()}}
    return doc


def load(items_path: Path, view: Path, store_root: Path, halves_path: Path) -> Items:
    """R5 items joined with the frozen val view (patient, offset_frac) and the store (union-box extent)."""
    from objdel.data import store
    from objdel.data.targets.frame import union_box

    rows = read_jsonl(items_path)
    with open(view / "lidc" / "val.jsonl") as f:
        recs = {r["item_id"]: r for r in map(json.loads, f)}
    extent = {}
    for it in store.load("lidc", "val", store_root):
        if it.item_id in recs:
            x0, y0, x1, y1 = union_box(it)
            extent[it.item_id] = max(x1 - x0, y1 - y0)
    ids = [r["item_id"] for r in rows]
    names = np.array([recs[i]["group_id"] for i in ids])
    off = np.array([recs[i]["view"]["offset_frac"] for i in ids])
    halves = read_json(halves_path)["halves"]
    items = Items(ids, np.array([r["iou_tc"] for r in rows], float), np.unique(names, return_inverse=True)[1],
                  cell_of(off, np.array([extent[i] for i in ids])), off > FAR, off < NEAR,
                  np.array([halves[g] == "select" for g in names]))
    for name, m in (("all", np.ones(len(ids), bool)), ("far", items.far), ("near", items.near)):
        if len(np.unique(names[m])) < MIN_PATIENTS:
            raise ValueError(f"the {name} subset has < {MIN_PATIENTS} patients: the gate would use cluster_t")
    return items


def table(doc: dict) -> str:
    out = [f"Gate B1 power: {doc['runs']} runs, effect {doc['effect']:.0%} of {doc['far_misses']} far misses; "
           f"far-miss point at the effect {doc['far_miss_point_at_effect']:.4%} against {FAR_MISS_MAX:.1%} "
           "(a quarter of an item under rule 2's threshold)"]
    for name, line in doc["lines"].items():
        p = line["p_pass"]
        out.append(f"{name:>10}: " + "  ".join(f"P({r}) {p[r]:.3f}" for r in B1_RULES)
                   + f"  joint P(pass) {line['joint']:.3f}  smallest effect with joint >= {POWER}: "
                   f"{line['smallest_effect']}  max empty-subset draws {line['nan_draws_max']}")
        out += [f"{'':>10}  Gate B3 report, further {s} of B1 misses: "
                + "  ".join(f"P({r}) {v:.3f}" for r, v in b3["p_pass"].items()) + f"  both {b3['joint']:.3f}"
                for s, b3 in line["b3"].items()]
    out.append("not simulated: " + ", ".join(doc["not_simulated"]) + "; rules are read from the c=0.077 line")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--items", type=Path, default=Path("runs/r5/eval/eval_val/items.jsonl"))
    ap.add_argument("--view", type=Path, default=Path("data/eval_views/lidc-val-shifted"))
    ap.add_argument("--store", type=Path, default=Path("data/store"))
    ap.add_argument("--halves", type=Path, default=Path("docs/reference/val_halves.json"))
    ap.add_argument("--runs", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--item-sd", type=float, help="A3's noise_floor.item_sd: adds a line, for the record only")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    items = load(a.items, a.view, a.store, a.halves)
    lines = dict(LINES)
    if a.item_sd is not None:
        c = c_from_item_sd(a.item_sd, items.mu)
        lines["A3 item_sd"] = (c, c)
    doc = report(items, weights(items.ids, items.groups), a.runs, a.seed, lines)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "gate_power.json").write_text(json.dumps(doc))
    print(table(doc))


if __name__ == "__main__":
    main()
