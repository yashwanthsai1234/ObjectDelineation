"""Task 4: replay group-advantage schemes on GRPO-like groups built from the fixtures.

Each group (one fixture item) has 8 outputs, built around the consensus (strict-majority) mask V
so that the known orderings are unambiguous:
  g1, g2  good: V, and V with sub-pixel boundary jitter
  o1, o2  slightly off: V shifted by (2, 1) px, V grown by 2 px (topology kept)
  b1, b2  topology broken, high IoU: g1 + speck, g2 + pinhole, both of min_area px (size filter)
  inf     inflated: rater union dilated by 2 tau (area-inflation hack)
  far     V shifted by max(5, 0.8 r_eq) px (low IoU)
Known-better pairs: g1>o1, g2>o2, g1>b1, g2>b2, g1>inf, g2>inf, g1>far, o1>far.
Terms: iou (vs vote), band (F5), topo (per-instance topology matches a rater), legal (1 unless a
scenario changes it); weights 0.3 / 0.3 / 0.3 / 0.1.

Dark-room scenarios: a term that should carry no information is made near-constant and the
largest change of any advantage it causes is measured (vs the same term exactly constant):
  D1 legality numeric noise: legal = 1 - U(0, 1e-3) for all outputs
  D2 one slightly illegal output: legal = 0.99 for one random output, 1 elsewhere
  D3 topology near-constant: topo = 1 - U(0, 1e-3) for all outputs
  D4 converged group: 8 copies of g1 with sub-pixel jitter (all terms nearly equal)

    ../.venv/bin/python scripts/advantage_replay.py
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from rewardcore import fixtures, perturb  # noqa: E402
from rewardcore.advantages import SCHEMES as V1_SCHEMES  # noqa: E402
from rewardcore.advantages import dm_loo_topology, reward_v2, reward_v3  # noqa: E402
from rewardcore.reward import ExpertCache, score, signature_id  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
NAMES = ("iou", "band", "topo", "legal")
W = np.array([0.3, 0.3, 0.3, 0.1])
FLOORS = np.array([0.05, 0.05, 0.25, 0.25])
LOW_VAR = np.array([False, False, True, True])
OUTPUTS = ("g1", "g2", "o1", "o2", "b1", "b2", "inf", "far")
PAIRS = [("g1", "o1"), ("g2", "o2"), ("g1", "b1"), ("g2", "b2"), ("g1", "inf"), ("g2", "inf"), ("g1", "far"), ("o1", "far")]
B = 256


def terms(m, cache):
    t = score(m, cache)
    # columns 0-3: v1 terms; 4-7: IoU_topocond, rater frequency of the topology, topology id, topo noise (v2)
    return [t.iou, t.band["F5"], t.topo_any, 1.0, t.iou_tc, t.topo_freq, float(signature_id(t.sig)), 0.0]


def build_group(args):
    key, stack, seed = args
    rng = np.random.default_rng(seed)
    if len(stack) < 2:
        return None
    cache = ExpertCache.build(stack)
    g1 = cache.vote
    g2 = perturb.jitter(g1, 0.3, rng, corr=4.0)
    r_eq = float(np.sqrt(g1.sum() / np.pi))
    A = cache.min_area  # defects at the signature size filter (iteration 3, R6)
    b1 = perturb.add_speck(g1, A, float(rng.uniform(3, 8)) + np.sqrt(A / np.pi), rng)
    b2 = perturb.punch_hole(g2, A, rng)
    if b1 is None or b2 is None:
        return None
    d = max(5, int(round(0.8 * r_eq)))
    masks = dict(g1=g1, g2=g2, o1=perturb.shift(g1, 2, 1), o2=perturb.grow(g1, 2.0), b1=b1, b2=b2,
                 inf=perturb.dilate(cache.union, 2 * cache.tau), far=perturb.shift(g1, d, d // 2))
    R = np.array([terms(masks[o], cache) for o in OUTPUTS])
    # good and slightly-off must keep a rater topology, broken must not
    if not (R[[0, 1, 2, 3], 2].all() and not R[[4, 5], 2].any()):
        return None
    conv = np.array([terms(perturb.jitter(g1, 0.3, rng, corr=4.0), cache) for _ in range(8)])
    return dict(key=key, R=R.tolist(), conv=conv.tolist())


W2, FLOORS2 = np.array([0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25])


def _v1(fn):
    return lambda X: fn(X[..., :4], W, floors=FLOORS, low_var=LOW_VAR, names=NAMES)


def _v2(c_topo, fn=reward_v2, debias=False):
    def run(X):
        valid = (X[..., 2] > 0) & (X[..., 4] >= 0.5)
        r_topo = dm_loo_topology(X[..., 6], X[..., 5], valid, lam=0.5, debias=debias) + X[..., 7]
        return fn(X[..., [4, 1, 3]], W2, FLOORS2, ("iou", "band", "legal"), r_topo, c_topo=c_topo)
    return run


SCHEMES = {name: _v1(fn) for name, fn in V1_SCHEMES.items()}
SCHEMES["v2 (c_topo=2)"] = _v2(2.0)
SCHEMES["v3 frozen (c_topo=2)"] = _v2(2.0, reward_v3, debias=True)


def run_schemes(R):
    out = {}
    for name, fn in SCHEMES.items():
        out[name] = fn(R)[0]
    return out


def ordering(adv):
    idx = {o: j for j, o in enumerate(OUTPUTS)}
    ok = np.array([[a[idx[b]] > a[idx[w]] for b, w in PAIRS] for a in adv])
    return ok.mean(), ok.mean(0)


def main():
    items = [(k, st) for ds, k, st in fixtures.load_all()]
    rng = np.random.default_rng(0)
    order = rng.permutation(len(items))
    args = [(items[j][0], items[j][1], int(j)) for j in order]
    with Pool(8) as pool:
        groups = []
        for g in pool.imap(build_group, args, chunksize=4):
            if g is not None:
                groups.append(g)
            if len(groups) >= B:
                break
    R = np.array([g["R"] for g in groups])  # (B, 8, 8)
    conv = np.array([g["conv"] for g in groups])
    (OUT / "advantage_groups.json").write_text(json.dumps(groups))
    print(f"groups: {len(R)} x {R.shape[1]} outputs; mean terms per output type:")
    for j, o in enumerate(OUTPUTS):
        print(f"  {o:4s} " + " ".join(f"{n}={R[:, j, i].mean():.3f}" for i, n in enumerate(NAMES))
              + f" iou_tc={R[:, j, 4].mean():.3f}")

    base = run_schemes(R)
    print("\n(i) known-better > known-worse ordering (share of pairs); per pair type: " + " ".join(f"{b}>{w}" for b, w in PAIRS))
    for name, a in base.items():
        tot, per = ordering(a)
        print(f"  {name:11s} {tot:6.1%}   " + " ".join(f"{p:6.1%}" for p in per))

    noise = np.random.default_rng(1)
    scen = {}
    R1 = R.copy()
    R1[..., 3] = 1 - noise.uniform(0, 1e-3, R1.shape[:2])
    scen["D1 legal noise 1e-3"] = (R1, 3)
    R2 = R.copy()
    R2[np.arange(len(R)), noise.integers(0, 8, len(R)), 3] = 0.99
    scen["D2 one output legal=0.99"] = (R2, 3)
    R3 = R.copy()
    u = noise.uniform(0, 1e-3, R3.shape[:2])
    R3[..., 2] = 1 - u  # v1: the topology term itself
    R3[..., 7] = -u  # v2: the same noise on the DM-loo topology reward
    scen["D3 topo near-constant"] = (R3, 2)
    print("\n(ii) dark room: max |change of an advantage| caused by a near-constant term"
          " (median over groups / max over batch), and ordering accuracy under it")
    print(f"  {'scheme':11s} " + " ".join(f"{s:>34s}" for s in scen))
    for name, fn in SCHEMES.items():
        cells = []
        for Rn, k in scen.values():
            Rc = Rn.copy()
            Rc[..., k] = 1.0  # the same term exactly constant
            Rc[..., 7] = 0.0
            an, ac = fn(Rn)[0], fn(Rc)[0]
            d = np.abs(an - ac).max(1)
            acc, _ = ordering(an)
            cells.append(f"{np.median(d):6.3f} / {d.max():6.3f} (ord {acc:5.1%})")
        print(f"  {name:11s} " + " ".join(f"{c:>34s}" for c in cells))
    half = len(R) // 2
    mixed = np.concatenate([R[:half], conv[half:]])  # converged groups inside a normal batch
    s = conv[half:, :, :4] @ W
    print("\n  D4 converged groups (8 near-identical outputs) mixed 50/50 into a normal batch:")
    print(f"    raw total reward of converged groups: median within-group std {np.median(s.std(1)):.2e}, "
          f"max range {np.max(np.ptp(s, 1)):.2e}")
    for name, fn in SCHEMES.items():
        a = fn(mixed)[0]
        ac, an = np.abs(a[half:]).max(1), np.abs(a[:half]).max(1)
        print(f"    {name:11s} converged groups: median max|A| {np.median(ac):.3f}, max {ac.max():.3f}"
              f" | normal groups: median max|A| {np.median(an):.3f}")


if __name__ == "__main__":
    main()
