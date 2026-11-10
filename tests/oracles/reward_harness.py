"""Gate harness for reward v4 (tests and experiments/phase1/reward_gates.py): fixture loaders, controlled mask
perturbations, realistic honest GRPO rollouts, the adversary suite, K9 pairs and cluster-bootstrap CIs.

Ported from the lab (archive/phase1_lab/t_reward: rewardcore/perturb.py, scripts/honest.py,
scripts/adversaries.py, scripts/iter4_adversaries.py, scripts/k9_discrimination.py, scripts/k9_v3.py) onto the
production objdel.reward. Self-contained (numpy, scipy, OpenCV, objdel) so it can be imported both as
`oracles.reward_harness` (pytest) and `tests.oracles.reward_harness` (experiments: one harness for the pytest gates
and the full-data gate tables, so they cannot drift; ISSUES I19).
"""
import csv
from pathlib import Path

import cv2
import numpy as np
from scipy import ndimage

from objdel.geom.topology import instance_topology, min_area, signature
from objdel.reward import ExpertCache, advantages_v4, score_mask, sig_id, voters
from objdel.reward.cache import edt  # exact EDT; 1e6 where no zero pixel exists

G = 8
# packed per-rollout terms (one row per rollout)
FIELDS = ("iou_tc", "f5", "g_area", "topo_any", "topo_freq", "sig_id", "area")
IOU, F5, GAREA, TANY, TFREQ, SID, AREA = range(len(FIELDS))
HOLES = ("holes_stamped", "minority_hole")  # adversaries judged against honest rollouts of their own topology


# ---------------------------------------------------------------- data

def load_fixture_items(fixtures_dir: Path) -> list[tuple[str, str, np.ndarray]]:
    """[(dataset, key, stack)] for RIGA ('<image>__disc' / '__cup') and LIDC fixture items with >= 2 raters."""
    out = []
    for ds in ("riga", "lidc"):
        with np.load(Path(fixtures_dir) / f"{ds}_experts.npz") as z:
            out += [(ds, k, z[k].astype(bool)) for k in sorted(z.files) if len(z[k]) >= 2]
    return out


def subset(ds: str, key: str) -> str:
    return f"{ds}/{key.split('__')[1]}"


def load_lidc_disagreement(fixtures_dir: Path) -> dict[str, np.ndarray]:
    """Full-LIDC train items (no duplicated reader, k >= 2) whose voting raters' size-filtered signatures differ
    (lab ITERATION3 R6 / M3 census: 251 train items; 239 once raters whose whole mask is below the size filter
    abstain, objdel.reward.voters: on the other 12 the abstainer was the only different answer)."""
    fixtures_dir = Path(fixtures_dir)
    with open(fixtures_dir / "lidc_full_meta.csv") as f:
        meta = {r["key"]: r for r in csv.DictReader(f)}
    out = {}
    with np.load(fixtures_dir / "lidc_full.npz") as z:
        for k in sorted(z.files):
            r = meta.get(k)
            if r is None or r["dup_reader_in_cluster"] == "True" or r["split"] != "train":
                continue
            st = z[k].astype(bool)
            a = min_area(int(st.any(0).sum()))
            sigs = [signature(m, a) for m in st]
            if len(st) >= 2 and len({sigs[i] for i in voters(sigs)}) > 1:
                out[k] = st
    return out


# ---------------------------------------------------------------- perturbations (lab rewardcore/perturb.py)

def signed_distance(m: np.ndarray) -> np.ndarray:
    m = np.asarray(m, bool)
    return np.where(m, edt(m), -edt(~m)).astype(np.float32)


def grow(m: np.ndarray, d: float) -> np.ndarray:
    """Euclidean dilation (d > 0) or erosion (d < 0) by |d| px."""
    m = np.asarray(m, bool)
    if d > 0:
        return m | (edt(~m) <= d)
    return edt(m) > -d if d < 0 else m.copy()


def shift(m: np.ndarray, dy: int, dx: int) -> np.ndarray:
    out = np.zeros_like(m)
    h, w = m.shape
    ys, yd = (slice(0, h - dy), slice(dy, h)) if dy >= 0 else (slice(-dy, h), slice(0, h + dy))
    xs, xd = (slice(0, w - dx), slice(dx, w)) if dx >= 0 else (slice(-dx, w), slice(0, w + dx))
    out[yd, xd] = m[ys, xs]
    return out


def disk_offsets(area: int) -> np.ndarray:
    """Offsets of the most compact digital disk with exactly `area` pixels."""
    r = int(np.ceil(np.sqrt(area / np.pi))) + 2
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    order = np.lexsort((xx.ravel(), yy.ravel(), np.hypot(yy, xx).ravel()))
    return np.stack([yy.ravel()[order][:area], xx.ravel()[order][:area]], 1)


def add_speck(m, area: int, gap: float, rng) -> np.ndarray | None:
    """A detached compact speck of `area` px about `gap` px outside the mask (None if no room)."""
    off = disk_offsets(area)
    reach = np.abs(off).max()
    dist = edt(~m)
    h, w = m.shape
    cand = np.argwhere((np.abs(dist - gap) <= 0.5) & (dist > reach + 1.5))
    cand = cand[(cand[:, 0] >= reach) & (cand[:, 0] < h - reach) & (cand[:, 1] >= reach) & (cand[:, 1] < w - reach)]
    if not len(cand):
        return None
    c = cand[rng.integers(len(cand))]
    out = m.copy()
    out[c[0] + off[:, 0], c[1] + off[:, 1]] = True
    return out


def punch_hole(m, area: int, rng) -> np.ndarray | None:
    """Remove a compact `area`-px region strictly inside the mask (None if no room)."""
    off = disk_offsets(area)
    reach = np.abs(off).max()
    depth = edt(m)
    if depth.max() <= reach + 1.5:
        return None
    cand = np.argwhere(depth > reach + 1.5)
    c = cand[rng.integers(len(cand))]
    out = m.copy()
    out[c[0] + off[:, 0], c[1] + off[:, 1]] = False
    return out


def bridge(m, a, b, width: int = 1) -> np.ndarray:
    """Join regions a and b of m with a straight strip between their (approximately) closest pixels."""
    pa, pb = np.argwhere(a), np.argwhere(b)
    ia = pa[np.linspace(0, len(pa) - 1, min(len(pa), 400)).astype(int)]
    ib = pb[np.linspace(0, len(pb) - 1, min(len(pb), 400)).astype(int)]
    i, j = np.unravel_index(np.argmin(((ia[:, None] - ib[None]) ** 2).sum(-1)), (len(ia), len(ib)))
    p, q = ia[i].astype(float), ib[j].astype(float)
    out = m.copy()
    for t in np.linspace(0, 1, int(np.ceil(np.abs(q - p).max())) * 2 + 1):
        y, x = np.round(p + t * (q - p)).astype(int)
        out[max(0, y - width // 2): y + (width + 1) // 2, max(0, x - width // 2): x + (width + 1) // 2] = True
    return out


def jitter(m, amp: float, rng, corr: float) -> np.ndarray:
    """Smooth random normal displacement of the boundary by about amp px (correlation length corr)."""
    f = ndimage.gaussian_filter(rng.standard_normal(m.shape), corr, mode="wrap")
    return signed_distance(m) + amp * f / (f.std() or 1.0) > 0.5


def topology_jitter(m, sig, min_area: int, amp: float, corr: float, rng) -> np.ndarray | None:
    """jitter(m, amp, rng, corr), redrawn up to 20 times until its size-filtered signature is `sig` (else None)."""
    for _ in range(20):
        x = jitter(m, amp, rng, corr)
        if signature(x, min_area) == sig:
            return x
    return None


def hull(m) -> np.ndarray:
    out = np.zeros(m.shape, np.uint8)
    cv2.fillConvexPoly(out, cv2.convexHull(cv2.findNonZero(m.astype(np.uint8))), 1)
    return out.astype(bool) | m


def holes_of(m, area: int) -> np.ndarray:
    """Pixels of the holes of m with >= area pixels."""
    h = ndimage.binary_fill_holes(m) & ~m
    lab, n = ndimage.label(h)
    return ((np.bincount(lab.ravel()) >= area)[lab] & h) if n else h


# ---------------------------------------------------------------- terms and advantages of packed rollouts

def pack(mask, cache: ExpertCache) -> tuple:
    return packed(score_mask(mask, cache))


def packed(t) -> tuple:
    """The FIELDS of score_mask's ObjectTerms."""
    return (t.iou_tc, t.f5, t.g_area, float(t.topo_freq > 0), t.topo_freq, sig_id((t.sig,)), float(t.area))


def advantages(T: np.ndarray, legal: np.ndarray | None = None) -> np.ndarray:
    """advantages_v4 of packed groups T (B, G, len(FIELDS)) -> (B, G)."""
    T = np.asarray(T, float)
    terms = {"iou_tc": T[..., IOU], "f5": T[..., F5], "g_area": T[..., GAREA], "topo_any": T[..., TANY],
             "topo_freq": T[..., TFREQ], "sig_id": T[..., SID],
             "legal": np.ones(T.shape[:2]) if legal is None else legal}
    return advantages_v4({k: v.reshape(-1) for k, v in terms.items()}, T.shape[1]).reshape(T.shape[:2])


def cluster_ci(ok, items, rng, n_boot: int = 2000) -> tuple[float, float]:
    """95 % percentile-bootstrap CI of a share, resampling items (clusters of groups)."""
    keys = sorted(set(items))
    pos = {k: i for i, k in enumerate(keys)}
    s = np.zeros((len(keys), 2))
    for o, k in zip(ok, items, strict=True):
        s[pos[k]] += (float(o), 1.0)
    idx = rng.integers(0, len(s), (n_boot, len(s)))
    b = s[idx, 0].sum(1) / s[idx, 1].sum(1)
    return tuple(np.percentile(b, [2.5, 97.5]))


def mean_ci(x, rng, n_boot: int = 4000) -> tuple[float, float, float]:
    """(mean, 95 % CI low, high) by a percentile bootstrap over the entries of x."""
    x = np.asarray(x, float)
    b = x[rng.integers(0, len(x), (n_boot, len(x)))].mean(1)
    lo, hi = np.percentile(b, [2.5, 97.5])
    return float(x.mean()), float(lo), float(hi)


# ---------------------------------------------------------------- honest rollouts (lab scripts/honest.py)

def honest_pools(stack, cache: ExpertCache, rng, mode: str, n: int, pack_fn=pack) -> list[list]:
    """pools[j] = n packed rollouts of rater j whose size-filtered topology is the rater's.

    mode 'real': with probability 0.25 every rater's exact mask, else a smooth boundary perturbation with
    amplitude ~ U(0, tau) and correlation length max(3, tau); 'jit1px': amplitude 1 px, correlation 3 px.
    One amplitude per pool entry is shared by all raters, so every topology class gets the same quality
    (lab ITERATION4 §1); up to 20 draws per rater and 30 amplitudes, then the exact masks.
    """
    k = len(stack)
    exact = [pack_fn(m, cache) for m in stack]
    pools = [[] for _ in range(k)]
    for _ in range(n):
        entry = None
        if not (mode == "real" and rng.random() < 0.25):
            for _amp in range(30):
                amp = 1.0 if mode == "jit1px" else float(rng.uniform(0, cache.tau))
                corr = 3.0 if mode == "jit1px" else max(3.0, cache.tau)
                out = []
                for j, m in enumerate(stack):
                    x = topology_jitter(m, cache.rater_sigs[j], cache.min_area, amp, corr, rng)
                    if x is None:
                        break
                    out.append(x)
                if len(out) == k:
                    entry = [pack_fn(x, cache) for x in out]
                    break
        for j in range(k):
            pools[j].append(entry[j] if entry is not None else exact[j])
    return pools


# ---------------------------------------------------------------- adversaries (lab scripts/adversaries.py)

def band_speck(cache: ExpertCache, rng) -> np.ndarray | None:
    """Vote + a detached compact speck of min_area px inside the band (U minus the dilated vote)."""
    allowed = cache.union & ~ndimage.binary_dilation(cache.vote, np.ones((3, 3), bool), iterations=2)
    off = disk_offsets(cache.min_area)
    cand = np.argwhere(allowed)
    rng.shuffle(cand)
    for c in cand[:400]:
        p = c + off
        if (p.min() >= 0 and p[:, 0].max() < allowed.shape[0] and p[:, 1].max() < allowed.shape[1]
                and allowed[p[:, 0], p[:, 1]].all()):
            out = cache.vote.copy()
            out[p[:, 0], p[:, 1]] = True
            return out
    return None


def adversaries(stack, cache: ExpertCache, rng) -> dict[str, np.ndarray]:
    """dilated_union (2 tau), hull, union, core, band_speck; holes_stamped (vote minus every rater's holes);
    minority_hole (vote with a minority rater's holes, matching that rater's topology)."""
    adv = {"dilated_union": grow(cache.union, 2 * cache.tau), "hull": hull(cache.union), "union": cache.union}
    if cache.core.any():
        adv["core"] = cache.core
    sp = band_speck(cache, rng)
    if sp is not None:
        adv["band_speck"] = sp
    holes = np.zeros_like(cache.vote)
    for m in stack:
        holes |= holes_of(m, cache.min_area)
    if (cache.vote & holes).any():
        adv["holes_stamped"] = cache.vote & ~holes
    vote_sig = signature(cache.vote, cache.min_area)
    for j, m in enumerate(stack):
        if cache.rater_sigs[j] != vote_sig and sum(cache.rater_sigs[j]) > sum(vote_sig):
            h = holes_of(m, cache.min_area)
            cand = cache.vote & ~h
            if (cache.vote & h).any() and signature(cand, cache.min_area) == cache.rater_sigs[j]:
                adv["minority_hole"] = cand
                break
    return adv


# ---------------------------------------------------------------- K9 pairs (lab scripts/k9_discrimination.py, k9_v3.py)

def composite(a: np.ndarray, b: np.ndarray, gap: int) -> np.ndarray:
    """Two-part item: raters paired by index, part a left of part b, `gap` px between the union boxes."""
    k = min(len(a), len(b))

    def box(st):
        ys, xs = np.nonzero(st.any(0))
        return st[:k, ys.min():ys.max() + 1, xs.min():xs.max() + 1]

    a, b = box(a), box(b)
    pad = 12
    out = np.zeros((k, max(a.shape[1], b.shape[1]) + 2 * pad, a.shape[2] + gap + b.shape[2] + 2 * pad), bool)
    out[:, pad:pad + a.shape[1], pad:pad + a.shape[2]] = a
    out[:, pad:pad + b.shape[1], pad + a.shape[2] + gap:pad + a.shape[2] + gap + b.shape[2]] = b
    return out


def correct_candidates(base, rng):
    for d in (-3, -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 3):
        yield grow(base, d)
    for dy in (-2, -1, 0, 1, 2):
        for dx in (-2, -1, 0, 1, 2):
            if dy or dx:
                yield shift(base, dy, dx)
    for amp in (0.5, 1.0, 1.5, 2.0, 3.0):
        for _ in range(2):
            yield jitter(base, amp, rng, 6.0)


def wrong_variants(base, cache: ExpertCache, rng, two_part: bool):
    """Topology defects at and above the size filter A: holes and specks of A and 2A px; bridges (composites)."""
    A = cache.min_area
    for a in (A, 2 * A):
        m = punch_hole(base, a, rng)
        if m is not None:
            yield m
    for a in (A, 2 * A):
        reach = int(np.ceil(np.sqrt(a / np.pi))) + 2
        m = add_speck(base, a, float(rng.uniform(reach + 2, reach + 12)), rng)
        if m is not None:
            yield m
    if two_part:
        lab, n = ndimage.label(base, np.ones((3, 3), bool))
        if n == 2 and len(instance_topology(base).areas) == 2:
            for w in (1, 2):
                yield bridge(base, lab == 1, lab == 2, w)


def vote_iou(m, cache: ExpertCache) -> float:
    inter = np.count_nonzero(m & cache.vote)
    union = np.count_nonzero(m | cache.vote)
    return inter / union if union else 1.0


def k9_groups(stack, rng, two_part: bool = False, pack_fn=pack) -> list[tuple[np.ndarray, bool]]:
    """(group, high_iou) for one item. A group packs 8 rollouts [correct, wrong, 0-3 other wrong variants,
    topology-correct fillers]; the correct partner is a topology-preserving variant whose IoU vs the vote is
    within 0.01 of the wrong one's (lab k9_v3.py); high_iou = both IoUs >= 0.5."""
    cache = ExpertCache.build(stack)
    out = []
    for base in (cache.vote, stack[int(rng.integers(len(stack)))]):
        if not cache.same_sig.get(signature(base, cache.min_area)):
            continue
        cands = [(t, vote_iou(m, cache)) for m in correct_candidates(base, rng) if (t := pack_fn(m, cache))[TANY] == 1]
        wrongs = [(t, vote_iou(m, cache)) for m in wrong_variants(base, cache, rng, two_part)
                  if (t := pack_fn(m, cache))[TANY] == 0]
        if len(cands) < G:
            continue
        for wi, (tw, iw) in enumerate(wrongs):
            near = [j for j, (_, ic) in enumerate(cands) if abs(ic - iw) <= 0.01]
            if not near:
                continue
            jc = near[int(rng.integers(len(near)))]
            others = [t for j, (t, _) in enumerate(wrongs) if j != wi]
            n_w = min(len(others), int(rng.integers(0, 4)))
            fill_w = [others[j] for j in rng.choice(len(others), n_w, replace=False)] if n_w else []
            pool = [j for j in range(len(cands)) if j != jc]
            fill_c = [cands[j][0] for j in rng.choice(pool, G - 2 - n_w, replace=False)]
            out.append((np.array([cands[jc][0], tw, *fill_w, *fill_c]), min(cands[jc][1], iw) >= 0.5))
    return out


# ---------------------------------------------------------------- gate computations (items -> numbers)

def k8_rows(ds: str, key: str, stack) -> list[dict]:
    """Leave-one-rater-out: each rater scored against the cache of the others (k >= 3 only: the remaining
    set must be RL-eligible). Area = mean rater area of the item (lab k8_band.py)."""
    if len(stack) < 3:
        return []
    area = float(stack.reshape(len(stack), -1).sum(1).mean())
    return [dict(subset=subset(ds, key), k=len(stack), area=area, key=key,
                 f5=score_mask(stack[j], ExpertCache.build(np.delete(stack, j, 0))).f5) for j in range(len(stack))]


def k8b_item(stack) -> dict:
    """F5 of the vote moved by c * tau (c = 0.25 .. 2, dilation and erosion) and monotonicity in |c|."""
    cache = ExpertCache.build(stack)
    out = {}
    for sign, name in ((1, "dil"), (-1, "ero")):
        seq = [score_mask(grow(cache.vote, sign * c * cache.tau), cache).f5 for c in (0.25, 0.5, 1.0, 1.5, 2.0)]
        out[name] = seq[-1]
        out[name + "_mono"] = all(b <= a for a, b in zip([1.0, *seq], seq, strict=False))
    return out


def r5_item(stack, seed: int) -> dict[int, list[float]]:
    """F5 of a held-out expert against m random reference raters (m = 2..5, 3 draws)."""
    rng = np.random.default_rng(seed)
    out = {}
    for m in (2, 3, 4, 5):
        out[m] = []
        for _ in range(3):
            pick = rng.permutation(len(stack))[:m + 1]
            out[m].append(score_mask(stack[pick[0]], ExpertCache.build(stack[pick[1:]])).f5)
    return out


def g_area_loo(stack) -> list[tuple[float, float]]:
    """AB4: each rater scored against the remaining raters' area range: (zero tolerance, tau*P tolerance)."""
    out = []
    for j in range(len(stack)):
        rest = np.delete(stack, j, 0)
        if len(rest) < 2:
            continue
        c = ExpertCache.build(rest)
        a = int(stack[j].sum())
        slack = c.tau * c.perimeter
        lo, hi = min(c.rater_areas), max(c.rater_areas)
        out.append((min(1.0, hi / a) * min(1.0, a / lo),
                    min(1.0, (hi + slack) / a) * (min(1.0, a / (lo - slack)) if lo > slack else 1.0)))
    return out


def item_pools(args) -> dict:
    """Per item: cache facts, honest pools (exact, jit1px, real) and the packed adversaries."""
    key, stack, seed, n_pool, pack_fn, with_adv = args
    rng = np.random.default_rng(seed)
    cache = ExpertCache.build(stack)
    sigs = list(cache.rater_sigs)
    vsig = signature(cache.vote, cache.min_area)
    counts = {s: len(r) for s, r in cache.same_sig.items()}
    top = max(counts.values())
    tied = sorted(s for s, c in counts.items() if c == top)
    mode = vsig if vsig in tied else tied[0]
    # minority: True / False for voting raters, None for abstaining ones (in no class of R4)
    out = dict(key=key, k=len(stack), minority=[(s != mode) if j in cache.voters else None for j, s in enumerate(sigs)],
               pools={"exact": [[pack_fn(m, cache)] for m in stack],
                      "jit1px": honest_pools(stack, cache, rng, "jit1px", n_pool, pack_fn),
                      "real": honest_pools(stack, cache, rng, "real", n_pool, pack_fn)})
    if with_adv:
        adv = adversaries(stack, cache, rng)
        out["adv"] = {a: pack_fn(x, cache) for a, x in adv.items()}
        out["same_class"] = {a: [j for j in range(len(stack)) if sigs[j] == signature(x, cache.min_area)]
                             for a, x in adv.items()}
    return out


def adversary_groups(data: list[dict], name: str, mode: str, rng, n_groups: int = 20):
    """(T (N, G, F), item ids, same-class honest positions) for one adversary: each group is the adversary
    plus 7 honest rollouts (raters drawn uniformly; for hole adversaries the first honest rollout is a rater
    of the adversary's own topology class)."""
    T, items, same = [], [], []
    for d in data:
        if name not in d.get("adv", {}):
            continue
        sc = d["same_class"][name]
        holes = name in HOLES
        if holes and not sc:
            continue  # a topology no rater drew is invalid and loses trivially
        for _ in range(n_groups):
            draws = rng.integers(0, d["k"], G - 1)
            if holes:
                draws[0] = sc[int(rng.integers(len(sc)))]
            pools = d["pools"][mode]
            T.append([d["adv"][name]] + [pools[r][rng.integers(len(pools[r]))] for r in draws])
            items.append(d["key"])
            same.append([1 + g for g, r in enumerate(draws) if r in sc])
    return np.array(T, float), items, same


def in_steps(adv_fn, T: np.ndarray, step: int = 32) -> np.ndarray:
    """Advantages of many groups, computed one GRPO step (32 groups) at a time."""
    return np.concatenate([adv_fn(T[s:s + step]) for s in range(0, len(T), step)]) if len(T) else np.zeros((0, G))


def neutrality(data: list[dict], mode: str, rng, adv_fn, n_batch: int = 300) -> np.ndarray:
    """Per item: mean advantage of minority-topology rollouts minus mode-class rollouts, groups of 8 drawn
    i.i.d. from the raters (honest rollouts of `mode`), one group per item per batch."""
    acc = [[[], []] for _ in data]
    for _ in range(n_batch):
        T, cls = [], []
        for d in data:
            draws = rng.integers(0, d["k"], G)
            pools = d["pools"][mode]
            T.append([pools[r][rng.integers(len(pools[r]))] for r in draws])
            cls.append([d["minority"][r] for r in draws])
        a = in_steps(adv_fn, np.array(T, float))
        for b, c in enumerate(cls):
            for g in range(G):
                if c[g] is not None:
                    acc[b][int(c[g])].append(a[b, g])
    return np.array([np.mean(x[1]) - np.mean(x[0]) for x in acc if x[0] and x[1]])


MOVES = (("dilate 1px", 1.0, 0.0), ("dilate 0.25tau", 0.0, 0.25), ("dilate 0.5tau", 0.0, 0.5),
         ("erode 1px", -1.0, 0.0), ("erode 0.25tau", 0.0, -0.25), ("erode 0.5tau", 0.0, -0.5))


def local_gradient_item(args) -> dict:
    """Groups of 8 realistic honest rollouts; rollout 0 is dilated / eroded (others fixed) and the change of
    its own advantage recorded per move. Moves that leave the mask unchanged are skipped."""
    key, stack, seed, n_groups = args
    rng = np.random.default_rng(seed)
    cache = ExpertCache.build(stack)
    masks = [[] for _ in stack]
    for j, m in enumerate(stack):  # realistic honest masks, as honest_pools 'real' (one amplitude per entry)
        for _ in range(6):
            amp = 0.0 if rng.random() < 0.25 else float(rng.uniform(0, cache.tau))
            x = topology_jitter(m, cache.rater_sigs[j], cache.min_area, amp, max(3.0, cache.tau), rng) if amp else m
            if x is not None:
                masks[j].append(x)
        masks[j] = masks[j] or [m]
    out = {name: [] for name, _, _ in MOVES}
    away = {name: [] for name, _, _ in MOVES}  # the move takes the area away from the median rater area
    med = np.median(cache.rater_areas)
    for _ in range(n_groups):
        group = [masks[r][rng.integers(len(masks[r]))] for r in rng.integers(0, len(stack), G)]
        base = np.array([[pack(m, cache) for m in group]])
        a0 = advantages(base)[0, 0]
        for name, px, frac in MOVES:
            moved = grow(group[0], px if px else frac * cache.tau)
            if np.array_equal(moved, group[0]):
                continue
            T = base.copy()
            T[0, 0] = pack(moved, cache)
            out[name].append(advantages(T)[0, 0] - a0)
            area = base[0, 0, AREA]
            away[name].append(area >= med if name.startswith("dilate") else area <= med)
    return dict(key=key, deltas=out, away=away)
