"""Synthetic probe items with exact ground truth (precision probe, A4 section 6; synthetic topology probe, H-A1-7).

Generator (deterministic per seed):
- equivalent diameter drawn from the LIDC train voting raters' masks at frame resolution (census.py);
- shape families on the 0-1000 grid of an S x S frame, rasterized by polycore (objdel.geometry.to_mask), so the
  ground truth raster has no tie error:
  * "fourier": star polygon r(theta) = r0 (1 + E e cos(2 theta + phi_2) + A sum_{k=3..12} c_k cos(k theta + phi_k)),
    e ~ U(0.5, 1), c_k ~ U(0, 1) / k, plus spicules (narrow radial bumps) with probability 0.3. (E, A) is calibrated
    (H-A4-7) to two real-LIDC targets: the centroid-disc IoU and the certified targets' vertex count (a one-parameter
    version, all harmonics on A, needed 16 vertices to reach the disc IoU, against 11 on real targets: run 1);
  * "lobulated": union of 2-4 discs around the centre;
  * "hole": a fourier outer with a star-shaped cavity carved inside (cavitation);
  * "two": two fourier parts separated by >= 2 px;
- frame: the targets.frame side rule (union box side + 2 margins, smallest divisor of 1000 >= need, >= 50 px) with
  the margin drawn from [8 px, 60 %] (scale jitter), then the object placed uniformly with >= 2 px clearance
  (translation jitter): the frame does not give the position away (A4 S5-S6);
- the certified target (what SFT teaches) is objdel.targets.build_target on the ground-truth raster, with the
  reward's size rule min_area(|GT|);
- every item's ground-truth signature at that rule must equal its family's intended topology, else it is redrawn.

Renderer: a real LIDC lung-window background patch (HU), the object composited with a partial-volume blur
(sigma 0.5-1 px) as solid (-100..+100 HU) or ground-glass (-700..-300 HU), CT noise sigma 10-30 HU, then the frame's
lung window [-1000, 400] and the Pillow bicubic upscale to 512 px, as objdel.targets.frame does for real items.

Metrics (frame pixels): boundary = foreground pixels with a 4-adjacent background pixel (outside the frame counts as
background): MONAI get_mask_edges (binary_erosion XOR mask) and A4 S2's convention. Distances are exact Euclidean
(scipy distance_transform_edt) from each boundary pixel to the other boundary set.
- ASSD: pooled mean over both directed sets, (sum d(a, B) + sum d(b, A)) / (|A| + |B|) (Yeghiazaryan & Voiculescu
  2018; MONAI's symmetric average surface distance concatenates both directions before the mean). DeepMind's
  surface-distance returns the two directed means instead, surfel-area weighted; the two differ by implementation
  (Podobnik & Vrtovec 2024), so the convention is stated with every number.
- HD95: max of the two directed 95th percentiles (np.percentile, linear interpolation), as MONAI and Metrics
  Reloaded; NSD@tau: share of both boundary sets within tau px of the other (Nikolov et al.).
- IoU, centroid error (pixel-centre centroids), relative equivalent-diameter error, topology = size-filtered
  signature at min_area(|GT|) equal to the ground truth's.
"""
import numpy as np
from PIL import Image
from scipy import ndimage

from objdel.checker import signature, topology
from objdel.geometry import GRID, canonicalize, parse, serialize, to_mask
from objdel.targets import build_target
from objdel.targets.frame import DIVISORS, LUNG_WINDOW

SIG = {"fourier": (0,), "lobulated": (0,), "hole": (1,), "two": (0, 0)}
CROSS = ndimage.generate_binary_structure(2, 1)


# ---------------------------------------------------------------- generator
def _star(cx, cy, r0, shape, rng, n=96, spikes=True):
    """Star polygon; shape = (E, A): E scales the 2nd harmonic (elongation), A the 3rd-12th (irregularity)."""
    E, A = shape
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    rr = 1 + E * rng.uniform(0.5, 1) * np.cos(2 * th + rng.uniform(0, 2 * np.pi))
    for k in range(3, 13):
        rr += A * rng.uniform(0, 1) / k * np.cos(k * th + rng.uniform(0, 2 * np.pi))
    if spikes and rng.random() < 0.3:
        for _ in range(rng.integers(1, 5)):
            c = rng.uniform(0, 2 * np.pi)
            dth = np.angle(np.exp(1j * (th - c)))
            rr += rng.uniform(0.2, 0.6) * max(A, E) * np.exp(-0.5 * (dth / 0.08) ** 2)
    rr = np.maximum(rr, 0.35)
    return np.c_[cx + r0 * rr * np.cos(th), cy + r0 * rr * np.sin(th)]


def _raster(poly_px, S):
    """Pixel-unit polygon -> grid ring (deduplicated) -> polycore raster on S x S; None if illegal."""
    g = np.rint(np.asarray(poly_px) * (GRID / S)).astype(np.int64)
    g = np.clip(g, 0, GRID)
    g = g[np.any(g != np.roll(g, 1, 0), 1)]
    if len(g) < 3:
        return None
    parts = canonicalize({"x": [(g, [])]})["x"]
    geo = to_mask(parts, S, S)
    return geo.mask if geo.legal else None


def _prims(family, d, shape, rng):
    """The object as primitives [(+1 union | -1 carve, pixel-unit polygon centred on 0)]; shape = (E, A)."""
    E, A = shape
    r0 = d / 2
    if family == "fourier":
        return [(1, _star(0, 0, r0, (E, A), rng))]
    if family == "lobulated":
        out = []
        for _ in range(rng.integers(2, 5)):
            a, rr, off = rng.uniform(0, 2 * np.pi), r0 * rng.uniform(0.5, 0.85), r0 * rng.uniform(0.2, 0.55)
            out.append((1, _star(off * np.cos(a), off * np.sin(a), rr, (E * 0.5, A * 0.3), rng, spikes=False)))
        return out
    if family == "hole":
        hr = r0 * rng.uniform(0.25, 0.45)
        return [(1, _star(0, 0, r0, (E, A * 0.5), rng, spikes=False)),
                (-1, _star(rng.uniform(-.15, .15) * r0, rng.uniform(-.15, .15) * r0, hr, (0.1, 0.2), rng, n=48,
                           spikes=False))]
    if family == "two":
        d1 = d * rng.uniform(0.55, 0.8)
        d2 = np.sqrt(max(d * d - d1 * d1, 9.0))
        a = rng.uniform(0, 2 * np.pi)
        sep = (d1 + d2) / 2 * (1 + 0.5 * max(A, E)) + rng.uniform(2.5, 6)
        return [(1, _star(-sep / 2 * np.cos(a), -sep / 2 * np.sin(a), d1 / 2, (E, A), rng, spikes=False)),
                (1, _star(sep / 2 * np.cos(a), sep / 2 * np.sin(a), d2 / 2, (E, A), rng, spikes=False))]
    raise ValueError(family)


def _compose(prims, S, dx, dy):
    m = np.zeros((S, S), bool)
    for sign, poly in prims:
        x = _raster(poly + [dx, dy], S)
        if x is None:
            return None
        m = (m | x) if sign > 0 else (m & ~x)
    return m


def _frame_side(side, rng):
    margin = max(8, int(np.ceil(rng.uniform(0.0, 0.60) * side)))
    need = max(side + 2 * margin, 50)
    return next(x for x in DIVISORS if x >= need)


def make_item(family, d, shape, rng, tries=50):
    """One probe item: frame side S, ground-truth raster, its certified target text, family; shape = (E, A)."""
    for _ in range(tries):
        prims = _prims(family, d, shape, rng)
        pts = np.concatenate([p for s, p in prims if s > 0])
        lo, hi = pts.min(0), pts.max(0)
        side = int(np.ceil((hi - lo).max())) + 1
        S = _frame_side(side, rng)
        if side + 6 > S or S > 500:
            continue
        # translation jitter: the polygon box lands uniformly inside the frame with >= 2 px clearance
        x0, y0 = rng.uniform(2, S - 2 - (hi - lo)[0]), rng.uniform(2, S - 2 - (hi - lo)[1])
        gt = _compose(prims, S, x0 - lo[0], y0 - lo[1])
        if gt is None or not gt.any():
            continue
        ma = topology.min_area(int(gt.sum()))
        if signature(gt, ma) != SIG[family] or signature(gt) != SIG[family]:
            continue
        try:
            t = build_target({"nodule": gt}, ("nodule",), min_area={"nodule": ma})
        except ValueError:
            continue
        return {"family": family, "S": S, "gt": gt, "min_area": ma, "text": t.text, "objs": t.objs["nodule"],
                "diam": float(2 * np.sqrt(gt.sum() / np.pi)), "d_drawn": float(d)}
    return None


# ---------------------------------------------------------------- metrics
def boundary(m):
    return m & ~ndimage.binary_erosion(m, structure=CROSS, border_value=0)


def surface(gt, pred):
    """ASSD (pooled), HD95, NSD@1/@2 between the boundaries; None when either mask is empty."""
    bg, bp = boundary(gt), boundary(pred)
    if not bg.any() or not bp.any():
        return None
    dg = ndimage.distance_transform_edt(~bp)[bg]  # each GT boundary pixel to the predicted boundary
    dp = ndimage.distance_transform_edt(~bg)[bp]
    d = np.concatenate([dg, dp])
    return {"assd": float(d.mean()), "hd95": float(max(np.percentile(dg, 95), np.percentile(dp, 95))),
            "nsd1": float((d <= 1).mean()), "nsd2": float((d <= 2).mean())}


def surface_brute(gt, pred):
    a = np.argwhere(boundary(gt)).astype(float)
    b = np.argwhere(boundary(pred)).astype(float)
    D = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))
    dg, dp = D.min(1), D.min(0)
    d = np.concatenate([dg, dp])
    return {"assd": float(d.mean()), "hd95": float(max(np.percentile(dg, 95), np.percentile(dp, 95)))}


def centroid(m):
    ys, xs = np.nonzero(m)
    return np.array([xs.mean() + 0.5, ys.mean() + 0.5])


def score(item, text):
    """Probe metrics of one answer text against the item's ground truth."""
    p = parse(text, ("nodule",))
    if not p.ok:
        return {"parse": False, "legal": False}
    geo = to_mask(p.objs["nodule"], item["S"], item["S"])
    out = {"parse": True, "legal": bool(geo.legal)}
    if not geo.legal:
        return out
    gt, pr = item["gt"], geo.mask
    u = np.count_nonzero(gt | pr)
    out["iou"] = np.count_nonzero(gt & pr) / u
    out["topology"] = signature(pr, item["min_area"]) == signature(gt, item["min_area"])
    if not pr.any():
        out.update(assd=np.inf, hd95=np.inf, nsd1=0.0, nsd2=0.0, centroid=np.inf, rel_diam=-1.0)
        return out
    out.update(surface(gt, pr))
    out["centroid"] = float(np.linalg.norm(centroid(pr) - centroid(gt)))
    out["rel_diam"] = float(np.sqrt(pr.sum() / gt.sum()) - 1)
    return out


def pass_rule(rows, rule, bins=(6.0, 16.0)):
    """A4 pass rule on the items with ground-truth diameter in [lo, hi): legal share over all of them; medians of
    ASSD / IoU / centroid error over legal answers; topology over the legal hole and two-part answers (an illegal
    answer counts as a topology miss)."""
    rs = [r for r in rows if bins[0] <= r["diam"] < bins[1]]
    lg = [r for r in rs if r["legal"]]
    topo = [r for r in rs if r["family"] in ("hole", "two")]
    m = {"n": len(rs), "legal": float(np.mean([r["legal"] for r in rs])) if rs else np.nan,
         "assd": float(np.median([r["assd"] for r in lg])) if lg else np.inf,
         "iou": float(np.median([r["iou"] for r in lg])) if lg else 0.0,
         "centroid": float(np.median([r["centroid"] for r in lg])) if lg else np.inf,
         "topology": float(np.mean([r.get("topology", False) for r in topo])) if topo else np.nan,
         "n_topology": len(topo)}
    ok = {"legal": m["legal"] >= rule["legal"], "assd": m["assd"] <= rule["assd"], "iou": m["iou"] >= rule["iou"],
          "centroid": m["centroid"] <= rule["centroid"],
          "topology": (m["topology"] >= rule["topology"]) if topo else True}
    m["checks"] = ok
    m["pass"] = all(ok.values())
    return m


# ---------------------------------------------------------------- synthetic models (answer texts)
def _text(parts):
    return serialize(canonicalize({"nodule": parts}), ("nodule",))


def model_exact(item, rng):
    return item["text"]


def model_noise(sigma):
    def f(item, rng):
        k = GRID / item["S"]
        parts = [(np.rint(o + rng.normal(0, sigma * k, o.shape)).astype(np.int64),
                  [np.rint(h + rng.normal(0, sigma * k, h.shape)).astype(np.int64) for h in hs])
                 for o, hs in item["objs"]]
        parts = [(np.clip(o, 0, GRID), [np.clip(h, 0, GRID) for h in hs]) for o, hs in parts]
        return _text(parts)
    return f


def model_shift(item, rng):
    k = GRID // item["S"]
    dx, dy = [(k, 0), (-k, 0), (0, k), (0, -k)][rng.integers(4)]
    return _text([(o + [dx, dy], [h + [dx, dy] for h in hs]) for o, hs in item["objs"]])


def shift_text(item, dx, dy):
    k = GRID // item["S"]
    return _text([(o + [dx * k, dy * k], [h + [dx * k, dy * k] for h in hs]) for o, hs in item["objs"]])


def _disc_text(cx, cy, r, S, n=48):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    g = np.rint(np.c_[cx + r * np.cos(th), cy + r * np.sin(th)] * (GRID / S)).astype(np.int64)
    g = np.clip(g, 0, GRID)
    g = g[np.any(g != np.roll(g, 1, 0), 1)]
    return _text([(g, [])])


def model_centroid_disc(item, rng):
    cx, cy = centroid(item["gt"])
    return _disc_text(cx, cy, np.sqrt(item["gt"].sum() / np.pi), item["S"])


def model_centred_prior(prior_diam):
    def f(item, rng):
        S = item["S"]
        d = prior_diam.get(S) or prior_diam[min(prior_diam, key=lambda s: abs(s - S))]
        return _disc_text(S / 2, S / 2, d / 2, S)
    return f


def model_a0(item, rng):
    """A0: holes emptied (outer rings only), parts kept."""
    return _text([(o, []) for o, _ in item["objs"]])


# ---------------------------------------------------------------- renderer
def render(item, bg_hu, rng, image_side=512):
    """uint8 RGB model image (image_side) and the noise-free HU composite (frame resolution)."""
    S, gt = item["S"], item["gt"]
    y0, x0 = rng.integers(0, bg_hu.shape[0] - S + 1), rng.integers(0, bg_hu.shape[1] - S + 1)
    bg = bg_hu[y0:y0 + S, x0:x0 + S].astype(float)
    hu = rng.uniform(-100, 100) if rng.random() < 0.7 else rng.uniform(-700, -300)
    alpha = ndimage.gaussian_filter(gt.astype(float), rng.uniform(0.5, 1.0))
    clean = bg * (1 - alpha) + hu * alpha
    noisy = clean + rng.normal(0, rng.uniform(10, 30), clean.shape)
    lo, hi = LUNG_WINDOW
    u8 = np.clip(np.rint((noisy - lo) / (hi - lo) * 255), 0, 255).astype(np.uint8)
    img = Image.fromarray(np.repeat(u8[..., None], 3, -1)).resize((image_side, image_side), Image.Resampling.BICUBIC)
    return np.asarray(img), clean, hu, float(np.median(bg))
