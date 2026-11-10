"""Data quality gates of the LIDC store and the scheduled RL rows (Implementation.md B1.3, B1.7; gates DG1-DG18).

Gates (B.6.1; `gates(report)` returns the 21 in this order, `failures` the hard ones that fail):
- DG1, DG2, DG3 `leakage(store_root, ...)`: patient, series and item overlap between splits; SHA-256 of every int16
  slice and of the 32x32 crop at the union-box centre shared across patients or splits; near-duplicate crops across
  splits (`near_duplicates`), each needing a review in `docs/reference/neardup_review.json`.
- DG1b `rows_leakage`: every row a train item of the store, no val or test patient, no val or test view image, no
  excluded item (and, on the Mac, the manifest's exclusions SHA-256 equal to the file's).
- DG4, DG5, DG6, DG13, DG18 `window_checks(rows_root, store_root)`: each row's window box against its item (doctor
  pixels inside, clearance, other outlined nodules, padding) and a fresh render (box, side, image SHA-256 of the
  render and of the PNG on disk).
- DG7-DG11, DG14, DG16 `schedule_report`: coverage of the off-centre rows against the train shifted law (the pool's
  draw-0 shifted windows) on every block prefix, the full mix (DG9b), policy-side share, exact block counts, the
  report numbers and the val coverage with its train-patient null.
- DG12 `tokens`, DG14b `agreement`, DG15 `frozen`, DG17 `dg17` (contour-to-slice check, below).

Contour-to-slice check (B1.7, DG17). `lidc.slice_index` puts each contour on the scan slice nearest to its z, so a
series with an ambiguous z can put a contour on the wrong slice. Per series of the store, from the XML and the fetch
sidecars only (no DICOM):
- (a) ``shared_z``: one z cited by two or more SOP UIDs (two images at one position);
- (b) ``far_rois``: a contour z farther than half the slice spacing (median gap of the scan's z values) from every
  slice;
- (d) ``sidecar``: the fetch sidecar has fewer or more kept slices than cited SOP UIDs, or matched slices by z
  (``lidc.ct_slices`` drops a slice key that is not within 0.01 mm of a scan slice, so the contour loses its image).
Every series with a hit goes to step 2 (the DICOM header check, `dicom_z`). Train hits must be dropped
(`train_exclusions.json`) or mapped (fixed, store rebuilt, re-run lists no hit); val and test hits keep their items
in the frozen views and are listed as known label issues.

No torch and no transformers here, not even lazily: callers build the tokenizer and pass it in.

References:
- LIDC-IDRI known issues ("Annotation and Markup Issues" item 6: LIDC-IDRI-0396, 139.xml, SOP Instance UID at
  position 1420, fixed 2018-06-28): https://www.cancerimagingarchive.net/collection/lidc-idri/
- Third-party list of duplicate-z series (0085, 0146, 0418, 0572, 0979):
  https://huggingface.co/datasets/YongchengYAO/LIDC-IDRI-Lite
- DICOM ImagePositionPatient (0020,0032): https://dicom.innolitics.com/ciods/ct-image/image-plane/00200032
- Two-sample Kolmogorov-Smirnov: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html
"""
import hashlib
import json
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from math import ceil
from pathlib import Path

import cv2
import numpy as np
import pyarrow.parquet as pq
from scipy import ndimage as ndi

from objdel.data import lidc, store, windows
from objdel.data.targets.frame import CLEARANCE, crowded, frame, hu_window, union_box
from objdel.hashing import file_sha256, pixels_sha256
from objdel.io import read_json

from .coverage import cell_coverage, mix_coverage
from .store import SPLITS
from .windows import KINDS, Window, cell_shares, offset_of

KNOWN_0396 = ("LIDC-IDRI-0396", 1420.0)  # TCIA known issue 6: patient, contour z of the corrected SOP UID


def _sidecar_hit(path: Path) -> dict | None:
    if not path.exists():
        return None
    j = read_json(path)
    s = {"slices": j["slices"], "annotated_sops": j["annotated_sops"], "matched_by_z": j.get("matched_by_z", 0)}
    return s if s["slices"] != s["annotated_sops"] or s["matched_by_z"] > 0 else None


def z_report(store_root, xml_root, sidecars, scans=None, previous=None) -> dict:
    """One entry per LIDC series of the store at ``store_root`` (split, patient, item count, hits (a), (b), (d)).

    ``step2``: the series with a hit; ``known_label_issues``: the item ids of val and test hit series; ``mapped``: the
    train hit series of ``previous`` (an older report) that list no hit now; ``check_0396``: the SOP UIDs that
    LIDC-IDRI-0396's XML cites at z 1420.0 and whether the fetch found every cited SOP UID (by UID, none by z)."""
    store_root, sidecars = Path(store_root), Path(sidecars)
    scans = lidc.scans() if scans is None else scans
    sop_z = lidc.sop_z(xml_root)
    series, item_ids = {}, {}
    for split in SPLITS:
        for row in pq.read_table(store_root / "lidc" / f"{split}.parquet", columns=["item_id", "group_id"]).to_pylist():
            s = row["item_id"].split("/")[1]
            series.setdefault(s, {"split": split, "patient": row["group_id"], "items": 0})["items"] += 1
            item_ids.setdefault(s, []).append(row["item_id"])
    for s, e in series.items():
        by_z = {}
        for sop, z in sop_z.get(s, {}).items():
            by_z.setdefault(z, []).append(sop)
        e["shared_z"] = {str(z): sorted(sops) for z, sops in sorted(by_z.items()) if len(sops) > 1}
        zvals = scans[s].zvals
        half = float(np.median(np.diff(zvals))) / 2
        e["far_rois"] = [{"sop": sop, "z": z, "distance_mm": round(d, 2), "spacing_mm": round(2 * half, 2)}
                         for sop, z in sorted(sop_z.get(s, {}).items(), key=lambda kv: (kv[1], kv[0]))
                         if (d := float(np.abs(zvals - z).min())) > half]
        e["sidecar"] = _sidecar_hit(sidecars / f"{s}.json")
    hits = sorted(s for s, e in series.items() if e["shared_z"] or e["far_rois"] or e["sidecar"])
    for s in hits:
        series[s]["item_ids"] = sorted(item_ids[s])
    old = [s for s in (previous or {}).get("step2", []) if previous["series"][s]["split"] == "train"]
    patient, z0396 = KNOWN_0396
    c = next((s for s, sc in scans.items() if sc.patient == patient), None)
    side = read_json(sidecars / f"{c}.json") if c else {}
    return {"train_sha256": file_sha256(store_root / "lidc" / "train.parquet"), "series": series, "step2": hits,
            "mapped": sorted(s for s in old if s not in hits),
            "known_label_issues": {sp: sorted(i for s in hits if series[s]["split"] == sp for i in item_ids[s])
                                   for sp in ("val", "test")},
            "check_0396": c and {"series": c,
                                 "sop_uids_at_1420": sorted(u for u, z in sop_z.get(c, {}).items() if z == z0396),
                                 "all_cited_sops_fetched": side.get("slices") == side.get("annotated_sops")
                                 and side.get("matched_by_z", 0) == 0}}


def train_hits_not_excluded(report, exclusions) -> dict[str, list[str]]:
    """Train hit series -> its item ids that are not keys of ``exclusions`` (series with none left out)."""
    out = {}
    for s in report["step2"]:
        if report["series"][s]["split"] == "train" and (m := [i for i in report["series"][s]["item_ids"]
                                                              if i not in exclusions]):
            out[s] = m
    return out


def dg17(report, exclusions, train_sha256) -> list[str]:
    """Gate DG17: failure lines, [] passes. A report on another train store is stale; a train hit is resolved only
    when every store item of its series is a key of ``exclusions`` (a mapped series lists no hit any more). Val and
    test hits never fail it (they are known label issues)."""
    if report["train_sha256"] != train_sha256:
        return [f"stale report: train store {train_sha256}, report {report['train_sha256']}"]
    return [f"{s}: {len(m)} of {report['series'][s]['items']} items not in train_exclusions.json"
            for s, m in train_hits_not_excluded(report, exclusions).items()]


def dicom_z(folder, zvals) -> dict:
    """Step 2 on the DICOM headers in ``folder`` (``*.dcm``, untrusted: run under ``python -I``): z positions
    (ImagePositionPatient[2], 0.01 mm) shared by two or more SOP Instance UIDs, and the scan z values that
    `lidc.slice_index` maps contours to but no DICOM has."""
    import pydicom  # [data] extra

    by_z = {}
    for p in sorted(Path(folder).glob("*.dcm")):
        ds = pydicom.dcmread(p, stop_before_pixels=True)
        by_z.setdefault(round(float(ds.ImagePositionPatient[2]), 2), []).append(str(ds.SOPInstanceUID))
    return {"duplicate_z": {str(z): sorted(u) for z, u in sorted(by_z.items()) if len(u) > 1},
            "slice_z_without_dicom": [float(z) for z in zvals if not any(abs(z - d) < 0.01 for d in by_z)]}


# --- B1.3: the data gates -----------------------------------------------------------------------------------------
NEAR = 0.98  # DG3 score limit (Pearson correlation of two 32x32 HU crops)
CROSS = (("train", "val"), ("train", "test"), ("val", "test"))
FROZEN = ("lidc-test-shifted", "lidc-val-edge", "lidc-val-shifted")  # docs/reference/frozen_views.json (B1.6)
SCHEDULE_ARGS = ("per_step", "block_steps", "shares", "disagreement_weight")  # manifest.schedule keys of the report
LOW_CONTRAST = 5  # grey levels between the doctors' union and the ring around it (DG14)
NULL_DRAWS = 500  # DG16 null: random train-patient subsamples


@dataclass(frozen=True)
class Gate:
    """One quality gate: its measured value, the limit as text, whether it is hard and whether it passed."""
    name: str
    value: float
    limit: str
    hard: bool
    ok: bool


def records(rows_root) -> list[dict]:
    """The rows root's records (lidc/train.jsonl), in schedule order."""
    return [json.loads(x) for x in Path(rows_root, "lidc", "train.jsonl").read_text().splitlines()]


def crop32(image: np.ndarray, ub) -> np.ndarray:
    """int16 (32, 32) HU crop at the union-box centre, padded with -1,000 HU below and right at the slice edge (the
    b1_pixels.py rule)."""
    cx, cy = (ub[0] + ub[2]) // 2, (ub[1] + ub[3]) // 2
    c = image[max(cy - 16, 0):cy + 16, max(cx - 16, 0):cx + 16]
    return np.pad(c, ((0, 32 - c.shape[0]), (0, 32 - c.shape[1])), constant_values=-1000).astype(np.int16)


def unit(crop: np.ndarray) -> np.ndarray:
    """The crop minus its mean over its norm (+1e-6), float32: the dot product of two is the DG3 score."""
    x = crop.astype(np.float32).reshape(-1)
    x = x - x.mean()
    return x / (np.linalg.norm(x) + 1e-6)


def near_duplicates(meta: Sequence[tuple[str, str, str]], units: np.ndarray) -> dict:
    """DG3 on (item_id, split, group_id) rows and their unit crops: every cross-split pair (train-val, train-test,
    val-test) with score >= NEAR, and (information only) the count of train items with such a train partner of
    another patient."""
    split, group = np.array([m[1] for m in meta]), np.array([m[2] for m in meta])
    pairs = []
    for a, b in CROSS:
        ia, ib = np.flatnonzero(split == a), np.flatnonzero(split == b)
        for i, j in np.argwhere(units[ia] @ units[ib].T >= NEAR):
            x, y = ia[i], ib[j]
            pairs.append({"a": meta[x][0], "b": meta[y][0], "splits": f"{a}-{b}",
                          "score": float(units[x] @ units[y])})
    it, within = np.flatnonzero(split == "train"), 0
    for k in range(0, len(it), 1024):
        c = units[it[k:k + 1024]] @ units[it].T
        c[group[it[k:k + 1024]][:, None] == group[it][None, :]] = -1
        within += int((c.max(1) >= NEAR).sum())
    return {"pairs": pairs, "within_train": within}


def reviewed(pair: dict, review: Mapping[str, str], exclusions: Mapping) -> bool:
    """A near pair is reviewed by the owner's OK, or by "drop train item" with its train item in the exclusions."""
    v = review.get(f"{pair['a']} {pair['b']}", "")
    return v.startswith("owner OK") or (v == "drop train item" and pair["splits"].startswith("train")
                                        and pair["a"] in exclusions)


def _png(path: Path, a: np.ndarray, b: np.ndarray) -> None:
    img = np.hstack([hu_window(a), np.full((32, 4), 255), hu_window(b)]).astype(np.uint8)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not cv2.imwrite(str(path), cv2.resize(img, None, fx=4, fy=4, interpolation=cv2.INTER_NEAREST)):
        raise OSError(f"{path}: cv2.imwrite failed")


def _shared(holders: Mapping[str, set]) -> int:
    """Hashes held by two (split, patient) pairs of different patients or splits."""
    return sum(len({s for s, _ in v}) > 1 or len({g for _, g in v}) > 1 for v in holders.values())


def leakage(store_root, review: Mapping[str, str], exclusions: Mapping, out) -> dict:
    """DG1, DG2, DG3 on every LIDC item of the store: counts of patients, series and items in two splits; of slice
    and crop SHA-256 values held by two patients or two splits; the cross-split near pairs (each with a side-by-side
    PNG under `out`/neardup and its review state). A store without a split raises (store.load), never passes."""
    meta, slices, crops, units, raw = [], {}, {}, [], []
    for split in SPLITS:
        for it in store.load("lidc", split, store_root):
            c = crop32(it.image, union_box(it))
            meta.append((it.item_id, split, it.group_id))
            slices.setdefault(pixels_sha256(it.image), set()).add((split, it.group_id))
            crops.setdefault(pixels_sha256(c), set()).add((split, it.group_id))
            units.append(unit(c))
            raw.append(c)
    keys = {"patients": lambda m: m[2], "series": lambda m: m[0].split("/")[1], "items": lambda m: m[0]}
    out_doc = {}
    for name, key in keys.items():
        where = {}
        for m in meta:
            where.setdefault(key(m), set()).add(m[1])
        out_doc[name] = sum(len(s) > 1 for s in where.values())
    near = near_duplicates(meta, np.stack(units))
    index = {m[0]: i for i, m in enumerate(meta)}
    for n, p in enumerate(near["pairs"]):
        p["png"] = f"neardup/{n:03d}.png"
        _png(Path(out, p["png"]), raw[index[p["a"]]], raw[index[p["b"]]])
        p["reviewed"] = reviewed(p, review, exclusions)
    return out_doc | {"slice_sha": _shared(slices), "crop_sha": _shared(crops), "near": near["pairs"],
                      "within_train": near["within_train"], "train_items": sum(m[1] == "train" for m in meta)}


def view_image_hashes(eval_views_root) -> set[str]:
    """image_sha256 of every record of every val and test view under `eval_views_root` (DG1b)."""
    return {json.loads(x)["image_sha256"] for p in sorted(Path(eval_views_root).glob("*/lidc/*.jsonl"))
            if p.stem in ("val", "test") for x in p.read_text().splitlines()}


def rows_leakage(rows_root, store_root, view_hashes: set[str], excluded: frozenset[str]) -> dict:
    """DG1b: the window ids of rows whose source item is not a train item of the store, whose patient is a val or
    test patient, whose image equals a val or test view image, or whose source item is in `excluded`."""
    ids = {}
    for split in SPLITS:
        t = pq.read_table(Path(store_root, "lidc", f"{split}.parquet"), columns=["item_id", "group_id"]).to_pylist()
        ids[split] = ({r["item_id"] for r in t}, {r["group_id"] for r in t})
    recs = records(rows_root)
    return {"not_train": [r["item_id"] for r in recs if r["source_item_id"] not in ids["train"][0]],
            "val_test_patient": [r["item_id"] for r in recs if r["group_id"] in ids["val"][1] | ids["test"][1]],
            "view_image": [r["item_id"] for r in recs if r["image_sha256"] in view_hashes],
            "excluded": [r["item_id"] for r in recs if r["source_item_id"] in excluded]}


def _contrast(it) -> float:
    """Grey levels (lung window) between the doctors' union and the ring 1-4 px around it (on the union box + 5 px:
    the same pixels as on the whole slice, 14 -> 1 ms)."""
    x0, y0, x1, y1 = union_box(it)
    box = slice(max(y0 - 5, 0), y1 + 5), slice(max(x0 - 5, 0), x1 + 5)
    m = np.logical_or.reduce(it.masks["nodule"], axis=0)[box]
    ring = ndi.binary_dilation(m, iterations=4) & ~ndi.binary_dilation(m, iterations=1)
    g = hu_window(it.image[box])
    return float(g[m].mean() - g[ring].mean()) if ring.any() else float("nan")


def window_checks(rows_root, store_root) -> dict:
    """DG4, DG5, DG6, DG13, DG18 over every row: its box against its source item (all doctor pixels inside;
    clearance >= CLEARANCE and no other outlined nodule, the policy frame included; inside the slice), its view flags,
    and a fresh render (frame) whose box, side and image SHA-256 must equal the record's and the PNG's; plus the
    window ids of low-contrast items (DG14)."""
    rows_root, recs = Path(rows_root), records(rows_root)
    by_src = {}
    for r in recs:
        by_src.setdefault(r["source_item_id"], []).append(r)
    out = dict.fromkeys(("inside", "clear", "fallback", "kept_side", "padded"), 0) | {"windows": len(recs)}
    bad, low = set(), []
    for it in store.load("lidc", "train", store_root):
        if it.item_id not in by_src:
            continue
        m = it.masks["nodule"]
        H, W = m.shape[-2:]
        ux0, uy0, ux1, uy1 = union_box(it)
        dark = _contrast(it) < LOW_CONTRAST
        for r in by_src.pop(it.item_id):
            x0, y0, x1, y1 = box = r["box"]
            v = r["view"]
            inside = m[:, max(y0, 0):max(y1, 0), max(x0, 0):max(x1, 0)].sum((1, 2))
            out["inside"] += bool((inside == m.sum((1, 2))).all())
            out["clear"] += min(ux0 - x0, uy0 - y0, x1 - ux1, y1 - uy1) >= CLEARANCE and not crowded(it, box)
            out["padded"] += x0 < 0 or y0 < 0 or x1 > W or y1 > H
            out["fallback"] += bool(v["fallback"])
            out["kept_side"] += bool(v["kept_side"])
            low += [r["item_id"]] if dark else []
            try:
                f = frame(it, offset=v["offset"], side=v["side"])
            except ValueError:
                bad.add(r["item_id"])
                continue
            png = cv2.imread(str(rows_root / "images" / f"{r['item_id']}.png"), cv2.IMREAD_COLOR)
            sha = pixels_sha256(f.image)
            if (list(f.box) != box or f.W != v["side"] or sha != r["image_sha256"] or png is None
                    or pixels_sha256(cv2.cvtColor(png, cv2.COLOR_BGR2RGB)) != sha):
                bad.add(r["item_id"])
    bad |= {r["item_id"] for rs in by_src.values() for r in rs}  # rows whose source item is not a train item
    return out | {"mismatch": sorted(bad), "low_contrast": sorted(set(low))}


def view_windows(view_dir, store_root) -> list[Window]:
    """The windows of a shifted eval view as Window objects (geometry only; DG16): kind "shifted", draw 0."""
    view_dir = Path(view_dir)
    split = read_json(view_dir / "manifest.json")["view"].split("-")[1]
    recs = {r["item_id"]: r for r in map(json.loads, (view_dir / "lidc" / f"{split}.jsonl").read_text().splitlines())}
    out = []
    for it in store.load("lidc", split, store_root):
        if r := recs.get(it.item_id):
            bare, ub = replace(it, image=None), union_box(it)
            s0, v = frame(bare).W, r["view"]
            out.append(Window(it.item_id, it.group_id, it.item_id.rsplit("/", 1)[0], "shifted", 0, *v["offset"],
                              v["side"], s0, offset_of(ub, r["box"]), max(ub[2] - ub[0], ub[3] - ub[1]),
                              r["disagreement"]))
    return out


def _law(windows) -> tuple[np.ndarray, list[float]]:
    return cell_shares(windows), [w.offset for w in windows]


def _counts(scheduled, pool, nb: int, shares, w: float) -> list[str]:
    """DG11: per block of nb rows, the kind counts, each off-centre kind's cell counts and each (kind, cell)
    disagreement count against the largest-remainder counts (windows.schedule steps 1-3)."""
    law = [x for x in pool.windows if x.kind == "shifted" and x.draw == 0]
    t = cell_shares(law).reshape(-1)
    n = Counter((x.kind, x.cell, x.disagreement) for x in pool.windows)
    bad = []
    for b in range(0, len(scheduled), nb):
        rows = scheduled[b:b + nb]
        kinds = Counter(x.kind for x in rows)
        want = windows._lrr(np.array([shares[k] for k in KINDS]) * nb, nb)
        if [kinds[k] for k in KINDS] != list(want):
            bad.append(f"block {b // nb}: kind counts {[kinds[k] for k in KINDS]}, want {list(want)}")
        for k in KINDS[1:]:
            got = Counter(x.cell for x in rows if x.kind == k)
            m = windows._lrr(t * kinds[k], kinds[k])
            for c, mc in zip(windows.CELLS, m, strict=True):
                if got[c] != mc:
                    bad.append(f"block {b // nb}: {k} cell {c} has {got[c]} rows, want {mc}")
                elif mc:
                    p = n[k, c, True] / max(n[k, c, True] + n[k, c, False], 1)
                    q = w * p / (1 + (w - 1) * p)
                    want_dis = int(windows._lrr(np.array([q, 1 - q]) * mc, int(mc))[0])
                    if (d := sum(x.disagreement for x in rows if x.kind == k and x.cell == c)) != want_dis:
                        bad.append(f"block {b // nb}: {k} cell {c} has {d} disagreement rows, want {want_dis}")
    return bad


def _null(scheduled, law, patients: int, seed: int = 0) -> dict:
    """Quantiles (5, 50, 95 %) of the schedule's coverage against the law of `patients` random train patients."""
    groups = {}
    for x in law:
        groups.setdefault(x.group_id, []).append(x)
    names, rng, got = sorted(groups), np.random.default_rng(seed), []
    for _ in range(NULL_DRAWS):
        sub = [x for g in rng.choice(names, min(patients, len(names)), replace=False) for x in groups[g]]
        got.append(cell_coverage(scheduled, *_law(sub)))
    return {k: [float(np.percentile([g[k] for g in got], q)) for q in (5, 50, 95)] for k in ("worst_ratio", "tv", "ks")}


def schedule_report(scheduled: Sequence[Window], pool, val: Mapping[str, Sequence[Window]], *, per_step: int,
                    block_steps: int, shares: Mapping[str, float], disagreement_weight: float) -> dict:
    """DG7-DG11, DG9b, DG14, DG16 of the scheduled windows (module docstring). The target is the train shifted law:
    the pool's draw-0 shifted windows. `val`: name -> windows (DG16, with its null from as many train patients)."""
    law = [x for x in pool.windows if x.kind == "shifted" and x.draw == 0]
    target, nb = _law(law), per_step * block_steps
    off = [x for x in scheduled if x.kind != "centred"]
    nodules = Counter(x.nodule for x in scheduled).most_common()
    dis = Counter(x.window_id for x in scheduled if x.disagreement)
    resc = {}
    for x in scheduled:
        if x.kind == "rescaled":
            resc.setdefault(f"{x.cell[0]}{x.cell[1]}", Counter())[x.item_id] += 1
    return {
        "prefixes": [cell_coverage(scheduled[:k], *target) for k in range(nb, len(scheduled) + nb, nb)],
        "mix": mix_coverage(list(scheduled), *target),
        "s0_share": sum(x.side == x.s0 for x in off) / len(off) if off else 0.0,
        "count_mismatch": _counts(list(scheduled), pool, nb, shares, disagreement_weight), "unknown": 0,
        "kinds": {k: sum(x.kind == k for x in scheduled) for k in KINDS},
        "disagreement_share": sum(dis.values()) / len(scheduled) if scheduled else 0.0,
        "max_disagreement_reuse": max(dis.values(), default=0),
        "top10_nodule_share": sum(n for _, n in nodules[:ceil(len(nodules) / 10)]) / max(len(scheduled), 1),
        "cells": {k: {f"{a}{b}": n for (a, b), n in sorted(Counter(x.cell for x in scheduled if x.kind == k).items())}
                  for k in KINDS},
        "dropped": dict(pool.dropped),
        "rescaled": {c: {"items": len(n), "max_reuse": max(n.values())} for c, n in sorted(resc.items())},
        "val": {name: {"patients": len({x.group_id for x in v}), "coverage": cell_coverage(scheduled, *_law(v)),
                       "null": _null(scheduled, law, len({x.group_id for x in v}))} for name, v in val.items()},
    }


def tokens(rows_root, tokenizer, limit: int) -> dict:
    """DG12: token counts of every certified target (no special tokens; the end token is the gate's + 1) and the
    count of target errors."""
    recs = records(rows_root)
    n = np.array([len(x) for x in tokenizer([t for r in recs for t in r["targets"] if t],
                                            add_special_tokens=False)["input_ids"]])
    return {"limit": limit, "targets": len(n), "max": int(n.max(initial=0)),
            "p99_9": float(np.percentile(n, 99.9)) if len(n) else 0.0, "over_1200": int((n > 1200).sum()),
            "over_limit_minus_100": int((n > limit - 100).sum()), "errors": sum(len(r["errors"]) for r in recs)}


def agreement(rows_root, low_contrast) -> dict:
    """DG14b (report): per RL row, the mean doctor leave-one-out IoU_tc (objdel.verify.ceiling.loo_caches, each left-out
    doctor scored with score_mask against the other doctors); per cell, overall, the share below 0.3 and that share
    jointly with low contrast (window ids in `low_contrast`)."""
    from objdel.verify.ceiling import loo_caches  # A2
    from objdel.verify.reward import score
    from objdel.verify.reward.cache import ExpertCache

    vals, low = {}, set(low_contrast)
    for r in records(rows_root):
        if r["rl"]:
            c = ExpertCache.load(Path(rows_root, "caches", r["item_id"], "nodule.npz"), r["caches"]["nodule"])
            v = float(np.mean([score.score_mask(c.raters[x.index], x.cache).iou_tc for x in loo_caches(c)]))
            vals[r["item_id"]] = (r["cell"], v)
    v = np.array([x for _, x in vals.values()])
    cells = {}
    for c, x in vals.values():
        cells.setdefault(c, []).append(x)
    return {"windows": len(v), "mean": float(v.mean()) if len(v) else 0.0,
            "share_lt_0.3": float((v < 0.3).mean()) if len(v) else 0.0,
            "share_lt_0.3_low_contrast": sum(x < 0.3 and i in low for i, (_, x) in vals.items()) / max(len(v), 1),
            "cells": {c: {"n": len(x), "mean": float(np.mean(x)), "share_lt_0.3": float(np.mean(np.array(x) < 0.3))}
                      for c, x in sorted(cells.items())}}


def view_content_sha256(view_dir) -> str:
    """B1.6 content hash: SHA-256 of the newline-joined JSON tuples [item_id, view.offset, view.side, image_sha256]
    of the view's LIDC records in file order (a prepared root without windows has null offset and side). Manifest
    versions and the reward fingerprint do not enter it."""
    recs = [json.loads(x) for f in sorted(Path(view_dir, "lidc").glob("*.jsonl")) for x in f.read_text().splitlines()]
    lines = [json.dumps([r["item_id"], (r.get("view") or {}).get("offset"), (r.get("view") or {}).get("side"),
                         r["image_sha256"]]) for r in recs]
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def target_set_sha256(view_dir) -> str:
    """The `target_set_sha256` of the one set listed in the view directory's `manifest.json`."""
    (s,) = read_json(Path(view_dir, "manifest.json"))["sets"].values()
    return s["target_set_sha256"]


def frozen(eval_views_root, entries: Mapping[str, dict]) -> dict[str, str]:
    """DG15: frozen view -> why it differs from its `frozen_views.json` entry (missing entry or view, content or
    target-set hash); {} passes."""
    bad = {}
    for v in FROZEN:
        d, e = Path(eval_views_root, v), entries.get(v)
        if e is None or not d.is_dir():
            bad[v] = "no frozen entry" if e is None else f"no view {d}"
        elif view_content_sha256(d) != e["content_sha256"]:
            bad[v] = "content_sha256 differs"
        elif target_set_sha256(d) != e["target_set_sha256"]:
            bad[v] = "target_set_sha256 differs"
    return bad


def _hard(name, value, limit, ok) -> Gate:
    return Gate(name, float(value), limit, True, bool(ok))


def _report(name, value, what) -> Gate:
    return Gate(name, float(value), f"report ({what})", False, True)


def gates(report: Mapping) -> list[Gate]:
    """The 21 data gates of B.6.1 in order DG1, DG1b, DG2..DG9, DG9b, DG10..DG14, DG14b, DG15..DG18. KeyError when a
    report part is missing (never a pass)."""
    lk, rw, wn, sc = report["leakage"], report["rows"], report["windows"], report["schedule"]
    tk, ag, fz, z = report["tokens"], report["agreement"], report["frozen"], report["z"]
    n, pre = wn["windows"], sc["prefixes"]
    share = lambda k: wn[k] / n if n else 0.0  # noqa: E731
    dg1 = lk["patients"] + lk["series"] + lk["items"]
    dg1b = sum(map(len, (rw["not_train"], rw["val_test_patient"], rw["view_image"], rw["excluded"])))
    dg1b += bool(rw["exclusions_stale"])
    worst = min((c["worst_ratio"] for c in pre), default=0.0)
    tv, ks = max((c["tv"] for c in pre), default=1.0), max((c["ks"] for c in pre), default=1.0)
    if not all(np.isfinite([c[k] for c in pre for k in ("worst_ratio", "tv", "ks")])):
        worst, tv, ks = 0.0, 1.0, 1.0  # min/max drop a NaN that is not first: a NaN prefix fails DG7-DG9
    val = sc["val"]["val"]["coverage"]
    return [
        _hard("DG1", dg1, "== 0", dg1 == 0),
        _hard("DG1b", dg1b, "== 0", dg1b == 0),
        _hard("DG2", lk["slice_sha"] + lk["crop_sha"], "== 0", lk["slice_sha"] + lk["crop_sha"] == 0),
        _hard("DG3", sum(not p["reviewed"] for p in lk["near"]), "== 0 unreviewed",
              all(p["reviewed"] for p in lk["near"])),
        _hard("DG4", share("inside"), "== 1", n and wn["inside"] == n),
        _hard("DG5", share("clear"), "== 1", n and wn["clear"] == n),
        _hard("DG6", wn["fallback"] + wn["kept_side"], "== 0", wn["fallback"] + wn["kept_side"] == 0),
        _hard("DG7", worst, ">= 0.8", worst >= 0.8),
        _hard("DG8", tv, "<= 0.10", tv <= 0.10),
        _hard("DG9", ks, "<= 0.06", ks <= 0.06),
        _report("DG9b", sc["mix"]["tv"], "full mix TV"),
        _hard("DG10", sc["s0_share"], ">= 0.75", sc["s0_share"] >= 0.75),
        _hard("DG11", len(sc["count_mismatch"]) + sc["unknown"], "== 0",
              not sc["count_mismatch"] and not sc["unknown"]),
        _hard("DG12", tk["max"] + 1, f"< {tk['limit'] - 100}; 0 target errors",
              tk["max"] + 1 < tk["limit"] - 100 and tk["errors"] == 0),
        _hard("DG13", len(wn["mismatch"]), "== 0", n and not wn["mismatch"]),
        _report("DG14", sc["disagreement_share"], "disagreement share"),
        _report("DG14b", ag["share_lt_0.3"], "doctor LOO IoU_tc < 0.3 share"),
        _hard("DG15", len(fz), "== 0", not fz),
        _report("DG16", val["tv"], "TV against the whole val"),
        _hard("DG17", len(z), "== 0", not z),
        _hard("DG18", wn["padded"], "== 0", wn["padded"] == 0),
    ]


def failures(gs: Sequence[Gate]) -> list[Gate]:
    """The hard gates that failed."""
    return [g for g in gs if g.hard and not g.ok]
