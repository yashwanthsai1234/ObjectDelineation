"""Near-duplicate / leakage audit across RIGA subsets.

Input: grey 576x384 thumbnails of each image's reference (unannotated) image,
cached by thumbs.py from riga_extract.py's refs/ (prime image, or the 6-copy
reconstruction for BinRushed1). Signals per pair:
  - pHash / dHash / wHash (imagehash, 64 bit) Hamming distance;
  - NCC: shift-tolerant (+-12 px at 1/8 scale) normalised cross-correlation of
    a high-pass central crop. Global hashes and plain thumbnail correlation are
    dominated by the shared fundus layout / FOV border (thousands of unrelated
    pairs within pHash <= 10); the high-pass central crop measures vessel
    texture, which identifies the eye.
Near-duplicate := NCC >= NCC_DUP (0.45, set in the empty gap between same-eye
pairs, >= 0.49, and unrelated pairs, <= 0.33; checked by eye).

    python -I dedup.py <thumbs.npz> <out_dir>
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import imagehash
import numpy as np
from PIL import Image
from scipy import ndimage

NCC_DUP = 0.45
SHIFT = 12


def group_of(i):
    return "test" if i.startswith("Magrabia") else "train"


def desc(g):
    c = g[96:288, 144:432].astype(np.float32).reshape(96, 2, 144, 2).mean((1, 3))
    hp = c - ndimage.gaussian_filter(c, 2)
    hp -= hp.mean()
    return hp / (np.linalg.norm(hp) + 1e-6)


def pair_ncc(D, iu, chunk=4000):
    F = np.fft.rfft2(D)
    out = np.empty(len(iu[0]), np.float32)
    S = SHIFT
    for s in range(0, len(out), chunk):
        a, b = iu[0][s:s + chunk], iu[1][s:s + chunk]
        cc = np.fft.irfft2(F[a] * np.conj(F[b]), s=D.shape[1:])
        w = np.concatenate([cc[:, :S + 1], cc[:, -S:]], 1)
        w = np.concatenate([w[:, :, :S + 1], w[:, :, -S:]], 2)
        out[s:s + chunk] = w.reshape(len(a), -1).max(1)
    return out


def main(thumbs, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    z = np.load(thumbs)
    T, ids = z["thumbs"], [str(i).replace("_image", "/image") for i in z["ids"]]
    n = len(ids)
    split = np.array([group_of(i) for i in ids])
    iu = np.triu_indices(n, 1)
    cross = split[iu[0]] != split[iu[1]]
    ncc = pair_ncc(np.stack([desc(g) for g in T]), iu)
    H = {}
    for name, fn in (("phash", imagehash.phash), ("dhash", imagehash.dhash), ("whash", imagehash.whash)):
        h = np.stack([fn(Image.fromarray(g)).hash.ravel() for g in T]).astype(np.uint8)
        H[name] = np.concatenate([(h[iu[0][s:s + 20000]] != h[iu[1][s:s + 20000]]).sum(1) for s in range(0, len(iu[0]), 20000)])
    flat = T.reshape(n, -1)
    exact = defaultdict(list)
    for i in range(n):
        exact[flat[i].tobytes().__hash__()].append(ids[i])
    bins = [-1, .1, .2, .3, .45, .6, .8, .95, .99, 1.01]
    s = {"images": n, "by_split": {k: int((split == k).sum()) for k in ("train", "test")}, "pairs": int(len(iu[0])),
         "cross_split_pairs": int(cross.sum()), "ncc_threshold": NCC_DUP,
         "ncc": {"max_cross_split": float(ncc[cross].max()), "bins": bins,
                 "hist_cross_split": np.histogram(ncc[cross], bins)[0].tolist(),
                 "hist_within_split": np.histogram(ncc[~cross], bins)[0].tolist()},
         "identical_thumbnail_groups": sum(len(v) > 1 for v in exact.values())}
    for k, d in H.items():
        s[k] = {"min_cross_split": int(d[cross].min()), "pairs_le10_cross": int((d[cross] <= 10).sum()),
                "pairs_le10_within": int((d[~cross] <= 10).sum()),
                "max_over_ncc_dups": int(d[ncc >= NCC_DUP].max())}
    dup = ncc >= NCC_DUP
    pairs = [{"a": ids[a], "b": ids[b], "ncc": round(float(c), 4), "phash": int(H["phash"][k]), "dhash": int(H["dhash"][k]),
              "kind": "same_photo" if c >= 0.95 else "same_eye", "cross_split": bool(split[a] != split[b]),
              "cross_subset": ids[a].split("/image")[0] != ids[b].split("/image")[0]}
             for k, (a, b, c) in enumerate(zip(iu[0], iu[1], ncc)) if dup[k]]
    pairs.sort(key=lambda p: -p["ncc"])
    border = [{"a": ids[a], "b": ids[b], "ncc": round(float(c), 4)} for a, b, c in zip(iu[0], iu[1], ncc) if 0.25 <= c < NCC_DUP]
    s["near_duplicate_pairs"] = len(pairs)
    s["near_duplicate_pairs_by_kind"] = {k: sum(p["kind"] == k for p in pairs) for k in ("same_photo", "same_eye")}
    s["near_duplicate_pairs_cross_split"] = sum(p["cross_split"] for p in pairs)
    s["near_duplicate_pairs_cross_subset"] = sum(p["cross_subset"] for p in pairs)
    s["borderline_0.25_0.45"] = sorted(border, key=lambda p: -p["ncc"])
    o = np.argsort(-np.where(cross, ncc, -9))[:5]
    s["closest_cross_split"] = [{"a": ids[iu[0][k]], "b": ids[iu[1][k]], "ncc": round(float(ncc[k]), 4), "phash": int(H["phash"][k])} for k in o]
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for p in pairs:
        ra, rb = find(p["a"]), find(p["b"])
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    groups = defaultdict(list)
    for i in ids:
        groups[find(i)].append(i)
    dup_groups = {k: sorted(v) for k, v in sorted(groups.items()) if len(v) > 1}
    s["near_duplicate_groups"] = len(dup_groups)
    s["images_in_dup_groups"] = sum(len(v) for v in dup_groups.values())
    (out / "pairs.json").write_text(json.dumps(pairs, indent=1))
    (out / "dup_groups.json").write_text(json.dumps(dup_groups, indent=1, sort_keys=True))
    (out / "summary.json").write_text(json.dumps(s, indent=1))
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:3])
