"""Group-aware, seeded, deterministic splits, near-duplicate grouping and leak checks (K11).

- LIDC: the unit is the patient (both scans of the 8 two-scan patients stay together), 70 / 10 / 20. Patients
  without reader-topology disagreement are dealt per stratum of multi-reader slice count; patients with disagreement
  are dealt by weighted greedy on their disagreement-slice count, so those (rare, concentrated) slices split about
  70 / 10 / 20 too.
- RIGA: the literature split (MRNet, TAB, MSE-Nets): test = Magrabia (95); train = BinRushed + MESSIDOR, with a 10 %
  validation carve-out stratified by subset. Near-duplicate photographs form one group; a group that touches the
  test set goes to test.
- Assignment: within a stratum, groups are ordered by SHA-256(seed, group id) and dealt by a largest-remainder
  schedule, so the result depends only on the seed and the population (no RNG state, no library version).
- Near-duplicates: a pair is a duplicate when either test fires.
  (1) Shift-tolerant normalised cross-correlation of vessel-texture descriptors (riga.descriptor), NCC >= 0.45. The
  lab measured same-eye pairs >= 0.49 and unrelated pairs <= 0.33 (cross-split maximum 0.21), while global
  perceptual hashes put thousands of unrelated fundus pairs within pHash distance 10.
  (2) Keypoint geometric verification (B3-8): ORB matches inside the field of view (riga.keypoints, 500 features)
  after Lowe's ratio test, then the RANSAC homography inlier count, >= 22. It tolerates the rotation, scale and
  shift of a reframed photograph of the same eye, which (1) does not (BinRushed4/image002 vs BinRushed3/image004:
  NCC 0.27, 154 inliers). Calibration over all 280,875 RIGA pairs: every pair with >= 8 inliers (512) was labelled
  independently by the correlation of high-pass vessel maps after a 1,500-feature similarity alignment. Same eye
  (correlation 0.30-0.95, 51 pairs, including all 37 NCC pairs): >= 29 inliers. Different eyes (correlation
  <= 0.07): <= 16 inliers.

References:
- RIGA literature split: Ji et al., MRNet, CVPR 2021, https://github.com/jiwei0921/MRNet ;
  TAB, https://arxiv.org/abs/2306.01340
- Patient-level LIDC splits "to reduce bias ... from the similarity between adjacent lung nodule slices":
  ProSeg, https://arxiv.org/abs/2512.00748 ; two-timepoint patients: https://www.cancerimagingarchive.net/collection/lidc-idri/
- Slice/image-level splits inflate CT segmentation scores: CTSCAN, https://arxiv.org/abs/2604.15561
- Lab evidence: archive/phase1_lab/t_data/riga/dedup.py, splits/build_splits.py (11 tests).
- Lowe, "Distinctive image features from scale-invariant keypoints", IJCV 2004 (ratio test 0.75 here):
  https://www.cs.ubc.ca/~lowe/papers/ijcv04.pdf ; spatial verification by RANSAC inliers for near-duplicate
  retrieval: Philbin et al., CVPR 2007, https://www.robots.ox.ac.uk/~vgg/publications/2007/Philbin07/philbin07.pdf
- OpenCV feature matching + findHomography (RANSAC): https://docs.opencv.org/4.x/d1/de0/tutorial_py_feature_homography.html
- Fundus pairs of one eye under rotation / scale / shift, registered from keypoint correspondences: FIRE
  (Hernandez-Matas et al., Modeling and Artificial Intelligence in Ophthalmology 1(4), 2017),
  https://doi.org/10.35119/maio.v1i4.42 ; https://projects.ics.forth.gr/cvrl/fire/
"""
import hashlib
from collections import Counter
from dataclasses import dataclass

import numpy as np

LIDC_FRACTIONS = (("train", 0.7), ("val", 0.1), ("test", 0.2))
RIGA_VAL_FRACTION = 0.1
NCC_DUP = 0.45
SHIFT = 12  # descriptor shift searched, px (about 48 px of a 2,000-px photograph)
ORB_DUP = 22  # RANSAC inliers: between different eyes (<= 16) and the same eye (>= 29), all RIGA pairs
RATIO = 0.75  # Lowe's ratio test


@dataclass(frozen=True)
class Patient:
    slices: int  # multi-reader slices
    disagree: int  # of which readers' topologies disagree


def _key(seed, gid) -> str:
    return hashlib.sha256(f"{seed}:{gid}".encode()).hexdigest()


def deal(groups, fractions, seed) -> dict:
    """Largest remainder: the k-th group (seeded hash order) goes to the split furthest below its target k * f."""
    counts, out = Counter(), {}
    for k, g in enumerate(sorted(groups, key=lambda g: _key(seed, g)), 1):
        name = max(fractions, key=lambda nf: (k * nf[1] - counts[nf[0]], nf[1]))[0]
        counts[name] += 1
        out[g] = name
    return out


def deal_weighted(groups, weight, fractions, seed) -> dict:
    """Heaviest group first (ties by seeded hash) to the split furthest below its share of the total weight."""
    total = sum(weight[g] for g in groups) or 1
    have, out = Counter(), {}
    for g in sorted(groups, key=lambda g: (-weight[g], _key(seed, g))):
        name = max(fractions, key=lambda nf: ((nf[1] * total - have[nf[0]]) / nf[1], nf[1]))[0]
        have[name] += weight[g]
        out[g] = name
    return out


def lidc_stratum(p: Patient) -> str:
    if p.disagree:
        return "dis=1"
    s = p.slices
    return "dis=0|slices=" + ("0" if s == 0 else "1-10" if s <= 10 else "11-40" if s <= 40 else "41+")


def lidc_split(patients: dict[str, Patient], seed=0) -> dict[str, str]:
    """Patient id -> split."""
    strata = {}
    for pid, p in patients.items():
        strata.setdefault(lidc_stratum(p), []).append(pid)
    out = {}
    for st in sorted(strata):
        if st == "dis=1":
            out.update(deal_weighted(strata[st], {p: patients[p].disagree for p in strata[st]}, LIDC_FRACTIONS, seed))
        else:
            out.update(deal(strata[st], LIDC_FRACTIONS, seed))
    return out


def riga_split(subsets: dict[str, str], groups: dict[str, str], seed=0) -> dict[str, str]:
    """Image id -> split, from each image's subset and near-duplicate group id."""
    test = {groups[i] for i, s in subsets.items() if s.startswith("Magrabia")}
    members = {}
    for i in sorted(subsets):
        if groups[i] not in test:
            members.setdefault(groups[i], []).append(subsets[i])
    strata = {}
    for g, subs in members.items():
        strata.setdefault(min(subs), []).append(g)
    fractions = (("train", 1 - RIGA_VAL_FRACTION), ("val", RIGA_VAL_FRACTION))
    assign = {}
    for st in sorted(strata):
        assign.update(deal(strata[st], fractions, seed))
    return {i: "test" if groups[i] in test else assign[groups[i]] for i in subsets}


def near_duplicates(descriptors, threshold=NCC_DUP) -> list[tuple[int, int, float]]:
    """(i, j, ncc) for i < j with max circular cross-correlation over shifts within +-SHIFT >= threshold."""
    f = np.fft.rfft2(descriptors)
    iu, ju = np.triu_indices(len(descriptors), 1)
    ncc = np.empty(iu.size)
    near = np.r_[0:SHIFT + 1, -SHIFT:0]
    for s in range(0, iu.size, 1000):  # 1,000 pairs per FFT batch bounds the memory
        a, b = iu[s:s + 1000], ju[s:s + 1000]
        cc = np.fft.irfft2(f[a] * np.conj(f[b]), s=descriptors.shape[1:])
        ncc[s:s + 1000] = cc[:, near][:, :, near].reshape(len(a), -1).max(1)
    hit = ncc >= threshold
    return [(int(i), int(j), float(c)) for i, j, c in zip(iu[hit], ju[hit], ncc[hit], strict=True)]


def _inliers(a, b) -> int:
    import cv2
    (pa, da), (pb, db) = a, b
    if da is None or db is None or len(da) < 2 or len(db) < 2:
        return 0
    good = [m[0] for m in cv2.BFMatcher(cv2.NORM_HAMMING).knnMatch(da, db, k=2)
            if len(m) == 2 and m[0].distance < RATIO * m[1].distance]
    if len(good) < 4:  # a homography needs 4 correspondences
        return 0
    src, dst = pa[[m.queryIdx for m in good]], pb[[m.trainIdx for m in good]]
    _, inl = cv2.findHomography(src, dst, cv2.RANSAC, 5.0)
    return 0 if inl is None else int(inl.sum())


def keypoint_duplicates(features, threshold=ORB_DUP, pairs=None) -> list[tuple[int, int, int]]:
    """(i, j, inliers) for the pairs (default: every i < j, in order) whose RANSAC homography inlier count over
    ratio-tested ORB matches of features[i] -> features[j] is >= threshold. features: riga.keypoints outputs."""
    pairs = pairs if pairs is not None else [(i, j) for i in range(len(features)) for j in range(i + 1, len(features))]
    out = [(int(i), int(j), _inliers(features[i], features[j])) for i, j in pairs]
    return [p for p in out if p[2] >= threshold]


def group_ids(ids, pairs) -> dict[str, str]:
    """Id -> group id: connected components of the duplicate pairs, named by their smallest member."""
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in pairs:
        ra, rb = find(a), find(b)
        parent[max(ra, rb)] = min(ra, rb)
    return {i: find(i) for i in ids}


def leaks(split_of: dict, group_of: dict) -> list:
    """Groups whose members fall in more than one split."""
    seen = {}
    for m, g in group_of.items():
        seen.setdefault(g, set()).add(split_of[m])
    return sorted(g for g, s in seen.items() if len(s) > 1)


def crossing(pairs, split_of) -> list:
    """Pairs whose two members are in different splits."""
    return [p for p in pairs if split_of[p[0]] != split_of[p[1]]]
