"""Group-aware, seeded, deterministic splits and leak checks.

- LIDC: the unit is the patient (both scans of the 8 two-scan patients stay together), 70 / 10 / 20. Patients
  without reader-topology disagreement are dealt per stratum of multi-reader slice count; patients with disagreement
  are dealt by weighted greedy on their disagreement-slice count, so those (rare, concentrated) slices split about
  70 / 10 / 20 too.
- Val halves: the val patients are dealt into "shape" (tunes reward v5) and "select"
  (judges it) per lidc_stratum, with the patient's item and disagreement-item counts in the frozen val view as its
  Patient. Seeds "val_half:0", "val_half:1", ... are tried; the first that balances items (within 4 %) and
  disagreement items (within 4) is kept. The table is stored in reference/val_halves.json; val_half reads it.
- Assignment: within a stratum, groups are ordered by SHA-256(seed, group id) and dealt by a largest-remainder
  schedule, so the result depends only on the seed and the population (no RNG state, no library version).

References:
- Patient-level LIDC splits "to reduce bias ... from the similarity between adjacent lung nodule slices":
  ProSeg, https://arxiv.org/abs/2512.00748 ; two-timepoint patients: https://www.cancerimagingarchive.net/collection/lidc-idri/
- Slice/image-level splits inflate CT segmentation scores: CTSCAN, https://arxiv.org/abs/2604.15561
"""
import functools
import hashlib
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from objdel.core.io import read_json

LIDC_FRACTIONS = (("train", 0.7), ("val", 0.1), ("test", 0.2))
VAL_HALVES = (("shape", 0.5), ("select", 0.5))
VAL_HALF_RULE = ("deal per lidc_stratum (dis=1 included), seed 'val_half:<s>', first s with |items(shape) - "
                 "items(select)| <= 4 % of items and |disagreement items(shape) - (select)| <= 4")


@dataclass(frozen=True)
class Patient:
    """Per-patient counts for the LIDC split strata: multi-reader slices and how many of them have disagreeing
    topologies.
    """
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
    """The stratum of a patient: `dis=1` when any slice disagrees, else a bin of the multi-reader slice count."""
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


def val_half_table(records, max_s=1000) -> dict:
    """The val halves from the val view's records (dicts with group_id and disagreement; nothing else is
    read): the first seed s < max_s that balances items and disagreement items (module docstring)."""
    if not records:
        raise ValueError("val_half_table: no records")
    items = Counter(r["group_id"] for r in records)
    dis = Counter(r["group_id"] for r in records if r["disagreement"])
    strata = {}
    for g in sorted(items):
        strata.setdefault(lidc_stratum(Patient(items[g], dis[g])), []).append(g)
    for s in range(max_s):
        halves = {}
        for st in sorted(strata):
            halves |= deal(strata[st], VAL_HALVES, seed=f"val_half:{s}")
        n = {h: sum(items[g] for g in halves if halves[g] == h) for h, _ in VAL_HALVES}
        d = {h: sum(dis[g] for g in halves if halves[g] == h) for h, _ in VAL_HALVES}
        if abs(n["shape"] - n["select"]) <= 0.04 * len(records) and abs(d["shape"] - d["select"]) <= 4:
            return {"rule": VAL_HALF_RULE, "patient_form": "items,disagreement_items", "s": s,
                    "halves": dict(sorted(halves.items())), "items": n, "disagreement_items": d}
    raise ValueError(f"no seed s < {max_s} balances the val halves")


@functools.cache
def _val_halves() -> dict:
    return read_json(Path(__file__).resolve().parents[3] / "reference/val_halves.json")["halves"]


def val_half(group_id: str) -> str:
    """"shape" or "select" for a val patient id (no prefix: "LIDC-IDRI-0026"); KeyError for any other id."""
    return _val_halves()[group_id]


def leaks(split_of: dict, group_of: dict) -> list:
    """Groups whose members fall in more than one split."""
    seen = {}
    for m, g in group_of.items():
        seen.setdefault(g, set()).add(split_of[m])
    return sorted(g for g, s in seen.items() if len(s) > 1)

