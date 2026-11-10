"""Expert-mask fixtures: keys '<item>__<structure>' -> bool stacks (k raters, h, w)."""
from __future__ import annotations

from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parents[2]
FIXTURES = {"riga": LAB / "fixtures" / "riga_experts.npz", "lidc": LAB / "fixtures" / "lidc_experts.npz"}


def load(name: str) -> dict[str, np.ndarray]:
    with np.load(FIXTURES[name]) as z:
        return {k: z[k].astype(bool) for k in sorted(z.files)}


def load_all() -> list[tuple[str, str, np.ndarray]]:
    """[(dataset, key, stack)] for both datasets, deterministic order."""
    return [(ds, k, v) for ds in ("riga", "lidc") for k, v in load(ds).items()]


def subset(dataset: str, key: str) -> str:
    """'riga/disc', 'riga/cup' or 'lidc/nodule'."""
    return f"{dataset}/{key.split('__')[1]}"


LIDC_FULL = LAB / "fixtures" / "lidc_full.npz"
LIDC_FULL_META = LAB / "fixtures" / "lidc_full_meta.csv"


def lidc_full_meta() -> list[dict]:
    import csv

    with open(LIDC_FULL_META) as f:
        return list(csv.DictReader(f))


def load_lidc_full(split: str | None = "train", disagree_only: bool = False, min_raters: int = 2) -> dict[str, np.ndarray]:
    """Full LIDC build (data agent): {key: bool (k, h, w)}. Rows whose cluster has a duplicated
    reader are dropped. disagree_only keeps items whose raters' SIZE-FILTERED per-instance
    signatures differ (recomputed here with reward.signature_min_area, not the meta flag)."""
    from .reward import signature_min_area
    from .topology import filtered_signature

    meta = {r["key"]: r for r in lidc_full_meta()}
    out = {}
    with np.load(LIDC_FULL) as z:
        for k in sorted(z.files):
            r = meta.get(k)
            if r is None or r["dup_reader_in_cluster"] == "True" or (split and r["split"] != split):
                continue
            st = z[k].astype(bool)
            if len(st) < min_raters:
                continue
            if disagree_only:
                a = signature_min_area(int(st.any(0).sum()))
                if len({filtered_signature(m, a) for m in st}) < 2:
                    continue
            out[k] = st
    return out
