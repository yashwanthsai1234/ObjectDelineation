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
