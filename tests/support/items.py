"""Synthetic store items."""
import numpy as np

from tests.support.shapes import annulus


def lidc_item(item_id, discs, *, other=(), group="LIDC-IDRI-0000", n=200, image=True):
    """A synthetic LIDC slice: one doctor mask per disc (cy, cx, r[, hole]), other outlined nodules as boxes."""
    from objdel.data.store import Item
    hu = np.random.default_rng(len(item_id)).integers(-1000, 400, (n, n), dtype=np.int16) if image else None
    return Item(item_id, "lidc", group, "train", hu, {"nodule": np.stack([annulus(n, n, *d) for d in discs])},
                [str(k + 1) for k in range(len(discs))], {"other_nodules": [list(b) for b in other]})


def speck_item():
    """A LIDC item of three raters whose third mask is a 4-pixel speck: below min_area 5, so it abstains from the
    topology vote."""
    from objdel.data.store import Item
    rng = np.random.default_rng(1)
    speck = np.zeros((100, 100), bool)
    speck[50:52, 50:52] = True
    return Item("lidc/1.2.9/n000/z001", "lidc", "LIDC-IDRI-0009", "train",
                rng.integers(-1000, 400, (100, 100), dtype=np.int16),
                {"nodule": np.stack([annulus(100, 100, 50, 50, 12), annulus(100, 100, 50, 50, 11), speck])},
                ["1", "2", "3"], {})
