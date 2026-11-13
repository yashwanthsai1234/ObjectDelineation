"""A tiny prepared view root (one item, one object) and answers for it, built in tmp_path."""

import numpy as np
import pytest

from objdel.data.targets.build import build_target
from objdel.verify.reward.cache import ExpertCache
from tests.support.prepared import write_prepared_root
from tests.support.shapes import disc as disc_hw

SIDE = 50
KEY = "lidc/1.2.0/n000/z000"


def disk(cy, cx, r):
    return disc_hw(SIDE, SIDE, cy, cx, r)


def build_root(root, items):
    """A prepared view root at `root`: {item key: list of rater masks} -> one lidc/val record and cache per key."""
    return write_prepared_root(root, {"lidc/val": [{"item_id": key, "rl": True, "caches": {"nodule": np.stack(masks)}}
                                                   for key, masks in items.items()]})


@pytest.fixture
def view_root(tmp_path, clean_store):
    return build_root(tmp_path, {KEY: [disk(25, 25, 10), disk(25, 26, 11), disk(24, 25, 9)]})


@pytest.fixture
def answers(view_root):
    cache = ExpertCache.load(view_root / "caches" / KEY / "nodule.npz")
    return [build_target({"nodule": m}, ("nodule",), min_area={"nodule": cache.min_area}).text
            for m in cache.raters]
