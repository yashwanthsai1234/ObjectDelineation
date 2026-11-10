"""A tiny prepared view root (one item, one object) and answers for it, built in tmp_path."""
import hashlib
import json

import numpy as np
import pytest

from objdel.data.targets.build import build_target
from objdel.verify.reward.cache import ExpertCache

SIDE = 50
KEY = "lidc/1.2.0/n000/z000"


def disk(cy, cx, r):
    yy, xx = np.ogrid[:SIDE, :SIDE]
    return (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r


def build_root(root, items):
    """A prepared view root at `root`: {item key: list of rater masks} -> one lidc/val record and cache per key."""
    lines = []
    for key, masks in items.items():
        cache = ExpertCache.build(np.stack(masks))
        (root / "caches" / key).mkdir(parents=True)
        lines.append(json.dumps({"item_id": key, "rl": True,
                                 "caches": {"nodule": cache.save(root / "caches" / key / "nodule.npz")}}))
    records = root / "lidc" / "val.jsonl"
    records.parent.mkdir(parents=True)
    records.write_text("\n".join(lines) + "\n")
    (root / "manifest.json").write_text(json.dumps({"files": {"lidc/val.jsonl": {
        "sha256": hashlib.sha256(records.read_bytes()).hexdigest()}}}))
    return root


@pytest.fixture
def view_root(tmp_path, clean_store):
    return build_root(tmp_path, {KEY: [disk(25, 25, 10), disk(25, 26, 11), disk(24, 25, 9)]})


@pytest.fixture
def answers(view_root):
    cache = ExpertCache.load(view_root / "caches" / KEY / "nodule.npz")
    return [build_target({"nodule": m}, ("nodule",), min_area={"nodule": cache.min_area}).text
            for m in cache.raters]
