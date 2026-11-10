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


@pytest.fixture
def view_root(tmp_path, clean_store):
    cache = ExpertCache.build(np.stack([disk(25, 25, 10), disk(25, 26, 11), disk(24, 25, 9)]))
    (tmp_path / "caches" / KEY).mkdir(parents=True)
    digest = cache.save(tmp_path / "caches" / KEY / "nodule.npz")
    records = tmp_path / "lidc" / "val.jsonl"
    records.parent.mkdir()
    records.write_text(json.dumps({"item_id": KEY, "rl": True, "caches": {"nodule": digest}}) + "\n")
    (tmp_path / "manifest.json").write_text(json.dumps({"files": {"lidc/val.jsonl": {
        "sha256": hashlib.sha256(records.read_bytes()).hexdigest()}}}))
    return tmp_path


@pytest.fixture
def answers(view_root):
    cache = ExpertCache.load(view_root / "caches" / KEY / "nodule.npz")
    return [build_target({"nodule": m}, ("nodule",), min_area={"nodule": cache.min_area}).text
            for m in cache.raters]
