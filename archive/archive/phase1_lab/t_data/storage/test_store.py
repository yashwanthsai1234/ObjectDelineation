"""Store checks on the built prototype (run build_store.py first): manifest SHA-256, bit-exact masks and images."""
import hashlib
import json
from pathlib import Path

import pytest

from build_store import verify

STORE = Path(__file__).parent / "store"


@pytest.mark.parametrize("codec", ["poly", "rle"])
def test_manifest_sha(codec):
    m = json.loads((STORE / codec / "manifest.json").read_text())
    data = (STORE / codec / "data.parquet").read_bytes()
    assert hashlib.sha256(data).hexdigest() == m["files"]["data.parquet"]


@pytest.mark.parametrize("codec", ["poly", "rle"])
def test_roundtrip_bit_exact(codec):
    v = verify(STORE / codec / "data.parquet", codec)
    assert v["masks"] > 0 and v["mask_mismatch"] == 0 and v["image_mismatch"] == 0


def test_no_group_in_two_splits():
    from datasets import Dataset
    ds = Dataset.from_parquet(str(STORE / "poly" / "data.parquet"), columns=["group_id", "split"])
    seen = {}
    for g, s in zip(ds["group_id"], ds["split"]):
        assert seen.setdefault(g, s) == s
