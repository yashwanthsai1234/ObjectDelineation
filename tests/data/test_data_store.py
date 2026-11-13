"""Parquet store: exact round trip, deterministic bytes (K10), manifest, contract checks (tampered files:
tests/data/test_breakit_data.py::test_store_tampering_is_detected)."""
import json

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

from objdel.core.hashing import file_sha256
from objdel.data import store
from objdel.data.store import Item


def riga_item(i, split="train", h=37, w=45):
    rng = np.random.default_rng(i)
    return Item(item_id=f"riga/MESSIDOR/image{i:03d}", dataset="riga", group_id=f"riga:MESSIDOR/image{i:03d}",
                split=split, image=rng.integers(0, 256, (h, w, 3), dtype=np.uint8),
                masks={"disc": rng.random((6, h, w)) < 0.4, "cup": rng.random((6, h, w)) < 0.2},
                raters=[str(e) for e in range(1, 7)], meta={"subset": "MESSIDOR", "index": i, "stages": {"1": "A"}})


def lidc_item(i, split="val", image=True, k=3):
    rng = np.random.default_rng(100 + i)
    return Item(item_id=f"lidc/1.2.3/n000/z{i:03d}", dataset="lidc", group_id="LIDC-IDRI-0001", split=split,
                image=rng.integers(-3024, 3000, (64, 64), dtype=np.int16) if image else None,
                masks={"nodule": rng.random((k, 64, 64)) < 0.1}, raters=[str(r) for r in range(1, k + 1)],
                meta={"series": "1.2.3", "slice": i})


def items():
    return [riga_item(1), riga_item(2), riga_item(3, "test"), lidc_item(4), lidc_item(5, image=False), lidc_item(6)]


def same(a, b):
    assert (a.item_id, a.dataset, a.group_id, a.split, a.raters, a.meta) == (
        b.item_id, b.dataset, b.group_id, b.split, b.raters, b.meta)
    assert (a.image is None and b.image is None) or (
        a.image.dtype == b.image.dtype and np.array_equal(a.image, b.image))
    assert a.masks.keys() == b.masks.keys()
    for k in a.masks:
        assert b.masks[k].dtype == bool and np.array_equal(a.masks[k], b.masks[k])


def test_round_trip_exact(tmp_path):
    src = items()
    store.write(src, tmp_path)
    got = [*store.load("riga", "train", tmp_path), *store.load("riga", "test", tmp_path),
           *store.load("lidc", "val", tmp_path)]
    assert len(got) == len(src)
    for a, b in zip(src, got, strict=True):
        same(a, b)


@settings(max_examples=40, deadline=None)
@given(arrays(bool, st.tuples(st.integers(1, 4), st.integers(1, 20), st.integers(1, 20))))
def test_mask_codec_round_trip(m):
    assert np.array_equal(store.decode_masks(store.encode_masks(m), m.shape[2]), m)


def test_same_items_same_bytes(tmp_path):  # K10
    a, b = store.write(items(), tmp_path / "a"), store.write(items(), tmp_path / "b")
    assert a == b
    for rel in a["files"]:
        assert (tmp_path / "a" / rel).read_bytes() == (tmp_path / "b" / rel).read_bytes()


def test_manifest_lists_hashes_and_rows(tmp_path):
    m = store.write(items(), tmp_path)
    assert json.loads((tmp_path / "manifest.json").read_text()) == m
    assert {k: v["rows"] for k, v in m["files"].items()} == {
        "lidc/val.parquet": 3, "riga/test.parquet": 1, "riga/train.parquet": 2}
    for rel, f in m["files"].items():
        assert file_sha256(tmp_path / rel) == f["sha256"]
    assert store.verify(tmp_path) == m


def test_missing_split_file_is_an_error(tmp_path):
    store.write(items(), tmp_path)
    with pytest.raises(KeyError):
        next(store.load("lidc", "test", tmp_path))


@pytest.mark.parametrize("bad", [
    lambda it: setattr(it, "raters", it.raters[:-1]),  # k != number of raters
    lambda it: setattr(it, "split", "holdout"),
    lambda it: setattr(it, "dataset", "other"),
    lambda it: it.masks.update(disc=it.masks["disc"][:, :-1]),  # mask frame != image frame
    lambda it: it.masks.update(disc=it.masks["disc"].astype(np.uint8)),
])
def test_write_rejects_contract_violations(tmp_path, bad):
    it = riga_item(1)
    bad(it)
    with pytest.raises(ValueError):
        store.write([it], tmp_path)


def test_write_requires_sorted_unique_ids(tmp_path):
    with pytest.raises(ValueError, match="order"):
        store.write([riga_item(2), riga_item(1)], tmp_path / "a")
    with pytest.raises(ValueError, match="order"):
        store.write([riga_item(1), riga_item(1)], tmp_path / "b")


def test_write_rejects_one_item_id_in_two_splits(tmp_path):  # B3-6: ids are unique across the whole store
    with pytest.raises(ValueError, match="two splits"):
        store.write([riga_item(1), riga_item(1, "test")], tmp_path)
    with pytest.raises(ValueError, match="two splits"):
        store.write([riga_item(1, "val"), riga_item(2), riga_item(1, "train")], tmp_path)


@pytest.mark.parametrize("dataset, image", [
    ("lidc", np.zeros((64, 64), np.uint16)),  # HU + 1024 as raw DICOM stores it
    ("lidc", np.zeros((64, 64), np.uint8)),  # already windowed
    ("lidc", np.zeros((64, 64), np.float32)),
    ("lidc", np.zeros((64, 64, 3), np.int16)),
    ("riga", np.zeros((37, 45), np.uint8)),  # grey
    ("riga", np.zeros((37, 45, 4), np.uint8)),  # RGBA
    ("riga", np.zeros((37, 45, 3), np.int16)),
])
def test_write_enforces_image_dtype_and_layout(tmp_path, dataset, image):  # B3-7
    it = lidc_item(1) if dataset == "lidc" else riga_item(1)
    it.image = image
    with pytest.raises(ValueError, match="image"):
        store.write([it], tmp_path)


def test_load_enforces_the_item_contract(tmp_path, monkeypatch):  # B3-7: a file written around the check is refused
    it = lidc_item(1)
    it.image = it.image.astype(np.uint16)
    monkeypatch.setattr(store, "_check", lambda it: it.masks["nodule"].shape[1:])
    store.write([it], tmp_path)
    monkeypatch.undo()
    with pytest.raises(ValueError, match="image"):
        next(store.load("lidc", "val", tmp_path))


@pytest.mark.parametrize("src, dst", [("riga/train.parquet", "riga/val.parquet"),  # wrong split
                                      ("lidc/val.parquet", "riga/val.parquet")])  # wrong dataset
def test_load_refuses_a_row_of_another_file(tmp_path, src, dst):  # review #12: split leakage is K11
    store.write(items(), tmp_path)
    (tmp_path / src).rename(tmp_path / dst)
    m = json.loads((tmp_path / "manifest.json").read_text())
    m["files"][dst] = m["files"].pop(src)  # a manifest that agrees with the moved file
    (tmp_path / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(ValueError, match=f"{dst}: row .* belongs to"):
        next(store.load(*dst.removesuffix(".parquet").split("/"), tmp_path))


def test_load_streams_one_row_at_a_time(tmp_path):  # I20: peak memory is about one row, not one batch of 16
    import tracemalloc
    rng = np.random.default_rng(7)
    src = [riga_item(i, h=600, w=800) for i in range(12)]
    for it in src:
        it.image = rng.integers(0, 256, it.image.shape, dtype=np.uint8)  # incompressible, 1.44 MB
    store.write(src, tmp_path)
    tracemalloc.start()
    try:
        for _ in store.load("riga", "train", tmp_path):
            pass
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    row = 600 * 800 * 3 + 2 * 6 * 600 * 800  # image + both decoded bool mask stacks
    assert peak < 4 * row  # 16-row batches held all 12 encoded rows at once


def test_small_row_groups_round_trip(tmp_path, monkeypatch):
    monkeypatch.setattr(store, "ROW_GROUP_BYTES", 1)  # one row per row group
    src = [lidc_item(i) for i in range(5)]
    store.write(src, tmp_path)
    import pyarrow.parquet as pq
    assert pq.ParquetFile(tmp_path / "lidc/val.parquet").metadata.num_row_groups == 5
    for a, b in zip(src, store.load("lidc", "val", tmp_path), strict=True):
        same(a, b)


# ---- data-backed -----------------------------------------------------------------------------------------------

@pytest.mark.data
def test_store_verifies(store_dir):
    m = store.verify(store_dir)
    assert {"riga/train.parquet", "riga/val.parquet", "riga/test.parquet", "lidc/train.parquet", "lidc/val.parquet",
            "lidc/test.parquet"} == set(m["files"])


@pytest.mark.data
def test_store_rewrite_gives_identical_bytes(store_dir, tmp_path):  # K10 on real data: load -> write is the identity
    m = store.write(store.load("lidc", "val", store_dir), tmp_path)
    assert m["files"]["lidc/val.parquet"] == store.verify(store_dir)["files"]["lidc/val.parquet"]
