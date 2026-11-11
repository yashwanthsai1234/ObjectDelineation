"""The item store: one Parquet file per (dataset, split), exact images and masks, and a SHA-256 manifest.

Layout: ``<root>/<dataset>/<split>.parquet`` plus ``<root>/manifest.json`` (SHA-256, rows and bytes per file, and
the library versions that wrote it). Rows are sorted by ``item_id`` and written with fixed settings, so the same
items give the same bytes (K10).

Encodings (chosen in the M4 sandbox, experiments/phase1/impl/M4_data.md):
- image: the decoded pixel array as NPY bytes (uint8 RGB for RIGA, int16 Hounsfield units for LIDC), compressed by
  Parquet zstd. Lossless for any dtype, decoder-independent (no JPEG/PNG library in the loop), 25x faster to decode
  than PNG; about 1.3x larger than PNG on RIGA photographs.
- masks: per structure, the bool stack (k, H, W) bit-packed along W (``np.packbits``) as NPY bytes. Exact, 8x
  smaller in memory than bool, and after zstd within a few KB of run-length encoding, which needs custom code.
- row groups close at about ROW_GROUP_BYTES of uncompressed payload; `load` decodes one row at a time without
  pre-buffering the row group (I20: a RIGA val pass peaks at 455 MB RSS instead of 1,543 MB with 16-row batches,
  at the same speed).

Contract checked on write and on load: item ids unique across the whole store (one id in two splits is a leak, B3-6);
every row's dataset and split are those of its file (on load); images uint8 (H, W, 3) RGB for RIGA and int16 (H, W)
Hounsfield units for LIDC, or None (B3-7: a uint16 HU + 1024 or a pre-windowed uint8 slice would silently get the
wrong lung window downstream).

References:
- Apache Parquet configuration (row group and page sizes): https://parquet.apache.org/docs/file-format/configurations/
- pyarrow.parquet.ParquetWriter / write_table options:
  https://arrow.apache.org/docs/python/generated/pyarrow.parquet.write_table.html
- Hugging Face datasets also closes Parquet row groups by uncompressed size (``config.MAX_ROW_GROUP_SIZE = "100MB"``
  in 5.x; 2.x used 100 rows for image data): https://github.com/huggingface/datasets/blob/main/src/datasets/config.py
- NPY format (self-describing dtype and shape, no pickle): https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html
- np.packbits: https://numpy.org/doc/stable/reference/generated/numpy.packbits.html
- pyarrow.parquet.ParquetFile (pre_buffer reads whole column chunks of a row group ahead) and iter_batches:
  https://arrow.apache.org/docs/python/generated/pyarrow.parquet.ParquetFile.html
- COCO RLE, the alternative considered: https://github.com/cocodataset/cocoapi/blob/master/PythonAPI/pycocotools/mask.py
"""
import functools
import io
import json
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from objdel.core.hashing import file_sha256
from objdel.core.io import read_json, write_json

DATASETS = ("riga", "lidc")
SPLITS = ("train", "val", "test")
ROW_GROUP_BYTES = 128 << 20
BATCH_ROWS = 1  # I20: rows are up to 16 MB (RIGA); one at a time bounds the reader's memory
IMAGE = {"riga": (np.dtype(np.uint8), 3), "lidc": (np.dtype(np.int16), 2)}  # dtype, ndim (B3-7)
SCHEMA = pa.schema([
    ("item_id", pa.string()), ("dataset", pa.string()), ("group_id", pa.string()), ("split", pa.string()),
    ("raters", pa.list_(pa.string())), ("height", pa.int32()), ("width", pa.int32()), ("image", pa.binary()),
    ("masks", pa.list_(pa.struct([("name", pa.string()), ("bits", pa.binary())]))), ("meta", pa.string()),
])
WRITE_OPTIONS = {"compression": "zstd", "compression_level": 6, "use_dictionary": False, "write_statistics": False}


@dataclass
class Item:
    """One prompt's data: an image and every rater's exact mask per structure.

    image: uint8 (H, W, 3) RGB for RIGA, int16 (H, W) Hounsfield units for LIDC, or None when there is no model
    image (LIDC CT pixels not available locally; the BinRushed1 reconstruction, B3-1). masks: structure name -> bool
    (k, H, W), one layer per rater, in the image frame. raters: k rater ids, in layer order. meta: JSON-serializable
    provenance (see objdel.data.riga / .lidc builders).
    """

    item_id: str
    dataset: str
    group_id: str
    split: str
    image: np.ndarray | None
    masks: dict[str, np.ndarray]
    raters: list[str]
    meta: dict


def _npy(a: np.ndarray) -> bytes:
    buf = io.BytesIO()
    np.save(buf, a, allow_pickle=False)
    return buf.getvalue()


def _from_npy(b: bytes) -> np.ndarray:
    return np.load(io.BytesIO(b), allow_pickle=False)


def encode_masks(stack: np.ndarray) -> bytes:
    """Pack a bool mask stack along its last axis into `.npy` bytes (`np.packbits`)."""
    return _npy(np.packbits(stack, axis=-1))


def decode_masks(bits: bytes, width: int) -> np.ndarray:
    """Unpack `encode_masks` bytes into a bool stack whose last axis has `width` entries."""
    return np.unpackbits(_from_npy(bits), axis=-1, count=width).view(bool)


def _check(it: Item) -> tuple[int, int]:
    """The (H, W) frame of a valid item; ValueError on any contract violation."""
    if it.dataset not in DATASETS or it.split not in SPLITS:
        raise ValueError(f"{it.item_id}: dataset {it.dataset!r} / split {it.split!r}")
    if not it.masks or {m.dtype for m in it.masks.values()} != {np.dtype(bool)}:
        raise ValueError(f"{it.item_id}: masks must be non-empty bool stacks")
    shapes = {m.shape for m in it.masks.values()}
    if len(shapes) != 1 or len(shape := shapes.pop()) != 3 or shape[0] != len(it.raters):
        raise ValueError(f"{it.item_id}: masks must share one (k, H, W) shape with k = len(raters)")
    if it.image is not None:
        dtype, ndim = IMAGE[it.dataset]
        if it.image.dtype != dtype or it.image.ndim != ndim or it.image.shape[:2] != shape[1:] or (
                ndim == 3 and it.image.shape[2] != 3):
            raise ValueError(f"{it.item_id}: image {it.image.dtype} {it.image.shape}, expected {dtype} "
                             f"{'(H, W, 3)' if ndim == 3 else '(H, W)'} on the mask frame {shape[1:]}")
    return shape[1:]


def _row(it: Item) -> dict:
    h, w = _check(it)
    if json.loads(json.dumps(it.meta)) != it.meta:
        raise ValueError(f"{it.item_id}: meta must be JSON-native")
    return {"item_id": it.item_id, "dataset": it.dataset, "group_id": it.group_id, "split": it.split,
            "raters": list(it.raters), "height": h, "width": w,
            "image": None if it.image is None else _npy(np.ascontiguousarray(it.image)),
            "masks": [{"name": k, "bits": encode_masks(v)} for k, v in it.masks.items()],
            "meta": json.dumps(it.meta, sort_keys=True)}


def write(items: Iterable[Item], root) -> dict:
    """Write items (ascending item_id within each dataset/split) and the manifest; returns the manifest."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    (root / "manifest.json").unlink(missing_ok=True)
    writers, buffers, last, rows, split_of = {}, {}, {}, {}, {}

    def flush(rel):
        rs, _ = buffers.pop(rel)
        if rel not in writers:
            (root / rel).parent.mkdir(parents=True, exist_ok=True)
            writers[rel] = pq.ParquetWriter(root / rel, SCHEMA, **WRITE_OPTIONS)
        writers[rel].write_table(pa.Table.from_pylist(rs, schema=SCHEMA))

    try:
        for it in items:
            rel = f"{it.dataset}/{it.split}.parquet"
            if split_of.setdefault(it.item_id, it.split) != it.split:  # B3-6
                raise ValueError(f"{it.item_id}: in two splits ({split_of[it.item_id]}, {it.split})")
            if rel in last and it.item_id <= last[rel]:
                raise ValueError(f"{rel}: item ids out of order at {it.item_id}")
            last[rel] = it.item_id
            row = _row(it)
            rs, size = buffers.get(rel, ([], 0))
            rs.append(row)
            buffers[rel] = (rs, size + sum(len(m["bits"]) for m in row["masks"]) + len(row["image"] or b""))
            rows[rel] = rows.get(rel, 0) + 1
            if buffers[rel][1] >= ROW_GROUP_BYTES:
                flush(rel)
        for rel in list(buffers):
            flush(rel)
    finally:
        for w in writers.values():
            w.close()
    manifest = {"format": 1, "versions": {"numpy": np.__version__, "pyarrow": pa.__version__},
                "files": {rel: {"sha256": file_sha256(root / rel), "rows": rows[rel],
                                "bytes": (root / rel).stat().st_size} for rel in sorted(rows)}}
    write_json(root / "manifest.json", manifest)
    return manifest


@functools.cache
def _verified_stat(path: str, size: int, mtime_ns: int, sha256: str, rows: int) -> None:
    """Hash a file against its manifest entry once per process: a failure raises and is not cached."""
    if file_sha256(path) != sha256 or pq.ParquetFile(path).metadata.num_rows != rows:
        raise ValueError(f"{path}: SHA-256 or row count differs from the manifest")


def _verified(root: Path, rel: str, entry: dict) -> None:
    path = (root / rel).resolve()
    st = path.stat()
    _verified_stat(str(path), st.st_size, st.st_mtime_ns, entry["sha256"], entry["rows"])


def verify(root) -> dict:
    """Check every file against the manifest (SHA-256, rows); returns the manifest, ValueError on mismatch."""
    root = Path(root)
    manifest = read_json(root / "manifest.json")
    for rel, entry in manifest["files"].items():
        _verified(root, rel, entry)
    return manifest


def load(dataset: str, split: str, root) -> Iterator[Item]:
    """Stream one dataset/split in item_id order, after checking the file against the manifest."""
    root = Path(root)
    rel = f"{dataset}/{split}.parquet"
    _verified(root, rel, read_json(root / "manifest.json")["files"][rel])
    for batch in pq.ParquetFile(root / rel, pre_buffer=False).iter_batches(batch_size=BATCH_ROWS):
        for r in batch.to_pylist():
            it = Item(item_id=r["item_id"], dataset=r["dataset"], group_id=r["group_id"], split=r["split"],
                      image=None if r["image"] is None else _from_npy(r["image"]),
                      masks={m["name"]: decode_masks(m["bits"], r["width"]) for m in r["masks"]},
                      raters=r["raters"], meta=json.loads(r["meta"]))
            if (it.dataset, it.split) != (dataset, split):  # review #12: a row from another file is a split leak (K11)
                raise ValueError(f"{rel}: row {it.item_id} belongs to {it.dataset}/{it.split}")
            _check(it)
            yield it
