"""Break-it B3 (data, targets, prepared set): reproducers of the confirmed breakdowns, now fixed (ids B3-<n>, see
experiments/phase1/impl/breakit/B3_data.md) and the robustness checks that held (plain tests).

Data-backed tests read data/store, data/prepared (skipped when absent) and only touch a few
rows each, so they stay fast."""
import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pyarrow")
pytest.importorskip("datasets")

import pyarrow.parquet as pq  # noqa: E402
from scipy import ndimage  # noqa: E402

from objdel.data import prepare as data_prepare  # noqa: E402
from objdel.data import store  # noqa: E402
from objdel.data.store import Item  # noqa: E402
from objdel.data.targets.frame import frame  # noqa: E402
from objdel.train import dataset as D  # noqa: E402
from objdel.verify import store as verify_store  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PREPARED = ROOT / "data" / "prepared"


def need(path: Path) -> Path:
    if not path.exists():
        pytest.skip(f"{path} not found")
    return path


def disc(h, w, cy, cx, r):
    y, x = np.ogrid[:h, :w]
    return (y - cy) ** 2 + (x - cx) ** 2 <= r * r


def prepared(tmp_path) -> Path:
    """A synthetic store (one RIGA item, two LIDC items of one patient with different images and different nodule
    sizes, so their frames and caches differ) run through data.prepare.prepare."""
    rng = np.random.default_rng(0)
    hu = lambda: rng.integers(-1000, 400, (100, 100), dtype=np.int16)  # noqa: E731
    two = lambda c: np.stack([disc(100, 100, 40 + 20 * c, 50, 8 + c + d) for d in (0, 1)])  # noqa: E731
    riga = Item("riga/MESSIDOR/image001", "riga", "riga:MESSIDOR/image001", "train",
                rng.integers(0, 256, (80, 80, 3), dtype=np.uint8),
                {"disc": np.stack([disc(80, 80, 40, 40, 20)] * 2), "cup": np.stack([disc(80, 80, 40, 40, 8)] * 2)},
                ["1", "2"], {})  # prepare needs every dataset's split file
    store.write([riga, *(Item(f"lidc/1.2.3/n00{c}/z001", "lidc", "P1", "train", hu(), {"nodule": two(c)}, ["1", "2"],
                              {}) for c in (0, 1))], tmp_path / "store")
    data_prepare.prepare(tmp_path / "store", tmp_path / "out", splits=("train",), procs=1, prompts=D.PROMPTS)
    return tmp_path / "out"


def swap(a: Path, b: Path) -> None:
    ab = a.read_bytes()
    a.write_bytes(b.read_bytes())
    b.write_bytes(ab)


# ---------------------------------------------------------------- confirmed breakdowns (now regression tests)
# Pinned elsewhere, on their break-it items where there is one: B3-2, B3-3, B3-4 and the LIDC item contract by
# tests/test_data_lidc.py (crowded frame left out, copied reader merged, every item has pixels); B3-6 and B3-7 by
# tests/test_data_store.py (one id in two splits, image dtype on write and load); B3-8 by tests/test_data_splits.py
# (keypoint K11).

@pytest.mark.data
def test_b3_1_binrushed1_inputs_carry_no_rater_ink():
    need(PREPARED / "riga" / "train.jsonl")
    import cv2
    inked, caches = [], verify_store.cache_loader(PREPARED, maxsize=1)
    for r in D.records(PREPARED, "riga", "train"):
        if "BinRushed1/" not in r["item_id"]:
            continue
        img = cv2.imread(str(PREPARED / "images" / (r["item_id"] + ".png")))[..., ::-1].astype(int)
        edge = np.zeros((r["H"], r["W"]), bool)
        for n in ("disc", "cup"):
            for m in caches(r["item_id"])[n].raters:
                edge |= m & ~ndimage.binary_erosion(m)
        band = ndimage.binary_dilation(edge, iterations=2)
        dark = (img[..., 0] < 80) & (img[..., 1] < 80)  # the pen is near-black; fundus near the disc is not
        if (dark & band).sum() >= 30:
            inked.append(r["item_id"])
    # Measured: 18 / 50 BinRushed1 items >= 30 dark pixels on the rater edges (0-1 per other subset).
    assert inked == []


def test_b3_5_tampered_prepared_image_is_detected(tmp_path):
    out = prepared(tmp_path)
    a, b = (out / "images" / f"lidc/1.2.3/n00{c}/z001.png" for c in (0, 1))
    swap(a, b)  # the two model inputs
    from PIL import Image as PILImage  # the swapped file really holds other pixels than its record pins
    pixels = np.asarray(PILImage.open(a).convert("RGB"))
    assert D.pixels_sha256(pixels) != D.records(out, "lidc", "train")[0]["image_sha256"]
    with pytest.raises(ValueError, match="image_sha256"):  # building rows checks the pinned image hash
        D.grpo_rows(out, "lidc", "train")
    with pytest.raises(ValueError, match="image_sha256"):
        D.sft_rows(out, "lidc", "train")


# ---------------------------------------------------------------- robustness checks that held

def test_store_load_is_independent_of_the_reader_batch_size(tmp_path, monkeypatch):
    rng = np.random.default_rng(1)
    items = [Item(f"lidc/1.2.3/n000/z{i:03d}", "lidc", "P1", "train",
                  rng.integers(-1024, 3000, (32, 48), dtype=np.int16),
                  {"nodule": rng.random((3, 32, 48)) < 0.3}, ["1", "2", "4"], {"i": i}) for i in range(37)]
    store.write(items, tmp_path)
    got = {}
    for b in (1, 5, 10 ** 6):
        monkeypatch.setattr(store, "BATCH_ROWS", b)
        got[b] = list(store.load("lidc", "train", tmp_path))
    for b in got:
        assert [x.item_id for x in got[b]] == [x.item_id for x in items]
        for x, y in zip(got[b], items, strict=True):
            assert np.array_equal(x.image, y.image) and x.image.dtype == np.int16
            assert np.array_equal(x.masks["nodule"], y.masks["nodule"]) and x.raters == y.raters and x.meta == y.meta


@pytest.mark.parametrize("attack", ["truncate", "append", "swap", "flip_last", "flip_middle"])
def test_store_tampering_is_detected(tmp_path, attack):
    rng = np.random.default_rng(2)
    mk = lambda s, i: Item(f"lidc/1.2.3/n000/z{i:03d}", "lidc", "P1", s,  # noqa: E731
                           rng.integers(-1024, 3000, (16, 16), dtype=np.int16),
                           {"nodule": rng.random((2, 16, 16)) < 0.5}, ["1", "2"], {})
    store.write([mk("train", 1), mk("train", 2), mk("val", 3)], tmp_path)
    tr, va = tmp_path / "lidc" / "train.parquet", tmp_path / "lidc" / "val.parquet"
    data = tr.read_bytes()
    if attack == "truncate":
        tr.write_bytes(data[:-100])
    elif attack == "append":
        tr.write_bytes(data + b"\0")
    elif attack == "swap":
        swap(tr, va)
    else:
        i = len(data) - 1 if attack == "flip_last" else len(data) // 2
        tr.write_bytes(data[:i] + bytes([data[i] ^ 1]) + data[i + 1:])
    with pytest.raises(ValueError, match="lidc/train.parquet"):
        store.verify(tmp_path)
    with pytest.raises(ValueError, match="lidc/train.parquet"):
        next(store.load("lidc", "train", tmp_path))


def test_frame_scale_gt_1_maps_grid_to_item_pixels_at_the_border():
    """f > 1 with a window larger than the image (padded on both sides): the frame mask is the block majority of
    the item mask on the box, and every rater pixel's block lands inside the frame."""
    h = w = 900
    m = disc(h, w, 450, 450, 440)
    it = Item("riga/MESSIDOR/image001", "riga", "riga:x", "train", np.full((h, w, 3), 200, np.uint8),
              {"disc": np.stack([m, m]), "cup": np.stack([disc(h, w, 450, 450, 100)] * 2)}, ["1", "2"], {})
    f = frame(it)
    x0, y0, x1, y1 = f.box
    assert f.scale > 1 and x1 - x0 == f.scale * f.W and x0 < 0 and x1 > w
    pad = np.zeros((x1 - x0, x1 - x0), bool)
    pad[-y0:-y0 + h, -x0:-x0 + w] = m
    s = f.scale
    vote = pad.reshape(f.H, s, f.W, s).sum((1, 3)) * 2 >= s * s
    assert np.array_equal(f.masks["disc"][0], vote)
    assert f.image[0, 0].tolist() == [0, 0, 0] and f.image[f.H // 2, f.W // 2].tolist() == [200, 200, 200]


@pytest.mark.data
def test_prepared_rows_match_records_and_oversampling(store_dir):
    """Every prepared record has its PNG and caches; sft rows carry one image placeholder; grpo oversampling counts
    match; record keys are unique and in the same store split (no leakage introduced by prepare)."""
    need(PREPARED / "manifest.json")
    m = json.loads((PREPARED / "manifest.json").read_text())
    for ds in ("riga", "lidc"):
        for split in ("train", "val"):
            recs = D.records(PREPARED, ds, split)
            ids = [r["item_id"] for r in recs]
            assert len(ids) == len(set(ids)) == m["sets"][f"{ds}/{split}"]["items"]
            col = pq.read_table(store_dir / ds / f"{split}.parquet", columns=["item_id"]).column(0).to_pylist()
            assert set(ids) <= set(col)
            assert all((PREPARED / "images" / (i + ".png")).is_file() for i in ids)
            dis = sum(r["disagreement"] for r in recs)
            g = D.grpo_rows(PREPARED, ds, split, oversample=D.DISAGREEMENT_OVERSAMPLE)
            assert len(g) == sum(r["rl"] for r in recs) + (D.DISAGREEMENT_OVERSAMPLE - 1) * dis
            assert dis == m["sets"][f"{ds}/{split}"]["disagreement_items"]
    s = D.sft_rows(PREPARED, "lidc", "val")
    row = s[0]
    assert sum(c["type"] == "image" for c in row["messages"][0]["content"]) == 1 and row["image"] is not None


@pytest.mark.parametrize("side", [6, 40, 90])
def test_lidc_upscaled_image_lines_up_with_the_mask(side):
    """Frame image (lung window, bicubic to 512) vs frame mask: a bright square drawn exactly on the mask spans the
    same model-image pixels as the mask scaled by 512 / S (half-pixel centres), within one output pixel."""
    hu = np.full((512, 512), -1000, np.int16)
    m = np.zeros((512, 512), bool)
    m[200:200 + side, 300:300 + side] = True
    hu[m] = 400
    it = Item("lidc/1.2.3/n000/z001", "lidc", "P1", "train", hu, {"nodule": np.stack([m, m])}, ["1", "2"], {})
    f = frame(it)
    k = 512 / f.W
    ys, xs = np.nonzero(f.masks["nodule"][0])
    bright = f.image[..., 0] >= 128
    by, bx = np.nonzero(bright)
    spans = ((ys.min(), ys.max() + 1, by.min(), by.max() + 1), (xs.min(), xs.max() + 1, bx.min(), bx.max() + 1))
    for lo, hi, blo, bhi in spans:
        assert abs(blo - lo * k) <= 1 and abs(bhi - hi * k) <= 1


def test_cache_loader_detects_swapped_caches(tmp_path):
    out = prepared(tmp_path)
    swap(*(out / "caches" / f"lidc/1.2.3/n00{c}/z001" / "nodule.npz" for c in (0, 1)))  # two items' caches
    with pytest.raises(ValueError, match="digests differ"):
        verify_store.cache_loader(out)("lidc/1.2.3/n000/z001")
