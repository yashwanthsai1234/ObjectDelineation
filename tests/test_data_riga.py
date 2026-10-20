"""RIGA staged extractor: every stage on synthetic pen drawings, plus data-backed counts and Phase 0 parity."""
from pathlib import Path

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays
from PIL import Image
from scipy import ndimage

from objdel.data import riga, store

H, W = 240, 280
CY, CX = 120, 140


def fundus(seed=0):
    rng = np.random.default_rng(seed)
    img = np.empty((H, W, 3), np.uint8)
    img[..., 0], img[..., 1], img[..., 2] = 170, 70, 30
    return (img.astype(np.int16) + rng.integers(-3, 4, img.shape)).clip(0, 255).astype(np.uint8)


def ring(r, cy=CY, cx=CX, width=3.0):
    yy, xx = np.mgrid[:H, :W]
    d = np.hypot(yy - cy, xx - cx)
    return (d >= r - width / 2) & (d < r + width / 2)


def pen(img, stroke):
    out = img.copy()
    out[stroke, 0], out[stroke, 1] = 45, 49  # the pen keeps B
    return out


def filled(stroke):
    return ndimage.binary_fill_holes(stroke)


def write_photo(root, subset, idx, prime, copies, prime_name=None):
    d = Path(root) / subset
    for s in riga.SUBSETS:
        (Path(root) / s).mkdir(parents=True, exist_ok=True)
    if prime is not None:
        Image.fromarray(prime).save(d / (prime_name or f"image{idx}prime.tif"))
    for e, c in copies.items():
        Image.fromarray(c).save(d / f"image{idx}-{e}.tif")
    return next(p for p in riga.discover(root) if p.subset == subset and p.index == idx)


def drawings(prime, n=6):
    """n raters, disc and cup radii 4 px apart per rater, so no pixel is inked by every rater."""
    return {e: pen(prime, ring(60 + 4 * e) | ring(26 + 4 * e)) for e in range(1, n + 1)}


def iou(a, b):
    return (a & b).sum() / (a | b).sum()


def test_discover_names(tmp_path):
    for s in riga.SUBSETS:
        (tmp_path / s).mkdir(parents=True)
    img = fundus()
    Image.fromarray(img).save(tmp_path / "Magrabia/MagrabiFemale/image48.tif")  # F2: no "prime" in the name
    Image.fromarray(img).save(tmp_path / "Magrabia/MagrabiFemale/image48-1.tif")
    Image.fromarray(img).save(tmp_path / "BinRushed/BinRushed1/image7-3.jpg")  # F1: no reference file
    Image.fromarray(img).save(tmp_path / "MESSIDOR/image12prime.tif")
    Image.fromarray(img).save(tmp_path / "MESSIDOR/image3-2.tif")
    (tmp_path / "MESSIDOR/Thumbs.db").write_bytes(b"x")
    photos = riga.discover(tmp_path)
    assert [p.image_id for p in photos] == ["BinRushed/BinRushed1/image007", "MESSIDOR/image003",
                                            "MESSIDOR/image012", "Magrabia/MagrabiFemale/image048"]
    assert photos[0].prime is None and list(photos[0].raters) == [3]
    assert photos[3].prime.name == "image48.tif" and list(photos[3].raters) == [1]


def test_max_rg_recovers_reference_where_some_copy_is_clean():  # F1
    ref = fundus()
    copies = [pen(ref, ring(60)), pen(ref, ring(60) | ring(30)), pen(ref, ring(30))]
    assert np.array_equal(riga.max_rg(copies), ref)


@settings(max_examples=60, deadline=None)
@given(arrays(np.int16, st.tuples(st.integers(1, 12), st.integers(1, 12)), elements=st.integers(0, 100)),
       st.integers(0, 60), st.integers(0, 60))
def test_hysteresis_matches_scikit_image(d, low, high):  # F3
    ski = pytest.importorskip("skimage.filters")
    low = min(low, high)
    assert np.array_equal(riga.hysteresis(d, low, high), ski.apply_hysteresis_threshold(d, low, high))


def test_hysteresis_keeps_faint_pixels_linked_to_strong_ones():
    d = np.array([[0, 30, 30, 70, 30, 0, 30, 30]])
    assert riga.hysteresis(d, 25, 60).tolist() == [[False, True, True, True, True, False, False, False]]


def test_stage_a_clean_drawings(tmp_path):
    prime = fundus()
    photo = write_photo(tmp_path, "MESSIDOR", 1, prime, drawings(prime))
    ref, kind, out = riga.extract(photo)
    assert kind == "prime" and np.array_equal(ref, prime)
    assert sorted(out) == [1, 2, 3, 4, 5, 6]
    for e, (disc, cup, stage) in out.items():
        assert stage == "A"
        assert np.array_equal(disc, filled(ring(60 + 4 * e))) and np.array_equal(cup, filled(ring(26 + 4 * e)))


def test_stage_b_bridges_a_pen_gap(tmp_path):  # F4
    prime = fundus()
    copies = drawings(prime)
    gap = ring(64)
    gap[:, CX - 1:CX + 2] = False  # a 3-px gap at the top and bottom of rater 1's disc outline
    copies[1] = pen(prime, gap | ring(30))
    _, _, out = riga.extract(write_photo(tmp_path, "MESSIDOR", 1, prime, copies))
    disc, cup, stage = out[1]
    assert stage.startswith("B:prime/b") and stage != "B:prime/b0"
    assert iou(disc, filled(ring(64))) > 0.98 and np.array_equal(cup, filled(ring(30)))
    assert all(out[e][2] == "A" for e in range(2, 7))


def test_stage_c_splits_touching_disc_and_cup(tmp_path):  # F5
    prime = fundus()
    copies = drawings(prime)
    cup = ring(30, cy=CY - 38)  # the cup outline touches the disc outline from inside
    copies[2] = pen(prime, ring(68) | cup)
    _, _, out = riga.extract(write_photo(tmp_path, "MESSIDOR", 1, prime, copies))
    disc, got_cup, stage = out[2]
    assert stage == "C:prime/b0"
    assert np.all(disc[got_cup]) and iou(disc, filled(ring(68))) > 0.98 and iou(got_cup, filled(cup)) > 0.85


def test_blank_copy_is_a_missing_rater(tmp_path):  # F6
    prime = fundus()
    copies = drawings(prime)
    copies[5] = prime.copy()
    _, kind, out = riga.extract(write_photo(tmp_path, "MESSIDOR", 1, prime, copies))
    assert kind == "prime" and sorted(out) == [1, 2, 3, 4, 6]


def test_reference_that_is_a_rater_copy_is_rebuilt(tmp_path):  # F7 (Magrabia Female image48)
    clean = fundus()
    copies = drawings(clean)
    photo = write_photo(tmp_path, "Magrabia/MagrabiFemale", 48, copies[1], copies, prime_name="image48.tif")
    ref, kind, out = riga.extract(photo)
    assert kind == "maxrg6:prime_is_rater1_copy" and np.array_equal(ref, clean)
    assert sorted(out) == [1, 2, 3, 4, 5, 6] and np.array_equal(out[1][0], filled(ring(64)))


def test_no_reference_file_uses_the_reconstruction(tmp_path):  # F1 end to end (BinRushed1)
    clean = fundus()
    ref, kind, out = riga.extract(write_photo(tmp_path, "BinRushed/BinRushed1", 3, None, drawings(clean)))
    assert kind == "maxrg6" and np.array_equal(ref, clean) and len(out) == 6


@pytest.mark.parametrize("subset, prime", [("BinRushed/BinRushed1", False), ("Magrabia/MagrabiFemale", True)])
def test_binrushed1_reconstruction_is_not_a_model_image(tmp_path, subset, prime):  # B3-1
    """The BinRushed1 max-R+G reconstruction keeps the raters' strokes where every copy is inked (and a seam on every
    rater edge): the item keeps its masks but has no image. The F7 reconstruction (lossless copies) stays."""
    clean = fundus()
    copies = drawings(clean)
    photo = write_photo(tmp_path, subset, 48, copies[1] if prime else None, copies,
                        prime_name="image48.tif" if prime else None)
    ref, kind, masks = riga.extract(photo)
    it = riga.item(photo, ref, kind, masks, "train", photo.image_id, [])
    assert kind.startswith("maxrg6") and it.meta["reference"] == kind
    assert it.image is None if not prime else np.array_equal(it.image, clean)
    assert it.masks["disc"].shape == (6, H, W) and it.masks["disc"].any((1, 2)).all()


def test_certified_requires_nesting_and_simple_curves():
    disc, cup = filled(ring(60)), filled(ring(25))
    assert riga.certified(disc, cup)
    assert not riga.certified(disc, filled(ring(25, cx=CX + 70)))  # cup outside the disc
    assert not riga.certified(disc & ~filled(ring(10, cx=CX + 40)), cup)  # disc with a hole
    assert not riga.certified(disc | filled(ring(5, cy=10, cx=10)), cup)  # two disc components


# ---- data-backed -----------------------------------------------------------------------------------------------

PHASE0 = Path(__file__).resolve().parents[1] / "data" / "riga_true"


@pytest.mark.data
def test_store_riga_counts(store_dir):
    n, masks, images, per_split, test_groups = 0, 0, 0, {}, set()
    for split in ("train", "val", "test"):
        for item in store.load("riga", split, store_dir):
            n += 1
            per_split[split] = per_split.get(split, 0) + 1
            assert item.masks["disc"].shape == item.masks["cup"].shape and len(item.masks["disc"]) == len(item.raters)
            # B3-1: the BinRushed1 reconstruction (F1) is not a model image
            assert (item.image is None) == (item.meta["reference"] == "maxrg6") == (item.meta["subset"].endswith("1"))
            if item.image is not None:
                assert item.image.shape == (*item.masks["disc"].shape[1:], 3) and item.image.dtype == np.uint8
            images += item.image is not None
            masks += len(item.raters)
            if split == "test":
                test_groups.add(item.group_id)
    assert n == 750 and masks == 4500 and images == 700  # the 50 BinRushed1 items have no image
    assert per_split["test"] == 95 and per_split["train"] + per_split["val"] == 655
    # Magrabia: 95 files, 88 unique photographs by NCC; keypoints (B3-8) join two more same-eye pairs
    # (MagrabiFemale image014 + 038, image023 + 045): 86 eyes
    assert len(test_groups) == 86


@pytest.mark.data
def test_store_riga_stage_a_masks_equal_phase0(store_dir):
    if not PHASE0.is_dir():
        pytest.skip("Phase 0 masks absent")
    checked = 0
    for item in store.load("riga", "val", store_dir):
        sub, idx = item.meta["subset"], item.meta["index"]
        for k, e in enumerate(item.raters):
            f = PHASE0 / f"{sub.replace('/', '_')}_image{idx}-{e}.png"
            if item.meta["stages"][e] != "A" or not f.exists():
                continue
            v = np.asarray(Image.open(f))
            assert np.array_equal(item.masks["disc"][k], v >= 100) and np.array_equal(item.masks["cup"][k], v >= 200)
            checked += 1
    assert checked >= 300


def test_item_format(tmp_path):
    prime = fundus()
    photo = write_photo(tmp_path, "BinRushed/BinRushed2", 7, prime, drawings(prime, n=2))
    ref, kind, masks = riga.extract(photo)
    it = riga.item(photo, ref, kind, masks, "val", "BinRushed/BinRushed2/image007", ["MESSIDOR/image010"])
    assert (it.item_id, it.dataset, it.group_id, it.split, it.raters) == (
        "riga/BinRushed/BinRushed2/image007", "riga", "riga:BinRushed/BinRushed2/image007", "val", ["1", "2"])
    assert it.masks["disc"].shape == it.masks["cup"].shape == (2, H, W) and it.image.dtype == np.uint8
    assert np.array_equal(it.masks["cup"][1], filled(ring(34)))
    assert it.meta == {"subset": "BinRushed/BinRushed2", "source": "BinRushed", "index": 7, "reference": "prime",
                       "stages": {"1": "A", "2": "A"}, "duplicates": ["MESSIDOR/image010"]}
