"""c_speed_check --sync: the pass rules of G-SY1..G-SY3 (the box run itself needs the real engine)."""

import pytest

torch = pytest.importorskip("torch")
from experiments.phase2 import c_speed_check as C  # noqa: E402

MB = 1 << 20


def test_sync_passes_a_fast_flat_series():
    v = C.sync_verdict([0.05] * 20, [100 * MB] * 50, equal=[True] * 3)
    assert v["passed"] and v["median_s"] == 0.05


@pytest.mark.parametrize(
    "kw",
    [
        {"seconds": [0.3] * 20},  # G-SY2: median over 0.2 s
        {"mem": [100 * MB] * 49 + [102 * MB]},  # G-SY3: a leak of 2 MB
        {"equal": [True, False, True]},  # G-SY1: one push differs from TRL's
        {"seconds": [0.05] * 19 + [0.5]},  # G-SY3: seconds not stable
    ],
)
def test_sync_fails_each_rule(kw):
    args = {"seconds": [0.05] * 20, "mem": [100 * MB] * 50, "equal": [True] * 3} | kw
    assert not C.sync_verdict(**args)["passed"]


def test_hash_params_sees_one_changed_value():
    a = {"w": torch.ones(4, dtype=torch.bfloat16)}
    b = {"w": torch.ones(4, dtype=torch.bfloat16)}
    assert C.hash_params(a.items()) == C.hash_params(b.items())
    b["w"][2] = 2
    assert C.hash_params(a.items()) != C.hash_params(b.items())


def test_parse_takes_sync_and_the_train_rl_flags():
    a = C.parse(["--sync", "--model", "qwen3vl2b", "--ckpt", "c", "--out", "o", "--train", "t", "--steps", "1"])
    assert a.sync and a.fast_sync is False


# --- C1: the five batches and the pixel check of `--pixels` (U-CB1, G-EQ1) ----------------------------------------

import numpy as np  # noqa: E402
import PIL.Image  # noqa: E402

pytest.importorskip("trl")
from datasets import Dataset, Image  # noqa: E402

from objdel.hashing import pixels_sha256  # noqa: E402
from objdel.train import dataset  # noqa: E402

G = 8


def bundle(n=170, side=32):
    """A fixture rows set: n train rows with one random RGB image each; (OpenOnce over them, key -> pixel SHA-256)."""
    imgs = [PIL.Image.fromarray(np.random.default_rng(i).integers(0, 256, (side, side, 3), dtype=np.uint8))
            for i in range(n)]
    rows = Dataset.from_list([{"key": f"item{i}", "image": im} for i, im in enumerate(imgs)])
    rows = rows.cast_column("image", Image())
    return dataset.OpenOnce(rows), {f"item{i}": pixels_sha256(np.asarray(im)) for i, im in enumerate(imgs)}


def fake_generate(seen):
    """The trainer's job without an engine: a record in the --record-batches format for the 256 rows it is given."""
    def generate(inputs, step):
        seen.append((step, inputs))
        n = len(inputs)
        ids = {k: torch.zeros(n, 4, dtype=torch.long) for k in ("prompt_ids", "completion_ids")}
        mask = {k: torch.ones(n, 4, dtype=torch.long) for k in ("prompt_mask", "completion_mask")}
        first = {id(r["image"]): i for i, r in reversed(list(enumerate(inputs)))}
        meta = {"image_id": torch.tensor([sorted(first.values()).index(first[id(r["image"])]) for r in inputs]),
                "group": torch.arange(n) // G}
        return {**ids, **mask, "image_grid_thw": torch.tensor([[1, 2, 2]] * n), "num_images": torch.ones(n).long(),
                "meta": meta}
    return generate


def test_the_draw_is_seeded_distinct_train_rows_in_five_batches_of_32():
    a = C.draw_batches(8320)
    assert a.shape == (5, 32) and len(set(a.ravel())) == 160 and a.max() < 8320
    assert (C.draw_batches(8320) == a).all() and (C.draw_batches(8320, seed=1) != a).any()
    with pytest.raises(ValueError, match="160"):
        C.draw_batches(159)


@pytest.mark.parametrize("path", ["data/eval_views/lidc-val-shifted", "/d/lidc-val-edge", "x/val", "eval_test"])
def test_a_val_view_path_raises(path):
    with pytest.raises(ValueError, match="val"):
        C.refuse_val(path)
    from pathlib import Path

    for ok in ("/workspace/data/lidc/rows/lidc-train-b1", "data/prepared", Path("runs/interval/x")):
        C.refuse_val(ok)


def test_main_refuses_a_val_root_before_building_anything():
    with pytest.raises(ValueError, match="val"):
        C.main(["--pixels", "--model", "m", "--ckpt", "c", "--out", "o", "--steps", "20",
                "--train", "data/eval_views/lidc-val-shifted"])


def test_the_advantages_are_seeded_normal_values_centred_per_group_of_8():
    a = C.centred_advantages(5)
    assert a.shape == (5, 256) and torch.allclose(a.view(5, 32, G).sum(-1), torch.zeros(5, 32), atol=1e-5)
    assert torch.equal(a, C.centred_advantages(5)) and (a[0] != a[1]).any()
    raw = torch.randn(5 * 256, generator=torch.Generator().manual_seed(0)).view(5, 32, G)
    assert torch.allclose(a.view(5, 32, G), raw - raw.mean(-1, keepdim=True))  # seed 0, as the plan says


def test_the_batch_maker_writes_exactly_five_files_that_weights_only_loads(tmp_path):
    """U-CB1: b0..b4 in the --record-batches format plus advantages, image keys and SHA-256 values."""
    rows, sha = bundle()
    seen = []
    paths = C.make_batches(fake_generate(seen), rows, len(rows), tmp_path / "batches")
    assert sorted(p.name for p in (tmp_path / "batches").iterdir()) == [f"b{i}.pt" for i in range(5)]
    assert [p.name for p in paths] == [f"b{i}.pt" for i in range(5)]
    assert [s for s, _ in seen] == [0, 1, 2, 3, 4] and all(len(x) == 256 for _, x in seen)
    for _, inputs in seen:  # groups of 8 identical rows share one image object, 32 distinct images a batch
        assert all(inputs[i]["image"] is inputs[i // G * G]["image"] for i in range(256))
        assert len({id(r["image"]) for r in inputs}) == 32
    chosen = [r["key"] for _, x in seen for r in x[::G]]
    assert len(set(chosen)) == 160  # five batches of 32 distinct rows, no row twice
    for i, p in enumerate(paths):
        rec = torch.load(p, weights_only=True)
        assert {"prompt_ids", "completion_ids", "meta", "advantages", "image_key", "image_sha256"} <= set(rec)
        assert rec["image_key"] == chosen[i * 32:(i + 1) * 32]
        assert rec["image_sha256"] == [sha[k] for k in rec["image_key"]]  # the fixture bundle's own hashes
        assert rec["advantages"].shape == (256,)
        assert torch.allclose(rec["advantages"].view(32, G).sum(-1), torch.zeros(32), atol=1e-5)
        assert torch.equal(rec["advantages"], C.centred_advantages(5)[i])


def test_the_batch_maker_leaves_a_non_empty_folder_alone(tmp_path):
    rows, _ = bundle()
    (tmp_path / "b").mkdir()
    (tmp_path / "b/old").write_text("x")
    with pytest.raises(FileExistsError):
        C.make_batches(fake_generate([]), rows, len(rows), tmp_path / "b")


def test_pixel_compare_counts_mismatches_and_names_the_first():
    """G-EQ1's comparison on CPU: the torch backend equals the HF processor on every view; a planted processor that
    changes one image shows up as one mismatch at its index."""
    from transformers import Qwen2VLImageProcessor

    ip = Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)
    imgs = [PIL.Image.fromarray(np.random.default_rng(i).integers(0, 256, (64, 64, 3), dtype=np.uint8))
            for i in range(70)]
    ok = C.compare_pixels(imgs, ip, "cpu", "torch", chunk=32)
    assert (ok["views"], ok["mismatches"], ok["first"]) == (70, 0, None) and ok["median_ms"] >= 0

    class Off(Qwen2VLImageProcessor):
        def __call__(self, images, *a, **k):
            out = super().__call__(images, *a, **k)
            if len(images) > 40 or len(images) == 6:  # chunk 32 -> 32, 32, 6: only the last chunk is touched
                out["pixel_values"][0, 0] += 1
            return out

    bad = Off(patch_size=16, temporal_patch_size=2, merge_size=2)
    r = C.compare_pixels(imgs, bad, "cpu", "torch", chunk=32)
    assert (r["mismatches"], r["first"]) == (1, 64)  # the chunk starting at view 64: its first image differs
