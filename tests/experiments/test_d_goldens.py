"""Phase D goldens written on the Mac CPU build by experiments/phase2/perf/make_d_goldens.py.

GPU02 (gpu, the box half of D-G3 and D-G10): every engine entry on the shape corpora, the replay set and 1,000 caches
of the R-D2 payload; the inputs' SHA-256 first ("input differs"), then dtype, shape and SHA-256 of every output equal
tests/goldens/d_gpu02.json. A mismatching entry's arrays go to $OBJDEL_RUN_DIR/d/gpu02/ (pull.sh brings them back).
D-R01 (always-run copy, LIDC only): 50 LIDC items' masks rebuild their pre-D target texts. D-R04 (no marker: the CPU
build here, the CUDA build on the box): score_polygons and score_masks of 200 items equal the golden floats bit for bit.
Reads no data/ except the gpu tests, which read the payload under OBJDEL_DATA_ROOT."""
import json
import os
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage

from experiments.phase2.perf import bench_checker as B
from experiments.phase2.perf import make_d_goldens as M
from objdel.data.targets.build import build_target
from objdel.geom.shape import to_mask
from objdel.geom.topology import filtered
from tests.geom.test_cuda_checker_gpu import _gpus

GOLDENS = Path(__file__).resolve().parents[1] / "goldens"
BOX_INPUTS = ("checker_eval/masks.npz", "checker_eval/polygons.npz", "prepared/lidc/val.jsonl", "prepared/caches/lidc",
              "tests/goldens/d_gpu02.json")


def _gpu02() -> dict:
    return json.loads((GOLDENS / "d_gpu02.json").read_text())


def test_gpu02_golden_lists_every_entry_and_no_riga():
    g = _gpu02()
    assert sorted(g["entries"]) == sorted(M.entries())
    assert all(e["inputs"] and e["outputs"] for e in g["entries"].values())
    assert "riga" not in json.dumps(g).lower()


@pytest.mark.parametrize("name", list(M.entries()))
@pytest.mark.gpu
@pytest.mark.skipif(_gpus() == 0, reason="no CUDA device (nvidia-smi -L)")
def test_gpu02_entry_equals_the_cpu_build(name):
    want = _gpu02()["entries"][name]
    inputs, r = M.entries()[name](B.data_root())
    assert inputs == want["inputs"], f"{name}: input differs"
    got = M.outputs(r)
    if got != want["outputs"]:
        out = Path(os.environ.get("OBJDEL_RUN_DIR", "."), "d/gpu02")
        out.mkdir(parents=True, exist_ok=True)
        np.savez(out / f"{name.replace('/', '__')}.npz", **r)
    assert got == want["outputs"], [k for k in want["outputs"] if got.get(k) != want["outputs"][k]]


@pytest.mark.gpu
@pytest.mark.skipif(_gpus() == 0, reason="no CUDA device (nvidia-smi -L)")
def test_gpu02_edt_equals_scipy():
    """D-G10 box half: the device EDT of the 1,000 caches' maps equals scipy bit for bit."""
    _, cs = M.caches(B.data_root())
    for i, m in enumerate(M.edt_maps(cs)):
        want = ndimage.distance_transform_edt(m).astype(np.float32)
        assert M.gc.edt(m.view(np.uint8)).tobytes() == want.tobytes(), i


def test_d_r01_lidc_items_rebuild_their_pre_d_texts():
    items = json.loads((GOLDENS / "d_r01_lidc.json").read_text())["items"]
    assert len(items) == 50
    assert all(it["item_id"].startswith("lidc/") and "riga" not in it["item_id"].lower() for it in items)
    for it in items:
        for r in it["raters"]:
            m = M.unrle(r["nodule"], it["H"], it["W"])
            assert build_target({"nodule": m}, ("nodule",), None, it["min_area"]).text == r["text"], it["item_id"]


def test_rle_round_trip():
    for m in (np.zeros((3, 4), bool), np.ones((2, 2), bool), np.eye(5, dtype=bool)):
        assert np.array_equal(M.unrle(M.rle(m), *m.shape), m)


def test_d_r04_scorer_floats_equal_the_golden():
    items = M.d_r04_items()
    assert len(items) == 200
    raw = [to_mask(parts, 40, 40).mask for _, parts, _ in items]
    kept = filtered(raw, [c.min_area for _, _, c in items])
    for (kind, _, _), m, (f, _, _) in zip(items, raw, kept, strict=True):  # what the size filter does to each answer
        assert {"nothing": np.array_equal(f, m), "small": 0 < f.sum() < m.sum(), "hole": f.sum() > m.sum(),
                "everything": m.any() and not f.any()}[kind], kind
    want = json.loads((GOLDENS / "d_r04_scorer.json").read_text())
    assert {k: M.floats(v) for k, v in M.d_r04_scores(items).items()} == want
