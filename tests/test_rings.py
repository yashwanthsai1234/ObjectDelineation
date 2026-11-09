"""The polygonizer seam (B1.4a, D3): objdel.checker.rings.trace and polygonize (the engine,
objdel._gpucheck.check_masks; the CPU build on the Mac) give the numpy oracles' rings and bit-quads
(tests/oracles/polygonize_oracle.py, bitquad_oracle.py): same order, corners and signed area; ring parts are the
first-pixel ranks of cv2's 8-components and b0, b1 those of cv2 and scipy. The oracle guard keeps every ring reference
numpy, so a ring test can still fail after the switch."""
import ast
import dataclasses
import importlib
import inspect
import json
import shutil
import sys
from pathlib import Path

import cv2
import numpy as np
import pytest
import scipy.ndimage as ndi

import objdel.checker
from experiments.cuda_checker import compare as C
from objdel.checker import certify, rings
from tests.oracles import bitquad_oracle as BQ
from tests.oracles import polygonize_oracle as PO
from tests.test_cuda_checker import k1_masks

ROOT = Path(__file__).resolve().parents[1]
MASKS = ROOT / "data" / "checker_eval" / "masks.npz"


def _flat(rs):
    return [(r.pts.tolist(), r.area) for r in rs]


def _canon(rs):
    return [(r.pts.tolist(), r.area, r.part) for r in rs]


def _families(seed: int):
    """The D-G1 families: C.fuzz_masks (edge cases + the ten random_masks families), K1's generators (noise 0.05-0.95
    at sides 31/32/33, blobs with holes, checkerboard, nested squares, comb, spiral, diagonals, stripes at four sizes,
    1xN, Nx1, empty), both checkerboard phases at 16/17 and nested squares of 1..10 levels."""
    out = [m for _, m in C.fuzz_masks(400, seed=seed)] + k1_masks(seed)
    out += [np.indices((n, n)).sum(0) % 2 == p for n in (16, 17) for p in (0, 1)]
    for levels in range(1, 11):
        n = 4 * levels + 1
        y, x = np.indices((n, n))
        out.append(np.minimum(np.minimum(x, y), np.minimum(n - 1 - x, n - 1 - y)) % 4 < 2)
    return out


def _first_pixel_rank_parts(m, rs):
    """Per ring, the first-pixel rank of the cv2 8-component right of its first edge (C.ring_parts)."""
    return list(C.ring_parts(m, [(r.pts, r.area) for r in rs]))


def test_trace_equals_the_oracles_in_one_engine_call(monkeypatch):
    """P01/P02 at the seam: rings (corners from the canonical start, order, signed area) and quads equal the numpy
    oracles; part = first-pixel rank of the cv2 8-component right of the first edge; b0 = cv2 8-components, b1 =
    scipy 4-connected bounded background components, b0 - b1 = chi8 of the quads; one check_masks call."""
    masks = _families(seed=5)
    calls, real = [], rings._gpucheck.check_masks
    monkeypatch.setattr(rings._gpucheck, "check_masks", lambda *a, **k: calls.append(1) or real(*a, **k))
    got = rings.trace(masks)
    assert len(calls) == 1 and len(got) == len(masks) > 480
    for i, (m, g) in enumerate(zip(masks, got, strict=True)):
        want, q = PO.polygonize(m), BQ.bitquads(m)
        assert _flat(g.rings) == _flat(want), i
        assert g.quads == (q.n1, q.n2, q.nD, q.n3, q.n4), i
        assert [r.part for r in g.rings] == _first_pixel_rank_parts(m, want), i
        b0 = cv2.connectedComponents(m.astype(np.uint8), connectivity=8)[0] - 1 if m.size else 0
        b1 = ndi.label(~np.pad(m, 1))[1] - 1
        assert (g.b0, g.b1) == (b0, b1) == (sum(not r.hole for r in g.rings), sum(r.hole for r in g.rings)), i
        assert b0 - b1 == q.chi8, i


def test_trace_empty_masks_between_others():
    """U01/U02 at the seam: empty and 0-sized masks give no rings, zero quads and (0, 0); neighbours equal their solo
    trace; a single pixel gives one 4-corner ring of area 1, part 0, quads (4, 0, 0, 0, 0)."""
    px = np.zeros((5, 5), bool)
    px[2, 3] = True
    empty = [np.zeros(s, bool) for s in ((0, 0), (0, 7), (7, 0), (1, 0), (1, 1), (3, 3))]
    masks = [px, *empty[:3], px.T, *empty[3:], np.ones((1, 1), bool)]
    got = rings.trace(masks)
    for m, g in zip(masks, got, strict=True):
        if not m.any():
            assert (g.rings, g.quads, g.b0, g.b1) == ([], (0, 0, 0, 0, 0), 0, 0)
        else:
            solo = rings.trace([m])[0]
            assert (_canon(g.rings), g.quads, g.b0, g.b1) == (_canon(solo.rings), solo.quads, solo.b0, solo.b1)
    assert _canon(got[0].rings) == [([[3, 2], [4, 2], [4, 3], [3, 3]], 1, 0)] and got[0].quads == (4, 0, 0, 0, 0)
    assert rings.trace([]) == []


def test_trace_offsets():
    """U07: offsets move every corner (int64, no wrap); area, part and order stay; polygonize is trace of one."""
    m = np.zeros((6, 7), bool)
    m[1:5, 1:3] = m[3:5, 1:6] = True
    m[0, 6] = True
    base = PO.polygonize(m)
    offs = [(0, 0), (5, -3), (-7, 100), (4096, 4096), (2 ** 20, 2 ** 20), (-1000, -1000)]
    got = rings.trace([m] * len(offs), offs)
    for o, g in zip(offs, got, strict=True):
        assert [(r.pts.tolist(), r.area) for r in g.rings] == [((b.pts + o).tolist(), b.area) for b in base]
        assert all(r.pts.dtype == np.int64 for r in g.rings)
        assert [r.part for r in g.rings] == [0, 1]
        assert _canon(rings.polygonize(m, o)) == _canon(g.rings)
    assert _canon(rings.trace([m, m.T])[1].rings) == _canon(rings.trace([m, m.T], [(0, 0)] * 2)[1].rings)
    with pytest.raises(ValueError, match="offsets"):
        rings.trace([m, m], [(0, 0)])


def test_trace_dtype_and_layout():
    """U08: any nonzero value is foreground, any layout reads as the C-contiguous mask; the input is not changed;
    1-D, 3-D and 0-D raise ValueError."""
    rng = np.random.default_rng(3)
    m = rng.random((9, 12)) < .5
    want = _canon(rings.trace([m])[0].rings)
    ro = m.copy()
    ro.flags.writeable = False
    variants = [m.astype(np.uint8), m.astype(np.uint8) * 255, np.where(m, rng.choice([-3, 7], m.shape), 0),
                m.astype(np.float32) * .5, m.tolist(), np.asfortranarray(m), ro]
    before = [np.array(v, copy=True) for v in variants]
    for v in variants:
        assert _canon(rings.trace([v])[0].rings) == want
    for v, b in zip(variants, before, strict=True):
        assert np.array_equal(np.asarray(v), b)
    assert not ro.flags.writeable
    big = rng.random((20, 40)) < .5
    for view in (big.T, big[::2, ::3]):
        assert _canon(rings.trace([view])[0].rings) == _canon(rings.trace([np.ascontiguousarray(view)])[0].rings)
    for bad in (np.ones(4, bool), np.ones((2, 2, 2), bool), np.array(True)):
        with pytest.raises(ValueError):
            rings.trace([bad])


def test_ring_and_mask_rings_values():
    """U09: Ring is frozen, hole == area < 0, pts int64 (c, 2); quads a 5-tuple; b0, b1, area, part Python ints."""
    m = np.ones((4, 4), bool)
    m[1:3, 1:3] = False
    g = rings.trace([m])[0]
    assert type(g) is rings.MaskRings and [type(r) for r in g.rings] == [rings.Ring] * 2
    assert [r.hole for r in g.rings] == [False, True] and [r.area for r in g.rings] == [16, -4]
    assert all(r.pts.dtype == np.int64 and r.pts.shape == (4, 2) for r in g.rings)
    assert len(g.quads) == 5
    assert all(type(v) is int for v in (*g.quads, g.b0, g.b1, *(r.area for r in g.rings), *(r.part for r in g.rings)))
    with pytest.raises(dataclasses.FrozenInstanceError):
        g.rings[0].area = 3
    with pytest.raises(dataclasses.FrozenInstanceError):
        g.b0 = 2


@pytest.mark.parametrize("mask, offset, expected", [
    (np.zeros((5, 7), bool), (0, 0), []),
    (np.zeros((0, 0), bool), (3, 4), []),
    (np.ones((1, 1), bool), (0, 0), [([[0, 0], [1, 0], [1, 1], [0, 1]], 1)]),
    (np.ones((1, 1), bool), (10, 20), [([[10, 20], [11, 20], [11, 21], [10, 21]], 1)]),
    (np.ones((500, 500), bool), (0, 0), [([[0, 0], [500, 0], [500, 500], [0, 500]], 250000)]),
    (np.pad(np.zeros((1, 1), bool), 1, constant_values=True), (0, 0),
     [([[0, 0], [3, 0], [3, 3], [0, 3]], 9), ([[1, 1], [1, 2], [2, 2], [2, 1]], -1)]),
])
def test_literal_rings(mask, offset, expected):
    """Clockwise on screen (y down) from the minimum edge key; holes counter-clockwise with negative area."""
    assert _flat(rings.polygonize(mask, offset)) == expected


def test_any_nonzero_value_is_foreground():
    """I1 at the seam: True bytes 2 or 255 and integer masks read as the canonical mask."""
    m = np.zeros((9, 9), bool)
    m[1:8, 1:8], m[4, 4] = True, False
    want = _flat(rings.polygonize(m))
    for byte in (2, 255):
        assert _flat(rings.polygonize((m.view(np.uint8) * np.uint8(byte)).view(bool))) == want
    assert _flat(rings.polygonize(m.astype(np.int32) * 7)) == want
    assert _flat(rings.polygonize(np.asfortranarray(m))) == want


@pytest.mark.slow
@pytest.mark.skipif(not MASKS.is_file(), reason=f"{MASKS} not found")
def test_trace_equals_the_oracle_on_the_lidc_crops():
    """R08 at the seam: the 37,093 crops in one trace call give 38,002 rings and 1,476,508 corners; every crop
    equals the numpy oracle, parts equal the cv2 first-pixel ranks."""
    z = dict(np.load(MASKS, allow_pickle=False))  # one read per array, not per access
    assert len(z["hw"]) == 37093
    masks = [z["bits"][z["start"][k]:z["start"][k + 1]].reshape(z["hw"][k]).astype(bool) for k in range(37093)]
    got = rings.trace(masks)
    assert sum(len(g.rings) for g in got) == 38002 and sum(len(r.pts) for g in got for r in g.rings) == 1476508
    for k, (m, g) in enumerate(zip(masks, got, strict=True)):
        want = PO.polygonize(m)
        assert _flat(g.rings) == _flat(want), k
        assert [r.part for r in g.rings] == _first_pixel_rank_parts(m, want), k


# ------------------------------------------------------------------ the seam is production's polygonizer
def test_production_callers_use_the_engine():
    from objdel.targets import build
    certify_module = importlib.import_module("objdel.checker.certify")  # the package's `certify` is the function
    assert objdel.checker.polygonize is rings.polygonize and objdel.checker.trace is rings.trace
    assert certify_module.trace is rings.trace and build.trace is rings.trace
    assert certify(np.ones((3, 3), bool)).certified


def test_the_numpy_oracles_left_the_package():
    """S01, static half (D3): no numpy polygonizer or bit-quads in src, no src import of the test oracles, the
    package exports the seam and no bitquads."""
    for name in ("polygonize.py", "bitquad.py"):
        assert not (ROOT / "src" / "objdel" / "checker" / name).exists()
    for path in (ROOT / "src").rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            names = [node.module or ""] if isinstance(node, ast.ImportFrom) else \
                [a.name for a in node.names] if isinstance(node, ast.Import) else []
            assert not any(n.startswith("tests") or "oracle" in n for n in names), path
    with pytest.raises(ImportError):
        from objdel.checker import bitquads  # noqa: F401
    for name in ("trace", "polygonize", "Ring", "MaskRings", "certify", "certify_many", "as_mask", "betti",
                 "instance_topology", "min_area", "signature", "size_filter"):
        assert name in objdel.checker.__all__ and hasattr(objdel.checker, name), name


# ------------------------------------------------------------------ oracle guard (final round 12, qa-1)
def test_mask_reference_imports_the_numpy_oracles():
    tree = ast.parse(inspect.getsource(sys.modules[C.mask_reference.__module__]))
    sources = {(a.name, node.module) for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
               for a in node.names if a.name in ("polygonize", "bitquads")}
    assert sources == {("polygonize", "tests.oracles.polygonize_oracle"), ("bitquads", "tests.oracles.bitquad_oracle")}


def test_importing_the_oracle_keeps_the_package_function():
    import tests.oracles.polygonize_oracle  # noqa: F401
    assert inspect.isfunction(objdel.checker.polygonize)
    assert objdel.checker.polygonize is rings.polygonize


def test_the_reference_does_not_come_from_the_seam(monkeypatch):
    m = np.zeros((5, 9), bool)
    m[1:4, 1:4] = m[1:4, 5:8] = True
    want = C.mask_reference(m)
    assert want["n_rings"] == 2
    for owner in (rings, objdel.checker):
        orig = owner.polygonize
        monkeypatch.setattr(owner, "polygonize", lambda mask, offset=(0, 0), f=orig: f(mask, offset)[:-1])
    assert len(objdel.checker.polygonize(m)) == 1  # the patch bites
    got = C.mask_reference(m)
    assert (got["n_rings"], got["ring_sha"]) == (2, want["ring_sha"])


# ------------------------------------------------------------------ the B1.4a hash gate (verify_target_hashes)
@pytest.fixture(scope="module")
def built(tmp_path_factory, synthetic_items):
    """A synthetic store, its prepared root and val views, and a phase1-layout target manifest whose hashes come from
    the prepared records (every synthetic LIDC item but z002 has pixels; RIGA has one per split)."""
    pytest.importorskip("datasets")
    from experiments.phase2 import build_eval_views as BE
    from objdel.data import store
    from objdel.targets import target_set_hash
    from objdel.train import dataset as D
    root = tmp_path_factory.mktemp("hashes")
    items = [it for it in synthetic_items(("train", "val")) if it.image is not None]
    store.write(items, root / "store")
    D.prepare(root / "store", root / "prepared", procs=1)
    BE.build(root / "store", root / "views", ("lidc", "riga"), "val", seed=0, procs=1)
    prep = json.loads((root / "prepared" / "manifest.json").read_text())
    every = [(r["item_id"], k, t) for ds, split in (s.split("/") for s in prep["sets"])
             for r in D.records(root / "prepared", ds, split) for k, t in zip(r["raters"], r["targets"], strict=True)
             if t]
    (root / "targets").mkdir()
    (root / "targets" / "manifest.json").write_text(json.dumps({
        "store_manifest_sha256": prep["store_manifest_sha256"], "target_set_sha256": target_set_hash(every),
        "files": {f"{k}.jsonl": {"sha256": s["target_set_sha256"]} for k, s in prep["sets"].items()}}))
    return root


def _verify(root, manifests):
    from experiments.phase2 import verify_target_hashes as V
    return V.main([*map(str, manifests), "--store", str(root / "store"), "--procs", "2"])


def test_hash_gate_passes_on_every_recorded_set(built, capsys):
    ms = [built / "prepared" / "manifest.json", built / "targets" / "manifest.json",
          *sorted((built / "views").glob("*/manifest.json"))]
    assert len(ms) == 7 and _verify(built, ms) == 0
    out = capsys.readouterr().out
    assert sum(line.startswith("ok ") for line in out.splitlines()) == 4 + 5 + 5 and "MISMATCH" not in out


@pytest.mark.parametrize("edit", ["hash", "store"])
def test_hash_gate_fails_on_a_changed_hash_or_store(built, tmp_path, capsys, edit):
    shutil.copytree(built / "views" / "lidc-val-rotated", tmp_path / "v")
    p = tmp_path / "v" / "manifest.json"
    m = json.loads(p.read_text())
    if edit == "hash":
        s = m["sets"]["lidc/val"]
        s["target_set_sha256"] = ("0" if s["target_set_sha256"][0] != "0" else "1") + s["target_set_sha256"][1:]
    else:
        m["store_manifest_sha256"] = "0" * 64
    p.write_text(json.dumps(m))
    assert _verify(built, [p]) == 1
    assert "MISMATCH" in capsys.readouterr().out
