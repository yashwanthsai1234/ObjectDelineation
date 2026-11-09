"""B1.4a polygonizer seam: objdel.checker.rings.polygonize (the engine, objdel._gpucheck.check_masks; the CPU build on
the Mac) gives the numpy oracle's rings (tests/oracles/polygonize_oracle.py): same order, corners and signed area.
The oracle guard keeps every ring reference numpy, so a ring test can still fail after the switch."""
import ast
import importlib
import inspect
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

import objdel.checker
from experiments.cuda_checker import compare as C
from objdel.checker import certify, rings
from tests.oracles import polygonize_oracle as PO

ROOT = Path(__file__).resolve().parents[1]
MASKS = ROOT / "data" / "checker_eval" / "masks.npz"


def _flat(rs):
    return [(r.pts.tolist(), r.area) for r in rs]


def _families(seed: int):
    """C.fuzz_masks (edge cases + the ten random_masks families) plus sides 31/32/33 noise, 1xN, Nx1, nested squares."""
    rng = np.random.default_rng(seed)
    out = [m for _, m in C.fuzz_masks(400, seed=seed)]
    out += [rng.random((s, s)) < p for s in (31, 32, 33) for p in (.05, .5, .95)]
    out += [rng.random((1, 77)) < .5, rng.random((77, 1)) < .5]
    nested = np.zeros((21, 21), bool)
    for k in range(0, 10, 2):
        nested[k:21 - k, k:21 - k] = k % 4 == 0
    return out + [nested]


def test_polygonize_equals_the_oracle_on_the_fuzz_families():
    rng = np.random.default_rng(1)
    for i, m in enumerate(_families(seed=5)):
        off = tuple(int(v) for v in rng.integers(-50, 600, 2))
        got = rings.polygonize(m, off)
        assert _flat(got) == _flat(PO.polygonize(m, off)), i
        assert all(type(r) is rings.Ring and r.pts.dtype == np.int64 for r in got)


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
def test_polygonize_equals_the_oracle_on_the_lidc_crops():
    z = dict(np.load(MASKS, allow_pickle=False))  # one read per array, not per access
    assert len(z["hw"]) == 37093
    for k in range(len(z["hw"])):
        m = z["bits"][z["start"][k]:z["start"][k + 1]].reshape(z["hw"][k]).astype(bool)
        assert _flat(rings.polygonize(m)) == _flat(PO.polygonize(m)), k


# ------------------------------------------------------------------ the seam is production's polygonizer
def test_production_callers_use_the_engine():
    from objdel.targets import build
    certify_module = importlib.import_module("objdel.checker.certify")  # the package's `certify` is the function
    assert objdel.checker.polygonize is rings.polygonize
    assert certify_module.polygonize is rings.polygonize and build.polygonize is rings.polygonize
    assert certify(np.ones((3, 3), bool)).certified


def test_the_numpy_polygonizer_left_the_package():
    assert not (ROOT / "src" / "objdel" / "checker" / "polygonize.py").exists()
    for path in (ROOT / "src").rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            names = [node.module or ""] if isinstance(node, ast.ImportFrom) else \
                [a.name for a in node.names] if isinstance(node, ast.Import) else []
            assert not any(n.startswith("tests") or "oracle" in n for n in names), path


# ------------------------------------------------------------------ oracle guard (final round 12, qa-1)
def test_mask_reference_imports_the_numpy_oracle():
    tree = ast.parse(inspect.getsource(sys.modules[C.mask_reference.__module__]))
    sources = {node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
               for a in node.names if a.name == "polygonize"}
    assert sources == {"tests.oracles.polygonize_oracle"}


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
