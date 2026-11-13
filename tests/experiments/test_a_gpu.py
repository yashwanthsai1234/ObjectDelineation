"""Phase A tests that need a GPU box or the Mac's data (Implementation.md A.9.6). One marker line per test: `gpu`
tests run on the box's CUDA build (make test on a GPU box) and skip on the CPU build; `data` tests read data/ and
runs/ on the Mac."""

import json
import os
import resource
import subprocess
import sys
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from objdel.core.hashing import file_sha256
from objdel.geom import engine as _gpucheck
from objdel.verify import ceiling
from objdel.verify.reward import score
from objdel.verify.reward.cache import ExpertCache
from tests.goldens import make_a_goldens as GOLD
from tests.support.paths import ROOT

PREPARED = Path(os.environ.get("OBJDEL_PREPARED", ROOT / "data" / "prepared"))  # the box's prepared root
CUDA = pytest.mark.skipif(_gpucheck.build() != "cuda", reason="CPU build: the CUDA build runs on the box")


def _golden() -> dict:
    gold = json.loads((ROOT / "tests" / "goldens" / "a_g02b.json").read_text())  # a missing golden fails
    assert GOLD.texts_sha256(gold["texts"]) == gold["texts_sha256"]
    return gold


@pytest.mark.gpu
@CUDA
def test_g02b_builds():
    """G02b: the CUDA build's doctor and model of ceiling.item equal the Mac CPU build's golden to the bit."""
    gold = _golden()
    assert GOLD.values(PREPARED, gold["texts"]) == gold["values"]


@pytest.mark.gpu
@CUDA
def test_g03_device_bytes(monkeypatch):
    """G03: the ceiling scores on the leave-one-out caches' device tables, and only per-object numbers come back
    (no returned array as large as a mask)."""
    from objdel.policy import rows as dataset
    from objdel.verify import store as verify_store

    called, sizes = set(), []
    real_device = ExpertCache.device

    def device(self):
        called.add(self.digest)
        return real_device(self)

    def spy(f):
        def g(*a):
            r = f(*a)
            sizes.append((len(a[-1]), max(np.asarray(v).size for v in r.values())))  # (objects, largest field)
            return r
        return g

    proxy = SimpleNamespace(**{n: getattr(_gpucheck, n) for n in dir(_gpucheck) if not n.startswith("__")})
    proxy.score_masks, proxy.score_polygons = spy(_gpucheck.score_masks), spy(_gpucheck.score_polygons)
    monkeypatch.setattr(ExpertCache, "device", device)
    monkeypatch.setattr(score, "_gpucheck", proxy)
    caches, gold = verify_store.cache_loader(PREPARED), _golden()
    recs = {r["item_id"]: r for r in dataset.records(PREPARED, "lidc", "val")}
    for item_id, texts in list(gold["texts"].items())[:5]:
        cache, r = caches(item_id)["nodule"], recs[item_id]
        ceiling.item(cache, texts, [False] * 8, r["W"], r["H"])
        assert {lo.cache.digest for lo in ceiling.loo_caches(cache)} | {cache.digest} <= called
    assert sizes and all(big <= 64 * n for n, big in sizes)


@pytest.mark.data
@pytest.mark.parametrize("split, legacy, found, n_found", [("val", .723, .734, 947), ("test", .732, .740, 1962)])
def test_g04_legacy_ceiling(split, legacy, found, n_found, tmp_path, clean_store):
    """G04 (A-G2): the legacy LOO d_max values of the ceiling CLI on the shifted view are the old ceiling's, on all
    items and on R5's found items (K-mean IoU_tc >= 0.1)."""
    from objdel.cli.eval import ceiling as CLI
    from objdel.verify import stats as signal

    view = ROOT / "data" / "eval_views" / f"lidc-{split}-shifted"
    if not (view / "manifest.json").is_file():
        pytest.skip(f"no {view.name} here (the test view is built on the box, A3)")
    CLI.main(["--view", str(view), "--split", split, "--out", str(tmp_path)])
    s = json.loads((tmp_path / "ceiling.json").read_text())
    a = s["all"]
    if split == "val":
        assert (a["n"], s["k=2"]["n"], s["k>=3"]["n"]) == (1084, 281, 803)
    assert min(abs(a["legacy_loo_dmax_item"][0] - legacy), abs(a["legacy_loo_dmax_pooled"] - legacy)) <= 0.005
    r5 = ROOT / "runs" / "r5" / "eval" / f"eval_{split}" / "items.jsonl"
    if not r5.is_file():
        pytest.skip(f"found-stratum line: no {r5.relative_to(ROOT)}")
    with open(r5) as f:
        hit = {d["item_id"] for d in map(json.loads, f) if d["iou_tc"] >= signal.MISS}
    with open(tmp_path / "ceiling_items.jsonl") as f:
        dmax = [d["legacy_loo_dmax_item"] for d in map(json.loads, f) if d["item_id"] in hit]
    assert len(dmax) == n_found and abs(np.mean(dmax) - found) <= 0.005


A3 = Path(os.environ.get("OBJDEL_A3_RUN", ROOT / "runs" / "a3"))  # pulled: eval_val, eval_val_edge, eval_test
R5 = ROOT / "runs" / "r5" / "eval"
VIEWS = ROOT / "data" / "eval_views"


def _need(path: Path) -> Path:
    if not path.exists():
        pytest.skip(f"no {path} (A3's box run is pulled after rental V1 / V3)")
    return path


def _jsonl(path: Path) -> list[dict]:
    with open(path) as f:
        return [json.loads(line) for line in f]


@pytest.mark.data
@pytest.mark.parametrize("folder", ["eval_val", "eval_val_edge", "eval_test"])
def test_g02a_rescore(folder, tmp_path, clean_store):
    """G02a (A-G3): final_eval --rescore of a pulled A3 folder on the Mac's view reproduces the box's per-sample
    legal, iou_tc, f5 and g_area to the bit."""

    run = _need(A3 / folder)
    doc = json.loads((run / "final_eval.json").read_text())
    view = _need(VIEWS / Path(doc["view"]).name)
    from objdel.cli.eval import run as CLI

    CLI.main(["--rescore", str(run), "--view", str(view), "--split", doc["split"], "--out", str(tmp_path)])
    for a, b in zip(_jsonl(run / "metrics.jsonl"), _jsonl(tmp_path / "metrics.jsonl"), strict=True):
        assert a["item_id"] == b["item_id"]
        for k in ("legal", "iou_tc", "f5", "g_area"):
            assert a["samples"][k] == b["samples"][k], (a["item_id"], k)


@pytest.mark.data
def test_g05_outline_lengths(tmp_path):
    """G05 (A8, A-G8): the CLI over the real LIDC train roots: 7,341 records, 23,469 targets, longest outline 1,182
    tokens (with the end token), none above 1,200; the target count equals the independent sum over the records; the
    shifted and rescaled view roots carry their manifest's SHA-256; under 120 s of CPU
    time (the child's user + system time, so a busy host does not fail it)."""
    prep = ROOT / "data" / "prepared"
    views = [ROOT / "data" / "views" / n for n in ("lidc-train-shifted", "lidc-train-rescaled")]
    for p in (prep, *views):
        _need(p / "manifest.json")
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    subprocess.run([sys.executable, "-I", "-m", "tools.model.outline_lengths", "--train", prep, "--views",
                    *views, "--model", "qwen3-vl-2b", "--datasets", "lidc", "--splits", "train", "--out", tmp_path],
                   cwd=ROOT, check=True, capture_output=True)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    assert after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime < 120
    roots = json.loads((tmp_path / "outline_lengths.json").read_text())["roots"]
    row = roots[str(prep)]["lidc"]["train"]
    recs = _jsonl(prep / "lidc" / "train.jsonl")
    assert (row["n_records"], row["n_targets"], row["max"], row["above"]["1200"]) == (7341, 23469, 1182, 0)
    assert row["n_targets"] == sum(t is not None for r in recs for t in r["targets"])
    for v in views:
        assert roots[str(v)]["manifest_sha256"] == file_sha256(v / "manifest.json")


@pytest.mark.data
def test_g06_nsys_baseline():
    found = sorted((ROOT / "runs").glob("*/nsys/nsys_baseline.json"))
    if not found:
        pytest.skip("missing runs/<V2 run>/nsys/nsys_baseline.json (C0 baseline arm of V2 item 3)")
    for path in found:
        doc = json.loads(path.read_text())
        assert [s["step"] for s in doc["steps"]] == [16, 17]
        for s in doc["steps"]:
            assert s["coverage"] >= 0.90 and {"step", "loader", "micro", "optimizer", "log", "generate"} <= set(
                s["per_range_s"])  # every step generates once (steps_per_generation = GA)
            assert {"gpu_idle", "graph_rows", "host_syncs", "host_sync_parts"} <= set(s)
        assert "overhead_affected" in doc["overhead"]
        assert (path.parent / "c0_base_tables.txt").exists() and (path.parent / "version.txt").exists()
        rl = path.parents[1] / "c0_base" / "rl.json"
        if rl.exists():
            steps = [h for h in json.loads(rl.read_text())["log"] if "time/step_s" in h]
            keys = [f"sampling/importance_sampling_ratio/{k}_share" for k in ("zero", "hi", "lo")]
            assert steps and all(k in steps[-1] for k in keys)


@pytest.mark.gpu
@CUDA
def test_g07_repeat(tmp_path):
    """G07 (A3; optional box test, `OBJDEL_G07_CKPT` names the checkpoint folder): two final_eval runs with --seed 11
    on 16 items write answers of one schema, and scoring the same answers twice (--rescore) gives byte-identical
    metrics. Texts that differ between the runs are reported, not failed (the vLLM build state can differ)."""
    from objdel.core.io import read_jsonl

    ckpt = os.environ.get("OBJDEL_G07_CKPT")
    assert ckpt, "set OBJDEL_G07_CKPT to the checkpoint folder"
    final_eval = [sys.executable, "-I", "-m", "objdel.cli.eval.run", "--view", PREPARED]
    runs = [tmp_path / n for n in "ab"]
    for out in runs:
        subprocess.run([*final_eval, "--ckpt", ckpt, "--seed", "11", "--limit", "16",
                        "--out", out], cwd=ROOT, check=True)
    a, b = (read_jsonl(r / "answers.jsonl") for r in runs)
    assert len(a) == len(b) == 16 and [sorted(x) for x in a] == [sorted(x) for x in b]
    if differ := [x["item_id"] for x, y in zip(a, b, strict=True) if x["samples"] != y["samples"]]:
        warnings.warn(f"{len(differ)} of 16 items sampled other texts on the second run: {differ}", stacklevel=2)
    scored = []
    for n in "cd":
        subprocess.run([*final_eval, "--rescore", runs[0], "--out", tmp_path / n], cwd=ROOT, check=True)
        scored.append((tmp_path / n / "metrics.jsonl").read_bytes())
    assert scored[0] == scored[1]


@pytest.mark.gpu
@CUDA
def test_g08_biou_device():
    """G08 (A5, D7): the build's boundary-IoU bands (gc.eval_sets) equal the host reference to the bit on 320 random
    masks (10 sets of 16 samples against 16 raters on a 512 frame), a full-frame sample among them, at d = 2 and 16."""
    from tests.verify.test_eval_device import _cache, _call, _expect, _polys

    rng = np.random.default_rng(8)
    sets = []
    for i in range(10):
        objs = _polys(rng, 16, 512, 512, "random")
        if i == 0:  # touches every frame edge
            objs[0] = ([(np.array([[0, 0], [1000, 0], [1000, 1000], [0, 1000]]), [])], 512, 512)
        sets.append({"objs": objs, "cache": _cache(rng, 16, 512, 512), "legal": np.ones(16, bool)})
    for d in (2, 16):
        got = _call(sets, 16, d=d)
        for i, s in enumerate(sets):
            for key, v in _expect(s, 16, d=d).items():
                assert np.array_equal(got[key][i], v), (d, i, key)
