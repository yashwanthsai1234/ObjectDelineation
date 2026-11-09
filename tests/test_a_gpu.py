"""Phase A tests that need a GPU box or the Mac's data (Implementation.md A.9.6). One marker line per test: `gpu`
tests run on the box's CUDA build (d_gpu_test.yaml step 2) and skip on the CPU build; `data` tests read data/ and
runs/ on the Mac."""

import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from experiments.phase2.perf import make_a_goldens as GOLD
from objdel import _gpucheck
from objdel.eval import ceiling
from objdel.reward import score
from objdel.reward.cache import ExpertCache

ROOT = Path(__file__).resolve().parents[1]
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
    from objdel.train import dataset

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
    caches, gold = dataset.cache_loader(PREPARED), _golden()
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
    from experiments.phase2 import ceiling as CLI
    from objdel.eval import signal

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
