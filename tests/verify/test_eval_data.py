"""Data-backed R10 tests on real validation items of the store: the raters' own targets score perfectly and are not
flagged; every R3 adversary, as a deterministic 64-sample model, is flagged by its per-output audit check."""
import json
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from oracles import reward_harness as H

from objdel.core import datasets
from objdel.data import store
from objdel.data.audit import hack_audit
from objdel.data.targets.build import build_target, min_areas
from objdel.data.targets.frame import frame
from objdel.geom.topology import signature
from objdel.verify.metrics import metrics
from objdel.verify.reward.cache import ExpertCache

pytestmark = pytest.mark.data

EXPECTED = {"dilated_union": "inflation", "hull": "inflation", "union": "inflation", "core": "deflation",
            "band_speck": "fragmentation", "holes_stamped": "topology_collapse", "minority_hole": "topology_collapse"}


@pytest.fixture(scope="module")
def frames(store_dir):
    """Every 16th RIGA val item; every 25th LIDC val item plus those with a rater hole (the hole adversaries)."""
    out = []
    for ds, every in (("riga", 16), ("lidc", 25)):
        for i, it in enumerate(store.load(ds, "val", store_dir)):
            fr = frame(replace(it, image=None))
            holes = ds == "lidc" and any(sum(signature(m, 5)) for m in fr.masks["nodule"])
            if i % every == 0 or holes:
                out.append((it.item_id, ds, fr.masks, fr.W))
    return out


@pytest.fixture(scope="module")
def models(frames):
    """{model: [metrics item]} for the raters' own targets and each pure adversary (64 copies)."""
    out = {"raters": []} | {a: [] for a in EXPECTED}
    for item_id, ds, masks, S in frames:
        names = datasets.get(ds).names
        caches = {n: ExpertCache.build(masks[n]) for n in names}
        ma = min_areas({n: masks[n] for n in names})
        k = caches[names[0]].k
        texts = [build_target({n: masks[n][r] for n in names}, names, min_area=ma).text for r in range(k)]
        out["raters"].append({"item_id": item_id, "dataset": ds, **metrics(texts, masks, S, S, names, caches=caches)})
        adv = {n: H.adversaries(masks[n], caches[n], np.random.default_rng(len(out["raters"]))) for n in names}
        for a in EXPECTED:
            if any(a in adv[n] for n in names):
                text = build_target({n: adv[n].get(a, caches[n].vote) for n in names}, names, min_area=ma).text
                out[a].append({"item_id": item_id, "dataset": ds,
                               **metrics([text] * 64, masks, S, S, names, caches=caches)})
    return out


def test_raters_own_targets_score_perfectly_on_real_items(models):
    items = models["raters"]
    assert len(items) >= 40 and {it["dataset"] for it in items} == {"riga", "lidc"}
    assert sum(it["modes"] > 1 for it in items) >= 10
    for it in items:
        s = it["samples"]
        assert min(s["format_ok"]) == min(s["legal"]) == min(s["topo_any"]) == min(s["passed"]) == 1.0, it["item_id"]
        for o in it["objects"].values():
            assert min(o["iou_tc"]) >= 0.99 and min(o["d_max"]) >= 0.99 and o["soft_dice"] >= 0.97, it["item_id"]
        assert it["tv"] == it["js"] == 0.0 and it["min_mode"] == it["coverage"] == 1.0, it["item_id"]


def test_audit_does_not_flag_the_raters_own_targets_on_real_items(models):
    a = hack_audit(models["raters"])
    assert a["fired"] == [], a["fired"]
    assert all(v["rhr"] == 0.0 for v in a["datasets"].values())


@pytest.mark.parametrize("name", list(EXPECTED))
def test_audit_flags_each_r3_adversary_on_real_items(models, name):
    items = models[name]
    assert len(items) >= 10, f"{len(items)} items with {name}"
    a = hack_audit(items, baseline=hack_audit(models["raters"]))
    assert EXPECTED[name] in a["fired"], (name, a["fired"], a["structures"])


def test_eval_run_cli_end_to_end(frames, store_dir, tmp_path):
    """objdel.cli.eval.score in a fresh -I interpreter: JSONL of the raters' own targets -> report files."""
    picked = [f for f in frames if f[1] == "lidc"][:3]
    rows = []
    for item_id, ds, masks, _ in picked:
        names = datasets.get(ds).names
        ma = min_areas({n: masks[n] for n in names})
        k = len(masks[names[0]])
        texts = [build_target({n: masks[n][r] for n in names}, names, min_area=ma).text for r in range(k)]
        rows.append({"item_id": item_id, "samples": texts, "truncated": [False] * len(texts)})
    (tmp_path / "c.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
    args = ["--completions", tmp_path / "c.jsonl", "--out", tmp_path / "out", "--store", store_dir, "--procs", "1"]
    subprocess.run([sys.executable, "-I", "-m", "objdel.cli.eval.score", *map(str, args)],
                   cwd=Path(__file__).resolve().parents[2], check=True, capture_output=True)
    report = json.loads((tmp_path / "out" / "report.json").read_text())
    assert report["audit"]["fired"] == [] and set(report["summary"]["datasets"]) == {"lidc"}
    for d in report["summary"]["datasets"].values():
        assert d["format"] == d["legal"] == d["pass_rate"] == 1.0
    assert (tmp_path / "out" / "report.md").read_text().startswith("# Evaluation of c.jsonl")
    assert len((tmp_path / "out" / "items.jsonl").read_text().splitlines()) == len(rows)
