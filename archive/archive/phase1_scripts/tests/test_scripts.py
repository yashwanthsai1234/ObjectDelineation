"""Tests of the four archived phase 1 scripts, removed from tests/ in CLEAN5. Not collected (archive/ is outside
testpaths) and not runnable here: the scripts are gone from experiments/phase1. Run from the tag pre-clean-2026-11-08."""

import json

import pytest

@pytest.mark.parametrize("cuda, supported, expect", [(False, True, False), (True, False, False), (True, True, True)])
def test_sample_eval_uses_bf16_only_on_cuda_cards_that_support_it(cuda, supported, expect, monkeypatch):
    """A pre-Ampere card (no bf16) must not get bf16=True (TrainingArguments raises); the Mac (no CUDA) samples in
    fp32, as F3 verified the sampler."""
    torch = pytest.importorskip("torch")
    from experiments.phase1 import sample_eval
    monkeypatch.setattr(torch.cuda, "is_available", lambda: cuda)
    monkeypatch.setattr(torch.cuda, "is_bf16_supported", lambda *a, **k: supported)
    assert sample_eval.bf16() is expect

def test_smoke_grpo_group_stats_use_the_gate_definitions():
    """REVIEW #10: smoke.json's GRPO-signal numbers are the Phase 2 gates' (objdel.eval: duplicates = equal rasters,
    sigma of the legality-gated geometric terms), not text duplicates and ungated sigma."""
    pytest.importorskip("trl")
    import numpy as np

    from experiments.phase1 import smoke_grpo
    from objdel.eval import metrics, summarize
    from objdel.reward import ExpertCache

    def square(x0, side=400, dx=0):  # grid units on a 50 px frame (20 per px); dx nudges x within a pixel
        return json.dumps({"nodule": [{"outer": [[x0 + dx, x0], [x0 + dx + side, x0], [x0 + dx + side, x0 + side],
                                                 [x0 + dx, x0 + side]], "holes": []}]}, separators=(",", ":"))

    m = np.zeros((3, 50, 50), bool)
    for i, (a, b) in enumerate([(8, 28), (8, 28), (8, 30)]):
        m[i, a:b, a:b] = True
    caches = {"lidc/a": {"nodule": ExpertCache.build(m)}, "lidc/b": {"nodule": ExpertCache.build(m)}}
    same = [square(160)] * 4 + [square(160, dx=d) for d in (1, 2, 3, 4)]  # 5 texts, one raster
    spread = [square(160 + 20 * i) for i in range(8)]  # 8 rasters
    records = [("lidc/a", t, False, 50, 50, ["nodule"]) for t in same]
    records += [("lidc/b", t, False, 50, 50, ["nodule"]) for t in spread]
    adv = np.r_[np.zeros(8), np.linspace(-1, 1, 8)]
    s = smoke_grpo.group_stats(records, caches.__getitem__, adv)
    assert s["duplicate_rate"] == pytest.approx((7 / 8 + 0) / 2)  # text duplicates would give (3 / 8 + 0) / 2
    items = [{"item_id": k, "dataset": "lidc", **metrics(ts, {"nodule": m}, 50, 50, ["nodule"], [False] * 8, caches[k])}
             for k, ts in (("lidc/a", same), ("lidc/b", spread))]
    gate = summarize(items)["datasets"]["lidc"]
    assert s["duplicate_rate"] == gate["dup_rate"] and s["sigma_above_floor_any_term"] == gate["signal_share"]
    assert s["groups"] == 2 and s["nonzero_adv_groups"] == 1 and s["adv_abs_max"] == 1.0

def test_eval_demo_builds_targets_with_the_prepare_size_rule(riga, monkeypatch):
    """REVIEW #11: the demo's rater, honest and adversary targets use targets.min_areas (with the RIGA rim), the rule
    build_targets and train.dataset.prepare use, so the "raters" model is the rater's own SFT target."""
    from experiments.phase1 import eval_demo
    from objdel.targets import min_areas
    raters, texts = riga
    seen = []

    def spy(masks, names, cfg=None, min_area=None):
        seen.append(min_area)
        return build_target(masks, names, cfg, min_area)

    monkeypatch.setattr(eval_demo, "build_target", spy)
    _, out = eval_demo._texts(("riga/x", {"dataset": "riga", "masks": raters}))
    assert seen and all(m == min_areas(raters) for m in seen) and "rim" in seen[0]
    assert out["raters"] == texts
