"""Controller tests moved from tests/test_eval.py (CLEAN6, from commit f531f86). Archived: never collected, linted or
run. To run them, restore the Controller from the tag pre-clean-2026-11-08 (README.md)."""
from test_eval import S, ExpertCache, hack_audit, metrics, nod, text_of  # live helpers (tests/test_eval.py)

from audit import alarms_for  # archive/controller/audit.py


def test_alarms_for_training_evaluations(nod):
    """B4-4: alarms_for audits what a training evaluation holds (per completion: key, text, truncated, W, H, names;
    a cache loader) exactly as hack_audit audits eval.metrics items; the baseline is a previous alarms_for result."""
    raters, texts = nod
    caches = {"nodule": ExpertCache.build(raters["nodule"])}
    union = text_of({"nodule": caches["nodule"].union}, ("nodule",), raters)
    base = alarms_for([("lidc/a", t, False, S, S, ["nodule"]) for t in texts * 2], {"lidc/a": caches}.__getitem__)
    assert base == hack_audit([{"item_id": "lidc/a", "dataset": "lidc", **metrics(texts * 2, raters, S, S)}])
    assert base["fired"] == []
    a = alarms_for([("lidc/a", union, False, S, S, ["nodule"])] * 8, {"lidc/a": caches}.__getitem__, baseline=base)
    assert {"duplicates", "inflation"} <= set(a["fired"])
