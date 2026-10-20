"""Evaluation harness (plan §4.7, judge R10): the training checker, parser, rasterizer and reward applied to K
sampled completions per item.

    items = [{"item_id": ..., "dataset": ..., **metrics(texts, frame.masks, frame.W, frame.H)} for ...]
    report = {"summary": summarize(items), "audit": hack_audit(items, baseline)}
    gates.phase2_entry(report)                    # per dataset; datasets=("riga",) for a model dropped from LIDC
    audit = alarms_for([(key, text, truncated, W, H, names), ...], cache_loader, baseline)  # in training
"""
from . import gates
from .audit import alarms_for, hack_audit
from .metrics import metrics, summarize

__all__ = ["alarms_for", "gates", "hack_audit", "metrics", "summarize"]
