"""Evaluation harness (plan §4.7, judge R10): the training checker, parser, rasterizer and reward applied to K
sampled completions per item.

    items = [{"item_id": ..., "dataset": ..., **metrics(texts, frame.masks, frame.W, frame.H)} for ...]
    view = {"summary": summarize(items), "audit": hack_audit(items, baseline)}       # one view (eval_run)
    gates.phase2_entry({"views": {"shifted": view, "rescaled": ...}, "prior": ..., "probes": ...,
                        "vision": vision.check(model, processor, images, items, reference=base)}, datasets)
    # plan §5.1, schema in objdel.verify.gates; datasets=("riga",) for a model dropped from LIDC
"""
