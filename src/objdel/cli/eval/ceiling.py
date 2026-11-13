"""Fair ceiling of a view: each doctor and the model against the same k - 1 doctors (objdel.verify.ceiling, plan A2).

    python -I -m objdel.cli.eval.ceiling --view data/eval_views/lidc-val-shifted --split val \
        [--answers runs/a3/eval_val/answers.jsonl] --out runs/a3/ceiling_val

Every RL item of the view's LIDC <split> records, on the view's caches (checked against the records' digests).
Without --answers only the doctor side (model keys null). With it, the K answers of each item (final_eval's
answers.jsonl, joined by item_id; ValueError naming the missing items, before any output) are scored against every
leave-one-out cache and against all doctors. Writes <out>/ceiling.json (ceiling.summarize with patient draws from the
records' group_id, plus view, split, answers and seconds) and <out>/ceiling_items.jsonl (one row per item).
Items are scored in this process, one after another (D7: no forked worker, one device context). Reads untrusted
data: run with python -I.
"""

import argparse
import json
import time
from pathlib import Path

from objdel.core.io import write_jsonl
from objdel.policy import rows as dataset
from objdel.verify import ceiling
from objdel.verify import store as serve


def _item(job) -> dict:
    """One item's ceiling row on the preloaded caches."""
    item_id, root, texts, truncated, W, H = job
    return ceiling.item(serve.CacheRef(item_id, root)["nodule"], texts, truncated, W, H)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--view", type=Path, required=True, help="a view directory (data/prepared layout)")
    ap.add_argument("--split", default="val")
    ap.add_argument("--answers", type=Path, help="final_eval's answers.jsonl of this view")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)

    t0 = time.time()
    recs = [r for r in dataset.records(a.view, "lidc", a.split) if r["rl"]]
    ans = {}
    if a.answers:
        with open(a.answers) as f:
            ans = {d["item_id"]: (d["samples"], d["truncated"]) for d in map(json.loads, f)}
        lost = [r["item_id"] for r in recs if r["item_id"] not in ans]
        if lost:
            raise ValueError(f"{len(lost)} items of the view have no answers in {a.answers}: {lost}")
    root = str(a.view.resolve())
    serve.preload(root, ["lidc"])
    jobs = [(r["item_id"], root, *ans.get(r["item_id"], (None, None)), r["W"], r["H"]) for r in recs]
    rows = [{"item_id": r["item_id"], "group_id": r["group_id"]} | _item(j) for r, j in zip(recs, jobs, strict=True)]
    doc = {"view": str(a.view), "split": a.split, "answers": a.answers and str(a.answers),
           **ceiling.summarize(rows, {r["item_id"]: r["group_id"] for r in recs}), "seconds": time.time() - t0}
    a.out.mkdir(parents=True, exist_ok=True)
    write_jsonl(a.out / "ceiling_items.jsonl", rows, sort_keys=False)
    (a.out / "ceiling.json").write_text(json.dumps(doc, indent=1))
    s = doc["all"]
    print("CEILING", a.split, "items", s["n"], "doctor", round(s["doctor"][0], 4), "legacy item",
          round(s["legacy_loo_dmax_item"][0], 4), "pooled", round(s["legacy_loo_dmax_pooled"], 4),
          f"s {doc['seconds']:.0f}")


if __name__ == "__main__":
    main()
