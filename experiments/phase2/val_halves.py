"""Write reference/val_halves.json (Implementation.md B1.0): the val halves of objdel.data.splits.val_half_table
on the frozen val shifted view, plus the 128 held-out item ids that chose R5 (B.6.2 Limitations).

    .venv/bin/python -I -m experiments.phase2.val_halves --view data/eval_views/lidc-val-shifted/lidc/val.jsonl \\
        --heldout runs/rl/eval/final/heldout.json --out reference/val_halves.json

Test surface: kept by the CLEAN16 decision (2026-11-10); only tests import it, so a change here is checked there.
"""
import argparse
import json
from pathlib import Path

from objdel.core.io import read_json
from objdel.data.splits import val_half_table


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--view", type=Path, required=True, help="the val view's lidc/val.jsonl")
    ap.add_argument("--heldout", type=Path, required=True, help="heldout.json of the run that chose R5")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    records = [json.loads(x) for x in a.view.read_text().splitlines()]
    chose = sorted(read_json(a.heldout)["per_item_iou"])
    if missing := set(chose) - {r["item_id"] for r in records}:
        raise SystemExit(f"{len(missing)} held-out items are not in the val view, e.g. {min(missing)}")
    table = val_half_table(records) | {"r5_selection_items": chose}
    a.out.write_text(json.dumps(table, indent=1) + "\n")


if __name__ == "__main__":
    main()
