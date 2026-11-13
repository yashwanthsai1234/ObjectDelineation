"""Copy the headline test and validation metrics of final_eval.json files into reference/headline.json.

    python -I -m objdel.cli.eval.headline --run r5 --model-ref svenk029/objdel-sft@b1-rows
    python -I -m objdel.cli.eval.headline --test T/final_eval.json --val V/final_eval.json --model-ref REF

--run ID reads runs/ID/eval/eval_{test,val}/final_eval.json. Values are copied at full precision, never typed.
The run id is ID, or the folder under runs/ of --test (else empty).
Schema 1 (no `schema` key) holds `ci95` (item bootstrap); schema 2 names its gate CI in `gate_ci` (patient).
Refuses (exit 1, writes nothing) on a missing file, a wrong split, no view, mixed schemas or a null gate CI.
"""

import argparse
import json
import sys
from pathlib import Path

from objdel.core.io import read_json

METRICS = ("iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote", "legal")


def section(path: Path, split: str) -> tuple[dict, str, str]:
    """One final_eval.json -> (headline section, ci_unit, schema key); ValueError when it may not be a headline."""
    if not path.is_file():
        raise ValueError(f"{path} missing")
    f = read_json(path)
    if f.get("split") != split:
        raise ValueError(f"{path} has split {f.get('split')!r}, expected {split!r}")
    if not f.get("view"):
        raise ValueError(f"{path} has no view")
    key = f.get("gate_ci", "ci95") if f.get("schema") == 2 else "ci95"
    metrics = {}
    for m in METRICS:
        if (ci := f[m].get(key)) is None:
            raise ValueError(f"{path}: {m} has no {key}")
        metrics[m] = {"mean": f[m]["mean"], "ci95": ci}
    unit = "item" if key == "ci95" else key.removeprefix("ci95_")
    sec = {"source": str(path), "view": Path(f["view"]).name, "items": f["items"], "metrics": metrics}
    return sec | {"_run": {k: f[k] for k in ("ckpt", "seed", "K")}}, unit, f.get("schema", 1)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run")
    ap.add_argument("--test", type=Path)
    ap.add_argument("--val", type=Path)
    ap.add_argument("--model-ref", required=True)
    ap.add_argument("--out", type=Path, default=Path("reference/headline.json"))
    a = ap.parse_args(argv)
    if (bool(a.run), bool(a.test) + bool(a.val)) not in ((True, 0), (False, 2)):
        print("headline: give --run ID, or both --test and --val", file=sys.stderr)
        return 2
    if a.run:
        a.test, a.val = (Path(f"runs/{a.run}/eval/eval_{s}/final_eval.json") for s in ("test", "val"))
    try:
        (test, unit, schema), (val, unit_v, schema_v) = section(a.test, "test"), section(a.val, "val")
        if (unit, schema) != (unit_v, schema_v):
            raise ValueError("the test and val files differ in schema or gate CI")
    except ValueError as e:
        print(f"headline: {e}", file=sys.stderr)
        return 1
    run = test.pop("_run")
    val.pop("_run")
    run_id = a.run or (a.test.parts[1] if a.test.parts[:1] == ("runs",) else "")  # runs/<id>/...
    doc = {"run": run_id, "model_ref": a.model_ref, **run, "ci_unit": unit, "test": test,
           "val": val}  # fmt: skip
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(doc, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
