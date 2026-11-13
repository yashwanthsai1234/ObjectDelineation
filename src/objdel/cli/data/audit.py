"""The B data gates DG1-DG18 (21 gates, objdel.data.quality) of a rows root: prints the gate table, exit 1 on a hard
failure (report gates DG9b, DG14, DG14b, DG16 never change the exit code). Implementation.md B1.3.

    .venv/bin/python -I -m objdel.cli.data.audit --rows data/views/lidc-train-b1-780 [--store data/store]
        [--val-view data/eval_views/lidc-val-shifted] [--json out.json]

The rows root's manifest.schedule gives the pool (draws, seed, excluded_ids) and the schedule arguments; the target
is the train shifted law of that pool. The val view feeds only the DG16 report (val-shape, val-select, whole val).
On the Mac the manifest's exclusions_sha256 must equal the SHA-256 of the current exclusions file (DG1b). DG14b needs
objdel.verify.ceiling (A2). Reads untrusted data: run with python -I.
"""
import argparse
import json
from dataclasses import asdict
from pathlib import Path

from objdel.core.provenance import repo_root
from objdel.data import quality as Q

ROOT = repo_root()
REF = ROOT / "reference"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--rows", type=Path, required=True)
    ap.add_argument("--store", type=Path, default=ROOT / "data/store")
    ap.add_argument("--eval-views", type=Path, default=ROOT / "data/eval_views", help="val and test views (DG1b, DG15)")
    ap.add_argument("--val-view", type=Path, default=None, help="DG16 only (default <eval-views>/lidc-val-shifted)")
    ap.add_argument("--exclusions", type=Path, default=REF / "train_exclusions.json")
    ap.add_argument("--review", type=Path, default=REF / "neardup_review.json")
    ap.add_argument("--frozen", type=Path, default=REF / "frozen_views.json")
    ap.add_argument("--series-z", type=Path, default=ROOT / "runs/b1/series_z.json")
    ap.add_argument("--model", default="qwen3-vl-2b", help="tokenizer and token limit of DG12")
    ap.add_argument("--out", type=Path, default=ROOT / "runs/b1/audit", help="near-duplicate PNGs (DG3)")
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)
    from objdel.policy import models

    report = Q.audit(a.rows, *models.tokenizer(a.model), a.store, a.eval_views, a.val_view, a.exclusions, a.review,
                     a.frozen, a.series_z, a.out)
    gs = Q.gates(report)
    for g in gs:
        state = "report" if not g.hard else "ok" if g.ok else "FAIL"
        print(f"{g.name:6s} {state:6s} {g.value:10.4f}  {g.limit}")
    bad = Q.failures(gs)
    print(f"AUDIT {'FAILED: ' + ', '.join(g.name for g in bad) if bad else 'PASSED'} ({a.rows})")
    if a.json:
        a.json.write_text(json.dumps({"gates": [asdict(g) for g in gs], "report": report}, indent=1) + "\n")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
