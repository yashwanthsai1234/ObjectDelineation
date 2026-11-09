"""The B data gates DG1-DG18 (21 gates, objdel.data.quality) of a rows root: prints the gate table, exit 1 on a hard
failure (report gates DG9b, DG14, DG14b, DG16 never change the exit code). Implementation.md B1.3.

    .venv/bin/python -I -m experiments.phase2.audit_data --rows data/views/lidc-train-b1 [--store data/store]
        [--val-view data/eval_views/lidc-val-shifted] [--json out.json]

The rows root's manifest.schedule gives the pool (draws, seed, excluded_ids) and the schedule arguments; the target
is the train shifted law of that pool. The val view feeds only the DG16 report (val-shape, val-select, whole val).
On the Mac the manifest's exclusions_sha256 must equal the SHA-256 of the current exclusions file (DG1b). DG14b needs
objdel.eval.ceiling (A2). Reads untrusted data: run with python -I.
"""
import argparse
import json
from dataclasses import asdict
from pathlib import Path

from objdel.data import quality as Q
from objdel.data import splits, store, windows

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "docs" / "reference"


def tokenizer(model: str):
    """(tokenizer, max completion tokens) of a train config model, from the local HF cache."""
    from transformers import AutoTokenizer

    from objdel.train.configs import MAX_COMPLETION, MODELS

    return AutoTokenizer.from_pretrained(MODELS[model]), MAX_COMPLETION[model]


def _json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.is_file() else {}


def audit(rows, store_root=ROOT / "data/store", eval_views=ROOT / "data/eval_views", val_view=None,
          exclusions=REF / "train_exclusions.json", review=REF / "neardup_review.json",
          frozen=REF / "frozen_views.json", series_z=ROOT / "runs/b1/series_z.json", model="qwen3-vl-2b",
          out=ROOT / "runs/b1/audit", plan=None) -> dict:
    """The report that quality.gates reads, for the rows root `rows` (module docstring). `plan`: the schedule dict
    (build_rl_rows passes it before the manifest exists); default the manifest's `schedule`."""
    rows, store_root, eval_views, exclusions = Path(rows), Path(store_root), Path(eval_views), Path(exclusions)
    plan = plan or json.loads((rows / "manifest.json").read_text())["schedule"]
    ex = _json(exclusions)
    pool = windows.pool(store.load("lidc", "train", store_root), plan["draws"], plan["seed"],
                        frozenset(plan["excluded_ids"]))
    by_id = {w.window_id: w for w in pool.windows}
    recs = Q.records(rows)
    scheduled = [by_id[r["item_id"]] for r in recs if r["item_id"] in by_id]
    val = Q.view_windows(val_view or eval_views / "lidc-val-shifted", store_root)
    halves = {f"val-{h}": [w for w in val if splits.val_half(w.group_id) == h] for h in ("shape", "select")}
    sched = Q.schedule_report(scheduled, pool, halves | {"val": val}, **{k: plan[k] for k in Q.SCHEDULE_ARGS})
    wn = Q.window_checks(rows, store_root)
    tok, limit = tokenizer(model)
    z = Path(series_z)
    return {
        "leakage": Q.leakage(store_root, _json(Path(review)), ex, out),
        "rows": Q.rows_leakage(rows, store_root, Q.view_image_hashes(eval_views), frozenset(ex))
        | {"exclusions_stale": plan["exclusions_sha256"] != Q.file_sha256(exclusions)},
        "windows": wn,
        "schedule": sched | {"unknown": len(recs) - len(scheduled)},
        "tokens": Q.tokens(rows, tok, limit),
        "agreement": Q.agreement(rows, wn["low_contrast"]),
        "frozen": Q.frozen(eval_views, _json(Path(frozen))),
        "z": Q.dg17(json.loads(z.read_text()), ex, Q.file_sha256(store_root / "lidc/train.parquet"))
        if z.is_file() else [f"no contour-to-slice report {z}: run check_series_z"],
    }


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
    report = audit(a.rows, a.store, a.eval_views, a.val_view, a.exclusions, a.review, a.frozen, a.series_z, a.model,
                   a.out)
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
