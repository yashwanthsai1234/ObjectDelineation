"""Evaluate sampled completions against the raters of the data store; write a JSON and a Markdown report.

    .venv/bin/python -I -m experiments.phase1.eval_run --completions runs/sft_val.jsonl --out data/eval/sft \
        [--store data/store] [--split val] [--baseline data/eval/sft0/report.json] [--datasets riga] [--procs 3]

Input: JSONL, one line per item, {"item_id": "<dataset>/...", "samples": [text, ...], "truncated": [bool, ...]}
("truncated" optional, default all False). Every item is framed exactly as for training (objdel.targets.frame) and
scored by objdel.eval.metrics on the frame (W = H = S) with the dataset's schema. The report holds the per-dataset
summary, the hack audit (against the baseline report's audit, normally the SFT checkpoint, for the component alarm)
and the Phase 2 entry gates, decided per dataset (--datasets: the datasets the model is judged on, default those in
the completions; e.g. riga alone for a model the precision probe dropped from LIDC); per-item results go to
items.jsonl. Reads untrusted data: run with python -I.

References:
- Metrics, audit and gates: objdel.eval (plan §4.7, §7, §9; experiments/phase1/impl/M8_eval.md).
- Training frame: objdel.targets.frame (EXECUTION.md §3, M5).
- Process pool: https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool ; JSON without NaN
  (allow_nan=False): https://docs.python.org/3/library/json.html#json.dump
- python -I (isolated mode) for untrusted data: https://docs.python.org/3/using/cmdline.html#cmdoption-I
"""
import argparse
import json
from dataclasses import replace
from multiprocessing import Pool
from pathlib import Path

from objdel.data import store
from objdel.eval import gates, hack_audit, metrics, summarize
from objdel.geometry import SCHEMAS
from objdel.targets import frame


def load_frames(root, split: str, ids=None, every: dict | None = None) -> dict:
    """{item_id: {"dataset", "masks", "W", "H"}}: the training frame of every requested item (all when ids is None;
    every = {dataset: n} keeps every n-th item). Streams the store and keeps only the framed masks (a RIGA item's
    full-size rater stacks are ~40 MB, its frame ~3 MB). KeyError if a requested id is missing."""
    want = None if ids is None else set(ids)
    datasets = store.DATASETS if want is None else sorted({i.split("/")[0] for i in want})
    out = {}
    for ds in datasets:
        for i, it in enumerate(store.load(ds, split, root)):
            if (want is None and i % (every or {}).get(ds, 1) == 0) or (want is not None and it.item_id in want):
                fr = frame(replace(it, image=None))  # the image is not needed
                out[it.item_id] = {"dataset": ds, "masks": fr.masks, "W": fr.W, "H": fr.H}
    missing = (want or set()) - out.keys()
    if missing:
        raise KeyError(f"{len(missing)} item ids not in the {split} store, e.g. {sorted(missing)[0]}")
    return out


def _job(args) -> dict:
    row, fr = args
    return {"item_id": row["item_id"], "dataset": fr["dataset"],
            **metrics(row["samples"], fr["masks"], fr["W"], fr["H"], SCHEMAS[fr["dataset"]], row.get("truncated"))}


def evaluate(rows, frames: dict, procs: int = 1) -> list[dict]:
    """metrics of every completions row (in order) on its item's frame."""
    jobs = [(r, frames[r["item_id"]]) for r in rows]
    if procs <= 1:
        return [_job(j) for j in jobs]
    with Pool(procs) as pool:
        return pool.map(_job, jobs, chunksize=4)


def build_report(results: list[dict], baseline: dict | None = None, datasets=None) -> dict:
    report = {"summary": summarize(results), "audit": hack_audit(results, baseline and baseline["audit"])}
    return report | {"gates": gates.phase2_entry(report, datasets)}


def _f(v, pct=False):
    if v is None:
        return "–"
    return f"{v:.1%}" if pct else f"{v:.3f}"


def markdown(report: dict, title: str) -> str:
    """The report's key numbers as Markdown tables."""
    L = [f"# {title}", ""]
    L += ["| dataset | items | samples | format | legal | topology-any | topo-any (unanimous) | pass | pass@k (k: value) |",
          "|---|---|---|---|---|---|---|---|---|"]
    for ds, d in report["summary"]["datasets"].items():
        pk = ", ".join(f"{k}: {v:.3f}" for k, v in d["pass_at_k"].items())
        L.append(f"| {ds} | {d['items']} | {d['samples']} | {_f(d['format'], 1)} | {_f(d['legal'], 1)} | "
                 f"{_f(d['topo_any'], 1)} | {_f(d['topo_any_unanimous'], 1)} | {_f(d['pass_rate'], 1)} | {pk} |")
    L += ["", "| dataset | TV | JS | coverage | disagreement items | TV (dis.) | JS (dis.) | min-mode (dis.) | "
              "group | duplicates | groups with σ > floor |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for ds, d in report["summary"]["datasets"].items():
        q = d["disagreement"]
        L.append(f"| {ds} | {_f(d['tv'])} | {_f(d['js'])} | {_f(d['coverage'])} | {q['items']} | {_f(q['tv'])} | "
                 f"{_f(q['js'])} | {_f(q['min_mode'])} | {d['group']} | {_f(d['dup_rate'], 1)} | "
                 f"{_f(d['signal_share'], 1)} |")
    L += ["", "| structure | IoU_tc median | IoU_tc mean | F5 median | IoU vote | Dice vote | D_max | g_area | "
              "soft Dice | GED |", "|---|---|---|---|---|---|---|---|---|---|"]
    for ds, d in report["summary"]["datasets"].items():
        for n, o in d["objects"].items():
            L.append(f"| {ds}/{n} | {_f(o['iou_tc_median'])} | {_f(o['iou_tc_mean'])} | {_f(o['f5_median'])} | "
                     f"{_f(o['iou_vote_mean'])} | {_f(o['dice_vote_mean'])} | {_f(o['d_max_mean'])} | "
                     f"{_f(o['g_area_mean'])} | {_f(o['soft_dice'])} | {_f(o['ged'])} |")
    a = report["audit"]
    L += ["", "| structure | outputs | median A/A_max | median A/A_min | median \\|M\\|/\\|E\\| | > 1.3 A_max | "
              "above envelope | below envelope | extra components | median components |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for s, v in a["structures"].items():
        if v:
            L.append(f"| {s} | {v['outputs']} | {_f(v['a_amax_median'])} | {_f(v['a_amin_median'])} | "
                     f"{_f(v['m_e_median'])} | {_f(v['above_tail'], 1)} | {_f(v['inflated'], 1)} | "
                     f"{_f(v['deflated'], 1)} | {_f(v['fragmented'], 1)} | {v['components_median']:g} |")
    L += ["", "| dataset | RHR | NRFR | rewarded outputs | newly rewarded failures |", "|---|---|---|---|---|"]
    for ds, v in a["datasets"].items():
        L.append(f"| {ds} | {_f(v['rhr'], 1)} | {_f(v['nrfr'], 1)} | {v['rewarded']} | {v['nrf']} |")
    state = lambda x: "not measured" if x["fired"] is None else ("**FIRED**" if x["fired"] else "ok")  # noqa: E731
    L += ["", "Alarms (plan §9): " + ", ".join(f"{k} {state(v)}" for k, v in a["alarms"].items()) + ".",
          "Audit flags: " + ", ".join(f"{k} {state(v)}" for k, v in a["flags"].items()) + ".", "",
          "| Phase 2 gate | requirement | values | pass |", "|---|---|---|---|"]
    for k, g in report["gates"]["gates"].items():
        vals = ", ".join(f"{n} {_f(v)}" for n, v in g["values"].items())
        L.append(f"| {k} | {g['requirement']} | {vals} | {'PASS' if g['pass'] else 'FAIL'} |")
    per = ", ".join(f"{ds} {'PASS' if ok else 'FAIL'}" for ds, ok in report["gates"]["datasets"].items())
    verdict = "PASS" if report["gates"]["pass"] else "FAIL"
    L.append(f"\nPhase 2 entry gates measurable here: **{verdict}** ({per or 'no dataset'}).")
    return "\n".join(L) + "\n"


def write(report: dict, results: list[dict], out: Path, title: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(report, indent=1, allow_nan=False) + "\n")
    (out / "report.md").write_text(markdown(report, title))
    with open(out / "items.jsonl", "w") as f:
        for r in results:
            f.write(json.dumps(r, allow_nan=False) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--completions", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--store", type=Path, default=Path("data/store"))
    ap.add_argument("--split", default="val")
    ap.add_argument("--baseline", type=Path, help="report.json of the SFT checkpoint (component alarm)")
    ap.add_argument("--datasets", nargs="+", help="datasets the model is judged on (default: those evaluated)")
    ap.add_argument("--procs", type=int, default=3)
    a = ap.parse_args()
    rows = [json.loads(line) for line in a.completions.read_text().splitlines() if line.strip()]
    results = evaluate(rows, load_frames(a.store, a.split, [r["item_id"] for r in rows]), a.procs)
    baseline = json.loads(a.baseline.read_text()) if a.baseline else None
    report = build_report(results, baseline, a.datasets)
    write(report, results, a.out, f"Evaluation of {a.completions.name} ({a.split})")
    print(markdown(report, a.completions.name))


if __name__ == "__main__":
    main()
