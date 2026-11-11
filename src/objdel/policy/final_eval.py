"""Final evaluation of a tuned model on a whole split: K samples per item, scored by the checker. The CLI is
experiments/phase2/final_eval.py; the flags are an objdel.core.spec.EvalStage.

Samples every RL-eligible item of the view's <split> records with objdel.policy.sampler.VllmSampler (T = 1, the
checkpoint's stop ids, `--adapter` served unmerged, each request seeded by request_seed(seed, item_id)), then scores
in this process through objdel.policy.evaluate. Write order per view, so a crash in scoring keeps the paid sampling:
1. <out>/answers.jsonl and a partial <out>/final_eval.json ("partial": true, the run keys, seconds.load_and_sample);
2. after scoring, <out>/metrics.jsonl (every metric per item), <out>/items.jsonl (the per-item K-means) and the full
   <out>/final_eval.json (schema 2, "partial": false). Every file goes through a temp file.
Run keys: ckpt, view, split, K, seed, sampler, view_manifest_sha256, view_content_sha256, code_sha, checker_build.

The summary holds the K-sample mean per item, then the mean over items with three 95 % percentile bootstrap intervals
(B = 2,000, seed 0): `ci95_patient` (the gate CI, named by `gate_ci`), `ci95_nodule` and `ci95_item`. Illegal samples
score 0 on every geometric term. Keys: the METRICS of objdel.policy.evaluate (iou_tc is IoU against the doctors'
band, the reward's term), miss_rate (items whose K-mean IoU_tc is below 0.1), iou_found (the mean over the other
items), sel_* (one answer chosen per item with no doctor masks, objdel.verify.select) and noise_floor (K-answer
sampling noise from two answer halves; printed beside gains, never a rule).

Other forms, none loads a model:
- --rescore DIR: score DIR/answers.jsonl again and write the same four files to --out (refused, exit 2, when the
  view's manifest.json is not the one DIR/final_eval.json was sampled on);
- --hash-view VIEW: the view's identity record (content and target-set hashes, the last commit of each file that
  renders it);
- --resummarize DIR: the schema-2 summary of a finished run's items.jsonl, group_id joined from the LIDC store.
A frozen view (objdel.data.quality.FROZEN) runs only when its content and target-set hashes equal its entry in
reference/frozen_views.json, and --split must name a set of the view (the test view needs --split test). Both checks
run before any model load or rescoring; a failed one raises Refusal. Reads untrusted data: run with python -I.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

from objdel.core.hashing import file_sha256
from objdel.core.io import code_sha, read_json, read_jsonl
from objdel.core.io import write_atomic as _write
from objdel.core.provenance import repo_root
from objdel.data import quality as Q
from objdel.geom import engine as _gpucheck
from objdel.policy.errors import Refusal
from objdel.policy.evaluate import METRICS, evaluate, per_item, score_items, summarize
from objdel.verify import device
from objdel.verify import stats as signal

FROZEN_VIEWS = repo_root() / "reference/frozen_views.json"  # /workspace/src on the box

RUN_KEYS = ("ckpt", "view", "split", "K", "seed")
SAMPLER_LABEL = "VllmSampler; sampler not box-verified (E-G10 open)"  # written into every result


def noise_floor(results: list[dict], groups: dict[str, str], draws: int = 2000, seed: int = 0) -> dict:
    """K-sample noise from one run's answers: d_i = mean of answers 0..K/2-1 minus mean of the rest.
    item_sd = sd(d) / sqrt 2 and mean_sd = the patient-bootstrap SD of mean(d) / sqrt 2. Both are the SD of a
    difference of two independent K-answer means (a retest difference), not of one mean: one K-answer mean has
    item_sd / sqrt 2 (math-stats §2, §11 item 1). Assumes the answers of one image are independent draws. A gate
    prints mean_sd beside a gain; no gate passes or fails on it."""
    ids = [r["item_id"] for r in results]
    iou = np.array([r["samples"]["iou_tc"] for r in results], float)
    if iou.shape[1] < 2:
        raise ValueError(f"noise_floor needs K >= 2 answers per item, got {iou.shape[1]}")
    h = iou.shape[1] // 2
    d = iou[:, :h].mean(1) - iou[:, h:].mean(1)
    D = signal.draws(ids, draws, seed, key=groups.__getitem__)
    return {"item_sd": float(np.std(d, ddof=1) / np.sqrt(2)),
            "mean_sd": float(np.std([d[x].mean() for x in D], ddof=1) / np.sqrt(2))}


def limit_rows(rows: list[dict], limit: int | None) -> list[dict]:
    """--limit: at most `limit` rows, evenly spaced through the view in record order (a smoke slice that spans the
    patients, not the first patient); all rows without a limit. Deterministic: no seed."""
    if limit is None or limit >= len(rows):
        return rows
    return rows[:: -(-len(rows) // limit)][:limit]


def view_groups(view: Path, split: str) -> dict[str, str]:
    """{item_id: group_id} of a view's LIDC records (never from the prompt rows: TRL forwards their columns)."""
    from objdel.policy import rows as dataset

    return {r["item_id"]: r["group_id"] for r in dataset.records(view, "lidc", split)}


def view_manifest_sha256(view: Path) -> str:
    """SHA-256 of the view's manifest.json bytes (it changes with harmless rebuild details: versions, fingerprint)."""
    return file_sha256(Path(view) / "manifest.json")


def _view_records(view: Path) -> list[dict]:
    """Every LIDC record of the view, file by file in name order."""
    return [r for f in sorted(Path(view, "lidc").glob("*.jsonl")) for r in read_jsonl(f)]


def frozen_entry(view_name: str) -> dict | None:
    """The view's frozen_views.json entry, or None (no entry, or no file)."""
    return read_json(FROZEN_VIEWS).get(view_name) if FROZEN_VIEWS.is_file() else None


def check_view(view: Path, split: str) -> None:
    """Refusal unless `split` is a set of the view and a frozen view's hashes equal its frozen_views.json entry."""
    sets = read_json(view / "manifest.json")["sets"]
    if f"lidc/{split}" not in sets:
        raise Refusal(f"final_eval: {view} holds {sorted(sets)}, not lidc/{split} (the test view needs --split test)")
    name = view.resolve().name  # an alias (symlink) of a frozen view is that view
    if name not in Q.FROZEN:
        return
    if (e := frozen_entry(name)) is None:
        raise Refusal(f"final_eval: {view.name} is a frozen view with no frozen_views.json entry "
                      f"(looked in {FROZEN_VIEWS})")
    for key, got in (("content_sha256", Q.view_content_sha256(view)), ("target_set_sha256", Q.target_set_sha256(view))):
        if got != e[key]:
            raise Refusal(f"final_eval: {view.name} {key} {got} differs from frozen_views.json {e[key]}: "
                          "refusing to score")


def hash_view(view: Path, out: Path, commits: dict[str, str]) -> None:
    """--hash-view: the view's identity record (the pre-rental record); `commits`: the last commit of each
    file that renders it (freeze_views.render_files)."""
    sets = read_json(view / "manifest.json")["sets"]
    (target_set,) = [s["target_set_sha256"] for k, s in sets.items() if k.startswith("lidc/")]
    _write(out, json.dumps({
        "view": str(view), "items": len(_view_records(view)), "content_sha256": Q.view_content_sha256(view),
        "target_set_sha256": target_set, "generator_commit": commits["experiments/phase2/build_eval_views.py"],
        "render_files": commits, "r5_commit": None,
        "note": "r5_commit not recorded: no file under runs/r5/ holds a code commit. Equal to the R5 box view by "
                "code identity, not measured"}, indent=1))


def _lines(rows: list[dict]) -> str:
    """Rows as JSON lines."""
    return "".join(json.dumps(r) + "\n" for r in rows)


def resummarize(run: Path, out: Path, store: Path) -> None:
    """Schema-2 final_eval.json of a finished run (module docstring); exit 2 when an item has no group_id."""
    import pyarrow.parquet as pq

    doc = read_json(run / "final_eval.json")
    items = read_jsonl(run / "items.jsonl")
    table = pq.read_table(store / "lidc" / f"{doc['split']}.parquet", columns=["item_id", "group_id"])
    found = {i: g for i, g in zip(*table.to_pydict().values(), strict=True) if g}
    lost = [it["item_id"] for it in items if it["item_id"] not in found]
    if lost:
        raise Refusal(f"resummarize: {len(lost)} items have no group_id in {store}, first {lost[0]}", 2)
    s = summarize(items, found)
    res = {k: doc[k] for k in RUN_KEYS} | s | {"metrics_absent": [m for m in METRICS if m not in s],
                                               "resummarized_from": str(run)}
    out.mkdir(parents=True, exist_ok=True)
    (out / "final_eval.json").write_text(json.dumps(res, indent=1))


def load_sampler(a):
    """(VllmSampler of the stage's ckpt, built once per process, serving a.adapter unmerged when given; params_for's
    dict at K = a.k and seed a.seed with the checkpoint's stop ids)."""
    from objdel.policy import sampler as sample

    s = sample.VllmSampler(Path(a.ckpt), a.model, adapter=a.adapter and Path(a.adapter))
    stop = sorted(sample.stop_ids(Path(a.ckpt), s.processor))
    return s, sample.params_for(a.model, n=a.k, seed=a.seed, stop_token_ids=stop)


def engine_record() -> dict:
    """The view's eval_sets wall time, the scratch cap the loop passed and the engine's memory counters, read
    after its scoring (the counters were reset just before it)."""
    s, torch = _gpucheck.stats(), sys.modules.get("torch")  # a rescore imports no torch
    return {"eval_sets_s": device.eval_sets_s, "max_scratch_bytes": device.MAX_SCRATCH_BYTES,
            **{k: s[k] for k in ("arena_bytes", "table_bytes", "eval_scratch_peak", "live_device_bytes")},
            "max_memory_reserved": torch.cuda.max_memory_reserved() if torch else 0}


def start_scoring() -> None:
    """Reset the engine counters before one view's scoring."""
    _gpucheck.reset_stats()
    device.eval_sets_s = 0.0


def write_scored(view: Path, split: str, out: Path, results: list[dict], run: dict, t0: float):
    """Writes metrics.jsonl, items.jsonl and the full final_eval.json (schema 2, with the view's engine record)."""
    items, groups = [per_item(r) for r in results], view_groups(view, split)
    doc = run | summarize(items, groups) | {"noise_floor": noise_floor(results, groups), "partial": False,
                                            "engine": engine_record()}
    doc["seconds"] = run["seconds"] | {"total": time.time() - t0}
    _write(out / "metrics.jsonl", _lines(results))
    _write(out / "items.jsonl", _lines(items))
    _write(out / "final_eval.json", json.dumps(doc, indent=1, allow_nan=False))
    print("FINAL_EVAL", split, json.dumps({m: round(doc[m]["mean"], 4) for m in METRICS if m in doc}),
          "items", doc["items"], f"s {doc['seconds']['total']:.0f}")  # fmt: skip


def _identity(view: Path, split: str) -> dict:
    """The run keys that name the view and the code (module docstring, Run keys)."""
    return {"view": str(view), "split": split, "view_manifest_sha256": view_manifest_sha256(view),
            "view_content_sha256": Q.view_content_sha256(view), "code_sha": code_sha(repo_root(), required=False),
            "checker_build": _gpucheck.build()}


def rescore(src: Path, view: Path, split: str, out: Path) -> None:
    """--rescore: score src/answers.jsonl again on the view, with no model; exit 2 when the view's manifest is not
    the one src was sampled on."""
    from objdel.policy import rows as dataset
    from objdel.verify import store as serve

    f = src / "final_eval.json"
    old = read_json(f) if f.exists() else {}
    if "view_manifest_sha256" not in old:
        raise Refusal(f"rescore: {f} is missing or has no view_manifest_sha256")
    if old["view_manifest_sha256"] != view_manifest_sha256(view):
        raise Refusal(f"rescore: {view}/manifest.json is not the manifest {src} was sampled on", 2)
    t0 = time.time()
    answers = read_jsonl(src / "answers.jsonl")
    by_key = {r["key"]: r for r in dataset.prompts(view, "lidc", split)}
    rows = [by_key[d["item_id"]] for d in answers]
    run = ({k: old[k] for k in ("ckpt", "K", "seed", "sampler") if k in old} | _identity(view, split)
           | {"rescored_from": str(src), "seconds": {}})
    serve.preload(view, ["lidc"])
    _write(out / "answers.jsonl", _lines(answers))
    start_scoring()
    results, _ = score_items(rows, {d["item_id"]: d for d in answers}, view.resolve())
    write_scored(view, split, out, results, run, t0)


def run(a) -> None:
    """--rescore, or sample and score every view of the EvalStage `a` (module docstring); exit 2 on a bad flag."""
    def error(msg: str):
        raise Refusal(f"final_eval: error: {msg}", 2)

    for name in ("adapter_config.json", "adapter_model.safetensors") if a.adapter else ():
        if not (Path(a.adapter) / name).is_file():
            error(f"--adapter {a.adapter}: {name} is missing")
    views, splits = [Path(v) for v in a.view], list(a.split) or ["val"] * len(a.view)
    outs = [Path(o) for o in a.out] if len(a.out) == len(views) else [Path(a.out[0]) / v.resolve().name for v in views]
    for view, split in zip(views, splits, strict=True):
        check_view(view, split)
    if len({o.resolve() for o in outs}) != len(outs):
        error("two views write the same --out folder: give distinct --view directory names or one --out each")
    if a.rescore:
        if outs[0].resolve() == Path(a.rescore).resolve():
            error("--out must not be the --rescore folder")
        if len(views) != 1:
            error("--rescore takes one --view")
        return rescore(Path(a.rescore), views[0], splits[0], outs[0])

    import torch

    from objdel.policy import rows as dataset
    from objdel.verify import store as serve

    t0 = time.time()
    for root in dict.fromkeys(v.resolve() for v in views):
        serve.preload(root, ["lidc"])
    sampler, params = load_sampler(a)
    for view, split, out in zip(views, splits, outs, strict=True):
        rows = limit_rows(dataset.prompts(view, "lidc", split), a.limit)
        groups = view_groups(view, split)
        run = {"ckpt": str(a.ckpt), "K": a.k, "seed": a.seed, "sampler": SAMPLER_LABEL} | (
            {"adapter": str(a.adapter), "adapter_sha256": file_sha256(Path(a.adapter) / "adapter_model.safetensors")}
            if a.adapter else {}) | _identity(view, split) | ({"limit": a.limit} if a.limit else {})

        def keep(answers, run=run, out=out, view=view, groups=groups, rows=rows, t0=t0):
            """answers.jsonl and the partial final_eval.json before scoring (a crash keeps the sampling)."""
            lines = [{"item_id": d["item_id"], "group_id": groups[d["item_id"]], "W": r["W"], "H": r["H"],
                      **{k: d[k] for k in ("samples", "truncated", "tokens")}}
                     for r, d in zip(rows, answers, strict=True)]
            run["seconds"] = {"load_and_sample": time.time() - t0}
            _write(out / "answers.jsonl", _lines(lines))
            _write(out / "final_eval.json", json.dumps(run | {"partial": True}, indent=1))
            start_scoring()

        torch.manual_seed(a.seed)
        results = evaluate(sampler, rows, params, view.resolve(), keep)["results"]
        write_scored(view, split, out, results, run, t0)
        t0 = time.time()
