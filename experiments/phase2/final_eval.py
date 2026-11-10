"""Final evaluation of a tuned model on a whole split: K samples per item, scored by the checker.

    python -I -m experiments.phase2.final_eval --ckpt DIR --view data/eval_views/lidc-val-shifted --split val \
        --out runs/eval/val [-k 8] [--seed 11] [--adapter LORA_DIR]
    python -I -m experiments.phase2.final_eval --ckpt DIR --view V1 --split val --out O1 \
        --view V2 --split val --out O2          # several views, one model load (one --out: O/<view name>/)
    python -I -m experiments.phase2.final_eval --ckpt DIR --view V --split val --out O \
        --limit 32                              # a GPU smoke: 32 evenly spaced items, not a result
    python -I -m experiments.phase2.final_eval --rescore runs/a3/eval_val --view V --split val --out O
    python -I -m experiments.phase2.final_eval --hash-view data/eval_views/lidc-test-shifted \
        --out runs/a3/pre_rental.json
    python -I -m experiments.phase2.final_eval --resummarize runs/r5/eval/eval_test --out runs/r5/eval_a1/eval_test \
        [--store data/store]
    python -I -m experiments.phase2.final_eval --spec specs/b1_rows.yaml --stage eval [--rd RD] [--data DATA]

The model forms read their flags as an objdel.core.spec.EvalStage (cli_stage: kebab flags, K as -k, or a spec stage);
--hash-view and --resummarize load no model and are read before it.

Every RL-eligible item of the view's <split> records (the G-probe uses a 128-item subset of the val view; this uses
all of them, and runs on the test split's view built by build_eval_views --split test). Samples with
objdel.train.sample.VllmSampler, no trainer (E4: vllm_native image only; one engine per process, built from --ckpt,
serving --adapter unmerged when given; params_for's sampling at T = 1 with the checkpoint's stop ids; torch seeded
before each view, each request seeded request_seed(seed, item_id)) and scores through objdel.train.evaluate.evaluate
(score_items in this process: D7, one CUDA context, after sampling). Until V4 passes every result carries
"sampler": SAMPLER_LABEL (E-G10 open). Write order per view, so a crash in scoring keeps the paid sampling:
1. <out>/answers.jsonl (item_id, group_id, W, H, samples, truncated, tokens per item) and a partial
   <out>/final_eval.json ("partial": true, the run keys, seconds.load_and_sample);
2. after scoring, <out>/metrics.jsonl (the whole metrics() result per item), <out>/items.jsonl (per-item metrics) and
   the full <out>/final_eval.json (summary, schema 2, "partial": false). Every file goes through a temp file.
Run keys: ckpt, view, split, K, seed, sampler, view_manifest_sha256 (manifest.json bytes), view_content_sha256 (B1.6's
content hash, view_content_sha256), code_sha (<repo>/CODE_SHA, written by vast sync_up on the box; else git HEAD;
else null), checker_build ("cuda" or "cpu").
Per item the K-sample mean; then the mean over items with three 95 % percentile bootstrap intervals (B = 2,000,
seed 0): `ci95_patient` resamples patients (`group_id`; the gate CI, named by `gate_ci`), `ci95_nodule` resamples
nodules, `ci95_item` single items (the schema-1 `ci95`). Illegal samples score 0 on every geometric term (the
reward's rule).
- iou_tc: IoU against the doctors' band (the reward's term); f5: the boundary score (calibrated so a new doctor
  scores 0.95); precision / recall / dice: pixels against the doctors' strict-majority outline;
- nsd, boundary_iou (best doctor) and boundary_iou_vote (the majority outline); soft_dice, hm_iou (HM-IoU8) and ged
  (lower is better) per item: objdel.eval.metrics defines them;
- miss_rate: the share of items whose K-mean IoU_tc is below 0.1 (not the best of K); iou_found: the mean IoU_tc of
  the other items, its CIs drawn on that subset;
- sel_*: one answer chosen per item with no doctor masks (objdel.eval.select; `per_item`), and `paired`: the per-item
  differences medoid - mean of K and vote - mean of K with their patient and nodule CIs. `gate_a` decides Gate A
  from the val patient CI. Run gates keep sel_mean_of_k (= iou_tc);
- noise_floor: K-answer sampling noise from the two answer halves (`noise_floor`); printed beside gains, never a rule.
--rescore DIR: score DIR/answers.jsonl again with no model and write the same four files to --out (exit 2, no output,
when the view's manifest.json is not the one DIR/final_eval.json was sampled on). --hash-view VIEW: the view's
identity record (content and target-set hashes, the last commit of each file that renders it) and exit.
`tost` and `view_identity` are A-G3's paired checks against an older run (tests/test_a_gpu.py G01).
--resummarize DIR: the schema-2 summary of a finished run's items.jsonl (no model, no scoring), with group_id joined
from the LIDC store of the run's split; written to --out, never into DIR; exit 2 when an item has no group_id.
A frozen view (objdel.data.quality.FROZEN) runs only when its content hash and target-set hash equal its entry in
docs/reference/frozen_views.json (B1.6; shipped to the box by vast.policy.code_files), and --split must name a set of
the view (the test view is read only with --split test). Both checks run before any model load or rescoring.
Reads untrusted data: run with python -I.
"""

import argparse
import contextlib
import json
import os
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np

from experiments.phase2.freeze_views import render_files
from objdel.data import quality as Q
from objdel.eval import device, signal
from objdel.geom import engine as _gpucheck
from objdel.hashing import file_sha256
from objdel.io import git, read_jsonl
from objdel.io import write_atomic as _write
from objdel.train.evaluate import METRICS, evaluate, per_item, score_items, summarize

ROOT = Path(__file__).resolve().parents[2]
FROZEN_VIEWS = ROOT / "docs/reference/frozen_views.json"  # /workspace/src on the box

RUN_KEYS = ("ckpt", "view", "split", "K", "seed")
TOPOLOGY_PERIOD_S = 5.0  # D-G17 (a): the sampler thread's period during scoring
GATE_A_GAIN = 0.015  # Gate A: the val mean gain over sel_mean_of_k a method needs (Direction.md §4)
SAMPLER_LABEL = "VllmSampler; sampler not box-verified (E-G10 open)"  # until V4's deterministic check passes


def gate_a(val: dict, test: dict | None) -> dict:
    """Gate A (A4): the medoid becomes the default answer choice when its val mean gain over sel_mean_of_k is
    >= GATE_A_GAIN and the lower end of its paired gate CI (patients) is > 0. The vote passing the same rule is
    reported, but it is no default before Phase D ships the CUDA polygonizer. `test_same_sign` (the medoid's test
    difference has the val sign; None without a test summary) is a direction line only: test never chooses."""
    def passes(arm: str) -> bool:
        p = val["paired"][arm]
        if p is None or val.get("patients", 2) in (0, 1):  # one patient's point interval proves nothing
            return False
        ci = p[val["gate_ci"]]
        if ci is None:
            raise ValueError(f"gate_a: {arm}.{val['gate_ci']} is null (a summary built without groups)")
        return round(p["mean"], signal.DIGITS) >= GATE_A_GAIN and ci[0] > 0

    medoid, vote = passes("medoid_minus_mean"), passes("vote_minus_mean")
    reason = "medoid passes on val" if medoid else "medoid does not pass on val"
    if vote:
        reason += "; the vote passes on val but needs the Phase D CUDA polygonizer to be an answer"
    same = None if test is None else bool(
        np.sign(test["paired"]["medoid_minus_mean"]["mean"]) == np.sign(val["paired"]["medoid_minus_mean"]["mean"]))
    return {"method": "medoid" if medoid else "mean_of_k", "change_default": medoid, "test_same_sign": same,
            "reason": reason}


def noise_floor(results: list[dict], groups: dict[str, str], draws: int = 2000, seed: int = 0) -> dict:
    """K-sample noise from one run's answers (A3 step 7): d_i = mean of answers 0..K/2-1 minus mean of the rest.
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


def _paired(new: dict[str, float], old: dict[str, float], ids: list[str]) -> np.ndarray:
    """The old run's values in `ids` order; ValueError unless both runs hold the same items."""
    lost = [i for i in ids if i not in old] + [i for i in old if i not in new]
    if lost:
        raise ValueError(f"{len(lost)} items are not in both runs, first {lost[0]}")
    return np.array([old[i] for i in ids], float)


def tost(new: dict[str, float], old: dict[str, float], groups: dict[str, str], bound: float = 0.01,
         draws: int = 2000, seed: int = 0) -> dict:
    """A-G3 equivalence of two runs' per-item values ({item_id: value}): the 90 % patient-bootstrap CI of the mean
    paired difference new - old lies strictly inside +-bound. A wider CI is "not shown", never a pass."""
    ids = list(new)
    d = np.array([new[i] for i in ids]) - _paired(new, old, ids)
    lo, hi = np.quantile([d[x].mean() for x in signal.draws(ids, draws, seed, key=groups.__getitem__)], [.05, .95])
    return {"mean": float(d.mean()), "ci90_patient": [float(lo), float(hi)], "bound": bound, "items": len(ids),
            "patients": len({groups[i] for i in ids}), "pass": bool(-bound < lo and hi < bound)}


def view_identity(results: list[dict], old: dict[str, float], groups: dict[str, str], draws: int = 2000,
                  seed: int = 0) -> dict:
    """A-G3 view identity: is another run's per-item IoU_tc ({item_id: K-mean}) as correlated with this run's K-mean
    as a second K-answer sampling of the same view would be? r_obs = Pearson(this K-mean, old); r_sh = 2 r_h / (1 + r_h)
    (Spearman-Brown, the K-answer reliability), r_h = Pearson of the two half means. The rule bootstraps r_obs - r_sh
    on the same patient draws and fails only when the upper end of its 95 % CI is < 0 (math-stats §11 item 7: the
    plan's "r_obs >= lower end of r_sh's CI" ignores r_obs's own error). A failure means another view."""
    ids = [r["item_id"] for r in results]
    iou = np.array([r["samples"]["iou_tc"] for r in results], float)
    h = iou.shape[1] // 2
    a, b, new, o = iou[:, :h].mean(1), iou[:, h:].mean(1), iou.mean(1), _paired(dict.fromkeys(ids), old, ids)

    def stats(x):
        r_h = np.corrcoef(a[x], b[x])[0, 1]
        return np.corrcoef(new[x], o[x])[0, 1], 2 * r_h / (1 + r_h)

    r_obs, r_sh = stats(np.arange(len(ids)))
    reps = np.array([stats(x) for x in signal.draws(ids, draws, seed, key=groups.__getitem__)])
    lo, hi = np.nanquantile(reps[:, 0] - reps[:, 1], [.025, .975])
    ok = bool(hi >= 0)
    msg = (f"view identity {'holds' if ok else 'fails'}: r_obs {r_obs:.4f}, r_sh {r_sh:.4f}, 95 % patient CI of "
           f"r_obs - r_sh [{lo:.4f}, {hi:.4f}] (fails when its upper end < 0)")
    return {"r_obs": float(r_obs), "r_sh": float(r_sh), "diff_ci95_patient": [float(lo), float(hi)],
            "r_sh_ci95_patient": np.nanquantile(reps[:, 1], [.025, .975]).tolist(), "pass": ok, "message": msg}


def limit_rows(rows: list[dict], limit: int | None) -> list[dict]:
    """--limit: at most `limit` rows, evenly spaced through the view in record order (a smoke slice that spans the
    patients, not the first patient); all rows without a limit. Deterministic: no seed."""
    if limit is None or limit >= len(rows):
        return rows
    return rows[:: -(-len(rows) // limit)][:limit]


def view_groups(view: Path, split: str) -> dict[str, str]:
    """{item_id: group_id} of a view's LIDC records (never from the prompt rows: TRL forwards their columns)."""
    from objdel.train import dataset

    return {r["item_id"]: r["group_id"] for r in dataset.records(view, "lidc", split)}


def view_manifest_sha256(view: Path) -> str:
    """SHA-256 of the view's manifest.json bytes (it changes with harmless rebuild details: versions, fingerprint)."""
    return file_sha256(Path(view) / "manifest.json")


def _view_records(view: Path) -> list[dict]:
    """Every LIDC record of the view, file by file in name order."""
    return [json.loads(line) for f in sorted(Path(view, "lidc").glob("*.jsonl")) for line in f.read_text().splitlines()]


def frozen_entry(view_name: str) -> dict | None:
    """The view's frozen_views.json entry, or None (no entry, or no file)."""
    return json.loads(FROZEN_VIEWS.read_text()).get(view_name) if FROZEN_VIEWS.is_file() else None


def check_view(view: Path, split: str) -> None:
    """Exit unless `split` is a set of the view and a frozen view's hashes equal its frozen_views.json entry."""
    sets = json.loads((view / "manifest.json").read_text())["sets"]
    if f"lidc/{split}" not in sets:
        sys.exit(f"final_eval: {view} holds {sorted(sets)}, not lidc/{split} (the test view needs --split test)")
    name = view.resolve().name  # an alias (symlink) of a frozen view is that view
    if name not in Q.FROZEN:
        return
    if (e := frozen_entry(name)) is None:
        sys.exit(f"final_eval: {view.name} is a frozen view with no frozen_views.json entry (looked in {FROZEN_VIEWS})")
    for key, got in (("content_sha256", Q.view_content_sha256(view)), ("target_set_sha256", Q.target_set_sha256(view))):
        if got != e[key]:
            sys.exit(f"final_eval: {view.name} {key} {got} differs from frozen_views.json {e[key]}: refusing to score")


def code_sha(root: Path = ROOT) -> str | None:
    """The code identity of a run: <root>/CODE_SHA stripped (written by vast sync_up on the box; empty = absent), else
    git rev-parse HEAD, else None."""
    f = Path(root) / "CODE_SHA"
    text = f.read_text().strip() if f.exists() else ""
    return text or git("rev-parse", "HEAD", root=root, required=False)


def hash_view(view: Path, out: Path) -> None:
    """--hash-view: the view's identity record (the pre-rental record of A3 step 1)."""
    sets = json.loads((view / "manifest.json").read_text())["sets"]
    (target_set,) = [s["target_set_sha256"] for k, s in sets.items() if k.startswith("lidc/")]
    commits = render_files()
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

    doc = json.loads((run / "final_eval.json").read_text())
    items = read_jsonl(run / "items.jsonl")
    table = pq.read_table(store / "lidc" / f"{doc['split']}.parquet", columns=["item_id", "group_id"])
    found = {i: g for i, g in zip(*table.to_pydict().values(), strict=True) if g}
    lost = [it["item_id"] for it in items if it["item_id"] not in found]
    if lost:
        print(f"resummarize: {len(lost)} items have no group_id in {store}, first {lost[0]}", file=sys.stderr)
        sys.exit(2)
    s = summarize(items, found)
    res = {k: doc[k] for k in RUN_KEYS} | s | {"metrics_absent": [m for m in METRICS if m not in s],
                                               "resummarized_from": str(run)}
    out.mkdir(parents=True, exist_ok=True)
    (out / "final_eval.json").write_text(json.dumps(res, indent=1))


def load_sampler(a):
    """(VllmSampler of the stage's ckpt, built once per process, serving a.adapter unmerged when given; params_for's
    dict at K = a.k and seed a.seed with the checkpoint's stop ids)."""
    from objdel.train import sample

    s = sample.VllmSampler(Path(a.ckpt), a.model, adapter=a.adapter and Path(a.adapter))
    stop = sorted(sample.stop_ids(Path(a.ckpt), s.processor))
    return s, sample.params_for(a.model, n=a.k, seed=a.seed, stop_token_ids=stop)


def _smi(query: str) -> list[str]:
    """The non-empty CSV rows of one nvidia-smi query."""
    out = subprocess.run(["nvidia-smi", query, "--format=csv,noheader"], capture_output=True, text=True, check=True)
    return [x for x in out.stdout.splitlines() if x.strip()]


def _proc_count() -> int:
    """Processes in /proc that hold /dev/nvidia<N> open: the container's own view (D-G17 (a))."""
    n = 0
    for p in Path("/proc").glob("[0-9]*"):
        with contextlib.suppress(OSError):
            n += any(re.fullmatch(r"/dev/nvidia\d+", os.readlink(f)) for f in (p / "fd").iterdir())
    return n


def _sample(method: str) -> tuple[int, int]:
    """(device processes by the method, memory.used in MiB summed over the GPUs)."""
    used = sum(int(x.split()[0]) for x in _smi("--query-gpu=memory.used"))
    return (len(_smi("--query-compute-apps=pid,used_memory")) if method == "nvidia-smi" else _proc_count()), used


@contextlib.contextmanager
def topology(name: str):
    """D-G17 (a), (c) around one view's scoring, only on a CUDA build with OBJDEL_RUN_DIR set: the device process
    count before scoring (nvidia-smi's compute-apps rows; when there are none, as in many containers, the /proc
    count), a sampler thread every TOPOLOGY_PERIOD_S during it and one sample at its end. Adds the view's method,
    count before, largest count during and peak memory.used (MiB) to $OBJDEL_RUN_DIR/d/topology.json. In this
    process because "before scoring" is a point inside final_eval, after sampling."""
    run_dir = os.environ.get("OBJDEL_RUN_DIR")
    if not run_dir or _gpucheck.build() != "cuda":
        yield
        return
    method = "nvidia-smi" if _smi("--query-compute-apps=pid,used_memory") else "proc"
    seen, stop = [_sample(method)], threading.Event()

    def run():
        while not stop.wait(TOPOLOGY_PERIOD_S):
            seen.append(_sample(method))

    thread = threading.Thread(target=run, daemon=True)
    thread.start()
    try:
        yield
    finally:
        stop.set()
        thread.join()
    seen.append(_sample(method))
    path = Path(run_dir) / "d" / "topology.json"
    doc = json.loads(path.read_text()) if path.exists() else {}
    doc[name] = {"method": method, "before": seen[0][0], "max_during": max(c for c, _ in seen[1:]),
                 "peak_memory_used_mib": max(u for _, u in seen)}
    path.parent.mkdir(parents=True, exist_ok=True)
    _write(path, json.dumps(doc, indent=1))


def engine_record() -> dict:
    """D-G14: the view's eval_sets wall time, the scratch cap the loop passed and the engine's memory counters, read
    after its scoring (the counters were reset just before it)."""
    s, torch = _gpucheck.stats(), sys.modules.get("torch")  # a rescore imports no torch
    return {"eval_sets_s": device.eval_sets_s, "max_scratch_bytes": device.MAX_SCRATCH_BYTES,
            **{k: s[k] for k in ("arena_bytes", "table_bytes", "eval_scratch_peak", "live_device_bytes")},
            "max_memory_reserved": torch.cuda.max_memory_reserved() if torch else 0}


@contextlib.contextmanager
def scoring(view: Path):
    """Around one view's scoring: the engine counters reset first, the topology sampler around it."""
    _gpucheck.reset_stats()
    device.eval_sets_s = 0.0
    with topology(view.name):
        yield


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
            "view_content_sha256": Q.view_content_sha256(view), "code_sha": code_sha(),
            "checker_build": _gpucheck.build()}


def rescore(src: Path, view: Path, split: str, out: Path) -> None:
    """--rescore: score src/answers.jsonl again on the view, with no model; exit 2 when the view's manifest is not
    the one src was sampled on."""
    from objdel.train import dataset, serve

    f = src / "final_eval.json"
    old = json.loads(f.read_text()) if f.exists() else {}
    if "view_manifest_sha256" not in old:
        sys.exit(f"rescore: {f} is missing or has no view_manifest_sha256")
    if old["view_manifest_sha256"] != view_manifest_sha256(view):
        print(f"rescore: {view}/manifest.json is not the manifest {src} was sampled on", file=sys.stderr)
        sys.exit(2)
    t0 = time.time()
    answers = read_jsonl(src / "answers.jsonl")
    by_key = {r["key"]: r for r in dataset.prompts(view, "lidc", split)}
    rows = [by_key[d["item_id"]] for d in answers]
    run = {k: old[k] for k in ("ckpt", "K", "seed", "sampler") if k in old} | _identity(view, split) | {"rescored_from": str(src),
                                                                                 "seconds": {}}
    serve.preload(view, ["lidc"])
    _write(out / "answers.jsonl", _lines(answers))
    with scoring(view):
        results, _ = score_items(rows, {d["item_id"]: d for d in answers}, view.resolve())
    write_scored(view, split, out, results, run, t0)


def no_model(argv: list[str]) -> None:
    """--hash-view VIEW --out FILE and --resummarize DIR --out DIR [--store S]: no model, no EvalStage."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--hash-view", type=Path, metavar="VIEW", help="write VIEW's identity record to --out and exit")
    mode.add_argument("--resummarize", type=Path, metavar="DIR", help="a finished run: summarize its items.jsonl again")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--store", type=Path, default=Path("data/store"), help="the item store (group_id per item)")
    a = ap.parse_args(argv)
    if a.hash_view:
        return hash_view(a.hash_view, a.out)
    if a.out.resolve() == a.resummarize.resolve():
        ap.error("--out must not be the --resummarize folder")
    return resummarize(a.resummarize, a.out, a.store)


def main(argv=None) -> None:
    """The no-model forms, --rescore, or sample and score every view (module docstring)."""
    from objdel.core.spec import EvalStage, cli_stage

    argv = sys.argv[1:] if argv is None else list(argv)
    if {"--hash-view", "--resummarize"} & set(argv):
        return no_model(argv)
    a = cli_stage(EvalStage, argv)

    def error(msg: str):
        print(f"final_eval: error: {msg}", file=sys.stderr)
        sys.exit(2)

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

    from objdel.train import dataset, serve

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
        with contextlib.ExitStack() as stack:

            def keep(answers, run=run, stack=stack, out=out, view=view, groups=groups):
                """answers.jsonl and the partial final_eval.json before scoring (a crash keeps the sampling)."""
                lines = [{"item_id": d["item_id"], "group_id": groups[d["item_id"]], "W": r["W"], "H": r["H"],
                          **{k: d[k] for k in ("samples", "truncated", "tokens")}}
                         for r, d in zip(rows, answers, strict=True)]
                run["seconds"] = {"load_and_sample": time.time() - t0}
                _write(out / "answers.jsonl", _lines(lines))
                _write(out / "final_eval.json", json.dumps(run | {"partial": True}, indent=1))
                stack.enter_context(scoring(view))

            torch.manual_seed(a.seed)
            results = evaluate(sampler, rows, params, view.resolve(), keep)["results"]
        write_scored(view, split, out, results, run, t0)
        t0 = time.time()

if __name__ == "__main__":
    main()
