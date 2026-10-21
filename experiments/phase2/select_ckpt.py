"""D12 checkpoint selection over a train_sft run's WSD branches (plan D12, D13, §3.2 ②, §5.1); never touches CUDA.

    python -I -m experiments.phase2.select_ckpt --run runs/r1 --shifted data/eval_views/lidc-val-shifted \
        --vision '{"prepared": root, "reference": R0 base value}' | --no-vision \
        [--sampler vllm|hf] [--items 256] [--rescaled V] [--rotated V] [--prepared P] [--probes DIR]

Order (D13): refuse unless train_sft has exited (manifest status "done", pid gone on this host); preload the eval
views' caches (objdel.train.serve.preload) and fork the reward pool before anything imports CUDA; per branch, export
the prompts (objdel.train.sample.export_prompts), launch the sampler child (vLLM's `sample_vllm.py` under
`/opt/vllm/bin/python`, or `experiments.phase2.hf_sample`; same JSONL bridge) after asserting
torch.cuda.is_initialized() is False here, with a timeout, at most TRIES attempts, and a completeness check (every
prompt id has its n samples; a partial file is a failed attempt); completions whose `.meta.json` records the same
prompts' SHA-256 are reused on a rerun. Then import the completions (truncation by TRL's rule), score them on the pool
(objdel.eval.metrics with the preloaded CacheRef), build the B2 signal panel (objdel.eval.signal) and G-vision on the
CPU (objdel.eval.vision). Eligible (D12) = gates.eligible: phase2_entry's G-format, G-signal (groups of 8) and
G-vision on the shifted view, one rule with the gate table, and an unmeasured G-vision fails (`--vision` is required; `--no-vision`
makes every branch ineligible). Rule R-B (prereg.select_rb) picks the branch; none eligible -> the deepest is gated
(§6.3). `--items N` is a seeded subset (SUBSET_SEED) for branch selection only (D13 HF path).

The selected (or gated) branch then gets, always on the full views: shifted (re-sampled when --items was used) and
rescaled at K = 8 (gates), shifted at K = 16 (HM-IoU16 / GED16, report-only), G-prior (b) (--prepared), both
probes on every probe record (--probes; probe records are rl: false, so not grpo_rows), and the free attribution
check (--rotated: rotated vs identity shifted view, paired median IoU_tc, report-only; plan §6.1). Writes
`<out>/selected.json`, `<out>/report.json` and `<out>/gates.json` (prereg.entry). Samples stay in
`<out>/<branch>/<view>/` for re-scoring on the Mac.

References: PLAN D12 (R-B), D13 (sampler sequencing), §5.1 gates, §6.1 R1; vLLM offline inference
https://docs.vllm.ai/en/stable/serving/offline_inference.html
"""

import argparse
import hashlib
import json
import os
import socket
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np

from objdel.eval import gates, metrics  # imported before the pool forks: a forked child must not import (macOS objc)
from objdel.train import serve

ROOT = Path(__file__).resolve().parents[2]
VLLM_PYTHON = "/opt/vllm/bin/python"
SUBSET_SEED = 0  # --items: one fixed subset for every branch and rung
TRIES = 2  # sampler child attempts


def exited(run: Path) -> dict:
    """train_sft's manifest, after checking that it finished and its process is gone. A pid recorded on another host
    (manifest "host") cannot be checked here; PermissionError means the pid exists (another user's process)."""
    man = json.loads((run / "manifest.json").read_text())
    alive = False
    if man.get("host", socket.gethostname()) == socket.gethostname() and man["pid"] != os.getpid():
        try:
            os.kill(man["pid"], 0)
            alive = True
        except PermissionError:
            alive = True
        except ProcessLookupError:
            pass
    if man.get("status") != "done" or alive:
        raise RuntimeError(f"{run}: train_sft status {man.get('status')!r}, pid alive {alive}: it must exit first")
    return man


def cuda_untouched() -> None:
    torch = sys.modules.get("torch")
    if torch is not None and torch.cuda.is_initialized():
        raise RuntimeError("CUDA is initialized in select_ckpt: the sampler child would share the GPU (D13)")


def sampler_cmd(kind: str, prompts: Path, model: Path, out: Path) -> list[str]:
    head = [VLLM_PYTHON, "-I", str(ROOT / "experiments/phase2/sample_vllm.py")]
    head = head if kind == "vllm" else [sys.executable, "-I", "-m", "experiments.phase2.hf_sample"]
    return head + ["--prompts", str(prompts), "--model", str(model), "--out", str(out)]


def complete(comps: Path, want: Counter) -> bool:
    """Every prompt id has exactly its n completions (a cut or partial file is not complete)."""
    try:
        return Counter(json.loads(x)["id"] for x in comps.read_text().splitlines()) == want
    except (OSError, ValueError, KeyError):
        return False


def run_sampler(kind: str, prompts: Path, model: Path, comps: Path, timeout: float) -> None:
    """The sampler child, bounded: reuse matching completions, else TRIES attempts of `timeout` s each."""
    want = Counter({p["id"]: p["n"] for p in map(json.loads, prompts.read_text().splitlines())})
    sha, meta = hashlib.sha256(prompts.read_bytes()).hexdigest(), comps.with_suffix(".meta.json")
    if meta.is_file() and json.loads(meta.read_text()).get("prompts_sha256") == sha and complete(comps, want):
        return
    meta.unlink(missing_ok=True)
    for attempt in range(1, TRIES + 1):
        cuda_untouched()
        try:
            subprocess.run(sampler_cmd(kind, prompts, model, comps), check=True, timeout=timeout)
            err = None if complete(comps, want) else "partial output (an id without its n samples)"
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            err = str(e)
        if err is None:
            meta.write_text(json.dumps({"prompts_sha256": sha}) + "\n")
            return
        print(f"select_ckpt: sampler attempt {attempt}/{TRIES} on {model}: {err}", file=sys.stderr)
    raise RuntimeError(f"sampler failed {TRIES} times on {prompts}: {err}")


def _score(job) -> dict:
    """One item's metrics on the preloaded caches (runs in a pool worker)."""
    item_id, texts, truncated, W, H, names, root = job
    return {
        "item_id": item_id,
        "dataset": item_id.split("/")[0],
        **metrics(texts, {}, W, H, names, truncated, caches=serve.CacheRef(item_id, root)),
    }


def view_rows(view: Path, items: int | None = None, rl_only: bool = True) -> list[dict]:
    """Prompt rows of an eval view's lidc/val records (images checked against their records): the scored items (rl,
    with caches) or, with rl_only=False, every record (probe records are rl: false); `items` = a seeded subset."""
    from objdel.train import dataset as D

    recs = [r for r in D.records(view, "lidc", "val") if r["rl"] or not rl_only]
    if items and items < len(recs):
        recs = [recs[i] for i in sorted(np.random.default_rng(SUBSET_SEED).choice(len(recs), items, replace=False))]
    names = list(D.SCHEMAS["lidc"])
    return [
        {"prompt": [D._user("lidc")], "image": D._image(view, r["item_id"]), "key": r["item_id"]}
        | {"W": r["W"], "H": r["H"], "names": names}
        for r in D._verified(view, recs)
    ]


def sample_view(a, branch: Path, rows: list[dict], out: Path, k: int) -> dict:
    """The sampler child on `rows` -> {key: completions row}."""
    from types import SimpleNamespace

    from transformers import AutoProcessor, GenerationConfig

    from objdel.train import rewards, sample

    out.mkdir(parents=True, exist_ok=True)
    proc = AutoProcessor.from_pretrained(branch)  # saved by train_sft with the run's processor kwargs
    stop = rewards.stop_ids(SimpleNamespace(generation_config=GenerationConfig.from_pretrained(branch)), proc)
    prompts, comps = out / "prompts.jsonl", out / "completions.jsonl"
    sample.export_prompts(rows, prompts, a.model, "lidc", k, processor=proc, stop_ids=stop, seed=a.seed)
    run_sampler(a.sampler, prompts, branch, comps, 600 + a.timeout_s * len(rows) * k)
    return {r["item_id"]: r for r in sample.import_completions(comps, stop, n=k, ids=[r["key"] for r in rows])}


def score_view(a, branch: Path, view: Path, out: Path, pool, k: int) -> dict:
    """sample_view on the view (a.items: subset), then metrics on the pool -> summarize (eval_run's report form), the
    signal panel and the per-item K-sample mean IoU_tc per structure (G-prior)."""
    from objdel.eval import signal, summarize

    rows = view_rows(view, a.items)
    done = sample_view(a, branch, rows, out, k)
    root = str(Path(view).resolve())
    jobs = [(m["key"], done[m["key"]]["samples"], done[m["key"]]["truncated"], m["W"], m["H"], m["names"], root)
            for m in rows]  # fmt: skip
    results = list(pool.map(_score, jobs, chunksize=4))
    (out / "items.jsonl").write_text("".join(json.dumps(r) + "\n" for r in results))
    groups = [signal.group(r["item_id"], done[r["item_id"]]["samples"], r) for r in results]
    iou = {r["item_id"]: {n: float(np.mean(o["iou_tc"])) for n, o in r["objects"].items()} for r in results}
    return {"summary": summarize(results), "panel": signal.panel(groups), "iou": iou}


def prior_ci(model_iou: dict, prior_rows: list[dict]) -> dict:
    """G-prior (b): paired per-item (model - image-blind prior) IoU_tc per structure, nodule-cluster bootstrap."""
    from objdel.eval import signal

    prior = {r["item_id"]: r["iou_tc"] for r in prior_rows}
    ids = [i for i in model_iou if i in prior]
    out, D = {}, signal.draws(ids)
    for n in model_iou[ids[0]] if ids else ():
        mean, lo, hi = signal.ci(np.array([model_iou[i][n] - prior[i][n] for i in ids]), D)
        out[n] = {"diff": mean, "ci": [lo, hi], "items": len(ids)}
    return {"lidc": out}


def judge(a, sel: dict, out: Path, pool) -> dict:
    """The §5.1 report of the selected (or gated) branch, on the full views (--items never applies here, D13)."""
    from experiments.phase2 import baselines
    from objdel.eval import probes, signal

    best, whole = a.run / "branches" / sel["name"], argparse.Namespace(**vars(a) | {"items": None})

    def full(view, sub, k=a.k):
        return score_view(whole, best, view, out / sub, pool, k)

    scored = {"shifted": full(a.shifted, "shifted-full") if a.items else sel}
    scored |= {"rescaled": full(a.rescaled, "rescaled")} if a.rescaled else {}
    report = {"views": {v: {"summary": r["summary"]} for v, r in scored.items()}, "vision": sel["vision"]}
    report |= {"signal": scored["shifted"]["panel"], "k16": full(a.shifted, "shifted-k16", 16)["summary"]}
    if a.rotated:  # free attribution check: rotated ~ identity means the dihedral training transferred
        rot = full(a.rotated, "rotated")
        diff = signal.paired(rot["panel"], scored["shifted"]["panel"], "iou_median")
        report["attribution"] = {"summary": rot["summary"], "iou_median_rotated_minus_identity": diff}
    if a.prepared:
        fit = baselines.fit_prior(a.prepared, "lidc", "train")
        views = {"shifted": a.shifted, "rescaled": a.rescaled}
        report["prior"] = {
            v: prior_ci(r["iou"], baselines.score_prior(views[v], "lidc", "val", fit)) for v, r in scored.items()
        }
    if a.probes:
        report["probes"] = {}
        for p, scorer in (("precision", probes.score_precision), ("topology", probes.score_topology)):
            d = a.probes / f"probe-{p}"
            done = sample_view(whole, best, view_rows(d, rl_only=False), out / f"probe-{p}", a.k)
            report["probes"][p] = scorer(list(done.values()), d)
    return report


def vision_config(a) -> dict | None:
    if a.no_vision:
        return None
    if a.vision is None:
        raise SystemExit("select_ckpt: --vision is required (G-vision decides D12 eligibility); --no-vision opts out")
    cfg = json.loads(a.vision)
    ref = cfg.get("reference")
    ok = Path(ref).is_file() if isinstance(ref, str) else isinstance(ref, int | float) and ref > 0
    if not isinstance(cfg.get("prepared"), str) or not ok:
        raise SystemExit(f"select_ckpt: --vision needs prepared (a path) and reference (> 0, or R0's JSON): {cfg}")
    return cfg


def vision_check(branch: Path, cfg: dict | None) -> dict:
    """G-vision on the CPU in float32 (B1 protocol); not measured (pass None, so it fails) under --no-vision."""
    if cfg is None:
        return {"pass": None, "measured": False}
    from transformers import AutoModelForImageTextToText, AutoProcessor

    from experiments.phase2.train_sft import vision_check as check

    model = AutoModelForImageTextToText.from_pretrained(branch, dtype="float32", device_map="cpu")
    return check(model, AutoProcessor.from_pretrained(branch), cfg)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", type=Path, required=True, help="train_sft output_dir")
    ap.add_argument("--shifted", type=Path, required=True, help="shifted-val eval view (selection)")
    ap.add_argument("--rescaled", type=Path, help="rescaled-val eval view (selected branch only)")
    ap.add_argument("--rotated", type=Path, help="rotated shifted-val eval view (attribution check, report-only)")
    ap.add_argument("--prepared", type=Path, help="prepared root: the image-blind prior is fitted on lidc/train")
    ap.add_argument("--probes", type=Path, help="dir with probe-precision/ and probe-topology/ (selected branch)")
    ap.add_argument("--out", type=Path, help="default <run>/select")
    ap.add_argument("--sampler", choices=("vllm", "hf"), default="vllm")
    ap.add_argument("--model", help="configs.MODELS key for export_prompts (default: the run's)")
    ap.add_argument("--items", type=int, help="seeded subset size for branch selection only (HF path: 256)")
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0, help="sampler seed")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--timeout_s", type=float, default=2.0, help="sampler child timeout per sample, plus 600 s")
    vis = ap.add_mutually_exclusive_group()
    vis.add_argument("--vision", help='G-vision JSON {"prepared": root, "reference": R0 value or JSON} (required)')
    vis.add_argument("--no-vision", action="store_true", help="G-vision not measured: every branch ineligible")
    a = ap.parse_args(argv)
    vision_cfg = vision_config(a)
    man = exited(a.run)
    a.model = a.model or man["config"]["model"]
    vision_cfg = vision_cfg and {"model": a.model} | vision_cfg  # train_sft.reference checks R0's JSON model
    out = a.out or a.run / "select"
    for v in filter(None, (a.shifted, a.rescaled, a.rotated)):
        serve.preload(v, ["lidc"])  # every scored root before the fork: workers see STORE as of fork time
    pool = serve.fork_pool(a.procs)  # before torch / CUDA (D14)
    branches = []
    for b in sorted((a.run / "branches").iterdir(), key=lambda p: float(p.name[1:])):
        res = score_view(a, b, a.shifted, out / b.name / "shifted", pool, a.k)
        res["vision"] = vision_check(b, vision_cfg)
        res["eligible"] = gates.eligible(res["summary"], res["vision"])  # D12: G-format, G-signal, G-vision
        branches.append({"name": b.name, "depth": float(b.name[1:]), **res})
    from experiments.phase2 import prereg

    name, why = prereg.select_rb(branches)
    gate, judged = None, name or (branches[-1]["name"] if branches else None)  # none eligible: the deepest (§6.3)
    if judged:
        report = judge(a, next(x for x in branches if x["name"] == judged), out / judged, pool)
        gate = prereg.entry(report, ["lidc"])  # phase2_entry with the pre-registered G-prior (a) slot
        (out / "report.json").write_text(json.dumps(report, indent=1, default=str) + "\n")
        (out / "gates.json").write_text(json.dumps(gate, indent=1, default=str) + "\n")
    pool.shutdown()
    torch = sys.modules.get("torch")
    doc = {
        "cuda_initialized": bool(torch and torch.cuda.is_initialized()),
        "run": str(a.run),  # g_probe --select reads the branches from here
        "selected": name,
        "reasons": why,
        "gates_branch": judged,
        "gates_pass": gate and gate["pass"],
        "sampler": a.sampler,
        "k": a.k,
        "items": a.items,
        "branches": {
            x["name"]: {k: x[k] for k in ("depth", "eligible", "vision")}
            | {"summary": x["summary"]["datasets"]["lidc"], "panel": x["panel"]["summary"]}
            for x in branches
        },
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "selected.json").write_text(json.dumps(doc, indent=1, default=str) + "\n")
    print(json.dumps({"selected": name, "eligible": why["eligible"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
