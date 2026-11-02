"""D12 checkpoint selection over a train_sft run's WSD branches (PLAN D12, D13, §3.2 ②, §5.1); never touches CUDA.

    python -I -m experiments.phase2.select_ckpt --run runs/r1 --shifted data/eval_views/lidc-val-shifted \
        --vision '{"prepared": root, "reference": R0 base value}' | --no-vision \
        [--sampler vllm|hf] [--items 256] [--rescaled V] [--rotated V] [--prepared P] [--probes DIR]
    python -I -m experiments.phase2.select_ckpt --run runs/c2 --probes DIR --only topology [--sampler hf]

Order (D13): refuse until train_sft has exited (manifest "done", pid gone on this host); preload the views' caches
and fork the reward pool before anything imports CUDA; per branch, export the prompts and run the sampler child
(`sample_vllm.py` under /opt/vllm, or `hf_sample`; one JSONL bridge) with CUDA asserted untouched here, a timeout,
TRIES attempts and a completeness check (every prompt id has its n samples). Completions whose `.meta.json` stamp
(prompts SHA-256, child command) matches are reused on a rerun. Score them on the pool (objdel.eval.metrics), build
the B2 signal panel and run G-vision on the CPU. Eligible (D12) = gates.eligible: G-format, G-signal (groups of 8),
G-vision; an unmeasured G-vision fails (--vision is required, --no-vision makes every branch ineligible). Rule R-B
(prereg.select_rb) picks; none eligible -> the deepest is gated (§6.3). --items N: a seeded subset (SUBSET_SEED)
for selection only (D13 HF path).

The selected (or gated) branch is then judged on the full views: shifted and rescaled at K = 8 (gates), shifted at
K = 16 (report-only), G-prior (b) (--prepared), both probes on every probe record (--probes) and the rotated-vs-
identity attribution check (--rotated, report-only, §6.1). Writes <out>/selected.json ("gates_pass" /
"gates_failed" are what g_probe --select reads), report.json and gates.json (prereg.entry); samples stay in
<out>/<branch>/<view>/ for re-scoring on the Mac.

--only topology (C2, PLAN §6.1 R3-R5): the topology probe alone on the run's single branch at K = --k into
<out>/report.json {"only", "branch", "sampler", "k", "probes": {"topology"}} (c2_analyze reads probes.topology).
No selection, other views, G-vision, prior or gates; their flags are ignored.

References: PLAN D12 (R-B), D13, §5.1, §6.1 R1; https://docs.vllm.ai/en/stable/serving/offline_inference.html
"""

import argparse
import hashlib
import json
import os
import signal
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
    man, host = json.loads((run / "manifest.json").read_text()), socket.gethostname()
    alive = man.get("host", host) == host and man["pid"] != os.getpid()
    try:
        alive and os.kill(man["pid"], 0)
    except ProcessLookupError:
        alive = False
    except PermissionError:
        pass
    if man.get("status") != "done" or alive:
        raise RuntimeError(f"{run}: train_sft status {man.get('status')!r}, pid alive {alive}: it must exit first")
    return man


def cuda_untouched() -> None:
    torch = sys.modules.get("torch")
    if torch is not None and torch.cuda.is_initialized():
        raise RuntimeError("CUDA is initialized in select_ckpt: the sampler child would share the GPU (D13)")


def sampler_cmd(kind: str, prompts: Path, model: Path, out: Path, extra=()) -> list[str]:
    head = [VLLM_PYTHON, "-I", str(ROOT / "experiments/phase2/sample_vllm.py")]
    head = head if kind == "vllm" else [sys.executable, "-I", "-m", "experiments.phase2.hf_sample"]
    return head + ["--prompts", str(prompts), "--model", str(model), "--out", str(out), *extra]


def complete(comps: Path, want: Counter) -> bool:
    """Every prompt id has exactly its n completions (a cut or partial file is not complete)."""
    try:
        return Counter(json.loads(x)["id"] for x in comps.read_text().splitlines()) == want
    except (OSError, ValueError, KeyError):
        return False


def run_group(cmd: list[str], timeout: float) -> None:
    """subprocess.run(check=True, timeout=...) in a new session, killing the whole group on timeout: vLLM's engine is
    a child process that holds the GPU, and the next job in the chain needs it."""
    with subprocess.Popen(cmd, start_new_session=True) as p:
        try:
            if p.wait(timeout):
                raise subprocess.CalledProcessError(p.returncode, cmd)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGKILL)
            raise


def run_sampler(kind: str, prompts: Path, model: Path, comps: Path, timeout: float, extra=()) -> None:
    """The sampler child (extra: more flags, e.g. --logprobs), bounded: reuse completions of the same prompts and
    child command (sampler, checkpoint, flags), else TRIES attempts of `timeout` s each."""
    want = Counter({p["id"]: p["n"] for p in map(json.loads, prompts.read_text().splitlines())})
    cmd, meta = sampler_cmd(kind, prompts, model, comps, extra), comps.with_suffix(".meta.json")
    stamp = {"prompts_sha256": hashlib.sha256(prompts.read_bytes()).hexdigest(), "cmd": cmd}
    if meta.is_file() and json.loads(meta.read_text()) == stamp and complete(comps, want):
        return
    meta.unlink(missing_ok=True)
    for attempt in range(1, TRIES + 1):
        cuda_untouched()
        try:
            run_group(cmd, timeout)
            err = None if complete(comps, want) else "partial output (an id without its n samples)"
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            err = str(e)
        if err is None:
            meta.write_text(json.dumps(stamp) + "\n")
            return
        print(f"select_ckpt: sampler attempt {attempt}/{TRIES} on {model}: {err}", file=sys.stderr)
    raise RuntimeError(f"sampler failed {TRIES} times on {prompts}: {err}")


def _score(job) -> dict:
    """One item's metrics on the preloaded caches (runs in a pool worker)."""
    item_id, texts, truncated, W, H, names, root = job
    m = metrics(texts, {}, W, H, names, truncated, caches=serve.CacheRef(item_id, root))
    return {"item_id": item_id, "dataset": item_id.split("/")[0], **m}


def view_rows(view: Path, items: int | None = None, rl_only: bool = True) -> list[dict]:
    """Prompt rows of an eval view's lidc/val records (images checked against their records): the scored items (rl,
    with caches) or, with rl_only=False, every record (probe records are rl: false); `items` = a seeded subset."""
    from objdel.train import dataset as D

    recs = [r for r in D.records(view, "lidc", "val") if r["rl"] or not rl_only]
    if items and items < len(recs):
        recs = [recs[i] for i in sorted(np.random.default_rng(SUBSET_SEED).choice(len(recs), items, replace=False))]
    return [{"prompt": [D._user("lidc")], "image": D._image(view, r["item_id"]), "key": r["item_id"], "W": r["W"],
             "H": r["H"], "names": list(D.SCHEMAS["lidc"])} for r in D._verified(view, recs)]


def stop_ids(path: Path, proc) -> set[int]:
    """rewards.stop_ids of a checkpoint directory. Base weights may have no generation_config.json (Qwen3.5-0.8B,
    whose eos is only text_config's): then the model config's, as from_pretrained falls back for a loaded model."""
    from types import SimpleNamespace

    from transformers import AutoConfig, GenerationConfig

    from objdel.train import rewards

    try:
        gen = GenerationConfig.from_pretrained(path)
    except OSError:
        gen = GenerationConfig.from_model_config(AutoConfig.from_pretrained(path))  # reads text_config's eos too
    return rewards.stop_ids(SimpleNamespace(generation_config=gen), proc)


def sample_view(a, branch: Path, rows: list[dict], out: Path, k: int, extra=()) -> dict:
    """The sampler child on `rows` (`out`/prompts.jsonl -> `out`/completions.jsonl) -> {key: completions row}."""
    from transformers import AutoProcessor

    from objdel.train import sample

    out.mkdir(parents=True, exist_ok=True)
    proc = AutoProcessor.from_pretrained(branch)  # saved by train_sft with the run's processor kwargs
    stop = stop_ids(branch, proc)
    prompts, comps = out / "prompts.jsonl", out / "completions.jsonl"
    sample.export_prompts(rows, prompts, a.model, "lidc", k, processor=proc, stop_ids=stop, seed=a.seed)
    run_sampler(a.sampler, prompts, branch, comps, 600 + a.timeout_s * len(rows) * k, extra)
    return {r["item_id"]: r for r in sample.import_completions(comps, stop, n=k, ids=[r["key"] for r in rows])}


def score_items(pool, rows: list[dict], done: dict, root: Path) -> tuple[list[dict], list[dict]]:
    """_score on the pool for every row (done: {key: completions row}) -> (per-item metrics, signal groups)."""
    from objdel.eval import signal

    jobs = [(r["key"], done[r["key"]]["samples"], done[r["key"]]["truncated"], r["W"], r["H"], r["names"], str(root))
            for r in rows]  # fmt: skip
    results = list(pool.map(_score, jobs, chunksize=4))
    return results, [signal.group(r["item_id"], done[r["item_id"]]["samples"], r) for r in results]


def score_view(a, branch: Path, view: Path, out: Path, pool, k: int) -> dict:
    """sample_view on the view (a.items: subset), then metrics on the pool -> summarize (eval_run's report form), the
    signal panel and the per-item K-sample mean IoU_tc per structure (G-prior)."""
    from objdel.eval import signal, summarize

    rows = view_rows(view, a.items)
    results, groups = score_items(pool, rows, sample_view(a, branch, rows, out, k), view.resolve())
    (out / "items.jsonl").write_text("".join(json.dumps(r) + "\n" for r in results))
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


def score_probe(a, branch: Path, p: str, out: Path) -> dict:
    """Probe `p` ("precision" / "topology", every record of a.probes/probe-<p>) at K = a.k -> its scorer's result."""
    from objdel.eval import probes

    d = a.probes / f"probe-{p}"
    done = sample_view(a, branch, view_rows(d, rl_only=False), out / f"probe-{p}", a.k)
    return getattr(probes, f"score_{p}")(list(done.values()), d)


def judge(a, sel: dict, out: Path, pool) -> dict:
    """The §5.1 report of the selected (or gated) branch, on the full views (--items never applies here, D13)."""
    from experiments.phase2 import baselines
    from objdel.eval import signal

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
        report["prior"] = {v: prior_ci(r["iou"], baselines.score_prior(views[v], "lidc", "val", fit))
                           for v, r in scored.items()}  # fmt: skip
    if a.probes:
        report["probes"] = {p: score_probe(whole, best, p, out) for p in ("precision", "topology")}
    return report


def dump(path: Path, doc: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=1, default=str) + "\n")


def topology_only(a, out: Path) -> int:
    """--only topology (C2): the topology probe on the run's single branch -> <out>/report.json."""
    names = sorted(p.name for p in (a.run / "branches").iterdir())
    if a.probes is None or len(names) != 1:
        raise SystemExit(f"select_ckpt --only topology: needs --probes and one branch, found {names}")
    topo = score_probe(a, a.run / "branches" / names[0], "topology", out / names[0])
    dump(out / "report.json", {"only": "topology", "branch": names[0], "sampler": a.sampler, "k": a.k,
                               "probes": {"topology": topo}})
    print(json.dumps({"branch": names[0], "pass": topo.get("pass")}))
    return 0


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
    ap.add_argument("--shifted", type=Path, help="shifted-val eval view (selection; required unless --only)")
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
    ap.add_argument("--only", choices=("topology",), help="C2: only the topology probe of a single-branch run")
    a = ap.parse_args(argv)
    if a.shifted is None and not a.only:
        ap.error("--shifted is required (unless --only topology)")
    vision_cfg = None if a.only else vision_config(a)
    a.model, out = a.model or exited(a.run)["config"]["model"], a.out or a.run / "select"
    if a.only:
        return topology_only(a, out)
    vision_cfg = vision_cfg and {"model": a.model} | vision_cfg  # train_sft.reference checks R0's JSON model
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
        dump(out / "report.json", report)
        dump(out / "gates.json", gate)
    pool.shutdown()
    torch = sys.modules.get("torch")
    dump(out / "selected.json", {
        "cuda_initialized": bool(torch and torch.cuda.is_initialized()),
        "run": str(a.run),  # g_probe --select reads the branches from here
        "selected": name, "reasons": why, "gates_branch": judged, "gates_pass": gate and gate["pass"],
        "gates_failed": [n for n, g in gate["gates"].items() if not g["pass"]] if gate else [],
        "sampler": a.sampler, "k": a.k, "items": a.items,
        "branches": {x["name"]: {k: x[k] for k in ("depth", "eligible", "vision")}
                     | {"summary": x["summary"]["datasets"]["lidc"], "panel": x["panel"]["summary"]} for x in branches},
    })  # fmt: skip
    print(json.dumps({"selected": name, "eligible": why["eligible"]}))
    return 0

if __name__ == "__main__":
    sys.exit(main())
