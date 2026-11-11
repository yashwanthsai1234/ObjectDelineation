"""Checkpoint selection over a train_sft run's WSD branches; never touches CUDA.

`run(stage, prereg, baselines)` runs one objdel.core.spec.SelectStage: the sampler children, scoring per view and
probe, G-vision, the HF check, the judge and the gates. experiments/phase2/select_ckpt.py is its CLI and passes the
pre-registered record (prereg: rule R-B, the entry gates) and the experiment baselines, which src never imports.

Order: refuse until train_sft has exited (manifest "done", pid gone on this host); preload the views' caches; per
branch, export the prompts and run the sampler child (`python -I -m objdel.policy.sampler --kind vllm|hf`) with CUDA
asserted untouched here, a timeout of 600 s + --timeout-s per item-sample and a completeness check (every prompt id
has its n samples). The completions' `.meta.json` stamp (prompts SHA-256, child command) is written before the child
starts: a matching stamp reuses the file when complete and resumes it otherwise, another stamp empties it; an attempt
that finishes no new item counts against TRIES. SIGTERM / SIGINT kill the child's process group before this exits.
Scoring runs in this process; G-vision runs on the CPU. A branch is eligible when it passes gates.eligible (G-format,
G-signal, G-vision; an unmeasured G-vision fails, so --no-vision makes every branch ineligible); rule R-B
(prereg.select_rb) picks, and with none eligible the deepest is gated. --items N: a seeded subset for selection only.

The selected (or gated) branch is judged on the full views: shifted and rescaled at K = 8 (gates), shifted at K = 16
(report-only), the prior (--prepared), both probes on every probe record (--probes) and the rotated-vs-identity check
(--rotated, report-only). With --sampler vllm, `hf_check` re-samples that branch with HF on HF_CHECK_ITEMS shifted-val
items; a vLLM - HF metric difference outside its 95 % CI fails the gates. Writes <out>/selected.json ("gates_pass" and
"gates_failed" are what g_probe --select reads), report.json and gates.json; samples stay in <out>/<branch>/<view>/.

--only topology: the topology probe alone on the run's single branch at K = --k into <out>/report.json; no selection,
other views, G-vision, prior or gates, and --vision / --no-vision are refused.
"""

import argparse
import contextlib
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

from objdel.core.io import read_json, read_jsonl
from objdel.policy import evaluate
from objdel.policy.errors import Refusal

TRIES = 2  # sampler child attempts that finish no new item


# sampler child timeout per item-sample, plus 600 s: 3x the measured time (EV-5), so a hung engine is killed in
# minutes, not hours. HF (~1 s per answer) resumes from its whole batches. vLLM writes once at the end, so a retry
# starts over: its time is scaled to answers at the 1,536-token cap (0.09 s at ~580 tokens in R0 v3 -> 0.25 s)
TIMEOUT_S = {"vllm": 0.75, "hf": 3.0}


HF_CHECK_ITEMS = 64  # HF re-samples the vLLM-judged branch on this seeded shifted-val subset


HF_CHECK_KEYS = ("format_ok", "legal", "iou_tc", "f5", "g_area")


def exited(run: Path) -> dict:
    """train_sft's manifest, after checking that it finished and its process is gone. A pid recorded on another host
    (manifest "host") cannot be checked here; PermissionError means the pid exists (another user's process)."""
    man, host = read_json(run / "manifest.json"), socket.gethostname()
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
    """Refuse when this process initialised CUDA: the sampler child needs the whole GPU."""
    torch = sys.modules.get("torch")
    if torch is not None and torch.cuda.is_initialized():
        raise RuntimeError("CUDA is initialized in select_ckpt: the sampler child would share the GPU (D13)")


def sampler_cmd(kind: str, prompts: Path, model: Path, out: Path, extra=()) -> list[str]:
    """The sampler child's argv: objdel.policy.sampler in this venv (vllm_native: one venv for vLLM and objdel)."""
    return [sys.executable, "-I", "-m", "objdel.policy.sampler", "--kind", kind, "--prompts", str(prompts), "--model",
            str(model), "--out", str(out), *extra]


def complete(comps: Path, want: Counter) -> bool:
    """Every prompt id has exactly its n completions (a cut or partial file is not complete)."""
    try:
        return Counter(r["id"] for r in read_jsonl(comps)) == want
    except (OSError, ValueError, KeyError):
        return False


def finished(comps: Path, want: Counter) -> int:
    """Prompt ids with all their n completions, cut lines skipped (progress of a resumable child)."""
    got = Counter()
    for x in comps.read_text().splitlines() if comps.exists() else []:
        try:
            got[json.loads(x)["id"]] += 1
        except (ValueError, KeyError, TypeError):
            continue
    return sum(got[i] == n for i, n in want.items())


def run_group(cmd: list[str], timeout: float) -> None:
    """subprocess.run(check=True, timeout=...) in a new session, killing the whole group on timeout and on SIGTERM /
    SIGINT to this process (the run cap's `timeout`): vLLM's engine is a child process that holds the GPU, and the
    next job in the chain needs it."""
    with subprocess.Popen(cmd, start_new_session=True) as p:

        def kill(*_):
            with contextlib.suppress(ProcessLookupError):
                os.killpg(p.pid, signal.SIGKILL)

        def stop(sig, frame):
            kill()
            raise SystemExit(128 + sig)

        old = {s: signal.signal(s, stop) for s in (signal.SIGTERM, signal.SIGINT)}
        try:
            if p.wait(timeout):
                raise subprocess.CalledProcessError(p.returncode, cmd)
        except subprocess.TimeoutExpired:
            kill()
            raise
        finally:
            for s, h in old.items():
                signal.signal(s, h)


def run_sampler(kind: str, prompts: Path, model: Path, comps: Path, timeout: float, extra=()) -> None:
    """The sampler child (extra: more flags, e.g. --logprobs), bounded: completions of the same prompts and child
    command (sampler, checkpoint, flags) are reused when complete and resumed otherwise; attempts of `timeout` s
    until TRIES of them finish no new item."""
    want = Counter({p["id"]: p["n"] for p in read_jsonl(prompts)})
    cmd, meta = sampler_cmd(kind, prompts, model, comps, extra), comps.with_suffix(".meta.json")
    stamp = {"prompts_sha256": hashlib.sha256(prompts.read_bytes()).hexdigest(), "cmd": cmd}
    if not (meta.is_file() and read_json(meta) == stamp):
        comps.unlink(missing_ok=True)  # another sampling's completions: never resumed
        meta.write_text(json.dumps(stamp) + "\n")
    stalled = 0
    while not complete(comps, want):
        before = finished(comps, want)
        cuda_untouched()
        try:
            run_group(cmd, timeout)
            err = "partial output (an id without its n samples)"
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            err = str(e)
        if complete(comps, want):
            return
        stalled += finished(comps, want) <= before
        print(f"select_ckpt: sampler on {model} ({stalled}/{TRIES} attempts without progress): {err}",
              file=sys.stderr)
        if stalled >= TRIES:
            raise RuntimeError(f"sampler failed {TRIES} times on {prompts}: {err}")


def sampler_timeout(a, samples: int) -> float:
    """600 s (engine start) + the per-sample time (--timeout-s, else TIMEOUT_S of the sampler) for `samples`."""
    return 600 + (a.timeout_s or TIMEOUT_S[a.sampler]) * samples


def sample_view(a, branch: Path, rows: list[dict], out: Path, k: int, extra=()) -> dict:
    """The sampler child on `rows` (`out`/prompts.jsonl -> `out`/completions.jsonl) -> {key: completions row}."""
    from transformers import AutoProcessor

    from objdel.policy import sampler as sample

    out.mkdir(parents=True, exist_ok=True)
    proc = AutoProcessor.from_pretrained(branch)  # saved by train_sft with the run's processor kwargs
    stop = sample.stop_ids(branch, proc)
    prompts, comps = out / "prompts.jsonl", out / "completions.jsonl"
    sample.export_prompts(rows, prompts, a.model, "lidc", k, processor=proc, stop_ids=stop, seed=a.seed,
                          max_tokens=a.max_tokens)
    run_sampler(a.sampler, prompts, branch, comps, sampler_timeout(a, len(rows) * k), extra)
    return {r["item_id"]: r for r in sample.import_completions(comps, stop, n=k, ids=[r["key"] for r in rows])}


def score_view(a, branch: Path, view: Path, out: Path, k: int) -> dict:
    """sample_view on the view (a.items: subset), then metrics -> summarize (eval_run's report form), the
    signal panel and the per-item K-sample mean IoU_tc per structure (G-prior)."""
    from objdel.verify import stats as signal
    from objdel.verify.metrics import summarize

    rows = evaluate.eval_rows(view, a.items)
    results, groups = evaluate.score_items(rows, sample_view(a, branch, rows, out, k), view.resolve())
    (out / "items.jsonl").write_text("".join(json.dumps(r) + "\n" for r in results))
    iou = {r["item_id"]: {n: float(np.mean(o["iou_tc"])) for n, o in r["objects"].items()} for r in results}
    return {"summary": summarize(results), "panel": signal.panel(groups), "iou": iou}


def prior_ci(model_iou: dict, prior_rows: list[dict]) -> dict:
    """G-prior (b): paired per-item (model - image-blind prior) IoU_tc per structure, nodule-cluster bootstrap."""
    from objdel.verify import stats as signal

    prior = {r["item_id"]: r["iou_tc"] for r in prior_rows}
    ids = [i for i in model_iou if i in prior]
    out, D = {}, signal.point_or_draws(ids)
    for n in model_iou[ids[0]] if ids else ():
        mean, lo, hi = signal.ci(np.array([model_iou[i][n] - prior[i][n] for i in ids]), D)
        out[n] = {"diff": mean, "ci": [lo, hi], "items": len(ids)}
    return {"lidc": out}


def score_probe(a, branch: Path, p: str, out: Path) -> dict:
    """Probe `p` ("precision" / "topology", every record of a.probes/probe-<p>) at K = a.k -> its scorer's result
    and the sampled item ids ("ids")."""
    from objdel.data import probes

    d = a.probes / f"probe-{p}"
    done = sample_view(a, branch, evaluate.eval_rows(d, rl_only=False), out / f"probe-{p}", a.k)
    return getattr(probes, f"score_{p}")(list(done.values()), d) | {"ids": list(done)}  # ids pair the seeds


def hf_check(a, name: str, out: Path) -> dict:
    """With --sampler vllm, do vLLM's answers score like HF's on the trained model? HF samples
    branch `name` on HF_CHECK_ITEMS seeded shifted-val items at K; vLLM's answers are that branch's shifted samples
    (<out>/<name>/shifted/items.jsonl). Per metric, the per-item mean difference vLLM - HF with the nodule-cluster
    bootstrap 95 % CI (every difference inside its CI); pass = 0 inside every CI."""
    from objdel.verify import stats as signal

    rows = evaluate.eval_rows(a.shifted, HF_CHECK_ITEMS)
    hf = argparse.Namespace(**vars(a) | {"sampler": "hf"})
    done = sample_view(hf, a.run / "branches" / name, rows, out / name / "hf-check", a.k)
    got = {"hf": {r["item_id"]: r["samples"] for r in evaluate.score_items(rows, done, a.shifted.resolve())[0]}}
    got["vllm"] = {r["item_id"]: r["samples"] for r in read_jsonl(out / name / "shifted" / "items.jsonl")}
    ids = [i for i in got["hf"] if i in got["vllm"]]
    D = signal.point_or_draws(ids)
    diff = {k: signal.ci(np.array([np.mean(got["vllm"][i][k]) - np.mean(got["hf"][i][k]) for i in ids]), D)
            for k in HF_CHECK_KEYS}  # fmt: skip
    return {"branch": name, "items": len(ids), "k": a.k, "vllm_minus_hf": diff,
            "pass": len(ids) == len(rows) and all(lo <= 0 <= hi for _, lo, hi in diff.values())}


def vision_check(branch: Path, cfg: dict | None) -> dict:
    """G-vision on the CPU in float32 ; not measured (pass None, so it fails) under --no-vision."""
    if cfg is None:
        return {"pass": None, "measured": False}
    from transformers import AutoModelForImageTextToText, AutoProcessor

    from objdel.sft.train import vision_check as check

    model = AutoModelForImageTextToText.from_pretrained(branch, dtype="float32", device_map="cpu")
    return check(model, AutoProcessor.from_pretrained(branch), cfg)


def judge(a, sel: dict, out: Path, baselines) -> dict:
    """The report of the selected (or gated) branch, on the full views (--items never applies here).
    baselines (the experiment module, passed in by the CLI) fits and scores G-prior (b) under --prepared."""
    from objdel.verify import stats

    best, whole = a.run / "branches" / sel["name"], argparse.Namespace(**vars(a) | {"items": None})

    def full(view, sub, k=a.k):
        return score_view(whole, best, view, out / sub, k)

    scored = {"shifted": full(a.shifted, "shifted-full") if a.items else sel}
    scored |= {"rescaled": full(a.rescaled, "rescaled")} if a.rescaled else {}
    report = {"views": {v: {"summary": r["summary"]} for v, r in scored.items()}, "vision": sel["vision"]}
    report |= {"signal": scored["shifted"]["panel"], "k16": full(a.shifted, "shifted-k16", 16)["summary"]}
    if a.rotated:  # free attribution check: rotated ~ identity means the dihedral training transferred
        rot = full(a.rotated, "rotated")
        diff = stats.paired(rot["panel"], scored["shifted"]["panel"], "iou_median")
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
    """Indented JSON to path (parents created)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=1, default=str) + "\n")


def topology_only(a, out: Path) -> int:
    """--only topology: the topology probe on the run's single branch -> <out>/report.json."""
    names = sorted(p.name for p in (a.run / "branches").iterdir())
    if a.probes is None or len(names) != 1:
        raise Refusal(f"select_ckpt --only topology: needs --probes and one branch, found {names}")
    from objdel.geom import engine as _gpucheck

    topo = score_probe(a, a.run / "branches" / names[0], "topology", out / names[0])
    dump(out / "report.json", {"only": "topology", "branch": names[0], "sampler": a.sampler,
        "checker": _gpucheck.build(), "k": a.k,
                               "probes": {"topology": topo}})
    print(json.dumps({"branch": names[0], "pass": topo.get("pass")}))
    return 0


def vision_config(a) -> dict | None:
    """The G-vision config of the stage (None with --no-vision); Refusal when prepared or reference is bad."""
    if a.no_vision:
        return None
    cfg = a.vision
    ref = cfg.get("reference")
    ok = Path(ref).is_file() if isinstance(ref, str) else isinstance(ref, int | float) and ref > 0
    if not isinstance(cfg.get("prepared"), str) or not ok:
        raise Refusal(f"select_ckpt: --vision needs prepared (a path) and reference (> 0, or R0's JSON): {cfg}")
    return cfg


PATHS = ("run", "shifted", "rescaled", "rotated", "prepared", "probes", "out")


def run(st, prereg, baselines) -> int:
    """The SelectStage `st` -> score every branch, judge, write <out>/selected.json. prereg (rule R-B and the
    phase-2 entry gates) and baselines are the experiment modules, passed in by the CLI (src never imports them)."""
    from objdel.geom import engine as _gpucheck
    from objdel.verify import gates
    from objdel.verify import store as serve

    a = argparse.Namespace(**{k: Path(v) if k in PATHS and v is not None else v
                              for k, v in st.model_dump(exclude={"kind", "name", "soft"}).items()})
    a.sampler = st.sampler or "vllm"
    vision_cfg = None if a.only else vision_config(a)
    a.model, out = a.model or exited(a.run)["config"]["model"], a.out or a.run / "select"
    if a.only:
        return topology_only(a, out)
    vision_cfg = vision_cfg and {"model": a.model} | vision_cfg  # sft.train.reference checks R0's JSON model
    for v in filter(None, (a.shifted, a.rescaled, a.rotated)):
        serve.preload(v, ["lidc"])
    branches = []
    for b in sorted((a.run / "branches").iterdir(), key=lambda p: float(p.name[1:])):
        res = score_view(a, b, a.shifted, out / b.name / "shifted", a.k)
        res["vision"] = vision_check(b, vision_cfg)
        res["eligible"] = gates.eligible(res["summary"], res["vision"])  # G-format, G-signal, G-vision
        branches.append({"name": b.name, "depth": float(b.name[1:]), **res})
    name, why = prereg.select_rb(branches)
    gate, judged = None, name or (branches[-1]["name"] if branches else None)  # none eligible: the deepest (§6.3)
    if judged:
        report = judge(a, next(x for x in branches if x["name"] == judged), out / judged, baselines)
        gate = prereg.entry(report, ["lidc"])  # phase2_entry with the pre-registered G-prior (a) slot
        if a.sampler == "vllm":  # D13: the gates stand only if vLLM's answers score like HF's
            report["hf_check"] = check = hf_check(a, judged, out)
            gate["gates"]["G-hf-check"] = {"requirement": "every vLLM - HF metric difference: 0 inside its 95 % CI",
                                           "values": check["vllm_minus_hf"], "pass": check["pass"]}
            gate["passed"]["G-hf-check"], gate["pass"] = check["pass"], gate["pass"] and check["pass"]
        dump(out / "report.json", report)
        dump(out / "gates.json", gate)
    torch = sys.modules.get("torch")
    dump(out / "selected.json", {
        "cuda_initialized": bool(torch and torch.cuda.is_initialized()),
        "run": str(a.run),  # g_probe --select reads the branches from here
        "selected": name, "reasons": why, "gates_branch": judged, "gates_pass": gate and gate["pass"],
        "gates_failed": [n for n, g in gate["gates"].items() if not g["pass"]] if gate else [],
        "sampler": a.sampler, "checker": _gpucheck.build(), "k": a.k, "items": a.items,
        "branches": {x["name"]: {k: x[k] for k in ("depth", "eligible", "vision")}
                     | {"summary": x["summary"]["datasets"]["lidc"], "panel": x["panel"]["summary"]} for x in branches},
    })  # fmt: skip
    print(json.dumps({"selected": name, "eligible": why["eligible"]}))
    return 0
