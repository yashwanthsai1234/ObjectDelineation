"""D12 checkpoint selection over a train_sft run's WSD branches (PLAN D12, D13, §3.2 ②, §5.1); never touches CUDA.

The sampler children, scoring per view and probe, G-vision and the HF check; experiments/phase2/select_ckpt.py is
the CLI (it keeps judge, which reads the experiment baselines, and the pre-registered rule R-B).


Order (D13): refuse until train_sft has exited (manifest "done", pid gone on this host); preload the views' caches;
per branch, export the prompts and run the sampler child
(`python -I -m objdel.train.sample --kind vllm|hf`; one JSONL bridge) with CUDA asserted untouched here, a timeout
of 600 s + --timeout-s per item-sample, and a completeness check (every prompt id has its n samples). The
completions' `.meta.json` stamp (prompts SHA-256, child command) is written before the child starts: a matching
stamp means the file is reused when complete and resumed otherwise (the HF sampler keeps whole finished batches), on a
retry or a rerun alike; another stamp empties it. An attempt that finishes no new item counts against TRIES.
SIGTERM / SIGINT (the run cap's `timeout`) kill the child's process group before select_ckpt exits. Score them in
this process (objdel.eval.metrics: the device work runs in one process, D7), build the B2 signal panel and run G-vision
on the CPU. Eligible (D12) = gates.eligible: G-format, G-signal (groups of 8),
G-vision; an unmeasured G-vision fails (--vision is required, --no-vision makes every branch ineligible). Rule R-B
(prereg.select_rb) picks; none eligible -> the deepest is gated (§6.3). --items N: a seeded subset (SUBSET_SEED)
for selection only (D13 HF path).

The selected (or gated) branch is then judged on the full views: shifted and rescaled at K = 8 (gates), shifted at
K = 16 (report-only), G-prior (b) (--prepared), both probes on every probe record (--probes) and the rotated-vs-
identity attribution check (--rotated, report-only, §6.1). With --sampler vllm, `hf_check` then re-samples that
branch with HF on HF_CHECK_ITEMS shifted-val items; a vLLM - HF metric difference outside its 95 % CI fails the gates
(gate "G-hf-check", D13, owner 2026-11-05). Writes <out>/selected.json ("gates_pass" /
"gates_failed" are what g_probe --select reads), report.json and gates.json (prereg.entry); samples stay in
<out>/<branch>/<view>/ for re-scoring on the Mac.

--only topology (C2, PLAN §6.1 R3-R5): the topology probe alone on the run's single branch at K = --k into
<out>/report.json {"only", "branch", "sampler", "k", "probes": {"topology"}} (archive/phase2_runs/c2_analyze.py read
probes.topology).
No selection, other views, G-vision, prior or gates; --vision and --no-vision are refused, the other flags ignored.

References: PLAN D12 (R-B), D13, §5.1, §6.1 R1; https://docs.vllm.ai/en/stable/serving/offline_inference.html
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

from objdel.io import read_json
from objdel.policy import evaluate

TRIES = 2  # sampler child attempts that finish no new item


# sampler child timeout per item-sample, plus 600 s: 3x the measured time (EV-5), so a hung engine is killed in
# minutes, not hours. HF (~1 s per answer) resumes from its whole batches. vLLM writes once at the end, so a retry
# starts over: its time is scaled to answers at the 1,536-token cap (0.09 s at ~580 tokens in R0 v3 -> 0.25 s)
TIMEOUT_S = {"vllm": 0.75, "hf": 3.0}


HF_CHECK_ITEMS = 64  # D13 (owner 2026-11-05): HF re-samples the vLLM-judged branch on this seeded shifted-val subset


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
    """Refuse when this process initialised CUDA: the sampler child needs the whole GPU (D13)."""
    torch = sys.modules.get("torch")
    if torch is not None and torch.cuda.is_initialized():
        raise RuntimeError("CUDA is initialized in select_ckpt: the sampler child would share the GPU (D13)")


def sampler_cmd(kind: str, prompts: Path, model: Path, out: Path, extra=()) -> list[str]:
    """The sampler child's argv: objdel.train.sample in this venv (vllm_native: vLLM and objdel share it, ADR 0023)."""
    return [sys.executable, "-I", "-m", "objdel.train.sample", "--kind", kind, "--prompts", str(prompts), "--model",
            str(model), "--out", str(out), *extra]


def complete(comps: Path, want: Counter) -> bool:
    """Every prompt id has exactly its n completions (a cut or partial file is not complete)."""
    try:
        return Counter(json.loads(x)["id"] for x in comps.read_text().splitlines()) == want
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
    want = Counter({p["id"]: p["n"] for p in map(json.loads, prompts.read_text().splitlines())})
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

    from objdel.train import sample

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
    from objdel.eval import signal, summarize

    rows = evaluate.eval_rows(view, a.items)
    results, groups = evaluate.score_items(rows, sample_view(a, branch, rows, out, k), view.resolve())
    (out / "items.jsonl").write_text("".join(json.dumps(r) + "\n" for r in results))
    iou = {r["item_id"]: {n: float(np.mean(o["iou_tc"])) for n, o in r["objects"].items()} for r in results}
    return {"summary": summarize(results), "panel": signal.panel(groups), "iou": iou}


def prior_ci(model_iou: dict, prior_rows: list[dict]) -> dict:
    """G-prior (b): paired per-item (model - image-blind prior) IoU_tc per structure, nodule-cluster bootstrap."""
    from objdel.eval import signal

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
    from objdel.eval import probes

    d = a.probes / f"probe-{p}"
    done = sample_view(a, branch, evaluate.eval_rows(d, rl_only=False), out / f"probe-{p}", a.k)
    return getattr(probes, f"score_{p}")(list(done.values()), d) | {"ids": list(done)}  # ids pair the seeds (C2)


def hf_check(a, name: str, out: Path) -> dict:
    """D13 (owner 2026-11-05): with --sampler vllm, do vLLM's answers score like HF's on the trained model? HF samples
    branch `name` on HF_CHECK_ITEMS seeded shifted-val items at K; vLLM's answers are that branch's shifted samples
    (<out>/<name>/shifted/items.jsonl). Per metric, the per-item mean difference vLLM - HF with the nodule-cluster
    bootstrap 95 % CI (PLAN G-parity: every difference inside its CI); pass = 0 inside every CI."""
    from objdel.eval import signal

    rows = evaluate.eval_rows(a.shifted, HF_CHECK_ITEMS)
    hf = argparse.Namespace(**vars(a) | {"sampler": "hf"})
    done = sample_view(hf, a.run / "branches" / name, rows, out / name / "hf-check", a.k)
    got = {"hf": {r["item_id"]: r["samples"] for r in evaluate.score_items(rows, done, a.shifted.resolve())[0]}}
    lines = (out / name / "shifted" / "items.jsonl").read_text().splitlines()
    got["vllm"] = {r["item_id"]: r["samples"] for r in map(json.loads, lines)}
    ids = [i for i in got["hf"] if i in got["vllm"]]
    D = signal.point_or_draws(ids)
    diff = {k: signal.ci(np.array([np.mean(got["vllm"][i][k]) - np.mean(got["hf"][i][k]) for i in ids]), D)
            for k in HF_CHECK_KEYS}  # fmt: skip
    return {"branch": name, "items": len(ids), "k": a.k, "vllm_minus_hf": diff,
            "pass": len(ids) == len(rows) and all(lo <= 0 <= hi for _, lo, hi in diff.values())}


def vision_check(branch: Path, cfg: dict | None) -> dict:
    """G-vision on the CPU in float32 (B1 protocol); not measured (pass None, so it fails) under --no-vision."""
    if cfg is None:
        return {"pass": None, "measured": False}
    from transformers import AutoModelForImageTextToText, AutoProcessor

    from objdel.sft.train import vision_check as check

    model = AutoModelForImageTextToText.from_pretrained(branch, dtype="float32", device_map="cpu")
    return check(model, AutoProcessor.from_pretrained(branch), cfg)
