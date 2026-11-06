"""GPU profile of the SFT step before the next training run: R1's settings and speed arms, one train_sft timing run
(STEPS optimizer steps, `branches: []`, no G-vision) per arm, each in its own process, torch.profiler on one step
after warm-up (train_sft `profile`), nvidia-smi sampled alongside.

    run.sh p0 experiments.phase2.profile_sft                  # on a Vast box ($OBJDEL_RUN_DIR/profile)
    python -I -m experiments.phase2.profile_sft --out DIR [--arms base no_ckpt]

Writes DIR/<arm>/{config.yaml, run/manifest.json, trace.json.gz, ops.txt, smi.csv, log.txt} and DIR/profile.json:
per arm the exit code, s/step and peak CUDA GiB (train_sft timing), and from the traced step its wall time, the GPU
busy share (union of kernel, copy and memset intervals / step time), kernel launches, CUDA runtime API time on the
CPU, host-sync calls (.item(), nonzero, stream syncs), the top kernels, and the nvidia-smi mean utilization and power. An arm that fails (out of memory) is recorded with
its exit code; the next arm runs.

Arms (owner question 2026-11-06: why is the GPU at ~60 % in R1; audit of the training path, same day):
- base: R1 (micro-batch 4 x 8, gradient checkpointing on).
- mb8, mb16: micro-batch 8 x 4, 16 x 2 (the CPU launches the same kernels per micro-batch at any batch size; TRL 1.16
  scores only label positions in chunks, so R1's "152k-vocab logits ~12 GB" limit does not apply).
- mb8_lean: mb8 without the per-micro-batch isnan/isinf sync (logging_nan_inf_filter) and with non-blocking copies.
- no_ckpt, no_ckpt_mb8: gradient checkpointing off (frozen base: no weight gradients, ~4N instead of ~6N FLOPs/token);
  the mb8 one may run out of memory (recorded).
- no_ckpt_compile: no_ckpt with torch.compile (Inductor; recompiles on new sequence lengths show in s/step).
"""

import argparse
import gzip
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R1 = ROOT / "experiments/phase2/configs/r1_qwen3vl2b.yaml"
STEPS = 16  # optimizer steps per arm: WARM (3) excluded from s/step, the profiler traces step 5
MB8, MB16 = ({"per_device_train_batch_size": b, "gradient_accumulation_steps": 32 // b} for b in (8, 16))
NO_CKPT = {"gradient_checkpointing": False}
ARMS = {
    "base": {},
    "mb8": MB8,
    "mb16": MB16,
    "mb8_lean": MB8 | {"logging_nan_inf_filter": False, "accelerator_config": {"non_blocking": True}},
    "no_ckpt": NO_CKPT,
    "no_ckpt_mb8": NO_CKPT | MB8,
    "no_ckpt_compile": NO_CKPT | {"torch_compile": True},
}
GPU_CATS = ("kernel", "gpu_memcpy", "gpu_memset")
SYNCS = ("aten::_local_scalar_dense", "aten::nonzero", "cudaStreamSynchronize", "cudaDeviceSynchronize")  # host waits
SMI = ["nvidia-smi", "--query-gpu=timestamp,utilization.gpu,power.draw,memory.used", "--format=csv,noheader,nounits",
       "-lms", "500"]


def arm_config(arm: str, out: Path) -> dict:
    """R1's config as a timing run of this arm: branches [], no G-vision, STEPS steps, profile to out/<arm>."""
    import yaml

    c = yaml.safe_load(R1.read_text())
    for k in ("sampler", "stages", "vision", "branches"):
        c.pop(k, None)
    return c | {"output_dir": str(out / arm / "run"), "branches": [], "vision": "off", "steps_per_epoch": STEPS,
                "rolling_minutes": 1e3, "profile": str(out / arm), "sft": c["sft"] | {"logging_steps": 1} | ARMS[arm]}


def busy(trace: dict) -> dict:
    """The traced optimizer step (ProfilerStep#N span): wall ms, GPU busy share (union of GPU intervals in the span),
    kernel launches, CUDA runtime API ms, host-sync calls by name, top kernels by total ms."""
    ev = [e for e in trace["traceEvents"] if e.get("ph") == "X"]
    steps = [e for e in ev if e.get("name", "").startswith("ProfilerStep#")]
    if not steps:
        raise ValueError("no ProfilerStep span in the trace")
    t0 = min(e["ts"] for e in steps)
    t1 = max(e["ts"] + e["dur"] for e in steps)
    gpu = sorted((max(e["ts"], t0), min(e["ts"] + e["dur"], t1)) for e in ev if e.get("cat") in GPU_CATS)
    covered, end = 0.0, t0
    for a, b in gpu:
        if b > end:
            covered += b - max(a, end)
            end = b
    kernels = [e for e in ev if e.get("cat") == "kernel"]
    top: dict[str, float] = {}
    for e in kernels:
        top[e["name"]] = top.get(e["name"], 0.0) + e["dur"] / 1e3
    api = sum(e["dur"] for e in ev if e.get("cat") == "cuda_runtime") / 1e3
    syncs = {n: sum(e.get("name") == n for e in ev) for n in SYNCS}
    return {"step_ms": (t1 - t0) / 1e3, "gpu_busy": covered / (t1 - t0), "kernels": len(kernels),
            "cuda_api_ms": api, "syncs": syncs, "top_kernels_ms": dict(sorted(top.items(), key=lambda x: -x[1])[:10])}


def smi_mean(path: Path) -> dict:
    rows = [r.split(", ") for r in path.read_text().splitlines() if r.count(",") == 3]
    vals = [(float(u), float(p)) for _, u, p, _ in rows if u.strip().isdigit()]
    if not vals:
        return {}
    return {"util_mean": sum(u for u, _ in vals) / len(vals), "power_w_mean": sum(p for _, p in vals) / len(vals)}


def run_arm(arm: str, out: Path) -> dict:
    import yaml

    d = out / arm
    d.mkdir(parents=True, exist_ok=True)
    (d / "config.yaml").write_text(yaml.safe_dump(arm_config(arm, out)))
    cmd = [sys.executable, "-I", "-m", "experiments.phase2.train_sft", "--config", str(d / "config.yaml")]
    t = time.time()
    with open(d / "smi.csv", "w") as smi_f, open(d / "log.txt", "w") as log:
        smi = subprocess.Popen(SMI, stdout=smi_f, stderr=subprocess.DEVNULL)
        try:
            code = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT).returncode
        finally:
            smi.terminate()
            smi.wait()
    res = {"exit": code, "wall_s": time.time() - t, "sft": ARMS[arm]}
    man = d / "run" / "manifest.json"
    if man.exists():
        tm = json.loads(man.read_text()).get("timing") or {}
        res |= {k: tm.get(k) for k in ("s_per_step", "samples_per_s", "max_memory_gib")}
    if (d / "trace.json.gz").exists():
        with gzip.open(d / "trace.json.gz", "rt") as f:
            res["traced_step"] = busy(json.load(f))
    return res | {"smi": smi_mean(d / "smi.csv")}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path(os.environ.get("OBJDEL_RUN_DIR", "runs")) / "profile")
    ap.add_argument("--arms", nargs="+", choices=list(ARMS), default=list(ARMS))
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    report = {}
    for arm in a.arms:
        report[arm] = run_arm(arm, a.out)
        (a.out / "profile.json").write_text(json.dumps(report, indent=1))
        print(arm, json.dumps({k: report[arm].get(k) for k in ("exit", "s_per_step", "max_memory_gib")}), flush=True)
    return 0 if all(r["exit"] == 0 for r in report.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
