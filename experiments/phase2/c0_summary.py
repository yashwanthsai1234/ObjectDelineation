"""C0 steps 5 and 6 (Mac, offline): c0_summary.json, the Phase C baseline per captured step, from the exports
nsys_baseline reads (no second capture):

    python -I -m experiments.phase2.c0_summary --stats RD/nsys/c0_base_stats --sqlite RD/nsys/c0_base.sqlite
        --log RD/c0_base/rl.json --window 16 17 --lengths LENGTHS.json --vision-kernel REGEX
        --out RD/nsys/c0_summary.json

Per step (FIELDS and SYNC_FIELDS; check() refuses a document missing one):
- wall_s: seconds by NVTX range name (nsys_baseline's per_range_s); gpu_active: per range name, the union of the
  device busy rows (gpu_active.sqlite_rows: kernels, memcpy, memset and graph runs) inside its ranges over their
  summed length ("step": the whole step);
- syncs_by_range: host syncs (nsys_baseline.is_sync) starting inside each range name, on the range's own thread;
- syncs: both counts of C0 step 1, torch's (the rl.json perf/host_syncs/<phase> keys and perf/host_syncs_vllm of
  the step; empty and None without --count-syncs) and nsys's: nsys_trainer_thread (calls on the `step` range's
  thread outside the vLLM ranges; the reward thread and the checker are not counted) and nsys_vllm (calls inside a
  vLLM range, any thread);
- h2d_bytes, d2h_bytes: cuda_gpu_trace memcpy bytes (Bytes (MB) x 1e6), "vllm" when the copy starts inside a vLLM
  range, else "trainer";
- 6a input_processing_s: from the step's first `generate` start to the first kernel inside it whose name matches
  --vision-kernel (the vision encoder's; None when none does);
- 6b decode_tail_s: the time with fewer than 8 sequences decoding, 4.3 ms x (longest answer - 8th-longest answer)
  from the step's answer lengths in --lengths ({"<step>": [tokens of each answer]}); check it against the trace's
  `generate/tail` range in wall_s;
- 6c generate_minus_handovers_s (each `generate` range less the sync, wake and sleep ranges inside it on its
  thread) and longest_answer.
Top level: window and graph_tables (the GRAPH tables gpu_active found by name).
"""

import argparse
import json
import re
import sys
from pathlib import Path

from experiments.phase2 import gpu_active
from experiments.phase2 import nsys_baseline as nb
from experiments.phase2.profile import covered
from objdel.core.io import read_json

SECONDS_PER_TOKEN = 0.0043  # decode seconds per token of the step's batch (B_gpu_audit.md:40, as runaway_cut)
TAIL_BELOW = 8  # the decode tail: fewer than this many sequences still decoding
HANDOVERS = tuple(n for n in nb.VLLM_RANGES if n != "generate")
H2D, D2H = re.compile(r"memcpy (Host-to-Device|HtoD)"), re.compile(r"memcpy (Device-to-Host|DtoH)")
API, GPU = ("Start(ns)", "Name", "Tid"), ("Start(ns)", "Bytes(MB)", "Name")
FIELDS = ("step", "wall_s", "gpu_active", "syncs_by_range", "syncs", "h2d_bytes", "d2h_bytes", "input_processing_s",
          "decode_tail_s", "generate_minus_handovers_s", "longest_answer")  # fmt: skip
SYNC_FIELDS = ("torch", "torch_vllm", "nsys_trainer_thread", "nsys_vllm")


def decode_tail_s(lengths) -> float:
    """4.3 ms x (longest - 8th-longest) of one step's answer lengths (C0 step 6b)."""
    top = sorted(lengths, reverse=True)
    if len(top) < TAIL_BELOW:
        raise ValueError(f"{len(top)} answers: the decode tail needs at least {TAIL_BELOW}")
    return SECONDS_PER_TOKEN * (top[0] - top[TAIL_BELOW - 1])


def summary(stats: str, db, log: list[dict], window: tuple[int, int], lengths: dict, vision: str) -> dict:
    base = nb.baseline(stats, db, log, window)
    ranges = [(int(r["Start(ns)"]), int(r["End(ns)"]), r["Name"], r["TID"]) for r in nb.read(stats, "nvtx_pushpop_trace")]
    syncs = [(int(r["Start(ns)"]), r["Tid"]) for r in nb.read(stats, "cuda_api_trace", API)
             if nb.is_sync(nb.api_name(r["Name"]))]  # fmt: skip
    gpu = [(int(r["Start(ns)"]), r["Name"], r["Bytes(MB)"]) for r in nb.read(stats, "cuda_gpu_trace", GPU)]
    busy, graph_tables, _ = gpu_active.sqlite_rows(db)
    steps = sorted(r for r in ranges if r[2] == "step" and r[1] > r[0])  # nsys_baseline's steps, in its order
    out = []
    for row, (t0, t1, _, tid) in zip(base["steps"], steps, strict=True):
        n = row["step"]
        if str(n) not in lengths:
            raise ValueError(f"no answer lengths for step {n}")
        lens = lengths[str(n)]
        inside = [r for r in ranges if t0 <= r[0] and r[1] <= t1]
        vllm = [(a, b) for a, b, name, _ in inside if name in nb.VLLM_RANGES]

        def in_vllm(t, vllm=vllm):
            return any(a <= t < b for a, b in vllm)

        step_syncs = [(t, th) for t, th in syncs if t0 <= t < t1]
        gen = sorted(r for r in inside if r[2] == "generate")
        kernels = [t for t, name, _ in gpu if gen and gen[0][0] <= t < gen[0][1] and re.search(vision, name)]
        handover = sum(b - a for g in gen for a, b, name, th in inside
                       if name in HANDOVERS and th == g[3] and g[0] <= a and b <= g[1])  # fmt: skip
        logged = {k: v for h in log if h.get("step") == n for k, v in h.items()}
        out.append({
            "step": n, "wall_s": row["per_range_s"], "gpu_active": active(inside, busy),
            "syncs_by_range": {name: sum(th == r[3] and r[0] <= t < r[1] for r in inside if r[2] == name
                                         for t, th in step_syncs) for name in dict.fromkeys(r[2] for r in inside)},
            "syncs": {"torch": {k.removeprefix("perf/host_syncs/"): v for k, v in logged.items()
                                if k.startswith("perf/host_syncs/")},
                      "torch_vllm": logged.get("perf/host_syncs_vllm"),
                      "nsys_trainer_thread": sum(th == tid and not in_vllm(t) for t, th in step_syncs),
                      "nsys_vllm": sum(in_vllm(t) for t, _ in step_syncs)},
            "h2d_bytes": copied(gpu, H2D, t0, t1, in_vllm), "d2h_bytes": copied(gpu, D2H, t0, t1, in_vllm),
            "input_processing_s": (min(kernels) - gen[0][0]) / 1e9 if kernels else None,
            "decode_tail_s": decode_tail_s(lens),
            "generate_minus_handovers_s": (sum(b - a for a, b, *_ in gen) - handover) / 1e9,
            "longest_answer": max(lens),
        })  # fmt: skip
    return {"window": list(window), "graph_tables": graph_tables, "steps": out}


def active(inside: list[tuple], busy: list) -> dict:
    """{range name: busy length inside its ranges / their summed length} (zero-length ranges skipped)."""
    total: dict[str, int] = {}
    hit: dict[str, float] = {}
    for a, b, name, _ in inside:
        if b > a:
            total[name] = total.get(name, 0) + b - a
            hit[name] = hit.get(name, 0) + covered(busy, a, b)
    return {name: hit[name] / total[name] for name in total}


def copied(gpu: list[tuple], kind: re.Pattern, t0: int, t1: int, in_vllm) -> dict:
    """{"trainer": bytes, "vllm": bytes} of the memcpy rows of `kind` starting in [t0, t1)."""
    out = {"trainer": 0, "vllm": 0}
    for t, name, mb in gpu:
        if t0 <= t < t1 and kind.search(name):
            out["vllm" if in_vllm(t) else "trainer"] += round(float(mb) * 1e6)
    return out


def check(doc: dict) -> None:
    """Raise ValueError naming every missing field (U-MS5)."""
    miss = [k for k in ("window", "graph_tables", "steps") if k not in doc]
    for i, s in enumerate(doc.get("steps", [])):
        lost = [f for f in FIELDS if f not in s] + [f"syncs.{f}" for f in SYNC_FIELDS if f not in s.get("syncs", {})]
        miss += [f"steps[{i}]: {f}" for f in lost]
    if miss:
        raise ValueError(f"c0_summary misses {', '.join(miss)}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stats", required=True, help="nsys stats --output prefix (PREFIX_<report>.csv)")
    ap.add_argument("--sqlite", type=Path, required=True, help="nsys export --type sqlite of the same capture")
    ap.add_argument("--log", type=Path, required=True, help="the arm's rl.json (time/step_s, perf/host_syncs/*)")
    ap.add_argument("--window", type=int, nargs=2, required=True, metavar=("S", "E"))
    ap.add_argument("--lengths", type=Path, required=True, help='{"<step>": [answer tokens]} of the window steps')
    ap.add_argument("--vision-kernel", required=True, help="regex of the vision encoder's first kernel name (6a)")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    try:
        doc = summary(a.stats, a.sqlite, read_json(a.log)["log"], tuple(a.window), read_json(a.lengths),
                      a.vision_kernel)  # fmt: skip
        check(doc)
    except (OSError, ValueError, KeyError) as e:
        sys.exit(f"c0_summary: {e}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(doc, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
