"""Phase C baseline from an nsys capture (plan A6 step 9; Mac, offline, stdlib + profile.covered):

    python -I -m tools.profiling.nsys_baseline --stats RD/nsys/c0_base_stats --sqlite RD/nsys/c0_base.sqlite
        --log RD/c0_base/rl.json --window 16 17 --out RD/nsys/nsys_baseline.json

Inputs are nsys_capture.sh's outputs: `nsys stats --output PREFIX` writes PREFIX_<report>.csv (time columns in ns)
for nvtx_pushpop_trace, nvtx_gpu_proj_trace, cuda_api_trace and cuda_gpu_trace; the SQLite export holds the GPU rows.
One entry per captured `step` range (NvtxSteps: the capture holds exactly steps S..E), numbered from S:
- coverage: the summed duration of the ranges directly inside `step` on its thread / the step's duration (nesting
  by time containment; a range that straddles the step is not inside; identical intervals are siblings; no clamp);
- gpu_busy: profile.covered of the kernel, memcpy and memset rows plus every GRAPH_TABLES table present (graph rows
  with start == end == 0 have no time and are dropped) / the step's duration; gpu_idle = 1 - gpu_busy. A graph row
  spans the whole graph run, so gaps between its nodes count as busy: on graph segments gpu_busy is an upper bound.
  gpu_busy_no_graph leaves the graph rows out. graph_rows: "present" (a GRAPH_TABLES table exists), "absent" (the
  step launches a CUDA graph but no table holds graph runs: gpu_busy_is_lower_bound, and Phase C must not read
  gpu_idle as its baseline; fix: add the table named in NAME_tables.txt to GRAPH_TABLES and re-run) or "none";
- host_syncs: calls by name (the `_v<digits>` suffix removed) of the stream, device/context and event syncs and the
  synchronous (non-Async) memcpys, runtime (cuda*) and driver (cu*) names, starting inside the step;
  host_sync_parts: "vllm" when the call starts inside a generate, sync, wake or sleep range of the step (any
  thread), else "trainer"; driver_api_calls: every cu* call in the step (does `-t cuda` list driver calls at all?);
- per_range_s: total duration by range name of the ranges inside the step (any thread); gpu_proj_s: projected GPU
  time by range name (nvtx_gpu_proj_trace, by the range's original start).
Top level: graph_rows_in_gpu_trace (the GPU trace CSV lists more rows in the window than the kernel, memcpy and
memset tables: it lists graph runs too) and overhead: time/step_s of the window steps against the 3 steps before,
overhead_affected when more than 10 % slower (the arm ends at E: there is no step after the window).
"""

import argparse
import csv
import json
import re
import sqlite3
import sys
from collections import Counter
from pathlib import Path

from objdel.core.io import read_json
from tools.profiling.profile import covered

GPU_TABLES = ("CUPTI_ACTIVITY_KIND_KERNEL", "CUPTI_ACTIVITY_KIND_MEMCPY", "CUPTI_ACTIVITY_KIND_MEMSET")
GRAPH_TABLES = ("CUPTI_ACTIVITY_KIND_GRAPH_TRACE",)  # the CUPTI record name; the nsys table name is unverified
VLLM_RANGES = ("generate", "sync", "wake", "wake_weights", "wake_kv", "sleep")
SYNCS = re.compile(r"(cuda|cu)(StreamSynchronize|DeviceSynchronize|CtxSynchronize|EventSynchronize|Memcpy\w*)$")
COLS = {  # the columns read, spaces removed (nsys has written both "Start (ns)" and "Start(ns)")
    "nvtx_pushpop_trace": ("Start(ns)", "End(ns)", "Name", "TID"),
    "nvtx_gpu_proj_trace": ("Name", "ProjectedDuration(ns)", "OrigStart(ns)"),
    "cuda_api_trace": ("Start(ns)", "Name"),
    "cuda_gpu_trace": ("Start(ns)",),
}


def read(prefix: str, report: str, cols: tuple = ()) -> list[dict]:
    """PREFIX_<report>.csv as rows keyed by the space-free column names `cols` (default COLS[report])."""
    cols, path = cols or COLS[report], Path(f"{prefix}_{report}.csv")
    with open(path, newline="") as f:  # a missing file raises FileNotFoundError naming it
        rows = csv.reader(f)
        head = [h.replace(" ", "") for h in next(rows, [])]
        miss = [c for c in cols if c not in head]
        if miss:
            raise ValueError(f"{path.name}: no column {', '.join(c.replace('(', ' (') for c in miss)} in {head}")
        idx = [head.index(c) for c in cols]
        return [{c: r[i] for c, i in zip(cols, idx, strict=True)} for r in rows if r]


def api_name(name: str) -> str:
    return re.sub(r"_v\d+$", "", name)


def is_sync(name: str) -> bool:
    return bool(SYNCS.match(name)) and "Async" not in name


def gpu_rows(db) -> tuple[list, list, bool]:
    """(kernel + memcpy + memset intervals, graph-run intervals, whether a GRAPH_TABLES table exists)."""
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    have = {n for (n,) in con.execute("select name from sqlite_master where type = 'table'")}
    miss = [t for t in GPU_TABLES if t not in have]
    if miss:
        raise ValueError(f"{db}: no table {', '.join(miss)} (tables: {sorted(have)})")
    plain = [iv for t in GPU_TABLES for iv in con.execute(f"select start, end from {t}")]
    graphs = [iv for t in GRAPH_TABLES if t in have for iv in con.execute(f"select start, end from {t}")
              if iv != (0, 0)]  # fmt: skip
    return plain, graphs, any(t in have for t in GRAPH_TABLES)


def coverage(step: tuple, ranges: list[tuple]) -> float:
    """The summed duration of the ranges directly inside `step` on its thread / the step's duration."""
    t0, t1, _, tid = step
    inside = sorted(((a, b) for a, b, n, t in ranges if t == tid and t0 <= a and b <= t1 and (a, b, n) != step[:3]),
                    key=lambda r: (r[0], -r[1]))  # fmt: skip
    total, end, top = 0, t0, None
    for a, b in inside:
        if a >= end or (a, b) == top:  # after the last top-level range, or its identical sibling
            total += b - a
            end, top = max(end, b), (a, b)
    return total / (t1 - t0)


def overhead(log: list[dict], window: tuple[int, int]) -> dict:
    s, e = window
    t = {h["step"]: h["time/step_s"] for h in log if "time/step_s" in h and "step" in h}
    need = [i for i in range(s - 3, e + 1) if i not in t]
    if need:
        raise ValueError(f"the log has no time/step_s for steps {need}")
    win = sum(t[i] for i in range(s, e + 1)) / (e - s + 1)
    base = sum(t[i] for i in range(s - 3, s)) / 3
    return {"window_s": win, "before_s": base, "overhead_affected": win > 1.1 * base}


def baseline(stats: str, db, log: list[dict], window: tuple[int, int]) -> dict:
    ranges = [(int(r["Start(ns)"]), int(r["End(ns)"]), r["Name"], r["TID"]) for r in read(stats, "nvtx_pushpop_trace")]
    proj = [(int(r["OrigStart(ns)"]), r["Name"], int(r["ProjectedDuration(ns)"])) for r in
            read(stats, "nvtx_gpu_proj_trace")]  # fmt: skip
    api = [(int(r["Start(ns)"]), api_name(r["Name"])) for r in read(stats, "cuda_api_trace")]
    traced = [int(r["Start(ns)"]) for r in read(stats, "cuda_gpu_trace")]
    plain, graphs, graph_table = gpu_rows(db)
    steps = sorted(r for r in ranges if r[2] == "step" and r[1] > r[0])  # a zero-length step has no time to share
    out = []
    for k, step in enumerate(steps):
        t0, t1 = step[:2]
        dur = t1 - t0
        inside = [r for r in ranges if t0 <= r[0] and r[1] <= t1]
        vllm = [(a, b) for a, b, n, _ in inside if n in VLLM_RANGES]
        calls = [(t, n) for t, n in api if t0 <= t < t1]
        syncs = [(t, n) for t, n in calls if is_sync(n)]
        launch = any(n in ("cudaGraphLaunch", "cuGraphLaunch") for _, n in calls)
        busy = covered(plain + graphs, t0, t1) / dur
        state = "present" if graph_table else "absent" if launch else "none"
        n_vllm = sum(any(a <= t < b for a, b in vllm) for t, _ in syncs)
        out.append({
            "step": window[0] + k, "start_ns": t0, "dur_s": dur / 1e9, "coverage": coverage(step, ranges),
            "gpu_busy": busy, "gpu_idle": 1 - busy, "gpu_busy_no_graph": covered(plain, t0, t1) / dur,
            "graph_rows": state, "gpu_busy_is_lower_bound": state == "absent",
            "host_syncs": dict(Counter(n for _, n in syncs)),
            "host_sync_parts": {"vllm": n_vllm, "trainer": len(syncs) - n_vllm},
            "driver_api_calls": sum(n.startswith("cu") and not n.startswith("cuda") for _, n in calls),
            "per_range_s": _sum((n, b - a) for a, b, n, _ in inside),
            "gpu_proj_s": _sum((n, d) for t, n, d in proj if t0 <= t < t1),
        })  # fmt: skip
    t0, t1 = (steps[0][0], steps[-1][1]) if steps else (0, 0)
    n_plain = sum(t0 <= a < t1 for a, _ in plain)
    return {"window": list(window), "steps": out,
            "graph_rows_in_gpu_trace": sum(t0 <= t < t1 for t in traced) > n_plain,
            "overhead": overhead(log, window)}  # fmt: skip


def _sum(pairs) -> dict:
    """{name: total ns / 1e9} in first-seen order."""
    out: dict[str, float] = {}
    for n, d in pairs:
        out[n] = out.get(n, 0) + d
    return {n: d / 1e9 for n, d in out.items()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stats", required=True, help="nsys stats --output prefix (PREFIX_<report>.csv)")
    ap.add_argument("--sqlite", type=Path, required=True, help="nsys export --type sqlite of the same capture")
    ap.add_argument("--log", type=Path, required=True, help="the arm's rl.json (time/step_s per step)")
    ap.add_argument("--window", type=int, nargs=2, required=True, metavar=("S", "E"))
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    try:
        doc = baseline(a.stats, a.sqlite, read_json(a.log)["log"], tuple(a.window))
    except (OSError, ValueError, KeyError) as e:
        sys.exit(f"nsys_baseline: {e}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(doc, indent=1) + "\n")
    if not doc["steps"]:
        print(f"nsys_baseline: no `step` range in {a.stats}_nvtx_pushpop_trace.csv (was OBJDEL_NVTX=1 and "
              "--profile-window set?)", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
