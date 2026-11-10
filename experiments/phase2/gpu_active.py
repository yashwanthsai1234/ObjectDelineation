"""GPU-active share of a step (plan C0 step 4): the length of the union of every device busy interval (kernels,
memcpy, memset and whole CUDA-graph runs, all streams) inside [start, end), over end - start, from the SQLite file an
nsys export writes. vLLM decodes inside CUDA graphs, which nsys records as one graph-trace row per run; the graph
table is found by name (it contains GRAPH) because no nsys schema page names it.

    python -I -m experiments.phase2.gpu_active TRACE START_NS END_NS

TRACE is that SQLite file, an nsys `cuda_gpu_trace` CSV (`nsys stats --format csv`; every row is a busy interval, ns)
or the Chrome trace of torch.profiler (the nsys fallback; kernel, gpu_memcpy and gpu_memset events, us). A bad file
exits non-zero with the reason.
"""

import argparse
import csv
import json
import sqlite3
import sys
from pathlib import Path

import numpy as np

GPU_CATS = ("kernel", "gpu_memcpy", "gpu_memset")  # torch.profiler's categories of device work
TABLES = ("CUPTI_ACTIVITY_KIND_KERNEL", "CUPTI_ACTIVITY_KIND_MEMCPY", "CUPTI_ACTIVITY_KIND_MEMSET")


def share(intervals, start: int, end: int) -> float:
    """Union length of the (begin, end) intervals clipped to [start, end), over end - start; empty ones ignored."""
    if end <= start:
        raise ValueError(f"a zero-length range [{start}, {end}]")
    iv = np.clip(np.asarray(intervals, dtype=np.int64).reshape(-1, 2), start, end)
    iv = iv[iv[:, 1] > iv[:, 0]]
    iv = iv[np.argsort(iv[:, 0], kind="stable")]
    reach = np.maximum.accumulate(iv[:, 1])  # the furthest end so far: a new run starts past it
    covered = np.minimum(iv[:, 1], end) - np.maximum(iv[:, 0], np.concatenate([[start], reach[:-1]]))
    return float(np.clip(covered, 0, None).sum() / (end - start))


def from_sqlite(path, start: int, end: int) -> dict:
    """share() of the export's kernel, memcpy, memset and graph-run rows; names the graph tables used and the GRAPH
    tables without start/end columns (ignored)."""
    db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        names = [r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type = 'table'")]
        cols = {n: {r[1] for r in db.execute(f'PRAGMA table_info("{n}")')} for n in names}
        graph = [n for n in names if "GRAPH" in n.upper() and {"start", "end"} <= cols[n]]
        ignored = [n for n in names if "GRAPH" in n.upper() and n not in graph]
        rows = [r for n in [*TABLES, *graph] if n in cols for r in db.execute(f'SELECT start, end FROM "{n}"')]
    finally:
        db.close()
    return {"share": share(rows, start, end), "graph_table": graph or None, "ignored_tables": ignored}


def from_csv(path, start: int, end: int) -> dict:
    """share() of an nsys cuda_gpu_trace CSV: each row is [Start (ns), Start + Duration (ns)]."""
    with open(path, newline="") as f:
        head = next(csv.reader(f), None)
        if head is None:
            raise ValueError(f"{path}: empty file")
        if not {"Start (ns)", "Duration (ns)"} <= set(head):
            raise ValueError(f"{path}: no column Start (ns) and Duration (ns) in {head}")
        t0, d = head.index("Start (ns)"), head.index("Duration (ns)")
        rows = [(int(r[t0]), int(r[t0]) + int(r[d])) for r in csv.reader(f) if r]
    return {"share": share(rows, start, end)}


def from_chrome(path, start: int, end: int) -> dict:
    """share() of a torch.profiler Chrome trace: the complete events (ph X) of the device categories, us to ns."""
    text = Path(path).read_text()
    if not text.strip():
        raise ValueError(f"{path}: empty file")
    try:
        doc = json.loads(text)
    except json.JSONDecodeError as e:
        raise ValueError(f"{path}: not JSON ({e})") from e
    events = doc.get("traceEvents") if isinstance(doc, dict) else doc
    if not isinstance(events, list):
        raise ValueError(f"{path}: no traceEvents list")
    rows = [(round(e["ts"] * 1000), round((e["ts"] + e["dur"]) * 1000))
            for e in events if e.get("ph") == "X" and e.get("cat") in GPU_CATS]
    return {"share": share(rows, start, end)}


def read(path, start: int, end: int) -> dict:
    """The share of a SQLite export, a cuda_gpu_trace CSV or a Chrome trace, chosen by the file's suffix."""
    readers = {".sqlite": from_sqlite, ".csv": from_csv, ".json": from_chrome}
    if (reader := readers.get(Path(path).suffix)) is None:
        raise ValueError(f"{path}: unknown trace type (needs {', '.join(readers)})")
    return reader(path, start, end)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("trace")
    ap.add_argument("start", type=int)
    ap.add_argument("end", type=int)
    a = ap.parse_args(argv)
    try:
        print(json.dumps(read(a.trace, a.start, a.end)))
    except (ValueError, OSError, sqlite3.Error) as e:
        sys.exit(f"gpu_active: {e}")


if __name__ == "__main__":
    main()
