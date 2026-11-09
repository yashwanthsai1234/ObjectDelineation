"""GPU-active share of a step (plan C0 step 4): the length of the union of every device busy interval (kernels,
memcpy, memset and whole CUDA-graph runs, all streams) inside [start, end), over end - start, from the SQLite file an
nsys export writes. vLLM decodes inside CUDA graphs, which nsys records as one graph-trace row per run; the graph
table is found by name (it contains GRAPH) because no nsys schema page names it.

    python -I -m experiments.phase2.gpu_active TRACE.sqlite START_NS END_NS
"""

import argparse
import json
import sqlite3

import numpy as np

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


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("sqlite")
    ap.add_argument("start", type=int)
    ap.add_argument("end", type=int)
    a = ap.parse_args(argv)
    print(json.dumps(from_sqlite(a.sqlite, a.start, a.end)))


if __name__ == "__main__":
    main()
