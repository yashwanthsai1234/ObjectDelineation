"""tools.profiling.gpu_active (plan C0 step 4): the GPU-active share of a step from an nsys SQLite export."""

import json
import sqlite3
import time

import numpy as np
import pytest

from tools.profiling import gpu_active


def test_share_is_the_union_of_busy_intervals_over_the_range():
    """U-GA1."""
    iv = [(0, 40), (30, 60), (50, 80), (90, 100), (100, 110)]
    assert gpu_active.share(iv, 0, 200) == 0.5
    assert gpu_active.share(iv, 10, 50) == 1.0
    assert gpu_active.share(iv, 85, 95) == 0.5
    assert gpu_active.share([], 0, 200) == 0.0
    assert gpu_active.share([(10, 20), (0, 10), (0, 10), (2, 5), (30, 30), (40, 35)], 0, 100) == 0.2
    assert gpu_active.share([(-50, 10), (190, 300)], 0, 200) == 0.1
    with pytest.raises(ValueError, match="zero-length"):
        gpu_active.share(iv, 5, 5)


def test_a_million_intervals_in_under_two_seconds():
    iv = np.arange(1_000_000)[:, None] * 3 + [0, 2]  # an array, as from_sqlite's rows become one
    t = time.perf_counter()
    assert gpu_active.share(iv, 0, 3_000_000) == pytest.approx(2 / 3)
    assert time.perf_counter() - t < 2


def export(path, graph=True, bad_graph=False):
    db = sqlite3.connect(path)
    for t in ("CUPTI_ACTIVITY_KIND_KERNEL", "CUPTI_ACTIVITY_KIND_MEMCPY", "CUPTI_ACTIVITY_KIND_MEMSET"):
        db.execute(f"CREATE TABLE {t} (start INT, end INT, streamId INT)")
    db.executemany("INSERT INTO CUPTI_ACTIVITY_KIND_KERNEL VALUES (?, ?, 7)", [(0, 40), (60, 80)])
    db.execute("CREATE TABLE NVTX_EVENTS (start INT, end INT, text TEXT)")  # ignored
    db.execute("INSERT INTO NVTX_EVENTS VALUES (0, 200, 'step')")
    if graph:
        db.execute("CREATE TABLE CUPTI_ACTIVITY_KIND_GRAPH_TRACE (start INT, end INT, graphId INT)")
        db.execute("INSERT INTO CUPTI_ACTIVITY_KIND_GRAPH_TRACE VALUES (100, 160, 1)")
    if bad_graph:
        db.execute("CREATE TABLE GRAPH_NODES (id INT)")
    db.commit()
    db.close()
    return path


def test_graph_runs_count_as_busy(tmp_path):
    """U-GA3: a table named like the graph-trace export with start/end counts; one without them is named, ignored."""
    out = gpu_active.from_sqlite(export(tmp_path / "a.sqlite", bad_graph=True), 0, 200)
    assert out == {"share": 0.6, "graph_table": ["CUPTI_ACTIVITY_KIND_GRAPH_TRACE"], "ignored_tables": ["GRAPH_NODES"]}
    out = gpu_active.from_sqlite(export(tmp_path / "b.sqlite", graph=False), 0, 200)
    assert out == {"share": 0.3, "graph_table": None, "ignored_tables": []}


# the same busy intervals as nsys rows (ns) and as torch.profiler events (us)
NS = [(0, 40), (30, 60), (60, 80), (100, 160)]  # union 0-80 and 100-160 = 140 of 200


def write_csv(path, rows=NS, head="Start (ns),Duration (ns),CorrId,Strm,Name"):
    path.write_text(head + "\n" + "".join(f"{s},{e - s},1,7,k\n" for s, e in rows))
    return path


def write_chrome(path, rows=NS, extra=()):
    ev = [{"ph": "X", "cat": c, "name": "k", "ts": s / 1000, "dur": (e - s) / 1000, "tid": 7}
          for (s, e), c in zip(rows, ["kernel", "gpu_memcpy", "gpu_memset", "kernel"], strict=True)]
    ev += [{"ph": "X", "cat": "cpu_op", "name": "aten::add", "ts": 0, "dur": 200},  # CPU distractors
           {"ph": "X", "cat": "cuda_runtime", "name": "cudaLaunchKernel", "ts": 0, "dur": 200},
           {"ph": "i", "cat": "kernel", "name": "no duration", "ts": 5}, *extra]
    path.write_text(json.dumps({"traceEvents": ev}))
    return path


def test_csv_equals_chrome_trace(tmp_path):
    """U-GA2: nsys cuda_gpu_trace (ns) and the torch.profiler Chrome trace (us) give one share."""
    a = gpu_active.from_csv(write_csv(tmp_path / "a.csv"), 0, 200)
    b = gpu_active.from_chrome(write_chrome(tmp_path / "a.json"), 0, 200)
    assert a == b == {"share": 0.7}
    assert gpu_active.read(tmp_path / "a.csv", 0, 200) == a and gpu_active.read(tmp_path / "a.json", 0, 200) == b
    ev = [{"ph": "X", "cat": "kernel", "name": "k", "ts": 0.5, "dur": 0.25}]  # binary-fraction us: exact in ns
    (tmp_path / "f.json").write_text(json.dumps(ev))  # the bare list form
    assert gpu_active.from_chrome(tmp_path / "f.json", 0, 1000) == {"share": 0.25}


def test_a_bad_trace_file_raises_a_clear_error(tmp_path):
    for name, text, msg in [("empty.csv", "", "empty"), ("bad.csv", "a,b\n1,2\n", "Start \\(ns\\)"),
                            ("empty.json", "", "empty"), ("bad.json", "{not json", "JSON"),
                            ("none.json", '{"x": 1}', "traceEvents"), ("t.txt", "x", "unknown")]:
        (tmp_path / name).write_text(text)
        with pytest.raises(ValueError, match=msg):
            gpu_active.read(tmp_path / name, 0, 200)


def test_main_exits_non_zero_on_a_bad_file(tmp_path, capsys):
    (tmp_path / "bad.csv").write_text("a,b\n")
    with pytest.raises(SystemExit) as e:
        gpu_active.main([str(tmp_path / "bad.csv"), "0", "200"])
    assert e.value.code != 0 and "Start (ns)" in str(e.value.code)
    gpu_active.main([str(write_csv(tmp_path / "ok.csv")), "0", "200"])
    assert json.loads(capsys.readouterr().out) == {"share": 0.7}
