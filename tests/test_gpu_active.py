"""experiments.phase2.gpu_active (plan C0 step 4): the GPU-active share of a step from an nsys SQLite export."""

import sqlite3
import time

import numpy as np
import pytest

from experiments.phase2 import gpu_active


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
