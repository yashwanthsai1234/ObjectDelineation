"""A6: tools/profiling/nsys_capture.sh (S04: pins, one capture flag, a stub nsys run) and
tools.profiling.nsys_baseline on the hand-made export in tests/goldens/nsys_small (I12-I15; the SQLite file is
built here from gpu_rows.json, so no binary is committed). G06 (the pulled V2 capture) is in
tests/experiments/test_a_gpu.py."""

import json
import os
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tools/profiling/nsys_capture.sh"


# S04 -------------------------------------------------------------------------------------------------------------


def test_capture_script_pins_and_flags():
    s = SCRIPT.read_text()
    assert subprocess.run(["bash", "-n", str(SCRIPT)]).returncode == 0
    for text in ("https://developer.download.nvidia.com/compute/cuda/repos/debian13/x86_64/$DEB",
                 "nsight-systems-2026.5.2_2026.5.2.323-1_amd64.deb",
                 "24954cc4d4cb5f051dc91c6a44101d7ccf15c802ab8ea8a32c09dd22542f9031",
                 "e71440978244e3883fc769da91a44a03451129ce21988d0db656a6bf7c1cb56b",
                 "sha256sum -c", "dpkg-deb -x", 'timeout --signal=INT --kill-after=120 "${NSYS_TIMEOUT:-0}"',
                 "OBJDEL_NVTX=1", '--profile-window "$S,$E"', "--capture-range=cudaProfilerApi", "exit 3", "exit 4",
                 "nvtx_pushpop_trace,nvtx_gpu_proj_trace,cuda_gpu_trace,cuda_api_trace", "pkill"):  # fmt: skip
        assert text in s, text
    assert "--cuda-memory-usage" not in s and "--cuda-graph-trace=node" not in s
    assert "apt-get install" not in s  # one source per package: the pinned .deb files
    hits = [p for p in (ROOT / "experiments").rglob("*") if p.is_file() and p.suffix in (".py", ".sh", ".yaml")
            and "--profile-steps" in p.read_text(errors="ignore")]  # fmt: skip
    assert hits == []


@pytest.mark.skipif(shutil.which("shellcheck") is None, reason="shellcheck not installed")
def test_capture_script_passes_shellcheck():
    out = subprocess.run(["shellcheck", str(SCRIPT)], capture_output=True, text=True)
    assert out.returncode == 0, out.stdout


@pytest.mark.parametrize(("timeout", "want"), [(None, "0"), ("1800", "1800")])
def test_capture_script_runs_the_command_under_nsys_with_the_window(tmp_path, timeout, want):
    """A stub nsys in NSYS_HOME logs its arguments; a stub `timeout` (none on the Mac) logs its duration and drops
    its own flags. No default cut: unset NSYS_TIMEOUT gives GNU timeout's 0 (no limit), so a whole-run capture
    (--profile-window 1,N) is not killed; a caller that wants a limit sets NSYS_TIMEOUT."""
    home, bin_ = tmp_path / "nsys" / "target-linux-x64", tmp_path / "bin"
    home.mkdir(parents=True)
    bin_.mkdir()
    log, dur = tmp_path / "calls.txt", tmp_path / "timeout.txt"
    (home / "nsys").write_text(f'#!/bin/sh\necho "$@" >> {log}\n')
    (bin_ / "timeout").write_text(f'#!/bin/sh\necho "$3" > {dur}\nshift 3\nexec "$@"\n')
    for f in (home / "nsys", bin_ / "timeout"):
        f.chmod(0o755)
    env = {k: v for k, v in os.environ.items() if k != "NSYS_TIMEOUT"} | {
        "NSYS_HOME": str(tmp_path / "nsys"), "PATH": f"{bin_}:{os.environ['PATH']}", "PY": sys.executable}
    env |= {"NSYS_TIMEOUT": timeout} if timeout else {}
    out = subprocess.run(["bash", str(SCRIPT), str(tmp_path / "rd"), "x", "16", "17", "--", "echo", "hi"], env=env,
                         capture_output=True, text=True)  # fmt: skip
    assert out.returncode == 0, out.stderr
    calls = log.read_text().splitlines()
    assert calls[:2] == ["--version", "status -e"]
    profile = next(c for c in calls if c.startswith("profile"))
    assert profile.endswith("echo hi --profile-window 16,17") and "--capture-range=cudaProfilerApi" in profile
    assert any(c.startswith("export --type sqlite") for c in calls) and any(c.startswith("stats") for c in calls)
    assert (tmp_path / "rd/nsys/x_tables.txt").exists() and dur.read_text().strip() == want


# I12-I15: nsys_baseline on the hand-made export ------------------------------------------------------------------

GOLD = ROOT / "tests/goldens/nsys_small"
REPORTS = ("nvtx_pushpop_trace", "nvtx_gpu_proj_trace", "cuda_api_trace", "cuda_gpu_trace")
HEADER = {r: (GOLD / f"stats_{r}.csv").read_text().splitlines()[0] for r in REPORTS}


def export(tmp_path, csv=None, rows=None, drop=()):
    """(stats prefix, sqlite path) of a copy of the golden export: `csv` replaces a report's data rows (header
    kept), `rows` the GPU tables of gpu_rows.json, `drop` leaves tables out of the SQLite file."""
    for r in REPORTS:
        text = (GOLD / f"stats_{r}.csv").read_text()
        if csv and r in csv:
            text = "\n".join([HEADER[r], *csv[r]]) + "\n"
        (tmp_path / f"stats_{r}.csv").write_text(text)
    db = tmp_path / "x.sqlite"
    con = sqlite3.connect(db)
    for table, iv in (rows or json.loads((GOLD / "gpu_rows.json").read_text())).items():
        if table not in drop:
            con.execute(f"create table {table} (start integer, end integer, deviceId integer)")
            con.executemany(f"insert into {table} values (?, ?, 0)", iv)
    con.commit()
    con.close()
    return str(tmp_path / "stats"), db


def log(times):
    """rl.json log: one entry per step 1..n, with time/step_s from step 2 on (the first log has none)."""
    return [{"step": i, "loss": 0.0} | ({"time/step_s": t} if i > 1 else {}) for i, t in enumerate(times, 1)]


def run(tmp_path, csv=None, rows=None, drop=(), window=(16, 17)):
    from tools.profiling import nsys_baseline as NB

    stats, db = export(tmp_path, csv, rows, drop)
    return NB.baseline(stats, db, log([10.0] * 17), window)


def test_hand_values(tmp_path):
    out = run(tmp_path)
    s1, s2 = out["steps"]
    assert (s1["step"], s2["step"]) == (16, 17)
    assert s1["coverage"] == 0.85 and s2["coverage"] == 0.5
    assert s1["gpu_busy"] == 0.55 and s1["gpu_idle"] == 1 - 0.55
    assert s2["gpu_busy"] == 0.4 and s2["gpu_busy_no_graph"] == 0.2  # the graph run [150, 170] counts as busy
    assert s1["graph_rows"] == s2["graph_rows"] == "present" and not s2["gpu_busy_is_lower_bound"]
    assert s1["host_syncs"] == {"cudaStreamSynchronize": 2, "cudaDeviceSynchronize": 1, "cudaEventSynchronize": 1,
                                "cudaMemcpy": 1, "cudaMemcpyToSymbol": 1}  # fmt: skip
    assert s1["host_sync_parts"] == {"vllm": 2, "trainer": 4}  # cudaMemcpy in generate, cudaEventSynchronize in sync
    assert s2["host_syncs"] == {"cuCtxSynchronize": 1, "cuStreamSynchronize": 1}
    assert s2["host_sync_parts"] == {"vllm": 1, "trainer": 1}
    assert (s1["driver_api_calls"], s2["driver_api_calls"]) == (0, 4)
    assert s1["per_range_s"] == pytest.approx({"step": 100e-9, "generate": 40e-9, "batch": 25e-9, "sync": 15e-9,
                                               "rewards": 30e-9})  # fmt: skip
    assert s2["per_range_s"] == pytest.approx({"step": 100e-9, "forward": 50e-9, "generate": 40e-9})
    assert s1["gpu_proj_s"] == pytest.approx({"step": 90e-9, "generate": 30e-9, "rewards": 10e-9})
    assert s2["gpu_proj_s"] == pytest.approx({"forward": 20e-9})
    assert out["graph_rows_in_gpu_trace"] is False and out["overhead"]["overhead_affected"] is False


def test_without_a_graph_table_busy_is_a_lower_bound(tmp_path):
    s1, s2 = run(tmp_path, drop=("CUPTI_ACTIVITY_KIND_GRAPH_TRACE",))["steps"]
    assert s2["graph_rows"] == "absent" and s2["gpu_busy_is_lower_bound"] is True and s2["gpu_busy"] == 0.2
    assert s1["graph_rows"] == "none" and s1["gpu_busy_is_lower_bound"] is False  # no graph launch in step 16


def test_graph_runs_listed_in_the_gpu_trace_are_seen(tmp_path):
    extra = (GOLD / "stats_cuda_gpu_trace.csv").read_text().splitlines()[1:] + ["150,20,12,,,,,,,,,,,,,,GPU,1,,7,g"]
    assert run(tmp_path, csv={"cuda_gpu_trace": extra})["graph_rows_in_gpu_trace"] is True


def ranges(*rows):
    """pushpop CSV rows from (start, end, name, tid)."""
    return [f"{a},{b},{b - a},0,0,{n},7,{t},0,0,0,,," for a, b, n, t in rows]


EMPTY = {"cuda_api_trace": [], "nvtx_gpu_proj_trace": [], "cuda_gpu_trace": []}


@pytest.mark.parametrize("rows,want", [
    ([(0, 100, "step", 1), (60, 100, "log", 1), (0, 10, "loader", 1)], 0.5),  # ends exactly at the step end
    ([(0, 100, "step", 1), (90, 110, "micro", 1)], 0.0),  # straddles the end: excluded (pin)
    ([(0, 100, "step", 1), (10, 20, "micro", 1), (10, 20, "micro", 1)], 0.2),  # identical siblings both count
    ([(0, 100, "step", 1), (0, 100, "micro", 1), (0, 100, "loader", 1)], 2.0),  # not clamped (pin)
    ([(0, 100, "step", 1), (0, 100, "micro", 2)], 0.0),  # another thread's range
    ([(0, 100, "step", 1)], 0.0),  # no children
])  # fmt: skip
def test_coverage_edges(tmp_path, rows, want):
    out = run(tmp_path, csv=EMPTY | {"nvtx_pushpop_trace": ranges(*rows)}, window=(16, 16))
    assert out["steps"][0]["coverage"] == want


def test_zero_length_steps_are_skipped_and_ns_become_seconds(tmp_path):
    rows = ranges((5, 5, "step", 1), (10, 1_000_000_010, "step", 1))
    out = run(tmp_path, csv=EMPTY | {"nvtx_pushpop_trace": rows}, window=(16, 16))
    assert [s["per_range_s"] for s in out["steps"]] == [{"step": 1.0}]


def test_empty_gpu_tables_mean_idle(tmp_path):
    rows = {t: [] for t in ("CUPTI_ACTIVITY_KIND_KERNEL", "CUPTI_ACTIVITY_KIND_MEMCPY", "CUPTI_ACTIVITY_KIND_MEMSET")}
    s = run(tmp_path, csv=EMPTY | {"nvtx_pushpop_trace": ranges((0, 100, "step", 1))}, rows=rows, window=(16, 16))
    assert (s["steps"][0]["gpu_busy"], s["steps"][0]["gpu_idle"]) == (0.0, 1.0)


def test_clear_errors(tmp_path):
    from tools.profiling import nsys_baseline as NB

    (tmp_path / "a").mkdir()
    with pytest.raises(ValueError, match="CUPTI_ACTIVITY_KIND_KERNEL"):
        run(tmp_path / "a", drop=("CUPTI_ACTIVITY_KIND_KERNEL",))
    stats, db = export(tmp_path)
    (tmp_path / "rl.json").write_text(json.dumps({"log": log([10.0] * 17)}))
    (tmp_path / "stats_cuda_api_trace.csv").unlink()
    with pytest.raises(SystemExit, match="stats_cuda_api_trace.csv"):
        NB.main(["--stats", stats, "--sqlite", str(db), "--log", str(tmp_path / "rl.json"), "--window", "16", "17",
                 "--out", str(tmp_path / "o.json")])  # fmt: skip
    bad = HEADER["nvtx_pushpop_trace"].replace("Start (ns)", "Start (us)")
    (tmp_path / "stats_nvtx_pushpop_trace.csv").write_text(bad + "\n")
    with pytest.raises(ValueError, match=r"stats_nvtx_pushpop_trace.csv.*Start \(ns\)"):
        NB.baseline(stats, db, log([10.0] * 17), (16, 17))


def test_no_step_ranges_writes_an_empty_list_and_fails(tmp_path, capsys):
    from tools.profiling import nsys_baseline as NB

    stats, db = export(tmp_path, csv={"nvtx_pushpop_trace": ranges((0, 10, "micro", 1))})
    (tmp_path / "rl.json").write_text(json.dumps({"log": log([10.0] * 17)}))
    argv = ["--stats", stats, "--sqlite", str(db), "--log", str(tmp_path / "rl.json"), "--window", "16", "17",
            "--out", str(tmp_path / "o.json")]  # fmt: skip
    assert NB.main(argv) == 1
    assert json.loads((tmp_path / "o.json").read_text())["steps"] == [] and "no `step` range" in capsys.readouterr().err


def test_main_writes_the_baseline(tmp_path):
    from tools.profiling import nsys_baseline as NB

    stats, db = export(tmp_path)
    (tmp_path / "rl.json").write_text(json.dumps({"log": log([10.0] * 15 + [12.0, 12.0])}))
    argv = ["--stats", stats, "--sqlite", str(db), "--log", str(tmp_path / "rl.json"), "--window", "16", "17",
            "--out", str(tmp_path / "o.json")]  # fmt: skip
    assert NB.main(argv) == 0
    doc = json.loads((tmp_path / "o.json").read_text())
    assert len(doc["steps"]) == 2 and doc["overhead"]["overhead_affected"] is True


@pytest.mark.parametrize("win,want", [(12.0, True), (10.9, False), (11.0, False)])
def test_overhead_compares_the_window_with_the_three_steps_before(win, want):
    from tools.profiling import nsys_baseline as NB

    o = NB.overhead(log([10.0] * 15 + [win, win]), (16, 17))
    assert o["overhead_affected"] is want  # 11.0 > 11.000000000000002 is False: `win > 1.1 * base` in floats
    assert o["before_s"] == 10.0 and o["window_s"] == win


def test_overhead_needs_the_step_times():
    from tools.profiling import nsys_baseline as NB

    with pytest.raises(ValueError, match="time/step_s"):
        NB.overhead([{"step": i} for i in range(1, 18)], (16, 17))


def test_api_names_lose_only_a_trailing_version_suffix():
    from tools.profiling import nsys_baseline as NB

    assert [NB.api_name(n) for n in ("cudaMemcpy_v3020", "cuStreamSynchronize_v2", "cudaMemcpy_v2_foo")] == [
        "cudaMemcpy", "cuStreamSynchronize", "cudaMemcpy_v2_foo"]  # fmt: skip


def test_cli_runs_isolated(tmp_path):
    """S02: `python -I -m` ignores a planted json.py in the working directory."""
    (tmp_path / "json.py").write_text("raise SystemExit('planted json.py ran')\n")
    out = subprocess.run([sys.executable, "-I", "-m", "tools.profiling.nsys_baseline", "--help"], cwd=tmp_path,
                         capture_output=True, text=True, env=os.environ | {"PYTHONPATH": str(ROOT)})  # fmt: skip
    assert out.returncode == 0 and "--stats" in out.stdout, out.stderr
