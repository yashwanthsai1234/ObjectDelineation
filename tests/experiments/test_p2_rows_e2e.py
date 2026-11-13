"""objdel.cli.data.rl_rows end to end (Implementation.md B.9 E03, E06, E07, E08) on the test_data_quality
fixture store, through the `run` fixture of test_p2_rows (one 13-step block of 2 rows)."""
import json
import multiprocessing
import os
import shutil
import signal
import time
from dataclasses import replace

import pytest
from test_data_quality import (
    world,  # noqa: F401 (the module fixture)
    write_store,  # noqa: F401
)
from test_p2_rows import run  # noqa: F401 (the fixture)

from objdel.core.hashing import file_sha256
from objdel.data import prepare as data_prepare
from objdel.data import quality as Q


def _files(root):
    """Relative path -> SHA-256 of every file of a rows root but manifest.json (the root's name is in it)."""
    return {p.relative_to(root).as_posix(): file_sha256(p)
            for p in sorted(root.rglob("*")) if p.is_file() and p.name != "manifest.json"}


def _ids_sha(root):
    return json.loads((root / "manifest.json").read_text())["schedule"]["ids_sha256"]


def _die(job):
    """A worker that is SIGKILLed mid-build (E06)."""
    os.kill(os.getpid(), signal.SIGKILL)


# ------------------------------------------------------------------------------------------- E03 reproducible
def test_e03_same_args_and_procs_1_vs_4_give_the_same_root(run, tmp_path):  # noqa: F811
    a, b, c = tmp_path / "a", tmp_path / "b", tmp_path / "c"
    assert run.go("--out", str(a)) == 0 and run.go("--out", str(b)) == 0
    assert run.go("--out", str(c), "--procs", "4") == 0
    assert (a / "lidc/train.jsonl").read_bytes() == (b / "lidc/train.jsonl").read_bytes()
    assert _files(a) == _files(b) == _files(c)  # train.jsonl bytes, every image and cache
    assert _ids_sha(a) == _ids_sha(b) == _ids_sha(c)


def test_e03_another_seed_differs(run, tmp_path):  # noqa: F811
    a, b = tmp_path / "a", tmp_path / "b"
    run.go("--out", str(a))
    run.go("--out", str(b), "--seed", "7")  # a seed whose 26-row pool passes DG7-DG9
    assert (a / "lidc/train.jsonl").read_bytes() != (b / "lidc/train.jsonl").read_bytes()


def test_e03_an_interrupted_run_rerun_into_the_same_out_equals_a_clean_run(run, tmp_path):  # noqa: F811
    clean = tmp_path / "clean"
    assert run.go("--out", str(clean)) == 0
    shutil.copytree(clean, run.out)  # an interrupted build: no manifest, a truncated image, a stale file
    (run.out / "manifest.json").unlink()
    png = next((run.out / "images").rglob("*.png"))
    png.write_bytes(png.read_bytes()[:10])
    (run.out / "images" / "stale.png").write_bytes(b"x")
    assert run.go() == 0
    assert _files(run.out) == _files(clean)
    assert _ids_sha(run.out) == _ids_sha(clean)


# ------------------------------------------------------------------------------------------ E06 pool safety
def test_e06_a_zero_byte_image_fails_naming_the_window(run, world, tmp_path):  # noqa: F811
    bad = world.sched[0].item_id
    st = write_store(tmp_path / "store", [replace(it, image=None) if it.item_id == bad else it for it in world.train])
    with pytest.raises(SystemExit, match=f"window {bad}~"):
        run.go("--store", str(st), "--procs", "4")
    assert not (run.out / "manifest.json").exists()
    assert multiprocessing.active_children() == []  # no orphan worker


def test_e06_a_killed_worker_ends_the_build_within_60_s(run, monkeypatch):  # noqa: F811
    monkeypatch.setattr(run.R, "_record", _die)
    t0 = time.monotonic()
    with pytest.raises(SystemExit, match="BrokenProcessPool"):
        run.go("--procs", "4")
    assert time.monotonic() - t0 < 60
    assert not (run.out / "manifest.json").exists()
    assert multiprocessing.active_children() == []


# ------------------------------------------------------------------------------------------ E07 path safety
@pytest.mark.parametrize("bad", ["../x", "lidc/../../x", "/tmp/x", "lidc/a\0b"])
def test_e07_unsafe_ids_are_refused_before_any_write(world, tmp_path, bad):  # noqa: F811
    from objdel.data.targets.frame import frame
    it = world.train[0]
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="unsafe id"):
        data_prepare._prepare((str(out), bad, "lidc", it.group_id, it.raters, frame(it)))
    assert not tmp_path.joinpath("x").exists() and not any(out.rglob("*"))


def test_e07_ids_that_differ_only_by_case_are_refused(run, world, tmp_path):  # noqa: F811
    w = world.sched[0]
    twins = [w, replace(w, item_id=w.item_id.upper())]
    with pytest.raises(SystemExit, match="only by case"):
        run.R.materialize(world.store, tmp_path / "out", twins, 1)
    assert not (tmp_path / "out").exists()


def test_e07_every_written_path_resolves_under_out(run):  # noqa: F811
    assert run.go() == 0
    root = run.out.resolve()
    paths = [p.resolve() for p in run.out.rglob("*")]
    assert paths and all(p.is_relative_to(root) for p in paths)
    assert len({p.as_posix().casefold() for p in paths}) == len(paths)
    assert Q.records(run.out)
