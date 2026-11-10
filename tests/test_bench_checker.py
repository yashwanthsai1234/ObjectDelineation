"""experiments/phase2/perf/bench_checker.py (Phase D, D1): the checker benchmark (I16), its replay arm without cv2 or
scipy (S10), and its inputs read only under OBJDEL_DATA_ROOT (S15). CPU build; reads no data/."""
import builtins
import hashlib
import json
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from experiments.phase2.perf import bench_checker as B
from objdel import _gpucheck
from objdel.verify.reward.cache import ExpertCache
from objdel.verify.reward.score import score_many

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "experiments/phase2/perf/bench_checker.py"
MODULE = "experiments.phase2.perf.bench_checker"


def _disc(n, r, cy, cx):
    y, x = np.mgrid[:n, :n]
    return np.hypot(y - cy, x - cx) <= r


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


@pytest.fixture()
def root(tmp_path):
    """A fake data root: checker_eval/masks.npz (8 masks), prepared/lidc/val.jsonl (4 rows) and their 4 caches."""
    T = tmp_path / "root"
    (T / "checker_eval").mkdir(parents=True)
    masks = [_disc(n, n / 3, n / 2, n / 2).astype(np.uint8) for n in (5, 8, 13, 21, 34, 1, 2, 40)]
    hw = np.array([m.shape for m in masks], np.int32)
    start = np.cumsum([0, *(m.size for m in masks)], dtype=np.int64)
    np.savez(T / "checker_eval/masks.npz", bits=np.concatenate([m.ravel() for m in masks]), hw=hw, start=start)
    rows = []
    for k in range(4):
        item = f"lidc/1.2.{k}/n000/z00{k}"
        c = ExpertCache.build(np.stack([_disc(40, r, 20, 20 + k) for r in (8, 9, 10)]))
        path = T / "prepared/caches" / item / "nodule.npz"
        path.parent.mkdir(parents=True)
        digest = c.save(path)
        text = '{"nodule":[{"outer":[[300,300],[700,300],[700,700],[300,700]],"holes":[]}]}'
        rows.append({"item_id": item, "rl": True, "W": 40, "H": 40, "caches": {"nodule": digest}, "targets": [text]})
    (T / "prepared/lidc").mkdir(parents=True)
    (T / "prepared/lidc/val.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows[::-1]))  # unsorted
    return T


def _run(args, cwd, env=None, check=True, py=sys.executable):
    p = subprocess.run([py, "-I", "-m", MODULE, *map(str, args)], cwd=cwd, capture_output=True, text=True,
                       env={**os.environ, **(env or {})}, timeout=300)
    assert not check or p.returncode == 0, p.stderr
    return p


def test_help(tmp_path):
    assert "--dump-replay" in _run(["--help"], tmp_path).stdout


def test_tiny_run_writes_build_device_inputs_times_and_stats(root, tmp_path):
    out = tmp_path / "out/b.json"
    _run(["--engine", "current", "--workload", "crops", "--workload", "lesion512_b1", "--stages", "--warmup", "1",
          "--repeats", "3", "--out", out], tmp_path, {"OBJDEL_DATA_ROOT": str(root)})
    r = json.loads(out.read_text())
    assert r["build"] == _gpucheck.build() and r["device"] and r["engine"] == "current"
    masks = root / "checker_eval/masks.npz"
    assert r["root"] == str(root) and r["inputs"] == {"checker_eval/masks.npz": {"path": str(masks),
                                                                                "sha256": _sha(masks)}}
    assert set(r["workloads"]) == {"crops", "lesion512_b1"}
    for w in r["workloads"].values():
        assert w["repeats"] == 3 and w["median_ms"] > 0 and w["iqr_ms"] >= 0
        assert w["stats"]["waits"] > 0 and w["stats"]["drain_syncs"] == 2  # one call's counters
        assert list(w["stages_ms"]) == ["ccl", "rings", "trace"]
    assert r["workloads"]["crops"]["masks"] == 8 and r["workloads"]["lesion512_b1"]["masks"] == 1


def test_out_is_required_and_old_engine_needs_its_build(tmp_path):
    assert _run(["--workload", "lesion512_b1"], tmp_path, check=False).returncode == 2  # argparse: --out missing
    p = _run(["--engine", "old", "--old-python", tmp_path / "none/python", "--workload", "lesion512_b1", "--out",
              tmp_path / "b.json"], tmp_path, check=False)
    assert p.returncode == 1 and "old engine" in p.stderr and str(tmp_path / "none/python") in p.stderr


def test_old_engine_child_gets_the_same_data_root(root, tmp_path):
    stub = tmp_path / "old_python"
    stub.write_text("#!/bin/sh\necho ARGS \"$@\"\nenv\n")
    stub.chmod(0o755)
    p = _run(["--engine", "old", "--old-python", stub, "--workload", "crops", "--out", tmp_path / "b.json"], tmp_path,
             {"OBJDEL_DATA_ROOT": str(root)})
    assert f"OBJDEL_DATA_ROOT={root}" in p.stdout.splitlines()
    args = next(ln for ln in p.stdout.splitlines() if ln.startswith("ARGS")).split()
    assert args[1:3] == ["-I", str(SCRIPT)] and "--engine current" in " ".join(args) and "--old-python" not in args


def test_inputs_only_under_the_data_root(root, tmp_path, monkeypatch):
    """S15: no data/ in the working folder; every input path starts with the root; a missing input names its path."""
    work, out = tmp_path / "work", tmp_path / "out"
    work.mkdir()
    env = {"OBJDEL_DATA_ROOT": str(root)}
    _run(["--dump-replay", out / "replay_inputs.npz"], work, env)
    meta = json.loads(str(np.load(out / "replay_inputs.npz")["meta"]))
    assert meta["root"] == str(root) and all(v["path"].startswith(str(root)) for v in meta["inputs"].values())
    assert len(meta["item_ids"]) == 4 and meta["item_ids"] == sorted(meta["item_ids"])
    _run(["--engine", "current", "--workload", "crops", "--warmup", "1", "--repeats", "1", "--out", out / "b.json"],
         work, env)
    r = json.loads((out / "b.json").read_text())
    assert r["root"] == str(root) and all(v["path"].startswith(str(root)) for v in r["inputs"].values())

    seen = []
    real_open, real_load = builtins.open, np.load
    monkeypatch.setattr(builtins, "open", lambda f, *a, **k: seen.append(str(f)) or real_open(f, *a, **k))
    monkeypatch.setattr(np, "load", lambda f, *a, **k: seen.append(str(f)) or real_load(f, *a, **k))
    monkeypatch.setenv("OBJDEL_DATA_ROOT", str(root))
    monkeypatch.chdir(work)
    B.main(["--workload", "crops", "--warmup", "1", "--repeats", "1", "--out", str(out / "c.json")])
    B.main(["--dump-replay", str(out / "r2.npz")])
    assert seen and all(s.startswith((str(root), str(out), str(REPO))) for s in seen), seen

    empty = tmp_path / "empty"
    empty.mkdir()
    p = _run(["--workload", "crops", "--out", out / "d.json"], work, {"OBJDEL_DATA_ROOT": str(empty)}, check=False)
    assert p.returncode == 1 and str(empty / "checker_eval/masks.npz") in p.stderr


def test_no_data_path_literal_in_the_script():
    assert not re.search(r"""["']data/""", SCRIPT.read_text())


# ------------------------------------------------------------------ S10: the replay arm needs numpy and _gpucheck only

HOOK = textwrap.dedent("""
    import sys, importlib.abc
    class Only(importlib.abc.MetaPathFinder):  # the old venv: the standard library, numpy and objdel._gpucheck
        def find_spec(self, name, path=None, target=None):
            top = name.split(".")[0]
            if top in sys.stdlib_module_names or top == "numpy" or name in ("objdel", "objdel._gpucheck"):
                return None
            raise ImportError(f"blocked in the old venv: {name}")
    sys.meta_path.insert(0, Only())
""")


def _replay_child(code, cwd, extra_path=(), flags=("-I",)):
    pre = f"import sys; sys.path[:0] = {list(map(str, extra_path))!r}\n" if extra_path else ""
    return subprocess.run([sys.executable, *flags, "-c", pre + HOOK + code], cwd=cwd, capture_output=True, text=True,
                          timeout=300)


def test_replay_arm_runs_without_cv2_or_scipy(tmp_path):
    npz, out = tmp_path / "replay_inputs.npz", tmp_path / "replay.json"
    _run(["--dump-replay", npz, "--synthetic", "8"], tmp_path)
    assert npz.is_file() and (tmp_path / "replay_inputs.npz.sha256").read_text().split()[0] == _sha(npz)
    argv = [str(SCRIPT), "--engine", "current", "--replay", str(npz), "--warmup", "1", "--repeats", "2", "--out",
            str(out)]
    run = f"import runpy; sys.argv = {argv!r}; runpy.run_path({str(SCRIPT)!r}, run_name='__main__')"
    p = _replay_child(run, tmp_path)
    assert p.returncode == 0, p.stderr
    r = json.loads(out.read_text())
    jobs = [(t, False, ("nodule",), {"nodule": c}, c.union.shape[1], c.union.shape[0])
            for t, c in B.synthetic_replay(8)]
    want = [o["nodule"] for _, o in score_many(jobs)]
    res = r["replay"]["results"]
    assert res["iou_tc"] == [o.iou_tc for o in want] and res["f5"] == [o.f5 for o in want]
    assert res["area"] == [o.area for o in want] and r["replay"]["objects"] == 8
    assert r["replay"]["stats"]["waits"] == 4 and r["inputs"]["replay_inputs.npz"]["sha256"] == _sha(npz)

    p = _replay_child(f"import runpy; runpy.run_path({str(SCRIPT)!r})", tmp_path)  # import alone
    assert p.returncode == 0, p.stderr

    fake = tmp_path / "fake"  # the pre-E old tree: objdel with only _gpucheck, no objdel.geom or objdel.verify.reward
    (fake / "objdel").mkdir(parents=True)
    (fake / "objdel/__init__.py").write_text("")
    (fake / "objdel/_gpucheck.py").write_text(textwrap.dedent(f"""
        import importlib.util
        _s = importlib.util.spec_from_file_location("objdel_real._gpucheck", {_gpucheck.__file__!r})
        _m = importlib.util.module_from_spec(_s)
        _s.loader.exec_module(_m)
        globals().update({{k: v for k, v in vars(_m).items() if not k.startswith("__")}})
    """))
    site = Path(np.__file__).parents[1]
    argv[-1] = str(tmp_path / "replay_fake.json")
    run = f"import runpy; sys.argv = {argv!r}; runpy.run_path({str(SCRIPT)!r}, run_name='__main__')"
    p = _replay_child("import objdel; assert objdel.__file__.startswith(sys.path[0])\n" + run, tmp_path,
                      (fake, site), ("-I", "-S"))
    assert p.returncode == 0, p.stderr
    assert json.loads((tmp_path / "replay_fake.json").read_text())["replay"]["results"] == res


def test_replay_refuses_a_changed_file(tmp_path):
    npz = tmp_path / "replay_inputs.npz"
    _run(["--dump-replay", npz, "--synthetic", "2"], tmp_path)
    with open(npz, "ab") as f:
        f.write(b"x")
    p = _run(["--replay", npz, "--out", tmp_path / "r.json"], tmp_path, check=False)
    assert p.returncode == 1 and "SHA-256" in p.stderr


# ------------------------------------------------------------------ rungs (D1 step 3, input to the D6 grid decision)

RUNG = {"dataset", "item_id", "rater", "object", "vertices", "max_ring", "validity_ms"}


def _check_rungs(r: dict) -> None:
    assert {"build", "root", "inputs", "rungs", "share_ge512", "time_share_ge512", "largest"} <= set(r)
    assert r["rungs"] and all(set(g) == RUNG and g["vertices"] >= g["max_ring"] >= 0 for g in r["rungs"])
    big = [g for g in r["rungs"] if g["vertices"] >= 512]
    assert r["share_ge512"] == len(big) / len(r["rungs"])
    assert r["time_share_ge512"] == pytest.approx(sum(g["validity_ms"] for g in big)
                                                  / sum(g["validity_ms"] for g in r["rungs"]))


def test_rungs_on_a_fake_root(tmp_path):
    T = tmp_path / "root"
    sq = '[[100,100],[900,100],[900,900],[100,900]]'
    rows = {"lidc": [{"item_id": "lidc/a", "rater": "1", "W": 40, "H": 40,
                      "text": '{"nodule":[{"outer":' + sq + ',"holes":[]}]}'}],
            "riga": [{"item_id": "riga/b", "rater": "2", "W": 60, "H": 60,
                      "text": '{"disc":[{"outer":' + sq + ',"holes":[]}],"cup":[]}'}]}
    for ds, rs in rows.items():
        (T / "targets" / ds).mkdir(parents=True)
        (T / "targets" / ds / "train.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rs))
    _run(["--rungs", tmp_path / "rungs.json"], tmp_path, {"OBJDEL_DATA_ROOT": str(T)})
    r = json.loads((tmp_path / "rungs.json").read_text())
    _check_rungs(r)
    assert [(g["dataset"], g["object"], g["vertices"]) for g in r["rungs"]] == [("lidc", "nodule", 4),
                                                                                ("riga", "disc", 4), ("riga", "cup", 0)]


@pytest.mark.data
def test_recorded_rungs_schema():
    """R09: runs/d1/rungs.json (bench_checker --rungs on the Mac data) holds every rung and the >= 512 shares."""
    path = REPO / "runs/d1/rungs.json"
    if not path.is_file():
        pytest.skip(f"{path} not recorded")
    r = json.loads(path.read_text())
    _check_rungs(r)
    assert {g["dataset"] for g in r["rungs"]} == {"lidc", "riga"}
