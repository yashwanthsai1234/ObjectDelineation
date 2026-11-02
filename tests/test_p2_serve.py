"""objdel.train.serve (M3, plan D14): preloaded caches, CacheRef, the warmed fork pool, TRLReward on the pool, and
the K3/K4 harness. CI-safe: a mini prepared layout is written from the conftest synthetic items (no datasets, no
torch); the real-data K3 run is evidence in experiments/phase2/impl/M3.md."""

import gc
import hashlib
import json
import os
import pickle
import shutil
import signal
import subprocess
import sys
import time
import types
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

import pytest

from experiments.phase1.reward_gates import digest, make_jobs
from objdel.geometry import SCHEMAS
from objdel.reward import ExpertCache, rl_eligible, score_batch
from objdel.train import serve
from objdel.train.rewards import TRLReward

ROOT = Path(__file__).resolve().parents[1]
EOS = 0


def _write_prepared(out: Path, items) -> Path:
    """The data/prepared layout serve reads: caches/<key>/<name>.npz, <ds>/<split>.jsonl, manifest.json files."""
    sets = {}
    for it in items:
        if it.image is None:
            continue
        rl = rl_eligible(next(iter(it.masks.values())))
        caches = {}
        for name, stack in it.masks.items():
            if rl:
                (out / "caches" / it.item_id).mkdir(parents=True, exist_ok=True)
                caches[name] = ExpertCache.build(stack).save(out / "caches" / it.item_id / f"{name}.npz")
        h, w = next(iter(it.masks.values())).shape[1:]
        sets.setdefault(f"{it.dataset}/{it.split}", []).append(
            {"item_id": it.item_id, "rl": rl, "caches": caches, "W": int(w), "H": int(h)}
        )
    files = {}
    for rel, rows in sets.items():
        p = out / f"{rel}.jsonl"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("".join(json.dumps(r) + "\n" for r in rows))
        files[f"{rel}.jsonl"] = {"sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
    (out / "manifest.json").write_text(json.dumps({"files": files, "sets": {k: {} for k in sets}}))
    return out


@pytest.fixture(scope="module")
def prepared(tmp_path_factory, synthetic_items):
    return _write_prepared(tmp_path_factory.mktemp("prep"), synthetic_items(("train", "val"), unanimous=True))


pytestmark = pytest.mark.usefixtures("clean_store")  # conftest: an empty serve.STORE per test


@pytest.fixture
def pool():
    ex = serve.fork_pool(2)
    yield ex
    ex.shutdown(cancel_futures=True)


def _prompts(root, ds):
    keys = sorted(k for k in serve.STORE[str(Path(root).resolve())] if k.startswith(ds))
    return [(SCHEMAS[ds], serve.loader(root)(k)) for k in keys]


def _workloads(root) -> dict[str, list]:
    riga, lidc = _prompts(root, "riga"), _prompts(root, "lidc")
    return {"riga": make_jobs(riga * 2, 1), "mixed": make_jobs(riga + lidc, 2), "lidc": make_jobs(lidc, 3)}


def _trl_call(reward, jobs) -> dict:
    """TRLReward's 8 functions in TRL's order with TRL's kwargs (objdel.train.rewards)."""
    funcs, _ = reward.funcs()
    kw = {
        "prompts": [[{"role": "user", "content": "p"}]] * len(jobs),
        "completions": [[{"role": "assistant", "content": j[0]}] for j in jobs],
        "completion_ids": [[7, EOS] for _ in jobs],
        "key": [j[3].key for j in jobs],
        "W": [j[4] for j in jobs],
        "H": [j[5] for j in jobs],
        "names": [list(j[2]) for j in jobs],
        "trainer_state": None,
        "log_metric": lambda *_: None,
        "log_extra": lambda *_: None,
    }
    return {f.__name__: f(**kw) for f in funcs}


def test_preload_lidc_only_loads_no_riga_cache(prepared, monkeypatch):
    loaded = []
    real = ExpertCache.load
    monkeypatch.setattr(ExpertCache, "load", lambda path, *a, **k: loaded.append(Path(path)) or real(path, *a, **k))
    assert serve.preload(prepared, ("lidc",), threads=2) == 4  # z001 and z004 of train and val
    assert loaded and all(p.relative_to(prepared / "caches").parts[0] == "lidc" for p in loaded)
    assert not any(k.startswith("riga") for s in serve.STORE.values() for k in s)
    ref = serve.loader(prepared)("lidc/1.2.1/n000/z001")
    assert set(ref) == {"nodule"} and ref["nodule"].k == 2 and len(ref) == 1
    with pytest.raises(KeyError):
        serve.loader(prepared)("lidc/1.2.1/n000/z003")  # one rater: not RL, never loaded


def test_preload_checks_digests(prepared, tmp_path):
    root = Path(shutil.copytree(prepared, tmp_path / "p"))
    a, b = (root / "caches" / f"lidc/1.2.1/n000/z00{z}" / "nodule.npz" for z in (1, 4))
    a.write_bytes(b.read_bytes())  # a valid cache under the wrong item: the record's digest differs
    with pytest.raises(ValueError, match="digests differ"):
        serve.preload(root, ("lidc",))
    root = Path(shutil.copytree(prepared, tmp_path / "q"))
    rec = root / "lidc" / "train.jsonl"
    rec.write_text(rec.read_text().replace('"W": 100', '"W": 101', 1))  # records no longer the manifest's
    with pytest.raises(ValueError, match="manifest"):
        serve.preload(root, ("lidc",))
    with pytest.raises(ValueError, match="no foo records"):
        serve.preload(prepared, ("foo",))


def test_cache_ref_pickles_as_its_key_and_roots_stay_apart(prepared, tmp_path):
    other = Path(shutil.copytree(prepared, tmp_path / "view"))  # e.g. an eval view: same item ids, other caches
    serve.preload(prepared, ("riga", "lidc"))
    key = "riga/MESSIDOR/image001"
    ref = serve.loader(prepared)(key)
    assert dict(ref) == serve.STORE[str(prepared.resolve())][key]
    assert len(pickle.dumps(ref)) < 300 and pickle.loads(pickle.dumps(ref))["disc"] is ref["disc"]
    serve.preload(other, ("riga",))
    assert serve.loader(other)(key)["cup"] is serve.STORE[str(other.resolve())][key]["cup"]
    assert serve.loader(prepared)(key)["cup"] is not serve.loader(other)(key)["cup"]


def test_pool_hashes_equal_serial_on_all_workloads(prepared):
    serve.preload(prepared, ("riga", "lidc"))
    batches = _workloads(prepared)
    ex = serve.fork_pool(2)  # forked after the preload: the workers resolve CacheRefs in their copy of STORE
    try:
        for name, jobs in batches.items():
            assert digest(score_batch(jobs, ex)) == digest(score_batch(jobs)), name
    finally:
        ex.shutdown()


def test_trl_reward_one_dispatch_per_generation_batch(prepared):
    serve.preload(prepared, ("riga", "lidc"))
    batches = list(_workloads(prepared).values())

    class Counting:
        def __init__(self, ex):
            self.ex, self.n = ex, 0

        def map(self, *a, **k):
            self.n += 1
            return self.ex.map(*a, **k)

    ex = serve.fork_pool(2)
    try:
        cp = Counting(ex)
        par, ser = TRLReward(serve.loader(prepared), {EOS}, pool=cp), TRLReward(serve.loader(prepared), {EOS})
        pids = set(ex._processes)
        for jobs in batches:
            n0 = cp.n
            assert _trl_call(par, jobs) == _trl_call(ser, jobs)
            assert cp.n - n0 == 1
        assert set(ex._processes) == pids  # no respawn
    finally:
        ex.shutdown()


def test_fork_pool_forks_every_worker_now(pool):
    assert len(pool._processes) == 2 and gc.get_freeze_count() > 0
    assert set(pool.map(int, ["1", "2"])) == {1, 2}


def test_killed_worker_raises_broken_process_pool_within_1s(prepared):
    serve.preload(prepared, ("lidc",))  # before the fork (BK2-11)
    pool = serve.fork_pool(2)
    try:
        jobs = make_jobs(_prompts(prepared, "lidc"), 4)
        fut = pool.submit(time.sleep, 30)
        os.kill(next(iter(pool._processes)), signal.SIGKILL)
        t0 = time.perf_counter()
        with pytest.raises(BrokenProcessPool):
            fut.result(timeout=5)
        assert time.perf_counter() - t0 < 1.0
        with pytest.raises(BrokenProcessPool):
            score_batch(jobs, pool)
    finally:
        pool.shutdown(cancel_futures=True)


def test_preload_is_refused_while_a_forked_pool_is_alive(prepared):
    """BK2-11: workers keep STORE as it was at the fork, so a later preload (a new root, or rebuilt caches of the same
    root) would be scored on the old caches; it is refused until the pool is shut down."""
    serve.preload(prepared, ("lidc",))
    ex = serve.fork_pool(1)
    try:
        with pytest.raises(RuntimeError, match="fork_pool"):
            serve.preload(prepared, ("lidc",))
    finally:
        ex.shutdown()
    assert serve.preload(prepared, ("lidc",)) > 0


def test_fork_pool_imports_the_scorer_modules_before_forking(monkeypatch):
    """S1-11: on macOS a worker importing a module after the fork can crash the pool; fork_pool imports the scorer
    modules (and any extra `modules`) in the parent first, so the workers inherit them."""
    monkeypatch.delitem(sys.modules, "json.tool", raising=False)  # (scorers stay: re-importing them splits classes)
    ex = serve.fork_pool(1, modules=("json.tool",))
    try:
        assert all(name in sys.modules for name in (*serve.SCORERS, "json.tool"))
        assert ex.submit(_imported, ("json.tool", *serve.SCORERS)).result()
    finally:
        ex.shutdown()


def _imported(names):
    return all(n in sys.modules for n in names)


def test_fork_pool_refuses_after_cuda_init(monkeypatch):
    def fake(initialized):
        return types.SimpleNamespace(cuda=types.SimpleNamespace(is_initialized=lambda: initialized))

    monkeypatch.setitem(sys.modules, "torch", fake(True))
    with pytest.raises(RuntimeError, match="CUDA"):
        serve.fork_pool(1)
    monkeypatch.setitem(sys.modules, "torch", fake(False))
    ex = serve.fork_pool(1)
    ex.shutdown()


def test_k3_harness_help():
    r = subprocess.run(
        [sys.executable, "-I", "-m", "experiments.phase2.k3_harness", "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert r.returncode == 0 and "--workers" in r.stdout


def test_k3_harness_smoke(prepared, tmp_path):
    from experiments.phase2 import k3_harness

    out = tmp_path / "k3.json"
    k3_harness.main(
        [
            "--root",
            str(prepared),
            "--out",
            str(out),
            "--synth",
            "2",
            "--lidc",
            "2",
            "--batches",
            "1",
            "--threads",
            "2",
            "--workers",
            "2",
        ]
    )
    res = json.loads(out.read_text())
    assert res["gate"]["all_hashes_equal"] and res["lidc_items_preloaded"] == 4
    for w in ("riga", "mixed", "lidc"):
        assert {"serial", "threads2", "fork2", "per_rollout", "hashes"} <= set(res[w])
    assert {"K3_pass", "K4_pass", "K4_p99_ms", "d14_cpp_trigger"} <= set(res["gate"])
