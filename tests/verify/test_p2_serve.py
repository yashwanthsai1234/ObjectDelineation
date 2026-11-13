"""objdel.verify.store (M3, plan D14): preloaded caches, CacheRef, TRLReward on
preloaded caches, and the K3/K4 harness. CI-safe: a mini prepared layout is written from the conftest synthetic items
(no datasets, no torch); the real-data K3 run is evidence in docs/archive/experiments/phase2/impl/M3.md."""

import json
import pickle
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from harness.k3_harness import make_jobs

from objdel.core import datasets
from objdel.rl.reward import TRLReward
from objdel.sft import select as sft_select
from objdel.verify import store as serve
from objdel.verify.reward.advantages import advantages_v4
from objdel.verify.reward.cache import ExpertCache, rl_eligible
from objdel.verify.reward.score import score_batch
from tests.support.prepared import write_prepared_root

ROOT = Path(__file__).resolve().parents[2]
EOS = 0


def _write_prepared(out: Path, items) -> Path:
    """The data/prepared layout serve reads: caches/<key>/<name>.npz, <ds>/<split>.jsonl, manifest.json files."""
    sets = {}
    for it in items:
        if it.image is None:
            continue
        rl = rl_eligible(next(iter(it.masks.values())))
        h, w = next(iter(it.masks.values())).shape[1:]
        sets.setdefault(f"{it.dataset}/{it.split}", []).append(
            {"item_id": it.item_id, "rl": rl, "caches": it.masks if rl else {}, "W": int(w), "H": int(h)}
        )
    write_prepared_root(out, sets, manifest_extra={"sets": {k: {} for k in sets}})
    return out


@pytest.fixture(scope="module")
def prepared(tmp_path_factory, synthetic_items):
    return _write_prepared(tmp_path_factory.mktemp("prep"), synthetic_items(("train", "val"), unanimous=True))


pytestmark = pytest.mark.usefixtures("clean_store")  # conftest: an empty serve.STORE per test


def _prompts(root, ds):
    keys = sorted(k for k in serve.STORE[str(Path(root).resolve())] if k.startswith(ds))
    return [(datasets.get(ds).names, serve.loader(root)(k)) for k in keys]


def _workloads(root) -> dict[str, list]:
    riga, lidc = _prompts(root, "riga"), _prompts(root, "lidc")
    return {"riga": make_jobs(riga * 2, 1), "mixed": make_jobs(riga + lidc, 2), "lidc": make_jobs(lidc, 3)}


def _trl_call(reward, jobs) -> dict:
    """TRLReward's 8 functions in TRL's order with TRL's kwargs (objdel.rl.reward)."""
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


def test_trl_reward_on_preloaded_caches(prepared):
    serve.preload(prepared, ("riga", "lidc"))
    for jobs in _workloads(prepared).values():
        out = _trl_call(TRLReward(serve.loader(prepared), {EOS}), jobs)
        assert out["advantage_v4"] == advantages_v4(score_batch(jobs)).tolist()


def test_k3_harness_help():
    r = subprocess.run(
        [sys.executable, "-I", "-m", "tests.harness.k3_harness", "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert r.returncode == 0 and "--batches" in r.stdout


def test_k3_harness_smoke(prepared, tmp_path):
    from harness import k3_harness

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
        ]
    )
    res = json.loads(out.read_text())
    assert res["gate"]["all_hashes_equal"] and res["lidc_items_preloaded"] == 4
    for w in ("riga", "mixed", "lidc"):
        assert {"serial", "per_rollout", "hashes", "hash_equal"} <= set(res[w])
    assert {"K3_pass", "K4_pass", "K4_p99_ms", "all_hashes_equal"} <= set(res["gate"])


def test_sampler_timeout_is_three_times_the_measured_time():
    """EV-5: 1 s per answer let a hung vLLM engine hold the box ~2.5 h. Three times the measured time plus the 600 s
    engine start: HF ~1 s per answer (it resumes); vLLM 0.09 s at ~580 tokens, scaled to the 1,536-token cap, since a
    vLLM retry starts over."""
    from argparse import Namespace


    assert sft_select.sampler_timeout(Namespace(sampler="vllm", timeout_s=None), 1000) == pytest.approx(600 + 750)
    assert sft_select.sampler_timeout(Namespace(sampler="hf", timeout_s=None), 512) == pytest.approx(600 + 1536)
    assert sft_select.sampler_timeout(Namespace(sampler="vllm", timeout_s=2.0), 10) == pytest.approx(620)
