"""Break-it round 1, wiring: edge cases across the rows-root contract, vast policy, final_eval and preflight.

Each test states the behaviour the code should have; strict xfail means the bug is still there (BR1-wire-<n>)."""
import json
import math
import shutil
import subprocess
from types import SimpleNamespace

import pytest

from experiments.phase2 import final_eval as CLI
from experiments.phase2 import train_rl as TR
from objdel.data import quality as Q
from objdel.infra import preflight as PF
from objdel.infra.vast import __main__ as P
from objdel.infra.vast import ledger as LG
from objdel.infra.vast import remote as RM
from objdel.policy import final_eval as FE
from objdel.policy.errors import Refusal
from objdel.rl import run as rl_run
from tests.experiments.test_final_eval import ev, ev_root, jsonl, run_fe  # noqa: F401  (fixtures and helpers)

SCHEDULE = {"steps": 100, "per_step": 32, "block_steps": 20, "ids_sha256": "0" * 64}


def cfg(**kw):
    base = dict(generation_batch_size=256, num_generations=8, steps_per_generation=4, gradient_accumulation_steps=4,
                num_iterations=1)
    return SimpleNamespace(**(base | kw))


# ---------------------------------------------------------------- train_rl

@pytest.mark.parametrize("steps", [0, -20])
def test_check_sampler_refuses_non_positive_steps(steps):
    with pytest.raises(ValueError):
        rl_run.check_sampler(cfg(), SCHEDULE, steps)


@pytest.mark.parametrize("flag,val", [("--budget-usd", "nan"), ("--rate-usd-per-h", "nan")])
def test_parse_refuses_nan_budget(flag, val):
    argv = ["--model", "m", "--ckpt", "c", "--out", "o", "--train", "t", "--steps", "20", "--budget-usd", "5",
            "--session-start", "0", "--rate-usd-per-h", "1"]
    i = argv.index(flag)
    argv[i + 1] = val
    with pytest.raises(SystemExit):
        TR.parse(argv)


# ---------------------------------------------------------------- vast policy

OFFER = {"dph_total": 0.90, "inet_up_cost": 0.01, "inet_down_cost": 0.01}


def test_launch_refusal_refuses_a_nan_price(tmp_path, monkeypatch):
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    (tmp_path / "ledger.csv").write_text(",".join(LG.FIELDS) + "\n")
    assert P.launch_refusal(dict(OFFER, dph_total=math.nan), 2.4).startswith("refused")


def _repo(tmp_path, monkeypatch, files):
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    for f in ["tests/verify/test_gpu_scoring.py", *files]:
        (tmp_path / f).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / f).write_text("x\n")
    monkeypatch.setenv("OBJDEL_REPO", str(tmp_path))
    return set(RM.code_files())


@pytest.mark.xfail(strict=True, reason="BR1-wire-4: code_files drops RIGA tables only by a name pattern "
                                       "riga[^/]*.(csv|json|jsonl): a RIGA table in a riga/ folder, or as "
                                       ".parquet/.tsv/.npz, is rsynced to the rented box")
def test_code_files_keep_riga_derived_tables_off_the_box(tmp_path, monkeypatch):
    shipped = _repo(tmp_path, monkeypatch, ["experiments/phase2/riga/table.csv",
                                            "experiments/phase2/riga_cache.parquet",
                                            "experiments/phase2/riga_labels.tsv"])
    assert not {f for f in shipped if "riga" in f}, shipped


# ---------------------------------------------------------------- final_eval

def test_check_view_frozen_check_survives_an_alias(ev, tmp_path, monkeypatch):  # noqa: F811
    frozen = ev.parent / "lidc-val-frozenfix"
    ev.rename(frozen)
    monkeypatch.setattr(Q, "FROZEN", {frozen.name: {}})
    monkeypatch.setattr(Q, "view_content_sha256", lambda v: "got")
    monkeypatch.setattr(Q, "target_set_sha256", lambda v: "got")
    monkeypatch.setattr(FE, "FROZEN_VIEWS", tmp_path / "frozen.json")
    bad = {"content_sha256": "bad", "target_set_sha256": "bad"}
    (tmp_path / "frozen.json").write_text(json.dumps({frozen.name: bad}))
    with pytest.raises(Refusal):
        FE.check_view(frozen, "lidc", "val")  # control: the real name is refused
    alias = tmp_path / "alias"
    alias.symlink_to(frozen)
    with pytest.raises(Refusal):
        FE.check_view(alias, "lidc", "val")


def test_multi_view_one_out_does_not_collide(ev, tmp_path, monkeypatch):  # noqa: F811
    other = tmp_path / "elsewhere" / ev.name
    shutil.copytree(ev, other)
    monkeypatch.setattr(FE, "load_sampler",
                        __import__("tests.experiments.test_final_eval", fromlist=["x"]).fake_sampler(ev))
    out = tmp_path / "O"
    try:
        CLI.main(["--ckpt", "c", "--view", str(ev), "--split", "val", "--view", str(other),
                 "--split", "val", "--out", str(out)])
    except SystemExit:
        return  # refusing is fine
    assert len([d for d in out.iterdir() if d.is_dir()]) == 2, "two views, one output folder"


def test_rescore_into_its_own_source_is_refused(ev, tmp_path, monkeypatch):  # noqa: F811
    src = tmp_path / "src"
    run_fe(ev, src, monkeypatch)
    before = (src / "final_eval.json").read_bytes()
    try:
        CLI.main(["--rescore", str(src), "--view", str(ev), "--split", "val", "--out", str(src)])
    except SystemExit:
        pass
    assert (src / "final_eval.json").read_bytes() == before


def test_summarize_one_item_is_finite():
    it = {"item_id": "lidc/1.2.1/n000/z001", "iou_tc": 0.5, "sel_mean_of_k": 0.5}
    s = FE.summarize([it], {it["item_id"]: "p1"}, draws=50)
    json.dumps(s, allow_nan=False)


def test_summarize_zero_items_is_refused_or_finite():
    try:
        s = FE.summarize([], {}, draws=50)
    except (ValueError, IndexError):
        return
    json.dumps(s, allow_nan=False)


def test_summarize_refuses_duplicate_items():
    a = {"item_id": "lidc/1.2.1/n000/z001", "iou_tc": 0.5, "sel_mean_of_k": 0.5}
    b = {"item_id": "lidc/1.2.2/n000/z001", "iou_tc": 0.1, "sel_mean_of_k": 0.1}
    with pytest.raises(ValueError):
        FE.summarize([a, a, b], {a["item_id"]: "p1", b["item_id"]: "p2"}, draws=50)


# ---------------------------------------------------------------- preflight

def test_refusal_of_a_malformed_stamp_is_a_string(tmp_path, monkeypatch):
    c = tmp_path / "cfg.yaml"
    c.write_text("a: 1\n")
    monkeypatch.setattr(PF, "stamp_path", lambda cfg: tmp_path / "s.json")
    monkeypatch.setattr(PF, "hashes", lambda cfg: {"code": "1", "config": "2", "data": "3"})
    (tmp_path / "s.json").write_text(json.dumps({"pass": True, "hashes": {"code": "1"}}))
    assert isinstance(PF.refusal(c), str)
