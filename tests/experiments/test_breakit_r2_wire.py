"""Break-it round 2, wiring: edge cases across the rows-root contract, preflight stamp, vast policy and final_eval.

Each test states the behaviour the code should have; each guards a fixed BR2-wire finding."""

import hashlib
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from objdel.cli.eval import run as CLI
from objdel.cli.train import rl as TR
from objdel.infra import preflight as PF
from objdel.infra.vast import __main__ as P
from objdel.infra.vast import ledger as LG
from objdel.infra.vast import remote as RM
from objdel.policy import final_eval as FE
from tests.experiments.test_final_eval import ev, ev_root, fake_sampler, jsonl  # noqa: F401  (fixtures and helpers)

OFFER = {"dph_total": 0.90, "inet_up_cost": 0.01, "inet_down_cost": 0.01}


# ---------------------------------------------------------------- preflight stamp


@pytest.mark.parametrize(
    "stamp",
    [
        "null",
        "[]",
        '{"pass": true, "hashes": null}',
        '{"pass": false, "checks": {"DATA x": 1}}',
        '{"pass": false, "checks": null}',
    ],
)
def test_refusal_of_a_wrongly_typed_stamp_is_a_string(tmp_path, monkeypatch, stamp):
    """BR2-wire-1: a stamp of the wrong shape, for this spec's path, refuses with a reason (no traceback)."""
    c = tmp_path / "cfg.yaml"
    c.write_text("a: 1\n")
    monkeypatch.setattr(PF, "stamp_path", lambda cfg: tmp_path / "s.json")
    monkeypatch.setattr(PF, "hashes", lambda cfg: {"code": "1", "config": "2", "data": "3"})
    (tmp_path / "s.json").write_text(stamp.replace('{"pass"', f'{{"config": {json.dumps(str(c.resolve()))}, "pass"'))
    assert isinstance(PF.refusal(c), str)


# ---------------------------------------------------------------- train_rl


@pytest.fixture
def stub(tmp_path, monkeypatch):
    """train_rl.main in process: no preload, no model; build_trainer records its arguments."""
    import gc

    from objdel.rl import build
    from objdel.rl import rows as dataset
    from objdel.verify import store as serve

    seen = SimpleNamespace(rows_calls=[], trainers=[], trained=[])
    monkeypatch.setattr(serve, "preload", lambda root, ds: None)
    monkeypatch.setattr(gc, "freeze", lambda: None)
    monkeypatch.setattr(dataset, "grpo_rows", lambda *a, **k: seen.rows_calls.append(k) or {"key": ["k1", "k2"]})
    saver = SimpleNamespace(save_pretrained=lambda p: None)
    trainer = SimpleNamespace(
        train=lambda resume_from_checkpoint=None: seen.trained.append(resume_from_checkpoint),
        model=saver,
        accelerator=SimpleNamespace(unwrap_model=lambda m: m),
        state=SimpleNamespace(global_step=2, log_history=[]),
        args=SimpleNamespace(beta=0.0, mask_truncated_completions=True),
    )
    monkeypatch.setattr(build, "build_trainer", lambda *a, **k: seen.trainers.append(k) or (trainer, saver, 5e-5))
    root = tmp_path / "root"
    root.mkdir()
    seen.root, seen.out = root, tmp_path / "out"
    seen.argv = [
        "--model",
        "qwen3-vl-2b",
        "--ckpt",
        str(tmp_path / "c"),
        "--out",
        str(seen.out),
        "--train",
        str(root),
        "--steps",
        "20",
    ]
    return seen


def test_an_empty_schedule_is_not_an_unscheduled_root(stub):
    raw = json.dumps({"schedule": {}}).encode()
    (stub.root / "manifest.json").write_bytes(raw)
    try:
        TR.main([*stub.argv, "--rows-manifest-sha256", hashlib.sha256(raw).hexdigest(), "--rows-revision", "a" * 40])
    except (ValueError, KeyError, SystemExit):
        return  # refusing the empty schedule is right
    assert stub.trainers == [] or stub.trainers[0]["check"] is not None, "trained shuffled with no schedule guards"


def test_resume_past_the_requested_steps_is_refused(stub):
    (stub.root / "manifest.json").write_text('{"sets": {"lidc/train": {}}}')  # the train set check_root looks for
    (stub.out / "trainer" / "checkpoint-100").mkdir(parents=True)
    argv = [a if a != "20" else "60" for a in stub.argv]
    try:
        TR.main([*argv, "--resume"])
    except (ValueError, SystemExit):
        return
    pytest.fail(f"trained or merged from checkpoint-100 for --steps 60: {stub.trained}")


@pytest.mark.parametrize("extra", [["--save-steps", "0"], ["--save-steps", "-5"]])
def test_parse_refuses_a_non_positive_save_steps(extra):
    with pytest.raises(SystemExit):
        TR.parse(["--model", "m", "--ckpt", "c", "--out", "o", "--train", "t", "--steps", "20", *extra])


# ---------------------------------------------------------------- vast policy


@pytest.mark.parametrize(
    "offer,cap", [(dict(OFFER, dph_total=-1000.0), 2.4), (dict(OFFER, inet_up_cost=-5.0), 2.4), (OFFER, -50.0)]
)
def test_launch_refusal_refuses_negative_money_inputs(tmp_path, monkeypatch, offer, cap):
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    (tmp_path / "ledger.csv").write_text(",".join(LG.FIELDS) + "\n")
    assert P.launch_refusal(offer, cap).startswith("refused")


def _repo(tmp_path, monkeypatch, files):
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    # pyproject.toml: repo_root's OBJDEL_REPO check
    for f in ["tests/verify/test_gpu_scoring.py", "pyproject.toml", *files]:
        (tmp_path / f).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / f).write_text("x\n")
    monkeypatch.setenv("OBJDEL_REPO", str(tmp_path))
    return set(RM.code_files())


def test_code_files_keep_riga_tables_off_the_box_even_under_kept_folders(tmp_path, monkeypatch):
    shipped = _repo(
        tmp_path,
        monkeypatch,
        ["docs/implementation2.0/kernels/K1/riga_labels.csv", "docs/implementation2.0/kernels/K1/riga_scores.json"],
    )
    assert not {f for f in shipped if "riga" in f}, shipped


# ---------------------------------------------------------------- final_eval


def _results(k, n=6):
    return [{"item_id": f"lidc/p{i % 3}/n000/z{i}", "samples": {"iou_tc": [0.1 * (i + 1)] * k}} for i in range(n)]


def test_noise_floor_of_one_answer_is_refused_or_finite():
    res = _results(1)
    groups = {r["item_id"]: r["item_id"].split("/")[1] for r in res}
    try:
        nf = FE.noise_floor(res, groups, draws=50)
    except ValueError:
        return
    json.dumps(nf, allow_nan=False)


def test_view_path_ending_in_dotdot_stays_inside_out(ev, tmp_path, monkeypatch):  # noqa: F811
    (ev / "sub").mkdir()
    other = Path(f"{ev}/sub/..")
    monkeypatch.setattr(FE, "load_sampler", fake_sampler(ev))
    out = tmp_path / "O" / "deep"
    try:
        CLI.main(
            [
                "--ckpt",
                "c",
                "--view",
                str(ev),
                "--split",
                "val",
                "--view",
                str(other),
                "--split",
                "val",
                "--out",
                str(out),
            ]
        )
    except SystemExit:
        return
    wrote = {p.parent.resolve() for p in tmp_path.rglob("final_eval.json")}
    assert all(out.resolve() in w.parents for w in wrote), wrote
