"""Break-it round 2, wiring: edge cases across the rows-root contract, preflight stamp, vast policy and final_eval.

Each test states the behaviour the code should have; strict xfail means the bug is still there (BR2-wire-<n>)."""

import hashlib
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from experiments.phase2 import final_eval as FE
from experiments.phase2 import preflight as PF
from experiments.phase2 import train_rl as TR
from experiments.phase2.vast import policy as P
from tests.test_final_eval import ev, ev_root, fake_sampler, jsonl  # noqa: F401  (fixtures and helpers)

OFFER = {"dph_total": 0.90, "inet_up_cost": 0.01, "inet_down_cost": 0.01}


# ---------------------------------------------------------------- preflight stamp


@pytest.mark.xfail(
    strict=True,
    reason="BR2-wire-1: refusal() guards the stamp keys (BR1-wire-9), not their types: "
    "a stamp that is JSON null / a list, or has hashes: null or a check that is not "
    "a dict, raises AttributeError / TypeError out of refusal(), so cmd_launch "
    "crashes with a traceback instead of the 'preflight gate' die",
)
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
    c = tmp_path / "cfg.yaml"
    c.write_text("a: 1\n")
    monkeypatch.setattr(PF, "stamp_path", lambda cfg: tmp_path / "s.json")
    monkeypatch.setattr(PF, "hashes", lambda cfg: {"code": "1", "config": "2", "data": "3"})
    (tmp_path / "s.json").write_text(stamp)
    assert isinstance(PF.refusal(c), str)


# ---------------------------------------------------------------- train_rl


@pytest.fixture
def stub(tmp_path, monkeypatch):
    """train_rl.main in process: no preload, no model; make_trainer records its arguments."""
    import gc

    import experiments.phase2.g_probe as g_probe
    from objdel.train import dataset, serve

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
    monkeypatch.setattr(g_probe, "make_trainer", lambda *a, **k: seen.trainers.append(k) or (trainer, saver, 5e-5))
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


@pytest.mark.xfail(
    strict=True,
    reason="BR2-wire-2: main() tests `if schedule:`, so a manifest whose `schedule` is {} "
    "(pinned by SHA-256, so rows_schedule returns it) is treated as a root without "
    "one: DISAGREEMENT_OVERSAMPLE rows, shuffled, no check_rows / check_sampler",
)
def test_an_empty_schedule_is_not_an_unscheduled_root(stub):
    raw = json.dumps({"schedule": {}}).encode()
    (stub.root / "manifest.json").write_bytes(raw)
    try:
        TR.main([*stub.argv, "--rows-manifest-sha256", hashlib.sha256(raw).hexdigest()])
    except (ValueError, KeyError, SystemExit):
        return  # refusing the empty schedule is right
    assert stub.trainers == [] or stub.trainers[0]["check"] is not None, "trained shuffled with no schedule guards"


@pytest.mark.xfail(
    strict=True,
    reason="BR2-wire-3: --resume does not compare the last checkpoint's step with --steps: "
    "--resume --steps 60 over checkpoint-100 trains nothing and still merges and "
    "saves the 100-step adapter as 'final' of a 60-step run",
)
def test_resume_past_the_requested_steps_is_refused(stub):
    (stub.root / "manifest.json").write_text("{}")
    (stub.out / "trainer" / "checkpoint-100").mkdir(parents=True)
    argv = [a if a != "20" else "60" for a in stub.argv]
    try:
        TR.main([*argv, "--resume"])
    except (ValueError, SystemExit):
        return
    pytest.fail(f"trained or merged from checkpoint-100 for --steps 60: {stub.trained}")


@pytest.mark.xfail(
    strict=True,
    reason="BR2-wire-4: parse accepts --save-steps 0 and negative, and --heldout-every 0: "
    "transformers takes save_steps=0 / -5 without error and then never saves, so "
    "the run has no checkpoint for --resume / --merge-last after a lost box",
)
@pytest.mark.parametrize("extra", [["--save-steps", "0"], ["--save-steps", "-5"]])
def test_parse_refuses_a_non_positive_save_steps(extra):
    with pytest.raises(SystemExit):
        TR.parse(["--model", "m", "--ckpt", "c", "--out", "o", "--train", "t", "--steps", "20", *extra])


# ---------------------------------------------------------------- vast policy


@pytest.mark.xfail(
    strict=True,
    reason="BR2-wire-5: launch_refusal checks only 'committed + worst <= HARD_STOP', so a "
    "negative price, a negative bandwidth cost or a negative --cap-hours makes the "
    "worst case negative or ~0 and passes; the box then runs with a negative "
    "WATCHDOG_CAP_S or books a negative cost in the ledger",
)
@pytest.mark.parametrize(
    "offer,cap", [(dict(OFFER, dph_total=-1000.0), 2.4), (dict(OFFER, inet_up_cost=-5.0), 2.4), (OFFER, -50.0)]
)
def test_launch_refusal_refuses_negative_money_inputs(tmp_path, monkeypatch, offer, cap):
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    assert P.launch_refusal(offer, cap)


def _repo(tmp_path, monkeypatch, files):
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    for f in ["tests/oracles/o.py", *files]:
        (tmp_path / f).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / f).write_text("x\n")
    monkeypatch.setattr(P, "ROOT", tmp_path)
    return set(P.code_files())


@pytest.mark.xfail(
    strict=True,
    reason="BR2-wire-6: code_files lets `keep` win over `drop` (`keep.match(f) or not "
    "drop.search(f)`), so a riga*.csv/.json table under a kept kernel folder "
    "docs/implementation2.0/kernels/K1/ (or a K1/data or K1/sub/build file) is "
    "rsynced to the rented box although drop's riga rule names it",
)
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


@pytest.mark.xfail(
    strict=True,
    reason="BR2-wire-7: --k is not validated, and noise_floor with K=1 splits the answers "
    "into an empty half: d is NaN, item_sd / mean_sd are NaN, and final_eval.json "
    "is written with bare NaN (invalid JSON) after the paid sampling",
)
def test_noise_floor_of_one_answer_is_refused_or_finite():
    res = _results(1)
    groups = {r["item_id"]: r["item_id"].split("/")[1] for r in res}
    try:
        nf = FE.noise_floor(res, groups, draws=50)
    except ValueError:
        return
    json.dumps(nf, allow_nan=False)


@pytest.mark.xfail(
    strict=True,
    reason="BR2-wire-8: with several --view and one --out, outs = out / view.name, and a "
    "--view path ending in '..' has the name '..': its outputs go to out's PARENT "
    "folder, outside the folder the caller named (and a '.' view would write into "
    "out itself)",
)
def test_view_path_ending_in_dotdot_stays_inside_out(ev, tmp_path, monkeypatch):  # noqa: F811
    (ev / "sub").mkdir()
    other = Path(f"{ev}/sub/..")
    monkeypatch.setattr(FE, "sample_items", fake_sampler(ev))
    monkeypatch.setattr(FE, "load_trainer", lambda a: object())
    out = tmp_path / "O" / "deep"
    try:
        FE.main(
            [
                "--ckpt",
                "c",
                "--train",
                str(ev),
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
