"""The implementation-2.0 GPU test rental ($0 checks on the Mac): specs/impl2_gpu_test.yaml and its two-box split
(impl2_gpu_test_a.yaml, impl2_gpu_test_b.yaml) parse, every module, flag, node id and file they name exists, the hard
and soft steps are where the recovery sequence needs them, and the money block gives a cap that fits the $3 hard stop.
Like tests/test_b1_plans.py: a missing flag fails here and not on the paid box."""

import csv
import json
import shlex
import subprocess
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")
from test_b1_plans import parse_step  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "specs"
PLAN, B1 = "impl2_gpu_test.yaml", "b1_rows.yaml"
SPLIT = {"impl2_gpu_test_a.yaml": (0, 1, 2, 3), "impl2_gpu_test_b.yaml": (0, 1, *range(4, 12))}
RD, DATA = "/workspace/runs/r", "/workspace/data/lidc"
HARD = {0, 1, 6}  # build, env_check, the lost-run simulation; every other step is soft ("~")
SOFT_TAIL_USD = 0.35  # front 6.4 rule 5: BudgetStop = approved - the tail after the RL steps
HARD_STOP_ROOM_USD = 3.00  # the ledger closed at $16.90, HARD_STOP 19.90 (owner, 2026-11-09)


# The money block of the old plan file (front 6 rules and planning inputs); the spec holds only compute.cap_hours and
# the RL stages' budget_usd, and the tests below check both against it.
BUDGET = {
    "approved_usd": 2.20,
    "planning_dph_usd": 0.641,  # front 6: the offer limit, the highest rate the launch accepts
    "planning_gb_usd": 0.0065,  # front 6: bandwidth ceiling per GB, up and down
    "fixed_gb_down": 15,  # LG.DOWN_GB
    "gb_up_per_h": 6,  # LG.UP_GB_PER_H
    "cap_hours": 3.05,
    "budget_stop_usd": 1.85,
    "ledger_closed_usd": 16.90,  # measured: ledger.csv, 2026-11-08 (hard stop 19.90 = 16.90 + 3.00)
    "minutes": {"launch_to_first_run": 38, "build": 3, "env_check": 1, "gpu_check": 44, "gpu_smoke": 5,
                "eval_smoke_r5": 6.2, "rl_run1": 15.2, "merge_last": 2, "rl_resume": 10.6, "eval_smoke_final": 6.0},
}


def raw() -> dict:
    return {"budget": BUDGET}


def steps() -> list[list[str]]:
    from objdel.infra.vast import chain as CH

    return CH.chain(CONFIGS / PLAN, "r")


def words(step: list[str]) -> list[str]:
    return step[1:] if step[0] == "~" else step


def after(step: list[str], flag: str) -> str:
    return step[step.index(flag) + 1]


def module_steps(name: str) -> list[list[str]]:
    return [s for s in steps() if "-m" in words(s)[:1] and words(s)[1].startswith("experiments.phase2.")
            and words(s)[1].endswith(name)]


def floor_to(x: float, q: float = 0.05) -> float:
    return round(int(x / q + 1e-9) * q, 2)


# --- the plan ---------------------------------------------------------------------------------------------------


def test_the_plan_formats_and_has_the_twelve_steps_in_order():
    s = steps()
    assert len(s) == 12
    assert [x[0] == "~" for x in s] == [i not in HARD for i in range(12)]  # hard: build, env_check, crash simulation
    assert s[0][:2] == ["bash", "-c"] and "run_box.sh --build-only" in s[0][2]
    assert "assert g.build() == 'cuda'" in s[0][2] and s[0][2].endswith("tests/test_gpu_scoring.py")
    assert s[1][:2] == ["bash", "experiments/phase2/docker/env_check.sh"]
    assert s[2][1:3] == ["bash", "-c"] and "run_box.sh --gpu-check" in s[2][3]
    assert words(s[3])[:4] == ["-m", "pytest", "-m", "gpu"]
    assert [words(x)[1].rsplit(".", 1)[-1] for x in s[4:6]] == ["final_eval", "train_rl"]
    assert words(s[7])[1].endswith("train_rl") and "--merge-last" in s[7]
    assert words(s[9])[1].endswith("train_rl") and "--resume" in s[9]
    assert words(s[11])[1].endswith("final_eval") and after(s[11], "--ckpt") == f"{RD}/smoke-rl/final"


def test_gpu_check_step_follows_the_d8_form():
    """D8: OBJDEL_DATA_ROOT and PY on the box, one run_box.sh --gpu-check entry, its exit code kept for pull.sh."""
    cmd = steps()[2][3]
    assert cmd.startswith(f"OBJDEL_DATA_ROOT={DATA} PY=/opt/venv/bin/python bash experiments/cuda_checker/run_box.sh "
                          "--gpu-check;")
    assert f"echo $rc > {RD}/gpu_check.exit" in cmd and cmd.endswith("exit $rc")
    box = (ROOT / "experiments/cuda_checker/run_box.sh").read_text()
    assert "OBJDEL_RUN_DIR" in box  # $OUT = $OBJDEL_RUN_DIR/d, which pull.sh brings back


def test_bash_steps_are_valid_shell_and_their_python_snippets_compile():
    n = 0
    for s in steps():
        if words(s)[0] != "bash" or words(s)[1] != "-c":
            continue
        n += 1
        text = words(s)[2]
        ok = subprocess.run(["bash", "-n", "-c", text], capture_output=True, text=True)
        assert ok.returncode == 0, ok.stderr
        toks = shlex.split(text)
        for i, t in enumerate(toks):
            if t == "-c":
                compile(toks[i + 1], "<plan>", "exec")
    assert n == 5


def test_named_files_exist():
    for f in ("experiments/cuda_checker/run_box.sh", "experiments/phase2/docker/env_check.sh",
              "tests/test_gpu_scoring.py", "tests/test_train_fast.py"):
        assert (ROOT / f).is_file(), f
    box = (ROOT / "experiments/cuda_checker/run_box.sh").read_text()
    assert "--build-only" in box


def test_smoke_node_ids_are_the_b1_rows_ones():
    """test_b1_plans collects exactly these two ids; here they must stay the same two."""
    from objdel.infra.vast import chain as CH

    b1 = CH.chain(CONFIGS / B1, "r")[2]
    ids = [w for w in words(steps()[3]) if "::" in w]
    assert ids == [w for w in b1 if "::" in w] and len(ids) == 2


def test_every_output_lands_under_the_run_folder_and_only_val_is_read():
    for s in steps():
        if "--out" in s:
            assert after(s, "--out").startswith(f"{RD}/"), s
    evals = module_steps("final_eval")
    assert len(evals) == 2
    for e in evals:
        assert [e[i + 1].rsplit("/", 1)[1] for i, w in enumerate(e) if w == "--view"] == ["lidc-val-shifted"]
        assert after(e, "--split") == "val" and "--train" not in e
        assert 1 <= int(after(e, "--limit")) <= 32  # a smoke slice, not a number to quote
    assert "lidc-test" not in " ".join(" ".join(s) for s in steps())
    assert not [s for s in steps() if "--split" in s and after(s, "--split") != "val"]


# --- every module step parses with its module's own parser ------------------------------------------------------


def test_every_module_step_parses_with_its_module():
    pytest.importorskip("trl")
    found = [s for s in steps() if words(s)[0] == "-m" and words(s)[1].startswith("experiments.phase2.")]
    assert len(found) == 5  # final_eval x2, train_rl x3 (first run, --merge-last, --resume)
    for s in found:
        parse_step(s)


def test_final_eval_steps_name_the_limit_flag_that_exists():
    pytest.importorskip("trl")
    for e in module_steps("final_eval"):
        ns = parse_step(e)
        assert ns.limit == int(after(e, "--limit")) and ns.view and ns.ckpt


def test_a_missing_required_flag_fails_the_parse():
    pytest.importorskip("trl")
    run1 = module_steps("train_rl")[0]
    no_model = run1[: run1.index("--model")] + run1[run1.index("--model") + 2:]
    with pytest.raises(SystemExit) as e:
        parse_step(no_model)
    assert e.value.code == 2
    ev = module_steps("final_eval")[0]
    bad = ev[: ev.index("--limit")] + ev[ev.index("--limit"):][:2] + ["--bogus-flag"] + ev[ev.index("--limit") + 2:]
    with pytest.raises(SystemExit) as e:
        parse_step(bad)
    assert e.value.code == 2


# --- the RL smoke and its recovery ------------------------------------------------------------------------------


def test_rl_smoke_is_one_block_with_one_lost_run_recovery():
    pytest.importorskip("trl")
    run1, merge, resume = module_steps("train_rl")
    a, m, r = parse_step(run1), parse_step(merge), parse_step(resume)
    from objdel.infra.vast import chain as CH

    b1 = parse_step(CH.chain(CONFIGS / B1, "r")[3])
    assert a.steps == r.steps == 20 and a.steps % 20 == 0 and a.steps < b1.steps  # whole block, inside the schedule
    assert a.save_steps == r.save_steps == 10 and a.steps % a.save_steps == 0 and a.save_steps < a.steps
    assert a.shared_prefix and r.shared_prefix and a.lr == r.lr == b1.lr and not a.resume and r.resume
    assert m.merge_last and not a.merge_last and not r.merge_last
    assert a.out == m.out == r.out and a.ckpt == m.ckpt == r.ckpt == f"{RD}/ckpt"
    assert a.train == r.train == b1.train == f"{DATA}/rows/lidc-train-b1"
    assert (a.rows_revision, a.rows_manifest_sha256) == (r.rows_revision, r.rows_manifest_sha256) \
        == (b1.rows_revision, b1.rows_manifest_sha256)
    assert a.budget_usd == r.budget_usd == pytest.approx(raw()["budget"]["budget_stop_usd"])
    crash = steps()[6][2]  # keeps checkpoint-10, drops the 20-step save and everything run 1 wrote at its end
    assert "checkpoint-10" in crash and "rm -rf" in crash and "trainer/checkpoint-20" in crash
    assert "$d/final $d/final_adapter" in crash and f"{RD}/rl_run1.json" in crash
    check = steps()[8][3]
    assert "d['partial'] is True and d['steps'] == 10" in check and "final_partial10" in check
    after_resume = steps()[10][3]
    assert "d['steps'] == 20" in after_resume and "checkpoint-10" in after_resume
    assert "final_adapter" in after_resume


def test_budget_stop_has_the_run_py_clock_inputs():
    """--budget-usd needs OBJDEL_SESSION_START and OBJDEL_RATE_USD_PER_H, which run.sh exports (policy.cmd_run)."""
    from experiments.phase2.vast import policy

    src = Path(policy.__file__).read_text()
    for var in ("OBJDEL_SESSION_START", "OBJDEL_RATE_USD_PER_H", "OBJDEL_FIXED_USD"):
        assert var in src


# --- the rows root pins -----------------------------------------------------------------------------------------


def test_rows_pins_equal_the_b1_rows_pins_and_resolve_for_fetch_hf():
    from experiments.phase2.vast import policy
    from objdel.infra.vast import chain as CH

    def pins(name):
        s = next(s for s in CH.chain(CONFIGS / name, "r") if "experiments.phase2.train_rl" in s
                 and "--rows-revision" in s)
        return after(s, "--rows-revision"), after(s, "--rows-manifest-sha256")

    rev, sha = pins(PLAN)
    assert (rev, sha) == pins(B1), "fill the same revision and SHA-256 in both plans after the rows upload"
    assert len(rev) == 40 and len(sha) == 64
    assert policy.rows_pin(rev) == sha  # fetch_hf.sh --rows REV finds exactly one pin across the plans
    for s in steps():
        if "--rows-revision" in s:
            assert after(s, "--rows-revision") == rev and after(s, "--rows-manifest-sha256") == sha


def test_the_spec_holds_what_the_launch_config_held_and_the_money_block():
    """E6.6: impl2_gpu_launch.yaml is archived; each spec has b1's bundle_view, compute.cap_hours is the money block's
    cap, and every RL stage stops at its budget_stop_usd."""
    from objdel.core.spec import load_spec

    b1 = load_spec(CONFIGS / B1)
    for name in (PLAN, *SPLIT):
        spec = load_spec(CONFIGS / name)
        assert spec.bundle_view == b1.bundle_view and spec.compute == b1.compute.model_copy(
            update={"cap_hours": BUDGET["cap_hours"]})
        assert {s.budget_usd for s in spec.stages if s.kind == "rl" and not s.merge_last} <= {BUDGET["budget_stop_usd"]}


@pytest.mark.parametrize("name", SPLIT)
def test_the_two_box_split_runs_the_plan_steps_unchanged(name):
    """Owner, 2026-11-09: box A = steps 0-3 (build, env check, gpu-check, GPU smokes); box B = steps 0, 1, 4-11."""
    from objdel.infra.vast import chain as CH

    assert CH.chain(CONFIGS / name, "r") == [steps()[i] for i in SPLIT[name]]


# --- the money --------------------------------------------------------------------------------------------------


def test_budget_block_gives_a_cap_inside_the_three_dollars(tmp_path, monkeypatch):
    from experiments.phase2.vast import policy
    from objdel.infra.vast import ledger as LG
    from objdel.infra.vast import remote as RM

    b = raw()["budget"]
    fixed = b["fixed_gb_down"] * b["planning_gb_usd"]
    rate = b["planning_dph_usd"] + b["gb_up_per_h"] * b["planning_gb_usd"]
    assert (fixed, rate) == (pytest.approx(0.0975), pytest.approx(0.680))
    assert b["cap_hours"] == floor_to((b["approved_usd"] - fixed) / rate)  # front 6: floor_to_0.05((amount - b) / r_s)
    assert b["budget_stop_usd"] == pytest.approx(b["approved_usd"] - SOFT_TAIL_USD)
    assert (b["fixed_gb_down"], b["gb_up_per_h"]) == (LG.DOWN_GB, LG.UP_GB_PER_H)
    # the worst case the launch checks: the offer limit, the cap plus the 600 s pull window
    offer = {"dph_total": b["planning_dph_usd"], "inet_up_cost": b["planning_gb_usd"],
             "inet_down_cost": b["planning_gb_usd"]}
    run = {"dph": offer["dph_total"], "inet_up_cost": offer["inet_up_cost"], "inet_down_cost": offer["inet_down_cost"],
           "started": 0.0}
    worst = LG.accrued(run, b["cap_hours"] * 3600 + RM.PULL_S)
    assert worst <= b["approved_usd"] + 0.12 and worst + b["ledger_closed_usd"] <= LG.HARD_STOP
    assert LG.HARD_STOP - b["ledger_closed_usd"] <= HARD_STOP_ROOM_USD + 1e-9 and worst <= HARD_STOP_ROOM_USD
    # launch_refusal on a ledger that holds the closed spend: this cap passes, a cap two hours longer is refused
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    with (tmp_path / "ledger.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, LG.FIELDS)
        w.writeheader()
        w.writerow({"ts": "2026-11-08T17:11:15Z", "event": "teardown", "run_id": "r5",
                    "cost_usd": b["ledger_closed_usd"]})
    assert policy.launch_refusal(offer, b["cap_hours"]) == ""
    assert "hard stop" in policy.launch_refusal(offer, b["cap_hours"] + 2)


def test_minutes_fit_the_cap_with_room_and_the_expected_spend_is_under_the_approved_amount():
    b = raw()["budget"]
    total_h = sum(b["minutes"].values()) / 60
    rate = b["planning_dph_usd"] + b["gb_up_per_h"] * b["planning_gb_usd"]
    assert total_h <= b["cap_hours"] * 0.8  # the cap is a guard with margin, not the forecast
    assert b["fixed_gb_down"] * b["planning_gb_usd"] + total_h * rate < b["approved_usd"]
    for key in ("launch_to_first_run", "build", "env_check", "gpu_check", "gpu_smoke", "eval_smoke_r5", "rl_run1",
                "merge_last", "rl_resume", "eval_smoke_final"):
        assert key in b["minutes"]


def test_plan_pushes_nothing_and_rents_nothing():
    text = json.dumps(steps())
    for bad in ("push_hf", "hf_push", "upload", "vastai", "--split\", \"test"):
        assert bad not in text
