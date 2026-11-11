"""B1.5 run guards of experiments.phase2.train_rl and objdel.policy.callbacks.BudgetStop (U41-U46, U61, I16, I33,
I34, R04): every refusal happens before a model or a vLLM engine loads, and the budget stop measures session spend
as LG.accrued does."""

import hashlib
import json

import pytest

from objdel.rl import run as rl_run

H = 3600.0


class Ctl:
    should_training_stop = should_save = False


def budget(rate, cap, start=0.0, fixed=0.0):
    from objdel.policy.callbacks import BudgetStop

    now = [start]
    return BudgetStop(rate, cap, start, fixed_usd=fixed, clock=lambda: now[0]), now


def step(cb, now, t) -> Ctl:
    c = Ctl()
    now[0] = t
    cb.on_step_end(None, None, c)
    return c


def test_budget_stop_trips_exactly_at_the_cap():
    cb, now = budget(2.0, 1.0, start=100.0)
    c = step(cb, now, 100.0 + 0.4999 * H)
    assert not c.should_training_stop and not c.should_save
    c = step(cb, now, 100.0 + 0.5 * H)
    assert c.should_training_stop and c.should_save
    c = step(cb, now, 100.0 + 0.1 * H)  # tripped once: never reset, even when the clock reads earlier
    assert c.should_training_stop and c.should_save


def test_budget_stop_counts_from_the_session_start_not_the_process_start():
    """Launch 30 min before the process: the guard trips 30 min earlier than a process-clock guard would."""
    cb, now = budget(2.0, 1.0, start=0.0)
    assert not step(cb, now, 0.4999 * H).should_training_stop  # process started at 0.5 h: 0.0 h of its own clock
    cb2, now2 = budget(2.0, 1.0, start=0.5 * H)
    assert step(cb, now, 0.5 * H).should_training_stop and not step(cb2, now2, 0.5 * H).should_training_stop


@pytest.mark.parametrize("rate,cap,hours", [(0.44, 3.00, 6.8182), (1.0, 2.0, 2.0)])
def test_budget_stop_threshold(rate, cap, hours):
    cb, now = budget(rate, cap)
    assert not step(cb, now, (hours - 1e-4) * H).should_training_stop
    assert step(cb, now, hours * H).should_training_stop


def test_budget_stop_edges():
    cb, now = budget(1.0, 0.0, start=50.0)
    assert step(cb, now, 50.0).should_training_stop  # cap 0: the first step
    cb, now = budget(1.0, 1.0, start=1000.0)
    assert not step(cb, now, 0.0).should_training_stop  # a clock behind the session start never trips
    for rate in (0.0, -1.0):
        with pytest.raises(ValueError, match="rate"):
            budget(rate, 1.0)


def test_budget_stop_fixed_cost_equals_policy_accrued():
    """fixed 15 x inet_down_cost, rate dph + 6 x inet_up_cost (run.sh exports): (2.05 - 0.0975) / 0.68 h = 2.8713 h,
    the time at which LG.accrued of the same run reaches 2.05 (front 6.2)."""
    from objdel.infra.vast import ledger as LG

    run = {"dph": 0.641, "inet_up_cost": 0.0065, "inet_down_cost": 0.0065, "started": 0.0}
    rate, fixed = 0.641 + 6 * 0.0065, 15 * 0.0065
    assert (rate, fixed) == pytest.approx((0.68, 0.0975))
    t = (2.05 - fixed) / rate * H
    assert t / H == pytest.approx(2.8713, abs=1e-4)
    assert LG.accrued(run, t) == pytest.approx(2.05)
    cb, now = budget(rate, 2.05, fixed=fixed)
    assert not step(cb, now, t - 1.0).should_training_stop and step(cb, now, t + 1e-6).should_training_stop


def test_budget_stop_in_a_real_trainer_saves_once_and_stops(tmp_path):
    """I16: a real Trainer and TrainerControl: the trip stops training and gives exactly one checkpoint save."""
    import torch
    from transformers import Trainer, TrainingArguments

    from objdel.policy.callbacks import BudgetStop

    class M(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.w = torch.nn.Parameter(torch.zeros(1))

        def forward(self, x):
            return {"loss": ((self.w - x) ** 2).mean()}

    now = [0.0]

    class Tick(BudgetStop):
        def on_step_end(self, args, state, control, **kw):
            now[0] = state.global_step * H  # one hour per step: $1 a step
            return super().on_step_end(args, state, control, **kw)

    saves = []

    class T(Trainer):
        def _save_checkpoint(self, model, trial):
            saves.append(self.state.global_step)

    args = TrainingArguments(str(tmp_path), max_steps=10, save_strategy="no", report_to="none", use_cpu=True,
                             per_device_train_batch_size=1, logging_strategy="no")
    data = [{"x": torch.ones(1)} for _ in range(20)]
    t = T(model=M(), args=args, train_dataset=data, callbacks=[Tick(1.0, 3.0, 0.0, clock=lambda: now[0])])
    t.train()
    assert t.state.global_step == 3 and saves == [3]


# --- train_rl: rows flags, budget flags, --merge-last (U42, U45, I49 parse part) ----------------------------------

BASE = ["--model", "qwen3-vl-2b", "--ckpt", "c", "--out", "o"]
TRAIN = ["--train", "t", "--steps", "20"]


def parse(*argv, env=None, monkeypatch=None):
    from experiments.phase2 import train_rl

    for k in ("OBJDEL_SESSION_START", "OBJDEL_RATE_USD_PER_H", "OBJDEL_FIXED_USD"):
        monkeypatch.delenv(k, raising=False)
    for k, v in (env or {}).items():
        monkeypatch.setenv(k, v)
    return train_rl.parse([*BASE, *argv])


def refused(capsys, *argv, monkeypatch, env=None) -> str:
    with pytest.raises(SystemExit) as e:
        parse(*argv, env=env, monkeypatch=monkeypatch)
    assert e.value.code == 2
    return capsys.readouterr().err


def test_rows_flags_take_exact_hex(capsys, monkeypatch):
    a = parse(*TRAIN, "--rows-revision", "a" * 40, "--rows-manifest-sha256", "0" * 64, monkeypatch=monkeypatch)
    assert (a.rows_revision, a.rows_manifest_sha256) == ("a" * 40, "0" * 64)
    for bad in (["--rows-revision", "a" * 39], ["--rows-revision", "g" * 40], ["--rows-manifest-sha256", "0" * 63],
                ["--rows-manifest-sha256", "A" * 65]):
        refused(capsys, *TRAIN, *bad, monkeypatch=monkeypatch)


def test_train_views_is_gone(capsys, monkeypatch):
    """U42 (B1.5 side): the old multi-root path is deleted, the one-root path stays."""
    from objdel.rl import rows as rl_rows
    from objdel.verify import store

    assert "--train-views" in refused(capsys, *TRAIN, "--train-views", "v", monkeypatch=monkeypatch)
    assert not hasattr(rl_rows, "grpo_rows_views") and not hasattr(store, "VIEW_SEP")
    assert not hasattr(store, "_tagged")


def test_budget_flags_need_a_session_start_and_a_rate(monkeypatch):
    """U54: the check lives in session_budget (run time on the box), not in the stage: a spec loads on the Mac."""

    for argv, missing in ((["--rate-usd-per-h", "0.68"], "--session-start"),
                          (["--session-start", "100"], "--rate-usd-per-h")):
        a = parse(*TRAIN, "--budget-usd", "2.05", *argv, monkeypatch=monkeypatch)
        with pytest.raises(rl_run.BudgetFlagError, match=missing):
            rl_run.session_budget(a, {})
    a = parse(*TRAIN, "--budget-usd", "2.05", "--session-start", "100", "--rate-usd-per-h", "0.68",
              "--fixed-usd", "0.0975", monkeypatch=monkeypatch)
    assert (a.budget_usd, a.session_start, a.rate_usd_per_h, a.fixed_usd) == (2.05, 100.0, 0.68, 0.0975)
    assert parse(*TRAIN, monkeypatch=monkeypatch).budget_usd is None


def test_budget_flags_default_to_the_run_sh_environment(monkeypatch):
    """I42 (train_rl side), U54: run.sh exports the active run's start, rate and fixed cost; session_budget reads them
    when the flag is left out (the stage holds None), and a flag wins."""

    env = {"OBJDEL_SESSION_START": "1700000000", "OBJDEL_RATE_USD_PER_H": "0.68", "OBJDEL_FIXED_USD": "0.0975"}
    a = parse(*TRAIN, "--budget-usd", "2.05", env=env, monkeypatch=monkeypatch)
    assert a.session_start is None
    cb = rl_run.session_budget(a, env)
    assert (cb.rate, cb.cap, cb.start, cb.fixed) == (0.68, 2.05, 1700000000.0, 0.0975)
    a = parse(*TRAIN, "--budget-usd", "2.05", "--fixed-usd", "0", env=env, monkeypatch=monkeypatch)
    assert rl_run.session_budget(a, env).fixed == 0
    assert rl_run.session_budget(parse(*TRAIN, env=env, monkeypatch=monkeypatch), env) is None
    for bad in ("x", "inf", "-1"):
        with pytest.raises(rl_run.BudgetFlagError, match="OBJDEL_SESSION_START"):
            rl_run.session_budget(a, env | {"OBJDEL_SESSION_START": bad})


def test_merge_last_needs_neither_steps_nor_train(capsys, monkeypatch):
    assert parse("--merge-last", monkeypatch=monkeypatch).merge_last
    refused(capsys, "--train", "t", monkeypatch=monkeypatch)  # a normal run needs --steps
    refused(capsys, "--steps", "20", monkeypatch=monkeypatch)  # and --train


# --- train_rl: schedule guards (U41, U61, I34 train_rl side) ------------------------------------------------------

def make_root(tmp_path, steps=4, per_step=32, block_steps=2, schedule=True):
    ids = [f"i{n}" for n in range(steps * per_step)]
    m = {"files": {}, "sets": {}, "reward_fingerprint": "x", "versions": {}}
    if schedule:
        m["schedule"] = dict(steps=steps, per_step=per_step, block_steps=block_steps,
                             ids_sha256=hashlib.sha256("\n".join(ids).encode()).hexdigest())
    (tmp_path / "manifest.json").write_text(json.dumps(m))
    return tmp_path, ids, hashlib.sha256((tmp_path / "manifest.json").read_bytes()).hexdigest()


def test_schedule_root_needs_the_pinned_manifest(tmp_path):

    root, _, sha = make_root(tmp_path)
    assert rl_run.rows_schedule(root, sha)["per_step"] == 32
    with pytest.raises(ValueError, match="rows-manifest-sha256"):
        rl_run.rows_schedule(root, None)  # a schedule root without the flag
    with pytest.raises(ValueError, match="SHA-256"):
        rl_run.rows_schedule(root, "0" * 64)  # the box loaded a different upload


def test_root_without_schedule_takes_the_old_path(tmp_path):

    root, _, _ = make_root(tmp_path, schedule=False)
    assert rl_run.rows_schedule(root, None) is None


@pytest.mark.parametrize("mutate,match", [
    (lambda k: k[:-1], "rows"),
    (lambda k: k + ["i0"], "rows"),
    (lambda k: k[:-1] + ["x"], "ids_sha256"),  # one id swapped for another
    (lambda k: k[1:2] + k[:1] + k[2:], "ids_sha256"),  # same ids, two reordered
])
def test_rows_must_be_the_schedule_in_order(tmp_path, mutate, match):

    root, ids, sha = make_root(tmp_path)
    schedule = rl_run.rows_schedule(root, sha)
    rl_run.check_rows(ids, schedule)
    with pytest.raises(ValueError, match=match):
        rl_run.check_rows(mutate(list(ids)), schedule)


def sampler(tmp_path, **over):
    """The sampler fields of the real B1 GRPOConfig, as a plain namespace to mutate."""
    from types import SimpleNamespace

    from objdel.rl import config as configs

    cfg = configs.grpo_config("qwen3-vl-2b", str(tmp_path / "trainer"), use_vllm=False)
    keys = ("generation_batch_size", "num_generations", "steps_per_generation", "gradient_accumulation_steps",
            "num_iterations")
    return SimpleNamespace(**{k: getattr(cfg, k) for k in keys} | over)


@pytest.mark.parametrize("over,steps,match", [
    ({"generation_batch_size": 16 * 8}, 4, "per_step"),
    ({"generation_batch_size": 31 * 8}, 4, "per_step"),
    ({"steps_per_generation": 2, "gradient_accumulation_steps": 1}, 4, "steps_per_generation"),
    ({"num_iterations": 2}, 4, "num_iterations"),
    ({}, 6, "steps"),  # more steps than the schedule holds
    ({}, 3, "block"),  # not whole blocks (U61: --steps 30 with block_steps 20)
])
def test_sampler_shape_and_steps_are_refused(tmp_path, over, steps, match):

    root, _, sha = make_root(tmp_path)
    schedule = rl_run.rows_schedule(root, sha)
    rl_run.check_sampler(sampler(tmp_path), schedule, 4)  # the real config fits a 32-prompt schedule
    rl_run.check_sampler(sampler(tmp_path), schedule, 2)  # one whole block of two
    with pytest.raises(ValueError, match=match):
        rl_run.check_sampler(sampler(tmp_path, **over), schedule, steps)


def test_u61_thirty_steps_with_blocks_of_twenty(tmp_path):

    root, _, sha = make_root(tmp_path, steps=40, per_step=32, block_steps=20)
    with pytest.raises(ValueError, match="block"):
        rl_run.check_sampler(sampler(tmp_path), rl_run.rows_schedule(root, sha), 30)


# --- train_rl main on a real rows root: refusals come first (R04), kind / cell columns (U43) ----------------------

@pytest.fixture(scope="module")
def prepared(tmp_path_factory, synthetic_items):
    """The synthetic_items fixture, not `from conftest import`: tests/verify has its own conftest, and the module
    name `conftest` is whichever one pytest loaded last."""
    pytest.importorskip("datasets")
    from objdel.data import prepare as data_prepare
    from objdel.data import store
    from objdel.policy import prompt as policy_prompt

    root = tmp_path_factory.mktemp("rows")
    store.write(synthetic_items(("train", "val"), unanimous=True), root / "store")
    data_prepare.prepare(root / "store", root / "prepared", procs=1, prompts=policy_prompt.PROMPTS)
    return root / "prepared"


def schedule_root(prepared, out, steps=4, per_step=1, block_steps=2):
    """A rows root as build_rl_rows writes it, from a prepared root: lidc/train.jsonl holds steps x per_step window
    records in schedule order (cell and view.kind added; the LIDC train RL items repeated), the manifest lists only
    that file and has `schedule`. Returns (root, keys in order, manifest SHA-256)."""
    import shutil

    shutil.copytree(prepared, out)
    m = json.loads((out / "manifest.json").read_text())
    for rel in list(m["files"]):
        (out / rel).unlink()
    recs = [r for r in map(json.loads, (prepared / "lidc/train.jsonl").read_text().splitlines()) if r["rl"]]
    rows = [recs[i % len(recs)] | {"cell": f"{i % 3}{i % 2}", "view": {"kind": ("shifted", "rescaled")[i % 2]}}
            for i in range(steps * per_step)]
    text = "".join(json.dumps(r) + "\n" for r in rows).encode()
    (out / "lidc/train.jsonl").write_bytes(text)
    keys = [r["item_id"] for r in rows]
    m["files"] = {"lidc/train.jsonl": {"sha256": hashlib.sha256(text).hexdigest()}}
    m["schedule"] = {"steps": steps, "per_step": per_step, "block_steps": block_steps, "files_count": 0,
                     "ids_sha256": hashlib.sha256("\n".join(keys).encode()).hexdigest()}
    (out / "manifest.json").write_text(json.dumps(m, indent=1))
    return out, keys, hashlib.sha256((out / "manifest.json").read_bytes()).hexdigest()


SMOKE = {"per_device_train_batch_size": 8, "gradient_accumulation_steps": 1, "per_device_eval_batch_size": 8,
         "max_completion_length": 8, "use_cpu": True, "bf16": False, "beta": 0.0, "use_vllm": False}


@pytest.fixture
def no_load(monkeypatch):
    """Model and processor both raise: a refusal must come before each of them (TRL creates the vLLM engine inside
    the trainer, after the model; the Mac has no vllm)."""
    def boom(*a, **k):
        raise AssertionError("loaded before the guards ran")

    monkeypatch.setattr("transformers.AutoModelForImageTextToText.from_pretrained", boom)
    monkeypatch.setattr("transformers.AutoProcessor.from_pretrained", boom)


@pytest.mark.parametrize("case,match", [
    ("sha", "SHA-256"), ("noflag", "rows-manifest-sha256"), ("count", "rows"), ("block", "block"),
    ("iterations", "num_iterations"), ("per_step", "per_step"),
])
def test_refusals_come_before_any_load(prepared, tmp_path, no_load, case, match):
    """R04: every train_rl refusal raises ValueError before the model, the processor or a vLLM engine loads."""
    from experiments.phase2 import train_rl

    root, _, sha = schedule_root(prepared, tmp_path / "rows")
    if case == "count":  # the schedule claims 5 steps, the file holds 4
        m = json.loads((root / "manifest.json").read_text())
        m["schedule"]["steps"] = 5
        (root / "manifest.json").write_text(json.dumps(m))
        sha = hashlib.sha256((root / "manifest.json").read_bytes()).hexdigest()
    grpo = SMOKE | ({"num_iterations": 2} if case == "iterations" else {})
    grpo |= {"per_device_train_batch_size": 16} if case == "per_step" else {}
    argv = ["--model", "qwen3-vl-2b", "--ckpt", str(tmp_path / "none"), "--train", str(root), "--out",
            str(tmp_path / "out"), "--steps", "3" if case == "block" else "4",
            "--grpo", json.dumps(grpo)]
    argv += [] if case == "noflag" else ["--rows-manifest-sha256", "0" * 64 if case == "sha" else sha,
                                         "--rows-revision", "a" * 40]  # RlStage: both pins or neither
    with pytest.raises(ValueError, match=match):
        train_rl.main(argv)


def test_window_records_give_kind_and_cell_columns(prepared, tmp_path):
    """U43: rows-root records give kind and cell columns equal to the record; prepared records have none."""
    from objdel.rl import rows as rl_rows

    root, keys, _ = schedule_root(prepared, tmp_path / "rows")
    rows = rl_rows.grpo_rows(root, "lidc", "train")
    assert rows["key"] == keys  # file order kept: the schedule
    assert rows["kind"] == ["shifted", "rescaled", "shifted", "rescaled"] and rows["cell"] == ["00", "11", "20", "01"]
    plain = rl_rows.grpo_rows(prepared, "lidc", "train")
    assert not {"kind", "cell"} & set(plain.column_names)


def test_heldout_holds_val_shape_items_only(prepared, tmp_path, monkeypatch):
    """I44 (selection): the held-out rows come from val-shape patients only; a planted val-select patient is left
    out, and a view with no val-shape item is refused."""
    import shutil


    view = tmp_path / "view"
    shutil.copytree(prepared, view)
    recs = [json.loads(x) for x in (view / "lidc/val.jsonl").read_text().splitlines()]
    planted = next(r for r in recs if r["rl"] and r["item_id"].endswith("z004"))
    planted["group_id"] = "LIDC-IDRI-9999"
    (view / "lidc/val.jsonl").write_text("".join(json.dumps(r) + "\n" for r in recs))
    halves = {"LIDC-IDRI-0002": "shape", "LIDC-IDRI-9999": "select"}
    monkeypatch.setattr("objdel.data.splits.val_half", halves.__getitem__, raising=False)
    rows = rl_run.heldout_rows(view, 64, "lidc")
    assert rows["key"] == ["lidc/1.2.2/n000/z001"] and planted["item_id"] not in rows["key"]
    halves["LIDC-IDRI-0002"] = "select"
    with pytest.raises(ValueError, match="val-shape"):
        rl_run.heldout_rows(view, 64, "lidc")


def test_run_sh_exports_the_session_clock_and_budget_stop_agrees_with_policy(tmp_path, monkeypatch):
    """I42: run.sh exports the active run's start, dph + 6 x inet_up_cost and 15 x inet_down_cost; train_rl takes them
    as its flag defaults; BudgetStop on them gives the session spend that LG.spend() adds to the ledger total."""
    import re

    from objdel.infra.vast import ledger as LG
    from tests.infra.test_vast import Env

    E = Env(tmp_path)
    E.launched()
    run = E.active()["r1"]
    line = E.tmux("experiments.phase2.final_eval", "--out", "x.json")
    env = dict(re.findall(r"(OBJDEL_SESSION_START|OBJDEL_RATE_USD_PER_H|OBJDEL_FIXED_USD)=([0-9.e+-]+)", line))
    assert float(env["OBJDEL_SESSION_START"]) == run["started"]
    assert float(env["OBJDEL_RATE_USD_PER_H"]) == pytest.approx(run["dph"] + 6 * run["inet_up_cost"])
    assert float(env["OBJDEL_FIXED_USD"]) == pytest.approx(15 * run["inet_down_cost"])
    a = parse(*TRAIN, "--budget-usd", "2.05", env=env, monkeypatch=monkeypatch)
    cb = rl_run.session_budget(a, env)
    t = run["started"] + 1.7 * H
    cb.clock = lambda: t
    monkeypatch.setenv("VAST_HOME", E.env["VAST_HOME"])
    (E.home / "ledger.csv").unlink(missing_ok=True)
    LG.ledger(event="teardown", run_id="old", cost_usd=3.25)  # a closed session
    monkeypatch.setattr(LG.time, "time", lambda: t)
    assert cb.spend() == pytest.approx(LG.accrued(run, t)) == pytest.approx(LG.spend() - 3.25)


# --- B2.2: --numerics and --opsm-start-delta --------------------------------------------------------------------


def test_numerics_defaults_to_the_b1_control(monkeypatch):

    a = parse(*TRAIN, monkeypatch=monkeypatch)  # E6.4: an unset field is None, the B1 control resolved at use
    assert a.numerics is None and a.opsm_start_delta is None
    assert rl_run.numerics_overrides(a) == {"numerics": "bf16-seqmask"}


def test_numerics_takes_the_presets_that_need_only_a_config(capsys, monkeypatch):
    """The bf16 presets need only the flag; lora-seqmask turns --f-check on by default, so it also needs the F
    comparator's B1 adapter (--b1-adapter), or --no-f-check, and refuses without either (B2.2 F step 4b.2)."""
    for name in ("bf16-seqmask", "bf16-tis", "bf16-tis-opsm"):
        a = parse(*TRAIN, "--numerics", name, monkeypatch=monkeypatch)
        assert a.numerics == name and a.f_check is False
    a = parse(*TRAIN, "--numerics", "lora-seqmask", "--b1-adapter", "b1", monkeypatch=monkeypatch)
    assert (a.numerics, a.f_check, a.b1_adapter) == ("lora-seqmask", True, "b1")
    a = parse(*TRAIN, "--numerics", "lora-seqmask", "--no-f-check", monkeypatch=monkeypatch)
    assert (a.numerics, a.f_check) == ("lora-seqmask", False)
    assert "--b1-adapter" in refused(capsys, *TRAIN, "--numerics", "lora-seqmask", monkeypatch=monkeypatch)


def test_unknown_numerics_is_refused(capsys, monkeypatch):
    assert "bf16-fp8" in refused(capsys, *TRAIN, "--numerics", "bf16-fp8", monkeypatch=monkeypatch)




def test_fp16_seqmask_parses(monkeypatch):
    assert parse(*TRAIN, "--numerics", "fp16-seqmask", monkeypatch=monkeypatch).numerics == "fp16-seqmask"


def test_opsm_start_delta_only_with_the_opsm_preset(capsys, monkeypatch):
    a = parse(*TRAIN, "--numerics", "bf16-tis-opsm", "--opsm-start-delta", "0.8", monkeypatch=monkeypatch)
    assert a.opsm_start_delta == 0.8
    assert "bf16-tis-opsm" in refused(capsys, *TRAIN, "--opsm-start-delta", "0.8", monkeypatch=monkeypatch)
    assert "> 0" in refused(capsys, *TRAIN, "--numerics", "bf16-tis-opsm", "--opsm-start-delta", "0",
                            monkeypatch=monkeypatch)


def test_the_grpo_json_cannot_set_the_opsm_threshold(capsys, monkeypatch):
    err = refused(capsys, *TRAIN, "--numerics", "bf16-tis-opsm", "--grpo", '{"off_policy_mask_threshold": 0.3}',
                  monkeypatch=monkeypatch)
    assert "--opsm-start-delta" in err
    assert parse(*TRAIN, "--grpo", '{"beta": 0.01}', monkeypatch=monkeypatch).grpo == {"beta": 0.01}


def test_numerics_reach_the_grpo_config(tmp_path, monkeypatch):

    a = parse(*TRAIN, "--numerics", "bf16-tis-opsm", "--opsm-start-delta", "0.8", monkeypatch=monkeypatch)
    assert rl_run.numerics_overrides(a) == {"numerics": "bf16-tis-opsm", "off_policy_mask_threshold": 0.8}
    a = parse(*TRAIN, "--numerics", "bf16-tis-opsm", monkeypatch=monkeypatch)
    assert rl_run.numerics_overrides(a) == {"numerics": "bf16-tis-opsm"}  # the preset's 0.5 stays
