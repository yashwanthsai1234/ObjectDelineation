"""E6: the run spec (objdel.core.spec): strict loading, one definition of each flag, and a CLI round trip."""
import json
import sys

import pytest
import yaml
from pydantic import ValidationError

from objdel.core import spec as S
from tests.support.paths import ROOT

SPECS = sorted((ROOT / "specs").glob("*.yaml"))
RD, DATA = "/workspace/runs/x", "/workspace/data/lidc"
HEX40, HEX64 = "a" * 40, "b" * 64

RL = {"kind": "rl", "name": "train", "model": "qwen3-vl-2b", "ckpt": "{rd}/ckpt", "train": "{data}/rows/r",
      "out": "{rd}/o", "steps": 10}
CMD = {"kind": "cmd", "name": "c", "argv": ["echo", "{rd}"]}
EVAL = {"kind": "eval", "name": "e", "ckpt": "{rd}/c", "view": ["{data}/v"], "out": ["{rd}/o"]}


def doc(*stages, **over) -> dict:
    return {"schema_version": 1, "name": "t", "compute": {"gpu": "RTX_5090", "cap_hours": 1.0},
            "bundle_view": "views/b", "stages": list(stages)} | over


def load(tmp_path, d) -> S.RunSpec:
    f = tmp_path / "s.yaml"
    f.write_text(yaml.safe_dump(d))
    return S.load_spec(f)


@pytest.mark.parametrize("path", SPECS, ids=lambda p: p.stem)
def test_every_shipped_spec_loads(path):
    assert S.load_spec(path).name == path.stem


def test_a_valid_spec_loads_and_is_frozen(tmp_path):
    s = load(tmp_path, doc(RL, CMD))
    assert [x.kind for x in s.stages] == ["rl", "cmd"] and isinstance(s.stages[0].steps, int)
    with pytest.raises(ValidationError):
        s.stages[0].steps = 3


@pytest.mark.parametrize("bad", [
    doc(RL | {"typo": 1}),  # unknown key in a stage
    doc(RL | {"kind": "nope"}),  # unknown kind
    doc(RL, extra=1),  # unknown top-level key
    doc(RL, data_revision="x"),  # the rows pin lives on the rl stage
    doc(RL, compute={"gpu": "RTX_5090", "cap_hours": 1.0, "regions": ["EU"]}),  # the region is a fixed rule
    doc(schema_version=2) | {"stages": [RL]},
    doc(RL, schema_version="1"),  # no coercion
    doc(RL | {"steps": "10"}),  # no coercion
    doc(RL | {"steps": 10.0}),
    doc(RL | {"lr": "5e-5"}),  # YAML reads 5e-5 as a string: write 5.0e-5
    doc(),  # stages: []
    doc(RL, RL),  # names are unique
    doc(RL | {"name": ""}),
    doc(RL, bundle_view=None),  # an rl stage without an sft stage needs bundle_view
    doc(RL | {"rows_revision": HEX40}),  # both pins or none
    doc(RL | {"rows_revision": "z" * 40, "rows_manifest_sha256": HEX64}),
    doc(RL | {"rows_revision": HEX40, "rows_manifest_sha256": "b" * 63}),
    doc(RL | {"train": None}),  # train and steps are required unless merge_last
    doc({"kind": "probe", "name": "p", "model": "m", "train": "t", "probe": "p", "out": "o"}),  # ckpt or select
    doc({"kind": "probe", "name": "p", "model": "m", "ckpt": "c", "select": "s", "train": "t", "probe": "p",
         "out": "o"}),
    doc(EVAL | {"ckpt": None}),  # ckpt or rescore
    doc(EVAL | {"split": ["val", "val"]}),  # one split per view or none
    doc(EVAL | {"view": ["a,b"]}),  # a comma splits a repeated flag
    doc(CMD | {"argv": ["echo", "{x}"]}),  # unknown placeholder
    doc({"kind": "select", "name": "s", "run": "r", "prepared": "p", "probes": "q", "shifted": "x"}),  # no vision
])
def test_a_bad_spec_fails_before_any_rental(tmp_path, bad):
    with pytest.raises((ValidationError, ValueError)):
        load(tmp_path, bad)


def test_both_row_pins_pass_and_none_passes(tmp_path):
    load(tmp_path, doc(RL | {"rows_revision": HEX40, "rows_manifest_sha256": HEX64}))
    load(tmp_path, doc(RL))


def test_a_duplicate_yaml_key_fails(tmp_path):
    f = tmp_path / "s.yaml"
    f.write_text("schema_version: 1\nname: a\nname: b\n")
    with pytest.raises(yaml.constructor.ConstructorError, match="duplicate key 'name'"):
        S.load_spec(f)


def test_to_argv_fills_placeholders_and_keeps_literal_braces():
    cmd = S.CmdStage(name="c", argv=("bash", "-c", "{{ echo {rd}; }}", "{data}"))
    assert S.to_argv(cmd, RD, DATA) == ["bash", "-c", "{ echo /workspace/runs/x; }", DATA]
    assert S.to_argv(cmd.model_copy(update={"soft": True}), RD, DATA)[0] == "~"
    with pytest.raises(ValueError, match="unknown placeholder"):
        S.to_argv(S.CmdStage(name="c", argv=("{x}",)), RD, DATA)


def test_to_argv_writes_only_the_flags_the_spec_sets_in_plan_order():
    st = S.RlStage(name="t", model="m", ckpt="{rd}/c", train="{data}/t", out="{rd}/o", steps=260, save_steps=50,
                   shared_prefix=True, lr=5.0e-5, budget_usd=2.05, resume=False)
    assert S.to_argv(st, RD, DATA) == [
        "-m", "objdel.cli.train.rl", "--model", "m", "--ckpt", f"{RD}/c", "--train", f"{DATA}/t", "--out",
        f"{RD}/o", "--steps", "260", "--save-steps", "50", "--shared-prefix", "--lr", "5e-5", "--budget-usd", "2.05"]


def test_floats_are_written_as_plans_write_them():
    got = [S.num(v) for v in (5.0e-5, 5.0e-6, 2.05, 1.0, 1e-3, 1e22)]
    assert got == ["5e-5", "5e-6", "2.05", "1.0", "0.001", "1e22"]


def test_eval_views_interleave_with_their_split_and_out():
    st = S.EvalStage(name="e", ckpt="c", view=("v1", "v2"), split=("val", "test"), out=("o1", "o2"))
    assert S.to_argv(st, RD, DATA)[2:] == ["--ckpt", "c", "--view", "v1", "--split", "val", "--out", "o1", "--view",
                                           "v2", "--split", "test", "--out", "o2"]
    one = st.model_copy(update={"out": ("o",), "split": ()})
    assert S.to_argv(one, RD, DATA)[2:] == ["--ckpt", "c", "--view", "v1", "--view", "v2", "--out", "o"]


def test_flag_order_names_every_field_of_its_stage_once():
    base = set(S.StageBase.model_fields) | {"kind"}
    for kind, order in S.FLAG_ORDER.items():
        fields = set(S.STAGE_CLASSES[kind].model_fields) - base
        names = [f for n in order for f in ({"view", "split", "out"} if n == "views" else {n})]
        assert sorted(names) == sorted(fields), kind


def filled(v):
    """The same value with {rd} and {data} filled in strings, lists and dicts."""
    if isinstance(v, str):
        return v.format(rd=RD, data=DATA)
    if isinstance(v, list):
        return [filled(i) for i in v]
    if isinstance(v, dict):
        return {k: filled(i) if k not in ("sft", "grpo") else i for k, i in v.items()}
    return v


def round_trips(stage) -> bool:
    """to_argv of the stage (without its `~`, which is no CLI token), parsed by its CLI, equals the stage with its
    placeholders filled."""
    argv = S.to_argv(stage.model_copy(update={"soft": False}), RD, DATA)
    assert argv[:2] == ["-m", S.MODULE[stage.kind]]
    want = type(stage).model_validate(filled(stage.model_dump(mode="json")))
    return S.same_stage(S.cli_parse(type(stage), argv[2:]), want)


@pytest.mark.parametrize("path", SPECS, ids=lambda p: p.stem)
def test_every_typed_stage_round_trips_through_its_cli_parser(path):
    assert all(round_trips(st) for st in S.load_spec(path).stages if st.kind != "cmd")


def test_sft_select_probe_stages_round_trip():
    sft = S.SftStage(name="t", model="qwen3-vl-2b", view="views/v", output_dir="r1", seed=1,
                     branches=(0.5, 1.0, 2.0, 3.0), vision={"prepared": "prepared", "reference": "vision_ref/m.json"},
                     sft={"gradient_accumulation_steps": 8})
    sel = S.SelectStage(name="s", run="{rd}/r1", prepared="p", probes="q", shifted="s", sampler="vllm", k=4,
                        vision=S.VisionRef(prepared="{data}/prepared", reference="{rd}/vision_ref/m.json"))
    top = S.SelectStage(name="o", run="r", prepared="p", probes="q", only="topology")
    prb = S.ProbeStage(name="p", model="m", select="s", train="t", probe="q", out="o", shared_prefix=True,
                       grpo={"max_steps": 2})
    off = sft.model_copy(update={"vision": False})
    assert all(round_trips(st) for st in (sft, sel, top, prb, off))
    argv = S.to_argv(off, RD, DATA)
    assert argv[argv.index("--vision") + 1] == "off"


def test_same_stage_ignores_only_name_and_soft():
    a = S.RlStage(name="a", soft=True, model="m", ckpt="c", train="t", out="o", steps=1)
    assert S.same_stage(a, a.model_copy(update={"name": "", "soft": False}))
    assert not S.same_stage(a, a.model_copy(update={"steps": 2}))


CLI_STAGES = {"train_rl": S.RlStage, "g_probe": S.ProbeStage, "final_eval": S.EvalStage, "select_ckpt": S.SelectStage,
              "train_sft": S.SftStage}
# CLEAN8 (full fine-tuning retired): train_sft drops freeze_vision (its D7 repair); select_ckpt gains max_tokens (the
# CPU smoke's completion cap on a random Qwen3-VL). The four stages that read a dataset gain dataset (default lidc).
# The golden is recorded once, so the changes are listed here.
CLEAN8_DROPPED = {"train_sft": {"freeze_vision"}}
CLEAN8_ADDED = {"select_ckpt": {"max_tokens", "dataset"}, "train_rl": {"dataset"}, "g_probe": {"dataset"},
                "final_eval": {"dataset"}}


@pytest.mark.parametrize("cli", CLI_STAGES)
def test_cli_dest_names_equal_stage_fields(cli):
    """U52: every flag of a CLI before E6.4 is a field of its stage, but the dropped ones; a stage adds only the
    trainer flags its CLI did not have (the TrainerFlags mixin). A `negated` dest (argparse BooleanOptionalAction:
    --x and --no-x on one dest) is two bool fields, x and no_x."""
    g = json.loads((ROOT / "tests/goldens/phase_e/cli_dests.json").read_text())[cli]
    dests = set(g["dests"]) | {f"no_{d}" for d in g.get("negated", ())}
    dropped = set(g["dropped"]) | CLEAN8_DROPPED.get(cli, set())
    fields = set(CLI_STAGES[cli].model_fields) - {"kind", "name", "soft"}
    added = set(S.TrainerFlags.model_fields) - dests if cli in ("train_rl", "g_probe") else set()
    added |= CLEAN8_ADDED.get(cli, set())
    assert dests - dropped == fields - added and dropped <= dests and not dropped & fields


def test_each_cli_gains_the_trainer_flags_it_lacked():
    rl = S.cli_parse(S.RlStage, ["--model", "m", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "2",
                                 "--plain", "--compile", "default", "--profile-window", "2,3"])
    assert (rl.plain, rl.compile, rl.profile_window) == (True, "default", (2, 3))
    probe = S.cli_parse(S.ProbeStage, ["--model", "m", "--ckpt", "c", "--train", "t", "--probe", "p", "--out", "o",
                                       "--lr", "5e-5"])
    assert probe.lr == 5e-5
    assert S.RlStage(model="m", ckpt="c", train="t", out="o", steps=2, profile_window=(1, 3)).profile_window == (1, 3)
    for bad in ({"profile_window": (0, 3)}, {"profile_window": (3, 2)}, {"profile_window": (1, 3), "trace": "T"}):
        with pytest.raises(ValidationError, match="profile_window"):
            S.RlStage(model="m", ckpt="c", train="t", out="o", steps=2, **bad)


def test_cli_stage_reads_a_spec_stage_with_its_folders(capsys):
    spec = ROOT / "specs/b1_rows.yaml"
    st = S.cli_stage(S.RlStage, ["--spec", str(spec), "--stage", "train", "--rd", RD, "--data", DATA])
    want = next(s for s in S.load_spec(spec).stages if s.name == "train")
    assert S.same_stage(st, type(want).model_validate(filled(want.model_dump(mode="json"))))
    assert st.out == f"{RD}/b1-rows" and st.name == "" and not st.soft
    for argv in (["--spec", str(spec), "--stage", "merge"], ["--spec", str(spec), "--stage", "nope"]):
        with pytest.raises(SystemExit) as e:
            S.cli_stage(S.EvalStage, argv)  # merge is an rl stage
        assert e.value.code == 2
    with pytest.raises(SystemExit) as e:
        S.cli_stage(S.EvalStage, ["--ckpt", "c", "--rescore", "r", "--view", "v", "--out", "o"])
    assert e.value.code == 2 and "exactly one of ckpt and rescore" in capsys.readouterr().err


def test_a_flag_value_out_of_range_exits_2():
    for argv in (["--steps", "0"], ["--budget-usd", "inf"], ["--session-start", "-1"], ["--save-steps", "0"]):
        with pytest.raises(SystemExit) as e:
            S.cli_stage(S.RlStage, ["--model", "m", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "2",
                                    *argv])
        assert e.value.code == 2, argv


def test_a_bool_field_named_no_x_parses_its_bare_flag_as_true():
    """argparse's BooleanOptionalAction reads any option that starts with --no- as False, so the one parser form uses
    store_true flags (cli_implicit_flags="toggle"); that is the old store_true form while every stage bool defaults
    to False."""
    sel = S.cli_parse(S.SelectStage, ["--run", "r", "--shifted", "s", "--no-vision"])
    assert sel.no_vision is True and S.same_stage(S.cli_parse(S.SelectStage, S.to_argv(sel, "R", "D")[2:]), sel)
    assert all(i.default is False for c in S.STAGE_CLASSES.values() for i in c.model_fields.values()
               if i.annotation is bool)


def test_an_rl_stage_without_save_steps_saves_every_50():
    """E2.7: train_rl's CLI default (save_steps 50) is the stage default."""
    st = S.RlStage(model="m", ckpt="c", train="t", out="o", steps=260)
    assert st.save_steps == 50 and "save_steps" not in st.model_fields_set


def test_an_rl_stage_refuses_save_steps_past_its_steps():
    """E2.7: save_steps > steps would write no checkpoint, so the run could not resume; save_steps == steps is one at
    the end."""
    with pytest.raises(ValueError, match="save_steps"):
        S.RlStage(model="m", ckpt="c", train="t", out="o", steps=20, save_steps=1000)
    assert S.RlStage(model="m", ckpt="c", train="t", out="o", steps=20, save_steps=20).save_steps == 20


def test_run_level_settings_fill_the_stages_that_leave_them_unset(tmp_path):
    """RunSpec model, dataset, rows and compute.budget_usd reach each stage that has the field and leaves it unset,
    unless the value is the field's default (to_argv writes only set fields, so the argv is the same either way); the
    budget reaches the rl stages that train only, rows (train and both pins) those that also leave train unset."""
    rows = {"name": "r2", "revision": HEX40, "manifest_sha256": HEX64}
    bare = {k: v for k, v in RL.items() if k not in ("model", "train")}
    merge = {"kind": "rl", "name": "m", "merge_last": True, "ckpt": "{rd}/c", "out": "{rd}/o"}
    s = load(tmp_path, doc(bare, merge, EVAL, RL | {"name": "own", "model": "x"}, model="m2", dataset="d2", rows=rows,
                           compute={"gpu": "RTX_5090", "cap_hours": 1.0, "budget_usd": 3.0}))
    train, m, ev, own = s.stages
    assert (train.model, train.dataset, train.train, train.rows_revision, train.rows_manifest_sha256,
            train.budget_usd) == ("m2", "d2", "{data}/rows/r2", HEX40, HEX64, 3.0)
    assert (m.model, m.dataset, m.train, m.rows_revision, m.budget_usd) == ("m2", "d2", None, None, None)
    assert (ev.model, ev.dataset) == ("m2", "d2")
    assert (own.model, own.train, own.rows_revision, own.budget_usd) == ("x", "{data}/rows/r", None, 3.0)
    plain, same = (load(tmp_path, doc(RL, EVAL, **kw)) for kw in ({"model": "qwen3-vl-2b", "dataset": "lidc"}, {}))
    assert not {"model", "dataset"} & plain.stages[1].model_fields_set  # the eval stage's defaults: not copied
    assert [S.to_argv(x, RD, DATA) for x in plain.stages] == [S.to_argv(x, RD, DATA) for x in same.stages]


def test_cli_parse_exit_on_error_false_raises_a_value_error_not_a_system_exit(capsys):
    ok = ["--model", "m", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "2"]
    with pytest.raises(SystemExit) as e:  # the default: argparse's exit 2 on an unknown flag
        S.cli_parse(S.RlStage, [*ok, "--no-such-flag"])
    assert e.value.code == 2
    for bad in (["--no-such-flag"], ["--steps", "x"], ["--steps", "0"]):
        with pytest.raises(ValueError):  # pydantic-settings' SettingsError or a ValidationError
            S.cli_parse(S.RlStage, [*ok, *bad], exit_on_error=False)
    assert S.cli_parse(S.RlStage, ok, exit_on_error=False).steps == 2
    capsys.readouterr()


def test_cli_parse_and_cli_stage_read_sys_argv_without_argv(monkeypatch):
    ok = ["--model", "m", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "3"]
    monkeypatch.setattr(sys, "argv", ["prog", *ok])
    assert S.cli_parse(S.RlStage).steps == 3 and S.cli_stage(S.RlStage).steps == 3
