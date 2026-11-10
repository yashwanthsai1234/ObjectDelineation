"""E6: C.chain() on a run spec gives, string for string, what the old plan files gave (golden 5), and the
G-vision references and the SFT launch config map onto spec stages."""
import json
from pathlib import Path

import pytest
import yaml

from objdel.core import spec as S
from objdel.infra.vast import __main__ as policy
from objdel.infra.vast import chain as C
from tests.experiments.test_move_invariance import later_edits

ROOT = Path(__file__).resolve().parents[2]
CONFIGS, SPECS = ROOT / "experiments/phase2/configs", ROOT / "specs"
GOLDEN = json.loads((ROOT / "tests/goldens/phase_e/invariance.json").read_text())["chain_argv"]
for _g in GOLDEN.values():
    _g["argv"] = later_edits(_g["argv"])


@pytest.mark.parametrize("plan", ["b1_rows", "rl_full", "r5_shared"])
def test_spec_chain_equals_golden_5(plan):
    assert C.chain(SPECS / f"{plan}.yaml", "x") == GOLDEN[plan]["argv"]


def test_r5_shared_keeps_its_soft_steps():
    argv = C.chain(SPECS / "r5_shared.yaml", "x")
    assert [a[0] == "~" for a in argv] == [True, False, False, True]


PLAN_ARGV = json.loads((ROOT / "tests/goldens/phase_e/plan_argv.json").read_text())
PLAN_ARGV = {plan: later_edits(argvs) for plan, argvs in PLAN_ARGV.items()}
# (plan, step) pairs whose old text a spec cannot write string for string; the stage they parse to is the same
SAME_STAGE = {("b3_v5", 3): "grpo JSON spacing and 2.20",
              ("c_speed_c6", 1): "the rows pins stood before --out"}


@pytest.mark.parametrize("plan", sorted(PLAN_ARGV))
def test_spec_chain_equals_the_plan_it_replaced(plan):
    """E6.2: each ported spec runs the argv of the steps plan it replaced (recorded before E6.6 archived them)."""
    new, old = C.chain(SPECS / f"{plan}.yaml", "x"), PLAN_ARGV[plan]
    assert len(new) == len(old)
    for i, (a, b) in enumerate(zip(new, old, strict=True)):
        if (plan, i) not in SAME_STAGE:
            assert a == b, (plan, i)
            continue
        assert a != b, f"{plan} step {i} is equal now: drop it from SAME_STAGE"
        a, b = (x[x[0] == "~":] for x in (a, b))
        cls = next(S.STAGE_CLASSES[k] for k, m in S.MODULE.items() if a[1] == f"experiments.phase2.{m}")
        assert S.same_stage(S.cli_parse(cls, a[2:]), S.cli_parse(cls, b[2:])), (plan, i)


def test_every_spec_is_ported_or_checked_elsewhere():
    assert {p.stem for p in SPECS.glob("*.yaml")} == set(PLAN_ARGV) | {"r1_qwen3vl2b"}


def r1_spec() -> S.RunSpec:
    """r1_qwen3vl2b.yaml (the SFT launch config) and its chain as one spec (the form E6.2 ports)."""
    raw = yaml.safe_load((CONFIGS / "r1_qwen3vl2b.yaml").read_text())
    train = S.SftStage(name="train", **{k: tuple(v) if k == "branches" else v for k, v in raw.items()
                                        if k not in ("sampler", "stages")})
    ev = "{data}/eval_views/"
    select = S.SelectStage(
        name="select", run="{rd}/r1-qwen3vl2b", prepared="{data}/prepared", probes="{data}/eval_views", sampler="vllm",
        shifted=ev + "lidc-val-shifted", rescaled=ev + "lidc-val-rescaled", rotated=ev + "lidc-val-rotated",
        vision=S.VisionRef(prepared="{data}/prepared", reference="{rd}/vision_ref/qwen3-vl-2b.json"))
    probe = S.ProbeStage(name="probe", model="qwen3-vl-2b", select="{rd}/r1-qwen3vl2b/select",
                         train="{data}/prepared", probe=ev + "lidc-val-shifted", out="{rd}/r1-qwen3vl2b/g_probe")
    gpu = S.CmdStage(name="gpu_scorers", argv=tuple(C.chain(CONFIGS / "r1_qwen3vl2b.yaml", "x")[0]))
    return S.RunSpec(schema_version=1, name="r1_qwen3vl2b", compute=S.ComputeSpec(gpu="RTX_5090", cap_hours=1.7),
                     stages=(gpu, train, select, probe))


def test_the_sft_launch_config_maps_onto_an_sft_stage_and_the_chain_stages_equal_golden_5():
    spec = r1_spec()
    gold = GOLDEN["r1_qwen3vl2b"]["argv"]  # gpu scorers, train (--config <yaml>), select, probe
    argv = [S.to_argv(s, "/workspace/runs/x", C.DATA) for s in spec.stages]
    assert [argv[0], argv[2], argv[3]] == [gold[0], gold[2], gold[3]]
    raw = yaml.safe_load((CONFIGS / "r1_qwen3vl2b.yaml").read_text())
    got = spec.stages[1].model_dump()
    assert {k: got[k] for k in raw if k != "sampler"} == {k: tuple(v) if k == "branches" else v
                                                          for k, v in raw.items() if k != "sampler"}


def test_vision_refs_of_a_spec_equal_the_old_config(tmp_path):
    f = tmp_path / "r1.yaml"
    f.write_text(yaml.safe_dump(json.loads(r1_spec().model_dump_json(exclude={"stages": {"__all__": {"soft"}}}))))
    assert C.vision_refs(f, "x") == GOLDEN["r1_qwen3vl2b"]["vision_refs"]


def test_push_target_takes_the_last_rl_or_sft_stage():
    b1 = S.load_spec(SPECS / "b1_rows.yaml")
    assert C.push_target(b1, "x") == ("/workspace/runs/x/b1-rows", "b1-rows", True)
    sft = b1.model_copy(update={"stages": (r1_spec().stages[1],), "bundle_view": None})
    assert C.push_target(sft, "x") == ("/workspace/runs/x/r1-qwen3vl2b", "r1-qwen3vl2b", False)


def test_push_target_of_a_spec_without_a_training_stage_dies():
    spec = S.load_spec(SPECS / "a3_test.yaml")
    with pytest.raises(SystemExit):
        C.push_target(spec, "x")


def test_rows_pin_reads_the_specs():
    """B1.4 step 5 pin (rows/lidc-train-b1): the plans and specs/b1_rows.yaml agree on one manifest SHA-256."""
    assert policy.rows_pin("34213dceb068128f1a1abf9c67aa65dd59bc1220") == (
        "472738141f47bb0f1bdb1e98c679076cb73b979a4b8e256d7e2051c63e352fb4")
