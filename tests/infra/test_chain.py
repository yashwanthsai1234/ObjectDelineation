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
SPECS = ROOT / "specs"
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


# r1_qwen3vl2b.yaml, the SFT launch config that specs/r1_qwen3vl2b.yaml replaced (archived in E6.6, local only),
# minus its runner-only key `sampler`
OLD_R1 = {"model": "qwen3-vl-2b", "view": "views/lidc-bundle-v3", "output_dir": "r1-qwen3vl2b", "seed": 1,
          "branches": (0.5, 1.0, 2.0, 3.0),
          "vision": {"prepared": "prepared", "reference": "vision_ref/qwen3-vl-2b.json"},
          "sft": {"per_device_train_batch_size": 4, "gradient_accumulation_steps": 8, "dataloader_num_workers": 6,
                  "dataloader_persistent_workers": True}}


def r1_spec() -> S.RunSpec:
    return S.load_spec(SPECS / "r1_qwen3vl2b.yaml")


def test_the_sft_launch_config_maps_onto_an_sft_stage_and_the_chain_stages_equal_golden_5():
    spec = r1_spec()
    gold = GOLDEN["r1_qwen3vl2b"]["argv"]  # gpu scorers, train (--config <yaml>), select, probe
    argv = C.chain(SPECS / "r1_qwen3vl2b.yaml", "x")
    assert [argv[0], argv[2], argv[3]] == [gold[0], gold[2], gold[3]]
    got = spec.stages[1].model_dump()
    assert {k: got[k] for k in OLD_R1} == OLD_R1


def test_chain_refuses_an_old_schema_plan(tmp_path):
    """I29: an old `steps:` plan raises ValidationError and an old SFT launch config too, before any remote call."""
    from pydantic import ValidationError

    old = tmp_path / "r2_qwen3vl2b.yaml"
    old.write_text('steps:\n  - ["-m", "experiments.phase2.g_probe", "--model", "qwen3-vl-2b"]\n')
    with pytest.raises(ValidationError):
        C.chain(old, "x")
    cfg = tmp_path / "r1_qwen3vl2b.yaml"
    cfg.write_text(yaml.safe_dump({**OLD_R1, "branches": list(OLD_R1["branches"]), "sampler": "vllm"}))
    with pytest.raises(ValidationError):
        C.chain(cfg, "x")
    assert not list((ROOT / "experiments/phase2/configs").glob("*.yaml"))  # E6.6: the old plans left git
    elsewhere = tmp_path / "b1_rows.yaml"
    elsewhere.write_text((SPECS / "b1_rows.yaml").read_text())
    with pytest.raises(ValueError, match="only a run spec in specs/"):
        C.chain(elsewhere, "x")


# I29: every spec a front 6.3 rental names runs from specs/, never from the old plan folder
RENTAL_PLANS = ("a3_measure", "c_speed", "c_speed_b", "c_speed_c6", "d_gpu_test", "b2_fp16", "a3_test", "b3_v5")


@pytest.mark.parametrize("plan", RENTAL_PLANS)
def test_rental_plans_run_only_from_specs(plan):
    assert C.chain(SPECS / f"{plan}.yaml", "x")
    with pytest.raises(FileNotFoundError):
        C.chain(ROOT / f"experiments/phase2/configs/{plan}.yaml", "x")


def test_vision_refs_of_a_spec_equal_the_old_config():
    assert C.vision_refs(SPECS / "r1_qwen3vl2b.yaml", "x") == GOLDEN["r1_qwen3vl2b"]["vision_refs"]


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
    """B1.4 step 5 pin (rows/lidc-train-b1-780): the plans and specs/b1_rows.yaml agree on one manifest SHA-256."""
    assert policy.rows_pin("d2c8edb8e15534363ebbced8f839cecace49c495") == (
        "a0f2c1d212b8c685e2b85aaac1aa9fa3de9c004ab310af13412d5aa630d4a6b9")
