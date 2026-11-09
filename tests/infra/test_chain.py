"""E6: policy.chain() on a run spec gives, string for string, what the old plan files gave (golden 5), and the
G-vision references and the SFT launch config map onto spec stages."""
import json
from pathlib import Path

import pytest
import yaml

from experiments.phase2.vast import policy
from objdel.core import spec as S

ROOT = Path(__file__).resolve().parents[2]
CONFIGS, SPECS = ROOT / "experiments/phase2/configs", ROOT / "specs"
GOLDEN = json.loads((ROOT / "tests/goldens/phase_e/invariance.json").read_text())["chain_argv"]


@pytest.mark.parametrize("plan", ["b1_rows", "rl_full", "r5_shared"])
def test_spec_chain_equals_golden_5(plan):
    assert policy.chain(SPECS / f"{plan}.yaml", "x") == GOLDEN[plan]["argv"]


def test_r5_shared_keeps_its_soft_steps():
    argv = policy.chain(SPECS / "r5_shared.yaml", "x")
    assert [a[0] == "~" for a in argv] == [True, False, False, True]


@pytest.mark.parametrize("plan", ["a3_measure"])
def test_spec_chain_equals_the_old_plan(plan):
    assert policy.chain(SPECS / f"{plan}.yaml", "x") == policy.chain(CONFIGS / f"{plan}.yaml", "x")


def test_a3_test_differs_from_the_old_plan_only_in_where_the_dropped_train_flag_stands():
    new, old = (policy.chain(d / "a3_test.yaml", "x") for d in (SPECS, CONFIGS))
    assert new[0] == old[0] and sorted(new[1]) == sorted(old[1])


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
    gpu = S.CmdStage(name="gpu_scorers", argv=tuple(policy.chain(CONFIGS / "r1_qwen3vl2b.yaml", "x")[0]))
    return S.RunSpec(schema_version=1, name="r1_qwen3vl2b", compute=S.ComputeSpec(gpu="RTX_5090", cap_hours=1.7),
                     stages=(gpu, train, select, probe))


def test_the_sft_launch_config_maps_onto_an_sft_stage_and_the_chain_stages_equal_golden_5():
    spec = r1_spec()
    gold = GOLDEN["r1_qwen3vl2b"]["argv"]  # gpu scorers, train (--config <yaml>), select, probe
    argv = [S.to_argv(s, "/workspace/runs/x", policy.DATA) for s in spec.stages]
    assert [argv[0], argv[2], argv[3]] == [gold[0], gold[2], gold[3]]
    raw = yaml.safe_load((CONFIGS / "r1_qwen3vl2b.yaml").read_text())
    got = spec.stages[1].model_dump()
    assert {k: got[k] for k in raw if k != "sampler"} == {k: tuple(v) if k == "branches" else v
                                                          for k, v in raw.items() if k != "sampler"}


def test_vision_refs_of_a_spec_equal_the_old_config(tmp_path):
    f = tmp_path / "r1.yaml"
    f.write_text(yaml.safe_dump(json.loads(r1_spec().model_dump_json(exclude={"stages": {"__all__": {"soft"}}}))))
    assert policy.vision_refs(f, "x") == GOLDEN["r1_qwen3vl2b"]["vision_refs"]


def test_push_target_takes_the_last_rl_or_sft_stage():
    b1 = S.load_spec(SPECS / "b1_rows.yaml")
    assert policy.push_target(b1, "x") == ("/workspace/runs/x/b1-rows", "b1-rows", True)
    sft = b1.model_copy(update={"stages": (r1_spec().stages[1],), "bundle_view": None})
    assert policy.push_target(sft, "x") == ("/workspace/runs/x/r1-qwen3vl2b", "r1-qwen3vl2b", False)


def test_push_target_of_a_spec_without_a_training_stage_dies():
    spec = S.load_spec(SPECS / "a3_test.yaml")
    with pytest.raises(SystemExit):
        policy.push_target(spec, "x")


def test_rows_pin_reads_the_specs():
    zero40 = "0" * 40
    assert policy.rows_pin(zero40) == "0" * 64
