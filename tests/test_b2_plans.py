"""B2.3 plan file of the optional fp16 arm: the B1 recipe with fp16 numerics, 60 steps, no eval."""

from pathlib import Path

import pytest

pytest.importorskip("yaml")

CONFIGS = Path(__file__).resolve().parents[1] / "specs"


def steps(name: str) -> list[list[str]]:
    from experiments.phase2.vast import policy

    return policy.chain(CONFIGS / name, "r")


def after(step: list[str], flag: str) -> str:
    return step[step.index(flag) + 1]


def test_b2_fp16_is_the_b1_recipe_with_fp16_numerics():
    build, g08b, train = steps("b2_fp16.yaml")
    b1 = next(s for s in steps("b1_rows.yaml") if s[:3] == ["~", "-m", "experiments.phase2.train_rl"])
    assert build == steps("b1_rows.yaml")[0]  # the CUDA build of b1_rows.yaml step 0, else CPU scoring
    assert g08b[-1] == "tests/test_train_fast.py::test_g08b_fp16_dtype" and "gpu" in g08b
    assert train[:2] == ["-m", "experiments.phase2.train_rl"]
    assert after(train, "--numerics") == "fp16-seqmask" and "--numerics" not in b1
    for flag in ("--ckpt", "--train", "--lr", "--rows-revision", "--rows-manifest-sha256"):
        assert after(train, flag) == after(b1, flag)  # same start model, rows root and LR as B1
    assert (after(train, "--steps"), after(train, "--save-steps")) == ("60", "60") and "--shared-prefix" in train
    assert after(train, "--out") == "/workspace/runs/r/b2-fp16"
    assert not [s for s in steps("b2_fp16.yaml") if "experiments.phase2.final_eval" in s]  # no eval
    assert "--grpo" not in train and "--heldout" not in train


def test_the_train_step_parses_with_train_rl():
    pytest.importorskip("trl")
    from experiments.phase2 import train_rl

    a = train_rl.parse(steps("b2_fp16.yaml")[2][2:])
    assert a.numerics == "fp16-seqmask" and a.steps == 60 and a.save_steps == 60 and a.shared_prefix and a.lr == 5e-5
    assert a.log_by_prompt
