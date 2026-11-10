"""B3.9 config files (I22b): the V3 launch config and the b3_v5.yaml steps in order; every
`-m experiments.phase2.<module>` step parses with its module's parser on the Mac; against B1 the train_rl step
changes only the reward, the truncation mask, beta and adam_epsilon, and the output folder."""
import json

import pytest
from test_b1_plans import CONFIGS, after, module_steps, parse_step, steps

yaml = pytest.importorskip("yaml")


def test_b3_v5_spec_holds_what_b3_launch_held():
    """E6.6: b3_launch.yaml is archived; the spec has b1's bundle_view and pushes its rl out, revision b3-v5."""
    from experiments.phase2.vast import policy
    from objdel.core.spec import load_spec

    spec = load_spec(CONFIGS / "b3_v5.yaml")
    assert spec.bundle_view == load_spec(CONFIGS / "b1_rows.yaml").bundle_view
    assert policy.push_target(spec, "x")[1:] == ("b3-v5", True)


def test_b3_v5_holds_the_b3_9_steps_in_order():
    build, smoke, b1_eval, train, merge, evals, rc0 = steps("b3_v5.yaml")
    assert build == steps("b1_rows.yaml")[0]
    assert smoke[:2] == ["bash", "-c"] and smoke[2].startswith("OBJDEL_DATA_ROOT=")
    assert smoke[2].endswith("tests/test_train_fast.py::test_g08b_smoke")
    assert after(b1_eval, "--ckpt").endswith("/b1/final") and "eval_b1v3_val_shifted" in after(b1_eval, "--out")
    assert train[0] == "~" and merge[0] != "~" and "--merge-last" in merge and after(merge, "--out") == after(
        train, "--out")
    assert after(evals, "--ckpt") == after(train, "--out") + "/final"
    assert rc0[0] == "~" and rc0[-1] == "tests/test_c_gpu.py::test_g_rc0_per_row_caps"


def test_b3_train_step_differs_from_b1_only_where_the_plan_says():
    b1, b3 = steps("b1_rows.yaml")[3], steps("b3_v5.yaml")[3]
    assert after(b3, "--reward") == "v5" and "--numerics" not in b3 and "--runaway-cap" not in b3
    grpo = json.loads(after(b3, "--grpo"))
    assert grpo["mask_truncated_completions"] is False and grpo["beta"] * 1e-5 == pytest.approx(grpo["adam_epsilon"])
    drop = {"--reward", "--grpo", "--out", "--budget-usd"}
    strip = lambda s: [w for i, w in enumerate(s) if w not in drop and (i == 0 or s[i - 1] not in drop)]  # noqa: E731
    assert strip(b1) == strip(b3)


def test_every_b3_plan_step_parses_with_its_module():
    pytest.importorskip("trl")
    found = module_steps("b3_v5.yaml")
    assert len(found) == 4
    args = [parse_step(s) for s in found]
    assert args[1].reward == "v5" and args[1].steps == 260
