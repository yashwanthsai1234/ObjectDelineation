"""B1.5 config files (I22a) and the plan argv check (I47): the V1 launch config, the b1_rows.yaml steps in order, and
every `-m experiments.phase2.<module>` step of a V1 plan parsed with its module's own parser on the Mac, so a missing
or unknown flag fails here and not on the paid box."""

import subprocess
import sys
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

ROOT = Path(__file__).resolve().parents[1]
NEEDS_ARCHIVE = pytest.mark.skipif(not (ROOT / "archive").exists(), reason="archive/ is local-only (git-ignored)")
CONFIGS = ROOT / "specs"  # E6.2: the run specs replaced the plan files
V1_APPROVED_USD = 2.40  # front 6.4 rule 5: the owner's V1 amount ...
TAIL_USD = 0.35  # ... minus the tail after train_rl (evals, pull, push, teardown); no later V1 item
V1_PLANS = ("b1_rows.yaml", "a3_measure.yaml")


def steps(name: str) -> list[list[str]]:
    from objdel.infra.vast import chain as C

    return C.chain(CONFIGS / name, "r")


def after(step: list[str], flag: str) -> str:
    return step[step.index(flag) + 1]


def test_b1_rows_spec_holds_what_b1_launch_held():
    """E6.6: b1_launch.yaml (view, output_dir, rows_plan) is archived; the spec's bundle_view is the r1 view, its rl
    out is the push revision b1-rows, and preflight reads the rows pin from the rl stage itself."""
    from objdel.core.spec import load_spec
    from objdel.infra.vast import chain as C

    spec, r1 = load_spec(CONFIGS / "b1_rows.yaml"), load_spec(CONFIGS / "r1_qwen3vl2b.yaml")
    assert spec.bundle_view == next(s.view for s in r1.stages if s.kind == "sft")
    assert C.push_target(spec, "x")[1:] == ("b1-rows", True)
    train = next(s for s in spec.stages if s.kind == "rl" and not s.merge_last)
    assert (train.rows_revision, train.rows_manifest_sha256) == ("34213dceb068128f1a1abf9c67aa65dd59bc1220",
                                                                 "472738141f47bb0f1bdb1e98c679076cb73b979a4b8e256d7e2051c63e352fb4")


def test_b1_rows_holds_the_b1_5_steps_in_order():
    build, env, smoke, train, merge, evals = steps("b1_rows.yaml")
    assert build[:2] == ["bash", "-c"] and "run_box.sh --build-only" in build[2]
    assert "assert g.build() == 'cuda'" in build[2] and build[2].endswith("tests/test_gpu_scoring.py")
    assert env[:2] == ["bash", "experiments/phase2/docker/env_check.sh"]
    assert smoke[:4] == ["-m", "pytest", "-m", "gpu"] and "-k" not in smoke
    assert train[:3] == ["~", "-m", "experiments.phase2.train_rl"]  # soft: a crash goes on to --merge-last
    t = dict(zip(train[3::2], train[4::2], strict=False))
    assert t["--model"] == "qwen3-vl-2b" and t["--ckpt"] == "/workspace/runs/r/ckpt"
    assert after(train, "--train") == "/workspace/data/lidc/rows/lidc-train-b1"
    assert after(train, "--out") == "/workspace/runs/r/b1-rows" and after(train, "--steps") == "260"
    assert after(train, "--save-steps") == "50" and float(after(train, "--lr")) == 5e-5
    assert "--shared-prefix" in train and "--numerics" not in train and "--grpo" not in train
    assert len(after(train, "--rows-revision")) == 40 and len(after(train, "--rows-manifest-sha256")) == 64
    assert float(after(train, "--budget-usd")) == pytest.approx(V1_APPROVED_USD - TAIL_USD)
    from objdel.train import configs

    cfg = configs.grpo_config("qwen3-vl-2b", "o", use_vllm=False)  # no --grpo: beta and masking are the defaults
    assert cfg.beta == 1e-3 and cfg.mask_truncated_completions is True
    assert merge[:3] == ["-m", "experiments.phase2.train_rl", "--merge-last"]  # hard
    assert after(merge, "--out") == after(train, "--out") and after(merge, "--ckpt") == after(train, "--ckpt")
    assert evals[:2] == ["-m", "experiments.phase2.final_eval"]
    assert after(evals, "--ckpt") == "/workspace/runs/r/b1-rows/final"
    assert after(evals, "--train") == "/workspace/data/lidc/prepared"
    views = [evals[i + 1].rsplit("/", 1)[1] for i, w in enumerate(evals) if w == "--view"]
    assert views == ["lidc-val-shifted", "lidc-val-edge"]
    assert [evals[i + 1] for i, w in enumerate(evals) if w == "--split"] == ["val", "val"]
    assert not [s for s in steps("b1_rows.yaml") if "experiments.phase2.final_eval" in s
                and after(s, "--ckpt") == "/workspace/runs/r/ckpt"]  # R5 is a3_measure.yaml, not this plan


def test_smoke_node_ids_collect_exactly_two_tests():
    smoke = steps("b1_rows.yaml")[2]
    ids = [w for w in smoke if "::" in w]
    r = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider", *ids],
                       cwd=ROOT, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    assert [ln for ln in r.stdout.splitlines() if "::" in ln] == ids


@NEEDS_ARCHIVE
def test_rl_shifted_is_archived_with_its_readme_row():
    name = "rl_shifted" + ".yaml"  # built, so this file does not match its own scan below
    assert not (CONFIGS / name).exists()
    assert (ROOT / "archive/phase2_runs/configs" / name).exists()
    readme = (ROOT / "archive/phase2_runs/configs/README.md").read_text()
    assert "| Was | Replaced by | Commit |" in readme and f"| {name} | b1_rows.yaml |" in readme
    files = subprocess.run(["git", "ls-files", "experiments", "src", "tests", "Makefile"], cwd=ROOT,
                           capture_output=True, text=True, check=True).stdout.split()
    assert [f for f in files if (ROOT / f).is_file() and name in (ROOT / f).read_text(errors="ignore")] == []


# --- I47: plan argv parse ----------------------------------------------------------------------------------------


def parse_step(step: list[str]):
    """The step's argv as its module's stage: objdel.core.spec.cli_stage, the parser each module's main() runs first
    (exit 2 on a missing or unknown flag)."""
    from objdel.core.spec import MODULE, STAGE_CLASSES, cli_stage

    words = step[1:] if step[0] == "~" else step
    assert words[0] == "-m"
    if words[1] == "experiments.phase2.c_speed_check":  # its own --pixels/--sync, then a train_rl stage
        from experiments.phase2 import c_speed_check

        return c_speed_check.parse(words[2:])[1]
    kind = next(k for k, m in MODULE.items() if words[1] == f"experiments.phase2.{m}")
    return cli_stage(STAGE_CLASSES[kind], words[2:])


def module_steps(name: str) -> list[list[str]]:
    return [s for s in steps(name) if "-m" in s[:2] and s[s.index("-m") + 1].startswith("experiments.phase2.")]


@pytest.mark.parametrize("plan", V1_PLANS)
def test_every_plan_step_parses_with_its_module(plan):
    pytest.importorskip("trl")
    found = module_steps(plan)
    assert found
    for s in found:
        parse_step(s)


def test_a_missing_required_flag_fails_the_parse():
    pytest.importorskip("trl")
    train, merge, ev = (module_steps("b1_rows.yaml")[i] for i in (0, 1, 2))
    assert parse_step(merge).merge_last and parse_step(train).steps == 260
    no_model = train[:train.index("--model")] + train[train.index("--model") + 2:]
    no_train = ev[:ev.index("--train")] + ev[ev.index("--train") + 2:]
    with pytest.raises(SystemExit) as e:
        parse_step(no_model)
    assert e.value.code == 2
    from experiments.phase2 import final_eval  # --train is checked after parse_args (it is optional with --rescore)
    with pytest.raises(SystemExit) as e:
        final_eval.main(no_train[2:])
    assert e.value.code == 2


def test_every_spec_formats():
    """policy.chain fills {rd}/{data} in every stage: a bash group `{ ...; }` raises on the Mac, `{{ ... }}` is fine."""
    from objdel.infra.vast import chain as C

    specs = sorted(CONFIGS.glob("*.yaml"))
    assert {"b1_rows.yaml", "rl_full.yaml", "impl2_gpu_test.yaml"} <= {p.name for p in specs}
    for p in specs:
        C.chain(p, "r")
    bad = ROOT / "build" / "test_cfgs" / "group.yaml"
    bad.parent.mkdir(parents=True, exist_ok=True)
    bad.write_text("schema_version: 1\nname: g\ncompute: {gpu: RTX_5090, cap_hours: 1.0}\n"
                   'stages:\n- {kind: cmd, name: g, argv: [bash, -c, "{ true; }"]}\n')
    with pytest.raises(ValueError):
        C.chain(bad, "r")
