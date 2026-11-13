"""Full GRPO run CLI: the RlStage of the flags -> objdel.rl.run.run_rl (or merge_last with --merge-last).

    python -I -m objdel.cli.train.rl --model qwen3-vl-2b --ckpt runs/r1/branches/e2 --train data/prepared
        --out runs/rl --steps 500 --save-steps 50 [--resume]
    python -I -m objdel.cli.train.rl --model qwen3-vl-2b --ckpt {rd}/ckpt --train {data}/rows/lidc-train-b1
        --rows-revision <hex40> --rows-manifest-sha256 <hex64> --out {rd}/b1-rows --steps 260 --budget-usd 2.05
    python -I -m objdel.cli.train.rl --merge-last --model qwen3-vl-2b --ckpt {rd}/ckpt --out {rd}/b1-rows

The flags are an objdel.core.spec.RlStage (cli_stage: one flag per field, --grpo as JSON, --profile-window S,E with
1 <= S <= E: `1,STEPS` traces every step, S >= 2 with --trace); objdel.rl.run.checked adds the rules that need objdel.rl
(the numerics preset names, the --grpo keys train_rl owns). The run itself (rows, saves, budget, held-out curve,
rl.json): objdel.rl.run.
"""

import sys


def parse(argv=None):
    """The RlStage of the flags (or of --spec/--stage) through objdel.rl.run.checked; exit 2 on a broken one."""
    from objdel.core.spec import RlStage, cli_stage
    from objdel.rl.run import checked

    return checked(cli_stage(RlStage, sys.argv[1:] if argv is None else argv))


def main(argv=None) -> int:
    """CLI: the parsed RlStage -> run_rl (or merge_last); exit 2 on a budget without a session start or rate, 3 when
    step 1 differs from --expect-step1."""
    from objdel.rl.run import BudgetFlagError, Step1Mismatch, merge_last, run_rl

    a = parse(argv)
    if a.merge_last:
        return merge_last(a)
    try:
        run_rl(a)
    except BudgetFlagError as e:
        print(f"train_rl: error: {e}", file=sys.stderr)
        return 2
    except Step1Mismatch as e:
        print(e, file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
