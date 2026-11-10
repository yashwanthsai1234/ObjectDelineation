"""Full GRPO run CLI: the RlStage of the flags -> objdel.rl.run.run_rl (or merge_last with --merge-last).

    python -I -m experiments.phase2.train_rl --model qwen3-vl-2b --ckpt runs/r1/branches/e2 --train data/prepared
        --out runs/rl --steps 500 --save-steps 50 [--resume]
    python -I -m experiments.phase2.train_rl --model qwen3-vl-2b --ckpt {rd}/ckpt --train {data}/rows/lidc-train-b1
        --rows-revision <hex40> --rows-manifest-sha256 <hex64> --out {rd}/b1-rows --steps 260 --budget-usd 2.05
    python -I -m experiments.phase2.train_rl --merge-last --model qwen3-vl-2b --ckpt {rd}/ckpt --out {rd}/b1-rows

The flags are an objdel.core.spec.RlStage (cli_stage: one flag per field, --grpo as JSON, --profile-window S,E with
1 <= S <= E: `1,STEPS` traces every step, S >= 2 with --trace); parse adds the rules that need objdel.rl (the
numerics preset names, the --grpo keys train_rl owns). The run itself (rows, saves, budget, held-out curve,
rl.json): objdel.rl.run.
"""

import sys


def parse(argv=None):
    """The RlStage of the flags (or of --spec/--stage), and the rules that need objdel.rl; exit 2 on a broken one.
    The returned stage has f_check resolved: --f-check, --no-f-check, else on under the f_check.PRESETS."""
    from objdel.core.spec import RlStage, cli_stage
    from objdel.rl import config as configs
    from objdel.rl import f_check

    a = cli_stage(RlStage, sys.argv[1:] if argv is None else argv)

    def error(msg: str):
        print(f"train_rl: error: {msg}", file=sys.stderr)
        sys.exit(2)

    numerics = a.numerics or configs.DEFAULT_NUMERICS
    if numerics not in configs.NUMERICS:
        error(f"--numerics {a.numerics}: one of {sorted(configs.NUMERICS)}")
    if a.f_check and numerics not in f_check.PRESETS:
        error(f"--f-check: only under --numerics {' or '.join(f_check.PRESETS)}, not {numerics}")
    a = a.model_copy(update={"f_check": a.f_check or (numerics in f_check.PRESETS and not a.no_f_check)})
    if a.f_check and numerics == "lora-seqmask" and not a.merge_last and a.b1_adapter is None:
        error("--f-check under --numerics lora-seqmask needs --b1-adapter (the comparator's B1 adapter)")
    if a.opsm_start_delta is not None and (numerics != "bf16-tis-opsm" or a.opsm_start_delta <= 0):
        error("--opsm-start-delta is a float > 0, valid only with --numerics bf16-tis-opsm")
    if "off_policy_mask_threshold" in a.grpo:
        error("--grpo cannot set off_policy_mask_threshold: use --opsm-start-delta")
    if a.grpo.get("logging_steps", 1) != 1:  # seq_gaps and cut_seq_gaps need one block per step
        error("--grpo logging_steps: must stay 1 (seq_gaps.f32 holds one block per optimizer step)")
    return a


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
