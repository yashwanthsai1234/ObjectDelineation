"""G-probe (PLAN §5.1, §3.2 ④): can GRPO improve this checkpoint? 20 GRPO steps at the pre-registered probe LR.

    python -I -m objdel.cli.train.probe --model qwen3-vl-2b (--ckpt runs/r1/branches/e1 | --select runs/r1/select)
        --train data/prepared --probe data/eval_views/lidc-val-shifted --out runs/r1/g_probe
    python -I -m objdel.cli.train.probe --spec specs/r1_qwen3vl2b.yaml --stage probe [--rd RD] [--data DATA]

The flags are an objdel.core.spec.ProbeStage (cli_stage), with the trainer flags train_rl has (TrainerFlags: --lr,
--grpo JSON, --shared-prefix, --compile, --plain, --trace, the C flags); objdel.rl.build.build_trainer builds the
trainer from them.

The probe itself (trainer, step-0 null, verdict, the D12(c) repair, <out>/g_probe.json): objdel.rl.probe, which takes
the pre-registered record (prereg: G_PROBE, g_probe) from here.
"""

import sys

from objdel.cli.train import prereg


def main(argv=None) -> int:
    """CLI: the ProbeStage of the flags (or of --spec/--stage) -> objdel.rl.probe.run with the pre-registered record."""
    from objdel.core.spec import ProbeStage, cli_stage
    from objdel.rl.probe import run

    return run(cli_stage(ProbeStage, sys.argv[1:] if argv is None else argv), prereg)


if __name__ == "__main__":
    sys.exit(main())
