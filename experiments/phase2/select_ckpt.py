"""D12 checkpoint selection CLI (PLAN D12, D13, §3.2 ②, §5.1): objdel.sft.select.run scores the branches
and applies rule R-B and the gates with the experiment modules this script passes; never touches CUDA.

    python -I -m experiments.phase2.select_ckpt --run runs/r1 --shifted data/eval_views/lidc-val-shifted \
        --vision '{"prepared": root, "reference": R0 base value}' | --no-vision \
        [--sampler vllm|hf] [--items 256] [--rescaled V] [--rotated V] [--prepared P] [--probes DIR]
    python -I -m experiments.phase2.select_ckpt --run runs/c2 --probes DIR --only topology [--sampler hf]
    python -I -m experiments.phase2.select_ckpt --spec specs/r1_qwen3vl2b.yaml --stage select [--rd RD] [--data DATA]

The flags are an objdel.core.spec.SelectStage (cli_stage): K is -k and the sampler timeout --timeout-s; --max-tokens
caps every sampled completion below the model's budget (the CPU smoke on a random model, which never stops). The order of
work, the sampler child and the outputs: objdel.sft.select.
"""

import sys

from experiments.phase2 import baselines, prereg
from objdel.policy.errors import Refusal


def main(argv=None) -> int:
    """CLI: the SelectStage of the flags (or of --spec/--stage) -> objdel.sft.select.run with prereg and baselines."""
    from objdel.core.spec import SelectStage, cli_stage
    from objdel.sft.select import run

    try:
        return run(cli_stage(SelectStage, sys.argv[1:] if argv is None else argv), prereg, baselines)
    except Refusal as e:
        raise e.system_exit() from None


if __name__ == "__main__":
    sys.exit(main())
