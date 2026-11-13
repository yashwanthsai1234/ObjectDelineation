"""Phase 2 SFT CLI: the SftStage of the flags -> objdel.sft.train.train (the run, its schedule and manifest).

    python -I -m objdel.cli.train.sft --model qwen3-vl-2b --view views/lidc-bundle-v3 --output-dir r1 \
        --vision '{"prepared": "prepared", "reference": "vision_ref/qwen3-vl-2b.json"}' [--resume]
    python -I -m objdel.cli.train.sft --spec specs/r1_qwen3vl2b.yaml --stage train [--rd RD] [--data DATA]
    python -I -m objdel.cli.train.sft <same> --vision-reference-only   # R0: G-vision base ref

The flags are an objdel.core.spec.SftStage (cli_stage: one kebab flag per field, lists and dicts as JSON, a bare flag
for True).
"""

import sys

from objdel.policy.errors import exits_on_refusal
from objdel.sft.train import train


@exits_on_refusal
def main(argv=None) -> int:
    """CLI: the SftStage of the flags (or of --spec/--stage) -> train."""
    from objdel.core.spec import SftStage, cli_stage

    return train(cli_stage(SftStage, argv))


if __name__ == "__main__":
    sys.exit(main())
