"""Final evaluation of a tuned model on a whole split: K samples per item, scored by the checker.

    python -I -m experiments.phase2.final_eval --ckpt DIR --view data/eval_views/lidc-val-shifted --split val \
        --out runs/eval/val [-k 8] [--seed 11] [--adapter LORA_DIR]
    python -I -m experiments.phase2.final_eval --ckpt DIR --view V1 --split val --out O1 \
        --view V2 --split val --out O2          # several views, one model load (one --out: O/<view name>/)
    python -I -m experiments.phase2.final_eval --ckpt DIR --view V --split val --out O \
        --limit 32                              # a GPU smoke: 32 evenly spaced items, not a result
    python -I -m experiments.phase2.final_eval --rescore runs/a3/eval_val --view V --split val --out O
    python -I -m experiments.phase2.final_eval --hash-view data/eval_views/lidc-test-shifted \
        --out runs/a3/pre_rental.json
    python -I -m experiments.phase2.final_eval --resummarize runs/r5/eval/eval_test --out runs/r5/eval_a1/eval_test \
        [--store data/store]
    python -I -m experiments.phase2.final_eval --spec specs/b1_rows.yaml --stage eval [--rd RD] [--data DATA]

The model forms read their flags as an objdel.core.spec.EvalStage (cli_stage: kebab flags, K as -k, or a spec stage)
and call objdel.policy.final_eval.run; --hash-view and --resummarize load no model and are read before it. What is
sampled, scored and written: objdel.policy.final_eval. Reads untrusted data: run with python -I.
"""

import argparse
import sys
from pathlib import Path

from experiments.phase2.freeze_views import render_files
from objdel.policy import final_eval as E


def no_model(argv: list[str]) -> None:
    """--hash-view VIEW --out FILE and --resummarize DIR --out DIR [--store S]: no model, no EvalStage."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--hash-view", type=Path, metavar="VIEW", help="write VIEW's identity record to --out and exit")
    mode.add_argument("--resummarize", type=Path, metavar="DIR", help="a finished run: summarize its items.jsonl again")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--store", type=Path, default=Path("data/store"), help="the item store (group_id per item)")
    a = ap.parse_args(argv)
    if a.hash_view:
        return E.hash_view(a.hash_view, a.out, render_files())
    if a.out.resolve() == a.resummarize.resolve():
        ap.error("--out must not be the --resummarize folder")
    return E.resummarize(a.resummarize, a.out, a.store)


def main(argv=None) -> None:
    """The no-model forms, else the EvalStage of the flags to objdel.policy.final_eval.run."""
    from objdel.core.spec import EvalStage, cli_stage

    argv = sys.argv[1:] if argv is None else list(argv)
    if {"--hash-view", "--resummarize"} & set(argv):
        return no_model(argv)
    return E.run(cli_stage(EvalStage, argv))


if __name__ == "__main__":
    main()
