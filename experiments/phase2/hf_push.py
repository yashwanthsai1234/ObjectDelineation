"""Upload a run's selected branch to the owner's private HF model repo (on the box; the Mac runs push_hf.sh).

    python -I -m experiments.phase2.hf_push --run RUN_DIR --repo svenk029/objdel-sft --revision r1-qwen3vl2b < token
    python -I -m experiments.phase2.hf_push --rl TRAIN_RL_OUT --repo svenk029/objdel-sft --revision b1-rows < token
    python -I -m experiments.phase2.hf_push --card-only --headline docs/results/headline.json --repo ... --revision ...
    (on the Mac, after pull.sh: add --view data/views/<name>, the local copy of the box's view; its view_hash must
    be the run's)

The token arrives on stdin (never on argv, in the environment or on disk). The library is objdel.infra.hf (push,
card_only); its docstring lists what it uploads and what it refuses.
"""

import argparse
import sys
from pathlib import Path

from objdel.infra import hf


def main(argv=None) -> int:
    """The hf_push CLI (--run, --rl or --card-only, --repo, --revision, --view, --headline; the token on stdin)."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--card-only", action="store_true", help="replace the card's Evaluation and Limitations only")
    ap.add_argument("--headline", type=Path, help="with --card-only: docs/results/headline.json")
    src.add_argument("--run", type=Path, help="train_sft output dir (select/ inside)")
    src.add_argument("--rl", type=Path, help="train_rl output dir (final/ and final_adapter/ inside)")
    ap.add_argument("--repo", required=True)
    ap.add_argument("--revision", required=True, help="the repo branch, e.g. r1-qwen3vl2b")
    ap.add_argument("--view", type=Path, help="local copy of the training view (default: the manifest's path)")
    a = ap.parse_args(argv)
    token = sys.stdin.readline().strip() or sys.exit("hf_push: no token on stdin")
    if a.card_only:
        return hf.card_only(a.repo, a.revision, a.headline, token) if a.headline else sys.exit(
            "hf_push: --card-only needs --headline")
    return hf.push(a.run, a.rl, a.repo, a.revision, a.view, token)


if __name__ == "__main__":
    sys.exit(main())
