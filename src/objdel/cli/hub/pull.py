"""One revision of the private HF model repo onto the box (RL starts from an SFT run's selected branch, e.g. r1-qwen3vl2b).

    python -I -m objdel.cli.hub.pull --repo svenk029/objdel-sft --revision r1-qwen3vl2b --out DIR < token
    python -I -m objdel.cli.hub.pull --repo svenk029/objdel-sft --revision HEX40 --out /workspace/data/lidc \
        --allow "rows/lidc-train-b1-780/**" < token     # the B1 rows root lands at /workspace/data/lidc/rows/lidc-train-b1

The token arrives on stdin (fetch_hf.sh streams it from objdel.infra.secrets over ssh, as push_hf does: never on argv,
in the environment or on disk). The library is objdel.infra.hf.pull.
"""

import argparse
import sys

from objdel.infra import hf


def main(argv=None) -> int:
    """The hf_pull CLI (--repo, --revision, --out, --allow; the token on stdin)."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--repo", required=True)
    ap.add_argument("--revision", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--allow", help="only repo paths matching this pattern (kept under --out at their repo paths)")
    a = ap.parse_args(argv)
    token = sys.stdin.readline().strip() or sys.exit("hf_pull: no token on stdin")
    return hf.pull(a.repo, a.revision, a.out, a.allow, token)


if __name__ == "__main__":
    sys.exit(main())
