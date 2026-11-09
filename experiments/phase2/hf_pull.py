"""One revision of the private HF model repo onto the box (RL starts from R1's selected branch, e.g. r1-qwen3vl2b).

    python -I -m experiments.phase2.hf_pull --repo svenk029/objdel-sft --revision r1-qwen3vl2b --out DIR < token
    python -I -m experiments.phase2.hf_pull --repo svenk029/objdel-sft --revision HEX40 --out /workspace/data/lidc \
        --allow "rows/lidc-train-b1/**" < token     # the B1 rows root lands at /workspace/data/lidc/rows/lidc-train-b1

The token arrives on stdin (fetch_hf.sh streams it from objdel.secrets over ssh, as push_hf does: never on argv, in
the environment or on disk) and is passed to huggingface_hub directly, so nothing is cached on the box.
"""

import argparse
import sys

from huggingface_hub import snapshot_download


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--repo", required=True)
    ap.add_argument("--revision", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--allow", help="only repo paths matching this pattern (kept under --out at their repo paths)")
    a = ap.parse_args(argv)
    token = sys.stdin.readline().strip() or sys.exit("hf_pull: no token on stdin")
    print(snapshot_download(repo_id=a.repo, revision=a.revision, local_dir=a.out, token=token,
                            allow_patterns=a.allow))
    return 0


if __name__ == "__main__":
    sys.exit(main())
