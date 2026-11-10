"""Preflight: everything that can be proven on the Mac, before any GPU is rented (objdel.infra.preflight).

    python -I -m experiments.phase2.preflight --config specs/b1_rows.yaml     (make preflight CONFIG=specs/...)
"""

import argparse
import sys
from pathlib import Path

from objdel.infra import preflight


def main(argv=None) -> int:
    """Parse --config and run every check on it (objdel.infra.preflight.run)."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", type=Path, required=True, help="a run spec (specs/<name>.yaml)")
    return preflight.run(ap.parse_args(argv).config)


if __name__ == "__main__":
    sys.exit(main())
