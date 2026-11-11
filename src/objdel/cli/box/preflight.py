"""Preflight: everything that can be proven on the Mac, before any GPU is rented (objdel.infra.preflight).

    python -I -m objdel.cli.box.preflight --config specs/b1_rows.yaml [--bundle-dir DIR]
    (make preflight CONFIG=specs/... [BUNDLE_DIR=DIR]; DIR as sync_up's --bundle-dir, default build/bundles)
"""

import argparse
import sys
from pathlib import Path

from objdel.infra import preflight


def main(argv=None) -> int:
    """Parse --config and --bundle-dir and run every check on them (objdel.infra.preflight.run)."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", type=Path, required=True, help="a run spec (specs/<name>.yaml)")
    ap.add_argument("--bundle-dir", type=Path, default=preflight.BUNDLE_DIR,
                    help="the folder of bundle-lidc.tar, as sync_up's --bundle-dir (default build/bundles)")
    a = ap.parse_args(argv)
    return preflight.run(a.config, a.bundle_dir)


if __name__ == "__main__":
    sys.exit(main())
