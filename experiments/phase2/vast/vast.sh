#!/usr/bin/env bash
# PLAN §4.5: the one shim. offers.sh, launch.sh, sync_up.sh, run.sh, pull.sh, teardown.sh and sweep_orphans.sh are
# symlinks to it; the link's name is the policy.py command, and all logic and usage live in policy.py.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
exec "${PYTHON:-$ROOT/.venv/bin/python}" -I -m experiments.phase2.vast.policy "$(basename "$0" .sh)" "$@"
