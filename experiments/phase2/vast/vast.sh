#!/usr/bin/env bash
# PLAN §4.5: the one shim. offers.sh, launch.sh, sync_up.sh, run.sh, pull.sh, teardown.sh and sweep_orphans.sh are
# symlinks to it; the link's name is the command of objdel.infra.vast (src/objdel/infra/vast/__main__.py), where
# all logic and usage live.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
exec "${PYTHON:-$ROOT/.venv/bin/python}" -I -m objdel.infra.vast "$(basename "$0" .sh)" "$@"
