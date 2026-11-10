#!/usr/bin/env bash
# The one Nsight Systems mechanism of the project (plan A6 step 8; C0 calls it for each arm, on the box):
#
#   nsys_capture.sh RD NAME S E -- CMD...
#
# Installs the pinned Nsight Systems CLI once per box (dpkg-deb -x: no GUI dependencies), runs CMD under
# `nsys profile` with `--profile-window S,E` appended (train_rl: callbacks.NvtxSteps starts the cudaProfilerApi
# capture at the end of step S-1 and stops it at the end of step E; S=1 starts it at train begin, so `1 STEPS` traces
# the whole run), then writes RD/nsys/NAME.{nsys-rep,sqlite},
# NAME_tables.txt (the export schema differs by nsys version) and NAME_stats_<report>.csv for
# experiments.phase2.nsys_baseline (--stats RD/nsys/NAME_stats).
# NSYS_TIMEOUT (seconds) bounds the whole command; unset it is 0, GNU timeout's "no limit", so a whole-run capture is
# never cut short (the run's own stop rules end it). Exit 3: install or `nsys status -e` failed; exit 4: capture failed
# or timed out. The caller then runs the fallback
# (same arm with `--trace DIR --profile-window 16,17` and OBJDEL_NVTX unset: torch.profiler on step 16).
# Pins: the Debian 13 CUDA repo package, then the CLI-only devtools package (Ubuntu 24.04 build) as the second try,
# each checked by SHA-256 [measured: curl, 2026-11-08]. A manual retry may install the pinned package with apt after
# the repo's cuda-keyring_1.1-1_all.deb. No CUDA memory-usage tracing (NVIDIA: "may cause significant runtime overhead"),
# and --cuda-graph-trace stays at its default `graph` (one row per graph run; nsys_baseline reads them).
set -euo pipefail; RD=$1; NAME=$2; S=$3; E=$4; shift 5        # $5 is the literal "--"
mkdir -p "$RD/nsys"
DEB=nsight-systems-2026.5.2_2026.5.2.323-1_amd64.deb
URL=https://developer.download.nvidia.com/compute/cuda/repos/debian13/x86_64/$DEB
SHA=24954cc4d4cb5f051dc91c6a44101d7ccf15c802ab8ea8a32c09dd22542f9031
DEB2=NsightSystems-linux-cli-public-2026.6.1.218-39266636.deb   # second try: CLI-only devtools package
URL2=https://developer.download.nvidia.com/devtools/repos/ubuntu2404/amd64/$DEB2
SHA2=e71440978244e3883fc769da91a44a03451129ce21988d0db656a6bf7c1cb56b
H=${NSYS_HOME:-/opt/nsys}                          # a test points it at a tmp folder with a stub
get() { curl -fsSL -o "/tmp/$1" "$2" && echo "$3  /tmp/$1" | sha256sum -c - && dpkg-deb -x "/tmp/$1" "$H"; }
if [ ! -d "$H" ]; then
  get "$DEB" "$URL" "$SHA" || { rm -rf "$H"; get "$DEB2" "$URL2" "$SHA2"; } || exit 3
fi
NSYS=$(find "$H" -type f -name nsys -path '*target-linux-x64*' -print -quit)   # no pipe: pipefail
[ -n "$NSYS" ] || NSYS=$(find "$H" -type f -name nsys -print -quit)   # the CLI package layout is not checked
{ "$NSYS" --version > "$RD/nsys/version.txt" && "$NSYS" status -e >> "$RD/nsys/version.txt"; } || exit 3
timeout --signal=INT --kill-after=120 "${NSYS_TIMEOUT:-0}" env OBJDEL_NVTX=1 "$NSYS" profile -t cuda,nvtx,osrt \
  --capture-range=cudaProfilerApi --capture-range-end=stop -o "$RD/nsys/$NAME" "$@" --profile-window "$S,$E" \
  || exit 4
"$NSYS" export --type sqlite -o "$RD/nsys/$NAME.sqlite" "$RD/nsys/$NAME.nsys-rep"
"${PY:-/opt/venv/bin/python}" -I -c 'import sqlite3, sys
for (n,) in sqlite3.connect(sys.argv[1]).execute("select name from sqlite_master order by name"): print(n)' \
  "$RD/nsys/$NAME.sqlite" > "$RD/nsys/${NAME}_tables.txt"      # the export schema differs by nsys version
"$NSYS" stats --report nvtx_pushpop_trace,nvtx_gpu_proj_trace,cuda_gpu_trace,cuda_api_trace --format csv \
  --output "$RD/nsys/${NAME}_stats" "$RD/nsys/$NAME.nsys-rep"
pkill -f "$H" || true
