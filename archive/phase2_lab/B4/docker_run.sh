#!/usr/bin/env bash
# B4: run the B4 lab inside the linux/amd64 objdel-x86:deps image (experiments/phase1/bench/Dockerfile; emulated on
# Apple silicon: relative speed-ups and correctness only, not x86 timing).
#
#   [STEPS="pool_bench --batches 10;trl_pattern"] experiments/phase2/lab/B4/docker_run.sh <out dir>
#
# The source goes in as a tarball of tracked + untracked-not-ignored files under src/ cpp/ experiments/ tests/oracles
# (never data/); data/prepared is bind-mounted read-only at run time only (RIGA governance: no RIGA data in any image
# or layer). Limits as the M9 bench: --cpus 3 --memory 4g --memory-swap 4g --network none. Container: b4-x86 (--rm).
# References: https://docs.docker.com/engine/containers/resource_constraints/ ,
# https://docs.docker.com/engine/storage/bind-mounts/ (readonly bind mounts)
set -euo pipefail
OUT=$(cd "$1" && pwd)
ROOT=$(cd "$(dirname "$0")/../../../.." && pwd)
cd "$ROOT"
mkdir -p "$OUT/io"
git ls-files -z -co --exclude-standard -- src cpp experiments tests/oracles tests/__init__.py CMakeLists.txt pyproject.toml README.md \
    | tr '\0' '\n' | grep -v -E '^data/' | while IFS= read -r f; do if [ -f "$f" ]; then printf "%s\n" "$f"; fi; done > "$OUT/io/files.lst"
COPYFILE_DISABLE=1 tar --no-mac-metadata --no-xattrs -T "$OUT/io/files.lst" -czf "$OUT/io/src.tar.gz"
if tar -tzf "$OUT/io/src.tar.gz" | grep -q -E '^data/'; then echo "refusing: data/ in the tarball"; exit 2; fi
docker run --rm --name b4-x86 --platform linux/amd64 --cpus 3 --memory 4g --memory-swap 4g --network none \
    -v "$OUT/io:/io" -v "$ROOT/data/prepared:/work/data/prepared:ro" -e B4_SCRATCH=/io -e STEPS="${STEPS:-k3syn --reps 20;pool_bench --batches 10;trl_pattern;preload --max-riga 200 --steps 6 --threads 2 --pool 4;k3_harness --batches 5}" \
    objdel-x86:deps bash -c '
set -uo pipefail
tar -xzf /io/src.tar.gz -C /work && cd /work
{ grep -m1 "model name" /proc/cpuinfo; nproc; cat /sys/fs/cgroup/cpu.max; cat /sys/kernel/mm/transparent_hugepage/enabled;
  python -c "import cv2; print([l.strip() for l in cv2.getBuildInformation().splitlines() if \"Parallel\" in l])"; } > /io/env.txt 2>&1
uv pip install --python /opt/venv/bin/python --offline --no-deps --no-build-isolation -e . > /io/build.log 2>&1 || exit 3
IFS=";" read -ra steps <<< "$STEPS"
for step in "${steps[@]}"; do
  name=${step%% *}
  echo "== $step $(date +%T)"
  if [ "$name" = fs_main_probe ]; then
    for h in 0 1; do HEAVY=$h python -I experiments/phase2/lab/B4/fs_main_probe.py > /io/fs_main_probe_$h.log 2>&1; done
  else
    python -I -m experiments.phase2.lab.B4.$step --out /io > /io/$name.log 2>&1
  fi
  echo "rc=$? $(date +%T)"
done'
