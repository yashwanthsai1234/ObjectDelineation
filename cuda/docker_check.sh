#!/usr/bin/env bash
# nvcc build check of objdel._gpucheck without a GPU, on the Mac: runs the build half of the box command
# (experiments/cuda_checker/run_box.sh --build-only) in a throwaway linux/amd64 container that mimics the Phase 2 image
# (experiments/phase2/docker/Dockerfile.phase2): same python:3.12-slim-trixie base (pinned digest), build-essential,
# cmake, ninja, uv and a venv at /opt/venv. run_box.sh adds the CUDA 13.0 compiler from PyPI and runs `make cuda` for
# sm_80 (A100), sm_86 (RTX 30xx, A10), sm_89 (RTX 40xx, L4), sm_90 (H100) and sm_120 (RTX 50xx). Then: the module
# imports as the cuda build, sees no device, and a call fails with a clean CUDA error instead of crashing.
# Only the sources are copied in (no data/). Several minutes under amd64 emulation; ~1.5 GB of downloads, nothing kept.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BASE=python:3.12-slim-trixie@sha256:a6e34c598f2467ed0e9a8d349809fcd8b5c603269512df273a0bb1784edc11b1
docker run --rm --platform linux/amd64 -v "$ROOT":/src:ro "$BASE" bash -euo pipefail -c '
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq && apt-get install -y -qq --no-install-recommends build-essential cmake ninja-build >/dev/null
pip install -q --root-user-action=ignore uv
mkdir /work && cd /src
cp -r --parents CMakeLists.txt Makefile pyproject.toml README.md cpp cuda src experiments/cuda_checker/run_box.sh /work/
cd /work && uv venv -q /opt/venv && uv pip install -q --python /opt/venv/bin/python numpy
CUDA_ARCHS="80;86;89;90;120" bash experiments/cuda_checker/run_box.sh --build-only
/opt/venv/bin/python -I -c "
import numpy as np
from objdel import _gpucheck as g
try:
    g.check_masks(np.ones(4, np.uint8), np.array([[2, 2]], np.int32), np.array([0, 4], np.int64))
    raise SystemExit(\"a kernel ran without a GPU?\")
except RuntimeError as e:
    print(\"no GPU, call fails cleanly:\", e)
"
echo "docker_check: OK"
'
