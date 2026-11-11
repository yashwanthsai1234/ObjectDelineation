#!/usr/bin/env bash
# nvcc build check of objdel._gpucheck without a GPU, on the Mac: runs the build half of the box command
# (experiments/cuda_checker/run_box.sh --build-only) in a throwaway linux/arm64 container that mimics the Phase 2 image
# (experiments/phase2/docker/Dockerfile.phase2): same python:3.12-slim-trixie base (pinned digest), build-essential,
# cmake, ninja, uv and a venv at /opt/venv. run_box.sh adds the CUDA 13.0 compiler from PyPI and runs `make cuda` for
# sm_80 (A100), sm_86 (RTX 30xx, A10), sm_89 (RTX 40xx, L4), sm_90 (H100) and sm_120 (RTX 50xx). Then: the module
# imports as the cuda build, sees no device, and both calls fail with a clean CUDA error instead of crashing. Last, the
# engine source plus one function with two GC_RANGEs in one scope compiles for sm_86 and sm_120 (the NVTX macro names
# each range apart), and ptxas reports no register spill in any kernel. Its report (registers and spills per kernel and
# arch) is kept as build/ptxas.txt, with one summary line printed.
# Native arm64, not amd64 emulation: the PyPI CUDA wheels (nvcc, cccl, nvtx, sanitizer) have aarch64 builds, the base
# digest is multi-arch, and the device code is the same (docs/implementation2.0/notes/research/cuda-kernels.md
# section 3). Only the sources are copied in (no data/); ~1.5 GB of downloads, nothing kept.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BASE=python:3.12-slim-trixie@sha256:a6e34c598f2467ed0e9a8d349809fcd8b5c603269512df273a0bb1784edc11b1
mkdir -p "$ROOT/build"
docker run --rm --platform linux/arm64 -v "$ROOT":/src:ro -v "$ROOT/build":/out "$BASE" bash -euo pipefail -c '
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq && apt-get install -y -qq --no-install-recommends build-essential cmake ninja-build >/dev/null
pip install -q --root-user-action=ignore uv
mkdir /work && cd /src
cp -r --parents CMakeLists.txt Makefile pyproject.toml README.md cuda src experiments/cuda_checker/run_box.sh /work/
cd /work && uv venv -q /opt/venv && uv pip install -q --python /opt/venv/bin/python numpy
CUDA_ARCHS="80;86;89;90;120" bash experiments/cuda_checker/run_box.sh --build-only
/opt/venv/bin/python -I -c "
import numpy as np
from objdel import _gpucheck as g
sq = np.array([[0, 0], [9, 0], [9, 9], [0, 9]], np.int32)
calls = {\"check_masks\": lambda: g.check_masks(np.ones(4, np.uint8), np.array([[2, 2]], np.int32),
                                               np.array([0, 4], np.int64)),
         \"check_polygons\": lambda: g.check_polygons(sq, np.array([0, 4], np.int64), np.array([0, 1], np.int64),
                                                     None, np.array([5], np.int32), np.array([5], np.int32), 10)}
for name, call in calls.items():
    try:
        call()
        raise SystemExit(f\"{name}: a kernel ran without a GPU?\")
    except RuntimeError as e:
        print(f\"no GPU, {name} fails cleanly:\", e)
"
CU=$(/opt/venv/bin/python -c "import nvidia, os; print(os.path.join(nvidia.__path__[0], \"cu13\"))")
{ cat cuda/gpucheck.cpp; printf "void gc_two_ranges() {\n  GC_RANGE(\"a\");\n  GC_RANGE(\"b\");\n}\n"; } > /tmp/two.cu
"$CU/bin/nvcc" -std=c++17 -O3 --fmad=false -gencode arch=compute_86,code=sm_86 -gencode arch=compute_120,code=sm_120 \
    -Xptxas -v -Icuda -c /tmp/two.cu -o /tmp/two.o 2> /out/ptxas.txt || { cat /out/ptxas.txt; exit 1; }
echo "two GC_RANGEs in one scope: compiled (sm_86, sm_120)"
if grep -E "[1-9][0-9]* bytes spill" /out/ptxas.txt; then echo "register spills (above)"; exit 1; fi
awk -f /dev/stdin /out/ptxas.txt <<"AWK"
/Compiling entry function/ { arch = $NF; gsub(/[^a-z0-9_]/, "", arch) }
/Used [0-9]+ registers/ { n[arch]++; if ($5 + 0 > top[arch]) top[arch] = $5 + 0 }
END { for (a in n) s = s sprintf(" %s %d kernels, at most %d registers;", a, n[a], top[a])
      print "ptxas:" s " 0 bytes spilled (report: build/ptxas.txt)" }
AWK
echo "docker_check: OK"
'
