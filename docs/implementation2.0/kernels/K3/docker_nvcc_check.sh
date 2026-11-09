#!/usr/bin/env bash
# nvcc build check of the K3 kernels without a GPU (on the Mac), mirroring REPO/cuda/docker_check.sh: the Phase 2
# image's base (python:3.12-slim-trixie, pinned digest) + build-essential + the CUDA 13.0 compiler and CCCL from PyPI
# (cuda-toolkit[nvcc,cccl]==13.0.3, the version torch 2.14.1 pins; REPO/experiments/cuda_checker/run_box.sh), then
# `make all ARCH=sm_120` and the CPU tests inside the container. Several minutes under amd64 emulation; nothing kept.
#   bash docker_nvcc_check.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
BASE=python:3.12-slim-trixie@sha256:a6e34c598f2467ed0e9a8d349809fcd8b5c603269512df273a0bb1784edc11b1
docker run --rm --platform linux/amd64 -v "$HERE":/src:ro "$BASE" bash -euo pipefail -c '
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq && apt-get install -y -qq --no-install-recommends build-essential >/dev/null
pip install -q --root-user-action=ignore uv
uv venv -q /opt/venv
uv pip install -q --python /opt/venv/bin/python "cuda-toolkit[nvcc,cccl]==13.0.3"
CUDA_HOME=$(/opt/venv/bin/python -c "import nvidia, os; print(os.path.join(nvidia.__path__[0], \"cu13\"))")
[ -e "$CUDA_HOME/lib64" ] || ln -s lib "$CUDA_HOME/lib64"
ln -sf libcudart.so.13 "$CUDA_HOME/lib/libcudart.so"
export PATH=$CUDA_HOME/bin:$PATH
nvcc --version | tail -2
mkdir /work && cp /src/*.cu /src/*.cuh /src/*.h /src/*.cpp /src/Makefile /work/ && cd /work
make all ARCH=sm_120 NVFLAGS="-O3 -std=c++17 -arch=sm_120 --fmad=false -lineinfo -Xcompiler -ffp-contract=off -Xptxas -v -L$CUDA_HOME/lib" 2>&1 \
  | grep -E "error|warning|registers|spill|smem|Compiling entry" | sed "s/^ptxas info    : //" | head -80
./build/test_cpu
ls -la build/
'
