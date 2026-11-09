#!/usr/bin/env bash
# nvcc_docker_check.sh <REPO/cuda> -- compile and link the GPU binaries of K1 with nvcc (CUDA 13.0.3 from PyPI, the
# toolkit REPO/experiments/cuda_checker/run_box.sh uses) for sm_120 in a throwaway linux/amd64 container, without a
# GPU. Catches nvcc errors before a GPU is rented. Several minutes under emulation; ~1.5 GB downloaded, nothing kept.
set -euo pipefail
K1=$(cd "$(dirname "$0")" && pwd)
GC=$(cd "$1" && pwd)
BASE=python:3.12-slim-trixie@sha256:a6e34c598f2467ed0e9a8d349809fcd8b5c603269512df273a0bb1784edc11b1
docker run --rm --platform linux/amd64 -v "$K1":/k1src:ro -v "$GC":/gc:ro "$BASE" bash -euo pipefail -c '
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq && apt-get install -y -qq --no-install-recommends build-essential >/dev/null
pip install -q --root-user-action=ignore uv
uv venv -q /opt/venv
uv pip install -q --python /opt/venv/bin/python "cuda-toolkit[nvcc,cccl]==13.0.3"
CUDA_HOME=$(/opt/venv/bin/python -c "import nvidia, os; print(os.path.join(nvidia.__path__[0], \"cu13\"))")
[ -e "$CUDA_HOME/lib64" ] || ln -s lib "$CUDA_HOME/lib64"
ln -sf libcudart.so.13 "$CUDA_HOME/lib/libcudart.so"
mkdir -p /work && cp /k1src/*.cpp /k1src/*.h /k1src/build_gpu.sh /work/
"$CUDA_HOME/bin/nvcc" --version | tail -2
NVCC="$CUDA_HOME/bin/nvcc" bash /work/build_gpu.sh /gc sm_120
ls -la /work/build/gpu/
echo "nvcc_docker_check: OK"
'
