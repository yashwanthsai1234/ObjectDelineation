#!/usr/bin/env bash
# The CUDA checker's GPU run, one command on a CUDA box with the Phase 2 image (repo rsynced, data/checker_eval only):
#   bash experiments/cuda_checker/run_box.sh [--build-only]
# 1. refuse to run if any RIGA data is on the host (its licence forbids rented hosts);
# 2. add the CUDA 13.0 compiler torch 2.14.1 pins (the image has none) and pytest/hypothesis to $PY's venv;
# 3. make cuda: objdel._gpucheck rebuilt with nvcc (CUDA_ARCHS, default the box's GPU);
# 4. the GPU test group and the equivalence suite, then experiments/cuda_checker/eval.py even when a test failed (it
#    shows where the GPU build goes wrong); writes results/cuda.json and the GPU section of RESULTS.md: copy both back.
#    Exit status: the tests' (1 if any failed), else the eval's.
# --build-only stops after 3 (cuda/docker_check.sh runs it on the Mac, without a GPU).
# Env: PY (default /opt/venv/bin/python, the image's training venv), CUDA_ARCHS (e.g. "80;86;89;90;120"),
#      DATA (the eval set, default data/checker_eval; R0 passes the bundle's copy).
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=${PY:-/opt/venv/bin/python}
if find data -maxdepth 3 -iname '*riga*' 2>/dev/null | grep -q .; then
    echo "refused: RIGA data on this host (find data -iname '*riga*'); remove it first" >&2
    exit 2
fi
CONSTRAINTS=()
[ -f /opt/requirements.txt ] && CONSTRAINTS=(-c /opt/requirements.txt)  # the image's lock: nothing else may move
uv pip install --python "$PY" "${CONSTRAINTS[@]}" "cuda-toolkit[nvcc,cccl]==13.0.3" pytest hypothesis
CUDA_HOME=$("$PY" -c "import nvidia, os; print(os.path.join(nvidia.__path__[0], 'cu13'))")
# PyPI layout vs what nvcc and CMake's FindCUDAToolkit expect: libraries in lib64, an unversioned libcudart.so
[ -e "$CUDA_HOME/lib64" ] || ln -s lib "$CUDA_HOME/lib64"
ln -sf libcudart.so.13 "$CUDA_HOME/lib/libcudart.so"
export CUDACXX=$CUDA_HOME/bin/nvcc
make cuda PY="$PY" ${CUDA_ARCHS:+CUDA_ARCHS="$CUDA_ARCHS"}
[ "${1:-}" = "--build-only" ] && exit 0
nvidia-smi -L
tests=0
"$PY" -m pytest -q -p no:cacheprovider tests/test_cuda_checker_gpu.py tests/test_cuda_checker.py || tests=1
"$PY" -I -m experiments.cuda_checker.eval ${DATA:+--data "$DATA"}
exit $tests
