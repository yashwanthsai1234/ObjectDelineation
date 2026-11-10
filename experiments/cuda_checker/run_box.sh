#!/usr/bin/env bash
# The CUDA checker's GPU run, one command on a CUDA box with the Phase 2 image (repo rsynced, data/checker_eval only):
#   bash experiments/cuda_checker/run_box.sh [--build-only]
# 1. refuse to run if any RIGA data is on the host (its licence forbids rented hosts);
# 2. add the CUDA 13.0 compiler and NVTX headers (the image has none; its torch is 2.13.0+cu130,
#    experiments/phase2/docker/locks/vllm_native.txt), compute-sanitizer 13.4.92 (the toolkit extra's 13.0.85 wheel
#    lacks its TreeLauncher executables: K1, 2026-11-06) and pytest/hypothesis to $PY's venv;
# 3. make cuda: objdel._gpucheck rebuilt with nvcc (CUDA_ARCHS, default the box's GPU);
# 4. the GPU test group (with the arena-reuse and two-stream tests), the equivalence suite, the packed-mask and
#    one-checker-call-per-reward tests, then experiments/cuda_checker/eval.py even when a test failed (it shows where
#    the GPU build goes wrong); writes results/cuda.json and the GPU section of RESULTS.md: copy both back.
#    Exit status: the tests' (1 if any failed), else the eval's.
# --build-only stops after 3 (cuda/docker_check.sh runs it on the Mac, without a GPU).
# --gpu-check (make gpu-check, Implementation.md D8) replaces 4: pytest on $D_GPU_TESTS (determinism included; the step
#   fails on any skipped test), then the 4 compute-sanitizer tools on its sanitizer corpus. Outputs go to
#   $OUT = $OBJDEL_RUN_DIR/d (run.sh sets it; pull.sh brings it back): gpu.xml, sanitize.json.
#   Needs OBJDEL_DATA_ROOT (the LIDC bundle: checker_eval/masks.npz).
# --gpu-check --dry (Mac, make gpu-check DRY=1): prints every box step in order, then collects $D_GPU_TESTS; reads no
#   data/ and installs nothing (tests/experiments/test_gpu_chain.py runs it).
# Env: PY (default /opt/venv/bin/python, the image's training venv), CUDA_ARCHS (e.g. "80;86;89;90;120"),
#      DATA (the eval set, default data/checker_eval; R0 passes the bundle's copy).
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=${PY:-/opt/venv/bin/python}
# D's box test list: every box pytest call names it (tests/experiments/test_gpu_chain.py: each entry collects a test under $SEL).
D_GPU_TESTS="tests/geom/test_cuda_checker_gpu.py tests/geom/test_cuda_checker.py tests/verify/test_gpu_scoring.py \
tests/geom/test_gpu_determinism.py tests/geom/test_gpu_sanitize.py tests/experiments/test_rings.py tests/geom/test_certify_polygons.py \
tests/geom/test_topology_device.py tests/geom/test_tables_device.py tests/geom/test_engine_counters.py tests/experiments/test_scorer_sums.py"
SEL="gpu or not (slow or data)"  # the GPU tests and every fast test of the list; data/ tests are the Mac's (-I)
SAN_TOOLS="memcheck initcheck synccheck racecheck"  # memcheck without --leak-check full: the engine is never freed
DRY=; [ "${1:-}" = "--gpu-check" ] && [ "${2:-}" = "--dry" ] && DRY=1
step() { echo "+ $1"; [ -n "$DRY" ] || eval "$1"; }
if [ "${1:-}" = "--gpu-check" ]; then
    if [ -n "$DRY" ]; then OUT=$(mktemp -d); else
        OUT="${OBJDEL_RUN_DIR:?run under run.sh}/d"
        : "${OBJDEL_DATA_ROOT:?set by the run config}"
        mkdir -p "$OUT"
    fi
fi
if [ -z "$DRY" ] && find data -maxdepth 3 -iname '*riga*' 2>/dev/null | grep -q .; then
    echo "refused: RIGA data on this host (find data -iname '*riga*'); remove it first" >&2
    exit 2
fi
CONSTRAINTS=()
[ -f /opt/requirements.txt ] && CONSTRAINTS=(-c /opt/requirements.txt)  # the image's lock: nothing else may move
step 'uv pip install --python "$PY" "${CONSTRAINTS[@]}" "cuda-toolkit[nvcc,cccl,nvtx]==13.0.3" "nvidia-cuda-sanitizer-api==13.4.92" pytest hypothesis'
step 'CUDA_HOME=$("$PY" -c "import nvidia, os; print(os.path.join(nvidia.__path__[0], '"'cu13'"'))")'
# PyPI layout vs what nvcc and CMake's FindCUDAToolkit expect: libraries in lib64, an unversioned libcudart.so
step '[ -e "$CUDA_HOME/lib64" ] || ln -s lib "$CUDA_HOME/lib64"'
step 'ln -sf libcudart.so.13 "$CUDA_HOME/lib/libcudart.so"'
step 'export CUDACXX=$CUDA_HOME/bin/nvcc'
step 'make cuda PY="$PY" ${CUDA_ARCHS:+CUDA_ARCHS="$CUDA_ARCHS"}'
[ "${1:-}" = "--build-only" ] && exit 0
if [ "${1:-}" = "--gpu-check" ]; then
    rc=0
    step '"$PY" -I -c "import objdel._gpucheck as g; assert g.build() == '"'cuda'"'"'
    step '[ -f "$OBJDEL_DATA_ROOT/checker_eval/masks.npz" ] || { echo "box input missing: $OBJDEL_DATA_ROOT/checker_eval/masks.npz" >&2; exit 1; }'
    step '"$PY" -I -m pytest -q -p no:cacheprovider $D_GPU_TESTS -m "$SEL" --junitxml "$OUT/gpu.xml" || rc=1'
    step '"$PY" -I -m experiments.cuda_checker.junit_skips "$OUT/gpu.xml" || rc=1'
    for t in $SAN_TOOLS; do
        step 's=$SECONDS; e=0; PYTORCH_NO_CUDA_MEMORY_CACHING=1 timeout 300 "$CUDA_HOME/bin/compute-sanitizer" --tool '"$t"' --error-exitcode 99 --print-limit 20 "$PY" -I -m pytest -q -p no:cacheprovider $D_GPU_TESTS -m "gpu and sanitizer" || e=$?; echo "{\"tool\": \"'"$t"'\", \"exit\": $e, \"seconds\": $((SECONDS - s))}" >> "$OUT/sanitize.jsonl"; [ $e = 0 ] || rc=1'
    done
    if [ -n "$DRY" ]; then
        "$PY" -I -m pytest --collect-only -q -p no:cacheprovider $D_GPU_TESTS -m "$SEL"
        exit 0
    fi
    exit $rc
fi
nvidia-smi -L
tests=0
"$PY" -m pytest -q -p no:cacheprovider tests/geom/test_cuda_checker_gpu.py tests/geom/test_cuda_checker.py \
    tests/geom/test_checker_host.py tests/rl/test_reward_one_checker_call.py || tests=1
"$PY" -I -m experiments.cuda_checker.eval ${DATA:+--data "$DATA"}
exit $tests
