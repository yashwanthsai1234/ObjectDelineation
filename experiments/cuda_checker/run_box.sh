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
# --gpu-check (make gpu-check, Implementation.md D8) replaces 4. Step 0: old_engine.tar and its .sha256 in the repo root
#   (make old-engine-tar on the Mac, D-G7), else exit 1 before the toolkit install. Then pytest on $D_GPU_TESTS
#   (determinism included; the step fails on any skipped test), the 4 compute-sanitizer tools on its sanitizer corpus,
#   the pre-D2 engine built from the tar into /opt/old (step 5), and the old/current benches, the replay and their
#   D-G7/D-G8 table (step 6, bench_checker.py --compare). Outputs go to $OUT = $OBJDEL_RUN_DIR/d (run.sh sets it;
#   pull.sh brings it back): gpu.xml, sanitize.jsonl, bench_*.json, replay_*.json, replay_inputs.npz.
#   Needs OBJDEL_DATA_ROOT (the LIDC bundle: $R_D2_PAYLOAD).
# --gpu-check --dry (Mac, make gpu-check DRY=1): with the owner's tag pre-d2, checks the tar and builds the old tree on
#   the CPU in a temp folder (else prints "skipped"); prints every box step in order, collects $D_GPU_TESTS and runs
#   S08, S10 and S15 (box inputs ship, the replay arm needs no cv2/scipy, the bench reads only under its data root).
#   Reads no data/ (tests/experiments/test_gpu_chain.py runs it).
# Env: PY (default /opt/venv/bin/python, the image's training venv), CUDA_ARCHS (e.g. "80;86;89;90;120"),
#      DATA (the eval set, default data/checker_eval; R0 passes the bundle's copy).
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=${PY:-/opt/venv/bin/python}
# D's box test list: every box pytest call names it (tests/experiments/test_gpu_chain.py: each entry collects a test under $SEL).
D_GPU_TESTS="tests/geom/test_cuda_checker_gpu.py tests/geom/test_cuda_checker.py tests/verify/test_gpu_scoring.py \
tests/geom/test_gpu_determinism.py tests/geom/test_gpu_sanitize.py tests/experiments/test_rings.py tests/geom/test_certify_polygons.py \
tests/geom/test_topology_device.py tests/geom/test_tables_device.py tests/geom/test_engine_counters.py tests/experiments/test_scorer_sums.py"
SEL="gpu or not (slow or data or post_commit)"  # GPU + fast tests; data/ is the Mac's (-I), post_commit needs .git
SAN_TOOLS="memcheck initcheck synccheck racecheck"  # memcheck without --leak-check full: the engine is never freed
# R-D2 payload under $OBJDEL_DATA_ROOT (Implementation.md D8): step 3 checks it before pytest, so a missing input never
# shows up only as a skip
R_D2_PAYLOAD="checker_eval/masks.npz checker_eval/polygons.npz checker_eval/manifest.json prepared/lidc/val.jsonl \
prepared/caches/lidc"
# S08, S10, S15 (make gpu-check DRY=1 runs them; none reads data/)
DRY_TESTS="tests/infra/test_box_inputs.py tests/experiments/test_bench_checker.py::test_replay_arm_runs_without_cv2_or_scipy \
tests/experiments/test_bench_checker.py::test_inputs_only_under_the_data_root \
tests/experiments/test_bench_checker.py::test_old_engine_child_gets_the_same_data_root \
tests/experiments/test_bench_checker.py::test_no_data_path_literal_in_the_script"
CONSTRAINTS=()
[ -f /opt/requirements.txt ] && CONSTRAINTS=(-c /opt/requirements.txt)  # the image's lock: nothing else may move
DRY=; [ "${1:-}" = "--gpu-check" ] && [ "${2:-}" = "--dry" ] && DRY=1
step() { echo "+ $1"; [ -n "$DRY" ] || eval "$1"; }
need_tar() {  # step 0 (D-G7): the code sync ships the two repo-root files that make old-engine-tar writes on the Mac
    [ -f old_engine.tar ] && [ -f old_engine.tar.sha256 ] && return
    echo "old_engine.tar missing: run make old-engine-tar on the Mac, then sync_up.sh" >&2
    exit 1
}
dry_old_engine() {  # the box's step 5 in a temp folder on the Mac, CPU build (make build: no OBJDEL_CUDA)
    sha256sum -c old_engine.tar.sha256
    "$PY" -I -c "from objdel.infra.vast.remote import code_files; assert {'old_engine.tar', 'old_engine.tar.sha256'} <= set(code_files())"
    T=$(mktemp -d)
    mkdir -p "$T/old_src" && tar -xf old_engine.tar -C "$T/old_src"
    uv venv --python "$PY" "$T/old"
    uv pip install --python "$T/old/bin/python" numpy
    make -C "$T/old_src" build PY="$T/old/bin/python"
    "$T/old/bin/python" -I -c "import objdel._gpucheck as g; assert g.build() == 'cpu'"
}
if [ "${1:-}" = "--gpu-check" ]; then
    if [ -n "$DRY" ]; then
        OUT=$(mktemp -d)
        if git rev-parse -q --verify 'pre-d2^{commit}' >/dev/null 2>&1; then need_tar; dry_old_engine
        else echo "old-engine tar: skipped, no pre-d2 tag"; fi
    else
        OUT="${OBJDEL_RUN_DIR:?run under run.sh}/d"
        need_tar
        : "${OBJDEL_DATA_ROOT:?set by the run config}"
        mkdir -p "$OUT"
    fi
fi
if [ -z "$DRY" ] && find data -maxdepth 3 -iname '*riga*' 2>/dev/null | grep -q .; then
    echo "refused: RIGA data on this host (find data -iname '*riga*'); remove it first" >&2
    exit 2
fi
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
    step 'miss=; for f in $R_D2_PAYLOAD; do [ -e "$OBJDEL_DATA_ROOT/$f" ] || miss="$miss $OBJDEL_DATA_ROOT/$f"; done; [ -f tests/goldens/d_gpu02.json ] || miss="$miss tests/goldens/d_gpu02.json"; [ -z "$miss" ] || { echo "box inputs missing:$miss" >&2; exit 1; }'
    step '"$PY" -I -m pytest -q -p no:cacheprovider $D_GPU_TESTS -m "$SEL" --junitxml "$OUT/gpu.xml" || rc=1'
    step '"$PY" -I -m experiments.cuda_checker.junit_skips "$OUT/gpu.xml" || rc=1'
    for t in $SAN_TOOLS; do
        step 's=$SECONDS; e=0; PYTORCH_NO_CUDA_MEMORY_CACHING=1 timeout 300 "$CUDA_HOME/bin/compute-sanitizer" --tool '"$t"' --error-exitcode 99 --print-limit 20 "$PY" -I -m pytest -q -p no:cacheprovider $D_GPU_TESTS -m "gpu and sanitizer" || e=$?; echo "{\"tool\": \"'"$t"'\", \"exit\": $e, \"seconds\": $((SECONDS - s))}" >> "$OUT/sanitize.jsonl"; [ $e = 0 ] || rc=1'
    done
    # step 5 (D-G7): the pre-D2 engine beside the current one, from the root tar; same interpreter, numpy and CUDACXX
    step 'sha256sum -c old_engine.tar.sha256'
    step 'mkdir -p /opt/old_src && tar -xf old_engine.tar -C /opt/old_src'
    step 'uv venv --python "$PY" /opt/old'
    step 'uv pip install --python /opt/old/bin/python "${CONSTRAINTS[@]}" numpy'
    step 'make -C /opt/old_src cuda PY=/opt/old/bin/python ${CUDA_ARCHS:+CUDA_ARCHS="$CUDA_ARCHS"}'
    # step 6 (D-G7, D-G8): the replay inputs (main venv: cv2, scipy), each engine's bench and replay, the gate table
    step '"$PY" -I -m experiments.phase2.perf.bench_checker --dump-replay "$OUT/replay_inputs.npz"'
    for e in old current; do
        step '"$PY" -I -m experiments.phase2.perf.bench_checker --engine '"$e"' --out "$OUT/bench_'"$e"'.json"'
        step '"$PY" -I -m experiments.phase2.perf.bench_checker --engine '"$e"' --replay "$OUT/replay_inputs.npz" --out "$OUT/replay_'"$e"'.json"'
    done
    step '"$PY" -I -m experiments.phase2.perf.bench_checker --compare "$OUT" || rc=1'
    if [ -n "$DRY" ]; then
        "$PY" -I -m pytest --collect-only -q -p no:cacheprovider $D_GPU_TESTS -m "$SEL"
        "$PY" -I -m pytest -q -p no:cacheprovider $DRY_TESTS
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
