#!/usr/bin/env bash
# Runs INSIDE the linux/amd64 container started by x86.sh (do not run on the host). /io is build/x86/io on the host:
# it holds src.tar.gz (the repo without data/, .venv, build/) and receives every log and report.
#
# Steps, each timed and logged to /io/<step>.log, status to /io/status.tsv:
#   build    editable install from source with no build isolation (as `make setup` / CI): scikit-build-core 1.1.1 +
#            nanobind 3.1.0 drive CMake Release (-O3 -DNDEBUG), g++ 14, x86-64 baseline ISA (no -march)
#   lint     ruff check src tests experiments (= make lint)
#   fast     pytest -q -m "not slow and not data" (data/ is absent, so data tests would skip anyway)
#   hashes   experiments.phase1.bench.xarch_hashes -> /io/hashes-x86.json; hashes2: a second run (x86 run-to-run);
#            noipp: the same with OpenCV's Intel IPP HAL off (diagnostic, impl/M9_x86.md)
#   slow     the K1/K2 exhaustive polycore oracle gates (pytest -m slow tests/test_polycore_oracle.py), or reduced
#            with SLOW=0
#   fuzz     make fuzz (ASan + UBSan, FUZZ_N shapes)
#   wheel    pip wheel of the same source, then auditwheel show and abi3audit --strict
#
# References:
#   scikit-build-core editable installs / build dir: https://scikit-build-core.readthedocs.io/en/latest/configuration/index.html
#   nanobind CMake (STABLE_ABI, NOMINSIZE): https://nanobind.readthedocs.io/en/latest/api_cmake.html
#   uv pip --no-build-isolation: https://docs.astral.sh/uv/reference/cli/#uv-pip-install
#   auditwheel show: https://github.com/pypa/auditwheel#usage ; abi3audit: https://github.com/pypa/abi3audit
set -uo pipefail
IO=/io
cd /work
: > "$IO/status.tsv"
[ "$(uname -m)" = x86_64 ] || { echo "not x86_64: $(uname -m)"; exit 2; }

step() {  # step <name> <command...>: run, log, time, record the exit status; never aborts the script
    local name=$1; shift
    local t0=$SECONDS
    echo "== $name: $*"
    "$@" > "$IO/$name.log" 2>&1
    local rc=$?
    printf '%s\t%s\t%ss\n' "$name" "$rc" "$((SECONDS - t0))" | tee -a "$IO/status.tsv"
    tail -3 "$IO/$name.log"
}

{
    echo "uname: $(uname -a)"
    grep -m1 'model name' /proc/cpuinfo
    echo "flags: $(grep -m1 flags /proc/cpuinfo | tr ' ' '\n' | grep -E '^(sse4_2|avx|avx2|fma|avx512f|avx512bw)$' | tr '\n' ' ')"
    echo "nproc: $(nproc)  cgroup cpu.max: $(cat /sys/fs/cgroup/cpu.max 2>/dev/null)  memory.max: $(cat /sys/fs/cgroup/memory.max 2>/dev/null)"
    python -V; g++ --version | head -1; cmake --version | head -1
    uv pip freeze --python /opt/venv/bin/python
} > "$IO/env.txt" 2>&1

step build uv pip install --python /opt/venv/bin/python --offline --no-deps --no-build-isolation -e .
# the flags the extension was actually compiled with (CMake Release + nanobind)
grep -h -E '^ *(FLAGS|DEFINES) *=' build/*/build.ninja 2>/dev/null | sort -u > "$IO/compile_flags.txt" \
    || grep -rh -E 'CXX_FLAGS|CXX_DEFINES' build/*/CMakeFiles/_polycore.dir/flags.make > "$IO/compile_flags.txt" 2>/dev/null
grep -h -E '^CMAKE_(BUILD_TYPE|CXX_COMPILER|CXX_FLAGS_RELEASE|CXX_FLAGS):' build/*/CMakeCache.txt >> "$IO/compile_flags.txt"

# Under Rosetta (cpuinfo "VirtualApple") the RLIMIT_AS out-of-memory test kills the translator itself ("rosetta error:
# mmap_anonymous_rw mmap failed") at random: it is deselected from the gate and run separately as information.
OOM=tests/test_polycore.py::test_score_batch_out_of_memory_raises_in_python
DESELECT=()
if grep -q VirtualApple /proc/cpuinfo; then DESELECT=(--deselect "$OOM"); fi

step lint python -m ruff check src tests experiments
step fast python -m pytest -q -m "not slow and not data" -p no:cacheprovider --durations=15 \
    --junitxml="$IO/junit-fast.xml" "${DESELECT[@]}"
[ ${#DESELECT[@]} = 0 ] || step oom-info python -m pytest -q -p no:cacheprovider "$OOM"
step hashes python -I -m experiments.phase1.bench.xarch_hashes --out "$IO/hashes-x86.json"
# run-to-run determinism on x86 (the IPP distance transform used before F1 was not: OpenCV issue 24082), and the
# IPP-off diagnostic
step hashes2 python -I -m experiments.phase1.bench.xarch_hashes --out "$IO/hashes-x86-run2.json"
step noipp python -I -m experiments.phase1.bench.xarch_hashes --opencv no-ipp --out "$IO/hashes-x86-noipp.json"
if [ "${SLOW:-1}" = 1 ]; then
    step slow python -m pytest -q -s -m slow -p no:cacheprovider tests/test_polycore_oracle.py --durations=0
fi
step fuzz make fuzz "FUZZ_N=${FUZZ_N:-10000}" CXX=g++
mkdir -p "$IO/wheel"
step wheel bash -c "rm -f $IO/wheel/*.whl && \
    uv build --wheel --offline --no-build-isolation --python /opt/venv/bin/python --out-dir $IO/wheel . \
    && /opt/tools/bin/auditwheel show $IO/wheel/*.whl && /opt/tools/bin/abi3audit --strict --report $IO/wheel/*.whl"
echo "== done"
