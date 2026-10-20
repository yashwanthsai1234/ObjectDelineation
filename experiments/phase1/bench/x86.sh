#!/usr/bin/env bash
# One-command x86-64 (linux/amd64) check of objdel on an Apple-silicon Mac (Phase 1, M9; plan R9).
#
#   experiments/phase1/bench/x86.sh [--quick] [--no-wheels]
#     --quick      skip the slow K1/K2 10k-case oracle gates and use a 2,000-shape fuzz run
#     --no-wheels  skip the sdist + cibuildwheel (manylinux_2_28 x86_64) release-wheel check (needs network)
#
# What it does (outputs in build/x86/, git-ignored; summary in build/x86/summary.txt):
#   1. pins: `uv pip freeze` of the Mac's .venv -> the container installs exactly the same versions (torch from the
#      CPU index); the image objdel-x86:deps is rebuilt only when the pins change (bench/Dockerfile);
#   2. source: a tarball of the tracked and untracked-but-not-ignored files (git ls-files -co --exclude-standard), so
#      data/ (RIGA licence: never leaves the Mac), .venv/ and build/ are never copied; the run asserts it;
#   3. Mac (native arm64): experiments.phase1.bench.xarch_hashes -> build/x86/hashes-mac.json;
#   4. container (--platform linux/amd64 --cpus 3 --memory 4g --network none): bench/x86_inside.sh builds the
#      extension from source and runs lint, the fast suite, the hash script, K1/K2, the ASan/UBSan fuzz and a local
#      wheel + auditwheel + abi3audit;
#   5. release wheel: sdist on the host, then cibuildwheel 4.3.1 in the manylinux_2_28 x86_64 image (same caps), its
#      own test command (pyproject [tool.cibuildwheel]) and abi3audit; auditwheel show and the hash script against
#      the installed wheel in the x86 container;
#   6. compare: Mac vs x86 (raw split scores accepted, xarch_hashes.ACCEPTED) and x86 run 1 vs run 2 (strict);
#      exit status 1 on any other difference or a failed step.
#
# Docker Desktop runs linux/amd64 under Rosetta 2 (or QEMU) on Apple silicon: correct results, not x86 timing (x86
# timing is a Phase 2 entry gate on the first Vast box).
#
# References:
#   Docker multi-platform: https://docs.docker.com/build/building/multi-platform/
#   docker run resource limits: https://docs.docker.com/engine/containers/resource_constraints/
#   cibuildwheel options (container-engine create_args, archs, test-command):
#     https://cibuildwheel.pypa.io/en/stable/options/
#   auditwheel: https://github.com/pypa/auditwheel ; manylinux policies (PEP 600): https://peps.python.org/pep-0600/
#   abi3audit: https://github.com/pypa/abi3audit ; stable ABI (PEP 384): https://peps.python.org/pep-0384/
#   git ls-files --exclude-standard: https://git-scm.com/docs/git-ls-files
set -uo pipefail
QUICK=0; WHEELS=1
for a in "$@"; do
    case $a in
        --quick) QUICK=1 ;;
        --no-wheels) WHEELS=0 ;;
        -h|--help) sed -n '2,8p' "$0"; exit 0 ;;
        *) echo "unknown option $a"; exit 2 ;;
    esac
done

ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
OUT=$ROOT/build/x86
IO=$OUT/io
IMAGE=objdel-x86:deps
T0=$SECONDS
FAIL=0
cd "$ROOT" || exit 2
rm -rf "$IO" "$OUT/src" "$OUT/dist" "$OUT/wheelhouse"
mkdir -p "$IO" "$OUT/ctx"
note() { echo "$*" | tee -a "$OUT/summary.txt"; }
: > "$OUT/summary.txt"
note "x86 check started $(date -u +%FT%TZ) on $(uname -m) $(uname -s)"

# 1. pins and image -------------------------------------------------------------------------------------------------
uv pip freeze --python .venv/bin/python | grep -v -E '^(-e |objdel)' > "$OUT/ctx/requirements.txt.new"
if ! cmp -s "$OUT/ctx/requirements.txt.new" "$OUT/ctx/requirements.txt" || ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    mv "$OUT/ctx/requirements.txt.new" "$OUT/ctx/requirements.txt"
    t=$SECONDS
    docker build --platform linux/amd64 -f experiments/phase1/bench/Dockerfile -t "$IMAGE" "$OUT/ctx" \
        > "$OUT/image-build.log" 2>&1 || { note "image build FAILED (see build/x86/image-build.log)"; exit 1; }
    note "image rebuilt in $((SECONDS - t)) s"
else
    rm "$OUT/ctx/requirements.txt.new"
    note "image up to date (pins unchanged)"
fi

# 2. source snapshot ------------------------------------------------------------------------------------------------
git ls-files -z -co --exclude-standard | while IFS= read -r -d '' f; do [ -f "$f" ] && printf '%s\0' "$f"; done \
    > "$OUT/files.lst"
if tr '\0' '\n' < "$OUT/files.lst" | grep -q -E '^(data|\.venv|build)/'; then
    note "refusing: the snapshot would contain data/, .venv/ or build/"; exit 2
fi
# macOS tar would add AppleDouble ._* entries for extended attributes (com.apple.provenance): disable them
COPYFILE_DISABLE=1 tar --no-mac-metadata --no-xattrs --null -T "$OUT/files.lst" -czf "$IO/src.tar.gz"
if tar -tzf "$IO/src.tar.gz" | grep -q -E '(^|/)\._'; then note "refusing: AppleDouble ._* files in the snapshot"; exit 2; fi
manifest() { tr '\0' '\n' < "$OUT/files.lst" | grep -E '^(src|cpp|tests|experiments/phase1/bench)/|^(CMakeLists.txt|pyproject.toml)$' \
    | while IFS= read -r f; do [ -f "$f" ] && shasum -a 256 "$f"; done | shasum -a 256 | cut -c1-16; }
SNAP=$(manifest)
note "snapshot: $(tr -cd '\0' < "$OUT/files.lst" | wc -c | tr -d ' ') files, $(du -h "$IO/src.tar.gz" | cut -f1), source manifest $SNAP"
if [ cpp/polycore/polycore.hpp -nt build/.polycore-stamp ] || [ src/objdel/_polycore.cpp -nt build/.polycore-stamp ]; then
    note "WARNING: the Mac's objdel._polycore is older than its sources (run make build); Mac hashes may be stale"
fi

# 3. Mac hashes -----------------------------------------------------------------------------------------------------
t=$SECONDS
if .venv/bin/python -I -m experiments.phase1.bench.xarch_hashes --out "$OUT/hashes-mac.json" > "$OUT/hashes-mac.log" 2>&1; then
    note "mac hashes: $((SECONDS - t)) s"
else
    note "mac hashes FAILED (build/x86/hashes-mac.log)"; FAIL=1
fi
[ "$(manifest)" = "$SNAP" ] || note "WARNING: sources changed while the Mac hashes ran; re-run for a clean comparison"

# 4. container run --------------------------------------------------------------------------------------------------
t=$SECONDS
docker run --rm --platform linux/amd64 --cpus 3 --memory 4g --memory-swap 4g --network none \
    -e SLOW=$((1 - QUICK)) -e FUZZ_N=$([ $QUICK = 1 ] && echo 2000 || echo 10000) -v "$IO:/io" "$IMAGE" \
    bash -c 'tar -xzf /io/src.tar.gz -C /work && exec bash /work/experiments/phase1/bench/x86_inside.sh' \
    2>&1 | tee "$OUT/container.log"
note "container run: $((SECONDS - t)) s"
while IFS=$'\t' read -r name rc secs; do
    note "  x86 $(printf '%-7s' "$name") $([ "$rc" = 0 ] && echo PASS || echo "FAIL (exit $rc)")  $secs"
    [ "$rc" = 0 ] || [ "${name%-info}" != "$name" ] || FAIL=1  # *-info steps are informational (Rosetta artefacts)
done < "$IO/status.tsv"
grep -E ' passed| failed' "$IO/fast.log" | tail -1 | sed 's/^/  fast suite: /' | tee -a "$OUT/summary.txt"
grep -E '^FAILED ' "$IO/fast.log" | sed 's/^/    /' | tee -a "$OUT/summary.txt"
grep -E '^[^ ].*\.py:[0-9]+:[0-9]+|^ *--> ' "$IO/lint.log" | sed 's/^/    lint: /' | tee -a "$OUT/summary.txt"

# 5. release wheel --------------------------------------------------------------------------------------------------
if [ $WHEELS = 1 ]; then
    t=$SECONDS
    mkdir -p "$OUT/src" && tar -xzf "$IO/src.tar.gz" -C "$OUT/src"
    # pyproject's test-command, minus the out-of-memory test that Rosetta's translator breaks (x86_inside.sh)
    TEST_COMMAND="$(sed -n 's/^test-command = "\(.*\)"$/\1/p' "$OUT/src/pyproject.toml") \
--deselect tests/test_polycore.py::test_score_batch_out_of_memory_raises_in_python"
    if uv build --sdist --out-dir "$OUT/dist" "$OUT/src" > "$OUT/sdist.log" 2>&1 \
        && CIBW_ARCHS=x86_64 CIBW_CONTAINER_ENGINE="docker; create_args: --cpus 3 --memory 4g" \
           CIBW_TEST_COMMAND="$TEST_COMMAND" \
           uvx cibuildwheel==4.3.1 --platform linux --output-dir "$OUT/wheelhouse" "$OUT"/dist/objdel-*.tar.gz \
           > "$OUT/cibuildwheel.log" 2>&1; then
        mkdir -p "$IO/wheelhouse" && cp "$OUT"/wheelhouse/*.whl "$IO/wheelhouse/"
        note "cibuildwheel x86_64: $(basename "$OUT"/wheelhouse/*.whl) in $((SECONDS - t)) s"
        grep -E 'passed|failed' "$OUT/cibuildwheel.log" | tail -1 | sed 's/^/  cibuildwheel test-command: /' | tee -a "$OUT/summary.txt"
        docker run --rm --platform linux/amd64 --cpus 3 --memory 4g --network none -v "$IO:/io" "$IMAGE" bash -c '
            set -e; tar -xzf /io/src.tar.gz -C /work
            /opt/tools/bin/auditwheel show /io/wheelhouse/*.whl
            /opt/tools/bin/abi3audit --strict /io/wheelhouse/*.whl
            uv pip install --python /opt/venv/bin/python --offline --no-deps /io/wheelhouse/*.whl
            cd /work && python -c "import objdel._polycore as m, sys; print(m.__file__)"
            python -m experiments.phase1.bench.xarch_hashes --out /io/hashes-x86-wheel.json > /dev/null' \
            > "$OUT/wheel-check.log" 2>&1 && note "  release wheel: auditwheel show + abi3audit + hashes OK" \
            || { note "  release wheel check FAILED (build/x86/wheel-check.log)"; FAIL=1; }
        grep -i -E 'consistent with|manylinux' "$OUT/wheel-check.log" | head -2 | sed 's/^/  /' | tee -a "$OUT/summary.txt"
    else
        note "sdist / cibuildwheel FAILED (build/x86/sdist.log, cibuildwheel.log)"; FAIL=1
    fi
fi

# 6. compare --------------------------------------------------------------------------------------------------------
# a: reference, b: other; the IPP-off comparison is a diagnostic (impl/M9_x86.md) and does not set the exit status
for pair in "hashes-mac hashes-x86" "hashes-x86 hashes-x86-run2" "hashes-mac hashes-x86-wheel" "hashes-mac hashes-x86-noipp"; do
    set -- $pair
    a=$OUT/$1.json; [ -f "$a" ] || a=$IO/$1.json
    [ -f "$a" ] && [ -f "$IO/$2.json" ] || continue
    note "compare $1 vs $2:"
    strict=$([ "$2" = hashes-x86-run2 ] && echo --strict)
    .venv/bin/python -I -m experiments.phase1.bench.xarch_hashes --compare "$a" "$IO/$2.json" $strict \
        > "$OUT/compare-$2.txt"; rc=$?
    grep -v '^OK' "$OUT/compare-$2.txt" | sed 's/^/  /' | tee -a "$OUT/summary.txt"
    [ $rc = 0 ] || [ "$2" = hashes-x86-noipp ] || FAIL=1
done
note "total $((SECONDS - T0)) s; $([ $FAIL = 0 ] && echo 'ALL PASS' || echo 'SOME CHECKS FAILED') (details: build/x86/)"
exit $FAIL
