#!/usr/bin/env bash
# PLAN §4.3: build the Vast images for linux/amd64: phase2 (two venvs), its sampler stage, and vllm_native
# (training + vLLM in one venv, D15; RL, owner 2026-11-06).
#   build_images.sh --lock      (Mac) rewrite locks/ from pyproject.toml and the Mac .venv; commit the result
#   build_images.sh [--check] [--push] [--prune] phase2|vllm_native|sampler ...
#                               (Mac, or GitHub Actions: .github/workflows/image.yml) build from the committed locks
#     --check  lint the Dockerfiles only (docker buildx build --check)
#     --push   push to $IMAGE_REPO (private GHCR, e.g. ghcr.io/<owner>/objdel-gpu; `docker login ghcr.io` first)
#     --prune  after the build, empty the objdel-builder cache (V0: ~24 GB of cache filled the Mac's disk)
# locks/ (committed): train.txt = pyproject's base deps + [train] + [build-system] requires (the box runs
# `pip install -e . --no-build-isolation`), for x86_64-manylinux_2_40 / CPython 3.12, constrained to the Mac .venv's
# versions so the box runs what the Mac tested (V0 audit: the old lock was the .venv freeze itself and carried 31 dev and
# download packages); vllm.txt = vllm==0.32.0 for the same platform; vllm_native.txt + vllm_native.override (B5 Q1). Tag
# p2-<pinhash> / sampler-<pinhash> / vllm_native-<pinhash> (first 12 hex of the SHA-256 over the Dockerfile and its inputs), so the image is rebuilt only when the
# pins change; record the pushed digest in .vast-allowlist ("image": "<repo>:p2-<hash>@sha256:..."). On the Mac the
# builder is capped at 2 CPUs / 4 GB (B5); a CI runner (CI=true) gives it the whole machine.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd); CTX=$ROOT/build/images/ctx; LOCKS=$HERE/locks
USAGE="usage: build_images.sh --lock | [--check] [--push] [--prune] phase2|vllm_native|sampler ..."
LOCK=0; CHECK=0; PUSH=0; PRUNE=0; WHICH=()
for a in "$@"; do
    case $a in
        --lock) LOCK=1 ;; --check) CHECK=1 ;; --push) PUSH=1 ;; --prune) PRUNE=1 ;;
        -h|--help) sed -n '2,16p' "$0"; exit 0 ;;
        phase2|vllm_native|sampler) WHICH+=("$a") ;;
        *) echo "$USAGE" >&2; exit 2 ;;
    esac
done
mkdir -p "$CTX"
if [ "$LOCK" = 1 ]; then
    cd "$ROOT"  # relative paths, so the locks' "# via" notes name no Mac directories
    compile() { uv pip compile "$@" --python-platform x86_64-manylinux_2_40 --python-version 3.12 --no-header --quiet; }
    mkdir -p "$LOCKS"; c=build/images/ctx; l=experiments/phase2/docker/locks
    uv pip freeze --python .venv/bin/python | grep -v -E '^(-e |objdel)' > "$c/mac.txt"
    .venv/bin/python -I -c "import tomllib; print(*tomllib.load(open('pyproject.toml', 'rb'))['build-system']['requires'], sep='\n')" > "$c/build.in"
    compile pyproject.toml --extra train "$c/build.in" -c "$c/mac.txt" -o "$l/train.txt"
    echo "vllm==0.32.0" > "$c/vllm.in" && compile "$c/vllm.in" -o "$l/vllm.txt"
    sed -e 's/^torch==.*/torch==2.13.0/' -e 's/^torchvision==.*/torchvision==0.28.0/' \
        -e 's/^transformers==.*/transformers==5.19.0/' -e '/^setuptools==/d' "$c/mac.txt" > "$c/vllm_native.in"
    echo "vllm==0.32.0" >> "$c/vllm_native.in"
    printf 'numba==0.68.0\nnumpy==2.5.4\n' > "$l/vllm_native.override"
    compile "$c/vllm_native.in" --override "$l/vllm_native.override" -o "$l/vllm_native.txt"
    echo "locks written: $(cd "$l" && echo *)"
    exit 0
fi
[ ${#WHICH[@]} -gt 0 ] || { echo "$USAGE" >&2; exit 2; }
for f in train.txt vllm.txt; do [ -f "$LOCKS/$f" ] || { echo "no $LOCKS/$f: run build_images.sh --lock on the Mac" >&2; exit 2; }; done
cp "$LOCKS"/* "$HERE/exclude-vllm.txt" "$CTX/"
if [ "$CHECK" = 1 ]; then
    for w in "${WHICH[@]}"; do f=${w/sampler/phase2}; docker buildx build --check -f "$HERE/Dockerfile.$f" "$CTX"; done
    exit 0
fi
docker buildx inspect objdel-builder >/dev/null 2>&1 || docker buildx create --name objdel-builder --driver docker-container --bootstrap >/dev/null
[ -n "${CI:-}" ] || docker update --cpus 2 --memory 4g --memory-swap 4g buildx_buildkit_objdel-builder0 >/dev/null
for w in "${WHICH[@]}"; do
    f=${w/sampler/phase2}; target=(); [ "$w" = sampler ] && target=(--target sampler)
    inputs=(train.txt vllm.txt exclude-vllm.txt); [ "$f" = vllm_native ] && inputs=(vllm_native.txt vllm_native.override)
    hash=$(cd "$CTX" && cat "$HERE/Dockerfile.$f" "${inputs[@]}" | shasum -a 256 | cut -c1-12)
    tag="${IMAGE_REPO:-objdel-gpu}:${w/phase/p}-$hash"
    out=(--load); [ "$PUSH" = 1 ] && out=(--push)
    t=$SECONDS
    docker buildx build --builder objdel-builder --platform linux/amd64 "${out[@]}" -t "$tag" "${target[@]}" \
        --label "org.opencontainers.image.source=https://github.com/yashwanthsai1234/ObjectDelineation" \
        -f "$HERE/Dockerfile.$f" "$CTX"
    echo "$tag built in $((SECONDS - t)) s"
    [ -z "${GITHUB_OUTPUT:-}" ] || echo "tag=$tag" >> "$GITHUB_OUTPUT"
done
if [ "$PRUNE" = 1 ]; then docker buildx prune --builder objdel-builder --force; fi
