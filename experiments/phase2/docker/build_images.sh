#!/usr/bin/env bash
# PLAN §4.3: build the Phase 2 (and, at Phase 3 start, Phase 3) Vast images for linux/amd64 on the Mac.
#   build_images.sh [--check] [--push] phase2|phase3|sampler ...
#     --check  lint the Dockerfiles only (docker buildx build --check); no locks, no build
#     --push   push to $IMAGE_REPO (private GHCR, e.g. ghcr.io/<owner>/objdel-vast; `docker login ghcr.io` first)
# Locks in build/images/ctx: train.txt (uv pip freeze of .venv, as x86.sh), vllm.txt (vllm==0.32.0 for
# x86_64-manylinux_2_40 / CPython 3.12), p3.txt + p3.override (B5 Q1). Tag p2-<pinhash> / p3-<pinhash> (first 12 hex
# of the SHA-256 over the Dockerfile and its locks), so the image is rebuilt only when the pins change; record the
# pushed digest in .vast-allowlist ("image": "<repo>:p2-<hash>@sha256:..."). Builder capped at 2 CPUs / 4 GB (B5).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd); CTX=$ROOT/build/images/ctx
CHECK=0; PUSH=0; WHICH=()
for a in "$@"; do
    case $a in
        --check) CHECK=1 ;; --push) PUSH=1 ;;
        -h|--help) sed -n '2,10p' "$0"; exit 0 ;;
        phase2|phase3|sampler) WHICH+=("$a") ;;
        *) echo "usage: build_images.sh [--check] [--push] phase2|phase3|sampler ..." >&2; exit 2 ;;
    esac
done
[ ${#WHICH[@]} -gt 0 ] || { echo "usage: build_images.sh [--check] [--push] phase2|phase3|sampler ..." >&2; exit 2; }
mkdir -p "$CTX"
if [ "$CHECK" = 1 ]; then
    for w in "${WHICH[@]}"; do f=${w/sampler/phase2}; docker buildx build --check -f "$HERE/Dockerfile.$f" "$CTX"; done
    exit 0
fi
compile() { uv pip compile "$@" --python-platform x86_64-manylinux_2_40 --python-version 3.12 --no-header --quiet; }
uv pip freeze --python "$ROOT/.venv/bin/python" | grep -v -E '^(-e |objdel)' > "$CTX/train.txt"
echo "vllm==0.32.0" > "$CTX/vllm.in" && compile "$CTX/vllm.in" -o "$CTX/vllm.txt"
sed -e 's/^torch==.*/torch==2.13.0/' -e 's/^torchvision==.*/torchvision==0.28.0/' \
    -e 's/^transformers==.*/transformers==5.19.0/' -e '/^setuptools==/d' "$CTX/train.txt" > "$CTX/p3.in"
echo "vllm==0.32.0" >> "$CTX/p3.in"
printf 'numba==0.68.0\nnumpy==2.5.4\n' > "$CTX/p3.override"
compile "$CTX/p3.in" --override "$CTX/p3.override" -o "$CTX/p3.txt"
docker buildx inspect objdel-builder >/dev/null 2>&1 || docker buildx create --name objdel-builder --driver docker-container --bootstrap >/dev/null
docker update --cpus 2 --memory 4g --memory-swap 4g buildx_buildkit_objdel-builder0 >/dev/null
for w in "${WHICH[@]}"; do
    f=${w/sampler/phase2}; target=(); [ "$w" = sampler ] && target=(--target sampler)
    locks=(train.txt vllm.txt); [ "$f" = phase3 ] && locks=(p3.txt p3.override)
    hash=$(cd "$CTX" && cat "$HERE/Dockerfile.$f" "${locks[@]}" | shasum -a 256 | cut -c1-12)
    tag="${IMAGE_REPO:-objdel-vast}:${w/phase/p}-$hash"
    out=(--load); [ "$PUSH" = 1 ] && out=(--push)
    t=$SECONDS
    docker buildx build --builder objdel-builder --platform linux/amd64 "${out[@]}" -t "$tag" "${target[@]}" \
        -f "$HERE/Dockerfile.$f" "$CTX"
    echo "$tag built in $((SECONDS - t)) s"
done
