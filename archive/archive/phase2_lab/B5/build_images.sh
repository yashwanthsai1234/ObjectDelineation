#!/usr/bin/env bash
# B5: build the Phase 2 (two venvs) and Phase 3 (one re-pinned venv) image candidates for linux/amd64 on the Mac,
# on a dedicated buildx builder capped at 2 CPUs / 4 GB, and record build time and image sizes.
#   build_images.sh <out_dir> [phase2|phase3|sampler ...]   (sampler = Dockerfile.phase2 --target sampler)
# <out_dir>/ctx must hold train.txt, vllm.txt, p3.txt, p3.override (see REPORT.md, "Locks").
# Images: b5-objdel-phase2:test, b5-objdel-phase3:test. Builder: b5-builder. Delete all three afterwards
# (cleanup at the end of REPORT.md).
set -uo pipefail
OUT=$1; shift
HERE=$(cd "$(dirname "$0")" && pwd)
CTX=$OUT/ctx
docker buildx inspect b5-builder >/dev/null 2>&1 || {
    docker buildx create --name b5-builder --driver docker-container >/dev/null
    docker buildx inspect --bootstrap b5-builder >/dev/null
}
docker update --cpus 2 --memory 4g --memory-swap 4g buildx_buildkit_b5-builder0 >/dev/null
for which in "$@"; do
    t=$SECONDS
    file=$which; target=()
    [ "$which" = sampler ] && { file=phase2; target=(--target sampler); }
    docker buildx build --builder b5-builder --platform linux/amd64 --load -t "b5-objdel-$which:test" \
        "${target[@]}" -f "$HERE/Dockerfile.$file" "$CTX" > "$OUT/build-$which.log" 2>&1
    rc=$?
    echo "$which build rc=$rc in $((SECONDS - t)) s" | tee -a "$OUT/build.txt"
    docker image inspect "b5-objdel-$which:test" --format "$which image size {{.Size}} bytes" 2>&1 | tee -a "$OUT/build.txt"
    docker buildx prune --builder b5-builder -f > /dev/null 2>&1  # disk budget: no build cache kept between builds
done
