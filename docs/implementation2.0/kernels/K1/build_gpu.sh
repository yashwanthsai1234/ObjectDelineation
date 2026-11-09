#!/usr/bin/env bash
# build_gpu.sh <dir with gpucheck.h/.cpp (REPO/cuda)> [arch=sm_120]
# Builds into build/gpu/:
#   test_k1_gpu   K1 on the GPU vs ref_polygonize + gpucheck::check_masks (CPU build) -- exactness
#   k1_dump_gpu   K1 on the GPU on a raw mask file, compared in-process with gpucheck::check_masks (CPU build)
#   bench_k1_gpu  K1 on the GPU vs gpucheck::check_masks on the GPU -- speed (and equality)
# Env: NVCC (default nvcc), CXX (default g++).
set -euo pipefail
K1=$(cd "$(dirname "$0")" && pwd)
GC=$(cd "$1" && pwd)
ARCH=${2:-sm_120}
NVCC=${NVCC:-nvcc}
CXX=${CXX:-g++}
O=$K1/build/gpu
mkdir -p "$O"
NVFLAGS=(-std=c++17 -O3 "-arch=$ARCH" -lineinfo)
set -x
"$NVCC" "${NVFLAGS[@]}" -Xptxas -v -x cu -c "$K1/k1poly.cpp" -I"$K1" -o "$O/k1poly_gpu.o" 2> "$O/ptxas_k1.txt" \
  || { cat "$O/ptxas_k1.txt"; exit 1; }
# registers per kernel (entry name, then its "Used N registers" line); full occupancy on CC 12.0 needs <= 40 at 128
# threads per block (64 K registers / 1536 threads)
grep -E "Compiling entry|registers" "$O/ptxas_k1.txt" | sed -E "s/.*entry function '([^']+)'.*/\1/" | (c++filt || cat) \
  | paste - - | sed -E 's/void k1poly::\(anonymous namespace\)::kernel<k1poly::\(anonymous namespace\)::([A-Za-z]+)>.*ptxas info *: /\1 /' \
  | head -40 || true
"$NVCC" "${NVFLAGS[@]}" --fmad=false -x cu -c "$GC/gpucheck.cpp" -I"$GC" -o "$O/gpucheck_gpu.o"
"$CXX" -std=c++17 -O2 -ffp-contract=off -c "$GC/gpucheck.cpp" -I"$GC" -o "$O/gpucheck_cpu.o"
"$CXX" -std=c++17 -O2 -DWITH_GPUCHECK -I"$K1" -I"$GC" -c "$K1/test_k1.cpp" -o "$O/test_k1.o"
"$CXX" -std=c++17 -O2 -DWITH_GPUCHECK -I"$K1" -I"$GC" -c "$K1/k1_dump.cpp" -o "$O/k1_dump.o"
"$CXX" -std=c++17 -O2 -I"$K1" -I"$GC" -c "$K1/bench_k1.cpp" -o "$O/bench_k1.o"
"$NVCC" "-arch=$ARCH" "$O/test_k1.o" "$O/k1poly_gpu.o" "$O/gpucheck_cpu.o" -o "$O/test_k1_gpu"
"$NVCC" "-arch=$ARCH" "$O/k1_dump.o" "$O/k1poly_gpu.o" "$O/gpucheck_cpu.o" -o "$O/k1_dump_gpu"
"$NVCC" "-arch=$ARCH" "$O/bench_k1.o" "$O/k1poly_gpu.o" "$O/gpucheck_gpu.o" -o "$O/bench_k1_gpu"
set +x
echo "build_gpu: OK ($ARCH) -> $O"
