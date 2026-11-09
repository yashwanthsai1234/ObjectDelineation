#!/usr/bin/env bash
# build_cpu.sh <REPO/cuda> -- CPU builds (clang++ on macOS, g++ on Linux) into build/cpu/:
#   test_k1_cpu{,_o1,_o2}  exactness tests; _o1/_o2 run every step's items in reverse / permuted order (K1_ORDER)
#   k1_dump_cpu            raw-file runner for cross_numpy.py;  bench_k1_cpu  K1 vs check_masks, both CPU builds
set -euo pipefail
K1=$(cd "$(dirname "$0")" && pwd); GC=$(cd "$1" && pwd); CXX=${CXX:-c++}; O=$K1/build/cpu; mkdir -p "$O"
F=(-std=c++17 -O2 -Wall -Wextra -I"$K1" -I"$GC")
"$CXX" "${F[@]}" -DWITH_GPUCHECK "$K1/k1poly.cpp" "$K1/test_k1.cpp" "$GC/gpucheck.cpp" -o "$O/test_k1_cpu"
for o in 1 2; do
  "$CXX" "${F[@]}" -DK1_ORDER=$o -DWITH_GPUCHECK "$K1/k1poly.cpp" "$K1/test_k1.cpp" "$GC/gpucheck.cpp" -o "$O/test_k1_cpu_o$o"
done
"$CXX" "${F[@]}" -DWITH_GPUCHECK "$K1/k1poly.cpp" "$K1/k1_dump.cpp" "$GC/gpucheck.cpp" -o "$O/k1_dump_cpu"
"$CXX" "${F[@]}" "$K1/k1poly.cpp" "$K1/bench_k1.cpp" "$GC/gpucheck.cpp" -o "$O/bench_k1_cpu"
echo "build_cpu: OK -> $O"
