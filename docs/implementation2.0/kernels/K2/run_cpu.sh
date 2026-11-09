#!/usr/bin/env bash
# CPU build of K2 (macOS clang++ or Linux g++): tests at scale 1 and 4, the sanitizer build, a small benchmark.
# Expected: ALL PASS, DIGEST b158d8132b18be83 (seed 20261108, scale 1) and f92b9a3939cec10f (scale 4).
set -euo pipefail
cd "$(dirname "$0")"
CXX=${CXX:-c++}
mkdir -p build
FLAGS="-std=c++17 -O2 -ffp-contract=off -Wall -Wextra -Wno-unused-parameter"
$CXX $FLAGS k2.cpp k2ref.cpp k2test.cpp -o build/k2test_cpu
$CXX $FLAGS k2.cpp k2ref.cpp k2bench.cpp -o build/k2bench_cpu
./build/k2test_cpu 20261108 1 | tail -3
./build/k2test_cpu 20261108 4 | tail -3
$CXX -std=c++17 -O1 -g -ffp-contract=off -fsanitize=address,undefined -fno-sanitize-recover=undefined \
    k2.cpp k2ref.cpp k2test.cpp -o build/k2test_asan
./build/k2test_asan 20261108 1 | tail -2
./build/k2bench_cpu 8 1
