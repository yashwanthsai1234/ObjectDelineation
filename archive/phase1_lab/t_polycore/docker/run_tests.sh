#!/bin/sh
# Runs inside the container. /src = this sandbox (read-only), /common = lab oracle (read-only).
set -e
mkdir -p /w/pc && cd /src && tar --exclude=./build --exclude=./.hypothesis --exclude=./results -cf - . | tar -C /w/pc -xf -
cd /w/pc
echo "== $(uname -m) $(g++ --version | head -1) $(cmake --version | head -1) $(python --version)"
pip install --no-build-isolation -q . 2>&1 | tail -3
echo "== native fuzz (g++ -fsanitize=address,undefined)"
g++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all -Wall -Wextra -I cpp \
    cpp/polycore/tests/fuzz_main.cpp -o /tmp/polycore_fuzz && /tmp/polycore_fuzz ${FUZZ_N:-100000}
echo "== pytest"
cd /w/pc && POLYCORE_ORACLE_DIR=/common PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=/tmp/numba \
    python -m pytest -p no:cacheprovider -q tests "$@"
