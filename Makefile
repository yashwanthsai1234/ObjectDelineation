PY := .venv/bin/python
FUZZ_N ?= 20000
# Sources of the objdel._polycore extension: editing them makes `make test` rebuild it first.
EXT_SRC := CMakeLists.txt cpp/polycore/polycore.hpp src/objdel/_polycore.cpp cuda/gpucheck.h cuda/gpucheck.cpp cuda/gpucheck_py.cpp
EXT_STAMP := build/.polycore-stamp

.PHONY: setup build test test-all cov lint fmt fuzz wheels cuda preflight

setup:            ## create .venv and install the package (editable, all extras); builds the C++ extension
	uv venv --python 3.12
	uv pip install --python $(PY) -e ".[train,data,dev]"
	@mkdir -p build && touch $(EXT_STAMP)

build: $(EXT_STAMP) ## rebuild the C++ extension if cpp/ or its binding changed (incremental, build/<wheel tag>)

$(EXT_STAMP): $(EXT_SRC)
	uv pip install --python $(PY) --no-deps --reinstall-package objdel -e .
	@mkdir -p build && touch $@

test: $(EXT_STAMP) ## run the test suite without the slow exhaustive runs
	$(PY) -m pytest -q -m "not slow"

test-all: $(EXT_STAMP) ## run every test, including the slow K1/K2 10k-case gates
	$(PY) -m pytest -q

cov: $(EXT_STAMP) ## run the fast suite with a line-coverage report of src/objdel; fails below 90 % (K13, as CI)
	$(PY) -m pytest -q -m "not slow" --cov=objdel --cov-report=term-missing:skip-covered --cov-fail-under=90

CONFIG ?= experiments/phase2/configs/r1_lfm.yaml
preflight: $(EXT_STAMP) ## prove data, checker and code on the Mac before renting a GPU (CONFIG=...); launch needs its stamp
	$(PY) -I -m experiments.phase2.preflight --config $(CONFIG)

lint:             ## lint library, tests and experiments
	$(PY) -m ruff check src tests experiments

fmt:              ## auto-fix lint issues
	$(PY) -m ruff check src tests experiments --fix

fuzz:             ## native polycore fuzz under ASan/UBSan (FUZZ_N random shapes)
	@mkdir -p build
	$(CXX) -std=c++17 -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all -Wall -Wextra -Werror \
		-I cpp cpp/polycore/tests/fuzz_main.cpp -o build/polycore_fuzz
	build/polycore_fuzz $(FUZZ_N)

wheels:           ## build manylinux x86-64 + aarch64 wheels in Docker with cibuildwheel (build/wheelhouse)
	uvx cibuildwheel==4.3.1 --platform linux --output-dir build/wheelhouse

cuda:             ## rebuild objdel._gpucheck with nvcc (CUDA box; nvcc = $CUDACXX; CUDA_ARCHS="80;90" cross-compiles)
	uv pip install --python $(PY) --no-deps --reinstall-package objdel -e . -C cmake.define.OBJDEL_CUDA=ON \
		$(if $(CUDA_ARCHS),-C "cmake.define.CMAKE_CUDA_ARCHITECTURES=$(CUDA_ARCHS)")
	$(PY) -I -c "from objdel import _gpucheck as g; assert g.build() == 'cuda'; print('_gpucheck: cuda build,', g.device_count(), 'device(s)')"
	@mkdir -p build && touch $(EXT_STAMP)
