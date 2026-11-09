PY := .venv/bin/python
# Sources of the objdel._gpucheck extension: editing them makes `make test` rebuild it first.
EXT_SRC := CMakeLists.txt cuda/gpucheck.h cuda/gpucheck.cpp cuda/gpucheck_py.cpp
EXT_STAMP := build/.gpucheck-stamp

.PHONY: setup build test test-all cov lint dead fmt wheels cuda preflight docs docs-serve headline

setup:            ## create .venv and install the package (editable, all extras); builds the C++ extension
	uv venv --python 3.12
	uv pip install --python $(PY) -e ".[train,data,dev]"
	@mkdir -p build && touch $(EXT_STAMP)

build: $(EXT_STAMP) ## rebuild the C++ extension if cuda/ changed (incremental, build/<wheel tag>)

$(EXT_STAMP): $(EXT_SRC)
	uv pip install --python $(PY) --no-deps --reinstall-package objdel -e .
	@mkdir -p build && touch $@

test: $(EXT_STAMP) ## run the test suite without the slow exhaustive runs
	$(PY) -m pytest -q -m "not slow" -n auto

test-all: $(EXT_STAMP) ## run every test, including the slow K1/K2 10k-case gates
	$(PY) -m pytest -q -n auto

cov: $(EXT_STAMP) ## run the fast suite with a line-coverage report of src/objdel; fails below 90 % (K13, as CI)
	$(PY) -m pytest -q -m "not slow" --cov=objdel --cov-report=term-missing:skip-covered --cov-fail-under=90

CONFIG ?= experiments/phase2/configs/r1_qwen3vl2b.yaml
MODEL_REF ?= svenk029/objdel-sft@rl-full-3ep
preflight: $(EXT_STAMP) ## prove data, checker and code on the Mac before renting a GPU (CONFIG=...); launch needs its stamp
	$(PY) -I -m experiments.phase2.preflight --config $(CONFIG)

lint:             ## lint library, tests and experiments, then the import-layer contracts
	$(PY) -m ruff check src tests experiments
	.venv/bin/lint-imports --no-logo

dead:             ## list unused code (vulture, [tool.vulture] in pyproject.toml)
	$(PY) -m vulture

fmt:              ## auto-fix lint issues
	$(PY) -m ruff check src tests experiments --fix

wheels:           ## build manylinux x86-64 + aarch64 wheels in Docker with cibuildwheel (build/wheelhouse)
	uvx cibuildwheel==4.3.1 --platform linux --output-dir build/wheelhouse

cuda:             ## rebuild objdel._gpucheck with nvcc (CUDA box; nvcc = $CUDACXX; CUDA_ARCHS="80;90" cross-compiles)
	uv pip install --python $(PY) --no-deps --reinstall-package objdel -e . -C cmake.define.OBJDEL_CUDA=ON \
		$(if $(CUDA_ARCHS),-C "cmake.define.CMAKE_CUDA_ARCHITECTURES=$(CUDA_ARCHS)")
	$(PY) -I -c "from objdel import _gpucheck as g; assert g.build() == 'cuda'; print('_gpucheck: cuda build,', g.device_count(), 'device(s)')"
	@mkdir -p build && touch $(EXT_STAMP)

docs:             ## build the docs site into site/ (strict: any warning fails)
	.venv/bin/zensical build --strict

docs-serve:       ## serve the docs site locally
	.venv/bin/zensical serve

headline:         ## copy the headline test and validation metrics into docs/results/headline.json (RUN=r5, or TEST=... VAL=...)
	@test -n "$(RUN)$(TEST)" || { echo "headline: give RUN=<id>, or TEST=<final_eval.json> and VAL=<final_eval.json>"; exit 2; }
	@test -n "$(RUN)" -o -n "$(VAL)" || { echo "headline: TEST needs VAL"; exit 2; }
	$(PY) -I -m experiments.phase2.headline $(if $(RUN),--run $(RUN),--test $(TEST) --val $(VAL)) --model-ref $(MODEL_REF)
