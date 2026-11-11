PY := .venv/bin/python
# Sources of the objdel._gpucheck extension: editing them makes `make test` rebuild it first.
EXT_SRC := CMakeLists.txt cuda/gpucheck.h cuda/gpucheck.cpp cuda/gpucheck_py.cpp
EXT_STAMP := build/.gpucheck-stamp

.DEFAULT_GOAL := help
.PHONY: help smoke setup build test test-all test-data test-order cov lint dead sweep fmt wheels cuda gpu-check preflight docs docs-serve headline

help:             ## list the targets (the default: a bare `make` changes nothing)
	@grep -hE '^[a-zA-Z_-]+:.*## ' $(MAKEFILE_LIST) | awk -F':.*## ' '{printf "  %-14s %s\n", $$1, $$2}'

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

test-data: $(EXT_STAMP) ## run the data tests (read data/ and runs/, so under -I), slow ones included
	$(PY) -I -m pytest -q -m data

test-order: $(EXT_STAMP) ## build the CPU test variants and run them (K1 order, poison scan, mutfg4); the planted alias must fail
	$(PY) -I -m tests.variants

cov: $(EXT_STAMP) ## run the fast suite with a line-coverage report of src/objdel; fails below 90 % (K13, as CI)
	$(PY) -m pytest -q -m "not slow" --cov=objdel --cov-report=term-missing:skip-covered --cov-fail-under=90

CONFIG ?= specs/b1_rows.yaml
MODEL_REF ?= svenk029/objdel-sft@b1-rows
preflight: $(EXT_STAMP) ## prove data, checker and code on the Mac before renting a GPU (CONFIG=... [BUNDLE_DIR=...]); launch needs its stamp
	$(PY) -I -m objdel.cli.box.preflight --config $(CONFIG) $(if $(BUNDLE_DIR),--bundle-dir $(BUNDLE_DIR))

smoke:            ## the CPU smoke test: train -> select on a tiny model (preflight's last check)
	$(PY) -m pytest -q -p no:cacheprovider tests/experiments/test_p2_scripts.py::test_cpu_smoke_train_then_select

lint:             ## lint library, tests, tools and ops, then the import-layer contracts
	$(PY) -m ruff check src tests tools ops
	.venv/bin/lint-imports --no-logo

dead:             ## list unused code (vulture, [tool.vulture] in pyproject.toml)
	$(PY) -m vulture

sweep:            ## final sweep: stale names, orphans, doc links, docstrings, dead code (exit 1 on a hard finding)
	$(PY) -I -m tools.sweep

fmt:              ## auto-fix lint issues
	$(PY) -m ruff check src tests tools ops --fix

wheels:           ## build manylinux x86-64 + aarch64 wheels in Docker with cibuildwheel (build/wheelhouse)
	uvx cibuildwheel==4.3.1 --platform linux --output-dir build/wheelhouse

cuda:             ## rebuild objdel._gpucheck with nvcc (CUDA box; nvcc = $CUDACXX; CUDA_ARCHS="80;90" cross-compiles)
	uv pip install --python $(PY) --no-deps --reinstall-package objdel -e . -C cmake.define.OBJDEL_CUDA=ON \
		$(if $(CUDA_ARCHS),-C "cmake.define.CMAKE_CUDA_ARCHITECTURES=$(CUDA_ARCHS)")
	$(PY) -I -c "from objdel import _gpucheck as g; assert g.build() == 'cuda'; print('_gpucheck: cuda build,', g.device_count(), 'device(s)')"
	@mkdir -p build && touch $(EXT_STAMP)

gpu-check:        ## the GPU check chain on a CUDA box (run_box.sh --gpu-check); DRY=1 prints it and collects its tests
	PY=$(PY) bash ops/box/run_box.sh --gpu-check $(if $(DRY),--dry)

.PHONY: old-engine-tar
old-engine-tar:   ## the pre-D2 engine (owner's tag pre-d2) as old_engine.tar + .sha256 in the repo root, for D-G7 on the box
	@git rev-parse -q --verify 'pre-d2^{commit}' >/dev/null || { echo "owner: git tag pre-d2 <commit>, see D-G7" >&2; exit 1; }
	git archive -o old_engine.tar pre-d2 cuda src CMakeLists.txt pyproject.toml Makefile README.md
	sha256sum old_engine.tar > old_engine.tar.sha256

docs:             ## build the docs site into site/ (strict: any warning fails)
	.venv/bin/zensical build --strict

docs-serve:       ## serve the docs site locally
	.venv/bin/zensical serve

headline:         ## copy the headline test and validation metrics into reference/headline.json (RUN=r5, or TEST=... VAL=...)
	@test -n "$(RUN)$(TEST)" || { echo "headline: give RUN=<id>, or TEST=<final_eval.json> and VAL=<final_eval.json>"; exit 2; }
	@test -n "$(RUN)" -o -n "$(VAL)" || { echo "headline: TEST needs VAL"; exit 2; }
	$(PY) -I -m objdel.cli.eval.headline $(if $(RUN),--run $(RUN),--test $(TEST) --val $(VAL)) --model-ref $(MODEL_REF)
