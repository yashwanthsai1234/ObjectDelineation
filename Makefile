PY := .venv/bin/python

.PHONY: setup test lint fmt cuda

setup:            ## create .venv and install the package with all extras
	uv venv --python 3.12
	uv pip install --python $(PY) -e ".[models,data,dev]"

test:             ## run the Python test suite
	$(PY) -m pytest -q

lint:             ## lint library, tests and experiments
	$(PY) -m ruff check src tests experiments

fmt:              ## auto-fix lint issues
	$(PY) -m ruff check src tests experiments --fix

cuda:             ## build the C++/CUDA engine and its tests
	cmake -S cuda -B build && cmake --build build -j
