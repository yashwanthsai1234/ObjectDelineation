# Testing strategy

The test suite protects one claim: the reward is exact, and every part that touches it gives the same answer as a slow reference.
This page says what is tested, how, where, and what runs in which environment.
Facts below were checked on the collected suite: 1,484 tests, of which 1,453 run by default (`-m "not slow"`) and 31 are marked slow.
The count changes with each commit. Re-run `.venv/bin/python -m pytest --collect-only -q` to refresh it.

## 1. Principles

1. **Exactness over tolerance.** Geometry, raster, legality, topology and the scorer are compared for equality with exact references.
   Float terms are compared bit for bit, not with a tolerance.
2. **A slow oracle for every fast path.** An oracle is simple Python with no code shared with the checker. It is the specification.
3. **Record what you retire.** When polycore was retired, its results were recorded per case and kept as goldens.
4. **Break it first.** Each phase had an adversarial review. Every confirmed breakdown became a named regression test (`test_breakit_*`, `test_p2_breakit_*`).
5. **Tests run on a laptop.** The CPU build of the checker runs the same code as the GPU build, so the whole checker suite runs on a Mac.
6. **Skip, never fake.** A test that needs data, a GPU, or a heavy dependency skips itself when it is missing.

## 2. Layers of the pyramid

```text
                    /  GPU box tests  \         few, slow, run on a rented card
                   /  slow exhaustive  \        31 tests in 11 files, among them the 10,000-shape and 10,000-ring runs
                  /  data-backed tests  \       need data/, skip in CI
                 /  integration (TRL, runners) \  tiny offline model, CPU
                /   unit and property tests      \  most of the 1,453 default tests
```

| Layer | Examples | Environment |
|---|---|---|
| Unit | `test_geometry.py`, `test_topology.py`, `test_volume.py`, `test_reward.py`, `test_targets_frame.py` | Any machine |
| Property | `hypothesis` in `test_data_riga.py`, `test_data_store.py` | Any machine |
| Oracle equivalence | `test_checker_oracle.py`, `test_gpu_scoring.py`, `test_geometry_oracle.py`, `test_cuda_checker.py` | Any machine (CPU build), and a GPU box |
| Golden equivalence | `test_checker_host.py`, `test_cuda_checker.py`, `test_gpu_scoring.py` | Same |
| Integration | `test_train_trl.py` (tiny offline model), `test_train_fast.py`, `test_p2_scripts.py` (CPU train then select) | Needs the `train` extra |
| Data-backed | `test_data_lidc.py`, `test_data_riga.py`, `test_eval_data.py`, `test_reward_gates.py`, `test_p2_views.py` | Needs `data/`; marked `data` or skipped by fixtures |
| GPU | `test_cuda_checker_gpu.py` | Skips unless `nvidia-smi -L` lists a GPU |
| Operations | `test_vast.py`, `test_preflight.py`, `test_hf_pull.py`, `test_p2_hf_push.py`, `test_images.py` | Fake `vastai`, `ssh`, `rsync` on `PATH`. No real call |
| Documentation | `test_readme.py` | Executes the README quick start and checks `make` targets |

## 3. Fast and slow

| Command | What runs |
|---|---|
| `make test` | `pytest -q -m "not slow"`. Rebuilds the extension first if `cuda/` or `CMakeLists.txt` changed |
| `make test-all` | Every test, including the slow K1 and K2 gates |
| `make cov` | The fast suite with coverage of `src/objdel`. Fails under 90 % |
| `make lint` | `ruff check src tests experiments` |

Markers (`pyproject.toml`):

| Marker | Meaning |
|---|---|
| `slow` | Long runs. 31 tests in 11 files. Among them K1 (10,000 shapes against the raster oracle) and K2 (10,000 rings against the validity oracle) in `test_checker_oracle.py`. CI runs only that file in the `gates` job |
| `data` | Needs `data/fixtures` or `data/store`. Skipped when absent |

Coverage is measured across `multiprocessing` workers and threads (`[tool.coverage.run]`). The target is 90 % of `src/objdel` (gate K13).

## 4. Oracles (`tests/oracles/`)

| File | Reference for | How |
|---|---|---|
| `raster_oracle.py` | The even-odd raster of `check_polygons` | Pure Python, exact integers |
| `validity_oracle.py` | The 13 defects. Also checks `DEFECTS` names | Pure Python, exact integers |
| `topology_oracle.py` | Per-instance topology and the size rule | Independent labelling |
| `score_numpy.py` | The device scorer (IoU_tc, F5, g_area, signature) | The host scorer that the GPU scorer replaced, unchanged: NumPy and OpenCV, float64 sums |
| `reward_reference.py`, `reward_variants.py`, `reward_harness.py` | Reward v4 terms and advantages | Reference implementation and variants used by gates |

The scorer test demands bit-identical results from the device. This works because the sums are exact (see [CUDA engine, section 5](reference/cuda-engine.md#5-exactness-guarantees)).

## 5. Goldens (`tests/goldens/polycore.npz`)

The CPU checker polycore is archived in `archive/cpu_checker`. Nothing there is built or tested.
`archive/cpu_checker/record_goldens.py` recorded polycore's result for every case that the tests and `experiments/cuda_checker/eval.py` compare with: masks and the 13 defect counts.
The tests (`test_cuda_checker.py`, `test_checker_host.py`, `test_gpu_scoring.py`) demand equality with it.

A case that polycore never recorded takes its reference from the exact integer oracles (`test_unrecorded_cases_fall_back_to_the_exact_oracles`).
This happens when a generator rounds `math.cos` or `math.sin`, whose last bit differs between macOS ARM and Linux x86 (4 of about 7,600 cases on the first Linux box).

Do not regenerate goldens to make a test pass. A golden changes only with a reviewed reason.

## 6. What each area must test

| Area | Must test | Where |
|---|---|---|
| Parser (`geometry.text`) | Every error code. Near-misses 200,000 characters long (linear time). Equality with the lab parser | `test_geometry.py`, `test_breakit_geometry.py`, `test_geometry_oracle.py` (a slow differential oracle, not marked `slow`) |
| Canonical form | Orientation, start vertex, ordering. Input is not modified | `test_geometry.py` |
| Checker | 13 defects exact. Arena reuse over 200 calls. Two-stream chunks equal one call. One checker call per reward call | `test_checker_oracle.py`, `test_cuda_checker.py`, `test_cuda_checker_gpu.py`, `test_reward_one_checker_call.py` |
| Topology | Betti numbers, size rule idempotence, 3-D counts on worked shapes (hollow cube has one cavity, a ring one tunnel) | `test_topology.py`, `test_volume.py` |
| Reward | Terms, gates, advantages, determinism, TRL advantage equals `advantages_v4` | `test_reward.py`, `test_train_rewards.py`, `test_train_trl.py` |
| Targets | Every target parses, is canonical and legal, keeps the signature | `test_targets.py`, `test_targets_frame.py` |
| Data | Exact round trip, deterministic bytes, split leak checks | `test_data_store.py`, `test_data_splits.py`, `test_data_lidc.py`, `test_data_riga.py` |
| Training speed patches | Same log-probabilities, entropies and gradients as TRL (atol 1e-5 and 1e-4) | `test_train_fast.py`, `test_open_once.py`, `test_image_memo.py` |
| Vast operations | Offer filter, cost caps, watchdog, teardown always destroys | `test_vast.py`, `test_p2_breakit_ops.py` |
| Evaluation | Metrics, audit alarms, gates, signal panel, `final_eval` summary | `test_eval.py`, `test_p2_eval.py`, `test_p2_signal.py`, `test_final_eval.py` |

Edge cases to cover when you add a test: empty input, one-element input, maximum caps (1,000 vertices, 10 rings, 4,096 pixels), inputs on a pixel centre or an edge, integers near int32 limits, and unparseable or truncated text.

## 7. GPU tests and sanitizers

GPU tests run on a rented box through `experiments/cuda_checker/run_box.sh`, which runs:

1. `tests/test_cuda_checker_gpu.py`, `test_cuda_checker.py`, `test_checker_host.py`, `test_reward_one_checker_call.py`.
2. `experiments.cuda_checker.eval` even when a test failed (it shows where the GPU build goes wrong).
3. In the RL chains: `tests/test_gpu_scoring.py`.
4. compute-sanitizer memcheck, then initcheck, on the five checker and scorer test files, with `PYTORCH_NO_CUDA_MEMORY_CACHING=1`.

Recorded results:

| Run | Result | Source |
|---|---|---|
| RTX 5090 (P1) | 40 GPU tests pass, 2 skipped (they read local fixtures). 0 mismatches | `perf/p1_results.md` |
| RTX 3090 (K1) | memcheck and initcheck: 0 errors on 40 tests | `perf/k1_results.md` |
| RTX 5090 (R3, full RL chain) | 52 tests in 5 files: 0 errors | `docs/research/notes/issues_solved_2026-11-07.md`, `full_rl_2026-11-08.md` |

`test_profile.py` checks that the chains contain the sanitizer step and that it fails loudly on a missing binary.
`racecheck` and `synccheck` are not part of any chain.

The checker equivalence set `data/checker_eval` (LIDC only, 19 MB) holds 195,809 cases. It ships in the LIDC bundle. `run_box.sh` refuses to run if RIGA data is on the host.

## 8. Continuous integration (`.github/workflows/tests.yml`)

| Job | Runs on | What |
|---|---|---|
| `python` | Linux x86-64, Linux arm64, macOS arm64 | `uv pip install -e ".[dev]"` (builds the CPU extension), `ruff check`, `pytest -m "not slow"` |
| `train` | Linux | `.[dev,data,train]` with the CPU torch index. Fast suite with coverage. Fails under 90 % |
| `gates` | Linux | `pytest -m slow tests/test_checker_oracle.py`: the K1 and K2 exhaustive runs |
| `wheels` | Linux x86-64, arm64, macOS arm64 | `cibuildwheel` for `cp312`. The wheel test runs `test_checker_oracle.py`, `test_cuda_checker.py`, `test_checker_host.py` |

Third-party actions are pinned to full commit SHAs. The `image` workflow builds a Vast image by hand (`workflow_dispatch`). It is not a test job. No CI job has a GPU.

## 9. Gaps

| Gap | Effect |
|---|---|
| No GPU in CI | Kernel regressions show only on a rented box, through preflight's CPU build and the box chain |
| `racecheck`, `synccheck` not run | Races in the union-find and the ballot raster are covered only by exact-equality tests |
| The exam metrics read masks on the host. Few tests compare them with a device version | Acceptable today. They run a few times per run |
| RL on RIGA has no test | RL reads LIDC only |
| `docs/` is not tested | Only the README quick start is. Keep examples in these documents in sync by hand |
| The Controller is tested with a tiny model (`test_train_trl.py`) but no production script attaches it | See [Training guide](guides/training.md) |

## 10. Adding a test

1. Put it in `tests/test_<area>.py`. Name it for the behaviour.
2. Use the shared fixtures in `tests/conftest.py`: `fixtures_dir`, `store_dir` (skip when absent), the synthetic items and the tiny offline TRL model.
3. Compare with an oracle or a golden when the code is exact. Do not widen a tolerance.
4. Mark exhaustive runs `@pytest.mark.slow` and data-backed runs `@pytest.mark.data`.
5. Run `make test` and `make lint`. Both must pass before a commit.
