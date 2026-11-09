# Runbook: troubleshooting

This page lists failures that the repository has met, with the cause and the fix.
Each entry names its source. When a failure is new, add an entry here and a test next to the fix.
Method for a new failure: reproduce, isolate the component, form a hypothesis and test it, find the root cause, fix it, add a test.

## 1. Training on the box

| Symptom | Cause | Fix | Source |
|---|---|---|---|
| FastGRPO with CUDA Graphs (`--compile reduce-overhead`) crashes | Per-layer graph replay overwrites outputs. It conflicts with vLLM's sleep-mode allocator | Do not use `reduce-overhead`. Use no compile, or `--compile default` (47.9 s per step against 44.3 s, so not worth it) | `README.md`, P1 table; `experiments/phase2/perf/p1_results.md` |
| Out of memory in the update | Gradient checkpointing is off. The step needs more than the 31.45 GiB left beside vLLM | Keep `gradient_checkpointing=True` (the default in `grpo_config`). Micro-batch 16 x 16 beats 32 x 8 (44.3 s against 59.1 s per step) | `README.md`; `configs.py` comment |
| OOM on SFT with a long batch | Qwen3-VL's 152k-token vocabulary makes the logits large | Micro-batch 4, accumulation 8 (as `r1_qwen3vl2b.yaml`) | `r1_qwen3vl2b.yaml` |
| `RuntimeError: vLLM expanded the image to N pads, the trainer's processor to M` | vLLM's processor expands the image tokens differently from the trainer's. vLLM would sample from other image features than the trainer scores | A real mismatch. Stop. `patch_vllm_prompts` sends one pad per image and checks each prompt. Check the processor and vLLM versions | `fast.py` |
| Log-probability gap between vLLM and the trainer is large (about 25 nats), importance ratio mean 0.65 to 0.81 | TRL as shipped sends processor-expanded prompt ids to vLLM for Qwen3-VL | Use `FastGRPO` (`one_pad_per_image`). Healthy: gap mean near 0.01 nats | `perf/p1_results.md` |
| `ValueError: OpenOnce keeps the step's stack in the main process: dataloader_num_workers must be 0` | `OpenOnce` is main-process only | Set `dataloader_num_workers` to 0 for RL | `fast.py` |
| `ValueError: reward call must receive whole groups of 8` | `per_device_train_batch_size * steps_per_generation` or the eval batch is not a multiple of 8 | Use multiples of 8. `grpo_config` guarantees it | `train/rewards.py` |
| `ValueError: num_generations=... reward v4 scores groups of exactly 8` | Another `num_generations` was passed | Remove it. The reward cannot score other group sizes | `train/configs.py`, `callbacks.py` |
| `StaleCacheError: built under another reward fingerprint` | `REWARD_V4`, the size rule, or the numpy, scipy or OpenCV version changed | Re-run `prepare_rl`, or `ExpertCache.load(..., rebuild=True)` | `reward/cache.py` |
| `train_sft` exits with code 3 | G-vision abort: projector ratio below 0.5 and TV below 0.002. The vision features collapsed | Restart once with `freeze_vision: true` (the D7 repair). Check LoRA excludes the vision tower | `train_sft.py`, `configs.lora` |
| `train_sft` refuses to start: missing G-vision reference | R0 did not write it | Run `train_sft --config CONFIG --vision_reference_only true` first | `train_sft.py` |
| KeyError `'inner'` at the first checkpoint save | `ImageMemo` was deep-copied by TRL with the processor | Fixed in commit `a300197`: `__deepcopy__` returns the plain processor. If it returns, check that fix | `fast.py` |
| Single-token KL spikes in the logged mean (1020, 5.9) | vLLM sampled a token the trainer gives near-zero probability | Harmless at the step seen. Judge drift by the median KL (0.19 to 0.40). Gradient norm and rewards did not change | `docs/research/notes/full_rl_2026-11-08.md` |
| `CUDA Error: invalid argument at cumem_allocator.cpp:249` at exit | vLLM's sleep-mode allocator shutting down | Ignore if the exit code is 0 and `final/` loads | same |
| Base weights differ from the SFT checkpoint after RL | TRL merges the adapter into bf16 weights in place and unmerges. `W + d - d` is not `W` (TRL #6688) | `FastGRPO.sync_without_merge` (`lora_merged_params`). One process only | [ADR-0011](../architecture/decisions/0011-lora-weight-sync-without-merging.md) |
| RL IoU stops improving after epoch 1 | RL trains on centred frames. Exams use shifted frames | `train_rl --train-views` with the shifted and rescaled train windows (`rl_shifted.yaml`) | `docs/research/notes/rl_plateau_diagnosis_2026-11-08.md` |
| Run lost with the box | Box failure or cap | `--resume` with the last `trainer/checkpoint-N` (every 50 steps) | `train_rl.py`; [GPU run runbook](gpu-run.md) |
| `--resume: no checkpoint in ...; starting at step 0` | There is no checkpoint folder | Expected on a first run | `train_rl.py` |

## 2. Scoring pool and data

| Symptom | Cause | Fix | Source |
|---|---|---|---|
| `BrokenProcessPool` within a second | A scoring worker died (OOM kill). Pending calls fail fast | The run fails loudly. Reduce `--procs`, or free memory. Resume from a checkpoint | `train/serve.py` |
| `RuntimeError` when calling `serve.fork_pool` | CUDA is already initialised | Fork the pool before the model loads | `train/serve.py` |
| `RuntimeError` from `serve.preload` | A pool is alive | `preload` every root before the fork. Preload raises while a pool is alive | `train/serve.py` |
| A worker crashes on macOS after it imports a module | A module imported after the fork | `fork_pool` imports `SCORERS` first. Pass extra `modules` to it | `train/serve.py` |
| RIGA val pass uses over 1.5 GB | Row groups decoded in batches | `BATCH_ROWS = 1` (455 MB peak) | `data/store.py` |
| Mask has bytes other than 0 and 1 | A PIL "1" image stores 255. bool casts keep the bytes | Use `as_mask`. It casts to uint8 and views back | `checker/topology.py` |
| OpenCV threads in worker pools | `cv2.setNumThreads(1)` is a no-op on macOS GCD | Call `cv2.setNumThreads(0)` | `checker/topology.py` |
| Results differ between arm64 and x86 | OpenCV's x86 distance transform (Intel IPP) is not correctly rounded | The code uses scipy's exact EDT. If `reward.advantages` still differs on x86 (AVX-512 SVML), set `NPY_DISABLE_CPU_FEATURES="X86_V4"` for the reward process | `reward/cache.py`; `experiments/phase1/RESULTS.md` section 8 |
| A prepared image fails its SHA-256 | The PNG was changed after `prepare` | Rebuild the prepared set. Rows refuse to train on it | `train/dataset.py` |

## 3. CUDA checker and build

| Symptom | Cause | Fix | Source |
|---|---|---|---|
| `build()` returns `cpu` on a GPU box | The CPU extension is installed | Run `make cuda` | `Makefile` |
| `make cuda` cannot find `nvcc` or `libcudart.so` | The PyPI toolkit layout differs from what CMake expects | `run_box.sh` links `lib64` and an unversioned `libcudart.so`, and sets `CUDACXX`. Use it | `experiments/cuda_checker/run_box.sh` |
| `CMAKE_CUDA_ARCHITECTURES` empty after `make cuda` | A target copies the variable when it is created | Fixed: the arch is set on the target (`CMakeLists.txt`, R0 v3) | `CMakeLists.txt` |
| Kernels are half as fast as expected | nanobind builds with `-Os` | The module uses `NOMINSIZE` | `CMakeLists.txt` |
| Scores differ by one rounding | FMA contraction is on | Build with `-ffp-contract=off` and `--fmad=false` | `CMakeLists.txt` |
| compute-sanitizer: "Target application terminated before first instrumented API call" | The 13.0.85 wheel lacks its `TreeLauncher` executables | Pin `nvidia-cuda-sanitizer-api==13.4.92` (done in `run_box.sh`) | `perf/k1_results.md` |
| memcheck hides out-of-bounds accesses | torch's caching allocator pools memory | Run with `PYTORCH_NO_CUDA_MEMORY_CACHING=1` | `rl_full.yaml` |
| 4 of about 7,600 checker test cases have no recorded reference on a Linux box | A test generator rounds `math.cos` and `math.sin`. The last bit differs between macOS ARM and Linux x86 | Unrecorded cases fall back to the exact integer oracles | `tests/test_checker_oracle.py`; `issues_solved_2026-11-07.md` |
| `make test` fails on the GIL-release check on a busy runner | Contention | The test tolerates 0.28 of the solo rate. A held GIL gives about 0 | `tests/` (commit `4ad4158`) |
| nvcc build check on the Mac | No GPU | `bash cuda/docker_check.sh`. With no GPU the calls must fail with a clean CUDA error | `cuda/docker_check.sh` |
| A bad input reaches a kernel | It cannot. The binding checks caps before releasing the GIL | `ValueError` or `TypeError` from `_gpucheck` names the cap | `cuda/gpucheck_py.cpp` |
| Calls from threads wait | One mutex serialises calls | By design. The call releases the GIL | `cuda/gpucheck.cpp` |
| CUDA error from a failed call poisons the next | The runtime keeps the last error | Each call clears it first (`clear_error`) | `cuda/gpucheck.cpp` |

## 4. Vast AI and image

| Symptom | Cause | Fix | Source |
|---|---|---|---|
| Image pull times out on several hosts | A single 4.14 GB layer that Docker cannot resume. Four full-RL launches on four hosts timed out | Use the four-layer image `vllm_native-3ef6e9179051` (largest layer 1.54 GB). Pick another host if a pull still stalls | commits `0a4b93f`, `9f83c07`; `Dockerfile.vllm_native` |
| SSH refused right after launch | Some hosts write `authorized_keys` with a foreign owner and sshd refuses it | The on-start command fixes modes for 120 s. `teardown.sh` retries the pull three times, and always destroys | `policy.onstart_cmd` |
| `ssh host CMD` fails with "not a terminal" | Vast's login tmux | The on-start command touches `/root/.no_auto_tmux` | `policy.onstart_cmd` |
| Launch fails with "onstart is N characters, over Vast's 4048" | A long run id or a bigger watchdog | Shorten the run id or `watchdog.sh` | `policy.onstart_cmd` |
| A card is busy or power-capped | A bad host (machine 28360: 150 W limit, 7.4 of 165 bf16 TFLOPS) | The GPU check runs before the bundle upload. Add the machine to `BAD_MACHINES` | `policy.py`; commit `999969b` |
| A host in Thailand cannot pull the image | Region | Only US and CA hosts are admitted | `policy.REGION` |
| Triton `@jit` fails in `env_check` | Triton reads its kernel's source, so the check must run from a file, not stdin | `env_check.sh` runs checks from a file | commit `28c4c72` |
| `launch.sh` refuses: no preflight stamp or stale stamp | Code, data or config changed after preflight | `make preflight CONFIG=...` as the last step | `preflight.refusal` |
| `fetch_hf.sh` or `push_hf.sh` refuses the token | The HF token is too broad for a rented box | Make a fine-grained token with write access to `svenk029/objdel-sft` only | `hf_push.token_problems` |
| Teardown on a box that never got SSH raised `KeyError 'port'` | No port in the run record | Fixed in commit `74a0135`: no port means no pull, only destroy | `policy.cmd_teardown` |
| Cost looks lower than the invoice | Vast bills late | The ledger estimate is `dph x time` plus bandwidth. The invoice row says `pending` until it arrives | `policy.cmd_teardown` |

## 5. Tests and CI

| Symptom | Cause | Fix | Source |
|---|---|---|---|
| Collection errors on the light CI job | The job installs `.[dev]` only. torch and trl are missing | Tests that need torch or trl skip themselves (`importorskip`). Keep that pattern | commits `24d12fb`, `51c0473` |
| Wheel test cannot import `experiments` | The test runs outside the repository | `test-sources` copies `tests`, `experiments`, `pyproject.toml` and `CMakeLists.txt`. The test command is `python -m pytest` | `pyproject.toml`, commits `8f51441`, `610a17c` |
| CPU dry run with `use_vllm` refused | TRL 1.17 refuses `use_vllm` without vLLM installed | The CPU dry run runs without vLLM | commit `421eac3` |
| `tests/test_readme.py` fails | The README quick start no longer prints what its comment says, or `make` targets changed | Fix the README or the target. The test executes the Python block | `tests/test_readme.py` |
| `make lint` fails on `archive/phase1_lab` | It is excluded on purpose (frozen evidence) | Do not lint lab files. They are not in the ruff scope | `pyproject.toml` |
