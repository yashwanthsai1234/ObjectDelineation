# K1: the CUDA checker on an RTX 3090, with sanitizers (2026-11-06)

**Run:**
- Box: one RTX 3090 (sm_86), 23 CPUs (cgroup v1 quota 2,304,000 / 100,000 of `nproc` 192), driver 615.71 (CUDA 13.4), Oregon, US.
- Image: `vllm_native-6008baf45`.
- Sync: code only, with no model and no data bundle.
- Chain: `experiments/phase2/configs/checker_gpu.yaml`.
- Cost: $0.096 estimated at teardown, $0.030 billed at teardown (Vast bills late).
- Phase 2 test allocation: $1.26 (P1) + K1, within $1.50.

## Sanitizers

| Tool | Tests | Result |
|---|---|---|
| memcheck, `tests/test_cuda_checker_gpu.py` | 6 passed | 0 errors |
| initcheck, same file | 6 passed | 0 errors |
| memcheck, the full GPU group (`run_box.sh`'s four files) | 40 passed, 2 skipped | 0 errors |
| initcheck, the full GPU group | 40 passed, 2 skipped | 0 errors |

All four ran with `PYTORCH_NO_CUDA_MEMORY_CACHING=1`, so torch's pool does not hide out-of-bounds accesses. The two skipped tests read local data fixtures that are not shipped.

**Root cause of the failed memcheck step:**
- The chain's first attempt failed with "Target application terminated before first instrumented API call" and no output from the target. The cause is that `cuda-toolkit[sanitizer]==13.0.3` installs `nvidia-cuda-sanitizer-api` 13.0.85. That wheel ships the binary and libraries, but not the `TreeLauncherSubreaper` and `TreeLauncherTargetLdPreloadHelper` executables the tool needs to launch its target. Colocating its files did not help.
- The 13.4.92 wheel (SHA-256 `d36bc50b…5103`, matching the host driver's CUDA 13.4) has them. Installed with uv, its `nvidia/cu13/bin/compute-sanitizer` wrapper ran memcheck on a torch program with 0 errors.
- `run_box.sh` now pins `nvidia-cuda-sanitizer-api==13.4.92`.

## Equivalence on sm_86

P1 measured these on sm_120; this run repeats them on a second architecture:

- **Fuzz polygons:** 8,559 cases, 0 mismatches.
- **Fuzz masks:** 6,024 cases, 0 mismatches.
- **LIDC sets:** absent, since the run synced code only.
- **GPU tests:** 40 passed, 2 skipped. `tests/test_gpu_scoring.py`: 6 passed.

Polygons per second at 16 vertices and batch 1,024: 3.29M on the 3090, against P1's 7.74M on the 5090.
