# K3 on one rented RTX 5090: build, test, bench

Goal: prove the three K3 kernels on sm_120 and time them, in **under 20 minutes of GPU time** (budget $3; at a
typical $0.40-0.60/h for a 5090 that is about $0.15-0.20 [estimate]). Nothing here trains or downloads a model's
weights; no dataset is copied (the tests make their own inputs).

Assumes the Phase 2 image (python 3.12, torch 2.14.1 in `/opt/venv`, no nvcc), the REPO rsynced to the box, and the
Qwen3-VL-2B processor files in the box's HF cache (if not, drop `HF_HUB_OFFLINE=1`: `AutoProcessor` then fetches only
the tokenizer and processor configs, a few MB). Follow the repo's preflight before renting (North-America host,
`docs/runbooks/gpu-run.md`).

## 0. Copy (on the Mac, ~10 s)

```bash
REPO=$(git rev-parse --show-toplevel)   # run inside the clone
K3=$REPO/docs/implementation2.0/kernels/K3/
rsync -az --exclude build "$K3"/ root@$HOST:/root/K3/
ssh root@$HOST
```

## 1. Toolchain (~1-2 min): CUDA 13.0 nvcc + CCCL, as `experiments/cuda_checker/run_box.sh` does

```bash
PY=/opt/venv/bin/python
CONSTRAINTS=(); [ -f /opt/requirements.txt ] && CONSTRAINTS=(-c /opt/requirements.txt)
uv pip install --python "$PY" "${CONSTRAINTS[@]}" "cuda-toolkit[nvcc,cccl]==13.0.3" "nvidia-cuda-sanitizer-api==13.4.92"
export CUDA_HOME=$("$PY" -c "import nvidia, os; print(os.path.join(nvidia.__path__[0], 'cu13'))")
[ -e "$CUDA_HOME/lib64" ] || ln -s lib "$CUDA_HOME/lib64"
ln -sf libcudart.so.13 "$CUDA_HOME/lib/libcudart.so"
export PATH=$CUDA_HOME/bin:$PATH
nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total --format=csv   # expect 12.0
nvcc --version | tail -1
```

## 2. Build the standalone programs (~1-2 min)

```bash
cd /root/K3
make all ARCH=sm_120 2>&1 | tee build.log
```

`make` builds `build/test_cpu`, `build/k3_cpu_cli` (C++ only) and `build/test_gpu`, `build/bench_gpu` (nvcc,
`-arch=sm_120 --fmad=false`). For register / shared-memory use add `NVFLAGS+=" -Xptxas -v"` (the lse kernel should
show 0 spill bytes and 33,792 bytes of static shared memory).

## 3. Correctness (~2-4 min)

```bash
./build/test_cpu        | tee cpu.log     # seconds: formulas, exp/log accuracy, CPU reference
./build/test_gpu        | tee gpu.log     # GPU vs CPU reference; the real head shape (192 x 151,936 x 2,048)
                                          # needs ~10-40 s of multithreaded CPU reference
HF_HUB_OFFLINE=1 "$PY" -I test_hf_cpu.py --cli build/k3_cpu_cli --tmp build/hf_tmp | tee hf_cpu.log
```

Pass criteria (each program prints `all passed` and exits 0):
- (a) preprocessing: GPU fp32 and bf16 output **bit-equal** to the CPU reference (256 rows, HWC and CHW); the CPU
  reference bit-equal to the HF processor (`test_hf_cpu.py`, already passed on the Mac).
- (c) rope: GPU == host emulation == C++ port of HF `get_rope_index`; the three adversarial rows raise the error flag.
- (b) lse exact mode: **bit-equal** to the CPU reference (logp and entropy); tensor-core mode within 2e-3 of fp64.
- (b) is_ratio (portable exp): ratios and the 7 stats **bit-equal** to the CPU reference, all 6 mode/bound cases.

## 4. Against Hugging Face / TRL / objdel on the GPU (~3-5 min, of which ~1.5 min builds the torch extension)

```bash
TORCH_CUDA_ARCH_LIST=12.0 HF_HUB_OFFLINE=1 "$PY" -I test_torch_gpu.py --out build/torch_gpu.json | tee torch.log
```

Checks: (a) GPU output vs the HF processor's own `pixel_values` (CPU) bit-equal, fp32 / bf16 / 256-row layout;
(c) GPU vs HF `get_rope_index` on a real prompt, 256 rows; (b) `is_ratio` (libdevice exp) vs `objdel.train.fast.is_ratio`
on CUDA, bit-equal for the token modes; (b) lse errors vs fp64 next to TRL's `ChunkedLogProbFunction`.
Timings land in `build/torch_gpu.json` (HF CPU vs GPU for preprocessing and rope; K3 vs TRL for the head; K3 vs torch
for the IS ratio; extra memory of K3 vs TRL).

## 5. Micro-benchmarks (~1 min)

```bash
./build/bench_gpu | tee bench.log
```

## 6. Optional, if time is left (~2-4 min)

```bash
SAN=$(command -v compute-sanitizer || find / -name compute-sanitizer -type f 2>/dev/null | head -1)
[ -n "$SAN" ] && "$SAN" --tool memcheck  ./build/test_gpu --quick | tail -3
[ -n "$SAN" ] && "$SAN" --tool racecheck ./build/test_gpu --quick | tail -3
# Nsight Compute, only if `ncu` exists on the image:
command -v ncu && ncu --kernel-name regex:lse_kernel --launch-skip 2 --launch-count 1 --set full ./build/bench_gpu > ncu_lse.txt
```

## 7. Copy back (on the Mac) and stop the instance

```bash
mkdir -p "$K3"/results
scp root@$HOST:/root/K3/{build.log,cpu.log,gpu.log,hf_cpu.log,torch.log,bench.log} "$K3"/results/ 2>/dev/null || true
scp root@$HOST:/root/K3/build/torch_gpu.json "$K3"/results/
```

## If something fails

| Symptom | Likely cause | Action |
|---|---|---|
| `nvcc fatal: Unsupported gpu architecture 'compute_120'` | old nvcc | step 1 (CUDA 13.0 from PyPI) |
| link error `-lcudart_static` / `libcudart.so` | PyPI layout | the two `ln -s` lines of step 1 |
| torch extension: "C++20 or later compatible compiler is required" | a `-std=c++17` flag | `test_torch_gpu.py` passes none; do not add one |
| (a) differs from HF only on GPU | fast-math or FTZ somewhere | check `--fmad=false` and no `--use_fast_math` in the build line |
| (b) exact mode differs | implicit FMA contraction or FTZ | same; then diff `lse_ref` vs GPU per split (DESIGN.md, debugging plan) |
| (b) token-mode `is_ratio` differs from torch by 1 ulp | torch.exp not libdevice `expf` on this build | report the ulp count; keep the portable exp and a tolerance |
| lse kernel "too many resources requested" | register use | build with `-maxrregcount=128` or lower `LSE_PER_LANE` unroll |
