# K2 GPU run: one RTX 5090 (sm_120), budget under 20 GPU-minutes

What this run proves:
1. The GPU build compiles and runs.
2. The GPU build gives the same bits as the CPU build: the same `DIGEST`, and the exact-equality tests pass.
3. compute-sanitizer finds no memory, sync or init errors.
4. The benchmark gives the device time of each stage.

Already checked without a GPU, on 2026-11-08:
- macOS clang++ 21: tests and ASan/UBSan pass (`build/run_cpu_mac.txt`).
- Linux g++ (Debian trixie) CPU build: the same digest as macOS.
- nvcc 13.0.88 (`-arch=sm_120 --fmad=false`, PyPI `cuda-toolkit[nvcc,cccl]==13.0.3`) builds and links the GPU binaries.
  It has no warnings, no spills and a 0-byte stack frame for every kernel (`build/docker/ptxas.txt`).
- The GPU binary stops cleanly when there is no driver: "CUDA driver version is insufficient".

Expected digests (seed 20261108):

| Command | Expected digest |
|---|---|
| `k2test 20261108 1` | `b158d8132b18be83` |
| `k2test 20261108 4` | `f92b9a3939cec10f` |

## 0. Rental preflight (main agent)

- Use Vast AI, a North America host, one RTX 5090, a driver of R580 or newer (CUDA 13), and an image with g++.
  A CUDA 12.8+ devel image works. So does the Phase 2 image.
- Cost [estimate]: about 15 min at about $0.44/h comes to about $0.11. The 5090 price paid in R3 was $0.440/h
  (repo memory), so this is far below the $3 budget.
- Copy only this folder. Nothing in it reads `data/`.

```bash
REPO=$(git rev-parse --show-toplevel)   # run inside the clone
# from the Mac
rsync -av --exclude build $REPO/docs/implementation2.0/kernels/K2/ BOX:~/K2/
```

## 1. Toolchain (about 2 min)

```bash
nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv   # expect: NVIDIA GeForce RTX 5090, >=580, 12.0
which nvcc && nvcc --version                                           # needs CUDA >= 12.8 for sm_120
# No system nvcc: use the PyPI toolkit, pinned as in experiments/cuda_checker/run_box.sh.
python3 -m venv ~/cu && ~/cu/bin/pip install -q "cuda-toolkit[nvcc,cccl]==13.0.3" "nvidia-cuda-sanitizer-api==13.4.92"
export CUDA_HOME=$(~/cu/bin/python -c "import nvidia, os; print(os.path.join(nvidia.__path__[0], 'cu13'))")
[ -e "$CUDA_HOME/lib64" ] || ln -s lib "$CUDA_HOME/lib64"; ln -sf libcudart.so.13 "$CUDA_HOME/lib/libcudart.so"
export PATH=$CUDA_HOME/bin:$PATH
# compute-sanitizer from the pinned wheel (13.0.85 lacks TreeLauncher: docs/runbooks/troubleshooting.md):
export SAN=$(~/cu/bin/python -c "import nvidia, os, glob; print(glob.glob(os.path.join(nvidia.__path__[0], '**', 'compute-sanitizer'), recursive=True)[0])")
# With a system toolkit instead: export SAN=compute-sanitizer
```

## 2. Build (about 3 min)

```bash
cd ~/K2 && mkdir -p build
nvcc -std=c++17 -O3 -arch=sm_120 --fmad=false -Xcompiler -ffp-contract=off -lineinfo --resource-usage \
     -x cu -c k2.cpp -o build/k2_gpu.o 2> build/ptxas.txt
for f in k2ref k2test k2bench; do g++ -std=c++17 -O2 -ffp-contract=off -c $f.cpp -o build/$f.o; done
nvcc -arch=sm_120 build/k2_gpu.o build/k2ref.o build/k2test.o  -o build/k2test_gpu
nvcc -arch=sm_120 build/k2_gpu.o build/k2ref.o build/k2bench.o -o build/k2bench_gpu
g++ -std=c++17 -O2 -ffp-contract=off k2.cpp k2ref.cpp k2test.cpp  -o build/k2test_cpu     # same box, CPU build
g++ -std=c++17 -O2 -ffp-contract=off k2.cpp k2ref.cpp k2bench.cpp -o build/k2bench_cpu
```

## 3. Exactness tests (about 1 min)

```bash
./build/k2test_gpu 20261108 1 | tee build/test_gpu_s1.txt | tail -3   # ALL PASS, DIGEST b158d8132b18be83
./build/k2test_cpu 20261108 1 | tail -3                               # the same DIGEST
./build/k2test_gpu 20261108 4 | tee build/test_gpu_s4.txt | tail -3   # ALL PASS, DIGEST f92b9a3939cec10f
./build/k2test_gpu 99 4 | tail -3; ./build/k2test_cpu 99 4 | tail -3  # a new seed: the two digests must be equal
```

Pass rule:
- Every run prints `ALL PASS`.
- The GPU and CPU digests are equal for each seed.
- The digests equal the table above.

A digest mismatch with `ALL PASS` means an output differs between builds, but each build still agrees with its own
reference. To find it, bisect by section: comment out the `test_*` calls in `main` and compare the digests again.

## 4. Sanitizers (about 5 min; scale 1 only)

```bash
$SAN --tool memcheck ./build/k2test_gpu 20261108 1 | tail -5   # ERROR SUMMARY: 0 errors
$SAN --tool synccheck ./build/k2test_gpu 20261108 1 | tail -3                     # 0 errors (warp collectives)
$SAN --tool initcheck ./build/k2test_gpu 20261108 1 | tail -3                     # 0 errors
$SAN --tool racecheck ./build/k2test_gpu 20261108 1 | tail -3                     # 0 hazards (no shared memory used)
```

## 5. Benchmark (about 2 min)

```bash
./build/k2bench_gpu 1 3 | tee build/bench_gpu.txt
./build/k2bench_cpu 8 1 | tee build/bench_cpu_box.txt      # CPU build on the box's CPU at 1/8 size, for scale
```

The workloads:
- 32 items, each with 8 samples and 4 raters on one 512 x 512 frame: EDT of 384 masks, boundary F and HD95 of
  1,024 pairs, and GED and HM-IoU of 32 sets.
- Typical validity: 4,096 one-ring shapes of 64 to 200 vertices. It runs once with the grid off and once with it on.
- Adversarial validity: 64 shapes of 10 rings with 1,000 vertices per ring. Half of them are touching strips, the
  ShapeTotals straggler.

Each line prints the old path against the new one.

## 6. Optional profile (about 3 min, only if time is left)

```bash
ncu --target-processes all -k regex:"warp_kernel|kernel" --launch-skip 0 --launch-count 40 \
    --section SpeedOfLight --section Occupancy --section MemoryWorkloadAnalysis \
    -o build/k2 ./build/k2bench_gpu 8 1
```

## 7. Copy back, then stop the instance

```bash
rsync -av BOX:~/K2/build/{ptxas.txt,test_gpu_s1.txt,test_gpu_s4.txt,bench_gpu.txt,bench_cpu_box.txt,k2.ncu-rep} ./K2_results/
```

Record the measured values in DESIGN.md, section 9, marked [measured].
