# K1 polygonizer: GPU run sheet (one RTX 5090, sm_120)

Who runs this: the main agent, on one rented RTX 5090. Use Vast AI, North America hosts only, and run the preflight
before renting. Budget: at most 20 minutes of GPU time (about $0.15 at the $0.440/h that R3 paid for an RTX 5090
[memory: r3-gpu-test]). The $3 cap is far away. Nothing here trains or downloads models.

What it proves:
1. Exactness: the K1 GPU output is identical to the reference polygonizer, `gpucheck::check_masks` and the 37,093 real
   LIDC masks.
2. Memory safety: compute-sanitizer finds no errors.
3. Speed: K1 on the GPU against `check_masks` on the GPU, same inputs, each result checked for equality.

## 0. Before renting (on the Mac, about 2 min, no cost)

Run inside the clone. `REPO` is the repo root.

```bash
REPO=$(git rev-parse --show-toplevel)
K1=$REPO/docs/implementation2.0/kernels/K1
cd $K1
./build_cpu.sh $REPO/cuda                                    # CPU builds (clang++)
for b in test_k1_cpu test_k1_cpu_o1 test_k1_cpu_o2; do ./build/cpu/$b 3000 1 | tail -1; done   # 3 x PASS
# the 37,093 real LIDC crops as raw files (also cross-checks K1 against the numpy polygonize; about 10 s)
$REPO/.venv/bin/python -I cross_numpy.py $K1/build/cpu/k1_dump_cpu $REPO/data/k1_bundle/lidc \
    --masks $REPO/data/checker_eval/masks.npz                # -> "cross_numpy: 37093 masks, 38002 rings, 0 mismatches -> PASS"
./nvcc_docker_check.sh $REPO/cuda                            # optional: nvcc 13.0 compile + link for sm_120, no GPU (~5 min)
# bundle: sources, scripts, gpucheck, raw masks (~46 MB)
mkdir -p /tmp/k1_bundle/gc /tmp/k1_bundle/lidc
cp $K1/*.cpp $K1/*.h $K1/*.sh /tmp/k1_bundle/ && cp $REPO/cuda/gpucheck.h $REPO/cuda/gpucheck.cpp /tmp/k1_bundle/gc/
cp $REPO/data/k1_bundle/lidc/in.bits $REPO/data/k1_bundle/lidc/in.hw $REPO/data/k1_bundle/lidc/in.start /tmp/k1_bundle/lidc/
tar -C /tmp -czf /tmp/k1_bundle.tgz k1_bundle
```

`data/checker_eval/masks.npz` is LIDC-IDRI only (CC BY 3.0, its manifest says "no RIGA"), so the bundle may leave the
Mac. `masks_gen.h` holds the synthetic masks.

## 1. On the box: toolkit and upload (about 3 min)

The image needs `nvcc` that knows sm_120. That means CUDA 12.8 or newer. The CUDA major version must not be newer
than the driver's (`nvidia-smi` prints the highest CUDA version it supports).

```bash
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv && nvidia-smi | head -4
nvcc --version || true
# Without nvcc (or with nvcc < 12.8), take the toolkit run_box.sh uses (CUDA 13.0.3; needs a driver for CUDA 13):
#   pip install "cuda-toolkit[nvcc,cccl]==13.0.3"
#   export CUDA_HOME=$(python -c "import nvidia, os; print(os.path.join(nvidia.__path__[0], 'cu13'))")
#   [ -e $CUDA_HOME/lib64 ] || ln -s lib $CUDA_HOME/lib64; export NVCC=$CUDA_HOME/bin/nvcc
# With a CUDA 12.x driver, use "cuda-toolkit[nvcc,cccl]==12.8.*" instead (cu12 folder).
# Upload from the Mac: scp -P <port> /tmp/k1_bundle.tgz root@<host>:/root/ ; then on the box:
cd /root && tar xzf k1_bundle.tgz && cd k1_bundle
```

## 2. Build (about 2 min)

```bash
./build_gpu.sh ./gc sm_120 2>&1 | tail -3        # -> "build_gpu: OK (sm_120)"
```

## 3. Exactness on the GPU (about 2 min)

```bash
./build/gpu/test_k1_gpu 3000 21 2>&1 | tee test_gpu.txt | tail -4
./build/gpu/test_k1_gpu 3000 22 2>&1 | tail -1
```

Pass: both last lines start with `PASS:` and report 0 failures. Each case runs K1 on the GPU 15 times (seed rows 0, 1,
2, 3 and 16; one chunk, 700-pixel chunks, and a 1-byte scratch limit that splits the batch down to one mask). Each
result is compared field by field with `ref_polygonize` (a literal port of `polygonize.py`) and with `check_masks`
(its CPU build). The GPU's thread schedule is the real test of order independence: union-find, atomicMin and atomicOr.
The Mac runs the same items forward, reversed and permuted.

```bash
./build/gpu/k1_dump_gpu lidc/in out16.bin 16 2>&1 | tee dump_gpu.txt
./build/gpu/k1_dump_gpu lidc/in out0.bin 0 2>&1 | tail -1
```

Pass: `37093 masks, 38002 rings, 1476508 corners` and `gpucheck::check_masks IDENTICAL`, for both seed-row settings.
These are the counts the Mac measured on the same file.

## 4. Memory checks (about 3 min, capped)

compute-sanitizer comes with the toolkit (pip: add `nvidia-cuda-sanitizer-api`, as run_box.sh does). The NVIDIA
catalog skill `compute-sanitizer-diagnose-cuda-bugs` covers this tool. It is recommended only, not installed.

```bash
timeout 300 compute-sanitizer --tool memcheck --leak-check full ./build/gpu/test_k1_gpu 25 5 2>&1 | tail -3 | tee san_mem.txt
timeout 300 compute-sanitizer --tool initcheck ./build/gpu/test_k1_gpu 10 6 2>&1 | tail -3 | tee san_init.txt
```

Pass: `ERROR SUMMARY: 0 errors` from both. racecheck and synccheck do not apply: K1 uses no shared memory and no
barriers. The cross-thread writes are atomics or happen in separate kernels.

## 5. Speed (about 2 min)

```bash
./build/gpu/bench_k1_gpu 10 lidc/in 2>&1 | tee bench_gpu.txt
```

There is one line per workload: the median wall time per call (end to end, from host bytes in to host results out),
masks/s, the speedup over `check_masks` (GPU build), `IDENTICAL`, and the K1 stage split from CUDA events (h2d, ccl,
rings, trace, d2h). The workloads:

| line | what | why |
|---|---|---|
| lesion512 b256 (seed rows 16, 0) | 256 frames of 512x512, 2 blobs each | the RL/eval frame size |
| noise512 b16 | 50 % noise, about 18 k rings per mask | many rings, many runs |
| comb512 b4 (16, 0) | one ring of about 262 k edges per mask | seed rows against one thread per ring |
| spiral / concentric / checker 512 b4 | deep union-find chains, 512 nested rings, 130 k rings per mask | adversarial |
| real (all) | the 37,093 LIDC crops in one call | production-like crops |

Record these. Pass/fail is only the `IDENTICAL` column. The numbers go into DESIGN.md "Measured".

## 6. Optional profile (at most 3 min, only if nsys exists)

```bash
timeout 180 nsys profile -o k1_prof --stats=true ./build/gpu/bench_k1_gpu 2 > nsys_stats.txt 2>&1 || true
grep -A25 "cuda_gpu_kern_sum" nsys_stats.txt | head -40
```

## 7. Bring back and stop

```bash
# from the Mac
scp -P <port> root@<host>:/root/k1_bundle/{test_gpu.txt,dump_gpu.txt,san_mem.txt,san_init.txt,bench_gpu.txt,nsys_stats.txt} $K1/results/
```

Destroy the instance right after the copy. Time: about 3 + 2 + 2 + 3 + 2 (+ 3) = 12 to 15 minutes.

## If something fails

- Build error: the docker check compiled this exact code with nvcc 13.0.88 for sm_120 (step 0). Check the toolkit
  version and flags first.
- `k1poly: invariant violated on the device`: a trace hit its step cap or a missing successor, a hole lacked its key,
  or a ring walk did not close. Rerun the failing case with `seed_rows=0` (`k1_dump_gpu ... 0`) to separate the seeds
  from the CCL. Then run memcheck on that case only.
- Mismatch only on the GPU, never on the Mac: suspect a missing ordering between kernels or a non-atomic
  read-modify-write. Every cross-thread write must be an atomic or happen in a different kernel (DESIGN.md
  "exactness").
