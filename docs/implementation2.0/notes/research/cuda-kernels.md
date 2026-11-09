# Research note: CUDA kernels for the geometry checker (phase D)

Date 2026-11-09. Scope: plan phase D (`Implementation.md` D.4, lines 9348-10256), kernel drafts K1/K2/K3
(`docs/implementation2.0/kernels/*/DESIGN.md`), `cuda/gpucheck.cpp`, `CMakeLists.txt`. Pins: CUDA 13.0
(PyPI `cuda-toolkit[nvcc,cccl,nvtx]==13.0.3` = nvcc 13.0.88, CCCL 3.0.1), `nvidia-cuda-sanitizer-api==13.4.92`,
sm_120 (RTX 5090) and sm_86, scikit-build-core 1.1.1, nanobind 3.1.0, CMake >= 3.26.
Tags: [measured: ...] = run in this job (docker on the owner's Mac, linux/arm64, no GPU); [web: URL] = read at that
URL; [code: path:line] = repo; [unverified] = not checked.

## 0. Skills

- `nvidia-skills:nvidia-skill-finder` loaded. Live catalog (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`,
  group "GPU Development", fetched 2026-11-09) is the same as `notes/nv_catalog.json`. Only relevant entry:
  `compute-sanitizer-diagnose-cuda-bugs`. Recommended, not installed (install needs owner approval:
  `npx skills add nvidia/skills --skill compute-sanitizer-diagnose-cuda-bugs --agent claude-code --global --yes`).
  No catalog skill covers CCL, contour tracing, EDT, NVTX or nsys.
- Key guidance of that skill [web: github.com/NVIDIA/skills/blob/main/skills/compute-sanitizer-diagnose-cuda-bugs/SKILL.md]:
  1. Build with `-lineinfo` (add `-Xcompiler -rdynamic` for host backtraces).
  2. Unknown bug class: memcheck first. racecheck, initcheck and synccheck assume a memcheck-clean program.
  3. racecheck/initcheck/synccheck need a vanilla binary; a `-fdevice-sanitize=memcheck` binary corrupts their results.
  4. Use `--error-exitcode N` in CI (the tool exits 0 on errors otherwise).
  5. A clean run proves only that tool's bug class. Global-memory races are out of scope (racecheck = shared memory).
  6. Imprecise error without file:line: rerun with `--force-blocking-launches`.
  7. Check `compute-sanitizer --version` against the release notes before calling a report a false positive.

## 1. Toolchain facts at the pins (D.4 "Common toolchain", D1, D2, D8)

| Fact | Source |
|---|---|
| CUDA 13.0 removed offline compilation for Maxwell, Pascal, Volta. Supported real/virtual targets run sm_75 .. sm_121a. Default `-arch` is sm_75. | [web: docs.nvidia.com/cuda/archive/13.0.2/cuda-toolkit-release-notes/index.html.md] [web: docs.nvidia.com/cuda/archive/13.0.2/cuda-compiler-driver-nvcc/index.html] |
| CUDA 13.x needs driver >= 580 (13.0 U2 Linux: >= 580.95.05). Preflight must check the Vast host driver. | release notes above |
| CUDA 13.0 removed deprecated `cudaDeviceProp` fields (`clockRate`, `memoryClockRate`, `deviceOverlap`, `computeMode`), deprecated `double4`/`long4`/`ulong4`/`longlong4`/`ulonglong4` (use `*_16a`/`*_32a`), removed `this_multi_grid`/`multi_grid_group` and `cudaLaunchCooperativeKernelMultiDevice`. K1's `device_info()` reads only `name`, `multiProcessorCount`, `l2CacheSize`, `major`, `minor`, `totalGlobalMem` [code: K1/k1poly.cpp:944-950] and compiles on 13.0 [measured]; new code must not add the removed fields. | release notes above |
| Family targets: "PTX for .target sm_XYf can be compiled to GPU targets sm_XZ, sm_XZf, sm_XZa where Z >= Y" (same family); `sm_XYa` PTX only to `sm_XYa`. `compute_120f` covers CC 12.0 and 12.1. Plain `sm_120` is enough for the box; `sm_120f` buys nothing for this code. | nvcc doc above; [web: docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/compute-capabilities.html] |
| `-arch=sm_100` == `--gpu-architecture=compute_100 --gpu-code=sm_100,compute_100` (SASS + PTX). | nvcc doc above |
| CMake `CUDA_ARCHITECTURES` "If no suffix is given then code is generated for both real and virtual architectures"; `native` since 3.24; env `CUDAARCHS` seeds `CMAKE_CUDA_ARCHITECTURES` on first configure only (3.20+). | [web: cmake.org/cmake/help/latest/prop_tgt/CUDA_ARCHITECTURES.html] [web: cmake.org/cmake/help/latest/envvar/CUDAARCHS.html] |
| CC 12.x and 8.6: 64 K regs/SM, 255 regs/thread, 100 KB smem/SM, 99 KB/block (>48 KB needs dynamic smem + opt-in), 48 resident warps/SM, unified L1/smem 128 KB. | compute-capabilities page above (guide v13.5) |
| nvcc 13.0.88, sm_120, sm_86 and sm_120f all compile `__ballot_sync`, `__reduce_add_sync`, `__reduce_min_sync`, nvtx3 and CUB `DeviceScan`. `gpucheck.cpp` compiles with `--fmad=false -x cu` for sm_120 and sm_86. | [measured: docker linux/arm64, python:3.12-slim, 2026-11-09] |
| CCCL in the wheel: `CCCL_VERSION 3000001`, `CUB_VERSION 300001` (3.0.1, not "13.0.85"; that is the wheel version). | [measured: `cu13/include/cccl/cuda/std/__cccl/version.h:17`] |
| D2 step 6 answered: `cub::DeviceScan::ExclusiveSum` is `template <typename InputIteratorT, typename OutputIteratorT, typename NumItemsT>` with `using OffsetT = detail::choose_offset_t<NumItemsT>`; `choose_offset` picks `uint32_t` when `sizeof(NumItemsT) <= 4`, else `unsigned long long`. So `i64` counts are legal. With `long long` the sm_120 scan kernel used 48 regs/31,632 B smem; with `int` 71 regs/33,808 B (K1 log). The `i32` output sums still need the `INT_MAX` guard. | [measured: `cub/device/device_scan.cuh:181-193` in the wheel] [web: raw.githubusercontent.com/NVIDIA/cccl/v3.0.1/cub/cub/detail/choose_offset.cuh:53-62] |
| D1 CMake check answered: with `CUDACXX` = the wheel's nvcc, `find_package(CUDAToolkit)` gives `CUDAToolkit_INCLUDE_DIRS=.../nvidia/cu13/include`, `CUDAToolkit_VERSION=13.0.88`; `find_file(... nvtx3/nvtx3.hpp PATHS ${CUDAToolkit_INCLUDE_DIRS} NO_DEFAULT_PATH)` finds `.../cu13/include/nvtx3/nvtx3.hpp`; `CUDA::nvtx3` has that include dir and `INTERFACE_LINK_LIBRARIES=dl`. Without the `nvtx` extra the compile fails with `nvtx3/nvtx3.hpp: No such file or directory`. | [measured: CMake 3.31.6, arm64 container] [web: cmake.org/cmake/help/latest/module/FindCUDAToolkit.html: `CUDA::nvtx3` added 3.25; `nvToolsExt` gone since CUDA 12.9] |
| `nvidia-cuda-sanitizer-api==13.4.92` installs `cu13/bin/compute-sanitizer` = "Version 2026.3.0.0", plus `TreeLauncherSubreaper`, `TreeLauncherTargetLdPreloadHelper`, `libsanitizer-collection.so`. | [measured] |
| `cuobjdump` is not in `cuda-toolkit[nvcc]`; it is `nvidia-cuda-cuobjdump` and needs `nvidia-cuda-nvdisasm` (no 13.0.x release on PyPI; 13.2+ only). Fine for reading SASS, not part of the build. | [measured; web: pypi.org/pypi/nvidia-cuda-cuobjdump/json] |
| nsys is not on PyPI (`nvidia-nsys-cli` 0.0.1.dev5 is "A fake package to warn the user"). Install `nsight-systems-cli` from `developer.download.nvidia.com/devtools/repos/<ubuntuXXXX>/<arch>/`. The box base is Debian trixie, not Ubuntu: repo fit [unverified]. | [web: pypi.org/pypi/nvidia-nsys-cli/json] [web: docs.nvidia.com/nsight-systems/InstallationGuide/index.html] |

## 2. Build flags (D1, D2, D8; CMakeLists.txt)

- Floating point, nvcc device code: `--fmad` default true; `--prec-div`, `--prec-sqrt` default true (IEEE RN);
  `--ftz` default false; `--use_fast_math` implies `--ftz=true --prec-div=false --prec-sqrt=false --fmad=true`
  [web: nvcc doc above]. Measured on sm_120 for `o = a*x + c` (double): default emits 1 `DFMA`; `--fmad=false` emits
  `DMUL` + add, 0 `DFMA` [measured: cuobjdump -sass]. Explicit `fma()`/`fmaf()`/`__fma_rn` stay fused under
  `--fmad=false` [web: cuda-programming-guide/05-appendices/mathematical-functions.html "Fused Multiply-Add Support"].
- **Finding: the CUDA build's host code is not covered by `-ffp-contract=off`.** `CMakeLists.txt:21-22` adds
  `-ffp-contract=off` only under `$<COMPILE_LANGUAGE:CXX>`, but `cuda/gpucheck.cpp` is `LANGUAGE CUDA` when
  `OBJDEL_CUDA=ON` [code: CMakeLists.txt:21-22, 29], so its host functions go through nvcc -> g++ with no contraction
  flag. GCC's default is "-ffp-contract=off for C in a standards compliant mode ..., -ffp-contract=fast otherwise"
  [web: gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html], i.e. `fast` for C++. Measured: nvcc 13.0.88 + g++ 14.2 on
  aarch64 contracts a host `a*x + c` into `fmadd`; with `-Xcompiler -ffp-contract=off` it does not [measured: objdump].
  On the x86-64 box the baseline ISA has no FMA, so it is latent until someone adds `-march=native`/`-mfma`, or builds
  on aarch64 (GH200/Grace, the Mac docker). Fix: add `$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-ffp-contract=off>` next
  to `--fmad=false`. K2/K3 Makefiles already pass `-Xcompiler -ffp-contract=off` [code: K3/Makefile].
- Other nvcc flags: `-lineinfo` (needed for sanitizer file:line, no perf cost), `-res-usage` = `--resource-usage`
  "implies --ptxas-options=--verbose" (registers/smem per kernel), `-t 0` = threads = CPU count, `-Werror all-warnings`
  [web: nvcc doc]. `-std` accepts up to `c++20`; gpucheck stays C++17 [code: CMakeLists.txt:17].
- nanobind flags never reach nvcc because the core is a static library linked into `nanobind_add_module`
  [code: CMakeLists.txt:10-11, 40-42]. Keep `NOMINSIZE` (nanobind's `-Os` halves kernel speed, per the repo comment).
- scikit-build-core 1.1.1 config: `-Ccmake.define.X=Y` (env `SKBUILD_CMAKE_DEFINE`), `-Cbuild-dir=...`
  (env `SKBUILD_BUILD_DIR`), `-Ccmake.build-type=...`; lists join with `;`. Prefer `cmake.define` over `cmake.args`
  [web: scikit-build-core.readthedocs.io/en/latest/configuration/index.html]. D2 step 7's
  `-C "cmake.define.OBJDEL_TEST_HOOKS=GC_POISON_SCAN;GC_ALIAS_SCAN_FOR_TEST"` matches this. For CUDA archs pass
  `-Ccmake.define.CMAKE_CUDA_ARCHITECTURES=120` or env `CUDAARCHS=120` (first configure only; a reused
  `build-dir` keeps the cached value) [web: CUDAARCHS page].
- nanobind 3.x device arrays (if a later task passes torch tensors in place of host copies): `nb::ndarray<T,
  nb::device::cuda, nb::c_contig, nb::shape<...>>`, `.data()`, `.device_id()`. DLPack `__dlpack__(stream=...)`
  "performs no synchronization": the caller must order the producer stream before the engine's lane stream
  [web: nanobind.readthedocs.io/en/latest/ndarray.html]. Today every entry takes host numpy, so this is future work.

## 3. Compile-only checks on the Mac (D1 step 3, D8 `make gpu-check DRY=1`, `cuda/docker_check.sh`)

- **Finding: run the docker checks as linux/arm64, not linux/amd64.** PyPI has manylinux aarch64 wheels for
  `nvidia-cuda-nvcc 13.0.88`, `nvidia-cuda-cccl 13.0.85`, `nvidia-nvtx 13.0.85`, `nvidia-cuda-sanitizer-api 13.4.92`
  [web: pypi.org JSON of each]. The pinned base `python:3.12-slim-trixie@sha256:a6e34c59...` is a multi-arch index
  that includes `linux/arm64/v8` [measured: `docker buildx imagetools inspect`]. Only `--platform` changes.
- Same device code on both hosts: K1's `k1poly.cpp` for sm_120 on arm64 gives the same `ptxas -v` register count for
  all 16 kernels as the amd64 log (Header 14, RingWalk 30, Trace 40, SeedEmit 30, SetCounts 12, SeedFlag 16, RingInit 25,
  Popc 10, RingMark 22, HoleKey 36, Flatten 14, Merge 26, RunEmit 40, RunCount 36, Quads 37, pack_kernel 22)
  [measured vs K1/results/nvcc_docker_sm120.log.txt]. Bit-identical cubins across hosts: [unverified].
- Speed: the K1 nvcc compile plus gpucheck for sm_120 and sm_86 took 9 s natively [measured]. K1's script says
  "Several minutes under emulation" [code: K1/nvcc_docker_check.sh:4].
- Host-side caveat: an arm64 container's g++ contracts FMAs (section 2), so the CPU tests inside an arm64 container
  are a stricter check of the `-ffp-contract=off` rule than the amd64 one.
- Alternative base: `nvidia/cuda:13.0.3-devel-ubuntu24.04` is amd64+arm64, about 3.7 GB for arm64 [web: Docker Hub tag
  API]. Its nvcc build may differ from 13.0.88 [unverified]; the PyPI route matches the box exactly, so keep it.

## 4. Bit-exactness and determinism rules (D-G6, D-G9, D-G10; all kernels)

1. Integer atomics (`atomicAdd`, `atomicMin`, `atomicOr`) commute: final values do not depend on order. Float atomics
   do not: "The ordering of arithmetic atomic operations depends on hardware scheduling which can vary between runs"
   [web: mathematical-functions.html "Deterministic Results"]. Never accumulate a float with atomics.
2. Slot order from `atomic_inc` is schedule-dependent (e.g. `Gather` fills `c.bg` by `atomic_inc(c.ncomp + o)`
   [code: gpucheck.cpp:1402]). It is safe only because `ObjTerms` sorts the list (`sort_descending`, :1455, :1493).
   Any new list built this way (K2 `Directed` d² lists) must be sorted or read by order-free functions before output.
   D-G6's 20-run byte check is the only GPU test for this: racecheck does not see global memory (section 9).
3. CUB: `DeviceScan` float results "may vary from run to run"; `DeviceReduce` is run-to-run deterministic for float
   only "on the same GPU device" and "may be inconsistent ... of a different compute-capability"
   [web: raw.githubusercontent.com/NVIDIA/cccl/v3.0.1/cub/cub/device/device_reduce.cuh:111-115, device_scan.cuh:107-109].
   So CPU = sm_86 = sm_120 needs integer scans/reductions only. K1/K2 use integer `ExclusiveSum` only.
4. Warp collectives on integers are exact in any lane order: `__reduce_add_sync`/`__reduce_min_sync` take 32-bit
   `int`/`unsigned` only (CC >= 8.0; compiled to `REDUX.SUM.S`/`REDUX.MIN.S` on sm_120 [measured]). K2's `wsum`/`wmin`
   on `i64` must stay as `__shfl_xor_sync` loops [code: K2/k2score.inc:74-91].
5. Float math by one thread, fixed order, `--fmad=false` + host `-ffp-contract=off` (section 2). `sqrt(double)` is
   0 ULP (correctly rounded); `sqrtf` is 0 ULP only without `--use_fast_math` [web: mathematical-functions.html table].
6. Never `--use_fast_math`, never `-ftz=true` anywhere in this target.

## 5. Connected-component labelling (D2 run CCL, D4b/D5 pixel union-find)

- Union with `atomicMin`, as in `gpucheck.cpp:1267-1282` and K1 `unite`: Allegretti et al. give
  `Find(L,i)`: follow `L[i]` until a root; `Union(L,a,b)`: loop { a = Find(a); b = Find(b); if a == b stop; if a > b
  swap; `old = atomicMin(&L[b], a)`; done = (old == b); b = old } (their version stores label+1)
  [web: prittt.ing.unimore.it/media/publications/pdfs/2018ipas_gpuccl.pdf, Alg. 1-4]. Playne and Hawick's
  contribution: skip the atomic when no other thread interfered, and decide locally which neighbours need a merge
  [same paper, §III; D. P. Playne, K. Hawick, IEEE TPDS 29(6):1217-1230, 2018, doi:10.1109/TPDS.2018.2799216].
- Pitfalls:
  - Parent loads in `Find` must not be cached across the loop: `volatile` (as gpucheck and K1 do) or `cuda::atomic_ref`
    relaxed loads. A plain load can be hoisted by the compiler and spin forever on a stale root [unverified on 13.0,
    standard memory-model reasoning].
  - Path compression during merge races with `atomicMin`; only compress in a separate flatten kernel. D5 step 4's
    path halving `atomicMin(par[x], par[par[x]])` is safe because parents only decrease (min-root union keeps
    `par[x] <= x`) [reasoning; matches the plan].
  - Root = minimum index makes the final labels independent of schedule (K1 DESIGN §4 "Order independence").
- Block-based alternatives (D-R2 fallback): BUF and BKE label 2x2 blocks because "any two foreground pixels in a 2 x 2
  block are connected" under 8-connectivity; BKE "is always the best-performing algorithm on 2D datasets" and is the
  one in OpenCV CUDA; BUF is usually best in 3D (relevant to the 3D-voxel direction). Run-based "performance
  significantly drops when foreground objects ... have elongated shapes on the y-direction" — K1's known risk
  (§10.1). [web: federicobolelli.it/media/publications/pdfs/2024tpds.pdf §2, §5; Allegretti, Bolelli, Grana, IEEE TPDS
  31(2):423-438, 2020]. BKE/BUF are 8-connected fg only; the engine also needs bg 4-connectivity, so K1's choice stands.
- Two labellers on purpose (D5 step 3, D-G16): the certificate's b0/b1 leg must come from the pixel chain, not from K1's
  run CCL. Nothing in the literature changes that.

## 6. Contour / ring tracing (D2)

K1-draft (CCL-derived ring keys + seeded crack tracing) has its own proofs and citations [code: K1/DESIGN.md §4, §11].
Not re-verified here: Suzuki-Abe 1985, Garcia-Molla 2022, G-Contour (ICCAD 2025). One citation could not be found:
"TRUCO, arXiv:2610.37297" returned nothing on web search (K1 also notes HF has no record) [unverified]. Port rules
that matter: keep polygonize's "first east edge" ring order (D2 step 8), keep `K1_ORDER` variants (D2 step 7), and keep
`pack_kernel`'s early exit warp-uniform: `__ballot_sync(0xffffffff, ...)` is valid only if every non-exited lane in the
mask reaches it; "Exited threads are ignored" [web: cpp-language-extensions.html "Warp __sync Intrinsic Constraints"].
K1's `if (w >= nw) return;` is uniform per warp because `w = t >> 5` [code: K1/k1poly.cpp:378-388].

## 7. Exact EDT (D7 `edt`, D7b tables, D-G10)

- scipy 1.18.1 (box and Mac) [measured: `.venv` and box lock :762]: `distance_transform_edt` computes the feature
  transform with Maurer's exact algorithm (`ni_morphology.c:1008` comment "Exact euclidean feature transform, as
  described in: C. R. Maurer") [web: raw.githubusercontent.com/scipy/scipy/v1.18.1/scipy/ndimage/src/ni_morphology.c],
  then `dt = (ft - indices).astype(float64)`, squares, `np.add.reduce(axis=0)`, `np.sqrt`
  [code: .venv scipy/ndimage/_morphology.py:2584-2600]. So scipy's value is `sqrt64(exact integer d²)` (exact while
  d² < 2^53). `cache.edt` then casts to float32 [code: src/objdel/reward/cache.py:130-135].
- Bit-exact device recipe: integer d² (min over sites; ties do not matter because only the value is returned), then
  `__double2float_rn(sqrt((double)d2))` or `(float)sqrt((double)d2)`. Do not use `sqrtf((float)d2)` unless you prove
  equality: float(d2) is exact for d2 < 2^24, and double rounding of sqrt from binary64 to binary32 is believed
  innocuous (Figueroa 1995) [unverified here]; the double route needs no proof.
- FH Algorithm 1 (1-D lower envelope): `s = ((f(q)+q²) - (f(v[k])+v[k]²)) / (2q - 2v[k])`; `if s <= z[k]` pop, else
  push `v[k]=q, z[k]=s, z[k+1]=+inf`; fill: `while z[k+1] < q: k++`; `D(q) = (q-v[k])² + f(v[k])`
  [web: cs.brown.edu/people/pfelzens/papers/dt-final.pdf, Theory of Computing 8:415-428, 2012, Alg. 1]. With float `s`
  and `f = +inf` for non-sites this is not exact (`inf - inf`). Use Meijster's integer separator
  `Sep(i,u) = (u² - i² + g(u)² - g(i)²) div (2(u - i))` with a finite sentinel g = H + W, as K2 does
  [code: K2/k2score.inc:456, 499, 510-512] [paper: Meijster, Roerdink, Hesselink 2000, doi:10.1007/0-306-47025-x_36;
  formula from K2 and memory, PDF not reachable: unverified wording].
- Pitfall: the numerator can be negative (g(u) < g(i)), and C++ `/` truncates toward zero. `div` must be floor
  division; K2 uses `floordiv<i64>` [code: K2/k2score.inc:512]. Keep it in the port; a plain `/` is off by one on some
  columns and only an adversarial test sees it.
- Pitfall: scipy has no defined "no zero pixel" output. `distance_transform_edt(np.ones((3,3), bool))` returns
  `[[1, 1.414, 2.236], [2, 2.236, 2.828], [3, 3.162, 3.606]]` [measured]. The D7 spec returns `_BIG` there; exclude
  that case from scipy-equality tests and test it against `_BIG`.
- K2 row pass: one warp per row, nearest site via `clz`/`ctz` on packed words (ballot idea of Manduhu and Jones, IEEE
  TIP 2019, doi:10.1109/TIP.2019.2916741) [code: K2/DESIGN.md §2.3; not re-verified].

## 8. Warp primitives and cooperative groups (D6 `ShapeTotalsWarp`, K1 pack, raster)

- Signatures: `unsigned __ballot_sync(unsigned mask, int pred)`, `int __any_sync/__all_sync(...)`,
  `T __shfl_sync(unsigned mask, T v, int srcLane, int width=warpSize)`, `__shfl_xor_sync(mask, v, laneMask, width)`,
  `T __reduce_{add,min,max}_sync(unsigned mask, T v)` (T 32-bit int/unsigned, CC >= 8.0), `__syncwarp(mask)`
  [web: docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/cpp-language-extensions.html §5.4.6].
- Rules (same page): each calling thread must be in `mask`; each non-calling, non-exited thread must not be; all
  named threads must reach the same call with the same mask, "In conditional code, all conditions must evaluate
  identically across all non-exited threads specified in the mask"; disjoint masks may run concurrently; violation =
  hang or UB. `__reduce_*_sync` "do not provide any memory ordering": use `__syncwarp()` before reading peers' global
  writes (K2 `wsync`) [same page].
- K2's `wballot`/`wsum` assume warp-uniform control flow up to the call [code: K2/k2score.inc:65]. Lane loops
  `for (i = lane; i < n; i += 32)` must call collectives outside the loop or with a uniform trip count (round n up to
  a multiple of 32 and guard the body, not the call). synccheck reports violations (section 9).
- Cooperative groups (optional): `cg::tiled_partition<32>(cg::this_thread_block())`, `cg::reduce(tile, v,
  cg::plus<int>())` (operators plus, less, greater, bit_and/or/xor; hardware accelerated on CC >= 8.0)
  [web: cuda-programming-guide/04-special-topics/cooperative-groups.html §4.4.6.1]. Not needed: raw intrinsics already
  work on both builds via the `HD` wrappers, and see the synccheck known issue below.

## 9. compute-sanitizer (D8 `sanitize`, D-G5, D5c leak test)

CLI at the box version (2026.3.0 from the 13.4.92 wheel) [measured `--help`; web:
docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html, doc v2026.4.0]:
- `--tool memcheck|racecheck|initcheck|synccheck` (default memcheck); `--error-exitcode N` (default 0);
  `--print-limit N` (default 100, 0 = no limit); `--leak-check full|no` (default no; reports at context destroy:
  exit, `cudaDeviceReset`, `cuCtxDestroy`); `--target-processes all` (default all); `--racecheck-report
  hazard|analysis|all` (default analysis); `--racecheck-detect-level info|warn|error` (default warn);
  `--track-unused-memory` (initcheck); `--track-stream-ordered-races all|use-before-alloc|use-after-free|no`
  (default all); `--kernel-name kns=<substr>|kne=<exact>|regex=<re>` (note `regex=`, not `regex:`);
  `-s/--launch-skip`, `-c/--launch-count`; `--padding N`; `--force-blocking-launches`; `--nvtx yes` (default);
  `--xml --save f`, `--log-file f` (`%p`, `%q{ENV}`); `--suppressions f` (CUDA API errors, initcheck, racecheck only);
  `--require-cuda-init yes` (default: error exit when the app makes no CUDA call).
- The plan's command form is right: `timeout 300 "$CUDA_HOME/bin/compute-sanitizer" --tool T --error-exitcode 99
  --print-limit 20 "$PY" -I -m pytest -q $D_GPU_TESTS -m "gpu and sanitizer"` with `PYTORCH_NO_CUDA_MEMORY_CACHING=1`.
  Add `-p no:cacheprovider`. `--require-cuda-init=yes` makes a run where every test skipped (CPU build) fail: good.
- Do not use `-fdevice-sanitize=memcheck`: it needs "The PTX compiler and compute-sanitizer ... from the same toolkit
  version" (nvcc 13.0.88 vs sanitizer 13.4.92 here), and it corrupts the other three tools [web: sanitizer doc; skill].
- Coverage gaps that matter for this engine (design facts, not tool bugs):
  1. memcheck sees allocation bounds only. The engine carves many arrays out of one `Block` with `Layout::put`
     (256-byte aligned offsets) [code: gpucheck.cpp:188-196], so an overrun from one array into the next one in the
     same arena is invisible, and `--padding` pads only whole allocations. Mitigation: a test-hook build
     (`OBJDEL_TEST_HOOKS`) where `Layout` adds a guard gap filled with a canary and the host checks it after `wait()`
     [proposal].
  2. initcheck tracks bytes ever written. Arenas are reused and never cleared [code: gpucheck.cpp:135-150], so a read of
     stale data from the previous call passes initcheck. The D5 `zero(changed, ...)` rule is still needed for the
     first call. Run each sanitizer test in a fresh process with the smallest corpus first, so the first call hits
     fresh memory.
  3. racecheck checks only shared memory: "Currently, this tool only supports detecting accesses to on-chip shared
     memory". The only shared buffer is `raster_kernel`'s `bufs` [code: gpucheck.cpp:418]; K1/K2 have none. Global
     races are covered only by D-G6 (20 runs byte-identical) and the `K1_ORDER` variants.
  4. synccheck covers `__syncthreads`/`__syncwarp` divergence and bad masks: relevant to raster, K1 `pack_kernel` and K2
     warp collectives.
- Known issues in the release notes [web: docs.nvidia.com/compute-sanitizer/ReleaseNotes/index.html]:
  - "The synccheck tool may incorrectly emit divergence in block errors in applications compiled with the 13.0 CUDA
    toolkit" when they use cooperative-group synchronizations. Our build is 13.0: do not introduce CG syncs, or expect
    false synccheck failures (another reason to keep raw `__syncwarp`).
  - "initcheck ... might report false positives for device-to-host cudaMemcpy operations on padded structs" initialized
    by a kernel; workaround `#pragma pack`. `ObjOut` is 2 x i64 + 4 x double = 48 B, no padding [code:
    gpucheck.cpp:1199-1202]; adding an `i32` field (or `RingInfo`-style mixed structs) creates padding. Keep D2H
    structs 8-byte homogeneous or `static_assert(sizeof(T) == sum of fields)`.
  - Driver 615.00-615.70 with apps built before 13.4 may print "CUDA API Error: Fabric handle support is dependent on
    IMEX channels". Our build is 13.0: if the rented host has such a driver, suppress that one API error or pick
    another host.
  - 2026.4 fixed "a device memory leak with cudaMallocAsync/cudaFreeAsync under initcheck"; we run 2026.3, so initcheck
    runs with the D5c pool may grow memory; keep the corpus small.
- Leaks: the plan correctly drops `--leak-check full` (engine never destroyed) and uses `live_device_bytes` instead.

## 10. NVTX (D1 step 1, Phase A `OBJDEL_NVTX`)

- C++ (header-only, `#include <nvtx3/nvtx3.hpp>`): `nvtx3::scoped_range r{"name"};` (alias of
  `scoped_range_in<>` global domain); domain: `struct gc_domain { static constexpr char const* name{"gpucheck"}; };`
  then `nvtx3::scoped_range_in<gc_domain>`; `NVTX3_FUNC_RANGE()`; `nvtx3::mark("msg")`;
  `nvtx3::registered_string_in<D>` for cheap repeated messages; `unique_range` only when start/end cross threads
  [web: nvidia.github.io/NVTX/doxygen-cpp/]. `#define NVTX_DISABLE` before the include compiles everything out;
  with no tool attached each call is "a handful of instructions"; `-ldl` needed only for glibc < 2.34 (trixie has
  2.41 [measured]) [web: github.com/NVIDIA/NVTX/blob/release-v3/README.md].
- Recommendation for D1: use a domain (`gpucheck`) so nsys can filter with `--nvtx-domain-include gpucheck`, and
  registered strings (or `NSYS_NVTX_PROFILER_REGISTER_ONLY=0`, below). The plan's `GC_RANGE(name)` macro with
  `__LINE__`-suffixed names is fine; make it `nvtx3::scoped_range_in<gc_domain>`.
- Python: `nvtx` (PyPI 0.2.17, aarch64 + x86_64 wheels): `nvtx.annotate(message=None, color=None, domain=None,
  category=None, payload=None)` (decorator/context manager), `push_range`/`pop_range(domain=None)`,
  `start_range(...) -> (int, int)`/`end_range(id)`, `get_domain(name)`; global-function messages are cached as
  registered strings [web: nvidia.github.io/NVTX/python/reference.html].
- torch (box 2.13.0 and Mac 2.14.1 identical): `torch.cuda.nvtx.range_push(msg)`, `range_pop()`,
  `range_start(msg) -> int`, `range_end(id)`, `mark(msg)`, `range(msg)` (context manager)
  [web: raw.githubusercontent.com/pytorch/pytorch/v2.13.0/torch/cuda/nvtx.py:27-116]. Pitfall: on a build without
  `torch._C._nvtx` (Mac CPU torch) every call raises `RuntimeError("NVTX functions not installed. Are you sure you
  have a CUDA build?")` [same file :8-22]. `OBJDEL_NVTX` code must check `torch.cuda.is_available()` or use the `nvtx`
  package. `torch.cuda.profiler.start()/stop()` call `cudaProfilerStart/Stop` [same tag, profiler.py:20-36].

## 11. nsys (D1 step 4, optional)

- `nsys profile -t cuda,nvtx,osrt` (default trace is `cuda,opengl,nvtx,osrt`); `-c cudaProfilerApi` (start/stop by
  `torch.cuda.profiler.start/stop`) or `-c nvtx -p "RANGE@DOMAIN"` (`range`, `range@domain`, `range@*`);
  `--capture-range-end stop|stop-shutdown(default)|repeat[:N]|none`; `-o path` (`%p`, `%h`, `%q{ENV}`); `-f true`;
  `--cuda-memory-usage true` (costly); `-s none --cpuctxsw none` in containers; `--stats true`
  [web: docs.nvidia.com/nsight-systems/UserGuide/index.html].
- NVTX capture matches registered strings only, unless `-e NSYS_NVTX_PROFILER_REGISTER_ONLY=0`: "If you do not use
  registered strings you will have to enable the full match by launching your application with
  NSYS_NVTX_PROFILER_REGISTER_ONLY=0" [same guide]. `nvtx3::scoped_range{"literal"}` is not registered.
- Containers: CPU sampling needs host `perf_event_paranoid <= 2` and `perf_event_open` (`--cap-add=SYS_ADMIN`);
  check with `nsys status --environment` [same guide]. Vast containers: use `-s none` [unverified per host].
- Reports: `nsys stats --report nvtx_sum,nvtx_gpu_proj_sum,cuda_gpu_kern_sum,cuda_api_sum --format csv` (names from
  [web: docs.nvidia.com/nsight-systems/AnalysisGuide/index.html]; `--report/--format` flag spelling [unverified],
  check `nsys stats --help`). `nvtx_gpu_proj_sum` projects NVTX ranges onto GPU time: the per-stage GPU numbers D1
  wants without extra events.
- The plan's CUDA-event stage timing (`timing=True`) stays the gate source; nsys is optional (D1 step 4).

## 12. Per-task checklist

| Task | Do | Section |
|---|---|---|
| D1 | `find_file` + `CUDA::nvtx3` work with the wheel layout (measured); add `nvtx` extra; NVTX domain + registered strings; torch nvtx raises on Mac; nsys from apt, `NSYS_NVTX_PROFILER_REGISTER_ONLY=0` | 1, 10, 11 |
| D1/D8 docker | `--platform linux/arm64` for `docker_check.sh`, K1/K3 docker scripts (same digest) | 3 |
| all | add `-Xcompiler=-ffp-contract=off` for the CUDA language; never fast-math | 2, 4 |
| D2 | `NumItemsT` may be i64; keep INT_MAX guard on i32 sums; volatile parents; uniform early exit before ballot | 1, 5, 6 |
| D4b/D5 | min-root `atomicMin` union; path halving only with `atomicMin`; label lists sorted before output | 4, 5 |
| D6 | `__reduce_add_sync` only for 32-bit; i64 via shuffles; uniform control flow to collectives | 4, 8 |
| D7/D7b | integer d², floor `div`, double sqrt then float32; all-ones input = `_BIG`, not scipy | 7 |
| D8 | sanitizer flags as planned + `-p no:cacheprovider`; no `-fdevice-sanitize`; arena overruns invisible to memcheck; initcheck padded-struct and CG-synccheck known issues; driver 615.00-615.70 caveat | 9 |
| box | driver >= 580 for CUDA 13; no removed `cudaDeviceProp` fields in new code | 1 |
