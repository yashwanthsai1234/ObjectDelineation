# CUDA C++ kernel optimization guidance applied to a custom polygon validation / rasterization kernel (nanobind, many small polygons)

Research date: 2026-11-06. Context: cuda/gpucheck.cpp is exact but measured 0.002-0.19x CPU speed on an RTX 4090 for batches of small LIDC nodule polygons (rasters ~64-512 px). Search tools returned page summaries (not raw pages) in several cases; items from forum posts are flagged as anecdotal.

## 1. CUDA C++ Best Practices Guide essentials (APOD, coalescing, shared memory, occupancy, divergence, atomics, launch config, instruction tips)

### Takeaway
The Best Practices Guide's high-priority rules are: minimize host-device transfers, coalesce global memory access, keep enough warps resident (blocks multiple of 32, 128-256 threads, many more blocks than SMs), and avoid intra-warp divergence. Each applies directly to a one-thread-per-pixel or one-block-per-polygon rasterizer, and the guide itself says a transfer-bound kernel should not be judged by kernel time alone.

### Cited Findings
- APOD (Assess, Parallelize, Optimize, Deploy) is described as a cyclical process in which early speedups can be deployed before further tuning begins (guide section 2.2 per the fetched text). — [CUDA C++ Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Sec 10.1 (High Priority): "Minimize data transfer between the host and the device." The guide contrasts 898 GB/s device memory (Tesla V100) with 16 GB/s PCIe x16 Gen3, and recommends minimizing transfers even when that means running kernels that do not beat the CPU. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Sec 10.2.1 (High Priority): "Ensure global memory accesses are coalesced whenever possible." For compute capability 6.0+, a warp's accesses are serviced as as many 32-byte transactions as needed; adjacent aligned 4-byte words need four; misaligned sequential access touches five 32-byte segments; stride 2 gives 50% load efficiency. Block sizes that are multiples of 32 help alignment. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Sec 10.2.3.1: shared memory has 32 banks (cc 5.x+), 32 bits per bank per clock; successive 32-bit words map to successive banks; conflicting addresses in one request serialize; same-address reads in a warp are broadcast without penalty. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Sec 11.1 Occupancy: ratio of active warps per SM to maximum; higher occupancy does not always mean higher performance, but low occupancy hurts latency hiding. On cc 7.0 (65,536 registers, 2048 threads/SM) full occupancy needs at most 32 registers/thread. Register allocation granularity matters (guide example: 37 regs/thread gives 75% occupancy with 128-thread blocks, 63% with 320-thread blocks). Tools: `--ptxas-options=-v`, Occupancy Calculator, `cudaOccupancyMaxActiveBlocksPerMultiprocessor`; `__launch_bounds__` for register control. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Sec 11.3 Execution configuration: threads per block a multiple of 32, minimum 64, 128-256 a good starting range; grid should have more blocks than multiprocessors (ideally thousands); a larger block size does not imply higher occupancy; several smaller blocks can beat one large block for kernels that call `__syncthreads()` often. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Sec 11.4: shared memory per block can limit occupancy. Sec 11.5: concurrent kernels need non-default streams. Sec 11.6: avoid multiple contexts per GPU in one app. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Sec 12.1.2 / 13.1 Divergence: divergent threads in a warp execute paths serially; write thread-ID-dependent conditions aligned to warp boundaries (e.g. `threadIdx / warpSize`); short branches may be predicated; since Volta, `__syncwarp()` can reconverge. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Sec 12.1 Instruction tips: prefer single precision; integer division/modulo can compile to up to 20 instructions (use shifts/masks for powers of two); use `rsqrtf()`; use `f` suffix on float literals; use `__syncthreads()` sparingly; signed loop counters; `-use_fast_math` trades accuracy for speed. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Atomics: warp aggregation (one lane performs the atomic for the whole warp) cuts atomic count by up to 32x; the same NVIDIA post says shared-memory atomics improved the filtering example but stayed within 1.5x of the original, and a later post says NVCC now warp-aggregates automatically in many cases and its output beat hand-written aggregation. — [NVIDIA blog, warp-aggregated atomics](https://developer.nvidia.com/blog/cuda-pro-tip-optimized-filtering-warp-aggregated-atomics/) (blog text as summarized by search; compiler claim is from a later post, same search result set)
- Forum anecdote (not authoritative): a developer found plain atomicAdd slightly faster than a hand-optimized warp-primitive version, with a reply suggesting the compiler already aggregates; advice was to check SASS. — [NVIDIA forum](https://forums.developer.nvidia.com/t/warp-aggregate-atomicadd/37130)
- The Best Practices Guide text I could read did not cover atomics (its sections 11-13 as fetched had no atomics section, and the last ~38k characters were not read). — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)

### Inferences
- A one-thread-per-pixel even-odd test over a small polygon is naturally coalesced along rows, but only if the raster is stored row-major with one thread per consecutive x; per-polygon launches of 64x64 = 4096 threads (32 blocks of 128) cannot fill 128 SMs, so "more blocks than SMs" fails unless polygons are batched into one launch.
- Edge loops with data-dependent early exits (crossing-number with branches on y-range) diverge within a warp when threads hold different pixels in the same row; branch-free formulations (compute crossing as a predicated add) are the standard fix.
- Per-polygon statistic accumulation via global atomics on one counter per polygon would serialize; use warp/block reduction then one atomic per block.

### Gaps
- Best Practices Guide section numbers are taken from the fetch tool's rendering of the current HTML; the final ~38k characters (sections after 13.x, possibly an atomics or compute-capability appendix) were not read. The guide's bandwidth numbers (V100, PCIe Gen3) are from that guide revision; it does not give RTX 4090 figures.
- I did not read the CUDA C++ Programming Guide directly this session (atomic function throughput, cooperative groups chapter).

## 2. Host-device transfer costs (PCIe 4.0 x16 / RTX 4090, pageable vs pinned, async, zero-copy, batching, per-call overhead)

### Takeaway
Each cudaMemcpy and each kernel launch carries a fixed cost on the order of several to tens of microseconds, so for rasters of 64-512 px (4 KB-256 KB at 1 byte/px) the fixed overhead, not bandwidth, dominates. Batch everything into one pinned buffer, one H2D, one launch, one D2H.

### Cited Findings
- PCIe 4.0 x16 theoretical ~32 GB/s per direction; measured pinned throughput on 4090-class systems reported at ~24-26 GB/s (forum reports), pageable ~6 GB/s on one user's system (explained as extra host-side staging copy). Anecdotal, system-dependent; a motherboard M.2 slot can drop the GPU slot to x8. — [NVIDIA forum](https://forums.developer.nvidia.com/t/bandwidth-test-of-pageable-memory-is-mush-different-in-2-computer/355919); [TechPowerUp RTX 4090 PCIe scaling](https://techpowerup.com/review/nvidia-geforce-rtx-4090-pci-express-scaling)
- Best Practices Guide 10.1.1: pinned (page-locked) memory attains highest host-device bandwidth (~12 GB/s on PCIe x16 Gen3 in the guide); do not overuse as it is scarce. 10.1.2: async copies require pinned host memory and non-default streams; overlap depends on `asyncEngineCount`. 10.1.3: zero-copy is low priority for discrete GPUs; mapped pinned memory should be read/written once, coalesced. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Batching many small transfers into one larger transfer "performs much better"; `cudaMemcpy2D/3D` help for strided layouts. — [NVIDIA blog, How to Optimize Data Transfers in CUDA C/C++ (Dec 2012)](https://developer.nvidia.com/blog/how-optimize-data-transfers-cuda-cc/)
- The same 2012 blog measured pinned vs pageable at 16 MB: GTX 680 desktop 6.19 vs 5.37 GB/s H2D (~15%), laptop NVS 4200M ~2.5x gain. Old hardware; shows gain depends on host system. — [same blog](https://developer.nvidia.com/blog/how-optimize-data-transfers-cuda-cc/)
- Kernel launch overhead (NVIDIA, Sept 2019, V100, CUDA 10.1): a 2.9 us kernel took 9.6 us per kernel when synchronizing after every launch (~6.7 us overhead), 3.8 us per kernel with overlapped launches, 3.4 us with a CUDA Graph (one cudaGraphLaunch per 20 kernels); graph creation+instantiate ~400 us; first launch ~33% slower. — [NVIDIA blog, Getting Started with CUDA Graphs](https://developer.nvidia.com/blog/cuda-graphs/)
- Per-call overhead anecdotes: smallest cuMemcpyDtoH ~18 us (old forum thread); ~5-7 us latency per cudaMemcpyAsync submission; forum respondents say tiny copies take ~10 us+; reported 25-67 us idle gaps after sync flag copies. Treat as indicative, driver/GPU dependent. — [forum](https://forums.developer.nvidia.com/t/cumemcpy-overhead-cumemcpy-call-has-unacceptable-overhead/4844); [forum](https://forums.developer.nvidia.com/t/delay-before-kernel-runs/304251); [forum](https://forums.developer.nvidia.com/t/massive-hiccups-when-transferring-flags-back-to-host-challenging-the-commonly-quoted-2-10us-latency/20503); [emergentmind summary](https://www.emergentmind.com/open-problems/breakdown-of-cudamemcpyasync-overhead)
- NVIDIA blog (via search) defines launch latency as API call to kernel start, ~20 us on V100 in its example; GPU-side gaps typically <1 us. — [NVIDIA blog](https://developer.nvidia.com/blog/?p=20916)

### Inferences
- Arithmetic: a 256x256 byte mask is 64 KB; at ~25 GB/s that is ~2.6 us of pure transfer, comparable to or smaller than the fixed ~5-20 us call cost. So per-polygon copies are overhead-bound by roughly an order of magnitude; one copy of N concatenated polygons amortizes it by N.
- If gpucheck.cpp calls cudaMalloc/cudaFree per call (not seen, I did not read the source), those are far more expensive than memcpy and are a prime suspect for 0.002x; the fix is a persistent device arena/stream-ordered allocator (`cudaMallocAsync`) and pinned staging buffers allocated once.
- Pageable numpy buffers passed through nanobind incur the hidden staging copy; allocate pinned output/input once and reuse.

### Gaps
- No controlled RTX 4090 pinned/pageable bandwidth or per-call latency benchmark from NVIDIA or a paper was found; figures are forum anecdotes. Recommend measuring with `bandwidthTest`/Nsight Systems on the actual host.
- I did not find a CUDA 12/13-specific published launch-latency number; the NVIDIA figure is V100/CUDA 10.1 (2019).

## 3. When small problems are faster on the CPU; amortization by batching

### Takeaway
No published formal break-even analysis for polygon workloads was found, but NVIDIA guidance (transfer minimization, graphs for short kernels, grid >> SM count) together with the fixed-cost numbers above explain the loss: work per call is below the cost of the call. The standard fix is a single launch that covers all polygons in the batch, one block (or warp) per polygon.

### Cited Findings
- Guide sec 10.1 explicitly accepts running kernels that are not faster than CPU to avoid transfers, i.e. it treats data movement as the dominant cost to avoid. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- CUDA Graphs cut per-kernel cost for short kernels from 9.6 to 3.4 us (2.8x) in NVIDIA's example, by launching 20 kernels through one CPU call. — [NVIDIA blog](https://developer.nvidia.com/blog/cuda-graphs/)
- Guide sec 11.3: the grid should have more blocks than multiprocessors, ideally thousands. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- cuSpatial vendor blog claims "up to 10,000x" over CPU and 30 billion points/s on GV100 for `point_in_polygon` (baseline unclear; large point counts); a later cuSpatial post reports 140x vs GeoPandas at 1e6 points with throughput saturating at about 1000 points. — [RAPIDS cuSpatial blog](https://medium.com/rapids-ai/releasing-cuspatial-to-accelerate-geospatial-and-spatiotemporal-processing-b686d8b32a9); [cuSpatial 22.12](https://medium.com/rapids-ai/good-things-come-in-twos-binary-predicate-support-and-more-in-cuspatial-22-12-1d89292e581d)

### Inferences
- GPU wins reported (e.g. 140x) are at 1e6 items; the user's 0.002-0.19x ratios plausibly correspond to batches where total pixels per launch are too few (or per-polygon launches) to cover launch+copy+sync cost, and where the CPU path is already a tight C++/numpy loop on KB-sized rasters that live in L1/L2.
- A back-of-envelope break-even: CPU time per 256x256 polygon raster is likely tens to hundreds of us; GPU fixed cost per synchronous round trip is ~20-60 us (launch + H2D + D2H + sync). GPU wins only when batch size B satisfies B * cpu_time > fixed + B * (gpu_marginal_time + transfer_time). This is my arithmetic, not a published result.
- Python/nanobind boundary adds per-call cost too (argument conversion, GIL, numpy allocation); do one call per batch, return one big array.

### Gaps
- No published break-even for small-polygon GPU vs CPU rasterization found; recommend an own sweep over batch size (1 to 10^5 polygons) and raster size, reporting launch/copy/kernel split from Nsight Systems.

## 4. GPU polygon rasterization / point-in-polygon algorithms and measured speedups

### Takeaway
Published GPU rasterizers win on large workloads (36x-290x) and use tile binning and edge pruning; the closest literature to ours is compute-shader tile renderers (Vello) and parity/winding formulations. No source benchmarked many tiny polygons on tiny rasters, which is the regime where we lose.

### Cited Findings
- Zhou et al., "A Parallel Scheme for Large-scale Polygon Rasterization on CUDA-enabled GPUs", Transactions in GIS, 2017: multi-GPU + MPI + OpenMP, ~1.3e8 polygons (50.32 GB), conversion time reduced from 25.43 h to 0.69 h (speedup ~36.9). — [abstract via Bohrium](https://www.bohrium.com/paper-details/a-parallel-scheme-for-large-scale-polygon-rasterization-on-cuda-enabled-gpus/814525717449539584-193)
- Siemens GPU rasterizer for computational lithography (SC25 poster / vendor blog, 2026): on H100, up to 290x for Manhattan and 45x for curvilinear geometry vs optimized CPU, <1% absolute error; computes fractional coverage; the poster's numbers are internally inconsistent (0.29 s vs 27.5 s in one table, 0.9 s vs 30 s in another), and only a poster/blog was found. — [SC25 poster](https://sc25.supercomputing.org/proceedings/posters/poster_files/post228s2-file2.pdf); [Siemens blog](https://blogs.sw.siemens.com/calibre/2026/04/13/why-a-gpu-rasterizer-matters-for-computational-lithography-in-both-performance-and-precision-at-scale/)
- NVIDIA Path Rendering (Kilgard and Bolz, SIGGRAPH Asia 2012): "stencil, then cover" two-step; stencil pass resolves per-sample coverage (even-odd or nonzero fill via the stencil buffer), cover pass draws conservative geometry; shipped in NV_path_rendering since 2011. No CPU speedup figure was retrieved in my search. — [NVIDIA GPU Accelerated Path Rendering](https://developer.nvidia.com/gpu-accelerated-path-rendering); [Kilgard annex](https://ghostscript.com/~robin/nvpr_annex.pdf)
- Vello (linebender): compute-shader 2D renderer with pipeline stages path counting, binning, tile allocation, backdrop, coarse, path tiling, fine (area/MSAA variants); stated rationale that fine rasterization stays on GPU since CPU-side raster would incur upload cost. Newer direction (RustWeek 2025 slides): sparse strips, vello_cpu and vello_hybrid; original design needs a modern GPU and has high unpredictable memory use. — [Vello shaders (docs.rs)](https://docs.rs/crate/vello/latest/source/src/shaders.rs); [Vello roadmap 2023](https://github.com/linebender/vello/blob/main/doc/roadmap_2023.md); [RustWeek 2025 slides](https://docs.google.com/presentation/d/1f_vKBJMaD68ifBO2j83lBly9Zdk-2bsvj_DIHXxvcuk/htmlpresent)
- AresDB geo-intersection: one thread per (point, edge) pair, atomic XOR to flip parity bit for crossings; odd count = inside. Design note, no benchmark. — [AresDB wiki](https://github.com/uber/aresdb/wiki/Geo-Intersection)
- GPUDay 2024 talk on polygon masks: winding-number contributions and segment distances can be computed in parallel over edges; only edges to the right of the point contribute, so per-tile edge lists are short. No break-even figure. — [GPUDay 2024 slides](https://indico.wigner.hu/event/1567/contributions/4024/attachments/2688/5442/DBerenyi_GPUDay2024_rzrlr.pdf)
- Crossing-number and winding-number tests differ for self-overlapping polygons (even-odd vs nonzero); Hormann-Agathos: angle-summing is costly, fast winding-number variant comparable to crossing count. — [Wikipedia, Point in polygon](https://en.wikipedia.org/wiki/Point_in_polygon); [Haines, Graphics Gems IV draft](https://erich.realtimerendering.com/ptinpoly/)
- Older GPU quadtree study (Zhang et al., 2014): up to 27x vs single-core, 2x vs multi-core CPU for polygon query modules. — [hgpu.org](https://hgpu.org/?p=12530)

### Inferences
- For our size (tens to a few hundred edges per polygon, 64-512 px), one-thread-per-pixel with a loop over edges in registers/shared memory (edge list loaded once per block into shared memory, broadcast reads are conflict-free) is the simplest; one-thread-per-edge needs atomics or reductions per pixel and wastes lanes at these edge counts. Scanline (one thread or warp per row, sort/span crossings, fill spans) reduces work from O(pixels x edges) to O(rows x crossings + pixels) and gives coalesced writes.
- Even-odd and nonzero can share one loop: accumulate signed crossing (winding) per pixel; even-odd = parity of crossing count, nonzero = winding != 0; cost is the same.
- Tile binning (Vello/NVPR style) pays off for large canvases with many paths, not for 512 px single-polygon rasters; it adds multi-kernel overhead we cannot afford at this size.

### Gaps
- No peer-reviewed paper with numbers for GPU point-in-polygon on many tiny polygons/rasters; all speedups above are for large inputs. No Vello CPU-vs-GPU per-scene numbers retrieved. Kilgard-Bolz paper figures not read.

## 5. Reduction patterns for per-polygon statistics (warp shuffle, cooperative groups, CUB)

### Takeaway
For per-polygon counts/areas, reduce within a warp or block (shuffle / CUB BlockReduce / cooperative groups), then write one value per polygon; use DeviceSegmentedReduce when per-polygon data lives in one flat array with offsets. Limited sourcing obtained.

### Cited Findings
- CUB `DeviceSegmentedReduce` takes begin/end offset iterators; for contiguous segments a single offsets array of length num_segments+1 serves both (end offsets = offsets+1); outputs must not overlap input or offset arrays; num_segments is 64-bit in current docs. — [CCCL docs, DeviceSegmentedReduce](https://nvidia.github.io/cccl/cub/api/structcub_1_1DeviceSegmentedReduce.html) (as reported by search; the direct fetch returned 404)
- Warp aggregation: shuffle-reduce within a warp then one atomic per warp reduces atomics up to 32x; Cooperative Groups and warp primitives (CUDA 9+) are the supported APIs. — [NVIDIA blog](https://developer.nvidia.com/blog/cuda-pro-tip-optimized-filtering-warp-aggregated-atomics/)

### Inferences
- With one block per polygon, `cub::BlockReduce` (or warp shuffle + one shared-memory stage) gives per-polygon counts with no global atomics and no second kernel; DeviceSegmentedReduce is a second pass with its own temp-storage query and launch, so only worth it if polygons are processed with different thread mappings.

### Gaps
- Could not confirm CUB's internal strategy for small segments or fixed-size-segment variants, nor any benchmark of BlockReduce vs DeviceSegmentedReduce for small segments; the CCCL source/docs pages for `cub/block` and `cub/warp` were not read.

## 6. CUDA tooling for kernel tuning (Nsight Compute, compute-sanitizer)

### Takeaway
Use Nsight Systems first to see whether time is in launch/copy/sync gaps (our suspected cause), then Nsight Compute per kernel for Launch Statistics, Occupancy, Memory Workload and Warp State sections.

### Cited Findings
- Nsight Compute Launch Statistics: grid/block size and resources required; Occupancy: achieved vs theoretical active warps, large gap usually indicates imbalanced workload; Memory Workload Analysis: memory unit utilization, bandwidth, cache hit rates; Warp State Statistics: cycles per instruction and stall states (focus on stalls only if schedulers do not issue every cycle); Source Counters: branch efficiency and sampled stall reasons. — [Nsight Compute Profiling Guide](https://docs.nvidia.com/nsight-compute/ProfilingGuide/index.html)
- Nsight Compute: first profiled kernel in each context has a one-time overhead; tool overhead is excluded from reported duration, so host timers/CUDA events do not match profiled durations. — [Nsight Compute Profiling Guide](https://docs.nvidia.com/nsight-compute/ProfilingGuide/index.html)
- Best Practices Guide recommends `--ptxas-options=-v` for registers and the Occupancy Calculator. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)

### Inferences
- Nsight Compute only measures kernel duration; our measured slowdown includes copies, allocation, sync and Python overhead, which Nsight Systems (timeline) exposes. Start there. Time with CUDA events around kernel only vs wall-clock around the whole nanobind call to split overhead from compute.
- Warm up before timing (first-launch/context creation costs; graph first launch ~33% slower in NVIDIA's example) and report steady-state.

### Gaps
- compute-sanitizer (memcheck/racecheck/synccheck/initcheck) documentation was not fetched this session; no cited facts. The remaining ~137k characters of the Nsight Compute guide (metric reference, occupancy section details) were not read.

## Likely causes and standard fixes for 0.002-0.19x (summary of inferences; verify with Nsight Systems)
1. Fixed per-call costs (launch ~5-20 us, memcpy ~5-20 us each, sync) exceed per-polygon CPU time: batch all polygons into one flat buffer with offsets, one H2D, one launch, one D2H, one sync.
2. Per-polygon or per-call device allocation: use a persistent arena / `cudaMallocAsync`; reuse pinned staging buffers (cudaMallocHost) across calls.
3. Underfilled GPU: grid below SM count (128 SMs on a 4090 is my own background knowledge, not sourced here); assign one block (or warp) per polygon in a 1-D grid of thousands of blocks.
4. Pageable host memory and large masks returned to host: pinned buffers, return bit-packed or per-polygon summary stats instead of full masks when only validation results are needed (reduce D2H volume).
5. Divergence and atomics in the inner loop: branch-free crossing accumulation; block-level reductions instead of global atomics.
6. Many small kernels per chain (build, check, score): fuse stages or capture in a CUDA Graph (2.8x per-kernel in NVIDIA's example).
