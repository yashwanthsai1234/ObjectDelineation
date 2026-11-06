# CUDA host-side design for a small Python-called library (persistent arena, pinned staging, packed output, multi-process)

Scope note: official sources fetched were NVIDIA blogs and API reference pages (blog posts date from CUDA 11.x/Pascal-Volta era; the runtime API reference is current). Several programming-guide chapters could not be fetched (landing page only), so some items are in Gaps. Items under "Inferences" that are not backed by a cited finding are engineering judgment and should be benchmarked on the 4090/5090.

## 1. Cost of cudaMalloc/cudaFree/managed/memset/sync; unified-memory faults vs explicit copy; prefetch/advise

### Takeaway
cudaFree on a non-stream-ordered pointer may implicitly synchronize the device, and "any CUDA API call may block or synchronize", so 20 managed allocs + free per call is both slow and serializing. On-demand page migration ran at roughly half of explicit-copy bandwidth in NVIDIA's own tests; prefetch matches memcpy but is not a free lunch. Explicit copies into a persistent device buffer are the straightforward replacement.

### Cited Findings
- `cudaFree` on memory still used by a running kernel "synchronizes the device before freeing the memory"; `cudaFreeAsync` needs no such sync. — [NVIDIA blog: stream-ordered allocator part 1](https://developer.nvidia.com/blog/using-cuda-stream-ordered-memory-allocator-part-1/)
- Runtime API reference: for pointers not from cudaMallocAsync/cudaMallocFromPoolAsync, "this API may perform implicit synchronization"; cudaFree "may also return error codes from previous, asynchronous launches." — [CUDA Runtime API: Memory Management](https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html)
- "Any CUDA API call may block or synchronize for various reasons"; this is subject to change and undocumented behavior should not be relied on. Kernel launches are async w.r.t. host. — [CUDA Runtime API: API synchronization behavior](https://docs.nvidia.com/cuda/cuda-runtime-api/api-sync-behavior.html)
- cudaMemset is async w.r.t. host except when target is pinned host memory; the Async versions are always async. Synchronous memcpy forms go through the default stream. — [API synchronization behavior](https://docs.nvidia.com/cuda/cuda-runtime-api/api-sync-behavior.html)
- Unified memory faults (Pascal/Volta): a fault locks the TLBs of affected SMs, new translations stall; one page drew 114 duplicate faults; a 64 MB profile showed ~5.98 ms of fault-group time vs ~2.57 ms of H2D transfer. — [NVIDIA blog: Maximizing Unified Memory Performance](https://developer.nvidia.com/blog/maximizing-unified-memory-performance-cuda/)
- H2D over PCIe (P100/V100, 16 GB/s peak): on-demand migration ~5.4 GB/s, cudaMemPrefetchAsync ~10.9 GB/s, cudaMemcpyAsync ~11.4 GB/s. Using one warp per 64 KB page cut faults and sped on-demand migration up to 2x, to within ~30% of explicit copy. — [Maximizing Unified Memory Performance](https://developer.nvidia.com/blog/maximizing-unified-memory-performance-cuda/)
- cudaMemPrefetchAsync updates page-table mappings (can limit concurrency); calls on a busy stream are deferred to a background thread, on an idle stream usually not; D2H prefetch blocks the CPU until complete. The article gives no concrete cudaMemAdvise recommendations. — [Maximizing Unified Memory Performance](https://developer.nvidia.com/blog/maximizing-unified-memory-performance-cuda/)
- Batching many small transfers into one larger one performs significantly better, even if packing into a contiguous buffer is needed; keep intermediates in device memory. — [CUDA C++ Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- cudaMalloc memory is aligned to at least 256 bytes. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)

### Inferences
- Replacing ~20 managed allocations with one arena sub-allocated by offset (align each sub-buffer to 256 B, or 128 B for coalescing) removes both the alloc cost and the cudaFree implicit sync; cudaMemsetAsync on the arena (or a zeroing kernel) replaces per-buffer cudaMemset.
- Packing all inputs (polygon vertices, offsets, params) into one contiguous pinned staging buffer yields a single H2D copy, consistent with the batching guidance.
- Managed memory plus prefetch could work, but the cited data show prefetch only equals memcpy bandwidth while adding page-table work; there is no cited advantage over explicit copies for a fixed per-call upload/download pattern.
- The numeric ratios above are from P100/V100-era hardware; on RTX 4090 (PCIe 4.0) / 5090 (PCIe 5.0) the absolute numbers differ and the fault-path costs on Ada/Blackwell consumer GPUs were not found.

### Gaps
- No per-call latency numbers (microseconds) for cudaMalloc, cudaMallocManaged, cudaFree, cudaDeviceSynchronize on Ada/Blackwell were found; must be measured (suggest a small microbenchmark).
- No cited guidance on cudaMemAdvise settings (SetReadMostly, SetPreferredLocation) beyond "hints can change migration behavior".
- Best Practices guide section 10.3 (allocation) was truncated in the fetch, so any explicit "avoid cudaMalloc in hot loops" statement is unverified here.

## 2. Stream-ordered allocation (cudaMallocAsync/pools) vs hand-rolled grow-only arena

### Takeaway
cudaMallocAsync/cudaFreeAsync (CUDA 11.2+) are stream-ordered, do not implicitly sync, and are fast when the pool's release threshold is raised from the default 0. For one library with a predictable per-call footprint, a single grow-only arena is simpler and equally fast; the pool is the better choice if lifetimes or sizes vary a lot.

### Cited Findings
- CUDA 11.2 added cudaMallocAsync/cudaFreeAsync; allocation/free errors (e.g., OOM) are returned immediately, "cannot fail asynchronously". — [NVIDIA blog: stream-ordered allocator part 1](https://developer.nvidia.com/blog/using-cuda-stream-ordered-memory-allocator-part-1/)
- Default pool per device via cudaDeviceGetDefaultMemPool; default `cudaMemPoolAttrReleaseThreshold` is 0, so unused pool memory is returned to the OS at every synchronization; setting it high (example UINT64_MAX) keeps memory across syncs and makes cudaMallocAsync cost independent of allocation size. The threshold is a hint; the driver can still release pool memory to satisfy unrelated allocations such as cudaMalloc in the same process. — [stream-ordered allocator part 1](https://developer.nvidia.com/blog/using-cuda-stream-ordered-memory-allocator-part-1/)
- `cudaMemPoolTrimTo` releases unused memory above bytesToKeep; reuse policies (follow-event-dependencies, opportunistic, internal-dependencies) can be toggled per pool; opportunistic reuse can add nondeterminism. — [stream-ordered allocator part 1](https://developer.nvidia.com/blog/using-cuda-stream-ordered-memory-allocator-part-1/)
- For cudaMallocAsync pointers, cudaFree performs no implicit sync and callers must ensure accesses completed; cudaFreeAsync is recommended for best reuse. — [Runtime API: Memory Management](https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html)

### Inferences
- Because the default threshold of 0 trims the pool at every sync, and the planned design syncs once per call, using the pool without setting the threshold would re-allocate every call. If the pool is used at all, set the threshold to a large value (or UINT64_MAX) at init.
- Recommended pattern (judgment): allocate one arena with cudaMalloc at first use sized to the high-water mark (e.g., batch 1024 x 512x512 rasters), grow-only with a geometric factor and a sync before reallocation (only on growth), bump-allocate sub-buffers per call, reset the offset each call. With one stream and one sync per call there are no lifetime hazards.
- Alternative of cudaMallocAsync of the whole per-call scratch plus a high release threshold is nearly equivalent and avoids arena bookkeeping; it fits better if multiple streams are added later.

### Gaps
- Part 2 of the stream-ordered allocator blog (benchmark numbers) was not fetched; no absolute allocation latencies were found.
- Behavior of the pool under cudaMalloc/cudaMallocAsync mixing and fragmentation on 32 GB (5090) was not found.

## 3. Pinned host memory: allocation vs registration, cost, safe amount, write-combined, async copy requirements

### Takeaway
cudaMemcpyAsync is only truly async with pinned host memory; pinning is expensive (allocation ~ms, free tens of ms in forum reports), so allocate one pinned staging buffer once per process and reuse it. cudaHostRegister on numpy buffers pins per call and should be avoided in the hot path.

### Cited Findings
- cudaMemcpyAsync "requires pinned host memory" to be non-blocking; overlap also needs non-default streams and asyncEngineCount > 0. — [Best Practices Guide 10.1.2](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Async copies between device and pageable host "might be synchronous with respect to host"; the driver "may synchronize with the stream and stage the copy into pinned memory." Sync H2D from pageable does a stream sync before the copy; D2H returns only after the copy completes. — [API synchronization behavior](https://docs.nvidia.com/cuda/cuda-runtime-api/api-sync-behavior.html)
- Pageable memory is copied by the driver into a temporary pinned buffer first; pinned transfers on a GTX 680 desktop reached ~6.2/6.7 GB/s vs ~5.4/5.6 GB/s pageable (laptop: 5.8/6.0 vs 2.3 GB/s). — [NVIDIA blog: How to Optimize Data Transfers](https://developer.nvidia.com/blog/how-optimize-data-transfers-cuda-cc/)
- "You should not over-allocate pinned memory"; it takes physical RAM from the OS and other programs; the safe limit is hard to predict, test on target and check return status. Pinned memory "should not be overused" and is heavyweight to allocate. — [How to Optimize Data Transfers](https://developer.nvidia.com/blog/how-optimize-data-transfers-cuda-cc/); [Best Practices Guide 10.1.1](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- Pinned transfers on PCIe x16 Gen3 reach roughly 12 GB/s of 16 GB/s theoretical. — [Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html)
- cudaHostRegister page-locks an existing range; default flags on UVA systems give mapped and portable; must be unregistered with cudaHostUnregister using the same base pointer; on systems where pageableMemoryAccessUsesHostPageTables is true the range may not actually be page-locked; supported only if cudaDevAttrHostRegisterSupported. — [Runtime API: Memory Management](https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html)
- cudaHostAllocWriteCombined can cross PCIe faster on some systems but "cannot be read efficiently by most CPUs"; suited to buffers the CPU writes and the device reads. cudaMallocHost: use sparingly for staging buffers. — [Runtime API: Memory Management](https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html)
- Forum anecdotes (not official): cudaMallocHost ~1.3 ms vs 0.1 ms for malloc in one report; cudaFreeHost ~18 ms on an A30; cudaHostRegister time grows with buffer size; advice is to pin once and reuse. — [NVIDIA forums: cudaMallocHost slow](https://forums.developer.nvidia.com/t/why-is-cudamallochost-so-slow/1540); [cudaFreeHost 20x slower](https://forums.developer.nvidia.com/t/cudafreehost-consistently-20x-slower-than-free-cudafree-full-runnable-example-code-available/217706); [cudaHostRegister time linear?](https://forums.developer.nvidia.com/t/cudahostregister-time-linear-in-buffer-size/57912)

### Inferences
- Use one pinned input staging buffer and one pinned output staging buffer (output is a 1-bit mask, so small), sized to the arena high-water mark, allocated once with cudaHostAlloc (flags: default or WriteCombined for the input buffer only if the CPU never reads it back; do not use WC for the output buffer).
- Memory budget: input and output staging for 1024 polygons are likely MBs; with 1 GPU process the pinned footprint is tiny relative to system RAM, so the "how much is safe" concern is negligible unless per-process buffers are multiplied by many workers.
- Writing from numpy into pinned memory costs one host memcpy; registering numpy arrays per call (page-locking each time) is likely slower for small payloads. Registration is only attractive for large, long-lived buffers.
- No cited number puts a cap on safe pinned memory; the usual rule of thumb (a fraction of RAM) is unsourced.

### Gaps
- No controlled benchmark of cudaHostRegister vs cudaHostAlloc cost on modern drivers/OS; forum data are anecdotal.
- No cited threshold for "safe" pinned memory on Linux (RLIMIT_MEMLOCK interaction not researched).

## 4. Returning results to Python with nanobind (ndarray ownership, copy vs expose, thread-safety, GIL)

### Takeaway
nb::ndarray records only metadata; ownership is expressed with an owner object or an nb::capsule deleter. For a persistent pinned staging buffer the safe pattern is to copy into a freshly allocated array (capsule-owned) rather than expose the shared buffer, since exposing it would let later calls overwrite data Python still holds. GIL release around GPU work is the standard approach but creating the return array needs the GIL.

### Cited Findings
- nb::ndarray behaves like a shared_ptr-style reference; creating it in C++ only records metadata (data, dtype, shape, owner), the Python object is built by the binding layer. — [nanobind ndarray docs](https://nanobind.readthedocs.io/en/latest/ndarray.html)
- The `owner` argument is a Python object whose continued existence keeps memory alive; for newly allocated memory use an `nb::capsule` with a destructor callback, and one capsule can back several arrays; with `rv_policy::reference` and no owner nanobind assumes memory "remains valid forever", so stale arrays give use-after-free. — [nanobind ndarray docs](https://nanobind.readthedocs.io/en/latest/ndarray.html)
- Default `automatic` policy copies when there is no owner and no existing Python object; `reference` never copies; `copy` always copies. `nb::numpy` tag returns a numpy.ndarray; `nb::device::cpu` / `cuda` tags constrain source device. — [nanobind ndarray docs](https://nanobind.readthedocs.io/en/latest/ndarray.html)
- Creating a new array from C++ and calling `cast()` require the GIL; read-only inspection of an ndarray can run without it; returning an array that wraps stack memory is unsafe. — [nanobind ndarray docs](https://nanobind.readthedocs.io/en/latest/ndarray.html)
- For exported device memory produced asynchronously, nanobind's array_api object ignores the `stream` argument in `__dlpack__`; code must sync itself. — [nanobind ndarray docs](https://nanobind.readthedocs.io/en/latest/ndarray.html)

### Inferences
- Pattern (judgment, based on the above): bind the function taking input `nb::ndarray<const float, nb::c_contig, nb::device::cpu>`; read pointers; `{ nb::gil_scoped_release r; run H2D, kernels, D2H, sync; }`; then allocate `new uint8_t[n]` (or numpy via a capsule with delete[]), memcpy from pinned staging, return `nb::ndarray<nb::numpy, uint8_t>(ptr, shape, capsule)`. The copy of a 1-bit mask is cheap (1/8 of the raster bytes).
- While the GIL is released, another Python thread could call the same function and race on the shared arena/staging; guard with a std::mutex (or one arena per thread) held across the GPU section. A mutex acquired after releasing the GIL avoids deadlock between the GIL and the lock.
- Input arrays must stay alive and not be mutated by other threads while the GIL is released; hold the nb::ndarray argument (which keeps a reference) for the duration.

### Gaps
- nb::gil_scoped_release, nb::call_guard<nb::gil_scoped_release>, and nb::mod_gil_not_used were not found in the fetched part of the nanobind API reference (page truncated at 100k of 272k characters); the usage described above is from general knowledge and should be verified in the nanobind docs.
- Free-threaded Python behavior for this module was not researched.

## 5. Bit-packed mask layout and warp-ballot packing

### Takeaway
numpy.unpackbits(bitorder="little") reads each byte LSB-first, so a row-major bitmap where pixel i of a row is bit (i mod 8) of byte i/8 round-trips. On the GPU, `__ballot_sync` returns a 32-bit mask over the warp's lanes, so one lane can write one 32-bit word per 32 pixels, which on a little-endian host equals four little-endian bytes in the right order.

### Cited Findings
- `np.unpackbits` default `bitorder='big'` puts the MSB first (3 -> [0,0,0,0,0,0,1,1]); `'little'` puts the LSB first (3 -> [1,1,0,0,0,0,0,0]); `axis` unpacks along that axis (expanding by 8); `count` trims padding (negative count trims from the end). — [NumPy unpackbits docs](https://numpy.org/doc/stable/reference/generated/numpy.unpackbits.html)
- `int __ballot_sync(unsigned mask, int predicate)` returns a mask of all warp threads with a true predicate; it is synchronizing; all threads in `mask` must reach it or the program may hang; use full mask 0xffffffff only if all threads execute it; do not use `__activemask()` as a collective membership mask. The article's examples imply bit i = lane i (e.g., `__popc(mask & ((1<<lane)-1))` counts lower lanes) but does not state it explicitly. — [NVIDIA blog: Using CUDA Warp-Level Primitives](https://developer.nvidia.com/blog/using-cuda-warp-level-primitives/)

### Inferences
- Layout recommendation: output shape (N, H, ceil(W/32)) of uint32 (pad rows to a multiple of 32 pixels, or W multiple of 8 for byte-packing), row-major; pixel x in row r is bit (x & 31) of word x >> 5. View as uint8 on a little-endian host, then `np.unpackbits(a, axis=-1, bitorder="little")[..., :W]` (or `count=W`) restores the raster. x86-64 and Apple/ARM Linux are little-endian, so the uint32 to uint8 view is safe; document the assumption.
- Kernel pattern: have each warp cover 32 consecutive x pixels of a row; every thread computes its pixel's predicate (all threads, including out-of-range ones with predicate false, must reach the ballot with full mask), then `unsigned w = __ballot_sync(0xffffffff, inside); if (lane == 0) out[row*words + wordIdx] = w;`. Lane 0 writes one coalesced 32-bit word per warp; consecutive warps write consecutive words.
- Bit i = lane i is the documented behavior in the CUDA Programming Guide warp vote functions, but this was not retrieved here; verify with a unit test using a known pattern against unpackbits.

### Gaps
- Programming Guide warp vote section not retrieved (landing page only), so the explicit lane-to-bit statement is unsourced.
- No measurements of packed vs unpacked write performance.

## 6. Multi-process use: fork, spawn, context cost, MPS, architecture

### Takeaway
A CUDA context/driver state created before fork() is not usable in the child (error "Cannot re-initialize CUDA in forked subprocess" in PyTorch). Use spawn/forkserver, or fork workers before anything touches CUDA. Each process that uses the GPU holds its own context (~600 MB in one measurement) and time-slices with others unless MPS is used; so the cleanest design is a single GPU-owning process.

### Cited Findings
- The error means CUDA was initialized in the parent before fork; the standard fix is `spawn` (or `forkserver` on Linux) and not touching CUDA in the parent; hidden triggers include torch.cuda.is_available() / device_count(); worker code must be importable under spawn. — [PyTorch forum / GitHub threads (search summary)](https://discuss.pytorch.org/t/cannot-re-initialize-cuda-in-forked-subprocess-on-network-to-device-operation/138090); [pytorch/pytorch#68861](https://github.com/pytorch/pytorch/issues/68861)
- NVIDIA forum: on CUDA 8.0 a child forked after the parent's cuInit() got CUDA_ERROR_NOT_INITIALIZED from its own cuInit(); an NVIDIA-side commenter advised not to initialize CUDA before fork() if the child needs CUDA. — [NVIDIA forums](https://forums.developer.nvidia.com/t/cuda8-0-bug-child-process-forked-after-cuinit-get-cuda-error-not-initialized-on-cuinit/45764)
- A 2026 project (PECOS PR #450) hit this with cuQuantum workers and changed messages to recommend `spawn`. — [PECOS PR #450](https://github.com/PECOS-packages/PECOS/pull/450)
- Measured on a bare CUDA context: ~600 MB per process under time-slicing (MIG 1g.12gb ~60 MB); MPS had a fixed ~600 MB total regardless of process count; under time-slicing kernels from different processes never execute simultaneously; context switching cost lowered bandwidth utilization (llama3: 35% time-slicing vs 52% MPS). — [arXiv 2604.08451](https://arxiv.org/pdf/2604.08451)
- MPS consists of nvidia-cuda-mps-control, nvidia-cuda-mps-server, and a client runtime in libcuda; it provides "improved GPU utilization and reduced GPU context switching". — [NVIDIA MPS docs](https://docs.nvidia.com/deploy/mps/index.html)
- Device memory of each process still adds up under time-slicing/MPS (Stack Overflow answer by Robert Crovella, via summary; secondary). — [kiwix mirror of StackOverflow](https://kiwix.lema.org/stackoverflow.com_en_all_2023-11/questions/74632058/memory-issue-in-running-multiple-processes-on-gpu)

### Inferences
- Current setup (fork-based ProcessPoolExecutor of 4): if the parent ever calls the checker (or anything that inits CUDA, including a CUDA build check) before creating the pool, forked children break; if each child initializes CUDA independently after fork, that is allowed but yields 4 contexts (~2.4 GB total, using ~600 MB each per the one measurement; the commonly quoted 300-500 MB varies by GPU/driver) and time-sliced execution without MPS.
- Recommended architecture: one dedicated GPU process (owning arena, pinned buffers, context) receiving batches from CPU workers via a queue/pipe/shared memory (multiprocessing with spawn context or a persistent server), or run the checker in the parent and use workers only for CPU-only work created via fork before CUDA init (or via forkserver/spawn). With a 1-bit mask output, IPC payloads shrink 8x.
- MPS helps only if several processes must keep submitting small kernels concurrently; for 4 workers submitting short batches, a single GPU process avoids both context overhead and the MPS daemon operations burden. MPS requires daemon setup and may not be possible on shared/rented (Vast) hosts (unverified).
- Lazy initialization inside the GPU process (first call, not at import) avoids accidental init in a parent; use `os.register_at_fork` or a pid check to detect use in a forked child and raise a clear error.

### Gaps
- Official CUDA documentation statement on fork (driver API behaviour) not found; only forum and framework sources.
- PyTorch multiprocessing notes page did not return content (redirect page); CuPy fork docs not fetched.
- The 300-500 MB per-context figure in the brief was not independently confirmed; the only measurement found was ~600 MB (arXiv 2604.08451, GPU model not captured in the summary).
- MPS client limits, "When to use MPS" and per-client memory were not retrieved from the NVIDIA page; the ORNL slides figure (~164-221 MB per client from thread stacks) is an older estimate from search snippet only ([ORNL MPS slides](https://www.olcf.ornl.gov/wp-content/uploads/2021/06/MPS_ORNL_20210817.pdf)).

## 7. CUDA Graphs for a short chain of small kernels

### Takeaway
Graphs reduce launch overhead for many short kernels and can be updated in place when topology is unchanged, but NVIDIA's own example shows modest gains (16-26%) on sections converted, and update cost scales with the number of changed nodes. For a short chain after removing allocations, they are probably not the first optimization.

### Cited Findings
- cudaGraphExecUpdate updates an executable graph when the new graph has the same topology (same node count/types/dependencies); if update fails, destroy and re-instantiate. Update cost is "roughly proportional to the number of changes made to CUDA graph nodes". — [NVIDIA blog: Employing CUDA Graphs in a Dynamic Environment](https://developer.nvidia.com/blog/employing-cuda-graphs-in-a-dynamic-environment/)
- In that test, kernels took 2-8 microseconds each on A100; graphs of 504 and 2520 kernels gave 16% and 26% speedup (update method) on the converted sections only. — [Employing CUDA Graphs in a Dynamic Environment](https://developer.nvidia.com/blog/employing-cuda-graphs-in-a-dynamic-environment/)
- The runtime API has a `cudaGraphKernelNodeUpdate` struct with a `gridDim` field, so grid dimensions can be updated via the update APIs. — [cudaGraphKernelNodeUpdate reference](https://docs.nvidia.com/cuda/cuda-runtime-api/structcudaGraphKernelNodeUpdate.html)

### Inferences
- Rasterize + reduce chain of a handful of kernels launched once per call: launch overhead is typically a few microseconds per launch (unsourced, general knowledge), small against a PCIe copy of MBs and the Python call overhead; graphs only pay off if the chain is dozens of kernels or latency per call is under ~100 us.
- Alternative to graph updates: capture with max sizes and have kernels early-exit on actual counts read from device memory (grid-stride loops with a device-side count), keeping the graph static. This is judgment, not cited.
- Worth revisiting only after measuring: if profiling (Nsight Systems) shows idle gaps between kernels.

### Gaps
- No benchmark of graph launch vs stream launch for 3-10 kernels on Ada/Blackwell, and no number for instantiate/update time.
- Conditional/device-side graph features in CUDA 12.x were not researched.

## 8. Error handling, lifetime, validation tooling

### Takeaway
Check errors from every API call, call cudaGetLastError after each launch, and expect asynchronous errors to surface at the next sync or at cudaFree. Free the arena explicitly at module unload with care about interpreter and driver shutdown ordering; validate with compute-sanitizer.

### Cited Findings
- cudaFree "may also return error codes from previous, asynchronous launches", so a late error can be reported from an unrelated call. — [Runtime API: Memory Management](https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html)
- Memory-pool allocation errors are returned immediately and cannot fail asynchronously. — [stream-ordered allocator part 1](https://developer.nvidia.com/blog/using-cuda-stream-ordered-memory-allocator-part-1/)
- cudaHostRegister memory must be unregistered with cudaHostUnregister using the same base address. — [Runtime API: Memory Management](https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html)

### Inferences
- Launch pattern (judgment): `kernel<<<...>>>(...); CUDA_CHECK(cudaGetLastError());` catches launch-configuration errors; execution errors appear at the cudaStreamSynchronize that ends the call, so check its return and convert to a Python exception. After a sticky error (illegal address) the context is unusable, so the module should mark itself failed.
- Lifetime: hold the arena in a singleton whose destructor calls cudaFree/cudaFreeHost, or register `atexit` from Python to release; avoid relying on static destructors after the CUDA runtime has been torn down (cudaErrorCudartUnloading); the OS/driver reclaims everything on process exit anyway, so a missing free at exit is benign. In a forked child never free a parent's arena.
- Validation: run the unit tests under `compute-sanitizer --tool memcheck` and `--tool racecheck` (also `initcheck`, `synccheck`); arena sub-allocation hides out-of-bounds writes between sub-buffers from memcheck because they stay inside one allocation, so add guard bands or a debug mode that uses separate cudaMalloc per buffer.

### Gaps
- No sources fetched for compute-sanitizer tools, cudartUnloading behavior, or atexit/static-destructor ordering; statements above are unsourced engineering judgment.
- No info found on CUDA 13.x changes relevant to these APIs (e.g., deprecations); versions of fetched sources: stream-ordered allocator blog is CUDA 11.2-era, unified memory blog is CUDA 8/9-era (Pascal/Volta), runtime API pages are the current docs (version not captured in fetch).
