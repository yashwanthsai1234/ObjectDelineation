# CUDA graphs, streams and GPU sharing for the RL run (RTX 5090)

Date: 2026-11-06. Scope: (A) gpucheck CUDA library, (B) LoRA GRPO update step, pinned H2D copies, colocated vLLM.
Research limits: this note was built from fetched doc pages plus web search. Claims without a fetched source are marked UNVERIFIED.
Numbers marked "typical" are from general CUDA experience, not from a fetched page. Measure them with nsys before relying on them.

## Sources

- [S1] NVIDIA CUDA Programming Guide, CUDA Graphs: https://docs.nvidia.com/cuda/cuda-programming-guide/04-special-topics/cuda-graphs.html
- [S2] NVIDIA CUDA Programming Guide, stream-ordered allocation: https://docs.nvidia.com/cuda/cuda-programming-guide/04-special-topics/stream-ordered-memory-allocation.html
- [S3] NVIDIA CUDA C++ Best Practices Guide: https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html
- [S4] PyTorch CUDA semantics (docs 2.15): https://docs.pytorch.org/docs/2.15/notes/cuda.html
- [S5] torch.cuda.make_graphed_callables (docs 2.15): https://docs.pytorch.org/docs/2.15/generated/torch.cuda.make_graphed_callables.html
- [S6] CUDAGraph Trees (docs 2.8 page; the 2.15 page returned empty): https://docs.pytorch.org/docs/2.8/torch.compiler_cudagraph_trees.html
- [S7] vLLM sleep mode: https://docs.vllm.ai/en/latest/features/sleep_mode.html
- [S8] vLLM and TRL colocate: https://docs.vllm.ai/en/latest/training/trl/ and https://huggingface.co/blog/vllm-colocate
- [S9] NVIDIA Megatron Bridge CUDA graphs doc: https://docs.nvidia.com/nemo/megatron-bridge/nightly/training/cuda-graphs.html
- [S10] Triton/ptxas sm_120 reports: https://discuss.pytorch.org/t/rtx-5070-ti-blackwell-pytorch-nightly-triton-still-getting-sm-120-is-not-defined-for-option-gpu-name-error/220460 and https://forums.developer.nvidia.com/t/running-llama-3-1-8b-fp4-get-triton-error-value-sm-121a-is-not-defined-for-option-gpu-name/348808
- [S11] PyTorch release 2.12 note on CUDA 13.2 and Blackwell: https://dev-discuss.pytorch.org/t/introducing-cuda-13-2-and-deprecating-cuda-12-8-release-2-12/3337.md
- [S12] Repo code read for this note: cuda/gpucheck.cpp, cuda/gpucheck_py.cpp.

Not fetched (no source): nanobind docs, CUDA fork behaviour page, PyTorch GitHub issues at the v2.13/v2.14 tags. Items depending on them are marked UNVERIFIED.

## 0. What the repo does today (A)

- Every `Buf` is a `cudaMallocManaged` allocation plus `cudaMemset` plus `cudaDeviceSynchronize()` in its constructor, and `cudaFree` in its destructor [S12, gpucheck.cpp lines 80-100].
- `for_each` launches 128-thread kernels with no stream argument (default stream) and `sync()` is a full `cudaDeviceSynchronize` [S12].
- Consequence: each Buf costs a device-wide sync, and `cudaMalloc`/`cudaFree` sync across all streams [S2]. Managed memory also adds page-migration cost. UNVERIFIED for the 5090 on Linux, but this is the first thing to remove.
- The binding already releases the GIL around the GPU call and validates inputs before release [S12, gpucheck_py.cpp lines 76-82]. Keep that.

## 1. Part A: is the plan optimal?

Verdict: the plan is sound. Biggest wins in order: (1) remove per-Buf malloc/sync, (2) bit-pack output, (3) one H2D and one D2H, (4) only then 2-stream chunking and graphs.

### 1.1 Pinned memory cost and reuse
- Pinned memory gives the highest host-device bandwidth, and async copies require it [S3].
- Do not overuse it: "Excessive use can reduce overall system performance" [S3]. PyTorch docs give the same warning [S4].
- `cudaMallocHost` is slow (it page-locks) and may synchronize. Typical cost is ms for tens of MB. UNVERIFIED number. So allocate once, grow only, never per call.
- Size to the worst case: input arrays are small (vertices), output is the big one. Packed 512x512 = 32 KB per mask, so 1024 masks = 32 MB per buffer. 4 buffers (2 pairs x in/out) is well under 200 MB. Fine.
- Reuse rule: a pinned buffer may be rewritten by the host only after the stream event for its last copy has completed. Use one `cudaEvent` per buffer pair, `cudaEventSynchronize` before reuse.

### 1.2 cudaMemcpyAsync ordering
- Operations in one stream run in order. Different streams may run in any relative order [S4].
- Copy/kernel overlap needs: pinned host memory, non-default streams, and a device with copy engines (`asyncEngineCount`) [S3]. The 5090 has copy engines (UNVERIFIED count). The legacy default stream serialises, so do not launch on stream 0 [S3].
- Order per chunk on its own stream: H2D copy, kernels, D2H copy, `cudaEventRecord`. Host waits with one `cudaStreamSynchronize` or event per chunk, not `cudaDeviceSynchronize`.
- Pageable-memory `cudaMemcpyAsync` silently falls back to staged, mostly synchronous behaviour. UNVERIFIED in the fetched text, but [S3] says the async version requires pinned memory.

### 1.3 When 2-stream chunking pays
- [S3]: each transfer has fixed overhead, so batching small transfers into one larger copy is "significantly better". Staging across `nStreams` helps when execution time and transfer time are comparable [S3].
- Numbers (typical, UNVERIFIED): a copy call costs about 5-10 us fixed; PCIe 5 x16 sustains about 50 GB/s pinned (theoretical 64 GB/s).
- Worked sizes, packed output: 8 masks = 256 KB = about 5 us of bandwidth, so latency dominates. Chunking is pure overhead. 1024 masks = 32 MB = about 0.65 ms. Chunking can hide at most that much, and only if kernels take comparable time.
- Unpacked output would be 8x more: 256 MB = about 5 ms. That is why bit-packing matters more than streams.
- Recommendation: one stream, one H2D, one D2H below about 128 polygons (or below about 4 MB output). Use 2 chunks above that. Make the threshold a measured constant, not a guess.
- Do not use more than 2 streams. One copy engine per direction is typical (UNVERIFIED), so H2D of chunk B and D2H of chunk A can overlap, but more streams add nothing.

### 1.4 cudaMallocAsync vs a custom arena
- Stream-ordered allocation avoids the cross-stream sync of `cudaMalloc`/`cudaFree` and caches physical pages [S2].
- By default the pool tries to minimise physical memory and releases the excess at the next sync. Set `cudaMemPoolAttrReleaseThreshold` to `UINT64_MAX` to keep memory [S2].
- Free memory is reusable in the same stream at once; after a CPU sync of that stream it is reusable in any stream [S2].
- Opportunistic reuse can cause run-to-run allocation variation [S2].
- Choice: a grow-only `cudaMalloc` arena is simplest, deterministic and one allocation. `cudaMallocAsync` with an infinite threshold is equivalent and less code to get wrong when sizes vary per call. Either is fine. Pick the arena if you want a fixed device address (needed for CUDA Graphs without updates, see 1.5). Pick `cudaMallocAsync` if you drop graphs.
- Growth must happen only when no GPU work reads the arena (after sync). With the mutex held this is easy.
- Cost to note for sharing the GPU with PyTorch: every CUDA process holds its own context (typical few hundred MB, UNVERIFIED). Budget it in vLLM `gpu_memory_utilization`.

### 1.5 CUDA Graphs for a small fixed kernel sequence
- Graphs cut CPU launch cost because setup is done in advance; this matters most for short kernels [S1]. [S1] gives no microsecond numbers.
- Typical (UNVERIFIED): about 2-5 us CPU per kernel launch, about 1-2 us per kernel node in a graph launch. A library call of 5-10 kernels saves tens of us.
- Update rules [S1]: `cudaGraphExecUpdate` needs identical topology. For few changes, update nodes one by one with `cudaGraphExecKernelNodeSetParams`. Memcpy nodes allow only 1D changes, with the same kind. Topology change needs re-instantiation.
- Sizes vary per call (8 to 1024 polygons). Options: (a) launch with a max grid and a device-side `n` guard, store `n` in a device scalar; (b) update kernel node params each call (host cost eats the gain); (c) cache one instantiated graph per size bucket (for example 8/32/128/512/1024).
- Recommendation: skip graphs at first. Per call, GPU time at 512x512 x 1024 masks is likely ms, so launch overhead is below 1-2%. Revisit only if nsys shows the small batches (8-32 polygons) are launch-bound. If so, use option (c). UNVERIFIED that the gain is worth it.
- The current code launches kernels in a loop with a device-wide sync between steps. Removing those syncs (use stream order) gives more than graphs.

### 1.6 __ballot_sync bit packing and coalescing
- Coalescing: a warp's loads and stores merge into as few 32-byte-segment transactions as possible on cc 6.0+ [S3]. Block sizes should be multiples of 32 [S3].
- Packing scheme: lane i tests pixel x = 32*w + i, `m = __ballot_sync(0xffffffff, inside)`, lane 0 writes `m` to `out[n][y][w]`. Bit i of word w is pixel 32w+i.
- All 32 lanes must reach the ballot, including out-of-range lanes. Set `inside=false` for x >= W. This also zeroes the padding bits of the last word, which `unpackbits` then slices off.
- numpy side: `out.view(np.uint8)` on a little-endian host puts byte 0 = bits 0-7 of word 0. `np.unpackbits(..., axis=-1, bitorder="little")` then yields pixel order x = 0..32*words-1. Slice `[..., :W]`. Correct on x86 and ARM little-endian. (numpy docs for `bitorder`: UNVERIFIED, not fetched.)
- Store pattern: one 4-byte store per warp, consecutive warps write consecutive words, so stores are contiguous across warps of a block. Good enough. A shared-memory staging of a row for 128-byte stores is optional.
- 512 wide = 16 words per row; use block = (32, 4 to 8) so one block covers several rows.
- Row stride in words is `ceil(W/32)`. If masks have different W/H in one batch, keep an offsets array (as `mask_start` does today).

### 1.7 Thread safety and fork safety
- A mutex around the whole GPU section is correct for the arena and the pinned buffers (shared mutable state). The CUDA runtime API itself is thread-safe (UNVERIFIED, not fetched). With one arena you serialise calls anyway, so the mutex costs nothing.
- Release the GIL before taking the mutex, not after, to avoid a Python thread holding the GIL while waiting on the mutex. Order: `nb::gil_scoped_release` then lock. nanobind pattern is already in use [S12]. nanobind docs not fetched. UNVERIFIED.
- Do not call back into Python while holding the mutex.
- Fork: a process that already initialised CUDA cannot use CUDA in a forked child. PyTorch raises "Cannot re-initialize CUDA in forked subprocess" (known behaviour; fetched docs do not quote it, UNVERIFIED). Rule: create the CUDA context lazily in the process that uses it, or use `spawn`. Never call `cuda_build()` paths that touch the runtime before a fork. `device_count()` calls `cudaGetDeviceCount`, which can initialise the driver. UNVERIFIED whether it poisons fork; avoid calling it before fork.
- Two processes on one GPU (trainer and checker) time-slice by default. Running the checker inside the trainer process avoids a second context. UNVERIFIED trade-off with the vLLM memory budget.

## 2. Part B: graphing the LoRA GRPO update

### 2.1 make_graphed_callables: hard requirements [S4, S5]
- `sample_args` must be Tensors only [S5].
- `requires_grad` of each sample tensor must match the real inputs [S4, S5]. Hidden states between decoder layers need `requires_grad=True` when the layer below is trainable. The first layer's input (embeddings from a frozen embedding) may be `requires_grad=False`. With LoRA in every layer, input to layer 0 usually does not require grad unless gradient checkpointing enabled `enable_input_require_grads`. Match what eager does.
- Only parameters may be trainable. Buffers must have `requires_grad=False` [S5]. LoRA adapter weights are parameters, so OK.
- You may not add or remove parameters or buffers after passing the module [S5]. So graph after `get_peft_model`, and do not merge/unmerge or swap adapters afterwards.
- No module hooks at the time of passing [S5]. Gradient checkpointing wrappers, accelerate hooks and PEFT hooks can break this. Check `_forward_hooks`/`_backward_hooks` on each layer before graphing. UNVERIFIED which hooks PEFT/TRL adds.
- No double backward [S5].
- Autocast must use `cache_enabled=False` [S5]. Reason: the cached fp32-to-bf16 weight casts are tensors created inside capture and would go stale. This matters here: fp32 LoRA adapters under autocast are exactly that case.
- Warmup: default 3 iterations; DDP needs 11 [S5]. Warmup runs on a side stream [S4].
- Unused inputs raise unless `allow_unused_input=True` [S5]. Attention mask and position ids are not Tensors-only friendly if they are None or Python ints. Pass only tensors, and bake constants in.
- Tuple form: pass layers in the exact order they run, only if they always run in that order and never concurrently [S4, S5]. For one layer called 28 times, you need 28 separate graphed callables (one per layer instance, since each has own parameters), or one graphed callable per distinct shape. With 7 buckets x 28 layers that is 196 graphs. UNVERIFIED memory cost; measure.
- Pool: `pool=` shares memory across graphs. Safe only if replays happen in the captured order and never concurrently. Replaying one graph can clobber another's outputs unless cloned [S4]. Forward of layer k feeds layer k+1, so outputs are read before the next replay overwrites the shared pool if the order is strict. Backward interleaving makes this risky. UNVERIFIED that `make_graphed_callables` with a shared pool across buckets is safe; test numerics against eager.
- Per-bucket graphs can share a pool only if buckets are never replayed in the same step in an interleaved order. They are replayed one bucket at a time, so a pool per bucket-group is plausible.

### 2.2 Pitfalls specific to this training loop
- No CPU sync inside the graphed region: no `.item()`, no `torch.nonzero` with data-dependent shape, no boolean mask indexing [S4]. Qwen3-VL code paths: the vision tower and the multimodal merge (placing image embeddings into text positions via masked scatter) are data-dependent. Graph only the text decoder layers, as planned.
- Dynamic shapes are prohibited [S4]. Bucketing to 7 shapes is the right fix. Batch size must also be fixed per bucket (pad the batch).
- Gradient accumulation across replays: the backward graph writes into static grad buffers and `.grad` handling is done by the autograd function wrapper. UNVERIFIED how it accumulates into existing `.grad` across replays. PyTorch source (`torch/cuda/graphs.py`, `Graphed.backward`) returns grads as outputs of the graph to autograd, which accumulates them into `.grad` as usual. Verify with a two-microbatch test against eager. Also call `optimizer.zero_grad(set_to_none=True)` outside the graph.
- Static outputs: graph outputs are reused on the next replay. Clone anything you keep across steps [S4].
- Parameters captured by address: the optimizer must update parameters in place. Never replace parameter tensors (`load_state_dict` is fine, `param.data = new` is not). A fused or foreach AdamW outside the graph is fine and is in-place.
- LoRA dropout: CUDA RNG is allowed in capture [S4]. With the default generator this works (PyTorch registers the generator state with the graph). With multiple generators use `graphsafe_get_state`/`graphsafe_set_state` [S4]. If dropout is 0, remove the dropout op entirely to be safe.
- Gradient checkpointing: `torch.utils.checkpoint` recomputes the forward in backward and is a different autograd path. Combining it with graphed callables is not covered by [S4, S5]. The plan to turn it off inside the graphed region is right. The memory cost: all activations of 28 layers at the largest bucket (1536 tokens x batch) must fit with the static buffers. NVIDIA's Megatron doc says graph buffers stay resident and a run near the memory limit can OOM, and that full recompute conflicts with scoped graphs [S9]. Same direction as our plan. Megatron-specific, so UNVERIFIED for plain PyTorch.
- Memory with vLLM colocated: static graph buffers stay allocated while vLLM sleeps. Budget them.
- Capture is slow for 196 graphs. Capture lazily per bucket on first use, or in a warmup phase before the first step. Estimate seconds per graph. UNVERIFIED.

### 2.3 CUDAGraph Trees (mode="reduce-overhead") [S6]
- Skipped (runs eager) when the function has CPU ops, multi-device ops, listed incompatible ops (for example `aten.multinomial`), or free unbacked symbols [S6].
- Input mutation: in-place writes to eager inputs cause "skipping due to mutated inputs". Mutating parameters and buffers, or outputs of other cudagraph functions, is supported. `torch._inductor.config.cudagraph_support_input_mutation = True` enables it [S6].
- Shapes: it re-records per unique input shape. Good for few fixed shapes, costly when shapes change often. Pad to a few shapes. `cudagraph_skip_dynamic_graphs=True` skips graphs for dynamic-shape functions [S6]. With 7 buckets this is the intended use case.
- 64 KB device memory per kernel launch in a graph up to CUDA 12.4 / driver 550 [S6]. Newer drivers: UNVERIFIED.
- One shared pool across all captures, so no extra memory versus eager [S6]. Without trees, no reuse and large memory regressions [S6].
- Iteration boundaries: `torch.compiler.cudagraph_mark_step_begin()` at the start of each step, or clone outputs. Accessing an overwritten output raises an error [S6]. Training starts a new iteration on each invocation unless a backward is pending [S6].
- Recompile limit: `torch._dynamo.config.cache_size_limit` (default 8, UNVERIFIED for 2.13/2.14) counts per code object. 7 buckets x dynamic=False fits just under 8 for one compiled module; if the same layer class is compiled per instance, each layer object may hit the limit separately. Raise `cache_size_limit` explicitly. UNVERIFIED.
- Trees plus LoRA plus fp32 adapters under autocast: not found in sources. UNVERIFIED. Test.

### 2.4 Regional compilation of repeated layers
- Idea: compile each decoder layer (same class, same graph) so Inductor compiles the layer graph once and reuses the code across 28 layers. PyTorch has an official recipe on regional compilation ("Reducing torch.compile cold start compilation time with regional compilation"). Not fetched, UNVERIFIED details.
- Reuse works only when the guards match: same shapes, dtypes, same Python class. Per-layer parameters are inputs, not constants, in recent versions (`inline_inbuilt_nn_modules`, default on since 2.5-2.6; UNVERIFIED for 2.13).
- Cold compile then scales with number of distinct shapes (7 buckets) rather than 28 layers. Use `fullgraph=True` on the layer to catch graph breaks. Set `dynamic=False`.
- Persist the cache: `TORCHINDUCTOR_CACHE_DIR` on a volume, so rented Vast instances skip recompiles. Mega-cache (`torch.compiler.save_cache_artifacts`) UNVERIFIED for 2.13.

### 2.5 Which option to choose
- Launch-bound update at 60% GPU busy means up to 40% headroom, so the best case is about 1.6x (Amdahl). NVIDIA says the gain depends on how launch-bound you are and the captured scope; compare replay iterations, never the capture step [S9]. No public number found for LoRA on small VLMs. Do not promise a speedup.
- Cheapest first: (1) regional `torch.compile` per layer with `mode=None` (fuses elementwise ops, fewer launches, no graph constraints). (2) Add `mode="reduce-overhead"` on the same layers. (3) Manual `make_graphed_callables` only if (2) is skipped or too fragile.
- Manual graphs give the most control but the most pitfalls (hooks, autocast cache, static inputs, 196 graphs).
- Tool: take an nsys or torch profiler trace first. Confirm idle gaps are CPU launch gaps and not syncs (`.item()`, `loss.item()`, `clip_grad_norm_` returning a tensor then `.item()`). Removing syncs may recover much of the 40% for free.

### 2.6 Blackwell sm_120, Triton and ptxas
- RTX 5090 is sm_120. PyTorch needs CUDA 12.8+ builds; older builds fail "not compatible with sm_120" [forum: https://discuss.pytorch.org/t/pytorch-support-for-sm-120/222119].
- Release 2.12 notes: CUDA 13.2 is added and CUDA 12.8 is deprecated; Blackwell SM 10.0, 11.0, 12.0 listed for CUDA 13.0.x/13.2.x (aarch64 line quoted; x86 matrix UNVERIFIED) [S11].
- Triton can fail with `ptxas fatal : Value 'sm_120' is not defined for option 'gpu-name'` when its bundled ptxas is too old. Same class of error for sm_121a. Workaround reported: `export TRITON_PTXAS_PATH=/usr/local/cuda/bin/ptxas` with a CUDA 12.8+ toolkit [S10]. The variable only redirects ptxas, not cuda.h lookup [S10].
- Whether cu13 wheels bundle ptxas for sm_120 in torch 2.13/2.14: not found. UNVERIFIED. Preflight check: run one tiny `torch.compile` on the rented box, and `python -c "import triton,os;print(triton.__file__)"` then list `triton/backends/nvidia/bin/`. If ptxas errors, install `cuda-nvcc` and set `TRITON_PTXAS_PATH`.
- Fallback seen in forums: Triton treating newer Blackwell as sm80, losing Blackwell features [forum above, NVIDIA thread in S10]. Not a fix.
- Attention backends on sm_120 (third-party reports only, UNVERIFIED): SDPA with the cuDNN backend is reported as the best path on RTX 5090 for non-causal BF16 forward; FlashAttention-4 reported blocked on sm_120; FlexAttention works via its Triton FA2-style kernel. Sources: search results (gist by solatticus; keenable write-up). Check `torch.backends.cuda` flags and benchmark `sdpa_kernel` choices on the actual host, with causal masks and your GQA shapes.
- Our model is Qwen3-VL-2B with padded buckets, so SDPA with an explicit mask or `is_causal` is simplest. FlexAttention adds compile risk for little gain at 1536 tokens. UNVERIFIED claim on gain; measure.

## 3. Pinned memory plus non_blocking H2D of large pixel tensors on a side stream

Correct pattern [S4]:

```python
copy_stream = torch.cuda.Stream()
pinned = torch.empty(shape, dtype=torch.uint8, pin_memory=True)   # allocate once, reuse
pinned.copy_(cpu_batch)                                           # fill on host
with torch.cuda.stream(copy_stream):
    dev = pinned.to("cuda", non_blocking=True)
    ready = torch.cuda.Event(); ready.record(copy_stream)
torch.cuda.current_stream().wait_event(ready)                     # consumer waits, host does not
dev.record_stream(torch.cuda.current_stream())                    # allocator keeps dev alive for consumer
```

- Allocation stream: `dev` is allocated by the caching allocator on `copy_stream`. If it is used on another stream, call `record_stream` so the block is not reused early [S4]. Alternative with no `record_stream`: allocate `dev` on the compute stream, then copy into it on the side stream after `copy_stream.wait_stream(current)`.
- `wait_stream` is needed before using memory that the other stream may still be using, because of caching-allocator reuse [S4]. Waiting right after launch kills the overlap [S4].
- Do not mutate or free the pinned source until the copy finished (`ready.synchronize()` or an event check). Reusing one pinned buffer for the next batch before then corrupts data silently. Keep 2 pinned buffers (double buffer).
- `.pin_memory()` on a fresh tensor each step is slow and allocates through the pinned caching allocator; keep a reusable buffer. `pin_memory=True` in DataLoader gives pinned batches [S4].
- `non_blocking=True` from pageable memory is not truly async. Source must be pinned [S3, S4].
- Size: 200 MB at about 50 GB/s is about 4 ms; 1.6 GB is about 32 ms (typical, UNVERIFIED). Pinning 1.6 GB is a real allocation (page-locking cost, host RAM); size it once for the max and slice it [S3 overuse warning].
- Send uint8 pixels and convert/normalise on GPU. It is 4x smaller than fp32.
- Backward streams: backward ops run on the stream of their forward op, and the autograd engine no longer syncs the default stream with all backward ops. If `backward()` ran in a different stream context, `wait_stream` before using grads [S4].
- Does the copy overlap with compute at all? Only if copy engines are free and the compute stream is not waiting on `ready` immediately. Prefetch batch N+1 while training on batch N.
- Peak GPU memory: with 1.6 GB tensors, check allocator fragmentation; set `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` (UNVERIFIED for interaction with cudagraph pools; test).

## 4. Missed opportunities for GPU utilisation in colocated RL

- Colocate means generation and training cannot overlap on that GPU; vLLM shares memory with the trainer [S8]. It is the simple, no-overlap setup. So utilisation gains come from within each phase and from overlapping CPU work (rewards) with GPU work.
- Sleep level: level 1 offloads weights to host RAM and discards KV cache; level 2 discards weights and KV cache and backs up little [S7]. For RL the weights change after every optimizer step, so level 2 is right: wake `weights` only, load the new weights, then wake `kv_cache`. This lowers peak memory [S7]. Level 1 wastes host RAM and a copy of weights that are immediately stale. With LoRA: sync only adapters (small) if vLLM supports LoRA serving, otherwise merge into base weights for the push (UNVERIFIED for vLLM LoRA plus sleep level 2: frozen base must be preserved via `sleep_preserve_parameter_names` [S7], or use level 1 so the frozen base is not reloaded).
  - Decision for LoRA with frozen base: level 1 keeps base on host and restores it by copy, adapters pushed separately. Level 2 needs the base reloaded from the trainer each step, which costs a full 2B bf16 (about 4 GB) copy per step, about 0.1 s over PCIe (typical). Either is cheap relative to generation. UNVERIFIED which has less bug risk; test both.
- Overlap reward scoring with the next generation: the GPU checker (A) and CPU rewards can run on a worker thread or process while vLLM generates the next batch. On one GPU this competes for SMs but the checker is short. The nanobind GIL release in [S12] already enables a Python thread to score while another thread works. Scoring chunk k while chunk k+1 generates is a pipeline with no GPU memory change. UNVERIFIED gain; needs profile.
- Async GRPO (off-policy by one step): needs a separate generation engine running concurrently with training. On one 32 GB GPU that means the trainer and vLLM both resident, no sleep, which contradicts colocate sleep. Not found in TRL docs as a one-GPU option [S8 search]. Axolotl's guide states colocate is slower than server mode since no overlap, and fits small models only (third-party). Verdict: not worth it on one GPU for a 2B model unless memory allows both resident. UNVERIFIED feasibility: 2B bf16 trainer plus activations plus a vLLM KV cache could fit in 32 GB, but SM contention removes most gain.
- Cheaper wins (UNVERIFIED sizes): larger generation batch per vLLM call (more sequences per wake); fewer sleep/wake cycles per optimizer step (generate for several update steps' worth of prompts, then train); drop `.item()` syncs; length-bucketing to cut padding; vLLM CUDA graphs stay on during generation (default).
- Keep CUDA graphs of vLLM and of the trainer in separate pools. The vLLM sleep mode frees its own graph and KV memory [S7]. Trainer graph pools stay allocated, so budget them in `gpu_memory_utilization`.

## 5. Concrete recommendations (numbers)

See the 15-line summary in the chat reply. Measure before and after each item with nsys. Do not trust any "typical" number above without a trace on the 5090.
