# CUDA execution model limits on PyTorch training speed: launch overhead, streams, CUDA Graphs, fusion, host-device sync

Scope note: notes written 2026-11-06. Research used about 17 tool calls. PyTorch docs now live at docs.pytorch.org with versions up to 2.15. Several pages were only partly readable, which is flagged under Gaps.

## What does a kernel launch cost (CUDA driver/runtime, PyTorch eager dispatcher)?

### Takeaway
No controlled, versioned per-op benchmark for Ada (RTX 4090) or Hopper on CUDA 12/13 turned up. Published figures range from about 3 us (best-case null-kernel GPU launch), to roughly 5-20 us per eager op (mostly blog or forum estimates), to about 30 us per kernel end-to-end in an H100 decode workload. NVIDIA's own CUDA Graph docs give layered estimates that put the Python and PyTorch dispatcher layers at 10-100 us and 5-20 us. A small kernel that runs for less than this per-op CPU cost leaves the GPU idle.

### Cited Findings
- NVIDIA's CUDA Graph basics page gives typical per-launch overhead layers: language transitions (Python-to-C++, marshalling, GIL) 10-100 us; runtime processing (operator dispatch, memory/workspace management, cuDNN/cuBLAS transitions, parameter checks) 5-20 us; driver operations (validation, command buffer, queue submission) 5-15 us; hardware submission (PCIe, GPU command processor) 1-5 us. The page says these "may vary" with hardware, driver, framework and workload. No GPU or CUDA version is stated. — [NVIDIA DL CUDA Graph basics](https://docs.nvidia.com/dl-cuda-graph/cuda-graph-basics/introduction.html)
- The same page says launch overhead hurts when the GPU is much faster than the CPU's preparation (it names older CPUs paired with modern GPUs), when the software stack has heavy dispatch overhead, and when there are many kernels that each run only a few microseconds. The GPU finishes its kernel and idles while the CPU prepares the next launch. — [NVIDIA DL CUDA Graph basics](https://docs.nvidia.com/dl-cuda-graph/cuda-graph-basics/introduction.html)
- Measured on V100 with CUDA 10.1 (Sep 2019, so older than 2023 and possibly outdated): a 2.9 us kernel costs 9.6 us per kernel with a sync after each launch, 3.8 us per kernel when launches overlap with one sync per timestep, and 3.4 us per kernel as a graph. Graph capture plus instantiation took about 400 us (paid once), and the first graph launch was about 33% slower. These are NVIDIA's own numbers for a 20-kernel chain of 2.9 us kernels. — [NVIDIA Technical Blog, Getting Started with CUDA Graphs](https://developer.nvidia.com/blog/cuda-graphs/)
- A forum post says null-kernel launch overhead appears to have dropped to about 3 us on recent hardware and software. The GPU model is unspecified and it is not an official figure. — found via [NVIDIA forum search results](https://forums.developer.nvidia.com/t/too-much-time-for-kernel-launch-latency/233645/3)
- A cluster paper measured launch latency on the RTX 4090 at about 5 us, but with clpeak through OpenCL, not CUDA. — [Dalek, arXiv 2508.10481](https://arxiv.org/pdf/2508.10481)
- A 2026 arXiv paper estimates about 30 us of CPU launch overhead per kernel in a batch-1 LLM decode workload on H100. It covers the framework and host path, not the raw API. CUDA Graphs improved decode latency 1.259x at that context length. — [arXiv 2605.30571](https://arxiv.org/pdf/2605.30571); I saw only a search summary, not the full paper.
- Eager per-op estimates from secondary sources: Python-to-C++ dispatch 1-5 us per op ([ZeroEntropy](https://www.zeroentropy.dev/concepts/pytorch-internals/)); 5-15 us of Python plus dispatcher cost per op ([blog](https://vishalkumarr2.github.io/Learning_from_building_systems/ml-compilers/weeks/week-02/day-12-eager-vs-graph-mode.html)); about 5 us in a noisy Jiterator measurement ([PyTorch dev-discuss](https://dev-discuss.pytorch.org/t/keeping-pytorchs-ops-maintainable-the-jiterator/468.md)); about 20 us per kernel in an anecdotal Reddit benchmark (36 kernels, 700 us per step). None cites a controlled test, so treat them as rough. A DTensor-specific figure (up to 580 us per op before optimization) does not transfer to plain tensors ([veScale](https://arxiv.org/pdf/2509.07003)).
- NVIDIA/PyTorch description of the symptom: at small batch, "CPU overhead can become larger than GPU run time", and the Nsight timeline shows the CPU at 100% load with the GPU idle and empty space between kernels. In Mask R-CNN (batch 1 per GPU) the backbone ran in 31 ms eager and 6 ms graphed, 5x on that part and about 1.7x overall (graph covered mostly the ResNet backbone). — [PyTorch blog, Accelerating PyTorch with CUDA Graphs](https://pytorch.org/blog/accelerating-pytorch-with-cuda-graphs/) (page updated 2024-11-15; the work is from PyTorch 1.10, 2021)
- PyTorch's own wording: eager execution has overheads from "Python, C++, and CUDA driver" layers, and a graph replay skips all of them with a single `cudaGraphLaunch` call. — [same blog](https://pytorch.org/blog/accelerating-pytorch-with-cuda-graphs/)

### Inferences
- Launch-bound test: if a kernel's GPU time is below the per-op CPU cost (assume roughly 10-30 us for eager PyTorch including autograd and module overhead; this is an estimate, not a measurement), the GPU idles between kernels. A 2B model at micro-batch 4 with LoRA has many small kernels (LoRA A/B matmuls, norms, elementwise ops, dropout, casts). That fits the observed 100% main-thread CPU with about 60% GPU use.
- Gradient checkpointing re-runs the forward, so the CPU issues roughly 1.3-1.5x as many ops per step. This is my estimate and has no source.
- A CPU-bound loop is sensitive to single-core speed. NVIDIA names an older CPU paired with a modern GPU as a risk. Rented hosts' CPU clocks vary, so check the CPU model.

### Gaps
- No published launch-latency number for RTX 4090 or H100 under CUDA 12.x/13.x with CUDA C++, measured directly. Measure it with a back-to-back empty-kernel microbenchmark on the target host.
- No controlled per-op PyTorch eager overhead number for a current release.
- CUDA C++ Best Practices Guide and Programming Guide launch-overhead statements were not retrieved (the Programming Guide fetch returned only the table of contents).

## What do CUDA streams do and NOT do?

### Takeaway
A stream is an in-order queue. Different streams may overlap on the GPU, but launching from one Python thread still costs the same CPU time per kernel, so streams do not cut launch overhead. They help only when the GPU has spare capacity and independent work (for example copy/compute overlap or independent kernels). They need explicit synchronization.

### Cited Findings
- PyTorch semantics: "By default, GPU operations are asynchronous." Operations inside one stream are serialized in creation order. With non-default streams, "it is the user's responsibility to ensure proper synchronization." Timings without synchronization are inaccurate. `CUDA_LAUNCH_BLOCKING=1` forces synchronous execution for debugging. — [PyTorch CUDA semantics (2.15)](https://docs.pytorch.org/docs/2.15/notes/cuda.html)
- A PyTorch forum thread: with a single Python thread on one core, kernel launches dominate when kernels are tiny. A responder said that launching onto several streams from several threads from Python helps little because of the GIL (a C++ extension would be needed), and that one batched kernel is much more efficient than many launches across streams. These are forum opinions. — [discuss.pytorch.org: How to minimize the gaps between kernels on GPU?](https://discuss.pytorch.org/t/how-to-minimize-the-gaps-between-kernels-on-gpu/212841)
- NVIDIA lists Python-to-C++ boundary crossing and GIL acquire/release as launch-cost sources, and recommends CUDA Graphs for short kernels. — [NVIDIA DL CUDA Graph basics](https://docs.nvidia.com/dl-cuda-graph/cuda-graph-basics/introduction.html)
- PyTorch tuning guide: `DataLoader(pin_memory=True)` enables faster asynchronous host-to-GPU copies. `num_workers>0` overlaps loading with training. The default is `num_workers=0`. No numbers given. — [PyTorch Tuning Guide](https://docs.pytorch.org/tutorials/recipes/recipes/tuning_guide.html)
- NVIDIA's graph blog shows that overlapping launches without per-kernel syncs (3.8 us per kernel) is already far better than syncing after each kernel (9.6 us), because the CPU runs ahead of the GPU. This is the benefit of asynchronous execution itself, not of extra streams. — [NVIDIA blog](https://developer.nvidia.com/blog/cuda-graphs/)

### Inferences
- In this setup streams offer little for a CPU-bound main thread. Extra streams add cross-stream sync and memory-allocator complexity (caching allocator is per-stream) without lowering the CPU cost per op.
- Streams are useful here only for overlapping host-to-device copies (pinned memory, `non_blocking=True`) with compute, which matters mainly for image-heavy VLM batches.

### Gaps
- No primary NVIDIA statement found that explicitly says "multiple streams do not reduce CPU launch cost". The conclusion rests on NVIDIA's cost breakdown plus forum evidence.
- Programming Guide streams section not retrieved.

## CUDA Graphs: how they work, capture rules, PyTorch APIs, CUDAGraph trees, limits

### Takeaway
A graph records the kernels and arguments once and replays them with one launch call, removing Python, dispatcher and driver costs. It requires static shapes, fixed addresses, no CPU syncs and no data-dependent control flow. In `torch.compile`, CUDAGraph trees re-record per distinct input shape, so dynamic sequence lengths need padding to a few buckets. Measured gains in published cases range from 1.12x (BERT at 4096 GPUs) to about 1.7x overall (Mask R-CNN).

### Cited Findings
- Mechanism and speedups: graphs launch many GPU operations through one CPU operation. Mask R-CNN gained 1.70x (272 GPUs) and BERT 1.12x (4096 GPUs) in MLPerf v1.0 training. To enlarge the graph scope NVIDIA removed CPU-GPU syncs, for example swapping `randperm`'s Thrust for CUB, and replaced dynamic-shape loss gathering with static tensors plus a mask. — [PyTorch blog](https://pytorch.org/blog/accelerating-pytorch-with-cuda-graphs/)
- PyTorch usage rules from the blog: warm up on a side stream, keep long-lived references to static input/output tensors, copy new data into them and replay. The whole step (forward, backward, optimizer) can be captured. `make_graphed_callables` graphs only the capture-safe parts when control flow, dynamic shapes or syncs prevent whole-network capture. Replay reuses the same memory addresses. The blog's advice is to try graphs when the network is graph-safe (static shapes and control flow) and "you suspect its runtime is at least somewhat CPU-limited". — [PyTorch blog](https://pytorch.org/blog/accelerating-pytorch-with-cuda-graphs/)
- PyTorch capture constraints (2.15 docs): capture must be on a non-default stream (raw API); only one capture at a time per process; ops that sync CPU with GPU (such as `.item()`) are prohibited; avoid `Generator.get_state/set_state` and use `graphsafe_get_state/graphsafe_set_state`; graphs may share a private memory pool only if always replayed in the same order, never concurrently. — [PyTorch CUDA semantics](https://docs.pytorch.org/docs/2.15/notes/cuda.html)
- CUDAGraph trees requirements (docs 2.8): same kernels with the same arguments (static sizes, addresses); same memory pattern at recording and replay; recordings invoked in a single order. — [CUDAGraph Trees (2.8)](https://docs.pytorch.org/docs/2.8/torch.compiler_cudagraph_trees.html)
- What breaks or skips them: "re-record CUDAGraph for every unique shape of an input tensor" (and frequent shape changes make it unprofitable, with padding to a few fixed shapes suggested); functions containing CPU operators are skipped; functions that mutate eager inputs in place are skipped ("skipping due to mutated inputs"), though mutating parameters and buffers is supported; multi-device operators and incompatible operators cause skips; control flow is not possible inside a graph; `.item()`-style syncs error; trees cannot persist memory between separate runs through the model. — [CUDAGraph Trees (2.8)](https://docs.pytorch.org/docs/2.8/torch.compiler_cudagraph_trees.html)
- Memory: all graphs share a single memory pool, so trees add no extra memory over eager in the common case. The legacy path (`cudagraph_trees=False`) can cause "large memory regressions". Each kernel launch inside a CUDA graph used 64 KB of device memory, which adds up with many re-recordings; the docs give this as true up to CUDA 12.4 and driver 550+ (a search summary of the docs, so confirm the exact version wording). — [CUDAGraph Trees (2.8)](https://docs.pytorch.org/docs/2.8/torch.compiler_cudagraph_trees.html); version detail via [search result summary of the stable docs](https://docs.pytorch.org/docs/stable/torch.compiler_cudagraph_trees.html)
- Input-mutation opt-in: `torch._inductor.config.cudagraph_support_input_mutation = True` for `reduce-overhead` mode (from the same docs, via search summary). — [CUDAGraph Trees docs](https://pytorch.org/docs/main/torch.compiler_cudagraph_trees)
- Tuning guide: enable via `torch.compile(m, "reduce-overhead")` or `"max-autotune"`; warns of higher device memory, possible compile failures and that support is still in development. — [PyTorch Tuning Guide](https://docs.pytorch.org/tutorials/recipes/recipes/tuning_guide.html)
- `reduce-overhead` mode "reduces framework overhead by a lot more" at the cost of a small amount of extra memory (PyTorch 2.0 launch page, 2023). — [PyTorch 2.x](https://pytorch.org/get-started/pytorch-2.0/)
- Graph updates (CUDA 11.3+ APIs, NVIDIA blog Aug 2022, A100): recapture plus `cudaGraphExecUpdate` gave 1.22x over no graph (580 ms vs 706 ms) and updating a stored kernel-node's parameters with `cudaGraphExecKernelNodeSetParams` gave 1.63x (433 ms). Reinstantiation is not needed if the topology is unchanged. — [NVIDIA blog](https://developer.nvidia.com/blog/constructing-cuda-graphs-with-dynamic-parameters/)
- Conditional nodes: IF and WHILE introduced in CUDA 12.4; IF-ELSE and SWITCH in CUDA 12.8 (post May 2024, updated Feb 2025). Body graphs may hold only kernel, empty, child-graph, memset, memcpy and conditional nodes; cooperative kernels are not supported in bodies. — [NVIDIA blog](https://developer.nvidia.com/blog/dynamic-control-flow-in-cuda-graphs-with-conditional-nodes/)
- Capture itself costs about 400 us to build a 20-kernel graph (V100, CUDA 10.1), amortized over replays. — [NVIDIA blog](https://developer.nvidia.com/blog/cuda-graphs/)

### Inferences
- Variable-length VLM batches (different image token counts, text lengths) make whole-step capture awkward. To use graphs you would pad to a handful of fixed shapes (for example sequence length buckets, fixed image resolution/token count) and accept padding waste. Each bucket costs one capture plus graph memory.
- Gradient checkpointing, the HF model's data-dependent Python (masks, `.item()`, `nonzero`, position-id construction in multimodal code) and TRL's loop are likely capture-hostile. Try partial graphing (`make_graphed_callables` on the decoder blocks, or `torch.compile` on the blocks) rather than the whole step. This is an inference; I found no source reporting CUDA Graph on a Qwen-VL LoRA run.
- Conditional nodes (CUDA 12.4/12.8) address GPU-side control flow in custom CUDA code. I found no evidence PyTorch exposes them for training loops.

### Gaps
- Official CUDA Programming Guide graph sections (capture restrictions, illegal calls, device graph launch) not retrieved.
- No source with a measured CUDA Graphs speedup for LoRA or VLM fine-tuning on RTX 4090.
- The docs.pytorch.org 2.15 CUDAGraph trees page returned empty; I used the 2.8 version, so newer behavior (dynamic shape handling, relaxed mutation rules) may differ.
- The exact 64 KB per-kernel memory statement and the `cudagraph_support_input_mutation` flag come from a search-result summary of the docs, not a full-page read.

## Kernel fusion: torch.compile, fused optimizers, Liger, FlashAttention

### Takeaway
Fusion reduces both kernel count (so CPU launch work) and memory traffic. torch.compile gave an average 43% training speedup on A100 across 163 models (21% fp32, 51% AMP), with smaller gains on desktop GPUs. Fused or foreach optimizers collapse the per-parameter loop. Liger reports about 20% throughput and 40-60% memory reductions on LLaMA-3-8B with 8x A100, not on LoRA or a 4090.

### Cited Findings
- PyTorch 2.0 launch benchmark (A100, training, 163 models: 46 HF, 61 TIMM, 56 TorchBench): works on 93% of models; average speedup 43% (weighted 0.75 AMP + 0.25 fp32), 21% in float32 and 51% in AMP. "Speedups are lower on desktop-class GPUs such as the NVIDIA 3090 than on server-class GPUs such as the A100." Inductor creates fusion groups. Dated 2023 and measured on PyTorch 2.0. — [PyTorch 2.x](https://pytorch.org/get-started/pytorch-2.0/)
- Hugging Face reports torch.compile gives "up to 30%" speedup for inference depending on model and GPU (for example Conditional-DETR batch 1 on A100 from about 31.2 ms to 18.0 ms). Vision models, inference. — [HF Transformers docs](https://huggingface.co/docs/transformers/perf_torch_compile)
- Optimizers have three implementations: for-loop, foreach (multi-tensor, horizontal fusion) and fused (a single kernel, vertical fusion on top). Speed order is fused > foreach > for-loop. Foreach is the default when available; fused is never the default. Fused supports Adam, AdamW, SGD, Adagrad (CPU only for Adagrad); fused AdamW is "stable" on CUDA. Fused is described as newer and less tested. — [PyTorch optim docs (2.8)](https://docs.pytorch.org/docs/2.8/optim.html)
- Liger-Kernel README: about 20% higher multi-GPU training throughput and about 60% lower memory; setup is LLaMA 3-8B, batch size 8, bf16, AdamW, gradient checkpointing, FSDP1 on 8x A100. Other examples: HF Trainer about 20% faster and over 40% memory reduction (4x A100); Lightning 15% and 40% (8x A100 ZeRO3); Medusa 40% and 80%. Kernels include RMSNorm, LayerNorm, SwiGLU, RoPE, CrossEntropy and FusedLinearCrossEntropy. Qwen2-VL and Qwen2.5-VL have `apply_liger_kernel_to_*` helpers. The README as fetched did not mention Qwen3-VL (it listed Qwen3.5 MoE multimodal) and did not mention LoRA. — [linkedin/Liger-Kernel](https://github.com/linkedin/Liger-Kernel)
- A secondary summary says Liger's fused linear cross entropy gives more than 2x speedup and more than 4x memory reduction for common vocab sizes, and that Liger scales more efficiently at larger batch sizes. Not independently verified. — [Liger-Kernel via search](https://github.com/linkedin/Liger-Kernel)
- "Chronicals" (arXiv 2601.02609) claims 3.51x over Unsloth for LLM fine-tuning. I could not read the body, so methodology, hardware and LoRA details are unknown. Treat as unverified single-author preprint. — [arXiv 2601.02609](https://arxiv.org/pdf/2601.02609)

### Inferences
- For LoRA, base weights are frozen, so the optimizer touches only adapter parameters; optimizer-launch savings come from fewer tensors, so they are modest versus full fine-tuning. The big per-step kernel counts come from transformer blocks (norms, rotary, activations, LoRA branches, casts), where Liger or torch.compile fusion would act. Savings are my reasoning, not measured.
- Qwen3-VL's text decoder follows the Qwen3 architecture; whether Liger's Qwen3-VL patch exists should be checked in the current repo before relying on it (README fetch did not list it).
- Gains reported on A100/8 GPUs at large batch (Liger) or in 2023 (compile) may not carry to a 4090 at micro-batch 4. The PyTorch note that desktop GPUs saw smaller compile gains supports caution.
- FlashAttention: no numbers collected here (see Gaps).

### Gaps
- FlashAttention paper/repo numbers were not retrieved. Use another researcher's notes or fetch the FlashAttention-2/3 papers.
- No LoRA-specific, 2B-VLM, or RTX 4090 fusion benchmark found. No Liger benchmark folder results read.
- torch.compile on gradient-checkpointed HF training with LoRA and recompilation behaviour with variable lengths: no source.
- No source for kernel-count reduction numbers (kernels per step before/after).

## Host-device synchronization points and how to find them

### Takeaway
Any operation that needs a GPU value on the CPU drains the queue and ends the CPU's run-ahead, so launch latency is exposed again after the sync. PyTorch's `torch.cuda.set_sync_debug_mode("warn" | "error")` flags many of these, with some gaps.

### Cited Findings
- PyTorch tuning guide lists sync triggers: `print(cuda_tensor)`, `.item()`, `.cpu()`/`.to(device)` copies, `nonzero()`, and Python branching on CUDA tensor values. It advises keeping the CPU ahead so the GPU queue stays full. — [PyTorch Tuning Guide](https://docs.pytorch.org/tutorials/recipes/recipes/tuning_guide.html)
- `torch.cuda.set_sync_debug_mode(debug_mode)`: "default"/0 ignores sync ops, "warn"/1 warns, "error"/2 raises. Experimental: not all synchronizing ops are covered, and `torch.distributed` and `torch.sparse` are not covered. — [torch.cuda.set_sync_debug_mode (2.15)](https://docs.pytorch.org/docs/2.15/generated/torch.cuda.set_sync_debug_mode.html)
- `optimizer.zero_grad(set_to_none=True)` avoids zeroing each gradient's memory and makes backward assign instead of accumulate; slightly different numerics. — [PyTorch Tuning Guide](https://docs.pytorch.org/tutorials/recipes/recipes/tuning_guide.html)
- Syncing after every launch took 9.6 us per kernel versus 3.8 us when launches overlap, on 2.9 us kernels (V100, CUDA 10.1). This is the cost of destroying CPU run-ahead. — [NVIDIA blog](https://developer.nvidia.com/blog/cuda-graphs/)
- NVIDIA removed syncs to widen graph scope: dynamic-shape tensors also force CPU-GPU sync via host-side memory management, and static shapes with masks avoid them. — [PyTorch blog](https://pytorch.org/blog/accelerating-pytorch-with-cuda-graphs/)
- Profiling cue: in an Nsight timeline the CPU at 100% with empty gaps between GPU kernels indicates launch-bound work. — [PyTorch blog](https://pytorch.org/blog/accelerating-pytorch-with-cuda-graphs/)

### Inferences
- A sync only hurts if the CPU is not already the bottleneck. When the loop is already CPU-bound (as here: main thread at 100%), the GPU queue is usually near-empty anyway, so removing syncs gives less than fixing the launch rate. Both matter; a sync after a long GPU phase stalls the CPU, while a sync in a CPU-bound loop mostly costs the tail.
- Likely sync candidates in a HF/TRL step (from general knowledge, not sourced here): `loss.item()` or logging each step, `clip_grad_norm_` when it returns a value that is inspected (the norm computation is async, but reading it is not), mask-dependent indexing and `nonzero` in multimodal code (placing image tokens via boolean masks), metric computation, and `torch.cuda.synchronize()` timers. Verify with `set_sync_debug_mode("warn")` plus the profiler.

### Gaps
- No source found specifically on `clip_grad_norm_` sync behavior or on HF Trainer/TRL logging cadence (`logging_steps`, `.item()` use).
- No numbers for the cost of a single sync in a training step.
- The PyTorch CUDA notes page was only read to character 100000 of about 164700; any `set_sync_debug_mode` text in the unread part was not seen.
