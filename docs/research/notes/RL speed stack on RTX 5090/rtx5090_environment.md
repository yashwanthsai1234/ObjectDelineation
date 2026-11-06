# Software environment for PyTorch training and vLLM inference on an RTX 5090 (sm_120) in a Docker image on Vast.ai

Research date: 2026-11-06. Most hardware figures come from secondary sources (the NVIDIA whitepaper itself was not fetched). Each claim is marked primary (NVIDIA, PyTorch, vLLM, Vast docs, GitHub) or secondary/community.

## 1. RTX 5090 specs and differences from datacenter Blackwell (sm_100)

### Takeaway
RTX 5090 (GB202, CC 12.0) is reported at 170 SMs, 21,760 CUDA cores, 680 5th-gen Tensor Cores, 32 GB GDDR7 on a 512-bit bus at about 1.79 TB/s. Per NVIDIA's CUDA programming guide, sm_120 has at most 100 KB shared memory per SM (vs 228 KB on sm_100) and is a separate compile target from sm_100. The sm_120 tensor cores accept FP8/FP6/FP4. Whether sm_120 lacks tcgen05/TMEM was only found in a community source.

### Cited Findings
- 21,760 CUDA cores, 170 SMs enabled of 192 on the die, 680 5th-gen Tensor Cores, 32 GB GDDR7, 512-bit, 1,792 GB/s, boost about 2.41 GHz (secondary spec sites agree) — [TechPowerUp](https://www.techpowerup.com/gpu-specs/geforce-rtx-5090.c4216), [Tom's Hardware review](https://www.tomshardware.com/pc-components/gpus/nvidia-geforce-rtx-5090-review)
- L2 cache: 96 MB per most sources; one listing says 88 MB and another 98 MB, so treat as unverified (secondary) — [TechPowerUp](https://www.techpowerup.com/gpu-specs/geforce-rtx-5090.c4216); conflicting: [PCGamesN](https://www.pcgamesn.com/nvidia/geforce-rtx-5090-release-date-price-specs-benchmarks)
- Peak FP4 tensor throughput from the RTX Blackwell whitepaper is quoted as 1676 dense / 3352 with sparsity (TFLOPS/TOPS); the 4090's "1321 AI TOPS" is an FP8 figure, so marketing numbers are not comparable across generations (secondary quoting the whitepaper) — [NeoGAF thread reproducing the whitepaper](https://www.neogaf.com/threads/rtx-blackwell-whitepaper-released.1680199/), [Tom's Hardware](https://www.tomshardware.com/pc-components/gpus/nvidia-announces-rtx-50-series-at-up-to-usd1-999)
- RTX Blackwell adds FP4 and FP6 Tensor Core operations and a second-generation FP8 Transformer Engine — [NeoGAF thread reproducing the whitepaper](https://www.neogaf.com/threads/rtx-blackwell-whitepaper-released.1680199/)
- NVIDIA CUDA guide, CC 12.x: unified data cache 128 KB with up to 100 KB shared memory per SM (CC 10.0: 256 KB cache, up to 228 KB shared); 99 KB opt-in per block; Tensor Core types include FP64, TF32, BF16, FP16, FP8, FP6, FP4, INT8; TMA, thread block clusters, distributed shared memory are listed as supported. The page does not mention tensor memory (TMEM) for 12.x. The fetch tool flagged table-column ambiguity, so verify against the rendered page — [NVIDIA CUDA compute capabilities](https://docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/compute-capabilities.html)
- Arch-specific target for 12.x is compute_120a; family target compute_120f covers 12.0 and 12.1 (NVIDIA guide, same URL above)
- The sm_100f / sm_120f family targets do not run on each other's devices; sm_100 cubins are not compatible with sm_120 (community guide) — [arnon.dk arch/gencode guide](https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/)
- Community claim (single blog, not NVIDIA): ptxas for the sm_120a target rejects tcgen05 instructions, so datacenter sm_100 kernels (e.g., FA4 as designed for B200) do not compile for sm_120 — [gist](https://gist.github.com/solatticus/aab6ec3a0436748b021cbbdd12e8c739)
- FlashAttention-4 status is contested. A September 2026 compilation says upstream flash-attention setup.py targets sm_120 with CUDA 12.8+, forward PR #2329 is merged, and an SM120 backward is committed; FA4's headline speedups are B200 numbers only — [keenable compilation, 2026-09-01](https://select.keenable.ai/r/sm120-attention-for-encoder-models/yWskIdfUzZDXx2FSsDhXIQ). A community Hugging Face package bundles six open upstream PRs for consumer Blackwell — [flash-attn-4-sm120](https://huggingface.co/SecondNatureComputing/flash-attn-4-sm120). Older (about March 2026) reports say FA4 was SM80/90/100/110 only — [SGLang discussion](https://github.com/sgl-project/sglang/discussions/10564), [gist](https://gist.github.com/solatticus/aab6ec3a0436748b021cbbdd12e8c739)
- PyTorch SDPA with the cuDNN backend reached 203.61 BF16 TFLOP/s on an RTX 5090 (one shape, B=1,H=8,Q=4096,KV=8192), described as 97.19% of a roofline; the roofline (about 209.5) is back-computed, and this is one forward shape only — [keenable compilation](https://select.keenable.ai/r/sm120-attention-for-encoder-models/yWskIdfUzZDXx2FSsDhXIQ)

### Inferences
- Kernels written for sm_100 (tcgen05, TMEM, 228 KB smem tiles) will not run on the 5090; sm_90/sm_80-style (mma.sync, cp.async, TMA, <=99 KB smem) code is what to expect. Our own CUDA checker should avoid assuming more than 99 KB dynamic smem per block.
- Official NVIDIA dense BF16/FP8 TFLOPS for the 5090 were not found; do not hard-code a roofline. Measure it with the matmul check (section 8).

### Gaps
- Official NVIDIA numbers for dense BF16/FP8 TFLOPS, L2 size and SM count (whitepaper PDF not fetched; spec sites are secondary).
- Primary confirmation that sm_120 has no tcgen05/TMEM (only the community ptxas observation above).

## 2. Driver and CUDA minimums; PyTorch wheels and sm_120; Triton

### Takeaway
sm_120 needs CUDA 12.8+ (driver >= 570.26 for 12.8 GA); CUDA 13.x needs driver >= 580. PyTorch 2.13.0 (released 2026-07-08) makes cu130 the default, dropped the 12.8/12.9 builds, and pins Triton 3.7.1. Whether the torch 2.13 cu130 wheel contains sm_120 SASS or PTX only was not found; check it with `torch.cuda.get_arch_list()`.

### Cited Findings
- CUDA 12.8 GA requires Linux driver >= 570.26 (update 1: >= 570.124.06); CUDA 13.x requires driver >= 580 for minor-version compatibility (CUDA 13.0 = R580 branch) — [NVIDIA CUDA release notes](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html)
- Blackwell GPUs need a PyTorch build made with CUDA 12.8 or newer (reply from ptrblck, 2025-10-08) — [PyTorch forum](https://discuss.pytorch.org/t/is-there-a-pytorch-build-that-supports-nvidia-rtx-5090-compute-capability-12-0-sm-120/223536.md)
- PyTorch 2.13.0 GA on 2026-07-08; CUDA 13.0 remains the default build; the CUDA 12.8/12.9 builds were removed; small wheels are always built for CUDA Linux; ptxas is no longer bundled in the cu13 binary; Triton pin advanced to 3.7.1; Inductor gains a CuTeDSL backend for GEMM/RMSNorm — [PyTorch 2.13 release blog](https://pytorch.org/blog/pytorch-2-13-release-blog/)
- The release blog does not mention sm_120 or give the arch list (same source); release-branch commits mention stripping "+PTX" from the CUDA arch list on release builds (search result summary, unverified) — [PyTorch 2.13 GA post](https://dev-discuss.pytorch.org/t/pytorch-2-13-0-general-availability/3412.md)
- CUDA 13.0 removed offline compilation support for Maxwell, Pascal, Volta (cuFFT 13.0 note says so for libraries) — [NVIDIA CUDA release notes](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html)
- A community build for torch 2.9.1 + CUDA 13 lists archs 7.5 to 12.1 including 12.0; this is not the official wheel — [search result, third-party](https://forums.developer.nvidia.com/t/request-for-official-cuda-and-pytorch-support-for-rtx-5080-5090-sm-120/327876)
- A cuBLAS known issue affecting CC 12.0/12.1 was introduced in CUDA 13.0 Update 2 (details not read) — [NVIDIA CUDA release notes](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html)

### Inferences
- cu128+ wheels have shipped sm_120 since 2025, and the 2.13 cu130 default almost certainly includes sm_120, but this is not confirmed here. The preflight should assert `"sm_120" in torch.cuda.get_arch_list()` (or a compute_120 PTX entry) rather than assume.
- "ptxas not bundled in cu13 binary" concerns PyTorch itself; Triton ships its own ptxas, but confirm on a 5090 that a `torch.compile` / Triton kernel compiles with Triton 3.7.1 (smoke test).
- Host driver must be >= 580 for cu130 wheels (pip cu13 runtime). Older Vast hosts on 570.x would work only for cu128-class runtimes, so filter on driver_version.

### Gaps
- Official arch list of torch 2.13.0+cu130 (SASS vs PTX for sm_120).
- Triton 3.7.1 sm_120 specific notes; no source found.
- Exact content of the cuBLAS CC 12.0 known issue.

## 3. vLLM on RTX 5090

### Takeaway
Prebuilt vLLM wheels failed on RTX 50-series ("no kernel image") until a fix closed 2026-05-15; the reasons are build arch lists. Backend choice on sm_120 is not documented in a priority table; community reports favor FlashInfer, with FA2 and Triton as fallbacks. vLLM 0.32.0 specifics for sm_120 were not found.

### Cited Findings
- Issue #35432 (opened 2026-02-26): prebuilt PyPI wheels and Docker images failed on RTX 50xx with "no kernel image is available"; a contributor said Blackwell needs CUDA 13.x builds; closed as completed 2026-05-15 via PR #35568 (content of the PR not read) — [vLLM #35432](https://github.com/vllm-project/vllm/issues/35432)
- Forum thread: vLLM 0.9.0+ added support but images must be built with torch_cuda_arch_list="12.0 12.1" (search summary) — [vLLM forum](https://discuss.vllm.ai/t/cuda-error-no-kernel-image-is-available-for-execution-on-the-device/1475)
- vLLM docs backend page has no sm_120 priority list. TRITON_ATTN supports any compute capability; FA2 requires >=8.0; FA4 requires >=10.0 and is default on "SM100+"; B12X backend (install `vllm[b12x]`, `--attention-backend b12x`) supports causal decoder attention on SM120/SM121, bf16 head sizes 128/192/256; no FlashInfer entry names 12.x; automatic selection picks the first compatible backend, and a forced `--attention-backend` is validated — [vLLM attention backends doc](https://docs.vllm.ai/en/latest/design/attention_backends/)
- 2025 working guide (vLLM source build, torch 2.9.0 nightly cu128, driver 575.64.03): `VLLM_FLASH_ATTN_VERSION=2` because FA3 is not supported on Blackwell; `TORCH_CUDA_ARCH_LIST="12.0"`; Qwen2.5-7B about 290 tok/s decode (user-reported, outdated) — [vLLM forum](https://discuss.vllm.ai/t/vllm-on-rtx5090-working-gpu-setup-with-torch-2-9-0-cu128/1492)
- Issue #37714: explicit `--attention-backend FLASH_ATTN` with `--kv-cache-dtype fp8` crashes at engine start on RTX 5090; users often pick FLASH_ATTN as a workaround for FlashInfer crashes (search summary, issue page not read) — [vLLM #37714](https://github.com/vllm-project/vllm/issues/37714)
- Issue #47749 (2026-07-06, vLLM 0.24.0, torch 2.11.0+cu130, flashinfer 0.6.12, driver 610.43.02): NVFP4 W4A16 on 5090 goes through Marlin with a misleading "no native FP4" warning; maintainers say Marlin for weight-only NVFP4 is by design, W4A4 uses native CUTLASS on SM120; FlashInfer autotune showed OOM fallback warnings on a 32 GB card; stale/open — [vLLM #47749](https://github.com/vllm-project/vllm/issues/47749)
- A 2025/2026 forum report of a segfault in CUDA graph capture and a deadlock under about 30 concurrent users on 5090 (vLLM 0.11.2 docker, driver 580.95.05, Qwen3-8B + LoRA with bnb); no confirmed fix, only an automated unverified reply — [vLLM forum, posted 2026-06-01](https://discuss.vllm.ai/t/bug-segfault-in-pythonsymnodeimpl-and-deadlock-on-rtx-5090-blackwell-with-vllm-0-11-2/2688)
- RTX PRO 6000 (sm_120) whole-host hard lock during vLLM CUDA-graph capture (vLLM 0.25.1, FlashInfer 0.6.13, Aug 2026); locking clocks (`nvidia-smi -lgc 180,2800`) or disabling CUDA graphs avoided it; different card from the 5090 — [NVIDIA forum](https://forums.developer.nvidia.com/t/rtx-pro-6000-blackwell-ws-whole-host-hard-lock-during-vllm-cuda-graph/379239)
- On RTX PRO 6000 (SM120), the FlashInfer attention backend is reported stable with FP8 KV cache, while some NVFP4 MoE CUTLASS paths gave bad output — [vLLM forum perf report](https://discuss.vllm.ai/t/sm120-rtx-pro-6000-nvfp4-moe-performance-report-qwen3-5-397b/2536)
- Vendor Docker description (not independent): standard flash-attn gave undefined-symbol errors on SM_120 so they use FlashInfer — [Docker Hub](https://hub.docker.com/r/malkaf/vllm-blackwell-optimizer)

### Inferences
- With `VLLM_USE_FLASHINFER_SAMPLER=0`, only the sampler avoids FlashInfer JIT. The FlashInfer attention backend on sm_120 may still JIT-compile (needs nvcc), unknown for 0.32.0. If the image has no nvcc, pin `--attention-backend` to FLASH_ATTN (FA2 path) or TRITON_ATTN and test, rather than rely on auto-selection. This is untested and not confirmed by any source.
- As a fallback for a crash in graph capture, use `enforce_eager=True`; cost is CPU launch overhead.
- FP8 KV cache + FLASH_ATTN is a reported crash combination; do not combine them.

### Gaps
- vLLM 0.32.0 release notes and sm_120 behavior (all found issues are on 0.11 to 0.25).
- Which backend vLLM 0.32.0 auto-selects on CC 12.0, and whether FlashInfer attention needs nvcc JIT there.
- Official list of recommended env vars for 5090.

## 4. Compiling our own CUDA C++ extension for sm_120 (pip nvcc, CMake)

### Takeaway
Target `compute_120`/`sm_120` (CUDA 12.8+); use 120a only if the kernel needs arch-specific features (not needed for ordinary mma/cp.async code); 120f (CUDA 12.9+) covers 12.0 and 12.1. The pip-nvcc and CMake specifics were not found in sources.

### Cited Findings
- 12.x targets: compute_120a (arch-specific), compute_120f (family; superset of baseline compute_120, covers 12.0 and 12.1) — [NVIDIA CUDA guide](https://docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/compute-capabilities.html)
- sm_120f was added in CUDA 12.9; sm_120a is not forward compatible; one forum user observed ptxas assuming different shared memory (49 kB vs 102 kB) for sm_120 vs sm_120a with no definitive answer — [NVIDIA forum](https://forums.developer.nvidia.com/t/choosing-the-right-arch-flag-for-nvcc/362941), [arnon.dk](https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/)
- LLVM's NVPTX table maps sm_120 to compute_120 with CUDA 12.8 as the minimum — [LLVM NVPTX target table, via search](https://mail-archive.com/cfe-commits@lists.llvm.org/msg723575.html)
- Community build setting for vLLM: `CMAKE_CUDA_ARCHITECTURES=120`, `TORCH_CUDA_ARCH_LIST="12.0+PTX"` (unconfirmed in thread) — [vLLM forum](https://discuss.vllm.ai/t/vllm-on-rtx5090-working-gpu-setup-with-torch-2-9-0-cu128/1492)

### Inferences
- For our checker, `-DCMAKE_CUDA_ARCHITECTURES=120` (SASS for sm_120) is the plain choice; adding `120-virtual` would add PTX. `native` requires a GPU visible at build time. CMake accepting the `120a` suffix was not verified (CMake 3.23+ documents the "a" suffix in general, from my knowledge, not from a fetched source).
- `nvidia-cuda-nvcc`/`cuda-toolkit[nvcc,cccl]` 13.0.x ships ptxas; CUDA 13 needs driver >= 580.

### Gaps
- Behavior of the pip `nvidia-cuda-nvcc` wheels (include paths, missing headers such as cuda_runtime.h, libcudart) for building extensions: no source found. Test with `nvcc --list-gpu-arch` and a compile of a trivial kernel in the image build.
- Official nvcc documentation for sm_120a vs sm_120 differences.

## 5. Vast.ai RTX 5090 hosts

### Takeaway
Prices are aggregator-reported at about $0.40 to $0.75 per hour (median about $0.60); there are no primary or current numbers for US/CA, and no source gives typical host driver versions. Use `driver_version`, `cuda_vers`, `geolocation`, `reliability`, `pcie_bw` filters.

### Cited Findings
- CLI query syntax: `vastai search offers 'FIELD=VALUE ...' -o 'dlperf_usd-'`; documented fields include `cuda_vers` (max CUDA supported by the driver; the docs list this, not `cuda_max_good`), `driver_version` (string like "535.86.05"), `geolocation` (two-letter country code, supports `in`/`notin`), `reliability`, `pcie_bw`, `cpu_arch`, `inet_down`, `disk_space`, `gpu_ram`, `dph` — [Vast CLI reference](https://docs.vast.ai/cli/reference/search-offers.md)
- Example filter: `'disk_space>146 duration>24 gpu_ram>10 cuda_vers>=12.1 direct_port_count>=2 driver_version >= 535.86.05'`; `gpu_name` example `RTX_4090` in the quick-start, with `"RTX 4090"` in lists — [Vast CLI reference](https://docs.vast.ai/cli/reference/search-offers.md), [Vast search docs](https://docs.vast.ai/documentation/instances/search)
- Machine tiers: Unverified (filtered by default), Verified, Secure Cloud; new machines start at 60% reliability — [Vast search docs](https://docs.vast.ai/documentation/instances/search)
- Price snapshots: $0.45 to $0.75 per hour, median $0.60, with a 3% commission in listed prices; some hosts use older Xeons that can bottleneck; seller-terminated instances can end abruptly (article dated 2026-02-27, aggregator) — [DeployBase](https://deploybase.ai/articles/rtx-5090-on-vastai-pricing-specs-how-to-rent). Other trackers range from $0.14 spot to $0.63 median and disagree — [Beam](https://www.beam.cloud/blog/nvidia-rtx-5090-pricing), [getdeploying](https://getdeploying.com/gpus/nvidia-rtx-5090), [gpuperhour](https://gpuperhour.com/providers/vastai/rtx-5090)
- Host CPUs seen: Intel i9, Ryzen 9 7950X, some older Xeons (aggregator) — [DeployBase](https://deploybase.ai/articles/rtx-5090-on-vastai-pricing-specs-how-to-rent)

### Inferences
- A suggested filter (not tested): `gpu_name=RTX_5090 num_gpus=1 geolocation in [US,CA] driver_version>=580.65.06 reliability>0.98 pcie_bw>=20 rentable=true verified=true`. The 580.65.06 value is my assumption (the R580 branch's first CUDA 13.0 driver); the docs only say >= 580. Verify field spelling with `vastai search offers --help`.
- A host with an older Xeon or PCIe 4.0 x8 hurts small-model workloads more (see section 7). Require `pcie_bw` and avoid hosts with low CPU clocks.

### Gaps
- Real driver / cuda_vers distribution of 5090 offers in the US/CA (needs a live `vastai search offers`).
- Whether Vast hosts commonly run PCIe 5.0 for the 5090; no source.
- Primary Vast pricing page for 2026.

## 6. Performance expectations vs the 4090 (small models)

### Takeaway
Published numbers (all secondary or vendor) show about 1.3 to 1.4x on small-model inference, tracking the bandwidth ratio (1792 vs 1008 GB/s), and about 1.3 to 1.5x on BF16 training. No source measured CPU launch overhead on a 5090 directly.

### Cited Findings
- Llama 3.1 8B FP16: 3,500 tok/s on 5090 vs 2,550 on 4090 (+37%) — [Spheron](https://www.spheron.network/blog/rtx-5090-vs-rtx-4090/); another vendor claims 7,198 vs about 4,200 to 4,800 with unclear setup (lower weight) — [Fluence](https://www.fluence.network/blog/rtx-5090-vs-4090/)
- Fleet data: 5090 about 36% faster than 4090 on average across 9 models — [Trooper](https://www.trooper.ai/benchmarks-compare-RTX-5090-with-RTX-4090-Pro)
- Training: about 1.5x for a translation job and about 1.3x for an 8B LLM on 4 GPUs (VNG Cloud); BF16 training +35 to 45% and FP8 about 1.7x vs BF16 on 5090 (Dev.to, self-reported) — [VNG Cloud](https://vngcloud.vn/blog/nvidia-rtx-5090-gpu-toi-uu-cho-inference-training-llm), [Dev.to](https://dev.to/thurmon_demich/rtx-5090-for-ai-in-2026-6-month-honest-retrospective-4e0g)
- Jarvislabs gives an estimate (not measured) of 1.3 to 1.7x — [Jarvislabs](https://jarvislabs.ai/ai-faqs/nvidia-rtx-5090-specs)

### Inferences
- For decode-heavy small-model work (batch small, kernels tiny), per-step time is bound by kernel launches and Python/CPU time, so the realized gain can fall below the bandwidth ratio. CUDA graphs and a fast host CPU matter more than the GPU swap. This is reasoning, not a measured result.
- Per-dollar: at about $0.60 per hour (5090) vs a typical 4090 price, a 1.3x speedup breaks even only if the 5090 premium is under about 30%. 4090 prices were not gathered, so compute this live.

### Gaps
- Rigorous benchmarks for <=3B models, vLLM with CUDA graphs, or CPU-launch-bound loops on a 5090.
- Independent confirmation of the benchmark figures above.

## 7. Known bugs on sm_120 summary (cross-reference)

### Takeaway
Reported problems are mostly about arch lists in builds, FA4/FA3 availability, NVFP4 paths, and CUDA graph capture hangs on SM120 hardware; none was confirmed for vLLM 0.32.0 + torch 2.13.

### Cited Findings
- See section 3 for vLLM issues #35432, #37714, #47749, the 0.11.2 segfault/deadlock thread, and the PRO 6000 CUDA-graph lock.
- "RuntimeError: Int8 not supported on SM120" was a forum topic dated 2025-11-19 (title only, not read) — [vLLM forum topic list](https://discuss.vllm.ai/t/bug-segfault-in-pythonsymnodeimpl-and-deadlock-on-rtx-5090-blackwell-with-vllm-0-11-2/2688)
- NVFP4 MoE on SM120 produced bad output with some CUTLASS paths, which the report says needs a CUTLASS fix — [vLLM forum](https://discuss.vllm.ai/t/sm120-rtx-pro-6000-nvfp4-moe-performance-report-qwen3-5-397b/2536)

### Inferences
- We use bf16 only, so NVFP4/FP8-KV/int8 bugs are outside our path. Keep to bf16 weights, bf16 KV, and no bitsandbytes.

### Gaps
- Open bugs against vLLM 0.32.0 on sm_120 (search the vLLM tracker with label/term "sm120" at rental time).

## 8. Start-of-rental checklist

### Takeaway
Run all of these before downloading the bundle; each maps to a failure found above. The commands are my proposals built from the cited requirements, not from a single source.

### Cited Findings
- Driver >= 580 is required for CUDA 13.x; sm_120 needs CUDA >= 12.8 — [NVIDIA CUDA release notes](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html)
- Blackwell consumer needs torch builds with CUDA >= 12.8 — [PyTorch forum](https://discuss.pytorch.org/t/is-there-a-pytorch-build-that-supports-nvidia-rtx-5090-compute-capability-12-0-sm-120/223536.md)
- Clock policy affected CUDA-graph stability on an SM120 card — [NVIDIA forum](https://forums.developer.nvidia.com/t/rtx-pro-6000-blackwell-ws-whole-host-hard-lock-during-vllm-cuda-graph/379239)

### Inferences
Proposed checks (in order):
1. `nvidia-smi --query-gpu=name,driver_version,memory.total,compute_cap,pcie.link.gen.current,pcie.link.width.current,clocks.max.sm,power.limit --format=csv`: expect RTX 5090, driver >= 580, about 32 GB (about 32607 MiB reported), compute_cap 12.0, PCIe gen 5 (or 4) x16 under load.
2. `python -I -c "import torch;print(torch.__version__,torch.version.cuda,torch.cuda.get_device_capability(),torch.cuda.get_arch_list())"`: expect (12, 0) and `sm_120` (or compute_120) in the list; fail early if absent.
3. bf16 matmul TFLOPS: 8192x8192 bf16 matmul, warm up, time with CUDA events, report TFLOPS. No official dense roofline was found; record the number and compare to a stored 4090/5090 baseline from the first rental. Also record `nvidia-smi -q -d CLOCK` during it.
4. Triton smoke: a `torch.compile` of a small function on cuda (confirms Triton 3.7.1 + its ptxas work on sm_120).
5. nvcc: `nvcc --version; nvcc --list-gpu-arch | grep 120`, then compile a trivial kernel with `-gencode arch=compute_120,code=sm_120` and run it.
6. vLLM: `import vllm`, then a tiny engine run (bf16, small model, `max_model_len` small) first with `enforce_eager=True`, then with CUDA graphs; log which attention backend was selected; keep `VLLM_USE_FLASHINFER_SAMPLER=0`.
7. Host: record `lscpu` model name and `nproc`; reject old low-clock Xeons for launch-bound jobs.
8. Optional: hold clocks (`nvidia-smi -lgc`) only if CUDA graph hangs appear (reported on PRO 6000, not shown for 5090).

### Gaps
- No primary source for expected TFLOPS thresholds; thresholds must come from our own first-run baseline.
