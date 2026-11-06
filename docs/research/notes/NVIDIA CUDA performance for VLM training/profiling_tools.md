# Profiling a PyTorch GPU training job and a vLLM inference job (Nsight Systems, Nsight Compute, torch.profiler, DCGM, nvidia-smi)

Doc versions seen: Nsight Systems v2026.6, Nsight Compute 2026.4.0, PyTorch docs 2.15 (tutorials need torch >= 2.3), vLLM `--profiler-config` needs >= 0.13.0 (so it applies to your 0.32). Not verified by me: exact Nsight versions that ship with CUDA 13 images.

## Nsight Systems (nsys): command lines, stats, overhead

### Takeaway
nsys is the first tool to run: it shows the timeline (GPU idle gaps, kernel vs API rows) and gives CLI summaries. Use `--capture-range=cudaProfilerApi` or `--delay/--duration` to capture only a few steps. File sizes and overhead numbers were not found in official docs.

### Cited Findings
- `--trace/-t` defaults to cuda, opengl, nvtx, osrt; other values include cudnn, cublas, nccl, python-gil; `none` disables API tracing — [Nsight Systems User Guide v2026.6](https://docs.nvidia.com/nsight-systems/UserGuide/index.html)
- `--capture-range` values: none (default), cudaProfilerApi, hotkey (graphics apps only), nvtx. `--capture-range-end` default is `stop-shutdown`; also none, stop, `repeat[:N]`, `repeat-shutdown:N` — [User Guide](https://docs.nvidia.com/nsight-systems/UserGuide/index.html)
- `--duration/-d` is collection seconds (>0); `--delay/-y` is start delay seconds (default 0) — [User Guide](https://docs.nvidia.com/nsight-systems/UserGuide/index.html)
- `--cuda-memory-usage` default false and "may add significant overhead" — [User Guide](https://docs.nvidia.com/nsight-systems/UserGuide/index.html)
- `--python-sampling` default false (x86 Linux supported), default 1000 Hz (1-2000) — [User Guide](https://docs.nvidia.com/nsight-systems/UserGuide/index.html)
- `--pytorch` adds automatic PyTorch annotations: autograd-nvtx, autograd-shapes-nvtx, functions-trace, functions-trace-shapes (default none) — [User Guide](https://docs.nvidia.com/nsight-systems/UserGuide/index.html)
- `--run-as` and `--ftrace` need root; some `--cpuctxsw` values may need admin/root; `--sample=xhv` needs root — [User Guide](https://docs.nvidia.com/nsight-systems/UserGuide/index.html)
- Built-in `nsys stats` reports: cuda_api_sum, cuda_api_trace, cuda_api_gpu_sum, cuda_gpu_kern_sum, cuda_gpu_kern_gb_sum (also groups by grid/block), cuda_gpu_mem_time_sum, cuda_gpu_mem_size_sum, cuda_gpu_sum, cuda_gpu_trace, cuda_kern_exec_sum (launch API vs queue vs execution time), cuda_kern_exec_trace, nvtx_sum, nvtx_gpu_proj_sum (NVTX ranges projected onto GPU work), osrt_sum. Kernel reports take suffixes `nvtx-name`, `base`, `mangled`. Formatters: column, table, csv, tsv, json, hdoc, htable — [Post-Collection Analysis Guide](https://docs.nvidia.com/nsight-systems/AnalysisGuide/index.html)
- `nsys analyze` runs expert-system rules: sync memcpy/memset, pageable-memory async memcpy, sync APIs; "GPU Starvation" finds idle gaps longer than 500 ms (adjustable); "GPU Low Utilization" flags time chunks under a threshold (time utilization, not resource utilization) — [Analysis Guide](https://docs.nvidia.com/nsight-systems/AnalysisGuide/index.html)
- `nsys stats --help-reports` lists reports for your installed version — [NVIDIA forum](https://forums.developer.nvidia.com/t/i-want-to-get-csv-file-about-kernel-launching-calling-to-kernel/192126)
- No standalone `cuda_gpu_gaps` report was found; one secondary source says canned reports cannot give idle percentage/gaps and need SQLite export — [tomevault summary](https://tomevault.io/tome/Snowfall99/nsys-report-skill) (low-quality secondary source)

### Inferences
- Suggested training command (my synthesis from the flags above, untested): `nsys profile -t cuda,nvtx,osrt,cudnn,cublas -s none --capture-range=cudaProfilerApi --capture-range-end=stop -o prof python -I -m ...`, with `torch.cuda.profiler.start()/stop()` around ~5 steps after warmup. `-s none` (CPU sampling off) is the safe choice without perf access (see the container section).
- Then run `nsys stats -r cuda_gpu_kern_sum,cuda_api_sum,cuda_kern_exec_sum,nvtx_sum,nvtx_gpu_proj_sum prof.nsys-rep` and `nsys analyze prof.nsys-rep`. For gap percentages, export to SQLite and query (no canned report found).
- Reading the timeline: big idle gaps in the GPU kernel row with a busy CPU thread = launch-bound (many tiny kernels; large `cudaLaunchKernel` share, large API-to-start latency in cuda_kern_exec_sum). Gaps preceded by a long `cudaStreamSynchronize`/`cudaMemcpy` (`.item()`, `.cpu()`) = synchronization-bound. Gaps while the main thread waits in DataLoader/`osrt` calls (poll/read/futex) = data-loader-bound. Dense kernels with no gaps = compute- or memory-bound; then use ncu on top kernels.

### Gaps
- File sizes and overhead percentages: not found in fetched official text (only a warning on `--cuda-memory-usage`). Treat as measure-yourself; keep captures to a few steps.
- Only the first ~100k characters of the 500k-character User Guide and Analysis Guide were read; flag details (e.g. `--stats`, `--force-export`, `nsys stats -r`) were not verified from the page. Check `nsys stats --help`.

## Running nsys / ncu in Docker without host root (Vast.ai)

### Takeaway
CPU sampling and counters depend on host settings; the documented container fixes (`--cap-add=SYS_ADMIN`, seccomp for `perf_event_open`) need docker-run control, which Vast docs do not show. nsys with CUDA/NVTX tracing and CPU sampling off is the likely-working fallback, but no source confirms it on Vast.

### Cited Findings
- Counter restriction (ERR_NVGPUCTRPERM) comes from NVIDIA's 2018 side-channel response (Linux driver 418.43+); `NVreg_RestrictProfilingToAdminUsers=0` in /etc/modprobe.d opens it to all users, `=1` restricts to CAP_SYS_ADMIN/CAP_PERFMON holders; check `RmProfilingAdminOnly` in `/proc/driver/nvidia/params`; may need reboot/module reload — [NVIDIA ERR_NVGPUCTRPERM](https://developer.nvidia.com/ERR_NVGPUCTRPERM)
- In a container, profiling needs host-side enablement or starting the container with `--cap-add=SYS_ADMIN` — [NVIDIA ERR_NVGPUCTRPERM](https://developer.nvidia.com/ERR_NVGPUCTRPERM)
- CPU sampling/scheduling needs perf paranoid level <= 2 on the host; a `--privileged` container with `sudo nsys` still failed at paranoid 4 — [NVIDIA forum](https://forums.developer.nvidia.com/t/cpu-sampling-in-privileged-docker-container-via-sudo-nsys/218361)
- NVIDIA's container guide suggests `--cap-add=SYS_ADMIN` or a seccomp profile allowing `perf_event_open`, and `nsys status -e` to check the environment — [NVIDIA blog (Korean)](https://developer.nvidia.com/ko-kr/blog/%EC%BB%A8%ED%85%8C%EC%9D%B4%EB%84%88-%EB%B0%8F-%ED%81%B4%EB%9D%BC%EC%9A%B0%EB%93%9C%EC%97%90%EC%84%9C-nvidia-nsight-%EC%8B%9C%EC%8A%A4%ED%85%9C-%EC%82%AC%EC%9A%A9%ED%95%98%EA%B8%B0)
- A user without sampling support in a non-privileged container used `nsys profile -s none` — [NVIDIA forum](https://forums.developer.nvidia.com/t/nsight-systems-initialization-errors/181287)
- ncu/nsys GPU metrics can conflict with DCGM running on the same host (shared counters; "Already under profiling") — [NVIDIA forum](https://forums.developer.nvidia.com/t/nsys-profile-gpu-metrics-devices-fails-with-already-under-profiling/333938)
- Vast docs describe template/CLI control of docker-create options (env, hostname, ports) with no mention of privileged/capabilities — [Vast docker environment](https://docs.vast.ai/documentation/instances/docker-environment)

### Inferences
- Plan: run `nsys status -e` and `cat /proc/driver/nvidia/params | grep RmProfilingAdminOnly` and `cat /proc/sys/kernel/perf_event_paranoid` in the preflight on the rented box. Use `nsys -s none` (and never rely on `--gpu-metrics-devices`, which uses counters). Treat ncu as likely unavailable; use torch.profiler plus nsys instead. Test once; I found no Vast-specific confirmation.

### Gaps
- Whether Vast allows `--cap-add=SYS_ADMIN`/`--privileged`, and whether its hosts set RestrictProfilingToAdminUsers: unanswered (Vast pages fetched returned 404 or lacked it; the NVIDIA English container blog URL 404'd).
- Whether nsys CUDA tracing works unprivileged: no authoritative statement found.

## Nsight Compute (ncu)

### Takeaway
ncu gives per-kernel SM/memory throughput, occupancy and roofline, but needs GPU counter permission (same ERR_NVGPUCTRPERM issue) and replays kernels, so profile only a few chosen kernels.

### Cited Findings
- Speed of Light section = high-level compute and memory throughput as % of peak (e.g. `sm__throughput`); Roofline sections exist for half/single/double/tensor precision; Occupancy = active warps / max warps, and a large theoretical-vs-achieved gap usually means imbalanced work — [Nsight Compute Profiling Guide 2026.4.0](https://docs.nvidia.com/nsight-compute/ProfilingGuide/index.html)
- No `--set/--section/--metrics` gives the `basic` set; `--set full` collects everything; `--list-sets`, `--list-sections`, `--query-metrics`, `--launch-count`, `--replay-mode`, `--cache-control`, `--clock-control` exist — [Profiling Guide](https://docs.nvidia.com/nsight-compute/ProfilingGuide/index.html)
- Replay: metrics that don't fit in one pass cause multiple passes (kernel replay saves/restores memory; application replay reruns the app); more metrics = more passes; software-patched metrics cost the most — [Profiling Guide](https://docs.nvidia.com/nsight-compute/ProfilingGuide/index.html)
- Permission error and fixes: see container section.

### Inferences
- Typical use (flags `--kernel-name`, `--launch-skip`, `-o` are standard ncu flags but were not in the excerpt I read): `ncu --set full -k regex:<name> -s <skip> -c 5 -o out python ...`. Reading: memory throughput high, SM throughput low = memory-bound; SM/tensor high = compute-bound; both low with low occupancy = latency/launch or small-problem issues.

### Gaps
- Exact kernel-filter flag docs, SM vs memory throughput metric names, and whether ncu works on Vast: not verified.

## torch.profiler and NVTX

### Takeaway
torch.profiler with a wait/warmup/active schedule gives CPU-op vs CUDA-kernel tables and a Chrome/Perfetto trace without special permissions (it uses CUPTI via Kineto; no counter privileges were documented as needed). Use NVTX ranges to label phases for nsys.

### Cited Findings
- Signature: `profile(activities, schedule, on_trace_ready, record_shapes=False, profile_memory=False, with_stack=False, with_flops=False, ...)`; activities default to CPU plus CUDA when available — [PyTorch 2.15 profiler docs](https://docs.pytorch.org/docs/2.15/profiler.html)
- `schedule(wait, warmup, active, repeat=0, skip_first=0, skip_first_wait=0)`; warmup traces but discards; `repeat=0` cycles forever — [PyTorch 2.15 profiler docs](https://docs.pytorch.org/docs/2.15/profiler.html), [Profiler recipe](https://docs.pytorch.org/tutorials/recipes/recipes/profiler_recipe.html)
- `record_shapes` and `with_stack` add overhead — [PyTorch 2.15 profiler docs](https://docs.pytorch.org/docs/2.15/profiler.html)
- `tensorboard_trace_handler(dir_name, worker_name=None, use_gzip=False)`; `export_chrome_trace(path)` (with Kineto only the last cycle is exported) — [PyTorch 2.15 profiler docs](https://docs.pytorch.org/docs/2.15/profiler.html)
- `prof.key_averages(group_by_input_shape=True).table(sort_by="self_cuda_time_total", row_limit=N)`; sort keys include cpu_time_total, self_cpu_time_total, self_cuda_time_total — [Profiler recipe](https://docs.pytorch.org/tutorials/recipes/recipes/profiler_recipe.html)
- `export_memory_timeline` is deprecated in favor of `torch.cuda.memory._record_memory_history` / `_export_memory_snapshot` — [PyTorch 2.15 profiler docs](https://docs.pytorch.org/docs/2.15/profiler.html)
- `torch.autograd.profiler.emit_nvtx(enabled=True, record_shapes=False)` emits an NVTX range per autograd op, with `seq=` tags linking backward to forward ops; the documented example is for nvprof — [PyTorch autograd docs](https://docs.pytorch.org/docs/2.15/autograd.html)

### Inferences
- Reading the table: self CPU time >> self CUDA time with many small ops = launch/CPU-bound; large `cudaStreamSynchronize`/`aten::item`/`aten::to` = sync-bound; in the Perfetto trace, gaps in the GPU stream rows with the `ProfilerStep` CPU row busy = launch-bound, and with `DataLoader` / `enumerate(DataLoader)` slices long = data-bound. `torch.cuda.nvtx.range_push/pop` and `torch.cuda.profiler.start/stop` are the standard APIs for labels and capture ranges, but I did not get a page confirming them.
- Use `--pytorch=autograd-nvtx` in nsys instead of `emit_nvtx` (documented in the nsys guide).

### Gaps
- `torch.cuda.nvtx` and `torch.cuda.profiler.start/stop` docs were not retrieved.
- Profiler overhead numbers: not found.

## nvidia-smi, DCGM, and why "GPU-Util" misleads

### Takeaway
nvidia-smi GPU-Util only says a kernel was running during the sample window. DCGM's profiling fields (SM active, occupancy, tensor active, DRAM active) are better, but DCGM docs do not confirm they work on GeForce (RTX 4090).

### Cited Findings
- NVML defines the gpu utilization field as the share of the past sample period in which one or more kernels were executing, with windows of 1 s to 1/6 s depending on product; it is a liveness signal, not an efficiency signal — via [PyTorch torch.cuda.utilization docs](https://docs.pytorch.org/docs/2.13/generated/torch.cuda.utilization.html) and [arthurchiao.art](https://arthurchiao.art/blog/understanding-gpu-performance/) (secondary; NVML page itself 404'd)
- DCGM profiling fields: 1001 GR_ENGINE_ACTIVE, 1002 SM_ACTIVE, 1003 SM_OCCUPANCY, 1004 PIPE_TENSOR_ACTIVE, 1005 DRAM_ACTIVE, also FP64/FP32/FP16 pipes 1006-1008, PCIe/NVLink bytes 1009-1012 — [DCGM feature overview](https://docs.nvidia.com/datacenter/dcgm/latest/user-guide/feature-overview.html)
- SM Active counts a warp as active even when it waits on memory; 0.8+ is necessary but not sufficient, below 0.5 likely ineffective; occupancy alone isn't an efficiency signal; values are interval averages — [DCGM feature overview](https://docs.nvidia.com/datacenter/dcgm/latest/user-guide/feature-overview.html)
- `dcgmi dmon -e 1001,1004,1005` (default 1 Hz, `-d` ms, min 100 ms); high rates may return zeros; `dcgmi profile --pause/--resume` frees counters for Nsight; profiling needs admin privileges (nv-hostengine as superuser) — [DCGM feature overview](https://docs.nvidia.com/datacenter/dcgm/latest/user-guide/feature-overview.html)

### Inferences
- `nvidia-smi dmon`/`pmon` fields (sm%, mem%, pwr, clocks, PCIe) were not verified from docs by me. Use them as coarse hints only: low power draw plus high "util" suggests tiny kernels or stalls. On a Vast container DCGM likely is not installed and counters may be restricted; prefer nsys/torch.profiler.

### Gaps
- nvidia-smi dmon/pmon column definitions: not retrieved.
- DCGM profiling support on GeForce 4090: not confirmed.

## Hugging Face Trainer / TRL integration

### Takeaway
No built-in Trainer profiler callback was found; write a small TrainerCallback that owns a `torch.profiler.profile` and calls `prof.step()` once per optimizer step.

### Cited Findings
- Oumi ships a `ProfilerStepCallback(TrainerCallback)` wrapping a torch profiler; it calls `profiler.step()` in `on_step_begin` and tags microsteps with `record_function` — [Oumi source](https://oumi.ai/docs/en/latest/_modules/oumi/core/callbacks/profiler_step_callback.html)
- A Trainer "step" is one optimizer update (n batches with gradient accumulation); hook order includes on_substep_end, on_pre_optimizer_step, on_optimizer_step, on_step_end — [Transformers callbacks](https://huggingface.co/docs/transformers/trainer_callbacks), [v4.1.1 callbacks](https://huggingface.co/transformers/v4.1.1/main_classes/callback.html)

### Inferences
- Sketch (mine, untested): in `on_train_begin` create `profile(activities=[CPU,CUDA], schedule=schedule(skip_first=3, wait=1, warmup=2, active=3, repeat=1), on_trace_ready=tensorboard_trace_handler(dir), record_shapes=True)` and `__enter__()`; in `on_step_end` call `prof.step()`; in `on_train_end` call `__exit__`. Pass via `SFTTrainer(callbacks=[...])`. Keep `with_stack` off first. Because a step spans gradient accumulation, active=3 covers 3 full optimizer steps.

### Gaps
- No official HF doc example for profiling; TRL-specific docs not checked.

## vLLM profiling

### Takeaway
In current vLLM use `--profiler-config` (not the old `VLLM_TORCH_PROFILER_DIR`) with `/start_profile` and `/stop_profile`, or wrap the server in nsys with cudaProfilerApi capture. Official docs warn it is for developers and slows inference.

### Cited Findings
- `--profiler-config` requires vLLM >= 0.13.0; nsys is low overhead, torch profiler is medium overhead and not for benchmarking — [vLLM profiling docs](https://docs.vllm.ai/en/latest/contributing/profiling/)
- Torch: `--profiler-config '{"profiler": "torch", "torch_profiler_dir": "./vllm_profile"}'`; options torch_profiler_record_shapes, _with_memory, _with_stack (default on), _with_flops, _use_gzip (default on), _dump_cuda_time_total (default on), _activities; each worker writes `rank0.<id>.pt.trace.json.gz` plus a summary table; view at ui.perfetto.dev — [vLLM profiling docs](https://docs.vllm.ai/en/latest/contributing/profiling/)
- Server flow: POST `/start_profile` (optional JSON: profile_prefix, delay_iterations, max_iterations, profiler_kwargs), send a few requests, POST `/stop_profile`; a second start before stop returns HTTP 409; stop can take minutes; keep requests few — [vLLM profiling docs](https://docs.vllm.ai/en/latest/contributing/profiling/)
- Offline example: `examples/features/profiling/simple_profiling_offline.py`; `vllm bench serve --profile` also works — [vLLM profiling docs](https://docs.vllm.ai/en/latest/contributing/profiling/)
- nsys with vLLM: set `VLLM_WORKER_MULTIPROC_METHOD=spawn`; flags `--trace-fork-before-exec=true --cuda-graph-trace=node`; server mode adds `--capture-range=cudaProfilerApi --capture-range-end repeat` and `vllm serve --profiler-config.profiler cuda`; analyze with `nsys stats` — [vLLM profiling docs](https://docs.vllm.ai/en/latest/contributing/profiling/)
- Also: Triton Proton profiler (`"profiler": "proton"`), `VLLM_GC_DEBUG=1` for GC cost — [vLLM profiling docs](https://docs.vllm.ai/en/latest/contributing/profiling/)

### Inferences
- With CUDA graphs on, `--cuda-graph-trace=node` is needed to see kernels inside graphs. For a VLM, check whether time sits in the vision encoder / image preprocessing (CPU) vs decode kernels; low GPU kernel density between steps = CPU scheduling/launch-bound, often fixed by larger batch.

### Gaps
- Whether `VLLM_TORCH_PROFILER_DIR` still works in 0.32 (docs show only the new flag; I did not check a deprecation note). vLLM 0.32 specifically was not checked against these docs ("latest" page).
