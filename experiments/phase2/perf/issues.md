# Speed and memory issues (2026-11-06)

What makes our GPU work slow, why, and where in the code. Read this first. The fixes and their tests are in
`fix_test.md`.

Sources:

- three read-only code audits (training path, evaluation path, CUDA checker), checked against the installed
  libraries: transformers 5.20.1, trl 1.16.0, peft 0.21.2, torch 2.14.1;
- CPU measurements on the Mac;
- the GPU numbers of R0 and R1 (RTX 4090);
- the research report `docs/research/reports/NVIDIA CUDA performance for VLM training.md`.

Numbers marked *estimate* are not measured. The profiling run (`fix_test.md` F7) measures them.

Next run: RL (R2, the G-probe on R1's e2), on an **RTX 5090** (owner, 2026-11-06). There are no more SFT runs, so the
issues are ordered by their effect on RL.

---

## 0. Background: how the CPU and the GPU work together

```
   CPU (the manager)                          GPU (the factory)
 ┌────────────────────┐                     ┌──────────────────────────┐
 │ one Python thread  │  orders (kernels)   │ thousands of small cores │
 │ writes the orders  │ ══════════════════► │ do the orders            │
 │                    │ ◄══════════════════ │                          │
 │ CPU memory         │  PCIe road          │ GPU memory (24 / 32 GB)  │
 └────────────────────┘                     └──────────────────────────┘
```

- **Kernel**: one job that the CPU gives to the GPU, for example one matrix multiply. Each launch costs the CPU about
  5–20 µs.
- **Sync**: the CPU stops and waits until the GPU has finished all its orders. Examples are `.item()`, `.tolist()`,
  printing a GPU number, and `cudaDeviceSynchronize`. During a sync the CPU writes no new orders.
- **Launch-bound**: the GPU does each order faster than the CPU can write the next one. The GPU then waits between
  orders.

```
GPU-bound (good): the GPU is always busy
CPU : [o][o][o][o][o][o][o][o]
GPU :   [■■■■■■][■■■■■■][■■■■■■][■■■■■■]

Launch-bound (R1): small orders, gaps between them
CPU : [o][o][o][o][o][o][o][o][o][o][o][o]
GPU :   [■] [■] [■] [■] [■] [■] [■] [■]          ← R1: GPU busy ~60 %, one CPU core at 100 %

A sync in the middle: the queue empties, so the GPU waits
CPU : [o][o][SYNC ......wait......][o][o]
GPU :   [■][■][■■■■]               [■][■]
                    ▲ GPU idle while the CPU starts again
```

Two tools fight launch cost:

- **torch.compile**: it merges many small kernels into fewer, bigger kernels, so the CPU sends fewer orders.
- **CUDA Graphs**: the GPU records a whole list of orders once and replays the list with one launch. The shapes and
  memory addresses must stay the same.

### Are we using them today?

| Path | torch.compile | CUDA Graphs |
|---|---|---|
| SFT training (`train_sft`) | **no** (eager) | **no** |
| RL G-probe: update step | **no** | **no** |
| RL G-probe: answer generation (HF `generate`, `use_vllm: False`) | **no** | **no** (dynamic KV cache) |
| Exam sampling (`sample_vllm`, vLLM 0.32) | **yes** (vLLM default, no `enforce_eager`) | **yes** (decode steps) |
| HF check in `select_ckpt` (`hf_sample`) | **no** | **no** |
| CUDA checker (`objdel._gpucheck`) | n/a | **no** |

Only vLLM uses them. That is one reason vLLM gives one answer in about 0.04 s and HF `generate` in about 1 s. It is
also part of vLLM's 95 s start: compiling the model and capturing the graphs.

---

## 1. RL issues (R2, the next run)

One GRPO step in our config (`configs.grpo_config`: G = 8 answers per prompt, 32 prompts, 256 answers):

```
one optimizer step
┌──────────────────────────┬──────────────┬────────────────────────────────────────────┐
│ 1 GENERATE 256 answers   │ 2 SCORE      │ 3 UPDATE                                   │
│ HF generate, one batch   │ CUDA checker │ ref forward (LoRA off) + 32 micro-batches  │
│ of 256 sequences         │ 32 groups    │ of 8 answers (forward + backward each)     │
└──────────────────────────┴──────────────┴────────────────────────────────────────────┘
```

### RL-1. Generation runs out of memory: 256 sequences in one HF `generate` call (blocker)

Evidence:
- `g_probe.py:82` sets `use_vllm: False` (HF generation, plan D13).
- `grpo_config`: `per_device_train_batch_size 8` × `gradient_accumulation_steps 32` → TRL's `steps_per_generation`
  is 32, so the generation batch is **256 sequences** (`grpo_trainer.py:1253-1275`, checked: `generation_batch_size`
  256).
- `grpo_trainer.py:1873`: one `unwrapped_model.generate(...)` call for the whole batch.

The KV cache is the memory where the model keeps its notes for each word so far. Qwen3-VL-2B needs **112 KB for each
token** (PLAN D13). HF `generate` keeps the notes for every sequence until the **longest** one ends:

```
KV cache = 256 sequences × (prompt ~330 + longest answer) tokens × 112 KB

longest answer   400 tokens   → 256 × 730  × 112 KB ≈ 20 GiB
longest answer   700 tokens   → 256 × 1030 × 112 KB ≈ 28 GiB   ← > 24 GB (4090), ≈ 32 GB (5090)
runaway to 1536 tokens        → 256 × 1866 × 112 KB ≈ 51 GiB   ← out of memory on any of our cards
+ model weights ~4.3 GiB + activations
```

About 1–2 % of R1's answers do not stop. With 256 answers per step, almost every step has one, so **R2 will very
likely stop with an out-of-memory error** on a 4090. A 5090 helps only for short steps. PLAN §4.8 line 552 already
says "G-probe memory is not measured". The plan's own HF limit is 64 sequences (8 items × 8).

### RL-2. HF generation is slow (~1 s per answer)

R0 measured HF `generate` at **0.973 s per answer** (batch 64), against **0.04 s** for vLLM. The reasons:

```
HF generate, one new token for all sequences:
CPU : [~600 small launches ...][SYNC "all done?"][~600 launches ...][SYNC]
GPU :   [■][■][■][■][■] ...                        [■][■][■] ...
          28 layers × (norm, q, k, v, rope, attn, o, mlp, + LoRA A/B for each linear)

- no CUDA Graph, no compile          → the CPU writes every order again for every token
- host check every token             → `unfinished_sequences.max() == 0` (generation/utils.py:385-396)
- DynamicCache torch.cat every token → copies the whole cache again (cache_utils.py:135-136)
- LoRA not merged during generation  → extra A and B matmuls on every linear layer
```

*Estimate* for the G-probe:
- 20 steps × 256 answers + the probe set (3 draws × 128 items × 8 answers) ≈ 8,200 answers;
- × ~1 s ≈ **2.3 h** of generation (PLAN budget: 0.5 h);
- with vLLM colocate at ~0.04–0.09 s per answer: **6–12 min**.

### RL-3. Update step: 32 small micro-batches

- `grpo_config`: micro-batch 8 answers × 32 accumulation steps.
- Each micro-batch costs the CPU the same number of launches, whatever its size (the same problem as R1, with 4× more
  micro-batches).
- Gradient checkpointing is on (the GRPOConfig default; checked: `gradient_checkpointing True`), so the forward runs
  twice.

```
today      : [8][8][8][8][8][8][8][8][8][8][8][8][8][8][8][8] ... 32 times
better     : [  32  ][  32  ][  32  ][  32  ][  32  ][  32  ][  32  ][  32  ]   8 times
```

### RL-4. Extra syncs in every micro-batch

These are shared with SFT. The audit counts 60–80 syncs per micro-batch:

- the `isnan/isinf` loss check (`trainer.py:1891-1894`, `logging_nan_inf_filter` default True; checked for GRPO);
- the TRL metrics `.item()`;
- Qwen3-VL `get_rope_index` `.tolist()` / `.item()` (`modeling_qwen3_vl.py:984-1021`);
- the vision `.tolist()` per layer (`:297-299`);
- FLOP counting walks ~2.5k PEFT modules per micro-batch (`trainer.py:1902, 4206-4209`);
- blocking host-to-GPU copies (`trainer_pt_utils.py:1518`).

### RL-5. Scoring: many small checker calls from 4 processes

- `rewards.TRLReward(..., pool=fork_pool(4))`: "one checker call per group on a pool" (`rewards.py:89-92`). That is
  32 calls of 8 answers per step.
- Each worker process opens its own GPU context. 4 contexts share the GPU by time-slicing, and they also use GPU memory
  next to the trainer.
- Each call has the ~8 ms fixed cost of issue CK-1 to CK-3 below.

---

## 2. CUDA checker issues (`cuda/gpucheck.cpp`)

The checker validates polygons (13 defect counts) and paints them into masks. It is **bit-exact** against the CPU code
(polycore), but on the RTX 4090 it ran at **0.002–0.19×** the CPU speed. The kernels (the GPU math) are not the problem.
The problem is the memory work around them.

### CK-1. New buffers and ~40 waits for every call

`gpucheck.cpp:85-103, 809-815, 731-736`. Each call does this about 20 times:
`cudaMallocManaged` → `cudaMemset` → `cudaDeviceSynchronize` → host write. At the end it does `cudaFree` 20 times.
`cudaFree` also waits.

```
today, one call:
alloc+wait alloc+wait alloc+wait ... (×20) │ kernels │ free+wait free+wait ... (×20)
└──────────────── ≈ 8 ms fixed cost, even for 1 polygon ────────────────────────┘

fixed: buffers made once, kept, reused (grown only when a batch is larger)
copy in │ kernels │ copy out │ wait       ≈ 50 µs (estimate)
```

### CK-2. Data moves by page faults, not by copies

There is no `cudaMemcpy` in the file. The buffers are "managed memory": a page moves between the CPU and the GPU only
when a side touches it. One page fault is one small trip on the PCIe road.

```
page faults (today, ~0.6 GB/s measured from RESULTS.md):
CPU ─pg─► GPU ─pg─► CPU ─pg─► GPU ─pg─► ...   thousands of small trips, in both directions

one copy (fix, ~25 GB/s on PCIe 4.0, ~50 GB/s on PCIe 5.0 / RTX 5090):
CPU ══[█████████ all inputs █████████]══► GPU
CPU ◄═[█████████ all masks  █████████]═══ GPU
```

Evidence from RESULTS.md: the cost per extra shape grows with the **bytes of output** (8 KB → 9.8 µs, 131 KB →
235 µs, 524 KB → 835 µs). It does not grow with the edge work: 30× more edges cost only 3.6× more.

### CK-3. Two full rasters go back, and one is not used

- `gpucheck.cpp:813, 823-824`; `gpucheck_py.cpp:33-36`; `shape.py:107`.
- The checker returns an even-odd mask and a nonzero mask, at 1 byte per pixel each.
- `to_masks` uses only even-odd. A mask needs 1 bit per pixel.
- The host also copies the mask again (`assign`), and the `cudaMemset` of the masks is useless because every pixel is
  written.

```
today :  even-odd [1 byte/px] + nonzero [1 byte/px]     = 16 bits per pixel
fixed :  even-odd [1 bit/px]                             =  1 bit per pixel  (16× less traffic)
```

### CK-4. Slow store pattern (second order, after CK-1 to CK-3)

- `gpucheck.cpp:221-236`: one thread paints one row. Neighbouring threads write bytes `W` apart, so the writes are not
  coalesced (32 separate memory sectors per warp store).
- The `buf[256]` array lives in slow local memory.
- Other second-order costs: emulated 64-bit division (`:147-151`) and O(E²) edge pairs with re-reads (`:373-388`).

### CK-5. Calls are too small

- Scoring calls the checker with one group of 8 answers (`score.py:159,173`; `metrics.py:246-253`).
- Fixed cost of 8 ms against <1 ms of real work: the GPU cannot win at 8.
- After CK-1 to CK-3 the fixed cost falls to about 50 µs. The GPU should then beat 8 CPU threads at batches of about
  32–64 polygons or more (*estimate*).

Note: RESULTS.md compares against `score_batch`, which returns areas only. The fair baseline for `to_masks` is the
polycore loop. Against the loop, the GPU is already 1.06× at 16 vertices and batch 1024.

---

## 3. Evaluation issues (`select_ckpt`, G-probe probe sets)

| # | Issue | Evidence | Cost *estimate* |
|---|---|---|---|
| EV-1 | One vLLM engine start (~95 s) for each sampler call, 8–9 per select run | `select_ckpt.py:136-163, 195-207` | ~14 min per select run |
| EV-2 | HF check: 512 answers at ~1 s each (issue RL-2) | `select_ckpt.py:284-302`; `hf_sample.py:62-79` | 8–9 min |
| EV-3 | Per-request seeds may make vLLM sample row by row (not verified) | `sample_vllm.py:44,52` | up to 1.5–3× decode |
| EV-4 | The GPU waits while scoring and the CPU fp32 vision check run | `select_ckpt.py:335-344, 379-383` | 2–8 min |
| EV-5 | Pool fixed at 4; sampler timeout ~2.5 h when an engine hangs | `select_ckpt.py:361, 58, 206` | cost only on failure |
| EV-6 | Resolved by F1 (2026-11-06): G-probe's probe sets go through the trainer's own `_generate`, which is vLLM colocate (`configs.grpo_config` `use_vllm` for every model but SmolVLM); P1's fast arm logged vLLM's IS metrics, and its step-0 sampling (2 × 8 items × 8) took 14.7 s | `g_probe.py:41-48`; `sample.py:49-71` | none |

```
select run today (one checkpoint):
GPU : [start 95 s][sample] [start 95 s][sample] [start 95 s][sample] ...
CPU :                     [score]                [score]          [vision]
better:
GPU : [start 95 s][sample all views of this checkpoint][next checkpoint ...]
CPU :                            [score ............... ][vision]
```

---

## 4. SFT issues (kept for reference; no SFT run is planned)

| # | Issue | Evidence | Gain *estimate* |
|---|---|---|---|
| SFT-1 | Micro-batch 4 × 8. The "152k-vocab logits ~12 GB" reason is out of date: TRL 1.16 scores only label positions, in chunks | `r1_qwen3vl2b.yaml:11-12`; `trl/trainer/utils.py:1548-1660` | −30 to −45 % |
| SFT-2 | Gradient checkpointing on, with 6.1 GiB of 24 GiB used (R1 manifest) | `configs.py:72,164` | −15 to −25 % |
| SFT-3 | Syncs (issue RL-4) | — | −5 to −15 % |
| SFT-4 | No torch.compile | — | −15 to −30 % more |

Measured on the Mac:
- the data path costs 19 ms per micro-batch, in 6 workers: **not the limit**;
- padding is 10 %;
- the loss mask covers the answer and `<|im_end|>`, which is correct.

---

## 5. The RTX 5090 (owner, 2026-11-06)

| | RTX 4090 (Ada, sm_89) | RTX 5090 (Blackwell, sm_120) |
|---|---|---|
| GPU memory | 24 GB | **32 GB** |
| bf16 dense peak (`image_check.sh` table) | 165.2 TFLOPS | 209.5 TFLOPS |
| Memory bandwidth | ~1.0 TB/s | ~1.8 TB/s |
| PCIe | 4.0 x16 | 5.0 x16 |
| Policy price ceiling (`policy.CARDS`) | $1.00/h (paid ~$0.40–0.44) | $1.50/h (quote not checked yet) |

What it changes:

- **Speed:** a faster GPU does **not** fix launch-bound work. The CPU still writes the same orders. Unless the fixes
  below land first, the 5090 waits even more (busy share falls).
- **Memory:** 8 GB more is good for RL-1. It is not enough for 256 sequences with a runaway answer (51 GiB).
- **Code support:**
  - the tooling already supports the card (`offers.sh RTX_5090`, `policy.CARDS`, the 5090 peak in `image_check.sh`);
  - the CUDA checker builds `native` on the box, and `cuda/docker_check.sh` builds sm_120;
  - the image's torch 2.14.1 is a cu130 build (Blackwell needs CUDA ≥ 12.8).
- **Gaps:**
  - `preflight.DPH_CEILING = 1.00` is the 4090 ceiling;
  - every timing in `project_r1` / preflight comes from the 4090;
  - vLLM 0.32 and the FlashInfer-off setting are not yet run on sm_120;
  - the CPU on Vast 5090 hosts is unknown, and a slow CPU makes launch-bound work worse.
