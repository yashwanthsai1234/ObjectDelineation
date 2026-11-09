# SD1: GPU-first throughput architecture, Phases A-E (system-design note)

Date: 2026-11-08. Author: fresh system designer (no prior opinion), lens = **GPU-first throughput**.
Scope: data -> SFT/RL loop on one RTX 5090 -> exact CUDA checker / scorer / polygonizer -> evaluation -> artefacts,
plus the growth path to 3D voxels and 2+ GPUs. Design note only: no repo edits, no git, no spend.

Tags: [measured: file] recomputed from a repo log or data file in this session; [code: file:line] read in code;
[paper/web: source]; [SP: note] taken from a scratch research note (not re-derived unless said); [estimate] my
arithmetic or judgment; [U] unverified, must be measured. REPO = the project root, SP = the scratch folder.

---

## 0. Skills used

| Skill | What it changed in this note |
|---|---|
| system-design-skills:system-design | The method: clarify -> estimate -> design -> trade-offs -> failure -> growth, the solves/worsens/when-to-change rule, the quality bar and Quick Diagnostic (§12). |
| system-design-skills:design | Ran the loop inline (I am already the orchestrator, so no sub-agent). It forced the coverage sweep (§11). |
| system-design-skills:requirements-scoping | §1: three lists (functional / non-functional with numbers / out of scope) and written assumptions in place of questions I cannot ask the owner. |
| system-design-skills:back-of-the-envelope | §2: every choice rests on a number: bytes per step, GPU-seconds per phase, memory on 32 GB, verifier load. It killed two ideas: a GPU-resident frame store (5.8 GB does not fit beside vLLM awake) and a cross-step vision-feature cache (~88 GB). |
| system-design-skills:data-storage | §5 store design: keyed by the access path (row id -> frame index, (item, object, view) -> doctor-table handle), preload, no growth in the step. |
| system-design-skills:caching | §5/§7: which caches pay (prompt ids + rope positions: 1 distinct value per run; vision features: per step only), versioned keys (reward version, data SHA), "no cache without a hit-rate number". |
| system-design-skills:resilience-failure | §8: one GPU = one process = one fault domain; containment by sanitizer gate + host-side input checks + checkpoint/resume; stale-beats-error does not apply to rewards (fail closed). |
| system-design-skills:observability | §9: NVTX phase ranges as the RED/USE frame of a training step, step-log alarms, "watch the watcher" (timing syncs are themselves syncs). |
| system-design-skills:architecture-diagram | Its rules applied to the ASCII diagrams (typed boxes, numbered arrows, dashed = async/off-path). I did not render the HTML file: the task asks for ASCII in a markdown note. |
| nvidia-skills:nvidia-skill-finder | Live catalog check (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`, fetched today). Relevant, recommend only, not installed: `compute-sanitizer-diagnose-cuda-bugs` (check 8), `nemo-mbridge-perf-cuda-graphs` (per-layer graphs, stage 1c), `nemo-mbridge-perf-sequence-packing` (token-budget packing), `dali-dynamic-mode` (GPU image pipeline alternative), `launch-nemo-rl` (stage 3), `tilegym-*` (cuTile kernels, not needed). No Nsight skill exists. |
| engineering:system-design | Requirements / high-level / deep dive / scale / trade-offs framing; "what to revisit as it grows" -> §10. |
| engineering:architecture | ADR shape for the two load-bearing decisions (§7.1 single process colocate; §7.2 in-place merged-weight push). |
| engineering:code-review | Correctness lens on each proposal: bit-exactness (fast-math, FMA contraction, float atomics), hidden syncs, memory spikes. It produced checks 2, 3, 5. |
| engineering:debug | "40 % idle" treated as a symptom: separate decode tail vs CPU phase vs handover with the R5 numbers before choosing fixes (§2.4). |
| engineering:documentation | Answer first (§0.1), then evidence; links to SP notes instead of copies. |
| huggingface-skills:huggingface-best | Used the main thread's live leaderboard lookup (docs/implementation2.0/notes/MAIN_skill_lookups.md): Qwen3-VL-4B RefCOCO 88.2 vs 84.8 for Qwen3-VL-2B and Qwen3.5-2B; Qwen3.5 is hybrid linear attention, which breaks the shared-prefix 4-D mask [SP: G_gpu_techniques §8]. Result: keep Qwen3-VL-2B through Phase E; Qwen3-VL-4B is the only drop-in scale point (same layers, same kernels). |

Web research done today (late 2026): vLLM weight-transfer API (IPC backend for colocated training)
[web: docs.vllm.ai/en/stable/training/weight_transfer/]; RTX 5090 = 32 GB GDDR7, 1,792 GB/s, ~209.5 TFLOPS dense
bf16 tensor (sources disagree on accumulate mode) [web: techpowerup c4216, tomshardware review]; nvImageCodec decodes
PNG only through a CPU fallback (GPU codecs are nvJPEG / nvJPEG2000 / nvTIFF) [web: github.com/NVIDIA/nvImageCodec];
vLLM accepts precomputed image embeddings for Qwen3-VL only with DeepStack features included, one known bug (#31708)
[web: discuss.vllm.ai/t/2294, /2203]; verl >= 0.8.0 ships fully async GRPO for Qwen3-VL (docs/implementation2.0/notes/D_rl_infra.md reports
v0.9.1) [web: newreleases.io verl v0.8.0]; FlexAttention on sm_120: no evidence found either way [U].

### 0.1 The answer in one screen

1. **The geometry is not the bottleneck; the plumbing around the model is.** The checker waits 5 ms per 23.4 s
   step [measured: runs/r5/shared_rl.json `time/rewards_s` 0.005]. The idle GPU time (~9-10 s) is the CPU phase
   (3.0-3.2 s), the weight handover (1.66 s), and the decode tail (4-5 s) [measured + SP: B_gpu_audit §1].
2. **Design rule: only three things cross PCIe per step**: 32 grey frames in (8.4 MB), answer token ids + sampling
   log-probs (<= 2.5 MB), numbers out (~25 KB). Today 1.61 GB of fp32 pixels is built twice on the CPU and copied
   once, blocking [code: fast.py:473; SP: B_gpu_audit §2 row 2].
3. **One process, one GPU, vLLM awake, in-place merged-weight push, token-budget micro-batches, GPU-resident doctor
   store, GPU pixel build from a uint8 frame store, CUDA polygonizer v2 everywhere offline.** Target 14 s/step with
   KL 0 and the cap at the measured target length; 16 s with KL kept [estimate, §2.4]. Stretch: ~11 s with
   group-streamed update during the decode tail (stage 1b) [estimate].
4. **Order of spend:** run the exact (bit-identical) speed changes and the NVTX capture in the same short GPU
   session as preflight, then run the data-first RL run (B1) once, on the faster code: B1 costs ~1.2 GPU-h instead of
   ~2.0 [estimate, §2.5]. The whole A-D GPU programme fits in about $1.0-1.4 of the $3.
5. **Growth:** 3D whole-nodule answers are ~6x longer (mean 5.8 slices per nodule [measured]) and the prompts carry
   ~6x more image tokens, so 3D is the point where one GPU stops being enough: disaggregate generation (GPU 0) and
   training + verifier (GPU 1), LoRA-only deltas over PCIe (35 MB), verl async.

---

## 1. Requirements (requirements-scoping)

**Restated in one sentence.** A 2B VLM writes a nodule outline as JSON polygon text for a 512x512 CT image; an exact
CUDA engine checks and scores it against the set of doctors; the loop SFT -> GRPO -> evaluation runs on one RTX 5090
and keeps the GPU busy, with the CPU only where no GPU path exists.

### 1.1 Functional (core four, ranked)

| # | Function | Notes |
|---|---|---|
| F1 | Build train / val / test data with certified targets and centred + shifted + rescaled train views, with quality gates | Owner rule: dataset quality first. 24,624 RL rows (8,208 x 3) [SP: BRIEF §2]. |
| F2 | Train: SFT (exists, reused), then GRPO with exact rewards | TRL 1.17.1 + FastGRPO + colocated vLLM [code: configs.py:183-193]. |
| F3 | Verify: legality (13 defect classes), IoU against the doctor set, boundary F, area, Betti topology; polygonize masks | One source, CPU and nvcc [code: cuda/gpucheck.h:1-6]. Polygonizer must be CUDA in production (owner). |
| F4 | Evaluate a checkpoint without a trainer, save all 8 answers, label-free selection (medoid), fair ceiling, nodule-level CIs, GED / HM-IoU / HD95 on the GPU | Direction A1-A5, D5. |
| F5 | Publish artefacts: adapter, answers, metrics, run record, profile | HF private repo (exists). |

### 1.2 Non-functional (numbers)

| Quality | Target | Today |
|---|---|---|
| Exact | CPU == GPU bit-for-bit for every geometric output; target hashes byte-identical | Yes for checker/scorer (195,809 golden cases) [SP: B_gpu_audit §3.2] |
| GPU busy | Fully idle < 1.5 s per step; prepare phase < 0.5 s | ~5.2-5.9 s fully idle + 4-5 s tail [SP: B_gpu_audit §1] |
| Step time | <= 16 s (KL kept), <= 14 s (KL 0 + cap), same step-1 reward | 23.43 s [measured: runs/r5/shared_rl.json] |
| Transfers | H2D < 50 MB/step, no pixel/mask D2H, host syncs < 200/step | 1.61 GB H2D, ~1,250 syncs [SP: B_gpu_audit §2] |
| Memory | Peak <= 30 GB of 32 GB under worst-case lengths | 17 GB peak with sleep mode [measured note: docs/research/notes/full_rl_2026-11-08.md:19] |
| Quality | Test IoU from 0.645 toward 0.70-0.73; far misses < 10 % | 0.645, F 0.816, legal 0.994 [measured: runs/r5/eval/eval_test/final_eval.json] |
| Money | Phases A-D GPU work <= $3; NA hosts only; preflight + owner OK per rental | rule in code [SP: BRIEF §3] |
| Reproducible | Every artefact names code SHA, checker build, reward version, data SHA-256 | partly |

### 1.3 Out of scope

CAD (parked by owner); multi-region / HA serving; a public inference API; Qwen3.5 swap (breaks shared-prefix
kernels [SP: G_gpu_techniques §8]); FP8 rollouts (raise the trainer-vLLM mismatch [SP: G §0 row 10]); grammar-
constrained decoding (masked sampler vs full-vocab trainer log-probs); any 2-GPU work before Phase C is measured.

### 1.4 Assumptions in place of owner questions (load-bearing; revisit if wrong)

| # | Question I would ask | Assumption | If wrong |
|---|---|---|---|
| Q1 | Synchronous on-policy, or 1-step-stale async allowed? | Synchronous on one GPU (exact on-policy); async only at 2+ GPUs | Allowing staleness on 1 GPU opens MPS async (§7.1 option B) |
| Q2 | KL weight 0? | Owner decision pending; design supports both (ref pass is a removable stage) | Keep beta 1e-3: +1.0-1.7 s/step |
| Q3 | Length cap from the true target maximum? | Cap = measured P100 of train target tokens + margin, set by a counted number, not by the unverified 1,182 note | Keep 1,536: tail stays 4-5 s |
| Q4 | Must vLLM see bit-identical pixels to the trainer? | Yes: one GPU pixel build feeds both (or n=8 requests so vLLM encodes 32 images) | - |
| Q5 | Host RAM on rented boxes >= 32 GB? | Yes (typical for 5090 offers) [U]; frame store needs 5.8 GB pinned | Fall back to memory-mapped shards, not pinned |

---

## 2. Estimates (back-of-the-envelope)

### 2.1 Fixed facts used

| Fact | Value | Source |
|---|---|---|
| Model text stack | 28 layers, hidden 2,048, 16 q heads / 8 kv heads, head_dim 128, FFN 6,144, vocab 151,936 | [code: runs/r3/ckpt/config.json] |
| Vision | depth 24, patch 16, merge 2, temporal patch 2, DeepStack layers 5/11/17 | same |
| Image tokens | 512x512 -> 32x32 patches = 1,024 patch rows of 1,536 values -> 256 merged tokens | arithmetic on config |
| Prompt length | ~325 tokens (256 image + ~69 text) | [SP: B_gpu_audit §2] |
| Answer length | mean 131 tokens; longest in a step mean 1,271; 38 % of R5 steps hit 1,536 | [SP: BRIEF §2] |
| KV cache per token | 28 x 2 x 8 x 128 x 2 B = 114,688 B ~ 112 KB | arithmetic |
| Parameters | 2.13 B -> 4.26 GB bf16; LoRA 17.4 M | [SP: B_gpu_audit §2 row 3, §4.6] |
| Frames | 512x512x3 uint8 PNG, R=G=B in all 5 files checked | [measured: data/views/lidc-train-shifted/images, 5 files] |
| Distinct train frames | 7,341 per view; 3 views | [measured: file counts in data/views/lidc-train-{shifted}] |
| Doctor tables | 14,682 shifted+rescaled caches = 1,022.5 MB as (k raters + uni + core + vote) u8 + dist f32; mean k 3.2; largest frame 200x200 | [measured: data/views/*/caches, this session] |
| Slices per nodule | 1,447 nodules, mean 5.8, median 4, p90 11, max 52 | [measured: data/prepared/images file layout] |
| RTX 5090 | 32 GB, 1,792 GB/s, ~209.5 TFLOPS dense bf16 | [web: techpowerup c4216] |
| Rental | $0.44-0.61/h | [SP: BRIEF; full_rl note: $0.610/h] |

### 2.2 Bytes per step (32 prompts x 8 answers)

| Path | Today | New design | Why |
|---|---|---|---|
| H2D pixels | 1.61 GB fp32 (1,024 x 1,536 x 4 B x 256 rows), pageable, blocking; built twice on CPU | 8.4 MB (32 x 512 x 512 x 1 B grey), pinned, async on a copy stream, prefetched one step ahead | GPU expands grey->3 channels, rescales, normalizes, patchifies |
| H2D ids / logps for the trainer | ~5.6 MB padded [SP: B §2 row 17] | <= 2.5 MB (256 x <= 1,200 x 8 B), mean 0.27 MB | vLLM returns these to the host anyway |
| H2D polygons for the checker | ~50 KB | ~50 KB | parse on CPU (§6.3) |
| H2D rope positions | 9.8 MB per heavy call | 0 | constant per run, computed once on GPU |
| D2H | rope 6.6 MB + many `.item()` | ~25 KB reward terms + ~2 KB metrics, once per step | numbers only |
| D2D weight push | 4.26 GB through ~700 `load_weights` calls + `empty_cache` + wake/sleep | 196 in-place slice writes, ~15 GB HBM traffic ~10 ms [SP: B §4.2] | no allocation, CUDA-graph addresses stay valid |
| D2D pixel tensor | 1.61 GB resident all step + 3.2 GB of compares | 0.2 GB (32 images fp32) or 0.1 GB bf16 | key rows by image index, not by pixel compare |
| Host syncs | ~1,250 | < 200 | §6.5 |

### 2.3 Memory on 32 GB (vLLM awake, single process)

| Item | GB | Basis |
|---|---|---|
| vLLM weights (bf16) | 4.3 | 2.13 B x 2 B |
| vLLM KV pool + graphs (util 0.35 = 11.2 GB total) | ~6.9 | 11.2 - 4.3; ~6 GB KV = ~53 k tokens |
| KV actually needed | 1.2 (32 prefixes x 325 tok) + ~3.8 (256 x 131 tok mean) | 112 KB/token; runaways are few |
| Trainer frozen base (bf16) | 4.3 | second copy (vLLM fuses qkv / gate_up, so sharing is not practical) |
| LoRA params + grads + AdamW (fp32) | 0.3 | 17.4 M x 16 B |
| Doctor store (full tables, all 3 views) | ~1.5 | 1.02 GB measured for 2 views, x1.5 [estimate]; 0.12 GB if bit planes + GPU EDT (§5.2) |
| Pixels + vision features for 32 images | 0.35 | 0.2 GB fp32 patches + 32 x 256 x 2,048 x 2 B x 4 (main + 3 DeepStack) |
| Activations, 8 k-token micro-batch, checkpointed | ~2-3 | ~115 KB/token boundaries + one layer recompute [estimate] |
| Allocator cache / fragmentation margin | 3-4 | [estimate] |
| **Total** | **~23-25** | headroom ~7 GB -> larger token budget or fewer checkpointed layers [U: measure] |

Two ideas killed by these numbers: a GPU-resident frame store (5.8 GB for 22,023 grey frames) would eat the headroom;
a cross-step vision-feature cache (22,023 x 4 MB with DeepStack ~ 88 GB) fits nowhere and has a low hit rate
(300 steps = 0.4 epoch, each frame seen < once).

### 2.4 Step time budget (seconds)

| Phase | Today [measured R5, steps 3-15] | New [estimate] | Change |
|---|---|---|---|
| Prompt prep (chat template, processor pass 1, 1.61 GB cat) | 0.3-0.8 | 0.02 | prompt ids memo (1 distinct prompt), GPU pixel build |
| Weight sync + wake + sleep | 1.27 + 0.11 + 0.29 | 0.05 | in-place merged push, vLLM awake |
| Decode (fixed ~3.4 + 4.3 ms x longest answer) | ~8.9 | ~7.5 | n=8 per request (-0.2-1.0), cap from measured targets (-0.6), `detokenize=False` (0-0.5) [U] |
| CPU phase (processor pass 2, 1.61 GB H2D, CPU rope, metrics) | 3.0-3.2 | 0.2 | GPU pixels, constant rope, device-side metrics |
| Reference pass | 1.4-1.7 | 0 (KL 0) or ~1.0 (64-row chunks) | owner decision |
| Update (fixed 4.5-5.0 + FLOP part ~2.8) | 7.67 | 4.0-5.0 | token-budget micro-batches (8 not 16: -2.2), sync cleanup (-0.3-0.8), partial checkpointing (-1) |
| Other (loader, optimizer, log) | 0.22 | 0.05 | frame gather is prefetched |
| **Total** | **23.43** | **~12-14 (KL 0) / ~13-15 (KL kept)** | |
| Stage 1b: group-streamed update in the tail | - | ~-2-3 | §7.3, only if NVTX shows tail > 2 s after the cap |

Debug check: the R5 fit `generate = 5.61 s + 4.3 ms x max_len` (r = 0.92) [SP: B §1] includes the 1.66 s handover and
~0.5 s prep, so pure decode is ~3.4 s + 4.3 ms per decode step. 4.3 ms per tail step vs a 2.4 ms weight read at
1,792 GB/s means the tail is half per-step overhead (scheduler, sampler, sync) [estimate].

### 2.5 GPU-seconds and dollars per phase of the programme

| Phase / job | GPU time | $ at $0.44-0.61/h | Basis |
|---|---|---|---|
| Preflight on box: build nvcc, golden files, sanitizer, memory worst case | ~15 min | 0.11-0.15 | [estimate] |
| A6 NVTX + C1-C5 exact changes, 2 x 15 steps (old vs new) + nsys capture of 2 steps | ~20 min | 0.15-0.20 | 15 steps x 23 s + 15 x 14 s + load |
| A3 re-run test eval saving answers (2,200 items x 8) | 10.5 min | 0.08-0.11 | [measured: final_eval.json `seconds.total` 633.6] |
| B1 RL shifted, 300 steps | 70 min at 14 s (117 min at 23.4 s) | 0.51-0.71 (0.86-1.19) | 300 x step time |
| B1 evals val + test | ~16 min | 0.12-0.16 | test 10.5 min + val pro-rata 1,084/2,200 |
| D6 GPU CI job (polygonizer v2 + EDT + sanitizer) | ~15 min | 0.11-0.15 | [estimate] |
| **Total A-D with C first** | **~2.5 h** | **~$1.1-1.5** | leaves ~$1.5 for B2/B3 or a re-run |

SFT is not repeated (R1 SFT exists: 3.40 s/step projected, 3,065 steps, billed $2.33 [SP: A_repo_history:199]).

### 2.6 Verifier load

256 answers/step, ~1 KB text each, 5 ms of GPU time [measured]. At 14 s/step: 18 answers/s and 0.04 % of the step.
Evaluation: 17,600 answers in ~630 s, i.e. 28 answers/s, sampling-bound [measured]. A batch of 256 answers is
~25 kernel launches; the engine is launch-bound, not throughput-bound. Conclusion: no service, queue, or HTTP is
justified on one GPU; an in-process function is the right interface (§7.4).

---

## 3. High-level architecture

### 3.1 Component diagram (one RTX 5090, one process; `==>` PCIe, `-->` on device, `..>` async/off path)

```
  OFFLINE (once per data version)                         ONLINE (each RL step)
 +--------------------------------------+   +-------------------------------------------------------------------+
 | [CPU] raw DICOM / RIGA readers (I/O) |   |                       RTX 5090 (32 GB)                            |
 |        |                             |   |                                                                   |
 |        v                             |   |  (2) GPU pixel build      (3) vision encode x32     policy cache   |
 | [CPU] PNG decode ONCE -> grey u8     |   |  u8 grey -> rescale,  -->  one batched call   -->  prompt ids +   |
 |        |                             |   |  normalize, patchify        features + DeepStack   rope (const)   |
 |        v                             |   |        ^                          |                     |         |
 | [GPU] view builder: shift / rescale  |   |        |                          v                     v         |
 |       frames + doctor masks          |   |        |              (4) vLLM in-process, awake, CUDA graphs     |
 |        |                             |   |        |                  32 requests x n=8, prefix cache        |
 |        v                             |   |        |                          |                               |
 | [GPU] CUDA polygonizer v2 -> targets |   |        |               ids+logps ==> [CPU] parse text -> ints      |
 | [GPU] checker: certify targets       |   |        |                          |        ==> polygons (50 KB)   |
 | [GPU] exact EDT -> doctor tables     |   |        |                          v              |                |
 | [GPU] gates: coverage, leakage,      |   |        |   (5) trainer: packed shared-prefix     v                |
 |       legality, length count         |   |        |       forward, token-budget  <-- (6) verifier (low-prio  |
 |        |                             |   |        |       micro-batches, loss,       stream): raster, legal,  |
 |        v                             |   |        |       backward, AdamW            IoU vs doctor set, F,    |
 | frame store: u8 shards (5.8 GB)  =====(1)=======+   |          |                       topology -> 256 x terms |
 | doctor store: tables (~1.5 GB)   =====(0) preload  |          v                              ^                 |
 | manifest + SHA-256               |   |        (7) in-place merged LoRA push --> vLLM params  |                 |
 +--------------------------------------+   |                                     doctor store (resident) -+     |
                                            +-------------------------------------------------------------------+
                                               ..> NVTX ranges, step metrics (~2 KB) ..> run record, checkpoints
                                               ..> every 50 steps: adapter (70 MB) -> disk -> HF private repo
```

### 3.2 Sequence of one RL step (numbered = order; `||` = concurrent)

```
 Loader     CopyStream   PixelBuild  VisionEnc   vLLM        CPU parse   Verifier     Trainer      WeightPush
   |            |            |           |          |            |           |            |              |
 1 |-idx[32]--->|            |           |          |            |           |            |              |
 2 |            |=8.4 MB H2D>|           |          |            |           |            |              |
   |  (prefetched during step k-1 update)           |            |           |            |              |
 3 |            |            |-patches-->|          |            |           |            |              |
 4 |            |            |           |-feats--->|  (or 32 images, n=8)   |            |              |
 5 |            |            |           |-feats-------------------------------------->|  kept for     |
   |            |            |           |          |            |           |            |  all passes   |
 6 |            |            |           |          |-decode 256 answers (graphs)       |              |
 7 |            |            |           |          |==ids,logps>|           |            |              |
 8 |            |            |           |          |            |=polys====>|            |              |
 9 |            |            |           |          |==ids,logps (pinned)================>|              |
10 |            |            |           |          |            |           |-256 terms->| (device)     |
11 |            |            |           |          |            |           |            |-advantage    |
12 |            |            |           |          |            |           |            |-fwd/bwd x8   |
13 |            |            |           |          |            |           |            |-AdamW step   |
14 |            |            |           |          |            |           |            |------------->|
15 |            |            |           |          |<-in-place write of 196 merged slices--------------|
16 |  || next step's frames already copied in step 12 (copy stream)                         | log ~2 KB D2H
```

Host waits per step in this design: one after decode (ids to host, inherent to vLLM), one for the verifier result,
one per micro-batch at most (loss is kept on device; read once in `log`), one in `log`. Target < 200 including vLLM's
per-decode-step sync [U for vLLM internals].

### 3.3 Offline data pipeline (data first)

```
raw LIDC / RIGA --[CPU, I/O: pydicom / PNG decode, isolated mode]--> grey u8 frames + doctor masks
      |
      v
group-aware split by patient (709 / 100 / 201 patients) ---- gate G1: zero shared patients across splits, views incl.
      |
      v
[GPU] views: centred, shifted, rescaled (frame + every doctor mask moved by the same integer transform)
      |
      +--> [GPU] polygonizer v2 -> targets -> [GPU] checker certifies -- gate G2: 100 % legal, hashes stable
      +--> [GPU] exact EDT + union / core / vote -> doctor tables     -- gate G3: bit-equal to scipy oracle
      +--> [GPU] coverage audit: offset / scale / size histograms, train vs test  -- gate G4: train covers test
      +--> [CPU] token count of every target                          -- gate G5: P100 length -> sets the cap
      |
      v
frame store shards (u8, safetensors) + doctor store (one file, preloaded) + manifest (SHA-256 per shard)
```

---

## 4. Component table

| # | Component | Layer (Architecture.md) | Runs on | Input -> output | Concrete interface | State |
|---|---|---|---|---|---|---|
| C1 | Readers (LIDC, RIGA) | data | CPU (file I/O, DICOM / PNG decode) | raw files -> grey u8 + masks | `read(item_id) -> (u8[H,W], masks u8[k,H,W])` | none |
| C2 | Frame store | data | pinned host -> GPU copy stream | row id -> u8 frame | `gather(idx: i64[32]) -> cuda u8[32,512,512]` (async, event) | 5.8 GB shards, read-only |
| C3 | View builder | data | GPU | frame + masks + (dx, dy, s) -> view frame + masks | `make_views(batch) -> batch` ; bit-equal to the current CPU builder or it stays CPU (gate) | offline |
| C4 | Polygonizer v2 | geom | GPU (CUDA) | bit-packed masks -> rings, quads | `polygonize_many(masks) -> (ring_start, pts_start, ring_area, pts, quads)` on device, 1 host wait | none |
| C5 | Checker / scorer | geom | GPU (CUDA) | polygons + table handles -> defects + terms | `score_polygons(Polygons, handles) -> Scores` [code: gpucheck.h] | pool, grow-free after start |
| C6 | Doctor store | verify | GPU resident | (item, object, view) -> handle | `preload(manifest) -> {key: handle}`; `handle(key) -> int` | ~1.5 GB (or 0.12 GB + EDT) |
| C7 | Verifier (reward vN) | verify | GPU + CPU parse | texts + keys -> rewards + terms + versions | `score_batch(spec, items[{id,key,text,truncated}]) -> results` (Architecture.md §4.2 schema, in process) | reward spec, versions |
| C8 | Policy cache | policy | GPU | image idx -> features; prompt -> ids, rope | `prompt_ids(text, HxW)`, `rope(len)`, `features(idx)` per step | 1 prompt per run |
| C9 | Sampler | policy | GPU (vLLM) | 32 prompts -> 256 answers | `generate(prompts, n=8, max_tokens, seed) -> ids, logps` | KV pool, CUDA graphs |
| C10 | Weight push | rl | GPU | LoRA A, B -> vLLM slices | `push(lora)`: `slice.copy_(W0); slice.add_((B@A*s).bf16())` | slice map built once |
| C11 | Trainer | rl | GPU | ids, logps, rewards, features -> LoRA update | TRL adapter + pure `forward(packed)` ; token budget per micro-batch | LoRA + AdamW |
| C12 | Evaluator | verify + policy | GPU | checkpoint + split -> metrics | `evaluate(ckpt, view, split, K=8) -> metrics.json, answers.jsonl` (no trainer) | none |
| C13 | Stats | verify | CPU (numbers only) | per-item scores -> CIs | nodule-level bootstrap | none |
| C14 | Run record + artefact push | infra | CPU | files -> HF repo | `record(run) -> json` with code SHA, checker hash, reward spec, data SHA | append-only |
| C15 | Preflight + cost guard | infra | CPU + short GPU test | config -> go / no-go + $ | `make preflight CONFIG=...`; NA hosts only | ledger |

IDs (sequencer sweep): answer id = `run/step/group/k` (deterministic, no service); item key =
`dataset/split/item::view` (exists: `lidc/train/0042::shift`, Architecture.md §4.2).

---

## 5. Data and stores (data-storage, caching)

### 5.1 Frame store

- **Key** = row index (int) -> frame index (many rows share a frame: disagreement rows x4 [SP: B §5 B4]).
- **Layout**: grey u8 [N, 512, 512] in safetensors shards of ~1 GB, memory-mapped; 22,023 frames = 5.8 GB
  [measured counts x arithmetic]. Pinned host if RAM allows (Q5), else mmap + a pinned bounce buffer of 2 x 8.4 MB.
- **Why not PNG per step**: PNG has no GPU decoder [web: nvImageCodec]; decoding once offline removes the only per-step
  CPU image work (~0.1 s [SP: B §3.1]) and the PIL objects TRL deep-copies.
- **Why not GPU-resident**: 5.8 GB vs ~7 GB headroom (§2.3). Change when on a 2nd GPU or a smaller crop.
- **Rule**: only if the 3 channels are equal for every frame (gate in build: assert R=G=B; 5 of 5 checked). If not,
  store 3 channels (17.3 GB, mmap only).

### 5.2 Doctor store

- **Key** = (item, object, view) -> device handle; preload all at start in one pooled allocation (fixes the epoch-1
  `cudaMalloc` + wait and the device-wide `cudaFree` stalls [code: gpucheck.cpp:138-149, 1741; SP: B §5 B7]).
- **Size**: 1.02 GB measured for 2 views as u8 masks + f32 distance; ~1.5 GB for 3 [estimate].
- **Phase D option (GPU-first and smaller)**: keep only bit planes (raters, uni, core, vote: (k+3)/8 B/px ~ 0.12 GB)
  and compute the two distance maps per step on the GPU for the 32 items in the batch (32 x <= 40,000 px, < 1 ms).
  Exactness: integer squared distances (separable exact EDT), then sqrt in float64 rounded to float32 = the
  correctly rounded value scipy produces [code: reward/cache.py docstring]. On the GPU use `__dsqrt_rn` then
  `__double2float_rn`, compile without fast-math. Gate G3 compares all tables to the scipy oracle.
- **Caches that pay** (caching skill: only with a hit-rate number): prompt ids + rope positions (1 distinct prompt per
  run for LIDC, hit rate ~100 % [SP: B TL;DR 5]); vision features within a step (each image used by 8 answers x 3
  passes). Caches that do not pay: cross-step features (§2.3). Versioned keys: tables carry the reward fingerprint and
  SHA-256 [code: reward/cache.py]; a new reward = a new spec name, never an in-place change.

### 5.3 Artefacts (blob-store sweep)

| Artefact | Size | Where |
|---|---|---|
| LoRA adapter (fp32) | ~70 MB | run dir -> HF private repo every 50 steps + final |
| Answers (all 8 per item, eval) | 17,600 x ~1 KB ~ 18 MB (zstd ~ 3 MB) | jsonl.zst beside metrics |
| Metrics + CIs + run record | < 1 MB | json |
| nsys report (2 steps) | ~100-500 MB [estimate] | local only |
| Step logs | ~2 KB x steps | jsonl (no log pipeline needed: one box; distributed-logging deferred) |

---

## 6. Deep dive: GPU vs CPU placement

### 6.1 Placement table

| Stage | Today | GPU-first placement | Why / "unavoidable" test |
|---|---|---|---|
| DICOM / PNG read and decode | CPU per step (32 PNG) | CPU, once, offline | No GPU PNG decoder [web]; file I/O is host-side |
| Patient split, manifest, hashes | CPU | CPU | metadata, KBs |
| View transform (shift / rescale) of frames + masks | CPU, 6 procs | GPU (offline) | Pixel work; gate: bit-equal to the current builder, else keep CPU and record why |
| Target polygonization | CPU numpy [code: targets/build.py:93] | GPU polygonizer v2 | Owner rule; v1 already bit-exact on 195,809 cases [SP: B §3.2] |
| Target certification | CPU numpy [code: checker/certify.py:46] | GPU checker | same engine |
| Doctor tables (EDT, union / core / vote) | CPU scipy | GPU exact EDT (Phase D) | pixels; exact by integer method + correct rounding |
| Coverage / leakage audit | CPU | GPU for pixel stats (offset, size from masks), CPU for joins | only numbers come back |
| Target length count | - | CPU tokenizer | text; runs once |
| Chat template / tokenization | CPU x512 per step | CPU, once per run (memo) | 1 distinct prompt |
| Pixel build (rescale, normalize, patchify) | CPU x2, fp32 1.61 GB | GPU from 32 u8 frames | elementwise; must equal the CPU processor bit-for-bit (test) or feed vLLM the same tensor |
| Rope positions | CPU Python loop | GPU constant per run (+ closed form for answer tokens) | integer, bit-exact by test |
| Vision encode | GPU, 16 calls of 2 images | GPU, 1 call of 32 | fewer launches, keyed by frame index (no hash) |
| Generation | GPU vLLM | GPU vLLM awake | - |
| Answer text parse | CPU, ~2 ms, overlapped | **CPU** | vLLM delivers token ids to the host; a GPU DFA would add an H2D and save ~0 s [SP: B §3.1]. Revisit only if vLLM keeps ids on device |
| Raster, legality, scoring | GPU | GPU | - |
| Advantage (256 numbers) | CPU numpy | GPU (fused after scoring) | removes a round trip; numbers only either way |
| Reference pass | GPU | removed (KL 0) or GPU 64-row chunks | owner |
| Update | GPU | GPU | - |
| Weight push | GPU via Python per-tensor | GPU in place | - |
| Evaluation metrics (IoU set, F, HD95, boundary IoU, GED, HM-IoU, medoid) | CPU pool for some | GPU; masks never leave device | memory rule: pixels never return to CPU |
| Bootstrap CIs | CPU | CPU | numbers only, 2,200 items |
| Logs, checkpoints, upload | CPU | CPU | I/O |

### 6.2 Pixel path (the largest CPU removal, -3 s)

32 frame indices -> copy stream gathers 8.4 MB pinned -> kernel: grey->RGB broadcast, `x * (1/255)`, normalize with
the processor's fused mean/std, temporal duplicate (temporal patch 2), reshape to [32 x 1,024, 1,536] -> vision
encoder once -> features + DeepStack features kept for every pass of the step. Rows map to images by an index list,
so the 256-row tensor is never built. vLLM gets the same 32 images with `n=8` (encodes 32, not 256) or the features
themselves via `--enable-mm-embeds` with DeepStack included [web: discuss.vllm.ai/t/2294; U for Qwen3-VL bugs #31708].
Feeding vLLM the trainer's features would make the image path identical in both by construction, which removes one
source of the trainer-vLLM mismatch (IS ratio min = 0 in 98.1 % of steps [SP: BRIEF §2]).

### 6.3 Generation

vLLM in process, awake (R4: awake = asleep in time and fits [SP: B §4.2]); `n=8` per prompt; prefix cache on; FULL +
PIECEWISE CUDA graphs (already on [SP: B §4.4]); `detokenize=False` (TRL decodes from ids) [U]; `max_tokens` = gate
G5 number. Optional (stage 1c): n-gram speculative decoding, lossless, aimed at the 1-2-sequence tail [SP: G §0 row 4;
U: acceptance on digit tokens].

### 6.4 Trainer

Packed shared-prefix forward (exists, R5 -24 %) with micro-batches sized by a token budget (e.g. ~8 k packed tokens
= ~6 groups of 1,373 tokens) instead of 16 rows: same dr_grpo update because the normaliser is constant
[SP: B §4.2 item 4]. Note the 2026-11-06 result "32x8 slower than 16x16" was on older code [SP: BRIEF §2]: re-measure,
do not assume. Attention: keep SDPA with the boolean mask now; FlexAttention block masks only after an sm_120 test
[U]. Partial checkpointing: drop checkpointing on k of 28 layers while peak <= 30 GB. DeepStack and fused-head indices
precomputed per forward (removes ~300 `nonzero` syncs [SP: B §2 rows 6, 9]).

### 6.5 Sync budget

Every `.item()`, `.tolist()`, boolean index, pageable `.to(device)` and `cudaFree` inside the step is a defect.
Allowed: vLLM's per-decode-step sync (inherent), 1 verifier result wait, 1 log read per step. The NVTX timer itself
must use events (`cudaEventRecord`), not `torch.cuda.synchronize` per phase, except in profiling mode.

---

## 7. Key decisions and trade-offs (solves / worsens / change when)

### 7.1 ADR-SD1-1: one process, colocated, synchronous, vLLM awake

| Option | Solves | Worsens | Change when |
|---|---|---|---|
| **A (chosen): one process, vLLM awake, sync** | Removes 1.66 s handover; zero policy lag; simplest debugging | One fault domain: a CUDA fault in any part kills all; memory pressure (~24 GB) | Peak > 30 GB in the worst-case test -> sleep level 1 for KV only |
| B: two processes on one GPU (MPS), async 1-step | Overlaps decode with update (~max(gen, train)) | NCCL cannot pair 2 ranks on 1 GPU, needs CUDA-IPC push; 1-step staleness; memory of two contexts | Owner allows staleness and stage 1b is not enough |
| C: partial rollouts | Kills the tail | Off-policy tokens inside an answer; not in TRL | Never on this budget |

### 7.2 ADR-SD1-2: in-place merged-weight push (not vLLM LoRA serving)

| Option | Solves | Worsens | Change when |
|---|---|---|---|
| **Merged slices written in place (chosen)** | ~20-50 ms vs 1.27 s; no allocation; decode at full bf16 speed; graphs stay valid | Needs a startup map of fused qkv / gate_up slices; must stay bit-equal to `lora_merged_params` (bf16 W + bf16 delta) | vLLM changes parameter layout -> contract test fails at start |
| vLLM native LoRA | Pushes only 35 MB; standard API | LoRA kernels slow decode 10-30 % [SP: B §4.6, estimate] | 2 GPUs over PCIe without P2P (stage 2) |
| vLLM weight-transfer IPC backend | Supported API for colocated RL [web: vLLM docs] | Full-weight path, more copies | It gains an in-place LoRA-merge mode |

### 7.3 Other decisions

| # | Decision | Solves | Worsens | Change when |
|---|---|---|---|---|
| T3 | u8 grey frame store, pinned host | -3 s CPU phase; 8.4 MB H2D | 5.8 GB host RAM; build step | Frames differ by channel -> 3-channel mmap |
| T4 | Doctor store preloaded, resident | No stalls; one handle lookup | ~1.5 GB of headroom | Memory tight -> bit planes + GPU EDT (0.12 GB) |
| T5 | Parse on CPU | Zero H2D of text; 2 ms hidden | One CPU stage remains (by necessity) | vLLM exposes device-side ids |
| T6 | Token-budget micro-batches | Halves the ~4.5-5 s fixed update cost | More activation memory | OOM -> smaller budget |
| T7 | KL 0 (owner) | -1.4-1.7 s, no reference model | Loses a drift brake (exact reward makes reward hacking unlikely) | KL spikes / collapse in entropy alarms |
| T8 | Cap from measured target length | Tail 4-5 s -> 1-2 s | Answers above cap become masked | Count shows targets near the cap |
| T9 | Verifier in process, versioned function | No serialization; 5 ms | Same fault domain as the trainer | A second consumer (verl, data filter at scale) -> HTTP service (Architecture.md §4.2) |
| T10 | TRL + adapter on 1 GPU; verl at 2+ | No migration cost now; 9 private names guarded by a start-up test | TRL upgrades are planned jobs | 2 GPUs, or TRL async gains images |
| T11 | Stage 1b: group-streamed update during the decode tail | Hides ~2-3 s: finished groups (all 8 done) run fwd/bwd while stragglers decode; same math (dr_grpo normaliser constant, gradients add) | Needs step-wise driving of vLLM's engine; float summation order changes (not bit-identical) | NVTX after the cap still shows tail > 2 s |
| T12 | No FP8 rollouts | Keeps the mismatch small | Leaves 1.3-1.6x decode on the table [estimate] | B2 has fixed the mismatch and decode dominates |
| T13 | Qwen3-VL-2B now, Qwen3-VL-4B as scale point | Same kernels and patches | 4B doubles weights (8.5 GB x 2 copies) -> memory forces sleep or 2 GPUs | Phase F |

---

## 8. Failure modes (resilience-failure)

One GPU = one process = one fault domain. Stale-beats-error does **not** apply to rewards: a wrong reward is worse
than a stopped run, so the verifier fails closed (error -> the step stops, run resumes from a checkpoint).

| # | Failure | Effect | Detection | Mitigation |
|---|---|---|---|---|
| FM1 | OOM after vLLM awake + bigger micro-batches | run dies mid-step | preflight worst case (256 x cap) | token budget from the measured peak; `expandable_segments`; fallback sleep level 1 |
| FM2 | Weight push writes to a stale or reallocated tensor; graphs read old weights | silent wrong policy | after push, assert data pointers unchanged; probe prompt log-probs vs trainer | in place only; GPU hash of vLLM slices vs fresh merge every 50 steps |
| FM3 | bf16 merge drift | silent drift of the base | same hash check | no merge/unmerge on the base (exists [SP: B §5 B10]) |
| FM4 | Checker kernel fault (illegal address) | CUDA context poisoned, whole process dies | exit code, heartbeat | host-side `check_*` before every launch [code: gpucheck.h]; sanitizer gate; checkpoint every <= 50 steps; tested resume |
| FM5 | Device-wide stall from `cudaMalloc` / `cudaFree` in the step | 10-100 ms stalls | allocator stats = 0 growth after step 1 | preload pools; caps from config |
| FM6 | CPU != GPU numbers (fast-math, FMA contraction, float atomics) | exactness claim breaks | golden files on the rented GPU | no `--use_fast_math`; integer-only or `--fmad=false` on float paths; no float atomics |
| FM7 | GPU pixel build != vLLM's pixels | larger IS mismatch | log-prob diff metric | one pixel tensor for both, or bit-equality test vs CPU processor |
| FM8 | Runaway decode tail | +4 s/step | share at cap > 5 % | cap from G5; stage 1b |
| FM9 | Trainer-vLLM mismatch masks most answers | lost signal | IS share at 0 > 0.5 % | B2 fix (token-level cut at 2, or fp16 both) |
| FM10 | Data leakage / train-test coverage gap | false gains or misses | gates G1, G4 | build fails closed |
| FM11 | TRL private name changes | silent wrong trainer | start-up contract test (9 names) | TRL pinned 1.17.1 |
| FM12 | Rented host dies / is reclaimed | lost time and $ | heartbeat | checkpoint 50 steps, HF push, `--resume`, cost cap per run |
| FM13 | Verifier stream steals time from the trainer | slower step | NVTX overlap | make the trainer stream high priority (negative); the checker's "lowest" = default today [SP: B §5 B6] |

---

## 9. Observability (what to measure)

- **NVTX ranges** (RED/USE frame of a step): `load`, `pixels`, `encode`, `push`, `decode`, `tail`, `parse`, `score`,
  `ref`, `update/mb[i]/{fwd,bwd}`, `optim`, `log`; in C++: one range per checker API call and stage.
- **Step log alarms** (from Architecture.md §8.2, kept): prepare > 0.5 s; IS share at 0 > 0.5 %; zero-spread groups
  > 25 %; share at cap > 5 %; entropy on shifted views < 0.03; far misses > 20 %; missing version -> stop.
- **GPU-first KPIs** (new): fully idle seconds per step (gaps > 1 ms in the nsys timeline) < 1.5; host syncs < 200;
  H2D < 50 MB; D2H < 1 MB; allocator growth after step 1 = 0; SM-active % during `tail` (nvidia-smi "100 %" lies at
  batch 1 [SP: B §1]).
- **Capture**: steps 4-5 only, `nsys --capture-range=cudaProfilerApi`, `--cuda-graph-trace=node`; ncu only for the
  polygonizer and EDT kernels.

---

## 10. Growth path

### 10.1 Stages

```
 Stage 1 (now): 1 GPU, one process         Stage 2: 2 GPUs (PCIe, no NVLink, P2P usually off)
 +-------------------------------------+   +------------------+  LoRA delta 35 MB  +----------------------+
 | vLLM | trainer | verifier | stores  |   | GPU 0: vLLM      | <================= | GPU 1: trainer +     |
 | in turns, 12-16 s/step [estimate]   |   | + base copy for  | =================> | verifier + doctor    |
 +-------------------------------------+   | in-place merge   |  ids+logps ~2.5 MB | store; frame store   |
                                           +------------------+  (via host)        +----------------------+
                                            verl async (CISPO/GSPO + IS correction), 1-step staleness

 Stage 3: 8+ GPUs                           k vLLM GPUs -> trainer (FSDP2 / Megatron-Bridge LoRA)
                                            verifier as HTTP batch service (same function, same bits)
```

| Stage | Trigger | Transfers per step | Notes |
|---|---|---|---|
| 1 | now | ~11 MB H2D, ~25 KB D2H | §2 |
| 1b | NVTX tail > 2 s after cap | same | group-streamed update |
| 1c | still decode-bound | same | n-gram speculative decode; per-layer CUDA graphs for the update [NVIDIA skill: nemo-mbridge-perf-cuda-graphs] |
| 2 | Phase C done and generation > 30 % of step; or 3D | LoRA bf16 35 MB GPU1 -> host -> GPU0 (~3-6 ms at PCIe 5 [estimate]) | GPU 0 keeps a frozen base copy and merges in place (no LoRA kernels in decode) |
| 3 | larger model or 3D at scale | NCCL broadcast | verl or NeMo-RL [SP: D_rl_infra §0] |

### 10.2 3D voxels (Phases G-H): what the numbers say

| Quantity | 2D now | 3D whole nodule | Basis |
|---|---|---|---|
| Slices per answer | 1 | mean 5.8, p90 11, max 52 | [measured] |
| Answer tokens | mean 131 | mean ~760, p90 ~1,440, max ~6,800 | 131 x slices [estimate] |
| Prompt image tokens (512x512 slices) | 256 | mean ~1,490, p90 ~2,820 | 256 x slices |
| KV for 32 prompts at p90 | 1.2 GB | ~10 GB prefixes alone | 112 KB/token |
| Verifier work | ~5 ms | 3D raster + run-based 3D labels + Euler (b0, b1, b2) on bit-packed volumes; max crop 200x200x52 ~ 260 KB bit-packed | still ms [estimate] |

Consequences: (a) on one GPU, 3D needs 8-16 prompts per step or smaller crops (a 256x256 slice crop = 64 image
tokens, 4x fewer); (b) compact codes (centre + 16 radii per slice, ~60 % fewer tokens, Direction F4/H3) matter more
than any kernel; (c) the decode share grows, so 3D is the natural trigger for stage 2. The engine design carries
over unchanged: bit words, run labels, Euler number by integer face counts, no marching cubes (Architecture.md D9),
doctor store as bit-packed volumes, EDT in 3D by the same separable exact method.

### 10.3 Physical AI (Phase I, after E): the same verifier with a mesh3d shape kind; out of scope here.

---

## 11. Coverage sweep

| Block | Status |
|---|---|
| blob-store | §5.3 (artefacts to HF private repo) |
| sequencer | deterministic ids `run/step/group/k`, keys `dataset/split/item::view` |
| dns / load-balancing / content-delivery | deferred: no clients, one box |
| observability | §9 |
| distributed-logging | deferred: one box, ~2 KB/step jsonl |
| distributed-search | not relevant |
| task-scheduling | the run spec's ordered steps (exists as YAML `steps:`); E6 typed spec |
| sharded-counters | not relevant |
| messaging-streaming | deferred to stage 3 (rollout queue between generators and trainer) |
| consistency-coordination | on-policy = strong; at stage 2 bounded staleness 1 with IS correction |

---

## 12. Score and diagnosis (quality bar, 0-5)

| Dimension | Score | What would raise it |
|---|---|---|
| Requirements | 4 | Owner answers to Q1-Q5 (sync vs async, KL, cap) |
| Numbers | 4 | Measured memory peak with vLLM awake + token budget; nsys split of the 4.65 s and 7.67 s phases |
| Building blocks under stress | 4 | vLLM internals (per-step sync, embeds path) are [U] |
| Trade-offs | 5 | - |
| Failure | 4 | One fault domain is accepted, not removed; MPS isolation would remove it at a cost |
| Growth | 4 | 3D token and KV numbers are estimates from slice counts, not from real 3D targets |
| **Weakest** | **Numbers on the trainer side** | The update's fixed 4.5-5 s and the memory headroom come from a two-point fit and nvidia-smi, not a trace. The first GPU session must be the NVTX capture. |

Current bottleneck after the design: decode (~7.5 s of ~13 s), half of it per-step overhead in the 1-2-sequence
tail. At the next order of magnitude (3D, or 10x answers per step) generation dominates and one GPU stops being the
right unit: disaggregate (stage 2). This design is a hypothesis; the NVTX capture can overturn §2.4.

---

## 13. Ten checks every implementation plan must pass

1. **Data gate before spend.** Patient-level leakage = 0 across all splits and views (G1); every RL row has >= 2
   doctors and a GPU-certified legal target (G2); the train offset / scale / size histograms cover the test histograms
   (G4: e.g. train share with offset > 0.3 >= test share of 37.5 %); the true target-length P100 is counted and sets
   the cap (G5). The build fails closed.
2. **Exactness on the real GPU.** Polygonizer v2, checker, scorer and GPU EDT are bit-identical to the numpy / scipy
   oracles on the 195,809 golden cases and on all LIDC and RIGA masks, run on sm_120; target hashes byte-identical;
   no fast-math, no float atomics, FMA contraction off on float paths; numpy code lives only in tests as the oracle.
3. **Transfer budget.** A per-step hook measures H2D < 50 MB, D2H < 1 MB, and asserts that no mask or pixel array
   reaches the host in training or evaluation (memory rule: only scores leave the GPU).
4. **Measured, not claimed.** An NVTX + nsys capture exists before and after each speed change; each change states
   its measured delta; fully idle < 1.5 s/step, host syncs < 200/step, prepare < 0.5 s.
5. **Equivalence.** Each speed change reproduces step-1 reward, loss and grad norm (bit-identical where claimed,
   else a stated tolerance); vLLM parameters equal a fresh merge (GPU hash) every 50 steps; GPU pixels equal the
   CPU processor's pixels or vLLM receives the trainer's tensor.
6. **Memory and allocator.** A worst-case preflight on the rented card (256 answers at the cap, vLLM awake, chosen
   token budget) peaks <= 30 GB; allocator growth after step 1 is 0; no `cudaMalloc` / `cudaFree` in the step.
7. **Signal health alarms.** The run logs and alarms on: IS share at 0 (< 0.5 %), zero-spread groups (< 25 %), share
   at the cap (< 5 %), far misses per offset bucket, entropy per view; a missing version field stops the run.
8. **Fault containment and resume.** compute-sanitizer memcheck / racecheck / initcheck / synccheck report 0 errors
   on every new kernel; host-side input checks precede every launch; checkpoints every <= 50 steps; a kill-and-resume
   test passes before the long run.
9. **Fair, reproducible evaluation.** Evaluation builds no trainer, saves all 8 answers' text, uses label-free medoid
   selection, nodule-level bootstrap CIs and the fair (k-1) ceiling; every artefact records code SHA, checker build
   hash, reward spec and data manifest SHA-256.
10. **Structure and spend.** Layers respect the import rule (geom < verify < data < policy < rl/sft < infra); old code
    is archived, dead code removed; the TRL 9-name contract test runs at start; `make test` and `make lint` pass; the
    plan's GPU total (from measured s/step) is <= $3, hosts are US/CA only, and every rental follows a passing
    preflight and the owner's approval.

---

Sources (web, checked 2026-11-08): [vLLM weight transfer](https://docs.vllm.ai/en/stable/training/weight_transfer/);
[RTX 5090 specs, TechPowerUp](https://techpowerup.com/gpu-specs/details.c4216);
[RTX 5090 review, Tom's Hardware](https://www.tomshardware.com/pc-components/gpus/nvidia-geforce-rtx-5090-review);
[nvImageCodec](https://www.github.com/NVIDIA/nvImageCodec);
[vLLM precomputed multimodal embeddings](https://discuss.vllm.ai/t/how-do-i-precompute-multimodal-embeddings/2294);
[vLLM image_embeds issue](https://discuss.vllm.ai/t/issue-unable-to-pass-precomputed-image-embeddings-to-vllm/2203);
[verl v0.8.0](https://newreleases.io/project/github/verl-project/verl/release/v0.8.0);
[NVIDIA skills catalog](https://raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json);
[NeMo-RL IPC weight sync](https://docs.nvidia.com/nemo/rl/nightly/apidocs/nemo_rl.weight_sync.ipc_weight_synchronizer.html).
