# H4: clean-code and algorithm review of ObjectDelineation (2026-11-08)

Scope: `src/objdel` (about 7.3k lines), `cuda/gpucheck.cpp` + `gpucheck.h` + `gpucheck_py.cpp`, `experiments/phase2/*.py`,
`experiments/phase2/configs/*.yaml`, `tests/` (71 files, 957 tests), `Makefile`, `CMakeLists.txt`, `pyproject.toml`, CI.
I only read code. I made no edits, changed no git state and spent nothing. Following the clean-context rule, I did not
read transcripts, memory, `.claude/`, or `docs/research/notes/*review*`. Branch: `gpu-scorer-open-once` at `be52b1c`.

---

## Skills used

| Skill | What it added to this review |
|---|---|
| `mattpocock-skills:codebase-design` | The words used for section 1 and section 4: module, interface, depth, seam, adapter, leverage, locality, and the deletion test. Shallow and deep modules are graded with it. Also the rule "one adapter is a hypothetical seam, two adapters make it real", which shaped the TRL-adapter proposal. |
| `engineering:tech-debt` | Debt categories (code, architecture, test, dependency, docs, infra) and the scoring formula Priority = (Impact + Risk) x (6 - Effort), used in the table in section 3.9. |
| `engineering:testing-strategy` | The testing pyramid and the "business-critical paths first" lens behind section 3.6: GPU kernels, the TRL contract and the statistics of the reported numbers. |
| `engineering:architecture` | ADR-style trade-offs for the target layout in section 4: options, consequences, and what to revisit later. |
| `engineering:code-review` | The review lenses (correctness, edge cases, concurrency, performance, maintainability) that order section 3. |
| `engineering:debug` | Reproduce, isolate, diagnose. I used it to check each correctness suspicion by tracing the real call path instead of guessing. For example, `final_eval` bootstraps items one by one while `signal.draws` bootstraps by cluster, and `train_rl` creates a fork pool that it never uses. |
| `engineering:system-design` | Requirements and data flow for the reward and eval path: what runs per step, what runs offline, and which host-device round trips matter. Used to rank the hot spots by how often they run. |
| `engineering:documentation` | "Keep it current" and "link, don't duplicate". This found doc drift (docstrings name TRL 1.16, but TRL 1.17.1 is pinned) and the weights that are copied by hand (section 3.8). |
| `system-design-skills:service-decomposition` | The modular-monolith verdict. Keep one deployable with enforced internal seams. Do not split the verifier into a network service until a second consumer exists (section 4). |
| `nvidia-skills:nvidia-skill-finder` | The main thread checked the catalog (`MAIN_skill_lookups.md`). Candidates are `nv-segment-ct` and `nv-segment-ct-finetune` as CT segmentation baselines, the `dicom-series-*` skills for data ingest, and `tao-generate-referring-expressions`. None of them changes the code findings. Recommendation only: nothing was installed. |
| `huggingface-skills:huggingface-best` | From the main-thread lookups: Qwen3.5-2B and 4B (Apache-2.0, RefCOCO avg 88.2 for 2B) are the natural next policy models. The design point that follows is in section 4: the policy adapter must not hard-code Qwen3-VL internals (`get_rope_index`, `visual.merger`, deepstack features), as `fast.py` does today. |

---

## 0. Top 12 findings

| # | Finding | Where | Severity |
|---|---|---|---|
| 1 | The training builder lives in an experiment script. `train_rl`, `final_eval` and the G-probe all build their trainer through `experiments/phase2/g_probe.py:74 make_trainer`. Its interface is a duck-typed `argparse.Namespace` read with `getattr` defaults. | g_probe.py:74-101, train_rl.py:41,53, final_eval.py:62,72 | High (architecture) |
| 2 | Evaluation builds a full GRPO trainer (LoRA wrap, optimizer config, vLLM colocate) and needs the train set (`--train` is required) just to sample. It then calls the private `trainer._generate`. | final_eval.py:51,69-75; sample.py:64 | High |
| 3 | `final_eval.summarize` bootstraps items one by one. Every other CI in the repo bootstraps by nodule cluster (`signal.draws`). The test CIs reported for the full RL run are therefore probably too narrow. | final_eval.py:33-43 vs signal.py:98-107 | High (statistics) |
| 4 | Nine TRL private names are overridden or patched, plus a module-global swap of `trl_grpo.shuffle_sequence_dict`. If TRL renames one, the override silently stops running. There is no startup check that the names exist. | fast.py:102,524,549,561-565,567,575,618,626,638; sample.py:64 | High (dependency) |
| 5 | The production RL learning rate comes from a pre-registration file of the G-probe experiment (`prereg.G_PROBE["lr_lora"] = 5e-5`). It overrides `configs.grpo_config`'s own 1e-5. So there are two sources of truth. | g_probe.py:86, prereg.py:61-65, configs.py:189 | Medium-High |
| 6 | `train_rl` forks up to 16 scoring workers holding copies of every preloaded cache, then never uses them (only `pool.shutdown()`). | train_rl.py:50,61 | Medium (waste) |
| 7 | Library dead weight. `callbacks.Controller` (about 130 lines) has no production caller, only tests. The GPU mask polygonizer `_gpucheck.check_masks` has no `src` caller, because the offline target builder uses the CPU numpy `checker.polygonize`. Two polygonizers are kept bit-identical by hand. | callbacks.py:81-209; targets/build.py:53,93; gpucheck.cpp:966-1058 | Medium |
| 8 | GPU-first gaps. Eval metrics pull masks back to the host (`to_masks`) and run IoU, NSD (scipy EDT), GED and HM-IoU on a CPU fork pool. The target builder runs one checker call per ladder rung per object, in sequence. | eval/metrics.py:249-270; targets/build.py:215-230 | Medium |
| 9 | The scorer always runs 4 union-find labelings plus 2 filters, even when the size filter changes nothing. The CPU path skips that work (`signature` fast path). | gpucheck.cpp:1612-1614 vs checker/topology.py:157-163 | Medium (perf) |
| 10 | No CI job builds or runs the nvcc path. The GPU raster (`raster_row_warp`, warp shuffles and ballots) is a different implementation from the CPU raster that CI tests. It is checked only on rented boxes, from config steps. | gpucheck.cpp:381-461; .github/workflows/tests.yml; configs/rl_full.yaml | Medium-High (test) |
| 11 | The weights are copied by hand: `signal.WEIGHTS` restates `REWARD_V4.weights`, and `grpo_config.reward_weights` restates the order of `TRLReward.funcs()`. | signal.py:31; configs.py:191; rewards.py:109 | Medium |
| 12 | Config sprawl: 16 YAML files (378 lines) that are really argv shell scripts. Hyperparameters live in CLI strings, with no schema, and most of each file is shared between files. | experiments/phase2/configs/*.yaml | Medium |

---

## 1. Module map

### 1.1 What exists

| Module | Lines | Job | Depth (interface vs behaviour) | Depends on |
|---|---|---|---|---|
| `cuda/gpucheck.{h,cpp}` | 150 + 1835 | Single-source C++/CUDA engine with five operations: polygon raster + 13 defect classes; mask -> rings (GPU polygonizer); resident expert tables; per-object scoring (IoU_tc, F5, g_area, GPU signature). | **Deep.** 8 functions in the header, about 1.8k lines behind them, one mutex. A good seam. | CUDA runtime only |
| `cuda/gpucheck_py.cpp` | 227 | nanobind adapter: validates the caps, releases the GIL. Accepts only `device::cpu` arrays. | Thin adapter (correct for its role) | gpucheck.h |
| `objdel.checker` | ~500 | CPU exact topology (OpenCV CCL), bit-quads, numpy polygonizer, certificate, 3-D Betti numbers. | Medium. `polygonize` duplicates the GPU engine. `betti3d` is used only by phase0. | numpy, cv2 |
| `objdel.geometry` | ~330 | Text grammar parse/serialize; canonical rings; packing for the engine; `to_masks`. | Medium, a good seam. Parsing is regex then `json.loads` (two passes, text.py:116-118). | `_gpucheck` |
| `objdel.reward` | ~600 | `ExpertCache` (offline rater statistics plus the device upload), `score_many/score_batch`, `advantages_v4`, `REWARD_V4`. | **Deep.** `score_batch(jobs) -> dict` hides parse, pack, the GPU call and the joint topology. The core verifier. | geometry, checker, `_gpucheck`, scipy |
| `objdel.targets` | ~820 | Frames, augmentation, simplification (numba), certified SFT target builder. | Deep, but it mixes CPU polygonize with GPU certification. | checker, geometry |
| `objdel.data` | ~1040 | LIDC/RIGA parsing, splits, Parquet store. | OK | pylidc, pyarrow |
| `objdel.eval` | ~1750 | Item metrics, signal/bootstrap, gates, audit, coverage, vision check, **synthetic probe generation (616 lines)**. | Mixed. `probes.py` generates synthetic CT data (`render`, `make_hole`, `build_precision`), which is data synthesis, not evaluation. | reward, geometry, scipy |
| `objdel.train.fast` | 732 | FastGRPO: windowed forwards, shared-prefix packing, VisionCache, RopeCache, ImageMemo, IS ratio, no-merge vLLM sync, vLLM prompt patch, timers. | **Shallow toward TRL.** Its interface is TRL's private method set, so every caller must know TRL internals. Seven separate concerns sit in one class. | TRL 1.17 privates, transformers, peft, Qwen3-VL internals |
| `objdel.train.rewards` | 160 | TRL reward adapter: weight-0 loggers plus the aggregator, prefetch thread. | Adapter. Fine, but it couples to `configs.reward_weights` by list order. | reward |
| `objdel.train.serve` | 177 | Global `STORE` of caches, `CacheRef`, fork pool, cgroup CPU math. | Process-global state is the interface (`STORE`, `_POOLS`). Hard to test in isolation. | reward |
| `objdel.train.dataset` | 392 | prepare(); SFT, GRPO and view rows; OpenOnce; cache_loader. | Medium. `cache_loader` duplicates `serve.preload`. | reward, targets, datasets |
| `objdel.train.callbacks` | 237 | Controller (audit, rollback, live switches), TraceStep. | Controller has no production caller. | TRL |
| `objdel.train.configs` | 196 | SFT/GRPO configs, LoRA, sampling names. | OK; the LR is overridden from outside. | TRL |
| `objdel.train.sample` | 158 | Sampling through `trainer._generate`; JSONL export/import. | Shallow, private-API coupled. | TRL privates |
| `objdel.secrets` | 38 | `.secrets` / keychain lookup. | Infra code in the library. | — |
| `experiments/phase2/g_probe.py` | 198 | G-probe **and the shared trainer factory** and `save_json`. | Hidden library | everything |
| `experiments/phase2/train_sft.py` | 440 | SFT orchestration: a 174-line `train()` with nested Trainer subclass, Rolling callback and closures. | Hidden library | TRL privates |
| `experiments/phase2/select_ckpt.py` | 426 | Sampler subprocesses, scoring pool, judge, HF check. Other scripts import `score_items` and `view_rows`. | Hidden library | train_sft, baselines, prereg |
| `experiments/phase2/vast/policy.py` | 628 | Cloud runner: offers, launch, sync, run chain, teardown, ledger. | Infra; fine in experiments, but tests import it 25 times. | hf_push, preflight |

### 1.2 The import graph shows a hidden library inside `experiments/`

```
train_rl ──> g_probe.make_trainer, g_probe.save_json
final_eval ─> g_probe.make_trainer, select_ckpt.score_items
g_probe ───> select_ckpt.{score_items, view_rows}, prereg.{G_PROBE, g_probe}
select_ckpt > train_sft.vision_check, baselines, prereg
preflight ─> train_sft, project_r1, vast.policy
project_r1 > train_sft;   hf_rate, r0_parity > select_ckpt;   k3_harness > phase1.reward_gates
vast.policy > hf_push, preflight
```

Tests import `experiments.*` in 35 of 71 files (79 imports of `experiments.phase2`). Experiment scripts are therefore
part of the tested API without being library code. The project rule ("library code lives in src/objdel and needs
tests") is broken in spirit.

### 1.3 Seams in the wrong place (deletion test)

1. **`make_trainer` (g_probe.py:74)**. If you delete it, its complexity comes back in 3 callers, so it earns its place,
   but it sits in the wrong package and has an implicit interface: `a.model, a.out, a.grpo, a.train` are required, and
   `a.lr, a.train_views, a.plain, a.trace, a.shared_prefix, a.compile` are optional through `getattr`. `final_eval`
   fakes it with `argparse.Namespace(model=..., out=..., grpo="{}", train=...)` (final_eval.py:71). It should move to
   `objdel.rl.build_trainer(policy: PolicySpec, data: RLData, run: RunSpec) -> Trainer`.
2. **Sampling exists in three implementations**: `sample.sample` through `trainer._generate` (sample.py:49-71),
   `experiments/phase2/sample_vllm.py` (a vLLM subprocess) and `hf_sample.py` (HF generate). One `Sampler` interface
   with vLLM and HF adapters is a real seam (two adapters). Evaluation should depend on it, not on a trainer.
3. **`serve.preload` vs `dataset.cache_loader`** both parse `manifest.json`, check SHA-256 and collect cache digests
   (serve.py:56-82, dataset.py:290-310), with different rules. `preload` globs `<ds>/*.jsonl` (serve.py:67) and raises
   on a stray file; `cache_loader` walks the manifest (the comment cites I45 "no glob"). Merge them into one
   `PreparedRoot.open(root) -> {key: digests}`.
4. **FastGRPO is one class with seven jobs** (fast.py:478-693): forward windowing, packing, the vision cache, the rope
   cache, the image memo, the IS ratio, vLLM sync and patching, and timers. Each is already a separate class or function
   except the trainer glue. Make the trainer a thin TRL adapter that wires independent, model-agnostic parts. Qwen3-VL
   knowledge (`get_rope_index`, `visual.merger`, `deepstack_features`, `mm_token_type_ids`) belongs in a policy adapter.
5. **Process-global `STORE`** (serve.py:51) is the interface between preload, the fork pool, `CacheRef` pickling and the
   reward. It works, but it forces the ordering rules (`preload` raises after `fork_pool`, serve.py:60) and makes
   parallel tests share state. A `CacheStore` object passed in explicitly keeps the fork copy-on-write trick (the
   object is still a module-level singleton at fork time) and gives tests a seam.

---

## 2. Algorithm hot spots

Frequency key: **S** = every RL step, **E** = every evaluation, **O** = offline (prepare or views).

| # | Where | Now | Better | Freq | Value |
|---|---|---|---|---|---|
| A1 | `gpucheck.cpp:1612-1614` score_words | Always: label fg, count, filter, label bg, count, filter, then relabel both, so 4 full union-find labelings per call. | Have Filter set a per-object "changed" flag; relabel only those crops (CPU `signature` already does this, topology.py:161). Most RL rollouts are clean blobs above min_area, so about 2 of the 4 labelings go away. | S | Medium |
| A2 | `gpucheck.cpp:1314-1334` LabelMerge, 1302 LabelInit, 1336 Flatten | One thread per **crop row**, serial along the row, atomics for union. Parallelism = total crop rows (about 256 objects x about 40 rows = about 10k threads), each looping over the width. | Pixel-per-thread or 2x2-block union-find (Playne-Hawick BUF / Allegretti BKE) on the bit-packed words, or one block per crop with shared memory. Crops are small, so one-block-per-object CCL in shared memory removes the global atomics. | S | Low-Medium (scoring overlaps forwards; see C9) |
| A3 | `gpucheck.cpp:1549-1640` score_words | 3 host round trips per call: boxes (sync 1574), terms (1626), signatures (1639). | Allocate the crop arena from the frame bound `(H + 2)(W + 2)` (frames are at most 512^2, so it is small), so trip 1 goes away. Copy signatures with a fixed cap (at most 10 rings means at most 10 components), so trip 3 goes away. One sync per call. | S | Low (about 30-100 us per call) |
| A4 | `gpucheck.cpp:585-623` EdgePairs | O(E^2) per shape; thread v loops over every later edge, so the load is triangular (the first edge does E work, the last does 0). Worst case at the caps is 10 rings x 1000 vertices = 10k edges, about 5e7 pairs per shape. | Flatten the pair space (one thread per edge pair, index into the triangle) or sweep by sorted x-interval. Fine for RL outputs (tens of vertices); it matters only for adversarial or long outputs. | S (tail) | Low |
| A5 | `gpucheck.cpp:1720-1763` add_tables | Host `pack_plane` loops pixel by pixel over k+3 planes (1167-1172, 1727); a synchronous `cudaMalloc` + copy + `cudaStreamSynchronize` per table, under the global mutex. Called lazily from `ExpertCache.device()` inside `score_many` (score.py:114). So in epoch 1 every new cache stalls the reward call (and the prefetch thread). | Batch-upload all preloaded caches once after `preload` (one arena, one copy). Pack the bits with numpy `packbits(bitorder="little")` or on the device. `drop_tables` -> `cudaFree` runs inside the GC finalizer (cache.py:256) and syncs the device at a random time during the step. Use a pooled arena and free at shutdown. | S (epoch 1) | Medium |
| A6 | `gpucheck.cpp:1808-1834` score_masks | Host per-pixel `pack_plane` (1826), then upload. | Pack on the device, or take already packed words. Used only by baselines and gates. | O | Low |
| A7 | `gpucheck.cpp:966-1058` check_mask_chunk | 3 host round trips plus a host `std::sort` of the rings (1032). | A device radix sort by (mask, key) and a device scan for offsets give one round trip. The kernel has no `src` caller today (finding 7). | — | Low |
| A8 | `targets/build.py:215-230` `_build` ladder | For each object, rung by rung: `_parts` -> `to_mask` (one checker call, full mask back to the host) -> CPU IoU -> CPU `kept()` (up to 2 OpenCV labelings per rule) -> next rung. Up to about 14 sequential GPU round trips per object, each with the mutex, a pinned copy and a sync. | Rungs do not depend on each other: build all rungs' parts, call `check_polygons` once for all rungs and all objects of a batch of items, compute IoU and signatures on the device (the GPU signature already exists inside score_words), and pick the first passing rung on the host. That is O(1) round trips per batch instead of O(items x rungs). | O | Medium (offline wall clock) |
| A9 | `checker/polygonize.py` (numpy) used by `targets/build.py:93` | CPU pointer jumping, O(E log E) work per mask, for every rater mask of every item and view epoch. | Use `_gpucheck.check_masks` (same output, already bit-identical, gpucheck.cpp:17-23) and delete the numpy copy, or keep it only as a test oracle under `tests/oracles/`. The owner rule ("CUDA for the polygonizer") then holds in the library. | O | Medium (rule + duplication) |
| A10 | `eval/metrics.py:249-270, 172-235` | Masks come back to the host (`to_masks`). Pairwise IoU as a float64 dense matmul over cropped masks; NSD with one scipy EDT per distinct sample plus one per rater; GED; HM-IoU. Runs in a CPU fork pool (serve.fork_pool). | Keep masks on the device: pairwise IoU = popcount(AND) on packed words (the same primitive as ScoreRows), NSD from the device EDT of the boundary (Felzenszwalb 1-D passes) and the corner-grid surfels. Only HM-IoU (Hungarian on K x k, tiny) stays on the host. Complies with "pixels never return to CPU". | E | Medium |
| A11 | `reward/cache.py:222-238` ExpertCache.__init__ | On every load or unpickle: 2 scipy EDTs, k signatures (OpenCV), vote, box. `preload` builds them all with 4 threads (serve.py:80-81). | Store the derived tables in the npz (they are deterministic; the digest already covers the inputs), or derive them on the device. Startup only. | O/startup | Low |
| A12 | `train/dataset.py:193-200, 206-207` sft_rows -> `_voting` | Builds a full `ExpertCache` (EDTs + signatures) for **every record**, only to read `voters`. | Store `voters` in the prepared record (it is fixed by the reward fingerprint), or compute signatures alone (`voters(signature(...))`). | O | Low-Medium |
| A13 | `train/fast.py:62-78` prompt_groups | `torch.stack(pix).flatten(1)` copies every row's pixel block per chunk per pass (old/ref/update) and compares element by element: O(rows x pixels) extra memory and traffic, although rows of one group share a single image object (OpenOnce). | Compare the VisionCache keys (already computed) or the `id()` of the image from ImageMemo plus the grid; or compare one representative row per run. | S | Low |
| A14 | `train/fast.py:174-204` packed_chunk | Builds the per-token index lists `gi, pi, ri, ci` with Python `+=`/`range` (O(tokens) Python ints per chunk; about 256 rows x about 130 answer tokens = about 33k ints, three passes per step), and a dense (c, W, W) bool mask (W can reach about 1.4k, so about 2 MB per group). | Build the indices with `torch.repeat_interleave` / `cumsum` from `lens` on the device. The dense mask is what makes SDPA fall back from flash attention. Use the FlexAttention block mask, or varlen flash-attn with per-segment `cu_seqlens` plus a prefix-sharing (cascade) kernel, so masked tiles are skipped. That is where most of the remaining forward time is. | S | Medium-High (perf) |
| A15 | `train/fast.py:277, 325, 339, 220` | Prefix sums as `[sum(x[:i]) for i in range(n)]`: O(n^2) Python (n up to 256). | `itertools.accumulate` / `np.cumsum`. Trivial. | S | Low |
| A16 | `eval/signal.py:102-107` draws | Cluster index lists built by `[j for j, x in enumerate(nod) if x == u]` for each cluster: O(n x U). | One pass with a `defaultdict(list)`. | E | Low |
| A17 | `geometry/text.py:116-118` parse | Regex `fullmatch`, then `json.loads`, then per-ring `np.fromiter`: two passes over every completion, and Python int objects per vertex. | A single-pass C++ tokenizer for the fixed grammar that writes the packed `xy/ring_start` arrays straight into the engine's input (it would also remove `pack_objects`' Python loops, shape.py:119-133). About 1-3 ms per 256 completions today, so low unless steps get much faster. | S | Low |
| A18 | `train/rewards.py:118-119, 141-143` prefetch | When the signature does not match, the discarded future still runs a full `score_batch` while the main thread scores again (serialized by the engine mutex, so the reward is computed twice). `ExpertCache.device()` check-then-set (cache.py:248) is not thread-safe between the prefetch thread and the main thread (a double upload, not a leak). | Cancel or ignore the future cheaply. Guard `device()` with a lock. | S (rare) | Low |
| A19 | `train/fast.py:642-650` timed | `torch.cuda.synchronize()` at the end of every timed phase. | Fine for profiling runs. Gate it behind a flag in production (it serializes the CPU-GPU overlap the prefetch design wants). | S | Low |

Already good, no change needed: `advantages_v4` is fully vectorized (advantages.py:95-105); `_pairwise_iou` dedups
masks (metrics.py:172-186); ScoreRows uses popcount over packed words (gpucheck.cpp:1411-1442); the polygon raster is
one warp per row with ballot-packed words (gpucheck.cpp:388-415); OpenOnce and ImageMemo remove repeated decoding.

---

## 3. Code-quality findings by severity

### 3.1 Correctness risks

1. **Item bootstrap in `final_eval`** (final_eval.py:33-43). `rng.integers(0, len(items), (draws, len(items)))`
   resamples items independently. Item ids are `lidc/<patient>/<nodule>/...`, and `signal.cluster` (signal.py:98-99)
   groups by nodule because slices of one nodule are correlated; `signal.draws` (102-107) is used for every other CI.
   The README-level test numbers (IoU 0.645 etc.) carry CIs from this function. Fix: `signal.ci(v, signal.draws(ids))`.
   Add a test that a split with duplicated-nodule items gives a wider CI than the item bootstrap.
2. **TRL overrides can silently stop running** (details in 3.2). If TRL renames `_get_per_token_logps_and_entropies`,
   `_prepare_inputs` or `_generate`, the FastGRPO methods become unused, with no error. Speed drops back, and C10
   (`skip_old_logps` turns off TRL's IS correction at fast.py:530-531 and expects the update forward to put it back
   through `_compute_loss`/`_get_per_token_logps_and_entropies`) could leave **no IS correction at all**. That is a
   silent off-policy bias, not just a slowdown. Add an import-time assertion that each overridden name exists on
   `GRPOTrainer` with the expected signature, and a test that `importance_sampling_ratio` reaches the loss.
3. **The VisionCache key is a hash, but the docstring calls it "exact"** (fast.py:242, 247, 262-272). It is a random
   linear combination of the bit patterns that wraps around in int64, and `place` multiplies again. If two images
   collide, the second gets the wrong vision embeddings with no error. The chance is tiny (about 2^-64 per pair), but the
   claim is wrong. Either key by exact bytes (like RopeCache, fast.py:321-327), or use a 128-bit digest and fix the
   docstring.
4. **The LR is defined in an experiment pre-registration** (g_probe.py:86 `G_PROBE["lr_lora"]`, prereg.py:64 = 5e-5)
   while `configs.grpo_config` says 1e-5 (configs.py:189) and its docstring says "GRPO learning rates: 1e-6 full,
   1e-5 LoRA". The full RL run used 5e-5. A reader of `configs.py` gets the wrong number.
5. **Duplicated reward constants**: `signal.WEIGHTS` (signal.py:31) restates `REWARD_V4.weights`, and
   `grpo_config.reward_weights` (configs.py:191) restates the function order of `TRLReward.funcs()` (rewards.py:109).
   `REWARD_V4` is called "frozen", but a change made in one place only would make the G-probe composite R and the
   training reward disagree silently. Derive both: `dict(zip(TERMS, REWARD_V4.weights))`, and have
   `funcs()` return the weights that configs consumes. (The current code also discards them: g_probe.py:91 `funcs, _ =`.)
6. **The identity cache in TRLReward** (rewards.py:130-131, `self._last[0] is completions`) depends on TRL passing the
   same list object to all 8 reward functions. If TRL ever copies the list, scoring runs 8 times (cost only). Key by
   `id()` plus the length plus the first and last elements, or by a step counter.
7. **`topology_valid` checks `topo_any > 0` and `topo_freq > 0`** (advantages.py:76). They are redundant (topo_any is
   `joint > 0`, score.py:140). Harmless, but it suggests two meanings.
8. **The GC finalizer frees device memory** (cache.py:256, `_drop_tables` -> `cudaFree`). `cudaFree` synchronizes the
   device, at any Python GC point, including inside the prefetch thread during a forward. It is not a correctness bug,
   but it is non-deterministic stalls; see A5.

### 3.2 TRL monkey-patching: what exists and how it breaks

Pinned to `trl==1.17.1` (pyproject.toml:27). Good. These are the items that break on an upgrade:

| Hook | Kind | File:line | What breaks |
|---|---|---|---|
| `_get_per_token_logps_and_entropies` (positional args `batch_size, compute_entropy, compute_aux_loss`) | override of a private method, called with super() using **positional** args | fast.py:102-107, 150-152, 575-579 | A rename leaves the override silently unused (speed, C10 IS ratio). Changed argument order means wrong arguments. |
| `fused_lm_head=True`, `out.log_probs`, `out.entropy`, `support_ids` | private forward protocol of TRL's patched model | fast.py:226, 232, 239 | AttributeError or wrong shapes |
| `_generate_and_score_completions` + toggling `self.vllm_importance_sampling_correction` | override + state flip | fast.py:524-538 | IS correction silently lost (3.1-2) |
| `_prepare_inputs` + swap of the module global `trl_grpo.shuffle_sequence_dict` | **module monkey-patch**, not thread-safe | fast.py:549-565 | If TRL calls `utils.shuffle_sequence_dict` or inlines it, the patch does nothing (ordering only, safe). The tests patch `_prepare_inputs` itself (test_train_fast.py:430), so they do not prove TRL calls the global. |
| `_compute_loss` | override | fast.py:567-573 | C10 path |
| `_generate` return tuple, read by **position** `out[3], out[1]` | private return protocol | fast.py:626; sample.py:64 | Silent wrong argument (completions vs ids) if the tuple changes |
| `_calculate_rewards`, `_metrics[mode][...]` | override, private dict | fast.py:638-640, 228, 609, 650 | KeyError |
| `vllm_generation._iter_named_params`, `_dist.is_fsdp`, `.llm.generate/.sleep/.wake_up`, `sync_weights` | instance monkey-patches of TRL's private vLLM helper | fast.py:498-504, 656-668, 675-690 | AttributeError at startup (loud, good) or a silent fall back to the drifting merge path |
| `trainer._generate` from sampling | private call | sample.py:64 | as above |
| `SFTTrainer._get_train_sampler`, `_load_from_checkpoint` | private overrides in a nested class | train_sft.py:322-345 | sampler order or resume silently changes |
| `processing_class.image_processor = ImageMemo(...)`, `visual.forward = ...`, forward pre-hook on the Qwen base | object surgery on transformers models | fast.py:430-441, 249-251, 301-303 | `__getattr__` forwarding hides missing attributes; `named_modules()` search by suffix `visual` + `merger` |

Mitigation, in order of value:
(a) one `objdel.rl.trl_adapter` module that owns every private touch point and checks at import time that each name
exists and that `inspect.signature` matches a recorded table, so an upgrade fails at once with a list of breaks;
(b) call super() with **keyword** arguments everywhere (fast.py:106, 150);
(c) replace `out[3], out[1]` with a named unpack helper used in one place;
(d) keep the existing TRL-equivalence tests (test_train_fast.py:93, 596, 612 are good contract tests) and add one
end-to-end test that `importance_sampling_ratio` is in the loss inputs when `skip_old_logps()` is true;
(e) longer term, upstream the windowing and shared-prefix packing to TRL, or own a small GRPO loop (the loss is about
50 lines with Dr.GRPO + clip-higher) and use TRL only for vLLM. Two adapters, TRL and a minimal own loop, would make
the seam real.

### 3.3 Dead or unused code and resources

- `callbacks.Controller` (callbacks.py:81-209): only tests construct it (tests/test_train_trl.py:190-395). `rewards.py`
  keeps hooks for it (`recorded`, rewards.py:97, 139-140; `trl_reward` attribute, 107). Either wire it into `train_rl`
  (it gives rollback and held-out audit, which a production run wants) or move it to `archive/`.
- `train_rl` fork pool (train_rl.py:50, 61): created, never used. Remove it (it saves 16 forked workers plus gc.freeze).
- `_gpucheck.check_masks` has no `src` caller; `checker.polygonize` has one (targets/build.py:93) plus `certify`. Keep
  one implementation (A9).
- `checker.volume.betti3d`: only `experiments/phase0/lidc_census.py:24` uses it. It can move with phase0.
- `dataset.prompt_rows` is used only by `grpo_rows` and one test; `grpo_rows` is a one-line alias
  (dataset.py:230-237). These are shallow modules; collapse them.
- `archive/cpu_checker` and `experiments/phase1/lab` (about 45k tracked lines, 24 byte-identical file pairs between
  `lab/t_*` and `lab/t_stack/vendor/*`). This is frozen evidence, excluded from lint (pyproject.toml:85). Move it to a
  tagged branch or a release artifact so the working tree is only live code.

### 3.4 Functions that are too big and names that are unclear

- `experiments/phase2/train_sft.py:263-437 train()` (174 lines) does validation, manifest and resume logic, a nested
  `Trainer` subclass, a nested `Rolling` callback, a nested factory, and a phase-by-phase checkpoint state machine.
  Split it into a `SFTRun` with `plan()`, `stable_phase()`, `decay_branch(name)`, `finish(dir)`, and a library
  `objdel.sft.Trainer`.
- `fast.py:164-239 packed_chunk` (76 lines, 13 parameters) and `fast.py:102-162` (61 lines, two code paths interleaved
  in one loop). Extract `pack_groups(...) -> PackedBatch` (pure, testable on the CPU) and `scatter_back(...)`.
- `eval/metrics.py:276-350 metrics()` (75 lines) mixes per-sample, per-object, per-item and per-group metrics.
- One-letter and abbreviated names on hot paths: `a` for the argparse namespace (g_probe.py:74), `T` for the terms
  dict (rewards.py:143), `P`, `W`, `c`, `k_`, `n_` (fast.py:174-204), `sl`, `v`, `L` (gpucheck.cpp). Short names are fine
  in kernels, but not as module interfaces (`make_trainer(a, ...)`).
- `fast.py:30-31`: `BUCKETS` vs `FINE`. Rename to `COMPILED_WIDTHS` / `EAGER_WIDTHS`.
- Repeated code: `by_length` and `by_group_length` each define the same `move` closure (fast.py:383-388, 399-404).
  `floating_point_ops -> 0` appears in both fast.py:692 and train_sft.py:319. There are three JSON writers:
  `g_probe.save_json` (184), `train_sft.write_json` (236), `select_ckpt.dump` (318). `_pixels_sha` exists twice
  (dataset.py:86, probes.py:369; the comment admits it).

### 3.5 Duplicated logic (locality)

| What | Copies | Keep one in |
|---|---|---|
| Polygonizer | `checker/polygonize.py`, `gpucheck.cpp:751-1058` | engine (GPU); numpy as a test oracle |
| Prepared-root manifest check and digests | serve.py:56-82, dataset.py:290-310 | `objdel.data.prepared` |
| Bootstrap CI | final_eval.py:33-43, signal.py:102-121 | `objdel.verify.stats` |
| Reward weights and function order | spec.py:28, signal.py:31, configs.py:191 | `REWARD_V4` |
| Sampling | sample.py:49, sample_vllm.py, hf_sample.py | `objdel.policy.sampler` with 2 adapters |
| Row reorder | fast.py:375-404 | one `reorder(batch, order)` |

### 3.6 Test gaps (critical paths first)

The suite is large (957 tests) and strong where it counts for exactness: the oracle K1/K2, TRL-equivalence of log-probs
and gradients for windowing and shared prefix, GPU-scoring equality. The gaps:

1. **No GPU build in CI** (tests.yml builds the CPU extension only). `raster_row_warp` (gpucheck.cpp:388-415), the
   atomics and `unite` races (1262-1282), and the two-lane split (1671-1693) run only in rented-box config steps
   (configs/rl_full.yaml, including compute-sanitizer). Add a self-hosted or spot GPU CI job that runs
   `tests/test_cuda_checker_gpu.py` and `tests/test_gpu_scoring.py` plus memcheck on each change to `cuda/`, or at least
   a required "box gate" artifact checked by preflight.
2. **No TRL contract test** (3.2). Nothing fails when a private name disappears. `set_is_ratio` and `skip_old_logps`
   are tested through monkeypatched trainers (test_train_fast.py:365-405), not through a real `GRPOTrainer.train()`
   step.
3. **The statistics of the headline numbers**: `final_eval.summarize` has no cluster test (3.1-1).
4. **Tests that grep the source** (16 places, for example test_train_fast.py:705, 768; test_p2_perf.py:18-57;
   test_p2_v1_select.py:220) assert literal source text such as
   `'trainer.shared_prefix = bool(getattr(a, "shared_prefix", False))' in src`. They pin the implementation, block
   refactoring, and prove nothing about behaviour. Replace them with behaviour tests on the factory once it is library
   code.
5. **The multi-root loader `serve._tagged`** (serve.py:120-122) has no direct test of a bad tag (`"x::id"` raises
   `ValueError` from `int()`, `"5::id"` raises `IndexError`). Only `grpo_rows_views` is covered (tests/test_p2_views.py).
6. **Table lifecycle under concurrency**: `drop_tables` from a finalizer while the prefetch thread scores; handle reuse
   (gpucheck.cpp:1752-1759). Only `resident()` counts are tested (tests/test_gpu_scoring.py:179-187).
7. (Checked, not a gap.) `tests/test_p2_train_rl.py:13-90` runs `train_rl` end to end on a tiny Qwen3-VL with the CPU
   checker (2 steps, resume, shared prefix, `--train-views`). What it does not cover is the vLLM path: `use_vllm` is
   off on CPU, so `patch_vllm_prompts`, `sync_without_merge` and C10 `skip_old_logps` run only on a GPU box.

### 3.7 Config sprawl (`experiments/phase2/configs/*.yaml`, 16 files, 378 lines)

- Each file is a **list of argv vectors** (`steps:`) plus a long comment of the manual shell commands. Hyperparameters
  (`--steps 300`, `--save-steps 50`, `--shared-prefix`, model name) sit inside CLI strings. There is no schema, no
  defaults and no inheritance. `rl_full.yaml` and `rl_shifted.yaml` differ by 20 lines (diff checked); `r5_shared.yaml`
  copies the box-probe step.
- Launch config and run config are separate concepts that the comments mix: every file says
  `launch.sh ... --config experiments/phase2/configs/r1_qwen3vl2b.yaml`, then `run.sh ... <this file>`.
- The real training hyperparameters live in three places: `configs.py` (code), `prereg.py` (experiment), and argv in
  YAML.
- Proposal: one typed `RunSpec` (a dataclass or pydantic model, saved with the run) for `{policy, data roots, rl:
  {steps, lr, shared_prefix, ...}, eval: [...]}`. Small per-run YAML files override it (`base: rl_full`). A tiny
  generic chain runner turns it into the step list. Keep the run-log comment as a `notes:` field. Retire c2_*, r0*, r1*,
  rx8, p1_rl, r3-r5 to `archive/` (completed runs whose results are recorded).

### 3.8 Doc drift

- rewards.py:1 and configs.py:1 say "TRL 1.16". The pin is 1.17.1, and fast.py cites 1.17 line numbers
  (fast.py:352, 672). Docstrings that quote line numbers of third-party files (`grpo_trainer.py 2776-2824`,
  `grpo_trainer.py 1844`) go stale on every upgrade. Put them in the contract table of 3.2 (a).
- Module docstrings are long design records (cache.py:1-75, metrics.py:1-95). They are valuable but mix interface and
  history. Move the history and evidence to `docs/research` and keep the interface facts in the code
  ("link, don't duplicate").

### 3.9 Debt priorities ((Impact + Risk) x (6 - Effort))

| Item | I | R | E | Score |
|---|---|---|---|---|
| Cluster bootstrap in final_eval (3.1-1) | 4 | 5 | 1 | 45 |
| Remove the unused pool in train_rl | 2 | 2 | 1 | 20 |
| TRL contract check + keyword super() calls (3.2 a-c) | 4 | 5 | 2 | 36 |
| Derive weights and order from REWARD_V4 / funcs (3.1-5) | 3 | 4 | 1 | 35 |
| One LR source (3.1-4) | 3 | 3 | 1 | 30 |
| Move make_trainer, score_items, view_rows, save_json to src (1.3-1) | 4 | 3 | 2 | 28 |
| Eval samples through a Sampler, not a trainer (1.3-2) | 4 | 3 | 3 | 21 |
| GPU CI job (3.6-1) | 3 | 4 | 3 | 21 |
| Skip the relabel when the filter is a no-op (A1) | 3 | 1 | 2 | 16 |
| Batch the target-ladder checker calls (A8) | 3 | 1 | 3 | 12 |
| Masked-attention kernel for shared prefix (A14) | 4 | 2 | 4 | 12 |
| Delete the numpy polygonizer from the library (A9) | 2 | 2 | 2 | 16 |
| Typed RunSpec replacing the argv YAML (3.7) | 3 | 2 | 3 | 15 |

---

## 4. Target package layout

Verdict (service-decomposition): a **modular monolith**. One wheel, enforced internal seams (import-linter
contracts), no network service yet. The verifier is the only part with a plausible second consumer (an external
reward or eval API). Give it a narrow in-process interface now, so that extracting it later is a packaging change.

```
src/objdel/
  geom/                 # exact geometry engine: the CUDA core + bindings (the deepest module)
    _engine (ext)       #   cuda/gpucheck.{h,cpp}, gpucheck_py.cpp  (unchanged interface)
    polygons.py         #   pack / parse-to-packed (geometry.text + shape), to_masks
    topology.py         #   CPU reference only (oracle), size rule, min_area
  verify/               # THE verifier: "completion(s) + item -> scores"; no torch, no TRL
    spec.py             #   REWARD_V4 (single source of weights, floors, order)
    experts.py          #   ExpertCache + CacheStore (replaces serve.STORE globals), device upload batching
    score.py            #   score_batch(jobs) -> Terms  (the one GPU call)
    advantages.py       #   advantages_v4, topology_reward
    metrics.py          #   eval metrics on device (A10), summarize
    stats.py            #   cluster bootstrap, Spearman, paired CIs (signal.py's statistics)
    service.py          #   Verifier(store).score(texts, keys) / .evaluate(samples) : the seam for an API later
  data/                 # store, LIDC/RIGA, splits, prepared roots (one manifest reader), views, augment,
    targets/            #   frame, simplify, build_target (batched ladder on geom)
    probes.py           #   synthetic probe *generation* (moved out of eval)
  policy/               # model + representation; no training logic
    spec.py             #   MODELS, processor kwargs, sampling(), prompt text, LoRA targets
    sampler.py          #   Sampler interface; adapters: VLLMSampler, HFSampler (two adapters = real seam)
    qwen_vl.py          #   VisionCache, RopeCache, ImageMemo, one_pad_per_image (model-specific adapters)
  rl/                   # algorithm + framework adapter
    build.py            #   build_trainer(PolicySpec, RLData, RunSpec)  (today's g_probe.make_trainer)
    trl_adapter.py      #   FastGRPO glue + the private-API contract table and import-time check
    forward.py          #   window(), pack_groups(), scatter_back(), is_ratio()  (pure, framework-free)
    reward_fn.py        #   TRLReward (adapter from verify.score to TRL's reward-func protocol)
    callbacks.py        #   Controller (wired into train_rl, or archived), TraceStep
  sft/                  # SFT run state machine (today's train_sft.train), sampler order, resume
  infra/                # cloud runners and ops (optional extra; not imported by the above)
    vast/               #   policy.py, bundle.py, watchdog
    runspec.py          #   typed RunSpec + chain runner (replaces argv YAML)
    secrets.py, hf.py   #   hf_push/pull
experiments/<phase>/    # thin CLIs: parse args -> RunSpec -> call the library; no cross-imports between scripts
```

Dependency rules (enforce with import-linter in `make lint`): `geom` <- `verify` <- `data` (targets) ;
`verify` never imports torch, transformers or trl ; `policy` <- `rl` and `sft` ; `rl` -> `verify` only through
`reward_fn` ; `infra` imports nothing from `rl` or `sft` (it runs CLIs) ; `experiments` import only `objdel.*`.

### Migration order (each step keeps `make test` green)

1. **Fixes with no moves** (one PR each): cluster bootstrap in final_eval plus a test; remove the train_rl pool;
   derive `signal.WEIGHTS` and `reward_weights`; one LR constant in `configs.py` that prereg reads; keyword `super()`
   calls; replace `out[3], out[1]` with a named helper. Add the TRL contract test now. It fails loudly on upgrade and
   guards the later moves.
2. **Promote the hidden library**: move `g_probe.make_trainer` -> `objdel/rl/build.py`, `save_json` -> `objdel/io.py`,
   `select_ckpt.score_items/view_rows` -> `objdel/verify/` and `objdel/data/`. Leave re-export shims in the old
   modules (`from objdel.rl.build import build_trainer as make_trainer`) so the 79 test imports keep working. Then
   rewrite the source-grepping tests as behaviour tests against `build_trainer`, and drop the shims.
3. **Make the `verify` seam**: create `objdel/verify/` by moving `reward/*` and the statistics from `eval/signal.py`.
   Keep `objdel.reward` and `objdel.eval` as re-export packages for one release. Replace `serve.STORE` with
   `CacheStore` (still module-level at fork time). `serve.preload` and `dataset.cache_loader` both call one
   `data.prepared.open_root`.
4. **Sampler seam**: add `policy/sampler.py` with `VLLMSampler` (from sample_vllm.py) and `HFSampler`. Move
   `final_eval` and the G-probe scoring onto it; `sample.sample(trainer, ...)` remains only inside RL. final_eval no
   longer needs `--train`.
5. **Split fast.py**: pure functions -> `rl/forward.py` (their tests already exist), Qwen caches ->
   `policy/qwen_vl.py`, trainer glue -> `rl/trl_adapter.py`. The equivalence tests in test_train_fast.py are the safety
   net; they import through `objdel.train.fast` re-exports until the end.
6. **Engine consolidation**: switch `targets/build` to `_gpucheck.check_masks`, batch the ladder (A8), and move the
   numpy polygonizer to `tests/oracles/`. Re-run the target-set hash check: it must be byte-identical (it is by the
   engine's bit-identity claim; preflight verifies it).
7. **Engine performance** (behind the existing GPU equality tests plus memcheck): A1 relabel skip, A3 one round trip,
   A5 batched table upload and a pooled free.
8. **Config**: introduce `RunSpec` + `base:` inheritance. Convert `rl_full`/`rl_shifted` first, archive finished run
   configs, and make `vast.policy.chain` read `RunSpec`.
9. **Infra extra**: move `experiments/phase2/vast` and `objdel.secrets` to `objdel.infra` (an `infra` optional extra).
   Archive `experiments/phase1/lab` and `archive/` to a tag.
10. **SFT**: extract `train_sft.train` into `objdel/sft/` as the last step (it has the most private overrides and the
   fewest reruns planned).

What to revisit later: if a second consumer of `verify.service` appears (an external grader, or a separate eval box),
expose it over gRPC with packed-array messages. The GPU-resident `CacheStore` makes a long-lived verifier process the
natural shape, but not before then (YAGNI).

---

### Checks behind the claims

- I traced the call sites with grep: `make_trainer` users, `check_masks` and `polygonize` callers, `Controller`
  constructors, and the use of `pool` in train_rl.
- I measured function lengths with an AST scan: train_sft.train 174, packed_chunk 76, metrics 75,
  `_get_per_token_logps_and_entropies` 61.
- Duplicate files: md5 over tracked files (24 identical pairs, all under `experiments/phase1/lab`).
- `rl_full.yaml` vs `rl_shifted.yaml`: diff checked.
- Not measured: the actual milliseconds of A1-A3 and A14 on a GPU (no GPU here). Their order of magnitude comes from
  the code structure only.
