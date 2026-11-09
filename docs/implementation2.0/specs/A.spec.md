# Phase A spec: measure fairly

Paths are relative to the repo root. Tests run on a Mac with no GPU unless marked [box]. The CUDA checker
(`objdel._gpucheck`) is available as its CPU build on the Mac and gives the same results as the GPU build.
Revision 2 (2026-11-08): matches `A.md` revision 2.

## 1. Goal

Make every evaluation number fair and repeatable:
- intervals that resample whole nodules and patients;
- a human ceiling where each doctor and the model are scored against the same k - 1 doctors;
- every answer saved, before any scoring;
- a label-free choice of one answer;
- the standard multi-rater metrics;
- NVTX time marks over each whole training step, and a per-step trace analysis;
- a start-up check of the private TRL names and the silent TRL couplings;
- a true count of outline token lengths (doctors on three train roots, and the model's answers).

## 2. Terms

- **item**: one image with its doctor outlines. Its id has the form `lidc/<series uid>/<nodule>/<slice>`, e.g.
  `lidc/1.3.6.../n000/z185`.
- **nodule cluster**: the item id without its last `/<slice>` part. A non-`lidc/` id is its own cluster
  (existing `objdel.eval.signal.cluster`).
- **patient cluster**: the record's `group_id` (e.g. `LIDC-IDRI-1001`).
- **K**: answers per item (8). **k**: doctors (raters) of the item (2 to 4 for LIDC).
- **legal**: the answer passes the checker's legality gate. Illegal and unparseable answers score 0 on every
  geometric term.
- **miss**: an item whose K-mean IoU_tc is < 0.1. **found**: not a miss.
- **step period**: the time from one `on_step_end` callback to the next.

## 3. Required behaviours and interfaces

### 3.1 Cluster bootstrap: `objdel.eval.signal`

`draws(ids: Sequence[str], boot: int = 2000, seed: int = 0, key: Callable[[str], str] = cluster) -> list[np.ndarray]`
- It returns `boot` index arrays. Each array concatenates the item indices of `C` clusters, drawn with replacement from
  the `C` distinct keys (in first-appearance order), with `np.random.default_rng(seed)`.
- With the default `key`, the output is bit-identical to the current implementation.
- `key` maps an item id to its cluster id.
- With one item per cluster, the arrays equal the rows of `np.random.default_rng(seed).integers(0, n, (boot, n))`.

### 3.2 Final evaluation: `experiments.phase2.final_eval`

`METRICS: dict[str, tuple[str, bool]]`: name -> (source, higher_is_better). It holds exactly these names:
`legal, iou_tc, f5, precision_vote, recall_vote, dice_vote, nsd, boundary_iou, boundary_iou_vote, soft_dice, hm_iou,
ged`. Only `ged` has higher_is_better = False.

`per_item(result: dict, obj: str = "nodule") -> dict`
- `result` is one `objdel.eval.metrics.metrics()` output with `item_id`.
- Returns `item_id`, every METRICS key and `miss`:
  - per-sample lists (`legal, iou_tc, f5` from `samples`; `precision_vote, recall_vote, dice_vote, nsd, boundary_iou,
    boundary_iou_vote` from `objects[obj]`) are averaged over the K answers;
  - item scalars (`soft_dice, hm_iou, ged`) are copied;
  - `miss = 1.0 if iou_tc < 0.1 else 0.0`.
- It also returns the selection keys of §3.5: `sel_mean_of_k, sel_legal_mean, sel_medoid, sel_medoid_legal,
  sel_oracle`, and `sel_vote` when available.

`summarize(items: list[dict], groups: Mapping[str, str] | None = None, draws: int = 2000, seed: int = 0) -> dict`
- Top level: `schema` = 2, `items`, `nodules`, `patients` (null without `groups`), `boot` = draws, `ci_level` = 0.95.
- Per METRICS key and per selection key: `{"mean", "ci95_nodule", "ci95_patient", "ci95_item", "higher_is_better"}`.
  There is no key named `ci95`.
  - `ci95_nodule`: the 2.5 % and 97.5 % quantiles of the item mean under nodule-cluster draws
    (`signal.draws(ids, draws, seed)`);
  - `ci95_patient`: the same under `key = groups[item_id]`; null without `groups`;
  - `ci95_item`: `np.percentile(v[idx].mean(1), [2.5, 97.5])` with
    `idx = np.random.default_rng(seed).integers(0, n, (draws, n))` (the previous behaviour, unchanged).
- `iou_tc` also has `median`, the median of the per-item IoU_tc.
- `miss_rate`: the mean of `miss`, with CIs. `iou_found`: the mean IoU_tc over found items, with a nodule CI on that
  subset.
- `paired`: `{"medoid_minus_mean": [d, lo, hi], "vote_minus_mean": [d, lo, hi] or null}`. These are the mean per-item
  differences with nodule-cluster CIs.

CLI (`python -I -m experiments.phase2.final_eval`):
- Sampling mode (existing arguments `--ckpt --view --split --train --out [--k 8] [--seed 11] [--per-call 64]`)
  [box]. Write order:
  1. after sampling and before any scoring: `<out>/answers.jsonl`, one line per item, `{"item_id": str,
     "group_id": str, "W": int, "H": int, "samples": [str] * K, "truncated": [bool] * K, "tokens": [int] * K}`; and
     `<out>/final_eval.json` with `"partial": true`, `ckpt, view, split, K, seed, view_manifest_sha256, code_sha,
     checker_build` and `seconds.load_and_sample`;
  2. after scoring: `<out>/metrics.jsonl` (one full `metrics()` result per item), `<out>/items.jsonl` (one
     `per_item` record per item), and `<out>/final_eval.json` replaced by the `summarize` output plus `ckpt, view,
     split, K, seed, seconds, view_manifest_sha256, code_sha, checker_build` and `"partial": false`.
  - If scoring raises, `answers.jsonl` and the partial `final_eval.json` remain on disk.
  - `view_manifest_sha256`: SHA-256 hex of `<view>/manifest.json`. `checker_build`: `"cuda"` or `"cpu"`.
  - `code_sha`: the stripped content of `<repo root>/CODE_SHA` when that file exists; else `git rev-parse HEAD` when
    the repo has `.git`; else null.
- Rescore mode: `--rescore DIR --view V --split S --out O`. It reads `DIR/answers.jsonl`, loads no model and writes the
  same four files to `O`. It accepts a partial or a full `DIR/final_eval.json`. It exits with code 2, writing no
  output, when the SHA-256 of `V/manifest.json` differs from `DIR/final_eval.json["view_manifest_sha256"]`.
- Rescoring the same answers gives per-sample `legal, iou_tc, f5, g_area` equal to the bit.

### 3.2a Code identity on the box: `experiments.phase2.vast.policy`

- `sync_up` (also with `--code-only`) writes `/workspace/src/CODE_SHA` on the box: the Mac's `git rev-parse HEAD`,
  followed by `-dirty-<12 hex chars>` when the working tree differs from HEAD (the first 12 hex characters of the
  SHA-256 of `git diff HEAD`).

### 3.3 Sampling: `objdel.train.sample.sample`

`sample(trainer, rows, k, items=16) -> list[dict]`: each dict has `item_id`, `samples` (K texts), `truncated`
(K bools) and `tokens`: K ints, each the length of that answer's completion id sequence as returned by
`trainer._generate`.

### 3.4 Metrics additions: `objdel.eval.metrics`

- `BIOU_D = {"nodule": 2, "disc": 9, "cup": 16}` (pixels).
- `band(m: np.ndarray, d: int) -> np.ndarray` (bool, same shape): `m` AND NOT `erode(m)`. `erode` pads `m` with one
  pixel of zeros on every side, applies `cv2.erode` with a 3x3 kernel of ones for `d` iterations, then removes the pad.
- `boundary_iou(a: np.ndarray, b: np.ndarray, d: int) -> float`: |band(a,d) AND band(b,d)| / |band(a,d) OR band(b,d)|.
  Two empty bands: 1.0. Exactly one empty band: 0.0.
- `metrics()` adds per object:
  - `boundary_iou`: list of K floats, the max over the item's raters of `boundary_iou(sample, rater, BIOU_D[name])`;
  - `boundary_iou_vote`: list of K floats, against the strict-majority vote mask;
  - `sample_iou`: K x K nested list, the IoU between answers' masks; illegal or failed answers count as empty masks;
    two empty masks have IoU 1.0.
- The `boundary_iou` values equal those computed on the full frame, also for masks that touch the frame edge.
- Illegal or failed answers get 0.0 for `boundary_iou` and `boundary_iou_vote`.
- Existing keys and values are unchanged.

### 3.5 Selection: `objdel.eval.select`

`medoid(sample_iou: np.ndarray, legal: Sequence[bool]) -> int`
- `sample_iou` is (K, K). Eligible = indices with `legal` true. If none, return 0.
- One eligible index: return it. Otherwise, for each eligible i, score_i = the mean of `sample_iou[i, j]` over
  eligible j != i. Return the argmax; on ties, the lowest index.
- It takes no rater input. Its result does not depend on any doctor mask.
- Legal empty answers have IoU 1.0 with each other, so a majority of legal empty answers wins.

`vote_mask(masks: np.ndarray, legal: Sequence[bool]) -> np.ndarray`
- `masks` is (K, H, W) bool. Returns an (H, W) bool mask: pixel true when 2 x (count of legal masks true there) >
  (number of legal masks). With no legal mask, all false.

Selection keys per item (in `per_item`; for several objects, `sample_iou` is averaged over the objects first):
- `sel_mean_of_k`: the mean IoU_tc over K;
- `sel_legal_mean`: the mean IoU_tc over legal answers (0.0 if none);
- `sel_medoid`: the IoU_tc of the medoid answer; `sel_medoid_legal`: 1.0 if that answer is legal, else 0.0;
- `sel_oracle`: the max IoU_tc over K;
- `sel_vote`: the IoU_tc of `vote_mask` scored against the item's cache.

### 3.6 Fair ceiling: `objdel.eval.ceiling` and `experiments.phase2.ceiling`

`LeaveOut = NamedTuple("LeaveOut", [("index", int), ("cache", ExpertCache), ("single", bool)])`

`loo_caches(cache: ExpertCache) -> list[LeaveOut]`
- There is one entry per voting rater r (`cache.voters`), in increasing r.
- k >= 3: `cache` = `ExpertCache.build` of the rater stack without r; `single` = False.
- k = 2: `cache` = `ExpertCache.build` of the stack `[other, other]` (the other rater twice); `single` = True. Scoring a
  mask against it gives IoU_tc equal to the plain IoU with the other rater.

`item(cache: ExpertCache, texts: Sequence[str] | None, truncated: Sequence[bool] | None, W: int, H: int,
name: str = "nodule") -> dict`. It returns:
- `k`, `single` (k == 2);
- `doctor`: the mean over the entries of `score_mask(cache.raters[r], loo.cache).iou_tc`;
- `legacy_dmax`: list, per voting r, of the max plain IoU of rater r with any other rater;
- `legacy_loo_dmax_item`: the mean of `legacy_dmax`;
- with `texts`:
  - `model`: the mean over entries and over K answers of the answer's IoU_tc against `loo.cache` (illegal or failed
    answers: 0.0);
  - `model_all_k`: the mean over K of IoU_tc against `cache`;
  - `diff` = `model - doctor`;
- without `texts`: `model`, `model_all_k` and `diff` are null.

`summarize(rows: Sequence[dict], groups: Mapping[str, str] | None = None, boot: int = 2000, seed: int = 0) -> dict`
- Strata `all`, `k>=3`, `k=2`. Each has `n` plus `doctor, model, model_all_k, diff, legacy_loo_dmax_item`, each given
  as `[mean, lo, hi]` under nodule-cluster draws, and `legacy_loo_dmax_pooled`: the mean of all `legacy_dmax` values
  of the stratum pooled together. `diff` is the CI of the per-item differences.
- `n(all) = n(k>=3) + n(k=2)`.
- When every item has the same k, `legacy_loo_dmax_pooled` equals the mean of `legacy_loo_dmax_item`.

CLI: `python -I -m experiments.phase2.ceiling --view V --split S [--answers FILE] --out O` writes `O/ceiling.json`
(summarize output plus `seconds`) and `O/ceiling_items.jsonl`.

Target text check: a target polygon text rasterized with `objdel.geometry.parse` + `to_masks` equals the rater mask it
was drawn from.

### 3.7 NVTX: `objdel.train.nvtx`, `objdel.train.callbacks.NvtxSteps`, FastGRPO counters

- Importing `objdel.train.nvtx` reads the environment variable `OBJDEL_NVTX` and calls no `torch.cuda` function.
- `enabled() -> bool`: true only when `OBJDEL_NVTX` equals `"1"` and `torch.cuda.is_available()`. The CUDA query
  runs at most once per process (the result is cached).
- `range(name: str)`: a context manager. When enabled, it calls `torch.cuda.nvtx.range_push(name)` on entry and
  `range_pop()` on exit (also on an exception). When disabled, it calls no `torch.cuda.nvtx` function.
- `push(name: str)`, `pop()`: `range_push` / `range_pop` when enabled, else nothing.
- `FastGRPO.timed(phase)` opens `range(phase)` and still appends the wall time to `time/<phase>_s`.
- `FastGRPO.training_step` runs inside range `micro`; `FastGRPO.get_batch_samples` runs inside range `loader`;
  `FastGRPO.log` runs inside range `log`. Their results are unchanged.
- Fixed range names: `step, loader, micro, batch, generate, sync, wake, sleep, rewards, checker, logps, forward,
  optimizer, log, vision_cache, rope_cache, image_memo`.
- `NvtxSteps(window: tuple[int, int] | None)`:
  - range `step` covers each step period: pushed at `on_train_begin`, popped and pushed again at each `on_step_end`,
    popped at `on_train_end`. Pushes and pops balance;
  - range `optimizer` from `on_pre_optimizer_step` to `on_optimizer_step`;
  - with `window = (s, e)` and CUDA available: `torch.cuda.profiler.start()` once, in `on_step_end` when
    `state.global_step == s - 1`, after the pop and before the next push; `torch.cuda.profiler.stop()` once, in
    `on_step_end` when `state.global_step == e`, after the pop;
  - with no CUDA: no profiler call.
- `g_probe.make_trainer` adds `NvtxSteps(tuple(a.profile_window))` when `a.profile_window` is set, and only then.
  `TraceStep` is added when `a.trace` is set, with `wait = a.profile_window[0] - 2` when a window is set, else 1.
- `train_rl --profile-window S E` and `train_rl --trace DIR` pass these to `make_trainer`. Without the flags, the
  behaviour is unchanged.
- `fast.raw_ratio(logps, sampling, mask, mode)` and `fast.bounded(r, mode, lo, hi)`: `bounded(raw_ratio(...))` equals
  `fast.is_ratio(...)` in all four modes (`sequence_mask, sequence_truncate, token_mask, token_truncate`).
- The step log gains, over the counted ratio entries (sequence rows in sequence modes, masked tokens otherwise):
  - `sampling/importance_sampling_ratio/zero_share`: the share whose final ratio equals 0;
  - `sampling/importance_sampling_ratio/hi_share`: the share whose raw ratio is above `clip_max`;
  - `sampling/importance_sampling_ratio/lo_share`: the share whose raw ratio is below `clip_min` or equals 0 (each
    entry counted once).
  These add no host read.
- The step log gains `time/step_s`: the wall time between consecutive `log` calls. It adds no device sync.

### 3.7a Trace analysis: `experiments.phase2.nsys_baseline`, `experiments.phase2.profile`

- `profile.covered(intervals: Iterable[tuple[float, float]], t0: float, t1: float) -> float`: the length of the union
  of the intervals clipped to [t0, t1]. `profile.busy` gives the same output as before.
- CLI: `python -I -m experiments.phase2.nsys_baseline --stats DIR --log RL_JSON --out FILE`. It reads
  `DIR/stats_nvtx_pushpop_trace.csv`, `stats_cuda_gpu_trace.csv`, `stats_cuda_api_trace.csv` and
  `stats_nvtx_gpu_proj_trace.csv`, and writes FILE (JSON) with one entry per `step` range instance:
  - `coverage`: the summed duration of the ranges directly inside `step` (nesting by time containment on the same
    thread) / the `step` duration;
  - `gpu_busy`: `covered` of the kernel, memcpy and memset intervals inside the step / the step duration;
    `gpu_idle` = 1 - `gpu_busy`;
  - `host_syncs`: counts by name of `cudaStreamSynchronize`, `cudaDeviceSynchronize`, `cudaEventSynchronize` and
    non-`Async` `cudaMemcpy*` calls inside the step, with any `_v<digits>` name suffix removed;
  - `per_range_s` (total duration by range name) and `gpu_proj_s` (projected GPU time by range name);
- and top-level `overhead`: `time/step_s` of the window steps and of the steps just before and after, and
  `overhead_affected` (true when the window steps are more than 10 % slower).

### 3.8 TRL check: `objdel.train.trl_adapter`

- `PIN = "1.17.1"`. `class TRLContractError(RuntimeError)`.
- `check() -> None`. It checks, in this order, and raises `TRLContractError` with the failing name in the message:
  1. `GRPOTrainer` has these methods, with exactly these parameter names in order:
     - `_generate_and_score_completions(self, inputs)`
     - `_prepare_inputs(self, generation_batch)`
     - `_compute_loss(self, model, inputs)`
     - `_get_per_token_logps_and_entropies(self, model, input_ids, attention_mask, logits_to_keep, batch_size,
       compute_entropy, compute_aux_loss, pixel_values, image_grid_thw, num_images, pixel_attention_mask,
       spatial_shapes, num_tiles, image_sizes, token_type_ids, mm_token_type_ids, image_position_ids, sampling_mask)`
     - `_generate(self, prompts)`
     - `_calculate_rewards(self, inputs, prompts, completions, completion_ids_list)`
  2. `trl.trainer.grpo_trainer.shuffle_sequence_dict` exists with one parameter `seq_dict`, and
     `trl.generation.vllm_generation.VLLMGeneration._iter_named_params` exists with parameters `(self)`;
  3. the source of `GRPOTrainer.__init__` contains `self._metrics`, and the source of `VLLMGeneration` contains
     `self._dist`;
  4. the source of `GRPOTrainer` contains `self.<name> =` for each of `use_vllm, loss_type, num_iterations,
     num_generations, top_entropy_quantile, vllm_importance_sampling_correction, vllm_importance_sampling_mode,
     vllm_importance_sampling_clip_min, vllm_importance_sampling_clip_max, vllm_generation`; the source of
     `VLLMGeneration` contains `self.llm`, `self.model` and `def sync_weights`;
  5. the source of `GRPOTrainer._generate_and_score_completions` contains the quoted keys `"sampling_per_token_logps"`,
     `"completion_mask"`, `"tool_mask"`, `"importance_sampling_ratio"`, and the substring
     `self.use_vllm and self.vllm_importance_sampling_correction`; the source of `GRPOTrainer._compute_loss` contains
     `"importance_sampling_ratio"`;
  6. `trl.__version__ == PIN`. On a mismatch with items 1-5 passing: it raises, unless the environment variable
     `OBJDEL_TRL_UNPINNED` equals `"1"`; then it emits one warning naming the found version and returns.
- It runs with no GPU and no vLLM installed or running.
- TRL's own importance-ratio block (cut from the source of `_generate_and_score_completions`, from the line above
  `per_token_logps_diff = (old_per_token_logps - sampling_per_token_logps) * mask` to the `else:` before
  `Unknown vLLM importance sampling level`) gives a ratio equal to `fast.is_ratio` in all four modes.
- `FastGRPO.__init__` calls `check()` before any other work.
- `experiments.phase2.preflight.CHECKS` has an entry named `"CODE TRL private names"`, called as `fn(a, cfg)`, that
  returns `[]` when `check()` passes and a one-element list with the error message when it raises. It is the first
  `CODE` entry.

### 3.9 Outline lengths: `experiments.phase2.outline_lengths`

`count(records: Iterable[dict], tokens_of: Callable[[str], int], voting: Callable[[dict], set[int]] | None = None)
-> dict`
- For each record and each non-None entry of `record["targets"]`, `tokens_of(text)` is the length.
- Returns `n_records, n_targets, n_voting_targets, mean, p50, p90, p99, p99_9, max`, and
  `above = {L: count of lengths > L}` for L in (256, 512, 800, 1000, 1200, 1536). Also `above_voting` (the same over
  voting targets only), and `longest`: up to 20 `{item_id, rater, tokens}` in decreasing tokens.

`answers(answer_rows: Iterable[dict], metric_rows: Iterable[dict]) -> dict`
- Joins by `item_id`. For each L in the same list: `share_above` (share of answers with `tokens` > L), and over those
  answers `legal_share`, `mean_iou_tc`, and `n_legal_good` (legal and IoU_tc >= 0.5).

- The CLI's `tokens_of(text)` = `len(tokenizer(text, add_special_tokens=False)["input_ids"]) + 1`. The tokenizer is
  `AutoProcessor.from_pretrained(configs.MODELS[model], **configs.processor_kwargs(model, dataset)).tokenizer`.
- `objdel.train.dataset.target_tokens` and `objdel.train.dataset.voting` replace `_tokens` and `_voting` (no aliases).
- CLI: `--train ROOT [--views ROOT ...] --model NAME --datasets ... --splits ... [--answers FILE ...] --out DIR` writes
  `DIR/outline_lengths.json`, keyed by root, dataset and split; each view root has its `manifest_sha256`; with
  `--answers`, a `model_answers` section per file.

## 4. Constraints

- Library code in `src/objdel`, with tests in `tests/`. Experiment scripts in `experiments/phase2/`, run as modules
  with `python -I`.
- `make test` and `make lint` pass. No new runtime dependency (numpy, scipy, OpenCV 5.0.0.93, torch 2.14.1,
  trl 1.17.1, transformers 5.20.1 as pinned).
- Mac tests need no GPU, no network, no vLLM engine and no file under `runs/` or `data/`. Fixtures live in
  `tests/goldens/` (e.g. `r5_test_items.jsonl.gz`: `item_id` and `iou_tc` of the 2,200 R5 test items). Tests that
  need a GPU are marked and skipped on the Mac.
- Pixels never leave the device on the GPU build for the checker terms (IoU_tc, F5, g_area, legality, topology).
  Only scores return.
- `final_eval.json` adds `schema: 2`; the CI keys are `ci95_nodule`, `ci95_patient`, `ci95_item`. Other existing keys
  keep their names.
- Selection never reads doctor masks. The selection method is chosen on val only.
- Data under `data/` is untrusted: read it with `python -I`.
- One GPU rental for the whole phase, North America host, after a passing preflight. Spend <= $1.00. The nsys step
  is the last plan step and stops after at most 900 s plus 120 s.

## 5. Acceptance gates

1. On `tests/goldens/r5_test_items.jsonl.gz` (2,200 items, 401 nodule clusters), IoU_tc:
   - `ci95_item` = [0.6334010210209712, 0.6561989221793691] exactly;
   - `ci95_nodule` = [0.6279017988891661, 0.6620867640136667] exactly (seed 0, 2,000 draws);
   - `miss_rate` = 0.10818181818181818.
   With one item per cluster, `ci95_nodule` equals `ci95_item` exactly.
2. With the default `key`, `signal.draws` output equals the current implementation (existing B2 tests pass).
3. Ceiling: a doctor equal to another doctor scores 1.0. A model answer equal to doctor r gives diff 0 for that r.
   For k = 2, IoU_tc equals the numpy IoU of the two masks. With all same signatures (k = 4), the doctor IoU_tc equals
   the legacy max IoU. An illegal answer scores 0. Strata counts add up. A target text rasterizes to its rater mask.
   On val, n(all) = 1,084, n(k=2) = 281. [box data] `legacy_loo_dmax_item` or `legacy_loo_dmax_pooled` is within
   ±0.005 of 0.723 (val) and 0.732 (test), both on all items; the found-stratum item value is within
   ±0.005 of 0.734 (val) and 0.740 (test); otherwise the gate fails (final round 8, data-1: 0.734 was the found value).
4. Selection: the hand-made medoid cases pass (tie, no legal sample, a majority of legal empty answers). Permuting
   rater masks never changes `medoid`.
5. Boundary IoU: a 5x5 square with d = 1 has a 16-pixel band. Identical masks give 1.0; disjoint masks give 0.0; two
   empty masks give 1.0. Values equal the full-frame reference formula (`cv2.copyMakeBorder` + `cv2.erode`) exactly on
   random masks, including masks that touch the frame edge.
6. Rescore: on a fixture, `--rescore` reproduces the direct scoring; a wrong manifest SHA exits with 2. If scoring
   raises, `answers.jsonl` and a partial `final_eval.json` exist. `code_sha` comes from `CODE_SHA` when present.
   [box] On the box run's answers, the Mac rescore matches per-sample `legal, iou_tc, f5, g_area` to the bit.
7. [box] `answers.jsonl` has 1,084 val and 2,200 test lines, with K = 8 entries in each list. Paired by `item_id` with
   `runs/r5/eval/eval_test/items.jsonl`, the nodule CI of the per-item IoU_tc difference contains 0 and the mean
   difference is within ±0.01.
8. Gate A [box data]: the medoid (or vote) becomes the default only if the val mean of `sel_medoid - sel_mean_of_k`
   (or `sel_vote - sel_mean_of_k`) is >= 0.015 and its nodule CI lower end is > 0. The test paired difference of the
   chosen method has the same sign as on val.
9. NVTX: importing `nvtx` with `OBJDEL_NVTX=1` calls no `torch.cuda` function; disabled means no NVTX call. `timed`
   still logs `time/<phase>_s`. `NvtxSteps` balances pushes and pops and starts and stops the profiler at the right
   steps. `make_trainer` adds `NvtxSteps` only with a window. The zero, hi and lo counts equal numpy counts on a toy
   tensor. `nsys_baseline` gives the hand-computed values on the CSV fixture. [box] `nsys_baseline.json` has the
   `step` instances of steps 4 and 5, each with `coverage` >= 0.90, and records `gpu_idle`, `host_syncs` and
   `overhead`; on the fallback path, `profile.busy` numbers for step 4.
10. TRL: `check()` passes on trl 1.17.1. It raises, naming the item, for a changed version, a missing
    `_prepare_inputs`, a changed `_compute_loss` signature, a missing `shuffle_sequence_dict` and a missing
    `"tool_mask"` key. With `OBJDEL_TRL_UNPINNED=1` and a changed version only, it passes with a warning. TRL's ratio
    block equals `fast.is_ratio` in all four modes. The preflight stamp shows `"CODE TRL private names"` passing.
11. Outline lengths: with a one-token-per-character stub, the counts match hand values. None targets are skipped and
    the +1 end token is applied. `answers` matches hand values on a fixture. On `data/prepared` LIDC train:
    7,341 records, 23,469 targets, max 1,182, 0 targets above 1,200. `outline_lengths.json` also has the
    `data/views/lidc-train-shifted` and `data/views/lidc-train-rescaled` rows with their manifest SHAs.

## Review log (final round 8)

Cross-file consistency fixer, round 8. Sources: `reviews/final_r8_{cuda,data,arch,qa}.md` and their `*_verdict.md`
files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-1 | blocking | Item 3: val reference 0.723 on all items; found-stratum line 0.734 (val) / 0.740 (test) | `notes/r8data/legacy.py`, `legacy_found.py` |
