# Phase 1 code-quality review (post simplify sweep S1–S4)

Independent, read-only review. In scope: `src/objdel/**`, `src/objdel/_polycore.cpp`, `cpp/polycore/polycore.hpp`,
`experiments/phase1/*.py`, `experiments/phase1/bench/*` and `tests/**`. `experiments/phase1/lab` is out of scope.
I edited no files except this one, ran no git commands that change state, and did not run the full suite or Docker.
Repros are small `.venv/bin/python -I` snippets. The only one that read `data/` was read-only.

There are 12 findings, ranked by severity: 1 high, 2 medium, 9 low. The data → targets → prepared → reward → TRL
→ controller → eval wiring held up under every probe I tried (see "Checked and not a defect" at the end). The one
high finding is a security-hygiene problem in a fetch script. The two medium findings are about what the gates and
logs report.

## CONFIRMED (reproduced, or certain from the code)

### 1. HIGH: `fetch_lidc_pixels` deletes and writes paths taken unchecked from the untrusted XML

- **Where:** `experiments/phase1/fetch_lidc_pixels.py:59-60` (`tmp = tmp_root / series; shutil.rmtree(tmp, ignore_errors=True)`),
  `:65` (`API + series`), `:88`, `:94`. The source is `src/objdel/data/lidc.py:252-260` (`sop_z`, whose keys are the
  `SeriesInstanceUid` text of every XML file).
- **Defect:** the series UID is used without validation as a path component and inside the URL. `rmtree` runs
  before any network access.
- **Evidence:** I put a planted XML in the scratchpad with `<SeriesInstanceUid>../../../victim_dir</SeriesInstanceUid>`.
  `sop_z` returned that key, and the `rmtree` target resolved to
  `/Users/sreeyashwanthsaivenkatesh/Workspace/projects/CUDA/victim_dir`. An absolute UID is worse:
  `Path('data/lidc_tmp') / '/Users/someone/work'` is `/Users/someone/work`. `main` adds every such series to `todo`
  automatically (`:123`), so one bad XML under `data/raw/lidc_xml` deletes an arbitrary directory. The same value
  can also inject URL query parameters. `build_data` is safe only because it asserts `set(files) == set(scans)`
  against pylidc's database.
- **Minimal fix:** in `main`, keep only series with `re.fullmatch(r"[0-9]+(\.[0-9]+)+", s)` and at most 64
  characters. A stronger alternative is to keep only `s in lidc.scans()`. Build the URL with `urllib.parse.quote(s)`.
  Put the same check in `lidc._series` so every consumer gets it.

### 2. MEDIUM: illegal answers count toward the Phase 2 geometry gates

- **Where:** `src/objdel/eval/metrics.py:202-203` (scores are zeroed only for `~ok`, i.e. unparsed answers) and
  `summarize` at `:256` (`iou_tc_pass`) and `:264-266` (`iou_tc_median`, `*_mean`, `f5_median`). The gates that
  read them are `geometry_median` and `geometry_pass` in `eval/gates.py`.
- **Defect:** a parsed but illegal rollout keeps its IoU_tc, F5 and g_area. The reward does not: legality is a hard
  gate in `advantages_v4`. Also, the metrics docstring (`:28-30`) says the geometry gate is the reward's pass
  criterion, but the gate actually reads `iou_tc_pass`, which ignores both legality and topology.
- **Evidence:** I scored 8 copies of a square with a collinear spike (`spike` defect) against two near-identical
  raters. Result: legal = 0.0, pass_rate = 0.0, yet `iou_tc_pass` = 1.0 and `iou_tc_median` = 1.0. With
  legal = 0.90 (the lowest the LEGAL gate accepts), all of the 10 % illegal rollouts can still count as geometric
  successes.
- **Minimal fix:** zero the `_ZEROED` scores where `~legal` (legal implies parsed), matching the reward's legality
  gate, and update docstring line 34. The audit is unaffected: it reads areas, not `_ZEROED`. If the current
  behaviour is intended, say so in the gate's REQUIREMENTS text.

### 3. MEDIUM: the controller's audit statistics never reach any log

- **Where:** `src/objdel/train/callbacks.py:152-174` (`on_evaluate` → `metrics.update(self._audit(...))`).
- **Defect:** transformers logs the evaluation metrics before it calls `on_evaluate`, so the `eval_audit/...`
  statistics added there never reach `state.log_history` or `report_to`. The `audit` event in `controller.jsonl`
  stores only `fired`, `new` and `completions`.
- **Evidence:** `.venv/.../transformers/trainer.py:2850` is `self.log(output.metrics)`, and `:2855` is
  `callback_handler.on_evaluate(...)`. During training the return value of `evaluate()` is used only for
  `metric_for_best_model`. So the plan §9 audit curves (area medians, envelope shares, NRFR/RHR) are computed and
  then discarded. The docstring's "the audit statistics enter the metrics dict" is literally true, but nothing
  persists them.
- **Minimal fix:** `self._event("audit", fired=..., new=new, completions=len(records), stats=_flat(audit))`. That
  is one argument and makes the jsonl the record. Alternatively call `self.trainer.log(_flat(audit))`.

### 4. LOW: the build report undercounts RIGA test images by one

- **Where:** `experiments/phase1/build_data.py:125`
  (`stats["riga_with_image"][...] += kind == "prime"`) versus `src/objdel/data/riga.py:333`
  (`image=None if kind == "maxrg6" else ref`).
- **Defect:** the F7 item (`kind == "maxrg6:prime_is_rater1_copy"`) keeps its image in the store but is counted as
  having none.
- **Evidence:** `data/store/build_report.json` says `"with_image": {"test": 94, ...}`. The store has 95 test rows
  with an image; the odd one is `riga/Magrabia/MagrabiFemale/image048`.
- **Minimal fix:** `+= kind != "maxrg6"`, or count `it.image is not None` on the yielded item.

### 5. LOW: stale numbers and descriptions in docstrings

- `src/objdel/geometry/text.py:20-21`: "the longest of the 31,922 SFT targets has 2,864". The current
  `data/targets` holds 30,843 targets (`data/prepared` 30,543). The maximum, 2,864 characters, is still right.
- `src/objdel/train/configs.py:15-18`: "Issue I5 (10 / 3,546 RIGA train answers above 1,536 Qwen tokens)".
  EXECUTION §3 gives the post-B3-1 figure, 1 / 3,282 prepared RIGA train targets.
- `experiments/phase1/bench/xarch_hashes.py:54-65`: `CAUSES` names cv2 `INTER_CUBIC` and the IPP
  `ippiTrueDistanceTransform` as the cause of a mismatch. Both were removed by F3 and S2/F1, so a future mismatch
  message would point at code that no longer exists. Reword these as "historical (fixed by F3 / S2)" or drop them.
- **Fix:** update the numbers, or cite the manifest instead of hard-coding counts.

### 6. LOW: citation hygiene

- **Missing `References:` list** (EXECUTION §1.6 requires one in every module docstring):
  - `src/objdel/checker/bitquad.py`: Gray 1971 is named but not linked.
  - `src/objdel/checker/certify.py`
  - `src/objdel/checker/polygonize.py`: Wyllie list ranking and pointer jumping are used but not cited.
  - `src/objdel/reward/spec.py`
  - `src/objdel/eval/gates.py`: should point at plan §7.
  - All 8 phase 1 scripts: `build_data`, `build_targets`, `census_sizerule`, `prepare_rl`, `smoke_grpo`,
    `eval_run`, `eval_demo`, `reward_gates`.
- **Spot checks with WebFetch:** arXiv 2601.05242 (GDPO), 2607.09492 (Yao et al.; RHR and NRFR defined, NRFR > RHR
  diagnostic), 2607.21273 (Dark Room), 2605.02435 (Mohri, Schneider, Zhong; Jensen bias of log-reward estimators)
  and 2504.07615 (VLM-R1) all exist and match their descriptions.
- **Two claims the abstracts do not support:**
  - 2609.31644 (MaD-RL) is "distribution matching with KL / JS divergences". The cited "gated reverse KL … Eq. 10,
    17" could not be confirmed from the abstract; check against the PDF.
  - VLM-R1's `min(1, L_gt / L_pred)` count guard (`reward/score.py:29`) is not in the abstract; cite the section.

### 7. LOW: two copies of the truncation rule

- **Where:** `src/objdel/train/rewards.py:108` (`not ids or ids[-1] not in self.stop`) and
  `src/objdel/train/sample.py:44,58`. `sample` rebuilds its stop set from `trainer.eos_token_ids` and the pad token,
  while `rewards.stop_ids(model, processing_class)` builds it from `generation_config`. Today they agree, and TRL's
  rule is `grpo_trainer.py:2301-2302`.
- **Minimal fix:** in `sample`, use `rewards.stop_ids(trainer.model, trainer.processing_class)` and share one
  `truncated(ids, stop)` helper.

### 8. LOW: other duplications

- `build_targets.TOKENIZERS` (`:29-33`) is identical to `train.configs.MODELS`. Import it.
- `smoke_grpo.LORA` (`:41`) is identical to `configs.lora("qwen3-vl-2b")`.
- `build_data.py:151-152` (`cross`, `kcross`) re-implements `splits.crossing`, which is now used only by tests.
  Either use it there or delete it.
- Two union-find implementations: `data/lidc.py:172` (`clusters.find`) and `data/splits.py:162`
  (`group_ids.find`).
- The foreground bounding-box crop is written five ways:
  - `checker/certify.crop`
  - `checker/topology._crop_pad`
  - `reward/cache._crop`
  - `eval/metrics._pairwise_iou` (rows/cols)
  - the ROI in `reward/score.score_mask` and `ExpertCache.box`

  Each is small, so this is low priority. A shared `bbox(mask) -> (y0, y1, x0, x1)` in `checker.topology` would
  remove about 20 lines.

### 9. LOW: polycore code the production path never uses

- **Where:** `polycore.score_batch`: the C++ thread pool and exception plumbing (`polycore.hpp:650-711`), the
  binding (`_polycore.cpp:134-180`) and the Python wrapper and `_mask` (`checker/polycore.py:73-79, 111-123`). It is
  called only from tests and `bench/xarch_hashes.py`.
- **Also:** the reward path (`geometry.to_mask`) fills the nonzero mask and the winding statistics on every call,
  then reads only `evenodd`.
- **Why it stays:** the names are §3 contract names, and the S-sweep rule forbids removing them while tests use
  them.
- **Suggestion:** decide at Phase 2 / CUDA time, and pass `nonzero=nullptr` from a `to_mask`-specific entry point if
  profiling shows it matters.

### 10. LOW: the smoke run's GRPO-signal numbers use different definitions from the gates

- **Where:** `experiments/phase1/smoke_grpo.py:49-57` (`group_stats`).
- **Defect:**
  - Duplicates are counted on raw texts. The eval and gate definition (B2-6, `eval/metrics.py:159-161`) counts
    equal rasters.
  - σ is taken over all four terms without the legality gate. `metrics` takes it over the three geometric terms
    multiplied by legality.
- **Consequence:** `smoke.json` reports a "Phase 2 GRPO-signal gate" number that differs from what
  `gates.grpo_sigma` and `grpo_duplicates` would compute on the same rollouts.
- **Minimal fix:** compute the groups with `objdel.eval.metrics` on the batch, or label the smoke numbers as text
  duplicates and ungated σ.

### 11. LOW: `eval_demo` builds the rater targets with a different size rule from `prepare`

- **Where:** `experiments/phase1/eval_demo.py:50` (`ma = {n: c.min_area ...}`) has no `"rim"` key, so
  `build_target` certifies the RIGA rim at every possible rule. `prepare` uses `targets.min_areas`, which includes
  the rim.
- **Consequence:** the "raters" model is documented as "the rater's own SFT target", but it is built under a
  different rule.
- **Evidence:** I compared 48 RIGA val targets built both ways and none differ, so the defect is latent.
- **Minimal fix:** `ma = min_areas(masks)`.

### 12. LOW: `store.load` does not check that a row belongs to the file it came from

- **Where:** `src/objdel/data/store.py:189-195`.
- **Defect:** a row whose `dataset` or `split` field differs from the requested `<dataset>/<split>.parquet` is
  still yielded. `_check` only checks that the values are valid names.
- **Why it is low:** the manifest SHA-256 makes accidental mismatch unlikely. But the store is untrusted data, and
  split leakage is the K11 property.
- **Minimal fix:** `if (it.dataset, it.split) != (dataset, split): raise ValueError(...)`.

## PLAUSIBLE (not reproduced; low)

- **`prepare` does not check the image write.** `dataset.py:102`: if `cv2.imwrite` returns False (disk full, bad
  path), the record still carries `image_sha256`. `_check` catches it later, at row-building time. `prepare` also
  never clears `out_root`, and `cache_loader` (`:227`) globs every `<ds>/*.jsonl`, so a leftover `test.jsonl` from
  an earlier run is loaded even though the manifest does not list it. Fix: check the return value, and build
  `cache_loader` from `manifest["files"]`.
- **`fetch_lidc_pixels` overwrites quietly on shared z keys.** `:81-85`: two kept DICOMs that map to the same z key
  overwrite each other in `slices[f"z{z:.2f}"]` with no warning. HU values go through float32 and then
  `astype(int16)`, which truncates (it does not round) when RescaleSlope is not an integer, and wraps outside the
  int16 range.
- **`eval_run` crashes on single-rater items.** `metrics` → `ExpertCache.build` raises ValueError for items with
  one rater, so one such row in the completions file kills the whole evaluation. `sample_eval` emits only RL items
  today, and the current val split has no k = 1 items. Fix: skip the item with a warning, or document it.

## Checked and not a defect

- Non-canonical bool bytes (255): `certify`, `polygonize` and `ExpertCache.build` / `score_mask` give the same
  terms and digest as canonical masks.
- `signature` of a 0×0 ROI (empty prediction with an empty union) returns `()`.
- TRL 1.16 attributes used by `Controller` and `sample` all exist under those names: `num_generations_eval`,
  `eos_token_ids`, `_tokenizer`, `_generate`'s return order, `reward_func_names`.
- The reward's truncation rule equals TRL's (`grpo_trainer.py:2302`).
- `log_extra` is called once per batch (identity cache), and `_pending_extra_logs` matches the completions table.
- Determinism: `instance_topology` and `size_filter` do not depend on label numbering. The JSON outputs sort their
  keys. `sig_id` is CRC32 of a repr. `np.load` calls use `allow_pickle=False` (or numpy's default False). The only
  pickle in scope (`reward_gates` K10) loads a file the script wrote itself, in a fresh `-I` interpreter.

## Summary

The pipeline is in good shape. No correctness bug turned up in the core geometry, reward, advantage or TRL wiring;
the break-it and fix waves closed the earlier ones.

- **Fix before Phase 2:**
  - Validate the LIDC series UID in `fetch_lidc_pixels` (high; directory deletion from a planted XML).
  - Gate the geometry metrics on legality, or document why not (medium; an all-illegal model scores 1.0 on
    `iou_tc_pass`).
  - Persist the controller's audit statistics (medium; they never reach a log).
- **Cleanups:** everything else is low: a report miscount, stale numbers, missing References lists, duplicate
  helpers, and contract-protected polycore code that production does not use.
