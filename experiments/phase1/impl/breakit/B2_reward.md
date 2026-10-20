# Break-it B2: reward v4, TRL wiring, evaluation

Scope: `src/objdel/reward/`, `src/objdel/train/{rewards,callbacks,dataset.cache_loader,configs}.py`, `src/objdel/eval/`,
and how they use geometry and topology. Reproducers: `tests/test_breakit_reward.py` (strict xfails `B2-<n>`, plus
22 passing robustness checks; no data, ~9 s). Scratch: `scratchpad/breakit/B2/`. No production code was changed.

**Totals: 39 attacks tried, 29 held, 10 broke** (1 high, 4 medium, 5 low). Known issues I13, I15, I16 were not
re-measured. B2-10 makes I1 worse than its log entry says.

## Findings

| id | severity | module | one line |
|---|---|---|---|
| B2-1 | **high** | train.callbacks | An alarm on an evaluation that also improved the metric "rolls back" to that same hacked checkpoint |
| B2-2 | medium | train.callbacks | An alarm (or a hand-written rollback) before any best checkpoint raises and stops training |
| B2-4 | medium | reward (cache / topology histogram) | A rater whose whole mask is below the size filter becomes a "nothing" mode; a sub-filter speck earns IoU_tc 1, F5 1 and r_topo |
| B2-5 | medium | train.configs / rewards | `num_generations ≠ 8` is accepted; TRLReward silently splits each TRL group into v4 groups of 8 |
| B2-6 | medium | eval (metrics / audit / gates) | The duplicate alarm and the Phase 2 duplicate gate hash raw text; 8 spellings of one mask count as 0 % duplicates |
| B2-3 | low | train.callbacks | A malformed or half-written `control.json` raises JSONDecodeError and stops training |
| B2-7 | low | reward.advantages | The r_topo validity gate reads the mean IoU_tc; a perfect disc plus a cup outside it is "valid" (eval says fail) |
| B2-8 | low | train.dataset.cache_loader | After a fingerprint change the loader reports tampering; `ExpertCache.load`'s rebuild path is never reached |
| B2-9 | low | reward.score | Rater counts across an item's caches are not checked: IndexError one way, raters silently dropped the other way |
| B2-10 | low (makes I1 worse) | reward.cache | `ExpertCache.build` keeps bool arrays whose True bytes are 255, so rater holes are lost, and the digest matches the correct cache's |

### B2-1 (high): rollback restores the hacked checkpoint
- **Reproduction:** `test_b2_1_alarm_rolls_back_to_the_checkpoint_before_the_hacked_evaluation`. Evaluation 1 has
  iou_tc 0.80 and no alarm, so it becomes best. The parameters then change. Evaluation 2 has iou_tc 0.90 and the alarm
  fires. The rollback goes to step 50, the hacked step itself. Events: `best 0.8`, `best 0.9`,
  `rollback to_step=50`.
- **Expected vs actual:** expected the parameters of the last checkpoint that did not alarm (0.80). Actual: nothing
  changes, `best` is now the hacked one, and every later honest checkpoint scoring below 0.90 can never become best.
- **Root cause:** `Controller.on_evaluate` saves and updates `best` before it calls `alarm(metrics)`.
- **Why it matters:** this is the plan §9 defence ("alarms fire and the best checkpoint is restored"). The hedges the
  alarms exist for score a **higher** IoU_tc than honest output. M8 demo, LIDC: union 0.946, core 0.936, hull 0.878
  vs honest 0.870. RIGA disc: 0.963 / 0.956 / 0.958 vs 0.928. So "metric improved and alarm fired" is the expected
  case, not a corner case.
- **Fix:** evaluate the alarm first. Update `best` only when no alarm fires. Optionally keep the last few
  non-alarmed checkpoints, as checkpoint-selection protocols do (Yao et al. 2026 select on a held-out set,
  https://arxiv.org/abs/2607.09492).

### B2-2 (medium): alarm without a best checkpoint stops the run
- **Reproduction:** `test_b2_2_…`. The metric is missing from the evaluation metrics (renamed, or `eval_strategy`
  logs under a different prefix) and the alarm fires: `ValueError: no best checkpoint to roll back to`, raised
  inside `on_evaluate`. The same happens when `control.json` asks for a rollback before the first evaluation.
- **Expected:** log the alarm, and either skip the rollback or raise beta / stop cleanly. Not an uncaught exception
  that kills a multi-hour run.
- **Fix:** in `on_evaluate`, if `best is None`, record an `alarm_no_best` event, and either set
  `control.should_training_stop = True` or apply a guard switch.

### B2-3 (low): malformed control.json stops the run
- **Reproduction:** `test_b2_3_…`. Write `{"reward_weights": {"g_area": 0.5`, which is what an editor can leave
  mid-save. `on_step_begin` raises JSONDecodeError. The mtime is recorded before parsing, so even a retry would
  skip the file.
- **Fix:** catch `json.JSONDecodeError` / `OSError`, log it, and retry at the next step without recording the
  mtime. Document writing the file atomically (write to a temp file, then `os.replace`,
  https://docs.python.org/3/library/os.html#os.replace).

### B2-4 (medium): sub-filter raters create a "nothing" topology mode that specks match
- **Data:** in `data/fixtures/lidc_full.npz`, 21 items have a rater whose size-filtered signature is `()`:
  every component of that rater's mask is below `min_area` (5 px), with areas of 2–11 px (the 9–11 px masks are fragmented).
  On **18 of the 446 disagreement items (4 %)**, that rater is the only source of the disagreement. Examples: `472765997306_c2_k107` (4 / 41 px),
  `194067808002_c0_k262` (4 / 5 / 95 px), `258447180130_c0_k106` (33 / 2 px). `prepare` flags these items as
  disagreement items, so `grpo_rows` repeats each one ×4 (`DISAGREEMENT_OVERSAMPLE`).
- **Reproduction:** `test_b2_4_…` (synthetic 4 px + 42 px item). On the real items
  (`scratchpad/breakit/B2/speck.py`), the rollout that copies the tiny rater gets `iou_tc 1, f5 1, g_area 1,
  topo_freq 1/k`, valid. Its r_topo is 1.55–1.93 and its advantage 1.68–2.77, against an honest median of 0.20–0.58
  in a group of the big rater's mask and its ±1 px variants.
- **Expected vs actual:** the size rule says features under `min_area` are noise, so "nothing" should not be a
  rater answer that GRPO pushes the policy toward on 1/k of its samples. Actual: frequency matching treats it as a
  real mode. A ≤ 4 px blob at that spot is rewarded as a perfect, valid answer. The eval harness counts it as a
  rater mode too (TV, coverage). `{"nodule":[]}` itself is not valid: IoU 0 against the tiny mask.
- **Root cause:** `ExpertCache` builds `rater_sigs` / `same_sig` from every rater. IoU_tc is computed on raw
  masks, so a sub-filter output can match a sub-filter rater.
- **Fix (pick one):** drop raters whose filtered signature is `()` from the topology histogram (keep them in the
  band and area range); or do not count an item as a disagreement item when the only differing answer is `()`; or
  exclude these items from RL. Also record the count in M4/M7 (I11's 308 / 313 disagreement-item counts include
  these).

### B2-5 (medium): groups of 16 are split silently
- **Reproduction:** `test_b2_5_…`. `configs.grpo_config(..., num_generations=16)` is accepted. TRLReward's check
  (blocks of 8 with one key) passes, because 16 equal keys form two such blocks. It returns `advantages_v4` on two
  independent halves. On 16 shifted squares, rollout 9's advantage is +0.725 under the halves and −0.023 under one
  group of 16. TRL then centres over 16. With a guard weight w > 0, the logger term is centred over 16 but v4 over
  8. The check cannot tell this apart from an oversampled disagreement item that drew two adjacent rows, which also
  gives 16 equal keys.
- **Fix:** `grpo_config` should refuse `num_generations` (and `num_generations_eval`) other than
  `REWARD_V4.group`. Or give TRLReward an explicit `group` and check it against `trainer.args` once, for example
  in the Controller's `on_train_begin`.

### B2-6 (medium): the duplicate alarm misses identical masks with different spellings
- **Reproduction:** `test_b2_6_…`. Eight completions of the same square either nudge one vertex by 0–7 grid units
  (sub-pixel at 20 units/px) or append 0–7 spaces (`parse` strips them). All eight masks are identical, so the
  reward sees one answer and drops the group as low-spread (advantage 0 everywhere). But `groups.dup = 0.0`, the
  `duplicates` alarm does not fire, and the Phase 2 `grpo_duplicates` gate passes.
- **Why it matters:** at temperature 1, a policy that has collapsed to one mode still produces sub-pixel vertex
  jitter. Measuring duplicates on text underestimates collapse, so the plan §9 duplicate alarm (a rollback trigger)
  becomes insensitive. The `grpo_sigma` gate does catch this group (signal False), but the audit alarm does not.
- **Fix:** hash the canonical answer instead of the text: the packed rasters of every object (for example
  `sha256(np.packbits(mask))`), or at least the canonicalized, re-serialized parse. Keep the text-level rate as a
  separate log line.

### B2-7 (low): RIGA validity gate on the mean IoU_tc
- **Reproduction:** `test_b2_7_…`. A perfect disc (IoU 1) and a cup outside the disc (IoU 0) average to 0.5, which
  passes the gate, so the rollout gets the full r_topo (1.0, the same as the honest group). `metrics` marks it
  `passed = 0`, because M8 deliberately chose a per-object rule.
- **Impact:** this is not a win over the honest median: the rollout still loses on IoU, F5 and g_area. But the
  2-unit topology penalty for missing the cup entirely disappears whenever the disc is good, and the reward and the
  eval disagree on what "valid" means. The rule is the lab rule (M6 "Multi-object completions"), but it is not
  listed in ISSUES.
- **Fix:** gate validity on min over objects of IoU_tc ≥ 0.5, matching eval. v4 is frozen, so this needs a decision.
  If it stays, record it as a documented residual.

### B2-8 (low): cache_loader after a fingerprint change
- **Reproduction:** `test_b2_8_…`, with the fingerprint monkeypatched. `ExpertCache.load` rebuilds from the stored
  masks as designed, but the rebuilt digest embeds the new fingerprint, so `cache_loader` raises "cache digests
  differ from the prepared records". That message reads as tampering. A numpy or OpenCV patch bump on the Vast box
  (both are in the fingerprint) therefore stops training at the first reward call with the wrong diagnosis.
- **Fix:** compare digests computed under the record's fingerprint (store it in the record or the manifest), or
  raise "reward fingerprint changed: re-run prepare_rl". Separately, scipy (`excursion_moment`: quad, norm) is not in
  the fingerprint, although T5 depends on it.

### B2-9 (low): rater count mismatch across objects
- **Reproduction:** `test_b2_9_…`. With `{"disc": k=2, "cup": k=3}`, `score_completion` returns a joint `topo_freq`
  over only 2 raters and gives no error. With the order swapped, it raises `IndexError: tuple index out of range`.
  The contract ("rater r of one object's cache is rater r of the others") is not enforced, and permuted rater order
  cannot be detected at all, because caches carry no rater ids. This cannot happen through `prepare`, which builds
  every object from the same item.
- **Fix:** in `score_completion` (or `cache_loader`), check that every cache of an item has the same `k`, and raise
  a `ValueError` naming them. Optionally store the rater ids in the cache file.

### B2-10 (low; makes I1 worse): the cache certifies a wrong topology
- **Reproduction:** `test_b2_10_…`. `ExpertCache.build` is given a bool stack whose True bytes are 255 (I1: for
  example `np.asarray(PIL "1" image, bool)`). It keeps those bytes (`np.ascontiguousarray(dtype=bool)` does not
  normalize them), so the rater signatures lose their holes: `((0,), (0,), (0,))` instead of `((1,), (1,), (0,))`.
  The digest is identical to the correct cache's, because `np.packbits` normalizes. After `save` / `load`, the
  bytes are still 255 and the signatures still wrong. So an item whose raters disagree on holes would be recorded as
  unanimous, and nothing downstream can detect it. The vote, core, τ and T5 are unaffected: numpy's sums normalize.
  The current store unpacks bits (`store.py`), so this cannot happen in today's pipeline.
- **Fix:** normalize with `raters.view(np.uint8) != 0` in `ExpertCache.__init__` and `score_mask`, as I1 proposes for
  topology.

## Attacks that held (29)

| # | attack | result |
|---|---|---|
| 1 | empty object list `{"nodule":[]}` | legal, IoU 0, g_area 0, topo invalid: worst legal score |
| 2 | one polygon over the whole crop | IoU 0.19, F5 0.01, g_area 0.19: loses to every honest rollout |
| 3–8 | duplicated part, nested outers, parts sharing an edge, hole = outer, hole wound clockwise, CCW outer, zero-area spike | all illegal (`ring_ring_touching`, `nested_outers`, `wrong_winding`, `zero_area`/`spike`), gated to the minimum |
| 9 | 10 parts of 0.5 px (no pixel centre) | legal empty raster, scores like the empty answer |
| 10 | 11 rings | `cap_rings` format failure |
| 11 | surrounding whitespace, including Unicode (U+3000, U+2028) | accepted; content scored normally (exploitable only via B2-6) |
| 12 | valid text flagged truncated | format failure (D10) |
| 13 | RIGA disc without cup, cup empty | cup g_area 0, IoU 0, topo invalid |
| 14 | groups: all illegal / all identical / all format failures / 1 valid among 7 illegal / 1 valid among 7 format failures / σ exactly at the floor / a lone minority answer | finite, sum 0 (≤ 1e-14), max \|A\| ≤ 6.7 |
| 15 | `advantages_v4` with G = 2, 3, 4, 16 | finite. G = 1 divides by zero, but G = 1 is unreachable (`metrics` needs K ≥ 2) |
| 16 | B·G not a multiple of G; partial or interleaved TRL groups | `ValueError` (existing tests) |
| 17 | determinism: serial, 8 threads, `cv2.setNumThreads(0)` + 8 threads, 3 spawn processes with chunk = 5, on a real 421 px RIGA fixture item | identical SHA-256 |
| 18 | `sig_id` CRC32 collision between topology answers | none among 144,273 distinct signatures (CRC32 detects all burst differences ≤ 32 bits; small tuple reprs differ in short bursts) |
| 19 | worst-case scoring time, 500×500 crop with 6 raters | honest 0.9 ms, whole crop 1.6 ms, 300-vertex comb 7.1 ms (5.7 k characters, above the token budget) |
| 20 | W/H not matching the cache | `ValueError` in `score_mask` |
| 21 | single-rater cache | refused by `build` |
| 22 | tampered cache file; cache from another item | SHA mismatch / digest mismatch `ValueError` |
| 23 | `cache_loader` LRU eviction | reloads bit-identical (test kept) |
| 24 | TRLReward without `trainer_state`; repeated calls in one step; two batches with one key | `**_` absorbs it; the identity key `self._last[0] is completions` holds a live reference, so ids cannot be reused |
| 25 | eval: K = 1, empty samples, all truncated, K varying across items, no baseline, empty item list, missing datasets | refused / TV 1, retention 0, advantage 0 / pass@k up to min K, `group` None / `components` None / empty audit / gate fails |
| 26 | empty rater masks (the empty answer would get IoU 1, F5 1 and a valid mode) | no empty rater mask in any of the three fixture files (11,923 stacks); the reachable variant is B2-4 |
| 27 | illegal rollouts sharing a valid rollout's signature | they lower its r_topo through π_loo (1.0 vs 2.32 among format failures); an illegal rollout's own advantage falls if it changes signature, so no exploitable incentive |
| 28 | NaN `best` in Controller | `best` would stick at NaN, but unreachable: no reward function returns None or NaN |
| 29 | `ExpertCache.load` of a uint8 / 0–2-valued raters array | `.astype(bool)` normalizes |

## Sources

- TRL v1.16.0 `grpo_trainer.py` (read in `.venv`): `_calculate_rewards`, `sum_then_normalize`, `nanmean` /
  `nan_to_num`, `mask_truncated_completions`, metric flush: https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/grpo_trainer.py
- Yao et al. 2026, multimodal reward hacking (RHR / NRFR; checkpoint selection on a held-out set): https://arxiv.org/abs/2607.09492
- DAPO (dynamic sampling; overlong filtering, the `mask_truncated_completions` analogue): https://arxiv.org/abs/2503.14476
- CRC error-detection properties: https://en.wikipedia.org/wiki/Cyclic_redundancy_check
- Atomic file replacement: https://docs.python.org/3/library/os.html#os.replace
- M8 demo IoU_tc of hedges vs honest output: `experiments/phase1/impl/M8_eval.md` (demo table)
