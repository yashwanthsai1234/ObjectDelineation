# Phase B test cases, set 1 (written from B.spec.md only)

## Skills used
- engineering:testing-strategy (loaded): pyramid split (unit > property > integration > GPU > end-to-end), one invariant per test, perf gates listed separately.
- The other named skills were not loaded for this job: it needed test design from the spec only (no design context, no web research).
- Repo style read: tests/test_reward.py (hypothesis, oracles, `disk()` helpers), tests/test_eval_coverage.py (small inline fixtures), reward `advantages.py` / `spec.py` (current v4 term names: iou_tc, f5, legal, g_area, topo_any, topo_freq, sig_id).

Conventions: [spec s.x] = section of B.spec.md. "Fixture" tests build masks with the `disk(n, r, cy, cx)` helper. Level GPU tests are skipped when CUDA is absent (as in test_cuda_checker_gpu.py). All numbers below were derived by hand from the spec formulas; items marked (derive) must be recomputed once in the test with an independent numpy oracle, not copied from here.

---------------------------------------------------------------------------------------------------------

## A. Cache fingerprint and val halves [spec 3.0]

**1. fingerprint_ignores_reward_weights** (unit)
Setup: monkeypatch/replace each of `alpha_loc, miss_iou, w_f5, w_area, c_topo, adv_scale, loc_width, area_gate, lam, eps, drop_spread, iou_gate` on RewardV5 and each RewardV4 weight/floor, then call `fingerprint()`.
Expect: digest identical to baseline for every one.
Why: if the fingerprint still covers weights, every reward tweak invalidates the whole store (silent mass rebuild or StaleCacheError mid-run).

**2. fingerprint_tracks_cache_inputs** (unit)
Setup: vary `CacheSpec.expert_score` (0.95 -> 0.94), `tau_floor` (1.0 -> 1.5), `tau_quantile` (95 -> 90), patch the size rule (`min_area`), patch `np.__version__`, `scipy.__version__`, `cv2.__version__` one at a time.
Expect: each change gives a digest different from baseline AND from every other change (7 distinct digests).
Why: a stale cache accepted after tau/T5 inputs change gives wrong rewards without any error.

**3. stored_values_survive_new_fingerprint** (regression, integration)
Setup: build an `ExpertCache` from the three-disk fixture (radius 18/19/20, n=60) on the old code path (golden saved before the change: tau, T5, ref values); rebuild after the change.
Expect: tau, T5 bit-equal to the golden; only the digest field differs.
Why: spec says cache values are unchanged; a refactor of CacheSpec that alters rounding would shift every reward.

**4. val_half_rule** (unit + property)
Input: 5000 random group_ids (uuid strings, empty string, unicode "Ünïcode", very long 10 000-char string).
Expect: `val_half(g) == "shape"` iff `hashlib.sha256(g.encode()).digest()[0] % 2 == 0`; result in {"shape","select"}; repeat call identical; two items with same group_id always same half.
Why: split leakage between the half used for targets and the half used for gates.

**5. val_half_balance** (property)
Input: 2000 distinct ids. Expect: shape share in [0.45, 0.55] (binomial 4 sigma = 0.045). Why: a degenerate hash (e.g. using the wrong byte) would put all patients in one half and empty the target set.

---------------------------------------------------------------------------------------------------------

## B. windows: cell, offset_of, Window [spec 3.1]

**6. cell_edges_go_up** (unit, table)
Input (extent, offset) -> expected:
(9,0.149)->(0,0); (10,0.15)->(1,1); (19,0.1999)->(1,1); (20,0.2)->(2,2); (39,0.2999)->(2,2); (40,0.3)->(3,3); (0,0)->(0,0); (10**6, 10.0)->(3,3); (10,0.0)->(1,0); (9,0.3)->(0,3).
Also assert float64 values exactly on the edge `0.15`, `0.2`, `0.3` (literal) go up, and np.float32(0.15) is documented (float32(0.15) < 0.15 in float64: expect bin 0 -- test pins the behaviour so a dtype change is noticed).
Why: `searchsorted(side="right")` vs "left" off-by-one is the classic bug; wrong bin shifts the whole target law.

**7. offset_of_known_geometry** (unit)
Input: union_box (40,40,60,60) (centre (50,50)), window_box (0,0,100,100) (centre (50,50), S=100) -> 0.0.
(60,40,80,60) centre (70,50) -> dist 20 / 100 = 0.2 -> cell offset bin 2.
Union (0,0,10,10) centre (5,5), window (0,0,100,100): dist = sqrt(45^2+45^2)=63.64 -> 0.6364 (derive).
Non-square rejection: window_box with width != height -> either ValueError or documented; test pins behaviour (S = width).
Half-open: box (0,0,1,1) single pixel: centre (0.5,0.5), extent 1 (not 0).
Why: half-open vs closed boxes shifts offsets by 0.5 px; S in the denominator must be the source-slice side, not the model input size.

**8. extent_single_pixel_and_thin** (unit)
Input: union of one pixel -> extent 1, size_bin 0. A 1 x 200 line -> extent 200, bin 3. Two distant doctors' masks (union box spans both) -> extent = union box, not per-doctor.
Why: extent must be of the UNION box [spec 2]; per-doctor or mask-sum extents change cells.

**9. window_id_format** (unit)
Input: item_id "lidc/LIDC-IDRI-0001/n3/s17" kind "shifted" draw 2.
Expect: window_id == item_id + "~shifted2"; the part after the last "/" of the item_id... precisely: `window_id.rsplit("/",1)[0] == item_id.rsplit("/",1)[0]` and `"/" not in window_id[len(item_id):]`. Also kind "centred" draw 0, "rescaled" draw 11.
Why: file paths/caches are derived from the id; a "/" in the suffix creates a nested folder and breaks path-uniqueness.

---------------------------------------------------------------------------------------------------------

## C. windows.pool [spec 3.1]

Fixture: synthetic Item set built in tmp store: 40 items, 300 x 300 .. 512 x 512 slices, nodules disks r = 3..40, 1..4 raters, some items with a second outlined nodule in the same slice, some with 1 rater.

**10. pool_not_rl_dropped** (unit)
Input: 10 items, 4 with one rater. Expect: `dropped["not_rl"] == 4`; no window has an item_id of those; 6 centred windows (draw 0 each).
**11. pool_one_centred_per_item** (property)
Expect: for every RL item exactly one window with kind "centred" and draw 0, for any draws in {0,1,5,16} (draws=0: still one centred; shifted/rescaled none, no drop counted).
**12. pool_at_most_one_per_kind_per_draw** (property)
Expect: count of (item, kind, draw) tuples is unique; shifted/rescaled count per item <= draws.
**13. pool_sides** (property)
Expect: every shifted `side == S0(item)`; every rescaled `side != S0`; every centred `side == S0`.
**14. pool_containment_and_clearance** (property, hypothesis over disks position/size)
For each window: all doctor-union pixels in [dx, dx+side) x [dy, dy+side), min distance from any doctor pixel to the window border >= 2 px (box edge pixel index: x - dx >= 2 and dx+side-1-x >= 2); no other outlined nodule's pixel in the window (check against the item's other-nodule masks).
Why: D4/D5 gates; a window that clips a doctor pixel teaches a truncated target.
**15. pool_nodule_at_slice_edge** (edge)
Input: nodule touching the slice border (x0 = 0). Expect: no window, because clearance 2 px outside the slice is impossible, OR a padded window if the spec allows padding (spec D14 mentions "padded share", so padding exists): test both branches: window box may exceed the slice, doctor pixels still inside with >= 2 px clearance; padded share computed. Pin whichever the implementation chooses and assert D4 holds either way.
**16. pool_nodule_larger_than_window** (edge)
Input: union extent bigger than S0 - 4. Expect: shifted draws all dropped, `dropped["shifted_none"] == draws`, centred window may still exist only if it satisfies D4/D5 (if not, that item raises or is dropped; test pins). No fallback window returned.
**17. pool_adjacent_nodule_blocks_window** (adversarial)
Input: second outlined nodule 3 px from the target; every window of side S0 that offsets toward it must be rejected. Expect: all windows exclude it; if none exist, `shifted_none` incremented.
**18. pool_rescaled_never_s0** (adversarial)
Input: item whose only admissible rescaled side equals S0 -> expect rescaled dropped (`rescaled_none` += draws) not a kept-side window. Also assert `fallback`/`kept_side` flags never present in pool windows.
**19. pool_matches_build_eval_views_draw0** (integration)
Input: same items, seed=0, call `experiments.phase2.build_eval_views` builder for shifted and rescaled; compare with pool(draws=1, seed=0).
Expect: for each item, draw-0 shifted (dx, dy, side) equal; draw-0 rescaled equal where the eval view is not a kept-side window; where it is, pool has no rescaled draw 0 and `rescaled_none` counts it.
Why: R5 / B1 comparison is only fair if the train law's draw 0 is the exam's recipe.
**20. pool_sorted_and_order_independent** (property)
Expect: `windows == tuple(sorted(windows, key=window_id))`; `pool(shuffle(items)) == pool(items)` for 20 random permutations; `pool(items) == pool(items)` (equal dataclasses and dropped dict).
**21. pool_seed_and_draw_independence** (property)
Expect: pool(draws=3, seed=0) windows with draw <3 are the same as in pool(draws=5, seed=0) (prefix stability: adding draws must not reshuffle earlier ones); different seeds give different shifted offsets for >= 90 % of items.
Why: a later "more draws" change silently changes the audited rows.
**22. pool_empty_inputs** (edge)
Input: items=[] -> `WindowPool((), {"not_rl":0,"shifted_none":0,"rescaled_none":0})`. Input: all not RL -> no windows, not_rl == len.
**23. pool_dropped_keys_always_present** (unit)
Expect: dropped has exactly the 3 keys, ints, even when 0.
**24. pool_disagreement_flag** (unit)
Expect: `disagreement` True iff raters disagree per the repo's existing rule (same value as in prepared records today); constant across a nodule's windows.

---------------------------------------------------------------------------------------------------------

## D. cell_shares and schedule [spec 3.1]

**25. cell_shares_basic** (unit)
Input: 3 windows in cell (0,0), 1 in (3,3). Expect: array shape (4,4), [0,0]=0.75, [3,3]=0.25, sum 1.0 within 1e-12, dtype float64. Empty list -> ValueError or zeros (pin; must not return NaN from 0/0).
**26. schedule_kind_counts_largest_remainder** (unit, table)
(steps, per_step) -> expected (centred, shifted, rescaled) with shares (0.15,0.75,0.10):
(10,32): total 320 -> (48,240,32).
(7,32): 224 -> raw 33.6,168,22.4 -> (34,168,22).
(1,32): 32 -> raw 4.8,24,3.2 -> floors 4,24,3 =31; largest remainder .8 -> (5,24,3).
(1,1): total 1 -> raw .15,.75,.10 -> (0,1,0).
(3,5): 15 -> 2.25,11.25,1.5 -> floors 2,11,1 =14; remainders .25,.25,.5 -> rescaled gets it: (2,11,2).
Tie case: shares (1/3,1/3,1/3), total 32 -> floors 10,10,10, 2 extra: tie broken deterministically (assert same output on 10 calls; the sum is 32; the chosen two kinds are documented).
Why: Python round() (banker's) or int() truncation gives totals != steps*per_step.
**27. schedule_cell_counts_largest_remainder** (unit)
Input: target uniform 1/16 over 16 cells, n = rows of kind. For n = 240: 15 each exactly. For n=241: 15 each + 1 extra cell (largest-remainder; all remainders equal -> deterministic tie). Peaked target (0.5, rest spread): check against an independent largest-remainder implementation in the test for 200 random targets (hypothesis) and random n.
**28. schedule_length_and_step_major** (property)
Expect: len == steps*per_step; rows `[per_step*s : per_step*(s+1)]` have no duplicate item_id (D11 and spec). steps=0 -> [] (or ValueError; pin). per_step=1 allowed.
**29. schedule_no_repeat_window_ids** (property)
Expect: all window_ids distinct over the run; also test pool just big enough (n = needed + 0): must succeed; pool one window short -> ValueError.
**30. schedule_error_messages** (unit)
- target shape (3,4) -> ValueError; target with -1e-3 entry -> ValueError; target summing 1.0000001 (within 1e-9? no: 1e-7 off) -> ValueError; sum off by 5e-10 -> accepted; NaN in target -> ValueError (NaN fails the sum test).
- shares missing "rescaled" key -> ValueError; extra key -> ValueError; negative; sum 0.99 -> ValueError.
**31. schedule_unfillable_names_kind_and_cell** (adversarial)
Input: pool with zero shifted windows in cell (2,3), target mass 0.1 there.
Expect: ValueError whose message contains "shifted" and "(2, 3)" (or "23"; test accepts both formats); no centred/rescaled window substituted (verify by catching the error, not by output). Also the case where cell has windows but all share items with other rows in each step (item-uniqueness impossible: only 2 items, per_step 3) -> ValueError.
**32. schedule_zero_target_cells** (edge)
Input: target with 0 mass in many cells: those cells get 0 rows even though the pool has windows there. Target with all mass in one cell: all off-centre rows in that cell; D7 ratio unaffected.
**33. schedule_determinism_and_seed** (property)
Expect: same args -> equal lists (also across processes: run in subprocess with PYTHONHASHSEED=0 and =123, compare ids); seed 0 vs 1 differ in order for a pool >> need.
**34. schedule_disagreement_weight_monotone** (property)
Input: pool with 30 % disagreement items; weights 1, 4, 16 with same seed (repeat 20 seeds).
Expect: mean disagreement share non-decreasing in weight across averaged seeds; weight 1.0 -> share approx pool share (within 3 sigma); weight 0 -> no disagreement rows if pool allows (or ValueError for negative; pin). Per-seed monotone is NOT guaranteed by spec, so assert on the 20-seed mean and also "(larger or equal)" for equal-seed pairs only as xfail-documented if it fails.
**35. schedule_does_not_mutate_inputs** (unit)
Expect: windows list, target array unchanged after the call (copy before, compare). `shares` default dict not mutated (mutable default bug).
**36. schedule_centred_cell_handling** (edge)
Expect: centred rows are scheduled without cell constraint (the per-cell rule applies to off-centre kinds only); a pool with centred windows in any cell works.
**37. schedule_regression_numbers** (regression)
Input: the real train pool at steps=100, per_step=32, seed 0 (integration, skipped if data absent). Expect: kind counts (480, 2400, 320) = 0.15/0.75/0.10 x 3200 exact; D7 >= 0.8, D8 <= 0.10, D9 <= 0.06 against val-shape target; D10 >= 0.75 -- exactly the spec gates, the tests print the values.

---------------------------------------------------------------------------------------------------------

## E. coverage and quality [spec 3.2]

**38. eval_coverage_moved** (unit)
Expect: `import objdel.data.coverage` works with the public functions props, distribution, gaps (existing tests in tests/test_eval_coverage.py pass after changing the import only); `importlib.import_module("objdel.eval.coverage")` raises ModuleNotFoundError; `objdel.eval.audit` still imports.
**39. cell_coverage_known** (unit)
Input: target uniform 1/16; scheduled 160 off-centre windows exactly 10 per cell; target_offsets spread uniformly in [0, 0.45].
Expect: worst_ratio == 1.0, tv == 0.0, ks small and < 0.06 (derive). Remove all rows of one cell: worst_ratio == 0.0, tv == 1/16 = 0.0625 (derive: 0.5 * (1/16 + 15 x (1/15-1/16))... recompute in test).
Skewed: 50 % of rows in cell (0,0): tv = 0.5*(0.5-1/16 + 15*(0.5/15-1/16))... pin via independent numpy.
Cells with target share < 1 % excluded from worst_ratio (target cell share 0.005, row share 0 -> ignored).
Centred rows in input must be ignored (add 100 centred rows; metrics unchanged).
Empty off-centre rows -> keys present, worst_ratio 0 (or error); D7/D8 must fail, not crash.
**40. gates_return_16_named** (unit)
Expect: `[g.name for g in gates(report)] == ["D1",...,"D16"]`; D14, D16 hard False; others True; `limit` is a str; `ok` bool.
**41. failures_only_hard_and_not_ok** (unit)
Input: report where D14 and D16 are "bad" and D5 fails. Expect: failures == [D5 gate]; with only D14/D16 bad -> [].
**42. D1_planted_cross_split_patient** (integration, tmp store)
Input: store where patient P in train and val. Expect: gate D1 ok False; clean store ok True; series overlap (different patient but same series uid) also fails; item overlap fails.
**43. D2_planted_hash_across_patients** (integration)
Input: two patients with byte-identical slice -> D2 fails. Within-patient duplicate (same patient, same crop SHA, two slices) -> D2 ok (the spec says across patients or splits). 32 x 32 nodule crop identical across patients -> fail.
**44. D3_near_duplicate_rate** (unit)
Input: descriptors with 4 of 7341 near-dups in train, cross-split rate 5/7341 -> fail; 4/7341 -> pass (<=); threshold cos 0.98 exact edge: pair at cos 0.98 counts as dup (>=).
**45. D4_cut_pixel** (unit)
Input: window with box clipping exactly one doctor pixel -> D4 fails, value < 100 %; value reported as fraction. No cut -> 1.0.
**46. D5_centred_with_other_nodule** (unit)
Input: centred (policy) window containing another outlined nodule's pixel -> D5 fails even though kind is centred; clearance 1 px -> fail; exactly 2 px -> pass.
**47. D6_fallback_and_kept_side** (unit)
Input: schedule containing a record/window with `fallback=True`, another with `kept_side=True` -> D6 fails with value >= 1 each; the "centred" kind is NOT a kept-side window (no false positive).
**48. D7_D8_no_offcentre_rows** (unit)
Input: schedule of only centred rows -> D7 and D8 fail (not NaN-pass, not ZeroDivisionError).
**49. D12_token_limit** (unit)
Input: tokenizer stub giving target length `limit - 100` -> D12 must pass? Spec: "target of more than `limit - 100` tokens fails". Pin: `limit-100` passes, `limit-99` fails; gate condition "max tokens + 1 < MAX_COMPLETION - 100" (cross-check with 5.1; resolve the off-by-one in the test with both spec sentences and report any inconsistency as a finding). Target error count > 0 -> fail.
**50. D10_policy_side_share** (unit)
Input: 100 off-centre rows, 75 at S0 -> pass (>=); 74 -> fail. Rows are shifted at S0 + rescaled non-S0.
**51. D11_exact_counts** (unit)
Input: schedule off by one in a single cell -> D11 fails; all correct -> passes; compare to independent largest-remainder.
**52. D13_fresh_render** (integration)
Input: 150 sampled windows; tamper one image byte -> exactly 1 mismatch -> gate fails; sample of 150 is deterministic (same sample for same seed); fewer than 150 windows -> samples all.
**53. D15_frozen_views** (integration)
Input: frozen_views.json; recompute content_sha256 per spec (newline-joined `(item_id, view.offset, view.side, image_sha256)` ordered tuples). Tamper one image sha -> D15 fails; change `versions` or `reward_fingerprint` in manifest -> D15 still passes.

---------------------------------------------------------------------------------------------------------

## F. Scripts [spec 3.3]

**54. build_rl_rows_layout_and_schema** (end-to-end, tiny store: 30 items, steps 3, per_step 4, draws 2)
Run `python -I -m experiments.phase2.build_rl_rows ...`. Expect:
- out has images/, caches/, lidc/train.jsonl, manifest.json;
- train.jsonl lines count == 12 and ordered == schedule order (`ids` NOT sorted unless by chance; assert not equal to sorted for a seed where schedule isn't sorted);
- every record: item_id == window_id, source_item_id, cell is str of 2 digits, group_id, `view == {kind, draw, offset:[dx,dy], side, fallback:false, kept_side:false, op:"id"}` (exact key set), rl true, plus the fields of today's prepared records (compare key sets with a record from the existing `prepare` output);
- all image paths unique and all cache folders unique, files exist;
- manifest: keys files, sets, reward_fingerprint, versions, schedule with exact keys {steps, per_step, seed, draws, shares, target_view, target_half, target_sha256, ids_sha256, dropped, gates}; target_half == "shape"; ids_sha256 == sha256("\n".join(ids)) (no trailing newline: pin).
**55. build_rl_rows_hard_gate_failure** (end-to-end)
Plant a leak; expect exit code 1, no manifest.json, (train.jsonl may exist), stderr names the gate.
**56. build_rl_rows_target_from_val_shape_only** (integration)
Input: VIEW_DIR where val-select items have an extreme offset law. Expect: the manifest target_sha256 equals sha of a target computed from val-shape only; changing val-select records leaves target_sha256 and the schedule unchanged. Training rows never include a val-select item (leakage property: set of group_ids of rows intersect val-select groups == empty). 
**57. build_rl_rows_reproducible_and_procs_invariant** (integration)
Same args run twice and with `--procs 1` vs `--procs 4`: identical train.jsonl bytes and identical image SHA-256s. Why: multiprocessing nondeterminism (set iteration, completion order).
**58. build_rl_rows_steps_too_large** (failure)
Pool too small for steps -> exit non-zero, message names kind and cell; no partial manifest.
**59. audit_data_exit_codes** (end-to-end)
Clean root -> prints 16 lines D1..D16, exit 0; root with hard failure -> exit non-zero; root with only D14/D16 "issues" -> exit 0; `--json FILE` writes parseable JSON of 16 gates.
**60. build_eval_views_kinds** (integration)
`--kinds shifted` builds only shifted (no rescaled/rotated records or images written); `--kinds shifted,rescaled` exact two; no `--kinds` -> three as today (regression: record counts equal to the stored lidc-test-shifted view); unknown kind -> argparse/ValueError. Shifted-only output equals the shifted records of the default build bit for bit (same image SHAs).
**61. final_eval_frozen_check** (end-to-end)
Tamper one `view.offset` in the test view -> non-zero exit with message mentioning content_sha256; tamper target set -> message mentions target_set_sha256; change manifest versions/reward_fingerprint only -> proceeds. A non-test view (val) is not checked.
**62. check_series_z** (integration, XML fixture)
Fixture XML: two SOP UIDs with same z 123.5 in series A, unique in B -> lists the shared z for A with both SOP UIDs; B empty; with no DICOM folder the script opens no `.dcm` (monkeypatch pydicom.dcmread to raise).
**63. fp16_vision_check_exit** (GPU, tiny stub model for unit)
Stub tower producing inf -> exit 1 and prints count > 0; clean stub -> exit 0 with count 0. Real run: Qwen3-VL-2B vision tower fp16 on 32 train images prints count (record it, expected 0 -> report).
**64. b2_attribution_output** (GPU)
Prints three means for (a),(b),(c); all finite non-negative; (c) with LoRA delta = 0 equals (b) within 1e-6 (merge of zeros is identity). Does not modify the checkpoint on disk (hash before/after).

---------------------------------------------------------------------------------------------------------

## G. train_rl guards and configs [spec 3.3, 3.4]

**65. train_rl_row_guards** (unit, no GPU; stub trainer)
Each of these must raise ValueError before any model load (assert model loader not called):
- row count 95 vs schedule 3 x 32 = 96; 97;
- same count but one id swapped -> ids_sha256 mismatch; same ids different order -> mismatch;
- `generation_batch_size // num_generations` = 16 vs per_step 32;
- `steps_per_generation != gradient_accumulation_steps`;
- `num_iterations == 2`;
- `--steps 101` with schedule.steps 100.
Passing case: all consistent -> proceeds; `shuffle_dataset is False` in the GRPOConfig. Without `schedule` in manifest -> old path (shuffle True, disagreement x4) until gate.
**66. train_rl_cli_after_b1** (unit)
`--train-views` -> argparse SystemExit(2); root without schedule -> refused; `dataset.grpo_rows_views` absent (`hasattr` False); serve key tagging absent. (xfail/skip until gate B1 flag says passed.)
**67. numerics_presets_exact** (unit)
Table check of NUMERICS keys and values: bf16-seqmask (sequence_mask, clip_max 3.0, clip_min None); bf16-tis (token_truncate, 2.0, None); bf16-geomask (sequence_geo_mask, 2.0, geo (0.99, 1.01)); fp16-tis = bf16-tis + fp16 True, bf16 False. `grpo_config(...)` for every preset: never both bf16 and fp16; unknown preset -> KeyError/ValueError; default equals `DEFAULT_NUMERICS`.
**68. grpo_config_group_check** (unit)
`grpo_config(..., reward=R, num_generations=4)` with R.group 8 -> ValueError; `num_generations_eval=4` -> ValueError; 8 and 8 pass.
**69. rl_json_additions** (integration with stub trainer)
`rl.json` has `rows={root, ids_sha256}`, `numerics`, `reward={name,fingerprint,c_topo,adv_scale,loc_width,area_gate}` for v5; for v4 reward the v5 fields are absent.
**70. reward_cli_choice** : `--reward v6` -> argparse error; v4 and v5 accepted.
**71. budget_flags** : `--budget-usd 3` without `--rate-usd-per-h` -> argparse error (or documented); both given -> BudgetStop built with the exact values.
**72. prompt_rows_kind_cell_columns** (unit)
Records with `view` -> dataset rows have columns `kind`, `cell` equal to record values; records without `view` -> columns absent (no None columns that break TRL).
**73. trlreward_logging** (unit)
Batch of 2 prompts x 8: prompt A kind shifted cell "21" all 8 misses; prompt B kind centred cell "00" with hits. Expect metrics `kind/shifted/iou_tc`, `kind/centred/iou_tc`, `cell/21/iou_tc`, `cell/21/allmiss_groups == 1.0` (share of groups in that cell that are all misses), `cell/00/allmiss_groups == 0.0`. No `kind/` keys when columns missing. 7 misses + 1 hit -> not all-miss (0.0). Group boundary: groups are 8 contiguous answers; a cell with 3 groups, one all-miss -> 1/3.

---------------------------------------------------------------------------------------------------------

## H. Number match: is_ratio, is_summary, FastGRPO [spec 3.4]

**74. is_ratio_geo_mask_values** (unit, CPU torch float32)
logps - sampling per token (mask all ones, 4 tokens, C = 2.0, geo (0.99, 1.01)):
- seq A diffs [0.005]*4: geo = exp(0.005)=1.00501 in range -> ratios all exp(0.005)=1.00501.
- seq B diffs [0.02]*4: exp(0.02)=1.0202 > 1.01 -> all zeros.
- seq C diffs [1.0, -1.0, 0.0, 0.0]: mean 0 -> geo 1.0 kept; ratios [min(e^1, 2)=2.0, e^-1=0.36788, 1, 1] (the clamp is truncate only above).
- seq D diffs [-0.02]*4: exp(-0.02)=0.9802 < 0.99 -> zeros.
- seq E diffs [0.0099]*4: exp(0.0099)=1.00995 -> inside (just below 1.01) kept; [0.0101]*4 -> 1.01015 -> out; exact edge exp(mean)=a and =b when mean = log(a) / log(b) in float64: pin inclusive "[a, b]".
Padding: sequence of length 6 with mask [1,1,1,0,0,0] and padded tokens having diff 5.0: geo uses valid tokens only (mean diff 0.005 -> kept), padded ratio positions not counted/zeroed per existing convention.
Empty sequence (all mask 0): no NaN/Inf (0/0), result finite; pin keep or mask.
Why: the mean must be over valid tokens, and a divide-by-zero must not propagate NaN into the loss.
**75. is_ratio_existing_modes_unchanged** (regression, property)
For 200 random (logps, sampling, mask) and modes sequence_mask, token_truncate, token_mask, sequence_truncate: output bit-equal to the previous implementation (copy the old function into tests/oracles as oracle). `geo=None` default must leave behaviour unchanged. `mode="sequence_geo_mask"` with `geo=None` -> ValueError.
**76. is_summary_values** (unit; hand-derived)
Inputs for the 5 sequences above (mode sequence_geo_mask): masked count = 3 of 5 (B, D, E-out) etc: `is_masked_share == 3/5 = 0.6` (counted per sequence), `masked_len_mean = mean valid length of masked`, `kept_len_mean` of kept = 0 when none kept (test all masked); clamped share for token_truncate: tokens with raw ratio outside [lo, hi] AND clamped equal to the bound; a token with raw ratio exactly = hi is not counted (raw not outside); lo None -> no lower clamping count. seq_logdiff = signed sum of (logps-sampling): seq C gives 0.0, seq A 0.02; `seq_logdiff_max` = max over sequences (derive); `tok_absdiff_gt_0.5_share` = share of valid tokens with |d| > 0.5, strictly (d = 0.5 exactly not counted; seq C contributes 2 of its 4 tokens). Output float32 on the input device (CUDA test: `.device == logps.device`, `.dtype == torch.float32`).
**77. is_summary_no_host_sync** (GPU unit)
Use `torch.cuda.set_sync_debug_mode("error")` around accumulation calls across 4 micro-batches; no error raised; `.item()` called only in the logging call. (Add a CPU variant that wraps tensors in a recording subclass counting `.item()`/`.cpu()`.)
**78. fastgrpo_existing_5_metrics_unchanged** (regression)
Same inputs through old and new code paths -> the 5 `sampling/*` metrics equal bit for bit (golden values captured before the patch). 7 new metrics present exactly once per logging step (not once per micro-batch), names exactly as in spec.
**79. fastgrpo_geo_mask_requires_own_ratio_path** (unit)
Construct with mode sequence_geo_mask while the ratio path is inactive -> ValueError at construction (not at step 1).
**80. vllm_dtype_context** (unit, stub)
Replace `trl.generation.vllm_generation.LLM` with a recording stub: inside `with vllm_dtype("float16")` a call receives `dtype="float16"`; after exit (also after an exception inside the with) original object restored (`is` identical); nested contexts restore in order; a `dtype` the caller passes explicitly -> pin (override or ValueError). Name missing in TRL (private-name check) -> clear error in start-up check.
**81. fp16_tis_dtype_guard** (unit/GPU)
fp16-tis with engine model dtype torch.bfloat16 -> raises at construction; torch.float16 -> ok; other presets don't check.
**82. fp16_overflow_adversarial** (GPU, optional)
Feed values > 65504 through the fp16 forward on stub -> inf is counted by fp16_vision_check (test 63), proving the check can fire.

---------------------------------------------------------------------------------------------------------

## I. CUDA scorer sum_x / sum_y [spec 3.5] (only if the scorer change is made)

**83. sums_vs_numpy_random** (unit CPU build, GPU build)
Input: 200 random masks (shape 8..512, density 0..1, including all-zero, all-one, single pixel at (0,0) and at (H-1,W-1), 1-pixel wide line).
Expect: sum_x == `np.nonzero(mask)[1].sum()`, sum_y == `np.nonzero(mask)[0].sum()` as int64, per object; empty -> 0, 0. Fixture masks from tests/goldens too.
**84. sums_overflow_range** (edge)
All-one 4096 x 4096 mask: sum_x = 4096 * (4095*4096/2) = 34 359 214 080 (> 2^32): fits int64 and equals numpy int64. Why: int32 accumulation overflow would wrap.
**85. sums_cpu_gpu_bit_equal** (GPU)
Same 200 masks via the CUDA build and the CPU build: all fields equal (np.array_equal), both new and old ones.
**86. old_fields_bit_identical** (regression)
Golden files saved from the previous build for the existing fixtures: every previous field `np.array_equal` (and dtypes equal). Result key set == old key set + {"sum_x","sum_y"}.
**87. sums_polygon_and_mask_agree** (integration)
For polygons rasterised by `score_polygons` and the same raster passed to `score_masks`: sums equal. Multi-object masks (two disjoint blobs): one value per object in the object order of the other fields; hole pixels excluded.

---------------------------------------------------------------------------------------------------------

## J. Reward v5 [spec 3.6, 5.2]

Terms builder for tests: `mk(iou, f5=0, legal=1, g_area=1, loc=0, topo_any=0, topo_freq=0, sig=0)` -> dict of arrays of length B*8.

**88. RewardV5_validation** (unit)
`RewardV5(c_topo=1, adv_scale=1)` ok (defaults alpha_loc 0.1 == miss_iou 0.1 allowed: not `>`); `alpha_loc=0.11` -> ValueError; `adv_scale=0` -> ValueError; negative -> ValueError; missing c_topo/adv_scale -> TypeError; frozen: assignment raises FrozenInstanceError; `group` default 8.
**89. location_formula** (unit)
`location(0,0,0,(1,1,1),1.0) == 0.0` (area 0). Centroids equal: sum_x=1000, sum_y=1000, area=100; ref_sum=(2000,2000,200) -> exactly 1.0. Shift dx=3: sum_x=1300 -> d2=9; radius2=36, width 1: exp(-9/72) = 0.88250 (derive). width 2: exp(-9/288)=0.96923. Dy=4 and dx=3 combined: d2=25 -> exp(-25/72)=0.70670. Strictly decreasing over d in 0..30 (31 values); always in (0,1]; huge d (1e6) -> value 0.0 allowed (underflow) but not NaN; radius2 tiny (1e-9) with d2 = 0 -> 1.0, with d2>0 -> 0.0 no ZeroDivisionError/NaN. radius2 == 0 -> pin (ValueError or 0).
**90. loc_term_in_score_completion** (integration, disk fixtures)
Three-disk cache (radii 18-20, n=60 crop), completion = a legal r=6 disk 25 px away -> iou_tc = 0 -> L = 1 * 1 * g_area * [topo_freq>0] * location. Expect `terms["loc"]` equals formula by hand (call `location` with ref_sum/radius2 from the cache) within 1e-12. Key "loc" is a number (float) and present in score_batch for every row.
- two objects (multi-nodule target), one miss with loc 0.8, one hit (iou 0.6) -> loc = mean(0.8*g, 0) = 0.4 g (derive).
- IoU exactly 0.1 -> L = 0 (indicator `< miss_iou` false). iou 0.0999999 -> (1 - 0.999999) ~ 1e-6 * location, continuous, not 1.
- Unparseable, truncated, empty polygon list -> loc = 0 (and no exception).
**91. T1_monotone_and_continuous** (property; hypothesis)
Fixed f5, g_area, legal=1, r_topo, area_gate False. R(iou, loc) = iou + 0.1*loc*[iou<0.1]*max(0, 1-iou/0.1) + const. Sweep iou in 2001 points of [0,1] and loc in {0, 0.37, 1}: R non-decreasing (tolerance 1e-12; at alpha_loc = miss_iou and loc = 1 slope is exactly 0 below 0.1, must not be negative); |R(0.1-1e-9) - R(0.1)| < 1e-8; Use the real pipeline (score_completion -> advantages inputs) at least for one point to link formula to implementation. Also with `alpha_loc=0.05`, `miss_iou=0.2` configs.
**92. T2_covers_union_big_area** (integration)
Legal outline covering union and larger than doctors' max area so g_area=0 -> `loc == 0` even if centroid equals the reference. Check boundary: area exactly = max area (g_area presumably 1) -> loc > 0 only if iou_tc < miss_iou: an outline covering the union has iou_tc high -> L = 0 either way (covers case is a hit). Use the version with a ring-shaped legal outline to get iou < 0.1 and g_area = 0: loc 0.
**93. T3_speck** (integration)
Legal speck below min doctor area at reference centroid -> g_area=0 -> loc=0. Speck of exactly min area (g_area 1) at centroid -> loc = 1 * (1 - iou/0.1) * 1 * 1 -> iou tiny -> ~1 (derive). Pins the gate boundary.
**94. T4_two_blobs_centroid_trick** (integration, adversarial)
Item with one doctor component; answer = two disjoint legal blobs each above the size filter, outside the union, joint centroid == reference centroid. Expect loc == 0: via `topo_freq == 0` (component count mismatch) per spec. Also repeat with the blobs placed symmetric around the reference at distances 5, 10, 20 px. Also the variant "one blob in the union + one outside" (hit; loc 0 by the iou indicator).
**95. T5_shift_away** (property)
Legal doctor-area disk shifted d = 15, 16, ..., 60 px from the nodule (all IoU_tc = 0 by construction: assert it). loc strictly decreasing in d, values in (0,1) (derive vs gaussian formula). With topo_freq > 0 required: pick shapes of the same topology (single disk).
**96. T6_illegal_perfect_centroid** (integration)
Illegal outline (self-intersecting bow-tie with centroid at reference) -> loc 0 and R 0 (via advantages terms: legal=0 so R = c_topo * r_topo with r_topo = 0 -> R exactly 0); unparseable text, truncated text -> R 0.
**97. advantages_v5_hand_computed** (unit)
Group of 8, legal 1, g_area 1, f5 0, topo_freq 0 (r_topo 0), iou_tc 0, loc = [0.9,0.8,0.7,0.5,0.4,0.3,0.2,0.1], defaults (alpha 0.1, w_area 0.1), adv_scale 1:
R_i = 0.1*loc_i + 0.1 = [0.19,0.18,0.17,0.15,0.14,0.13,0.12,0.11]; mean = 0.14875 (derive: sum 1.19/8).
A_i = (R_i - (S - R_i)/7) = 8/7 (R_i - mean): A_0 = 8/7 * 0.04125 = 0.047143; A_7 = 8/7 * (0.11-0.14875) = -0.044286; sum A = 0 (<1e-9); strictly ordered like loc.
**98. T7_zero_sum** (property; hypothesis)
Random terms (iou in [0,1], f5, loc, legal in {0,1}, g_area in {0,1}, topo values, sig_id ints); B in 1..8 groups: per-group sum A within 1e-9, including groups with zeroed A (spread rule), groups with illegal answers, c_topo large.
**99. T8_T13_scaling** (property)
All reward terms scaled by k in {0.01, 0.5, 3, 100}: not literally possible for legal and indicator terms; use scaling of the weights: `RewardV5(alpha_loc ...)` scale (alpha, w_f5, w_area, c_topo) and iou term... Spec says "all terms scaled": test on `R` directly by passing terms with iou_tc, f5, g_area scaled x k, with `legal`=1 and c_topo=0 and loc=0, area_gate False: A(k) == k A(1) within 1e-12 relative. Flat-group threshold interacts (spread < 0.01): choose spreads >= 1 after scaling so zero-rule not triggered (and test that tiny k does trigger zero-rule: documented). T13: `adv_scale = s` in {0.1, 2, 17.5} gives A = s * A(1) within 1e-12 for random groups. No std division: group with spread 0.02 gives A of order 0.02 not O(1) (assert max|A| < 0.05).
**100. T9_identical_group** (unit)
8 identical answers (same everything, same sig_id) -> A all exactly 0.0, no NaN/inf (np.isfinite). 8 identical R but different sig_id -> not flat by the sig rule: A from c_topo * r_topo differences may be non-zero; assert zero only if r_topo equal (derive). Spread = 0.0099 (<0.01) with one sig -> zero; spread = 0.0101 -> not zero (boundary).
**101. T11_hit_beats_miss** (property)
Pair of answers equal in f5, g_area, legal, topo terms; one with iou >= 0.1 and the other iou < 0.1 with the best loc 1 -> R_hit >= R_miss for every combination of hit iou in [0.1, 1] and loc in [0, 1] under default spec (equality allowed at iou=0.1 exactly with loc 1: R_hit = 0.1 vs R_miss = 0 + 0.1*1 = 0.1). With a non-default `alpha_loc` larger than miss_iou -> constructor already rejects.
**102. T12_closest_miss_largest** (unit)
Group of 8 misses at distances (px): [5, 5, 10, 10, 20, 20, 30, 40] with r=7 -> radius2 = 49 (derive; union area/pi): loc = exp(-d^2/98) = [0.7748, 0.7748, 0.3604, 0.3604, 0.0169, 0.0169, ~1e-4, ~0] (derive in test). Expect: argmax A is a d=5 answer, A_i > 0 for the closest, A strictly decreasing in d (pairs tie). Variant: some answers in group have g_area 0 -> loc 0, ranks below.
Spearman(A, -d/r) over this group > 0.9 (derive).
**103. area_gate_option** (unit)
`area_gate=True`: area term multiplied by min(1, iou/0.1): iou 0 -> area 0; iou 0.05 -> 0.5 * w_area * g_area; iou 0.3 -> full. False: unchanged. Check R values numerically through advantages.
**104. advantages_v5_input_validation** (unit)
Term arrays of different length -> ValueError; length not multiple of group (e.g. 12) -> ValueError; group arg 4 with length 12 passes; empty arrays (length 0) -> empty result (a multiple of group) or documented error; missing key "loc" -> KeyError/ValueError (pin; not silently 0).
**105. advantages_v5_leave_one_out_exactness** (property)
For random groups: A_i == adv_scale * (R_i - mean(R without i)) computed by an explicit O(G^2) loop in the test, tolerance 1e-12; with flat group exceptions; compare to `advantages_v4` NOT equal (different design; guards accidental alias).
**106. legal_gate_threshold** (unit)
legal = 0.5 -> legal (>=), 0.4999 -> illegal; illegal answer contributes R = c_topo * r_topo = 0 exactly (r_topo is 0 for illegal).
**107. ExpertCache_ref_sum_radius2** (unit)
Disks: three raters r=18,19,20 concentric at (29.5, 29.5) crop coords in a 60 x 60 image placed at known offset: vote mask (>=2 of 3) has pixels -> ref_sum = (sum x, sum y, count) of the vote mask (= r=19 disk) equals numpy; raters with no vote majority (three disjoint small disks) -> vote mask empty -> fall back to union mask sums; `radius2 == union_area / pi` (union = r=20 disk area count / pi). Coordinates in source-slice pixel frame (not crop frame): shifting the crop offset changes ref_sum by n*offset; single-pixel union -> radius2 = 1/pi (nonzero). Cache round-trip via pickle (existing tests use pickle) keeps the fields; a cache saved by the old code loads (derived on load) without StaleCacheError beyond the fingerprint rule.
**108. TRLReward_v5_integration** (unit with stub TRL)
`TRLReward(caches, stop_ids, spec=REWARD_V5)` logged terms include `loc`; the aggregator's returned advantages == `advantages_v5(terms)`; after TRL's advantage step the float32 result `np.allclose(.., rtol=1e-6)`; no std normalisation applied (scale_rewards off): check the GRPOConfig used has the matching option.
**109. v4_still_default_until_archive** (regression)
`score_completion` for v4 spec gives the same terms as the golden for all existing fixtures (extra "loc" key permitted but v4 advantages unchanged): existing tests/test_reward.py passes unmodified.

---------------------------------------------------------------------------------------------------------

## K. Run guards [spec 3.7]

**110. BudgetStop_threshold** (unit, fake clock)
rate 0.44 $/h, cap 3.00 -> trips at elapsed 6.8182 h. start_time = t0. Step at t0 + 6.8181 h: no flags; t0 + 6.8182 h (>=, use rate*elapsed exactly equal: choose rate=1, cap=2, elapsed exactly 2.0 h -> trips). Both `should_training_stop` and `should_save` set together. Never before: 10 000 steps with elapsed 0..cap-epsilon never trip. Flags are not reset by later calls. cap 0 -> trips at first step. Negative or zero rate -> ValueError (pin). Clock monotonic injection (time.time vs monotonic: use start_time as passed).
**111. LocGuard_slopes** (unit)
Window 50. Feed 49 logged steps with loc rising and iou falling -> no stop (never with fewer than 50). 50th step -> stop and save. Cases (loc slope, iou slope): (+,-) stop; (+,0) stop (<= 0); (+,+) no stop; (-,-) no stop; (0,-) no stop (loc slope must be strictly > 0); noise: loc = i*1e-6 + noise sd 1e-3 over 50 points -> slope sign from least-squares exactly as numpy.polyfit (compare to polyfit sign, not eyeball). Rolling: stops at step 60 when the last-50 window turns bad although the first 50 were fine. `rewards/loc/mean` never logged (v4 run) -> never stops even over 500 steps; logged only in some steps -> guard counts logged steps of loc, defined in the test (pin). NaN in loc history -> no crash; no stop (pin).
**112. guards_with_trainer_state** (integration with HF TrainerControl)
Callbacks set `control.should_training_stop = True` and `control.should_save = True` through `on_step_end`/`on_log`; stop leads to a checkpoint save (stub trainer counts saves == 1).

---------------------------------------------------------------------------------------------------------

## L. Cross-cutting, failure and performance gates

**113. determinism_end_to_end** (integration)
`build_rl_rows` twice -> identical manifest except timestamps/`versions` ordering (pin the fields); ids_sha256 equal; `train_rl --dry-run` (if exists) consumes rows in file order: first 3 batches' keys equal schedule[0:96].
**114. no_val_select_in_training_rows** (property, real store)
Every training row's group_id maps to `val_half`... train rows come from the TRAIN split, so the check is: no row's group_id appears in the val split at all, and the val-select group set is disjoint from all rows' group_ids (D1 reuse). Also `build_rl_rows` run never opens a val-select image (monkeypatch open/imread to record paths: none under val-select items).
**115. test_view_not_scored_in_phase** (guard)
Grep-style test: no phase-B script other than final_eval references the test view for scoring; `final_eval` with test view requires explicit flag (pin).
**116. masks_stay_on_gpu** (GPU)
During `score_batch`/`TRLReward` call, `torch.cuda` host-copy counter of mask-sized tensors is zero (use a hook on `Tensor.cpu/numpy` with size > 1 MiB; only score arrays returned). Why: spec constraint "masks never leave the GPU".
**117. perf_pool_and_schedule** (performance)
pool() on 7 341 items x draws 8 finishes in < 120 s on 4 procs (budget, to calibrate once); schedule() with 100 steps x 32 on the real pool < 2 s; memory < 2 GB. Gate: the build of 3 200 windows + images with `--procs 8` < 15 min (pin as [estimate], recalibrate).
**118. perf_advantages_v5** (performance)
advantages_v5 on 32 x 8 terms < 1 ms median (numpy vectorised; no Python loop over groups beyond 32) -- fails if a per-answer Python loop sneaks in (> 5 ms). 10 000-group batch < 1 s.
**119. perf_geo_mask_overhead** (GPU)
FastGRPO step time with sequence_geo_mask vs sequence_mask on the same micro-batch: overhead < 2 % (accumulate on device, no sync).
**120. perf_sums_kernel** (GPU)
score_masks on 256 masks of 512 x 512: time with sums <= 1.05 x time of the previous build (golden timing recorded).
**121. fail_fast_before_spend** (failure)
Every ValueError of tests 30, 65, 68 occurs before `AutoModel...from_pretrained` / vLLM engine creation (monkeypatch both to raise AssertionError("loaded")). Why: a bad row file must not cost GPU rent.

---------------------------------------------------------------------------------------------------------

# Pytest sketches (the 5 most important tests)

Importance: (1) advantages_v5 hand numbers + zero sum (reward core), (2) cell edges + largest-remainder schedule counts (data law), (3) sum_x/sum_y CPU vs numpy vs CUDA bit-equality (scorer exactness), (4) is_ratio geo-mask values (number match), (5) train_rl row guards before any model load (spend safety).

```python
# tests/test_phase_b_core.py
import hashlib
import math
import numpy as np
import pytest
from hypothesis import given, settings, strategies as st

# ---------------------------------------------------------------- 1. advantages_v5
def terms(iou, loc, g_area=1.0, legal=1.0, f5=0.0, topo=0.0, sig=0.0):
    n = len(iou)
    full = lambda v: np.full(n, v, np.float64) if np.isscalar(v) else np.asarray(v, np.float64)
    return {"iou_tc": full(iou), "f5": full(f5), "legal": full(legal), "g_area": full(g_area),
            "loc": full(loc), "topo_any": full(topo), "topo_freq": full(topo), "sig_id": full(sig)}

def test_advantages_v5_hand_numbers_and_zero_sum():
    from objdel.reward import REWARD_V5, advantages_v5
    loc = [0.9, 0.8, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1]
    A = advantages_v5(terms([0.0] * 8, loc))
    R = 0.1 * np.array(loc) + 0.1                   # alpha_loc * loc + w_area * g_area
    expect = REWARD_V5.adv_scale * (R - (R.sum() - R) / 7)
    np.testing.assert_allclose(A, expect, atol=1e-12)
    assert abs(A.sum()) < 1e-9 and np.all(np.diff(A) < 0)      # T7 and T12 ordering
    np.testing.assert_allclose(A[0], REWARD_V5.adv_scale * 8 / 7 * (0.19 - 0.14875), atol=1e-12)

@settings(max_examples=200, deadline=None)
@given(st.lists(st.tuples(st.floats(0, 1), st.floats(0, 1), st.booleans()), min_size=8, max_size=64)
       .filter(lambda x: len(x) % 8 == 0), st.floats(0.01, 20))
def test_advantages_v5_zero_sum_and_scale(rows, s):
    import dataclasses
    from objdel.reward import REWARD_V5, advantages_v5
    iou, loc, legal = map(np.array, zip(*rows))
    t = terms(iou, loc, legal=legal.astype(float))
    a1 = advantages_v5(t, spec=dataclasses.replace(REWARD_V5, adv_scale=1.0))
    as_ = advantages_v5(t, spec=dataclasses.replace(REWARD_V5, adv_scale=s))
    np.testing.assert_allclose(as_, s * a1, rtol=1e-9, atol=1e-12)          # T13
    assert np.all(np.abs(a1.reshape(-1, 8).sum(1)) < 1e-9)                    # T7

def test_identical_group_is_zero_not_nan():
    from objdel.reward import advantages_v5
    A = advantages_v5(terms([0.4] * 8, [0.0] * 8, sig=3.0))
    assert np.all(A == 0.0) and np.all(np.isfinite(A))                       # T9

# ---------------------------------------------------------------- 2. cell edges + schedule counts
@pytest.mark.parametrize("extent,off,want", [(9, .149, (0, 0)), (10, .15, (1, 1)), (19, .1999, (1, 1)),
                                             (20, .2, (2, 2)), (39, .2999, (2, 2)), (40, .3, (3, 3)),
                                             (0, 0.0, (0, 0)), (10**6, 10.0, (3, 3))])
def test_cell_edges_go_to_upper_bin(extent, off, want):
    from objdel.data.windows import cell
    assert cell(extent, off) == want

def largest_remainder(weights, n):
    raw = np.asarray(weights, float) * n
    out = np.floor(raw).astype(int)
    for i in np.argsort(-(raw - out), kind="stable")[: n - out.sum()]:
        out[i] += 1
    return out

@pytest.mark.parametrize("steps,per_step", [(10, 32), (7, 32), (1, 32), (3, 5)])
def test_schedule_kind_and_cell_counts(steps, per_step, big_pool, uniform_target):   # fixtures: pool >> need
    from collections import Counter
    from objdel.data.windows import KINDS, schedule
    rows = schedule(big_pool.windows, uniform_target, steps, per_step, seed=0)
    assert len(rows) == steps * per_step
    kinds = Counter(w.kind for w in rows)
    shares = (0.15, 0.75, 0.10)
    assert [kinds[k] for k in KINDS] == list(largest_remainder(shares, steps * per_step))
    for k in ("shifted", "rescaled"):
        sel = [w for w in rows if w.kind == k]
        got = np.zeros(16, int)
        for w in sel:
            got[w.cell[0] * 4 + w.cell[1]] += 1
        np.testing.assert_array_equal(got, largest_remainder(uniform_target.ravel(), len(sel)))
    for s in range(steps):
        ids = [w.item_id for w in rows[s * per_step:(s + 1) * per_step]]
        assert len(ids) == len(set(ids))
    assert len({w.window_id for w in rows}) == len(rows)

def test_schedule_unfillable_names_kind_and_cell(tiny_pool_without_shifted_cell_23, uniform_target):
    from objdel.data.windows import schedule
    with pytest.raises(ValueError, match=r"shifted.*(\(2, 3\)|23)"):
        schedule(tiny_pool_without_shifted_cell_23.windows, uniform_target, 2, 32)

# ---------------------------------------------------------------- 3. sum_x / sum_y exactness
def _masks(rng):
    yield np.zeros((64, 64), bool)
    yield np.ones((64, 64), bool)
    one = np.zeros((64, 64), bool); one[63, 63] = True; yield one
    for _ in range(200):
        h, w = rng.integers(8, 513, 2)
        yield rng.random((h, w)) < rng.random()

def test_sums_match_numpy_cpu():
    from objdel import _gpucheck as G
    rng = np.random.default_rng(0)
    for m in _masks(rng):
        r = G.score_masks(m[None])           # call shape per existing tests/test_cuda_checker.py (adapt)
        ys, xs = np.nonzero(m)
        assert r["sum_x"].dtype == np.int64 and r["sum_y"].dtype == np.int64
        assert int(r["sum_x"][0]) == int(xs.sum()) and int(r["sum_y"][0]) == int(ys.sum())

@pytest.mark.gpu
def test_sums_cpu_gpu_bit_equal(cpu_build, cuda_build):    # fixtures import the two builds
    rng = np.random.default_rng(1)
    for m in _masks(rng):
        a, b = cpu_build.score_masks(m[None]), cuda_build.score_masks(m[None])
        assert a.keys() == b.keys()
        for k in a:
            np.testing.assert_array_equal(a[k], b[k], err_msg=k)

# ---------------------------------------------------------------- 4. is_ratio geo mask
def test_is_ratio_geo_mask_values():
    import torch
    from objdel.train.fast import is_ratio
    d = torch.tensor([[0.005] * 4, [0.02] * 4, [1.0, -1.0, 0.0, 0.0], [-0.02] * 4, [0.005] * 3 + [5.0]],
                     dtype=torch.float64)
    mask = torch.ones_like(d); mask[4, 3] = 0                      # padded token must be ignored
    logps, sampling = d, torch.zeros_like(d)
    r = is_ratio(logps, sampling, mask, "sequence_geo_mask", None, 2.0, geo=(0.99, 1.01))
    assert torch.all(r[1] == 0) and torch.all(r[3] == 0)           # outside [a, b]
    torch.testing.assert_close(r[0], torch.full((4,), math.exp(0.005), dtype=torch.float64))
    torch.testing.assert_close(r[2], torch.tensor([2.0, math.exp(-1), 1.0, 1.0], dtype=torch.float64))
    torch.testing.assert_close(r[4, :3], torch.full((3,), math.exp(0.005), dtype=torch.float64))

def test_is_ratio_old_modes_unchanged(old_is_ratio_oracle):        # copy of pre-change function in tests/oracles
    import torch
    from objdel.train.fast import is_ratio
    g = torch.Generator().manual_seed(0)
    for mode, lo, hi in [("sequence_mask", None, 3.0), ("token_truncate", None, 2.0)]:
        lp, sp = torch.randn(8, 50, generator=g) * .1, torch.randn(8, 50, generator=g) * .1
        m = (torch.rand(8, 50, generator=g) > .2).float()
        assert torch.equal(is_ratio(lp, sp, m, mode, lo, hi), old_is_ratio_oracle(lp, sp, m, mode, lo, hi))

# ---------------------------------------------------------------- 5. train_rl guards before any model load
@pytest.fixture
def no_model_load(monkeypatch):
    def boom(*a, **k): raise AssertionError("model/engine loaded before the guards ran")
    monkeypatch.setattr("transformers.AutoModelForImageTextToText.from_pretrained", boom, raising=False)
    monkeypatch.setattr("vllm.LLM", boom, raising=False)

def _root(tmp_path, ids, steps, per_step, sha=None):
    import json
    (tmp_path / "lidc").mkdir()
    (tmp_path / "lidc/train.jsonl").write_text("\n".join(json.dumps({"item_id": i}) for i in ids))
    sha = sha or hashlib.sha256("\n".join(ids).encode()).hexdigest()
    (tmp_path / "manifest.json").write_text(json.dumps({"schedule": {"steps": steps, "per_step": per_step,
                                                                       "ids_sha256": sha}}))
    return tmp_path

@pytest.mark.parametrize("n_ids,tamper,kw", [
    (95, False, {}),                                   # row count != steps*per_step
    (96, True, {}),                                    # ids hash mismatch
    (96, False, {"generation_batch_size": 16 * 8}),    # gbs//G = 16 != 32
    (96, False, {"num_iterations": 2}),
    (96, False, {"steps_per_generation": 3}),          # != gradient_accumulation_steps
    (96, False, {"steps": 4}),                         # --steps > schedule.steps (3)
])
def test_train_rl_guards_raise_before_model_load(tmp_path, no_model_load, n_ids, tamper, kw):
    from experiments.phase2 import train_rl
    ids = [f"w{i}" for i in range(n_ids)]
    root = _root(tmp_path, ids, 3, 32, sha="0" * 64 if tamper else None)
    with pytest.raises(ValueError):
        train_rl.main(["--train", str(root), "--steps", str(kw.pop("steps", 3)), *sum(
            ([f"--{k.replace('_', '-')}", str(v)] for k, v in kw.items()), [])])
```

Notes for the implementer of the sketches: fixtures `big_pool`, `uniform_target`, `tiny_pool_without_shifted_cell_23`, `cpu_build`, `cuda_build`, `old_is_ratio_oracle` must be added to tests/conftest.py or tests/oracles; the `train_rl.main(argv)` signature and the `score_masks` call shape are assumed and must be adapted to the real entry points; the spec has no unit named for them.

---------------------------------------------------------------------------------------------------------

# Spec ambiguities found while writing tests (resolve before coding)
1. D12 wording: "more than `limit - 100` tokens fails" (3.2) vs "max tokens + 1 below `MAX_COMPLETION - 100`" (5.1) differ by one at the boundary (tests 49).
2. `cell` for float32 offsets exactly at an edge (test 6): dtype is not specified.
3. `alpha_loc == miss_iou` (the default) makes R flat (slope 0) under a miss when loc = 1; T1 holds only as non-strict monotone (test 91).
4. T8 "all terms scaled" cannot scale the `legal` indicator or `loc` in [0,1]; the test scales iou_tc, f5, g_area, weights and sets c_topo = 0 (test 99). The drop_spread rule (0.01) breaks proportionality for tiny k.
5. `schedule` with `steps=0`, `cell_shares([])`, `radius2 = 0`, rater masks empty: behaviour unspecified (pinned as test decisions).
6. Padded windows: D14 reports "padded share" but `pool` rules say only "doctor pixels inside with >= 2 px clearance"; whether a window may extend beyond the slice is unspecified (test 15).
7. `LocGuard` with sparse logging of `rewards/loc/mean`: the window counts logged steps; unspecified for a metric logged on a subset of steps (test 111).
8. T12 uses Direction.md s5 which test authors must not read; the test uses distances (5,10,20,30,40) px, r = 7 stated in the spec.
