# Review of phases/B.md (Phase B: quality)

Reviewer: fresh critical reviewer, 2026-11-08. I checked every code claim in B.md against REPO. Tags:
[code: file:line], [measured: file], [web: url], [estimate].

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | Checked B.md's two catalog picks (DICOM preflight, compute-sanitizer). They are fine as recommendations. The real gap is that B1.7 has no local DICOM to run on (item 9). |
| engineering:code-review | Line-by-line check of the files B.md names: fast.py, configs.py, dataset.py, serve.py, build_eval_views.py, gpucheck.cpp, advantages.py, spec.py, TRL 1.17.1. Found items 1-4, 6-8 and 12. |
| engineering:debug | Re-traced the evidence for "half the signal is lost" and "a miss gets no gradient". Neither holds as stated (items 13 and 14). |
| engineering:system-design | Traced the data flow from schedule to TRL sampler to reward. Found the row-drop break (item 1) and the batch-shape check gap (item 6). |
| engineering:architecture | Checked the ADR options. Found two confounds in the A/B design (items 5 and 11) and the circular coverage gates (item 10). |
| engineering:documentation | Flagged vague steps: the guard callback, test reads and freeze contents (items 21, 23 and 4). |
| huggingface-skills:huggingface-best | Keeping Qwen3-VL-2B for Phase B is fine. No leaderboard query needed. Web check on fp16 for Qwen-VL vision towers (item 15). |
| engineering:testing-strategy | Found a metric with a composition bias (item 16), missing paired tests (item 17) and missing assertion tests (items 1 and 15). |
| engineering:tech-debt | Found a duplicate of the existing coverage and audit modules (item 7). Found that the old path is deleted earlier than the v4 rule allows (item 20). |
| Web (Exa) | Confirmed arXiv 2510.26788: it also tests LoRA. Confirmed the SkyRL defaults. SkyRL recommends a geometric sequence mask first (item 18). Found the vLLM issue on Qwen2.5-VL fp16 vision overflow (item 15). |

## Findings (most severe first)

1. **blocking. B1.5, "rows = dataset.grpo_rows(...)": `prompts()` drops rows, so step boundaries shift.**
   Problem: `prompts()` keeps only `r["rl"]` records [code: src/objdel/train/dataset.py:221]. `rl` is
   `len(raters) >= 2` [code: src/objdel/reward/cache.py:103-105]. The pool in B1.1 never filters on it. Every
   single-rater window in the schedule is silently dropped. Every later row then moves into another step. That
   breaks "step s gets rows 32(s-1)..32s-1", "one item per step" and the kind shares per step.
   Fix: in `windows.pool`, skip items with `not rl_eligible(item.raters)` and count them. In `train_rl`, after
   `grpo_rows`, assert `len(rows) == schedule.steps * 32` and that `[r["key"] for r in rows]` equals the ordered
   window ids (compare their SHA-256 with `ids_sha256`). Add both checks to `tests/test_p2_rows.py`.

2. **blocking. B1.1 step 2 and 3: the shifted draw must pass `side=S0`, and the rescaled fallback is not a "fallback".**
   Problem: `augment.draw(item, rng)` with no `side` draws S from `sides(item)`, which is S0 plus its neighbours
   [code: src/objdel/targets/augment.py:85-94]. B1.md F4 notes this gives a 50/50 side mix
   [measured: notes/B1.md:101]. The exam view calls `augment.draw(bare, rng, side=S0)`
   [code: experiments/phase2/build_eval_views.py:60]. If B1.1 is coded as written, D10 (S0 share >= 0.75) fails and
   draw 0 does not equal the built views. Also, `_record` turns a failed `neighbour_side` into a shifted S0 window
   with `kept_side=True`. It is not a fallback. B1.1 step 3 ("no admissible draw is not added") does not say which
   case it means.
   Fix: write the exact calls. Shifted: `augment.draw(bare, rng, side=S0)`. Rescaled: `augment.neighbour_side(bare,
   S0, rng)`; on None, apply the B1.1 rule to `kept_side` windows explicitly. Either drop them, or keep them,
   label them `kind="shifted"` and count them under D6. Add a test that S0 windows are >= 0.75 of off-centre rows
   on the fixture.

3. **blocking. B1.6 vs B3.3: the frozen test-view hash breaks at the planned fingerprint change, and also across machines.**
   Problem: the view manifest stores `reward_fingerprint`, the Python version and the package versions
   [code: experiments/phase2/build_eval_views.py:89-92]. The jsonl stores cache digests. B3.3 changes the
   fingerprint, so after B3.3 `final_eval` refuses the frozen test view. "R5's GPU-box view must give the same
   hash" also fails for harmless reasons (another Python or OpenCV build).
   Fix: freeze a content hash, not the manifest SHA. Use SHA-256 of the ordered `(item_id, view.offset, view.side,
   image_sha256)` tuples plus `target_set_sha256`. Do the fingerprint split (B3.3) as task B1.0, before B1.4 and
   B1.6. Then the caches, views and the B1 rows root are built and uploaded once, not twice.

4. **major. B1.6 / D15: the freeze contents and their storage place are left open ("constant ... or a small JSON").**
   Fix: choose one place, `docs/reference/frozen_views.json`, with the fields from item 3 plus the `build_eval_views`
   arguments (seed, kinds, store manifest SHA). `final_eval` loads it, recomputes the hash and stops on a mismatch.

5. **major. B3.8 + B3.9: the Neyman tilt down-weights the cells that B3 targets, cannot pass gate D7, and confounds B3.**
   Problem: (a) Under v4, groups where all 8 answers miss have near-zero spread. So `cell/<ij>/group_sd` is lowest
   in the miss-heavy off-centre cells. `tilt = sd / mean(sd)` then gives those cells 0.5, which is the opposite of
   "tilt for hard cells". (b) The sd is measured on v4's reward input but used for a v5 run. (c) A tilted target
   with factors 0.5-2.0 cannot meet D7 (worst ratio >= 0.8) or D8 (TV <= 0.10) against the untilted val target, so
   the rebuild writes no manifest. (d) B3.9 says the A/B is clean, but B3 then changes both the reward and the data
   mix.
   Fix: remove B3.8 from Phase B. If it stays, (i) measure sd on a reward-independent signal (sd of IoU_tc), or on
   the shadow v5 R logged during B1. (ii) Run D7-D9 against the tilted target and report the untilted values.
   (iii) Run B3 with B1's untilted rows so that only the reward changes.

6. **major. B1.5 batch-shape check: the check tests the wrong quantity.**
   Problem: TRL's sampler takes `batch_size = generation_batch_size // num_generations` and `repeat_count =
   num_iterations * steps_per_generation` [code: .venv trl/trainer/grpo_trainer.py:1372-1377]. "32 prompts per
   optimizer step" holds only when `steps_per_generation == gradient_accumulation_steps` and `num_iterations == 1`.
   Fix: check `args.generation_batch_size // args.num_generations == schedule.per_step`,
   `args.steps_per_generation == args.gradient_accumulation_steps` and `args.num_iterations == 1`. Raise
   `ValueError` otherwise.

7. **major. B1.3: new `data/audit.py` and `coverage()` duplicate existing modules.**
   Problem: `src/objdel/eval/coverage.py` already holds "train-vs-exam coverage" (props, distribution, gaps).
   `src/objdel/eval/audit.py` is the hack audit and `src/objdel/eval/gates.py` holds `Gate`-style entry gates. B.md
   names none of them. Two coverage modules and two `audit.py` files break the owner rule (archive old code, no
   confusion).
   Fix: name the overlap in B1.3. Either move `eval/coverage.py` into `data/coverage.py` and extend it with the
   cells, TV and KS, or archive it with a README row. Name the new gate module `data/quality.py` (not `audit.py`).
   Reuse the gate record type from `eval/gates.py` if it fits.

8. **major. B1.4 step 3: the key and folder rule for windows is not stated.**
   Problem: `_prepare` writes `images/<item_id>.png` and `caches/<item_id>/` [code: dataset.py:102,110]. Many
   windows share one source item. If `window_record` passes the source `item_id`, windows overwrite each other's
   images and caches. The proposed `window_id = f"{item_id}/{kind}{draw}"` puts "/" in the key, which makes nested
   folders and a sub-folder inside the source item's cache folder.
   Fix: `window_record` must pass `window_id` as `_prepare`'s item id. Use a separator without "/" (for example
   `f"{item_id}~{kind}{draw}"`). Test that every window's PNG and cache path is unique.

9. **major. B1.7 cannot run as written: no LIDC DICOM on the Mac.**
   Problem: `data/raw` holds only `lidc_xml` and `riga` [measured: `ls data/raw`]. The `lidc_slices` JSONs keep no
   per-SOP z values. The script "read data/raw DICOM headers" has no input.
   Fix: first do a $0 check from the LIDC XML: `lidc.sop_z` already reads `imageSOP_UID` and `imageZposition`.
   List z values shared by different SOP UIDs in the 4 series. Download the 4 series headers from TCIA only if the
   XML shows duplicates. Write the download into the task, with the folder (`data/raw/lidc_dicom/<series>`) and
   `python -I`.

10. **major. D7-D9 are circular, and val both shapes and judges the run.**
    Problem: the scheduler draws cells from the val target, so D7/D8 pass by construction. They catch only
    implementation bugs and cell exhaustion. Val then judges gate B1. Section 7 adds a val split only after a
    val/test disagreement, which is too late to be useful. SD3 shows that train-shifted alone has TV 0.114 to test,
    and val 0.074 [measured: notes/SD3_data.md:282]. So the val target corrects real size differences, and the
    coupling is real.
    Fix: split val by patient now ($0). Use val-shape for `target` and val-select for gate B1 and B3 values. Label
    D7-D9 "implementation checks". Add one independent report: coverage of the schedule against the train-split
    shifted law (the same generator as the exam).

11. **major. B2.4 rule 3 and the B1 run: the numerics change confounds the data result, and B1 is blocked on B2.**
    Problem: R5 used `sequence_mask`. The B1 run uses the B2 winner, so gate B1 ("B1 - R5") measures data plus
    numerics. "Else stop and report (do not run B1 on sequence_mask)" puts the owner's top priority (data) behind
    a numerics experiment.
    Fix: if no arm passes gate B2, run B1 on `bf16-seqmask` (the control) and still log the B2.1 stats. If a
    winner is used, write in ADR-0014/0012 that the B1-vs-R5 gain includes the numerics change. The B2 probe gives
    the per-step reward delta that bounds it.

12. **major. B3.4 / B3.9: RLOO without std division changes the step size, so B3 vs B1 is not like for like.**
    Problem: v4 advantages are z-scores with weights summing to 1 plus c_topo = 2 [code: advantages.py:95-107], so
    their typical size is order 1. v5's raw advantages for IoU differences are about 0.05-0.3 (B3.4's own example:
    +0.057). With the same lr 1e-5 and the dr_grpo loss linear in A, B3 takes steps several times smaller. A B3
    failure then may come from the step size, not the reward.
    Fix: in the B3.7 lab round, measure the median |A| for v4 and v5 on non-dropped groups. Scale v5's A by that
    ratio (one constant in RewardV5, recorded in ADR-0013), or scale the lr. Add a test that the constant is applied.

13. **major. Section 1 item 5: "about half the signal is thrown away" is not supported by the cited numbers.**
    Problem: `importance_sampling_ratio/min = 0` in a step means that at least 1 of 256 sequences was masked. A mean
    ratio of 0.82 points to about 18 % masked if kept ratios are near 1 [estimate]. Neither number gives "half".
    Fix: write "an unmeasured share (at least 1 sequence in 98.1 % of steps; about 18 % if kept ratios are near 1
    [estimate])". B2.1 then measures it. Also add a $0-cost attribution to the B2 probe. At step 1 the weights are
    the same in trainer and vLLM, so log `sampling_logp_difference/mean` for (a) FastGRPO default, (b)
    `--shared-prefix` off and (c) merged LoRA sync on. The gap can come from the packed forward or from bf16
    rounding of the merged weights (`lora_merged_params`), not only from the dtype. Today B2 assumes the dtype.

14. **major. Section 1 item 6 and gate B3 "pre-run": v4 does give gradient to all-miss groups.**
    Problem: v4 z-scores g_area and F5 as well as IoU. Each sigma has a floor of 0.05 and is not dropped. A group
    is zeroed only when every term spans < 0.01 [code: advantages.py:100-106]. In an all-miss group, g_area still
    varies (doctor-sized versus not), so A is not zero. The pre-run gate ("higher share of all-miss groups with
    non-zero advantage under v5") can then fail or pass for reasons that have nothing to do with location.
    Fix: change the B3.7 metric and gate. In all-miss groups, use Spearman(A, -centroid distance / radius) under v4
    and under v5. Gate: v5 > 0.3 and v5 > v4 [estimate thresholds]. Also report the share of A variance from the
    `w_area * g_area` term against the `alpha_loc * loc` term (see item 19).

15. **major. B2.2 fp16 arm: known Qwen-VL fp16 vision overflow, a fallback that cannot be built, and no proof that the patch took effect.**
    Problem: vLLM has fixed fp16 overflow in the Qwen2.5-VL vision encoder before (PR 13968, issue 13817)
    [web: https://github.com/vllm-project/vllm/pull/13968]. Your pinned vLLM may or may not have a matching guard for
    Qwen3-VL. The fallback "vision bf16, language fp16" has no `LLM(dtype=...)` setting: vLLM takes one dtype per
    model. The wrapper on `trl.generation.vllm_generation.LLM` works only because TRL imports `LLM` at module level
    [code: .venv trl/generation/vllm_generation.py:43,356], and nothing checks that it was applied.
    Fix: (i) add a $0 Mac check before the rental: run one fp16 forward of the vision tower on 32 train images (CPU
    or MPS) and count inf/NaN. (ii) Drop the mixed-dtype fallback; on NaN/inf the arm fails and bf16-tis is used.
    (iii) After `FastGRPO.__init__`, assert `self.vllm_generation.llm.llm_engine.model_config.dtype == torch.float16`
    and test the assertion.

16. **major. Gate B3 "val IoU on found images falls by <= 0.01": composition bias.**
    Problem: when v5 turns misses into hits, the new hits are hard images with lower IoU. They join the "found" set
    and pull its mean down even if no image got worse. The gate then punishes success.
    Fix: compare IoU on images found by both runs (the paired set), with the nodule bootstrap CI. Gate: lower bound
    of (B3 - B1) > -0.01.

17. **minor. Gates B1/B3: CI rules are incomplete and unpaired.**
    Fix: state that all run comparisons use a paired nodule bootstrap on the same val items (R5 answers come from
    `runs/r5/eval/eval_val/items.jsonl`, 1,084 items [measured]). Give B3's miss-share gate a CI rule too (upper
    bound of the relative change <= -25 %, or point value plus lower bound < 0). For B2, compare arms by per-step
    paired differences on the same 640 rows, not by the control's sd.

18. **minor. B2 options: a length-fair mask is missing, and FastGRPO computes the ratio itself.**
    Problem: SkyRL's first recommendation is a geometric-mean sequence mask (ratio of the mean token gap, in
    [0.99, 1.01]) plus token TIS [web: https://docs.skyrl.ai/docs/algorithms/off_policy_correction]. TRL 1.17.1 has
    no such mode, but `FastGRPO.set_is_ratio` already computes the ratio through `fast.is_ratio`
    [code: src/objdel/train/fast.py:351-363, 581-591]. Adding one mode there is about 3 lines and fixes the
    long-answer bias that section 1 item 5 describes.
    Fix: put `sequence_geo_mask` (exp of the masked mean gap) into `is_ratio` as an option in the ADR-0012 table.
    Either add it as a 4th probe arm (+20 steps, about $0.08 [estimate: 20 x 23.65 s at $0.641/h]) or record why
    it is left out.

19. **minor. B3.3 "alpha_loc <= miss_iou keeps every hit above every miss" is false once g_area and F5 count.**
    Problem: a miss with g_area = 1 and full loc scores 0.1 + 0.1 x ramp, about 0.2. A hit at IoU 0.1 with
    g_area = 0 scores 0.1. T11 hides this with "all other terms equal". `w_area * g_area` pays for size alone, in
    any position.
    Fix: reword the claim to "for equal f5 and g_area". In the lab round, measure how often a miss outscores a hit
    in the same group. If it happens > 1 % of the time, gate `w_area * g_area` by `min(1, iou_tc / miss_iou)`.

20. **minor. B1.5: the old `--train-views` path is deleted before gate B1, but v4 is kept until gate B3.**
    Fix: use one rule for both. Delete `grpo_rows_views`, `VIEW_SEP`/`_tagged`
    [code: src/objdel/train/serve.py:108-123] and `--train-views` [code: experiments/phase2/train_rl.py:26,47] in
    the commit after gate B1 passes. Also drop the "without schedule: oversample 4 shuffled" branch after gate B1,
    so that only one training path is left.

21. **minor. Section 7, "stop (callbacks guard)": no task builds the guard.**
    Fix: add task B3.5b. Write a `callbacks.LocGuard(window=50)` that stops training when the 50-step slope of
    `rewards/loc/mean` is > 0 and the slope of `rewards/iou_tc/mean` is <= 0. Add a unit test with a synthetic
    log.

22. **minor. B1.2 step 3: substituting the other kind breaks D11's "exact counts".**
    Fix: give D11 a tolerance (substitutions <= 1 % of slots, reported), or make cells by kind fail loudly
    (`ValueError`, then a larger K). Also allocate cell counts by largest remainder, not by random draws. Then
    D7/D8 depend only on the pool, not on sampling noise.

23. **minor. Section 8: the cost inputs are optimistic, and the test view is read twice.**
    Problem: 23.65 s/step comes from a short R5 test on centred rows [measured: runs/r5/shared_rl.json]. The full run
    averaged 27.5 s against an earlier 36.3 s estimate [measured: docs/runbooks/gpu-run.md:136-137]. Shifted rows
    may raise the decode tail, because the step time follows the longest answer. Uploading the 1.2-1.6 GB rows
    root to the host is not costed. S1 and S2 both run the test eval, but section 7 says "test is read once".
    Fix: cost at 27.5 s/step as the base case. Upload the rows root to the private HF repo before the rental. Add a
    dollar stop in the run (callback on elapsed time x rate). Run the test eval only for the final candidate. This
    saves about 10.6 min per session [measured: 634 s].

24. **minor. B1.4 materialization uses the numpy polygonizer: conflict with the owner rule "CUDA for the polygonizer".**
    Problem: section 3.3 says why (the Mac has no NVIDIA GPU, and Phase D comes later), but the owner rule is
    binding, so the exception needs his sign-off.
    Fix: record the exception and its end date (Phase D K1 parity) in ADR-0014. Add a line to the owner-decision
    list in section 8.

25. **minor. B3.2: the location kernel may reach too short a distance for real misses.**
    Problem: `exp(-d^2 / (2 |U|/pi))` is almost 0 beyond about 3 radii. Nobody has measured how far typical misses
    are from the nodule. If most are farther, all-miss groups still get no signal. The example numbers come from
    Direction.md, not from a measurement.
    Fix: in B3.7, report the distribution of d / r for misses on val-shape. If the median is > 3, pick the width
    from that data (for example sigma = max(r, 0.25 S)) or a heavy-tailed kernel `1 / (1 + d^2 / r^2)`. Record the
    choice in ADR-0013. Option to check: the scorer already sums the exact distance-to-union of the spill pixels
    (ScoreRows `spill`, `spill_n`) [code: cuda/gpucheck.cpp:1428-1440]. If the EDT covers the whole frame, the mean
    distance to the union is available with no new row partial.

26. **minor. Gate D9 and section 1: the baseline uses another reference.**
    Problem: D9 "today 0.059" is train-shifted against test, but the gate is against the val target. 0.06 also sits
    right at today's value.
    Fix: give the today value against val-shape for all of D7-D9.

27. **minor. Gate B1 "val IoU >= 0.66" is not linked to the miss gate.**
    Problem: offset > 0.3 has 406 of 1,084 val images. Cutting misses there from 27.3 % to 20 % turns about 30
    misses into hits. That gives about +0.017 val IoU [estimate: 30 x 0.6 / 1,084]. The +0.042 to 0.66 needs gains
    beyond the miss gate, and nothing says where they come from.
    Fix: add the decomposition to section 6.2, or use the CI rule alone (lower bound > 0) and report 0.66 as a
    target, not a gate.

## Checked and correct (no action)

TRL `shuffle_dataset` is public, and `RepeatSampler(shuffle=...)` [code: grpo_config.py:509, grpo_trainer.py:1377].
`logprobs_mode="processed_logprobs"` and no dtype in `LLM(...)` [code: vllm_generation.py:356-372]. `_bf16` refuses
fp16 together with bf16 [code: configs.py:140-142]. Line refs for `is_ratio`, `set_is_ratio`, `log` and
`drop_lora_input_cast` are correct. `kPart = 5`, LSB = leftmost, and the ScoreRows layout match gpucheck.cpp:1144,
1187 and 1409-1441. The popcount column-sum identity is exact. RLOO gives a group-mean of 0 exactly. Cost arithmetic
in 8.1 is right for its inputs (189 min, $2.02). The SkyRL `token_tis_ratio_clip_high: 2.0` claim is confirmed.
arXiv 2510.26788 covers LoRA (B.md says "no VLM", which is also true).
