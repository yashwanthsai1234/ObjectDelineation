# Final review, round 1: RL for LLM/VLM (B2, B3, C, model choice)

Reviewer role: RL research scientist (GRPO, IS correction, reward design, TRL/vLLM internals).
Target: `REPO/docs/implementation2.0/Implementation.md` (11,117 lines). Earlier final reviews: none (`final_r*_*.md`
absent). Items already closed in the phase review logs (B.R 1-27, C.R 1-24) are not repeated.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No new NVIDIA skill fits RL numerics. The doc's single recommendation (`compute-sanitizer-diagnose-cuda-bugs`) stays. No item raised from it. |
| engineering:code-review | Every claim on B2/B3/C was checked against TRL 1.17.1 source in `.venv`, `src/objdel/train/{fast,configs}.py`, `src/objdel/reward/{spec,advantages}.py`. Items 2, 5, 6, 10 come from this. |
| engineering:debug | I re-split `runs/rl/train_rl/rl.json` by step range. This showed the mismatch grows with the adapter (item 1) and that the KL spike did not move the update (item 21). |
| engineering:system-design | I checked how the B2 probe, the B1 run and the B3 run hand settings to each other. Found that KL 0 has no owner (item 4) and that the front gate number differs from Phase C (item 15). |
| engineering:architecture | ADR-0012 options: added TRL's public OPSM (item 6) and vLLM LoRA serving (item 10) as options, with the facts that favour each. |
| engineering:documentation | Each item gives the section, the problem and the exact text or code change. |
| huggingface-skills:huggingface-best | HF API checks of the model-table facts: parameter counts and configs of the 2B, 4B, Qwen3.5 and MedMO models (items 3, 14, 17, 18, 20). |
| huggingface-skills:trl-training | Read TRL 1.17.1 `grpo_trainer.py` 2750-2830 (IS modes), 3066-3089 (OPSM), 3255-3290 (loss, dr_grpo normaliser), 1029-1036 (beta and PEFT ref), `rewards/other_rewards.py:83` (soft overlong punishment). |
| huggingface-skills:huggingface-papers | Checked arXiv 2601.13633 (EGM) tables and 2510.26788 (FP16) on the web. EGM numbers are mis-cited (item 17). The FP16 paper is LLM-only with LoRA tested; the doc states this correctly. |
| Web | SkyRL off-policy page (geo mask [0.99, 1.01], token TIS 2.0, "start with geometric sequence masking"): the doc is correct. EGM HTML Tables 4 and 9. FP16 paper markdown. |

## Measurements made for this review

- Full RL run, `runs/rl/train_rl/rl.json`, 770 steps, by step range [measured now]:

  | steps | median logp diff | median IS ratio mean | share of steps with IS min = 0 |
  |---|---|---|---|
  | 1-20 | 0.011 | 0.998 | 0.25 |
  | 21-100 | 0.028 | 0.827 | 1.00 |
  | 101-300 | 0.024 | 0.771 | 1.00 |
  | 301-500 | 0.019 | 0.810 | 1.00 |
  | 501-770 | 0.016 | 0.845 | 1.00 |

- R5 shared arm, 15 steps from a fresh adapter: logp diff rises step by step, 0.0095 to 0.0134. IS ratio mean stays
  between 0.96 and 1.04 [measured now: runs/r5/shared_rl.json].
- KL max 398,284 at step 746; grad_norm at that step 0.0041; median grad_norm 0.0048, max 0.022; median KL 0.32
  [measured now: rl.json].
- HF API: Qwen3-VL-2B 2.128 B, Qwen3-VL-4B 4.438 B, Qwen3.5-4B 4.660 B (layer types: linear_attention x3, then
  full_attention), MBZUAI/MedMO-4B 4.438 B (Qwen3VLForConditionalGeneration), LiquidAI/LFM3-VL-2B 1.918 B. Vision
  config of 2B and 4B: depth 24, hidden 1024, deepstack [5, 11, 17]; `out_hidden_size` 2048 vs 2560.

## Items

### Blocking

**1. [blocking] B (B2.3, B2.4, B.6.2 Gate B2). The probe measures the mismatch where it does not yet exist.**
Problem: the four arms start from `rl-full-3ep`. This is a merged model folder (`train_rl.py:58-59`
`merge_and_unload`), so `make_trainer` adds a fresh LoRA with B = 0. The full run shows the mismatch appears only
after the adapter grows. In steps 1-20 the median IS ratio mean was 0.998 and only 25 % of steps had a masked
sequence. In steps 21-100 the ratio mean was 0.827 and 100 % of steps had one (table above). In R5 the logp gap rose
every step. A 20-step probe is therefore in the benign regime. Gate B2 (ratio mean >= 0.97, masked < 0.5 %) will
pass for almost any arm, even the control. Rules 2 and 3 of B2.4 then decide on noise, and the probe costs about
84 step-equivalents of a budget that is already over (B.8.1).
Exact fix:
- (a) Replace the 4-arm trained probe with counterfactual logging (item 7) in the B1 run itself. B1 runs on the
  control `bf16-seqmask`, as B2.4 rule 4 already allows. Choose the preset for the B3 run from the B1 run's steps
  50-300, where the adapter has drifted.
- (b) Run B2.3b (no training, item 8) on the B1 run's step-100 and final adapters, not on the step-20 adapter.
- (c) If an fp16 arm is still wanted, run it as a >= 60-step arm. Compare steps 41-60, not 11-20.
- (d) Rewrite Gate B2 to read "over steps 41-60 (or B1 steps 50-300 for counterfactual stats)". Update B.8.1 and
  §6.3: the probe line goes from 84 to 0 steps, or to 60 for the fp16 arm.
- (e) Add one line to B.1 item 5: "The mismatch grows with the adapter: steps 1-20 ratio mean 0.998, steps 21-100
  0.827 [measured: rl.json]."

**2. [blocking] B (B3.4 Advantages v5). The drop rule removes the far-miss groups that the location term is for.**
Problem: v5 writes `A[group] = 0 where ptp(R) < drop_spread` on the weighted reward R. Here `drop_spread` = 0.01
and the location term enters as `alpha_loc * loc` = 0.1 x ramp x loc_raw. Consider an all-miss group whose
loc_raw values span less than 0.1 (every answer more than about 2 radii away: exp(-2) = 0.135, exp(-4.5) = 0.011
at 3 radii) and whose g_area and F5 are equal. It has ptp(R) < 0.01, so it gets A = 0. v4 tests the span of each
unweighted term instead (`np.ptp(R, axis=1).max(-1) < drop_spread` on the raw terms
[code: src/objdel/reward/advantages.py:103]). The v5 rule is therefore stricter than v4 exactly where v5 adds
signal. The text "one sig_id in the group" is also not code.
Exact fix: in B3.4 replace the drop line with
`flat = (np.ptp(np.stack([iou_tc, f5, g_area, loc, r_topo], -1), axis=1).max(-1) < spec.drop_spread) & (sig_id == sig_id[:, :1]).all(1)`.
This tests the span of each unweighted term, as v4 does, with `loc` = ramp x g_area x topo x loc_raw before
`alpha_loc`. Add test T14: an all-miss group with loc_raw = (0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01) and
equal other terms keeps a non-zero A. Add to B3.7: "share of all-miss groups dropped, v4 vs v5". Gate: v5 <= v4.

**3. [blocking] Model choice (§3 "Gate for the 4B arm"; owner decision 12). The 4B gate reads the test set.**
Problem: the gate is "test IoU on shifted views is still below about 0.70 and the miss rate is still above about
6 %". Phase B reads test once, on the final candidate only (B.6.2 "test is not read here"; B1.6 freeze). A
model-selection decision on test leaks test into the choice. The 0.70 also predates the fair ceiling of A2.
Exact fix: in §3 and owner decision 12, write: "run 4B only if, after the 2B B-runs, val-select IoU (paired
nodule-level CI, A1) is below about 0.70 and the val-select miss share is above about 6 %, and the found-image gap to
the A2 fair ceiling is <= 0.02. Test is read once, on the final candidate of whichever model wins."

### Major

**4. [major] C (C8), B (B2.4), C.6 G-C, owner decisions. KL weight 0 has no owner.**
Problem: C8 says "owner decision, run inside Phase B2". B2.4 says "keep beta 1e-3 in Phase B. Phase C owns KL 0".
G-C (<= 18.0 s) needs C8. The owner-decision list has no beta item. Nobody runs it, and G-C can never be measured.
Exact fix: add owner decision 8b: "beta for the B runs: [0.0 from the B1 run]". Evidence for the default: with LoRA,
TRL needs no reference copy at any beta (`grpo_trainer.py:1033-1036`), so beta can be raised again later. The median
KL is 0.32 and the 398 k spike did not move the update (item 21). The reference pass costs 1.4-1.7 s per step. Dr.
GRPO and DAPO run without KL. With beta = 0, item 5's constant goes too. In B2.4 replace the KL bullet with "beta
follows owner decision 8b; ADR-0012 records it". In C8 replace "run inside Phase B2" with "set by owner decision 8b
for every B run". If the owner keeps 1e-3, G-C is reported as a computed number, not gated.

**5. [major] B (B3.4 "Why adv_scale", B3.7 rule 5, ADR-0013). The reason for `adv_scale` is wrong under AdamW.**
Problem: B3.4 says that without `adv_scale` "B3 would take steps several times smaller". AdamW divides by
sqrt(v), so a constant gradient scale leaves the step almost unchanged. C7 of the same document says so: "AdamW is
invariant to a constant gradient scale except through eps and gradient clipping". Clipping is inactive (grad_norm
median 0.0048, max 0.022, against max_grad_norm 1.0 [measured: rl.json]). Per-parameter gradients are about 1e-6
(0.0048 / sqrt(17.4 M LoRA params) [estimate]), far above eps 1e-8. What `adv_scale` really changes is the weight of
the policy term against `beta * KL`, because the KL term is not scaled.
Exact fix: if owner decision 8b is beta = 0, delete `adv_scale`, B3.7 rule 5 and test T13. Write "Adam is
scale-invariant; v5 needs no scale constant". If beta > 0, keep the constant and rename the reason to "keeps v4's
policy-to-KL ratio". Fix the B.7 row "v5 step size differs from v4" the same way.

**6. [major] B (B2.1 `sequence_geo_mask`, B2.2 `bf16-geomask`, ADR-0012). The plan writes its own mask where TRL
1.17.1 has a public one (owner rule: prefer a library).**
Problem: TRL 1.17.1 ships `off_policy_mask_threshold` (OPSM, DeepSeek-V3.2 Eq. 9). It drops a sequence when the
advantage is < 0 and the length-normalised mean of (sampling logp - trainer logp) is > delta. It reads
`sampling_per_token_logps` [code: .venv trl/trainer/grpo_trainer.py:3066-3089, 3190-3203, 3262]. It is length-fair,
it is a public GRPOConfig key, and it works on both TRL's path and FastGRPO's ratio path (it reads `inputs`, not
the old-logp pass). The custom mode works only on FastGRPO's path, adds a private-name risk, and must be kept alive
through Phase C (B.8.4).
Exact fix: in ADR-0012 add option E, "token TIS (C = 2) + OPSM, delta = q99.5 of the B1 run's per-sequence mean
gap (item 7)". Make E the length-fair arm, not D. Keep D (geo mask) only if the counterfactual stats show that OPSM's
one-sided rule (negative advantage only) leaves > 0.5 % of tokens with |gap| > 0.5 unmasked. Drop the `is_ratio`
mode and its `ValueError` guard if D goes.

**7. [major] B (B2.1). Log the counterfactual effect of every mode in every arm. It is free and replaces most of the
probe.**
Problem: the IS correction does not change generation, only the weights in the loss. The masked, clamped and
down-weighted shares of all modes can be computed from the same `logps` and `sampling` tensors in one run. B2.1
logs only the active mode. B.1 item 5 also reads the median ratio mean 0.82 as a possible "18 % masked". But TRL masks
only ratios > 3 (`clip_min` None gives min_val = -inf [code: grpo_trainer.py:2803-2818]). A kept sequence with a sum
gap of -3 has weight 0.05, is silently down-weighted, and is counted nowhere.
Exact fix: in `is_summary`, add per step, for the logged tensors and regardless of the active mode:
`cf/seqmask3_share`, `cf/tok_clamp2_share`, `cf/geo_0.99_1.01_share`, `cf/opsm_share@{0.05,0.1,0.5}`,
`cf/seq_ratio_lt_0.1_share`, the quantiles q50/q99/q99.5 of the per-sequence mean gap, and the effective sample size
`(sum w)^2 / (N sum w^2)` of the active weights. All are device sums, read in the one existing `.tolist()`. Add a
unit test with planted values for each.

**8. [major] B (B2.3b). The attribution cannot separate dtype from the engine.**
Problem: B2.3b reads "only a gap present in (b) and (c) alike points to the dtype". A gap in both (b) (plain TRL
forward) and (c) (merged bf16 weights) is equally explained by vLLM-versus-HF kernels: attention backend, GEMM
tiling, CUDA graphs, and vLLM's own CPU preprocessing and vision encode (C §3.3 rows). None of the three ways
changes the dtype. Mean |gap| alone also hides the sign, and the sign decides what sequence masks do.
Exact fix: add
- (d) the plain trainer forward in fp32 (no autocast, base cast to fp32) on the same 256 answers, as the reference;
- (e) vLLM scoring of the same 256 prompt+answer sequences by teacher forcing (`SamplingParams(max_tokens=1,
  prompt_logprobs=0)`), in bf16 and in fp16 or fp32. The 2B model in fp32 is about 8.5 GB [estimate] and fits.

Report the signed mean gap, the mean |gap| and the gap by token position (bins of 64) for each pair. Reading:
- vLLM-bf16 vs HF-bf16 gap >> each one's gap to fp32: the engine kernels cause it.
- Both close to fp32 only in fp16 or fp32: the dtype causes it.
- (c) >> (b) on a drifted adapter: merge rounding causes it, and item 10 applies.

Run on the B1 adapters (item 1b).

**9. [major] B (B2.4 rule 3). The rule compares two different units.**
Problem: "the smaller `is_masked_share` + `is_clamped_share`". Under geo or OPSM masking, `is_masked_share` is a
share of sequences. Under TIS, `is_clamped_share` is a share of tokens. The geo arm also clamps tokens. The sum is
not comparable across arms.
Exact fix: define one number for every arm: `modified_token_share` = share of valid tokens whose final weight differs
from the raw ratio by more than 1 % (masked tokens included), plus the ESS of item 7. Rule 3: the arm with the higher
ESS wins. On a tie of 0.01 or less, take the one with the lower `modified_token_share`. On a second tie, take
bf16-tis.

**10. [major] C (C5 "Rejected ... vLLM native LoRA serving") and B2. The rejection rests on a false fact, and
it removes the fix for the most likely cause.**
Problem: C5 rejects vLLM LoRA serving because it would "change numerics vs the trainer's merged W + d". The trainer
never merges. PEFT runs W x + B A x unmerged, and `lora_merged_params` builds W + d only for vLLM, in W's bf16
[code: src/objdel/train/fast.py:696-717; `sync_without_merge` docstring "the base weights are never written"]. vLLM
LoRA serving would keep d apart from W, like the trainer does. It would remove the bf16 rounding of W + d, whose
growth with |d| fits item 1's trend.
Exact fix: in C5 replace the reason with: "Not adopted for speed: LoRA kernels may slow decode [estimate 10-30 %].
It is the numerics fix if B2.3b (c) >> (b) on a drifted adapter." Add to ADR-0012 option F, "vLLM `enable_lora` with
the adapter pushed every step", triggered by that B2.3b result. Add an optional G1 arm (15 steps) that measures its
decode cost.

**11. [major] B (B3.5b LocGuard). The guard will stop a healthy paid run.**
Problem: it fires when, over the last 50 logged steps, the least-squares slope of `rewards/loc/mean` is > 0 and the
slope of `rewards/iou_tc/mean` is <= 0. Both are noisy batch means: the cell mix changes every step, and only the
kind counts are fixed. When both true trends are near 0, each sign is a coin flip. The rule is checked every step for
about 250 windows, so a false stop is close to certain.
Exact fix: stop only when (a) the loc slope t-statistic is > 3 and (b) the upper 95 % bound of the IoU slope is < 0.
Compute both on per-step means over the off-centre rows only (`kind/shifted/*` logs), and require two disjoint
consecutive windows. Add a unit test: seeded noise with flat true trends (step sd as in R5) stops in < 1 % of 300-step
runs. The existing "rising loc, flat IoU" test must still stop within 100 steps.

**12. [major] B (B3.7 rule 1, Gate B3). The Spearman gate is not defined.**
Problem: "Spearman(A, -d / r) in all-miss groups". A is relative within a group, but d / r pooled across groups mixes
groups of different scale. An engineer can compute it in two ways that give different answers.
Exact fix: write "for each all-miss group with >= 3 distinct d values, the Spearman correlation of A_i and -d_i / r
over its 8 answers; the gate value is the mean over groups, with a nodule-level bootstrap 95 % CI (A1); gate: mean
> 0.3, CI lower bound > 0, and v5 mean - v4 mean > 0 paired by group".

**13. [major] C (C7) and B (B3). Runaway answers get no gradient, and a lower limit makes this worse.**
Problem: `mask_truncated_completions=True` [code: configs.py:188] masks every truncated answer. A runaway answer
(38 % of R5 steps hit 1,536; decode tail 4-5 s, BRIEF §2) then never gets a negative gradient. It gets only the
indirect effect through the group baseline. C7 at L = 1,200 would also turn the naturally ended 1,200-1,493-token
answers into masked rows: 6 of 15 R5 steps had one, and all are longer than any LIDC doctor outline. Those answers
carry a negative signal today. DAPO's overlong shaping exists for this case, and TRL ships it.
Exact fix: add to C7 (and to the B3 config as a measured option): "Option O: keep `mask_truncated_completions=True`
and add `trl.rewards.get_soft_overlong_punishment(max_completion_len=L, soft_punish_cache=200)` [code: .venv
trl/rewards/other_rewards.py:83]. With v5, add it to R as a term (weight 1, before RLOO), not as a separate TRL
reward (TRL would add it after our group-centred aggregator)." Log `completions/clipped_ratio` and the decode-tail
seconds (C0 6b) per run. Decide L only together with O.

**14. [major] Model choice (§3 table, owner decisions). The candidate models exceed a strict "up to 4B
parameters".**
Problem: the owner rule is "any model up to 4B parameters". Qwen3-VL-4B-Instruct has 4.438 B parameters,
Qwen3.5-4B 4.660 B and MedMO-4B 4.438 B [measured: HF API safetensors.total]. The doc shows 4.44 B but never says
that this is above 4.0 B.
Exact fix: add owner decision 12b: "'up to 4B' means the 4B size class (4.44 B allowed) or a hard 4.0 B limit (then
no 4B arm exists among the candidates, and the 4B gate is dropped)". Default: [size class].

**15. [major] 00_front (§1 summary, the §1 flow diagram, §6.3). The Gate C value differs from Phase C.**
Problem: the front says "Target 18 s per step or less, same rewards (Gate C)". The flow diagram gates the use of C
code in B on "<= 18 s/step". Phase C sets the must-pass gate G-SP1 at <= 20.0 s (mode A) or <= 20.7 s (mode B), and
keeps 18 s only for C1-C8 (G-C), which needs C7 and C8 inside B. An engineer following the front would keep C code
out of B even when G-SP1 passes.
Exact fix: in §1 and in the flow diagram write "Gate C (G-EQ1..6 and G-SP1: <= 20.0 s mode A / <= 20.7 s mode B)
decides whether B uses C code; 18 s (G-C) is reported after C7/C8".

**16. [major] C (C3, G-EQ4, G-EQ6). The bit-equality premise of C3 is fragile on the GPU.**
Problem: today `VisionCache` encodes the images new to each call: a varying subset per forward
[code: fast.py:280-291]. C3 encodes all 32 in one call. The varlen attention is per image, but the vision MLP and
QKV GEMMs see a different M (32 x 1,024 rows vs a few x 1,024). cuBLAS may pick another kernel and change bf16 bits.
That is not a bug. G-EQ4 (`torch.equal` on logprobs) and G-EQ6 (`grad_norm` to 4 digits) may then fail, which keeps
all C code out of Phase B.
Exact fix: in G-EQ4 compare against a baseline that encodes with the same grouping (call the tower once on the same
32 images through the old `pixel_values` path). Add a fallback rule: "if the two only differ through the encode
grouping (shown by the same-grouping test passing), accept max |delta logp| <= 2e-3 and step-1 metrics within the
R5 step-to-step spread". State this in C.1, the goal of "no gradient bit changes".

**17. [major] B (B2.2 fp16 plumbing step 0). The fp16 check covers only the vision tower.**
Problem: fp16 overflow can occur in the language model (residual-stream outliers above 65,504). The FP16 paper
tested no VLM and no Qwen3-VL [web: hf.co/papers/2510.26788]. Step 0 runs only the vision tower on 32 images.
Exact fix: step 0 runs the full model (vision + all 28 text layers) in fp16 on CPU. Use 32 prompts, each with its
saved A3 answers, including the longest (>= 1,200 tokens). Record the max |activation| per decoder layer and the
count of inf/NaN in logits. Drop the fp16 arm if any layer exceeds 0.5 x 65,504 or any inf/NaN appears. State in
ADR-0012 that vLLM's fp16 path is untested until the box.

### Minor

**18. [minor] Model choice (§3 reason 2). EGM numbers are mis-cited.**
Problem: "84.2 (2B) to 87.3 (4B) before training, and 89.6 vs 91.0 after SFT+RL; gap +3.1 to +1.4 [paper:
2601.13633 Table 2]". Table 2 of that paper is the COMPLEX-PROMPT subset. The 84.2 and 87.3 are Instruct base rows
of Table 4. The 89.6 and 91.0 are EGM on the Thinking variants, whose base is 83.6 and 87.2 (Table 9) [web:
arxiv.org/html/2601.13633].
Exact fix: "Qwen3-VL-Thinking RefCOCO-family average 83.6 (2B) vs 87.2 (4B) before, 89.6 vs 91.0 after SFT+RL: the
gap goes from +3.6 to +1.4 [paper: 2601.13633 Tables 4 and 9]".

**19. [minor] Model choice (§3 table row 2). "Vision tower is identical" is too strong.**
Problem: the blocks are the same (depth 24, hidden 1024, deepstack [5, 11, 17]). The merger output differs (2048 vs
2560) [measured: HF config.json], so cached features and the C3/C6b shapes differ.
Exact fix: "same vision blocks (depth 24, width 1024); merger output 2,048 vs 2,560".

**20. [minor] Model choice (§3 "Free check"). The check is already done, and its reading is wrong.**
Problem: the gap 0.740 - 0.722 = 0.018 is already measured (reason 1). "A bigger model cannot gain there" bounds only
found images. It says nothing about misses, which is what the 4B gate is about.
Exact fix: delete the paragraph, or reword: "found-image headroom is 0.018 [measured]; a bigger model could gain
only through fewer misses, which the 4B gate measures".

**21. [minor] B (B.7 row "KL spike", B2.4). The KL spike did not affect the update.**
Problem: the 398,284 KL at step 746 came with grad_norm 0.0041, against a run median of 0.0048 [measured: rl.json].
The risk row implies harm.
Exact fix: reword: "one logged KL spike (398 k, step 746) left grad_norm at 0.0041 (median 0.0048); logged as a
diagnostic, not a driver; median KL 0.32".

**22. [minor] 99_back (owner decision 13) vs §3. Two names for the medical model.**
Problem: owner decision 13 says `MedMO-4B-Next`, while §3 says `MBZUAI/MedMO-4B`. The HF id `MBZUAI/MedMO-4B`
exists (Qwen3VLForConditionalGeneration, 4.438 B) [measured: HF API].
Exact fix: use `MBZUAI/MedMO-4B` in both places, or name the exact repo id of "-Next" after checking it on the API.

**23. [minor] 00_front (§2.5 diagram). The diagram shows parts that are not adopted.**
Problem: the line "old logp (no-merge LoRA, shared-prefix forward, fused log-prob K3b)". K3b is "Not in Phase C"
(C.5), and the old-logp pass is skipped in the aligned case [code: fast.py `skip_old_logps`].
Exact fix: write "update forward computes logp and the IS ratio (no separate old-logp pass)". Remove "K3b".

**24. [minor] B (Gate B1, ADR-0014). The B1-vs-R5 gain mixes two causes.**
Problem: B1 = rl-full-3ep plus 300 more steps. R5 = rl-full-3ep. The gain is "shifted data" plus "more steps". No
centred-data control exists, and none fits the budget.
Exact fix: in ADR-0014, cite the plateau evidence (held-out IoU flat over the last epochs of the 770-step run,
`rl_plateau_diagnosis_2026-11-08.md`) as the reason more centred steps would add about 0. State that the per-cell
logs (offset < 0.15 cells flat, off-centre cells rising) are the within-run check.

## Counts

blocking 3, major 14, minor 7.
