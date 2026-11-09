# Final review, round 2: RL for LLM/VLM (B2, B3, C, model choice)

Reviewer role: RL research scientist (GRPO, IS correction, reward design, TRL/vLLM internals).
Target: `REPO/docs/implementation2.0/Implementation.md` (13,065 lines, round 2). Earlier file read:
`reviews/final_r1_rl.md` (24 items). I checked B.R2, C.R3, the front and back review logs. Items fixed there are not
repeated here. Where a round-1 fix created a new defect, the item says so.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill covers RL numerics or LoRA serving. No item comes from it. Nothing installed. |
| engineering:code-review | I read the TRL 1.17.1 source in `.venv` again: IS modes (`grpo_trainer.py:2776-2824`), OPSM (`:3066-3089`, `:3190-3203`), loss and `dr_grpo` normaliser (`:3255-3290`), truncation mask (`:2603-2611`). I also read `fast.py:350-363, 495-717`, `reward/advantages.py`, `reward/score.py:25-27` and `configs.py`. Items 2, 4, 9, 13 and 14 come from this. |
| engineering:debug | I tested the hypothesis "merge rounding causes the mismatch" on the saved adapters at $0 (item 1). I also split the runaway share by step range (item 3). |
| engineering:system-design | I checked how V1 hands work to B-S1 and A3. The V1 minutes leave out three items that Phase B needs (item 10). |
| engineering:architecture | ADR-0012/13/14 in B against the front ADR table (item 9). Option F is now first-class, backed by a library API (vLLM `load_inplace`) (item 1). |
| engineering:documentation | Each item has a severity, a section, the problem and an exact fix. |
| huggingface-skills:huggingface-best | HF API listing on 2026-11-08 (trending image-text-to-text, and Qwen org searches for "Qwen3.8" and "Qwen3.5"). No new model of 4B or less exists: the Qwen3.8 family is 27B and larger. The section 3 candidate list stands. Item 11 is about the 4B path, not the list. |
| huggingface-skills:trl-training | TRL facts: `mask_truncated_completions` zeroes only the loss mask, and the truncated row's reward stays in the group baseline (item 3). `get_soft_overlong_punishment` formula (item 7). OPSM docstring "expects a positive value (e.g., 0.5)" (item 15). TRL 1.17.1 has no vLLM LoRA path: `grep -rn "enable_lora\|LoRARequest" trl/` gives 0 hits (item 1). |
| huggingface-skills:huggingface-papers | I re-read DAPO (2503.14476) Eq. 13 and its overlong-filter rationale, used in item 3. I re-read the FP16 paper (2510.26788) mechanism, used in item 1: fewer mantissa bits are lost. |
| Web | vLLM PR 31326 "allow inplace loading lora" (`LoRARequest(load_inplace=True)`, merged 2026-01-20; docs "In-Place LoRA Reloading"). ms-swift PR 9077 (LoRA-only weight sync with in-place reload for GRPO). vLLM issue 20149 (TRL merges LoRA because vLLM has no tensor-update API). |

## Measurements made for this review (2026-11-08)

**M1. LoRA delta lost when vLLM receives W + d rounded to bf16.** This is what `lora_merged_params` pushes
[code: src/objdel/train/fast.py:717]. Script:
`docs/implementation2.0/notes/r2rl/merge_rounding.py` (`.venv/bin/python -I`). Inputs: the adapters in
`runs/rl/train_rl/trainer/checkpoint-{50..770}` (r 16, alpha 32). The W proxy is
`runs/rl/train_rl/final/model.safetensors`; |d| is much smaller than |W|, so the bf16 exponent of W is the same.
The sample is layers 0, 9, 18 and 27, modules q, v, down and up.

| step | share of delta elements with \|d\| < ulp(W)/2 (erased) | \|\|round(W+d) - W - d\|\| / \|\|d\|\| | same in fp16 |
|---|---|---|---|
| 50 | 0.692 | 0.605 | 0.128 |
| 100 | 0.643 | 0.531 | 0.102 |
| 300 | 0.557 | 0.426 | 0.071 |
| 500 | 0.512 | 0.379 | 0.060 |
| 770 | 0.464 | 0.339 | 0.051 |

At step 50, vLLM samples from a policy whose LoRA update is about 60 % wrong in norm. Two thirds of the delta
elements do not reach vLLM at all. The trainer computes W x + s·B A x unmerged, so it sees the whole delta. This fits
B.1 item 5. The logp gap rises while d is still below one ulp (steps 21-100: 0.028). It then falls as the relative
error falls (steps 501-770: 0.016).

**M2. Runaway share by step range** [measured now: runs/rl/train_rl/rl.json, `completions/clipped_ratio`,
`completions/max_length`]:

| steps | clipped share of answers | share of steps with a 1,536-token answer |
|---|---|---|
| 1-20 | 0.00195 | 0.400 |
| 21-100 | 0.00059 | 0.138 |
| 101-300 | 0.00182 | 0.255 |
| 301-500 | 0.00111 | 0.135 |
| 501-770 | 0.00373 | 0.315 |

**M3.** The full RL run and R5 used `lr` 5e-5 [measured: rl.json, runs/r5/shared_rl.json `lr`]. `configs.py` says
1e-5 for LoRA [code: configs.py grpo_config]. `train_rl` takes "the G-probe's LR unless --lr"
[code: experiments/phase2/train_rl.py:8, 33].

## Items

### Blocking

**1. [blocking] B (B2, B2.3b, B2.4, ADR-0012 option F, B.8.3 item 9) and C (C5, G-EQ5). Merge rounding is a
measured, dominant cause of the mismatch, but the plan only tests it after a paid run and makes the fix
conditional.**
Problem: M1 shows that 46-69 % of the LoRA delta elements are erased in the bf16 push. The relative error of the
pushed delta is 0.34-0.61. vLLM therefore samples from a lagging, shrunken policy at every step, by design. This is
a systematic off-policy bias, not noise. The plan does three things wrong here:
- It runs B2.3b (c) only at the end of the paid B-S1 (15 min, about $0.16).
- It triggers option F only "if B2.3b points to merge rounding".
- It keeps fp16 as an optional extra arm.
B2.4 rules 1-3 cannot choose F at all. C5 then fixes the lossy push in place ("keep the values exactly
`p.data + delta.to(p.dtype)`"). Also, TRL 1.17.1 has no vLLM LoRA path (0 grep hits), so F as written ("vLLM
`enable_lora`") has no plan for its code.
Exact fix:
- (a) New task **B2.0 ($0, Mac, before V1)**. Promote `merge_rounding.py` to
  `experiments/phase2/merge_rounding.py` with a test on a planted 2x2 case. Run it on all 28 layers and all 7
  modules of the saved checkpoints. Write the table into ADR-0013 (front numbering, item 9).
- (b) Move B2.3b ways (a)-(d) to the Mac CPU on `checkpoint-100` and `checkpoint-770` of the full run, at $0. Use 32
  of A3's saved answers, or any 32 saved R5 answers before A3. Keep only way (e) and the generation-time logps
  (item 5) for the box. Remove the 15 min from B-S1 (B.8.1) and from V1 (item 10).
- (c) ADR options: make F a first-class candidate for B3. F means vLLM `enable_lora=True, max_lora_rank=16`, and
  each step `LoRARequest("policy", 1, <tmpfs adapter dir>, load_inplace=True)` with
  `model.save_pretrained(dir, safe_serialization=True)` of the adapter only [web: vLLM PR 31326; docs "In-Place
  LoRA Reloading"]. Its patch points are the `LLM(...)` kwargs (the same wrapper as `fast.vllm_dtype`) and
  `lora_request=` in the `patch_vllm_prompts` wrapper. Add both names to the A7 start-up check. Keep A (fp16 both)
  as the second candidate (M1: relative error 0.13 → 0.05).
- (d) New B2.4 rule 0: "If B2.0 shows a relative error >= 0.2 at step 100, B3 uses F if its 15-step decode arm costs
  <= 15 % of step time; else A if the fp16 Mac check passes; else the control with OPSM (option E)".
- (e) C5: add "under F, the push is the adapter only (about 17 M parameters), not W + d; G-EQ5 then compares the
  vLLM LoRA tensors with the trainer's".
- (f) B.1 item 5: add the M1 table and the sentence "merge rounding erases 46-69 % of delta elements
  [measured: B2.0]".

**2. [blocking] B (B3.4 `advantages_v5` drop rule). The round-1 fix dropped `legal` from the span test, so
groups where only legality differs get A = 0.**
Problem: v5 writes
`flat = ptp(stack([iou_tc, f5, g_area, loc, r_topo])) ... < drop_spread`. A parsed but illegal answer keeps its raw
terms for logging [code: src/objdel/reward/score.py:25-27]. Take a group where 6 legal answers and 2 illegal answers
draw almost the same ring with IoU_tc 0.3 (illegal through a pinch or a spike, for example). The `iou_tc`, `f5`,
`g_area` and `loc` spans are below 0.01. r_topo is 0 for all, because IoU < `iou_gate` 0.5. So `flat` is True and
A = 0. Yet R is 0.3 + 0.1·f5 + area for the legal answers and 0 for the illegal ones. The illegal answers get no
negative gradient. v4 tests the span of `TERMS`, which includes `legal`
[code: src/objdel/reward/spec.py:26, advantages.py:136].
Exact fix: in B3.4 write
`flat = (np.ptp(np.stack([legal, iou_tc, f5, g_area, loc, r_topo], -1), axis=1).max(-1) < spec.drop_spread) & (sig_id == sig_id[:, :1]).all(1)`.
Add test **T15**: the group above (6 legal, 2 illegal, equal raw terms, IoU 0.3) keeps a non-zero A, and the
illegal answers get A < 0. Add T15 to B.9.13 and to Gate B3's "T1-T12, T14, T15".

**3. [blocking] B (B.1 item 6 wording, B3.7 rule 9, B.7 row "Runaway answers") and C (C7 step 5, C11 step 4).
Masking truncated answers is not neutral. It biases the update toward truncation, and C11 as written makes
runaways more likely.**
Problem: TRL zeroes only the loss mask of a truncated row. The row's reward (0: unparseable, so every term is 0)
stays in the group baseline [code: grpo_trainer.py:2603-2611; score.py:25]. For one row,
E[1{not truncated}·(R − b)·∇log π] = ∇E[R·1{nt}] − b·∇P(nt). Since R = 0 when truncated, the first term is ∇J. So the
masked update is ∇J − b·∇P(not truncated): it raises P(truncated) with weight b = baseline > 0. With RLOO, b is
independent of the row's own answer, so the algebra holds. The doc says a runaway "gets no direct negative
gradient. It acts only through the group mean". That reads as neutral, and it is not. M2 fits the bias: the runaway
share rose from 0.00059 (steps 21-100) to 0.00373 (steps 501-770). The share of steps with a 1,536-token answer rose
from 13.8 % to 31.5 %. That is the decode tail Phase C fights. C11 caps rows at `1.5·T_i + 64` and masks them, so
the same bias then pushes answers past cap_i. Eval has no cap, so more runaways reach test (scored 0). DAPO's reason
to filter (a cut answer may be correct) does not apply at L = 1,536: eval truncates at the same 1,536 and scores the
answer 0.
Exact fix:
- B3 / reward v5: set `mask_truncated_completions=False` for the B3 run at L = 1,536. A truncated row keeps R = 0,
  which is exactly its eval score, and gets A = −b̄. Keep `True` in B1, so Gate B1 stays a data-only test. Record it
  in ADR-0014 (v5) and in B3.5.
- C11: do not mask capped rows. G-RC1 already bounds found answers above cap_i at < 0.5 %. Optionally, zero the loss
  mask of the first min(len, T_i) tokens of a capped row, so only the tail past the doctor-length point gets the
  negative signal. Change C11 step 4 and the C.7 row.
- Replace the wording in B.1 item 6, B3.7 rule 9, B.7, C7 step 5 and C11 with: "masked truncated rows bias the
  update by −b·∇P(not truncated) (b = group baseline > 0); runaways rose 0.00059 → 0.00373 over the full run
  [measured: rl.json]".
- Add to every B run report: the slope of `completions/clipped_ratio` over steps. B3 gate line (report): the B3
  clipped share is no higher than B1's.
- Option O then becomes unnecessary for runaways (see item 7). Default: drop it.

**4. [blocking] B (B2.4 rule 3, Gate B2, B2.1 `ess`). The preset rule cannot be computed, and it is degenerate.**
Problem:
- (a) B2.1 logs `ess` only "of the active weights". The B1 run uses C, so the ESS of B and E is never logged, yet
  rule 3 picks "the one with the higher counterfactual ESS". Gate B2 compares "the chosen preset's ESS" with the
  control's.
- (b) The ESS of token weights (B, E) and of sequence weights (C) are different units.
- (c) Maximising ESS favours the least correction. Weight 1 everywhere (no IS) has ESS 1. E is B times a 0/1
  sequence mask, so ESS(E) <= ESS(B) almost always. E, the length-fair option that ADR-0012 prefers, can never win.
  The rule picks B by construction.
Exact fix:
- B2.1: add per-mode counterfactual ESS on one unit (token level; sequence weights broadcast to their tokens):
  `cf/ess_C`, `cf/ess_B`, `cf/ess_E@δ`. Report them; do not maximise them.
- B2.4 rule 3, rewrite as: "Default E (token TIS C = 2 + OPSM δ). Veto rules over B1 steps 50-300: if
  `cf/opsm_share@δ` > 2 % of sequences, raise δ to the q99.9 gap; if `cf/tok_clamp2_share` > 1 % of tokens, record
  it and keep E; fall back to C only if E's `cf/ess_E` < 0.5." Rule 0 (item 1) comes first.
- Gate B2: replace the ESS line with "the chosen preset and its veto values are recorded in ADR-0013; NaN/inf = 0".
- Update U51 and U54.

### Major

**5. [major] B (B2.3b way (e)). Teacher-forced vLLM scoring is not the quantity that enters the IS ratio.**
Problem: `SamplingParams(max_tokens=1, prompt_logprobs=0)` computes logprobs in a prefill pass. The IS ratio uses
the logprobs vLLM returned while decoding (`sampling_per_token_logps`). Prefill and decode use different attention
kernels and batch shapes, so "(e) bf16 vs (b)" can differ from the real gap. An fp32 vLLM engine also needs a
second `LLM` in the same process, beside the colocated bf16 one and the trainer. Its attention backend in fp32 is
untested for Qwen3-VL, so the "about 15 min" may fail.
Exact fix: add way **(f)**, the generation-time `sampling_per_token_logps` of the same 256 answers. This is the
"vLLM-bf16" line of every reading. Keep (e)-bf16 as a separate "prefill vs decode" line. Run (e)-fp32 only in a
fresh process after the trainer is freed, with `VLLM_ATTENTION_BACKEND` named. If it fails, record the failure and
do not retry. Readings: replace "vLLM-bf16 vs (b)" by "(f) vs (b)".

**6. [major] B (B3.7, Gate B3). No check that misses still get a usable share of the gradient once the std
division is gone.**
Problem: v4 z-scored each term with floors 0.05, so all-miss groups got unit-scale advantages. v5 is mean-centred
only. In an all-miss group A comes from `0.1·loc` and small `w_area` terms. B3.4's own worked check gives
+0.057 for the best answer. A group with hits and misses has A of about ±0.3-0.6 [estimate: IoU 0.7 vs 0]. So the
location signal is about 10x weaker per answer, and AdamW (shared second moments) does not undo a relative scale
between groups. Rule 1 tests only the order of A inside all-miss groups, not its weight in the update.
Exact fix: add B3.7 rule 10, **"Gradient share of all-miss groups"**:
`S = Σ_{all-miss groups} Σ_i |A_i| / Σ_{all groups} Σ_i |A_i|`, for v4 and v5, on val-shape answers. Report S
beside the all-miss group share of the batch. Gate (pre-run, in Gate B3): S_v5 >= 0.5·S_v4. If it fails, raise
`alpha_loc` and `miss_iou` together (keeping `alpha_loc <= miss_iou`), or set `loc_width` per rule 5, and re-measure.
Record S in ADR-0014.

**7. [major] C (C7 option O) and B (B3.7 rule 9, B.8.3 item 7). At L = 1,200, option O penalises correct
outlines of the largest nodules.**
Problem: TRL's penalty is ((L − cache) − |y|) / cache for L − cache < |y| <= L [code: trl/rewards/other_rewards.py:89-95].
With L = 1,200 and cache = 200, it starts at 1,000 tokens. LIDC doctor outlines reach 1,182 tokens [measured: A1 F9
via C7]. A correct 1,182-token answer then gets −0.91. That is about its whole IoU reward.
Exact fix: if option O is kept after item 3, write in C7 step 5: "require `L − soft_punish_cache >= max doctor
outline tokens + 32` (1,214 today). At L = 1,200 option O is off; at L = 1,536 cache <= 290". Add a config check
that raises otherwise, plus a unit test. Default: off (item 3 removes the need).

**8. [major] B (B3.7 rule 8, B3.3 `c_topo`). The topology weight has a lower bound only, so v5 can become a step
function at IoU 0.5.**
Problem: rule 8 sets c as the smallest c with `c·eps >= q95 margin`. For a valid answer with π ≈ f, r_topo ≈ 1, so
R jumps by about c at IoU_min = 0.5 (valid needs IoU_min >= `iou_gate`). Say the q95 margin is 0.2 [estimate]: then
c = 2 and the jump is 2, against IoU's whole range of 1. Above the gate, IoU then matters about one third as much as
passing it. The B.7 row names the risk, but no rule or gate holds it.
Exact fix: add to B3.7 rule 8: "report the R jump at the gate (c × median r_topo of valid answers) and the share of
A variance in hit groups from `c_topo·r_topo`. Gate: share <= 0.5. If the lower bound forces a larger c, the owner
decides between a smaller `eps` and accepting the share; record the choice in ADR-0014".

**9. [major] B (header, B.3.4, B2, B3), C (C5, C.8.3). ADR numbers still collide with the front table.**
Problem: front section 4 says 0012 = docs site (DOC10), 0013 = number match, 0014 = reward v5, 0015 = scheduled RL
rows. The back matter agrees ("0013 numerics and beta", "0015 polygonizer exception"). B still writes ADR-0012 =
number match, 0013 = reward v5 and 0014 = scheduled rows plus the polygonizer exception. C5 and C.8.3 cite
"ADR-0012 option F". An engineer would write two files numbered 0012.
Exact fix: in B and C replace ADR-0012 → 0013, ADR-0013 → 0014, ADR-0014 → 0015 (every occurrence, including B.0,
B.7, B.8.3, B.9 and C5). Delete the B header sentence "rename them here in the same commit" once done. Add a DOC4
drift test: every `ADR-00NN` in the doc matches the front table's title.

**10. [major] 00_front (section 6.2 V1 table). V1 leaves out three items that Phase B and A3 need, so the $2.20 cap
is too small.**
Problem: the front V1 has 30 + 137.5 + 6.9 + 28 = 202 min ($2.16). The following are missing:
- B-S1's `edge` view eval of the B1 model (6.9 min). Gate B1 has a hard rule on it.
- B2.3b (15 min).
- The R5 `edge` baseline. B.8.4 says A3 measures it, but A3's 17.5 min is val + test only (A.8).
That is +28.8 min, so V1 is about 231 min = $2.47 at 27.5 s, above the $2.20 cap. `BudgetStop` at $1.90 would cut
the B1 run or skip A3.
Exact fix: move B2.3b (a)-(d) to the Mac (item 1b), which saves 15 min. Add the two `edge` evals (13.8 min) to V1:
217 min = $2.32 at 27.5 s, $2.11 at 23.65 s. Set the cap at $2.35, or drop A3's test half from V1 (test is needed
only by C-G-EQ2/3 and D7, not by Gate B1). Update section 6.3's money left, A.8 and B.8.1.

**11. [major] Model choice (front section 3 "If the gate fires", back decision 17, open question 2). The 4B path
has no SFT step.**
Problem: steps 1-4 add the model id, a memory value, "a 15-step probe first" and an accept rule. RL needs a 4B
checkpoint that already writes the polygon text. The 2B had SFT first (`rl-full-3ep` comes from an SFT merge).
Reason 3 says the 4B "needs its own SFT, whose cost is unknown", but no step runs it, prices it or names its
checkpoint. A probe from raw `Qwen3-VL-4B-Instruct` would measure a model that cannot yet write the format.
Exact fix: insert step 0, "SFT the 4B on the B1 rows' targets with `sft_config` (LoRA r16, same schedule), priced by
a 50-step timing probe ([estimate] until measured). The RL probe and run start from that merged checkpoint." Add the
SFT minutes and dollars to the 4B gate text. The owner sees the full 4B cost (SFT + RL) before the gate decision.

### Minor

**12. [minor] B (B2.1 quantiles). The merge rule cannot hold quantiles.**
Problem: `set_is_ratio` runs once per micro-batch (16 per step) and merges `is_stats` by sums, min and max
[code: fast.py:581-602]. q50, q99 and q99.5 of the per-sequence gap cannot be merged that way.
Exact fix: keep a device buffer of the 256 per-sequence mean gaps (and advantages) per generation batch, made by
`torch.cat` in `set_is_ratio`. Compute the quantiles and the OPSM counts in `log`, before the one `.tolist()`. ESS
stays mergeable (Σw, Σw², N).

**13. [minor] B (B2.2 `NUMERICS`). The sketch is not valid Python.**
Problem: `"bf16-tis-opsm": {**NUMERICS["bf16-tis"], ...}` refers to `NUMERICS` inside its own literal.
Exact fix: define `_TIS = {...}` first, then `NUMERICS = {"bf16-tis": _TIS, "bf16-tis-opsm": {**_TIS, ...},
"fp16-tis": {**_TIS, ...}}`.

**14. [minor] B (B2.3 fp16 arm). The reward comparison is confounded by the IS mode.**
Problem: the fp16 arm runs `fp16-tis` (token TIS), and the B1 run runs `bf16-seqmask`. The (fp16 − B1)
`rewards/iou_tc` bound mixes dtype and IS mode.
Exact fix: run the arm as `fp16-seqmask` (control IS mode, fp16 numerics). Add that preset. The logp-gap rule is
unaffected.

**15. [minor] B (B2.4, ADR option E). OPSM δ from the B1 run is stale if B3 changes numerics.**
Problem: δ = q99.5 of the B1 run's gap, which is measured under bf16 merge rounding (M1). Under F or fp16 the gap
shrinks, so δ is too loose and OPSM masks almost nothing. TRL's docstring suggests 0.5.
Exact fix: "δ = q99.5 of the per-sequence mean gap over the first 20 steps of the B3 run under its own preset,
computed live from the counterfactual buffer (item 12), then frozen. Until step 20, δ = 0.5 (TRL docstring)."

**16. [minor] B (B3.5b LocGuard, B1.5 logging). The guard reads a log that nobody writes.**
Problem: LocGuard needs per-kind means of the location term. B1.5 logs only `kind/<k>/iou_tc`, and B3.5 adds `loc`
only to the global `LOGGED`.
Exact fix: in B3.5 add "`TRLReward` logs `kind/<k>/loc` (mean of the gated term before `alpha_loc`) for every kind".
U47/U58 read that key.

**17. [minor] 00_front (section 1 flow diagram, section 6.3, FR4). The "B2 probe" is stale.**
Problem: B removed the 4-arm probe (B.R2 item 4). The front still says "V3: R-D2 first, B2 probe, B3 run", "Then the
B2 probe, the B3 run" and FR4 "B2 probe with Gate B2".
Exact fix: write "optional fp16 arm (60 steps)". Change FR4 to "B2.0 merge-rounding table ($0) and B1-run
counterfactual stats; preset by B2.4 (F, A or E)".

**18. [minor] B (B1.5, B3.9 configs). The learning rate is not stated, and B2.0 depends on it.**
Problem: R5 and the full run used `lr` 5e-5, while `configs.py` says 1e-5 for LoRA (M3). At 1e-5 the deltas are
5x smaller and the share erased by merge rounding is larger (M1).
Exact fix: write `--lr 5e-5` (as R5) in `b1_rows.yaml` and `b3_v5.yaml`, and record `lr` in ADR-0013 beside the
B2.0 table. Fix the `configs.py` docstring or the default in CLEAN9 so there is one value.

## Counts

blocking 4, major 7, minor 7.
