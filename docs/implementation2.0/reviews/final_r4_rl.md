# Final review, round 4: RL for LLM/VLM (B2, B3, C, model choice)

Reviewer role: RL research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Target: `REPO/docs/implementation2.0/Implementation.md` (14,787 lines), read through its sources
`docs/implementation2.0/{00_front,B,C,99_back}.md` (B2, B3, C5-C11, C.6, front sections 2.5, 3, 6, back decisions 13-17).
Earlier files read: `reviews/final_r3_rl.md` (15 items) and `REPO/docs/implementation2.0/review-log.md`. The round-3
items are handled in the B text (rule (b), F sync order, `add_lora` once, packed check, `seq_kl` sign, live delta
callback, `adam_epsilon = 1e-8 / r`, LocGuard Newey-West, `truncated_grad_share`). None is repeated here. Item 1 is
the part of round-3 item 8 that B fixed but the front, C and the back matter did not follow. Item 2 is a gap that
99_back decision 17 itself names ("Front section 3 step 0 ... must add it") and that is still open, plus one new
defect in the same step.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers GRPO numerics, vLLM LoRA serving or IS correction. No item comes from it. Nothing installed. |
| engineering:code-review | Read TRL 1.17.1 in `.venv` (`grpo_trainer.py:882, 1029-1036, 2775-2824, 3060-3089, 3180-3195, 3286`; `grpo_config.py:348-352, 411`), `src/objdel/train/fast.py:351-363, 495-717`, `configs.py:73-196`, `experiments/phase2/train_sft.py:158-175, 285-292`, `experiments/phase2/vast/policy.py:422-452`, `experiments/phase2/configs/r1_qwen3vl2b.yaml`, `g_probe.py:88`. Items 2, 3 and 6 come from this. |
| engineering:debug | Traced the F timing rule through vLLM 0.32.0 (`lora/layers/base_linear.py:198-233`): a request with `lora_request=None` on an `enable_lora` engine still runs the punica path (item 4). Traced what feeds the F output-level check (item 5). Traced the 4B SFT launch to the first raise (item 2). |
| engineering:system-design | Traced the rental order V1 -> V2 -> V3 against ADR-0013 rule (b): which preset each run uses, and which C measurement carries to which run (items 1, 6). |
| engineering:architecture | ADR-0013 rule (b) is the rule in B; the front, C and 99_back still state the round-3 rule "F is the expected B3 preset" (item 1). |
| engineering:documentation | Each item names the source file and heading, quotes the doc line and the code line, and gives an exact fix. |
| huggingface-skills:huggingface-best | Re-checked the model-choice section against the HF facts it cites (Qwen3-VL-4B 4.438 B, Qwen3.5-4B 4.660 B). The candidate list stands. Item 2 is about the 4B SFT launch recipe, not the list. |
| huggingface-skills:trl-training | GRPO facts behind items 3, 6 and 7: `logging_steps` defaults to 1 in our runs (`g_probe.py:88`); TRL's OPSM reads `self.off_policy_mask_threshold` per loss call; the sequence IS ratio sums over every unmasked completion token; dr_grpo divides by `max_completion_length`. |
| huggingface-skills:huggingface-papers | Tried the EGM paper (2601.13633) through the HF papers API: not indexed on HF ("Paper not found"). The 2B/4B numbers of front section 3 reason 2 were not re-verified; no item depends on them. |
| Web | Downloaded vLLM v0.32.0 source (raw GitHub, tag `v0.32.0`, 2026-11-08): `lora/layers/base.py` (wrapper `load_weights` delegates to the base layer, so a mode-B base push into a LoRA-enabled engine is plausible; no item raised), `lora/layers/base_linear.py:198-233` (item 4), `lora/worker_manager.py:274-323` and `lora/model_manager.py:346-404, 1338-1405` (with `load_inplace=True`, `remove_adapter` deactivates and `activate_adapter` re-runs `set_lora` into the slots, so B2.2's re-add after a level-2 wake does restore the slots; no item raised), `model_executor/models/qwen2.py:453-521`. |

## Items

### Blocking

**1. [blocking] 00_front (section 1 diagram; section 2.5 "Sources" paragraph; section 6.3 rule "Numerics before
speed", V2 items 2b and 3, V3 items 2a and 2b), C (header "Numerics before speed"; C5 paragraph 2; C.6 protocol
"Numerics first"; C.8.2; C.8.3 row "Phase B2"), 99_back (decision 16) against B (B2.2 "Option F plumbing", B2.4
"Scope", B.8.1, B.8.3 item 9, B.8.4 row "Faster step", Gate B2). Two rules for the B3 preset and the F arm, and both
claim to be the rule.**

Problem. B now applies ADR-0013 rule (b):
- B2.4 "Scope": "The B3 run keeps B1's preset `bf16-seqmask`. ... B2.4 picks the preset for the runs after Gate B3."
- B2.2: "box check in the first run that uses F, after Gate B3 under rule (b)".
- B.8.4 row "Faster step": "Under ADR-0013 rule (b) the B3 run uses `bf16-seqmask` and C5's merged push, so a
  C-G-SP1 pass measured on the merged push carries over to B3, and no F arm is needed in V2."
- B.8.3 item 9: "Its F check is the first 15 steps of that run (0 extra minutes)".

The other parts still carry the round-3 rule:
- 00_front 2.5: "The seed value is 0.531 [measured: B.1 item 5], so F is the expected B3 preset."
- 00_front 6.3 V2 item 2b: "F arm (B2.4 rule 0): 15 steps under `lora-seqmask` ... | 13 | 0.15 | ... | F fails:
  items 3 and B3 use the preset that rule 0 falls back to". Item 3: "both under item 2b's preset".
- 00_front 6.3 V3 item 2a: "B2.3b box part, only if B2.4 rule 0 did not pick F or A"; item 2b: "If B then picks
  `fp16-seqmask`, C's G-EQ1, G-EQ4 and G-EQ5 re-run before item 4" (item 4 is the B3 run).
- C header: "C measures the weight path that B3 will run. ... The F check therefore runs before the C arms (V2 item
  2b, C.8.2), and the C arms run under the preset it picks". C5 and C.6 "Numerics first" say the same.
- 99_back 16: "its 15-step arm is rental V2 item 2b, before the C arms (front 6.3 prices it); if it fails, B3 uses
  the preset that rule 0 falls back to".

Consequence for an engineer who follows front section 6 (which says "If a phase number differs from this section,
this section wins"): V2 spends 13 min = $0.15 on an F arm that B says is not needed, then measures C-G-SP1 under
`lora-seqmask`. C itself says a speed result does not carry across the push type: "a speed result on the merged push
does not carry over to F" (C.6). So the C-G-SP1 pass that adopts C code for B3 is measured on a path B3 does not run
(B3 runs `bf16-seqmask` with the merged push, B3.9). Also V3 item 2a is skipped on a condition ("rule 0 picked F")
that under rule (b) is decided only after Gate B3, and V3 item 2b re-runs G-EQ gates before a B3 run that never uses
fp16. Code facts behind the push difference: merged push values `p.data + d.to(p.dtype)` [code:
src/objdel/train/fast.py:717]; F needs `enable_lora` and the punica path in every linear [web: vllm v0.32.0
lora/layers/base_linear.py:208-233].

Exact fix (adopt rule (b) everywhere, as B and review-log row 12 already decided):
- 00_front 6.3: delete V2 item 2b and the "Numerics before speed" rule. New rule: "C arms run under `bf16-seqmask`,
  the B3 preset (ADR-0013 rule (b)). The F check is the first 15 steps of the first run after Gate B3 (B2.4 rule
  0)." V2 item 3: "both under `bf16-seqmask`". Re-sum "Required for C adoption" to items 1, 2, 3 = 90 min = $1.12 at
  r_s, and "Items 1-5" to 138 min = $1.66 [estimate: the table rows minus 13 min / $0.15].
- 00_front 6.3 V3 item 2a: "B2.3b box part, only if the B2.0 step-100 relative error is < 0.2 (else the F check of
  the first run after Gate B3 needs the comparator of item 5 of this review)". V3 item 2b: "If B2.4 later picks
  `fp16-seqmask`, C's G-EQ1, G-EQ4 and G-EQ5 re-run inside the first run that uses it", not "before item 4".
- 00_front 2.5: "so F is the expected preset of the runs after Gate B3; B1 and B3 run `bf16-seqmask` (ADR-0013 rule
  (b))". Section 1 diagram: V2 line "c_speed_check, C arms, R-D1", without "F arm".
- C header, C5 paragraph 2, C.6 "Numerics first", C.8.2 and C.8.3 row "Phase B2": "C arms run under `bf16-seqmask`
  (the B3 preset, ADR-0013 rule (b)). Under `lora-seqmask` (runs after Gate B3), the first F run reports its step
  time against the B3 run as a report line; G-SP1 is not re-gated." Keep the F rows of C5 step 5(c) and G-EQ5 for
  that first F run.
- 99_back 16: "Option F [picked by B2.4 rule 0 for the runs after Gate B3 when the B2.0 step-100 relative error is
  >= 0.2; its check is the first 15 steps of that run; B1 and B3 run `bf16-seqmask`]".
- Add a DOC4 drift check: the strings "expected B3 preset" and "item 2b's preset" must not occur.

**2. [blocking] 00_front (section 3 "If the gate fires" step 0) against 99_back (decision 17). The 4B SFT launch
recipe fails at start-up on the paid box, and its "val stays unread" claim is false.**

Problem.
- Front step 0: "Use the 2B recipe unchanged: `experiments/phase2/configs/r1_qwen3vl2b.yaml` with only `model:
  qwen3-vl-4b` and `output_dir` changed."
- That file holds `vision: {prepared: prepared, reference: vision_ref/qwen3-vl-2b.json}` [code:
  experiments/phase2/configs/r1_qwen3vl2b.yaml:9]. `train_sft.main` calls `preload(args.vision, args.model)` before
  any training [code: experiments/phase2/train_sft.py:292]. `reference()` raises when the JSON's model differs from
  the run's: `if (doc["model"], doc["set_sha256"]) != (model or doc["model"], set_sha256): raise ValueError(...)`
  [code: train_sft.py:163-164]. The file says `"model": "qwen3-vl-2b"` (it is written with `args.model`, train_sft.py
  :193). So a 4B SFT launched as written stops with `ValueError` after launch, sync and model download: one session
  overhead (30 min = $0.44 at r_s, front 6 inputs) for nothing. Preflight does not check `vision.reference` [code:
  experiments/phase2/preflight.py: no vision check among its checks]. 99_back decision 17 already says "Front
  section 3 step 0 lists only `model` and `output_dir` and must add it". The front still does not, so the two parts
  disagree.
- Front step 0 also says: "Train only the e2 branch (`branches: [2.0]`), so no branch selection runs and val stays
  unread." With `branches` set, `policy.chain` runs the default stages `["train", "select", "probe"]` [code:
  experiments/phase2/vast/policy.py:450]. `select` runs `select_ckpt` on `lidc-val-shifted`, `lidc-val-rescaled`
  and `lidc-val-rotated`, and `probe` runs `g_probe` on `lidc-val-shifted` [code: policy.py:441-448]. So val is read
  and the session pays for both stages.

Exact fix: front step 0 and 99_back decision 17 read the same text:
"Copy `r1_qwen3vl2b.yaml` to `experiments/phase2/configs/r1_qwen3vl4b.yaml` and change exactly: `model:
qwen3-vl-4b`; `output_dir: r1-qwen3vl4b`; `branches: [2.0]`; `stages: [train]` (no select, no probe: val stays
unread); `vision.reference: vision_ref/qwen3-vl-4b.json`. Before the rental, write that reference with
`python -I -m experiments.phase2.train_sft --config experiments/phase2/configs/r1_qwen3vl4b.yaml
--vision_reference_only true` (CPU, float32; about 17.8 GB of host RAM for 4.438 B parameters [estimate: 4.438 B x
4 B]). Add a preflight check: the `vision.reference` JSON's `model` equals the config's `model`." Add a unit test:
`preflight` refuses `r1_qwen3vl4b.yaml` while the reference names `qwen3-vl-2b`.

### Major

**3. [major] B (B2.1 quantile buffer; B2.4 rules 3 and 4; U68). The per-sequence `seq_kl` buffer is cleared every
logged step, but rules 3 and 4 need it over 20 and 250 steps, and the B1 run (once, in V1) does not save it.**

Problem.
- B2.1: "`log` computes q50, q99 and q99.5 of `seq_kl` with one `torch.cat` + `torch.quantile`, before the existing
  single `.tolist()`, and clears the list. The same buffer gives the live OPSM delta of B2.4."
- Our RL runs log every step: `over |= {"logging_steps": 1, ...}` [code: experiments/phase2/g_probe.py:88]. So the
  buffer holds one step.
- B2.4 rule 4: "delta = q99.5 of `seq_kl` (B2.1) over the first 20 steps of the run that uses OPSM ... from the
  B2.1 buffer". Guard: "recompute delta once as the q99.5 of the last 20 steps". U68: "q99.5 of the buffered
  per-sequence mean gaps of steps 1-20". A q99.5 over 20 steps cannot be rebuilt from 20 per-step q50/q99/q99.5
  values.
- B2.4 rule 4 pricing: "Priced at $0 on the B1 run: apply the steps 1-20 q99.5 of the B1 buffer to steps 50-300 as a
  counterfactual and report the share it would mask." Rule 3: "raise the starting delta to the q99.9 of `seq_kl`"
  over B1 steps 50-300. Both need per-sequence values of the B1 run. B2.1 logs only fixed-delta counts
  (`cf/opsm_share@{0.05,0.1,0.5}`) and per-step quantiles. B1 runs once, in V1, on the `v1-b1` tag that holds B2.1
  (front order rule 2). After V1 the inputs of rules 3 and 4 do not exist, and the engineer must guess.

Exact fix:
- B2.1: "`log` also writes the step's per-sequence `seq_kl`, advantage and valid-token count (256 x 3 float32) to
  `<out>/seq_gaps.npy` (append; 300 steps x 256 x 3 x 4 B = 0.9 MB). They are concatenated into the same device
  tensor as the stats, so the single `.tolist()` stays the only host read." Keep the per-step quantiles as log keys.
- B2.4 rule 4: "the callback keeps its own device list of steps 1-20 (and, for the guard, a ring buffer of the last
  20 steps), separate from the per-step list that `log` clears." Rules 3 and 4 pricing read `seq_gaps.npy` of the
  B1 run.
- U68: the test runs with `logging_steps=1` and asserts that the step-20 delta equals `np.quantile` over all 20 x 256
  planted values. Add R01b: `seq_gaps.npy` has 256 rows per step and the B1 tag writes it (P13 tree).

**4. [major] B (B2.4 rule 0 "Timing rule"; Gate B2; U51) and C (C5 paragraph 2, "in-process LoRA/base
generate-time ratio"). The base side of the timing rule runs on the LoRA engine, so the gate does not measure F's
decode cost.**

Problem. Rule 0: "after the first F sync, generate the same 32 x 8 batch twice, with the `lora_request` and with
`lora_request=None` (the base, same engine), in the order ABAB, seeded; pass if the median LoRA / base
generate-time ratio is <= 1.15". The engine is built with `enable_lora=True` (B2.2 F step 1). In vLLM 0.32.0 every
LoRA-wrapped linear runs `self.punica_wrapper.add_lora_linear(output, x, self.lora_a_stacked, self.lora_b_stacked,
...)` on every forward, whatever the request's adapter [web: vllm v0.32.0 lora/layers/base_linear.py:208-233]. So
the "base" side already pays the LoRA wrapper and kernel launches. The ratio measures only the adapter's compute,
not the cost of LoRA serving against today's merged engine, which is the 10-30 % decode risk that C5 and C.6 name
("vLLM's LoRA kernels may slow decode by 10-30 %"). Two more defects: with prefix caching on, the first A and the
first B generation are cold and the second of each is warm (the cache key holds the LoRA name), and with 2 repeats
per side the "median" is a mean of one cold and one warm run.

Exact fix: rule 0 timing, Gate B2 and U51: "Timing is a report, not a pass rule. F is a numerics fix (B2 intro), so
rule 0 decides on NaN/inf = 0 and the equality checks. Report two numbers: (1) the in-process ratio, labelled
'adapter compute only', with `llm.reset_prefix_cache()` before every generation and 5 repeats per side (ABABABABAB),
median; (2) the decode slope of the first 15 F steps, `time/generate_s` = a + b x longest answer (least squares, as
the 4.3 ms per token fit [measured fit: B_gpu_audit.md:40]), against the same fit on the B3 run's steps 1-15 (cross
host, so a report line). The owner sees both before the second F run." If the owner wants a hard timing gate, it
needs a merged arm on the same host (one extra process, 6 min + 15 steps), priced in front 6.

**5. [major] B (B2.2 F step 4 "Output-level check"; B2.3b box part; B2.4 rule 0) and 00_front (6.3 V3 item 2a).
The F output-level check compares against "the B1 run's (f) against (b)", which is never measured on the expected
path.**

Problem. B2.2 F step 4: "Output-level check, independent of vLLM internals: under F, way (f) against way (b) of
B2.3b on the 256-answer batch has a mean |gap| at most half of the B1 run's (f) against (b). This is the number F
exists to fix." The B1 (f)-vs-(b) number comes only from the B2.3b box part, which runs "only if B2.4 rule 0 does
not already pick `lora-seqmask` or `fp16-seqmask`" (B2.3b) and as V3 item 2a "only if B2.4 rule 0 did not pick F
or A" (front 6.3). The B2.0 seed value is 0.531 >= 0.2, so rule 0 fires and the comparator is skipped exactly when
the check is needed. The F check then has a pass rule with no input.

Exact fix: B2.2 F step 4 and B2.4 rule 0: "The comparator is measured in the F run itself, before step 1, in the
same process and on the same host: load the B1 final adapter (HF revision `b1-rows`); push it merged (`values=
"merged"`, C5) and generate the 32 x 8 batch, keep (f_merged); then `sync_lora` the same adapter and generate the
same batch with the same seed, keep (f_lora); score both with way (b) on the trainer. Pass: mean |f_lora - b| <= 0.5
x mean |f_merged - b|." Cost: two generations and two scoring passes, about 2 min inside the run [estimate]. V3 item
2a then keeps its own purpose (attribution, ADR-0013) and is not an input of rule 0.

**6. [major] C (C7 "Files"; C.6 G-C row; C11 step 5) and 00_front (section 1 "C-G-C (<= 18 s) needs C7 and C8 and
is measured inside Phase B") against B (B3.9). Whether the B3 run uses the new outline limit L is not stated, and it
changes the B3 A/B.**

Problem. C7 "Files": "`experiments/phase2/configs/rl_shifted.yaml` (or the B config) adds `--grpo
'{"max_completion_length": L}'`". C.6 G-C: "With C7 and C8 added (owner decisions, inside Phase B): mean step <=
18.0 s". C11 step 5: "Arm (inside a Phase B run, as C7; owner approval)". The only B run after C adoption is B3 (V1
has no C code, front order rule 2). B3.9 lists every B3 change against B1: "the reward changes, plus the truncation
mask ..., plus `beta` and `adam_epsilon` ... The serving path and the IS mode do not change." L is not in the list,
and B1 runs at 1,536 [code: configs.py:73 `MAX_COMPLETION ... "qwen3-vl-2b": 1536`]. At L = 1,200 (Direction's
recommendation, C7 step 1), naturally ended 1,200-1,493-token answers become truncated rows: "6 of 15 steps had a
naturally ended answer longer than 1,200 tokens" [measured: C7 Facts, runs/r5/shared_rl.json]. Under B3's unmasked
truncation (B3.5), those rows get R = 0 and A = -b. dr_grpo also rescales the loss by 1,536 / L [code:
grpo_trainer.py:3286]. So an L change in B3 is a fifth change that B3.10 would credit to "v5 + unmasked truncation".

Exact fix: B3.9: "B3 keeps `max_completion_length` 1,536 and `--runaway-cap 0`, as B1." C7 "Files", C11 step 5,
C.6 G-C and front section 1: "L (C7) and the runaway cap (C11) apply from the first run after Gate B3; G-C and
G-RC2 are measured there, not in B3." Add a config test: `b3_v5.yaml` has no `max_completion_length` and no
`--runaway-cap`.

### Minor

**7. [minor] B (B3.5 "Truncated rows"; Gate B3 report lines). Under `bf16-seqmask`, the sequence IS ratio of an
unmasked 1,536-token row sums 1,536 token gaps, so the IS mask may drop part of the new negative signal.**
Problem. TRL's sequence ratio is `exp(sum_t (old - sampling) * mask)` and `sequence_mask` zeroes ratios above 3
[code: grpo_trainer.py:2779-2814]. In B1 truncated rows are masked, so their ratio is exp(0) = 1 and no data exists
on them. With independent token gaps of sd about 0.025 [estimate from the median |gap| 0.018], a 1,536-token sum has
sd about 1.0, so about 14 % of such rows exceed ln 3 and get weight 0, and others get weights far below 1 [estimate].
Runaway text may differ more. Fix: B2.1 adds `sampling/trunc_is_masked_share` and `sampling/trunc_weight_mean`
(truncated rows only); B3.5 prints them; Gate B3 reports them beside `truncated_grad_share`.

**8. [minor] B (B2.4 rules 3 and 4). Rule 3's raised starting delta is overwritten by rule 4 at step 20.**
Problem. Rule 3: "if `cf/opsm_share@0.5` > 2 % of sequences, raise the starting delta to the q99.9 of `seq_kl`".
Rule 4: "Until step 20, delta = 0.5 ... at step 20 sets `trainer.off_policy_mask_threshold = delta`" (q99.5 of the
run's own steps 1-20). So the veto acts for 20 steps only, and rule 4 may then set a lower value. Fix: rule 3 reads
"... raise the steps 1-20 delta to the B1 q99.9, and at step 20 use max(rule 4 value, B1 q99.9)"; U51 and U68 test
it.

**9. [minor] C (C.6 G-SP1 "Under F" thresholds; C.8.1 F rows). After item 1, the F thresholds of G-SP1 have no run
that gates on them.** Fix: mark the F rows of G-SP1 and C.8.1 "report only, first F run after Gate B3", or delete
them.

**10. [minor] 00_front (section 3 reason 2). The EGM numbers (84.2 / 87.3 / 89.6 / 91.0) could not be re-checked
through the HF papers API (2601.13633 is not indexed).** Fix: cite the arXiv HTML table URL beside the paper id so
DOC4 can re-check it.

## Counts

Blocking 2 (items 1, 2). Major 4 (items 3-6). Minor 4 (items 7-10).
