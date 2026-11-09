# Final review round 13: RL (B2, B3, C, model choice)

Reviewer role: RL for LLM/VLM research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Scope read in the sources: B.md B.1 items 5 and 8, B1.5, B2.0-B2.4, B3.1-B3.10, B.6.2 (Gates B2, B3), B.9 U54/U73;
C.md C7, C8, C9, C10, C11, C.6, C.9.1 (U-RC1..RC6), C.9.4 (G-RC0, G-GR1, "Model of the two vLLM box tests");
00_front sections 3 and 6.3 (V3 table, V3 amount rule); review-log.md round 12 (all phases) and
`reviews/final_r12_rl.md`. The round-12 change in my area (qa-2: G-RC0 as the last, soft step of `b3_v5.yaml`, its
model and U-RC6) was checked against the code and vLLM 0.32.0 and holds (below).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL GRPO, sequence-level IS or vLLM LoRA serving. No item depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:351-362` (`is_ratio`: sequence modes sum the token gaps), `:581-602` (`set_is_ratio`), `src/objdel/train/rewards.py:60-76` (truncation rule), `src/objdel/train/configs.py:1-30, 70-100, 180-196`. Gave item 1. |
| engineering:debug | Traced a truncated row through TRL 1.17.1 under `mask_truncated_completions=False` and `sequence_mask`: summed log gap over all 1,536 tokens -> ratio -> `per_token_loss * importance_sampling_ratio`. Gave item 1. |
| engineering:system-design | Checked that the B3 A/B changes only what B3.10 credits: the credited "unmasked truncation" acts through rows whose IS weight the doc never logs (item 1). |
| engineering:architecture | Checked front section 3 (model decision) against the code: the 2B base is `Qwen/Qwen3-VL-2B-Instruct` [code: configs.py:73-77], and the 4B arm is the Instruct base too; EGM's Instruct-base gap is quoted. No item. |
| engineering:documentation | Each item names the source file and heading and quotes the doc line and the code line or measurement. |
| huggingface-skills:huggingface-best | Re-checked the model table: round 12 re-ran the HF API on 2026-11-09 (no new Qwen VLM up to 4B); front section 3 still parks Qwen3.5-2B/4B for the hybrid-attention reason and the E5 guard. No model-choice item. |
| huggingface-skills:trl-training | Read TRL 1.17.1 in `.venv`: `grpo_trainer.py:2397-2399` (`clipped_ratio` logged whatever the mask flag), `:2603-2611` (truncation zeroes `completion_mask`), `:2776-2824` (IS modes; sequence modes sum the token gaps), `:2870-2910` (RLOO passes through `sum_then_normalize` + `scale_rewards="none"`), `:3062-3092` (OPSM uses the mean gap), `:3186-3190`, `:3262-3290` (IS weight on the policy term only; KL added after; dr_grpo normaliser). |
| huggingface-skills:huggingface-papers | No new paper claim in scope. The item-1 argument is the standard identity E_q[log p/q] = -KL(q||p) (no paper needed). |
| Web (WebSearch/WebFetch loaded; raw GitHub at tag v0.32.0) | `vllm/multimodal/encoder_budget.py:20-40, 150-170`: a small `max_model_len` only lowers the items per prompt (`max(1, min(mm_limit, max_model_len // max_tokens_per_item))`); no start-up error. So G-RC0's `vllm.LLM(..., max_model_len=256)` on the Qwen3-VL checkpoint is fine. |

## Verified, no item

- Round-12 qa-2: `[~, -m, pytest, ...]` runs as `/opt/venv/bin/python -I -m pytest ...` with `|| echo` (soft)
  [code: experiments/phase2/vast/policy.py:484-487]. The model folder `runs/rl/train_rl/final` holds a merged
  `model.safetensors` (4,255,140,312 B), `config.json` (`Qwen3VLForConditionalGeneration`), tokenizer, processor
  and chat template [measured: ls, 2026-11-09], so `vllm.LLM(model={rd}/ckpt)` loads it by path.
- V3 arithmetic: item 4 = 3 + 5 + 13 + 6 + 260 x 27.5 / 60 + 10.4 = 156.6 min; step 260 ends at minute 185.2;
  `BudgetStop` = 0.10 + 185.2 x 0.01133 = $2.20; tail $0.35 / $0.01133 = 30.9 min >= 10.4 + 2 + 3. Holds.
- Logp-gap numbers of B2.2 F step 4: step 1 0.009563 / 0.009515 / 0.009556; full run medians 0.0176 (240-280),
  0.0236 (101-300), 0.0185 (301-500); steps 1-20 0.0110, 21-100 0.0278; grad norm 0.0018-0.0221, median 0.00475
  [measured: runs/r5/base/rl.json, runs/r5/shared_rl.json, runs/rl/train_rl/rl.json, recomputed 2026-11-09].
- C11 "5 of R5 steps 3-15 hit the 1,536 limit": 5 [measured: runs/r5/shared_rl.json `completions/max_length`].
- TRL line citations used by B2/B3/C (882, 2823, 3078, 3187-3190, 3286, 3416-3420, 3425) match the installed 1.17.1.
- RLOO through TRL, the legality gate on `r_topo`, the Adam/eps scaling argument (g / r with eps / r; KL scaled by
  the same beta / r) and the dr_grpo 1,536 / L argument hold, as round 12 found.

## Items

### 1. MAJOR (B, C): under `sequence_mask`, an unmasked 1,536-token truncated row gets an IS weight near 0, so B3's "truncated row gets A = -b" signal is mostly re-masked, and `truncated_grad_share` ignores the weight

- Sources and quotes:
  - B.md "B.1. Goal and why" item 5: "A sequence ratio is exp(sum of token gaps), so long answers pass the cap
    first."
  - B.md "B3.5 Wire v5 into the TRL adapter": "Truncated rows: the B3 config sets
    `mask_truncated_completions=False` (B.1 item 8). A truncated row keeps R = 0, its eval score, and gets A = -b
    (minus the group baseline)." and "B3.5 logs `completions/truncated_grad_share` = sum over truncated rows of |A|
    x tokens / sum over all rows of |A| x tokens. ... If the median exceeds 10 % over steps 50-260, the owner may
    switch on C11's option".
  - B.md "B3.9 B3 run": B3 runs "B1's numerics preset `--numerics bf16-seqmask`"; B.6.2 Gate B3: "the gate credits
    'v5 + unmasked truncation'".
  - B.md B.9 U73: the hand value 0.626 uses |A| x tokens only.
  - C.md "C11" step 4: "A capped answer is cut mid-text, so it is unparseable and gets R = 0 and a negative
    advantage (-b). Masking it would bias the update toward more capping".
- Code:
  - `fast.py:355-358` (`is_ratio`): `diff = torch.nan_to_num((logps - sampling) * mask, nan=0.0)`; `if
    mode.startswith("sequence"): diff = diff.sum(-1, keepdim=True)`; `r = torch.exp(diff)`; mask mode zeroes only
    `r > hi` (clip_min None -> `lo = -inf`), so a very negative sum is kept with a tiny weight.
  - TRL `grpo_trainer.py:3268-3269`: `per_token_loss = per_token_loss * inputs["importance_sampling_ratio"]`
    (the policy term; the KL term is added after, unweighted).
  - TRL `grpo_trainer.py:2606-2608`: with the mask flag on, a truncated row's `completion_mask` is zeroed, so its
    summed gap is 0 and its ratio is 1. With the flag off (B3), the sum runs over all 1,536 tokens.
- Why the weight is small, not only random: tokens are sampled from vLLM, so E_vLLM[log pi_trainer - log pi_vLLM]
  = -KL(pi_vLLM || pi_trainer) <= 0 per token. The summed log weight of a row therefore drifts down linearly with
  its length. Measurement: the median step-mean sequence ratio of the full run is 0.823 at a mean answer length of
  105.5 tokens [measured: runs/rl/train_rl/rl.json `sampling/importance_sampling_ratio/mean`,
  `completions/mean_length`, 2026-11-09]. By Jensen, E[log w] <= log 0.823 = -0.195 at 105 tokens, about -0.0018
  per token, so a 1,536-token row has E[log w] <= about -2.8, w about 0.06 [estimate: linear in length]. The 12x
  token mass that B3.5 cites (1,536 / 131) is then cut by about 16x in the loss.
- Consequences: (a) the B3 mechanism ("a direct signal against runaways") largely does not reach the loss for the
  rows it targets, and capped rows under C11 are hit the same way (shorter, so less); (b) `truncated_grad_share`
  as defined overstates the truncated rows' share of the gradient by about the inverse weight, so the 10 % owner
  switch reads the wrong quantity; (c) ADR-0014's credit sentence names a mechanism nobody measured; (d) the runs
  after Gate B3 default to token-level presets (`bf16-tis-opsm`, or F), where a long row's per-token weights stay
  near 1, so B3's unmasking result does not carry over one to one. No planned check covers it: B1 runs with the
  mask flag on (ratio 1 for every truncated row, and TRL also drops those tokens from attention, so V1 cannot
  price it), U73 plants no weight, and Gate B3 reads only outcome lines (`clipped_ratio`).
- Fix:
  - B.md B3.5: "`completions/truncated_grad_share` = sum over truncated rows and their valid tokens of |A_i| x w_it
    / the same sum over all rows, where w is the `importance_sampling_ratio` the loss applies (a sequence weight
    broadcast to its tokens). B3.5 also logs `sampling/trunc_logw_median` (median summed log gap of the truncated
    rows of the step) and `sampling/trunc_w_lt_0.1_share` (truncated rows with w < 0.1, masked rows included). The
    run report prints their medians over steps 50-260." Update U73: plant w = 0.06 on the truncated row (hand
    value 0.3 x 0.06 x 1,536 / (0.3 x 0.06 x 1,536 + 7 x 0.3 x 131) = 0.091) and w = 1 elsewhere.
  - B.md B.6.2 Gate B3 "report only": add "the median `sampling/trunc_w_lt_0.1_share` of B3". ADR-0014 (B3.10
    "What is credited"): "v5 + unmasked truncation on `bf16-seqmask`, with the measured truncated-row weight
    (median `trunc_logw_median`); under a token-level preset the truncated rows keep weights near 1, so the first
    run after Gate B3 re-reads `truncated_grad_share`."
  - C.md C11 step 4: add one sentence: "Under a sequence-level IS mode a capped row's weight is exp(sum of its
    token gaps) and falls with its length (B3.5 logs it); the G-RC2 arm reports `trunc_w_lt_0.1_share` for capped
    rows."

### 2. MINOR (A, B, C): the three round-12 RL minors were neither applied nor logged

- Sources: review-log.md has no "rl-" row under "Review log (final round 12)" in any phase, and no
  `final_r12_rl_verdict.md` exists. The doc still reads: B.md "B2.2" F step 2.4 and C.md "C5" F sync step 4 "the
  prefix-cache key holds the LoRA name, not its weights" (vLLM 0.32.0 keys on name and path,
  `v1/core/kv_cache_utils.py:578-594`); A.md "A7" has no `_llm_weights_sleeping` entry (`grep` finds none) while C5
  writes that private TRL attribute (`trl/generation/vllm_generation.py:379, 503, 773`); B.md "B2.1" gives `ess`
  no unit, while `set_is_ratio` counts sequence ratios per sequence [code: fast.py:594].
- Problem: same as `reviews/final_r12_rl.md` items 1-3; none changes a result (U-SY9 and `make test` catch the
  first, the second's conclusion holds, the third is a report line).
- Fix: apply `final_r12_rl.md` items 1-3 as written, or add a round-12 log row that rejects them with a reason.

## Counts

blocking 0, major 1 (item 1), minor 1 (item 2).
