# Final review round 9: RL (B2, B3, C, model choice)

Reviewer role: RL for LLM/VLM research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Scope read: B.md B.1 items 5 and 8, B1.5, B2.0-B2.4, B3.1-B3.10, B.6.2 (Gates B2, B3); C.md C5-C8, C11, C.6, C.8.1;
00_front sections 1, 2.5 and 3; 99_back "RL settings" and "Model"; the B, C and front review logs of rounds 7 and 8.
Round 8's RL file had no verdict and was not a source of the round-8 fixer, so its three minors are still open in
the text. Items 4-6 carry them, with new evidence where I have it.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL/vLLM LoRA serving, GRPO or IS correction. No item depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:350-363` (`is_ratio`), `:407-415` (`one_pad_per_image`), `:500-620` (`skip_old_logps`, `_compute_loss`, `set_is_ratio`, `log`), `:618-622` (the patch goes in lazily in `_generate`), `:670-690` (`patch_vllm_prompts`), `:693-717` (`lora_merged_params`). Read `experiments/phase2/train_rl.py:52-60` (`merge_and_unload()` before `save_pretrained(final)`), `final_eval.py:48` (only `--ckpt`), and `src/objdel/reward/{advantages,spec,score}.py`. These reads gave items 1 and 2. |
| engineering:debug | Traced the F comparator before step 1. `patch_vllm_prompts` runs only inside `_generate`, so before step 1 `llm.generate` is vLLM's own function. The doc says the comparator bypasses the wrapper, and no other place converts the prompts to one pad per image (item 1). |
| engineering:system-design | Traced the life of an F-trained adapter: sampled unmerged in training, then merged in bf16 at save, then evaluated merged (item 2). Re-sized the fixed-length timing batch against vLLM's KV budget, using a measured bytes-per-token value (item 4). |
| engineering:architecture | Checked that B2.2 4b, C5 step 5(c), G-EQ5, U70 and U76 agree on the comparator. They agree on order and adapter names. None of them defines the prompts (item 1). Checked that B2.4 rules 3 and 4 both set the starting OPSM delta (item 8). |
| engineering:documentation | Each item has a source file, a heading, the quoted doc line, and a code line or measurement. |
| huggingface-skills:huggingface-best | Re-read front section 3 against the owner rule ("any model up to 4B") and the 4B SFT recipe in `r1_qwen3vl2b.yaml` and `configs.py:73-77`. Model ranking and size facts: no item. |
| huggingface-skills:trl-training | Read TRL 1.17.1 in the repo .venv: `grpo_config.py` field names (`vllm_importance_sampling_clip_max/min`, `off_policy_mask_threshold`, `mask_truncated_completions`), `grpo_trainer.py:882, 1029-1036, 1186-1188, 1203, 1838-1841, 2603-2611, 2750-2824, 2880-2924, 3066-3089, 3180-3195, 3255-3290`, and `generation/vllm_generation.py:356-379, 494-528, 718-776`. The presets build, OPSM's sign matches B2.1, the KL term comes after the IS weight, and `scale_rewards="none"` only centres. So RLOO passes through unchanged, as B3.4 says. |
| huggingface-skills:huggingface-papers | No new paper claim in scope. The front section 3 EGM numbers (2601.13633) were verified in round 7; I did not re-fetch them. |
| Web | Fetched vllm PR 13968: an fp16 overflow in the Qwen2.5-VL ViT. Confirmed as B2.2 fp16 step 0 cites it. Searched for vLLM LoRA-vs-base decode overhead on small models: no benchmark found, so item 5 stays minor. |

## Verified, no item

- Preset field names exist in TRL 1.17.1 `GRPOConfig` (`vllm_importance_sampling_clip_max` / `_clip_min` at
  grpo_config.py:952-965; `off_policy_mask_threshold` at :968). `_bf16()` plus the explicit `"bf16": False` gives
  fp16 without bf16 [code: configs.py:140-142].
- OPSM: `kl_div = sampling - per_token_logps.detach()`, keep if A >= 0 or mean <= delta [code: grpo_trainer.py:
  3078-3089]. TRL reads `self.off_policy_mask_threshold` at each loss call [:882, :3190]. FastGRPO's skip path leaves
  `sampling_per_token_logps` in `inputs` (set_is_ratio raises otherwise [code: fast.py:584-586]).
- RLOO identity `R_i - (S - R_i)/(G-1) = G/(G-1) (R_i - mean)`. The B3.4 worked check gives
  0.1 x (0.62 - 0.374/7) = 0.0567. `topology_valid` requires `legal` [code: advantages.py:79-82], so an illegal
  answer has r_topo = 0 and R = 0 under v5's `legal * (...) + c_topo * r_topo`.
- KL reference with PEFT is the adapter-off model at any beta [code: grpo_trainer.py:1029-1036]. `beta*per_token_kl`
  is added after the OPSM mask and the IS weight [:3266-3272]. So the (beta/r, eps/r) argument of B2.4 / B3.4 holds.
- `max_num_seqs` = 16 x 1 x 16 = 256 [code: grpo_trainer.py:1186-1188], as C6 step 2 says.

## Items

### 1. MAJOR (B; C follows): the F comparator calls vLLM's own `generate`, but the doc never says which 32 prompts it uses or how to cut them to one image pad. With the trainer's expanded ids, vLLM expands the image a second time.

- Sources: B.md "B2.2 Numerics presets", Option F step 4 and step 4b.3; C.md "C5", step 5(c) and G-EQ5 (they point to
  B2.2 4b); B.md B.9.1 U70/U76.
- Doc lines:
  - B2.2 step 4: "(i) It pushes the adapter merged ..., generates the 32 x 8 batch with a fixed seed and keeps the
    generation-time logps (f_merged)". "The 32 x 8 batch" is not defined anywhere in B2.2, B2.4 or C5. Grep for
    "32 x 8 batch" in B.md: lines 1245 and 1409 only, and neither defines the batch.
  - B2.2 4b.3: "The comparator calls `llm.generate` directly, not through the `patch_vllm_prompts` wrapper: arm
    (ii) with `LoRARequest("b1", 2, <b1 dir>)`, arm (i) with `lora_request=None`."
- Proof (code):
  - The wrapper is the only place that cuts each run of image pads to one pad before vLLM sees it. It is also the
    only place that checks that vLLM's own expansion gives the trainer's ids:
    `outputs = generate(one_pad_per_image(prompts, pad), ...)` and `raise RuntimeError("vLLM expanded the image to
    {n} pads, the trainer's processor to {m} ...")` [code: src/objdel/train/fast.py:680-688; helper :407-415].
  - The wrapper goes in only on the first `_generate` [code: fast.py:618-622]. The comparator runs before step 1, so
    at that point `llm.generate` is vLLM's own, unwrapped function.
  - The trainer scores answers on processor-expanded ids [code: fast.py:671-674 docstring, "TRL trains on its own
    processor-expanded ids (grpo_trainer.py 1844)"]. TRL 1.17.1 sends those expanded ids plus the image to vLLM
    unchanged [code: .venv trl/generation/vllm_generation.py:720-743]. vLLM 0.32 then replaces the first pad with N
    pads and keeps the other N-1 pads as text tokens. The source note behind the wrapper records this (N -> 2N-1
    pads; trl issue #6294 measured 855 -> 1672 tokens on SmolVLM) [code: commit f68b32d, the added research note
    §3.2].
  - So an engineer who follows 4b.3 and feeds the comparator the prompts the trainer scores gets vLLM sampling from
    a broken prompt in all four arms. Then f0, f_merged and f_lora measure prompt corruption, not merge rounding.
    The likely result is `informative: false` (an F "pass" with no evidence), or a crash on a placeholder count. A
    pass or fail of the F check is what decides the preset of every run after Gate B3 (B2.4 rule 0).
  - U70 and U76 use a fake engine, so neither test can catch this.
- Fix:
  - B2.2 4b.3: "The comparator calls vLLM's own `generate` (the wrapper goes in only at the first `_generate`,
    fast.py:618-622). It passes `one_pad_per_image(prompts, pad)` and runs the same expansion check. Move the body of
    `one_pad` into `fast.vllm_generate_checked(generate, prompts, pad, **kw)`; the wrapper and the comparator both
    call it. The only difference is `lora_request`: the wrapper adds `"policy"`, and the comparator passes `"b1"` or
    `None`."
  - Define the batch in B2.2 step 4: "the 32 prompts of step 1 of the run's rows schedule (rows 0-31 of the
    `--train` root, file order, B1.5), each repeated 8 times (`SamplingParams(n=1)`, 256 rows, seed s = 0 for
    arms (0), (i), (ii)). Prompt ids and images come from TRL's own `_tokenize_prompts` path (ImageMemo on), so the
    trainer's way (b) scoring reads the same expanded ids."
  - U76: a fake `generate` that raises on a prompt with more than one pad per image. The comparator must pass the
    check. A fake that returns N+1 pads must fail the F check with a `prompt_expansion` flag.
  - C5 step 5(c): add "(prompts as B2.2 4b.3)".

### 2. MAJOR (B; front 3 and 99_back 16 follow): an F-trained run is evaluated on its bf16-merged `final/`, which brings back the merge rounding F exists to remove

- Sources: B.md "B1.5 Train on the ordered rows" ("B1 adapters to V3"), "B2.4 Decide" rule 0 and the ADR-0013 text;
  99_back "RL settings" 16; 00_front section 3 (the 4B gate reads val IoU of "the last 2B run").
- Doc lines:
  - B1.5: "`hf_push --rl <train_rl out>`: it uploads `final/` (the merged model, the form `fetch_hf.sh` and
    `final_eval --ckpt` read) and `final_adapter/`".
  - B2.4 rule 0: "if B2.0 shows a bf16 relative error >= 0.2 at step 100 (the seed run gives 0.531), the first run
    after Gate B3 starts under `lora-seqmask`". No line in B, C, front or back says how that run is evaluated. Grep
    for "lora-seqmask" or "option F" together with eval, final or serve: 0 hits.
- Proof:
  - `train_rl` merges before it saves: `model = model.merge_and_unload() ...; model.save_pretrained(a.out / "final")`
    [code: experiments/phase2/train_rl.py:58-59]. PEFT's merge adds the delta into the bf16 base in place. This is the
    same `round(W + d)` that B2.0 measures.
  - `final_eval` takes only `--ckpt` [code: experiments/phase2/final_eval.py:48]. Nothing in the repo starts vLLM
    with `enable_lora` (grep `enable_lora` in src and experiments: 0 hits).
  - Measured size of the loss: bf16 relative error ||round(W+d) - W - d|| / ||d|| = 0.531 at step 100 and 0.426 at
    step 300 [measured: B.1 item 5 table, B2.0 seed script]. A 260-step F run would be evaluated on a policy whose
    delta is off by about 45 % of its norm [estimate: interpolation]. Training never sampled that policy, because
    under F vLLM samples W + BA unmerged.
  - Under the merged presets (B1, B3), training sampled `round(W + d)` too, so their evals match training. Only F
    breaks this. So the F run's val numbers, the 4B gate and the final test candidate would score a policy other
    than the one trained.
- Fix:
  - B2.4 (after rule 0) and ADR-0013: "A run trained under `lora-seqmask` is evaluated with vLLM LoRA serving:
    `final_eval --ckpt <its --ckpt base> --adapter <out>/final_adapter`. `final_eval` builds `LLM(...,
    enable_lora=True, max_lora_rank=16)` through `fast.vllm_kwargs` and passes `LoRARequest("eval", 1, <dir>)` to
    every generate (prompts via `vllm_generate_checked`, item 1). `final/` is still saved, for the merged report line
    only."
  - The F comparator already measures this gap on the B1 adapter (f_merged against f_lora). Write it into ADR-0013
    as the expected merged-eval loss.
  - Tests: an I-test in which `final_eval --adapter` on a fake engine passes `enable_lora` and one `LoRARequest` per
    generate; refuse `--adapter` without `adapter_config.json`.
  - 99_back 16: add "an F run is evaluated unmerged (B2.4)".

### 3. MINOR (B): B3 unmasks truncated rows, but the B3 preset's sequence IS mask is most likely to drop exactly those rows, and no log separates the two effects

- Sources: B.md "B.1 Goal and why" item 5 and item 8; "B3.5 Wire v5" (truncated rows, `truncated_grad_share`);
  "B3.9 B3 run" (`--numerics bf16-seqmask`); Gate B3 ("the gate credits 'v5 + unmasked truncation'").
- Doc lines:
  - B.1 item 5: "A sequence ratio is exp(sum of token gaps), so long answers pass the cap first."
  - B3.5: "A truncated row keeps R = 0, its eval score, and gets A = -b".
- Code: under `mask_truncated_completions=True` (B1), a truncated row has an all-zero mask, so its summed gap is 0
  and its ratio is exactly 1 [code: fast.py:355-358; grpo_trainer.py:2603-2611]. Under B3's `False`, the same row
  sums 1,536 token gaps, then `sequence_mask` zeroes it above 3 or shrinks it toward 0 [code: grpo_trainer.py:
  2784-2818]. `truncated_grad_share` (|A| x tokens) ignores that weight. No measurement exists of how many
  truncated rows the mask drops, so this item is minor.
- Fix:
  - B2.1: add `sampling/trunc_is_masked_share` (truncated rows with weight 0 / truncated rows) and
    `sampling/trunc_is_weight_mean` to the same stats tensor. One more truncated flag per row; no new host read.
  - B3.5: compute `truncated_grad_share` from |A| x IS weight x tokens.
  - Gate B3 report: print both. If more than half of the truncated rows are masked, ADR-0014 states that the
    unmasking reached fewer than half of the rows it targets.

### 4. MINOR (B, carried from round 8 rl-1, not processed): the fixed-length timing batch needs about 3x the KV cache that vLLM has, so it runs with preemption, not "at batch 256"

- Source: B.md "B2.2", Option F step 4b.4 ("Limit: a full-batch 512-token decode measures the kernel cost at batch
  256") and step 4 ("14 generations of about 7.8 s each").
- New measured anchor: 14.69 GiB of KV holds 137,488 tokens, so about 112 KiB per token [measured:
  `runs/*/env_check_vllm.log`, "Available KV cache memory: 14.69 GiB" and "GPU KV cache size: 137,488 tokens"].
- Arithmetic:
  - At `vllm_gpu_memory_utilization` 0.35 [code: configs.py:192] of 31.84 GiB, vLLM has 11.1 GiB in total. Take
    away 4.54 GiB of weights and non-torch memory, 0.78 GiB peak activation and 0.2 GiB CUDA graphs [measured:
    same logs]. That leaves about 5.6 GiB, or about 52k tokens [estimate].
  - One timed generate needs 32 x 325 + 256 x 512 = 141,472 tokens [prompt length measured: runs/r5/shared_rl.json].
  - Under F, the LoRA buffers take more memory still.
- Fix (as round 8):
  - Size T from the live engine: `kv = cache_config.num_gpu_blocks * block_size`;
    `T = min(512, 16 * floor((0.9 kv - 32 x prompt_len) / (256 x 16)))`; require T >= 64.
  - Record T, kv and the preemption delta (must be 0) in `f_check.json`.
  - Rewrite the "Limit" sentence and the cost line for T.

### 5. MINOR (B, C; carried from round 8 rl-2): the 1.15 timing bar is inside the doc's own LoRA decode estimate, so F can fail for speed alone

- Sources: B.md "B2.4 Decide" rule 0 ("pass if the median of the 6 LoRA / base generate-time ratios is <= 1.15
  [estimate threshold]"); C.md "C5" and "C.6" ("vLLM's LoRA kernels may slow decode by 10-30 %").
- New evidence [estimate]:
  - vLLM wraps 4 linear groups per layer (qkv, o, gate_up, down). Each needs a shrink and an expand kernel, so a
    decode step gets 224 extra kernels over 28 layers.
  - The measured decode cost is 4.3 ms per token of the longest answer [measured: B_gpu_audit.md:40].
  - At 5-15 us per small Triton kernel, that is 1.1-3.4 ms per step, which can reach 25-80 %.
  - No benchmark exists for this model (web search, item list above), so this item stays minor.
- Fix: a timing failure alone (every other F flag true) keeps F and writes the ratio into ADR-0013 as an owner line.
  If the owner wants a hard bar, he sets it from the same step-time budget as G-C (C.6), before the run.

### 6. MINOR (B, carried from round 8 rl-3): B.1 item 5 reads an IS ratio of 0 as "masked", but an fp32 underflow also gives 0

- Source: B.md "B.1 Goal and why" item 5: "`importance_sampling_ratio/min` is 0 in 98.1 % of 770 steps, so at least
  1 of 256 sequences was masked in those steps." A.md says that A6 splits the two cases.
- Code: `torch.exp(diff)` of a summed gap below about -103.97 is 0 in fp32 [code: fast.py:358].
- Fix: "... was masked (ratio > 3) or underflowed to 0 (A6 counts the two apart)".

### 7. MINOR (C): G-C's estimate leaves out F, but the run that measures G-C is expected to run under F

- Sources: C.md "C.6" G-C row ("estimate 16.3-19.4 s (A) or 17.0-20.1 s (B) (§8.1)"); "C.8.1" row "C1-C5 + C7 + C8"
  (built from the merged C1-C5 savings 2.25-4.7 s); C.6 "Numerics first" ("Under F that run reports its step time
  ... vLLM's LoRA kernels may slow decode by 10-30 %, that is 1.1-3.3 s").
- G-C is measured "in the first run after Gate B3", and B2.4 rule 0 expects that run under `lora-seqmask` (seed
  error 0.531). With C.8.1's own F row (1.47-3.97 s) and the 1.1-3.3 s decode cost, the estimate becomes about
  17.1-23.4 s [estimate: arithmetic]. So the 18 s bar is unlikely to be met. No action follows from a G-C fail.
- Fix:
  - C.8.1: add a row "C1-C5 (F) + C7 + C8 + measured `lora_over_base`".
  - G-C row: "Under F, G-C is a report line against that row."

### 8. MINOR (B): rules 3 and 4 of B2.4 both set the starting OPSM delta

- Source: B.md "B2.4 Decide" rule 3 ("raise the starting delta to the q99.9 of `seq_kl` ... B1 steps 50-260") and
  rule 4 ("Until step 20, delta = 0.5 (TRL docstring)").
- One reading keeps 0.5 for steps 1-20. The other uses rule 3's value. The effect is limited to 20 steps, because
  rule 4 replaces delta at step 20.
- Fix: rule 4 reads "Until step 20, delta = the start value: 0.5, or rule 3's raised value when its veto fires".
  U68 asserts the raised start value on a planted veto.

## Counts

blocking 0, major 2, minor 6.
