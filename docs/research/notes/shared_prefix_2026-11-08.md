# Shared-prefix forward, 2026-11-08

Item #1 of `rl_code_optimality_review_2026-11-08.md`: compute each prompt once per group of G = 8 answers in the
training forwards (old-logprob, reference, update). Code: `FastForward.packed_chunk` in `src/objdel/train/fast.py`,
turned on by `train_rl --shared-prefix` / `g_probe --shared-prefix`. Off by default until the GPU test (R5).

## Layout

TRL scores each answer as its own row, `[prompt | answer]`, so the prompt (about 325 real tokens, mostly image
tokens) runs 8 times per group in every forward and backward. The packed row of one group is

    [prompt minus its last token] [p_last a1 ...] [p_last a2 ...] ... [p_last a8 ...]

- The last prompt token is repeated at the head of each answer block, with its own position. Its hidden state is
  then the same as in TRL's row (same token, position and visible keys), and it predicts the answer's first token.
  TRL's fused LM head scores every block with its usual one-token shift; `support_ids` (the sampling mask) is placed
  on the same shift.
- A 4-D boolean mask from segment numbers (0 = prompt, k = block k, -1 = padding): the prompt attends causally, each
  block attends to the prompt and causally to itself, never to another answer. Padding rows see only themselves.
  SDPA and eager attention are both supported; flash attention is refused.
- Positions are each row's own 3-D mrope positions from `RopeCache`, sliced, so no position is recomputed from the
  4-D mask. A model without the RopeCache hook (a separate reference model) takes the per-row path.
- The vision tower gets one image per group (the 8 rows share it).
- Results are scattered back into TRL's `(rows, answer width)` layout; positions off the answers are 0.

## When it packs

- Update micro-batches must hold whole groups: `by_group_length` keeps each group's 8 rows together (groups longest
  first). Under `dr_grpo` with every micro-batch in one optimizer step the partition does not change the update
  (the `by_length` argument). With `top_entropy_quantile < 1` (a per-micro-batch cut) TRL's shuffle is kept.
- Cost check per chunk (`packing_wins`): SDPA's masked kernel skips no tiles, so packed attention is dense over the
  packed row. Cost per row of n tokens = n + n^2 / 12,000 (linear-layer FLOPs per token / attention FLOPs per
  query-key pair for Qwen3-VL-2B). At the full run's sizes (325-token prompt, 127-token answers, 16 rows = 2 groups)
  packing costs about 0.40 of TRL's rows. With 8 answers at the 1,536-token cap it would cost 1.5x, so it falls back.
- Off with `--compile` (packed widths are not bucketed, so each new width would recompile).
- The share of tokens is logged per chunk as `shared/token_ratio` (packed tokens / the rows' own prompt + answer
  tokens). Expected about 0.37 at the full run's sizes.

## Checks (Mac, tiny Qwen3-VL with GQA and deepstack)

- Log-probs and entropies equal TRL's per-row forward (atol 1e-5), batch sizes 8 and 4, SDPA and eager, with the
  sampling mask; gradients of every parameter equal TRL's (atol 1e-4).
- Fresh-context review (Opus) re-derived it against the trl / transformers source and ran its own edge cases: two
  images per row, zero-length answers, a group split across a chunk boundary, all answers empty. Max differences
  5e-7 (log-probs, entropies) and 1.1e-6 (gradients). Its findings (long-answer cost, compile recompiles, host
  syncs, missing hook on a separate reference model, entropy-cut ordering) are fixed and tested.
- `train_rl --shared-prefix` runs 2 steps of TRL's real loop on the tiny model and saves the merged model;
  `shared/token_ratio` is logged between 0 and 1 on every step.

## R5 on the GPU (2026-11-08)

RTX 5090 (Vast offer 52598969, $0.593/h), e2, `experiments/phase2/configs/r5_shared.yaml`: 15 steps without the flag,
then 15 with `--shared-prefix`, same seed.

| | base | shared-prefix |
|---|---|---|
| 15 steps, wall clock | 7 min 45 s | **5 min 54 s** |
| time per step (train_runtime / 15) | 31.04 s | **23.65 s** (-24 %) |
| last-step rate (tqdm) | 30.43 s/it | 22.24 s/it |
| `shared/token_ratio` (steps 3-15) | - | 0.361 (predicted 0.37) |
| step 1 IoU / legal / answer length / grad norm | 0.5467 / 0.8984 / 117.89 / 0.0033 | 0.5467 / 0.8984 / 117.89 / 0.0033 |
| mean training IoU, steps 3-15 | 0.534 | 0.537 |

- Step 1 samples the same answers in both arms (same weights, same seed) and gives the same gradient norm. The logged
  loss differs (0.0109 vs 0.0115) because it is the mean of the micro-batch losses and the micro-batches hold
  different rows (whole groups); under dr_grpo their sum, and so the update, is the same. KL differs in the 4th
  decimal (bf16 kernels on different shapes). From step 2 the arms sample different answers, as expected.
- The base arm matches R4 (29.9 s/step), which ran on the same host.
- Verdict: the packed path runs on the real model and cuts a step by about a quarter. A 770-step run would take
  about 5 h 3 min instead of 6 h 38 min at these step times.

### Box notes (R5)

- First offer (57086600) never finished pulling the image (timeout, destroyed by launch, $0.24).
- On 52598969 the direct SSH port reset every connection; Vast's proxy (`ssh2.vast.ai`, the instance's `ssh_port`)
  worked, so the run record was pointed at it.
- run.sh caps the chain at `deadline - now - margin`; after a slow sync through the proxy that left 1,310 s, which
  would have killed the second arm at step 3. The chain's `timeout` was killed by hand and the deadline extended.
- The Mac sweep (launchd `com.objdel.vast.sweep`) was not loaded; a local timer did its job for this run.
