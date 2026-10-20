# M7: TRL integration (`objdel.train`)

Plan: D7, D8, D9, D10, D12, §3.1–3.2, §4.6, §10; gates R2 (v4) and K12 (v4). Promotes `lab/t_stack` (`prod.py`
TRLReward, `r2_trl_equiv.py`, `k12_grpo.py`, `t3_grpo.py` controller, `t4_sft.py`) and the controller design of
`lab/t_dynamics`. TRL chat-row export moved here from M4 (EXECUTION.md §2).

## Result

| Gate / check | Result |
|---|---|
| **R2** TRL advantages = `advantages_v4` | max \|TRL − v4\| = **2.37e-7** (float32 rewards_per_func), 5 seeds / 5 group orders, real `GRPOTrainer` path (`_generate_and_score_completions`), 5 hand-set groups (illegal rollout, identical, minority topology, low-spread, answers no rater drew) |
| R2 dropped groups | identical and low-spread groups: advantage **exactly 0** |
| R2 group contiguity | single process: every reward call holds whole groups in order (checked inside `TRLReward`). World size 1/2/4 simulated with TRL's `RepeatSampler` and accelerate's `prepare_data_loader`: every process's first reward call holds whole groups, no prompt split across processes, for all three `grpo_config`s; a TRL-accepted config whose local batch is 4 (< G) is refused with `ValueError` |
| R1 on the prepared train/val sets | every target scored against its own item's caches: RIGA train 3,546 / val 384, LIDC train 8,157 / val 1,555: format 100 %, legal 100 %, topo_any 100 %, IoU_tc min 0.9900 (RIGA) / 1.0 (LIDC), F5 ≥ 0.95 100 % (min 0.988) |
| Live guard switches | `reward_weights` and `beta` mutated by the callback at step start take effect in that step's generation (advantage = v4 + w (g_area − group mean), ≤ 1e-5) and loss (`trainer.beta` read in `_compute_loss`: [1e-3, 0.5]); same via `control.json` |
| Rollback | restores the best evaluation's trainable parameters bit for bit, clears the optimizer state; held-out evaluation through TRL's `eval_strategy` reaches `on_evaluate` with `eval_rewards/<term>/mean` |
| **K12** (v4, Mac) | **PASS**: SmolVLM-256M, 40-step SFT, 4 GRPO steps with the production reward: every per-term log present, no NaN, 8 / 8 groups with non-zero advantages, duplicates 0 %, σ above the floor in 100 % of groups; 35 s/step, 340 s wall, peak footprint 6.2 GB |
| Tests | 40 in `tests/test_train_*.py` (18 s, 2 data-backed); `make lint` clean; `make test` 451 passed, 1 xfailed (M6's strict xfail), 9 deselected |

## What was built

| File | Lines | Contents |
|---|---|---|
| `src/objdel/train/dataset.py` | 197 | `PROMPTS`, `DISAGREEMENT_OVERSAMPLE`, `prepare`, `records`, `sft_rows`, `grpo_rows`, `cache_loader` |
| `src/objdel/train/rewards.py` | 113 | `LOGGED`, `AGGREGATOR`, `stop_ids`, `TRLReward` |
| `src/objdel/train/callbacks.py` | 111 | `Controller` (best checkpoint, rollback, alarm hook, live switches, control file, event log) |
| `src/objdel/train/configs.py` | 94 | `MODELS`, `MAX_COMPLETION`, `lora`, `processor_kwargs`, `sft_config`, `grpo_config` |
| `experiments/phase1/prepare_rl.py` | 48 | builds `data/prepared`, prints counts and the I5 token evidence |
| `experiments/phase1/smoke_grpo.py` | 144 | K12: short SFT, then GRPO with the production reward |
| `tests/test_train_dataset.py`, `test_train_rewards.py`, `test_train_trl.py`, `test_train_configs.py` | | below |

`objdel/train/__init__.py` re-exports nothing: `rewards` needs only numpy, `dataset` needs `datasets`, `configs` and
`callbacks` need torch/TRL/PEFT, so importing one does not pull in the others.

## Contract (EXECUTION.md §3, as delivered)

```python
from objdel.train import dataset, rewards, callbacks, configs
dataset.prepare(store_root, out_root, splits=("train", "val"), procs=3) -> manifest
dataset.records(root, dataset, split) -> list[dict]
dataset.sft_rows(root, dataset, split, tokenizer=None, max_tokens=None) -> datasets.Dataset   # messages, image, key, rater
dataset.grpo_rows(root, dataset, split, oversample=1) -> datasets.Dataset                     # prompt, image, key, W, H, names
dataset.cache_loader(root, maxsize=64) -> Callable[[str], dict[str, ExpertCache]]
rewards.stop_ids(model, processing_class) -> set[int]
rewards.TRLReward(caches, stop_ids, pool=None).funcs() -> (funcs, weights)   # LOGGED (w 0) + AGGREGATOR "advantage_v4" (w 1)
callbacks.Controller(trainer, metric="eval_rewards/iou_tc/mean", alarm=None)  # .switch(reward_weights=, beta=), .rollback(reason)
configs.sft_config(model, output_dir, **overrides) -> SFTConfig
configs.grpo_config(model, output_dir, **overrides) -> GRPOConfig
configs.MODELS, MAX_COMPLETION, lora(model), processor_kwargs(model, dataset)
```
*Superseded by F3 / F2 / S4 (EXECUTION.md §3 is current):* `Controller(trainer, metric=…, alarm=audit_alarm,
baseline=None, audit=True)` runs `eval.alarms_for` on held-out completions (B4-4) and refuses G ≠ 8 (B2-5);
`grpo_config` raises for G ≠ 8; `sft_rows` emits one row per *voting* rater (B2-4); `objdel.train.sample` added.

Changes from the §3 sketch: `TRLReward(provider)` became `TRLReward(caches, stop_ids, pool)`: the provider is
always "score the completions against the prompt's caches" (`objdel.reward.score_batch`), so only the cache lookup
and the stop ids vary; `sft_config` / `grpo_config` take the output directory and overrides.

Wiring (Phase 2):
```python
proc = AutoProcessor.from_pretrained(configs.MODELS[m], **configs.processor_kwargs(m, ds))
reward = rewards.TRLReward(dataset.cache_loader("data/prepared"), rewards.stop_ids(model, proc))
funcs, _ = reward.funcs()                       # weights are already in grpo_config
trainer = GRPOTrainer(model=model, reward_funcs=funcs, args=configs.grpo_config(m, out, eval_strategy="steps", eval_steps=50),
                      train_dataset=dataset.grpo_rows(root, ds, "train", dataset.DISAGREEMENT_OVERSAMPLE),
                      eval_dataset=dataset.grpo_rows(root, ds, "val").select(range(64)), processing_class=proc,
                      peft_config=configs.lora(m))
trainer.add_callback(callbacks.Controller(trainer, alarm=hack_audit_alarm))   # alarm: M8's plan-§9 lines
```

## Research (online, 2026-10-19) and what it decided

| Topic | Finding | Source |
|---|---|---|
| Reward function kwargs | prompts, completions, completion_ids, every dataset column except prompt, trainer_state, log_extra, log_metric (environments only with environment_factory) | https://huggingface.co/docs/trl/main/en/grpo_trainer |
| Truncation | no flag reaches reward functions. TRL cuts completion_ids after the first EOS (`_generate_single_turn`) and calls a completion truncated when `ids[-1] not in [*eos_token_ids, pad]`, with eos_token_ids = generation_config.eos_token_id + tokenizer.eos_token_id (grpo_trainer.py init). `rewards.stop_ids` reproduces the rule; SmolVLM: eos 49279 `<end_of_utterance>` + pad 2 | https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/grpo_trainer.py |
| Aggregation | `sum_then_normalize`: A = Σw r − group mean, divided by std only for scale_rewards ≠ "none"; `normalize_then_sum` = GDPO with a 1e-4 epsilon and a batch z (rejected in D7). reward_weights is a tensor read on every generation | same source; https://arxiv.org/abs/2601.05242 |
| Distribution over processes | `_calculate_rewards` runs each function on the local slice, then `gather`; advantages are computed globally and sliced back (`process_slice`). The train sampler is `RepeatSampler(mini_repeat_count=G, batch_size=generation_batch_size // G, repeat_count=num_iterations·steps_per_generation)`, the DataLoader batch is per_device·steps_per_generation and accelerate gives process p the p-th batch of every num_processes. TRL's own diagram shows a group split across GPUs (per_device 3, G 2): harmless for TRL's global advantages, fatal for an aggregator that computes per-group advantages locally. Hence local batch % G == 0 (configs) and a runtime check (TRLReward) | grpo_trainer.py `_get_train_sampler`, trl/trainer/utils.py `RepeatSampler`; https://huggingface.co/docs/accelerate/main/en/package_reference/torch_wrappers |
| Eval | `num_generations_eval` None → num_generations; GRPOConfig checks only the *global* eval batch % G, so a per-device eval batch below G would split groups; configs use 16. `GRPOTrainer.log` merges `eval_rewards/...` into the metrics dict passed to `on_evaluate` (verified) | grpo_config.py `__post_init__` |
| log_completions | writes `<output_dir>/completions/completions_<step>.parquet` with every reward column, `log_extra` columns and advantages; the hack audit (M8) can read it. On in grpo_config | grpo_trainer.py `log` |
| Live beta / weights | no TRL API (the PR for it was closed); `self.beta` is read when computing reference log-probs and the loss, `self.reward_weights` in `_generate_and_score_completions`: mutation at `on_step_begin` works (lab T3, re-verified here on 1.16.0). The KL reference model is re-loaded from the policy's `name_or_path` when beta > 0 without PEFT | grpo_trainer.py init, `_compute_loss` |
| SmolVLM + vLLM | issue #6294 open, fix PR #6718 open and unmerged (checked with `gh api`, 2026-10-19; latest release v1.16.0, 2026-10-15) → HF generation for SmolVLM | https://github.com/huggingface/trl/issues/6294 , https://github.com/huggingface/trl/pull/6718 |
| TRL VLM GRPO/SFT | GRPO tested on Qwen2/2.5-VL, SmolVLM2, Gemma3, LLaVA-NeXT; `image` column of PIL images with a conversational prompt; SFT needs `max_length=None` for VLMs and `assistant_only_loss` needs transformers ≥ 5.18 with vision data; TRL ships training chat templates (`idefics3_training.jinja`, `qwen3_vl_training.jinja`, `lfm2_2_5_vl.jinja`) | https://huggingface.co/docs/trl/main/en/grpo_trainer , https://huggingface.co/docs/trl/main/en/sft_trainer , `.venv/.../trl/chat_templates/` |
| transformers 5 | `apply_chat_template` returns a BatchEncoding; `dtype` replaces `torch_dtype`; SFTTrainer loads float32 when no dtype is given (unlike `from_pretrained`, which infers it in v5); SFTConfig/GRPOConfig default bf16=True (Mac runs override) | https://www.datacamp.com/tutorial/transformers-v5-tokenization , TRL SFT docs |
| SmolVLM processor (sandbox) | resizes every image so its longest edge equals `size["longest_edge"]`, up or down: 64 image tokens at 512, 320 at 1024 (5 tiles of 512), for 250, 500 and 512 px inputs alike. D12 kept: 1024 for RIGA, 512 for LIDC (`processor_kwargs`) | `scratchpad/impl/M7/proc.py`; https://huggingface.co/HuggingFaceTB/SmolVLM-256M-Instruct |
| Dr. GRPO / DAPO | constant length normaliser and no std scaling (dr_grpo, scale_rewards none); clip-higher 0.28; DAPO dynamic sampling ≈ v4's low-spread drop | https://arxiv.org/abs/2503.20783 , https://arxiv.org/abs/2503.14476 |

## Decisions

**Prepared set (`data/prepared`, built once).** Frames, targets and caches are computed once and hashed; rows only
reference them. Images are PNG files (lossless; `datasets.Image` decodes a path lazily, so a 3,546-row RIGA SFT set
does not hold 3,546 copies of six-times-repeated images); caches are `ExpertCache.save` files loaded through an LRU
(`cache_loader`, 64 items: a RIGA item's two caches hold ≈ 9 MB derived state, all 591 RIGA items would need
≈ 5 GB) whose digests are checked against the records. Items framed in batches of 12 per `pool.map`: main-process
peak 1.8 GB at 6, 2.4 GB at 32 on RIGA, of which 1.7 GB is `store.load`'s 16-row Parquet batches (M4); identical
manifests at both sizes (the prepared set below was built at 32).

**SFT rows: every rater is a row** (not a random rater per epoch; since B2-4 every *voting* rater, 10 / 23,469
final LIDC train targets dropped). An epoch then shows each image's full rater
distribution, which is what the RL topology term matches (π_SFT ≈ f), and it needs no sampler; items weigh by
their rater count (RIGA 6, LIDC 2–4), which only matters when datasets are mixed (a Phase 2 choice).

**I5: drop over-budget answers from SFT, keep the budget.** Measured on the prepared train set (answer tokens + 1
end-of-turn token):

| model (budget) | RIGA train over budget | max | LIDC train over budget | max |
|---|---|---|---|---|
| SmolVLM-256M (1,536) | 10 / 3,546 | 2,269 | 0 / 8,157 | 1,187 |
| LFM2.5-VL-450M (1,024) | 8 / 3,546 | 1,428 | 0 | 751 |
| Qwen3-VL-2B (1,536) | 10 / 3,546 | 2,262 | 0 | 1,182 |

*Superseded:* on the final prepared set (BinRushed1 has no model image since B3-1) 1 / 3,282 RIGA train targets is
over budget for each of the three models, 0 / 23,469 LIDC train (re-counted with each tokenizer, 2026-10-19).

The ten are 9 BinRushed1 + 1 BinRushed3 targets on 6 images (BinRushed1 image018, 030 ×3, 042, 045 ×3, 048; BinRushed3
image014), the jagged reconstructed pen strokes whose cup needs τ ≤ 1 (M5). RIGA p50 / p99 = 310 / 474 Qwen tokens.
Raising max_completion_length to 2,304 for 0.28 % of targets would make every runaway rollout 1.5× longer (HF
generation runs the batch until the longest sequence stops; at the start of RL many rollouts are truncated) and
would teach 2,000-token answers that mostly encode reconstruction artefacts. The items stay in RL: the reward needs
the caches, not the target. `sft_rows(..., tokenizer, MAX_COMPLETION[model])` applies it per model.

**I7: no augmentation in Phase 1, flips only through the masks.** Sandbox (`scratchpad/impl/M7/i7_flip.py`):
mirroring a target's vertices (x → 1000 − x, re-canonicalized) and rasterizing does not give the mirrored mask on
**7,078 / 7,092 RIGA objects (99.8 %, 121,810 px)** and on 0 / 8,157 LIDC objects: at S = 500 one pixel is 2 grid
units, so slanted edges pass through pixel centres where the top-left tie rule decides, and the rule is not
mirror-symmetric. A flip must flip image and rater masks and rebuild targets (`build_target`) and caches
(`ExpertCache.build`). Not built: Phase 1 runs no training; the note in `dataset.py` records the rule.

**Oversampling disagreement items ×4 (`DISAGREEMENT_OVERSAMPLE`).** TRL's `RepeatSampler` takes no weights, so a
disagreement item is a repeated row (no sampler subclass). Disagreement = the raters' joint size-filtered
signatures differ (from the caches). Counts: LIDC train 276 / 2,557 items with pixels (10.8 %), val 43 / 483; RIGA
0 (disc and cup always one component without holes). On the full LIDC set (after the CT fetch, issue I9) the rate
is ≈ 4.0 % (308 / 7,623 under the size rule, M4). With 32 prompts per step, uniform sampling gives 1.3 disagreement
groups per step at 4 % and no topology-matching signal at all in 27 % of steps (0.96^32); ×4 gives 4p / (1 + 3p) =
14.6 % of prompts, 4.7 groups per step, a step without one in 0.6 %, and unanimous items stay ≥ 85 % of prompts.
On the current pixel subset (10.8 %) ×4 gives 33 % of prompts; re-check the factor after the store rebuild.
*Final (I23, F4, B2-4):* LIDC train 289 / 7,341 disagreement items (3.9 %), val 44 / 1,084: ×4 gives 14.1 % of
prompts, 4.5 groups per step, a step without one in 0.8 %; the factor stands.

**Configs.** SmolVLM-256M and LFM2.5-VL-450M full fine-tuning, Qwen3-VL-2B LoRA r16 all-linear (24 GB card);
SFT lr 5e-5 full / 5e-4 LoRA, GRPO 1e-6 / 1e-5; one step = 256 completions (8 per device × 32 accumulation),
eval batch 16; vLLM colocate for LFM and Qwen, HF generation for SmolVLM; log_completions on. These are Phase 2
starting points; only D8's settings are gated (tests).

**Controller.** Held-out evaluation is TRL's own (`eval_strategy="steps"`); the callback only tracks the best
`eval_rewards/iou_tc/mean` (configurable), saves trainable parameters per rank to `<output_dir>/best/`, rolls back
when `alarm(metrics)` names a reason (M8's plan-§9 lines plug in here), and applies guard switches at the next
step, from code or from `<output_dir>/control.json` (D9: "guards switched by hand"). The typed alarm itself stays
deferred to Phase 5 (D9).

**Lab code not promoted:** the provider indirection, `history`, the v2 advantages, the lab parser and geometry
(M2/M6 own them), `normalize_then_sum` arms, `riga_samples` crops (M5 frames replace them).

## Tests (written first and run red: `objdel.train` did not exist)

- `test_train_rewards.py` (no model): names and weights; values equal `score_batch` + `advantages_v4` +
  `topology_reward` on a mixed batch (good, shifted, bow-tie, junk, truncated); one computation per batch;
  truncation from completion_ids (EOS, pad, cut); partial and interleaved groups refused; typed `error/` and
  `defect/` metrics and the `log_extra` column; `stop_ids` = TRL's rule.
- `test_train_trl.py` (tiny offline Llama + WordLevel tokenizer, CPU, ≈ 4 s): R2; world size 1/2/4 contiguity for
  the three configs; the split-group negative control; guard switches via API and control file; beta raise without
  a reference model refused; rollback; evaluation feeding the controller.
- `test_train_dataset.py`: synthetic store (RIGA-like item, a LIDC disagreement item, one without pixels, one
  single-rater item): counts, skips, layout, manifest hashes, determinism with 1 and 2 workers, cache digests and a
  swapped cache detected, SFT schema / placeholders / one row per rater / token budget, GRPO schema / RL-only /
  oversampling, row hashes stable, targets score as perfect answers; `@pytest.mark.data` R1 on the prepared
  RIGA and LIDC train sets.
- `test_train_configs.py`: D8 fields for each model, budgets, vLLM choice, whole-group batch sizes, SFT settings,
  LoRA choice, D12 processor kwargs.

## Prepared set (`.venv/bin/python -I -m experiments.phase1.prepare_rl`, 61 s with 3 workers)

Store manifest SHA-256 `20f315fa…d4475d`; reward fingerprint `0a6a4c09…8b5f5d`. 683 MB (images 655 MB, caches 22 MB).

| split | items | skipped (no pixels) | targets | RL items | disagreement | target-set SHA-256 | jsonl SHA-256 |
|---|---|---|---|---|---|---|---|
| riga/train | 591 | 0 | 3,546 | 591 | 0 | `370abe7c23e8fd5eae4d58f239534a2eb830cf8c12f2d8068cdcf819144b3f59` | `993806772ebee05d2e18c320e59aefcbece6e108160f103775e9f8979c16e060` |
| riga/val | 64 | 0 | 384 | 64 | 0 | `cf1a482aa57837d3cc146f27c83ac7d26cebe51cbfaaf5c4f72e1dec698af6b8` | `6499f26925004cc5d8191988e0c69924467d9e8692952943a1947781cd04240c` |
| lidc/train | 2,557 | 5,066 | 8,157 | 2,557 | 276 | `d87d119b94935616cf6989f72305357c1dd3d5056eb500a91a1b7e17496fea02` | `d5565717eaf42beb9a5b3cd9742909e85c9ac75e8fc033886f2fa6c0b1099df4` |
| lidc/val | 483 | 602 | 1,555 | 483 | 43 | `b72925d19fcfa006d6e0277b8e0d38d866d705c093f0c8b9c52bebcb797307ca` | `a8dd94158682ab1b69787ef347d7948ed9e7ac8133e4d519b4ad13fcd79965ec` |

**Superseded by the final rebuild** (`data/prepared/manifest.json`; store `37fe7665…`, reward fingerprint
`b1eb9f51…`): riga/train 547 items, 45 skipped (BinRushed1, B3-1), 3,282 targets, 0 disagreement; riga/val 58, 5,
348, 0; lidc/train 7,341, 0, 23,469, 289; lidc/val 1,084, 0, 3,444, 44; all items RL-eligible.

No target errors. The RIGA target-set hash differs from M5's train+val hash because M5 hashed train and val together
(and LIDC items without pixels); per (item, rater) the texts are built by the same `build_target` call.

## K12 smoke

`HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase1.smoke_grpo --grpo-steps 4` (Apple M3, MPS, fp32, load
average 5–10 from other agents). SFT: SmolVLM-256M, 48 LIDC train items (163 rater targets), 40 steps of 4 rows
(micro-batch 1 × 4), LoRA r16 lr 5e-4, merged. GRPO: 8 other LIDC train items (all disagreement items, ×4 rows),
`grpo_config` with 2 prompts × 8 completions per step (micro-batch 1 × 16), LoRA r16 lr 1e-5, β 1e-3, dr_grpo,
max_completion_length 1,536, HF generation, production `TRLReward` on the prepared caches.

| | value |
|---|---|
| SFT | loss 0.952 → 0.517 (≈ 1 epoch over the rows); 4.8 s/step; 192 s |
| GRPO | 4 steps, 35.4 s/step (generation + scoring 13–17 s, update 28–33 s), 142 s |
| Wall | 340 s total |
| Memory | peak footprint 6.21 GB (`/usr/bin/time -l`: RSS 2.17 GB + MPS driver ≤ 3.87 GB, cap `--mps-fraction 0.35` = 4.45 GB) |
| Per-term logs | `rewards/{format_ok, legal, iou_tc, f5, g_area, topo_any, r_topo, advantage_v4}/mean` at every step, plus `error/*`, `defect/*` |
| NaN | none in any log |
| Advantages | 8 / 8 groups non-zero (no group dropped, `frac_reward_zero_std` 0); max \|A\| 1.65–3.37 |
| Duplicate rate within groups | 0 % at every step (16 / 16 distinct completions) |
| σ above the floor | 100 % of groups for IoU_tc (σ > 0.05) and for any z-scored term |
| Terms (per step) | format_ok 0.94–1.0; legal 0.31–0.50; IoU_tc 0.11–0.15; F5 0.57–0.72; g_area 0.23–0.64; topo_any 0.75–0.94; r_topo 0.053 at step 1, then 0 |
| Errors | `whitespace` 1 / 16 at step 3; illegal shapes: self-crossings, pinches, spikes (`defect/*`) |
| Completion length | mean 85–95, max 256 tokens; clipped ratio 0 (no truncation) |
| KL | 0 → 2e-4 |

Reading: a 40-step SFT is far from usable geometry (IoU_tc ≈ 0.13, half the shapes illegal) but it is
non-collapsed: no duplicates and 100 % of groups above the σ floors, the Phase 2 GRPO-signal gate's direction
(lab K12 (a): a 50-step overfit gave 100 % dropped groups). r_topo is non-zero only when a legal rollout reaches
IoU_tc ≥ 0.5 with a rater topology, which happened once.

Memory work (all measured with `/usr/bin/time -l`): micro-batch 4 SFT took the MPS driver to 8.07 GB (footprint
10.1 GB with held-out evaluation); micro-batch 1 SFT 2.76 GB. The SFT trainer's tensors survive `del` (reference
cycles): `gc.collect()` before `torch.mps.empty_cache()`. GRPO micro-batch 2 still reached 4.8 GB of live MPS
tensors in `training_step`; micro-batch 1 fits under the 4.45 GB cap. TRL's held-out evaluation computes
log-probabilities for `per_device_eval_batch_size` ≥ G = 8 completions at once and OOMs under any cap that keeps
the footprint near 6 GB, so the final smoke has no evaluation. The evaluation path ran in the first, uncapped smoke
(same SFT, 3 GRPO steps, evaluation of one validation prompt every 2 steps, footprint 10.1 GB):
`eval_rewards/*/mean`, `eval_error/*` and `eval_defect/*` reached `Controller.on_evaluate`, which recorded
`best` = 0.286 (`eval_rewards/iou_tc/mean`) at step 2; 3 / 3 train generations had non-zero advantages in both
groups, duplicates 0 %, σ above the floor 100 %, 43–53 s/step.

## Open issues

- **The prepared LIDC set is the 214-series pixel subset** (issue I9). After the CT fetch and the store rebuild,
  re-run `prepare_rl` (no code change); counts, hashes and the oversampling rate change. *(Resolved: I23; final
  counts above.)*
- `store.load` holds 16 RIGA rows per Parquet batch (≈ 1.7 GB in the reading process); a smaller `BATCH_ROWS` for
  RIGA would cut `prepare`'s footprint (M4's file).
- K12 runs SFT and GRPO with LoRA and micro-batch 1 (memory), not the configured full fine-tune of SmolVLM; the full
  fine-tune path is exercised only by the tiny-model tests. Its footprint is 6.2 GB, 0.2 GB over the 6 GB target.
- On the Mac, TRL's held-out evaluation (log-probabilities of a whole group at once) does not fit the memory budget;
  on the Phase 2 GPU it is a non-issue, but per_device_eval_batch_size cannot go below G = 8 per process.
- `vllm` is not installed in the Mac venv: the vLLM settings for LFM / Qwen are configured but not run (Phase 2).
- Item weighting in mixed RIGA + LIDC SFT runs (rows per rater) and frame jitter (M5) are Phase 2 choices.
