# Training guide: SFT, then RL

This guide follows the path that produced the current model: Qwen3-VL-2B with LoRA, SFT on LIDC-IDRI, then GRPO.
It tells you which command runs each stage, what each stage writes, and how to resume.
Renting a box is in the [GPU run runbook](../runbooks/gpu-run.md). The reward is in the [reward reference](../reference/reward.md).

## 1. The path

```text
  data/store --prepare_rl--> data/prepared --build_views--> data/views/lidc-bundle-v3
                                   |                                |
                                   |                          train_sft   (stable run + decay branches)
                                   |                                |
                                   |                          select_ckpt (D12 rule, entry gates)
                                   |                                |
                                   |                          g_probe     (20 RL steps: does RL help?)
                                   |                                |
                                   +--------------------------> train_rl   (770 steps, LoRA merged at the end)
                                                                    |
                                                                final_eval (whole val and test splits)
```

| Stage | Module | Input | Output |
|---|---|---|---|
| Prepare | `experiments.phase1.prepare_rl` | `data/store` | `data/prepared` |
| Views | `experiments.phase2.build_views` | store, prepared | `data/views/<view>` |
| SFT | `experiments.phase2.train_sft` | a training config | `runs/<run>/{stable,branches}/e<epochs>`, `manifest.json` |
| Select | `experiments.phase2.select_ckpt` | the run, exam views | `selected.json`, `report.json`, `gates.json` |
| G-probe | `experiments.phase2.g_probe` | a checkpoint, the probe view | `g_probe.json` |
| RL | `experiments.phase2.train_rl` | an SFT model folder, `data/prepared` | `trainer/checkpoint-N`, `final/`, `rl.json` |
| Evaluate | `experiments.phase2.final_eval` | a model folder, an exam view | `final_eval.json`, `items.jsonl` |

All commands run as `.venv/bin/python -I -m MODULE` on the Mac, and as `/opt/venv/bin/python -I -m MODULE` on a box.
In practice the GPU stages run through a chain config (see [Configs](../reference/configs.md)).

## 2. Before you train

1. `make test` and `make lint` pass.
2. The prepared set and the training view exist and match their manifests.
3. `make preflight CONFIG=experiments/phase2/configs/r1_qwen3vl2b.yaml` passes. `launch.sh` refuses to rent without a fresh stamp.
4. The G-vision base reference exists under `run_dir`. R0 writes it (`train_sft --vision_reference_only true`). `train_sft` refuses to start without it.

## 3. SFT (`train_sft`)

```bash
.venv/bin/python -I -m experiments.phase2.train_sft --config experiments/phase2/configs/r1_qwen3vl2b.yaml [--resume true]
```

`TrlParser` reads the YAML into `RunArgs`. Any field of `RunArgs` may also be given as a flag.

What the run does:

- **Model.** Qwen3-VL-2B in bf16 with a frozen base. LoRA rank 16, alpha 32, on all linear layers of the language model. The vision tower and projector are excluded (`configs.lora`).
- **Learning rate.** 5e-4 with LoRA (5e-5 for full fine-tuning).
- **Schedule.** WSD: linear warmup over 3 % of the stable run's steps, then constant. For each branch point b in `branches` (0.5, 1, 2, 3 epochs), a linear decay to 0 over 15 % of b.
- **Stable run.** One run that stops at the deepest branch point and saves a full checkpoint at every b. Each branch is a new trainer resumed from it and decaying. The result is bit-identical to an uninterrupted WSD run that decays at b.
- **Rows.** One row per voting rater with a target, in view-epoch order, length-grouped within each epoch (`train_sampling_strategy="group_by_length"`). "b epochs" means the first b passes over fresh views.
- **Batch.** Micro-batch 4, accumulation 8: effective batch 32 (`r1_qwen3vl2b.yaml`). Gradient checkpointing is on for Qwen3-VL-2B.
- **Assistant-only loss.** `assistant_only_loss=True`, `max_length=None`, no packing.
- **G-vision.** Each rolling checkpoint is checked. Ratio below 0.5 and TV below 0.002 together abort the run (exit code 3, reason in `manifest.json`). One alone only warns.
- **Memory.** `PYTORCH_ALLOC_CONF=expandable_segments:True` unless you set it.

Outputs under `runs/<output_dir>`:

| Path | Content |
|---|---|
| `manifest.json` | Config, view hash, code SHA-256, step plan, G-vision log, pid, status, timing, peak memory |
| `stable/e<b>/`, `branches/e<b>/` | Model-only folders with the processor. LoRA runs are saved merged |
| `rolling/`, `decay/` | Rolling checkpoints with optimizer. Deleted once the stable and branch folders exist |

Resume: `--resume true` continues a crashed run whose manifest status is `running` with the same config and view. It uses the newest complete rolling checkpoint, then each unsaved branch.

## 4. Select and probe

```bash
.venv/bin/python -I -m experiments.phase2.select_ckpt --run runs/r1-qwen3vl2b \
    --shifted data/eval_views/lidc-val-shifted --rescaled data/eval_views/lidc-val-rescaled \
    --rotated data/eval_views/lidc-val-rotated --prepared data/prepared --probes data/eval_views \
    --vision '{"prepared": "data/prepared", "reference": "runs/vision_ref/qwen3-vl-2b.json"}'
```

`select_ckpt` waits until `train_sft` has exited. It samples K = 8 answers per item with vLLM (`--sampler vllm`, the default) or with HF (`--sampler hf`).
The rule R-B picks the deepest branch that passes G-format, G-signal and G-vision, with a veto on signal (see `docs/research/` and `experiments/phase2/PLAN.md` D12).
It then judges that branch on the full views and writes `selected.json`.
`g_probe` runs 20 GRPO steps on the selected branch and compares held-out reward before and after:

```bash
.venv/bin/python -I -m experiments.phase2.g_probe --model qwen3-vl-2b --select runs/r1-qwen3vl2b/select \
    --train data/prepared --probe data/eval_views/lidc-val-shifted --out runs/r1-qwen3vl2b/g_probe
```

`--ckpt DIR` probes one model folder instead of a selection. `--select` refuses to probe when the gates failed and writes `{"verdict": null, "skipped": ...}`.

## 5. RL (`train_rl`)

```bash
.venv/bin/python -I -m experiments.phase2.train_rl --model qwen3-vl-2b --ckpt runs/r1/branches/e2 \
    --train data/prepared --out runs/rl --steps 770 --save-steps 50 --resume
```

| Flag | Default | Meaning |
|---|---|---|
| `--model` | required | Key of `configs.MODELS` |
| `--ckpt` | required | The SFT model folder to tune |
| `--train` | required | Prepared root with `lidc/train.jsonl` and caches |
| `--train-views` | none | View directories of train windows with their own expert tables (shifted, rescaled). Trained beside the centred frames |
| `--out` | required | Output directory |
| `--steps` | required | Optimizer steps in total, a resumed run included |
| `--save-steps` | 50 | Checkpoint every N steps |
| `--shared-prefix` | off | Pack each group's prompt once per chunk ([ADR-0007](../architecture/decisions/0007-shared-prefix-forward.md)). 24 % faster at the full run's sizes |
| `--lr` | probe rate | Learning rate. Default is 5e-5 with LoRA (`G_PROBE["lr_lora"]`) |
| `--resume` | off | Continue from the last `<out>/trainer/checkpoint-N` |
| `--procs` | usable CPUs minus 1 | Scoring pool workers |
| `--grpo` | `{}` | JSON overrides of `configs.grpo_config` |

What the run does:

1. Preload every expert table of the roots into memory and fork the scoring pool before CUDA starts. The pool is for exam metrics. Training rewards score on the GPU.
2. Build GRPO rows: LIDC train items, disagreement items repeated 4 times (8,208 rows). With `--train-views` the rows of all roots mix. A key with `::` marks a view row.
3. Build `FastGRPO` through `g_probe.make_trainer`: Qwen3-VL-2B with LoRA, vLLM colocate with sleep mode, `TRLReward` as the reward functions, G = 8, 32 prompts per step.
4. `trainer.train(resume_from_checkpoint=last)`.
5. At the end, merge the LoRA into the base and save `final/` with its processor. This is a plain model folder that `select_ckpt` and vLLM read unchanged.
6. Write `rl.json`: steps, learning rate, timing, TRL's log.

Checkpoints: `<out>/trainer/checkpoint-N` holds the LoRA adapter, optimizer, scheduler, RNG and trainer state. A lost box loses at most 50 steps (about 30 minutes).
To continue on a new box, sync the old `train_rl/` folder up and run the same chain again. `--resume` picks the last checkpoint.

### Settings that matter (from `grpo_config`)

| Setting | Value |
|---|---|
| Answers per step | 256 (32 prompts x 8) |
| Loss | `dr_grpo`, clip 0.2 and 0.28 (clip-higher), beta 1e-3, truncated answers masked |
| Learning rate | Constant. 5e-5 for LoRA |
| Sampling | Temperature 1, `top_p` 1, `top_k` 0, up to 1,536 tokens |
| vLLM | Colocate, 35 % of GPU memory, sleep mode |
| Gradient checkpointing | On. Turning it off ran out of memory beside vLLM on a 32 GB card |

### Measured behaviour of the full run

Source: `docs/research/notes/full_rl_2026-11-08.md`, `final_eval_2026-11-08.md`.

| Quantity | Value |
|---|---|
| Steps | 770 (3 epochs of 8,208 rows) in 5 h 52 min, 27.5 s per step, about $4.56 on an RTX 5090 |
| Training reward parts, first to last 10-step mean | IoU 0.52 to 0.77, F5 0.72 to 0.94, topology 0.66 to 0.93, legal 0.91 to 0.99 |
| Held-out set (128 shifted-val items, 8 answers each) | IoU_tc 0.511 (SFT e2), 0.634 (step 250), 0.659 (step 770) |
| Test split (2,200 items) | IoU_tc 0.645, F5 0.816, precision 0.695, recall 0.748, legal 99.4 % |

Most of the gain comes in the first 50 steps. The cause of the later plateau is that RL trained on centred frames while the exams use shifted frames.
Median object offset in the RL rows is 0.018 of the frame side. In the shifted validation view it is 0.259 (`docs/research/notes/rl_plateau_diagnosis_2026-11-08.md`).
The fix is built and not yet run: `train_rl --train-views` with `rl_shifted.yaml`.

## 6. Watching a run

| Where | What |
|---|---|
| `time/generate_s`, `time/rewards_s`, `time/sync_s`, `time/wake_s`, `time/sleep_s`, `time/batch_s` in TRL's log | Phase timers. `FastGRPO.timed` synchronises once per phase |
| `rewards/<term>/mean` | Each reward term. Geometric terms are gated by legality |
| `error/<code>`, `defect/<class>` | Share of answers per parser error or checker defect |
| `sampling/sampling_logp_difference/mean` | Token log-probability gap between vLLM and the trainer. Near 0.01 nats is healthy |
| `sampling/importance_sampling_ratio/{min,mean,max}` | The importance ratio |
| `shared/token_ratio` | With `--shared-prefix`: packed tokens over the rows' own tokens. About 0.36 at the full run's sizes |
| KL, median per step | 0.19 to 0.40 in the full run. Single-token spikes in the mean are a vLLM token with near-zero trainer probability. Judge drift by the median |
| `<out>/controller.jsonl` | Events of the `Controller` callback: best, switch, rollback, audit |
| `<out>/control.json` | Write it (atomically, with `os.replace`) to change reward weights or beta, or to force a rollback, during a run |

The `Controller` callback adds a hack audit to held-out evaluation, keeps the best checkpoint by `eval_rewards/iou_tc/mean` (legality-gated), and rolls back when an alarm fires.
`train_rl` as written does not attach it. Attach it in a script with `trainer.add_callback(Controller(trainer))` and an `eval_dataset` (see `experiments/phase1/RESULTS.md`, TRL notes).

## 7. Changing the training

| You want to | Change | Care |
|---|---|---|
| Another RL schedule | `--steps`, `--save-steps`, `--lr`, `--grpo` | `num_generations` must stay 8 |
| Train on shifted windows | `--train-views` | Build the windows with `build_eval_views --split train`. Ship them in the bundle with `--view` |
| Another model | `configs.MODELS`, `MAX_COMPLETION`, `_LORA`, config YAML | Check memory first. Full fine-tune of Qwen3-VL-2B needs about 34 GB |
| A different reward weight | `REWARD_V4` | A new reward. All expert tables become stale. Needs a new lab round |
| A guard on one term | `Controller.switch(reward_weights={"g_area": 0.5})` | Adds `w * (term - group mean)` to every advantage |
