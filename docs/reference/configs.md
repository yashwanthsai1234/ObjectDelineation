# Experiment configs reference

The YAML files in `experiments/phase2/configs/` are of two kinds. A **training config** holds the settings of an SFT run (`train_sft`).
A **chain config** lists steps that `run.sh` executes on a box, one after the other.
Both kinds are read by code in `experiments/phase2/`. There is no JSON schema file. The schema is the code named below.

## 1. Training config (`RunArgs` in `train_sft.py`)

Used by `python -I -m experiments.phase2.train_sft --config FILE`. Also read by `preflight.load`, `policy.chain` and `policy.vision_refs`.

| Key | Type, default | Meaning |
|---|---|---|
| `model` | required | Key of `objdel.train.configs.MODELS`: `smolvlm-256m`, `lfm2.5-vl-450m`, `qwen3-vl-2b` |
| `view` | required | Training view directory, relative to `data_root`, for example `views/lidc-bundle-v3` |
| `output_dir` | required | Run directory, relative to `run_dir`. Also the Hugging Face revision name on push |
| `data_root` | `$OBJDEL_DATA_ROOT`, else `/workspace/data/lidc` on a box, else `data` | |
| `run_dir` | `$OBJDEL_RUN_DIR`, else `runs` | On a box `/workspace/runs/<run id>`, set by `run.sh` |
| `model_path` | `None` | Local weights. Default is `MODELS[model]` |
| `branches` | `[0.5, 1.0, 2.0]` | Epochs at which WSD decay branches are cut. `[]` means stable only, one epoch, nothing saved (timing) |
| `a0` | `false` | Empty every `holes` list in the targets (the A0 control of experiment C2) |
| `freeze_vision` | `false` | Freeze the vision tower and projector (repair of a G-vision abort) |
| `seed` | 1 | |
| `rolling_minutes` | 10.0 | Rolling checkpoints with a G-vision check at this interval |
| `vision` | `{}` | `{prepared, reference}`. `reference` is a number or a JSON path under `run_dir`. `off` for discarded timing runs. `{}` is refused |
| `vision_reference_only` | `false` | R0: write the G-vision base reference and exit |
| `resume` | `false` | Continue a crashed run in `output_dir` |
| `steps_per_epoch` | `None` | Override for smoke tests and timing |
| `sft` | `{}` | Overrides of `configs.sft_config` (batch size, accumulation, workers) |
| `deterministic` | `false` | `torch.use_deterministic_algorithms(True, warn_only=True)` for timing arms |
| `branches_from` | `None` | A `g_probe.json` under `run_dir`. Branches become the depth of the checkpoint with the final G-probe verdict |
| `profile` | `None` | `torch.profiler` output directory (timing runs only) |
| `stages` | `[]` | Runner only: subset of `train`, `select`, `topology`, `probe` (`policy.chain`) |
| `sampler` | `None` | Runner only: `hf` or `vllm`, passed to `select_ckpt` |

A wrong `vision` value fails in `__post_init__`, not at the first rolling checkpoint.

When `stages` is absent and the config has `branches` or `branches_from`, `policy.chain` uses `train`, `select`, `probe`.
A config without branches only trains. A chain that scores first builds the CUDA checker and runs `tests/test_gpu_scoring.py`, so a failed build stops the chain before training.

## 2. Chain config (`steps:`)

```yaml
steps:
  - ["~", bash, -c, "nvidia-smi -q > {rd}/box.txt 2>&1"]       # "~": a failure is logged, the chain goes on
  - [bash, experiments/phase2/docker/env_check.sh, "{rd}/ckpt", "{rd}/env_check.json"]
  - [-m, experiments.phase2.train_rl, --model, qwen3-vl-2b, --ckpt, "{rd}/ckpt", ...]
```

- Each step is a list of arguments. `-m module ...` runs as `/opt/venv/bin/python -I -m module ...`. `bash script ...` runs as is.
- `{rd}` is `/workspace/runs/<run id>`. `{data}` is `/workspace/data/lidc`. Literal braces in shell text are doubled: `{{`.
- The first list item `"~"` makes a step soft.
- The steps are joined with `&&`. A hard failure stops the chain. `run.sh` wraps the chain in `timeout` and writes `EXIT` and a heartbeat file.
- `run.sh RUN_ID [--cap-hours H] CONFIG` takes flags before the config. A trailing flag is refused.

## 3. Configs in this repository

| File | Kind | Purpose | Status |
|---|---|---|---|
| `r0.yaml` | chain | R0 box checks for Qwen3-VL-2B: image and card check, cross-architecture hashes, K3/K4 idle timing, sampler parity, HF sampling rate, G-vision base reference, a 30-step timing arm, the checker on the GPU | Record of the R0 chain. Its run date is not stated in the header |
| `r0_qwen3vl_arm.yaml` | training | R0 timing arm: 30 steps with the R1 settings, discarded. Gives s per step and peak memory | Record |
| `r1_lfm.yaml` | training | R1 rung 1: LFM2.5-VL-450M, full fine-tune, branches 0.5, 1, 2 epochs. Trained 2026-11-02 and failed the entry gates | Record. Not to be rerun |
| `r1_qwen3vl2b.yaml` | training | R1 rung 2: Qwen3-VL-2B, LM-only LoRA rank 16, view `lidc-bundle-v3`, branches 0.5, 1, 2, 3, effective batch 32 (4 x 8), sampler `vllm` | Current SFT recipe. Also the config that preflight, launch and the fetch scripts read for every RL run |
| `c2_a0_seed1.yaml` | training | C2: A0 control (holes emptied), seed 1, to the depth of the G-probe's checkpoint (`branches_from`), stages `train`, `topology` | Defined |
| `c2_a0p_seed2.yaml`, `c2_a0p_seed3.yaml` | training | C2: A0' (holes kept), seeds 2 and 3, same depth rule | Defined |
| `rx8.yaml` | training | R-x8: one A0' seed, non-simple rows repeated 8 times (`views/lidc-bundle-v2-x8`) | Defined. Needs owner approval |
| `checker_gpu.yaml` | chain | K1: the CUDA checker alone (no model, no data): build, memcheck, initcheck, then the equivalence and timing step | Run 2026-11-06, RTX 3090 |
| `p1_rl.yaml` | chain | P1: RL readiness on an RTX 5090. Env check, checker, sanitizers, profile arms (`experiments.phase2.profile`) | Run 2026-11-06 |
| `r2_qwen3vl2b.yaml` | chain | R2: the pre-registered 20-step G-probe on e2 | Run. PASS: reward 0.694 / 0.679 to 0.769 |
| `r3_gpu_test.yaml` | chain | R3: scorer on the GPU, OpenOnce and ImageMemo, 8 GRPO steps | Run 2026-11-07, billed $0.440 |
| `rl_full.yaml` | chain | The full RL run: 770 steps, checkpoint every 50, always `--resume` | Run 2026-11-08 |
| `r4_speed.yaml` | chain | R4: the full run's settings for 15 steps with vLLM asleep, then awake | Run 2026-11-08 |
| `r5_shared.yaml` | chain | R5: 15 steps without and with `--shared-prefix` | Run 2026-11-08 |
| `rl_shifted.yaml` | chain | RL on centred frames plus shifted and rescaled train windows, 300 steps, with `--shared-prefix`, then scoring of the whole val and test splits | Planned. Needs the owner's OK and budget |

The statuses come from the config headers and the notes in `docs/research/notes/`.

### Reading the RL chains

`rl_full.yaml` has five steps:

1. A hardware snapshot (soft).
2. `env_check.sh`: capability, bf16 peak, Triton, `torch.compile`, a vLLM engine on the checkpoint with CUDA graphs.
3. `run_box.sh`: build the checker for the card, GPU tests, equivalence suite, timing sweep. Then `tests/test_gpu_scoring.py`.
4. Sanitizers (soft step, but a missing binary or an error fails it): memcheck then initcheck.
5. `train_rl --steps 770 --save-steps 50 --resume`.

`rl_shifted.yaml` adds `--train-views {data}/views/lidc-train-shifted {data}/views/lidc-train-rescaled` and `--shared-prefix`. After training it runs `final_eval` on the validation and test views.

### Where model and sampling settings live

Config YAML does not hold GRPO settings. They are in code: `objdel.train.configs.grpo_config` and `sampling`.
Override them for one run with `train_rl --grpo '{"key": value}'` (JSON). The key facts:

| Setting | Value | Source |
|---|---|---|
| Group size `num_generations` | 8. Any other value raises | `grpo_config`, `Controller` |
| Batch | `per_device_train_batch_size` 16, `gradient_accumulation_steps` 16: 256 answers per step | `grpo_config` |
| Loss | `dr_grpo`, clip 0.2 and 0.28, beta 1e-3, truncated answers masked | `grpo_config` |
| Learning rate | 1e-5 for LoRA in `grpo_config`, but `make_trainer` sets 5e-5 for LoRA unless `--lr` is given (`G_PROBE["lr_lora"]`). Full fine-tune: 1e-6 and 5e-6 | `configs.py`, `prereg.py`, `g_probe.make_trainer` |
| Sampling | Temperature 1, no filter, `max_tokens` 1,536 (Qwen3-VL) | `configs.sampling` |
| vLLM | colocate, memory share 0.35, sleep mode, model length 4,096 | `grpo_config` |
| Precision | bf16 when `is_torch_bf16_gpu_available()`, else fp32 | `configs._bf16` |
