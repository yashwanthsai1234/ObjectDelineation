# B5: environments for Phase 2 and Phase 3 (H-A5-7, H-A5-1, H-A5-2, Vast images)

The pre-registered criteria are in `experiments/phase2/lab/B5/prereg.py`. They were fixed at 18:12 PDT on 2026-10-19, before the vLLM parity runs and the container tests. That file holds the text verbatim, because the harness refused a `.md` write.

## Pre-registered criteria (summary of prereg.py)

**Q1 / H-A5-7**
- **(a) one venv** passes only if all four hold:
  - a lock resolves with TRL 1.16, vLLM ≥ 0.30 and transformers ≥ 5.18, keeping numpy 2.5.4 and numba 0.68.0, with every override named;
  - the linux/amd64 image imports vLLM, `vllm.LLM`, `trl.trainer.grpo_trainer`, `trl.generation.vllm_generation` and `objdel.train.*`;
  - `pytest -m "not slow"` has 0 failures, with the same counts as the pinned venv;
  - `xarch_hashes --strict` is identical in every section.
- **(b) two environments:**
  - passes for Phase 2 if the two-venv image builds and imports;
  - passes for Phase 3 only if (i) the trainer needs no vLLM import, (ii) weight sync works on one GPU, and (iii) the TRL client matches the server's vLLM version. This is checked by reading the source.
- **(c) SGLang** is considered only if it resolves without overrides, and then only as a sampler.

**Q2 / H-A5-1:** 32 validation items per dataset × 3 models.
- vLLM's `prompt_token_ids` must equal the training env's GRPO `input_ids` in 100 % of cases.
- The SFT collator's token sequence must start with the generation prompt in 100 % of cases.

**Q2 / H-A5-2(a), CPU proxy:** 8 items × n = 8 samples at T = 1, top-p 1, top-k off, 96 tokens.
- fp32: median per-token |Δ| ≤ 1e-3 nats, and p99 per-sequence |ΣΔ| ≤ 0.5.
- bf16: median ≤ 0.02, and p99 per-sequence |ΣΔ| ≤ 0.5.
- Leak check: no generation_config default (temperature, top-k, top-p) may leak into sampling.

**Q3:** each image builds for linux/amd64 on a builder capped at 2 CPUs and 4 GB. Each venv must report CUDA 13.0, an arch list containing sm_80 and sm_120, and import without a GPU. Sizes are reported and everything is deleted afterwards.

## Results

### Q1: what resolves (`resolve_matrix.py`; `uv pip compile`, x86_64-manylinux_2_40, Python 3.12)

| Candidate | Resolves | Key versions or the reason it fails |
|---|---|---|
| Training stack as pinned | yes | torch 2.14.1, transformers 5.20.1, TRL 1.16.0, numpy 2.5.4, numba 0.68.0 |
| `vllm==0.32.0` alone (sampler) | yes | torch 2.13.0, transformers 5.19.0, numpy 2.3.5, numba 0.65.0, xgrammar 0.2.8, llguidance 1.7.6, flashinfer 0.7.0.post1 |
| vLLM 0.32 + our core pins | **no** | vLLM pins `numba==0.65.0`; we pin 0.68.0 |
| Re-pin (torch 2.13 / transformers 5.19 / TRL 1.16 / vLLM 0.32) | **no** | numba conflict |
| … plus an override of numba only | **no** | vLLM requires `mistral-common[image]`, which caps numpy below 2.4 on Python < 3.13 |
| … plus overrides of numba 0.68.0 **and** numpy 2.5.4 | **yes** | from the Mac's freeze only torch, torchvision, transformers and setuptools change (setuptools < 81 for vLLM), plus vLLM's added dependencies |
| Re-pin with numba 0.65 / numpy 2.3.5 (no overrides) | yes | hash cost below |
| `trl[vllm]==1.16.0` (TRL's declared range ≤ 0.31) plus both overrides | yes | picks **vLLM 0.30.0** with transformers 5.20.2: 0.30 has no transformers ceiling, while 0.31 caps < 5.18. Not functionally tested |
| vLLM 0.32 + torch 2.14.1 | **no** | `torch==2.13.0` pin |
| vLLM 0.32 + transformers 5.20.1 | **no** | `transformers<5.20` (and numba) |
| SGLang 0.5.22 + training pins | **no** without overrides | pins `numba==0.65.1`, `transformers==5.19.0`, `tokenizers==0.23.2`, manylinux_2_34 only. Resolves with both overrides, but has no TRL integration and no SmolVLM. Rejected |

**Cost of the clean re-pin (numpy 2.3.5, numba 0.65.0)**, from `xarch_hashes --strict` against the Mac's hashes:
- 24 of 26 sections are identical, including `targets` and `target_set_hash`;
- only `fingerprint` and `caches` differ, because the fingerprint contains the numpy version string;
- so that path means re-running prepare_rl to rebuild the caches. The overrides avoid this.

**What I checked in TRL 1.16 source**
- `is_vllm_available` only **warns** outside vLLM 0.22–0.31. The warning was confirmed in the image.
- `_init_vllm` raises `ImportError` when vLLM is missing, in **both** modes.
- `vllm_client.py` imports vLLM in the trainer: `get_ip`, `nccl_common.trainer_init` and `packed_nccl_broadcast_producer`.
- The client's `init_communicator` is NCCL-only. It raises `NotImplementedError` for a non-CUDA device and has no IPC path; I searched both the v1.16.0 tag and main.
- TRL's docs: "the vLLM server and the trainer must run on separate CUDA devices". NCCL's docs: "Using the same CUDA device multiple times as different ranks of the same NCCL communicator is not supported and may lead to hangs."
- TRL PR #7799 (open, 2026-10-19): vLLM 0.32 changed the logprob format of `/inference/v1/generate`, the endpoint multimodal prompts use. TRL 1.16 server mode therefore gets a `KeyError` on the first batch. Colocate mode is unaffected.
- So (b) fails Phase 3 criteria (i), (ii) and (iii).

**Option (a) image (`Dockerfile.phase3`, lock p3.txt plus p3.override)**

| Check | Result |
|---|---|
| Build | OK in 1,455 s |
| Image size | 4.34 GB; `/opt/venv` 8.6 GB unpacked |
| torch | 2.13.0+cu130, CUDA 13.0, arch sm_75 / 80 / 86 / 90 / 100 / 120 |
| Imports | vllm 0.32.0, xgrammar, llguidance 1.7.6, `vllm.LLM` / `SamplingParams` / `StructuredOutputsParams`, `trl.trainer.grpo_trainer`, `trl.trainer.sft_trainer`, `trl.generation.vllm_generation`, `trl.generation.vllm_client` (stateful-engine path True) and all `objdel.train.*` modules: all OK. TRL printed its vLLM-range `UserWarning` |
| numpy override | mistral-common 1.12.0 and `vllm.tokenizers.mistral` import fine with numpy 2.5.4 |
| Fast suite | same command as `x86_inside.sh`: 666 passed, 1 failed, 53 deselected. The failure is `test_geometry.py::test_length_cap`, a wall-clock limit of 1 ms: it measured 1.026 ms under emulation at load average about 29, and passed 3/3 when rerun. It touches no torch, transformers, numba or numpy path |
| Hashes | `xarch_hashes --compare build/x86/io/hashes-x86.json hashes-p3.json --strict`: **26 / 26 identical, PASS**, against the pinned Phase 1 x86 run |

**Online facts**
- vLLM main bumped to torch 2.14.0, torchvision 0.29.0 and triton 3.8.0 on 2026-10-14 (PR #52183). Main also allows transformers `>=5.16.1,<5.21.0` and still pins `numba == 0.65.0`.
- vLLM's transformers caps are precautionary: PR #59614 was prompted by breakages in other models (MiniMax-M3, Qwen3.8-Flash-Next), not ours.
- transformers 5.19 and 5.20 have the same `TrainingArguments` fields: WSD schedule, `train_sampling_strategy`, `batch_rebalance`. 5.19 also satisfies the ≥ 5.18 requirement for assistant-only loss.

### Q2: parity (vLLM 0.32.0 macOS CPU wheel, transformers 5.19.0 / torch 2.13.0, against the training env's transformers 5.20.1 / torch 2.14.1)

`export_prompts.py` builds the reference prompts in the training env, following TRL 1.16's GRPO `_tokenize_prompts` and its SFT collator.

| Model | LIDC prompt ids equal | RIGA prompt ids equal | Prompt lengths | SFT sequence starts with the generation prompt |
|---|---|---|---|---|
| LFM2.5-VL-450M | 32 / 32 | 32 / 32 | 332; 150 / 342 | 32 / 32, loss span equal to completion |
| Qwen3.5-0.8B | 32 / 32 † | 32 / 32 † | 329; 144 / 336 | 32 / 32 ‡ |
| SmolVLM-256M | 32 / 32 | 32 / 32 | 140; 417 | 32 / 32, loss span equal to completion |

† Checked through vLLM's renderer (`vllm_prompt_only.py`): the same `render_cmpl` → `process_for_engine` path that `LLM.generate` uses, without loading weights. The same path reproduced LFM's 64 / 64.

‡ For Qwen3.5 the token sequences are identical, but SFT's training template puts the 4 tokens `<think>\n\n</think>\n\n` inside the loss span. The generation prompt already contains them. This is benign (same sequence) and worth one line in the plan.

| Logprob parity (n = 8 at T = 1, processed logprobs) | Tokens | Median \|Δ\| | p99 \|Δ\| | Max \|Δ\| | p99 per-sequence \|ΣΔ\| | Verdict |
|---|---|---|---|---|---|---|
| LFM2.5-VL-450M, fp32 | 5,276 | 1.0e-7 | 4.2e-5 | 1.8e-4 | 2.4e-4 | **PASS** |
| SmolVLM-256M, fp32 (real items) | 4,092 | 4.0e-6 | 5.8e-5 | 1.5e-4 | 4.3e-4 | **PASS** |
| Qwen3.5-0.8B, bf16 | — | — | — | — | — | **inconclusive** |

- **Why Qwen3.5 is inconclusive:** the macOS CPU wheel's `_C` lacks `causal_conv1d_weight_pack`. I removed that wrapper at runtime (`vllm_parity.py`), but its `_C` also lacks `causal_conv1d_fwd_cpu`, which the Gated-DeltaNet forward pass needs. This is a CPU-only packaging gap; re-test H-A5-2 for Qwen3.5 on the first GPU box.
- **Leak check passed:** with `generation_config` left at "auto" and explicit `SamplingParams`, the fp32 parity at T = 1 rules out any leaked temperature, top-k or top-p.
- **Not tested:** H-A5-2(b), distribution parity at K = 64, needs an SFT checkpoint and a GPU.
- **CPU-run notes:**
  - run the engine in-process (`VLLM_ENABLE_V1_MULTIPROCESSING=0`, `VLLM_CPU_OMP_THREADS_BIND=0-1`), or vLLM CPU spawns extra processes and binds 7 cores;
  - cap the profiling dummy image with `limit_mm_per_prompt={"image": {"count": 1, "width": 512, "height": 512}}`, or Qwen profiles a 16k-token encoder pass;
  - `VLLM_CPU_KVCACHE_SPACE` must be an integer;
  - LFM's vLLM run peaked at 3.44 GB RSS.

### Q3: images (builder `b5-builder`, docker-container driver, capped at 2 CPUs and 4 GB)

| Image | Build | Size | Venv on disk | Checks |
|---|---|---|---|---|
| Phase 3 candidate (`Dockerfile.phase3`) | OK, 1,455 s | 4.34 GB | `/opt/venv` 8.6 GB | as in Q1 |
| Sampler stage (`Dockerfile.phase2 --target sampler`) | OK, 1,560 s | 4.14 GB | `/opt/vllm` 7.9 GB | torch 2.13.0+cu130, same arch list, vLLM 0.32.0, transformers 5.19.0, numpy 2.3.5, numba 0.65.0, imports in 109 s, no TRL or objdel |
| Phase 2 full (`--target phase2` = sampler + A3's `/opt/venv`) | not built (disk budget) | about 7.5 GB compressed, about 14.3 GB of venvs (estimated) | — | `/opt/venv` is A3's verified layer (3.61 GB, 6.4 GB venv) |

All `b5-*` images, the builder, its state volume and the `moby/buildkit` image it pulled are deleted. Only Docker objects I created were touched.

## What this changes in the plan

1. **Phase 2:** the image carries `/opt/venv` (as pinned, untouched) and `/opt/vllm` (`vllm==0.32.0` lock).
   - The sampler takes JSONL from the training env: chat text, image path, `mm_processor_kwargs` = `processor_kwargs` minus `padding_side`. It writes token ids and texts back; truncation, parsing and scoring stay in `/opt/venv`.
   - This is exactly the pattern tested here (`export_prompts.py` → `vllm_parity.py`).
   - No objdel in the sampler: it cannot share the numpy and numba pins.
2. **Phase 3:**
   - Colocate in one re-pinned venv, with the numba and numpy overrides written down in the lock.
   - Re-run the fast suite and `xarch_hashes --strict` on every re-pin.
   - Expect TRL's out-of-range warning until a TRL release declares vLLM ≥ 0.32 (PR #7799).
   - Re-evaluate when vLLM 0.33 is out: torch 2.14.0, transformers ≤ 5.20.x, so only a torch patch-level change for us. vLLM pins exactly 2.14.0, not our 2.14.1.
3. **Server mode** only on a 2-GPU box, and only with TRL ≥ the release that merges #7799. The trainer env would still need vLLM installed.
4. **First GPU session:** H-A5-2(a) for Qwen3.5-0.8B, plus H-A5-2(b) and H-A5-1 on CUDA (cheap, scripts ready).
5. **Rejected:**
   - TRL's in-range set, vLLM 0.30 + transformers 5.20: untested by vLLM above its own later ceilings;
   - the clean re-pin to numba 0.65 / numpy 2.3.5: caches have to be rebuilt, though targets are unchanged;
   - SGLang: no TRL, no SmolVLM, needs overrides.

## Lab files (`experiments/phase2/lab/B5/`, ruff-clean)
- `prereg.py`: the pre-registered criteria.
- `resolve_matrix.py`: the resolution matrix.
- `export_prompts.py`: prompts in the training env; run with `-I`.
- `vllm_parity.py`: the vLLM sampler, prompt ids and samples.
- `vllm_prompt_only.py`: prompt ids through the renderer.
- `hf_teacher.py`: HF teacher forcing.
- `Dockerfile.phase2`, `Dockerfile.phase3`, `build_images.sh`, `image_inside.sh`: images and in-container checks.

Outputs are in `<scratch>/phase2/B5/`: `data/*.json`, `logs/`, `resolve/`, `docker/ctx/{train,vllm,p3}.txt` with `p3.override`, and `docker/io/` (status, junit, hashes, freezes).

## Sources
- PyPI metadata: [vllm 0.32.0](https://pypi.org/pypi/vllm/0.32.0/json), [0.31.0](https://pypi.org/pypi/vllm/0.31.0/json), [0.30.0](https://pypi.org/pypi/vllm/0.30.0/json), [trl 1.16.0](https://pypi.org/pypi/trl/1.16.0/json), [sglang 0.5.22](https://pypi.org/pypi/sglang/0.5.22/json), [numba 0.65.0](https://pypi.org/pypi/numba/0.65.0/json).
- TRL v1.16.0 source: `trl/import_utils.py`, `trl/generation/vllm_generation.py` (`_init_vllm`), `trl/generation/vllm_client.py`; [TRL vLLM integration docs](https://huggingface.co/docs/trl/main/en/vllm_integration); [TRL PR #7799](https://github.com/huggingface/trl/pull/7799); [TRL PR #7567](https://github.com/huggingface/trl/pull/7567).
- vLLM: [PR #52183, torch 2.14.0 on main](https://github.com/vllm-project/vllm/pull/52183); [main requirements/cuda.txt](https://github.com/vllm-project/vllm/blob/main/requirements/cuda.txt); [main requirements/common.txt](https://github.com/vllm-project/vllm/blob/main/requirements/common.txt); [PR #59614](https://github.com/vllm-project/vllm/pull/59614); [issue #59562](https://github.com/vllm-project/vllm/issues/59562); [issue #59756](https://github.com/vllm-project/vllm/issues/59756).
- [NCCL "Creating a Communicator"](https://docs.nvidia.com/deeplearning/nccl/archives/nccl_2312/user-guide/docs/usage/communicators.html).
- [TorchRL VLLMWeightSyncScheme (NCCL vs CUDA IPC)](https://docs.pytorch.org/rl/main/reference/generated/torchrl.weight_update.llm.VLLMWeightSyncScheme.html).
- [uv dependency overrides](https://docs.astral.sh/uv/concepts/resolution/#dependency-overrides).
- [transformers v5.20.0 release notes](https://github.com/huggingface/transformers/releases/tag/v5.20.0).

(Report text returned by the B5 agent and saved by the coordinator; the harness blocked the agent's write.)
