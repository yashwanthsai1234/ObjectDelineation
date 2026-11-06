# Phase 2 Stage 1 issues

Break-it findings (wave 2). Format: id, severity, component, repro, expected vs actual, owner.

## Seeded by the coordinator from wave-1 reports

| id | sev | component | issue | owner |
|---|---|---|---|---|
| S1-1 | medium | M1 views | Under scale jitter, 11 of 24,701 LIDC rows (side ≥ 100) are certified but not pixel-exact against their mask (plan says 100 %); root cause unknown | V |
| S1-2 | low | M1 / M2 | `build_eval_views` calls the private `dataset._prepare(job)`; a signature change in T's file silently breaks it | V / T |
| S1-3 | low | M4b | Only 71 lung-background patches on the probe-test patient half (B3 used 200); topology probe draws on 19 val patients | P |
| S1-4 | low | M4b | MedSAM / SAM 2.1-tiny baseline paths never executed (no downloads in Stage 1) | P (Mac baseline step) |
| S1-5 | info | M7 | Unverified until V0: self-destroy endpoint, `VAST_API_KEY`, `destroy -y`, invoices format, onstart length limit | O (V0) |
| S1-6 | info | M1 | RIGA brightness/contrast/gamma jitter not implemented (Phase 2 runs are LIDC-only) | V (if RIGA permitted) |
| S1-7 | medium | M2 / M8 | `view_rows` has no `epoch` column; train_sft parses it from `row_id` (contract gap, `dataset.py:278`) | T |
| S1-8 | medium | M1 / M8 | `build_views` defaults to `--epochs 1`, R1 needs 3 (2-epoch branch decays into epoch 3); default and configs must agree | V / R |
| S1-9 | low | M4 / M8 | `vision.check` leaves the model in eval mode; caller must restore train mode | E |
| S1-10 | low | M2 / M8 | `configs.sampling` uses vLLM field names; g_probe maps to GRPOConfig names — one mapping should live in configs | T |
| S1-11 | medium | M3 | macOS: a pool worker importing a module after fork crashes the pool; imports must precede `fork_pool` (undocumented, untested) | S |
| S1-12 | medium | M8 | Not done: K = 16 + HM-IoU16 on the selected branch, rotated-view attribution check, WiSE-FT report, automatic D12(c) repair loop; `--items` takes the first N not a seeded subset | R |
| S1-13 | high | M2 | Qwen3.5 full FT loads bf16 base weights: no fp32 master weights for AdamW (updates below bf16 resolution are lost). Load fp32 + `bf16=True` autocast for full FT; bf16 base only for LoRA | T |
| S1-14 | medium | tests | `tests/test_breakit_reward.py::prepared()` writes no `manifest.json`; 2 tests fail after the I45 manifest-driven `cache_loader` | T |
| S1-15 | medium | M4 | `tests/test_p2_eval.py::test_nsd_equals_the_deepmind_reference` fails in the full suite (passes alone? order / randomness) | E |
| S1-16 | low | M5 | hf_teacher logprob parity (H-A5-2) not promoted; R0 needs it | T / R |

## Break-it round (wave 2): BK1 7, BK2 12, BK3 16, BK4 20 findings

Details, repro tests and suggested fixes: `impl/BK1.md` … `BK4.md`; strict-xfail tests in
`tests/test_p2_breakit_{data,scoring,runners,ops}.py`. High: BK2-1, BK2-2, BK3-1, BK3-3, BK3-4, BK4-1 … BK4-5.

## Fix round (wave 3): assignment by file owner

| Fixer | Files | Issues |
|---|---|---|
| FX-train | `train_sft.py`, `configs/*.yaml` | BK3-1 (yaml side), BK3-2, BK3-3, BK3-4 (yaml), BK3-5, BK3-8, BK3-10, BK3-12, BK3-13, BK3-14, BK4-5, BK4-6, S1-8 |
| FX-select | `select_ckpt.py`, `g_probe.py`, new `hf_sample.py` | BK2-1, BK2-2, BK2-5, BK2-6, BK3-4 (select), BK3-6, BK3-9, BK3-11, BK3-15, BK3-16, S1-12, BK3 simplicity cuts |
| FX-T | `train/{configs,dataset,sample}.py`, `sample_vllm.py`, T's tests, `test_breakit_reward.py` fixture | BK3-1 (sft_config), BK2-3 / BK3-7, BK1-3 (check side), BK1-4, BK1-5, BK4-1, S1-7, S1-10, S1-13, S1-14, S1-16 |
| FX-V | `targets/*`, `build_views.py`, `build_eval_views.py` | BK1-1, BK1-3 (write side), BK1-6, BK4-2, S1-2, S1-8 (default), BK1 simplicity cuts |
| FX-O | `vast/**`, `docker/**` | BK4-3 … BK4-19, BK1-2, BK1-7, BK4 simplicity cuts |
| FX-EGS | `eval/{gates,metrics,vision,signal}.py`, `train/serve.py`, `prereg.py`, `c2_analyze.py`, `eval/probes.py` (dedupe only), `tests/test_data_lidc.py` | BK2-4, BK2-7 … BK2-12, S1-9, S1-11, S1-15, BK4-20, c2/probes bootstrap dedupe |

## Status after the fix round (2026-10-20)

All 55 break-it findings and the seeded S1-* items are fixed, each with a test that failed first (no strict xfail
remains in `tests/test_p2_breakit_*.py`), except:

| id | state | reason / where it goes |
|---|---|---|
| S1-2 | open (low) | `build_eval_views` still calls `dataset._prepare` on one line; needs public `dataset.prepare_record` / `set_summary` helpers |
| WiSE-FT (S1-12 part) | deferred | report-only, runs only when the D12(c) repair fires; needs a run-chain stage after g_probe |
| S1-3, S1-4 | open | probe patch supply (71 vs 200); MedSAM / SAM 2.1-tiny paths first run in the Mac baseline step before V1 |
| S1-5, BK4-18, BK4-19 | UNVERIFIED until V0 | Vast self-destroy endpoint, `VAST_API_KEY`, invoice format, bandwidth field names, /dev/shm size |
| S1-6 | deferred | RIGA brightness / contrast / gamma jitter (LIDC-only until RIGA permission) |
| policy.py size | open (low) | 524 lines vs the ~380 target; the fixes added the safety checks |

Coverage config: `source_pkgs = ["objdel"]` added to `[tool.coverage.run]` (forked reward-pool workers otherwise traced
tests/ and experiments/, reporting 77 % for a library that is at 94 %).

## V0 (2026-11-01): RTX 3090 watchdog dry run, instance 58715600, offer 50693169 at $0.230/h

| Item | Result |
|---|---|
| Offer filter | Fixed (4f6159a): Vast sends `verification: "verified"` (no `verified` bool) and `cpu_ram` in MB, 1 GB = 1000 MB |
| vastai errors | Fixed (4f6159a): exit 0 with empty stdout, error JSON on stderr; offer IDs change between searches (`--pick top`) |
| Image | Rebuilt on GitHub Actions (a6d8fb3, `image.yml`): objdel-gpu:p2-99c6eb484c88, 5.88 GB compressed (was 7.57), causal-conv1d compiled, 2.58 GB deduped. Pull + start on the 3090: 5.7 min |
| Self-destroy (S1-5) | VERIFIED: cap 21:39:46 UTC, account empty at 21:40:50 with no destroy from the Mac; the sweep then closed the run |
| Invoice (BK4-18) | `show invoices-v1 --charges` lists one Instance row, $0.022 for the ~11 min; credit 4.660 -> 4.639 |
| On-box checks (BK4-19: /dev/shm, `CONTAINER_API_KEY`, watchdog.log, GPU imports) | NOT RUN: `~/.ssh/id_ed25519` has a passphrase and is not in the agent (the box accepted the key ~1 min after `running`). Owner runs `ssh-add --apple-use-keychain ~/.ssh/id_ed25519`, then a short re-check before V1 |
| Cap vs pull | The 10-min cap counts the image pull (5.7 min); V1 caps must allow ~6 min of pull |

## V1 day 1 (2026-11-02): RTX 4090, Texas instance 58762863, offer 41247264 at $0.440/h; paused by the owner

| Item | Result |
|---|---|
| Launch | 5 boxes. Thailand and Nevada (twice): layer `f5e58e423bc0` (5.61 GB of 5.88, both venvs) kept retrying, > 15-30 min; Quebec: pulled in ~5 min but ssh refused (`authorized_keys` modes, vast-cli #336); Texas (6.7 Gb/s down): pulled in ~5 min. Fixes: US/CA-only hosts as a fixed admissibility rule (owner rule), onstart key repair, teardown always destroys (f6f5eb1, 787c096). Open (low): split the venv layer if slow pulls recur |
| Box bugs | Fixed with tests (77ed632): Vast's login tmux broke `ssh host CMD` (`~/.no_auto_tmux`); `git ls-files` split paths with spaces; two-source rsync made the `.sha256` the destination; `run.sh` silently ignored flags after the config |
| R0 | image_check pass (165 of 168 bf16 TFLOPS, /dev/shm 50 GB); xarch 23/26 identical + accepted (the `eval` diff was a stale Oct-19 Mac reference, refreshed); K3 idle done; LFM timing 0.68 s/step (6 workers, det off), 16.4 GB peak; stop rule GO (R1 0.67 h vs 3.4 h) |
| Sampler parity (H-A5-2) | FAIL for both: prompt ids 32/32, median token diff LFM 1.4e-5 / Qwen3.5 1.8e-3 nats, but per-sequence p99 1.11 / 1.44 > 0.5. D13 fallback applied: `sampler: hf` in R1/C2 configs (b37917b). Owner question: is 0.5 nats on a summed bf16 sequence (~160 tokens) too strict? |
| GDN route (D15) | Qwen3.5 gated delta rule runs the PyTorch reference path although fla imports; causal-conv1d via hub kernels |
| Qwen3.5 arm | OOM at step 0 (full FT, batch 8, no checkpointing, > 23.4 GB): D6 says LoRA r16 for that rung; s/step not measured |
| R1 | Stable phase (rolling checkpoint 1,542) and branch e0.5 saved, then branch e1 OOM at step 0: the branch loop kept every earlier trainer (model, grads, fp32 AdamW ~7 GB) alive. Fixed with a test (trainers freed before each branch). Checkpoints (29 GB) not pulled: R1 reruns from scratch (~20 min); manifest kept in `runs/v1/r1-lfm-oom/` |
| Cost | Texas box $0.627 (billed); credit 4.59 -> 3.96 |

## Pre-rerun scan (2026-11-02): six reviewers (train_sft, HF sampling, Vast tooling, R0/analysis, CPU chain dry run, dead code), all fixed test-first

| Finding | Effect if unfixed | Fix |
|---|---|---|
| Branch decays and crash resumes reset the LFM vision tower to base (checkpoint keys `vision_tower.vision_model.*` vs live `vision_tower.*`; Trainer's non-strict load only warned) | every LFM branch = trained LM + untrained tower; G-vision blind on branches | `_load_from_checkpoint` via `from_pretrained` (key conversion), missing / unexpected keys raise; real-weights test |
| HF evaluation as chained: no `--items 256`, one item per `generate`, timeouts restart from zero | R1 eval ≈ 14 GPU-h (9-23) vs ~1 h planned: past the box cap and the credit | chain passes `--items 256` for `sampler: hf`; `hf_sample` batches 16 items, resumes whole batches; `sample.sample` batches items; SIGTERM kills the sampler group |
| HF sampling dtype: from_pretrained "auto" kept fp32 for trained branches, bf16 for base weights | R0's rate undercounts R1 by up to 2x | the GPU always samples in bf16 (D13's dtype) |
| `project_r1` eval allowance fixed at the vLLM 0.3 h | GO not shown for the HF path | R0 step `hf_rate` (32 items x K 8 per model); `project_r1 --hf_rate` from the chain's real item-samples; settings checked against the arm; exit 2 = bad input |
| Qwen3.5 configs `lora: false` after R0's OOM | Qwen3.5 R1 OOM at step 0 | `lora: true` (D6) in r1_qwen35 and r0_qwen35_arm |
| `rolling/` (full AdamW checkpoints) kept to the end; stable/e* copied | ~33 GB LFM, ~60 GB Qwen3.5 full FT on an 80 GB disk | stable/e* hard-linked; `rolling/` deleted once every stable/e* exists, `decay/<b>` once branches/<b> does |
| No `PYTORCH_ALLOC_CONF`; one peak for the whole run | fragmentation OOM (1.15 GB reserved unallocated at R1's crash); peaks hidden | `expandable_segments:True` (D4); `peak_memory_gib` per phase |
| Resume edge cases | lost stable copy found hours later; partial timing overwrote the first; an aborted decay saved as a branch | refused up front; first timing kept; no branch saved on abort |
| Vast: no excludes on pull; PULLED_OK stamped after rsync; region only a search filter; teardown missed Ctrl-C; `timeout` without `-k`; create/poll orphan window; onstart size unchecked; box with no job idles to the cap | GBs of checkpoint bandwidth per pull; results lost to the watchdog; a non-NA host launchable; boxes billing to the cap | `PULL_EXCLUDE` (rolling, decay, stable) shared with the watchdog; `PULLED_OK.next` → `PULLED_OK`; US/CA in `reasons()` (no `--region`); teardown retries then always destroys; `timeout -k 120`; launch cleanup covers create, poll tries 3; onstart gzip + length check, run_id ≤ 40; watchdog destroys a never-started box after 2 x idle |
| `k3_harness` exit 0 on a failed gate; `c2_analyze` read the wrong report key; rx8 / R-diag configs without `sampler: hf`; g_probe wrote bare NaN | gates that cannot fail; C2 crash; vLLM after its parity failure | fixed; R-diag views now buildable (`build_views --no-dihedral / --no-window-jitter`, same draws except that component) |
| Dead / stale code | — | removed `serve.ref_loader` and root-less `CacheRef`, `configs.DTYPE`, unbuilt C++ files; test-only probe answerers moved to `tests/probe_models.py`; one `dataset.prompts` for select_ckpt's rows; duplicate fixture and `.gitignore` merged; docs corrected |
| Open (low) | — | onstart gzip not yet seen on a real box (check in the first minute); spend estimate counts 15 GB download per launch (conservative); GHCR PAT on the Mac's `vastai create` argv |

## V1 rerun (2026-11-02, box 58836372, RTX 4090 Quebec, $0.401/h, billed $0.727)

- Launch: three offers vanished between search and create (no cost); a New York 5090 was still loading at the 15 min
  launch timeout and was destroyed (~$0.20 estimated). The gzip watchdog onstart ran on the real box.
- R0 rerun: G-vision reference (LFM) and hf_rate: 0.132 s/sample on base weights (5.9 % of samples hit 1024 tokens);
  project_r1 GO 2.70 h. Trained branches sampled at ~0.048 s/sample, so R1 took ~1.5 h, not 2.7.
- R1 training finished: 1,542 stable steps at 0.69 s/step; peaks stable 20.8, e0.5 17.8, e1 20.8, e2 20.9 GiB (no OOM);
  G-vision ok at 386 / 771 / 1542 (ratio 1.00 / 0.96 / 0.84).
- select_ckpt crashed at its CPU G-vision check: transformers 5.20 runs LFM2's short conv through the CUDA-only
  `causal_conv1d` on any device (R0's reference runs only the tower, the Mac has no such package). Fixed (024089f):
  select_ckpt blocks the package. Select + probe rerun from the saved branches.
- Selection (D12): e2 ineligible (legal 87.0 % < 90 %, though pass rate 0.69 vs e1 0.48 on 256 items); e1 selected.
- e1 gates: FAIL. G-format, G-signal, G-vision, G-prior-b pass; G-geometry fails (median IoU_tc shifted 0.495 < 0.60,
  rescaled 0.092, rescaled legal 71 %); G-precision and G-topology-a probes fail. R2 G-probe skipped (gates failed).
- Branches on the Mac: runs/v1/r1-lfm/branches/{e0.5,e1,e2}. Open then: why rescaled collapses (answered below); next rung (decided below).

## Rung 1 post-mortem and the revision 9 review (2026-11-02)

Cause of the rung-1 failure and the owner decisions: RESULTS §1 (correction note) and PLAN revision 9. The review found
these problems. "Fixed" means fixed in this change; owner names the file's owner.

| id | sev | where | finding | state |
|---|---|---|---|---|
| R9-1 | high | D10 / `augment.draw` | Train scale jitter did not cover the rescaled exam: nodules under 10 % of the crop were 0.5 % of train rows, 27 % of the exam. No check compared the training distribution with the exam views before renting | data work: `lidc-bundle-v2`, `eval/coverage.py` |
| R9-2 | high | PLAN §6.3 rule 1 | Hole oversampling x4 to x1 targets G-topology (a), but the control failures are 4.3 % illegal, 1.3 % unreadable, 0.7 % false holes | rule and its three configs dropped; `prereg.fallback_component` and its test removed |
| R9-3 | high | G-precision | PLAN §5.1 says "zero-shot fail, then the capacity run". No runner config or chain stage trains it; `gates.phase2_entry` counts a missing capacity result as a fail. Rung 1's G-precision verdict is zero-shot only | RESULTS correction note. Open: build the capacity stage only if rung 2 fails zero-shot (`probes --n-train 20000` and the scorer slot exist) |
| R9-4 | high | `hf_sample` default 16 items x 8 | Qwen3-VL-2B KV cache is 112 KB per token: 128 sequences x 1,856 tokens = 27 GB. One runaway answer (1,536 tokens) overflows a 24 GB card | fixed: default 8 items (13 GB) |
| R9-5 | medium | `r1_*.yaml` micro-batch 8 | With a 152k-token vocabulary the logits of one long batch are about 12 GB; LFM (65k) peaked at 20.8 GiB without them | fixed: Qwen3-VL-2B uses 4 x 8; R0 arm measures the peak on the longest batch (length grouping puts it first) |
| R9-6 | medium | `select_ckpt` chain | K = 16 and the rotated check are report-only but cost 26,016 of 62,944 samples (41 %): about 1.8 h of HF sampling for Qwen3-VL-2B | open: `select_ckpt` / `vast/policy.py` (PLAN §4.9) |
| R9-7 | medium | `project_r1` | Stop value was 2 x a hand-written high estimate, blind to the credit | fixed: explicit stop value per model (5.0 h for Qwen3-VL-2B) |
| R9-8 | medium | LoRA path | `train_sft` LoRA (branch decays, merge, G-vision on a PEFT model) was covered only by a tiny LFM with a flag, never by the real architecture | fixed: slow CPU smoke trains a tiny Qwen3-VL (LoRA, checkpointing, G-vision, merged branches) and runs G-probe on it; `hf_sample` batching is tested on it too (left padding, different image grids). `select_ckpt` on it is not smoke-tested (a random model never stops: 1,536 tokens x hours on a CPU) |
| R9-9 | medium | G-probe on Qwen3-VL-2B | GRPO micro-batch 8 (= G), completions up to 1,536 tokens, 152k vocabulary, no gradient checkpointing in `grpo_config`: memory not estimated | open: measure before R2 (only if the gates pass) |
| R9-10 | low | configs | `lora`, `force_lora`, `TEMPLATES`, `training_template`, `_BF16_BASE`, `qwen3.5-0.8b` entries, 4 timing arms, Qwen3.5 arm and rung configs, 3 R-diag configs: dead after the decisions | removed with their tests |
| R9-11 | low | `docker/` | Image holds Qwen3.5-only extras (`fla`, `causal-conv1d`, `kernels`) and `image_check.sh` records the GDN route; `build_views --no-dihedral / --no-window-jitter` and `test_leave_one_out_views_change_only_that_component` served only R-diag | fixed: extras, the `ccd` stage, `patch_arch.py`, the GDN record, `select_ckpt`'s causal_conv1d block and the R-diag flags removed with their tests |
| R9-12 | low | docs | PLAN header, D1, D6, D13, D15, D22, run graph, §4.7–4.9, §6.1–6.4, §8; RESULTS "Open"; V0_V1 Part 4 and 5.7; EXECUTION waves; ROADMAP Phase 2 / 8 / GPUs; README status | fixed |
| R9-13 | high | probe-precision | Probe-test backgrounds came from LIDC train patients (split by a patient hash): 43 of them are in SFT training (found by preflight) | fixed: probe-test patches from val patients (40 patches, 19 patients); R1's probe archived in `data/eval_views_archive/probe-precision-r1` |
| R9-14 | high | image (`Dockerfile.phase2`) | One 5.61 GB layer (both venvs, hardlinked) of the 5.88 GB image: 6 of 9 launches never finished the pull (Thailand, Nevada x2, Texas offer, New York, 2026-11-03 box 58931942 at $0.437/h, ledger $0.13) | fixed: each venv in four layers (cuBLAS + cuDNN, other CUDA wheels, torch, rest); CI image `objdel-gpu:p2-e392c8ccfdf9` (`sha256:e1b373e1…`): 18 layers, 7.32 GB, largest 1.45 GB (was 5.61 of 5.88) |
| R9-15 | medium | `vast/bundle.py` | The LIDC bundle shipped every `views/lidc-*` (the R1 view too, ~260 MB of upload on a paid box) | fixed: `--view NAME`; the bundle also carries `checker_eval/` for the CUDA checker step at the end of R0 |
| R9-16 | high | Vast host (machine 28360, Ontario, $0.464/h, box 58936618) | R0 `image_check`: 7.4 of 165 bf16 TFLOPS. The 4090 was at 100 % and 234 W with none of our processes, power limit 150 W. The check ran only after the 12 min bundle upload (billed $0.153 plus late charges). The new image pulled in under 2 min there | fixed: `sync_up` runs `image_check.sh` after the code and before the bundle; `policy.BAD_MACHINES` rejects 28360 |
| R9-17 | medium | `CMakeLists.txt` (CUDA build) | R0 v3: `make cuda` failed, "CUDA_ARCHITECTURES is empty for target gpucheck" (default `native` set after `add_library`); the docker check passed archs on the command line | fixed: property set on the target, test reads the file; rerun on the same box: 21 GPU tests pass, eval 195,809 of 195,809 exact. Speed fails (GPU 0.002 to 0.19 x the CPU path): `experiments/cuda_checker/RESULTS.md` |
| R9-18 | high | R1 cost (`project_r1`, R0 v3) | Qwen3-VL-2B R1 projects 18.87 h > stop 5.0 h: training 1,947 steps x 3.40 s = 1.84 h (peak 5.9 GiB), HF eval sampling 62,944 samples x 0.973 s = 17.0 h. vLLM parity fails on the sequence sum only (token median 4.6e-6, token p99 0.11 nats, sequence p99 1.79 > 0.5) | decided 2026-11-05: R1 samples with vLLM, see R9-20 |
| R9-19 | medium | scoring (`objdel.reward`, `objdel.eval`) | Owner 2026-11-04: the exam and RL scoring use the CUDA checker only (although scoring is 0.006 % of an answer's time: CPU score 0.06 ms vs HF generation 973 ms) | done: `geometry.to_masks` makes the masks of many objects in one `_gpucheck.check_polygons` call; `score_many` / `score_batch` / `score_completion` / `eval.metrics` / `TRLReward` all use it (no switch). `policy.chain` starts every scoring chain with `run_box.sh --build-only` and `tests/test_gpu_scoring.py` on the GPU build, so a failed build or a difference from polycore stops it before training. `selected.json` records `_gpucheck.build()`. polycore stays for building targets and as the tests' reference |
| R9-20 | high | R1 sampler (PLAN D13 revision 10) | Owner 2026-11-05 (option A), before any R1 result: R1 Qwen3-VL-2B samples with vLLM. The sequence p99 of G-parity becomes report-only; ids and token median stay gates. R0 v3 scored on the Mac: HF and vLLM base-weight answers both 0 on every metric (220 / 224 wrapped in a ```json fence, the rest truncated), tokens 581 vs 568 (CI [-69, 101]): no evidence either way | done: `select_ckpt.hf_check` re-samples the judged branch with HF (64 seeded shifted-val items x K 8) and gate G-hf-check needs 0 inside the 95 % CI of every vLLM - HF metric difference. `project_r1` uses the R0 v3 vLLM rate (0.09 s per answer, 95 s start; 256 answers in 23 s) plus the check (505 s): 4.04 h <= 5.0 h (was 18.87 h). `configs/gpu_scoring.yaml` removed (R1's first stage does that check) |
| R9-21 | medium | R1 depth | Owner 2026-11-05: train 3 epochs. The e3 branch decays to 3.45 epochs; `lidc-bundle-v2` holds 3 | done: `lidc-bundle-v3` = `build_views --epochs 4` (epochs 0-2 and every context byte-identical to v2; 98,804 rows, 558 MB); R1 branches 0.5 / 1 / 2 / 3, C2 configs on v3. Projected 3,065 steps x 3.40 s + vLLM exam 2.45 h = 5.34 h; the stop rule (`project_r1.R1`) rises from 5.0 to 6.0 h |
| R9-22 | high | G-vision on CUDA (`objdel.eval.vision.next_token_probs`) | R1 2026-11-05, box 59299927: crash at the first rolling checkpoint (step 170): `mm_token_type_ids` extended with CPU zeros (the target ids stayed on the CPU); G-vision had only run on the CPU (R0 reference, Mac tests) | fixed: target ids moved to the model's device; test on the `meta` device; R1 restarted from step 0 (lost about 25 min, ~$0.20) |
| R9-23 | high | vLLM sampler (`sample_vllm.py`) | R1: engine start failed: FlashInfer's top-k sampler JIT-compiles with nvcc, absent from /opt/vllm (R0 sampled with --logprobs, another path) | fixed: `VLLM_USE_FLASHINFER_SAMPLER=0` before vLLM is imported (PyTorch sampler, same distribution at top_k 0 / top_p 1); `select_ckpt` rerun alone without `--rotated` (report-only) to fit the box deadline. vLLM rate measured: ~0.04 s per answer (estimate 0.09) |
| R9-24 | medium | HF push (`hf_push.base_card`) | Qwen/Qwen3-VL-2B-Instruct has no LICENSE file (card: apache-2.0): the push failed at 404 | fixed: `experiments/phase2/licenses/Apache-2.0.txt` (apache.org, sha256 cfc7749b...) when the card declares apache-2.0; e2 pushed from the Mac to svenk029/objdel-sft@r1-qwen3vl2b. The HF token is user-wide (box guard refuses it); `objdel.secrets` reads <repo>/.secrets before the keychain |
| R9-25 | — | R1 result (C1 rung 2, Qwen3-VL-2B, 3 epochs, vLLM) | Selected e2 (eligible e2, e1; e0.5 format 97.8 %, e3 legal 89.1 %). Gates FAIL: G-geometry (shifted median IoU 0.583 < 0.60, pass share 61 % < 70 %; rescaled 0.458, 46 %), G-precision (probe IoU 0.64). Pass: G-format, G-prior-b (+0.42), G-topology-a, G-signal, G-vision, G-hf-check (vLLM - HF within +-0.035, 0 in every CI). pass@8 0.91 vs pass@1 0.61; hole recall 3-7 %. Box $2.33 billed (+ $0.16 failed launch) | open: owner decision (§6.3, §8.2) |
