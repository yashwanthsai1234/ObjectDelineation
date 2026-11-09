# Final review r1: fresh outside architect (whole Implementation.md, 11,117 lines)

Reviewer stance: no earlier opinion. There are no earlier `final_r*` files, so nothing below repeats a fixed item.
Phase reviews (`reviews/{A..CLEAN}_phase_review.md`) were not re-raised; every item here is cross-phase or
front/back-matter, or a defect that is still in the round-1 text.

## Skills used

| Skill | What it changed in this review |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | Read the live catalog (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`). `compute-sanitizer-diagnose-cuda-bugs` is confirmed. The catalog also lists `nv-segment-ct`, `nv-segment-ct-finetune` and `dicom-series-to-volume`, which the front matter misses (item 14). Nothing was installed. |
| `engineering:code-review` | Checked the code claims against the files: `configs.py:73-77,177-196`, `targets/build.py:53,93`, `certify.py:46`, `train_rl.py:50`, `run_box.sh:24`, TRL `grpo_config.py`, transformers `training_args.py:860` (`max_grad_norm` 1.0). |
| `engineering:debug` | Traced the `adv_scale` claim down to the optimiser. Grad norm over 770 steps is 0.0018 to 0.022, median 0.0048 [measured: runs/rl/train_rl/rl.json], so clipping never fires (item 3). |
| `engineering:system-design` | Followed the critical path from data to Gate B1 under the $3 limit. Found the ordering problem (item 1) and the unowned decode tail (item 2). |
| `engineering:architecture` | Checked the ADR numbers and the decision tables across phases. Found clashes (items 5, 7). |
| `engineering:documentation` | Checked front-matter facts against the phase sections. Found conflicts in gates, costs, the diagram and layer placement (items 8, 9, 10). |
| `huggingface-skills:huggingface-best` | Checked section 3 against the HF API: Qwen3-VL-4B-Instruct `safetensors.total` = 4,437,815,808 (matches 4.44 B). Checked the EGM table numbers against arXiv HTML (item 11). |
| `system-design-skills:system-design` | Ran the clarify, estimate, stress-test loop on the paid path: what B1 needs, and what G1 buys (item 1). Treated the train/eval shift law as an assumption to stress-test (item 4). |
| `engineering:tech-debt` | Listed section 4 "ADD" items that have no task in any phase. These would become instant docs debt (item 6). |

Web and package checks: PyPI `trl` latest 1.17.1, `vastai-sdk` 1.11.0, `import-linter` 2.16 (all match the doc).
arXiv 2601.13633 (EGM) and 2604.27720 were read.

## Verdict

The layered target (geom < verify < data < policy < rl/sft < infra), the one-source CUDA engine with a CPU oracle,
the move-only Phase E and the data-first diagnosis are sound. A top lab would build the same parts. A top lab would
change three things. First, it would spend the scarce GPU dollars on the one experiment that tests the diagnosis
before any infra validation. Second, it would attack the largest measured idle block, the decode tail. Third, it would
not tune the training shift law to one exam law without a check on a different law. The rest are consistency defects
that block a clean implementation.

## Items

### 1. BLOCKING. 00_front (sections 1, 6.2, 6.3) and 99_back (owner decision 1): the paid order puts infra validation before the data fix and breaks the "dataset quality first" rule

Problem. G1 is the first rental. It spends $1.62 (cap $1.74) of the $3. Its items 2-7 (K1, K2, K3, D `gpu-check`, C
speed arms; $1.04) do nothing for Gate B1. Item 8 (A3) is not needed for Gate B1 either. Gate B1 is paired against
R5's existing `runs/r5/eval/eval_val/items.jsonl` (B.6.2), and A3's saved answers feed only B3.7. Because G1 comes
first, section 6.3 concludes that "the planned 300-step B1 run does not fit". Owner decision 1 then recommends a
shortened B1 (160-190 steps). The highest-priority experiment gets cut to pay for kernel checks. By the doc's own
numbers, a B1 rental alone on today's code costs 30 + 300 x 27.5/60 + 6.9 = 174 min, which is $1.74 at $0.60/h or
$1.86 at $0.641/h (base case). At 23.65 s it costs 155 min, $1.55 to $1.66. B.R item 11 already allows B1 on the
control numerics (bf16-seqmask), so B1 does not need the probe first.

Exact fix.
- Section 6: rename the rentals. **R1 (first) = B1 run**: today's code with `--shared-prefix`, bf16-seqmask, 300
  steps, val eval, `BudgetStop` at $1.90. If the ledger allows, add A3 (`a3_measure.yaml`, 17.5 min) to the end of the
  same session; this saves a second launch (about 40 min, $0.40). **R2 = G1 without A3**, run with the money left.
  Its items are ordered by what the next run needs: C `c_speed_check`, then the C arms, then D `gpu-check`, then
  K1/K2/K3.
- Section 1 diagram: move "G1 validation rental" after "G2 B1 run". Change the text "C code is used by the B runs only
  after..." so that it applies to B2 and B3 only.
- 99_back decision 1: change the default to "(b') keep $3: B1 300 steps first (about $1.7-1.9), then whatever G1
  items fit; B2/B3 wait for a budget decision". Change decision 2 to match.
- Section 7 R1: change the mitigation to the new order.

### 2. MAJOR. C (C.8.3, C.R item 23), B, D: the largest measured idle block, the decode tail, has no task

Problem. The BRIEF measures 4-5 s of decode tail in a 23.4 s step. It also measures that 38 % of R5 steps hit 1,536
tokens; those answers are runaways and get masked. C hands the tail to "Phase B/D" (line 4633). B has no decode-tail
task beyond the C7 limit. D has nothing on decoding. A lab would treat this as the first speed item, because it is
bigger than all of C1-C5 together (2.25-4.7 s, C.8.1).

Exact fix. Add task **C11 "Runaway cut" (Mac first, $0 to size)**.
1. On the A3 saved answers (or on R5 `log_completions`), compute for every image `cap_i = min(L, ceil(1.5 x
   max doctor target tokens of image i) + 64)`. The doctor token counts are already in A8/B1.3 `tokens()`. Report the
   share of *found* answers (IoU >= 0.1) longer than `cap_i`. Gate: < 0.5 %.
2. If the gate passes, `rl.trl_adapter` passes one `SamplingParams` per prompt with `max_tokens = cap_i`. (vLLM
   `LLM.generate` accepts a sequence of `SamplingParams`. Test this on the pinned vLLM with the tiny model.)
   `mask_truncated_completions` stays True.
3. Add an NVTX sub-range `generate/tail` (A6) and the counter `completions/capped_share`.
4. Speed gate: tail seconds (C0 step 6b) fall by >= 50 %. Step-1 rewards stay equal up to the capped rows.
5. In C.8.1, record it as a separate arm next to C7. Remove the line "the measured tail goes to B (limit, grammar) and D".

### 3. MAJOR. B (B3.4, B3.7 rule 5, ADR-0013): the `adv_scale` reason is wrong under AdamW, and "RLOO" is a constant factor of mean-centring

Problem. B3.4 says that without `adv_scale`, "B3 would take steps several times smaller". The optimiser is AdamW
(`configs.py:100`; TRL default) with `max_grad_norm` 1.0 [code: transformers training_args.py:860]. The measured grad
norm never comes near the clip: 0.0018-0.022, median 0.0048 over 770 steps [measured: runs/rl/train_rl/rl.json].
AdamW's step does not depend on a constant gradient scale, except through `eps` = 1e-8. C7 step 2 of the same doc
says exactly this ("AdamW is invariant to a constant gradient scale except through eps and gradient clipping"). So
`adv_scale` changes only the weight of the KL term relative to the policy term. The effective beta is beta /
adv_scale, and with C8 (beta 0) `adv_scale` does nothing. Also, `R_i - (sum R - R_i)/(G-1) = G/(G-1) (R_i - mean
R)`. With G = 8, "RLOO" is mean-centring times 8/7. The real change from v4 is "no per-group std division", not the
estimator.

Exact fix.
- B3.4: remove `adv_scale`. Write `A_i = R_i - mean_g(R)` and name it "mean-centred, no std division (equals RLOO up
  to 8/7)".
- B3.3: drop `adv_scale` from `RewardV5`. Keep the constructor rule for `c_topo`.
- B3.7: replace rule 5 with "report median |A_v4| / median |A_v5| and the implied effective beta; if the owner wants
  v4's KL balance, set `beta = 1e-3 / ratio` in `b3_v5.yaml` (a config value, not a reward constant)".
- B.3.4 ADR-0013, the glossary row and the test "adv_scale multiplies A" (B3.6): update to match.

### 4. MAJOR. B (B1.4 step 1, B.6.2): the training shift law is fitted to one exam law, with no check on another law

Problem. B1.4 sets the target cell shares from the val-shape half of the *val shifted view*. That view is built by
`build_eval_views` with one shift law. All B gates are also scored on views from the same builder. A model can pass
Gate B1 by learning that builder's offset law, not translation robustness. B.2 cites Chen et al. ("the test risk when
the test law is the augmentation law"), which is exactly that coupling. A top lab would add a held-out law the
training never targets.

Exact fix.
- B1.6: add `--kinds uniform` to `build_eval_views`. This builds a val view with offsets drawn uniformly over the legal
  support (every doctor mask inside the window, same sides), with seed 1 and a different RNG stream. Freeze it like the
  test shifted view (D15).
- B.6.2 Gate B1: add a report-only line with a hard floor: "val-select IoU on the uniform view: CI lower bound of
  (B1 - R5) > -0.02".
- B.8.1: add the eval cost, about 2 x 7 min, R5 and B1 [estimate from 415 s].
- B1.4: state in one sentence that the val view is the target law *by choice*, and that the uniform view is the guard
  against fitting it.

### 5. MAJOR. 00_front (section 4, "EDT on the eval path") vs D (D.2, D7): opposite decisions on cuCIM

Problem. Section 4 says "REPLACE scipy EDT with cuCIM/CuPy ... Custom K2 EDT is the fallback". D.2 says cuCIM is "Not
adopted for production" (the Mac has no CUDA, it is a large dependency) and makes the K2 EDT the production path, with
cuCIM as an optional box-only oracle. An engineer gets two opposite instructions.

Exact fix. Change the section 4 row to "KEEP our K2 separable EDT (one source, CPU oracle on the Mac) for eval and
reward. cuCIM `distance_transform_edt(float64_distances=True)` = optional second oracle in one box test (D7). Trigger
to revisit: EDT time > 10 % of eval wall time". Change the 99_back line "cuCIM/CuPy EDT and labels for eval metrics
(D5)" (line 9875 region in CLEAN) to match.

### 6. MAJOR. 00_front (sections 2.1, 2.3, 2.4, 3, 4): ten architecture parts are "ADD" with no task in any phase

Problem. These parts appear only in the front matter (checked with grep: no task line after line 430):
- the view bank (`data/viewbank.py`; C calls it "optional T1");
- `uv.lock` with `uv sync --locked` in CI;
- pyright;
- `ruff format`;
- pre-commit with the `no-ai-attribution` hook;
- `events.jsonl`;
- `data_revision` with a private HF dataset repo;
- checkpoint `resume` (adapter + optimizer + step every K steps);
- the `linear_attention` guard in `FastGRPO` (section 3: "Guard now, because it is cheap");
- Trackio `report_to`.

Requirement 2.1 ("Pixels never come back to the CPU after the view bank is loaded") depends on a part nobody builds.
The view bank is also sized for the old fixed 22,023-view set (2.4: "~2.2 GB"), while B1 builds a new scheduled rows
root of 9,600 windows (1.2-1.6 GB, B1.4 step 5) and B1.8 retires the old views.

Exact fix. For each part, either assign an owner task with files and a test, or mark it "DEFERRED: not in this plan"
with its trigger:
- `linear_attention` guard → E5 (in `rl/trl_adapter.py`, unit test with a fake config);
- pre-commit and `no-ai-attribution` → DOC12;
- `uv.lock`, pyright, `ruff format` → CLEAN9 or E0;
- `resume` → B1.5 (the doc already relies on "HF resume skips consumed batches", line 3061);
- Trackio, `events.jsonl` → C0 or DEFERRED;
- `data_revision` and the HF dataset repo → B1.4 step 5 (which already uploads the rows root);
- view bank → DEFERRED until C1 measures PNG decode > 0.2 s per step. Until then, reword 2.1 as "after the images are
  decoded and uploaded (25 MB per step)".
Fix 2.4 and section 4 to name the B1 rows root, not the 22,023-view bank.

### 7. MAJOR. B, D, E, 99_back: ADR numbers clash

Problem. `docs/architecture/decisions/` ends at 0011. B assigns 0012 = number match, 0013 = reward v5, 0014 = scheduled
rows. E assigns "ADR 0012" = keep `objdel.infra` in the package (lines 6511, 7294, 7483; owner decision 19), "ADR
0013" = live runs use the vllm_native image / spec result (lines 7064, 7205), and "ADR 0015" (lines 7125, 7379). D's
G-3 also cites "ADR 0012" (line 6295). Two different files would get the same number.

Exact fix. Add a registry table to front-matter section 4 or 8:
- 0012 number match (B2);
- 0013 reward v5 (B3);
- 0014 scheduled RL rows (B1);
- 0015 sampler venv / vllm_native (E4);
- 0016 infra inside the package (E8);
- 0017 run spec result (E6);
- 0018 library keep list (CLEAN14);
- D's decision record (D9) gets 0019.
Then renumber every E, D and 99_back reference to match.

### 8. MAJOR. 00_front (2.2, 2.3) vs B (B1.3, B.8.4) vs E (E3.1, E3 move table): layer placement conflicts that break the exhaustive import contract

Problem.
- `coverage`: B1.3 moves `eval/coverage.py` to `src/objdel/data/coverage.py` (line 2574), and front 2.2 lists coverage
  in data. E3.1 (line 6641) puts coverage in verify, and the E3 move table maps `eval/coverage.py` -> `verify/coverage.py`
  (line 6946). After B, that source file no longer exists. A verify module that uses `data.windows` would break the
  contract.
- NVTX marks: front 2.2 has them in `core/marks.py`, but the 2.3 tree also lists "infra/ ... nvtx (marks)".
- 2.3 has both `verify/store.py` and `data/store.py`, with no rule for which loads what.

Exact fix.
- One placement table in section 2.3: `data/coverage.py` (B1.3 wins); `core/marks.py` only (delete "nvtx (marks)" from
  infra); `verify/store.py` = doctor-store preload/loader, `data/store.py` = source store.
- Change E3.1, the E3 move row (`data/coverage.py` stays; no move) and the import-linter contract test I1 to match.

### 9. MAJOR. 00_front (sections 1, 6.3, 7 R3) vs B (B.4 order, B.8.1) vs C (C.6, C.8): gates, order and cost basis disagree

Problem.
- Front Gate C is "<= 18 s/step", and 6.3 says "20 s (Gate C passed)". C's speed gate for C1-C5 is G-SP1 = 20.0 s;
  18 s is G-C, which needs C7 and C8, and those are owner decisions inside Phase B.
- Front 6.3 plans G2 = B1 alone and puts the B2 probe "outside the $3". B.4 orders session S1 = B2.3 probe *then* the
  B1 run, with the probe winner.
- Front uses $0.60/h and 20 or 23.65 s. Section 7 R3 says "Cost table uses 27.5 s as base case". B.8.1 uses 27.5 s and
  $0.641/h. C.8.2 uses $0.593/h.
An engineer cannot tell which rental to plan or which gate unlocks the faster code.

Exact fix.
- Front: "Gate C = G-SP1 (<= 20.0 s mode A / 20.7 s mode B) for use in B runs; G-C (<= 18 s) applies only after
  C7+C8 are approved".
- One cost table in 6.3, base 27.5 s and $0.641/h, with 23.65 s and 20 s as columns (the B.8.1 layout).
- State that B1 runs on control numerics without the probe (B.R item 11). Change B.4 "session S1 (B2.3 probe, B1 run)"
  to "S1 = B1 run; S2 = B2 probe + B3 run".
- C.8.2: re-cost at $0.641/h.

### 10. MINOR. 00_front (section 2.5): the step diagram shows default-off parts as the default

Problem. The diagram draws "32 x n=8 requests" (C6 is a "measured arm, default off") and "fused log-prob K3b". By the
section 5 rule, K3 lands only if the torch-op C1 misses its gate, and K3b has no landing task at all.

Exact fix. Draw 256 single requests (today), and mark n=8 as "[C6 arm]". Replace "fused log-prob K3b" with
"log-prob (torch; K3b only if a bench wins, no task yet)".

### 11. MINOR. 00_front (section 3): citation details

Problem. The EGM numbers are in Table 4, not Table 2 [web: arxiv.org/html/2601.13633]. Also, 89.6 / 91.0 are trained
from the Thinking bases (83.6 / 87.2), so the like-for-like gap goes 3.6 -> 1.4, not 3.1 -> 1.4. The abstract of
2604.27720 covers MedGemma and Lingshu (all systems 0.05-0.24 against a centre box at 0.10). It does not name
HuatuoGPT-Vision.

Exact fix. Write "Table 4"; write "84.2/87.3 Instruct, 83.6/87.2 Thinking bases; 89.6/91.0 after SFT+RL; gap 3.6 ->
1.4". Cite HuatuoGPT-Vision's exclusion to its own source (size / gating), or drop it from the 2604.27720 claim.

### 12. MINOR. B (B1.5, B.6.2): Gate B1 mixes "more steps" with "new data"

Problem. B1.5 does not state B1's start checkpoint. B2.3 and B3.9 imply rl-full-3ep, so B1 is R5 plus 300 more steps.
A paired win over R5 then includes the effect of extra training.

Exact fix. In B1.5, state `--ckpt rl-full-3ep` in `b1_rows.yaml`. In B.6.2, add one line: "Control for extra steps:
R5's held-out IoU plateau over the last 300 of 770 steps (runs/rl/train_rl/rl.json eval series) is the 'more steps on
old data' baseline; quote its slope next to the B1 result".

### 13. MINOR. A / B gates: the CIs cover data variance only, not training-seed variance

Problem. Every paired CI resamples nodules from one training run (seed 0). Run-to-run RL variance is unknown. A
+0.01-0.02 gain could be seed noise.

Exact fix. In A1 and B.6.2, state "CIs exclude training-seed variance". Add a $0 proxy: A3 runs the R5 val eval twice
with two sampling seeds, and the doc reports the paired difference as the sampling-noise floor. Gates B1/B3 must
exceed that floor.

### 14. MINOR. 00_front (Skills used row 1, 99_back open questions): the NVIDIA catalog has CT parts the doc says are absent

Problem. The front matter says only `compute-sanitizer-diagnose-cuda-bugs` fits. The live catalog (Vision AI) also
lists `nv-segment-ct` / `nv-segment-ct-finetune` and `dicom-series-to-volume`. These fit the stated next direction
(3D voxels) and could serve as an independent CT reference segmenter.

Exact fix. Add one open-question row: "3D direction: evaluate `dicom-series-to-volume` for the voxel loader and
`nv-segment-ct` as a reference segmenter (recommendation only; install needs owner approval)". Change the Skills row to
match.

## Counts

blocking 1 (item 1); major 8 (items 2-9); minor 5 (items 10-14).
