# ObjectDelineation: plan vs actual (repo history map)

## Skills used

Invoked on the coordinator's order. Each skill is a method guide; none ran code, installed anything or spent money. Read-only.

| Skill | What it contributed here |
|---|---|
| nvidia-skills:nvidia-skill-finder | Catalog-check method. I did not run the catalog check myself. I used the main thread's lookup (MAIN_skill_lookups.md) and only recommend, never install. Relevant catalog matches from that lookup: nv-segment-ct / nv-segment-ct-finetune / nv-segment-ctmr (CT segmentation baselines), dicom-series-to-volume and dicom-series-preflight (relevant to the unbuilt 3D / 2.5D arm, section 7), i4h-workflow-train-rl, tao-generate-referring-expressions. Not used by the repo today. |
| engineering:code-review | Lens (correctness, performance, maintainability) for sections 5-6: the debts listed there (CK-4 slow 1,000-vertex batches, host syncs, KL logging spikes, no-merge sync fix) are the review-style findings the docs already carry. Verdict style for the repo: reward and checker are well tested (bit-exact, sanitizers clean); the weak spot is that docs lag (RESULTS.md has no Qwen rung). |
| engineering:debug | Reproduce / isolate / diagnose / fix structure, applied to the RL plateau: expected (IoU keeps rising) vs actual (+0.006 per 100 steps after step 100); isolated by offset; root cause centred training crops; ruled out memorisation, dead signal, LR decay, divergence; fix built, not yet run (section 5). |
| engineering:system-design | Framework for section 4: requirements (reward exact, scoring on GPU, one 32 GB card), data flow (vLLM generate -> one checker call for 256 answers -> trainer passes), trade-offs (sleep mode on costs 1.6 s/step; shared-prefix packing wins only below a cost threshold; compile and CUDA graphs rejected). |
| engineering:architecture | ADR-style reading of the decisions in sections 2-3. Implied ADRs with committed evidence: (1) Qwen3-VL-2B over LFM2.5-VL-450M after the coverage-gap failure; (2) CUDA checker over the CPU polycore, with polycore archived; (3) all scoring on device; (4) FastGRPO over stock TRL (80.3 -> 44.3 s/step); (5) no-merge LoRA sync (TRL #6688). Consequences: faster and exact, but a custom trainer coupled to TRL 1.17.1 and vLLM 0.32. |
| engineering:documentation | Principle "keep it current, link do not duplicate": drove the doc-debt list in section 6 (RESULTS.md says "in progress" and covers only LFM; Qwen rung only in impl/ISSUES.md R9-25; PLAN.md revision 9 vs ROADMAP drift). Writer should link to the cited files, not copy their tables. |
| huggingface-skills:huggingface-best | Method for model choice. From the main-thread lookup (live HF leaderboards): the Qwen3.5 small family (0.8B, 2B, 4B, 9B; natively multimodal, Apache-2.0, released 2026-02-27/28) is newer than the repo's Qwen3-VL-2B-Instruct (2025-10-19). The repo's committed PLAN revision 9 explicitly SKIPPED Qwen3.5-0.8B (owner, 2026-11-02, and commit 734b1ec removed the Qwen3.5 configs and extras). Qwen3.5-2B card (unverified column parse): RefCOCO 88.2. 2B and 4B fit one RTX 5090 with LoRA RL; 9B is tight beside vLLM. This is an open option, not something the repo tested. Not checked by me against the cards. |

Skills I judged not needed: huggingface-papers, trl-training, hf-cli (no new external facts needed; TRL facts come from the repo's own notes), tech-debt and testing-strategy (section 6 already is the debt list from committed docs), codebase-design (inventory in section 4 covers layout).

Addendum from the lookups (clearly external, not from the repo): option for a future ablation is Qwen3.5-2B or 4B as the base model. The repo evidence says the failure that led away from the small-model rung was data coverage, not model family.

---------------------------------------------------------------------------------------------------------


Sources: committed docs, code, and `git log` only (87 commits, 2025-03-20 initial, real work 2026-10-18 to 2026-11-08, branch gpu-scorer-open-once, head be52b1c).
NOTE ON SOURCES: per the coordinator's late instruction, I do NOT cite the owner memory folder or `.claude/`. I had read the memory notes before the instruction; every owner decision below is re-confirmed from a committed doc, code comment, or commit message (cited inline). Decisions I could only find in memory are listed in section 3c as "NOT confirmed in committed sources" and left out of the claims.
One uncommitted-but-on-disk file is used as a cross-check only: `experiments/phase2/vast/ledger.csv` (git-ignored via .gitignore:26). Committed spend figures come from docs.

---------------------------------------------------------------------------------------------------------
## 1. Timeline

| Date | Phase / step | What happened | Key numbers | Cost / GPU | Source |
|---|---|---|---|---|---|
| 2026-10-18 | Phase 0a census (Mac) | checker v0 (`src/objdel/checker`), census of expert disagreement | 0 certificate failures on 26,820 RIGA + 40,189 LIDC masks (67,009); LIDC topology disagreement 114/1,877 nodules = 6.1 % filtered (12 % raw); every nodule with a hole has a reader who drew none (65/65); RIGA cup band 56 %, vCDR spread >= 0.2 in 29 %; RIGA topology disagreement 0 % (by construction of the hull mirror, later confirmed from raw outlines) | $0, Mac | experiments/phase0/RESULTS.md; commit 29868b9 |
| 2026-10-18 | Phase 0b (Mac, inference only) | MedSAM + zero-shot Qwen3-VL-2B baselines | MedSAM LIDC IoU median 0.33 vs reader-vs-readers 0.74; 7.7 % of IoU>=0.5 outputs have wrong topology; where readers disagree, single-reader GT penalizes 46 %, majority 25 %, set-valued 11 %; Qwen3-VL-2B zero-shot polygon: 1/30 (disc) and 2/30 (cup) valid -> "SFT first" | $0 | experiments/phase0/RESULTS_0b.md |
| 2026-10-18 | ROADMAP written, repo restructured | 9-phase plan, ~150 GPU-h A100/H100 total, Track A alone ~70 | see sec 2 | - | ROADMAP.md; faaacc4 |
| 2026-10-19 | Phase 1 (Mac, no GPU) | C++ polygon core ("polycore"), reward v4, SFT target builder, data store, TRL wiring, eval harness; 136 break-it attacks -> 33 breakdowns -> fixed | gates K1..K13, R1..R10 pass or re-registered (K1: 0/849,802,293 px differ on 10,000 shapes; K13 coverage 93.03 %); store: 10,625 LIDC items (7,341/1,084/2,200 by patient 709/100/201), RIGA 750 photos 592/63/95; K12 smoke: SmolVLM-256M LoRA SFT 40 steps + 3 GRPO steps on Mac | $0 | experiments/phase1/RESULTS.md; commits ab1b0b9, 1c54bc4 |
| 2026-10-20 | Phase 2 plan + Stage 1 build | PLAN revision 8 (4 judge rounds), augmented views, SFT runner, vLLM bridge, eval gates, Vast tooling | - | $0 | experiments/phase2/PLAN.md; 7516c9f, d8891ce |
| 2026-11-01 | V0/V1 (Vast, first rentals) | V0 watchdog dry run on RTX 3090, V1 day 1 on RTX 4090: several failed launches (ssh key mode, 5.9 GB image layer never pulled on 6 of 9 hosts, Thai host) | image `p2-99c6...` 5.88 GB; V0 $0.2466 | small | experiments/phase2/V0_V1.md; commits 787c096, f6159..., 4371df5 |
| 2026-11-02 | R1 rung 1: LFM2.5-VL-450M full FT (RTX 4090 Quebec $0.401/h) | 1,542 stable steps at 0.69 s/step; selected e1; FAILED entry gates | median IoU_tc shifted 0.495, rescaled 0.092 (gate >= 0.60); hole recall 0.5 %; simple controls 93.7 % correct (gate 95 %); G-precision probe IoU 0.61 | billed $0.727 (R0 rerun + R1) | experiments/phase2/RESULTS.md sec 1 |
| 2026-11-02 | "Phase 2 reset" | owner decisions, CUDA checker committed, preflight gate, data rebuilt (lidc-bundle-v2: nodule < 10 % of crop 0.5 % -> 10.3 % of rows) | see sec 3 | - | commit 734b1ec; PLAN revision 9 |
| 2026-11-03 | R0 v3 CUDA checker on RTX 4090 (BC) | first GPU run of cuda/gpucheck.cpp | 195,809/195,809 cases bit-identical to CPU polycore; but SLOW: 54,710 polygons/s at 16 vertices, batch 1024, batch 1 only 134/s (GPU 0.002 to 0.19 x CPU path) | $0.249 | experiments/cuda_checker/RESULTS.md; ISSUES R9-17; commit 87c0fb3 |
| 2026-11-04 | scoring switched to CUDA checker only | batched `to_masks`, GPU build + check before every scoring chain | CPU score 0.06 ms vs HF generate 973 ms/answer: the switch was not a speed win | - | commit b905e7e; impl/ISSUES.md R9-19 |
| 2026-11-05 | R1 rung 2: Qwen3-VL-2B + LoRA r16, 3 epochs (lidc-bundle-v3, 98,804 rows) | trained, selected e2; FAILED the entry gates | shifted median IoU_tc 0.583 (< 0.60, pass share 61 % < 70 %), rescaled 0.458; G-precision probe IoU 0.64; hole recall 3-7 %; pass@8 0.91 vs pass@1 0.61; G-format/G-prior-b (+0.42)/G-topology-a/G-signal/G-vision/G-hf-check pass; e3 legal 89.1 % ineligible. Crash at step 170 (G-vision CPU/GPU device bug) cost ~25 min / ~$0.20 | billed $2.327 (+ $0.16 failed launch) | impl/ISSUES.md R9-25, R9-20..R9-23 |
| 2026-11-06/07 | P1: RL readiness on RTX 5090 (US, $0.694/h) | FastGRPO profile arms, CUDA checker on sm_120, TRL/vLLM check | FastGRPO 44.3 s/step vs TRL-as-shipped 80.3 (1.81x), on-policy (logprob gap mean 0.009 nats vs TRL's up to 25 nats, IS ratio 0.65-0.81); checker 7.74 M polygons/s at 16 vertices batch 1024 (141x its own first build; 87x polycore loop); 0 mismatches on 8,559 fuzz + all 13 defect classes; compile 47.9 s, CUDA graphs crash, no-checkpointing OOM | billed $1.253 (allocation $1.50) | experiments/phase2/perf/p1_results.md; commit 647ec01 |
| 2026-11-06/07 | K1: checker sanitizers on RTX 3090 | memcheck + initcheck | 0 errors on 40 GPU tests (compute-sanitizer 13.4.92; 13.0.85 lacked launchers) | $0.030 | perf/k1_results.md; e0fc6c4 |
| 2026-11-07 | R2 = the G-probe on e2 (RTX 5090 $0.563/h) | 20 GRPO steps | composite reward 0.694 / 0.679 -> 0.769 (+0.082, CI 0.067-0.097; noise 0.018), legal 94.5 % -> 98.2 %, 56.6 s/step | billed $0.550 | docs/research/notes/issues_solved_2026-11-07.md; explanations_2026-11-07.md |
| 2026-11-07 | Scoring moved wholly onto the GPU; polycore archived; OpenOnce + ImageMemo | same scores bit-for-bit as old CPU scorer (goldens + numpy oracle); 4,096 -> 32 image opens/step; 2.12 -> 0.78 s image prep | commit 9c5dd2e, 6161da8 | - | issues_solved_2026-11-07.md |
| 2026-11-07 | R3: GPU test of all of the above (RTX 5090 $0.544/h) | 8 GRPO steps end to end | 40 + 45 GPU tests pass, 0 mismatches (144,133 LIDC polygons, 37,093 masks), memcheck/initcheck 0 errors on 52 tests; 56.6 -> 36.3 s/step (time outside trainer step 22.6 -> 0.3 s) | billed $0.440 (budget < $1) | issues_solved sec 8 |
| 2026-11-07/08 | Full RL ("rl"): 3 epochs, 770 steps (RTX 5090 Oregon $0.610/h, 04:08-11:20 UTC) | image build `vllm_native-3ef6e9179051` (4 layers; one-layer 4.14 GB image timed out on 4 hosts); earlier launch attempts aborted | 27.5 s/step mean, 5 h 52 min, 197,120 answers, GPU use 67-73 %, peak 17 of 32 GB; train IoU 0.52 -> 0.77, F5 0.72 -> 0.94, topo 0.66 -> 0.93, legal 0.91 -> 0.99 | ~$4.56 (ledger invoice $4.362 early); Phase 2 committed spend $14.94 | docs/research/notes/full_rl_2026-11-08.md; commit 1870328 |
| 2026-11-08 | Held-out check (128 shifted-val items, K=8) | e2 IoU 0.511 / reward 0.694 / legal 94.5 %; step 250 0.634 / 0.802 / 98.8 %; final 0.659 / 0.817 / 98.8 % | e2 -> final +0.148 [+0.115, +0.178], 112 better / 14 worse | - | full_rl_2026-11-08.md |
| 2026-11-08 | R4 speed test (RTX 5090 CA $0.593/h) | no-merge vLLM sync + 64-wide chunks; vLLM sleep vs awake | 32.66 -> 29.94 s/step (8 %), identical rewards; handovers 1.6 s/step (sync 1.28, wake 0.12, sleep 0.25); awake saves 0.05 s | $0.65; committed $15.59 of $17.50 | notes/r4_speed_2026-11-08.md; commits 255cb67, 3fd83d8 |
| 2026-11-08 | R5 shared-prefix forward test | each prompt computed once per group of 8 answers (`--shared-prefix`, off by default) | 31.04 -> 23.65 s/step (-24 %), token ratio 0.36, step 1 identical | ~ $0.60 plan; first offer timed out pulling image ($0.24) | notes/shared_prefix_2026-11-08.md; be52b1c |
| 2026-11-08 | Final eval on whole splits (R5's box) | K=8, T=1, seed 11, shifted crops | TEST 2,200: IoU_tc 0.645, F5 0.816, P 0.695, R 0.748, Dice 0.701, legal 99.4 %; VAL 1,084: IoU 0.618, F5 0.806, P 0.689, R 0.716, legal 99.5 % | - | notes/final_eval_2026-11-08.md; 302ca3a |
| 2026-11-08 | RL plateau diagnosis; fix built, NOT run on GPU | see sec 5 | - | rl_shifted.yaml ~ $1.6 needs owner budget | notes/rl_plateau_diagnosis_2026-11-08.md; be52b1c |

Total GPU rentals: all on Vast AI. Committed Phase-2 spend $15.59 after R4 (R5 and eval add a little more; the $17.50 hard stop is in code: `experiments/phase2/vast/policy.py:58`). Total GPU time is order of 20 hours (rough: R1 LFM ~1.5 h, R1 Qwen ~5 h, P1 1.8 h, R2 1 h, R3 1 h, RL 7.2 h, R4 ~1 h, R5 ~1 h); I did not find one committed GPU-hours total.

---------------------------------------------------------------------------------------------------------
## 2. Planned vs actual

Plan docs: ROADMAP.md (2026-10-18, 9 phases, ~150 GPU-h on A100/H100; Track A alone ~70), experiments/phase1/PLAN.md, experiments/phase2/PLAN.md (revision 9, 2026-11-02, 115 KB).

| Planned item (where) | What happened | Why changed |
|---|---|---|
| ROADMAP Phase 1: Douglas-Peucker targets, Python/PIL polygon check, reward terms from research | Superseded by phase1/PLAN.md: potrace-style targets, C++ polygon core (polycore), reward v4 (topology-conditional IoU, self-calibrated band, rater-frequency topology term, area guard). ROADMAP itself says "the plan wins" | online research + sandbox tests + independent review (ROADMAP Phase 1 status box) |
| Phase 1 gates: native x86 timing + first GitHub CI run as Phase 2 entry steps | never fully done for polycore on x86 box; polycore was later archived | CPU checker replaced by CUDA checker (see below) |
| ROADMAP Phase 2: SFT rung screening SmolVLM -> LFM -> Qwen, A0 vs A0' hole arms, ~6 GPU-h | rung 1 LFM2.5-VL-450M trained 2026-11-02 and failed; Qwen3.5-0.8B skipped; rung 2 Qwen3-VL-2B (LoRA) trained 2026-11-05 and ALSO failed G-geometry/G-precision (median 0.583 vs 0.60). C2 (A0 vs A0', 3 seeds, c2_*.yaml) configs exist but there is no recorded C2 verdict | owner reset 2026-11-02 (PLAN rev 9); after e2 failed gates the project went straight to RL (owner override of 5.1 for the G-probe, perf/implementation.md "Stage G") |
| PLAN rev 8: gates must pass before the G-probe | e2 failed G-geometry and G-precision; G-probe ran anyway as "owner override of 5.1" | owner decision (perf/implementation.md:521) |
| PLAN rev 8: sampler = HF generate (D13) | R1 Qwen3-VL-2B sampled with vLLM, HF re-check on the selected branch | owner 2026-11-05, option A (ISSUES R9-20; r1_qwen3vl2b.yaml:7). Reason: HF projection 18.87 h vs 5 h stop (R9-18) |
| PLAN: 2 epochs, branches 0.5/1/2 | 3 epochs, branches 0.5/1/2/3 on lidc-bundle-v3 | owner 2026-11-05 (ISSUES R9-21; PLAN 4.9) |
| PLAN rev 8 fallback rules (hole oversampling x4 -> x1 etc.) | dropped | evidence: control failures were 4.3 % illegal / 1.3 % unreadable / 0.7 % false holes (RESULTS correction note; PLAN 6.3) |
| PLAN rev 8: G-precision capacity run after zero-shot fail | never run | R9-3 open |
| ROADMAP Phase 3: GRPO arms A1-A9 (random, format-only, IoU-only, IoU+binary validity, clDice, buffer, A7 single-GT, A8 ours, A9 policy-band), best-of-8, DPO, 3 seeds, ~50 GPU-h | NOT run. Only one RL recipe was run: the full composite reward (the A8-style reward) with LoRA; no controls, no ablation arms, no H1/H2/H3 tests | RL work turned into speed engineering + one full run; owner drives development; budget hard stop $17.50 vs ~50 GPU-h plan |
| ROADMAP Phases 4-6 (Track B natural images, Level-2 controller, Track C), 1c census | not started. A Phase-1 `Controller` class (best checkpoint, rollback, audit alarm) exists in `src/objdel/train/callbacks.py` but was not exercised in the full RL run (not evidenced) | out of scope after the deviations |
| ROADMAP Phase 7: 3D checker, 2.5D input arm (added 2026-10-21 / commit 4ec9f87) | only the checker side exists: `betti3d`, `topology_any`, `stack` (commit b004ca3 2026-11-06); no model, no 3D reward training, no 2.5D arm, no neighbour-slice re-fetch | see sec 7 |
| ROADMAP Phase 8: CUDA checker "later with lithography"; Python-hash identical port of bit-quads/pointer jumping | Built NOW (owner 2026-11-02) and became the core of the reward. Design differs from the roadmap text: `cuda/gpucheck.cpp` is a polygon legality (13 defect classes) + even-odd/nonzero raster + scorer (IoU_tc, F5, area, topology) engine, bit-exact to the CPU polycore (not "hash-identical to Python"). Lithography: not started | owner decision; CUDA checker later also took over all scoring (owner rule "checker stays on the GPU") |
| Roadmap compute: A100/H100 80 GB, SmolVLM or 2B | one RTX 4090 for SFT (rung 1 and 2), then RTX 5090 for everything RL; total spend ~$15-16 not ~150 GPU-h | smallest-model rule + measured step times; RTX 5090 chosen by owner 2026-11-06 (perf/issues.md:16, sec 5) |
| Roadmap RL step estimate: ~25 s per step on A100/H100 | TRL as shipped 80.3 s -> FastGRPO 44.3 -> (OpenOnce, ImageMemo) 36.3 -> full run 27.5 mean -> R4 29.9 vs 32.7 -> shared-prefix 23.6 vs 31.0 | custom trainer `FastGRPO`, found by profiling; P1 trace GPU busy 70 %, 7,233 `.item()` reads/step |
| Reward scoring: CPU worker pool (Phase 1/2 design, EV-5 pool sizing) | all scoring on the CUDA checker (IoU/F/area/topology on device, expert tables resident, pixels never return to CPU) | owner rule (docs/research/notes/gpu_scoring.md line 3: "the checker stays on the GPU"); commits b905e7e, 9c5dd2e |
| HF TRL 1.16 / separate vLLM venv (PLAN D15: two venvs) | TRL 1.17.1, `vllm_native` single-venv image, 4 layers (largest 1.54 GB, was 4.14 GB) | TRL 1.17 gave vLLM 0.32 support (commit 96ed7b7); large image layer timed out on 4 hosts (commit 0a4b93f) |
| PLAN: TRL's merge/unmerge sync | replaced by no-merge sync (`fast.lora_merged_params`) after the full run | TRL #6688 bf16 base-weight drift (6.4 % of entries changed after 770 cycles on a Mac test); fix is verified on GPU in R4 |
| Budget: Phase 2 alarm $10 / hard stop $15 (PLAN G10) | HARD_STOP now $17.50, ALARM $10 (`policy.py:58`) | owner raised in steps: cf2f507 "$15 -> $16.60: owner approved up to $7 more ($9.56 spent) for the full RL run"; 74a0135 $17.00; edf3129 $17.50 for the 10 h cap |
| Planned RL length: "300 steps" (roadmap) / G-probe 20 steps (PLAN) | 770 steps (3 epochs of 8,208 RL rows / 32 per step) | owner chose 3 epochs, 10 h cap, up to $7 (rl_full.yaml header, full_rl note) |
| Phase 9 evaluation/release | partial: `final_eval.py` (whole-split eval, bootstrap CIs), HF private repo svenk029/objdel-sft; no deployment guard, no data release, no paper | scope |

What is the ROADMAP's "Done when" vs reality: Phase 2 "valid >= 90 %, A0' recovers holes": legal is 99.4 % on test after RL; holes recall after SFT only 3-7 %, A0 vs A0' never decided. Phase 3 decision on H1-H3: not made.

---------------------------------------------------------------------------------------------------------
## 3. Owner's fixed rules and decisions (each with committed evidence)

### 3a. Confirmed from committed sources
1. Vast AI is the GPU provider; Mac only for data/checker/analysis. ROADMAP "Compute": "GPUs are rented on Vast AI (decided 2026-10-18)".
2. Smallest VLM that passes, screened smallest first, stop at first passing rung. ROADMAP Phase 2 text; PLAN D1/C1. Qwen3-VL-2B is rung 2.
3. 2026-11-02 after R1 (LFM) failed: do not retrain LFM; skip Qwen3.5-0.8B; next rung Qwen3-VL-2B (LoRA, one RTX 4090). PLAN "Revision 9 ... owner decisions ... binding" item 1; r1_qwen3vl2b.yaml line 1 ("owner decision 2026-11-02").
4. CUDA checker is built NOW, in parallel with Phase 2. PLAN rev 9 item 2; ROADMAP Phase 8 ("built now (owner decision, 2026-11-02)").
5. Preflight gate: `make preflight` must pass before any GPU rental; launch.sh refuses without a matching stamp; owner sees GPU-hours and dollars first. PLAN rev 9 item 3 and item 7; commit 734b1ec; perf/implementation.md lines 16-17.
6. No automatic fallback run after a failed gate; owner decides the single next run. PLAN rev 9 item 5, sec 6.3.
7. North America only (US/CA) for Vast hosts, enforced in code with no flag: `experiments/phase2/vast/policy.py:75` "REGION = ("US","CA") # owner (2026-11-02): North American hosts only; V1: a Thailand host could not pull the image"; commit 787c096. RIGA data never on a rented host (PLAN rev 9 item 7, D16).
8. Money rules: hard stop and alarm in code (`policy.py:58` ALARM 10.0, HARD_STOP 17.5 "owner 2026-11-07: full RL with a 10 h cap"); commit messages cf2f507, 74a0135, edf3129; per-test allocations ($1.50 for P1 perf/p1_results.md, $1 for R3 r3_gpu_test.yaml, ~$0.50 R4, ~$0.60 R5 configs).
9. Scoring all on the GPU: "Owner's rule: the checker stays on the GPU. Every check and every score runs in the CUDA checker" (docs/research/notes/gpu_scoring.md line 3). Also ISSUES R9-19 (owner 2026-11-04: exam and RL scoring use the CUDA checker only).
10. No HF generation in RL: "No HF generation anywhere (owner 2026-11-06): vLLM generates; the checker scores on the GPU" (configs/p1_rl.yaml:9); sampler vLLM for R1 (owner 2026-11-05 option A).
11. RTX 5090 for RL (owner 2026-11-06, perf/issues.md sec 5); TRL 1.17.1 (owner decision 2026-11-06, perf/implementation.md B1); `vllm_native` image name (owner 2026-11-06).
12. R1 Qwen: 3 epochs (owner 2026-11-05, ISSUES R9-21); owner override of the 5.1 gates to run the G-probe on e2 (perf/implementation.md:521); full RL: 3 epochs, 10 h cap, up to $7 (rl_full.yaml header).
13. Repo git rules in committed CLAUDE.md: no Claude co-author/attribution; commits authored by the repo owner only; `make test` and `make lint` pass before commit; owner runs git (perf/implementation.md:16); explanations in ASD-STE100 style (CLAUDE.md). Confirmed: git log shows only owner author (git user yashwanthsai1234).
14. Clean/simple code with no dead code: PLAN EXECUTION.md rule 3 (simplicity pass), commits that delete the replaced path (e.g. 734b1ec "Removed Qwen3.5 configs, R-diag fallback rules and their dead code"; 9c5dd2e archived polycore). Direct owner quote not in committed docs.
15. HF: private repo `svenk029/objdel-sft`, one branch/revision per run (r1-lfm, r1-qwen3vl2b, rl-full-3ep); ISSUES R9-24 says the HF token is user-wide and the box guard refuses it; `objdel.secrets` reads `<repo>/.secrets` before keychain; RIGA-trained weights never pushed (D16 / licence).
16. Owner runs GitHub image builds (perf/implementation.md B2 "owner runs it").

### 3b. Decisions visible as a result of owner choice but whose wording is only in memory
"From 2026-11-07 the owner drives development; ignore phase-2 PLAN.md" is NOT stated in any committed doc. Evidence in git: after 2026-11-06 commits follow ad-hoc owner-sized tasks (GPU scoring, OpenOnce, full RL, R4, R5) rather than PLAN stages, and the committed ROADMAP/PLAN still describe SFT arms / C2 / GRPO arms that were never run. For the writer: say "the plans were not followed after the R1 Qwen3-VL-2B gate failure" (supported by RESULTS/ISSUES R9-25 "open: owner decision" and the later commits), not that a stated rule exists.

### 3c. NOT confirmed in committed sources (leave out or flag as owner statement)
- That Qwen3-VL-2B was chosen over LFM2.5-VL-450M explicitly because of "coverage gap data": committed RESULTS correction note says the cause was data coverage (train rows with nodule < 10 % of crop 0.5 % vs rescaled exam 27 %), and decision 4 says the LFM is not trained again. This part is confirmed. The sentence "no LFM retrain" = confirmed.
- Agent-model-choice (Opus vs Sonnet), "one word per thing", HF token scoping history: memory only.

---------------------------------------------------------------------------------------------------------
## 4. Current system inventory

Language mix: `cuda/gpucheck.cpp` 1,835 lines (+ .h 150, py binding 227); src/objdel 11.6 k lines Python incl. cuda; `archive/cpu_checker` (polycore.hpp 713 lines C++ header-only, nanobind wrapper, tests), 802 tracked files.

| Component | Path | Runs on | Language | Maturity / evidence |
|---|---|---|---|---|
| Data store | `src/objdel/data/{lidc,riga,splits,store}.py` | CPU (Mac, one-time) | Python (pandas/pyarrow, pylidc parity) | Done, hashed (manifest 37fe7665...), K11 0 leakage; RIGA stays local |
| Targets / polygonizer | `src/objdel/targets/{frame,build,simplify,augment}.py`, `checker/polygonize.py` | CPU, offline (numba kernels) | Python/numba | Done; K5 topology kept 100 %, IoU_min >= 0.99; frames, shifted/rescaled views via `experiments/phase2/build_views.py`, `build_eval_views.py` |
| Output spec / parser | `src/objdel/geometry/{text,shape}.py` | CPU parser; polygon->mask via the checker | Python | Strict regex parser, MAX_CHARS 16,384; `to_masks` batch call into the GPU checker |
| CPU checker (polycore) | `archive/cpu_checker/` | not built; tests use recorded goldens `tests/goldens/polycore.npz` + `tests/oracles/*` (numpy/integer oracles) | C++17 | ARCHIVED 2026-11-07 (9c5dd2e). Was exact: K1 0/849,802,293 px |
| CUDA checker + scorer | `cuda/gpucheck.cpp,h,gpucheck_py.cpp`, `CMakeLists.txt` | GPU (also a CPU build as its test path); 2 CUDA streams, grow-only arena, pinned staging, warp-ballot raster, union-find CC | C++/CUDA (nvcc, nanobind), sm_120/86/89 | Bit-exact vs polycore on two architectures (0 mismatches RTX 4090, 5090, 3090); memcheck/initcheck 0 errors on 52 tests; 7.74 M polygons/s (5090). Open: 1,000-vertex shapes at batch 1,024 slower than 8-thread CPU (CK-4); 1.6 GB pixel copy |
| Reward / caches | `src/objdel/reward/{score,cache,advantages,spec}.py`; `train/rewards.py` | expert tables resident on GPU; per-rollout terms on GPU; `advantages_v4` group-normalised; `TRLReward` | Python | Reward v4: format gate, ring legality, IoU_tc (topology-conditional), F5 band, rater-set topology (b0,b1), area guard; K8-K10 gates; R4 neutrality ACCEPTED with mild pro-minority lean ~0.04 advantage units |
| SFT | `experiments/phase2/train_sft.py`, `train/{dataset,configs}.py` | GPU (HF Trainer/TRL SFT, bf16, LoRA r16 LM-only for Qwen3-VL-2B) | Python | Ran twice (LFM, Qwen). Best SFT checkpoint e2 not meeting gates |
| RL trainer FastGRPO | `src/objdel/train/fast.py` (732 lines); `experiments/phase2/train_rl.py`, `g_probe.py` | GPU (trainer); CPU thread for prefetch; vLLM colocated | Python on TRL 1.17.1 | Production: per-chunk trimmed forwards (64-wide), VisionCache, RopeCache, no old-logprob pass, length-sorted micro-batches, OpenOnce, ImageMemo, shared-prefix packed forward (flag `--shared-prefix`, off by default, GPU-tested R5), no-merge LoRA sync |
| vLLM bridge | TRL colocate + `FastGRPO` one-pad-per-image check; `train/sample.py`, `experiments/phase2/sample_vllm.py` | GPU | Python | vLLM 0.32.0, sleep mode ON (costs 1.6 s/step total, R4) |
| Eval | `src/objdel/eval/{metrics,gates,probes,signal,vision,audit,coverage}.py`; `experiments/phase2/{select_ckpt,final_eval,baselines}.py` | scoring via GPU checker; exam metrics partly copy pictures to CPU (few times per run, not per step) | Python | `final_eval` whole-split K=8 + bootstrap CIs; precision/recall vs strict-majority outline |
| Vast tooling | `experiments/phase2/vast/{offers,launch,sync_up,run,pull,teardown,watchdog,sweep_orphans,fetch_hf,push_hf}.sh`, `policy.py`, `bundle.py`, launchd plists | Mac | bash + Python | Guards G1-G10: US/CA only, price ceilings (3090 $0.50, 4090 $1.00, 5090 $1.50), ALARM $10, HARD_STOP $17.50, deadline watchdog, orphan sweep, image allowlist |
| Preflight | `experiments/phase2/preflight.py`, `make preflight CONFIG=...` | Mac (CPU dry run) | Python | Data integrity, patient disjointness, exam coverage, prompt, CT window, target legality, checker tests, CPU dry run, cost projection; stamps keyed on code/data/config hashes (commit 734b1ec) |
| Image | GHCR `objdel-gpu`, `vllm_native` (Dockerfile in experiments/phase2/docker) | GPU box | Docker, uv locks | 4 layers, largest 1.54 GB; built in GitHub Actions |
| Tests | `tests/` ~65 files; CI light job (no torch), GPU tests on box | | pytest/hypothesis | K13 coverage 93.03 % at Phase 1; Phase 2 tests added; OpenOnce/ImageMemo tests skip without training stack |
| 3D volume | `src/objdel/checker/volume.py` (64 lines) | CPU (scipy.ndimage) | Python | see sec 7 |
| Docs | `docs/research/{reports,notes}`, `docs/figures/*.svg|png`, README | | | README has figures, gallery (rl_gallery.png) |

Where each stage runs in a step: vLLM generates (GPU, ~7-12 s) -> checker scores all 256 answers in ONE device call (rewards wait 0.003-0.01 s) -> trainer old/ref/update passes (GPU). Idle GPU (67-73 % mean utilisation) is inside generate/batch phases, not in handovers.

---------------------------------------------------------------------------------------------------------
## 5. Final metrics, ceiling, failure modes, RL plateau

Final model = full RL run, step 770, merged (HF svenk029/objdel-sft revision rl-full-3ep), LoRA r16 on Qwen3-VL-2B, shifted crops, K=8, T=1, seed 11 (notes/final_eval_2026-11-08.md).

| Metric | VAL (1,084) | TEST (2,200) |
|---|---|---|
| Legal | 99.5 % [99.2, 99.7] | 99.4 % [99.1, 99.6] |
| IoU_tc (vs doctors' band, best of k doctors) | 0.618 [0.601, 0.635] | 0.645 [0.633, 0.656] |
| IoU_tc median item | 0.737 | 0.756 |
| Boundary F5 (a new doctor scores 0.95) | 0.806 | 0.816 [0.805, 0.829] |
| Precision vs strict-majority outline | 0.689 | 0.695 |
| Recall | 0.716 | 0.748 |
| Dice | 0.680 | 0.701 |

Before/after on 128 shifted-val items: e2 0.511 -> step 250 0.634 -> final 0.659; reward 0.694 -> 0.802 -> 0.817. The e2 -> test-split comparison was never measured (only on the 128-item probe).
Doctor ceiling (plateau note): mean over doctors of best IoU against the other doctors (leave-one-out, ignores topology): whole test 0.732; by offset 0.71-0.74. Model without misses scores 0.722 vs ceiling 0.740; 10.8 % of items are misses (IoU < 0.1).
Reader-vs-readers IoU yardstick from Phase 0b: 0.74 (all slices), 0.66 on topology-disagreement slices.
MedSAM given the expert box (R1 baselines): shifted 0.77, rescaled 0.80 (this is an upper-ish reference with the box handing over the location; not a model comparison).
SFT-era failure modes: crop-size/small-in-frame detection failures (<10 % of crop -> IoU 0.00), holes never drawn (even 111-200 px holes), illegal 4.3 % + unreadable 1.3 % on simple controls (LFM rung).
Remaining failure modes after RL: off-centre small nodules missed; outline spilling into chest wall (README gallery 10th-percentile crop); recall > precision (outlines slightly larger than majority); hole recall (SFT stage 3-7 %; not re-measured after RL in committed docs).

RL plateau (notes/rl_plateau_diagnosis_2026-11-08.md): training IoU 0.641 (steps 1-50) -> 0.733 (51-100) -> 0.753 (251-300) -> 0.771 (751-770); slope +0.48 per 100 steps early, then +0.005 to +0.008. Entropy 0.271 -> 0.062; advantage sd 1.04 -> 0.34; groups with no spread max 2.7 %; LR constant 5e-5; KL median 0.27-0.40; PPO clip never active. ROOT CAUSE: RL `grpo_rows` used only centred policy frames (object offset median 0.018 of side, 0.1 % > 0.2) while val/test shifted views have median 0.26 (68-70 % > 0.2); SFT did train on shifted+rescaled windows. Test IoU by offset: <0.15: 0.76-0.77 (>= ceiling); 0.2-0.3: 0.670 (5.9 % misses); >0.3: 0.545 (23 % misses; 47 % for small off-centre objects). Ruled out: memorisation, dead signal, LR decay, divergence.
Not-yet-run fix: `train_rl --train-views` (dataset.grpo_rows_views, serve.loader(*roots) with `<i>::<id>` keys): 8,208 centred + 8,208 shifted + 8,208 rescaled = 24,624 rows, each window scored on its own rater masks. Config `experiments/phase2/configs/rl_shifted.yaml`: 300 steps from rl-full-3ep with `--shared-prefix`, then whole val and test eval; ~2.6 h, ~$1.6 at ~$0.60/h; headroom estimated ~ +0.09 IoU_tc. Needs owner budget: Phase 2 committed $15.59 of $17.50 (rl_shifted.yaml / plateau note say $0.60 left, written before the $17.50-15.59 = $1.91 reading; follow the note's number only with caution, rl_shifted needs ~$1.6 so it is tight). Tests written (Mac): tests/test_p2_views.py, test_p2_train_rl.py, test_profile.py. Train windows built on the Mac in 2 min 11 s.

---------------------------------------------------------------------------------------------------------
## 6. Open issues / known debts (committed docs)

- G-precision capacity run was never run (ISSUES R9-3); R1 rung 2 gate fail decision "open: owner decision" (R9-25); no C2 verdict (A0 vs A0'); `c2_*.yaml` configs unused.
- ISSUES R9-6: select_ckpt K=16 and rotated check cost 41 % of samples (open); R9-9 GRPO memory for Qwen3-VL-2B (later measured in P1: peak 19 GB).
- Held-out e2 baseline on the whole test split not measured, so the full RL gain on whole test is unmeasured.
- Phase-1 holdovers: first GitHub Actions run / native x86 timing "Phase 2 entry steps" (I34); R3(b) holes-stamped adversary wins 61.6 % (limit 60 %) ACCEPTED, monitored (I15/I16); R4 neutrality mild pro-minority lean accepted (I13).
- Checker: CK-4 1,000-vertex batch-1024 shapes slower than the 8-thread CPU path; 4 of ~7,600 test cases had Linux vs Mac cos/sin rounding differences, now falling back to exact integer oracles (commit 196ea3f).
- Trainer: shared-prefix is off by default (flag); TRL weight-sync drift fixed locally (TRL #6688, PR #7427 open); KL log shows single-token spikes from vLLM-vs-trainer sampling mismatch (judge by median); `CUDA Error: invalid argument at cumem_allocator.cpp:249` printed at exit by vLLM sleep allocator (exit code 0).
- Speed debts (notes/issues_solved_2026-11-07.md, cuda_streams_graphs_review, rl_code_optimality_review, grpo_utilization_research): vLLM's own image pass (pass 3) unchecked; 1.6 GB pixel copy; exam metrics still copy pictures to CPU; per-step host syncs (e.g. 384 `lengths.tolist()` in vision blocks, 112 loss `.item()`); idle 27-33 % of GPU inside generate/batch phases; beta=0 (drop reference pass) listed as an owner decision item; manual CUDA graph capture (F4) untried; regional torch.compile off; fixed cost per micro-batch dominates after packing (b ~ 3-5 s).
- Vast ops: direct SSH port can reset (use ssh2.vast.ai proxy port); big image layers time out (hosts cannot resume); launchd sweep job not loaded on Mac at one point.
- Data/licences: RIGA-derived tables and anything trained on RIGA stay local; RIGA never used for the RL/SFT in Phase 2 (LIDC only) - the project headline now mentions RIGA but the trained model is LIDC-only (check README claims: "Input: CT slice or fundus crop" vs Phase 2 PLAN D16 / S6 removal). RIGA was removed or triggered in Phase 2 (PLAN 6.4).
- Never reached: GRPO ablation arms, hacking alarms in a live run, Level-2 controller, natural-image tracks, 3D training, lithography.
- Doc debts: README results table cites SFT checkpoint stats only for LFM; Qwen3-VL-2B SFT (R1 rung 2) result lives only in impl/ISSUES.md R9-25, NOT in experiments/phase2/RESULTS.md (RESULTS.md still says "in progress", section 1 is LFM only).

---------------------------------------------------------------------------------------------------------
## 7. 3D voxel work: what exists

- Code: `src/objdel/checker/volume.py` (64 lines, commit b004ca3 2026-11-06 "3-D voxel topology: stack slices, betti3d (b0, b1 via Euler, b2), set-valued topology_any; Phase 0 census imports again"). CPU only (scipy.ndimage), no CUDA path.
  - `stack(slices)` places 2-D masks keyed by z into a voxel volume, missing slices empty.
  - `euler(vol)` = V - E + F - C of the closed unit-cube complex (each cell counted once when any cube containing it is present).
  - `betti3d(vol)` -> (b0, b1, b2): b0 = 26-connected foreground components, b2 = 6-connected background components not reaching the border, b1 = b0 + b2 - chi.
  - `topology_any(pred, readers)` = predicted (b0, b1, b2) equals some reader volume's.
- Tests: `tests/test_volume.py` (hollow cube one cavity, ring one tunnel, corner contact joins two parts, stacked slices of a hollow nodule).
- Used by: Phase 0 census (`experiments/phase0/lidc_census.py`, 3D clusters, `b0_3d_*`, `cavity_any_*` in clusters.csv): 3D component-count disagreement 16/1,877 nodules (0.9 %) filtered; slices with differing (b0,b1) 2.9 %. Phase 1 had removed `betti3d/FG8/BG4` (M3) and commit b004ca3 restored a 3-D version.
- NOT built: ROADMAP Phase 7 certificate (chi from 256 bit-cube patterns vs surface Euler number V-E+F; mesh volume = voxel count), 3D reward in training, whole-nodule outline over slices by the VLM, 2.5D input arm (z-1,z,z+1 as RGB) (ROADMAP 2026-10-21 / commit 4ec9f87) and the neighbour-slice re-fetch it needs, CUDA 3D checker, H6/H6b.
- README section "3D volumes" describes slice -> stacked voxel -> b0/b1/b2 -> `topology_any`. It is an analysis/reward-definition library, not part of the training loop.

---------------------------------------------------------------------------------------------------------
## 8. Quick caveats for the report writer
- Cost numbers: committed docs give RL ~$4.56, committed spend $14.94 (after full RL) and $15.59 (after R4). Local git-ignored ledger invoices list: v1 $0.627/$0.727, v2 $0.153, v3 $0.249, r1 $2.327, p1 $1.253, k1 $0.030, r2 $0.550, r3 $0.440, rl $4.362 (late accrual), r4 $0.519 - use committed numbers in the report.
- GPU-hours: no single committed total; steps x seconds numbers above are measured. RL alone: 770 x 27.5 s = 5 h 52 min.
- SFT step times: LFM 0.69 s/step; Qwen3-VL-2B 3.40 s/step projected (ISSUES R9-18), 3,065 steps, ~5.34 h projected (R9-21); actual R1 billed $2.327.
- Dates: Phase 0/1 2026-10-18/19; Phase 2 build 10-20; V1 11-01/02; LFM fail 11-02; CUDA checker GPU 11-03; Qwen SFT 11-05; RL-ready/P1 11-06/07; G-probe, GPU scoring, R3, full RL 11-07/08; R4, R5, final eval, plateau diagnosis 11-08.
