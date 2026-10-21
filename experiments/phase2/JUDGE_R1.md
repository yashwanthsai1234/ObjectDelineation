# Judge verdict, round 1: Phase 2 plan (draft C, 2026-10-20)

Independent review of `experiments/phase2/PLAN.md` against `RESEARCH.md`, notes A1–A6, lab reports B3/B4/B5 (and
their JSON), B1 `prereg.py`, `experiments/phase1/RESULTS.md` §2/§8, `experiments/phase1/PLAN.md`, `ROADMAP.md`, the
source (`targets/frame.py`, `eval/gates.py`, `train/configs.py`, `data/lidc.py`), and external spot checks (PyPI
metadata for trl 1.16.0 and vllm 0.32.0; the Qwen3.5-0.8B and LFM2.5-VL-450M cards; arXiv 2606.18487, 2510.01624,
2601.19228). B1/B2 are judged on their decision rules only.

## Overall verdict: APPROVE WITH CHANGES

| Dimension | Score (1–10) | One line |
|---|---|---|
| Scientific validity | 7 | Crop-leak fix and shifted-val gates are right; C2's synthetic probe is under-calibrated for the hole claim; seed handling and the GRPO-readiness KPI are under-specified |
| Evidence | 8 | Nearly every number traces to a note/JSON and checks out; one citation is reversed (Quagmires), one mislabelled (B3 P1 vs P2), one cost basis mismatched (2 vs ~4 epochs) |
| Engineering | 8 | Two-venv image, JSONL bridge, fork-before-CUDA pool, hashes and guardrails are tested and sound; one GPU-sequencing gap (trainer vs vLLM at 0.85) and one unnecessary preload |
| Feasibility / cost | 8 | ~50 GPU-h / ~$47 is credible and guarded; S0 replaces the estimates; the §4.9 unit arithmetic needs the WSD epoch count |
| KPIs | 6 | Gates are measurable and mostly justified, but "a model GRPO can improve" is only proxied (G-signal); G-topology (b) is vacuous; probe primary metric and seed pooling are not defined |
| Simplicity | 5 | 10 milestones, 9 A/Bs, 8 sessions and 4 new metrics to answer two questions (C1, C2); E2/E4/E7/E9, nnU-Net and K = 64 evals can be deferred without losing either verdict |

Scores are against the owner's bar (research-driven, test-evidenced, goal-driven, KPIs, sensible Vast use, C++ vs
Python justified). The plan clears that bar on method; the changes below are about measuring the two deliverables
directly and trimming scope.

## MUST-FIX (blocking approval)

1. **C2 is a hole-only claim on an uncalibrated synthetic probe (§1 C2, D20, §5.2).**
   Problem: A0 is "A0-multi" (parts kept, D2), so A0 and A0′ differ *only* in holes; "non-simple hit" paired
   difference is driven entirely by the 160 one-hole items. B3 P2 calibrated the precision generator's *outline*
   irregularity (disc IoU 0.822 vs 0.806) but nothing calibrates hole or part geometry to real LIDC non-simple
   targets, and A1 §1.4 says most real LIDC holes are tiny (the size rule is `max(5 px, 0.5 %·|U|)`). A probe
   with large, obvious holes would pass A0′ while saying nothing about real cavitations.
   Why it matters: C2 is the headline Phase 2 claim and feeds ROADMAP's "A0′ recovers holes and parts".
   Fix: (a) in M4, draw hole area / part-area ratio from the empirical distribution of the 414 non-simple train
   targets (census.json), with a pre-registered minimum at the size rule; (b) report non-simple hit stratified by
   hole-area bin (≤ 10 px, 10–30, > 30) and require the C2 CI on the pooled set *and* no bin below hole recall 0.30;
   (c) rewrite C2 and ROADMAP Phase 2 "done when" as "recovers holes (parts are already expressible in A0-multi)";
   (d) keep the 44-item real val slice descriptive as planned.

2. **No thresholded KPI for "a checkpoint GRPO can improve" (§5.1, D12).**
   Problem: G-signal (duplicates ≤ 10 %, σ-share ≥ 50 %) is a necessary proxy; the 20-step GRPO probe is only a
   selection input with no pass rule. A checkpoint can pass G-signal and still sit below the vanishing-gradient
   threshold of 2606.18487 (binary-reward advantage variance p(1−p)(g−1)/g).
   Fix: add **G-probe** to §5.1 for the selected checkpoint: 20 GRPO steps (G = 8, 32 prompts, HF generation in
   Phase 2), pass iff mean aggregate reward over steps 11–20 exceeds steps 1–5 (paired over prompts, CI > 0), legal
   rate at T = 1 does not fall by > 5 pts, and ≥ 50 % of groups keep σ above the floor at step 20. Pre-register it
   in `ab_prereg.py`; make D12(c)'s WiSE-FT repair fire on G-probe failure too.

3. **Statistics: seed pooling and primary-metric definitions are missing (§5, §5.2, §6).**
   Problem: every A/B and C2 uses "2 seeds per arm" and "paired CI over items", but the plan never says how two
   seeds per arm enter the paired bootstrap (per-item mean over seeds? seed as a cluster? one seed decides?). The
   topology probe's primary metric is also undefined across decodes ("K = 8 at T = 1 + greedy": greedy hit, any-of-8,
   or mean over 8 give different answers and different power).
   Fix: state once in §5: per-item statistic = mean over seeds and over the K samples (greedy reported separately);
   bootstrap resamples items (cluster = nodule where applicable) and, for 2 seeds, reports the CI both pooled and
   per seed, with the decision taken on the pooled CI and a "seed-discordant" flag when the two per-seed CIs
   disagree in sign. State that E1–E9 are decision rules without multiplicity correction and that C2 is the single
   confirmatory test.

4. **Baseline fairness on the shifted view (D17, §5.1 "near MedSAM").**
   Problem: the VLM is gated on *shifted* val; nnU-Net 2D is specified as "1 fold, 250 epochs" on frames with no
   view stated, and SAM 2.1/MedSAM rows are defined on centred frames with an oracle box. A supervised baseline
   trained and tested on centred frames inherits the exact leak S5 exposed (prior 0.904 / 0.783 / 0.588).
   Fix: specify that nnU-Net trains on the jittered train view (D10 shifts, scale jitter) and is evaluated on the
   same shifted val; add a shifted-val row for MedSAM/SAM 2.1 (box from the shifted frame's union, so they keep
   their oracle advantage but the view matches); label every baseline row with its view.

5. **GPU sequencing between trainer and vLLM sampler (§3.2, D13).**
   Problem: `sample_vllm` runs at `gpu_memory_utilization 0.85` "while the trainer process is stopped or idle". An
   idle trainer still holds model, optimizer and allocator cache on the GPU; vLLM at 0.85 will OOM or be refused.
   Fix: make `select_ckpt` a separate process that starts after `train_sft` has exited (branch checkpoints are
   model-only on disk), or have `train_sft` `del` the model/optimizer and `torch.cuda.empty_cache()` and assert
   `mem_get_info()` free ≥ 90 % before launching the sampler. Add this to the M8 CPU smoke (assert the sampler is
   launched only when `torch.cuda.is_initialized()` is False in the launching process).

6. **Evidence corrections.**
   (a) Quagmires (2510.01624) is cited backwards: A2 §4.10 and §6.1 say "one epoch on unique data beats more epochs on
   half the data"; the paper says "SFT training on unique examples for one epoch underperforms training on half
   examples for two epochs, either after SFT or SFT-then-RL". D5's "each epoch a fresh augmented view" and the
   augmentation rationale lean on this. The design survives (WSD branches at 0.5–3 epochs cover both regimes), but
   the plan must cite it correctly and drop "unique > repeated" as a reason. (b) §1 C2 and D20 cite "B3 P1" for
   "A0 scores exactly 0.50 [0.45, 0.55]"; that number is B3 P2's synthetic topology line (REPORT §3 P2, last bullet).
   (c) §4.9 unit costs come from `estimate_cost.py` at **2 epochs**, while D5 trains 3 stable epochs plus four 15 %
   decay branches (≈ 4 epoch-equivalents); recompute S3/S4 rows at 4 epochs or state the scaling. Also record B1
   AMEND 2 (150-step secondary budget) and AMEND 3 (Qwen3.5 bf16 base) in D1, which cites only AMEND 1.

## SHOULD-FIX

1. **Scope.** Promote only what C1 and C2 need to the critical path: S0, E3 (jitter), E1 (oversampling), screening,
   topology probe, baselines prior + MedSAM/SAM 2.1. Mark E2 (NTL), E4 (vision tower), E7 (determinism), E9,
   nnU-Net ×2, the K = 64 / I35 recalibration and HD95/ASSD/centroid as "after verdicts, budget permitting". The plan
   itself notes E3's outcome does not change the recipe ("jitter stays on" either way); consider dropping E3 to one
   seed as a sanity check rather than a 4-unit A/B.
2. **G-topology (b) ≥ 0.50 is vacuous** for any model that handles two-part items (A0 scores 0.50 by construction).
   Replace with "hole recall ≥ 0.30 on one-hole items" (the C2 fail line) or fold it into C2.
3. **Shifted val uses translation only**; training uses translation + scale jitter (D10). A model can still read size
   from S. Either add a scale-jittered val view (one fixed draw per item, S moved to a neighbouring divisor where
   the invariants allow) or state why translation-only suffices (G-prior partially covers it).
4. **Preload only what the process needs.** `train_sft` (LIDC-only screening) preloads *every* RL cache (~4.5 GB,
   mostly RIGA) and forks a reward pool it never uses during SFT. Preload per dataset and only in the process that
   scores (`select_ckpt` / GRPO probe). Keeps `cpu_ram ≥ 48` honest and shortens start-up.
5. **RIGA key delivery (D16, §4.4)** is unspecified ("key never on disk"). State: age identity piped over the direct
   SSH session's stdin into the decrypt command, never an env var, template field or file; ledger records the
   identity's public half only.
6. **Precision-probe patient disjointness (D19).** The capacity run trains on 20k items on train patches and tests
   on 1,000 held-out items; say the held-out items use patches from patients not used by the 20k (B3's probe uses
   val patches for topology and train patches for precision; keep the split explicit and tested in M4).
7. **Checker-filtered self-SFT as a cheap bridge.** The research does not consider rejection-sampling fine-tuning
   (sample K from the SFT checkpoint, keep legal + IoU_tc ≥ τ by the exact checker, one short SFT pass). It is the
   cheapest known repair for the legality knife-edge (A4 S4: legal ≥ 90 % needs σ ≤ 0.6 px) and a natural
   pre-GRPO step given a certified verifier; ROADMAP Phase 3 already has best-of-8 as a control. Pre-register it as
   an optional E10 triggered by "format passes, legal 85–90 % at every branch".
8. **D12 selection bias.** "Highest σ-above-floor share" prefers the least-trained branch that passes format; add the
   tie-break order explicitly as a lexicographic rule with minimum deltas (σ-share differences < 3 pts count as
   ties) so a 0.5-epoch branch does not win over a 1-epoch branch on noise.
9. **G-format reporting.** Add truncation share and mean/p95 completion tokens next to parse/legal (truncation is the
   runaway-cost failure A5 §7 flags), and held-out NLL per branch (D12 uses it as a tie-break but §5.1 never
   records it).
10. **Qwen3.5-0.8B sampler risk.** vLLM GDN kernel bug #34948 and the CPU-only parity gap mean the Qwen3.5 rung may
    have to use the HF sampler for every K = 8 branch eval; budget that (HF static batches with a 1,024-token
    runaway tail) in §4.9 or pre-decide "HF only, K = 8, 1,084 items" cost for that rung.
11. **Phase 3 env re-pin is already a Phase 2 deliverable (S7)** but §4.9 gives it 3 h with E8; keep S7 to the
    4-step colocate smoke only; E8 (coordinate tokens) is Phase 3 material unless S0 shows generation > 60 % of the
    GRPO step, as the plan says.
12. **Minor:** `experiments/phase2/lab/B1/common.py` `_shift` has 0 px clearance and no `other_nodules` check; the
    plan's M1 spec (≥ 2 px, no other nodule) supersedes it, so say the B1 shifted numbers were taken under the looser
    draw. `eval/gates.py` still encodes the topology-on-unanimous gate; M4 must replace, not add, it.

## Strengths

- The centring leak (A4 S5–S6, I43) was found before any GPU was rented, and the fix is complete across training
  (D10), evaluation (shifted val) and the null model (G-prior); the prior's collapse under shift (0.00 / 0.32 / 0.07,
  `null_prior_jitter.json`) is reproduced from code I read.
- The precision and topology probes are exact, brute-force-checked and calibrated (B3 P2; `probe_rule.json`), and
  the plan correctly replaces the vacuous "≤ 2 px" and "unanimous items" gates (0 unanimous val items,
  `census.json`).
- Every environment claim is tested, not asserted: trl 1.16.0 (`vllm<=0.31.0`) and vllm 0.32.0 (`torch==2.13.0`,
  `transformers<5.20`, `numba==0.65.0`) match PyPI; the two-venv image, prompt-id parity (32/32 × 3 models) and
  the numpy/numba override re-pin with 26/26 strict hashes are in B5.
- Reward serving is evidence-driven and minimal: preload + fork-after-preload pool, `CacheRef` (211 KB vs 96 MB),
  `BrokenProcessPool` semantics, C++ scorer deferred behind a measurable trigger (B4 §4–5). This is the right C++
  vs Python call.
- Checkpoint selection by GRPO signal with WSD branches follows the best current evidence (2606.18487 rank
  inversion, ICML 2026) and costs one run.
- Governance (RIGA permission AND Secure Cloud, bundle refusal tested by H-GOV, LIDC-only screening) and cost
  guardrails G1–G10 are concrete, cheap and testable.
- Unconstrained gate sampling with the constrained arm as diagnostic only (A5 §5–6, measured −18 to −111 nats
  IS mismatch) avoids a subtle RL bug.
- Five GPUs are chosen from a dated snapshot with roles, ceilings and fallbacks; money is correctly identified as
  non-binding and S0 replaces every estimate.

## What would make me APPROVE

MUST-FIX 1–6 applied in the plan text (and the C2 wording in ROADMAP), plus SHOULD-FIX 1 (a critical-path list
with the deferred experiments marked conditional) and SHOULD-FIX 2 (G-topology (b) replaced). With those, the
gates measure the two deliverables directly, the claim is scoped to what the probe can show, the statistics are
decidable before the first run, and the plan is roughly a third smaller on the critical path without losing a
verdict.
