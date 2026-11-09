# B2: SFT depth, GRPO signal and interpolation

Wave B task B2 of `experiments/phase2/RESEARCH.md`. Mac M3, 16 GB, MPS. Run 2026-10-20 in about 50 min (CPU analysis plus 2 MPS evals). Code: `experiments/phase2/lab/B2/`. Outputs: `<scratch>/phase2/B2/` (`scored/`, `runs/lfm_a0.75`, `runs/lfm_a0.5`, `analysis.json`, `posthoc.json`, `logs/`).

## 0. Hypotheses

- **H-B2-1 (A2's H-A2-3, Mac form).** Deeper SFT lowers the GRPO learning signal: within-group reward spread, headroom and exploration.
- **H-B2-2 (A2's H-A2-7, WiSE-FT, Mac form).** Scaling the LoRA delta toward the base (α < 1) recovers that signal without losing legality.
  - Label note: the brief calls interpolation H-A2-8. In `A2_sft_stack.md` and PLAN D12, H-A2-7 is WiSE-FT and H-A2-8 is DFT. DFT needs training and is not tested here.
  - This is a base↔late interpolation, the only one B1's kept adapters allow (no early LFM adapter was saved). A2's WiSE-FT interpolates an early and a late SFT checkpoint.
- **H-A1-8 (context).** Legality at T = 1 for each interpolated checkpoint.

## 1. Pre-registered design and criteria

Fixed in `prereg.py` before any B2 measurement. AMEND 1 applied the pre-registered 25-minute cut.

**Data**
- No new training. B1's stored K = 8 samples at T = 1 on 64 LIDC val items, centred and shifted views.
- Depth pairs (seed 0, 36 vs 150 steps): LFM2.5-VL-450M, SmolVLM-256M, Qwen3.5-0.8B (bf16 base).
- Seed noise: seed 1 vs seed 0 at 36 steps (LFM; SmolVLM reported).
- Interpolation: the LFM 150-step adapter (vision intact) with every `lora_B` scaled by α ∈ {0.75, 0.5}, then merged, then B1's exact eval path (`train_eval.evaluate` → `common.sample_k`, `torch.manual_seed(0)` per view).
- Path check before the evals: α = 1 reproduced B1's stored samples exactly (32 of 32 identical on 4 items).

**Scoring**
- Every sample goes through `objdel.eval.metrics`, i.e. the frozen reward v4: legality gate, `advantages_v4`, REWARD_V4 floors.
- Per-sample R is the G-probe's raw composite: R = 0.3 IoU_tc + 0.3 F5 + 0.1 legal + 0.3 g_area. Illegal or unparseable samples score 0.

**Metrics per group of G = 8**

| id | name | definition |
|---|---|---|
| M1 | r_std | std of R |
| M2 | σ-share | some gated geometric term has std > 0.05: the K12 / RESULTS §8 floor, `REWARD_V4.floors` |
| M3 | mixed | ≥ 1 legal and ≥ 1 illegal sample |
| M4 | headroom | max R − mean R |
| M5 | distinct | distinct texts / 8; raster-distinct = 1 − duplicate share |
| M6 | legal_wrong | legal with IoU_tc < 0.5 |
| M7 | active | v4 advantages not all 0; mean \|A\| also reported |

**Statistics**
- Nodule-cluster percentile bootstrap (PLAN §5.0), B = 2,000, 52 nodules across the 64 items.
- Paired differences resample the same nodules. "Significant" means the CI excludes 0.

**Rules**
- **H-B2-1:** primary statistics are M1 and M4 on the shifted view, 150 − 36.
  - ACCEPT iff, for ≥ 2 of 3 arms, both fall significantly **and** each fall exceeds LFM's seed noise.
  - REJECT iff no significant fall and ≥ 2 arms rise significantly.
  - Otherwise INCONCLUSIVE.
  - Sub-claim H-B2-1b: distinct falls by ≥ 5 points in ≥ 2 arms.
- **H-B2-2:** ACCEPT iff some α meets all three conditions:
  - (i) legal ≥ α=1's − 1 pt on both views;
  - (ii) M1 or M4 significantly above α = 1 on the shifted view;
  - (iii) σ-share ≥ the 36-step run's − 3 pts.
  - REJECT iff every α loses legality significantly on some view, or has M1 and M4 both significantly lower.
- **D12:** rule R-A (PLAN's signal-first lexicographic rule) is kept iff H-B2-1 is ACCEPTED; otherwise R-B (deepest eligible branch, with a signal veto at ≥ 15 % relative). The repair is interpolation iff H-B2-2 is ACCEPTED, otherwise the next-earlier branch.
- Stated in advance: σ-share is expected to be saturated (B1: 97–100 %), so M1 and M4 are primary.

**Amendments and post-hoc work**
- **AMEND 1:** α = 0.75 took 23.6 min for 64 items (9.9–12.2 s/item, against 3.8 s/item at α = 1, because answers get longer). α = 0.5 therefore ran on 32 items, a subset of the 64 (same seeded permutation); it took 19.6 min.
- **Post-hoc, exploratory (`posthoc.py`, no verdict attached):** because illegal R = 0, the group variance splits exactly into two parts:
  - Var(R) = p(1−p)·μ_L² + p·Var_L(R)
  - The first term is the legality part, the second the geometric part (spread among legal samples). sd_geom = √(p·Var_L).

## 2. Results

### 2.1 Depth (zero-training), shifted view (primary)

Point values are 36 → 150 steps. Δ is the paired difference 150 − 36 with its 95 % nodule-cluster CI.

| arm | legal | r_std (Δ) | headroom (Δ) | mixed (Δ) | legal_wrong (Δ) | mean \|A\| (Δ) | σ-share |
|---|---|---|---|---|---|---|---|
| LFM | 0.555 → 0.709, +0.154 [0.099, 0.213] | 0.193 → 0.164, **−0.028 [−0.044, −0.014]** | 0.262 → 0.215, **−0.047 [−0.070, −0.025]** | 1.00 → 0.86, −0.141 [−0.238, −0.061] | +0.144 [0.084, 0.204] | −0.055 [−0.104, −0.001] | 1.00 / 1.00 |
| SmolVLM | 0.459 → 0.668, +0.209 [0.159, 0.259] | 0.184 → 0.180, −0.004 [−0.014, +0.004] | 0.289 → 0.230, **−0.059 [−0.081, −0.037]** | −0.016 [−0.051, 0] | +0.203 [0.152, 0.254] | −0.026 [−0.059, +0.007] | 1.00 / 1.00 |
| Qwen3.5 | 0.822 → 0.834, +0.012 [−0.030, +0.056] | 0.143 → 0.148, +0.005 [−0.011, +0.019] | 0.171 → 0.163, −0.008 [−0.028, +0.012] | +0.031 [−0.102, +0.155] | +0.004 [−0.040, +0.049] | +0.050 [−0.002, +0.102] | 0.98 / 1.00 |

- **Seed noise** (LFM, s1 − s0, shifted): r_std −0.008 [−0.018, +0.002]; headroom +0.017 [−0.007, +0.041]; legal −0.092. LFM's depth falls (−0.028, −0.047) exceed it.
- **Saturated in every run, both views, all depths:** distinct texts 1.000, active groups 1.000, raster-distinct ≥ 0.98, duplicates ≤ 1.4 %.

**Centred view, paired 150 − 36**

| arm | r_std | headroom | mean \|A\| | other |
|---|---|---|---|---|
| LFM | −0.011 [−0.033, +0.011] | −0.117 [−0.155, −0.074] | **+0.424 [0.327, 0.515]** | IoU_tc mean +0.227 |
| SmolVLM | +0.007 [−0.006, +0.019] | −0.058 [−0.093, −0.023] | +0.228 [0.140, 0.316] | — |
| Qwen3.5 | **+0.023 [0.001, 0.043]** | +0.009 [−0.019, +0.037] | +0.368 [0.273, 0.451] | legal_wrong −0.166 [−0.218, −0.106] |

**Variance decomposition (post-hoc).** Paired Δ of sd_geom, 150 − 36:

| arm | centred | shifted |
|---|---|---|
| LFM | **+0.023 [0.009, 0.036]** | +0.005 [−0.004, +0.014] |
| SmolVLM | **+0.035 [0.023, 0.048]** | **+0.021 [0.013, 0.028]** |
| Qwen3.5 | +0.013 [−0.001, +0.029] (var_geom +0.0037 [0.0002, 0.0075]) | +0.001, not significant |

- LFM's shifted legality part, by contrast, falls significantly: −0.009 [−0.015, −0.003].
- Geometric share of the group variance, 36 → 150 steps:
  - LFM: 0.11 → 0.28 shifted, 0.17 → 0.35 centred.
  - SmolVLM: 0.08 → 0.17 shifted, 0.11 → 0.21 centred.
  - Qwen3.5: about 0.35–0.45 at both depths.

**Reading:**
- The shifted-view fall in r_std and headroom is legality saturation: the Bernoulli p(1−p) term that 2606.18487 identifies.
- It is not a loss of geometric exploration. The spread among legal samples, the part GRPO uses to learn geometry, grows with depth, and \|A\| on centred crops rises 0.23–0.42.
- Headroom = max − mean falls mechanically whenever the mean rises.

### 2.2 Interpolation toward the base (LFM 150-step adapter)

Paired Δ vs α = 1 unless noted.

| | α = 1 (B1) | α = 0.75 (64 items) | α = 0.5 (32 items) | 36 steps |
|---|---|---|---|---|
| parse, shifted / centred | 0.986 / 0.994 | 0.967 / 0.979 | **0.797 / 0.918** | 0.940 / 0.951 |
| legal, shifted | 0.709 | 0.623, **−0.086 [−0.150, −0.029]** | 0.375, **−0.324 [−0.415, −0.234]** | 0.555 |
| legal, centred | 0.832 | 0.762, **−0.070 [−0.109, −0.030]** | 0.598, **−0.234 [−0.315, −0.148]** | 0.588 |
| IoU_tc mean, centred | 0.414 | 0.356, −0.058 [−0.083, −0.031] | −0.209 [−0.248, −0.170] | 0.187 |
| r_std, shifted | 0.164 | +0.005 [−0.007, +0.019] | +0.000 [−0.028, +0.025] | 0.193 |
| headroom, shifted | 0.215 | +0.021 [−0.000, +0.044] | +0.095 [0.059, 0.132] (mean R −0.130) | 0.262 |
| sd_geom, shifted (post-hoc) | — | −0.001 [−0.011, +0.008] | **−0.016 [−0.030, −0.003]** | — |
| duplicates, shifted | 0.000 | 0.004 | **0.105** (raster-distinct −0.105 [−0.152, −0.062]) | 0.010 |
| σ-share, shifted | 1.00 | 0.98 | 1.00 | 1.00 |
| tokens mean, shifted; s/item | 91; 3.8 | 108; 12.2 | 141; 19.2 (truncation 1.6–2.0 %) | 80 |

- **α = 0.75 vs the 36-step run:** legal shifted +0.068 [0.002, 0.139], centred +0.174 [0.098, 0.250]; r_std shifted −0.023 [−0.037, −0.010]; centred IoU +0.169.
  - α = 0.75 behaves like an intermediate-depth checkpoint, between 36 and 150 steps on every metric.
  - It does not dominate the shallower branch.
- **α = 0.5** is below the 36-step run on shifted legality (−0.191 [−0.270, −0.117]) and fails G-signal's duplicate limit (10.5 % > 10 %).
- The base cannot write the D10 format (B1: 0 % parse zero-shot), so moving toward it un-learns the format first. Diversity does not come back; duplicates rise.
- Peak MPS memory 3.6 GB.

## 3. Verdicts

| hypothesis | verdict | basis |
|---|---|---|
| **H-B2-1** (deeper SFT lowers GRPO signal) | **INCONCLUSIVE** (pre-registered rule) | LFM: both primaries fall beyond seed noise. SmolVLM: headroom only. Qwen3.5: neither (its legality had already plateaued at 36 steps). Post-hoc: every fall is in the legality part; the geometric part is flat or rises. At Mac depths, deeper SFT is not over-training. |
| H-B2-1b (exploration narrows) | **not shown** | distinct texts 100 % and duplicates ≤ 1.4 % at every depth |
| **H-B2-2** (base interpolation recovers signal, keeps legality) | **REJECT** | both α lose legality significantly on both views; α = 0.75 raises no primary significantly; α = 0.5 lowers geometric spread and doubles answer length |
| H-A1-8 at α < 1 | fail | legal 0.62 / 0.38 shifted (< 80 %); α = 0.5 duplicates 10.5 % |
| PLAN D12(a) / (b) (does signal predict GRPO gain?) | **not testable on the Mac** | no GRPO run; see §5 |

## 4. D12 rule text for PLAN.md

> **D12 Checkpoint selection (B2 filled; rule R-B).**
>
> **Eligible:** WSD branches (0.5 / 1 / 2 epochs) that pass G-format (parse ≥ 98 %, legal ≥ 90 %) and G-signal (duplicates ≤ 10 %, σ-share ≥ 50 %) at T = 1 on the shifted val, K = 8.
>
> **Select the deepest eligible branch.**
>
> **Signal veto:** replace it with the next shallower eligible branch only if both hold:
> - that branch's **geometric reward spread** sd_geom = √(p·Var_legal(R)) on the shifted val (R = the G-probe raw composite, per item, then averaged) is higher, with the paired nodule-cluster CI excluding 0, **and** by ≥ 15 % relative;
> - its shifted-val median IoU_tc is not significantly lower.
>
> **Ties stay deep.**
>
> **Never select by:** greedy legal or IoU; σ-share, active-group share or distinct-text share (saturated at 98–100 % in B2); total R std or best-of-8 − mean. The last two fall mechanically as legality rises (B2: the entire 36→150 fall is the p(1−p) legality term).
>
> **Reported per branch, not used for selection:** r_std, headroom, mixed-legality share, legal-but-wrong share, mean |A_v4|, TRL entropy, held-out NLL, pass@8.
>
> **The selected branch then runs G-probe (§5.1).**
>
> **D12(c) repair** (the selected branch fails G-signal or G-probe): the next-earlier eligible branch.
> - **LoRA-delta / base interpolation is not a repair:** B2 H-B2-2 REJECTED (α = 0.75: −7 to −9 legal pts; α = 0.5: 38 % legal, 10.5 % duplicates).
> - WiSE-FT between two SFT branches (0.5-epoch ↔ selected, full FT) stays untested. When the repair fires, it may be evaluated alongside the next-earlier branch at inference cost only (merge + shifted val K = 8, same rule); otherwise it is not used.
>
> **D12(a) / (b)** (does the signal branch out-gain the greedy branch under GRPO?) are decided in R2 by the folded measurement below, not on the Mac.

Note on the veto statistic: sd_geom came out of B2's post-hoc decomposition. B2's pre-registered veto statistics were r_std and headroom, and B2 shows those penalise legality gains. Inside the eligible set (legal ≥ 90 %), the legality part is capped at 0.09·μ_L², so sd_geom and r_std mostly agree there anyway. If the plan prefers only pre-registered statistics, use r_std and keep the rule otherwise unchanged.

## 5. What folds into the Vast runs (no extra A/B; PLAN §6.2 additions)

| Measurement | Folded into | How | Cost |
|---|---|---|---|
| B2 signal panel per branch: r_std, headroom, mixed, legal_wrong, \|A_v4\|, sd_geom / legality-part split, σ-share, duplicates | R1 `select_ckpt` | CPU scoring of the shifted-val K = 8 samples it already draws per branch (promote `lab/B2/score.py` + `posthoc.parts` into `select_ckpt`) | 0 GPU; seconds of CPU per branch |
| Depth curve with 3 points | R1 | the 0.5 / 1 / 2-epoch decay branches already saved; the stable-phase checkpoints at those points are kept too, at no compute cost | 0 |
| Does the signal predict GRPO gain (H-A2-3 mechanism) | R2 G-probe | per item (128 items): Spearman of the step-0 sd_geom (from the step-0 double sample) and the mixed flag against the per-item G-probe gain R₂₀ − R₀, with a bootstrap CI | 0 (uses G-probe's own samples) |
| Shallow vs deep under GRPO (D12(a)/(b)) | R2, only when the plan already probes a second branch (selected not PASS) | report both branches' paired gains side by side | 0 extra |
| Entropy proxy | R1 SFT logs | TRL `entropy` and held-out NLL at each branch (already logged) | 0 |
| WiSE-FT (0.5-ep ↔ selected) | R1, only if the D12(c) repair fires | merge at α = 0.5, shifted val K = 8 | about 0.05 h inference |

## 6. What the Mac cannot answer

- **Whether any signal metric predicts the GRPO gain.** B2 ran no GRPO, so H-A2-3's decisive half (20-step gain) is open.
- **The over-trained regime.** Mac checkpoints are under-trained: legal 46–85 %, LoRA r16, ≤ 150 steps. The rank inversion of 2606.18487 appears only once SFT has compressed the distribution, and the WSD 2-epoch full-FT branches may reach that.
- **Image-conditioned geometric signal.** On shifted crops no B1 checkpoint uses the image (median IoU_tc 0), so shifted-view geometric spread here is image-blind.
- **WiSE-FT between two SFT checkpoints.** No early LFM adapter was kept, and the Vast recipe is full FT.
- **Seed variance at 150 steps**, and Qwen3.5 / SmolVLM with a frozen vision tower at depth.

## 7. Limits

- B1 checkpoints only: LoRA all-linear at lr 5e-4. Qwen3.5 and SmolVLM have collapsed vision (B1 §2.5); LFM's is intact. Single seed at 150 steps.
- α = 0.5 ran on 32 items (AMEND 1). Its CIs are wider and its 36-step comparison uses those 32 items.
- Base↔late interpolation is a harsher test than A2's early↔late WiSE-FT, because the base cannot write the format. The REJECT applies to the former only.
- The variance decomposition and the sd_geom veto statistic are post-hoc.
- σ-share and distinct share were saturated, so any branch differences in them are below this test's resolution.
- Multiplicity: 3 arms × 2 views × 12 metrics, with no correction. Only the shifted-view M1/M4 rules carry verdicts.
- The nodule-cluster bootstrap gives slightly wider CIs than B1's item bootstrap.

## 8. Reproduce

Run from the repo root. Scripts read untrusted `data/`, hence `-I`.

```bash
HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B2.interp --alpha 1 --check 4          # path check: 32/32 identical to B1
HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B2.interp --alpha 0.75                 # 23.6 min, MPS 3.6 GB
HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B2.interp --alpha 0.5 --eval-items 32  # 19.6 min (AMEND 1)
.venv/bin/python -I -m experiments.phase2.lab.B2.score                                   # B1 depth/seed runs
.venv/bin/python -I -m experiments.phase2.lab.B2.score --runs lfm_a0.75 lfm_a0.5
.venv/bin/python -I -m experiments.phase2.lab.B2.analyze      # -> <scratch>/phase2/B2/analysis.json (verdicts)
.venv/bin/python -I -m experiments.phase2.lab.B2.posthoc      # -> <scratch>/phase2/B2/posthoc.json (exploratory)
```

## 9. Sources

- Rank inversion: pre-RL pass@1 rises with SFT depth while peak GRPO pass@10 falls from 0.806 to 0.481; pre-RL entropy ρ = +0.69; binary-reward advantage variance p(1−p)(g−1)/g. Abstract re-read 2026-10-20: https://arxiv.org/abs/2606.18487
- WiSE-FT for reasoning ("interpolating the weights of the latest SFT checkpoint with an early checkpoint … almost completely recovers Pass@k while also improving Pass@1"), re-read 2026-10-20: https://arxiv.org/abs/2504.10478 ; original WiSE-FT: https://arxiv.org/abs/2109.01903
- Quagmires in SFT→RL (more SFT epochs, smaller post-RL gains; held-out loss and pass@k as proxies): https://arxiv.org/abs/2510.01624
- DAPO (equal-reward groups give zero advantage; dynamic sampling): https://arxiv.org/abs/2503.14476
- Entropy mechanism of RL for reasoning models: https://arxiv.org/abs/2505.22617
- PEFT LoRA (scaling, `merge_and_unload`): https://huggingface.co/docs/peft/main/en/package_reference/lora
- TRL GRPOTrainer HF generation path (B1's `sample_k`): https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/grpo_trainer.py
- Repo:
  - `src/objdel/reward/spec.py` (REWARD_V4 weights and σ floors 0.05 / 0.05 / 0.25 / 0.05);
  - `src/objdel/eval/metrics.py` (groups.signal, duplicate hash);
  - `experiments/phase1/RESULTS.md` K12 (σ above floor in 100 % of groups);
  - `experiments/phase2/lab/B1/REPORT.md`;
  - `docs/research/notes/Phase 2 SFT/A2_sft_stack.md` §4.10, H-A2-3 / 7 / 8;
  - `A1_models_tokens.md` H-A1-8;
  - PLAN §5.0, §5.1, §6.2, D12.
