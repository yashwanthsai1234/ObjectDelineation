"""B2 pre-registration: fixed on 2026-10-20 BEFORE any B2 measurement (no B1 sample has been re-scored and no
interpolated checkpoint has been sampled). signal.py / analyze.py import these constants; a criterion is never changed
after a run (a later re-specification gets a new name and an AMEND note here and in the report).

Question (PLAN D12): how should Phase 2 pick the SFT checkpoint handed to GRPO? Does deeper SFT reduce the GRPO
learning signal (H-A2-3, Mac form), and does scaling the LoRA delta toward the base recover it without losing
legality (H-A2-7 WiSE-FT, Mac form: base <-> 150-step interpolation, the only one B1's kept adapters allow)?

Inputs (no new training; owner cost directive):
- Zero-training: B1's per-sample outputs (<scratch>/phase2/B1/runs/*/samples.json), 64 LIDC val items x K = 8 at
  T = 1, centred and shifted views, LoRA r16 all-linear, lr 5e-4. Depth pairs (seed 0): 36 vs 150 steps for
  LFM2.5-VL-450M, SmolVLM-256M and Qwen3.5-0.8B (bf16 base). Seed noise: seed 1 vs seed 0 at 36 steps for LFM and
  SmolVLM.
- Interpolation (MPS, <= 2 evals, <= 7 GB): LFM 150-step adapter (vision intact, B1 §2.5) with every lora_B scaled by
  alpha in {0.5, 0.75} (W = W0 + alpha * (B A) * s, then merged), same 64 items x K = 8 at T = 1, both views, B1's
  common.sample_k path and torch.manual_seed(0) per view. Reference arms: alpha = 1 (B1's 150-step run) and the
  36-step seed-0 run. If one eval would exceed ~25 min, cut to the first 32 items (declared).

Scoring: every sample is re-scored with the frozen Phase 1 reward through objdel.eval.metrics (score_objects,
legality gate, advantages_v4, REWARD_V4 floors) on the item's own view (B1's common.views). Per-sample scalar reward R
= the raw composite of the v4 weights on the legality-gated terms, R = 0.3 IoU_tc + 0.3 F5 + 0.1 legal + 0.3 g_area
(PLAN §5.1 G-probe's score; illegal / unparseable samples score 0).

Metrics per GRPO group (one group of G = 8 per item):
  M1 r_std      population std of R over the group (the GRPO advantage scale before normalisation);
  M2 sigma      some legality-gated geometric term (IoU_tc, F5, g_area) has std > its REWARD_V4 floor (0.05): the
                K12 smoke / RESULTS §8 G-signal statistic (objdel.eval.metrics groups.signal);
  M3 mixed      >= 1 legal and >= 1 illegal sample in the group;
  M4 headroom   max R - mean R (best-of-8 minus mean);
  M5 distinct   distinct completion texts / 8 (and raster-distinct = 1 - duplicate share, reported);
  M6 legal_wrong share of samples that are legal with IoU_tc < 0.5 (the reward's "legal but wrong" mass);
  M7 active     share of groups whose v4 advantages are not all 0 (DAPO-style drop not triggered); mean |A| reported.
  Context: legal, mean R, mean IoU_tc.

Statistics: 95 % percentile bootstrap, B = 2,000, resampling LIDC NODULES (PLAN §5.0; item id up to the slice), the
same clusters for both arms of a paired difference. "Significant" = the paired CI excludes 0. Seed noise for a metric
= |seed 1 - seed 0| at 36 steps (LFM; SmolVLM reported), point estimate.

H-B2-1 (H-A2-3, Mac form): deeper SFT lowers the GRPO signal. Primary: SHIFTED view, M1 (r_std) and M4 (headroom),
  paired 150 - 36 per arm. ACCEPT iff for >= 2 of 3 arms both M1 and M4 fall significantly AND each fall exceeds
  LFM's seed noise on that metric. REJECT iff no arm has a significant fall in M1 or M4 on the shifted view and
  >= 2 arms have M1 or M4 significantly higher at 150 steps. Else INCONCLUSIVE. The centred view is reported.
  Sub-claim H-B2-1b (exploration): M5 distinct falls by >= 5 pts (paired CI below 0) for >= 2 arms -> "deeper SFT
  narrows exploration" accepted; else not shown.

H-B2-2 (H-A2-7, Mac form: interpolation toward the base recovers signal without losing legality). For each
  alpha in {0.5, 0.75}, paired against alpha = 1 (LFM 150):
  (i)   legal does not fall by > 1 pt on either view (point estimate >= alpha=1's - 0.01: A2's "legal >= late - 1 pt");
  (ii)  M1 or M4 significantly higher than alpha = 1 on the shifted view;
  (iii) M2 sigma share >= the 36-step run's - 0.03 on the shifted view (A2's "sigma share >= early's").
  ACCEPT (interpolation becomes a recommended D12 repair) iff some alpha meets (i)-(iii). REJECT iff every alpha
  loses legality significantly on some view (paired legal CI upper bound < 0) or every alpha's M1 and M4 are
  significantly LOWER than alpha = 1. Else INCONCLUSIVE.
  Secondary: an alpha that beats the 36-step run on legal (significantly) at no significant loss in M1 / M4 ->
  "interpolation dominates the shallower branch" (would favour interpolation over the next-earlier-branch repair).

D12 rule chosen by the results (written into the report's D12 text):
  R-A (signal-first, PLAN's lexicographic rule) is kept iff H-B2-1 is ACCEPTED.
  R-B ("deepest eligible, signal veto") otherwise: among branches that pass G-format and G-signal at T = 1 on the
      shifted val, take the deepest; replace it by the next shallower eligible branch only if that branch's paired
      r_std (M1) AND headroom (M4) are both significantly higher (nodule-cluster CI excludes 0) AND each by
      >= VETO_REL (15 %) of the deeper branch's value. Ties never move to the shallower branch.
  Repair when the selected branch fails G-signal or G-probe: interpolation (alpha = B2's best) iff H-B2-2 ACCEPTED;
      else the next-earlier branch (PLAN D12(c) default).

Also stated before measuring: SIGMA (M2) is expected to be ~saturated (B1: 97-100 % of groups at every depth), so it
cannot rank branches at Mac depths; M1 / M4 are continuous and are the primary statistics.
"""

WEIGHTS = (0.3, 0.3, 0.1, 0.3)  # R = sum of the v4 weights on the gated IoU_tc, F5, legal, g_area (REWARD_V4.weights)
IOU_OK = 0.5  # legal_wrong: legal and IoU_tc < this (REWARD_V4.iou_gate)
G, BOOT, CI = 8, 2000, 0.95
ARMS = ("lfm2.5-vl-450m", "smolvlm-256m", "qwen3.5-0.8b")
DEPTH_RUNS = {a: (f"{a}_s0", f"{a}_s0_n150") for a in ARMS}  # (36 steps, 150 steps), seed 0
SEED_RUNS = {a: (f"{a}_s0", f"{a}_s1") for a in ("lfm2.5-vl-450m", "smolvlm-256m")}  # seed noise at 36 steps
PRIMARY_VIEW = "shifted"
PRIMARY = ("r_std", "headroom")
METRICS = (
    "r_std",
    "sigma",
    "mixed",
    "headroom",
    "distinct",
    "raster_distinct",
    "legal_wrong",
    "active",
    "abs_adv",
    "legal",
    "r_mean",
    "iou_mean",
)
H1_ARMS_NEEDED = 2
H1B_DISTINCT_DROP = 0.05
ALPHAS = (0.5, 0.75)
INTERP_RUN = "lfm2.5-vl-450m_s0_n150"
INTERP_LEGAL_TOL = 0.01
INTERP_SIGMA_TOL = 0.03
VETO_REL = 0.15
MAX_EVAL_MIN = 25  # cut to 32 items if one alpha eval would exceed this

# AMEND 1 (2026-10-20, during the alpha = 0.75 eval, before any alpha = 0.5 sample; applies MAX_EVAL_MIN as written):
# alpha = 0.75 took 10.5 min for the centred view alone (9.9 s/item vs 3.8 s/item at alpha = 1: longer answers), so
# the full eval exceeds ~25 min. alpha = 0.75 keeps its 64 items; alpha = 0.5 runs on 32 items (eval_records("lidc",
# 32) is the first 32 of the same seeded permutation, a subset of the 64), and its paired comparisons use those 32.
ALPHA_ITEMS = {0.75: 64, 0.5: 32}
