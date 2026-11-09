"""B1 pre-registration: criteria fixed on 2026-10-19 BEFORE any B1 measurement (zero-shot or SFT eval). analyze.py
imports these constants; a criterion is never changed after a run (a later re-specification gets a new name and an
AMEND note here and in the report).

Common. Arms SmolVLM-256M-Instruct, LFM2.5-VL-450M, Qwen3.5-0.8B (non-thinking default). Eval items: 64 LIDC val
items, seeded random subset (seed 0) of the RL-eligible ones; test 1 also all 58 RL-eligible RIGA val items. Views
(coordinator, after A4 findings 1-3): centred (targets.frame) and shifted (same window moved by a uniform offset that
keeps the raters' union inside: A4 null_prior_jitter procedure, one draw per item). Prior reference row: A4's
image-blind centred disc (train median diameter per frame side S), scored with reward.score_mask on both views.
Statistics: 95 % percentile bootstrap over items (B = 2,000); paired comparisons resample the same items; "significant"
= the paired 95 % CI excludes 0.

T1  H-A1-1 zero-shot box prior. Greedy, <= 96 new tokens; native prompt per model card (Qwen3.5 bbox_2d 0-1000;
    LFM2.5-VL bbox in [0, 1]) and a generic {"bbox": [...]} 0-1000 prompt. Strict parse: key + 4 numbers in range,
    x1 < x2, y1 < y2. Box IoU vs the box of the raters' union (RIGA disc, LIDC nodule); unparsed = 0.
    ACCEPT (Qwen3.5 kept for its prior): shifted RIGA disc, native prompt: parse >= 90 % AND median box IoU >= 0.6
    AND paired mean (model - prior box) IoU CI lower bound > 0. REJECT: median box IoU < 0.4 on the shifted crops of
    both datasets. Else inconclusive. The centred version (A1's original rule) is reported, not decisive.
T2  H-A1-2 format learnability at equal budget. Same 1,000 LIDC train items (seed 0), every voting rater's target,
    configs.sft_config + configs.lora (r16, alpha 32, all linear), lr 5e-4 constant after 10 warmup steps, 8 rows per
    step, fp32, MPS <= 7.0 GB, N equal steps (N from a timing pilot so the slowest arm trains <= ~45 min). Eval:
    64 LIDC val items x 2 views, K = 8 at T = 1 (top_p 1, top_k 0), TRL GRPOTrainer HF generation, objdel.eval.
    Judged on the SHIFTED view (centred reported):
    - learns the format: parse >= 95 % AND legal >= 90 %;
    - uses the image: paired mean (item-mean IoU_tc - prior IoU_tc) CI lower bound > 0;
    - ranking: pooled median IoU_tc (illegal = 0); A > B iff the paired CI of the difference excludes 0, else tie;
    - H-A1-2 "LFM >= Qwen3.5 > SmolVLM": ACCEPT iff LFM not significantly below Qwen3.5 AND (Qwen3.5 significantly
      above SmolVLM OR SmolVLM fails "learns the format" while Qwen3.5 passes); REJECT iff Qwen3.5 significantly
      above LFM or SmolVLM significantly above Qwen3.5; else inconclusive;
    - A1's SmolVLM rule: drop SmolVLM from the Vast screen if its shifted median IoU_tc is > 0.05 below the best arm.
T3  H-A1-8 sampling legality at T = 1: per arm, legal >= 80 % AND duplicates <= 10 % on both views -> no constrained
    decoding needed at GRPO start; else adopt the pattern() grammar (D16 / A5). LFM zero-led number tokens ('042'):
    material failure mode iff >= 1 % of LFM samples fail with parse error leading_zero.

AMEND 1 (2026-10-19, after the timing pilots, before any eval): micro-batch 2 OOMs at the 6.5 GiB cap for every arm
without gradient checkpointing and for Qwen3.5 with it, so micro-batch 1 x 8 accumulation (still 8 rows per step);
Qwen3.5 runs at 79 s/step on MPS (reference PyTorch GDN / causal-conv1d kernels; fla is CUDA/Triton only), so
N = 36 steps for every arm (Qwen3.5 ~48 min); max_completion_length 1,024 (the longest LIDC train target is 881
tokens on the digit tokenizers, so 512 was not "above every target").
"""

H1_PARSE, H1_IOU, H1_FAIL_IOU = 0.90, 0.60, 0.40
FORMAT_PARSE, FORMAT_LEGAL = 0.95, 0.90
SMOL_GAP = 0.05
H8_LEGAL, H8_DUP, LEADING_ZERO_MATERIAL = 0.80, 0.10, 0.01
STEPS, ROWS_PER_STEP, LR, MAX_COMPLETION = 36, 8, 5e-4, 1024  # AMEND 1
EVAL_ITEMS, K, BOOT = 64, 8, 2000

# AMEND 2 (2026-10-19, after the 36-step SmolVLM and LFM evals, before Qwen3.5's): at 36 steps no arm met "learns the
# format" (legal 46-59 %), so the 36-step ranking can only speak to the earliest phase. Added, as an EXPLORATORY
# secondary budget with the same criteria: 150 steps for LFM2.5-VL-450M and SmolVLM-256M (Qwen3.5 at 150 steps would
# take ~3.3 h of MPS training, outside the budget), then seed 1 at 36 steps for LFM and SmolVLM (seed variance).
STEPS_SECONDARY = 150

# AMEND 3 (2026-10-19, before any Qwen3.5 eval): the fp32 Qwen3.5 run hit the 6.51 GiB MPS cap at step 13 (a long
# LIDC row). Qwen3.5 is re-run with its frozen base weights in bf16 (LoRA adapters fp32, PEFT's default upcast); the
# merged model samples in bf16. SmolVLM and LFM stay fp32. A precision asymmetry, declared as a limit.
QWEN_BASE_DTYPE = "bfloat16"

# AMEND 4 (2026-10-19, before any Qwen3.5 eval): with bf16 base weights Qwen3.5 trains at 34 s/step, so the
# exploratory 150-step tier includes Qwen3.5 too (~86 min of MPS training, over the ~45 min guideline; declared).
# Queue: LFM 150, SmolVLM 150, Qwen3.5 150, then seed 1 at 36 steps for LFM and SmolVLM.

# AMEND 5 (2026-10-20, POST-HOC and EXPLORATORY, motivated by results): the post-hoc diagnostics (image_use.py,
# vision_feats.py) show the configured adapter (LoRA on all linear layers, vision tower and projector included,
# lr 5e-4) collapses the visual features of Qwen3.5 and SmolVLM after SFT on centred LIDC crops. Added arms: the same
# 36-step recipe with LoRA on the language model only (vision tower + projector frozen; common.lora("lm")) for
# Qwen3.5 and SmolVLM. Read with image_use / vision_feats and the T2 metrics; no verdict is attached.
