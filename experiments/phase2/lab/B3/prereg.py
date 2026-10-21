"""B3 pre-registration: pass criteria fixed on 2026-10-19 BEFORE any B3 test was run. Scripts import these constants
and print PASS / FAIL against them; a criterion is never changed after a run (a later re-specification gets a new
name and a RE-REG note in the report).

P1. Topology probe (H-A1-7). Frame resolution, the reward's size rule min_area(|union|), voting raters as in
    reward.score.rater_answers. Non-simple item: >= 1 voting rater's size-filtered signature != (0,). Unanimous
    non-simple: every voting rater has the same non-simple signature. Leakage unit: the LIDC patient (group_id).
    1. 95 % CI half-width of the primary metric <= 6 points at true rate 0.5, one model, by a CLUSTER bootstrap over
       nodules (adjacent slices are correlated), not the binomial formula; also reported at 0.3 / 0.7 and for the
       paired A0 vs A0' difference.
    2. 0 probe patients among SFT-train patients.
    3. The probe leaves >= 50 % of the train split's non-simple targets for SFT (H-A1-6 oversamples them).
    4. Expected (A1): the 48-item val slice alone fails 1.
P2. Precision probe (A4 section 6, H-A4-7, H-A4-8).
    1. Exact on the ground-truth raster: ASSD = HD95 = 0, IoU = 1, NSD@1 = 1, centroid error 0, topology correct on
       100 % of >= 1,000 items; report the share of certified target texts that rasterize pixel-exactly.
    2. ASSD / HD95 equal a brute-force implementation to 1e-9 on >= 1,000 random pairs.
    3. A 1-px translation (4 directions): ASSD > 0 and centroid error = 1 px on 100 % of items; the pass rule fails it.
    4. Calibration (H-A4-7): centroid equal-area disc median IoU = 0.806 +- 0.03; vertex noise sigma 0.5 / 1 px gives
       median IoU within +- 0.05 and legality within +- 5 points of A4 S4 (0.90 / 95.7 %, 0.79 / 83.7 %).
    5. Discrimination (H-A4-8), A4 pass rule on 6-16 px objects (legal >= 95 %, median ASSD <= 0.6 px, median IoU >=
       0.80, median centroid error <= 0.5 px, topology >= 90 % on hole / two-part items): exact and sigma <= 0.5 pass;
       sigma >= 1, centroid disc, centred prior and 1-px shift fail.
P3. Offline dihedral targets (H-A6-6, H-A4-5): every train target (LIDC 23,469; RIGA 3,552) x 8 D4 ops rebuilt with
    objdel.targets.build_target from transformed masks.
    1. 100 % build and pass the text-only re-check (parse, legal, size-filtered signature at min_areas, IoU >= 0.99).
    2. Inverse-transformed variant raster: IoU >= 0.99 and same signature vs the original rater mask, 100 %; F5 >= 0.95
       and IoU_tc >= 0.99 against the original item's reward cache (D4 maps the lattice onto itself; EDT is
       D4-invariant).
    3. Identity variant per-file hashes = data/targets (LIDC train b4814506..., RIGA train 54bc88bf...).
    4. Two runs with different process counts give identical per-variant target_set_hash.
    5. <= 20 CPU-minutes for all 8 variants; wall <= 15 min at <= 2 processes.
    6. Tokens (SmolVLM-256M, LFM2.5-VL-450M, Qwen3-VL-2B, Qwen3.5-0.8B): no variant over MAX_COMPLETION unless its
       identity is; median relative length change within +- 5 %.
P4. Crop before trace (H-A6-9), sandbox monkeypatch of objdel.targets.build._rings (src/ untouched):
    1. every per-file and the overall target-set hash unchanged (1af819f4...);
    2. RIGA median per-target build time >= 15 % lower, paired per target, interleaved in one process.
"""

P1_HALF_WIDTH = 0.06
P1_MIN_TRAIN_NONSIMPLE_KEPT = 0.50

P2_EXACT_N = 1000
P2_BRUTE_TOL = 1e-9
P2_DISC_IOU = (0.806, 0.03)
P2_NOISE_REF = {0.5: (0.90, 0.957), 1.0: (0.79, 0.837)}  # sigma -> (median IoU, legal), A4 S4
P2_NOISE_TOL = (0.05, 0.05)
PASS_RULE = {"legal": 0.95, "assd": 0.6, "iou": 0.80, "centroid": 0.5, "topology": 0.90}
PASS_BINS = (6.0, 16.0)  # equivalent diameter, frame px
P2_SHOULD_PASS = ("exact", "sigma_0.25", "sigma_0.5")
P2_SHOULD_FAIL = ("sigma_1", "sigma_2", "centroid_disc", "centred_prior", "shift_1px")

P3_IOU = 0.99
P3_F5 = 0.95
P3_CPU_MIN = 20.0
P3_WALL_MIN = 15.0
P3_TOKEN_MEDIAN_CHANGE = 0.05
IDENTITY_HASH = {"lidc/train": "b481450688d3b8920562be3516ddd77401379288e4925f55c10954bd99c2c4c6",
                 "riga/train": "54bc88bfce0f7ee80d572ac88989442bb8556e3fcf6c5f0b00f2b579123a86ef"}

P4_HASH = "1af819f474536dffe5a4a33888ed343bebd3855383fc162b3ca00e0ab82e91d3"
P4_FILE_HASH = {**IDENTITY_HASH,
                "lidc/val": "138290299b5afcf72c10f45ff32171473ba8cc56c2b3b2c2853aae568d0e6c40",
                "riga/val": "39b2f17131b59c8ae9a3cbadefaa770d7de8bd83a3851a104f7a1e4a51015190"}
P4_RIGA_SPEEDUP = 0.15

# ---------------------------------------------------------------------------------------------------------------
# P5 (addendum, judge MUST-FIX 1 + 3; written 2026-10-19 BEFORE hole_census / topo_probe_v3 were run).
# P5.1 Empirical hole / part geometry of the real non-simple LIDC targets (all splits, voting raters, after the size
#      rule): hole count, hole area (px, fraction of the filled object), hole offset (centroid distance / object
#      equivalent radius), wall (px of object between hole and outside); parts: count, area ratio (2nd / 1st), gap
#      (px). Descriptive; the generator draws from the TRAIN records only (val/test reported for comparison).
# P5.2 Generator fidelity (generator-2 outlines, E = 0.2, A = 0.15; hole / part geometry resampled jointly from the
#      train records; minimum at the size rule, i.e. every hole >= min_area(|GT|)): generated vs train empirical
#      two-sample KS p >= 0.01 for hole area, hole fraction, hole offset, part ratio and gap; hole-area bin shares
#      (bins below) within +- 5 points of the train shares.
# P5.3 Probe design: hole-area bins (judge) <= 10, 11-30, > 30 px; equal allocation, 120 one-hole items per bin
#      (360), pooled hole recall = train-share-weighted mean of the bins (post-stratified), stratified bootstrap.
#      Pass: pooled 95 % CI half-width <= 6 points at true recall 0.5 under the recommended metric (P5.5);
#      per-bin half-widths reported (expected about +- 9: the per-bin floor "hole recall >= 0.30" is a
#      descriptive floor checked on point estimates, not a powered test); also report the floor's false-fail
#      probability at true bin recall 0.5 and its false-pass probability at true 0.15.
# P5.4 Sanity on the recalibrated probe: exact -> hole recall 1.0 and non-simple hit 1.0 in every bin, false
#      non-simple 0 on controls; A0 (holes emptied) -> hole recall exactly 0 in every bin; centroid disc -> 0;
#      sigma = 0.5 px vertex noise -> reported per bin (no criterion; expectation: recall rises with hole size).
#      Part metrics reported separately, descriptive only (A0 keeps parts).
# P5.5 Primary metric across decodes, chosen by a pre-stated rule: any-of-K is excluded (it rewards hedging, i.e.
#      emitting holes at random, and saturates); between greedy and mean-over-K = 8 at T = 1, the one with the
#      smaller simulated paired (A0' - A0) CI half-width at n = 360 is primary; if the two are within 10 % of each
#      other, mean-over-K (the quantity GRPO samples) is primary. Seeds: per-item statistic = mean over the 2 seeds
#      (and the K samples); bootstrap resamples items (nodules on real data), pooled CI decides, per-seed CIs and a
#      seed-discordant flag (per-seed CIs on opposite sides of 0) are reported.
# ---------------------------------------------------------------------------------------------------------------
P5_BINS = (10, 30)  # hole area px: <= 10 | 11-30 | > 30
P5_PER_BIN = 120
P5_KS_P = 0.01
P5_SHARE_TOL = 0.05
P5_HALF_WIDTH = 0.06
P5_FLOOR = 0.30
