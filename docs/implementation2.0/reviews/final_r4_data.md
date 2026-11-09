# Final review, round 4: data quality and statistics (Implementation.md)

Reviewer role: medical imaging data scientist and statistician (LIDC, multi-rater, leakage, cluster bootstrap).
Scope: the whole doc, with A and B1 first. I read the sources `phases/{00_front,A,B}.md` section by section, plus the
V1 parts of C. I checked every item against the review logs (`docs/implementation2.0/review-log.md`: A parts 1-4 and
the consistency pass, B.R to B.R4 and the consistency pass, front rounds 1-3). No item below repeats a fixed one.
Counts: **1 blocking, 4 major, 6 minor (11 items).**

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill covers split design, cluster bootstrap or LIDC label QA. The doc's `dicom-series-preflight` recommendation (B1.7 step 2) still fits item 8. Nothing installed or newly recommended. |
| engineering:code-review | Checked the plan's commands against the argparse code (`train_rl.py:24-27`, `final_eval.py:48-51`; item 5), the scheduler text against `augment.neighbour_side` (item 7), and test I22 against B1.5 (item 1). |
| engineering:debug | Reproduced before judging. Ran the B1.2 sampling rule as a simulation (item 3). Ran B1.9's power rule with its noise model and with a measured one (item 2). Ran B1.7 check (a) on the local XML (item 8). |
| engineering:system-design | Traced each V1 command (b1_rows.yaml steps 3-4) from the plan text to the box (items 1, 5). |
| engineering:architecture | Treated "what noise the gate sees" and "what share of disagreement rows B1 trains on" as single decisions with one owner each (items 2, 3). |
| engineering:documentation | Every fix names the source file, the heading, the text to change and the new text. |
| huggingface-skills:huggingface-best | No model choice in scope. Used only to confirm that selection stays on val and test is read once (it does). |
| huggingface-skills:huggingface-datasets | Re-read how B1.7 uses the third-party card `YongchengYAO/LIDC-IDRI-Lite`. The card is context; it is not a measurement on our store (item 8). |

Web checks (2026-11-08): TCIA LIDC-IDRI collection page, "Annotation and Markup Issues" items 1-8
[web: https://www.cancerimagingarchive.net/collection/lidc-idri/]. Items 4 (LIDC-IDRI-0101), 6 (0396) and 7 (0510)
and the 8 two-timepoint patients are already handled (B1.7, splits by patient). No new issue affects A or B1.

New measurements (read-only, `.venv/bin/python -I`, nothing written under `data/` or the repo). Scripts are in
`(planning probe, not kept: r4data)`:

- **Test-retest noise of the R5 model** [measured: `r4data/retest.py`]. `runs/rl/eval/final/heldout.json`
  (`train_rl/final` = `rl-full-3ep`, K 8, seed 11, rl rental) and `runs/r5/eval/eval_val/items.jsonl` (same model,
  K 8, seed 11, R5 box) share 128 `lidc-val-shifted` items. Per-item difference: mean -0.0076, SD 0.0569, SE 0.0050.
  Pearson 0.9712. Variance of the difference 0.00324, which is **7.7 %** of the two-arm bound
  mean(2 mu (1 - mu) / 8) = 0.0422.
- **B1.9 power, approximate** [estimate: `r4data/power2.py`, 300 runs, B = 400 patient draws, percentile CIs, R5 val
  items with geom.json offsets]. Effect: 25 % of far misses converted, chosen by patient. P(pass) for IoU, far-miss,
  near, DiD and joint:
  - noise on the B1 arm only, at the bound (B1.9 as written): 0.81 / 0.43 / 0.50 / 0.75 / **0.13**;
  - both arms at the bound: 0.61 / 0.42 / 0.30 / 0.57 / 0.03;
  - both arms at the measured level (c = 0.077): 0.997 / 0.92 / 1.00 / 0.99 / **0.92**;
  - measured level with no effect: IoU 0.03, DiD 0.027, far-miss 0, joint 0.
- **Disagreement share under B1.2** [estimate: `r4data/disagree.py`]. Weight 4, no window used twice in a run,
  7,200 shifted rows in 15 blocks: see item 3. Inputs: 289 / 7,341 train disagreement items, 44 / 1,084 val
  [measured: `data/prepared/lidc/train.jsonl`, `data/eval_views/lidc-val-shifted/lidc/val.jsonl`].
- **Val subsets** [measured: R5 val items joined with geom.json]. Far (> 0.3): 406 items, 75 patients, miss 27.3 %.
  Near (< 0.15): 206 items, 52 patients, IoU 0.751. Far misses average IoU_tc 0.012. Far found items average 0.688,
  all found items 0.706. Padded val shifted items: 5.
- **B1.7 check (a) on the local XML** [measured: `lidc.sop_z("data/raw/lidc_xml")` joined with the store]. Store
  series with a z value shared by different `imageSOP_UID`: 380786 (train, LIDC-IDRI-0013, 1 z), 228428 (train,
  -0017, 5 z), 074972 (test, -0659, 11 z). No val series. The pylidc `zvals` hold no duplicate z for any of the 800
  store series. The sidecars of LIDC-IDRI-0085 and -0146 have `slices == annotated_sops` (9/9, 15/15).

---

## Items

### 1. [blocking] [B] Test I22 asserts a `b1_rows.yaml` that B1.5 and front 6.2 forbid

- Problem. Two texts both claim to define `b1_rows.yaml`, and `make test` must pass before the commit.
  - B.md, B.9.3, row I22: "`b1_rows.yaml` has exactly the four steps of B1.5 in that order (`train_rl` with ...
    `--budget-usd 1.90`, ...; B1 on val shifted; B1 on val edge; **R5 on val edge**)".
  - B.md, B1.5: "`--budget-usd 2.05` ... is the front 6.4 rule 5 value ($2.40 - $0.35)", and "R5 on val edge is
    not a `b1_rows.yaml` step: it is front item 4, in A3's `a3_measure.yaml`". The four B1.5 steps are
    `env_check`, the GPU smoke, `train_rl`, one `final_eval`.
  - 00_front.md, 6.2 "Plan files": `(3) train_rl, item 2, with --out {rd}/b1-rows --budget-usd 2.05`; 6.4 rule 5:
    "V1: $2.40 - $0.35 = $2.05".
  - An engineer who writes the plan from B1.5 fails I22. One who writes it to pass I22 runs R5 on val edge twice
    (6.9 extra min) and sets the stop at $1.90. At $1.90 the stop trips at session minute (1.90 - 0.10) / 0.01133
    = 158.9, that is (158.9 - 30) x 60 / 27.5 = 281 steps at the base step time [estimate: front 6.2 inputs]. 281
    is not a multiple of 20, so the run ends inside a block and the B1.2 "every block prefix" law does not hold for
    the last block.
- Fix. B.md B.9.3, I22, replace the `b1_rows.yaml` clause with: "`b1_rows.yaml` has exactly the four steps of B1.5
  in this order: `env_check.sh`; `pytest -m gpu tests/test_train_fast.py -k "G02 or G08"`; `train_rl` with
  `--model qwen3-vl-2b`, `--ckpt {rd}/ckpt`, `--numerics bf16-seqmask`, beta 1e-3, `--lr 5e-5`, `--rows-revision`,
  `--rows-manifest-sha256`, `--budget-usd 2.05`, steps a multiple of 20, `--save-steps 50`,
  `mask_truncated_completions` true; one `final_eval` with `--ckpt {rd}/b1-rows/final`, `--train
  {data}/prepared` and the views `lidc-val-shifted` then `lidc-val-edge`. No step runs R5, and no step has
  `--split test`." Add to I22: "`--budget-usd` equals front 6.4 rule 5 (read from one constant or from the front
  table), not a literal in the test."

### 2. [major] [B, A] B1.9 prices Gate B1 with a noise model that is not a bound on the gate's noise and is 13x the measured per-arm noise

- Problem.
  - B.md, B1.9 "Noise model": "Each simulated B1 value is the simulated mean plus noise with SD
    `sqrt(mu_i (1 - mu_i) / 8)` ... This is an exact upper bound ... The baseline values stay fixed (they exist).
    ... The rules are fixed from the bound run." And: "If joint P(pass) < 0.5 for the 25 % effect, the owner decides
    before the rental whether to change a rule."
  - It is not a bound on what the gate sees. The primary baseline is A3's fresh K = 8 sampling of R5 (B.6.2
    "Paired baseline": "Primary: A3's R5 val items"). So the gate's per-item difference has two sampling noises.
    Even in the fallback (R5's file), the simulated difference cancels R5's realized noise, which the real gate
    keeps. The one-arm model can give half the true variance.
  - In practice the bound is far too large. Two independent K = 8 samplings of the same model on 128 val-shifted
    items differ with variance 0.00324, which is 7.7 % of the two-arm bound and 15 % of B1.9's one-arm bound
    [measured: `r4data/retest.py` on `runs/rl/eval/final/heldout.json` and `runs/r5/eval/eval_val/items.jsonl`].
  - Effect on the owner's choice: B1.9 as written gives joint P(pass) about 0.13. Measured noise on both arms gives
    about 0.92. The false-pass rate stays low (IoU 0.03, DiD 0.027 with no effect) [estimate: `r4data/power2.py`].
    So B1.9 as written triggers the "< 0.5: change a rule" path and invites the owner to loosen rules that already
    have power. B1.9 itself says "Nothing from A3 is needed", but it does not use the measured test-retest data
    that is already on the Mac.
- Fix. B.md B1.9, replace the "Noise model" bullet with:
  "Noise model (Mac, $0). Both arms get noise, because the baseline (A3) is a fresh sampling too. Item i gets
  independent noise on each arm with variance `c x mu_i (1 - mu_i) / 8`. `c` is measured: two independent K = 8
  samplings of `rl-full-3ep` on 128 `lidc-val-shifted` items (`runs/rl/eval/final/heldout.json` against
  `runs/r5/eval/eval_val/items.jsonl`, same seed, two boxes) have a per-item difference variance of 0.00324 =
  0.077 x mean(2 mu (1 - mu) / 8) [measured]. The rules are fixed from the run with c = 0.077. The run with c = 1
  (the Bernoulli bound on both arms) is printed beside it as a stress line. After B-S1, re-run with A3's split-half
  `item_sd` (A3 step 7), for the record only."
  Delete "This is an exact upper bound: Var(X) <= mu (1 - mu) for X in [0, 1]". Add a B1.9 test: with c = 0 on the
  baseline arm and c = 1 on the B1 arm, the result equals the old one-arm run. Add the 128-item retest numbers to
  B.6.2 "Noise floor" as the pre-V1 reference: paired per-item SD 0.0569 between two samplings
  (per-arm variance 0.00162 against the bound's mean 0.0211, a factor of 13).

### 3. [major] [B] The disagreement weight with "no window twice in a run" does not give 14 %, and it front-loads disagreement rows

- Problem.
  - B.md, B1.2 step 3: "draw a window of the slot's kind with weight `disagreement_weight` for disagreement items,
    else 1 (replaces `DISAGREEMENT_OVERSAMPLE = 4` row repetition)". Step 4: "Reject a window ... that is already
    used in the run." "Weight 4 inside each cell gives about 4p / (1 + 3p) = 14 % of rows [estimate: arithmetic]
    ... Default: 4, so B1 changes only the geometry against R5." DG14: "Disagreement share (expected about 14 % at
    weight 4)".
  - R5's rows repeat each of the 289 disagreement items 4 times and shuffle: a stationary 4 x 289 / 8,208 = 14.1 %
    in every step [code: src/objdel/train/dataset.py:81, 235-237; measured: 289 / 7,341 in train.jsonl].
  - B1.2 draws without replacement. Heavy items run out first. For the 7,200 shifted rows (15 blocks of 480)
    [estimate: `r4data/disagree.py`]:

    | `draws` | overall share | block 1 | block 8 | block 15 |
    |---|---|---|---|---|
    | 1 | 4.0 % | 17.7 % | 3.3 % | 0.0 % |
    | 2 | 7.5 % | 14.6 % | 5.6 % | 2.7 % |
    | 4 | 10.8 % | 12.3 % | 9.8 % | 7.7 % |
    | 8 | 12.0 % | 14.2 % | 9.8 % | 8.1 % |

  - So at the small `draws` that B1.4 step 2 prefers ("take the smallest value that passes DG7-DG9"), B1 trains on
    about half of R5's disagreement share, mostly in its first blocks. B1 then differs from R5 in geometry **and**
    in the disagreement exposure over time. That breaks the stated reason for the default ("B1 changes only the
    geometry"), weakens Gate B1's data-effect reading, and makes the owner choose weight 4 vs 1 (B.8.3 item 3) on a
    wrong number. DG14 is a report on the whole schedule, so it cannot see the per-block drift.
- Fix.
  - B.md B1.2: replace step 3 with a quota. "Per block and per (kind, cell) slot set, the disagreement slot count is
    the largest-remainder rounding of `q x slots`, q = 4p / (1 + 3p) with p the disagreement share of that (kind,
    cell) pool (weight 1: q = p). Fill the disagreement slots from the disagreement windows of that cell, the other
    slots from the rest, each uniform without replacement. When a cell's disagreement windows run out, a
    disagreement window may be used again, at most 4 times in the run (R5 repeated each disagreement row 4 times);
    beyond that raise `ValueError` naming (kind, cell, 'disagreement')."
  - B1.4 step 2: choose `draws` so that DG7-DG9 pass and no disagreement window is used more than once (record the
    largest reuse in the manifest).
  - DG14: split out "DG14a, hard: the disagreement share of every 20-step block is within ±1.5 points of q's
    block value". Change "expected about 14 %" to "q per cell; 14.1 % overall at weight 4 (R5's share)".
  - B.9 property test: the disagreement share of every block prefix equals the quota counts exactly.

### 4. [major] [B] Gate B1's far-miss point rule has two values

- Problem. B.md, B.6.2, Gate B1 header: "'R5' below means the paired baseline of the common rules: A3's value, with
  R5's 0.618 for reference". Same bullet list: "val misses at offset > 0.3: point <= 20.5 % (25 % below R5's
  27.3 %)". Under the header's definition, "25 % below R5" means 0.75 x A3's far-miss share. The bullet also fixes
  20.5 %. The two differ whenever A3's far-miss share differs from R5's file, and also after the
  `known_label_issues` exclusion changes the item set. The measured cross-host shift of the same model is -0.0076
  in mean IoU on 128 items [measured: `r4data/retest.py`], so A3's far-miss share will not equal 27.3 % exactly. An
  engineer must pick one. Near the threshold the choice decides pass or fail. B1.9 must simulate the same rule.
- Fix. B.md B.6.2 Gate B1, replace the bullet with: "val misses at offset > 0.3: point: B1's far-miss share <= 0.75 x
  the paired baseline's far-miss share on the same items (after the `known_label_issues` exclusion); with R5's file
  that is 0.75 x 27.3 % = 20.5 %; and the CI upper bound of (B1 - baseline) far-miss share < 0." Make the same edit
  in B1.9 ("apply every rule exactly as B.6.2 writes it") and in the Gate B3 analogue if it has a fixed number (it
  does not today: "falls by >= 25 % relative" is already relative).

### 5. [major] [B] B1.5's V1 commands omit required flags, so `train_rl` and `final_eval` exit 2 on the box

- Problem.
  - B.md, B1.5 step 3: "`train_rl --ckpt {rd}/ckpt` ... `--train {data}/rows/lidc-train-b1 --rows-revision <hex40>
    --rows-manifest-sha256 <hex64> --out {rd}/b1-rows --numerics bf16-seqmask --lr 5e-5 --steps <planned_steps>
    --save-steps 50 --shared-prefix --budget-usd 2.05`". No `--model`. Code: `ap.add_argument("--model",
    required=True, help="configs.MODELS key")` [code: experiments/phase2/train_rl.py:24].
  - B1.5 step 4: "`--ckpt {rd}/b1-rows/final --view {data}/eval_views/lidc-val-shifted --split val --out
    {rd}/eval_b1_val_shifted --view {data}/eval_views/lidc-val-edge --split val --out {rd}/eval_b1_val_edge`". No
    `--train`. Code: `ap.add_argument("--train", type=Path, required=True, ...)` [code:
    experiments/phase2/final_eval.py:50]. E removes `--train` only in E4, after V1 (E.md: "E4 (before E6) deletes
    final_eval's --train").
  - B1.5 is the only text that specifies these two steps (front 6.2: "B1.5 writes the commands"). A3's plan for
    the R5 process lists `--train` explicitly, so the gap is only on the B side. With step 4 as written, front item 3
    makes no output, and Gate B1 has no B1 input until a later rental re-runs it. No Mac test parses the plan's argv:
    I22 checks flag names only, and only after item 1's fix.
- Fix. B.md B1.5 step 3: add `--model qwen3-vl-2b` first. Step 4: add `--train {data}/prepared` (the same root as
  `a3_measure.yaml`). Add a test to I22 (Mac, no GPU): "for every `-m` step of `b1_rows.yaml` and
  `a3_measure.yaml`, build the module's argparse parser and parse the step's argv with `{rd}` and `{data}`
  replaced by temp paths; no `SystemExit`".

### 6. [minor] [B] B1.9 names the wrong decision interval for Gate B1

- Problem. B.md B1.9: "The decision interval of the simulation is the CR1 interval of the item-weighted mean
  (B.6.2), the same as the gate's." B.6.2 "Pairing and CI": "the percentile patient interval when the gate set has
  >= 50 patients (Gate B1, whole val, 78 patients), and `signal.cluster_t` when it has < 50". The far subset has 75
  patients and the near subset 52 [measured], so every Gate B1 rule uses the percentile interval. The DiD rule has
  no CR1 formula in A1 at all. The earlier B1.9 sentence ("apply every Gate B1 and Gate B3 rule exactly as B.6.2
  writes it ... the A1 decision interval") already gives the right answer, so this is a wording error.
- Fix. B.md B1.9, replace the sentence with: "The decision interval of the simulation is the one B.6.2 names for
  each rule: the percentile patient interval for a set with >= 50 patients (all Gate B1 rules today: whole val 78,
  far 75, near 52), and `signal.cluster_t` below 50 (Gate B3 on val-select)."

### 7. [minor] [B] A B1.1 pool test cannot pass: S0 share >= 0.75 is a schedule property

- Problem. B.md B1.1 tests: "S0 windows are >= 0.75 of off-centre windows". The pool holds one shifted window
  (side S0) and one rescaled window (side != S0) per item and draw. `neighbour_side` tries only the LADDER
  neighbours of S0 [code: src/objdel/targets/augment.py:97-104, `_neighbours`]. So the pool's S0 share of off-centre
  windows is about 0.5 by construction. The 0.75 limit belongs to the schedule (DG10: shares 0.75 / (0.75 + 0.10) =
  0.88).
- Fix. Move the bullet to the B1.2 tests: "the schedule's S0 share of off-centre rows is >= 0.75 (DG10)". In the
  B1.1 tests write: "for every draw, an item that has a rescaled window also has a shifted window, unless
  `shifted_none` counts it".

### 8. [minor] [B] B1.7 and DG17 cite two "flagged" val patients that our own checks do not flag

- Problem. B.md B1.7: "The third-party duplicate-z list puts two patients in val: LIDC-IDRI-0085 (2 items) and
  LIDC-IDRI-0146 (15 items), 1.6 % of val". DG17 "Today": "the third-party duplicate-z list names 2 val patients
  (17 items)". The doc does not say whether these 17 items go to `known_label_issues`. That choice changes the Gate
  B1 item set. Measured today: check (a) on the local XML flags 3 store series, all already in the sidecar list
  (380786 train 1 z, 228428 train 5 z, 074972 test 11 z) and no val series. The pylidc `zvals` hold no duplicate z
  for any of the 800 store series. The two val sidecars match (9/9, 15/15 slices) [measured: this review]. So by the
  B1.7 steps, 0085 and 0146 are not hits.
- Fix. B.md B1.7, after the third-party sentence, add: "The list is context only. A series is a hit only when
  steps 1(a)-(d) flag it. Today (a) flags 380786, 228428 and 074972 (1, 5 and 11 shared z values), all also (d)
  hits; 0085 and 0146 are clean on (a), (b) source data and (d), so their items stay in every Gate B1 rule." DG17
  "Today": replace the third-party clause with "check (a): 3 series (2 train, 1 test), 0 val".

### 9. [minor] [A] A3's TOST power is said to be unknown; it can be estimated now, and the risk is a host shift, not noise

- Problem. A.md A3 step 5: "The power of this test is not measured [estimate: the paired CI is narrower than the
  unpaired one, by how much is unknown]." The 128-item retest gives a paired per-item SD of 0.0569 and a cross-host
  mean shift of -0.0076 (SE 0.0050) [measured: `r4data/retest.py`]. On 1,084 items the 90 % half-width is about
  1.645 x 0.0569 / sqrt(1,084) = 0.0028 before patient clustering [estimate]. So the TOST passes when the host
  shift is below about 0.007, and the shift measured between two earlier boxes is at that edge. A failure is then
  likely a host or engine difference, which A3 step 5 already routes (Gate B1 still pairs with A3's items).
- Fix. A.md A3 step 5, replace the power sentence with the numbers above and: "If the TOST fails, report the mean
  shift with its CI. Gate B1 is not affected (same box, same view as B1)." Same line in A-G3.

### 10. [minor] [B] Gate B1's +0.018 decomposition uses the wrong found IoU and treats a miss as 0

- Problem. B.md B.6.2 Gate B1: "at the found IoU 0.706 that adds about 28 x 0.706 / 1,084 = +0.018". 0.706 is the
  found IoU of all val items. Far found items average 0.688, and far misses already score 0.012 on average
  [measured]. The gain is 28 x (0.688 - 0.012) / 1,084 = +0.017. B1.9's simulation draws from "the found-image
  distribution of the same cell", so it is consistent; only the prose is off.
- Fix. Replace the arithmetic with "28 x (0.688 - 0.012) / 1,084 = +0.017 (far found IoU 0.688; far misses 0.012)
  [estimate: arithmetic on measured values]".

### 11. [minor] [A] A3's header comment bundles after it builds the test view, so front 6.1 check 4c fails

- Problem. A.md A3 header comments: "build_eval_views --split test --kinds shifted --out data/eval_views ..." then
  "... vast.bundle lidc --view lidc-bundle-v3". The bundle ships every `data/eval_views/lidc-*` folder [code:
  experiments/phase2/vast/bundle.py:53-56], and A3's own text says so ("`lidc-test-shifted` reaches the box when it
  exists on the Mac at bundle time"). 00_front.md 6.1 step 4c requires "no `lidc-test-*` folder". The check catches
  it, but the listed order makes it fail every time.
- Fix. A.md A3 header comments, insert before the bundle line: "`mv data/eval_views/lidc-test-shifted
  data/frozen_views/` (front 6.1 step 4c; move it back before V3 item 3)", or name the bundle `--view` filter if B
  adds one.

---

Not repeated (fixed in earlier rounds, checked in the current text): the one V1 list and the A3 plan split (r3
data 1), the primary baseline and its fallback (r3 data 2), the val-half seed search (r3 data 3; re-run of
`valhalf.py` logic gives s = 78, 533 / 551 items), CR1 `cluster_t` (r3 data 4), the DiD informative check (r3
data 5), val label hits (r3 data 6), R5 selected on val (r3 data 7), the LIDC-IDRI-0659 decision (r3 data 8), the
patient-CI golden (r3 data 9), the near-duplicate score (r3 data 10), stop rule 7 (r3 data 12). The A2 voter rule
and the k = 2 stack were re-checked against `cache.py:103-111, 203-207` and hold. The R5 stamp chain was re-checked
against `preflight.py:59-69` and holds.
