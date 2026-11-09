# Final review round 6: fresh outside architect

Scope: the whole plan (`docs/implementation2.0/Implementation.md`, built from `phases/*.md`), read by section. Focus:
is this the best end-to-end architecture up to Phase E with GPU everywhere, and what would a top lab do
differently. I checked the round 1-5 review logs and did not repeat any item they already fixed (the money path
under (b'), V4's eval path, E6's SftStage, spec inputs and push target are all fixed in the round-5 sources).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | This review needs no NVIDIA catalog skill. The front already names `compute-sanitizer-diagnose-cuda-bugs`, `nv-segment-ct` and `dicom-series-to-volume` as recommendations only. I did not install anything. I used a web check (nvImageCodec) for the front's claim that no GPU PNG decoder exists. |
| engineering:code-review | Every quote below was re-read in code: `final_eval.py:49-84` (one view, fork pool), `select_ckpt.py:223-231`, `serve.py:132-161`, `policy.py:425-536, 595-610` (chain, run cap, push/fetch, CLI defaults), `runs/r5/shared_rl.json` (23.653 s/step, 15 steps), `ledger.csv` (133 rows, $16.9025), and the docker locks (box transformers 5.19.0, Mac 5.20.1; C.md already records the split). |
| engineering:debug | I traced the default-budget path (b') to its end: V1 runs, then nothing more is rented. Then I asked which planned commits can still run. That trace gives item 1. |
| engineering:system-design | I judged the end-to-end flow (data, V1, Mac gates, C/D/B3 code, E). The weak point is not a box or a number. It is the missing end state for code whose gate needs a rental that the default never funds. |
| engineering:architecture | I checked the rule hierarchy (front order rules 1-5, "this section wins", E's "When E starts") to find whether a rule already decides item 1. None does. |
| engineering:documentation | I checked that each quoted line exists word for word in the named source file and heading. I left out review-log text. |
| huggingface-skills:huggingface-best | No model-choice change. Section 3's Qwen3-VL-2B decision and the 4B gate hold. No Hub leaderboard covers polygon text on CT. |
| system-design-skills:system-design | I used the "state the assumption, then the end state" check. Front 6.3 states the (b') assumption, but no phase states the end state it leads to (item 1). Item 7 also comes from this check: one pre-registered conditional read instead of a later standalone rental. |
| engineering:tech-debt | Applied to the default end state. Two dual paths stay behind temporary flags with no planned removal under (b') (item 1). `e_v4.yaml` would be a typed spec outside `specs/` (item 6). |

Web checks (2026-11-08):
- nvImageCodec lists only nvJPEG, nvJPEG2000 and nvTIFF as GPU codecs. PNG goes through the CPU fallback, so the
  front's "no GPU PNG decoder" holds [web: https://developer.nvidia.com/nvimagecodec].
- `cupyx.scipy.ndimage.distance_transform_edt` is an exact GPU EDT (Parallel Banding Algorithm) with a scipy-style
  API and `float64_distances=True` by default [web:
  https://docs.cupy.dev/en/latest/reference/generated/cupyx.scipy.ndimage.distance_transform_edt.html]. It feeds
  item 8.

## Overall judgement

The architecture is sound for this scale. Its parts:
- One exact CUDA engine source that builds for both the CPU and nvcc (the CPU build is the oracle).
- Data first, with frozen hashed views and patient-cluster paired gates.
- One paid experiment (V1) that changes one variable.
- Typed specs, import-linter layers and a single TRL adapter.

A top lab would make the same core choices. Where it would differ:
- It defines what the repo looks like when the money stops (item 1).
- It folds a pre-registered conditional test read into the one funded rental instead of paying for a second launch
  (item 7).
- It takes a library where the plan writes its own kernel (item 8).

Counts: blocking 1, major 0, minor 7.

## Items

### 1. BLOCKING. Under the default (b'), no part says what happens to code whose gate needs V2 or V3. E's start is undefined, and old and new code stay side by side with no end.

Sources: 00_front, B, C, E (and 99_back decision 1).

Problem. Round 5 made (b') the default: no rental follows V1. The time line still makes every engineer land the
C, B3 and D code, and still starts E only after the "last" A-D commit.
- 00_front, section 1, time line step 5: "C1-C5, B2/B3 and D2-D8 code land on the Mac, on pre-E paths."
- Step 6: "V2 and V3 run only after the owner raises the budget (owner decision 1, section 6.3). Under the default
  (b') neither runs."
- Step 7: "E moves files after the last A-D commit that edits `fast.py`, `dataset.py`, `reward/` or `checker/`."
- Order rule 1: "Every A-D task path is pre-E. ... There is no exception".
- E.md, E.4 "Order of execution / When E starts" says the same: "E0 starts after the last A-D commit that edits
  `fast.py`, `dataset.py`, `reward/` or `checker/`."

Two planned A-D commits edit those paths and can run only after a funded gate:
- B.md B3.10 "Clean after gate B3": "Pass: v5 becomes the only reward. Move `RewardV4`, `advantages_v4` ... to
  `archive/reward_v4/` ... Fail: keep v4, archive v5 code to `archive/reward_v5_try/`". B1.8 adds: "The same rule
  applies to v4 (B3.10): old code goes only after the gate that replaces it passes." Gate B3 is read after the B3
  run, which is V3 item 4 (front 6.3). Under (b'), V3 does not run.
- C.md C10: "Files and actions (after each gate passes): delete from `src/objdel/train/fast.py` ..." and "remove the
  temporary flags `--gpu-pixels`, `--encoder-once`, `--host-meta`, `--fast-sync`, `--prompt-memo` once their path is
  the default." C1 says: "`train_rl --gpu-pixels {auto,off}` (temporary; default `auto` after the gate)". C.6 says:
  "**Gate C-adopt** ...: G-EQ1..6 and G-SP1 pass in V2". C.8.2 confirms that V2 runs only after a raise:
  "under the default (b') it does not run."

Code facts (verified): today `fast.py` and `reward/` are the live training path. The C flags and the v4/v5 split do
not exist yet, so the plan alone decides their fate.

What an engineer gets when following the doc exactly under (b'):
- (a) E0 can never start, because B3.10 and C10 stay "pending" forever. Or E starts, and B3.10 and C10 later land
  on post-E paths. Order rule 1 forbids that ("no exception"). E.7's escape covers only "A B/C/D fix that lands
  after E", not a planned cleanup task.
- (b) The repo keeps reward v4 beside v5, and the CPU pixel path beside the GPU pixel path behind temporary flags,
  with no planned removal. This breaks the owner rule "Old code must be ARCHIVED (not left beside new code)"
  (BRIEF section 3).
- (c) The GPU pixel path never becomes the default, so the "GPU everywhere" end state is never reached.

Two parts that both claim to be the rule (the time line with order rule 1, and the (b') default) cannot both be
met. That is the blocking definition.

Exact fix:
1. 00_front, section 1: add a paragraph "End state under (b')" after time line step 8. Under (b'), these land and
   finish on the Mac:
   - A1-A9;
   - B1.0-B1.9 (B1.8 on its rows);
   - B2.0 and B2.1;
   - the $0 part of B3 (B3.1-B3.4, B3.6 and the B3.7 lab round);
   - D1-D9 Mac parts (D-G4, D-G9, D-G11 and D-G16 are Mac gates);
   - E, DOC and CLEAN.

   These are written only on an owner-approved branch, or not at all, until a raise funds their gate: C1-C5, C10,
   B2.2-B2.4, B3.5 wiring, and B3.9-B3.10. So `main` never holds a temporary flag or a second reward.
   The other option: land them, and have C10 and B3.10 archive the unadopted new code (not the old code) at E0 with
   a README row "not adopted: no budget". Either way, state one rule.
2. 00_front, order rule 1, and E.md "When E starts": replace "after the last A-D commit" with a named commit list
   (the list in fix 1). Add: "B1.8 (when its row waits), B3.10, C7, C8, C10 and C11, when they run after E, use
   E's move tables as the path map (E.7 row)." Change E.7's row text from "fix" to "fix or deferred cleanup task".
3. 99_back decision 1, option (b'): add one line that points to the new end-state paragraph, so the owner sees
   which Mac work the default skips.
4. C.8.3 "Phase C blocks" and B.8.4: one sentence each that points to the end-state paragraph.

### 2. MINOR. The front 2.4 diagram row count disagrees with section 0.

Sources: 00_front.

00_front 2.4 diagram: "-> B1 rows root (9,600 rows)". Section 0, "rows root": "B1: 8,320 rows = 260 steps x 32".
B.md B1.2 says: "the B1 rows root holds 8,320 rows. The 9,600-window sizes and times in B1 are the 300-step upper
bound." Section 0 wins, and B1.5's `len(rows) == schedule.steps * schedule.per_step` assert catches a wrong build,
so the impact is low.
Fix: in 2.4, change "(9,600 rows)" to "(8,320 rows)".

### 3. MINOR. B.8.3 item 11 gives the wrong step count for the $2.20 fallback.

Sources: B (vs 00_front, 99_back).

B.md B.8.3 item 11: "Keeping the older $2.20 means a 260-step B1 run (front 6.2 notes)." Front 6.2 notes: "If the
owner keeps $2.20 anyway: the cap is 3.05 h = 183 min, `--budget-usd` is $1.85, and the run plans 220 steps".
99_back decision 2 says: "The older $2.20 would cut the B1 run to 220 steps". Front 6 says "this section wins", and
$2.20 is not the default, so this is minor.
Fix: B.8.3 item 11, change "a 260-step B1 run" to "a 220-step B1 run with `--budget-usd 1.85`".

### 4. MINOR. B-S2's item list has three spellings.

Sources: B, 00_front.

- B.md B.8.1: "**B-S2 = rental V3 items 0, 2a, 2b and 4, front 6.3**".
- 00_front 6.3 bullet: "B-S2 in B.8.1 is V3 items 0, 2a and 4 (plus 5 when the gate is closed)".
- The 6.3 table row: "**B-S2 = items 0 + 4**".

Item 2b (fp16 arm) is optional (99_back decision 16, default off).
Fix: one definition in 00_front 6.3: "B-S2 = V3 items 0 and 4; optional 2a, 2b; 5 only when the 4B gate is
closed". B.8.1 points to that line only.

### 5. MINOR. Under (b'), the final candidate and the test read are implied, not stated.

Sources: 00_front (1.1, 6.3), 99_back (decisions 2b, 17).

Front 1.1 row "Pass": final candidate = "Best of B1 and B3 on val-select by the Gate B3 rules". Under (b') no B3
exists. Front 6.3 V3 item 5: "The one test eval of the final candidate, only if the owner has closed the 4B gate".
Under (b') V3 never runs. The suggested 4B gate (val IoU_tc < about 0.68 and misses > about 7 %, against R5 val
0.618 / 12.6 %) is likely to fire, so test then "waits for the 4B session", which also needs a raise. An engineer
can infer that, under the default, no test number for the data fix is ever produced and the headline stays R5's
0.645. The plan should say so, because DOC's `results/headline.json` and DOC11 "Known issues" must carry it.
Fix: add one row note under the 1.1 table: "Under (b'): the final candidate is B1 if Gate B1 passes, else R5. No
test eval runs (V3 item 5 waits for a raise). The headline stays R5's test numbers, and B1 is reported on val only,
labelled 'val, single run'." DOC11 adds the same line.

### 6. MINOR. `e_v4.yaml` is a `RunSpec` outside `specs/`.

Sources: E, 00_front (2.3).

00_front 2.3 tree: "specs/  typed run specs (replace the old YAML plans)". E.md E6 step 8b: "E writes
`experiments/phase2/configs/e_v4.yaml` ... so the file is a `RunSpec` at the path front 6.3 names". E6 step 7's
test is "every `specs/*.yaml` loads", so it does not cover the one spec that a rental would run. U51 covers it
separately, and DOC4's glob includes `experiments/**/configs/*.yaml`, so nothing breaks. But two homes for typed
specs is the confusion the owner's archive rule aims to remove.
Fix: E writes `specs/e_v4.yaml`. Front 6.3 V4, 99_back 3b, E.8 and U51 name that path. Front section 0 "plan file"
gets one line: after E6, plan files are `specs/*.yaml`.

### 7. MINOR (top-lab design). Pre-register a conditional test read inside V1 instead of paying a second launch later.

Sources: 00_front (6.2, 6.3, section 3), 99_back (decisions 2, 17).

Today a test read of B1 needs a separate rental. Its fixed part alone is about 39 min plus b (front 6.3 V3 item 0:
$0.54), and the read itself is 10.6 min ($0.12). Front 6.3 already says, for V3: "The 4B gate numbers are set before
the spend, so the owner can close it in the session from the val evals."

The same pre-registration fits V1:
- Add a step 5 to `b1_rows.yaml`: a small `$0` script reads `eval_b1_val_shifted/final_eval.json`.
- Only if Gate B1's point rules pass and the owner's 4B thresholds say "do not run 4B" (val IoU_tc >= 0.68 and val
  misses <= 7 %), the script runs `final_eval` on `lidc-test-shifted`.
- Arithmetic [estimate]: the V1 amount goes $2.40 -> $2.52, the worst case $2.61, plus the $0.25 reserve = $2.86
  <= $3.00. `--budget-usd` stays $2.05, because rule 5 subtracts later items, and the extra 10.6 min sits in the
  tail.

This keeps "test read once, on the final candidate, after the gate closes": the gate closes by a rule fixed before
the spend. It saves about $0.54 of launch overhead in the case where B1 wins outright.
Fix (owner option, not default): add it as option (b'') in 99_back decision 2, with the numbers above. If the owner
takes it, the 6.1 step 4c test-view exclusion becomes "ships `lidc-test-shifted`, frozen hash checked" for V1.

### 8. MINOR (owner library rule). Name `cupyx.scipy.ndimage.distance_transform_edt` beside K2's own EDT.

Sources: 00_front (section 4, EDT row), D (D7).

Front section 4: "KEEP our K2 separable EDT on the device ... cuCIM is not adopted; it is an optional second oracle".
The owner rule says: "Prefer a good library (even with a little customization) over our own code." CuPy ships an
exact GPU EDT with a scipy-compatible call and float64 distances by default (web check above). The Mac oracle stays
scipy either way, so the "one source for CPU and nvcc" argument covers only the polygonizer and validity totals,
not EDT. CuPy is not in either box lock today [code: grep of `experiments/phase2/docker/locks/*.txt`], so adopting
it adds a dependency. That is why this is minor.
Fix: in the section 4 EDT row and D7, add "Considered: `cupyx.scipy.ndimage.distance_transform_edt` (exact,
scipy-compatible). Trigger to switch: K2 EDT fails D-G10 bit identity, or K2 EDT costs more than 1 engineer-day to
keep." Then CLEAN14's keep-list ADR (0021) records the choice.
