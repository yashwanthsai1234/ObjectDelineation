# Final review round 10: fresh outside architect

Scope: the whole doc `docs/implementation2.0/Implementation.md` (round 10), read by sections, with the focus on the
end-to-end architecture up to Phase E (front sections 1-8, the money path of 6.2-6.4, the D and E design and gates).
I checked every proof against the code or a measurement. Items already in `review-log.md` (rounds 1-9) are not
repeated.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item depends on an NVIDIA product choice. The front already lists the catalog skills for the 3D direction. I did no new catalog lookup and installed nothing. |
| engineering:code-review | I read `experiments/phase2/vast/policy.py:420-455` (`chain`: a scoring chain gets the CUDA build plus `tests/test_gpu_scoring.py` first), `experiments/cuda_checker/run_box.sh:1-40` (`--build-only` exists), `cuda/gpucheck_py.cpp:113` (`build()`), `experiments/phase2/train_rl.py:30-65`, `experiments/phase2/g_probe.py:74-100` and TRL's `vllm_generation.py:600-760` (colocate mode sends `n=1`, so 256 single requests today). These confirm the front claims I relied on. |
| engineering:debug | I traced two money paths to the end. (1) Under the default (b'), no rental after V1 reads test (item 1). (2) Under a raise, the V3 launch lines run items 1 and 3 before the B3 run, but the decision-1 total holds only items 0 and 4, so `BudgetStop` cuts the B3 run (item 2). I re-did both sums with the front's r_s = $0.01133/min and b = $0.10. |
| engineering:system-design | I checked the whole flow end to end: Mac data gates, then V1, then the Gate B1 decision table, then the branch `after-raise`, then E. The weak point is not the design of any one box. It is what the default path delivers: no test number for the B1 model (item 1), and branch code written on pre-E paths that must then be ported (item 5). |
| engineering:architecture | I compared how box-only gates are treated. C code whose gate needs a box stays off `main`. D code whose GPU gates (D-G3 box half, D-G5, D-G6) need a box lands on `main` (item 4). I wrote the trade-off and the fix in ADR form (option, consequence). |
| engineering:documentation | Every blocking or major item quotes the source file and heading, plus the code line or arithmetic that proves it. |
| huggingface-skills:huggingface-best | No item is about the model choice. I checked the HF private-storage limit because the rows root and the checkpoints live in the private repo: 100 GB on a free account [web: huggingface.co/docs/hub/storage-limits]. The plan needs far less (rows root about 1-2 GB, a few 2B checkpoints), so there is no item. |
| system-design-skills:system-design | I applied "quantify before choosing" to the leftover money. A separate rental pays about 38 min of launch before any work, which is $0.54 at r_s. So the cheap evidence (a test eval) belongs inside V1, not in a later rental (item 1). |
| engineering:tech-debt | I looked at the debt that the default (b') creates on purpose: C1-C11 and B2.2-B3 written on pre-E paths, on a branch that never merges before E rewrites the paths and the plans (item 5). |
| WebSearch | HF storage limits (above). No other external claim changes an item. |

## Items

### 1. major: under the default (b') the B1 model never gets a test number, although the $3 can buy it inside V1

Sources: 00_front, B, A, 99_back.

Problem.
- 00_front "6.3 V2, V3 and V4": "(b') keep $3: V1 only. V2, V3 and V4 wait for a budget raise. The money left
  ($0.46-0.65) stays unspent as a second reserve."
- The only test eval of a final candidate is V3 item 5 (00_front 6.3 V3 table): "The one test eval of the final
  candidate, only if the owner has closed the 4B gate". V3 never runs under (b').
- 00_front 6.1 step 4c forbids a test view in V1 ("no `lidc-test-*` folder"), and step 4d forbids `--split test`
  in the V1 plan files.
- 00_front 1.1 row "Pass": the final candidate is "Best of B1 and B3". Under (b') B3 never runs, so B1 is the final
  candidate. But no rental is left that reads test for it.
- A separate test rental does not fit the money left. It costs launch to first `run` 38.2 min [measured, 6
  inputs] + build 3 + test eval 10.6 + imports 2.6 = 54.4 min, which is $0.10 + 54.4 x $0.01133 = $0.72
  [arithmetic]. The money left is $0.46-0.65.

Result. An engineer who follows the doc exactly under the default ends Phase E with a B1 model that passed Gate B1
on val and has no test IoU. The README headline (DOC6) stays R5's 0.645. The plan's main question (does the data
fix close the shifted-test misses?) has no test answer. A top lab would buy the cheap evidence inside the rental that
is already paid for.

Proof that it fits inside V1. Item 3's process already loads the B1 final model (00_front 6.2, item 3). Adding
`lidc-test-shifted` (2,200 items) to that same process costs 2,200 x 0.193 s = 425 s = 7.1 min = $0.08 (the front's
per-item rate; no second fixed load). The view already exists and is hashed before the rental (A3:
`build_eval_views --split test --kinds shifted`, `runs/a3/pre_rental.json`). With V1 amount $2.48 and `BudgetStop`
kept at $2.05:
- cap = floor_to_0.05((2.48 - 0.10) / 0.680) = 3.50 h = 210 min;
- eval tail = 26 + 7.1 = 33.1 min, and 172.3 + 33.1 = 205.4 min <= 210 min;
- worst case = 0.10 + 0.680 x (3.50 + 0.167) = $2.59; `launch_refusal`: $16.90 + $2.59 = $19.49 <= $19.90;
- with the $0.25 reserve: $2.59 + $0.25 = $2.84 <= $3.00.

Exact fix (00_front 6.1 step 4c and 4d, 6.2, 6.3 totals; B B1.5 step 4; A3; 99_back decisions 1 and 2):
- Add V1 item 3b: in item 3's process, `--view lidc-test-shifted --out {rd}/sealed_test_b1`. Make the rule
  "sealed": `pull.sh` brings the folder back, and the Mac records only its SHA-256 in `runs/<V1>/sealed.json`.
  Nobody opens it until the decision-table row of 1.1 names B1 the final candidate and the owner closes the 4B
  gate. If R5 or B3 becomes the final candidate, the folder is deleted unread. This keeps "test is read once".
- 6.1 step 4c: allow exactly `lidc-test-shifted` in the V1 bundle, with its hash equal to `pre_rental.json`. Step
  4d: the grep test allows exactly one test view, only in `b1_rows.yaml` step 4, only with `--out
  {rd}/sealed_test_b1`.
- 6.2: V1 approved amount $2.48, cap 3.50 h, worst case $2.59. 6.3: money left is $0.38-0.57. Decision 2 changes to
  match.
- If the owner rejects the sealed eval, state it in "End state under (b')" and in decision 1: "Under (b') the B1
  model has no test number; the headline stays R5's; a test-only rental costs about $0.72 and does not fit."

### 2. major: the V3 price (decision 1 (a)) omits V3 items 1 and 3, which the launch lines run before the B3 run, so `BudgetStop` cuts B3 to about 86 of 260 steps

Sources: 00_front, B, 99_back.

Problem.
- 00_front "6.3", V3 table: "**B-S2 = items 0 + 4** | **182.6 = $2.17**". "Totals for owner decision 1": "(a) B1
  then B3: V1 + B-S2 + reserve = $2.29 + $2.17 + $0.25 = $4.71 at 27.5 s". 99_back decision 1 (a) points at these
  totals.
- The V3 launch lines (00_front 6.3) run item 1 (`d_gpu_test.yaml`, 27-46 min, $0.31-0.52) and item 3
  (`a3_test.yaml`, 13.6 min, $0.15) before item 4 (`b3_v5.yaml`). Neither line says "only if approved" (only item
  2b does).
- 00_front 6.4 rule 2 stops a skippable item only when "Session spend reaches 90 % of the amount the owner
  approved". Rule 5: the B3 `BudgetStop` = "approved amount - $0.35 - the $ of every later item".
- There are three definitions of B-S2. 00_front 6.3 table: items 0 + 4. 00_front 6.3 bullet: "B-S2 in B.8.1 is V3
  items 0, 2a and 4". B B.8.1: "B-S2 = rental V3 items 0, 2a, 2b and 4".

Arithmetic (27.5 s, front r_s). Take V3 approved = $2.17, the decision-1 (a) B-S2. Then item 0 $0.54 + item 1 $0.52 +
item 3 $0.15 = $1.21. All three run, because $1.21 < 90 % x $2.17 = $1.95. `BudgetStop` = $2.17 - $0.35 = $1.82.
The B3 build, smoke and start-up take 14 min, so spend reaches $1.37 before step 1. That leaves $0.45 = 39.7 min =
about 86 steps at 27.5 s. B3.9 requires "the same `--steps` as the B1 run (260, front 6.2)". Gate B3 then compares
an 86-step B3 with a 260-step B1. That is a wrong result, and the B3 money is spent for nothing.

Exact fix (00_front 6.3; B B.8.1; 99_back decision 1):
- Add a row "**V3 as launched (items 0, 1, 3, 4)**": 39 + 46 + 13.6 + 143.6 = 242.2 min = $2.84 at 27.5 s
  (39 + 46 + 13.6 + 126.9 = 225.5 min = $2.66 at 23.65 s); with item 1 at 27 min: $2.63 / $2.44. Decision 1 (a)
  then reads V1 + V3 + reserve = $2.29 + $2.84 + $0.25 = $5.38 at 27.5 s.
- Add a launch rule to 6.4: before the V3 launch, the approved amount must be at least the sum of every item that
  has a launch line, plus $0.35. Otherwise drop items 1 and 3 first (both skippable), so the B3 run keeps its 260
  steps.
- Use one B-S2 definition: "items 0 and 4 (+ 2a, 2b, 5 when approved)", the same in the 6.3 bullet and in B.8.1.

### 3. minor: diagram 2.4 says the B1 rows root has 9,600 rows; section 0 and B1.2 say 8,320

Sources: 00_front.
00_front "2.4" diagram: "-> B1 rows root (9,600 rows)". 00_front section 0 "rows root": "B1: 8,320 rows = 260 steps
x 32". B B1.2: "the B1 rows root holds 8,320 rows. The 9,600-window sizes and times in B1 are the 300-step upper
bound." Fix: change the diagram to 8,320 rows.

### 4. minor (architecture): the doc gates C and D box code by two different rules, and D-R1's rollback cannot be followed on the time-line order

Sources: 00_front, D.
- 00_front "End state under (b')": "Code whose adoption gate needs V2 or V3 therefore never lands on `main`". Yet
  the same paragraph lands "D1-D9 (D's box gates stay 'open until a rental')". D's GPU gates D-G3 (box half), D-G5,
  D-G6 and D-G7 all run in R-D2 = V3 item 1 (D.6). So under (b'), `main` ends with a work-optimal CUDA polygonizer,
  device certify, warp totals and device EDT that have never run on a GPU. C1 is held back for the same kind of
  reason.
- D "D.7" D-R1: "Rollback: D2 is not merged; `check_masks` and numpy stay." D.4: "D2-D8 ... do not wait for R-D1",
  and each merge "archives the old code in the same commit". So when the first GPU run of K1 happens (R-D2, V3), D2
  is already merged and the old code is archived.
- Why this is minor: a V3 B3 run is guarded. `b3_v5.yaml` step 0 repeats `b1_rows.yaml` step 0, which runs
  `tests/test_gpu_scoring.py` on the CUDA build [code: policy.py:451-454 does the same for scoring chains]. Step 1
  runs G08b. A scorer mismatch stops the chain before training.
- Fix: D-R1 rollback reads "reinstall the `pre-d2` engine (D-G7's old tree) and revert the D2-D7 commits". Add a
  front risk FR16: "under (b') the D CUDA paths are merged with no GPU run; the first rental after D runs R-D2
  before any training item". Reword the (b') sentence to "code whose Mac gate cannot pass without a box".

### 5. minor (top-lab sequencing): under (b') the branch work is written on pre-E paths and then must be ported through E

Sources: 00_front, E.
00_front time line steps 5 and 7: C1-C5 and B2/B3 code "waits on the owner branch `after-raise`", then "E moves
files after the last commit on `main`". "When a raise comes ... If E6 has also landed, the merge is more than a path
map: E R-4 'After E6' ports the branch plans and launch configs to `specs/` and adds the post-E0 flags as stage
fields." So under the default, the plan writes about 25 tasks (C0-C11, B2.2-B2.4, B3.1-B3.10) on paths and plan
formats that E then replaces, and the port is paid later. A top lab would order it so the branch work starts on the
post-E layout: under (b'), run E (and CLEAN waves 2-3) right after D and Gate B1. Then write C and B3 only when a
raise is in view, directly on `specs/` and the layer paths. Fix: add this as the (b') order in time line step 5, or
as owner decision 1's consequence: "under (b'), C/B2.2-B3 branch work starts after E".

### 6. minor: the doc's size works against its own "one source" rule

Sources: 00_front, all phases.
The plan is 16,352 lines for a $3 spend (one rental). Section 6 says it is the "only place" for rental numbers, but
item 2 above shows that three copies of the B-S2 definition still drift. A top lab would generate the rental tables
(minutes, dollars, caps, `BudgetStop`) from one small YAML with a script and a test, and keep the prose short. Fix
(nice to have): move the 6.2-6.3 numbers into `docs/reference/rentals.yaml`. A DOC4 check then recomputes every
total and every `BudgetStop`, and fails on a mismatch.

## Counts

blocking 0, major 2, minor 4.
