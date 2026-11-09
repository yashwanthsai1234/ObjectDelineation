# Final review round 3: QA lead and technical writer (ASD-STE100)

Target: `docs/implementation2.0/Implementation.md` (14,506 lines, 1,480,096 bytes, round 3), the code it names, and the
earlier reviews `final_r1_qa.md` and `final_r2_qa.md`. I checked each round-2 QA item first. Items that round 3 fixed
are not repeated. No REPO file was edited. No git write. No rental. No Python read `data/`.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No new catalog lookup. The plan names NVIDIA skills only as recommendations. The NVIDIA toolchain pins (nsys .deb, `cuda-toolkit 13.0.3`, sanitizer 13.4.92) were web-checked in round 2 and did not change. No item. |
| engineering:code-review | I read the code behind each launch claim: `experiments/phase2/preflight.py:49-56` (`load` builds `train_sft.RunArgs(**raw)`), `train_sft.py:77-101` (no `steps` field), `vast/policy.py:308-329` (launch reads the stamp of `--config`), `policy.py:80` (`PULL_EXCLUDE`), `policy.py:511-524` (`push_hf` reads `output_dir` from a RunArgs config), `pyproject.toml:65-71`, `Makefile:19-28`. Items 1, 2, 7, 9. |
| engineering:debug | For each gate I asked: which input does it read, which rental item makes it, and does that item run in the base case (27.5 s, where V1 item 5 is skipped)? Items 3, 4, 5, 6, 8. |
| engineering:system-design | I put V1, V2, V3 and the tag `v1-b1` on one line with the code each gate needs. Items 4, 5, 7, 8. |
| engineering:architecture | ADR numbers now match the front table in B, C, D, E and DOC (grep of every `ADR[- ]00NN`). No item. |
| engineering:documentation | One word per thing, one number per fact, cross-references between phase tables. Items 13, 14, 16, 17, 18, 19. |
| huggingface-skills:huggingface-best | Model ids and sizes in sections 3 and 99_back still agree. The HF part of this review is the adapter transport to V3 (item 7). |
| engineering:testing-strategy | I matched each Tests section to its rental and its runner: which command collects each `gpu` test, and on which box. Items 9, 10, 12. |
| mattpocock-skills:tdd | I looked for tests that cannot fail and gates that fail by construction. Round 2's either-outcome tests are fixed (U12, U32, U-CF2). New: DOC-G13 fails by construction (item 11). |
| Web | No new external claim needed a web check in round 3; the round-2 checks (PyPI, CUDA Debian repo) cover the pins that the doc still cites. |

Measured during this review (2026-11-08):
- `docs/implementation2.0/Implementation.md` is 1,480,096 bytes (`ls -l`).
- Prose sentences over 25 words, my script (join wrapped lines, list items split apart, code and tables skipped):
  1,165 of 7,034. This rule splits list items, so it is not comparable with the round-2 count (993 of 5,203).
- Round-2 QA items checked as fixed: 1 (G1a, nsys in V2), 3 (one interval rule, no BCa), 4 (ADR numbers), 5 (gate
  prefix), 6 (`view_content_sha256`), 7 (`cov`, `test-all`), 8 (`box` marker), 9 (either-outcome tests), 10 (G-SP1
  band), 11 (G06 3 % clause), 12 (`S[P]/` paths), 13 (D order), 14 (V4), 15 (S06, S08 run conditions), 16 (V1 configs
  named; but see item 1), 17-24, 26.

Severity: blocking = wrong, unsafe or breaks an owner rule; major = missing or vague, so an engineer cannot
implement it; minor = clarity or consistency.
Count: **5 blocking, 6 major, 10 minor**.

---

## Blocking

### 1. [blocking] [00_front, B] Preflight cannot stamp `b1_rows.yaml` or `a3_measure.yaml`, so V1 can never be approved

Problem.
- Front order rule 2 and 6.1 step 4 require `make preflight CONFIG=experiments/phase2/configs/b1_rows.yaml` and
  `make preflight CONFIG=experiments/phase2/configs/a3_measure.yaml` on the tag, with "both stamps have no failed
  check".
- Front 6.2 says "items 3-4 are steps of that config (B must add them)". So `b1_rows.yaml` is a `steps:` plan, like
  `rl_shifted.yaml`, which it replaces (B1.5). `a3_measure.yaml` is a `steps:` plan (A3 code block).
- `preflight.load` builds `train_sft.RunArgs(**raw)` [code: experiments/phase2/preflight.py:49-56]. `RunArgs` has
  no `steps` field [code: experiments/phase2/train_sft.py:77-101]. CLEAN9 step 3 measured the result:
  `TypeError: RunArgs.__init__() got an unexpected keyword argument 'steps'`.
- A3 itself says "a `steps:` plan file cannot get a stamp" and "the preflight stamp belongs to the launch config".
  So A contradicts the front.
- CLEAN9 step 3 flagged this for the front and B in round 2. Neither changed. B1.3 says "call `quality.failures` on
  every rows root named in the run config", but no task says how preflight reads a steps plan.
- Result: the 6.1 step 4 pass rule cannot pass. `launch.sh` refuses any config without a passing stamp
  [code: experiments/phase2/vast/policy.py:308-316].

Exact fix (pick one and write it in front order rule 2, 6.1 step 4, 6.2, B1.3, B1.5 and A3):
- (a) Keep today's mechanism. Launch V1 with an SFT-style launch config, `experiments/phase2/configs/r1_qwen3vl2b.yaml`
  (as R5 and `rl_shifted.yaml:17-19` do). Run `make preflight` for that file only. Delete the `a3_measure.yaml`
  preflight from order rule 2 and 6.1 step 4. B1.5: `b1_rows.yaml` stays a steps plan. B1.3: add a new preflight
  check "ROWS root" that reads the rows-root path and `rows_manifest_sha256` from a new launch-config key
  `rows_plan: experiments/phase2/configs/b1_rows.yaml`, then runs `quality.failures` on it. Add a test: preflight on
  the launch config fails when the planted rows root fails DG1b.
- (b) Teach `preflight.load` to accept a steps plan: when the YAML has `steps:`, read the `train_rl` stage for the
  rows root and skip the SFT dry run. Then add a unit test: `preflight.load(b1_rows.yaml)` returns without
  `TypeError`. Assign it to B1.3 by name.

### 2. [blocking] [A, B] Phases A and B still plan a different V1 than the front, and both versions exceed the V1 cap

Problem.
- The front 6.2 table (the owner of V1 minutes) has 5 items: B1 run; B1 on val shifted and val edge; R5 on val
  edge (a step of `b1_rows.yaml`); A3-val = R5 on val shifted only (`a3_measure.yaml` "reduced to val (A must
  follow)"). No test step. No nsys window. Total $2.08 at 27.5 s, cap $2.20.
- A did not follow:
  - `a3_measure.yaml` still has three `final_eval` steps: val shifted, val edge (R5) and test shifted.
  - A3 title, A3 step 3, A.7 row "V1 runs short of money before A3's test step" and A.8 ("Test step (last; may move
    to V3): 634 s"; "The front's V1 table must list these three steps") still plan the test step in V1.
  - A-G3 and G01 expect R5 val-edge answers from A3, but the front puts that eval in `b1_rows.yaml` (item 4).
  - So R5 on val edge runs twice (item 4 and inside A3), and the test eval (10.6 min) runs in V1.
- B did not follow:
  - B.8.1 says "B-S1 is the front's rental V1 ... must list the same items", but its table has 6 items. Item 6 adds
    "the A6 nsys window (about 10.5)". The front says "V1 has no nsys capture".
  - B.8.1 total: 205.6 min = $2.20. Front total: 195.1 min = $2.08. B says "raise the cap to $2.35". The front and
    decision 2 say $2.30.
  - B.8.1 item 6 skip rule: "only if the ledger plus its cost stays within the cap". The ledger total is already
    $16.90, so the rule never lets item 6 run. The front rule uses session spend above $1.98.
  - B.6.2 calls the R5 edge eval "B-S1 item 5"; the front calls it item 4.
- With A's plan file at 27.5 s: 195.1 + 6.9 (second R5 val edge) + 10.6 (test) = 212.6 min = $2.27 at $0.641/h
  [estimate: front 6 inputs], above the $2.20 cap. The cap then kills the box before the plan ends.

Exact fix.
- A3 `a3_measure.yaml`: keep only the `env_check` step and the val-shifted `final_eval` step. Delete the val-edge and
  test steps. Title: "Save every answer; re-evaluate R5 on val (V1 item 5); test half in V3".
- A3 step 3, A.7 row 6, A.8 GPU list, A-G3, G01: "R5 val-edge answers come from `b1_rows.yaml` (front V1 item 4,
  output `eval_r5_val_edge/`). The test half is V3 item 3." Delete "test step is last in `a3_measure.yaml`".
- A.8: "A3 in V1: 6.9 min ($0.07). Test half in V3: 10.6 min ($0.11)."
- B.8.1: replace the B-S1 table by a copy of front 6.2 (5 rows, 195.1 / 175.9 min, item 5 skip rule on session
  spend above $1.98). Delete the nsys window and "$2.35". B.6.2: "R5 on the edge view is V1 item 4".
- B1.5: list the `b1_rows.yaml` steps: `train_rl` (300 steps), `final_eval` B1 final on `lidc-val-shifted`, on
  `lidc-val-edge`, and R5 (`{rd}/ckpt`) on `lidc-val-edge`, with their `--out` folders.

### 3. [blocking] [00_front, B] Gate B1 has two different primary baselines, and its hard view-identity rule has no input in the base case

Problem.
- Front 6.2 note: "The main rule pairs item 3 against R5's `runs/r5/eval/eval_val/items.jsonl`; item 5, when it
  runs, is a same-view R5 baseline."
- B.6.2: "The baseline of Gate B1 is A3's R5 val items ... R5's `runs/r5/eval/eval_val/items.jsonl` is a second
  line." So one gate has two decision baselines. B1 can pass against one and fail against the other.
- B.6.2 and U49 make the view-identity rule hard for Gate B1. It needs A3 answers (`r_h` from answers 0-3 and 4-7).
- Front 6.2: at 27.5 s (the base case) item 5 (A3-val) is skipped. Then the hard rule has no input. B gives a
  fallback for the baseline (hash check), but none for the view-identity rule.
- The fallback hash check compares "the B1.0 pre-rebuild content hash" with "the frozen `lidc-val-shifted` entry".
  DG15 freezes that entry from the same B1.0 hash, so the check is always true. It does not compare against the view
  that R5 scored.

Exact fix.
- Pick one primary baseline and write it in front 6.2 and B.6.2:
  - Default: R5's `runs/r5/eval/eval_val/items.jsonl` is the primary baseline (it always exists). A3's items, when
    item 5 ran, are the second line.
  - View identity when item 5 is skipped: compare the frozen `lidc-val-shifted` `target_set_sha256` with the R5
    bundle manifest of the same view (the preflight `EXAMS` hash covers `lidc-val-shifted` [code: preflight.py:28]).
    Mark the rule "hash only, correlation not run" in the report.
- U49: add a case "item 5 skipped": the gate decides on the hash rule, and the report line says so.
- Or the owner raises the V1 cap to $2.30 (decision 2), and item 5 becomes "not skippable".

### 4. [blocking] [00_front, B] The tag `v1-b1` does not hold B2.1, but Gate B2 and B2.4 read B2.1 stats from the B1 run

Problem.
- Order rule 2: `v1-b1` = be52b1c + A3, A6, B1 (B1.5) + CLEAN wave 1. "B1 does not wait for any B2 work".
- B2 intro: "B2 measures it inside the B1 run". B2.1 adds `is_summary` and the `cf/*` counterfactual stats to
  `fast.py`.
- Gate B2: "B1 run: ... the counterfactual stats are present for steps 50-300". B2.4 rules 1-4 "use the B1 run's
  counterfactual stats over steps 50-300". FR4 says the same.
- If B2.1 is not on the tag, the B1 run logs no `cf/*` keys. Gate B2 then has no input, and the B3 preset can only be
  chosen with a second paid run.
- The tag list also names A3 and A6 only. A3 writes schema-2 `final_eval.json` (A1 keys) and `items.jsonl` with the
  A5 keys, and preflight runs the A7 check ("CODE TRL private names", A-G7). The list does not say A1, A5 and A7 are on
  the tag.

Exact fix.
- Order rule 2: "`v1-b1` = be52b1c + every Mac task of A (A1-A8) + B1.0-B1.9 + B2.0 + B2.1 (logging only; the active
  preset stays bf16-seqmask) + CLEAN wave 1. No C, D, B3 or kernel code."
- Replace "B1 does not wait for any B2 work" with "B1 waits for B2.1 logging only, not for a B2 decision".
- Add to the order-rule-2 test list: "`tests/test_train_fast_numbers.py` (U35, P13) passes, and the E08 dry run log
  has every `cf/*` key once per logging step (R01)".

### 5. [blocking] [00_front, A, C] When V1 skips A3-val, Gate C-adopt deadlocks: V2 needs A3 answers that V2 itself makes

Problem.
- C.6 and C.8.2: "The C items enter V2 only after every C1-C5 Mac test passes and the Mac half of G-EQ2/G-EQ3 passes
  on A3's answers from V1". Front time line step 4 says the same.
- Front 6.2: at 27.5 s, V1 item 5 (A3-val) is skipped and "moves to V2 item 3b".
- V2 item 3b comes after items 2-3 (the C items). So in the base case the C items need A3 answers that only V2 makes.
  The C items can then never enter V2, or V2 runs them against its own rule.

Exact fix (pick one; write it in front 6.3, C.6, C.8.2, A3 step 6):
- (a) Gate C-adopt matters only for the B runs in V3. Change the rule to: "The C items enter V2 after the C1-C5 Mac
  tests pass. The Mac half of G-EQ2/G-EQ3 runs on A3's val answers after V1 (or after V2 if V1 skipped item 5),
  and must pass before V3." Keep item 3b first in V2 (before item 2).
- (b) Make V1 item 5 "not skippable" with the $2.30 cap (owner decision 2).

---

## Major

### 6. [major] [00_front, B] `BudgetStop` measures from the training start, but the plan uses it as a session-spend stop

Problem.
- Front 6.2 item 2 and decision 2: "`BudgetStop` at session spend $1.90". Session spend starts at the launch.
- B1.5: `callbacks.BudgetStop(rate_usd_per_h, cap_usd, start_time)` "stops training ... when elapsed time x rate
  reaches the cap". `train_rl --budget-usd X --rate-usd-per-h Y`. No text says what `start_time` is. The box has no
  ledger, so `train_rl` cannot read `policy.spend()`.
- If `start_time` is the process start, the stop fires at $1.90 + item 1 (30 min = $0.32) = $2.22 session spend.
  Items 3-4 ($0.22) then push V1 to $2.44, above the $2.20 cap. The cap deadline kills the Gate B1 evals.
- B.8.1 computes "$1.90 buys 177.8 min, and items 1-2 need 167.5", which assumes the launch time.

Exact fix.
- B1.5: "`start_time` = the launch time, read from `/workspace/runs/<run>/LAUNCHED_AT` (written by `cmd_launch` or
  `cmd_sync_up`). `--budget-usd` = the session-spend stop ($1.90 for V1)." Or: "`--budget-usd` = $1.90 minus $0.32
  for item 1 = $1.58, measured from the `train_rl` start."
- U46: add a case where `start_time` is 30 min before the process start: the stop fires 30 min earlier.

### 7. [major] [00_front, B] No step carries the B1 adapters from V1 to V3, and B2.3b names other adapters than the front

Problem.
- Front 6.2 item 2: "save the step-100 and final adapters (B2.3b reads them in V3)".
- V1 ends with teardown. `pull.sh` copies `runs/<id>/` to the Mac except `rolling/`, `decay/`, `stable/`
  [code: policy.py:80, 499-508]. Nothing puts the adapters on the V3 box. `fetch_hf.sh` reads a revision of the HF
  repo; `push_hf.sh` reads `output_dir` from a RunArgs config [code: policy.py:511-524], which a steps plan has not.
- B2.3b uses `checkpoint-100` and `checkpoint-770` of the full run for the Mac part, and "the B1 final adapter" for
  the box part. G04 also says "B1 final adapter". The step-100 B1 adapter of the front is not used by B2.3b.

Exact fix.
- Front 6.2 item 2 and B1.5: "After the B1 run, a plan step runs `hf_push --run {rd}/train_rl --repo
  svenk029/objdel-sft --revision b1-rows` (adapters of steps 100 and final). V3 starts with `fetch_hf.sh <run>
  b1-rows`." Add a unit test that the plan step exists (B I22).
- Make front 6.2 and B2.3b name the same adapters: either "B1 final only" or "B1 step 100 and final".

### 8. [major] [B] B2.4 is "after B-S1, $0", but its rules read V3 arms, and its preset deletion is a code edit in the middle of V3

Problem.
- B2.4 header: "Decide and clean (after B-S1, $0)".
- Rule 0 needs "its 15-step arm in B-S2" (F arm). Rule 1 needs the fp16 arm (B2.3, V3 item 2). B-S2 is V3.
- V3 item 4 is the B3 run, which needs the preset that B2.4 picks. So the decision falls inside V3.
- B2.4 also says "Make the chosen preset the `grpo_config` default. Delete the other presets and their code paths".
  A code edit inside V3 makes the preflight stamp stale (`objdel_and_script_sha256`) and needs a new sync.

Exact fix.
- B2.4 header: "Decide (rule 0 and rule 1 inside V3, after items 2-3; rules 2-4 after V1, $0). Clean after V3."
- Write: "Inside V3 the B3 run picks the preset with `--numerics <name>` (no code edit). The preset deletion and the
  `grpo_config` default change land in the first commit after V3."
- U51: add a case "rule 0 and rule 1 inputs absent (arms skipped)": the rule returns E (or bf16-seqmask by veto).

### 9. [major] [C, D, A, B, E] `pytest -m gpu` in R-D2 collects every phase's box tests, so D-G3 cannot pass as written

Problem.
- One marker `gpu` for every phase (front section 4, C.9.0). D-G3: "On the RTX 5090: `pytest -m gpu` passes". D8:
  `make gpu-check` runs `pytest -m gpu`. R-D2 is V3 item 1, priced "`pytest -m gpu` about 5" min.
- By V3, C and B code has landed. `pytest -m gpu` then also collects:
  - C G-PX1 (2,000 views incl. `lidc-train-shifted`, the box HF processor), G-EN1 (real Qwen3-VL-2B), G-SY1 (vLLM
    engine, 3 pushes), G-RP2 (recorded batches from V2);
  - A G06 (reads `nsys_baseline.json` from V2), G07 (two vLLM samplings);
  - B G03 (whole Qwen3-VL-2B in fp16), G04 (`b2_attribution` box part), G05 (vLLM engine).
- These need inputs that R-D2 does not ship, and they take much more than 5 min. A failure there fails D-G3 for a
  reason outside D.

Exact fix.
- D-G3 and D8: "`make gpu-check` runs `pytest -m gpu` on D's file list only: `tests/test_cuda_checker_gpu.py`,
  `tests/test_gpu_scoring.py`, `tests/test_gpu_determinism.py`, `tests/test_gpu_sanitize.py`, `test_rings.py`,
  `test_certify_polygons.py`, `test_topology_device.py`, `test_tables_device.py`."
- Each box test of A, B, C and E names its rental item and its command, for example "C G-PX1: V2 item 2,
  `c_speed_check --pixels`". A box test that needs a run artifact skips with a reason when the artifact is absent.
- Front section 4 Tests row: add "`pytest -m gpu` without a file list is never run on a box".

### 10. [major] [A] `tests/test_a_gpu.py` has no per-test markers, and G06 has no runner

Problem.
- A.9.0: `tests/test_a_gpu.py` holds G01-G08 with "markers `gpu` and `data`". A.9.6 gives a marker for G01 (`data`),
  G06 (`gpu`) and G07 (`gpu`) only. G02-G05 and G08 have none.
- G01 (A-G3 val part) runs on the Mac on pulled answers. If the file uses a module-level `pytestmark` with both
  markers and a no-CUDA skip, G01 never runs on the Mac, and A-G3 has no test.
- G06 (A-G6 box half) reads `nsys_baseline.json` of the C0 baseline arm. No V2 plan step (`c_speed.yaml`) runs
  pytest. After `pull.sh` the file is on the Mac, where a `gpu` test skips. So G06 runs nowhere.

Exact fix.
- A.9.6: one marker line per test: G01 `data` (Mac, after `pull.sh`); G02 `data` for the rescore half, `gpu` for the
  "both builds" half (split into G02a/G02b); G03 `gpu`; G04 `data`; G05 `data`; G06 `data` (Mac, reads the pulled
  `nsys_baseline.json`; the 3-step toy capture becomes G06b, `gpu`, run by C0 in V2 item 3); G07 `gpu`; G08 `gpu`.
- A.9.0: "no module-level `pytestmark` in `tests/test_a_gpu.py`".
- C0: add a `c_speed.yaml` step that runs G06b, or drop the 3-step toy from G06.

### 11. [major] [DOC] DOC-G13 fails by construction: the plan itself is over 1 MB

Problem.
- DOC-G13: "no tracked file over 1 MB under `docs/implementation2.0/`".
- DOC-G1 requires `docs/implementation2.0/Implementation.md` to be tracked.
- The file is 1,480,096 bytes [measured: `ls -l`, 2026-11-08]. It is over 1 MB with either meaning of "MB".
  So DOC-G13 cannot pass, and the DOC phase gate cannot pass.

Exact fix.
- DOC-G13: "no tracked file over 1 MB under `docs/implementation2.0/` except `Implementation.md`" and add a separate
  size line for it, or
- DOC1: split `Implementation.md` into one file per phase (`docs/implementation2.0/{00_front,A,B,C,D,E,DOC,CLEAN,
  99_back}.md`), each under 1 MB. Then DOC2's site navigation and DOC4 check 14 read the folder.

---

## Minor

### 12. [minor] [C, A, D, E] Three different GPU skip rules, and C says A adds one that A does not add

C.9.0 says "A registers `gpu` with `--strict-markers` and a `conftest.py` skip when `not torch.cuda.is_available()`
(A.9)". A.9 registers the marker but writes no `conftest.py` skip. D8 uses `skipif(_gpus() == 0)` per file. E uses
`pytest.skip(...)` when `geom.engine.build() != "cuda"`. GP1 uses `torch.cuda.is_available()`.
Fix: A.9 adds one root `tests/conftest.py` hook: a `gpu` test skips unless `objdel._gpucheck.build() == "cuda"` (no
torch import, so the verifier tests keep their no-torch rule). C.9.0, D8 and E.9 point at that hook and drop their
own skip lines.

### 13. [minor] [A, B, C, 99_back] "Outline length" has two definitions, so one config check has two limits

B header: "Outline length = certified target tokens without the end token" (1,181). C header: "A doctor outline
length here counts the end token" (1,182). A says 1,182 with the end token. B3.7 rule 9 writes
`L - soft_punish_cache >= 1,181 + 32` (= 1,213). C7 step 5 and decision 13 write `1,182 + 32 = 1,214`.
Fix: front section 0 adds one row: "outline length = tokens without the end token; longest LIDC train target 1,181".
A, C and decision 13 then write "1,182 with the end token" only where the count includes it, and B3.7, C7 and
decision 13 use one limit: `L - soft_punish_cache >= 1,182 + 32` (end token counted, because the answer cap counts it).

### 14. [minor] [C] C.8.2 and C.8.3 keep stale front numbers

C.8.2: "The money left after V1 is $0.59 (V1 at 27.5 s) or $0.79 (V1 at 23.65 s) (front 6.3)". Front 6.3: $0.74 and
$0.87. C.8.2 also says "the front's own V2 table still omits the start-up rows", but front 6.3 now has them. C.8.3
row "Phase A3 saved answers": "front 6.2 item 4"; A3-val is item 5.
Fix: write $0.74 / $0.87; delete the start-up sentence; write "front 6.2 item 5 (or V2 item 3b)".

### 15. [minor] [D] D-G14 and D-G17 read the V3 test eval, which runs only if the owner closes the 4B gate

Front 6.3 V3 item 5: "If the gate fires, test waits for the 4B session". D-G14 ("Same box, test eval") and D-G17
("Box, V3 test eval") then have no input. Tying a D gate to the one test read also mixes an engineering check with
the final result.
Fix: D-G14 and D-G17 read the V3 val evals of the B3 run (`final_eval` on `lidc-val-shifted`), not the test eval.

### 16. [minor] [00_front, B, E, 99_back] The same rental items have different minutes and dollars

- V4: front 6.3 and decision 3b: 39 min = $0.42. E.8 and E-Q7: 26-40 min, $0.28-0.43.
- B-S2: decision 1: $1.98 at 27.5 s. B.8.1: $1.94.
- B2.3b box part: front 6.2: 15 min. B.8.1 and B2.3b: 10 min.
Fix: the front values win. E.8: "about 39 min = $0.42 (front 6.3)". Decision 1: $1.94 (B.8.1), or B.8.1 re-prices
and both write one number. B2.3b and B.8.1: 15 min, or the front writes 10 min.

### 17. [minor] [B] "Ledger" still stands for session spend in Gate B1

Gate B1 report line: "rescaled and rotated ... (when the ledger allows, B.8.1)". The ledger total is $16.90, so the
phrase has no meaning. (B.8.1's skip rule is in item 2.)
Fix: "when session spend allows (front section 0)".

### 18. [minor] [B] Numerics option letters A-F reuse the phase letters

B2.4: "The B1 run used C, so Gate B1 has no numerics confound"; "fall back to C"; "Else A (`fp16-seqmask`)". In the
same doc, "C" and "A" are phases. A reader can read "The B1 run used C" as "used Phase C code", which V1 forbids.
Fix: name the options by preset: `bf16-seqmask` (was C), `fp16-seqmask` (A), `bf16-tis-opsm` (E), `lora-seqmask`
(F). Keep the letters only in the ADR-0013 options table.

### 19. [minor] [B, A, 00_front] Heading order does not match section 8

B.9 runs 9.11, 9.13, 9.14, then 9.12 (lines 5154, 5185, 5208, 5225). A has review-log parts 2 and 3 after A.9, but
section 8 lists the review log before part 9.
Fix: move B.9.12 before B.9.13. Section 8: "Review logs of the final rounds may follow part 9."

### 20. [minor] [00_front] The 6.1 step 0 unit test does not use the launch rule it guards

Step 0 test: "`spend()` + V1 cap <= `HARD_STOP`". `cmd_launch` refuses when
`spend(committed=True) + accrued(run, deadline) > HARD_STOP`, where the deadline adds `PULL_S` = 600 s to the cap
[code: experiments/phase2/vast/policy.py:73, 316-326].
Fix: the test calls the same expression as `cmd_launch` (factor it into `policy.launch_refusal(offer, cap_hours)`),
on a copy of the ledger, with the V1 offer price and cap.

### 21. [minor] [all phases] Long sentences remain

My count (rule above): 1,165 of 7,034 prose sentences have more than 25 words. Most are in the review logs, the gate
cells and the B and C task steps. Examples: G-SP4 (C.6) and the A3 `code_sha` step.
Fix: split the gate cells and task steps first: one rule per sentence, the reason on its own "Why:" line. DOC.6
reports the count with the rule written beside it, so the next count is comparable.
