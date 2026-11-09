# Final review r3: fresh outside architect (whole Implementation.md, round 3, 14,506 lines)

Stance: no earlier opinion. I read `final_r1_*.md` and `final_r2_*.md` (all six roles) and checked the round-2
fresh-arch items against the round-3 text. These are fixed and are not raised again:
- the budget guard (front 6.1 step 0, FR15, decision 4);
- the V1 item list in the front;
- G1a in A (9 hits remain, all in review logs);
- C's rental mapping and rate;
- the ADR renumbering and ADR 0022;
- the RL-trainer row;
- the audit-callback path;
- the `v1-b1` tag;
- the 2.5 labels.

Every item below is new, or is a round-2 fix that a phase section did not follow.

## Skills used

| Skill | What it changed in this review |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | The plan's NVIDIA parts (CUDA 13 wheels, compute-sanitizer 13.4.92, NVTX/nsys, CT skills) were checked against the catalog rows that round 1 saved. No new catalog match changes an item. Nothing was installed. |
| `engineering:code-review` | Checked these citations against the code: `fast.py:540-547` (`skip_old_logps`), `fast.py:696-717` (merged push), `final_eval.py:60-95` (one `--ckpt`, one `--view`, full trainer per call), `policy.py:196-200` (`spend()` counts live runs), `policy.py:324-329`, `policy.chain` (no spend-based skip), `gpucheck.cpp:58-200` (own streams, raw `cudaMalloc`), `callbacks.py:1-11` (held-out eval lives in `Controller`, which CLEAN6 archives), and TRL 1.17.1 `grpo_config.py:489` (`num_generations_eval`). All match the doc, except the claims in items 1 and 2. |
| `engineering:debug` | Traced the V1 money guard from `--budget-usd` to `BudgetStop(rate, cap, start_time)`. Its clock starts when the training process starts, not at launch (item 2). |
| `engineering:system-design` | Walked the decision path after V1. No branch exists for "Gate B1 fails" (item 4). The measured value 0.531 already fires B2.4 rule 0, so the weight path that C measures is probably not the one B3 runs (item 3). |
| `engineering:architecture` | Compared the end-state diagram (2.5) and ADR 0022 with ADR-0013 rule 0. The plan holds two opposite end states for the weight path (item 3). Checked 3D readiness of the long-lived E interfaces (item 7). |
| `engineering:documentation` | Compared the B-S1 table (B.8.1) row by row with the V1 table (front 6.2). They differ in 5 places, one of them an unsafe rule (item 1). |
| `huggingface-skills:huggingface-best` | PyPI and the HF API, 2026-11-08: TRL 1.17.1 is the latest release. vLLM 0.33.0 is out (the plan pins 0.32.0; no item, because the pin is deliberate). `vastai-sdk` 1.11.0 says "DEPRECATED, use `pip install vastai`". import-linter 2.16 exists. Section 3 (model choice) stands. |
| `system-design-skills:system-design` | Stress-tested V1 against slow steps (FR3). Above about 30.7 s per step, the cap ends the box before the Gate B1 evals (item 2). Did a bottom-up estimate of eval cost: about 205 s of fixed cost per `final_eval` process (item 5). |
| `engineering:tech-debt` | Looked for round-2 fixes that a phase did not carry over. B.8.1 still holds the round-2 draft of V1 (item 1). CLEAN6 removes the only held-out-curve path with no replacement (item 6). |

Web checks:
- PyPI JSON for trl, vllm, vastai, vastai-sdk, import-linter, cuda-toolkit and nvidia-cuda-sanitizer-api (versions
  above; `cuda-toolkit` 13.0.3 and the sanitizer API 13.4.92 both exist).
- TRL GRPO docs: `GRPOTrainer(eval_dataset=...)`, `eval_strategy`, `eval_steps`, and
  `num_generations_eval` ("fewer generations during evaluation to save computation")
  [web: huggingface.co/docs/trl/main/en/grpo_trainer]. The installed 1.17.1 has the same field
  [code: .venv trl/trainer/grpo_config.py:489, grpo_trainer.py:1137].

## Verdict

The target architecture is sound. A top lab would keep its main parts:
- data first;
- one CUDA engine source with a CPU oracle;
- a verifier with no torch;
- layers checked by import-linter;
- one TRL adapter with a start-up check;
- move-only refactoring checked by goldens;
- every CPU item listed with its reason.

"GPU everywhere" is honest: the CPU items left are PNG decode, the vLLM lane, the tokenizer and text parsing, and
each one has a measured reason.

A top lab would do four things differently:
1. Choose numerics before speed. The plan's own measured merge-rounding value already picks vLLM LoRA serving (F)
   for B3. C still spends its V2 money measuring the merged push (item 3).
2. Write the decision tree for the one paid experiment: what happens when Gate B1 fails (item 4).
3. Make the money guards measure the right interval on the box (item 2). Make the V1 tables agree (item 1).
4. Pay the model-load cost once per checkpoint, not once per view. That alone lets A3-val fit V1 at 27.5 s (item 5).

## Items

### 1. BLOCKING. B (B.8.1, B.6.2), 00_front (6.2): B-S1 still holds the round-2 draft of V1; its skip rule reads the ledger total and would always skip A3-val

Problem. Front 6.2 says "B.8.1 must make B-S1 exactly this table". The two tables differ in five places:

| Field | Front 6.2 (V1) | B.8.1 (B-S1) |
|---|---|---|
| A6 nsys window in V1 | none ("V1 has no nsys capture"; A.8: "adds no V1 minutes") | item 6 = A3-val 6.9 min + "the A6 nsys window (about 10.5)" = 17.4 min |
| Totals | 195.1 min = $2.08 (27.5 s); 175.9 min = $1.88 | 205.6 min = $2.20; 186.4 min = $1.99 |
| Skip rule for A3-val | "session spend above $1.98 at its start" | "each part runs only if **the ledger** plus its cost stays within the cap" |
| Cap raise option | $2.30 | $2.35 |
| Item numbers | R5 on edge = item 4, A3-val = item 5 | R5 on edge = item 5, A3-val = item 6. B.6.2 cites "B-S1 item 5" for the edge baseline |

The B skip rule is wrong as written. The ledger total is already $16.90 [measured: ledger.csv, 133 rows]. So "ledger
+ cost <= $2.20" is false at every point, and an engineer who follows B skips A3-val in every case. A3-val produces:
- the Gate B1 paired baseline (B.6.2 "the baseline of Gate B1 is A3's R5 val items");
- the Gate A input;
- the Mac half of C-G-EQ2/3;
- the C11 sizing.

Section 0 defined "session spend" for exactly this reason. B did not carry the definition over.

Exact fix.
- B.8.1 B-S1 table: copy front 6.2 items 1-5 word for word.
- Delete the nsys part of item 6 and the "10.5 min" row. Delete the line "The A6 nsys window about 10.5 min
  [estimate: A.8]" from the B.8.1 inputs.
- Totals: 195.1 / 175.9 min.
- Skip rule: "A3-val is skipped when session spend is above $1.98 at its start (front 6.4 rule 2)".
- Cap raise: $2.30.
- Replace "If the owner keeps $2.20, item 6 runs only near 23.65 s; else raise the cap to $2.35" with the front
  sentence (26.9 s threshold, $2.30).
- B.6.2: "(R5 on the edge view is B-S1 item 4)".
- Add a DOC4 drift check: the item names and minutes of the B-S1 table equal the front V1 table.

### 2. MAJOR. B (B1.5 `BudgetStop`, U45, U46), 00_front (6.2 item 2, 6.4 rules 2 and 5): the money guards measure the wrong interval on the box, so a slow B1 run loses every Gate B1 eval

Problem.
- The front's stop rules use session spend: `policy.spend()` now minus at launch. `spend()` counts live runs, so
  the Mac can read it [code: policy.py:196-200].
- The box guard is `callbacks.BudgetStop(rate_usd_per_h, cap_usd, start_time)`. It "stops training and saves when
  elapsed time x rate reaches the cap". `train_rl` gets only `--budget-usd` and `--rate-usd-per-h` (B1.5; U45). No
  flag passes the launch time, so `start_time` is the start of the training process.
- Launch to first sync took 17-31 min [measured: ledger.csv]. The rows-root download also comes before
  `train_rl` starts. So BudgetStop's "$1.90" is about $0.21-0.32 later than the front's "session spend $1.90"
  [estimate: 20-30 min at $0.641/h].
- Failure case. `--cap-hours` = floor_to_0.05(2.20 / 0.641) = 3.40 h = 204 min. The Gate B1 evals (items 3-4)
  need 20.7 min after training. With session overhead 30 min, the evals do not fit once 30 + 5 s + 20.7 > 204,
  that is s > 30.7 s per step [estimate: arithmetic].
  - With a session-clock BudgetStop at $1.90 (177.8 min), training stops early and saves, and the evals still run.
  - With the process clock, BudgetStop trips only at about 30 + 177.8 = 208 session minutes. The cap destroys the
    box first. The last checkpoint is up to 49 steps old (`save_steps 50`). Items 3-4 never run, so Gate B1 cannot
    be decided.
  FR3 rates "shifted rows lengthen the decode tail" Medium likelihood.
- The same gap applies to the A3-val skip rule. `policy.chain` has no step that is skipped by spend [code:
  policy.py `chain`, "~" only lets a step fail]. So the rule can only be applied by the main agent on the Mac,
  between two `run.sh` calls. The plan does not say so.

Exact fix.
- B1.5: `BudgetStop(rate_usd_per_h, cap_usd, session_start)`. `train_rl` gets `--session-start <unix epoch>`
  (required with `--budget-usd`). `run.sh` exports it from the active-run record that `policy.cmd_launch` writes
  (launch time), and `--rate-usd-per-h` from that record's `dph_total`. The budget is session spend: (now -
  session_start) x rate.
- U45: `--budget-usd` without `--session-start` is an argparse error.
- U46: a case where `session_start` is 30 min before the process start trips 30 min earlier than a process-clock
  guard.
- Front 6.2 and 6.4 rule 2: "The A3-val skip rule is applied by the main agent on the Mac. It reads
  `policy.spend()` after `run.sh <run> b1_rows.yaml` ends, and starts `run.sh <run> a3_measure.yaml` only if session
  spend is <= $1.98."
- Add a unit test on a fake ledger and an active run.

### 3. MAJOR. 00_front (2.5, section 4 ADR 0022, 6.3), C (C5, C.6 Gate C-adopt, C.8.2), B (B2.4 rule 0): the plan's own measured number picks vLLM LoRA serving for B3, but the end state and the V2 speed gate are built on the merged push

Problem.
- B2.4 rule 0: if the B2.0 bf16 relative error at step 100 is >= 0.2, B3 uses option F (vLLM serves the adapter
  unmerged), if its 15-step arm passes.
- The seed measurement is already 0.531 at step 100 [measured: B.1 item 5 table,
  `docs/implementation2.0/notes/r2rl/merge_rounding.py`, re-run 2026-11-08]. So rule 0 will almost surely fire. F is
  the expected B3 numerics.
- But the end state does not show F:
  - front 2.5 shows "push merged weights W + BA to vLLM (all tensors, one call, C5)";
  - ADR 0022 is titled "... one-call weight push (LoRA serving not adopted)";
  - Gate C-adopt (G-SP1 <= 20.0 s) is measured in V2 on the merged push only.
- C5 itself estimates that vLLM LoRA kernels slow decode by 10-30 % [estimate: B_gpu_audit §4.6]. On the 11.1 s
  generate phase, that is 1.1-3.3 s per step. C1-C5 save 2.25-4.7 s [estimate: C.8.1].
- So a G-SP1 pass measured on the merged push does not carry over to the B3 configuration. The $0.95 that the
  required V2 items cost can validate a path that B3 then drops. The F decode arm is planned "in V3 before the B3
  run, only if B picks F" (C.8.2), after the C decision has been made.
- A top lab chooses numerics (correctness) first and then makes the chosen path fast.

Exact fix.
- Move the 15-step F arm (B2.4 rule 0; about 13 min with start-up = $0.14 [estimate: C.8.2]) into V2, as item 2b,
  before the C arms. It needs only B2.2's F code and the B2.0 Mac table, which are both $0 before V2.
- C.6: "The fast arm runs under the numerics that B2.4 rule 0 picks. Under F, C5 is the adapter-only push (C5
  'Under option F'), and G-SP1 is re-based on a baseline arm under F."
- Front 2.5: draw the push as "C5 one call: adapter only (F, expected per B2.0) or merged W + BA (if rule 0 does not
  fire)".
- ADR 0022 title: "GPU pixel backend, K3, one-call weight push (merged or adapter-only per ADR 0013)". Drop "LoRA
  serving not adopted" from the title. Keep "not adopted for speed" in the body.
- Re-price V2 items 1-3 with item 2b: 102 min = $1.09 [estimate].
- Decision 1 in 99_back and front 6.3 must follow.

### 4. MAJOR. 00_front (section 1, FR2), B (B.4 order, B1.8, B3.9, B.7), 99_back (decisions 1 and 17): no decision tree for the outcome of the one paid experiment

Problem. V1 is the only experiment that fits $3. The plan defines what happens when Gate B1 passes. It does not
define what happens when Gate B1 fails:
- Section 1's time line goes from Gate B1 straight to C/B2/B3/D code. B.4 step 3 is "Gate B1, B1.8, B2.4".
- B1.8 deletes the old path "in the commit after gate B1 passes". Nothing says what happens to the B1 rows, the
  schedule code, or the `--train-views` path on a fail.
- B3.9 trains v5 on "the same B1 rows root" from `rl-full-3ep`. It is not stated whether B3 still runs on B1 rows
  after B1 failed its gate, or whether B3 then pairs against B1 or against R5.
- FR2 ("data fix does not close the miss gap") lists tools (paired CIs, the 4B gate, B3), not a branch.
- The 4B gate (decision 17) fires "after the 2B B runs". With the $3 budget (decision 1 (b')), B3 may never run, so
  the trigger time is undefined.
- "The single final candidate" that reads test is not defined among R5, B1, B3 and 4B.

Grep, 2026-11-08: no "Gate B1 fails", "if Gate B1" or "B1 fails" anywhere in the doc.

Exact fix. Add front section 1.1, "After V1: decision table" (owner approves it before V1, so it cannot be bent
after the result):

| Gate B1 result | Next paid step | Code | Final-candidate rule |
|---|---|---|---|
| Pass | B-S2 (B3 on B1 rows) after the budget decision | B1.8 runs | best of B1 / B3 on val-select by Gate B3 rules; ties go to the simpler model (B1) |
| Fail on far misses only (IoU CI > 0, far-miss rule fails) | B3 on B1 rows (the location term targets the misses that are left) | B1.8 runs | as above |
| Fail on IoU or difference in differences (no data effect) | no B3 spend; Mac diagnosis on A3 and B1 saved answers (per-cell `allmiss_groups`, `kind/*` logs); the 4B gate is decided now on B1's val numbers | B1.8 waits; rows code stays | R5 stays the final candidate unless the 4B arm wins on val |
| Fail on the edge guard or legal | stop; owner review | B1.8 waits | R5 |

Also:
- Decision 17: "the 4B gate is evaluated after the last 2B run that the budget allows (B1 if B-S2 is not funded)".
- C/D Mac work does not depend on the row taken. State that, so the speed work is not blocked by a fail.

### 5. MAJOR. 00_front (6.2, 6.3 item 3b), A (A3, `a3_measure.yaml`), B (B.8.1): each val eval pays about 205 s of fixed load; V1 pays it four times

Problem.
- `final_eval` takes one `--ckpt` and one `--view`. For each call it preloads the train and view stores, forks the
  CPU pool and builds a full GRPO trainer with vLLM [code: experiments/phase2/final_eval.py:48-56, 64-73].
- Measured totals: val 1,084 items in 414.7 s; test 2,200 items in 630.5 s (`load_and_sample`) [measured:
  runs/r5/eval/eval_{val,test}/final_eval.json].
- A two-point linear fit gives 0.193 s per item plus about 205 s fixed per process [estimate: (630.5 - 414.7) /
  (2,200 - 1,084) = 0.1934; 414.7 - 1,084 x 0.1934 = 205.1 s].
- V1 runs four eval processes: B1 on shifted, B1 on edge, R5 on edge, R5 on shifted (A3-val). That is 13.7 min
  of fixed load.
- With one process per checkpoint (two views each), V1 saves 2 x 3.42 = 6.8 min = $0.07 [estimate].
  - At 27.5 s, the spend at the start of the R5 process is 30 + 137.5 + 10.4 = 177.9 min = $1.90, which is below
    $1.98. So A3-val runs.
  - Total V1 is 188.3 min = $2.01.
  - The step-time limit for A3-val rises from 26.9 s to about 29.0 s [estimate: (30 + 5 s + 10.4) min x
    $0.641/h <= $1.98].
- Today that limit is the reason decision 2 offers a cap raise.

Exact fix.
- A3 (pre-E path): `final_eval --view` accepts several views (`action="append"`), with one trainer build and one
  store preload. It writes one output folder per view (`<out>/<view name>/`), with each view's own `seconds`
  record. Add a unit test: two views give the same per-view items as two separate calls (same seed per view).
- Front 6.2: item 3 = one process (B1 on shifted + edge, 10.4 min). Item 4 + item 5 = one R5 process (edge first,
  then shifted with saved answers, 10.4 min). The skip decision for the shifted view is made before the R5 process
  starts.
- Re-price front 6.2, B.8.1, A.8 and decision 2 from this table. The $2.30 option can then be dropped.
- After E4 (sampler without a trainer), the same rule holds for `policy.evaluate`.

### 6. MINOR. B (B1.5, B.6.2 "report only" line), CLEAN (CLEAN6): the B1 run logs no held-out curve, and the only code path for one is archived

Problem.
- B.6.2 says the held-out series does not exist: "`rl.json` has no eval key in its 771 log entries".
- The repo's held-out path is TRL's own eval (`eval_strategy="steps"`), read by `Controller` [code:
  src/objdel/train/callbacks.py:1-11]. CLEAN6 archives `Controller` by default.
- Without a curve:
  - the B-S2 step count (300 or "short 200", B.8.1) is chosen blind;
  - the extra-steps confound (+0.014 for 300 steps [estimate, ADR-0015]) has no in-run evidence apart from the
    training per-cell logs.
- TRL 1.17.1 supports `eval_dataset` with `num_generations_eval` [code: .venv trl/trainer/grpo_config.py:489; web:
  TRL GRPO docs]. This is a library feature, as the owner prefers.

Exact fix.
- B1.5: an optional `--heldout <view> --heldout-items 128 --heldout-every 50`.
  - Use `GRPOConfig(eval_strategy="steps", eval_steps=50, num_generations_eval=4)`.
  - The eval set is 128 items from **val-shape** only, so val-select stays clean.
  - Log `eval/iou_tc` and `eval/allmiss_share` to `rl.json`. Report only: never a stop rule.
- Cost: 6 x 512 answers, about 2 min [estimate: 2 generation batches x about 11 s each]. This is less than item 5
  saves.
- Before the `v1-b1` tag: the tiny TRL loop test (`tests/test_p2_train_rl.py`) runs one held-out eval with
  FastGRPO's patches on.
- If the owner keeps it off, B.6.2 says so and B-S2's step count stays an owner guess.

### 7. MINOR. E (E7 `verify.api`, ADR 0020), 00_front (section 4 "3D" rows): the long-lived interfaces have no 3D note, though 3D voxels are the owner's next direction

Problem.
- E is the last phase before 3D voxels (BRIEF section 3). E freezes `verify.api` at `API_VERSION = "1.0"` with
  `ScoreItem(id, key, text, truncated)`: a 2D polygon-text grammar on a 0-1000 grid, keyed into a 2D doctor store
  (E7 step 1).
- `geom/volume.py` (`betti3d` via `scipy.ndimage.label`) is the only 3D code. D defers its GPU form "to the voxel
  phase" (D.3.3).
- No ADR says which seams survive the 3D move and which do not:
  - The verifier API, the engine, the doctor store and the layer rules can survive.
  - The polygon-text output, the strict parser and the 0-1000 grid do not extend to voxels.
- So the first 3D task will reopen ADR 0020 and ADR 0010 with no recorded trade-off. A top lab writes this
  paragraph while E freezes the interfaces. It costs nothing to build.

Exact fix. Add a "3D extension" paragraph to ADR 0020 (E7), with no code:
- (a) A minor version may add optional `ScoreItem.depth_key` and `spacing` fields. A 3D output format is a major
  version.
- (b) `geom` gains a 3D run CCL (6/26 connectivity) and 3D Betti on the device, through the same engine source and
  CPU oracle rule.
- (c) The output representation (polygon text per slice, a stack of per-slice outlines, or a voxel/mesh token
  format) is the first decision of the 3D phase. Trigger: the 3D direction starts after E. Inputs:
  `docs/implementation2.0/notes/H3_shape_representation.md`.

Front section 4: link the two "3D" rows to this paragraph.

## Counts

blocking 1 (item 1); major 4 (items 2-5); minor 2 (items 6-7).
