# Final review r2: fresh outside architect (whole Implementation.md, round 2, 13,065 lines)

Stance: no earlier opinion. I read `final_r1_fresh-arch.md` and checked each of its 14 items against the round-2
text. All 14 are fixed in the section that owned them: V1 = B1 first, C11 runaway cut, `adv_scale` removed, the
`edge` guard law, the cuCIM row, deferred rows with owners, the ADR table in the front, the placement table, the gate
and rate table, the 2.5 diagram, EGM Table 4, the start checkpoint, the noise floor, and the NVIDIA CT skills. None of
them is raised again here. Every item below is new, or is a fix that was made in the front but never carried into a
phase section.

## Skills used

| Skill | What it changed in this review |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | The plan's NVIDIA content (CUDA 13 toolchain, compute-sanitizer, NVTX/nsys, CT skills) was checked against the catalog rows that round 1 recorded. I found no new catalog match that would change an item. Nothing was installed. |
| `engineering:code-review` | Checked these citations against the code: `final_eval.py:26-30`, `fast.py:540-547`, `configs.py:73-77`, `run_box.sh:24`, `gpucheck.cpp:1612-1615, 1720-1763`, `augment.py:58-71`, TRL `vllm_generation.py:520-522`, `policy.py:58, 196-200, 318-330, 580-595`, `serve.py:85-125`, `callbacks.py:155-166`, `g_probe.py:93`. All match the doc except the claims in items 1 and 7. |
| `engineering:debug` | Traced the V1 launch path through `policy.cmd_launch` to `spend(committed=True)` and `HARD_STOP`. The ledger sum is $16.90, so the V1 launch is refused (item 1). |
| `engineering:system-design` | Walked the critical path from the B data gates through V1 to Gate B1, and listed every input Gate B1 needs. Two of them are produced by no rental (item 2). |
| `engineering:architecture` | Checked the ADR table against each phase. B and C still use the old numbers, and C10's ADR has no number (item 5). Checked the stack table for a missing decision: the RL framework (item 6). |
| `engineering:documentation` | Checked the cross-references to rentals left over from the round-1 reorder. A still says "G1a", and C still says "validation rental = section 6.2" and uses $0.60/h (items 3, 4). |
| `huggingface-skills:huggingface-best` | HF API, 2026-11-08: no new Qwen VL model of 4 B or less (the Qwen3.8 family is 27 B or larger; trending image-text-to-text models are all above 4 B). Section 3 stands. |
| `system-design-skills:system-design` | Stress-tested the paid path: "what if the ledger guard fires", and "what if A3 is skipped by its own $1.90 rule". Both break Gate B1 (items 1, 2). |
| `engineering:tech-debt` | Looked for stale text left by the round-1 fixes (rental names, ADR numbers, rates, task labels). Items 3, 4, 5 and 9 are documentation debt that would mislead an implementer. |

Web checks: vLLM `LLM.generate(sampling_params: SamplingParams | Sequence[SamplingParams])`, where a list must have
the same length as the prompts [web: docs.vllm.ai/en/latest/api/vllm/entrypoints/llm/]. This confirms C11 step 3 as
written. verl ships a Qwen3-VL GRPO example with vLLM rollout [web: github.com/verl-project/verl
examples/grpo_trainer/run_qwen3_vl_8b_fsdp.sh]. verl also runs LoRA with `lora.merge=False`, which sends only the
adapter deltas to vLLM [web: verl.readthedocs.io/en/latest/advance/ppo_lora.html].

## Verdict

The target architecture is sound, and a top lab would build the same parts:
- data first;
- one CUDA engine source, with a CPU oracle;
- a verifier with no torch;
- layers enforced by import-linter;
- a single TRL adapter;
- move-only refactoring with goldens.

The "GPU everywhere" claims hold up when checked against the code. Each remaining CPU item has a stated reason, and
D.1 item 9 honestly says that D is not an RL speed item.

A top lab would change three things:
1. Make the first paid rental launchable. Today the code's own budget guard refuses it (item 1).
2. Make the first rental produce every number its gate needs (item 2).
3. Write down the RL framework decision, because the plan rebuilds features that a maintained library already has
   (item 6).

The other items are stale text that the round-1 reorder left in A, B and C.

## Items

### 1. BLOCKING. 00_front (6.1, 6.4), 99_back (decision 4), DOC (U19): V1 cannot launch, because the code's Phase-2 hard stop has $0.60 left

Problem.
- `policy.py:58` sets `ALARM, HARD_STOP = 10.0, 17.5`.
- The closed ledger spend is $16.90 [measured: sum of `cost_usd` in `experiments/phase2/vast/ledger.csv`, 133 rows,
  2026-11-08].
- `cmd_launch` dies when `spend(committed=True) + accrued(run, deadline) > HARD_STOP` [code: policy.py:324-326]. For
  V1, that check is 16.90 + a cap of about 2.20 + `PULL_S` = about 19.1 > 17.5. The launch is refused.
- If the cap were cut to make the launch pass, `cmd_sweep_orphans` would destroy every instance once `spend() >=
  HARD_STOP` [code: policy.py:588-592]. `com.objdel.vast.sweep.plist` runs that sweep every 3,600 s. So the box would
  die about 56 min into V1, before the B1 run ends.
- Spend is also above `ALARM` (10.0), so every launch needs `--owner-approved` [code: policy.py:327-329]. Section 6.1
  does not mention that flag.
- Decision 4 in 99_back offers "keep 17.5 and rely on BudgetStop" as an option. That option cannot work.
- DOC test U19 asserts `HARD_STOP == 17.5`, and DOC line 10152 says "the agent does not change it". So the fix would
  also break a planned test.
- Section 6.2's "`BudgetStop` at $1.90" and "skip if the ledger is above $1.90" are ambiguous. The ledger total is
  already $16.90.

Exact fix.
- 99_back decision 4: remove the "keep 17.5" option. New default: before V1, the owner approves setting `HARD_STOP` to
  (ledger spend at plan start, $16.90 measured) + $3.00 = 19.90. The `ALARM` stays at 10.0.
- Add a CLEAN wave-1 task, or a step 0 of front 6.1: edit `policy.py:58` to that number. Add one unit test:
  `spend() + 2.20 <= HARD_STOP` on a copy of the ledger.
- Change DOC U19 and the DOC9 sentence to read the constant from the module, not the literal 17.5.
- Front 6.1: add step 4b, "`policy.spend()` + V1 cap <= `HARD_STOP`; the launch passes `--owner-approved` (spend >
  ALARM)".
- Front 6.2: replace "ledger above $1.90" and "BudgetStop at $1.90" with "session spend (`spend()` minus the
  spend at launch) above $1.90".

### 2. BLOCKING. 00_front (6.2), 99_back (decision 2), B (B.6.2, B.8.1, B.8.4), A (A3 step 7): V1 does not produce the inputs Gate B1 needs, and the V1 table does not match B-S1

Problem. Gate B1 (B.6.2) needs four inputs:
- (a) the B1 model on the val shifted view;
- (b) the B1 model on the val `edge` view (a hard floor);
- (c) R5 on the `edge` view, "from A3" (B.8.4);
- (d) a noise floor, which B.6.2 defines as "Phase A3 evaluates R5 on val twice with two sampling seeds".

What the plan produces:
- Front V1 has only (a) (item 3, 6.9 min).
- A3 (front item 4, A.4) evaluates R5 on val and test. It does not evaluate the edge view.
- A3 step 7 defines the noise floor as the two K = 4 halves of one sampling, and says that "a second val sampling...
  is not planned". This contradicts B.6.2.
- A3 also has a skip rule at $1.90. If it is skipped, input (c) does not exist, and Gate B1 cannot be decided.
- B.8.1 B-S1 costs a different session: 2 val evals + B2.3b (15 min), no A3, 196 min. 99_back decision 2 copies the
  front list.

Cost, with the missing items added to the front table: 202 + 13.8 + 15 = 231 min = $2.47 at 27.5 s. That is above the
$2.20 cap.

Exact fix.
- One noise-floor method: the A3 step 7 halves ($0). Change B.6.2 "Noise floor" and the B.8.4 A3 row to "A3 step 7
  (two K = 4 halves of one sampling, scaled by 1/sqrt 2)". Delete "R5 val eval with a second sampling seed".
- New V1 table in front 6.2. The A3 part is val only:
  1. setup, 30 min;
  2. B1 run, 300 steps;
  3. B1 on val shifted, 6.9 min;
  4. B1 on val `edge`, 6.9 min;
  5. A3-val: R5 val with saved answers (6.9) + R5 on `edge` (6.9) + load + the 6-step nsys window (about 10.5) =
     24.3 min.

  Totals: 205.6 min = $2.20 at 27.5 s; 186.4 min = $1.99 at 23.65 s.
- Why the A3 test re-sampling (10.6 min) can leave V1:
  - A-G3 can pair against the existing `runs/r5/eval/eval_val/items.jsonl`, so it does not need a test sampling.
  - Every user of A3's answers (Gate A, C11, B3.7, the Mac half of C-G-EQ2/3) reads val only.
  - Move the A3 test re-sampling and D7's test answers to V3.
- Move B2.3b (15 min) to V3. Gate B2 picks B3's preset, and B3 is in V3.
- Skip order if the session spend passes $1.90 after step 3: drop the nsys window first, then R5-val; never drop the
  two `edge` evals.
- B.8.1: B-S1 = this table. B-S2 gets B2.3b and the A3 test re-sampling.
- 99_back decision 2: copy the new item list.
- Bundle: B1.6's frozen `lidc-val-edge` view goes into the V1 bundle. Add a front 6.1 check that the bundle holds it.

### 3. BLOCKING. A (A.3.1, A.4 header, A3 steps 3/6, A6 step 8, A.6 A-G6, A.8, A.9.6, G06): Phase A still plans the removed rental G1a, with Gate A before the B spend

Problem. The front removed G1 in round 1. Phase A still says G1a in about 20 places (lines 742-755, 809-812, 961,
1014-1035, 1097, 1492, 1513-1551, 2079, 2156, 2178). The text says:
- "A3 = first G1a item after setup";
- "the A6 Nsight baseline is the C0 baseline arm capture of the same rental" (17 steps, window 16-17);
- "Gate A is decided on the Mac after the G1a rental, before the owner approves the B spend".

The front says otherwise:
- A3 is V1 item 4, after the B1 run.
- A3 has its own 6-step window (capture steps 4-5).
- The C0 arms are in V2.
- Gate A is decided after V1.

An engineer who follows A would rent an A3 + C0 box before the B1 spend. That breaks the binding data-first order.
Also, A3 step 3 runs `make preflight CONFIG=experiments/phase2/configs/r1_qwen3vl2b.yaml`, which is the SFT plan
(E6, line 8631). It is not the V1 RL/eval configs.

Exact fix.
- Replace every "G1a" in A with "V1".
- A.3.1 box diagram: replace the "C0 baseline arm (a C item of G1a)" block with "A6 nsys window: 6 steps of the
  `rl-full-3ep` config, `--profile-window 4 5` (V1 item 5)".
- A.4 header: "A3 is V1's last item, after the B1 run and its evals (front 6.2)". Drop the reason "a cap slip must
  not cut it" (item 2 above sets the skip order).
- Line 1097: "Gate A is decided on the Mac after V1. It does not gate the B1 spend."
- A3 step 3: `make preflight CONFIG=experiments/phase2/configs/b1_rows.yaml` and `CONFIG=.../a3_measure.yaml`.
- A-G6 and G06: take the box half from the A3 window. The C0 capture in V2 is C's own baseline.
- A.8 "Setup ... belong to G1a items 1 and 7": change to "V1 item 1 (setup) and item 5 (A3-val)".

### 4. MAJOR. C (C.8.2, C.8.3, C.R3 item 15, C11 step 3, C.4 header line 4952, C0 line 4976, C5 line 5301): stale rental mapping and rate after the round-1 reorder

Problem.
- C.8.2 says the C items are "inside the front's validation rental, section 6.2". Section 6.2 is now V1 (the B1 run).
  Validation is V2 (6.3).
- C.8.2 says "Rate: the front's section 6 planning rate, $0.60/h". Front section 6 sets $0.641/h. C.R3 item 15 repeats
  $0.60.
- C.8.3 says "A3 runs in the same rental as the C items". A3 is in V1 and the C items are in V2.
- The payback line uses $0.60/h.
- C11 step 3 runs the vLLM check "inside the first B rental". V1 is the first B rental, but it runs today's code, and
  C11 is built only after V1 (its sizing needs A3's answers). So the check can only run in V3.

An engineer cannot tell which rental takes the C items, or which rate the C caps use.

Exact fix.
- C.8.2 title and first lines: "inside V2 (front 6.3)"; rate $0.641/h. Remove the $0.60 column and keep one column.
- C.8.3 A3 row: "A3 runs in V1. The Mac half of G-EQ2/G-EQ3 runs after V1, before V2".
- C.R3 item 15: "$0.641/h (front 6)".
- C11 step 3: "box, first rental after C11 lands (V3), 1 min, before the B3 run".
- Lines 4952, 4976 and 5301: "validation rental (V2, front 6.3)".

### 5. MAJOR. B (B.0, B.3.4, B2, B.7, B.8.3-4, B.R, B.9, B.R2), C (C5, C.8.3, C.R3, C10): ADR numbers were not renumbered to the front table, and C10's ADR has no number

Problem.
- The front table (section 4) sets: 0012 = docs site (DOC10), 0013 = number match (B2), 0014 = reward v5 (B3),
  0015 = scheduled rows + numpy-polygonizer exception (B1).
- DOC, D, E, CLEAN and 99_back follow it.
- B and C still use the old numbers: ADR-0012 = number match, 0013 = reward v5, 0014 = scheduled rows / polygonizer
  exception. Lines: 2747, 2763, 2934, 2945, 2966, 2975, 3031, 3339-3464, 3494, 3508, 3601-3644, 3658, 3754, 3761,
  3824, 3847, 3850, 3871-3884, 3984, 4560, 4645-4691, 5300, 5688, 6573.
- Result: ADR-0012 names two different files (DOC's docs site and B's number match). B.8.3 item 2 says "ADR-0014" for
  the polygonizer exception, while 99_back decision 9 says "ADR 0015".
- C10 asks for an ADR ("GPU pixel backend and K3; one-call full push; LoRA serving not adopted") "with the number from
  the front's ADR allocation table". The table has no Phase C row.

Exact fix.
- In B and C only, renumber in descending order so that no number is replaced twice: 0014 -> 0015, then 0013 -> 0014,
  then 0012 -> 0013.
- Add the row "0022 | GPU pixel backend, K3, one-call weight push | C10" to the front table.
- Add one DOC4 drift check: every `ADR[- ]00NN` in `docs/implementation2.0/Implementation.md` maps to exactly one title
  in the front table.

### 6. MAJOR. 00_front (section 4; ADR 0019): no decision row for the RL training framework, although the plan rebuilds features that a maintained library has

Problem.
- The owner rule is "prefer a good library (even with a little customization) over our own code".
- FastGRPO overrides six private TRL methods and patches a module global (E.1 item 7).
- Phase C then builds, on top of those patches: GPU pixels (C1), encoder-once (C3), host metadata (C4), a one-call
  weight push (C5), and per-row caps (C11).
- verl has a Qwen3-VL GRPO recipe with vLLM rollout and dynamic token batching. It also runs LoRA with
  `lora.merge=False`, sending only the adapter deltas to vLLM (see the web checks). That last feature is the fix for
  the merge-rounding mismatch that B2 option F looks for.
- Section 4 has rows for tracking, profiling, orchestration, the spec and data. It has no "RL trainer" row. ADR 0019
  ("TRL adapter and start-up check") records how to live with TRL, not whether to.

A top lab would keep TRL for B1 (one variable at a time, $3). It would also write the decision down with a trigger.

Exact fix.
- Add a section 4 row: "RL trainer | KEEP TRL 1.17.1 + FastGRPO behind `rl.trl_adapter` | E5 (ADR 0019) | E.1, web
  verl docs | Considered: verl (Qwen3-VL GRPO, vLLM rollout, LoRA `merge=False`); not now: B1 must run on today's
  code, and Ray + FSDP on one 32 GB card is unmeasured. Trigger: a TRL bump fails `check_private_names()` on more
  than 2 names, OR B2.3b points to merge rounding, OR the 3D phase needs a second GPU."
- ADR 0019: add "Options considered: TRL + adapter (chosen), verl, own minimal loop", with the same trigger.
- Optional, only when a trigger fires: a $0 Mac spike that runs verl's reward-function interface with
  `verify.api.score` on 8 rows.

### 7. MINOR. D (D.1 item 9, D.R, D-G14, D7 step 5, D.0 system-design row): the "audit callback inside RL" path does not exist in live runs, and CLEAN6 archives it before D

Problem.
- D.1 item 9 says: "the audit callback calls `alarms_for -> metrics()` [callbacks.py:162, eval/audit.py:154]". It uses
  this to put D7 inside the RL process and to justify D-G14 ("with the audit callback on").
- But the `Controller` is never attached in live runs: `g_probe.make_trainer` builds `calls = [callbacks.TraceStep(...)]`
  only [code: experiments/phase2/g_probe.py:93]. No `experiments/phase2/*.py` names `Controller`.
- CLEAN6 (wave 1, before D) archives `Controller` and `alarms_for` by default.
- So the D-G14 precondition cannot be met, and its test would be skipped or fail.

Exact fix.
- D.1 item 9: delete the audit sentence. Keep the reward-path reason (`add_tables`, `drop_tables`, arenas).
- D-G14: "First 6 steps of the first B rental after D8: `arena_bytes` and `max_memory_reserved()` from the reward path.
  Add the eval path (`eval_sets`) in `final_eval` on the same box." Delete "with the audit callback on" there, in D7
  step 5 and in D.0.
- If the owner picks CLEAN6's alternative branch (`--audit`), D-G14 adds the audit leg back.

### 8. MINOR. 00_front (section 1 order rule 2, 6.2 header): "today's code" for V1 is not a defined tree

Problem.
- V1 is said to run "today's code with `--shared-prefix`".
- Before V1, several changes land on the Mac:
  - B1.5 edits `train_rl` (rows root, `rows_revision` / SHA check) and adds the `BudgetStop` callback (`BudgetStop`
    does not exist in `src/` today [measured: grep]);
  - CLEAN2 edits `train_rl` (idle pool);
  - CLEAN6 edits `train/rewards.py`;
  - CLEAN1 edits the bundle manifest.
- So the V1 tree differs from be52b1c. The plan has no tag and no single gate that ties the paid run to a tested tree.

Exact fix.
- Order rule 2: "V1 runs owner tag `v1-b1`: be52b1c + B1 + CLEAN wave 1. Before the tag: B E08 dry run, the reward
  goldens, `tests/test_p2_train_rl.py` (tiny TRL loop) and `make preflight` on `b1_rows.yaml` and `a3_measure.yaml`
  all pass on the tagged tree."
- Front 6.1 step 4: name those two configs.

### 9. MINOR. 00_front (section 2.5): two labels in the step diagram are wrong

Problem.
- (a) "no separate old-logp pass, C10": in this doc, C10 is "Archive, delete, document". The old-logp skip is the code
  name `skip_old_logps` (its "C10" comes from an older plan) [code: fast.py:540-547].
- (b) "push LoRA weights to vLLM (one call, C5)": C5 pushes the merged W + BA weights. That is every one of about 700
  tensors (1.41 B LoRA-touched parameters) in one call, not the LoRA weights. vLLM LoRA serving is "not adopted for
  speed" [C5, lines 5282-5300].

Exact fix.
- (a) Write "no separate old-logp pass (`FastGRPO.skip_old_logps`, existing)".
- (b) Write "push merged weights W + BA to vLLM (all tensors, one call, C5)".

## Counts

blocking 3 (items 1-3); major 3 (items 4-6); minor 3 (items 7-9).
