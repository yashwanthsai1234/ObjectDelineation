# Final review round 7: fresh outside architect

Scope: the whole `docs/implementation2.0/Implementation.md` (round 7), read by section, with the sources in
`impl/phases/*.md`. I focused on the end-to-end architecture up to Phase E (GPU everywhere, the (b') default, the
V1-V4 rentals, how the phases join). I checked code and kernel claims in the repo (read-only). I checked the review
logs first, so no item below repeats a fixed item.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item needs an NVIDIA catalog skill. The front already lists `compute-sanitizer-diagnose-cuda-bugs`, `nv-segment-ct` and `dicom-series-to-volume`. I did no catalog lookup and installed nothing. |
| engineering:code-review | I read the code lines that the items depend on: `train_rl.py:23-70` (`--train-views` at 26 and 47, the unscheduled path), `serve.py:108-123`, `dataset.py:240-249`, `final_eval.py:20-75`, `configs.py:170-195`, `bundle.py:20-72`, `preflight.py:20-60`, `policy.py:50-80, 180-200, 300-335, 370-500, 600-606`, `fast.py:535-550, 694-718`, `gpucheck.cpp:175-186, 1284-1290, 1574-1802`, the kernel build scripts. `ls src/objdel/train/`: no `gpu_pixels.py`, `vllm_sync.py` or `vllm_requests.py` exist. Most code citations in the doc are exact. |
| engineering:debug | I traced one path end to end: V1, then Gate B1 "Pass", then B1.8 on `main`, then a raise, then V2. The trace found item 1: V2's C arms call a flag and a code path that B1.8 deleted. |
| engineering:system-design | I judged the end state of each owner option. Under (b') the target architecture of front section 2 is never reached, and part of E and of section 5 still assume it (items 2, 3, 5). |
| engineering:architecture | I checked each finding against the governing rule (front order rule 1, "End state under (b')", section 5 kernel rule, the owner archive rule) and named which text must change. |
| engineering:documentation | Each blocking or major item quotes the source file and heading, word for word. |
| huggingface-skills:huggingface-best | No model item. Section 3 (keep Qwen3-VL-2B, 4B gate on val) is sound. No Hub leaderboard covers polygon text on CT. |
| system-design-skills:system-design | Stress test: "what if Gate B1 passes and the owner funds V2" (item 1) and "what if Gate B1 shows no data effect" (item 4). Both paths had no clean end state. |
| engineering:tech-debt | Item 4 (two training data paths left on `main`) and item 3 (kernel rule without a bench) are debt the plan would create. |
| Web (WebSearch) | TRL GRPOTrainer adapter-only vLLM sync: only the experimental async trainers have it (also checked in `.venv` TRL, item 10). CuPy/cuCIM `distance_transform_edt` is exact PBA+; the doc already rejects it with a reason, so no item. |

## Verdict in one paragraph

The core architecture is right: one CUDA engine source for CPU and GPU, the verifier with no torch, data first, one
paid experiment (V1) with a decision table fixed before the result, and a relative speed gate. A top lab would keep
these. It would do three things differently. (1) It would make the C speed arms and B3 run on the exact data path of
the B runs (the scheduled rows root). Today they run on a path that B1.8 deletes (item 1). (2) It would plan the
default end state first. Under (b') most of C, B2.2-B3 and V2-V4 never land. Front sections 2 and 5 and Phase E
still describe the funded state as if it comes (items 2, 3, 5). (3) It would shrink the V1 tree to what V1 needs and
treat the rest of this 15,750-line doc as an appendix (items 8, 9).

## Items

### 1. BLOCKING. V2's C arms train with `--train-views`, which B1.8 deletes on the Gate B1 "Pass" row

Sources: 00_front, C, B.

Problem.
- 00_front, section 1.1 (decision table), row "Pass": "B1.8 runs". Row "Fails only the far-miss rule": "B1.8 runs".
- B.md, "#### B1.8 Clean after gate B1": "In the commit after Gate B1 lands on one of those two rows: delete
  `dataset.grpo_rows_views`, `serve.VIEW_SEP` and `_tagged` ... `--train-views` [code:
  experiments/phase2/train_rl.py:26,47] ... and the "without schedule" branch of `train_rl` (one training path
  left)."
- 00_front, section 6.3 (V2), "Own batches": "The check arm of item 3 reads the same two views (`--train-views`)."
  And under "V2 launch config": "It has no `rows_plan`: the C arms read the bundle's train views, not the B1 rows
  root".
- C.md, C2 tests: "The check steps train on the centred + shifted + rescaled train views (`--train-views`, as
  `rl_shifted.yaml`)". C.md C.8.2: "Check arm = steps 18-32 of the fast arm's process (...), on the centred +
  shifted + rescaled train views". The check arm is the same process as the fast arm, so the fast arm trains on that
  path too.
- Code: `--train-views` and the unscheduled path are `train_rl.py:26` (`ap.add_argument("--train-views", ...)`) and
  `train_rl.py:47` (`roots = [a.train, *a.train_views]`), with `dataset.py:240` (`grpo_rows_views`) and
  `serve.py:108` (`VIEW_SEP`). These are exactly what B1.8 deletes.
- No Mac test catches it: B's I47 parses only `b1_rows.yaml` and `a3_measure.yaml` (B.9, I47 row). C's only
  `c_speed.yaml` test is `preflight.load(c_speed.yaml)` raises `TypeError` (C0 files).

Failure. Time line: V1, then Gate B1 on the Mac. On "Pass" (the row where a raise for V2/V3 is most likely),
B1.8 lands on `main`. The owner branch `after-raise` was created "from `main` after V1" and holds `c_speed.yaml`.
When a raise merges it, `train_rl` has neither `--train-views` nor an unscheduled path. The engineer launches V2,
pays item 1 (42 min, $0.58, front 6.3), and item 3 (C0 baseline, fast and check arms) dies with argparse exit 2 or
the schedule check. C-G-SP1 is never measured. Also, even without the crash, a speed result measured on the
shuffled, oversample-4 `--train-views` path does not describe the B3 run, which trains on the scheduled rows root.
That breaks the front claim "a C-G-SP1 pass carries over to the B3 run" (00_front 6.3, "One preset for the C arms").

Fix.
- 00_front 6.3 and C.md (C0 files, C.8.2 arms): every C arm (baseline, fast, check, 3c, 3d) trains on the B1 rows
  root through the schedule path: `--train {data}/rows/lidc-train-b1 --rows-revision <hex40>
  --rows-manifest-sha256 <hex64>`, with `--steps` a multiple of `schedule.block_steps` (17 and 32 are not: use 20
  and 40, or let `train_rl` accept a step count below one block for arms marked `--timing-arm`, and say which).
  `c_launch.yaml` gets `rows_plan: experiments/phase2/configs/c_speed.yaml`. The V2 launch lines add
  `fetch_hf.sh <run> --rows <rows_revision> lidc-train-b1`. The V2 bundle no longer needs `lidc-train-shifted`
  and `lidc-train-rescaled` (drop them from the bundle command and from 6.1 step 4 and 6.4 rule 7).
  `c_speed_check` draws its 5 batches from the rows root too.
- C.md: a new test (I47 form) parses every `-m` step of `c_speed.yaml` with the module's own parser on the tree
  that V2 will run. Front 6.1 step 4 (V2) lists it.
- B.md B1.8: one line, "after B1.8 no plan may name `--train-views`; C's arms use the rows root (front 6.3)".

### 2. MAJOR. Phase E still says B2, B3 and C land on `main` before E, with "no post-E exception"

Source: E (with 00_front).

Problem.
- 00_front, section 1, order rule 1: "The one exception is the `after-raise` branch merged after E0 ("End state
  under (b')"): it maps its paths through E's move tables (E R-4)." And "End state under (b')": "Code whose adoption
  gate needs V2 or V3 therefore never lands on `main`."
- E.md, "## E.8. Cost and time; dependencies", "Dependencies on other phases": "B2 (trainer-vLLM match) and B3
  (reward v5) land before E on pre-E paths (front order rule 1; no post-E exception)." Next bullet: "Phase C lands
  on pre-E paths; E moves `train/gpu_pixels.py` -> `policy/gpu_pixels.py`, `train/vllm_sync.py` ->
  `rl/vllm_sync.py` ... C's removal of host syncs is in the E-G11 baseline".
- E.md, "### E.3.3 Where each part runs": "The runtime placement stays as Phases C and D left it (they land before
  E)." and "recorded on the E0 commit (after C has removed syncs)".
- E.md, E5 move table (no "if landed" on these rows, unlike `vllm_requests.py`): "`rl/vllm_sync.py` |
  `train/vllm_sync.py` (C5)", "`policy/gpu_pixels.py` | `train/gpu_pixels.py` (C1)", "`policy/marks.py` |
  `train/nvtx.py` (A6, with C0's `count_syncs()`)". E3b row: "`advantages_v5` (B3.4 ...) also stays in
  `verify/reward/advantages.py`". E5 step 1: "move the code exactly as in the table".
- Code: `ls src/objdel/train/` today = `callbacks.py configs.py dataset.py fast.py rewards.py sample.py serve.py`.
  Under (b') none of the C files are ever created on `main`.
- C0's `count_syncs()` (C.md C0, a change to A6's `train/nvtx.py`) is in neither list of front "End state under
  (b')": "Lands on `main`" does not name C0, and "Waits on ... `after-raise`" names only "C0's V2-only tools
  (`c_speed_check`, `c_speed.yaml`, `c_launch.yaml`)".

Failure. Round 6 fixed "When E starts", R-4 and the E.7 row, but these four places were left. An engineer under the
default (b') reads E.8 and E.3.3 as the order rule and waits for B2/B3/C, or follows the E5 move commit "exactly as
in the table" and hits `git mv` on files that do not exist, and must guess whether `count_syncs` belongs on `main`.

Fix.
- E.md E.8 "Dependencies": replace the two bullets with: "B2.0, B2.1 and the Mac part of B2.3b land before E.
  B2.2-B2.4, B3 and C1-C11 land before E only under a raise that funds V2 or V3 before E; under (b') they wait on
  `after-raise` (front section 1 "End state under (b')"), and a merge after E0 maps paths through the move tables
  (R-4)." Delete "no post-E exception".
- E.md E.3.3: "as Phases C and D left it (they land before E)" -> "as Phase D, and Phase C if a raise landed it,
  left it". The E-G11 baseline is recorded on the E0 commit, whatever C state that commit holds.
- E.md E5 and E3b tables: add "(only if landed on the E0 commit; else the row is skipped and its names are absent
  from `old_names.json`)" to the C1, C5, C0 `count_syncs` and B3.4 `advantages_v5` rows. E5 step 1: "exactly as in
  the table, rows whose source is absent on the E0 commit skipped".
- 00_front "End state under (b')": put C0's non-V2 parts (`count_syncs()`, `OBJDEL_NVTX_SYNC`, the Dockerfile nsys
  layer) in one list: either "Lands on `main`" or "Waits on `after-raise`".

### 3. MAJOR. The kernel rule of section 5 says a kernel stays only after a bench; under (b') D's kernels land with no bench

Sources: 00_front, D.

Problem.
- 00_front, "## 5. Custom kernels index": "Rule from K3 section 9: a custom kernel stays only if it beats the
  torch-op or library version in a measured bench. Otherwise the draft is dropped."
- 00_front, section 1, "End state under (b')": "Lands on `main` under (b'): ... D1-D9 (D's box gates stay "open
  until a rental")".
- D.md, D.6, D-G7 "Polygonizer speed (same rental)": "Where: R-D2". 00_front 6.3 V3 item 1 is R-D2, and under (b')
  "no rental follows V1" (00_front 6.3). R-D1 is dropped on the time-line order (00_front time line step 5).
- So under (b') K1 (D2-D3) and K2 (D4-D7) merge into `cuda/gpucheck.cpp` and replace the numpy polygonizer
  (D-G11: "no `src/objdel/checker/polygonize.py` or `bitquad.py`") with no GPU run and no bench.

Failure. Two rules that both claim to decide whether K1/K2 stay give different answers under the default. One
engineer merges D2-D8 (front end state). Another holds them (section 5 rule: no bench, so "dropped"). The outcome
differs: the CUDA polygonizer in production, or not (an owner "must" item).

Fix. 00_front section 5: scope the rule. "K3 (C1 fallback): stays only if it beats the torch-op or HF version in a
measured bench. K1 and K2: they replace our own CPU code, not a library. They merge on CPU-build exactness (D-G3 Mac
half, D-G4, D-G16). D-G7 and D-G8 are speed checks: when R-D2 runs, a failure reverts the kernel (D-R1 rollback);
until then they show "open until a rental"." D.md D.6 header: the same sentence.

### 4. MAJOR. On the "no data effect" row, the old `--train-views` path stays beside the rows path on `main` with no plan and no owner

Sources: 00_front, B.

Problem.
- 00_front section 1.1, row "Fails the IoU rule or the difference-in-differences rule": "B1.8 waits; the rows code
  stays". B.md B1.8: "B1.8 runs only on the B.4.1 rows "Pass" and "Fail on far misses only". On the other rows it
  waits, and the rows code stays."
- B.md B1.5: "Archive `experiments/phase2/configs/rl_shifted.yaml` ... From B1.5 on, the live RL plan is
  `b1_rows.yaml`". `rl_shifted.yaml:26` is the only config that names `--train-views` [code: grep of
  `experiments/phase2/configs/*.yaml`]. So after B1.5 the `--train-views` code (`train_rl.py:26, 47`,
  `dataset.py:240`, `serve.py:108-123`, `g_probe.py:90`) has no live plan.
- No later task removes it on that row: E keeps it ("`train_views` ... only if that flag exists on the E0
  commit", E2.1 `RlStage`); CLEAN names B1.8 only (CT-33, CT-35 "with and without B1.8").
- BRIEF section 3 owner rule: "Old code must be ARCHIVED (not left beside new code) so nobody is confused. Remove
  dead code."

Failure. "No data effect" is a likely row (FR2: Medium / High). An engineer who follows the doc then leaves two
training data paths on `main` through E, DOC and CLEAN waves 2-3, one of them with no plan. That breaks the owner
rule. (Row 4, "Stop. Owner review", may wait; row 3 has no reason to: the effect test failed, not the rows code.)

Fix. B.md B1.8 and 00_front 1.1: "B1.8 runs on rows 1-3. Row 3 failed the effect, not the path: the scheduled rows
root is a superset of `--train-views` (same three kinds, plus the schedule and the data gates). On row 4 the owner
decides in the review." If item 1's fix lands, `--train-views` has no user at all after B1.5. Then the simpler fix is
to delete it in B1.5 itself.

### 5. MINOR. Front section 2 describes the funded end state, not the default

Source: 00_front.

"### 2.3 Target tree (after E and CLEAN)" lists "`archive/ ... reward_v4 (never imported)`", "`policy/ ...
gpu_pixels.py`", "`rl/ ... vllm_sync.py vllm_requests.py`". "### 2.4 GPU-everywhere data path (after B1, C and
D)". Section 2.1: "Pixels do not come back to the CPU after the 32 images of a step are decoded and uploaded".
Under the default (b'): "End state on `main` under (b'): reward v4 only; today's CPU pixel path only". No task
archives `reward_v4` under (b') (no CLEAN or DOC text names it), so this is a reader trap, not a wrong action. Fix:
add one line under 2.3 and 2.4: "This is the state after a raise funds V2 and V3. Under (b'): reward v4 stays live
in `verify/reward/`, there is no `gpu_pixels.py` or `vllm_*.py`, and 2.4's pixel path is today's CPU path (front
section 1 "End state under (b')")."

### 6. MINOR. Two row counts for the B1 rows root

Source: 00_front. Section 0, "rows root": "B1: 8,320 rows = 260 steps x 32". Section 2.4 diagram: "-> B1 rows root
(9,600 rows)". B.md explains 9,600 as the 300-step upper bound. Fix: diagram "8,320 rows (260 x 32)".

### 7. MINOR. The offer limits are a manual check; the code allows a dearer host that cuts the B1 run to about 140 steps

Source: 00_front (section 6 rules, 6.1 step 5).

The doc: "Offer limit: `dph_total` <= $0.641/h and `inet_up_cost`, `inet_down_cost` <= $0.0065/GB. A dearer offer
needs a new plan". Code: `INET_CEIL = 0.03` and the RTX 5090 `ceiling` is 1.50 [code: policy.py:52, 61], so
`reasons()` admits a $0.03/GB host. `BudgetStop` reads the real rate from the run record (B1.5), so it does not
overspend, but at $0.03/GB: b = 15 x 0.03 = $0.45 and r_s = 0.641 + 0.18 = $0.821/h, so `BudgetStop` at $2.05 trips at
minute (2.05 - 0.45) / 0.821 x 60 = 117, and with a 53-min start the B1 run stops near step (117 - 53) x 60 / 27.5
= 139 [estimate: arithmetic], about half the plan. Fix: 6.1 step 4b adds "`policy.launch_refusal` also refuses
`dph_total` > 0.641 or either `inet_*_cost` > 0.0065 for a config with `rows_plan`", with one unit test. That turns
a manual rule into a code gate on the one paid run.

### 8. MINOR (top-lab difference). The V1 tree holds much more than V1 uses

Source: 00_front (order rule 2), CLEAN.

Order rule 2: the tag holds "every Mac task of A (A1-A8)" and "CLEAN wave 1". V1 reads from A only A3 (saved
answers, repeated `--view`) and the budget guard. 6.2 says "V1 has no nsys capture", so A6 adds nothing to V1.
CLEAN2 edits `train_rl.py` (idle pool removed) and CLEAN6 edits `rewards.py` on the V1 tree (CLEAN.md "Order"). Each
extra commit is risk on the only paid run and delays the one measurement that decides B3, the 4B gate and the test
read. A top lab freezes the smallest tree that answers the question: be52b1c + B1.0-B1.7, B1.9 + A3 + B2.0/B2.1
logging + the 6.1 step 0 guard. A1, A2, A4-A8 and CLEAN wave 1 then run on `main` in parallel with V1, because they
read V1's outputs or need no box. Not a defect: the current rule is consistent. It is the main schedule lever.

### 9. MINOR (top-lab difference). The doc is now a hazard in itself

Source: all.

The doc has 15,750 lines and six review rounds, and each round still finds contradictions between the copies (items 1,
2 and 3 are of that kind). Under the default (b') about half of the text (C, B2.2-B3.10, V2-V4, most of the C and B
tests) describes work that waits on a branch that never merges. A top lab would keep two executable artifacts:
(a) the V1 runbook plus the decision table (front 1.1, 6.1, 6.2, B1.5, A3), and (b) the Mac work that lands under
(b') (A, B1, D, E, DOC, CLEAN). It would move the funded-only parts into a "funded appendix" with one line on top
of each: "runs only after a raise (owner decision 1)". DOC4 could then check that no non-appendix section depends
on an appendix task. This also lowers the cost of the next review round.

### 10. MINOR. Library candidate for adapter-only vLLM sync is not named

Source: 00_front (section 4, "RL trainer" row), B (B2.2 option F), C (C5).

Option F builds its own `sync_lora(adapter_name, slot)` path (B2.2 F step 4b). TRL 1.17.1 already ships an
adapter-only sync in its experimental async trainers: `self.vllm_client.load_lora_adapter(f"{self._lora_name}-v{version}",
adapter_dir)` [code: .venv trl/experimental/async_grpo/async_grpo_trainer.py:1936] over vLLM's
`/v1/load_lora_adapter` [code: .venv trl/experimental/async_grpo/vllm_client.py:99-110]. It is server mode, not
colocate, so it does not drop in today. The owner rule prefers a library. Fix: the section 4 "RL trainer" row names
it beside verl as the library option for F, with the trigger "F is picked by B2.4 rule 0 and colocate memory allows
a server-mode vLLM on the same card".

## Counts

blocking 1 (item 1), major 3 (items 2, 3, 4), minor 6 (items 5-10).
