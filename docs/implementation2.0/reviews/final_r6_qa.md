# Final review round 6: QA lead and technical writer (ASD-STE100)

Target: `docs/implementation2.0/Implementation.md` (round 6 build; 1,542,206 bytes, 15,440 lines [measured: `ls -l`,
`wc -l`]), built from `docs/implementation2.0/*.md` by `assemble.py`. I read front sections 1, 1.1 and 6 in full, the
round-5 changes in every phase (from `review-log.md`, "final round 5" parts), the Tests sections A.9 (box part),
B.9 (U49, U51, I22a/b, I40/I40b, I47, R01-R01c, U70, U76, E10-E10c, G02-G08b, the B.9.8 sketch), C.9.0/C.9.5,
D.9.6, E.9 (conventions, U51, GP1-GP5), and the E.8 V4 plan, plus 99_back decisions 1-4. I checked each item
against `review-log.md` (rounds 1-5); fixed items are not repeated. No REPO file was edited. No git write. No
rental. No Python read `data/`.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No new NVIDIA product claim in the round-5 text, so no catalog lookup. It pointed me at the box toolchain path again (`sync_up` installs the CPU build; `run_box.sh --build-only` makes the CUDA one), which I used to re-check that `e_v4.yaml` step 0 now exists (it does) and then what V4 still lacks (items 1, 3). |
| engineering:code-review | Read the code behind each launch and test claim: `policy.py:376-384` (`code_files` = `git ls-files --exclude-standard`), `policy.py:422-492` (`chain`, `cmd_run`, `OBJDEL_RUN_DIR`), `policy.py:308-316, 527-536` (stamp check, `cmd_fetch_hf`), `.gitignore:23` (`/runs/`), `env_check.sh:3-9`, `fast.py:567-590` (`set_is_ratio` only under `use_vllm`), `tests/test_p2_scripts.py:392-393` (`GRPO_SMOKE`), `train_rl.py:18-88`, `final_eval.py:46-84`. Items 1, 3, 5. |
| engineering:debug | Traced V4 step by step on a fresh box: what `sync_up` ships, what step 1 needs, what GP5 reads at step 2, what exit code a skipped GP5 gives (0, measured). Items 1, 3. |
| engineering:system-design | Put B2.4 rule 0, the B2.4 "Missing inputs" rule, B2.2 F step 5, front 6.3 ("V2 has no F arm"), U51 and the B.9.8 sketch side by side: the in-run F check makes rule 0 unreachable (item 2). |
| engineering:architecture | Checked the launch-config role (front section 0 "launch config", 6.1 step 4, `cmd_launch`) for every rental: only V1 (and V4 after E6) has one (item 4). |
| engineering:documentation | One rule per place: A3 step 5 still says "stop" for the correlation (item 7); 99_back 3b price text (item 8); `e_v4.yaml` outside `specs/` (item 9); test id S09 used twice (item 11). |
| huggingface-skills:huggingface-best | No model choice changed. Used only for the HF part of V4: the checkpoint must come from the private repo by `fetch_hf.sh` from the Mac (token on ssh stdin), not by a plan step (item 3). |
| engineering:testing-strategy | Matched each box test to its input, its step order and its rental (GP5 needs step 3's output at step 2, item 1); matched each B2.1 buffer test to the loop that can feed it (item 5). |
| mattpocock-skills:tdd | Looked for tests that pass by construction or cannot fail: GP5 skips, and a skip exits 0, so V4 "passes" with no TOST (item 1); the R01b/R01c tiny loop cannot fill the buffers it asserts on (item 5). |
| Web | pytest exit codes page: 5 = no tests collected; an all-skipped run is not listed [web: https://docs.pytest.org/en/stable/reference/exit-codes.html], so I measured it: one skipped test exits 0 [measured: `.venv/bin/python -I -m pytest` on a one-test file in `(planning probe, not kept: r6skip)`, 2026-11-08]. |

Measured or read during this review (2026-11-08):
- `git check-ignore -v runs/r5/eval/eval_val/items.jsonl` -> `.gitignore:23:/runs/`. `code_files()` lists files with
  `git ls-files -z -co --exclude-standard` [code: experiments/phase2/vast/policy.py:376-384], so nothing under `runs/`
  reaches the box by `sync_up`.
- An all-skipped pytest run exits 0 (above).
- `GRPO_SMOKE = {"per_device_train_batch_size": 8, "gradient_accumulation_steps": 1, ..., "use_vllm": False}`
  [code: tests/test_p2_scripts.py:392-393]; `tests/test_p2_train_rl.py` passes it to `train_rl --grpo`.
- `FastGRPO._compute_loss` sets `loss_inputs` only `if ... self.use_vllm and self.vllm_importance_sampling_correction`
  [code: src/objdel/train/fast.py:568]; `set_is_ratio` runs only when `loss_inputs` is set [code: fast.py:576-578].
- `env_check.sh MODEL_DIR [OUT_JSON]`: `MODEL=${1:?usage ...}` and a `vllm_engine` check that loads Qwen3-VL-2B at
  MODEL_DIR; a crashed check is a failed check, exit 1 [code: experiments/phase2/docker/env_check.sh:3-9, 17-24].
- `grep -n "_launch.yaml"`: the only launch config named in the doc is `b1_launch.yaml`.

Severity follows the task rules. Count: **2 blocking, 3 major, 7 minor**.

---

## Blocking

### 1. [blocking] [E, 00_front] V4 cannot decide the E-G10 TOST: GP5 runs before the sampler arm, its trainer-arm file never reaches the box, and its skip exits 0

Problem. Two parts both claim to be the rule, and the V4 rental ends with E4 still "not done".
- E "## E.8", V4 steps: "2. the E.9 box test list by file (`pytest -m gpu tests/verify/test_verifier_gpu.py
  tests/infra/test_box_smoke.py tests/policy/test_sampler_box.py`) ... 3. the E4 sampler arm: `final_eval` on
  `lidc-val-shifted`, K 8, seed as R5, on the CUDA build with pool size 1". 00_front "### 6.3", V4: "Then come
  `env_check.sh` (E-G13), the E.9 box test list, and the sampler-arm `final_eval` (pool size 1, as A3)". So the
  sampler arm runs after the box tests.
- E "### E.9.5", GP5 (in `tests/policy/test_sampler_box.py`, so in step 2): "trainer arm read from
  `runs/r5/eval/eval_val/items.jsonl`, sampler arm run on the box on the same checkpoint, view, K and seed; the 90 %
  patient-level ... CI ... lies inside +/-0.01".
- E "## E.9. Tests", conventions: "A box test that needs a run artifact skips with a reason when the artifact is
  absent: ... GP5 without the sampler-arm output or `runs/r5/eval/eval_val/items.jsonl`."
- E "## E.6", E-G10: "box (V4; open until V4 runs): ... TOST ±0.01 IoU_tc on 1,084 val-shifted items, sampler arm
  on the box against the trainer arm on disk (`runs/r5/eval/eval_val/items.jsonl`) ... E4 is done only after the box
  half". E4 gate: "**E4 is not done until E-G10's box part passes.**"
- Code: `runs/` is git-ignored [code: .gitignore:23], and `sync_up` ships only `git ls-files --exclude-standard`
  [code: experiments/phase2/vast/policy.py:376-384, 390-392]. No V4 step uploads the R5 file. A skipped test exits 0
  [measured: (planning probe, not kept: r6skip,) exit 0], so step 2 does not stop the `&&` chain [code: policy.py:490].

Result. At step 2 GP5 has neither input: the sampler arm is made later (step 3), and the trainer-arm file is not on
the box. GP5 skips, step 2 exits 0, step 3 runs the 1,084 x 8 eval, and no step or Mac test computes the TOST from
it. The V4 rental ($0.58, front 6.3, after a budget raise) ends with E-G10's TOST half open, so E4 stays "not done"
and the old samplers are never archived. GP5's text "sampler arm run on the box" also reads as a second 7-min eval
inside pytest.

Exact fix (keep the sampler arm on the box; compute the TOST on the Mac, as A does for G01):
- E.9.5 GP5: rename to `test_e_g10_tost` in `tests/policy/test_sampler_tost.py`, marker `data` (no `gpu`), run place
  "Mac, after `pull.sh` of V4", command `.venv/bin/python -I -m pytest -q
  tests/policy/test_sampler_tost.py::test_e_g10_tost`. Inputs: `runs/<v4 run>/eval_v4_sampler/items.jsonl` (step 3
  output) and `runs/r5/eval/eval_val/items.jsonl`. Same CI rule and +/-0.01 bound. It reads no GPU.
- E.8 V4 step 3: name the output `--out {rd}/eval_v4_sampler`. Add after the step list: "After `pull.sh`, the Mac
  runs GP5 (TOST). E-G10's box half passes when GP4 passes on the box and GP5 passes on the Mac."
- E.9 conventions: the "one file list" keeps GP1-GP4; replace the GP5 skip line with "GP5 is a Mac test after the V4
  pull (`data` marker); it fails, never skips, when `eval_v4_sampler/items.jsonl` is absent after V4."
- E-G10 row and E4 gate: replace "sampler arm on the box against the trainer arm on disk" with "sampler arm made
  on the box (V4 step 3), TOST on the Mac after the pull (GP5)".
- 00_front 6.3 V4: add "The TOST is computed on the Mac after the pull (E.9 GP5)."
- U51 (E): add "step 3 writes `--out {rd}/eval_v4_sampler`".

### 2. [blocking] [B] B2.4 rule 0 can never pick `lora-seqmask`: the F check runs only inside a `lora-seqmask` run, but the rule returns `bf16-tis-opsm` when the F check has not run

Problem. Two rules of B2.4 both claim to decide the preset of the first run after Gate B3, and the test encodes the
one that blocks rule 0.
- B "#### B2.4", rule 0: "If B2.0 shows a bf16 relative error >= 0.2 at step 100 (the seed run gives 0.531), the
  next run uses `lora-seqmask` if its F check passes. The F check is the first 15 steps of that run itself, in one
  process on one host, run by `train_rl --f-check` in the training process (B2.2 F step 5; no separate pytest)".
- B "#### B2.4", end of the rule list: "Missing inputs: when the F check or the fp16 arm did not run, its rule is
  skipped, and the rule returns `bf16-tis-opsm` (or `bf16-seqmask` by the rule 3 veto)."
- B "#### B2.2", F step 5: "`train_rl --f-check` (on by default under `lora-seqmask` and `fp16-seqmask`; setting it
  under another preset is an argparse error)". 00_front "### 6.3": "V2 has no F arm. The F check is the first 15
  steps of the first run after Gate B3 (B2.4 rule 0)."
- B "### B.9.1", U51: "Inputs absent (F check and fp16 arm not run) -> bf16-tis-opsm". B "### B.9.8" sketch:
  `(0.531, None, None, "bf16-tis-opsm"),   # inputs absent (checks not run)`.
- Code: today no F check exists anywhere; `train_rl` has no `--numerics` or `--f-check` flag [code:
  experiments/phase2/train_rl.py:23-36], so the only source of F-check values is the in-run check of B2.2 F step 5.

Result. `decide_numerics` is called before the first run after Gate B3 (it writes ADR-0013 and picks that run's
`--numerics`). At that time the F check cannot have run, because it runs only inside a `lora-seqmask` run. With the
step-100 error 0.531, the function built to U51 and the sketch returns `bf16-tis-opsm`, so no run ever starts under
`lora-seqmask`, the F check never runs, and rule 0 (the fix for the measured 0.531 rounding error) never fires.
ADR-0013 then records the wrong preset.

Exact fix (two-stage rule; keep the in-run check):
- B2.4 rule 0: "Before the run: B2.0 error >= 0.2 -> start the run under `lora-seqmask` with `--f-check`. In the
  run: the F check (B2.2 F step 5) passes -> keep it; it fails -> `train_rl` raises and the run restarts under
  `fp16-seqmask` (rule 1 passed) or `bf16-tis-opsm`."
- B2.4 "Missing inputs": "The fp16 arm not run: rule 1 is skipped. The F check is never a missing input before the
  run: rule 0 reads only the B2.0 error; the F result is read only by the restart choice."
- `decide_numerics`: split into `start_preset(cf, rel_err_100, fp16, b1)` and `restart_preset(f_check, fp16, b1)`.
- U51 and the sketch: replace `(0.531, None, None, "bf16-tis-opsm")` with `start_preset(rel_err=0.531, fp16=None)
  == "lora-seqmask"`; move the `F(...)` cases to `restart_preset` (F passed -> `lora-seqmask` kept; F failed and fp16
  half -> `fp16-seqmask`; F failed, no fp16 -> `bf16-tis-opsm`); keep `(0.1999, ...) -> "bf16-tis-opsm"`.

---

## Major

### 3. [major] [00_front, E] V4 has no launch lines and no checkpoint fetch: `env_check.sh` step 1 has no model folder, and GP4 and the sampler arm name no checkpoint

Problem.
- 00_front "### 6.3", V4 paragraph, and E "## E.8" V4 steps 0-3 list the plan steps only. Neither gives launch
  lines or a `fetch_hf.sh` call. E.8 step 1: "`env_check.sh` (E-G13: the `trl_private_names` check and
  `IMPORT_CHECK`)", with no argument. E.8 step 3: "`final_eval` on `lidc-val-shifted`, K 8, seed as R5" (no
  `--ckpt`). GP4: "On 64 val-shifted rows of one checkpoint" (no path). E U51 checks only the step kinds.
- 00_front "### 6.2", code-change list: "The token reaches the box only on ssh stdin, so a plan step cannot fetch the
  private repo."
- Code: `env_check.sh` needs `MODEL_DIR` (`MODEL=${1:?usage: env_check.sh MODEL_DIR [OUT_JSON]}`) and its
  `vllm_engine` check loads Qwen3-VL-2B there; a failed check gives exit 1 [code:
  experiments/phase2/docker/env_check.sh:3-9]. `cmd_fetch_hf` is the only way a revision reaches `{rd}/ckpt`
  [code: policy.py:527-536].

Result. An engineer who runs `launch.sh`, `sync_up.sh` and `run.sh <run> e_v4.yaml` gets a failed step 1 (no
`{rd}/ckpt`), and the `&&` chain stops before the box tests: the V4 money is spent for nothing. The TOST also needs
the R5 checkpoint (`rl-full-3ep`, the `ckpt` of the trainer-arm file, E4 gate), which the doc never names for V4.

Exact fix.
- 00_front 6.3 V4: add launch lines in the V1 form: `launch.sh OFFER <run> --config
  experiments/phase2/configs/e_v4.yaml --cap-hours <floor_to_0.05((approved - b) / r_s)> --owner-approved`;
  `sync_up.sh <run> && fetch_hf.sh <run> rl-full-3ep`; `run.sh <run> experiments/phase2/configs/e_v4.yaml`;
  `pull.sh <run> && teardown.sh <run> && sweep_orphans.sh`.
- E.8 V4 step 1: `bash experiments/phase2/docker/env_check.sh {rd}/ckpt {rd}/env_check.json`; step 3: `--ckpt
  {rd}/ckpt` (`rl-full-3ep`, the trainer arm's checkpoint). GP3 reads `$OBJDEL_RUN_DIR/env_check.json`; GP4 reads
  `$OBJDEL_RUN_DIR/ckpt`.
- E U51: add "step 1 argv is `env_check.sh {rd}/ckpt {rd}/env_check.json`; step 3 has `--ckpt {rd}/ckpt`".

### 4. [major] [00_front, B, C] V2 and V3 have no launch config, so `make preflight` and `launch.sh --config` have no file to name, and the B3 push revision is undefined

Problem.
- 00_front "## 6": "This section is the only place that lists rental items, their order, minutes, dollars, caps,
  skip rules and the launch configs that `make preflight` stamps." 00_front "## 0. Words", launch config: "`make
  preflight` stamps it and `launch.sh --config` names it. V1: `b1_launch.yaml`." 00_front 6.1 step 4 stamps only
  `b1_launch.yaml`. 6.3 lists V2 plan files (`c_speed.yaml`, `r_d1.yaml`) and V3 plan files (`b3_v5.yaml`,
  `a3_test.yaml`, `d_gpu_test.yaml`) but no launch config and no launch lines for either rental. D "### D.9.6"/D1
  step 0: "Preflight runs it [S09] when the plan is `r_d1.yaml`", but a plan is never stamped (front 6.1 step 4).
- Code: `cmd_launch` refuses without a passing stamp of its `--config` file [code: experiments/phase2/vast/policy.py:313-316];
  `preflight.load` builds `RunArgs(**raw)`, so a plan file cannot be stamped (front 6.2 "Why a new field"); the
  push revision is the launch config's `output_dir` [code: policy.py:517-521].

Result. For V2 and V3 the engineer must invent a launch config. Reusing `b1_launch.yaml` stamps the B1 plan (its
`rows_plan` is `b1_rows.yaml`), not `b3_v5.yaml`, so the "ROWS root" check never sees the B3 stage, and
`push_hf.sh` would push `/workspace/runs/<run>/b1-rows`, which does not exist on V3. The B3 final model then has no
HF revision (it leaves the box only by `pull.sh`).

Exact fix.
- 00_front 6.3: add "Launch config" lines. V2: `experiments/phase2/configs/c_launch.yaml` (RunArgs: `model:
  qwen3-vl-2b`, `view: views/lidc-bundle-v3`, `output_dir: c-speed`, `rows_plan: experiments/phase2/configs/c_speed.yaml`).
  V3: `experiments/phase2/configs/b3_launch.yaml` (same form, `output_dir: b3-v5`, `rows_plan:
  experiments/phase2/configs/b3_v5.yaml`). Add launch lines in the V1 form, with `push_hf.sh <run>
  experiments/phase2/configs/b3_launch.yaml` before teardown, and `b3_v5.yaml` writing `--out {rd}/b3-v5`.
- D1 step 0 / D.9.6 S09: replace "Preflight runs it when the plan is `r_d1.yaml`" with "6.1 step 4 runs it before
  V2 (it is a `make test` test)".
- B.9.3 I22b: add the `b3_launch.yaml` checks (RunArgs loads; `rows_plan` names `b3_v5.yaml`; `output_dir` equals
  the basename of step 3's `--out`).

### 5. [major] [B, 00_front] R01b and R01c (and R01) use the tiny TRL loop, which never fills the B2.1 buffers: it runs with `use_vllm: False` and 8 rows per step

Problem.
- B "### B.9.16", R01b: "a 3-step tiny TRL loop with `logging_steps=1` writes `<out>/seq_gaps.f32` of exactly 3 x
  256 x 3 float32 values in step order ...; the B1 tag's `train_rl` writes the file (P13 tree)". B "### B.9.17",
  R01c: "the R01b tiny TRL loop ... writes `seq_gaps.f32` of exactly 2 x 256 x 3 float32 values". B "### B.9.7",
  R01: "The 5 old `sampling/*` metrics equal goldens ...; every new key of B2.1 (... the `cf/*` keys) appears
  exactly once per logging step". 00_front order rule 2 puts "`tests/test_p2_train_rl.py` (tiny TRL loop); B2.1's
  P13 and R01" in the pre-tag list.
- Code: the tiny loop passes `GRPO_SMOKE` = `per_device_train_batch_size 8`, `gradient_accumulation_steps 1`,
  `use_vllm: False` [code: tests/test_p2_scripts.py:392-393; tests/test_p2_train_rl.py:24-27], so one step has 8
  rows, not 256. `FastGRPO` feeds `set_is_ratio` (the source of `is_stats`, `seq_gaps` and the `cf/*` keys) only
  when `self.use_vllm and self.vllm_importance_sampling_correction` [code: src/objdel/train/fast.py:568, 576-578].
  TRL logs `sampling/*` only under the same condition [code: .venv trl/trainer/grpo_trainer.py:2979].

Result. On the tiny loop as it exists, `seq_gaps.f32` stays empty and no `cf/*` or `sampling/*` key appears, so
R01b, R01c and R01 fail for a correct B2.1, and the 3 x 256 x 3 size cannot occur. R01 is a pre-tag test (order
rule 2), so the owner cannot tag `v1-b1` until the engineer invents a fake vLLM sampling path and a 256-row batch.
The B2.4 rule 3 and 4 pricing (read only from `seq_gaps.f32`) then rests on an untested writer.

Exact fix.
- B.9.0: define one fixture `tiny_loop_vllm_fake` in `tests/conftest.py`: the tiny loop with `use_vllm: True`
  routed to a fake `VLLMGeneration` that returns the HF-sampled ids plus planted `sampling_per_token_logps` (seeded),
  `vllm_importance_sampling_correction: True`, `vllm_importance_sampling_mode: "sequence_mask"`.
- R01, R01b, R01c: run on that fixture. Replace "3 x 256 x 3" and "2 x 256 x 3" with "steps x rows_per_step x 3,
  rows_per_step = `per_device_train_batch_size x gradient_accumulation_steps` of the fixture (8 with GRPO_SMOKE)",
  and assert the row count from the fixture's config, not a literal 256.
- 00_front order rule 2: name the fixture beside R01.

---

## Minor

### 6. [minor] [B, 00_front] Fallback case: the frozen-view hash rule names an A3 file that does not exist

B.6.2 Gate B1: "frozen-view hash: each B1 and A3 `final_eval.json` `view_content_sha256` equals its
`frozen_views.json` entry"; front 1.1 row 4 makes a failure "Stop". In the front 6.2 fallback ("only if the R5
process (items 4-5) gives no output") there is no A3 `final_eval.json`. U49 handles the fallback with the hash chain
but does not say how the A3 half of this rule reads. Fix: B.6.2 and front 1.1 row 4: "each B1 and A3 (A3 only when
item 5 ran) ..."; U49: "fallback: the A3 half of the frozen-view hash rule is not evaluated; the B1 half is".

### 7. [minor] [A] A3 step 5 still says "stop" when the correlation fails

A "#### A3" step 5, "View identity (hard, val ...)": "A lower correlation means a different view: stop, and pair
only on rebuilt, hashed views". The next bullet, B.6.2 and front 1.1 say Gate B1 continues on A3's items. Fix:
replace "stop, and pair only on rebuilt, hashed views" with "the R5 `items.jsonl` line is dropped; Gate B1 still
pairs with A3's items (next bullet)".

### 8. [minor] [99_back] Decision 3b gives the old V4 work time

99_back 3b: "V4 is its own rental: 30 min overhead + about 9 min of work". Front 6.3 (round 5): "30 min overhead + 3
+ 9 = about 42 min = $0.58". Fix: "30 min overhead + 3 min CUDA build + about 9 min of work (front 6.3 holds the
price)".

### 9. [minor] [E] `e_v4.yaml` is a `RunSpec` outside `specs/`

E6 step 2 and E "Hand-off to DOC": "After E6 the live plans are `specs/*.yaml`"; `tests/core/test_spec.py` loads
"every `specs/*.yaml`". E6 step 8b: `e_v4.yaml` is "a `RunSpec` at the path front 6.3 names"
(`experiments/phase2/configs/e_v4.yaml`). Fix: put it at `specs/e_v4.yaml` (front 6.3, E.8, E U51, 99_back 3b), so
the spec-load test, DOC4 `config_files` and `docs/reference/configs.md` cover it.

### 10. [minor] [D] D1 step 0 "changes" `tests/test_box_inputs.py`, which D8 creates later

D "D1" Files: "change `tests/test_box_inputs.py` (step 0): add S09". The file does not exist today [measured: `ls
tests/`], and S08 (task D8) is its first test. D1 runs before D8 and S09 must pass before V2 (front 6.1 step 4).
Fix: "create `tests/test_box_inputs.py` (step 0) with S09; D8 adds S08 to it".

### 11. [minor] [A, D] Test id S09 names two tests

A.9: "**S09 golden file integrity** (S1 X10). A1." D.9.6: "S09 | R-D1 and K3 inputs ship". Front 6.1 step 4 says
"D's box-input test S09", but A's coverage rows (A1, A-G1) also list S09. Fix: rename D's test S10 (D.9.6, D1 step
0, D.9 coverage, front 6.1 step 4, front section 5).

### 12. [minor] [E] E.9 names a V4 item "E checks" that E.8 does not have

E.9 conventions: "All five run at V4 item "E checks" (E.8)". E.8 numbers the V4 steps 0-3 and has no item of that
name. Fix: "run at V4 step 2 (E.8)" (and "GP1-GP4" after item 1).

---

Totals: 2 blocking (items 1, 2), 3 major (items 3, 4, 5), 7 minor (items 6-12).
