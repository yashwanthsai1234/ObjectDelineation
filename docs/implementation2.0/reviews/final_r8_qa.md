# Final review round 8: QA lead and technical writer (ASD-STE100)

Target: `docs/implementation2.0/Implementation.md` (round 8 build: 1,608,894 bytes, 15,987 lines), built from
`docs/implementation2.0/*.md` by `assemble.py`. First I read the round-7 fix rows of every phase in `review-log.md`. Then I
read the text that those fixes changed and the Tests sections that touch rentals: front 6.1-6.4; A3 (`a3_measure.yaml`,
`a3_test.yaml`), A6 step 8, A.6, A.9.6; B2.2 F steps 4-5, B2.3, B2.4, B3.9, B.9.1 U51/U76, B.9.5, B.9.8; C0 files and
steps, C.6; D7 steps 5-6, D8, D.6, D.9 box rows; E.8, E.9 intro, E.9.5; DOC.6; CLEAN.6. I do not repeat items that
rounds 1-7 fixed. I also do not repeat the round-7 minors that are still open (r7 qa 3-9).

No REPO file was edited. No git write. No rental. No `assemble.py` run. No Python read `data/`.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | The round-7 text adds no new NVIDIA product claim, so I did no catalog lookup. It sent me to the box CUDA-build path: which plan file builds `_gpucheck` with nvcc, and in which order the plans run on each rental. That gave items 1 and 5. |
| engineering:code-review | I read the code behind each claim: `policy.py:420-490` (`chain`, `cmd_run`, `job`: a `steps:` plan runs as written, and a non-`bash` argv runs as `python -I`), `final_eval.py:46-56` (`--k` 8 and `--seed` 11 are the defaults, so `a3_measure.yaml` is right on K and seed), `checker_gpu.yaml`. Items 1, 2, 5. |
| engineering:debug | I traced each rental in launch-line order and asked "which build is in `/opt/venv` now?". V2 item 1b runs before the only V2 build step (item 1). V3 item 2b has no build when item 1 is skipped (item 5). |
| engineering:system-design | I matched each box gate to the plan step that produces its record. D-G14, D-G15 and D-G17 have no producer (item 3). G08b has no plan step (item 4). |
| engineering:architecture | "One owner per rule": A6 still writes a full argv for the C0 baseline arm that C0 owns, and the two argvs differ (item 2). |
| engineering:documentation | I checked that each test names its runner, its input files and its rental step (the B.9.5 and A.9.6 rule). G02b has an input that no part names (item 6). I also found stale numbers (items 8, 9). |
| huggingface-skills:huggingface-best | No model choice is in scope this round. No change. |
| engineering:testing-strategy | For each box test: runner, inputs, rental minutes. Items 3, 4, 6 come from this check. U-CF4 only parses `-m` steps, so it cannot see the nsys-wrapped arms (item 2). |
| mattpocock-skills:tdd | I looked for gates whose test cannot fail. B04/B06 have no runner, so D-G14/D-G15/D-G17 cannot fail in V3 (item 3). |
| Web | Loaded WebSearch, WebFetch and Exa. I needed no web fact: every item turns on repo text and code only. |

Measured or read during this review (2026-11-08):
- `cmd_run` and `job`: a `steps:` plan returns its steps as written, with `{rd}`/`{data}` filled in [code:
  experiments/phase2/vast/policy.py:434-435]. A step that does not start with `bash` runs as
  `/opt/venv/bin/python -I <argv>` [code: policy.py:486]. Steps are joined with `&&` [code: policy.py:490].
- `final_eval` defaults: `--k` 8, `--seed` 11, `--split val`; `--ckpt`, `--view`, `--train`, `--out` are required
  [code: experiments/phase2/final_eval.py:48-55]. `runs/r5/eval/eval_val/final_eval.json` has `K` 8, `seed` 11
  [measured: `python3 -c json.load`].
- `Implementation.md` is 1,608,894 bytes [measured: `ls -l`].
- grep over every phase source: no task outside D writes `gc.stats()`, `sync_alloc`, `pinned_grows`,
  `max_memory_reserved` or `runs/d/topology.json`. `G08b` / `test_g08b_smoke` appears only in B.9.5 (two rows) and in
  the coverage maps.

Severity follows the task rules. Count: **1 blocking, 5 major, 4 minor**.

---

## Blocking

### 1. [blocking] [00_front, A] V2 item 1b (`a3_measure.yaml`) runs before the only V2 CUDA build, and its second step asserts the CUDA build, so item 1b always fails

Problem.
- 00_front "### 6.3", V2 table, item 1: "... `env_check`, toolchain and CUDA build (`c_speed.yaml` step 0 =
  `b1_rows.yaml` step 0, 6.2) ...". So the V2 build is step 0 of `c_speed.yaml`.
- 00_front "### 6.3", "V2 launch lines":
  ```
  run.sh <run> experiments/phase2/configs/a3_measure.yaml    # item 1b, only under the 6.2 fallback
  run.sh <run> experiments/phase2/configs/c_speed.yaml       # items 2-3 (3c, 3d only if approved); wait for EXIT
  ```
  So item 1b runs before `c_speed.yaml`, that is, before the build.
- A "#### A3", plan file `a3_measure.yaml`: step 2 is
  `[-c, "import objdel._gpucheck as g; assert g.build() == 'cuda', 'CPU build: b1_rows.yaml step 0 did not run'"]`,
  and the text says "`env_check.sh` and the CUDA build already ran in `b1_rows.yaml` on the same box (front 6.2), so
  `a3_measure.yaml` does not repeat them." That holds in V1 only.
- Code: `sync_up` installs the package with `uv pip install ... -e .` [code: policy.py:399], and CMake builds the CPU
  `_gpucheck` unless `OBJDEL_CUDA=ON` [code: CMakeLists.txt:15]. The front itself says this in 6.2 "Why step 0" and
  in the V4 text ("Without it `sync_up` leaves the CPU build").

Result. Under the 6.2 fallback (V1's R5 process gave no output), an engineer who follows the V2 launch lines runs
`a3_measure.yaml` on the CPU build. Step 2 raises `AssertionError: CPU build: b1_rows.yaml step 0 did not run`, the
`&&` chain stops [code: policy.py:490], and item 1b gives no A3 answers. Item 1b exists only to make the Gate B1
paired baseline, the A-G3 val answers, the Mac half of C-G-EQ2/3 and A8's answer table (A.7 row "The R5 process of V1
... gives no output"). Item 1's row says the build is part of item 1, but no item-1 command builds. Two parts both
claim the order, and the engineer cannot proceed with item 1b.

Exact fix (pick one; the first is smaller).
- A3: make `a3_measure.yaml` step 2 the full build, as `a3_test.yaml` step 1 already does:
  `[bash, -c, "bash experiments/cuda_checker/run_box.sh --build-only && /opt/venv/bin/python -I -c 'import
  objdel._gpucheck as g; assert g.build() == \"cuda\"'"]`. In V1 it adds about 3 min (front 6.2 item 1 "build 3"), so
  front 6.2 item 4 becomes 9.9 min / $0.11 and the V1 total 196.0 min. Or keep the assert in V1 and add the build
  step only to a V2 copy `a3_measure_v2.yaml`.
- Or front 6.3: move the `a3_measure.yaml` line after a build-only call. The build must not depend on `c_speed.yaml`.
  Add `run.sh <run> experiments/phase2/configs/box_build.yaml` (one step = `b1_rows.yaml` step 0) right after the sync
  line. Item 1 then names that file, not "`c_speed.yaml` step 0".
- Either way, add to front 6.1 step 4d (the plan-file grep test): in every rental, the first plan that runs on the
  box builds the CUDA checker, or runs after a plan that does.

---

## Major

### 2. [major] [A, C, 00_front] A6 still gives the C0 baseline arm a full argv on the centred `prepared` root with `--steps 17`; C0 (round 7) needs the B1 rows root and `--steps 20`; U-CF4 does not parse nsys-wrapped arms

Problem.
- A "#### A6" step 8 (the `nsys_capture.sh` text): "The A6 baseline is the C0 baseline arm of V2 item 3 (R5 code and
  settings, `--shared-prefix`, 17 steps, ...): `nsys_capture.sh {rd} c0_base 16 17 -- /opt/venv/bin/python -I -m
  experiments.phase2.train_rl --model qwen3-vl-2b --ckpt {rd}/ckpt --train {data}/prepared --out {rd}/c0_base
  --steps 17 --save-steps 1000 --shared-prefix`. C0 owns the arm and its minutes".
- C "### C0" Files: "**Every C arm trains on the B1 rows root** ... `--train {data}/rows/lidc-train-b1
  --rows-revision <hex40> --rows-manifest-sha256 <hex64>` ... B1.8 deletes the 'without schedule' branch of
  `train_rl` ... `train_rl` trains only whole blocks (B1.3: `--steps % schedule.block_steps == 0`, `block_steps`
  20). So each 17-step arm (baseline, 3c, 3d, C9) runs `--steps 20`". C0 step 3 then says: "The one command is A6's
  script (A6 step 8): `nsys_capture.sh {rd} <arm> 16 17 -- /opt/venv/bin/python -I -m experiments.phase2.train_rl
  ...`". The only full argv for the baseline arm is A6's, so a reader who follows C0 step 3 copies it.
- C0 Files, U-CF4: "parses every `-m experiments.phase2.<module>` step of `c_speed.yaml` with that module's own
  argparse parser". C0 step 3 says "Both the baseline and the fast arm are profiled the same way", so both are `bash`
  steps (`nsys_capture.sh ... -- ... train_rl ...`), not `-m` steps. A `bash` step runs as written [code:
  policy.py:486]. So a check limited to `-m` steps does not parse the two main arms.
- A "### A.9.6" G06 and the A.3 diagram (A.md line 185) also say "17 steps", "train_rl 17 steps".

Result. The A6 argv leaves out both rows flags and uses `--steps 17`. On the V2 tree (B1.8 landed) `train_rl` refuses
it: there is no branch without a schedule, and 17 is not a multiple of 20. The baseline arm fails on the paid box
(V2 item 3, $0.49 for the item), and C-G-SP1 is relative to that arm, so it cannot be computed. U-CF4 as written
cannot catch it. The front 6.1 step 4 wording ("every `train_rl` step names the rows root") may catch it, but only
if the test reads `bash` strings too. The text does not say so.

Exact fix.
- A6 step 8: replace the argv with
  `nsys_capture.sh {rd} c0_base 16 17 -- /opt/venv/bin/python -I -m experiments.phase2.train_rl --model qwen3-vl-2b
  --ckpt {rd}/ckpt --train {data}/rows/lidc-train-b1 --rows-revision <hex40> --rows-manifest-sha256 <hex64> --out
  {rd}/c0_base --steps 20 --save-steps 1000 --shared-prefix --profile-window 16 17` (values of `b1_rows.yaml`), or
  replace it with "C0 files hold the argv". Change "17 steps" to "17 steps read, `--steps 20`" in A6 step 8, in the
  `timeout` line ("model load + 20 steps"), in G06 and in the A.3 diagram.
- C0 U-CF4 and front 6.1 step 4: "every step whose argv holds `experiments.phase2.train_rl` (also inside a `bash -c`
  string after `nsys_capture.sh ... --`, split with `shlex.split`) is parsed with `train_rl`'s parser and names the
  rows root, both rows flags and a `--steps` that is a multiple of 20". Add a planted case: an nsys-wrapped arm with
  `--train {data}/prepared` fails the test.

### 3. [major] [D, B] D-G14, D-G15 and D-G17 (box tests B04, B06) have no runner: no plan step and no code task writes their records in V3

Problem.
- D "## D.6": D-G14 "V3, first 6 steps of the B3 run: the reward path's `arena_bytes` and
  `torch.cuda.max_memory_reserved()` go into the run JSON ... the B3 run's val eval ... records the `eval_sets` time,
  `arena_bytes` and peak memory". D-G15 "Box, V3 B3 run, steps 2-6: the `gc.stats()` delta of `sync_alloc` and of
  `sync_free` is 0 ... `pinned_grows == 0` over the whole run". D-G17 "(a) ... the count, sampled every 5 s, stays
  equal to it ... `runs/d/topology.json` records which method ran. (b) Each pool worker ... logs
  `torch.cuda.is_initialized()` and `gc.stats()["contexts"]`".
- D "### D.9" box rows: B04 "B3 run, first 6 steps: reward-path `arena_bytes` and `max_memory_reserved` logged ...";
  B06 "... result written to `runs/d/topology.json`". Neither row names a file, a command or a plan step. The B.9.5
  and A.9.6 rule is that "Each box test names its rental item and its command". B04 and B06 do not.
- D "#### D7" step 5 says only "In V3, record ..." and names no writer.
- B "#### B3.9" and B "### B.9.3" I22b list every element of `b3_v5.yaml` ("same start checkpoint ..., LocGuard and
  BudgetStop on, no test eval, `train_rl --out {rd}/b3-v5`"). They hold no counter logging, no `nvidia-smi` sampler
  and no `topology.json` step. D.8 bills the gates at "0 extra" minutes inside the B3 run.
- grep over all phase sources: no task in A, B, C or E makes `train_rl` or `final_eval` log `gc.stats()`,
  `sync_alloc`, `pinned_grows` or `max_memory_reserved`. No task writes `runs/d/topology.json`.

Result. In V3 the B3 run and its val eval run, and none of these records exist after `pull.sh`. The box halves of
three D gates cannot pass or fail, and an engineer must invent the logging code and a sampler during the paid
session. D-G17 is the only box check of D7 step 6 ("one process for device eval"), the change that stops up to 16
CUDA contexts beside vLLM. D-G15 is the only box check of D5c's pool.

Exact fix.
- D7 step 5 (or D5c): name the code. `train_rl` gains `--engine-stats N`: for steps 1..N it logs
  `engine/arena_bytes`, `engine/sync_alloc`, `engine/sync_free`, `engine/pool_allocs`, `engine/pool_frees` (deltas of
  `gc.stats()` per step) and `torch.cuda.max_memory_reserved()` into `rl.json`, and it logs `engine/pinned_grows` at
  the end. `final_eval` writes `eval_sets` seconds, `arena_bytes` and peak memory into `final_eval.json["engine"]`.
- D-G17 (a): a script `experiments/phase2/gpu_topology.py --out {rd}/d/topology.json --every 5 -- <final_eval
  argv>` starts the sampler, runs the eval as a child, and writes the counts and the method. Or the B3 eval step of
  `b3_v5.yaml` is wrapped with it.
- B3.9 and I22b: `b3_v5.yaml` adds `--engine-stats 6` to its `train_rl` step, and its eval step uses the wrapper. I22b
  asserts both. B04/B06 rows: "Run place: V3 item 4; Mac check after `pull.sh`:
  `tests/test_d_records.py::test_b04`, `::test_b06` read `runs/<v3>/b3-v5/rl.json`, `eval_*/final_eval.json` and
  `d/topology.json`".

### 4. [major] [B, D, 00_front] G08b (the smoke gate before the B3 rewards use B3.1's `sum_x`/`sum_y`) has no plan step in any V3 plan

Problem.
- B "### B.9.5": G08b is the "Smoke gate of the first rental after B3.1 merges (V3, `test_g08b_smoke`): G01 parity on
  1,000 windows (with `sum_x`, `sum_y`) and, only for the optional fp16 arm, the `vllm_dtype` engine check". Its
  runner row: "V3, the step after the CUDA build of the first plan after B3.1 merges | `python -I -m pytest -m gpu
  tests/test_train_fast.py::test_g08b_smoke`".
- The V3 plans, in launch order (00_front "### 6.3"): `d_gpu_test.yaml` (item 1, "Skippable"), `b2_fp16.yaml` (2b),
  `a3_test.yaml` (3), `b3_v5.yaml` (4).
  - D "### D8" and S05: "`d_gpu_test.yaml` has exactly 2 steps, both start with `bash`, `-c`; step 2 names
    `tests/test_a_gpu.py::test_g02b_builds` and `::test_g03_device_bytes`". There is no room for G08b.
  - A "#### A3": `a3_test.yaml` has exactly a build step and a `final_eval` step.
  - B "#### B2.3" / I22b: `b2_fp16.yaml` has no build and no smoke step.
  - B "#### B3.9" / I22b list `b3_v5.yaml`'s steps: "step 0 equal to `b1_rows.yaml` step 0" and then `train_rl`, with
    no smoke step.
- 00_front 6.3 V3 item 4 minutes are "3 + 6 + 260 x step time + 10.4", with no smoke minutes.
- The B.9.5 runner also writes `python -I -m pytest`, but a plan step must start with `-m, pytest` [code:
  policy.py:486].

Result. The engineer must guess where G08b goes. If it never runs, the B3 run (the largest V3 item, $1.57) trains on
reward v5's centroid sums with no box parity check of `sum_x`/`sum_y` against the CPU build. That breaks the "no spend
before scan" rule that G08a holds for V1. G01 runs only in R-D2 (item 1), which is skippable.

Exact fix.
- B3.9 and I22b: `b3_v5.yaml` step 1 (after the step-0 build) is `[-m, pytest, -q, -p, no:cacheprovider, -m, gpu,
  tests/test_train_fast.py::test_g08b_smoke]`. I22b asserts it and checks that the node id collects 1 item. Name
  G08b's 1,000-window input (for example, the first 1,000 rows of `{data}/rows/lidc-train-b1/lidc/train.jsonl` and
  their caches) and say that it is in the V3 payload.
- B.9.5 runner row: "V3 item 4, `b3_v5.yaml` step 1" and the `-m, pytest` argv. If the owner funds 2b, `b2_fp16.yaml`
  gets the same step for the `vllm_dtype` part (see item 5).
- 00_front 6.3 V3 item 4: add the smoke minutes (< 5 [estimate], as V1 item 1) and recompute B-S2.

### 5. [major] [B, 00_front] `b2_fp16.yaml` (V3 item 2b) has no CUDA build step, so when item 1 is skipped, the fp16 arm's rewards score on the CPU checker

Problem.
- 00_front "### 6.3" V3 table: item 1 (R-D2, the only V3 item before 2b that builds) has stop rule "Any failure:
  stop; bisect on the Mac. Skippable".
- A "#### A3", `a3_test.yaml`: "its first step is the CUDA build of `b1_rows.yaml` step 0, because no earlier V3 item
  is sure to build it". The same holds for item 2b, which runs before item 3.
- B "#### B2.3": "Config `experiments/phase2/configs/b2_fp16.yaml`: start `rl-full-3ep`, the B1 rows root, seed 0,
  `--numerics fp16-seqmask --lr 5e-5 --steps 60 --save-steps 60 --shared-prefix`. No eval." I22b lists its contents
  with no step 0. `b3_v5.yaml` does get one ("Step 0 of `b3_v5.yaml` is the CUDA build of `b1_rows.yaml` step 0").
- Code: `sync_up` leaves the CPU build [code: policy.py:399; CMakeLists.txt:15]. A `steps:` plan runs only what it
  lists [code: policy.py:434-435].

Result. If item 1 is skipped or fails (both are allowed), the paid 60-step fp16 arm computes every reward on the CPU
checker. That breaks the owner rule "CUDA is a must for the checker, the scorer AND the polygonizer" (BRIEF section
3). It also slows the CPU phase that B2.3 prints as s/step. G08b's `vllm_dtype` check for this arm has no step either
(item 4).

Exact fix.
- B2.3 and I22b: `b2_fp16.yaml` step 0 = `b1_rows.yaml` step 0 (build, assert, `tests/test_gpu_scoring.py`), and
  I22b asserts "step 0 equal to `b1_rows.yaml` step 0" for `b2_fp16.yaml` as it does for `b3_v5.yaml`.
- 00_front 6.3 item 2b: "60 x step time + 6 + 3 (build)": 36.5 / 32.7 min, $0.41.

### 6. [major] [A, D] G02b's golden scores `model` of `ceiling.item` for 50 val items, but no part names the model texts, and the saved answers are not on the V3 item 1 box

Problem.
- A "### A.9.6", G02b: "`make_a_goldens.py` ... writes `tests/goldens/a_g02b.json`: for 50 val items, `doctor` and
  `model` of `ceiling.item`; ... The box test scores the same items on the CUDA build and compares with `==`. ...
  the 50 val items come from `prepared/lidc/val.jsonl` in the bundle."
- A "#### A2": "`item(cache, texts, truncated, W, H, name)` ... With `texts=None` it returns the doctor side only."
  So `model` needs K answer texts per item.
- D "### D8", R-D2 payload: "Not in the payload: the saved val answers". D8 step 2 gives G02b's inputs as only
  "`prepared/lidc/val.jsonl` and `prepared/caches/lidc/` (R-D2 payload item 2) and the tracked golden
  `tests/goldens/a_g02b.json`".
- No part says which 50 items are used (first 50? a seed?) or which texts feed `model`.

Result. The natural reading is A3's answers (`runs/a3/eval_val/answers.jsonl`, the ceiling CLI's `--answers`). Those
are on the Mac but not on the V3 box, so the box half fails on a missing file in step 2 of `d_gpu_test.yaml`. A
different choice (target texts, synthetic answers, a different 50) is a guess, and the box and the golden must make
the same guess to compare with `==`.

Exact fix. A.9.6 G02b: "the 50 items are the first 50 rows of `prepared/lidc/val.jsonl` sorted by `item_id`; the
`model` side uses K = 8 texts per item stored in the golden itself (`texts` field: the item's target texts, cycled
to 8, plus 2 adversarial texts from the 500-polygon set), so the box reads no `runs/` file". `make_a_goldens.py`
writes the texts and their SHA-256. The box test asserts the SHA-256 before it scores.

---

## Minor

### 7. [minor] [C] G-EQ6 reference values come from R5's centred rows, but the round-7 arms train on the B1 rows root

C "## C.6", G-EQ6: "Reference values, R5 shared arm step 1: IoU_tc 0.5467, legal 0.8984, mean length 117.89, grad
norm 0.003264, logp difference 0.009515". C0 now says "Every arm reads the same first block of rows" of the B1 rows
root, so step 1 sees different prompts from R5. The gate (fast vs baseline, same host) is still right. A reader may
compare against these values and think the arm is broken. G-SP1's "host not comparable to R5" note (23.43 s ± 10 %)
has the same issue: shifted rows may lengthen the decode tail (front 6 "Step time, base").
Fix: G-EQ6: "Reference values (R5, centred rows; for information only, the arms read the B1 rows root, so step 1
differs)". G-SP1: add "(rows root: R5 is a different batch)" to the note.

### 8. [minor] [00_front, B] V3 item 2a (B2.3b box part) has no launch line and no full argv

00_front "### 6.3" V3 launch lines: "`run.sh <run> experiments/phase2/configs/b2_fp16.yaml  # item 2b, only if
approved (item 2a uses the B2.3b box command)`". B.9.5 G04 runner: "`python -I -m experiments.phase2.b2_attribution
--box`". No line gives the adapter path (`{rd}/b1/final_adapter`), `--ckpt {rd}/ckpt` or the output folder, and no
`run.sh` line exists for item 2a. `cmd_run` accepts the module form [code: policy.py:476-477], so the fix is text only.
Fix: front 6.3 adds `run.sh <run> experiments.phase2.b2_attribution --box --ckpt {rd}/ckpt --adapter
{rd}/b1/final_adapter --out {rd}/b2_attr   # item 2a, only if B2.0 error < 0.2`, with the box paths written out
(`run.sh` fills `{rd}` only in config files), and B2.3b names the same flags.

### 9. [minor] [DOC] DOC-G13's "today" cell gives a stale size for `Implementation.md`

DOC "## DOC.6", DOC-G13: "`Implementation.md` 1,480,096 bytes [measured: `ls -l`, 2026-11-08]". The round-8 build is
1,608,894 bytes [measured: `ls -l`]. The 2 MB limit still holds, with 391 KB of margin. Fix: write the current size,
or write "measured at DOC3 step 3".

### 10. [minor] [A] A6's `timeout` budget is written for 17 steps

A "#### A6" step 8: "`timeout` (20 min = model load + 17 steps at <= 31 s + nsys flush [estimate])". The baseline arm
now trains 20 steps (C0). At 31 s that is 10.3 min of steps instead of 8.8. Add about 6 min of start-up (front 6) and
the nsys flush, and the 20-min budget is tight. Fix: "20 min = start-up about 6 + 20 steps at <= 31 s (10.3) + flush",
or `NSYS_TIMEOUT=1500` for the baseline arm as C0 sets 1800 for the fast arm.

---

Checked and found correct this round (no item): the round-7 qa-1 fix (B2.4 rule 0, U51 and the B.9.8 sketch now take
`preset` and `passed` from a full `f_check.json` record; the ratios in the sketch are right: (0.0137 - 0.0096) /
(0.0176 - 0.0096) = 0.51, (0.0140 - 0.0096) / 0.008 = 0.55, (0.0130 - 0.0096) / 0.008 = 0.425, 0.0115 - 0.0096 =
0.0019 < 0.002). The round-7 qa-2 fix: D8 `d_gpu_test.yaml` step 2 is a `bash -c` step with a build assert, so
`cmd_run`/`job` run it as written [code: policy.py:486], and `{data}` is filled by `str.format` [code: policy.py:435].
The E.8 V4 step 2 argv starts with `-m, pytest`. `a3_measure.yaml` K 8 and seed 11 come from the `final_eval`
defaults [code: final_eval.py:54-55]. The V2 required total of 98 min = 42 + 13 + 43 = $0.10 + 98 x $0.01133 = $1.21.
