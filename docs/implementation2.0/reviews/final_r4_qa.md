# Final review round 4: QA lead and technical writer (ASD-STE100)

Target: `docs/implementation2.0/Implementation.md` (round 4; 1,457,224 bytes, 14,787 lines [measured: `ls -l`, `wc -l`]),
built from `docs/implementation2.0/*.md` by `assemble.py`. I read the Tests sections (A.9, B.9, C.9, D.9, E.9, DOC.9,
CLEAN.9), every gate table, front sections 1, 2.5 and 6, and the code each claim names. I checked every item against
`review-log.md` (rounds 1-3 and the round 4 prep pass); fixed items are not repeated. No REPO file was edited. No git
write. No rental. No Python read `data/`.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No new catalog lookup was needed: no new NVIDIA product claim. It pointed me at the box toolchain: the image has no nvcc (`run_box.sh` header), which is the root of item 2. |
| engineering:code-review | Read the code behind each launch and test claim: `policy.chain`/`cmd_run` (`&&` chain, policy.py:425-490), `preflight.load/hashes/refusal` (preflight.py:49-95), `bundle.members` (`LIDC_EVAL`, bundle.py:26, 53-62), `train_rl.py:10, 59` (merged `final/`), CMake `OBJDEL_CUDA` default OFF. Items 2, 3, 5, 6, 8. |
| engineering:debug | For each box test I traced: which build is on the box, which hook skips it, which exit code the plan step gets, and whether `&&` stops the chain. Two pytest exit codes measured (items 2, 3). |
| engineering:system-design | Put V1, V2, V3, the tag `v1-b1` and the numerics preset of each run on one line; found the F/B3 contradiction (item 1). |
| engineering:architecture | ADR-0013 rule (b) in B against the front 2.5/6.3 text and 99_back decision 16 (item 1). |
| engineering:documentation | One word per thing ("final adapter" vs merged `final/`, item 8); long-sentence count (item 15). |
| huggingface-skills:huggingface-best | Model choice did not change. The HF part of this review is what `hf_push --rl` uploads to revision `b1-rows` (item 8). |
| engineering:testing-strategy | Matched each Tests row to its runner, file, marker and rental (items 3, 4, 7, 10-13). |
| mattpocock-skills:tdd | Looked for tests that pass by construction or cannot run: the V1 GPU smoke passes vacuously when every test skips (item 2); I22 asserts the stale plan (item 4). |
| Web | No new external claim needed a web check: the pytest exit codes were measured locally instead (below). |

Measured during this review (2026-11-08, scratch dir `(planning probe, not kept: qa4probe)`):
- `pytest -k "G02 or G08"` on a file whose test names do not contain the ids: `1 deselected`, **exit 5**.
- `pytest -k "G02 or G08"` on a file where the only selected test skips: `1 skipped`, **exit 0**.
- `grep -rn sum_x cuda/ src/objdel`: 0 hits (B3.1's sums do not exist before B3.1).
- `runs/r5/eval/eval_test/` holds `final_eval.json`, `items.jsonl` (scores only), `trainer/`; no answers file.
- Long sentences, rule of DOC.6 (join wrapped lines, split list items, skip code and tables, > 25 words): 1,314 of
  7,204 sentences (18.2 %), my script `(planning probe, not kept: qa4probe)/count.py`.

Severity follows the task rules. Count: **4 blocking, 4 major, 7 minor**.

---

## Blocking

### 1. [blocking] [00_front, B, C, 99_back] Two rules for the B3 preset and for the V2 F arm

Problem. Two parts both claim to be the rule.
- 00_front "### 2.5 One RL step after all phases": "The seed value is 0.531 [measured: B.1 item 5], so F is the
  expected B3 preset."
- 00_front "### 6.3 V2 ...": "B2.4 rule 0 is expected to pick F (section 2.5). So the F arm (item 2b) runs before
  the C arms, and the C0 baseline arm and the fast arm run under the preset that item 2b leaves." Item 2b stop rule:
  "F fails: items 3 and B3 use the preset that rule 0 falls back to".
- 99_back "### RL settings" decision 16: "its 15-step arm is rental V2 item 2b ... if it fails, B3 uses the preset
  that rule 0 falls back to".
- C "## C.6. Acceptance gates" (C.md:1059): "V2 item 2b (B's F check) runs before the arms ... The arms then use
  `lora-seqmask` if the F check passed".
- Against that, B (round 3 fix, ADR-0013 rule (b)):
  - B "#### B2.4 Decide": "The B3 run keeps B1's preset `bf16-seqmask` ... B2.4 picks the preset for the runs after
    Gate B3."
  - B "#### B2.2 Numerics presets": "box check in the first run that uses F, after Gate B3 under rule (b), B2.4".
  - B "### B.8.4 Dependencies on other phases", row "Faster step": "Under ADR-0013 rule (b) the B3 run uses
    `bf16-seqmask` and C5's merged push, so a C-G-SP1 pass measured on the merged push carries over to B3, and no F
    arm is needed in V2."
  - B "### B.6.2 Run gates", Gate B3: "the B3 run used B1's preset `bf16-seqmask` (ADR-0013 rule (b))"; U51: "the
    B3 config is `bf16-seqmask` whatever the rule returns".
- Result. An engineer cannot tell (a) whether V2 spends 13 min ($0.15) on item 2b, (b) which preset the C arms use,
  so whether C-G-SP1 is measured on the push that B3 uses, (c) which `--numerics` the V3 B3 run gets. If the C arms
  run under F (front) and B3 runs `bf16-seqmask` (B), the C-G-SP1 result does not carry over to B3: GPU money for
  nothing.

Exact fix (B's rule (b) is the newer decision and keeps Gate B3 clean; make the other parts follow it):
- 00_front 2.5: replace "so F is the expected B3 preset" with "B2.4 picks the preset for the runs after Gate B3;
  F is the expected preset there. The B3 run keeps `bf16-seqmask` (ADR-0013 rule (b))."
- 00_front 6.3: delete "Numerics before speed" bullet and row 2b, or move row 2b to V3 as "after Gate B3, the first
  15 steps of the next run". The C arms run under `bf16-seqmask`. Recompute "Required for C adoption" (103 - 13 =
  90 min, $1.11 [estimate: $1.26 - $0.15, front 6.3 rows]).
- 99_back decision 16: "F is decided by B2.4 rule 0 for the runs after Gate B3; B3 uses `bf16-seqmask`."
- C.6 (C.md:1059), C intro (C.md:34-37), C5 "Under option F" and C.8.2: the arms use `bf16-seqmask`; delete "V2
  item 2b ... runs before the arms".

### 2. [blocking] [00_front, B, A, C] The V1 plan has no CUDA build, so V1 scores on the CPU build and the GPU smoke passes by construction

Problem.
- The box builds `objdel._gpucheck` for the CPU. `sync_up` runs `uv pip install ... -e .` [code:
  experiments/phase2/vast/policy.py:399]. CMake builds the CPU version unless `OBJDEL_CUDA=ON`
  [code: CMakeLists.txt:15 `option(OBJDEL_CUDA ... OFF)`]. Only `run_box.sh` installs nvcc, "(the image has none)"
  [code: experiments/cuda_checker/run_box.sh:5], and runs `make cuda`. D says the same: D "### D8": "Why one shell
  entry: the box image has no nvcc."
- The 770-step full run had that step: `rl_full.yaml:25` runs `run_box.sh` before training [code]. The V1 plan
  does not. 00_front "### 6.2 V1": "`b1_rows.yaml` (B1.5 writes the commands): (1) `env_check.sh`; (2) B's GPU
  smoke ...; (3) `train_rl` ...; (4) one `final_eval` process". B "#### B1.5 Train on the ordered rows" lists the
  same four steps. `a3_measure.yaml` (A "### A3") has no build step either.
- Consequence 1, owner rule (BRIEF section 3, "CUDA is a must for the checker, the scorer"): every reward and every
  V1 eval is scored by the CPU build on the host [code: src/objdel/reward/score.py:75, 115 call `_gpucheck`].
- Consequence 2, a gate that cannot fail. A "## A.9. Tests": "`pytest_collection_modifyitems` marks every `gpu`
  test skipped unless `objdel._gpucheck.build() == "cuda"`". On the V1 box `build()` is `"cpu"`, so G02 and G08
  skip. pytest then exits 0 [measured: all-skipped run, exit 0], and the chain goes on. B "### B.9.5 GPU tests":
  "Every new GPU path runs before money is spent" is then false.
- The same skip hits C's V2 box pytest step (C "### C.9.4": "`"$PY" -I -m pytest -m gpu tests/test_c_gpu.py`, a step
  of `c_speed.yaml`"; no build step in `c_speed.yaml`) and E's GP1/GP2 in V4. Only D8's `--gpu-check` asserts the
  build and fails on skipped tests.

Exact fix.
- B1.5 and front 6.2: make step 1 of `b1_rows.yaml`
  `[bash, -c, "bash experiments/cuda_checker/run_box.sh --build-only"]` (the build lands in `/opt/venv`, so later
  steps and `a3_measure.yaml` use it). Add its minutes to item 1 (build about 3 min [estimate: D "## D.8", R-D2
  row]) and recompute the V1 totals.
- Same first step in `c_speed.yaml` (C.8.2 / front 6.3 item 1), the B3 plan `b3_v5.yaml`, and the V4 plan (E.8).
- A.9 root hook: on a box (`OBJDEL_RUN_DIR` under `/workspace/runs`) a `gpu` test on a CPU build fails with
  "CPU build on a GPU box: run run_box.sh --build-only first" instead of skipping. Keep the skip on the Mac.
- `env_check.sh`: add a check `gpucheck_build == "cuda"`, so `train_rl` never starts on the CPU build.
- Tests: B I22 asserts the build step is step 1 of `b1_rows.yaml`; a Mac unit test of the hook with a faked box
  path asserts "fail", not "skip".

### 3. [blocking] [B, 00_front] The V1 GPU smoke step either selects no test or fails on code that is not on the tag

Problem. 00_front "### 6.2 V1" and B "#### B1.5": step 2 is `pytest -m gpu tests/test_train_fast.py -k "G02 or G08"`.
The `steps:` entries are joined by `&&` [code: policy.py `jobs = " && ".join(...)` in `cmd_run`], so a non-zero exit
stops `train_rl` and the evals of item 3.
- (a) File and name. B "### B.9.8 Pytest sketches" puts B2.1's tests in `# tests/test_train_fast_numbers.py (U35,
  P13; B2.1)`, and B's test names carry no id (`test_is_ratio_geo_mask_values`, `test_sums_cpu_gpu_bit_equal`).
  B "#### B2.1" says "Tests in `tests/test_train_fast.py`". No rule says the G02/G08 test names contain "G02"/"G08".
  If they do not, `-k` deselects all and pytest exits 5 [measured]. V1 then pays item 1 ($0.44) and gets no B1 run.
- (b) G08 content. B "### B.9.5": "G08 | Preflight smoke gate (< 5 min): G01 parity on 1,000 windows, G02 no-sync
  run ..." and "G01 | CPU build and CUDA build give bit-equal results on every field ..., including `sum_x`,
  `sum_y`". `sum_x`/`sum_y` come from B3.1. 00_front "## 1", order rule 2: the tag "holds no C, D, B2.2-B2.4, B3 or
  kernel code". `sum_x` has 0 hits today [measured: grep]. So, once item 2 gives the box a CUDA build, G08 fails on
  `v1-b1` and stops the chain.

Exact fix.
- B.9.5: "G08 on `v1-b1` = G02 plus the `score_masks` old-field parity (CPU vs CUDA, 1,000 windows, goldens of the
  previous build). The `sum_x`/`sum_y` part of G01 joins G08 only after B3.1." Split G08 into G08a (V1) and G08b
  (after B3.1).
- Name the node ids, not `-k`: `python -I -m pytest tests/test_train_fast.py::test_g02_is_summary_no_sync
  tests/test_train_fast.py::test_g08a_smoke`. Put G02, G08a and B2.1's U35/P13 in one file (`tests/test_train_fast.py`;
  rename the sketch header) and add the B.9.0 rule "box test names start with `test_<id>_`" (A already does this:
  `test_g01_answers`).
- Add to B I22: "the step-2 node ids exist in the file
  (`pytest --collect-only` on the Mac lists exactly 2 items)". A collect-only check catches exit 5 at $0.

### 4. [blocking] [B] Test I22 asserts a stale `b1_rows.yaml` ($1.90 stop, R5 val-edge step)

Problem.
- B "### B.9.3 Integration tests", I22: "`b1_rows.yaml` has exactly the four steps of B1.5 in that order
  (`train_rl` with ... `--budget-usd 1.90`, ...; B1 on val shifted; B1 on val edge; R5 on val edge)".
- B "#### B1.5": the four steps are `env_check`, the GPU smoke, `train_rl ... --budget-usd 2.05`, and one
  `final_eval` process with two views. "R5 on val edge is not a `b1_rows.yaml` step: it is front item 4, in A3's
  `a3_measure.yaml`." 00_front "### 6.4 Stop rules" rule 5: "V1: $2.40 - $0.35 = $2.05".
- An engineer who makes I22 pass writes the wrong plan. At $1.90 the stop trips at session minute
  (1.90 - 0.10) / 0.01133 = 158.9, so at 27.5 s the run ends near step (158.9 - 30) x 60 / 27.5 = 281, not 300
  [estimate: front 6.2 formula]. R5 on val edge would also run twice (6.9 min). An engineer who follows B1.5 gets a
  red I22, and `make test` blocks the tag (order rule 2).

Exact fix. I22 text: "`b1_rows.yaml` has exactly the steps of B1.5 in that order: (1) the CUDA build (item 2 of this
review), (2) `env_check`, (3) the GPU smoke node ids, (4) `train_rl` with `--ckpt {rd}/ckpt`, `--numerics
bf16-seqmask`, beta 1e-3, `--lr 5e-5`, `--rows-revision`, `--rows-manifest-sha256`, `--budget-usd 2.05`, steps a
multiple of 20, `--save-steps 50`; (5) one `final_eval` with `--view lidc-val-shifted` and `--view lidc-val-edge`.
No R5 step; no `--split test`." Read the budget from one constant (front 6.4 rule 5) in the test, not a literal.

---

## Major

### 5. [major] [B] Preflight's "ROWS root" check cannot find the local rows root from the plan's box path

Problem.
- B "#### B1.3 Coverage and quality gates", change 2: "`preflight` reads that plan, takes the `train_rl` stage's
  `--train` root and its `--rows-manifest-sha256`, requires the local root's `manifest.json` SHA-256 to equal the
  flag". A missing manifest "is a failed check".
- B "#### B1.5", step 3: `--train {data}/rows/lidc-train-b1` (a box path; `{data}` = `/workspace/data/lidc`
  [code: policy.py `chain`]).
- B "#### B1.4": the local root is `--out data/views/lidc-train-b1`; HF holds it under `rows/lidc-train-b1/`.
- No text maps `{data}/rows/lidc-train-b1` to `data/views/lidc-train-b1`. A literal read finds no manifest, so the
  V1 stamp fails; a guessed map can check a different folder than the one uploaded. I19 does not pin the map.

Exact fix. B1.3 change 2: "preflight maps the box path `{data}/rows/<name>` to the local `data/views/<name>`
(one function `preflight.local_rows_root`)." Or move B1.4's `--out` to `data/rows/lidc-train-b1` and map `{data}` to
`data`. Add an I19 case: the plan's literal `--train {data}/rows/lidc-train-b1` resolves to the fixture root.

### 6. [major] [00_front] Stop rule 7 tells the engineer to move out views that the bundle and the stamp require

Problem.
- 00_front "### 6.4 Stop rules" rule 7: "Nothing from `data/` leaves the Mac except ... the LIDC eval views named in
  the rental's table ... so a view that is not in the table is moved out first (6.1 step 4c)." The V1 table names
  only `lidc-val-shifted` and `lidc-val-edge`.
- `bundle.members` dies when any of `LIDC_EVAL = ("lidc-val-shifted", "lidc-val-rescaled", "probe-precision",
  "probe-topology")` is missing [code: experiments/phase2/vast/bundle.py:26, 58-60].
- `preflight.hashes` reads the manifests of `lidc-val-rescaled`, `lidc-val-rotated` and both probes
  [code: preflight.py:28, 66], and `refusal` recomputes them at launch [code: preflight.py:84]. A moved view gives
  "cannot check ... against the current code, data and config" and the launch is refused.
- So the literal rule stops V1 at $0, and the engineer must guess which views may stay.

Exact fix. Rule 7: "... the LIDC eval views named in the rental's table, plus the views the bundle and the stamp
always hold (`bundle.LIDC_EVAL` and `preflight.EXAMS`: `lidc-val-rescaled`, `lidc-val-rotated`, `probe-precision`,
`probe-topology`; LIDC val only). Only `lidc-test-*` is moved out (6.1 step 4c)." Add to 6.1 step 4c: "the bundle
list equals `LIDC_EVAL` + `lidc-val-edge`".

### 7. [major] [B] G05 has no runner that can see the engine it checks

Problem.
- B "### B.9.5", runner table: "G05 | the first run that uses `fp16-seqmask` or `lora-seqmask` (after Gate B3),
  inside its first 15 steps | `python -I -m pytest -m gpu tests/test_train_configs.py -k G05`".
- G05 checks the live engine: "after the first sync and again after one sleep/wake cycle, ... the vLLM slice of the
  packed `lora_a` equals `A.to(bf16)`". A separate pytest process has no access to the engine inside a running
  `train_rl`, and a second engine on the same 32 GB card while training runs is not planned.
- B "#### B2.4" rule 0 puts the same check inside the run: "The F check is the first 15 steps of that run itself, in
  one process on one host". Gate B2 lists G05 (B "### B.9.10": "Gate B2 (run) | ... G02, G05"). So the gate's test
  has no defined place to run.

Exact fix. B.9.5: G05 becomes two parts. G05a (pytest, box, before training, on the tiny model and the pinned vLLM;
the dtype and packed-layout assertions on a fresh engine). G05b = the in-process F check of B2.4 rule 0, written by
`train_rl` to `<out>/f_check.json` at step 15 (NaN count, per-module `torch.equal` result, ABAB ratio); a Mac
`data` test reads the pulled file. Gate B2 cites G05a and G05b. Put the run in V3 after Gate B3 (item 1 of this review).

### 8. [major] [B, 00_front] V3's B2.3b box part needs the B1 LoRA adapter, but only the merged model reaches HF

Problem.
- B "#### B1.5", "B1 adapters to V3": "B1.5 adds `hf_push --rl <train_rl out>`: it uploads `final/` (the merged
  model, the form `fetch_hf.sh` and `final_eval --ckpt` read) to revision `b1-rows`." `train_rl` merges the
  adapter into the base before it saves `final/` [code: experiments/phase2/train_rl.py:10, 59].
- B "#### B2.3b", box part: "Generate ... from the B1 final adapter ... (a)-(d) again on the same 256 answers".
  Way (b) is the plain trainer forward with LoRA, way (c) "the trainer forward with the LoRA delta merged". With a
  merged model there is no LoRA delta, so (c) - (b) is 0 by construction, and the merge-rounding reading ("(c) >> (b)
  on a drifted adapter") cannot appear. G04 ("on the B1 final adapter") then passes vacuously.
- 00_front "### 6.2", item 2 also says "final adapter pushed to HF (revision `b1-rows`)".

Exact fix. B1.5: `hf_push --rl` uploads `final/` and the last adapter checkpoint (`trainer/checkpoint-<last>/adapter_*`,
about 17.4 M parameters [measured: B2.2, checkpoint-100 shapes]) to `b1-rows`. B2.3b box part: load the base
`rl-full-3ep` plus that adapter. Front 6.2 item 2 and G04: "B1 final adapter (`trainer/checkpoint-<last>`) and
merged model". Add a unit test that `hf_push --rl` lists both folders.

---

## Minor

### 9. [minor] [A, B, D] `make test-data` has two definitions

A "## A.9. Tests": "new `test-data` ... runs `$(PY) -I -m pytest -q -m "data and not slow"`" and a separate
`test-all-data` with `-m data`. B "### B.9.0" and D "### D8" say `make test-data` is `$(PY) -I -m pytest -q -m data`.
Fix: B.9.0 and D8 copy A's two targets ("`test-data` = `data and not slow`; `test-all-data` = `data`") and say which one
each `data` test of B and D needs (B F01, F02, F06 are performance tests: name the target).

### 10. [minor] [C] Stale text about A.9 and V1

C "### C.9.0 Conventions": "A.9 does not have this hook yet (checked 2026-11-08)". A.9 now defines it. C "## C.0.
Skills used", system-design row: "V1 at 27.5 s skips A3-val". Front 6.2 now says item 5 is "Not skippable". Fix:
delete both clauses.

### 11. [minor] [B] B.9.5 runner table names a short D8 file list

B "### B.9.5": "G01, G06 | R-D2 ... | `make gpu-check` (D8 file list: `tests/test_gpu_scoring.py`,
`tests/test_cuda_checker_gpu.py`)". G01 lives in `tests/test_scorer_sums.py` (B.9.8 sketch header), which D8 adds
"once B3.1 lands". Fix: "D8 file list (`D_GPU_TESTS`), which includes `tests/test_scorer_sums.py` after B3.1".

### 12. [minor] [E, 00_front] GP1 and GP2 are V4 tests, but V4 does not list them

E "### E.9.5": GP1 and GP2 are "(`gpu`; ...; V4)". 00_front "### 6.3", V4: "E-G10 box part ..., E4 deterministic
check and E-G13. Work: about 9 min". Fix: add GP1/GP2 to the V4 list with their minutes, or move them to R-D2's
`D_GPU_TESTS` after E3.

### 13. [minor] [E] X2 needs saved test answers that do not exist until V3 item 3

E "### E.9.6", X2: "Re-score the saved `rl-full-3ep` test completions". `runs/r5/eval/eval_test/` has scores only
[measured: ls; first line of `items.jsonl`]. The answers exist only after V3 item 3 (`a3_test.yaml`). Fix: X2 reads
`runs/a3/eval_test/answers.jsonl`, skips with that file name as the reason, and the row says "scores only; it reads
no new test label".

### 14. [minor] [B, 00_front] "Final adapter" names two things

B1.5 calls `final/` "the merged model"; front 6.2 item 2, B2.3b and G04 call it "the B1 final adapter". Section 0
rule: one word per thing. Fix: "merged model" for `final/`, "adapter" for `trainer/checkpoint-*` (see item 8).

### 15. [minor] [all phases] Long sentences remain

1,314 of 7,204 prose sentences have more than 25 words (18.2 %), with DOC.6's rule [measured: `(planning probe, not kept: qa4probe)/count.py`;
my splitter differs from round 3's, so compare shares, not counts]. Fix: as DOC.6 says, split gate cells and task steps
first; new text keeps one rule per sentence.
