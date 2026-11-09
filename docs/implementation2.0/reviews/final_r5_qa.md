# Final review round 5: QA lead and technical writer (ASD-STE100)

Target: `docs/implementation2.0/Implementation.md` (round 5; 1,501,724 bytes, 15,116 lines [measured: `ls -l`, `wc -l`]),
built from `docs/implementation2.0/*.md` by `assemble.py`. I read the Tests sections (A.9, B.9, C.9, E.9; D.9, DOC.9 and
CLEAN.9 by grep), the gate tables (A.6, B.6, C.6, D.6, DOC.6), front sections 1, 1.1 and 6, 99_back decisions 1-3,
and the code each claim names. I checked every item against `review-log.md` (rounds 1-4 and the round 4 prep pass);
fixed items are not repeated. No REPO file was edited. No git write. No rental. No Python read `data/`.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No new NVIDIA product claim, so no catalog lookup. It pointed me at the box toolchain again (`run_box.sh` installs nvcc; the image has none), which is the root of item 3 (V4 has no build). |
| engineering:code-review | Read the code behind each launch and test claim: `train_rl.py:23-36` (no `--numerics`), `policy.py:399, 424-494, 527-536` (`sync_up`, `chain`, `cmd_run`, `cmd_fetch_hf`), `CMakeLists.txt:15`, `gpucheck_py.cpp:113`, `bundle.py:26, 50-64`, `preflight.py:28`, `final_eval.py:48-56`, `hf_pull.py`, `run_box.sh`, `runs/k1/checker_gpu.log`. Items 1, 3, 4, 6. |
| engineering:debug | Traced each V1 plan step on the `v1-b1` tree: which flags exist on that tree, which test parses them (I47), and what exit code stops the `&&` chain. Item 1. |
| engineering:system-design | Put the after-V1 decision table (front 1.1, 99_back 2b) beside the Gate B1 rules of B.6.2 and A-G3; found the view-identity contradiction (item 2). |
| engineering:architecture | Checked the order rule 2 tag contents against the phase that owns each flag (B2.2 owns `--numerics`). Item 1. |
| engineering:documentation | One rule per place: `make test-data` has two definitions (item 9); stale text about the A.9 hook (items 7, 8). |
| huggingface-skills:huggingface-best | No model choice changed. The HF part is where `fetch_hf.sh` puts each revision on the V3 box (item 4). |
| engineering:testing-strategy | Matched each box test to its rental, its build step and its selection form (node id or `-k`); matched each money guard to a test that cannot skip. Items 3, 5, 6, 11. |
| mattpocock-skills:tdd | Looked for tests that pass by construction: the `BudgetStop` sketch skips whenever the signature of B1.5 is followed (item 5); two "CPU build equals CUDA build" tests have no second build to compare (item 6). |
| Web | pytest exit code 5 = "No tests were collected" [web: https://docs.pytest.org/en/stable/reference/exit-codes.html]. pytest `-k` matching is case-insensitive since pytest 5.4/6.0 (issue #6316) [web: github.com/pytest-dev/pytest commit 81da5dac], so item 11 is minor, not major. |

Measured or read during this review (2026-11-08):
- `grep add_argument experiments/phase2/train_rl.py`: flags `--model --ckpt --train --train-views --out --steps
  --save-steps --shared-prefix --lr --resume --procs --grpo`. No `--numerics`.
- `.venv/bin/python -I -c "import experiments, tests"` from `/tmp` works (the `_objdel_repo.pth` adds the repo root),
  so `-I -m pytest` on the box can import `experiments` and `tests.oracles`.
- `runs/k1/checker_gpu.log`: the CUDA build step took about 20 s (wheels 3.1 s, build 15.7 s), so the 3-min build
  estimate in front 6.2 is safe.
- Tables with a `$` header: only lines 688, 755, 777 (all in front section 6), so DOC check 16 stays green.

Severity follows the task rules. Count: **3 blocking, 3 major, 8 minor**.

---

## Blocking

### 1. [blocking] [B, 00_front] V1 step 3 passes `--numerics`, a B2.2 flag that order rule 2 keeps off the `v1-b1` tag

Problem. Two parts both claim to be the rule.
- B "#### B1.5 Train on the ordered rows", step 3: "`train_rl --model qwen3-vl-2b --ckpt {rd}/ckpt` ...
  `--out {rd}/b1-rows --numerics bf16-seqmask --lr 5e-5 --steps 260 ...`".
- B "### B.9.3 Integration tests", I22: "step 3 `train_rl` with ... `bf16-seqmask` ..."; and in the same test:
  "`b2_fp16.yaml` (`fp16-seqmask`, ...) and `b3_v5.yaml` (... `--numerics bf16-seqmask`, `--reward v5`, ...) parse
  and name existing scripts and flags".
- The flag is B2.2 work: B "#### B2.2 Numerics presets": "- `train_rl --numerics {bf16-seqmask,bf16-tis,
  bf16-tis-opsm,fp16-seqmask,lora-seqmask}`; `rl.json` records it." B1.5 adds only `--rows-revision` and
  `--rows-manifest-sha256` ("New flags ...").
- 00_front "## 1", order rule 2: the tag "holds no C, D, B2.2-B2.4, B3 or kernel code", and "B2.1 as logging only".
- Code today: `train_rl` has no `--numerics` [code: experiments/phase2/train_rl.py:23-36, measured by grep above]
  and no `--reward`.
- B "### B.9.16", I47: "for every `-m experiments.phase2.<module>` step of `b1_rows.yaml` ... build that module's
  argparse parser ... and parse the step's argv: no `SystemExit`".

Result. On the tag, I47 raises `SystemExit(2)` on `--numerics`, so `make test` is red and the owner cannot tag
(order rule 2 and front 6.1 step 1). I22 also fails on the tag: `--reward v5` and `fp16-seqmask` do not exist before
B2.2 and B3.5. If an engineer skips I47, `run.sh` adds only `-I` [code: policy.py:486], and argparse exits 2 at step 3
on the paid box after item 1 (53 min, $0.70, front 6.2), and the `&&` chain [code: policy.py:490] also skips step 4.
If the engineer adds B2.2 to the tag instead, order rule 2 is broken.

Exact fix (keep order rule 2; today's defaults already are `bf16-seqmask`, `configs.py:188` and the TRL defaults
that B2.2 `_SEQMASK` writes out):
- B1.5 step 3: delete `--numerics bf16-seqmask`. Add one sentence: "No `--numerics` on the tag: `train_rl` without the
  flag runs today's `sequence_mask`, cap 3 (B2.2 `_SEQMASK` writes the same values out). `rl.json` records
  `numerics: "bf16-seqmask"` (B1.5 writes the literal until B2.2 adds the flag)."
- B.9.3 I22: split it. I22a (on the tag, `make test`): `b1_launch.yaml` and `b1_rows.yaml` checks as written, with
  "no `--numerics` key in step 3" instead of "`bf16-seqmask`", plus `--shared-prefix` (item 12). I22b (lands with
  B2.3 and B3.9): the `b2_fp16.yaml` and `b3_v5.yaml` checks. Add to B.9.0: "a test that names a config file is
  written in the commit that creates that file."
- I47: add the case "every flag of a `b1_rows.yaml` step exists in the parser of the `v1-b1` tree" (it already does
  this by parsing; state that it runs in the order rule 2 pre-tag list).
- 00_front order rule 2 pre-tag list: add "I22a and I47".

### 2. [blocking] [00_front, 99_back, B, A] View identity: the decision table stops, the Gate B1 rules continue

Problem. Two parts both claim to be the rule, and the table says the result cannot bend it.
- 00_front "### 1.1 After V1: decision table": "The rules named are Gate B1's hard rules (B.6.2)." Row: "Fails the
  offset < 0.15 guard, the edge guard, legal or view identity | Stop. Owner review | B1.8 waits | R5".
- 99_back "## Owner decisions needed", 2b: "The offset < 0.15 guard, the edge guard, legal or view identity fails:
  stop, owner review; B1.8 waits; R5."
- B "### B.6.2 Run gates", "View identity (from A-G3)": "If this correlation rule fails, Gate B1 still pairs with
  A3's items (same frozen view, same box); only the R5 report line is dropped. If A-G3's TOST fails, Gate B1 is still
  decided on A3's items". View identity is not in the Gate B1 hard-rule list of B.6.2 (val IoU, far misses, near IoU,
  difference in differences, edge, legal).
- A "## A.6. Acceptance gates", A-G3: "A val failure does not block Gate B1: B.6.2 then pairs with A3's items and
  drops only the R5 second line".
- B "### B.9.1", U49 holds both: "B1 passes only when all hard rules hold (...; when item 5 ran, view-identity
  correlation >= the `r_sh` lower bound)" and "item 5 ran but the correlation rule fails: the gate still pairs with
  A3's items and drops the R5 line".

Result. With all hard rules passed and the correlation rule failed, the table says "stop, R5 is final"; B.6.2, A-G3
and the second half of U49 say "pass, B-S2 next". The next paid step (B-S2, $2.11 front 6.3) depends on which part
the engineer reads, and U49 cannot be written (its two clauses disagree).

Exact fix (B.6.2 is the reasoned rule: the primary pair is A3 vs B1 on the same frozen view and box, so the R5
correlation only guards the R5 report line):
- 00_front 1.1 row 4 and 99_back 2b: replace "or view identity" with "or the frozen-view hash (DG15:
  `view_content_sha256` of the B1 and A3 eval folders equal `frozen_views.json`)". Add under the table: "The
  view-identity correlation of A-G3 is not a Gate B1 rule: when it fails, Gate B1 pairs with A3's items and drops the
  R5 report line (B.6.2)."
- B U49: delete "when item 5 ran, view-identity correlation >= the `r_sh` lower bound" from the hard-rule list; add
  "a planted `view_content_sha256` mismatch fails Gate B1 with 'stop'".

### 3. [blocking] [E, 00_front] V4 has no CUDA build step: GP1/GP2 fail on the CPU build and the E4 sampler arm scores on the host

Problem. The round 4 fix (review-log "cuda-2, qa-2": build step 0 in `b1_rows`, `a3_measure`, `c_speed`, `a3_test`,
`b3_v5`) did not reach V4. The r4 QA item 2 asked for "Same first step in ... the V4 plan (E.8)"; the E part of the
round 4 log has no qa-2 row.
- 00_front "### 6.3 V2 ...", V4 paragraph: "E-G10 box part (prompt ids and SamplingParams equal; TOST sampler arm on
  val), E4 deterministic check and E-G13. Work: about 9 min ... 30 min overhead + 9 min". No build.
- E "## E.8. Cost and time; dependencies", V4: "It runs E-G13 ..., the E4 deterministic check ... and the E4 sampler
  arm. ... the box runs one eval only: 1,084 items x K 8". No build.
- E "## E.9. Tests": "All five run at V4 item 'E checks' (E.8) with one file list: `pytest -m gpu
  tests/verify/test_verifier_gpu.py tests/infra/test_box_smoke.py tests/policy/test_sampler_box.py`", and "One skip
  rule ... it skips unless `objdel._gpucheck.build() == "cuda"`". GP1: "Skips through the root conftest hook (A.9)".
- A "## A.9. Tests": "On a box the hook fails instead of skipping: when `OBJDEL_RUN_DIR` is set ... and the build is
  `"cpu"`, every `gpu` test fails".
- Code: `sync_up` installs with `uv pip install ... -e .` [code: policy.py:399]; CMake builds the CPU extension unless
  `OBJDEL_CUDA=ON` [code: CMakeLists.txt:15]; a `steps:` plan runs only its listed steps [code: policy.py:432-433].

Result. On V4 ($0.54, front 6.3) GP1-GP5 fail by the A.9 hook, so the E-G3/E-G11 box halves cannot pass, and the E-G10
sampler arm (1,084 x 8 answers) is scored by the CPU build, which breaks the owner rule "CUDA is a must for the
checker, the scorer" (BRIEF section 3). The rental is spent for nothing. E.9's "skips" text also contradicts A.9's
"fails on a box".

Exact fix.
- E.8 V4: name the plan file `experiments/phase2/configs/e_v4.yaml` (E creates it; front 6.3 names it). Step 0 =
  `b1_rows.yaml` step 0 (CUDA build, `build() == 'cuda'` assert, `tests/test_gpu_scoring.py`); step 1
  `env_check.sh` (E-G13); step 2 the E.9 pytest list by file; step 3 the sampler arm (`final_eval` sampler path, CUDA
  build, pool size 1 as A3).
- 00_front 6.3 V4: "plan `e_v4.yaml`, step 0 the CUDA build (3 min [estimate: front 6.2]); 30 + 3 + 9 = 42 min =
  $0.58 at r_s" (recompute the V4 line and the "(a)" total "V4 adds").
- E.9 "One skip rule": add "On a box (`OBJDEL_RUN_DIR` set) the hook fails a `gpu` test on a CPU build (A.9)". GP1:
  "Skips on the Mac; fails on a box with a CPU build (A.9)."
- Test: extend B I22a's step-0 check (or a new E test) to "`e_v4.yaml` step 0 equals `b1_rows.yaml` step 0".

---

## Major

### 4. [major] [00_front, B] V3 item 0 relies on an output-folder option of `fetch_hf.sh` that B1.5 never defines

Problem.
- 00_front "### 6.3", V3 item 0: "`fetch_hf.sh <run> b1-rows` and `fetch_hf.sh <run> rl-full-3ep` (each into its own
  folder: B1.5 gives `fetch_hf.sh` an output folder name, because today every revision lands in `{rd}/ckpt` [code:
  policy.py:527-536])".
- B "#### B1.5", "B1 adapters to V3": "V3 starts with `fetch_hf.sh <run> b1-rows` when B2.3b's box part runs." B1.5
  defines only `--rows REV ROOT` ("Rows-root download"). B "#### B2.3b", box part: "The B1 final adapter reaches the
  box by `fetch_hf.sh <run> b1-rows` (B1.5). Load `final_adapter/` on `rl-full-3ep`".
- B "### B.9.16", I40 tests only `--rows`.
- Code: `cmd_fetch_hf` always runs `hf_pull ... --out /workspace/runs/{a.run_id}/ckpt` [code: policy.py:534-535], and
  `hf_pull` calls `snapshot_download(local_dir=a.out)` [code: experiments/phase2/hf_pull.py:22].

Result. The engineer must invent the option name, its box path and the paths that `b3_v5.yaml` (B3.9: "the same start
checkpoint as B1 (`rl-full-3ep`)") and B2.3b read. Without it both revisions land in `{rd}/ckpt`: the B1 `final/`
and `final_adapter/` folders mix with the R5 model files, and any file both revisions hold (model card, LICENSE that
`hf_push` writes) is overwritten by the second fetch.

Exact fix.
- B1.5: "`fetch_hf.sh <run> REV [--as NAME]` downloads to `{rd}/<NAME>` (default `ckpt`, today's path). V3 item 0 runs
  `fetch_hf.sh <run> rl-full-3ep` (to `{rd}/ckpt`) and `fetch_hf.sh <run> b1-rows --as b1` (to `{rd}/b1`). B2.3b reads
  `{rd}/b1/final_adapter`; `b3_v5.yaml` `--ckpt {rd}/ckpt`."
- 00_front 6.3 item 0: replace "B1.5 gives `fetch_hf.sh` an output folder name" with the two exact lines above.
- B.9: new I40b: `cmd_fetch_hf` with `--as b1` passes `--out /workspace/runs/<id>/b1`; without `--as` it passes
  `.../ckpt`; `NAME` with `/` or `..` raises.

### 5. [major] [B] The `BudgetStop` boundary test skips whenever B1.5's signature is followed

Problem.
- B "#### B1.5": "Run guard `callbacks.BudgetStop(rate_usd_per_h, cap_usd, session_start, fixed_usd=0.0)` ... Unit
  test with a fake clock." The signature has no clock parameter.
- B "### B.9.8 Pytest sketches", `tests/test_train_rl_guards.py`:
  `if "clock" not in BudgetStop.__init__.__code__.co_varnames: pytest.skip("BudgetStop has no injectable clock; patch time.time instead")`.
- B "### B.9.1", U46 needs exact trip times with a fake clock (0.4999 h vs 0.5 h; 6.8182 h; 2.8713 h with
  `fixed_usd`), and the sketch has no `fixed_usd` case.

Result. An engineer who writes `BudgetStop` exactly as B1.5 says gets a skipped test, and `make test` is green. The
guard that stops the paid V1 run at $2.05 (front 6.4 rule 5) is then untested. U46 must guess how the clock enters.

Exact fix.
- B1.5: signature `BudgetStop(rate_usd_per_h, cap_usd, session_start, fixed_usd=0.0, clock=time.time)`.
- B.9.8 sketch: delete the `pytest.skip` branch; add the U46 case `fixed_usd=0.0975, rate 0.68, cap 2.05` trips at
  2.8713 h and not at 2.8712 h.
- U46: add "`clock` is a keyword argument (the test fails, never skips, without it)".

### 6. [major] [A, E] Two box tests compare the CPU build with the CUDA build in one process; only one build exists there

Problem.
- A "### A.9.6", G02b: "Marker `gpu`. Run place: box, V3 item 1 (R-D2) ... `ceiling.item` on 50 val items:
  `doctor`, `model` equal with `==` under both builds. 500 random items ...: the four terms equal per sample."
- E "### E.9" table, GP1: "scored with `checker_build="cuda"` and `"cpu"`: reward, terms, defects equal by hex".
  `checker_build` is a read-only stamp: E "E7.1": "`checker_build` = `engine.build()` ("cuda" or "cpu")".
- Code: one extension `objdel._gpucheck` per venv; `build()` returns `"cuda"` or `"cpu"` for the one compiled build
  [code: cuda/gpucheck_py.cpp:113]; the option is a CMake switch at install time [code: CMakeLists.txt:15]. On the
  box, step 0 replaces the CPU build in `/opt/venv` (front 6.2 "Why step 0").
- D solved the same problem with Mac goldens: D "## D.8" (round 4 log): "the tracked digest golden
  `tests/goldens/d_gpu02.json`, written on the Mac by `experiments/phase2/perf/make_d_goldens.py`". A and E name no
  golden and no second venv.

Result. Neither test can be written as specified. The engineer must guess the reference (Mac golden, numpy host
oracle, or a second CPU venv), and each guess checks a different thing. G02b feeds A-G2's coverage row; GP1 is the
box half of E-G3.

Exact fix.
- A G02b: "On the Mac (CPU build) `experiments/phase2/perf/make_a_goldens.py` writes `tests/goldens/a_g02b.json`
  (per item: `doctor`, `model`; per random item: the four terms as float hex). The box test loads the golden and
  compares with `==`; a missing golden fails (not skips). The golden and the val item ids ship in the code sync
  (`tests/` ships [code: policy.py:376-383])."
- E GP1: the same pattern, `tests/goldens/e_gp1.json` from a Mac script, or reuse `d_gpu02.json` if E3a keeps the
  engine API. Delete "scored with `checker_build="cuda"` and `"cpu"`".

---

## Minor

### 7. [minor] [C] C.9.0 says A.9 has no hook

C "### C.9.0 Conventions": "A.9 does not have this hook yet (checked 2026-11-08); C writes no skip line of its own and
relies on it." A.9 now has the hook (A "## A.9. Tests", "One GPU skip rule for all phases"). Fix: delete "A.9 does not
have this hook yet (checked 2026-11-08);".

### 8. [minor] [B] B.9.0 still defines G tests by the old skip

B "### B.9.0": "**G** GPU (rented host or local CUDA; skipped without CUDA, as `tests/test_cuda_checker_gpu.py`)". That
file skips on `nvidia-smi -L` [code: tests/test_cuda_checker_gpu.py:27], and D8 deletes that line. Fix: "**G** GPU
(marker `gpu`; skipped on the Mac, failed on a box with a CPU build: the A.9 root hook)".

### 9. [minor] [A, B, D, CLEAN] `make test-data` has two definitions

A "## A.9. Tests": "new `test-data: $(EXT_STAMP)` runs `$(PY) -I -m pytest -q -m "data and not slow"`" and
"`test-all-data` ... `-m data`". B "### B.9.0", D "### D8" and CLEAN (CLEAN0 text) say "`make test-data`, which is
`$(PY) -I -m pytest -q -m data`". Fix: in B.9.0, D8 and CLEAN write "`make test-data` (A.9: `-m "data and not slow"`;
`make test-all-data` adds the slow `data` tests)".

### 10. [minor] [C] F leftovers in C after rule (b)

C "## C.6", G-SP1 row and C "### C.9.5" E-SP1 keep F margins ("Under F: `- 2.7 s` (A) or `- 2.1 s` (B)"), while C.6
says "Under F (runs after Gate B3) no per-parameter TRL path exists, so there is no F baseline arm" and "G-SP1 is not
re-gated". C "### C.9.4" table: "G-SY1, G-SY2, G-SY3 | V2 item 2 | `c_speed_check --sync` (G-SY1 runs C5 step 5(c)
under F)", while front 6.3 says "V2 has no F arm" and C.6 says G-EQ5 re-runs in the first run that uses F. Fix: delete
the F margins from G-SP1 and E-SP1 (keep them in C.8.1 as a report line); C.9.4: "G-SY1 merged part in V2 item 2; the
F part (C5 step 5(c)) runs in the first run after Gate B3 that uses `lora-seqmask`".

### 11. [minor] [B, C] `-k` box selections against the B.9.0 rule

B "### B.9.0": "A box command names node ids, never `-k`". B.9.5 G05: "`python -I -m pytest -m gpu
tests/test_train_configs.py -k G05`"; C.9.4: "G-GR1 ... pytest `-k gr1`", "G-RC0 ... pytest `-k rc0`". pytest `-k` is
case-insensitive [web, see Skills used], so `test_g05_...` matches, but an empty selection exits 5 [web: pytest exit
codes] and stops a `&&` chain. Fix: name node ids (`tests/test_train_configs.py::test_g05_vllm_dtype`,
`tests/test_c_gpu.py::test_g_gr1_...`, `::test_g_rc0_...`) and add them to the collect-only check of front 6.1 step 4d
for their plans.

### 12. [minor] [B] I22 does not check `--shared-prefix` in step 3

00_front "### 6.2": "It runs the owner tag `v1-b1` ... with `--shared-prefix`". B1.5 step 3 passes it. I22's step-3
list has no `--shared-prefix`. Without the flag R5 ran 31.0 s/step against 23.6 s (README, commit be52b1c), so 260
steps would not fit before `BudgetStop` (front 6.2 "Why 260 steps": the stop trips above about 27.5 s). Fix: add
"`--shared-prefix`" to I22(a)'s step-3 list.

### 13. [minor] [DOC] DOC-G13 baseline size is stale

DOC "## DOC.6", DOC-G13: "`Implementation.md` 1,480,096 bytes [measured: `ls -l`, 2026-11-08]". Today 1,501,724 bytes
[measured: `ls -l`]. Fix: update the number (still under 2 MB).

### 14. [minor] [00_front] The eval-tail margin does not count `train_rl`'s final save

00_front "### 6.2", Notes: "`BudgetStop` at $2.05 ends training at session minute ... 172.3 at the latest, and 172.3
+ 26 = 198.3 <= 201." After the stop, `train_rl` saves the checkpoint, saves `final_adapter/`, merges and writes the
4.0 GB `final/` [measured: `du -sh runs/rl/train_rl/final` = 4.0G; code: train_rl.py:58-59]. That time is not in the
26-min tail, and the margin is 2.7 min. No measurement shows it is over. Fix: add "plus the final save (about 1 min
[estimate]; the full run's `final/` is 4.0 GB)" to the tail and re-check `198.3 + 1 <= 201`; or log the save time in
`rl.json` (`timing_s.save`) on V1 for the next plan.
