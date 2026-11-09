# Final review round 15: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc reviewed: the built
`docs/implementation2.0/Implementation.md` (17,105 lines, 1,760,475 bytes, built 2026-11-09 03:28, the same minute as the
newest source `phases/B.md` [measured: `wc -l`, `ls -l`, `stat`]). Findings are written against the sources in
`impl/phases/`. I read `review-log.md` (round 14 sections of 00_front, B, D, E, CLEAN, 99_back; round 13 for D and
CLEAN), `final_r14_arch.md` (no verdict file) and `final_r14_fresh-arch.md` with its verdict. Focus this round: the
round-14 changes as they meet E, CLEAN and DOC (fresh-arch-1: old samplers archived in the E4 step 3 commit; rl-2:
`train_rl --merge-last` step 3b in B3.9; rl-3: `train_rl --opsm-start-delta`; qa-1: `tests/goldens/d_r01_lidc.json`).
No item repeats an item that a log marks as handled. The six round-14 arch minors (strict docs build with 3 `] [`
autorefs, `mkdocs.yml (DOC10)`, stale DOC1 size, `test_hooks` / `tests/variants.py` exception, `variants.py` test
paths, README status) are still open in the sources [measured: A.md line 84, 00_front.md line 358, DOC.md line 325,
no `test_hooks` in E.md] and are not filed again. I edited no repo file and ran no git command that changes anything
(only `git ls-files`, `git check-ignore`).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item depends on an NVIDIA product claim (this round's changes are sampler archiving, CLI flags and test plumbing). No catalog lookup, no install. |
| engineering:code-review | Re-grepped every user of `sample_vllm` / `hf_sample` in `tests/`, `src/`, `experiments/` and checked E4 step 1's list against it (items 2, 3). Read `select_ckpt.py:55,88-96` (second-venv helper) and `train_rl.py:23-36` against `RlStage` (item 4). |
| engineering:debug | Walked the V4-failure branch step by step: restore the two files from the tag, then what calls them, what `make test` says (U2 "no V4 exception"), and which eval path is live (item 1). |
| engineering:system-design | Order walk E4 step 1 -> E4 step 3 -> E6.0 spike -> E6.4 -> E8 gate for each new or moved CLI and flag (items 2, 4, 5). |
| engineering:architecture | Rollback design: a rollback must restore a working path, not only files. Compared E4 step 1's rollback with E.7's rollback row (item 1). |
| engineering:documentation | Found the stale E.7 mitigation "samplers archived only after V4" left by the round-14 fix (item 1). |
| huggingface-skills:huggingface-best | No model choice in scope. Nothing changed. |
| system-design-skills:system-design | Same order walk, framed as "which artifact does each step break, and which planned check sees it before spend". |
| system-design-skills:service-decomposition | One package with an exhaustive import-linter `layers` contract (a modular monolith) is still the right form; no split proposed. |
| mattpocock-skills:codebase-design | Interface check of `RlStage` as the seam of `train_rl`: its required fields (`train`, `steps`) do not fit the `--merge-last` mode that B added, so the seam is too narrow for one of its two callers' modes (item 4). Deletion test on `vllm_python` / `VLLM_VENV` after the sampler switch (item 3). |
| WebSearch / WebFetch | Loaded. Not needed for a finding: every item rests on repo text and code. |

## Checked and found correct (no item)

- E4 step 1 code citations: `select_ckpt.py:93-96` is `sampler_cmd` and the only live caller of both old samplers;
  `parity.py:5,71,89` and `r0_parity.py:7` name `sample_vllm.py --logprobs`; `r0_parity.py:44` passes
  `extra=("--logprobs",)`; `sample_vllm.py:47-49,62-63,97` hold the logprobs request, record and engine mode;
  `final_eval.py:51` is `--train`, `:73-75` is `make_trainer` + `sample.sample`; `g_probe.py:44-48` is the K = 8 draw
  [code, read 2026-11-09].
- `pyproject.toml:67` `testpaths = ["tests"]`, so tests moved to `archive/phase2_samplers/tests/` are not collected.
- No live doc outside the frozen sets names `sample_vllm` or `hf_sample` (only `docs/research/reports/Code performance
  audit 2026-11-06.md`, a frozen report), so DOC4 check 2 stays green on the archive commit.
- `experiments/phase2/vast/code_files.txt` lists both samplers, but it is generated and git-ignored
  (`.gitignore:30`), so the archive commit needs no edit there.
- `tests/goldens/d_r01_lidc.json` (round 14 qa-1) has no `riga` in its name, so CLEAN's CT-89 name rule allows it.

## Items

### 1. major: the V4-failure rollback restores two files that nothing calls, turns U2 red, and leaves `final_eval` on the failed sampler; E.7 gives a different rollback and a stale mitigation.

Sources: E (E4 step 1, E4 step 3, E-G14, U2, E.7), 00_front (6.3 V4 paragraph).

- E.md "### E4" step 1: "Rollback is the tag `pre-e4-eval` (E4 step 3). Owner command: `git checkout pre-e4-eval --
  experiments/phase2/sample_vllm.py experiments/phase2/hf_sample.py`. V4 still closes E-G10's box half. If V4 fails,
  the owner restores them from the tag." 00_front.md 6.3: "If V4 fails, the owner restores the old samplers from that
  tag (E4 step 1)."
- E.md "### E4" step 1 also: "`sft.select.sampler_cmd` now builds `[sys.executable, "-I", "-m",
  "objdel.policy.sampler", "--kind", kind, ...]`", and step 3: "`final_eval` calls `VllmSampler` + `evaluate`. It no
  longer calls `build_trainer` or `trainer._generate`." The two restored scripts had one caller,
  `sampler_cmd` [code: experiments/phase2/select_ckpt.py:94-95, `head = [vllm_python(), "-I", str(ROOT /
  "experiments/phase2/sample_vllm.py")]` / `"-m", "experiments.phase2.hf_sample"`]. After the restore nothing calls
  them, and `final_eval` still samples through `VllmSampler`, the path V4 just failed.
- E.md E-G14: "`sample_vllm.py` and `hf_sample.py` exist only under `archive/phase2_samplers/` from the E4 step 3
  commit on; no exception". U2: "the two samplers from the E4 step 3 commit on, with no V4 exception". So the restore
  makes `make test` fail (R-3), and the engineer has to undo it.
- E.md "## E.7" row "VllmSampler differs from TRL's generation path": Mitigation "...TOST ±0.01 on V4; samplers
  archived only after V4" (stale since round 14); Rollback "`final_eval` keeps the trainer path in git history (tag
  `pre-e4-eval`, E4 step 3); revert E4 only". "Revert E4 only" is not mechanical after E5/E6/E8, which delete
  `make_trainer`, `fast.py` and `--train` (E-G14; E6.1 EvalStage "no `train`").
- Effect: on the post-E raise path, the documented V4-failure action does nothing to the eval path, so later evals
  (B3's both arms, any test re-eval) keep the sampler whose prompt ids or SamplingParams V4 found different. No $0
  check sits before V4 on this branch (V4 is the check). The engineer must guess the real action.

Fix (E.md E4 step 1, E.7 row; 00_front 6.3 V4 paragraph):
1. Replace "If V4 fails, the owner restores them from the tag." with: "If V4 fails, the old scripts are NOT restored
   (they have no caller after step 1, and U2 forbids them outside `archive/`). The deterministic part failing is a
   bug in `policy.sampler`: fix it in a change commit, add the failing row to `tests/policy/test_sampler.py`, and
   re-run the V4 deterministic test on the next owner-approved rental. The TOST failing alone is the owner's call
   (E4 gate, 'not shown'). Until V4 passes, every post-E eval result carries the label 'sampler not box-verified
   (E-G10 open)'. The tag `pre-e4-eval` stays as the reference copy for reading only."
2. E.7 row: Mitigation "...TOST ±0.01 on V4; samplers archived in the E4 step 3 commit (final round 14)". Rollback
   "Fix forward in `policy.sampler` (E4 step 1, V4-failure rule); the trainer path stays readable at tag
   `pre-e4-eval`".
3. 00_front 6.3: "If V4 fails, E4 step 1's V4-failure rule applies (fix forward; no restore)."

### 2. minor: E4 step 1's test list misses one test, and the switch of `sampler_cmd` (step 1) and the test fixes (step 3 commit) are in different commits, so the step-1 commit fails R-3.

Sources: E (E4 step 1).

- E.md E4 step 1 lists the tests to fix: "`tests/test_p2_bridge.py`, `test_p2_breakit_scoring.py` and
  `test_p2_v1_select.py` (import `sample_vllm`), `test_p2_hf_eval.py` and `test_p2_breakit_runners.py` (import
  `hf_sample`), and `test_p2_scripts.py:115,280` (path asserts)". Missed: `tests/test_p2_scripts.py:35`
  `@pytest.mark.parametrize("script", [..., "hf_sample"])` in `test_help`, which runs `python -I -m
  experiments.phase2.hf_sample --help` [code]. `tests/test_profile.py:113` names `hf_sample` too, but as a negative
  assert, so it still passes.
- E.md E4 step 1: "This step switches it to `objdel.policy.sampler` (below)" and "The same commit fixes every test
  that imports or names them" where "the same commit" is "the E4 step 3 commit". `tests/test_p2_scripts.py:115`
  (`cmd[2].endswith("sample_vllm.py")`) and `tests/test_p2_breakit_runners.py:359` (`cmd[2:4] == ["-m",
  "experiments.phase2.hf_sample"]`) fail on the commit that switches `sampler_cmd` if that is a step-1 commit.
- Caught-before-spend: R-3's `make test` fails at $0 on that commit and names the test. So minor.

Fix (E4 step 1): add "`test_p2_scripts.py:35` (`test_help` parameter `hf_sample`: drop the parameter; the new
`tests/policy/test_sampler.py` runs `python -I -m objdel.policy.sampler --help`)". Replace "This step switches it"
with "The E4 step 3 commit switches it", so the switch, the archive and the test fixes are one commit.

### 3. minor: after the switch, `VLLM_VENV` and `vllm_python` are dead code with three tests, and the HF carry-over list is not written.

Sources: E (E2 name table `sft/select.py`, E4 step 1), CLEAN (CLEAN16 step 2).

- Dead helper: E.md E4 step 1 makes `sampler_cmd` use `sys.executable`. Today `vllm_python()` is called only by
  `sampler_cmd` [code: experiments/phase2/select_ckpt.py:88-90, 94], and `VLLM_VENV` only by `vllm_python`
  [code: select_ckpt.py:55, 89]. E2 still moves both to `sft/select.py` (E.md name table row `sft/select.py`). Their
  tests: `tests/test_p2_scripts.py:118-124` (`test_sampler_uses_the_second_venv_when_the_image_has_one`) and
  `tests/test_vast.py:1248-1263` (two tests; line 1255 also asserts `sampler_cmd(...)[0]` is the second-venv python, so it fails on the switch commit and is a second test E4 step 1 misses) [code]. Owner rule: "Remove dead code." CLEAN12's vulture scans `src, experiments,
  tests`, so test use hides it; CLEAN16 step 2 (vulture without `tests`, advisory) lists it at $0 and forces a
  decision. So minor.
- HF carry-over: step 1 lists the `sample_vllm` items it carries [code: sample_vllm.py:17,34-101] but for
  `hf_sample.py` says only "`HfSampler` ... (body of hf_sample.py)" with `device: str = "cpu"`. `hf_sample.py` has
  `--batch_items` (default 8), `--timing` (read by `hf_rate.py` until CLEAN8 archives it), `todo` (resume of whole
  batches), `load` (bf16 on the GPU) [code: hf_sample.py:34, 82, 97-103]. The HF path stays live after E: `select_ckpt`'s
  G-hf-check forces `sampler: hf` on the selected branch [code: select_ckpt.py:305-306], and `r1_qwen3vl2b` is one of
  the four ported specs. "A case that step 1 carried over runs on `objdel.policy.sampler`; the rest moves to archive"
  then leaves the engineer to decide which `test_p2_hf_eval.py` cases (batching, resume, bf16-on-GPU) are "carried".
  No funded rental runs SFT selection after E, so minor.

Fix (E4 step 1): (a) "The E4 step 3 commit deletes `VLLM_VENV` and `vllm_python` from `sft/select.py` and the two
tests that only test them (`test_p2_scripts.py::test_sampler_uses_the_second_venv_when_the_image_has_one`,
`test_vast.py::test_sampler_python_uses_the_vllm_venv_when_present` and `::test_sampler_python_else_the_running_python`); the E-G16 record lists them as deleted." (b) Add an HF carry-over sentence:
"From `hf_sample.py` [code: hf_sample.py:34-103]: `load` (bf16 on CUDA), `records`, `generate`, `todo` (whole-batch
resume) and the flags `--batch-items` and `--timing`; `HfSampler(device=...)` defaults to `cuda` when
`torch.cuda.is_available()`. The `test_p2_hf_eval.py` cases for batching, first-stop cut, resume and bf16 load move to
`tests/policy/test_sampler.py`."

### 4. minor: `RlStage` cannot express B's `train_rl --merge-last` step, and the `--opsm-start-delta` cross-flag rule has no home after E6.4.

Sources: E (E2.1 `RlStage`, E6.0 spike, E6.1 `FLAG_ORDER`, E2.3 `run_rl`), B (B1.5 step 3b, B3.9, B2.2).

- B.md B1.5: "3b. `train_rl --merge-last --model qwen3-vl-2b --ckpt {rd}/ckpt --out {rd}/b1-rows` (a hard step) ...
  In `--merge-last` mode `train_rl` reads no rows and trains" nothing. B3.9 (round 14 rl-2) adds the same step for
  `b3_v5.yaml`. The step has no `--train` and no `--steps`.
- E.md E2.1: "`class RlStage(StageBase, TrainerFlags ...)`: ... `train: str` ... `steps: int = Field(gt=0)`" (both
  required) and "every flag that A-D add to `train_rl` ... becomes a `TrainerFlags` or `RlStage` field". E.md has no
  `merge_last` and no text on `--merge-last` [measured: grep]. So `load_spec` on the ported `b1_rows` raises
  "Field required" on step 3b, and after E6.4 `cli_parse(RlStage, ["--merge-last", ...])` exits non-zero, which breaks
  the crash-recovery step of every later B run. E2.3 `run_rl` is "resume ..., train, merge LoRA, save `<out>/final`",
  with no merge-last branch.
- Round 14 rl-3: B2.2 `--opsm-start-delta X` is "only with `bf16-tis-opsm`, argparse error otherwise or for X <= 0".
  After E6.4 the CLI is a pydantic model, so that argparse check is gone unless it becomes a validator; E names none.
- Caught-before-spend: the E6.0 spike ("must parse every golden-5 argv list of the `rl` steps"), E-G9 and B's I49 /
  U45 fail at $0 and point at the step. So minor.

Fix (E.md E2.1 `RlStage`, E2.3, E6.1): add `merge_last: bool = False`; make `train: str | None = None` and
`steps: int | None = Field(None, gt=0)`; extend `_save`: "with `merge_last`, `train`, `steps`, `resume` and
`rows_*` must be unset; without it, `train` and `steps` are required". Add `opsm_start_delta: float | None =
Field(None, gt=0)` beside `numerics`, with the validator "set only when `numerics == 'bf16-tis-opsm'`" (pure
intra-stage, so `load_spec` may hold it, unlike the session-spend rule). `FLAG_ORDER["rl"]` places `merge_last` and
`opsm_start_delta` where B1.5 / B2.2 write them. E2.3: "`run_rl` has two modes: `merge_last` calls
`merge_last_checkpoint(stage)` (B1.5's body, I49) and returns; otherwise the train path". Add both cases to U49 and
golden 5 (`b1_rows` step 3b string for string).

### 5. minor: `LIVE_CLIS` is a fixed list of 14 CLIs, but the E8 gate requires every live `-m experiments.*` string to be in it, and A-D add more.

Sources: E (E8 step 1, step 2, Gate).

- E.md E8 step 2: "The live CLIs: build_data, build_targets, fetch_lidc_pixels, prepare_rl, build_views,
  build_eval_views, train_sft, select_ckpt, train_rl, g_probe, final_eval, preflight, hf_push/pull; `LIVE_CLIS` (step
  1) lists their module names". E8 Gate: "every `-m experiments.*` module string in `src/`, `specs/`, `Makefile` and
  live `experiments/` resolves (`find_spec`) and is in `LIVE_CLIS`".
- A-D add live CLIs run as modules: `experiments.phase2.ceiling` (A.md A1 CLI), `check_series_z`, `nsys_baseline`,
  `outline_lengths` (A), `fp16_check`, `b2_attribution`, `b2_summary` (B), `perf.bench_checker` (D: called by the
  `d_gpu_test` box script, a live plan while V3 has not run) and `perf.make_d_goldens` (D, round 14).
  `experiments.phase1.bench.xarch_hashes` is also allow-listed in `vast/policy.py:476`.
- Caught-before-spend: the E8 gate and `test_cli_smoke` fail at $0 and name the module. So minor.

Fix (E8 step 2): "`LIVE_CLIS` = the 14 above plus every `experiments.*` module that a live spec, the Makefile, a
live box script or a live CLI names on the E8 commit (re-grep then; at plan time: `ceiling`, `check_series_z`,
`nsys_baseline`, `outline_lengths`, `fp16_check`, `b2_attribution`, `b2_summary`, `perf.bench_checker`,
`perf.make_b_goldens`, `perf.make_d_goldens`). A module whose plan is finished is archived instead (CLEAN8 rule)."

## Counts

blocking 0, major 1 (item 1), minor 4 (items 2-5).
