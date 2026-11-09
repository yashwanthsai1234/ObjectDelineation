# Final review round 6: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc: `docs/implementation2.0/Implementation.md`
(built from `docs/implementation2.0/*.md`). Repo at be52b1c. I checked every item against the round 1-5 entries in
`docs/implementation2.0/review-log.md` and `reviews/final_r5_arch*.md`. No item below repeats a fixed one. Item 10
collects round 5 minors that are still open.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item needs an NVIDIA product skill (this round makes no new CUDA, NIM or driver claim). I recommend none and install none. |
| engineering:code-review | I opened every code line I quote: `train_sft.py:65-111,150-157,257-305`, `policy.py:295-330,420-536`, `final_eval.py:44-80`, `env_check.sh:1-12`, `bundle.py:40-75`, `configs/r1_qwen3vl2b.yaml`, `configs/r5_shared.yaml`, `tests/test_p2_scripts.py:300-370`, `runs/r5/eval/eval_val/final_eval.json`, `docs/research/notes/final_eval_2026-11-08.md:3`. |
| engineering:debug | I traced the post-E paths by hand, then ran them in a scratch venv (pydantic 2.14.3, pydantic-settings 2.16.2, the box-lock version): spec -> `to_argv` -> `CliApp.run` (item 2); `train` -> `resolve` on a frozen stage (item 3); V4 launch -> `env_check.sh` -> `final_eval` -> push -> teardown (item 4). |
| engineering:system-design | I checked the run-spec seam end to end: every ported plan step must have a stage type that can hold its argv. That found item 1 (A3's repeated `--view`). |
| engineering:architecture | I checked layer placement: `core/spec.py` cannot import `sft` for the `data_root` default (item 7). I checked the E.3.1 and E.3.2 diagrams against the E6.4 text (item 9). |
| engineering:documentation | I ran DOC4 check 16 and the check 14 scan on the built doc. Check 16 passes. Check 14 fails on the plan's own text (item 5). |
| huggingface-skills:huggingface-best | This review makes no model choice, so nothing changed. |
| system-design-skills:system-design | I stress-tested the V4 rental as a run: inputs, outputs, failure path. A spec with nothing to push stops the `&&` line before teardown (item 4). |
| system-design-skills:service-decomposition | I confirmed the modular-monolith call: layers in one package, an optional E7 HTTP face. No item. |
| mattpocock-skills:codebase-design | I tested the stage interface as a seam: the CLI and the spec are two adapters of one stage model, so the round trip must exclude fields that belong to only one adapter (`name`, `soft`). That gives item 2's fix. |
| WebFetch / local probe | I checked pydantic-settings behaviour by running it, not from the docs: `CliSuppress`, `CliSettingsSource(..., cli_kebab_case=True, cli_implicit_flags=True)`, repeated flags into tuples, and the `CliApp.run` keyword form. Probes are in `(planning probe, not kept: r6arch)/probe*.py`. |

## Counts

blocking 1, major 4, minor 5.

## Items

### 1. blocking. EvalStage holds one view; the ported `b1_rows.yaml` step 4 has two views, and golden 5 must match it string for string

Sources: E, B, A.

- E.md "### E6", step 1 code: `class EvalStage(...)`: `kind: Literal["eval"]; ckpt: str; view: str; split:
  Literal["val", "test"]; out: str`, and "**Flag order** ... `to_argv` emits fields in `FLAG_ORDER[stage.kind]`".
- E.md "### E6", step 2: "port the live RL plan `b1_rows.yaml` (B1.5 ...)". Step 7 and I13: "`to_argv` on the ported
  `b1_rows`, `rl_full` and `r5_shared` specs equals, string for string, golden 5".
- B.md "#### B1.5", `b1_rows.yaml` step 4: "One `final_eval` process (front item 3), with A3's repeated `--view` and
  `--out`: `--ckpt {rd}/b1-rows/final --train {data}/prepared --view {data}/eval_views/lidc-val-shifted --split val
  --out {rd}/eval_b1_val_shifted --view {data}/eval_views/lidc-val-edge --split val --out {rd}/eval_b1_val_edge`".
- A.md A3: "`--view` may be given more than once (`action="append"`), with `--split` per view in the same order ...
  `--out` may repeat too, one per `--view`".
- Code today: `final_eval.py:47-50` has one `--view`, one `--split` and one `--out` (no `append`). E6.1 copied that
  form ("[code: final_eval.py:48-56]"), but A3 changes it before E.
- Measured ((planning probe, not kept: r6arch)/probe5.py): pydantic-settings parses the interleaved argv
  `--view a --split val --out oa --view b --split val --out ob` into `view=('a','b') split=('val','val')
  out=('oa','ob')` when the three fields are tuples. With `view: str` the second view cannot be held at all.

Problem: two parts that both claim to be the rule disagree. E6.1 defines a single-view `EvalStage` and a field-order
`to_argv`. E6.7 and I13 require the ported `b1_rows` spec to reproduce a step with two interleaved
`--view/--split/--out` groups. A field-by-field `to_argv` writes all views, then all splits, so it cannot give the
interleaved B1.5 order. The engineer cannot port `b1_rows.yaml` (the live RL plan) and keep golden 5 equal.

Fix (E.md E6 step 1 code, the Flag order paragraph, and U49):
- `view: tuple[str, ...] = Field(min_length=1)`, `split: tuple[Literal["val", "test"], ...]`,
  `out: tuple[str, ...]`, with a validator: `len(split) == len(view)`, and `len(out)` is 1 or `len(view)` (A3: one
  `--out` gives `<out>/<view name>/`).
- New rule in the Flag order paragraph: "`FLAG_ORDER["eval"]` has one entry `views` where the three tuples go. At that
  place `to_argv` writes, for each index i, `--view view[i] --split split[i]`, then `--out out[i]` when `out` has the
  same length (else one `--out out[0]` after the last group)." The `--train` token stays in golden 5 through E4's
  dropped-flag map.
- `final_eval` (E4 step 3) keeps A3's loop over views: one `VllmSampler` build, one folder per view.
- Add a test to E6.7: `to_argv` of the `b1_rows` eval stage equals B1.5 step 4 string for string (minus `--train`),
  and `CliApp.run(EvalStage, ...)` on that argv gives the same stage back.

### 2. major. The CLI round trip cannot pass: `kind` and `name` are required, but `to_argv` never writes them

Sources: E.

- E.md "### E2", step 1: `class StageBase(...)`: `kind: str`, `name: str   # unique per spec (E6 rule)` (no default);
  E6 step 1: `class SftStage(StageBase, ...)`: `kind: Literal["sft"]; ...` (no default).
- E.md E6 step 1 Flag order, and U49: "each `FLAG_ORDER[kind]` names every non-base field of its stage exactly once".
  The SFT argv rule: "`FLAG_ORDER["sft"]` = the field order of `RunArgs`". So `kind`, `name` and `soft` are never
  written.
- E.md E6 step 0 and I14: "requires `CliApp.run(SftStage, cli_args=to_argv(stage)[2:]) == stage`"; I14: "old flags
  give stage `S1`, `--spec P --stage NAME` gives `S2`, `S1 == S2`". E6 step 4: `main` becomes
  `CliApp.run(SftStage, cli_args=argv)`, and the box chain uses the flag form.
- Measured ((planning probe, not kept: r6arch)/probe.py, pydantic-settings 2.16.2): `CliApp.run(SftStage, cli_args=[--model ... --vision ...])`
  exits 2 with "the following arguments are required: --kind, --name". The same happens on the box for the
  `train_sft` stage that `to_argv` builds. (`RlStage` has `kind="rl"` as its default, but `name` is still required.)
- Measured ((planning probe, not kept: r6arch)/probe2.py): the form the spike text names, `CliApp.run(Model, cli_args=..., cli_kebab_case=True,
  cli_implicit_flags=True)`, passes the two keywords to the model as init data. On a frozen `extra="forbid"` stage this
  raises "cli_kebab_case: Extra inputs are not permitted".

Problem: as written, the spike, I14 and the preflight dry run (`test_cpu_smoke_train_then_select`, which now passes
`to_argv(stage)[2:]`) all fail. The engineer must invent how a stage gets `kind` and `name` from a CLI, and how to
compare a CLI stage with a spec stage. The spike's fallback ("one generator that builds argparse from the model
fields") has the same gap.

Fix (E2 step 1, E6 steps 0, 4 and 7, I14):
- `StageBase`: `kind: CliSuppress[str]`, `name: CliSuppress[str] = ""`, `soft: CliSuppress[bool] = False`. Each stage
  sets its literal as the default (`kind: CliSuppress[Literal["sft"]] = "sft"`, and the same for every kind).
- New `RunSpec._rules` rule: every stage of a spec sets a non-empty `name` (checked with `model_fields_set`). Add the
  case to U21.
- CLIs build the parser with `CliApp.run(Stage, cli_args=argv, cli_settings_source=CliSettingsSource(Stage,
  cli_kebab_case=True, cli_implicit_flags=True))`, not with keyword arguments.
- Every round trip (spike, I14, the E6.7 r1 train stage test, `S1 == S2`) compares
  `a.model_dump(exclude={"name", "soft"}) == b.model_dump(exclude={"name", "soft"})`.
- Measured: this form round-trips ((planning probe, not kept: r6arch)/probe3.py prints `True`).

### 3. major. E6.4's edit list leaves `train` dropping the stage that `resolve` returns

Sources: E.

- E.md "### E6", step 4, "**train_sft after E6**": "`resolve(stage) -> SftStage` returns
  `stage.model_copy(update={...})` with the resolved paths (`model_copy` builds a new instance ... so the frozen input
  never changes)"; "`train` rebinds `args = args.model_copy(update={"branches": ...})`". The list is introduced as
  "So the same change commit (R-1) makes these edits".
- Code: `train_sft.py:272` is `resolve(args)`, a bare statement that drops the return value. Today this works only
  because `resolve` assigns to `args` in place (`:152`, `:155`).
- Measured ((planning probe, not kept: r6arch)/probe4.py): with a frozen model and the `model_copy` form, `resolve(s)` without a rebind leaves
  `s.view == "views/x"`, and `s = resolve(s)` gives `/data/views/x`.

Problem: if the engineer makes exactly the listed edits, `train` keeps the relative `view`, `output_dir` and
`vision.reference`. On the box it reads `views/lidc-bundle-v3/manifest.json` relative to `/workspace/src` and
crashes. In the CPU smoke, `vision_reference` writes to `vision_ref/<model>.json` under the cwd, not under `run_dir`,
and the test's `ref.read_text()` fails, so preflight refuses every rental. The Mac test catches it, but the edit list
claims to be complete.

Fix (E6 step 4 list, new bullet after the `resolve` bullet): "`train` rebinds the result: `args = resolve(args)`
[code: train_sft.py:272 today drops it]". Extend the E6.7 `tests/sft/test_train.py` case: after `train` starts on a
frozen stage with relative `view` and `output_dir`, the manifest's `view.dir` is `Path(data_root, view)`.

### 4. major. V4 has no launch lines: no checkpoint fetch, no checkpoint path, and the V1 push line would skip teardown

Sources: 00_front, E.

- 00_front.md "### 6.3 V2, V3 and V4: only after a budget raise", V4 paragraph: "Plan file:
  `experiments/phase2/configs/e_v4.yaml` ... Then come `env_check.sh` (E-G13), the E.9 box test list, and the
  sampler-arm `final_eval` (pool size 1, as A3)." It gives no launch, sync, fetch or teardown lines (V1 and V3 have
  them: 6.2 launch lines, V3 item 0 "`fetch_hf.sh <run> rl-full-3ep` (to `{rd}/ckpt`)").
- E.md "## E.8", V4 list: "1. `env_check.sh` (E-G13 ...)"; "3. the E4 sampler arm: `final_eval` on
  `lidc-val-shifted`, K 8, seed as R5". E4 gate: "checkpoint = that file's `ckpt` key, the R5 final model". U51 checks
  neither the checkpoint nor its path.
- Code: `env_check.sh:9` `MODEL=${1:?usage: env_check.sh MODEL_DIR [OUT_JSON]}` (required). `EvalStage.ckpt: str` is
  required (E6.1). `runs/r5/eval/eval_val/final_eval.json` `ckpt` = `/workspace/runs/r5/final_model`, a path on the
  old R5 box. `docs/research/notes/final_eval_2026-11-08.md:3`: "Model: the full RL run's final model (HF
  `svenk029/objdel-sft`, revision `rl-full-3ep`)". `cmd_fetch_hf` pulls a revision to `/workspace/runs/<id>/ckpt`
  [code: policy.py:534-535]. Nothing in V4 fetches it.
- E.md "### E6", step 8b: "A spec with no such stage dies with "nothing to push"". `e_v4.yaml` has only `cmd` and
  `eval` stages. Step 8b itself notes that in front 6.2's line `pull.sh && push_hf.sh && teardown.sh`, a failing push
  skips teardown.

Problem: the engineer must guess where the TOST checkpoint comes from (the doc names a dead box path) and which path
`env_check.sh` and the eval read. A wrong guess (another revision) makes the TOST compare two different checkpoints,
which is a wrong result. A V4 run with no fetch stops at `env_check.sh` after the build, on a paid box. Copying V1's
last line gives `push_hf.sh <run> e_v4.yaml`, which dies by design, and `&&` then skips teardown.

Fix:
- 00_front 6.3 V4: add launch lines, as 6.2 does:
  `launch.sh OFFER <run> --config experiments/phase2/configs/e_v4.yaml --cap-hours <...> --owner-approved`;
  `sync_up.sh <run> && fetch_hf.sh <run> rl-full-3ep` (to `{rd}/ckpt`);
  `run.sh <run> experiments/phase2/configs/e_v4.yaml`;
  `pull.sh <run> && teardown.sh <run> && sweep_orphans.sh` (no `push_hf.sh`: V4 trains nothing).
- E.8 V4 list: step 1 = `bash experiments/phase2/docker/env_check.sh {rd}/ckpt {rd}/env_check.json`. Step 3 =
  `EvalStage(ckpt="{rd}/ckpt", view=("{data}/eval_views/lidc-val-shifted",), split=("val",),
  out=("{rd}/eval_val_sampler",), k=8, seed=11)`.
- E4 gate: "checkpoint = HF `svenk029/objdel-sft` revision `rl-full-3ep` (the R5 final model,
  `docs/research/notes/final_eval_2026-11-08.md:3`); the trainer arm scored it at `/workspace/runs/r5/final_model`".
- U51: assert the eval stage's `ckpt == "{rd}/ckpt"`, the `env_check.sh` argv, and that no stage has kind `rl` or
  `sft`. A second assert: `push_target` on `e_v4.yaml` dies, so the V4 lines must not call `push_hf.sh`.

### 5. major. DOC4 check 14 fails on the plan's own text, so `make test` fails at the DOC4 commit

Sources: 00_front, D, DOC.

- DOC.md "### DOC4", check 14: "`test_no_scratch_or_local_paths`: 0 matches of `/[p]rivate/tmp`, `/[U]sers/`, ... in
  `doc_files` (plan files included, frozen files excluded) and in `docs/implementation2.0/**/*.md`".
- DOC.md "### DOC1", step 0: "replace each absolute scratchpad prefix and each short scratchpad prefix (`S[P]/impl/`,
  `S[P]/`) ... (This page writes those patterns with brackets, as regex classes, so that it does not trip its own
  scan and check 14.)" Only DOC.md uses the brackets.
- Measured (`grep -nE '/[p]rivate/tmp|/[U]sers/|docs/implementation2.0'` on the built files):
  `Implementation.md:67` (00_front "## 0. Words": "an absolute `/[p]rivate/tmp` or `/[U]sers/` path"), `:544` (00_front
  "## 5. Custom kernels index": "contains `/[p]rivate/tmp`"), `:889` (00_front "## 7. Risk register", FR13), `:8498`
  (D.md "### D1", "contains `/[p]rivate/tmp`"), `:9254` (D test S06). `review-log.md:30, 77, 886, 917`
  (`/[p]rivate/tmp`) and `:565` (`(planning probe, not kept: b_r4)`) are under `docs/implementation2.0/**/*.md` too.
- These mentions name the string. They are not scratchpad prefixes, so DOC1's replace rule does not remove them.

Problem: an engineer who follows DOC1 and then DOC4 gets a red check 14, so `make test` fails (owner rule: it must
pass before a commit). Nothing in the doc says how to fix it.

Fix: in 00_front (sections 0, 5, 7), D.md (D1 step and test S06) and the source review-log rows, write the patterns
with brackets as DOC.md does (`/[p]rivate/tmp`, `/[U]sers/`, `S[P]/impl`). Add to DOC1 step 0: "also bracket every
literal mention of these patterns in the phase pages and their review logs; re-run the check 14 scan on the built
`Implementation.md` and `review-log.md` until it prints nothing".

### 6. minor. `SelectStage.vision` has no stated path form, and the select `--vision` JSON spacing differs from the SFT rule

Source: E (E6 step 1, step 8b). Today `chain()` builds `{"prepared": f"{data}/{v['prepared']}", "reference":
f"{rd}/{v['reference']}"}` with default `json.dumps` spacing [code: policy.py:437-440]. E6.8b's `vision_refs` adds
`/workspace/runs/<run_id>/` to `reference`, so it reads a relative reference. But `to_argv` fills only `{rd}`/`{data}`
placeholders, and the SFT rule writes compact JSON. Fix: in E6.1, say that `VisionRef` holds paths relative to the
data root and the run folder (as `SftStage.vision` does), and that `to_argv(SelectStage)` adds the `{data}/` and
`{rd}/` prefixes and writes `--vision` with default `json.dumps` spacing (golden 5 form).

### 7. minor. `SftStage.data_root`'s default needs `BOX_DATA`, which E2 moves to `sft/train.py`

Source: E (E6 step 1, E2.6 table). "same names and defaults" as `RunArgs`. The default factory reads `BOX_DATA`
[code: train_sft.py:71,82-84]. E2.6 puts `BOX_DATA` in `sft/train.py`, and `core` imports nothing from objdel. Fix:
`core/spec.py` defines `BOX_DATA = Path("/workspace/data/lidc")`, and `sft/train.py` imports it from `core.spec`.

### 8. minor. `to_argv` drops `vision: False`

Source: E (E6 step 1 SFT argv). "a `True` boolean as the bare flag ..., a `False` boolean never". `SftStage.vision`
may be `False` (`vision: off`). If it is dropped, the default `{}` applies, and `train` raises "vision: {} skips
G-vision" [code: train_sft.py:282-283]. No ported spec has `vision: off` today; the archived timing arm `r0_qwen3vl_arm.yaml:9` did. Fix: "`vision`
`False` is written as `--vision off`", plus a U22 case.

### 9. minor. Names in E.3.1 and E.3.2 differ from E6.4

Source: E. E.3.1 says "`train.py run_sft(stage)`", but E6.4 says "`sft.train.train(stage)`" (R-5 keeps `train`).
E.3.2 shows "`box: experiments/phase2/train_rl.py --spec ... --stage train`", but E6.4 says "The box chain uses the
flag form (`to_argv`)". Fix: `train.py train(stage)`, and `train_rl.py <flags from to_argv>` in the diagram.

### 10. minor. Round 5 minors still open

Sources: 00_front, E. These are unchanged since round 5:
- r5 arch-6: 00_front 2.2/2.3 still list `rl ... callbacks` and `infra/ paths.py` (00_front.md lines 254, 309, 312).
- r5 arch-7: E5's `fast.py` row still lacks `FINE`, `ROW_KEYS`, `OTHER_VISION`, `PAIR` and `INT`, and still turns
  `packed_chunk` into a function inside a move commit (E.md line 694).
- r5 arch-8: E.8 decision 1 still reads "Keep `Controller` + `eval/audit.py`" (E.md line 1645).
- r5 arch-10: `CliVariadicArg` is still not named (0 hits).
- r5 arch-11: the `cmd_run` strings at `policy.py:473-474,479` are still not covered.

Apply the round 5 fixes as written in `final_r5_arch.md` items 6, 7, 8, 10 and 11.
