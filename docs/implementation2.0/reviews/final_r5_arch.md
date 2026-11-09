# Final review round 5: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc: `docs/implementation2.0/Implementation.md`
(built from `docs/implementation2.0/*.md`). Every item below was checked against the repo code at be52b1c and against the
round 1-4 review logs (`docs/implementation2.0/review-log.md`), so no item repeats a fixed one.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item needs an NVIDIA product skill (no CUDA claim is new in this round). Nothing recommended. |
| engineering:code-review | Opened every code line I quote: `preflight.py:49-56,121-190,244-251`, `policy.py:308-330,422-461,515-524`, `train_sft.py:76-105,152-155,257-260,276`, `fast.py:30-33,81,242`, `select_ckpt.py:170-262`, `final_eval.py:1-80`, `configs.py:108-114`, `tests/test_p2_scripts.py:306-345`, `runs/r5/eval/eval_val/*`. |
| engineering:debug | Traced each post-E path by hand: spec -> `cmd_launch` -> `refusal` -> DATA checks; spec -> `cmd_push_hf`; spec -> `vision_refs`; `SftStage` -> `train_sft.train`; CLEAN8 smoke -> select sampler. That found items 2-5. |
| engineering:system-design | Checked the run-spec seam end to end (who reads each spec field, which old reader still expects a `RunArgs` YAML). Items 2 and 3 come from this. |
| engineering:architecture | Checked layer placement against the E5/E2 tables and the front 2.2/2.3 copies (item 6), and decision consistency between front 6.3, 99_back 3b and E4 (item 1). |
| engineering:documentation | Checked "one source" rules: the front tree claims to follow E.3.1 but drifted (item 6); E.8 owner decision 1 drifted from CLEAN6 (item 8). |
| huggingface-skills:huggingface-best | No model choice in this review; no change. |
| system-design-skills:system-design | Stress-tested the spec as an interface: every consumer of the old config must have a named new input (items 2, 3). |
| system-design-skills:service-decomposition | Confirmed the modular-monolith call (layers in one package, E7 HTTP face optional) is right; no item. |
| mattpocock-skills:codebase-design | Deletion test on the `SftStage` seam: the frozen model cannot replace the mutable `RunArgs` without a body change (item 3). Seam check on `policy.sampler` (item 5). |
| WebFetch | pydantic-settings docs: `CliApp.run` accepts a plain `BaseModel`; list fields parse as repeated flag, JSON or comma form, and space-separated values need `CliVariadicArg` (minor item 10). |

## Counts

blocking 1, major 4, minor 6.

## Items

### 1. blocking. Front 6.3 and 99_back 3b keep the trainer eval path until V4; E4 deletes it before V4

Sources: 00_front, 99_back, E.

- 00_front.md, "### 6.3 V2: validation with the money left; V3 after the budget decision; V4 after E", V4 paragraph:
  "Until V4 passes, E4 keeps the trainer eval path."
- 99_back.md, decision 3b: "Until V4 passes, the box halves of E-G10 and E-G13 stay "open until V4", and E4 keeps the
  trainer eval path and the old samplers."
- E.md, "### E4", step 3: "`final_eval` calls `VllmSampler` + `evaluate`. It no longer calls `build_trainer` or
  `trainer._generate`. The same change commit deletes its `--train` flag". E.6 E-G10: "The E phase gate passes on the
  Mac half; E4 is done only after the box half". E.9 I22: "AST: `final_eval` imports no `build_trainer` and has no
  `_generate`." E.7: "`final_eval` keeps the trainer path in git history". E6.1 `EvalStage`: "no `train`: E4 (before
  E6) deletes final_eval's --train".
- Code: today the trainer path is `final_eval.py:73-75` (`g_probe.make_trainer(args, a.ckpt, train, 1)` then
  `sample.sample(trainer, ...)`).
- History: review-log front round 2 says "E4 keeps the trainer path until V4 passes" and lists "E4 trainer path until
  V4" as owned by E; E never made the edit.

Problem: two parts that both state the rule disagree. E removes the trainer path from `final_eval` on the Mac, before
V4, and E6 builds `EvalStage` on that removal. The front and back say it stays until V4 passes. An engineer cannot
follow both.

Fix (pick one; option A matches every other E task and costs nothing):
- A (recommended). In 00_front 6.3 V4 paragraph replace "Until V4 passes, E4 keeps the trainer eval path." with
  "Until V4 passes, E4 keeps the old samplers (`sample_vllm.py`, `hf_sample.py`); `final_eval` already uses the new
  sampler from E4 step 3, and the trainer path is in git history (the owner tags the last pre-E4 commit
  `pre-e4-eval`)." In 99_back 3b replace "and E4 keeps the trainer eval path and the old samplers" with "and E4 keeps
  the old samplers; the trainer eval path is in git history (tag `pre-e4-eval`)". In E4 step 3 add one sentence: "The
  owner tags the commit before this change `pre-e4-eval`".
- B. Keep the front rule: move E4 step 3 (switch and `--train` delete) into the V4 follow-up commit, give `EvalStage`
  a `train: str | None` field until then, and move I22's AST assert to the follow-up commit. (More edits in E4, E6,
  E.9 and E-G9.)

### 2. major. E6.8 makes preflight take a RunSpec but no DATA check has an input from an RL spec

Sources: E (E6 step 8).

- E.md, "### E6", step 8: "`preflight` takes a `RunSpec` (via `load_spec`). Data and view hashes come from the spec
  stages (the `train`, `train_views`, `view`, `probes` roots)." Tests: "preflight passes on `specs/b1_rows.yaml` with
  fixture data".
- Code: every DATA check takes a `RunArgs` and reads `a.data_root` and `a.view` (an SFT bundle with `rows.jsonl`):
  `view_integrity` (`preflight.py:101-104`), `patients` (`:126-128`), `coverage` (`:154-155`, reads `r["epoch"]`),
  `prompts` (`:170`), `ct_window` (`:181`), `targets` (`:194`); `CHECKS` at `:244-251`. `RlStage` (E2.1) has `train`,
  `train_views`, `out`, but no `data_root` and no `view`. The paths hold `{data}` placeholders that `chain()` fills
  with the box path `/workspace/data/lidc` (`policy.py:433`), not a Mac path.
- B's pre-E answer is a separate RunArgs launch config (`b1_launch.yaml` with `rows_plan`, B.md B1.3/B1.5, I19). E
  never names it (grep of E.md for `b1_launch`, `rows_plan`: 0 hits).

Problem: after E6.8 the only V1-style input (a RunArgs YAML) is refused, and for `specs/b1_rows.yaml` (one RL stage,
one eval stage) the six DATA checks have no `data_root`, no `view` and no Mac value for `{data}`. The engineer must
invent which checks run per stage kind and how `{data}` resolves on the Mac. A wrong guess either drops the leakage
and coverage checks (owner rule: dataset quality first) or makes every spec fail preflight, so no rental can launch.

Fix: add to E6 step 8 a table "check -> input per stage kind": for `rl` stages run `patients`, `coverage`, `prompts`,
`ct_window`, `targets` on the rows root `local_root(stage.train)` and each `train_views` root, using B1.3's
`preflight.local_rows_root` rule for `{data}` (`{data}/rows/<x>` -> `data/views/<x>`; any other form raises), plus
B1.3's "ROWS root" check against `rows_manifest_sha256`; for `sft` stages run them on `data_root/view` as today; for
`eval`/`select`/`probe` stages hash `view`/`probes` only. State `{data}` on the Mac = `OBJDEL_DATA_ROOT` or `data/lidc`.
Add a test: a spec with one `rl` stage runs exactly these checks on the fixture rows root and refuses a planted
patient overlap.

### 3. major. `SftStage` replaces `RunArgs`, but its argv form is undefined and the frozen model breaks `train()`

Sources: E (E2 step 6, E6 steps 0, 1, 7), E.9 (U22).

- E.md, "### E2", step 6: "`RunArgs` becomes `core.spec.SftStage` in E6." E6.1: `class SftStage(StageBase,
  frozen=True, extra="forbid"): kind: Literal["sft"]; ...  # fields of train_sft.RunArgs`. E6.0 spike parses "the
  `eval`, `select` and `probe` argv lists" (no `sft`). E6.1 `FLAG_ORDER`: "one tuple per kind: the order of the flags
  in today's plans". E6.7: "the train stage is compared by value (today it passes `--config <yaml path>` ...,
  and that file is archived)". E2.6 table: `main` stays as `TrlParser([RunArgs]) -> train`.
- Code: today the SFT stage argv is `["-m", "experiments.phase2.train_sft", "--config", <yaml>]` (`policy.py:443`);
  `main` is `TrlParser([RunArgs]).parse_args_and_config(argv)` (`train_sft.py:257-260`), which needs a dataclass.
  `train()` assigns to its argument: `args.view, args.output_dir = ...` (`train_sft.py:152`), `args.vision = {...}`
  (`:155`), `args.branches = [...]` (`:276`), and `RunArgs.__post_init__` sets `self.vision = False` (`:105`). A
  `frozen=True` pydantic model raises `ValidationError` on each of these assignments (U20 asserts it).

Problem: no text says what `to_argv(SftStage)` emits after the YAML is archived (flags in which spelling? a temp
YAML?), no spike or golden covers it (E6.7 compares values only), and the frozen stage cannot be passed to the
unchanged `train()` body. An engineer who follows E6 gets a crash at `train_sft.py:152` on the first SFT stage, or an
untested argv that fails on the box.

Fix: in E6 step 1 add "`to_argv(SftStage)` emits `-m experiments.phase2.train_sft` then each set field as
`--<field_name> <value>` in `FLAG_ORDER["sft"]` (snake case, the TrlParser spelling; lists and dicts as JSON)";
add `sft` to the E6.0 spike and a golden-5 row for it. In E6 step 4 add "train_sft's `main` builds `SftStage` with
`CliApp.run` and passes `sft.train.train(stage.model_copy())`; `train()` builds a local mutable `Run` dataclass (or
uses `model_copy(update=...)`) instead of assigning to its argument; lines 105, 152, 155, 276 change in that change
commit (R-1)". Also fix U22: `set(SftStage.model_fields) - {"kind", "name", "soft"}` (StageBase has three fields, E2.1),
and drop `stages` and `sampler` ("Vast runner only (policy.chain)", `train_sft.py:101-102`) from `SftStage`, since the
spec holds the stage list.

### 4. major. After E6, launch, push and the vision guard still read RunArgs keys; `b1_launch.yaml` has no fate

Sources: E (E6 steps 1, 6, 8; E8 name table).

- E.md, "### E6", step 8: "The stamp is keyed on the spec path and the spec's sha256. `policy.launch` calls
  `refusal(spec_path)`." E6.1: `class ComputeSpec(...): gpu: str; cap_hours: float`. E8.1 table moves `cmd_push_hf`
  to `infra/vast/__main__.py` and `vision_refs` to `infra/vast/chain.py`, unchanged.
- Code: `cmd_push_hf` reads `yaml.safe_load((ROOT / a.config).read_text())["output_dir"]` (`policy.py:517`); a
  RunSpec has no top-level `output_dir` (`extra="forbid"`), so it raises `KeyError`. `vision_refs` reads the
  top-level `vision` key (`policy.py:460`); in a spec it sits in `SelectStage.vision`, so the guard list is silently
  empty and `cmd_run` loses its "missing G-vision reference" refusal (`policy.py:472-475`). `cmd_launch` takes the
  cap from `--cap-hours` (`policy.py:321`); E names no reader of `compute.cap_hours` or `compute.gpu`.
- `b1_launch.yaml` (B1.5, the V1 launch config) is neither ported nor in the E6.6 archive list (grep E.md: 0 hits).
- Front 6.2 runs `pull.sh <run> && push_hf.sh <run> <config> && teardown.sh <run>`, so a crash in `push_hf`
  skips the teardown in that line.

Problem: the engineer must guess (a) the HF revision source for a spec, (b) where the vision guard reads, (c) which
cap wins (CLI flag or spec), (d) whether `b1_launch.yaml` is archived. Today's code gives a crash for (a), a silent
loss of a safety check for (b), and two cap sources for (c).

Fix: add E6 step 8b: "`cmd_push_hf` takes the revision from the spec's last `rl` or `sft` stage `out` (`output_dir`)
basename; `vision_refs(spec)` returns the `reference` of every `SelectStage.vision` that is a string; `cmd_launch`
reads `compute.cap_hours` and `compute.gpu` from the spec, and `--cap-hours` may only lower it (die if higher);
`b1_launch.yaml` joins the E6.6 archive list." Tests in `tests/infra/test_chain.py`: `vision_refs(specs/r1_qwen3vl2b)`
equals the old `vision_refs(r1_qwen3vl2b.yaml)` from golden 5's E0 record; `push_hf` revision for `specs/b1_rows`
is `b1-rows`; `--cap-hours` above `compute.cap_hours` dies.

### 5. major. CLEAN8 step 1 needs a token cap on the select sampler, but no stage or sampler has one

Sources: CLEAN (CLEAN8 step 1), E (E4 step 1, E6 step 1).

- CLEAN.md, "### CLEAN8", step 1: "If `select_ckpt` samples outside `build_trainer` on the post-E tree, the same
  commit gives its stage the same override and a test. Switch the `_tiny_vlm` default and the
  `test_cpu_smoke_train_then_select` parameters to `"qwen3-vl-2b"` only".
- After E4, select samples outside the trainer: `sft.select.sampler_cmd` runs `-m objdel.policy.sampler` (E4 step 1).
  `params_for(model, *, n, seed, stop_token_ids)` takes `sampling(model, "vllm")`, whose `max_tokens` is
  `MAX_COMPLETION[model]` = 1536 for Qwen3-VL (`configs.py:112-113`); it has no cap argument. `SelectStage` (E6.1)
  has no cap field and `extra="forbid"`.
- Code: `tests/test_p2_scripts.py:333-334`: "a random Qwen3-VL never stops, so its 1,536-token samples would take hours
  on a CPU". This test is preflight's `DRY_RUN` check (`preflight.py:34`).

Problem: the condition in CLEAN8 step 1 is always true after E4, and the step names no field, no CLI flag and no path
to the sampler child. The engineer must design the seam; without it the qwen-only smoke (and so preflight's dry run)
runs for hours.

Fix: in E4 step 1 add `max_tokens: int | None = None` to `params_for` (overrides `sampling`'s value when set) and a
`--max-tokens` flag to `policy.sampler.main`; in E6.1 add `max_tokens: int | None = None` to `SelectStage` (emitted by
`to_argv` only when set) and pass it through `sft.select.sampler_cmd`. Then CLEAN8 step 1 reads: "the smoke's select
stage sets `max_tokens=32`; its probe stage sets `grpo={"max_completion_length": 32}`".

### 6. minor. Front 2.2 and 2.3 drifted from E's move tables

Sources: 00_front, E.

- 00_front.md "### 2.2": `rl ... callbacks`; `core spec (RunSpec, pydantic)`; `infra paths, vast, ...`. "### 2.3":
  `core/ spec.py`; `rl/ ... callbacks.py`; `infra/ paths.py ...`; `data/` without `audit.py`.
- E.md: E5 table "`policy/callbacks.py` | `TraceStep` ... Not `rl/`"; E8.1 "There is no `infra/paths.py`"; E2.5
  creates `core/provenance.py`; E3c row `data/audit.py`.

Front order rule 1 says E's move tables are the only source of post-E paths, so an engineer is not misled, but DOC13
step 2 rewrites `package.md` "per layer" and the front is the readable copy. Fix: regenerate 2.2 and 2.3 from E.3.1
(callbacks under policy, `core/provenance.py`, `data/audit.py`, no `infra/paths.py`).

### 7. minor. E5's fast.py row misses five constants; `packed_chunk` changes shape in a move commit

Source: E (E5 table).

- E5 lists `BUCKETS, window, per_row, ...` and says "The table also covers `configs.py`, `fast.py`, ...". Code:
  `fast.py:31-33` `FINE`, `ROW_KEYS`, `OTHER_VISION` (used by `FastForward`, lines 99, 104, 139), `:81` `PAIR` (used by
  `packing_wins`), `:242` `INT` (used by `VisionCache`, line 268) have no row; `test_move_complete` fails on them.
- "`packed_chunk` as a function `packed_chunk(model, ids, ..., metrics: dict)`": it is a method today (`fast.py:164`),
  so this is a logic change inside the move commit, against R-1.

Fix: add `FINE, ROW_KEYS, OTHER_VISION` -> `rl/trl_adapter.py`, `PAIR` -> `rl/forward.py`, `INT` -> `policy/caches.py`;
move `packed_chunk` as a method in the move commit and make it a function in a separate change commit.

### 8. minor. E.8 owner decision 1 still archives `eval/audit.py`

Source: E (E.8). "Keep `Controller` + `eval/audit.py` (-> `rl/controller.py`) or archive them? [archive ...]" vs E3b/E3c
rows and CLEAN6 step 2: `hack_audit` stays and moves to `data/audit.py`. Fix: "Keep `Controller` and `alarms_for` (->
`rl/controller.py`) or archive them? [archive; CLEAN6; `hack_audit` stays and moves to `data/audit.py` in E3c]".

### 9. minor. V4 trainer arm: the `ckpt` key is in `final_eval.json`, not in `items.jsonl`

Source: E (E4 gate). "checkpoint = that file's `ckpt` key" where the file is `runs/r5/eval/eval_val/items.jsonl`.
Measured: `items.jsonl` lines hold `item_id, legal, iou_tc, f5, precision_vote, recall_vote, dice_vote` only;
`final_eval.json` holds `ckpt: /workspace/runs/r5/final_model`, `K` 8, `seed` 11, `items` 1084. Fix: "checkpoint =
the `ckpt` key of `runs/r5/eval/eval_val/final_eval.json`".

### 10. minor. Name the pydantic-settings option for `--train-views a b`

Source: E (E6 step 0). Today `train_rl --train-views` is `nargs="*"` (`train_rl.py:26-28`), so plans write space-separated
values. pydantic-settings parses list fields as repeated flags, JSON or comma form; space-separated values need
`CliVariadicArg` [web: pydantic.dev pydantic-settings docs, "CliVariadicArg"]. The spike would find it, but naming it
saves the fallback: add "`train_views: CliVariadicArg[tuple[str, ...]]`, or an alias" to E6.0.

### 11. minor. `cmd_run` messages still name the archived `r0.yaml`

Source: E (E6.6, E8.1). `policy.py:473-474` tells the user to run `experiments/phase2/configs/r0.yaml first`, and
`:479` says "module must be ... experiments/phase2/configs/<name>.yaml". E6.6 archives `r0.yaml` and moves plans to
`specs/`. Fix: in the E6.6 commit change both strings to `specs/<name>.yaml` and point the vision-reference message to
front section 3 step 0's `--vision_reference_only true` form.
