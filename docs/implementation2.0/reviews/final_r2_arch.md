# Final review, round 2: architecture, code organization, layering, APIs, DOC, CLEAN

Reviewer role: principal software architect for production ML systems.
Document: `REPO/docs/implementation2.0/Implementation.md` (13,065 lines, round 2).
Earlier file read: `reviews/final_r1_arch.md`. Its items 1, 2, 5-13, 15-22 are fixed in the text and are not
repeated. Item 4 (ADR numbers) is fixed in front, D, E, DOC and CLEAN, but not in B and C: item 3 below covers only
what is still open. Item 14 (scratchpad) is fixed except in C: item 15 below covers only that rest.
Clean context kept: no transcripts, no memory folder, no `REPO/.claude/`, no git write, no rental, no REPO edit.

## Skills used

| Skill | What it changed in this review |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | No NVIDIA catalog skill covers package layering, run specs or ADR hygiene. No install, no recommendation. It only confirmed that the K3 build home (E.5) stays a pointer to `tilegym-monkey-patch-kernels-to-transformers`. Item 14 notes that E.0 calls K3 "Phase D's"; it is C1's fallback. |
| `engineering:code-review` | I checked each claim against the code at the working tree: `vast/policy.py:46-47,84-85,169-199,422-453,484-486`, `.gitignore:25-35`, `train/dataset.py:26-32,73-365`, `train/serve.py:49-52,290-302`, `g_probe.py:139`, `select_ckpt.py:374`, `final_eval.py:48-56`, the 3 RL plans. This found items 2, 7 and 12. |
| `engineering:debug` | I traced what fails at run time after a move, not at import: `home()` resolving to the package folder (item 2), the soft-fail `"~"` step that no stage model can express (item 7), and the A7 negative case that stops raising after E5 (item 6). |
| `engineering:system-design` | I checked the front time line against each gate's needs. E-G10 and E-G13 need a rental after E; the time line has none (item 10). |
| `engineering:architecture` | I mapped every ADR number in the text (grep of `ADR[- ]00NN`). B and C still use the old numbers (item 3). Item 17 covers who writes ADRs 0017-0020. |
| `engineering:documentation` | "One source": two names and two places for the rows pin (item 8), two archive folders for old plans (item 1), and two owners for E's ADRs (item 17). |
| `huggingface-skills:huggingface-best` | No leaderboard query: model choice is outside this review. I checked only that "a new model = one row" still holds after CLEAN8. It does not: CLEAN8 deletes `_LORA`, which `rl_lr` and the 4B steps read (item 9). |
| `system-design-skills:system-design` | The phase order is the main hypothesis. I checked that each E move table matches what A-D leave on the E0 commit. Items 1, 4 and 5 show three places where it does not. |
| `system-design-skills:service-decomposition` | A modular monolith with enforced layers is still the right shape. No item asks for a service split. Item 12 moves one torch-free loader into `verify`, so that the verifier owns its own data. |
| `mattpocock-skills:codebase-design` | Deletion test and seam check: `cache_loader` and `verify.store.loader` are two adapters for one job (item 12). The `rl.trl_adapter` interface changes in E5 without its callers' tests changing (item 6). |
| Web (WebFetch) | import-linter layers page: siblings joined by `:` are non-independent and may import each other; `|` and `:` cannot mix on one line [web: import-linter.readthedocs.io/en/stable/contract_types/layers/]. Contract types page: `ignore_imports` supports `*` and `**` wildcards, and `unmatched_ignore_imports_alerting` defaults to `error` [web: import-linter.readthedocs.io/en/stable/contract_types/]. Item 4 uses both. |

## Counts

blocking 3, major 8, minor 9 (20 items).

## Items

### 1. blocking | E (also B) | E ports and pins `rl_shifted.yaml`, but B1.5 archives it before E starts

Problem. B1.5 archives `experiments/phase2/configs/rl_shifted.yaml` to `archive/phase2_configs/` ("replaced by
b1_rows.yaml") (`Implementation.md:3240`, B.9 I22 line 4031, B.4 row line 3851). B runs before E (front order rule 2).
E still treats the file as a live plan:
- E0 golden 5 records its argv, and golden 6 records "the `train_rl` stage of `rl_shifted.yaml`" (lines 8186-8191).
- E3.2 starts at `specs/rl_shifted.yaml` (line 7982).
- E6.2 ports it, E6.7 and E-G9 compare it, and E6.8 sets `CONFIG ?= specs/rl_shifted.yaml` (lines 8630, 8651, 8660,
  8831).
- The E.9 sketch reads `ROOT / "specs" / "rl_shifted.yaml"` (line 9272).
On the E0 commit the file is missing, so goldens 5 and 6 cannot be recorded, and the Makefile default names a file
that does not exist. There are also two archive folders for old plans: B uses `archive/phase2_configs/`, while E6.6
and CLEAN8 use `archive/phase2_runs/configs/`.

Fix. In E, replace the live RL plan `rl_shifted` with `b1_rows` everywhere: E0 goldens 5 and 6(a), E3.2, E6.2 (port
list `b1_rows.yaml`, `rl_full.yaml`, `r5_shared.yaml`, `r1_qwen3vl2b.yaml`), E6.6 (the 12 archived names stay; add
"`rl_shifted.yaml` was archived by B1.5"), E6.7, E6.8 (`CONFIG ?= specs/b1_rows.yaml`), E-G9 and the E.9 sketch.
Use one archive folder: change B1.5, B.9 I22 and the B.4 row from `archive/phase2_configs/` to
`archive/phase2_runs/configs/`. B1.5 creates the README there with the four headings (Was / Replaced by / Commit),
and E6.6 appends to it.

### 2. blocking | E | E8 moves `vast/policy.py` into the package; the spend ledger, active runs and image allowlist then resolve to an empty folder

Problem. `home()` returns `VAST_HOME` or else `HERE`, the folder that holds the module [code:
experiments/phase2/vast/policy.py:46,84-85]. All Vast state files live under `home()`: `ledger.csv` (line 169),
`active_runs.json` (line 186), `.vast-allowlist`, `offers/<id>.json` (line 179 `load`), `known_hosts` (line 237) and
`code_files.txt` (line 395). The module also reads `HERE / "watchdog.sh"` (line 271). E8.1 moves the module to
`src/objdel/infra/vast/` and fixes only `ROOT` (lines 8673-8681). After the move, `home()` is
`src/objdel/infra/vast/`. There:
- `spend()` returns 0.0 when `ledger.csv` is missing [code: policy.py:198-199]. The `HARD_STOP` refusal and front
  stop rule 2 (90 % of the approved amount) then see no spend.
- `active()` returns `{}`, so `sweep_orphans` cannot see a live instance.
- `load(".vast-allowlist")` returns `{}` with no error [code: policy.py:178-181].
- `watchdog.sh` is not found.
`.gitignore:25-35` and both launchd plists also name `experiments/phase2/vast/` paths. No E test checks `home()`:
`conftest.py` fakes only `OBJDEL_SECRETS`. This is unsafe, because the budget and teardown rules depend on these files.

Fix. Add to E8.1, under "Repo root, one function": "`infra.vast` state stays where it is.
`home() = Path(os.environ.get("VAST_HOME", repo_root() / "experiments/phase2/vast"))`. `watchdog.sh`, the plists,
`.vast-allowlist` and the shell shims stay in `experiments/phase2/vast/`. The code reads `watchdog.sh` through
`home()`, not `HERE`. `.gitignore` is unchanged." Add three tests in `tests/infra/test_paths.py`:
(a) with `VAST_HOME` unset, `home() == repo_root() / "experiments/phase2/vast"`;
(b) `load(".vast-allowlist")` holds the key `image`;
(c) a launch with an empty allowlist or a missing `ledger.csv` refuses (today a missing ledger counts as $0 spend).
Add (b) to E-G15.

### 3. blocking | B, C (still open from round 1, item 4) | B and C still use the colliding ADR numbers

Problem. The front table assigns these numbers (lines 391-404): 0012 docs site (DOC10), 0013 number match (B2),
0014 reward v5 (B3), 0015 scheduled rows + numpy exception (B1). B still defines ADR-0014 = scheduled rows,
ADR-0012 = number match and ADR-0013 = reward v5 (B.3.4, lines 2945, 2966, 2975). B uses these numbers about 40 times
(lines 2747-4691). C cites "ADR-0012 option F" for the B2 number-match decision (lines 5300, 5688, 6573). DOC10
creates `0012-docs-site-with-zensical.md` (line 10219). DOC4 check 4 fails on two `0012-*.md` files. "Record it in
ADR-0013" is ambiguous.

Fix. In B (lines 2739-4700), replace in this order: `ADR-0014` -> `ADR-0015`, then `ADR-0013` -> `ADR-0014`, then
`ADR-0012` -> `ADR-0013`. Also change the combined forms on lines 2747, 2763, 3761, 3850 and 3871 (for example
`ADR-0012/0013/0014` -> `ADR-0013/0014/0015`). In C, lines 5300, 5688 and 6573: `ADR-0012` -> `ADR-0013`. Then
`grep -nE 'ADR-?0012' Implementation.md` must hit only DOC lines.

### 4. major | E (also B) | E0's interim import contract cannot pass on the E0 commit: B1's data modules create package cycles

Problem. E0 writes `layers = ["train", "eval", "reward | targets", "data", "checker | geometry"]` with one ignore
(line 8113). It also says E0 "re-derives the graph on the E0 commit". B1 lands before E0 and adds three upward edges
from `data`:
- `data/windows.py` and `data/quality.py` import `targets.frame` and `targets.augment` (B1.1 steps 2-3; B1.3 "reuse
  `augment` and `frame`").
- They import `reward.cache.rl_eligible` (B1.1 step 1).
- They import `objdel.eval.ceiling.loo_caches` (B1.3 DG14, line 3150).
`targets -> data` and `eval -> data` (`eval/probes.py`) already exist. So `data <-> targets` and `data <-> eval` are
cycles. No ordering of a `layers` list can hold a cycle. `|` siblings are independent, so `reward | targets` fails
too. "Derive from the graph" gives the engineer no rule that passes. The E3c merge (targets into data) and the E3b
move (ceiling into verify, below data) end these cycles. Until then E0 needs an explicit escape.

Fix. In E0, after "the list below is the be52b1c graph", add: "B1 adds the edges `data -> targets`, `data -> reward`
and `data -> eval.ceiling`. The E0 contract therefore lists them as named ignores:
`ignore_imports = ["objdel.eval.vision -> objdel.train.dataset", "objdel.data.windows -> objdel.targets.**",
"objdel.data.quality -> objdel.targets.**", "objdel.data.quality -> objdel.eval.ceiling",
"objdel.data.* -> objdel.reward.**"]` (wildcards are supported; an unmatched ignore is an error by default).
E3b removes the `eval.ceiling` ignore (ceiling moves to verify, below data). E3c removes the `targets` ignores
(targets merges into data). The `reward` ignore goes when `reward` becomes `verify` (E3b)." Also add negative test
(e): a planted extra `objdel.data.lidc -> objdel.targets` import (not on the ignore list) makes `layers-interim` fail.

### 5. major | E (also A, B, C) | The E5 `dataset.py` table is stale: A8 renames two names, B1.4 adds one, and one name gets a home that breaks `rl | sft`

Problem. A8 renames `dataset._tokens` -> `target_tokens` and `_voting` -> `voting`, with "no aliases kept" (lines
1413, 1852). B1.4 adds a public `dataset.window_record(out, item, window)` that calls `_prepare` (line 3181). E5 still
places `_voting` in `policy/rows.py` and `_tokens` in `sft/rows.py` (lines 8399-8418), and it has no row for
`window_record`. `old_names.json` is recorded on the E0 commit, so `test_move_complete` (E-G16) fails on three names.
There is a second fault. C11 builds each row's cap from `dataset.target_tokens` (line 5481) on the RL row path. If
`target_tokens` sits in `sft/rows.py`, `rl` must import `sft`, and the `rl | sft` contract forbids that.

Fix. In the E5 move table: `policy/rows.py` takes `voting` (not `_voting`) and `target_tokens` (not `_tokens`, and
not in `sft/rows.py`); `data/prepare.py` (E3c row) takes `window_record` beside `_prepare`. Change "all 29 top-level
names ... at be52b1c" to "the top-level names on the E0 commit (29 at be52b1c, plus `window_record`, with A8's two
renames)". Add to E5: "`sft/rows.py` and `rl/rows.py` import `target_tokens` from `policy.rows`."

### 6. major | E (also A) | E5 changes A7's TRL contract, but A's count test, an A7 negative case, the coupling page and E-G5 are not updated

Problem. E5 step 3 removes `shuffle_sequence_dict` from A7's `GLOBALS` (line 8467). E5 step 2 adds the attributes
`_step`, `_buffered_inputs`, `num_iterations` and the `SOURCE_SHA256` entries. A's S07 asserts exactly 29 items, with
`GLOBALS (1)`, and that `docs/reference/trl-coupling.md` has 29 rows (lines 2119-2121, 1408). One of A7's 5 negative
cases is "remove `shuffle_sequence_dict` -> raises" (A7 step 5). After E5 that case no longer raises, so the test
fails. E-G5 counts ">= 7 negative cases (A7's 5 ... plus the source-hash and attribute cases)" (line 8828), but only
4 of A7's 5 cases stay valid.

Fix. Add to E5 step 3, in the same change commit:
(a) replace A7's negative case "remove `shuffle_sequence_dict`" with "change TRL's `_prepare_inputs` source ->
raises with `_prepare_inputs` (SOURCE_SHA256)";
(b) S07 reads its expected total from `len(trl_adapter.ALL_ITEMS)` (a tuple that E5 adds), not the literal 29, and
checks that `trl-coupling.md` has that many rows;
(c) the E row of DOC9 adds "`reference/trl-coupling.md` (new item count)".
Change E-G5 to ">= 7 negative cases: A7's 4 that still apply, the `_prepare_inputs` source case, and the `_step` and
`_buffered_inputs` attribute cases".

### 7. major | E | No stage model can express the soft-fail `"~"` step, so golden 5 cannot round-trip

Problem. `chain()` passes `steps:` lists through literally, and `cmd_run` treats a leading `"~"` as "a failure is
logged and the chain goes on" [code: policy.py:424,484-485]. `r5_shared.yaml` uses it on a `train_rl` step and on a
`bash` step [code: configs/r5_shared.yaml]. A3's `a3_measure.yaml` and `checker_gpu.yaml` (the base of D's
`d_gpu_test.yaml`) use it too (line 998). None of `RlStage`, `EvalStage` or `CmdStage` has a field for it
(lines 8475-8490, 8591-8612). So `to_argv` cannot give back golden 5 for `r5_shared`, and E-G9 fails.

Fix. In E6.1, add `soft: bool = False` to a shared base `StageBase(BaseModel, frozen=True, extra="forbid")` with
`kind` and `name`. Every stage inherits it (`TrainerFlags` stays a mixin). `to_argv` puts `"~"` first when `soft` is
true, then the stage argv. Add to E6.7: "`r5_shared` stage 4 has `soft: true`, and its argv equals golden 5, `"~"`
included".

### 8. major | E (also B) | The rows pin has two names, two places and two checkers: `data_revision` (E6) vs `rows_revision` + `rows_manifest_sha256` (B1.4/B1.5)

Problem. B1.4 step 5 records `rows_revision` (the HF commit in the model repo `svenk029/objdel-sft`, under `rows/`) and
`rows_manifest_sha256` in `b1_rows.yaml`. `train_rl` asserts the manifest SHA on the box (B1.5, lines 3194, 3207-3215).
Front section 4 and 99_back decision 10 agree. E6 instead adds a spec-level `RunSpec.data_revision`, described as
"HF dataset commit", and only `infra.preflight` checks it on the Mac (lines 8610, 8621, 8640, 8651). The spec has no
`rows_manifest_sha256`, so `to_argv` cannot pass the value that `train_rl` must assert. E.R item arch 10 (line 9519)
used option 1 of round 1. The front had already switched to `rows_revision`.

Fix. In E6.1, delete `RunSpec.data_revision`. Add to `RlStage`: `rows_revision: str | None = None` (40 hex) and
`rows_manifest_sha256: str | None = None` (64 hex), with a validator (both set or both None). `to_argv` emits them as
the `train_rl` flags that B1.5 reads. `run_spec.json` records them. `infra.preflight` refuses when the local rows-root
manifest SHA differs. Change E6.5, E6.7, U21 and E.R (arch 10) to match. Replace "HF dataset commit" with "commit of
`svenk029/objdel-sft` that holds `rows/<root>/`".

### 9. major | CLEAN, E, 00_front | CLEAN8 (Wave 3, after E) deletes names that E1, E5 and front section 3 depend on, and uses pre-E names

Problem. CLEAN8 step 3 deletes `_LORA`, `_GRAD_CKPT`, `param_groups` and the float32 branch of `dtype` (line 11429).
E1.1 defines `rl_lr(model) = 5e-5 if model in _LORA else 5e-6`, and `tests/experiments/test_prereg_lr.py` loops
`MODELS` with `_LORA` (lines 8228-8239). E5 moves `_LORA`, `_GRAD_CKPT` and `param_groups` to `policy/models.py`
(line 8398). Front section 3, "If the gate fires", step 1 adds the 4B model to `_LORA` and `_GRAD_CKPT` (line 349).
CLEAN8 runs after E2 (header, line 11419) but names pre-E places ("Library: `MODELS = ...`", `preflight.py:35`,
`build_trainer(spec, ...)`, `RlSpec.max_completion_length`). E2 defines no `RlSpec` and no such field (line 8494).

Fix. Rewrite CLEAN8 steps 1 and 3 on post-E paths:
- step 1: the test cap is `flags.grpo={"max_completion_length": 32}` on `RlStage` / `build_trainer(..., flags=...)`;
- step 3: in `policy/models.py`, keep `_LORA` (as `frozenset(MODELS)`) and `_GRAD_CKPT`. Delete only the
  full-fine-tune branches. `rl.config.rl_lr` and `test_prereg_lr` keep working. `G_PROBE["lr_full"]` stays, because
  it is the frozen pre-registered record.
Change front section 3, step 1, to "`policy/models.py` (after E; `configs.py:73-77` before E)". Change CLEAN.9 line
12307 the same way as step 1.

### 10. major | E, 00_front | E-G10 (box part) and E-G13 need a rental after E, and the time line has none

Problem. E.8 says that the E4 box checks and E-G13 "piggyback on the next owner-approved Phase B run" (line 8866).
Front section 6.3 says "about $0.09 extra if its box checks ride a B run" (line 498). In the front time line, V1 and
V3 (all B runs) come before E (steps 3-7, lines 115-122). So no B run follows E, unless V3 slips. Yet E4 archives
`sample_vllm.py` and the trainer eval path in the same commit (line 8541). The only eval path left then has not run on
a box. Its binding deterministic gate (prompt ids, SamplingParams) cannot run on the Mac, because it needs vLLM.

Fix. Add a row to front section 6.3: "V4 (optional, after E): E-G13 + E4 deterministic check + E4 TOST sampler arm,
about 15 min, about $0.16 at $0.641/h; only with owner approval; it counts inside the $3 unless the cap is raised".
Add the same item to 99_back decision 1. In E4, add: "E4 is not done until E-G10's box part passes. Until then, keep
the trainer path, and archive `sample_vllm.py` / `hf_sample.py` in a follow-up commit after the box check."

### 11. major | E | E8.1 splits `vast/policy.py`, `preflight.py`, `bundle.py` and `hf_push/pull.py` with no name table, and E8.2 lists the same files as live CLIs that stay

Problem. E3 and E5 give a name-level move table and an E-G16 test. E8.1 only says "split ... by function: offer
terms/reasons with REGION; ledger/accrued; ssh/rsync/onstart; chain/code_files" for a 628-line module with 9
`cmd_*` subcommands and `main` [code: policy.py:289-627]. No home is given for:
- `ALARM`, `HARD_STOP` (DOC4 check 8 imports them through `POLICY_MODULE`);
- `IMPORT_CHECK`, `CPU_RAM`, `vision_refs`, `home`, `load`, `active`, `spend`;
- the `cmd_*` dispatcher.
The shell shims would call `python -I -m objdel.infra.vast`, which needs a `vast/__main__.py`. E8.1 moves
`preflight.py` and `hf.py (hf_push + hf_pull)` into `infra`, while E8.2 lists `preflight` and `hf_push/pull` as live
CLIs that stay in `experiments/phase2/`. E-G16's `old_names.json` covers `src` modules only.

Fix. Add an E8.1 move table:

| New module | Names |
|---|---|
| `infra/vast/offers.py` | `terms`, `query`, `reasons`, `REGION`, `CPU_RAM` |
| `infra/vast/ledger.py` | `FIELDS`, `SECRET_RE`, `ledger`, `load`, `active`, `accrued`, `spend`, `home`, `ALARM`, `HARD_STOP` |
| `infra/vast/remote.py` | ssh/rsync/onstart, `code_files`, `IMPORT_CHECK`, `LIVE_CLIS` |
| `infra/vast/chain.py` | `chain`, `vision_refs` |
| `infra/vast/__main__.py` | `cmd_*` and `main` |

State the rule "library body in `infra/`, thin CLI stays" for `preflight`, `hf_push` and `hf_pull`: each
`experiments/phase2/<x>.py` keeps only argument parsing and calls `objdel.infra.<x>`. Add
`experiments/phase2/vast/policy.py`, `preflight.py`, `bundle.py`, `hf_push.py` and `hf_pull.py` to `old_names.json`,
so E-G16 covers them. Set `POLICY_MODULE = "objdel.infra.vast.ledger"` in DOC13 step 3.

### 12. minor | E | `cache_loader` is a torch-free doctor-store loader, but E5 puts it in `policy`, beside a second loader in `verify`

Problem. `cache_loader` returns `ExpertCache` objects, checked against the manifest digests. It uses no torch or TRL
[code: dataset.py:290-302]. Torch-free code calls it too: `experiments/phase1/eval_run.py:90`, A2's ceiling CLI
(line 935), B1.4's tests, and E7's `Verifier` (which "loads the doctor store"). E5 places it in `policy/rows.py`
(line 8401), and E3b places `serve.loader` / `preload` in `verify/store.py`. That gives two adapters for one job, in
two layers. Verify-level code must reach up into `policy` for it.

Fix. In the E5 table, move `cache_loader`, `_caches` and `_voting`/`voting` from `policy/rows.py` to
`verify/store.py`. `policy/rows.py` and `rl/reward.py` import them from there. Note in the E3b `verify/store.py`
row: "the one doctor-store loader module (`preload`, `loader`, `cache_loader`)".

### 13. minor | E | E3c's prompt injection names a module that exists only after E5

Problem. E3c says that the CLI `prepare_rl` passes `policy.prompt.PROMPTS` (line 8343). E3 runs before E5, and
`policy/prompt.py` is created in E5.

Fix. In E3c, write: "`prepare_rl` passes `objdel.train.dataset.PROMPTS` until E5; the E5 move commit changes it to
`objdel.policy.prompt.PROMPTS`. The manifest bytes stay equal (golden 2 and U39)."

### 14. minor | 00_front, E, CLEAN | Stale names left after round 1

Problem. Five stale names remain:
- Front section 2.2, seam table (line 215), and FR5 (line 529) say `check_private_names()`. A7 and E define
  `check()` and `TRLContractError` (E.R line 9542).
- E6.2 and the E log list `b2_probe` (lines 8632, 9552). B renamed it `b2_fp16.yaml` (line 4595).
- E.0 calls K3 "Phase D's K3 kernels" (line 7780). K3a is C1's fallback (E.5).
- CLEAN8 step 1 and CLEAN.9 (lines 11423-11426, 12307) still say `RlSpec.max_completion_length` and
  `build_trainer(spec, ...)`.

Fix. Front: `check()` (A7) in both places. E: `b2_fp16`. E.0: "Phase C's K3a fallback". CLEAN: see item 9.

### 15. minor | DOC, C | DOC1 step 0 cannot map every cited scratchpad file; C still cites `S[P]/`

Problem. C still has about 20 `S[P]/` citations, plus a header that says the notes "are not in the repo yet" (lines
4703-4708, 4796-4806, 5147, 5385, 5431, 5445, 5476, 5701, 5735-5736, 6561). C.R3 line 6582 says that these paths were
removed. DOC1 step 0's copy table (lines 9763-9771) has no row for `impl/c_rev/` (`C.v1.md`, `C.spec.v1.md`,
`pass2_text_bench.py`, `render_bench.py`). C cites these files as evidence 4 times. The mechanical prefix rewrite
therefore has no target for them. E.0 cites `docs/implementation2.0/notes/nv_catalog.json` (line 7778), but DOC1 copies
that file to `docs/implementation2.0/notes/doc/` (99_back uses that path).

Fix. Add a DOC1 step 0 row: `impl/c_rev/{C.v1.md,C.spec.v1.md,pass2_text_bench.py,render_bench.py}` ->
`docs/implementation2.0/notes/c_rev/`. In C, replace `docs/implementation2.0/` and `S[P]/` with the tracked paths, and delete the
header lines 4707-4708. In E.0, write `docs/implementation2.0/notes/doc/nv_catalog.json`. DOC4 check 14 then reports 0.

### 16. minor | DOC, E | DOC4 checks 6 and 7 read `experiments/**/configs/*.yaml`; after E6 the plans live in `specs/`

Problem. `config_files(root)` globs only `experiments/**/configs/*.yaml` (line 9994). Check 6 asserts that each
config is documented. Check 7 scans the same glob for version strings. After E6, the live plans are `specs/*.yaml`.
Check 6 then passes with almost nothing to check, and the specs are never compared with `docs/reference/configs.md`.
DOC13 step 3 lists `POLICY_MODULE` and `FROZEN_GLOBS`, but not this glob.

Fix. In DOC4, set `config_files` to `specs/*.yaml` plus `experiments/**/configs/*.yaml`, and assert that the list is
not empty. Add `specs/*.yaml` to check 7's scan list. Add to DOC13 step 3: "after E6, `docs/reference/configs.md`
documents every `specs/*.yaml`".

### 17. minor | E, DOC | Two owners for ADRs 0017-0020

Problem. E.8 says "After E8, DOC writes ... ADRs 0017-0020" (line 8913), and E8.5 is a "hand-off to DOC" (line 8777).
E6.0 says "Record the result in ADR 0018" during E6 (line 8585). DOC9 says "each phase writes its own pages ... DOC
only owns the table and the checks", and its E row lists no ADRs (lines 10203-10214).

Fix. Make E the owner, as for B (B2 and B3 rows list their ADRs). Add to the DOC9 E row: "ADRs 0017-0020". In E.8 and
E8.5, write: "E writes ADRs 0017-0020 in the commit of the task named in the front table; DOC4 checks 4 and 5 hold
them."

### 18. minor | E | The E4 TOST plans a trainer-path eval that already exists

Problem. E.8 counts "one extra val-shifted eval, because both the trainer path and the sampler path must score the
same checkpoint" (lines 8866-8872). The trainer-path arm on `rl-full-3ep`, `lidc-val-shifted`, 1,084 items, K 8 is
already on disk: `runs/r5/eval/eval_val/items.jsonl` (made by `final_eval` through `make_trainer` +
`sample.sample` [code: final_eval.py:73-75]; front 6.2 uses it for Gate B1).

Fix. In E4's statistical gate, write: "trainer arm = `runs/r5/eval/eval_val/items.jsonl` (rl-full-3ep); the box runs
only the sampler arm on `rl-full-3ep`". In E.8, the added box time becomes about 2 min (engine start) plus 7 min
(one eval), with no second eval.

### 19. minor | 00_front | The target tree (2.3) leaves out what DOC and E8 keep

Problem. Front 2.3 lists `docs/` as "site, architecture/decisions (ADRs), reference, guides, runbooks,
implementation/" and lists `experiments/phase2/` only. DOC1 and DOC6 add `docs/plan/` (the live plan, moved from
`docs/research/research2.0/`), `docs/results/headline.json`, `mkdocs.yml` and `CITATION.cff`. E8.2's default keeps the
data CLIs in `experiments/phase1/`. E8.1 keeps `experiments/phase2/vast/` (shims and state, item 2).

Fix. In 2.3, add `docs/plan/`, `docs/results/`, `mkdocs.yml`, `CITATION.cff`,
`experiments/phase1/ (data builders, thin CLIs)` and `experiments/phase2/vast/ (shell shims, state: ledger,
allowlist, never shipped)`.

### 20. minor | B, E | `data/quality.py` starts from a seed that imports `transformers` at the top

Problem. B1.3 says "Start `quality.py` from ... `b1_tokens.py`". That seed starts with
`from transformers import AutoTokenizer` [measured: line 4 of the B1 seed script]. `quality.tokens(rows_root,
tokenizer, limit)` already takes the tokenizer as an argument. But a copied top-level import puts `transformers` under
`objdel.data`. That breaks `pure-core` (sources include `objdel.data`) and E-G4.

Fix. Add to B1.3: "`data/quality.py` imports no `transformers` or torch, not even lazily. `audit_data.py` and
`preflight` build the tokenizer and pass it in." Add a B1.3 test:
`python -I -c "import objdel.data.quality, sys; assert 'transformers' not in sys.modules"`.
