# Final review round 8: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc: `docs/implementation2.0/Implementation.md`
(built 23:39 from the phase sources, 1,608,894 bytes). Repo at be52b1c. I checked every item against the round 1-7
rows in `docs/implementation2.0/review-log.md` and `reviews/final_r7_arch*.md`. The round 7 fixes are in place:
arch-2 (the `SelectStage`/`ProbeStage` fields, U52), arch-3 (the V4 step 2 `-m, pytest` argv, U51), arch-4
(`EvalStage.ckpt` optional, `--hash-view` read before `cli_parse`, U53) and arch-5 (the kernel-folder scrub in DOC1).
No item below repeats a fixed item. Scratchpad path patterns in this file are written with brackets (`S[P]/`), so
this file does not trip DOC4 check 14 after DOC1 step 0 copies it.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item makes a new NVIDIA product, CUDA or driver claim. I recommend no catalog skill and install none. |
| engineering:code-review | I opened every code line I quote: `select_ckpt.py:359-390`, `g_probe.py:127-157`, `train_rl.py:23-50`, `train_sft.py:68-112, 258-260, 388-390`, `final_eval.py:46-77`, `policy.py:422-490`, `tests/test_p2_scripts.py:306-386`, `tests/test_p2_train_rl.py:26,64,91`, TRL 1.17.1 `trl/scripts/utils.py:293-326` (`TrlParser`), `pyproject.toml`. |
| engineering:debug | I traced two paths by hand: V4's `run.sh <run> experiments/phase2/configs/e_v4.yaml` through `chain()` after E6 (item 1), and the pool names from D7 step 6 through E3b's R-2 branch to the `--procs` defaults (item 3). |
| engineering:system-design | I walked the V4 rental in time order (launch, sync, fetch, `run.sh`): the first code that reads the plan path with the E6 chain runs after about 38 min of paid box time, so a path rule in `chain()` costs money (item 1). |
| engineering:architecture | I compared front section 2 (target layers and tree) with E's move tables, as two texts that both describe the end state (item 4). |
| engineering:documentation | I ran the DOC4 check 14 patterns on the files that DOC1 step 0 copies (`reviews/`, `notes/`, `tests/`, `*.spec.md`): 118 lines in 50 files, and 35 lines in 27 files have no row in the step 0 table (item 6). |
| huggingface-skills:huggingface-best | No model choice in this round; nothing changed. |
| system-design-skills:system-design | Stress test of each seam on the failure path: a test-only cap that has no field (item 5), a field whose default calls a deleted function (item 3). |
| system-design-skills:service-decomposition | I re-confirmed the modular-monolith call (one package, layers enforced by import-linter, optional E7 HTTP face). The front tree places `TraceStep` in `rl`, which would break the `rl | sft` independence rule (item 4). |
| mattpocock-skills:codebase-design | I tested each stage model as the interface of its CLI, with two adapters (CLI and spec). U52 asserts an exact name match, but the `TrainerFlags` mixin gives some stages fields their CLI never had (item 2). |
| WebSearch/WebFetch | Loaded. Not needed: every claim below is checked against local code (TRL 1.17.1 is installed in `.venv`). |

## Counts

blocking 1, major 5, minor 6.

## Items

### 1. blocking. I29 and E-G9 make `chain()` read only `specs/*.yaml`, but V4's `RunSpec` lives in `experiments/phase2/configs/`

Sources: E, 00_front.

- E.md "## E.9" table, I29: "`policy.chain()` on an old-schema YAML ... raises a `ValidationError` ...; the 12 old
  plans are not in `experiments/phase2/configs/`; `chain()` reads only `specs/*.yaml`."
- E.md "## E.6", E-G9: "the box chain reads only `specs/*.yaml` (the 4 ported plans plus the live A-D plans)".
- E.md "### E6" step 8b, "V4's file": "E writes `experiments/phase2/configs/e_v4.yaml` (E.8, front 6.3) after E6.
  After E6, `chain()` reads only specs (step 3), so the file is a `RunSpec` at the path front 6.3 names". E.8 V4 and
  U51 name the same path.
- 00_front.md "### 6.3", V4 launch lines: `run.sh <run> experiments/phase2/configs/e_v4.yaml # steps 0-3 of E.8`.
- Code: `cmd_run` calls `chain(ROOT / a.module, a.run_id)` on the Mac [code: experiments/phase2/vast/policy.py:469-470],
  after `launch.sh`, `sync_up.sh` and `fetch_hf.sh`. `make preflight` reads the spec with `load_spec` (E6.8), and U51
  reads it with `load_spec` and `to_argv`; neither calls `chain()`.

Problem: two rules both claim the box chain. To make I29 pass ("`chain()` reads only `specs/*.yaml`"), `chain()`
must refuse a path outside `specs/`. Then V4's `run.sh` dies after launch, sync and the checkpoint download, which
front 6.3 prices at up to 38.2 min of the 54 min ($0.71) rental. No Mac test catches it, because preflight and U51
do not call `chain()`. The rental is paid and gives no E-G10 or E-G13 result.

Fix (E6 step 8b "V4's file", E.8 V4, U51, front 6.3 V4 launch config and launch lines, 99_back decision 3b): move the
file to `specs/e_v4.yaml`, the one folder for live plans after E6 ("the box chain reads only `specs/*.yaml`").
Launch lines: `make preflight CONFIG=specs/e_v4.yaml`, `launch.sh OFFER <run> --config specs/e_v4.yaml ...`,
`run.sh <run> specs/e_v4.yaml`. I29 adds: "`chain(specs/e_v4.yaml)` gives the U51 argv; `chain()` on a path outside
`specs/` raises before any `remote` call".

### 2. major. U52 asserts "argparse dest names = stage fields", but E2.1 and E6.1 define fields that the CLIs do not have, and drop dests the CLIs do have

Sources: E.

- E.md "### E.9.1", U52: "For each of `train_rl`, `g_probe`, `final_eval`, `select_ckpt`, `train_sft`: the set of
  argparse `dest` names on the E0 commit ... equals the stage fields minus `kind`, `name`, `soft`, minus a listed
  dropped set (`train` of `final_eval`, E4; `hash_view` of `final_eval`, read before `cli_parse`)."
- E.md "### E2" step 1: `class TrainerFlags`: `lr`, `shared_prefix`, `compile`, `plain`, `trace`, `grpo`;
  `class RlStage(StageBase, TrainerFlags, ...)`. E.md "### E6" step 1: `class ProbeStage(StageBase, TrainerFlags,
  ...)`; `class SftStage`: "the fields of train_sft.RunArgs ... minus the runner-only `stages`, `sampler` ... and
  B1.3's `rows_plan`".
- Code: `train_rl` has no `--compile` and no `--plain` [code: experiments/phase2/train_rl.py:23-36]; A6 adds only
  `--profile-window` and `--trace` to it (A.md A6 step 3), and no phase adds `--compile` or `--plain` (grep of A-D).
  `g_probe` has no `--lr` [code: experiments/phase2/g_probe.py:131-145]. `train_sft` parses with
  `TrlParser([RunArgs])` [code: train_sft.py:260], whose dests are the `RunArgs` field names, `stages` and `sampler`
  included [code: train_sft.py:101-102]; `TrlParser` adds no `config` dest (it pops `--config` from argv [code:
  .venv trl/scripts/utils.py:321-325]).

Problem: U52 fails on E's own models: for `train_rl` the fields hold `compile`, `plain` (not dests); for `g_probe`
they hold `lr` (not a dest); for `train_sft` the dests hold `stages`, `sampler`, `rows_plan` (not fields). The
engineer must choose: drop the mixin fields (and change `build_trainer` and P3), or widen the test. Each choice
changes the stage interface. U52 is an E-G9 test, so E6 cannot close green as written.

Fix (U52): replace the equality by two listed sets per CLI. (a) `dests - dropped == fields - added`, with dropped =
`{train, hash_view}` for `final_eval` and `{stages, sampler, rows_plan}` for `train_sft`, and added = `{compile,
plain}` for `train_rl` (minus any of them a phase adds as a flag) and `{lr}` for `g_probe`. (b) Each added field is
a new CLI flag on purpose; state it in E6 step 4 "Flag names that change at E6.4" ("`train_rl` gains `--compile`,
`--plain`; `g_probe` gains `--lr`, from the shared `TrainerFlags`"). CLEAN16 step 1 then adds `procs` to the
dropped set of `select_ckpt` and `g_probe`, as today.

### 3. major. After D7 step 6 no `fork_pool` caller is left, so E3b's R-2 branch archives `pool_size`, which the `--procs` defaults and E6.1 still call

Sources: E, CLEAN (D owns D7).

- D.md "### D7" step 6: "`final_eval`, `select_ckpt` and `g_probe` stop creating the pool for metrics"; the ceiling
  CLI "also stops using the pool in this commit"; "The pool code itself (`serve._Pool`, `fork_pool`, `pool_size`,
  `cpus`, `quota`) is deleted by CLEAN16 step 1".
- E.md "### E3", `verify/pool.py` row: "`_Pool`, `_POOLS`, `CGROUP`, `POOL_CAP`, `quota`, `cpus`, `pool_size`,
  `fork_pool` ... If no `fork_pool` caller is left on the E0 commit, archive these names instead of moving them
  (R-2)." E.3.3 says the same.
- E.md "### E6" step 1: `SelectStage.procs: int | None = None  # None -> serve.pool_size() in the CLI`;
  `ProbeStage.procs ... # None -> serve.pool_size() (g_probe.py:139)`.
- CLEAN.md "### CLEAN16" step 1: "This step is the single owner of the pool deletion".
- Code: the only live `fork_pool` callers are `final_eval.py:69`, `g_probe.py:157`, `select_ckpt.py:390` and
  `train_rl.py:50` [measured: grep of src and experiments]; CLEAN2 (wave 1) removes the `train_rl` one, D7 step 6
  the other three, both before E0. `select_ckpt.py:374` and `g_probe.py:139` keep `default=serve.pool_size()`.

Problem: on the E0 commit the R-2 condition holds by the doc's own order, so E3b archives `pool_size`. Then
`select_ckpt.main` and `g_probe.main` raise `AttributeError` at `ap.add_argument(..., default=serve.pool_size())`,
the preflight CPU dry run (`--procs 1` to both [code: tests/test_p2_scripts.py:369, 384]) fails, and preflight
refuses every rental. If the engineer reads the branch the other way (tests count as callers), E6 adds `procs`
fields that feed nothing: after D7 no CLI creates a pool, so `--procs` is dead code that E then types and ports
(owner rule: remove dead code). Three texts (E3b, E6.1, CLEAN16) disagree on who deletes what, and when.

Fix (one owner, earliest point): D7 step 6, in the same commit, deletes `--procs` from `select_ckpt` and `g_probe`
and `"--procs", "1"` from `tests/test_p2_scripts.py:369, 384` (nothing reads the flag after the switch). Then: E6.1
drops `SelectStage.procs`, `ProbeStage.procs` and `procs` from `FLAG_ORDER["select"]`/`["probe"]` and from U52's
`cli_parse(... "--procs", "1")` cases; E3b's R-2 branch fires (archive the pool names; U33's pool part is skipped,
as it already allows); CLEAN16 step 1 becomes a check ("`fork_pool`, `pool_size` and `--procs` are gone") and drops
its `SelectStage.procs` edit and its U52 note. If the owner keeps the pool until CLEAN16 instead, delete the R-2
branch from E3b and E.3.3 ("E3b always moves the pool names to `verify/pool.py`").

### 4. major. Front section 2 gives a target tree that E's move tables forbid (`infra/paths.py`, `rl/callbacks`, no `core/provenance.py`)

Sources: 00_front, E.

- 00_front.md "### 2.2 Layers": `| infra   paths, vast (offers, ledger, remote, chain), bundle, preflight, hf, |`;
  rl box: "build, trl_adapter ..., forward, rows, reward, config, run, callbacks, vllm_sync, vllm_requests"; policy
  box has no `callbacks`; `| core    spec (RunSpec, pydantic) |`.
- 00_front.md "### 2.3 Target tree": `core/ spec.py`; `rl/ ... run.py callbacks.py`; `infra/ paths.py vast/
  bundle.py preflight.py hf.py secrets.py`.
- E.md "### E8" step 1: "The one function is `objdel.core.provenance.repo_root()` (E2.5). It is not in `infra`,
  because `sft` and `rl` need it and must not import `infra`. There is no `infra/paths.py` and no re-export."
- E.md "### E5" table: "`policy/callbacks.py` | `callbacks.py`: `TraceStep` ... Not `rl/`: `train_sft.py:390` (->
  `sft/train.py`) and `g_probe.py:93` (-> `rl/build.py`) both use it". E.md "### E2" step 5 creates
  `core/provenance.py`.
- Code: `train_sft.py:390` builds `callbacks.TraceStep(args.profile, WARM)`; `g_probe.py:93` builds
  `callbacks.TraceStep(a.trace, wait=1)`.

Problem: front section 2 is headed "Target architecture and code organization" and says "the phases follow it".
An engineer who builds the tree from it puts `TraceStep` in `rl/callbacks.py`; `sft/train.py` then imports `rl`, and
the `rl | sft` independence contract (E0, E8.4) fails. An `infra/paths.py` re-export is the second root source that
E8 rules out, and a tree without `core/provenance.py` loses the code-hash fix of E2.5.

Fix (00_front 2.2 and 2.3): infra box "vast (offers, ledger, remote, chain), bundle, preflight, hf (+ licenses/),
secrets; repo root from `core.provenance`"; rl box without `callbacks`; policy box adds "callbacks (`TraceStep`)";
core box "spec (RunSpec, pydantic), provenance (repo_root, code_hashes)". Tree: `core/ spec.py provenance.py`,
`policy/ ... callbacks.py`, `rl/` without `callbacks.py`, `infra/ vast/ bundle.py preflight.py hf.py licenses/
secrets.py`, `data/` adds `audit.py`. Add one line under 2.3: "E's move tables win where this tree differs."

### 5. major. CLEAN8 step 1's token cap for `select_ckpt` has no field: after E4 `select_ckpt` always samples outside `build_trainer`, and `SelectStage` forbids extra keys

Sources: CLEAN, E.

- CLEAN.md "### CLEAN8" step 1: "Add a test-only cap through the trainer flags: the smoke test's stage sets
  `flags.grpo={"max_completion_length": 32}` ... If `select_ckpt` samples outside `build_trainer` on the post-E tree,
  the same commit gives its stage the same override and a test. Switch the `_tiny_vlm` default and the
  `test_cpu_smoke_train_then_select` parameters to `"qwen3-vl-2b"` only".
- E.md "### E4" step 1: the select bridge in `sft/select.py` runs the child `[sys.executable, "-I", "-m",
  "objdel.policy.sampler", "--kind", kind, ...]` with `params_for(model, ...)` = `policy.models.sampling(model,
  "vllm") | {...}`. E.md "### E6" step 1: `class SelectStage(StageBase, frozen=True, extra="forbid")` has no `grpo`,
  no `max_tokens` and no cap field; `--kind` and the child flags hold no cap either.
- Code: the dry run's docstring: "a random Qwen3-VL never stops, so its 1,536-token samples would take hours on a
  CPU" [code: tests/test_p2_scripts.py:331-333]; the dry run is preflight's `DRY_RUN` [code: preflight.py:35].

Problem: the "if" in CLEAN8 step 1 is always true on the post-E tree (E4), so the step needs a cap that reaches the
sampler child, and no stage field, child flag or `params_for` argument carries it. The engineer must design one (a
`SelectStage` field, a `--max-tokens` child flag, a `params_for` override) or drop `select` from the Qwen dry run;
the parameter switch "to `qwen3-vl-2b` only" does not say whether the Qwen case keeps `select=True`. With `select`
and no cap, the preflight dry run takes hours; without `select`, the dry run no longer covers `select_ckpt`, and
CLEAN16's edit of "the preflight dry run's `select_ckpt` call" names code that is gone.

Fix (CLEAN8 step 1; E6.1 and U52 follow): "The Qwen case keeps `select=True`. `SelectStage` gains `max_tokens: int |
None = None` (test-only; `None` keeps `MAX_COMPLETION`); `sft.select.sampler_cmd` passes `--max-tokens N` to the
child, and `params_for(model, ..., max_tokens=None)` overrides `max_tokens` when set. The smoke test passes
`--max-tokens 32` to `select_ckpt` and `flags.grpo={"max_completion_length": 32}` to `g_probe`. Test: the child argv
holds `--max-tokens 32`, and `params_for(..., max_tokens=32)["max_tokens"] == 32`." Add `max_tokens` to U52's added
set for `select_ckpt`.

### 6. major. DOC1 step 0's replace rule has no target for 35 scratchpad paths in the copied files, so check 14 cannot be made green by the rule

Sources: DOC.

- DOC.md "### DOC1" step 0: "in every copied file and in `Implementation.md`, replace each absolute scratchpad
  prefix and each short scratchpad prefix (`S[P]/impl/`, `S[P]/`) with its repo path from the table"; "Then run the
  check 14 scan on ... the copied files ... until it prints nothing; only then is check 14 green at the DOC4 commit."
  The table maps `impl/phases/*.spec.md`, `impl/notes/`, `impl/c_rev/`, `impl/r2rl/`, `impl/tests/`, `impl/reviews/*.md`,
  `impl/docfix/`, `impl/doc_mk/`, `impl/nv_catalog.json` and the deep findings.
- DOC.md "### DOC4" check 14 scans "`docs/implementation2.0/**/*.md`" (the copied reviews, notes, tests and specs).
- Measured: `grep -nE "S[P]/impl/(phases|r[4-7]data|r6skip|b_r4|nv_catalog_live)"` over the scratchpad `reviews/*.md`,
  `notes/*.md`, `tests/*.md`, `phases/*.spec.md` (minus `phases/*.spec.md` targets): 35 lines in 27 files, for
  example `reviews/A_phase_review.md:1` (`S[P]/impl/phases/A.md`), `reviews/final_r4_data.md:27` (`S[P]/impl/r4data/`),
  `reviews/final_r5_data.md:34,146` (`S[P]/impl/r5data/...`), `reviews/final_r6_data.md:7` (`S[P]/impl/r6data/`),
  `reviews/final_r6_qa.md:63` (`S[P]/impl/r6skip`), `reviews/final_r6_arch.md:176` (`S[P]/impl/b_r4/valhalf.py`),
  `reviews/final_r1_cuda.md:11` (`S[P]/impl/nv_catalog_live.json`). All check 14 patterns: 118 lines in 50 files.
  Each later review round adds more.

Problem: `impl/phases/<X>.md` (the phase sources, not copied: `Implementation.md` replaces them), `impl/r4data/`,
`r5data/`, `r6data/`, `r7data/`, `r6skip/`, `b_r4/` and `nv_catalog_live.json` have no row, so "its repo path from
the table" does not exist. The engineer must invent a mapping or delete text in review records. Round 7 arch-5 was
the same defect for the kernel folders.

Fix (DOC1 step 0 table and replace rule): add rows: `impl/phases/<X>.md` -> `docs/implementation2.0/Implementation.md`
(text form "phase <X> of Implementation.md"); `impl/{r4data,r5data,r6data,r7data,r6skip,b_r4}/` ->
`docs/implementation2.0/notes/review_probes/<folder>/` (copy only `*.py`, `*.json`, `*.txt` under 1 MB, as for notes);
`impl/nv_catalog_live.json` -> `docs/implementation2.0/notes/doc/nv_catalog_live.json`. Then: "any scratchpad path still
left after the table is replaced by `(scratch work file, not kept)` and listed in `docs/implementation2.0/README.md`".

### 7. minor. `to_argv(SelectStage)` still has no rule for the `--vision` JSON (prefixes and spacing); open since round 6

Sources: E. Round 6 arch-6 (minor, not fixed). Today `chain()` writes `json.dumps({"prepared": f"{data}/{v['prepared']}",
"reference": f"{rd}/{ref}"})` with default spacing [code: policy.py:437-440]. E6.1 types `vision: VisionRef | None`,
the SFT rule writes compact JSON, and "JSON that `to_argv` builds from fields ... is never formatted", so `{data}`
and `{rd}` inside the JSON are never filled. `vision_refs` (E6.8b) prefixes `/workspace/runs/<run_id>/`, so it reads a
relative reference. Golden 5 (r1_qwen3vl2b select stage, string for string) catches it, so no wrong result. Fix:
E6.1: "`VisionRef` holds paths relative to the data root (`prepared`) and the run folder (`reference`, when a
string); `to_argv(SelectStage)` adds `{data}/` and `{rd}/`, fills them, and writes `--vision` with default
`json.dumps` spacing (the chain form)."

### 8. minor. DOC1 step 0 does not copy evidence files the plan cites

Sources: DOC. The plan cites `docs/implementation2.0/reviews/r3arch_probe/field_order.py` (E6 "Flag order"
measurement), but the table copies `impl/reviews/*.md` only; `r4data/disagree.py` and `r4data/retest.py` (B, three
`[measured]`/`[estimate]` tags) have no row; `docs/implementation2.0/notes/nv_catalog.json` and
`docs/implementation2.0/notes/reviews/final_r2_{rl,data}.md`, `.../final_r2_data_check.py` name folders the table does
not create (it copies `nv_catalog.json` to `notes/doc/`, reviews to `reviews/`). Check 2 skips plan files, so no test
fails, but the evidence is lost when the scratchpad is cleared. Fix: copy `impl/reviews/**` (minus `docrev/`), add the
`r4data` row of item 6, and rewrite the three wrong citations to `notes/doc/nv_catalog.json` and `reviews/...`.

### 9. minor. DOC stale size and nav

Sources: DOC. DOC1 step 2 says `Implementation.md` is "1,480,096 bytes"; it is 1,608,894 now (still under DOC-G13's
2 MB, but it grows about 130 KB per round). DOC3's `nav` lists `implementation/README.md` and
`implementation/Implementation.md`, not `implementation/review-log.md`, which is not excluded, so check 13 reports it
at DOC3 step 3. Fix: write "about 1.6 MB, re-measured at the DOC1 commit" and add `implementation/review-log.md` to
the Plan section of `nav`.

### 10. minor. CLEAN16 step 1 deletes `--procs` from `final_eval`, which has none

Sources: CLEAN. "delete ... `--procs` from `select_ckpt`, `g_probe`, `final_eval`". `final_eval` has no `--procs`; it
calls `serve.fork_pool(serve.pool_size())` [code: final_eval.py:46-56, 69], and `EvalStage` has no `procs`. Fix:
"`--procs` from `select_ckpt` and `g_probe`" (or nothing, if item 3's fix moves the deletion to D7).

### 11. minor. E.3.2 shows the box chain with `--spec ... --stage`, against E6.4

Sources: E. E.3.2: "box: experiments/phase2/train_rl.py --spec ... --stage train" and "final_eval.py --spec ... --stage
eval_test". E6 step 4: "The box chain uses the flag form (`to_argv`), which golden 5 pins." Fix: E.3.2 shows
`train_rl.py <to_argv flags>` and notes that `--spec/--stage` is the Mac form.

### 12. minor. E6.4 calls `to_argv(stage)[2:]` with one argument

Sources: E. E6 step 4 (train_sft): "passes `to_argv(stage)[2:]` of an `SftStage`" and step 0 `to_argv(stage)[2:]`;
the signature is `to_argv(stage: Stage, rd: str, data: str)` (E6.1), and E6.7 uses `to_argv(stage, rd, data)[2:]`. Fix:
write `to_argv(stage, rd, data)[2:]` everywhere (the dry run passes its tmp folders).
