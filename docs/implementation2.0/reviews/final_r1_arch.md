# Final review, round 1: architecture, code organization, layering, APIs, DOC, CLEAN

Reviewer role: principal software architect for production ML systems.
Document: `REPO/docs/implementation2.0/Implementation.md` (11,117 lines, round 1).
Earlier `final_r*_*.md` files: none existed, so no items are repeated from them.
Clean context kept: no transcripts, no memory folder, no `REPO/.claude/`, no git write, no rental, no REPO edit.

## Skills used

| Skill | What it changed in this review |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | No NVIDIA catalog skill covers package layering or ADR hygiene. It shaped item 22 only: the K3 torch extension needs a build home, and the doc's own pointer (`tilegym-monkey-patch-kernels-to-transformers`) is a recommendation, not an install. |
| `engineering:code-review` | Each claim checked against code at the working tree (be52b1c + untracked docs): `vision.py:69`, `callbacks.py:63`, `policy.py:69,75`, `g_probe.py:86`, `configs.py:189`, `signal.py:31`, `fast.py:561-565`, `serve.py:52,172`, `train_rl.py:35,50`. All hold. Found item 5 (lost `gc.freeze`) this way. |
| `engineering:debug` | Traced what breaks at run time, not at import rewrite: the A6 import edge into the reward package (item 2), the CLEAN11 archive of a module that D3 keeps live (item 3), pytest import mode after the test move (item 16). |
| `engineering:system-design` | Compared the front matter target architecture (sections 2.2-2.5) with the phase designs box by box (item 9). |
| `engineering:architecture` | Checked every ADR number the phases allocate (item 4) and every tech-stack row for an owning task (item 10). |
| `engineering:documentation` | "Link, do not duplicate": found the scratchpad dependency of a tracked doc (item 14) and duplicate ownership of the same edit (items 15, 17). |
| `huggingface-skills:huggingface-best` | No model query: the model choice is not in this review's scope. Checked only that a new model row needs no code outside `policy/models.py` (E.0 claim): it holds. |
| `system-design-skills:system-design` | Treated the phase order as the main hypothesis. Found that the order is stated three different ways (item 1). |
| `system-design-skills:service-decomposition` | Modular monolith is the right default (one package, enforced layers, verifier HTTP face optional). No item asks for a service split. Checked that the in-process `verify.api` and the HTTP face share one scoring path (E7.3): it does. |
| `mattpocock-skills:codebase-design` | Checked each seam for "two adapters = real seam". The polygonizer seam has three conflicting designs (item 3). The sampler interface differs between the front matter and E4 (item 9). |
| Web (WebFetch, WebSearch) | import-linter layers and forbidden pages: `exhaustive_ignores`, `(layer)` optional, `|` independent, and "Indirect imports will also be checked" (forbidden contract) [web: import-linter.readthedocs.io/en/stable/contract_types/{layers,forbidden}/]. The docs do not say if function-level imports are scanned, so item 2 asks for a negative test. |

## Counts

blocking 5, major 9, minor 8 (22 items).

## Items

### 1. blocking | 00_front, E (also B, C, D) | The phase order is stated three different ways

Problem. Front section 1, order rule 2 and the flow diagram say: E moves files after B1, C and D land
(`Implementation.md:109-110`, diagram lines 88-100). C.8.3 says "land C first, then E moves the new modules unchanged"
(line 4637). D.8 and D.7 R10 say E's moves come after D3 and D4b (lines 6052, 6088). But E says the opposite: "E5
before C coding starts, so C edits the split files" (E.7, line 7428), "B2 and B3 land after E3/E5: v5 in
`verify/reward/`" (line 7463), "C edits the split files after E5" (line 7465), and "Phase D ... turns
`geom/polygonize.py` into an oracle after E3" (line 7468). Every A-D task names pre-E paths (`src/objdel/train/...`,
`checker/...`, `reward/...`). An engineer cannot tell which tree a C or D task edits.

Fix. Adopt the front order (E after A-D, CLEAN wave 1 before E), because all A-D task texts already use pre-E paths.
Edit E: (a) E.7 row 3, mitigation -> "E runs after the last A-D commit that edits `fast.py`, `dataset.py`, `reward/`,
`checker/`; B1 runs on the tagged pre-E tree"; (b) E.8 dependency bullets for B, C, D -> "B2/B3, C and D land on
pre-E paths; E moves their modules (see item 6)"; (c) E.3.3 and E.5 -> the seam is D3's `rings.polygonize`, which
E3a moves (item 3). Add one line to front section 1: "Every A-D task path is pre-E; E's move tables are the only
place where post-E paths are defined."

### 2. blocking | A, E | A6's NVTX module breaks both import-linter contracts that E0 adds

Problem. A6 creates `src/objdel/train/nvtx.py` and adds a `checker` range in `src/objdel/reward/score.py`
`score_many` (lines 869-874, 905-906; test sketch line 2065 `from objdel.train import nvtx`). So the reward package
imports `objdel.train` (today the graph is `reward -> {_gpucheck, checker, geometry}`, E.1 item 10). A lands first, so
E0's interim contract `layers = ["train", "eval", "reward | targets", ...]` (line 6771) is broken on day one: a lower
layer imports a higher one. `nvtx.py` imports torch inside `enabled()` (line 893). The forbidden contract checks
indirect imports [web: import-linter forbidden page], so `pure-core` (source `objdel.reward`) also sees
`reward -> train.nvtx -> torch` if grimp counts the function-level import. After E, `core/marks.py` (line 6650)
holds the same lazy torch import, and E8.4 adds `objdel.core` to `pure-core` sources (line 7325). E's G4 checks only
`sys.modules` after import, so it passes while lint fails. C0 adds `count_syncs()` to the same file (line 4051).

Fix. (a) Move the `checker` range out of the verifier side: open it in `train/rewards.py` `TRLReward.__call__` around
the `score_batch` call (same time span, one call per step), and drop the `reward/score.py` edit from A6. (b) Keep the
NVTX module in the ML layers: E5 move row `train/nvtx.py -> policy/marks.py` (not `core/marks.py`); delete `marks.py`
from the `core` box in E.3.1, front 2.2 and 2.3, and from section 4 row "Profiling". (c) Add E0 negative test (d): a
planted function-level `import torch` inside a function in `objdel/data/__init__.py` must break `pure-core`; if it
does not, record that grimp skips function bodies and keep G4 as the guard.

### 3. blocking | D, E, CLEAN | Three conflicting designs for the polygonizer seam; CLEAN11 archives live code

Problem. D3 creates `checker/rings.py` with `trace` and `polygonize(mask, offset)` that always call
`_gpucheck.check_masks` (also on the CPU build), moves `checker/polygonize.py` and `checker/bitquad.py` to
`tests/oracles/`, and keeps `checker/certify.py` as a device-backed module that `data/riga.py:226` calls through
`certify_many` (lines 5715-5741). E3a moves `checker/polygonize.py` to `geom/polygonize.py` and adds a seam in
`geom/__init__.py` that calls the numpy code; E.3.3 says D later uses the engine only "when
`engine.build() == "cuda"`" (lines 6693-6696, 6926-6928). That keeps numpy on every CPU build, against the owner rule
"the polygonizer must be CUDA" and against D3. CLEAN11 (lines 9762-9786) moves `polygonize.py` a third time, and
archives `certify.py` and `bitquad.py` "if phase D's device flag word covers all three methods". After D3, `certify`
is live (RIGA uses `certify_many`), so CLEAN11 step 2 breaks RIGA. CLEAN11 also names pre-E paths and
`experiments/cuda_checker/compare.py`, which E3b moves to `tests/oracles/cuda_parity.py`. Its gate
(`grep -rnw polygonize src/objdel` = 0 import lines) fails on any caller of the E3a seam.

Fix. One design, D3's. (a) E3a move table: `checker/rings.py -> geom/rings.py`, `checker/certify.py ->
geom/certify.py`, `checker/{topology,volume}.py -> geom/`; delete the rows for `polygonize.py` and `bitquad.py` (already
in `tests/oracles/` after D3); the seam is `geom.polygonize = rings.polygonize`, re-exported from `geom/__init__.py`.
(b) Delete E.3.3 paragraph 2 and the E3a "Change commit" paragraph; replace E.5 K1 bullet with "Production reaches K1
through `geom.rings.trace` / `geom.polygonize` (D3)". (c) Replace CLEAN11 with a check-only task: "after D3, assert
`tests/oracles/{polygonize,bitquad}_oracle.py` exist, `src/` has no numpy ring tracer, and `certify` uses
`rings.trace`; nothing is archived". (d) Front 2.3 tree: `geom/ engine.py text.py shape.py rings.py topology.py
volume.py certify.py`.

### 4. blocking | B, D, E, DOC, CLEAN | ADR numbers collide

Problem. ADR 0012 is allocated four times: B (number match, line 2456), D
(`0012-cuda-polygonizer-in-production.md`, line 6000), E (`0012-layers-and-import-linter.md`, line 7329) and DOC10
(`0012-docs-site-with-zensical.md`, line 8611). 0013 is B (reward v5, line 2461) and E (run spec, line 7330). 0014 is B
(scheduled rows, line 2450) and E (TRL adapter, line 7331). CLEAN14 takes "next free number". The last ADR on disk is
0011 [measured: `ls docs/architecture/decisions`]. DOC4 check 4 (ADR index and headings) fails on duplicates, and
cross-references such as "record it in ADR-0013" are ambiguous.

Fix. Put one allocation table in front section 4 and use it everywhere: 0012 docs site (DOC10, lands first), 0013
trainer-vLLM number match (B2), 0014 reward v5 (B3), 0015 scheduled RL rows and the numpy-polygonizer exception (B1),
0016 CUDA polygonizer in production (D9), 0017 layers and import-linter + wheel (E), 0018 typed run spec + CliApp +
vllm_native image (E), 0019 TRL adapter and start-up check (E), 0020 verifier API versioning + sampler on vllm_native
(E), 0021 library keep list (CLEAN14). Replace every `ADR-00xx` / `ADR 00xx` in B, D, E, DOC, CLEAN and 99_back items 7
and 19 with the new number. Add a DOC4 check: ADR file numbers are unique.

### 5. major | E, CLEAN | E1.3 and CLEAN2 give different code for the same change

Problem. E1.3 removes `--procs`, `serve.fork_pool` and `pool.shutdown()` from `train_rl.py` (line 6893). CLEAN2
replaces the `fork_pool` line with `gc.freeze()`, because `fork_pool` freezes the GC after the preload
(lines 9494-9500). `fork_pool` does call `gc.freeze()` [code: src/objdel/train/serve.py:172]. E.8 says "E1.3 is the
same item as CLEAN item 2, so do it once" (line 7470) but does not say which text wins. If E1.3 lands, the freeze is
lost and step time can change for a reason nobody measured.

Fix. Delete E1.3's steps and write: "E1.3 = CLEAN2 (CLEAN wave 1 runs before E). E1 only checks that `--procs` is
gone from `train_rl` and that `gc.freeze()` is called once after the last `serve.preload`."

### 6. major | E (also A, C, D, CLEAN) | E's move tables have no home for modules that A-D and CLEAN wave 1 create

Problem. E's tables are built from be52b1c. Modules created before E (front order) are missing: `eval/ceiling.py`,
`eval/select.py` (A, lines 449, 787), `train/nvtx.py` (A6), `train/gpu_pixels.py` (C1, line 4098),
`train/vllm_sync.py` (C5, line 4338), `train/vllm_requests.py` (C6, line 4379), `eval/device.py` (D7, line 5894),
`checker/rings.py` (D3), `geometry/box.py` (CLEAN9, line 9717). G14 requires `src/objdel/{checker,geometry,eval,train}`
to be gone (line 7416). G16 records old names "on the base commit" (line 7000), so these modules escape G16 and fail
only at G14, with no destination given.

Fix. Add rows to E3a/E3b/E5: `eval/ceiling.py -> verify/ceiling.py`, `eval/select.py -> verify/select.py`,
`eval/device.py -> verify/device.py`, `geometry/box.py -> geom/box.py`, `checker/rings.py -> geom/rings.py`,
`train/gpu_pixels.py -> policy/gpu_pixels.py`, `train/vllm_sync.py -> rl/vllm_sync.py`,
`train/vllm_requests.py -> rl/vllm_requests.py`, `train/nvtx.py -> policy/marks.py` (item 2). Change G16 and
`old_names.json`: record on the E0 commit (after A-D and CLEAN wave 1), not on be52b1c. Add a line: "a module created
after E0 and before E5 must add its row to this table in the same commit".

### 7. major | B, E, 00_front | `eval/coverage.py` has two destinations

Problem. B1.3 moves `eval/coverage.py` to `data/coverage.py` and extends it (lines 2299, 2570-2571); front 2.2 and
2.3 agree (`data ... coverage`). E3b moves `eval/coverage.py` to `verify/coverage.py` (line 6946) and the E.3.1 box
lists `coverage` under verify. B1 runs first, so at E3b the source file no longer exists, and the E3.1 diagram is
wrong.

Fix. Delete the `verify/coverage.py` row in E3b and `coverage` from the E.3.1 verify box; add `coverage, windows,
quality` to the E.3.1 data box; add `data/coverage.py` (B1.3) to the list of files E3c leaves in place.

### 8. minor | B, E, 00_front | RLOO advantages: verify or rl?

Problem. B3.4 puts `advantages_v5` in `src/objdel/reward/advantages.py` (line 2865), which E moves to
`verify/reward/advantages.py`; front 2.2 lists `verify reward (... advantages)`. E3b's note and E.8 say "the RLOO
baseline in `rl/`" (lines 6942, 7464). The function is numpy and fingerprinted with `REWARD_V5`, so verify is a
valid home.

Fix. In E3b's note and E.8, write: "advantages_v5 (B3.4, numpy, fingerprinted) stays in `verify/reward/advantages.py`;
`rl/reward.py` only calls it."

### 9. major | 00_front | Sections 2.2-2.5 do not match the phase designs

Problem. Front 2.2 says "Source: phase E section 3.1", but differs from it: verify lacks `store`, `pool`, `http`;
data lacks `probes`, `prepare`; data lists `viewbank` and front 2.3 lists `data/viewbank.py` and `data/hashing.py`,
which no task in A-CLEAN creates (`grep viewbank` hits only lines 140, 170, 296); 2.3 lists `infra/ ... nvtx (marks)`
beside `core/marks.py`, and lists `geom/polygonize.py certify.py` after CLEAN11 removes them. Front 2.4 draws the view
bank and "GPU patchify (K3a)" as the target path; C.5 makes K3a fallback number 3 and C.8.3 makes the view bank
optional. Front 2.5 routes RL generation through the "sampler seam" with "32 x n=8 requests" and "fused log-prob K3b";
E4 makes the sampler eval-only (training stays on TRL's `VLLMGeneration`), C6 is "default off", and C.5 says K3b is
not in Phase C. The seam table says `sample(rows, k) -> texts`; E4 defines `sample(rows, params) -> list[Sample]`.

Fix. Regenerate 2.2 and 2.3 from E.3.1 plus items 2, 3, 6, 7 (no `viewbank.py`, no `hashing.py` unless a task is
added; no `infra/nvtx`). In 2.4, label the view bank "optional, no task in A-E (T1)" and replace "GPU patchify (K3a)"
with "GPU patchify (`gpu_pixels.patchify`, torch ops; K3a fallback)". In 2.5, remove the sampler seam column (training
generation is TRL `VLLMGeneration`), write "32 prompts x 8 answers (n=8 requests only if C6 passes)", and replace
"fused log-prob K3b" with "log-prob (torch; K3b not planned)". Seam table row: `sample(rows, params) -> list[Sample]`.

### 10. major | 00_front | Section 4 stack rows with no owning task, and three rows that contradict the phases

Problem. Rows with no task in any phase: `uv.lock` + `uv sync --locked`, pyright, `ruff format`, pre-commit with the
`no-ai-attribution` commit-msg hook, `events.jsonl`, Trackio `report_to`, view bank, `data_revision` in the run spec
(E6's `RunSpec` has no such field, line 7229), checkpoint resume every K steps [measured: grep, 1-4 hits each, all in
section 4 or 99_back]. Rows that contradict a phase: "REPLACE scipy EDT with cuCIM/CuPy" versus D.2 "cuCIM ... not
adopted for production" (line 5523) and K2 EDT; "Move transport to `vastai-sdk`" versus CLEAN13 ("deprecated in favour
of `vastai`; keep the subprocess CLI", lines 9214, 9851); "ADD pytest-xdist" versus D8 "pytest-xdist is not a
dependency" (line 5978).

Fix. Add a column "Owning task" to section 4. Fill it, or write "deferred, not in A-E" for each row with no task.
Change three rows: EDT -> "K2 separable EDT on the device (D7); cuCIM not adopted (D.2)"; GPU orchestration -> "keep the
`vastai` CLI 1.11.0 as a `uv tool` (CLEAN13); no SDK"; Tests -> "pytest-xdist deferred (D8 runs without `-n`)". Either
add `data_revision: str | None` to `RunSpec` in E6 with a test, or drop it from the row.

### 11. blocking | E, 99_back | The default for E8.2 breaks the repo rule "experiments/<phase>/"

Problem. CLAUDE.md and BRIEF section 3 bind "experiment scripts in `experiments/<phase>/`, run as modules". E8.2
moves live CLIs to `experiments/{data,sft,rl,verify,infra}/` with the default "[yes]" (line 7478; 99_back item 15).
DOC13 only updates the run-command example in CLAUDE.md (line 8673), not the rule.

Fix. Change the default to "[no]: live CLIs stay in `experiments/phase2/`; finished one-offs go to `archive/`". Keep
the layer grouping as the alternative, with the precondition "owner edits the CLAUDE.md Project rule first". Change
golden 5's rename map, G15 and `test_cli_smoke` text to work with either choice (module names from one constant).

### 12. major | E, CLEAN | Nobody owns the archive of the 12 old YAML plans; wrong CLEAN cross-references

Problem. E6.6 says "CLEAN item 11 does the move" (line 7251); CLEAN11 is the polygonizer task. CLEAN8 moves only the
LFM configs "in the same commit as phase E6" (line 9707). No CLEAN task archives the other plans. E8.2 says finished
one-off scripts go to archive by "CLEAN items 5, 7, 12" (line 7301); CLEAN12 is the hygiene test.

Fix. E6.6: "E6 moves the 12 old plans to `archive/phase2_runs/configs/` with a README (Was / Replaced by / Commit) in
the commit that switches `chain()` to specs; CLEAN8's LFM configs go in the same commit." E8.2: "CLEAN items 5 and 7".
Add the 12 file names to E6.6 (from `ls experiments/phase2/configs`, 16 files [measured] minus the 4 ported).

### 13. major | DOC | DOC1 leaves `docs/implementation2.0/` untracked; its counts are stale

Problem. DOC1's `git add` list (lines 8256-8257) has no `docs/implementation2.0`. That folder holds this plan and the
kernel drafts: 348 files [measured: `find docs/implementation2.0 -type f | wc -l`]. Files outside `docs/research/` are
now 385, not 37 [measured]. DOC1 step 3 (tracked count equals on-disk count) fails. DOC3's strict-build trial
(2.1-2.7 s, "No issues found") ran before this folder existed, and check 13 (nav complete) will report it.

Fix. Add `docs/implementation2.0` to the DOC1 `git add` line, after the pre-publication scan covers it (it holds
absolute local paths, item 14). Re-measure the counts in DOC.1 item 1 and DOC.3.3. In DOC3, add
`implementation/kernels/**` to `exclude_docs` (source files, not pages) and add `implementation/Implementation.md` to
the Plan section of `nav`; re-run `zensical build --strict` and record the result.

### 14. major | 00_front, B, C, D, E, 99_back | A tracked doc depends on an ephemeral scratchpad

Problem. The doc cites `(planning probe, not kept: ...)` 112 times and holds 10 absolute paths under `/[p]rivate/tmp/...` and
`/[U]sers/...` (lines 431, 2970-2977, 3824-3825, 5430-5431, 6477-6478) [measured: grep -c]. Some are inputs, not only
citations: E5.2 "write it as specified in docs/implementation2.0/notes/A1.md R7" (line 7010); B1.3 starts `quality.py` from
`docs/implementation2.0/seeds/b1_scripts/` (line 2580); E2 test contracts `docs/implementation2.0/specs/E.spec.md` (line 6479). The scratchpad is
session storage and will go. The repo is public (DOC1 step 1), so the absolute paths also expose the local user path.

Fix. Copy every cited note, spec and script into `docs/implementation2.0/notes/`, `docs/implementation2.0/specs/` and
`docs/implementation2.0/seeds/` (as was done for the kernels), replace `docs/implementation2.0/` with `docs/implementation2.0/`, and delete
the REPO/SP header lines. Add a DOC4 check: 0 matches of `/[p]rivate/tmp|/[U]sers/|docs/implementation2.0` in tracked Markdown.

### 15. minor | A, E, DOC | The "TRL 1.16" text fix has three owners and two texts

Problem. A7 test S06 (line 1742), E1.4 (line 6896, "-> 1.17.1") and DOC2 (line 8282, "-> TRL (version pinned in
pyproject.toml)") edit the same docstrings.

Fix. DOC2 owns it (DOC runs first). Delete E1.4; change A7 S06 to assert that no live file contains "TRL 1.16"
(text-agnostic).

### 16. major | E | Test file locations disagree, and the test move has no import-mode rule

Problem. E.4 creates new tests in `tests/core/`, `tests/verify/`, `tests/experiments/` from E0 on (lines 6751, 6883,
6892). E.9 says new files "sit in `tests/phase_e/` and move in the E8 test move", and its sketches use
`tests/phase_e/` with `parents[2]` (lines 7548, 7688-7695). `tests/` has no `__init__.py` [measured], so pytest's
default prepend import mode needs unique basenames across all folders. E8.3 moves 63 files into 8 folders and adds
new ones (`test_api.py`, `test_paths.py`, `test_spec.py`, ...) with no rule. Helper imports `from oracles import ...`
(13 uses [measured]: grep of tests/*.py) work only because `tests/conftest.py` puts `tests/` on `sys.path`.

Fix. Pick `tests/<layer>/` from E0 (delete the `tests/phase_e/` sentence in E.9 and fix the sketch paths to
`parents[2]` from `tests/<layer>/`). Add to E8.3: "test basenames stay unique across `tests/**`; `tests/infra/
test_unique_basenames.py` asserts it; `tests/conftest.py` stays at the root so `from oracles import` resolves".

### 17. minor | A, E, CLEAN, D | Duplicate files and duplicate ownership of the same check or deletion

Problem. A creates `tests/test_repo_hygiene.py` (S01-S11, line 1291) and CLEAN12 also creates "new
`tests/test_repo_hygiene.py`" (line 9791). E0 `test_src_never_imports_tests` and CLEAN12
`test_src_never_imports_experiments_or_tests` test the same edge. E8.1's `.claude/`/`archive/` drop repeats CLEAN1.
E1.5 repeats CLEAN9 step 3 (Makefile default). The CPU pool deletion has three owners: "Phase D deletes `pool.py`"
(E.3.3, E.8), D5, and CLEAN16 step 1.

Fix. A owns `test_repo_hygiene.py`; CLEAN12 "extends" it. Delete E0 `test_src_never_imports_tests` (CLEAN12's test
runs first; E8 only moves it to `tests/infra/`). Delete E8.1's drop-pattern bullet and E1.5 (point to CLEAN1, CLEAN9.3).
CLEAN16 step 1 owns the pool deletion; E and D say "after D5, CLEAN16 deletes `verify/pool.py`".

### 18. minor | C, E | C names E homes that E does not define

Problem. C.0 and C.8.3 say C's new modules land at `rl/caches.py` and `rl/sync.py` (lines 3853, 4637). E puts the
caches in `policy/caches.py` and has no `rl/sync.py`.

Fix. In C.0 and C.8.3 write `policy/caches.py`, `policy/gpu_pixels.py`, `rl/vllm_sync.py`, `rl/vllm_requests.py`
(same names as item 6).

### 19. minor | CLEAN, E | CLEAN1's bundle test hard-codes module names that E8 renames

Problem. `test_bundle_imports_standalone` imports `experiments.phase1.reward_gates`, `tests.oracles.reward_harness`
and four `experiments.phase2.*` CLIs by literal name (line 9478). E8 changes `IMPORT_CHECK` and moves k3_harness to
`tests/harness/` and may move the CLIs. E8 does not list this test.

Fix. Make the test read its module list from the `IMPORT_CHECK` constant plus a `LIVE_CLIS` constant in the vast
policy module; add `tests/test_vast.py::test_bundle_imports_standalone` to E8.2's list of strings to rewrite.

### 20. minor | E | `dataset.py` name count disagrees inside E

Problem. E.0 revision row says 31 top-level names (line 6508); E5 and E.R say 29 (lines 6996, 7494). An AST count
gives 29 [measured: ast over `src/objdel/train/dataset.py`].

Fix. Change 31 to 29 in the E.0 revision row.

### 21. minor | 99_back | Owner decision labels point at the wrong E tasks

Problem. Items 14-19 are labelled "E1:" to "E6:", but they are E.8 owner decisions 1-6 (E1 is "small debts", E3 is
the lower layers). Item 19 cites "ADR 0012" (collides, item 4). Item 11 (archive 51,400 PNGs after the view bank)
depends on a view bank that no task builds (item 9).

Fix. Relabel items 14-19 "E decision 1" to "E decision 6"; use the new ADR number from item 4; mark item 11 "only if
a view-bank task is added".

### 22. minor | 00_front, C, E | No build home for the K3 torch extension

Problem. Front 2.3 says K1-K3 merge into `cuda/`. E.5 says K3 is a torch extension that belongs in `policy/` and
`rl/`. `cuda/` today builds one nanobind extension, `objdel._gpucheck`, with scikit-build-core [code: pyproject.toml
`wheel.packages`, `CMakeLists.txt`]. No task names the K3 build target, its extra, or an import ban.

Fix. Add to C1 step 7 (the only K3 use, fallback number 3): "if K3a is adopted, its source goes to `cuda/k3/`, built
as a separate torch extension `objdel._k3` only with the `train` extra; only `policy/gpu_pixels.py` imports it; add a
TID251 ban `objdel._k3` with that one per-file exception". In front 2.3, write "K1, K2 merged into `gpucheck.cpp`
(D); K3 only if adopted (C1)".
