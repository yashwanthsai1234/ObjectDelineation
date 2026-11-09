# Final review round 16: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc reviewed: the built
`docs/implementation2.0/Implementation.md` (17,141 lines, 1,765,919 bytes, built 2026-11-09 03:49 [measured: `wc -l`,
`ls -l`]). Findings are written against the sources in `impl/phases/`. I read the round-15 sections of
`review-log.md` (00_front, B, D, E, CLEAN), `final_r15_arch.md` with its verdict, and the round-15 verdicts of cuda,
data and fresh-arch. Focus this round: the round-15 changes as they meet E, CLEAN and DOC (arch-1: the V4-failure
rule in E4 step 1, E.7, E owner decision 7, front 6.3, 99_back 3b; data-2: `docs/reference/val_halves.json` and the
widened keep regex in B1.0, B1.6, D1, CLEAN1; cuda-1: the new engine counters as they meet E4's `final_eval`
rewrite). No item repeats an item that a log marks as handled. The four round-15 arch minors (arch-2 to arch-5) and
the six round-14 arch minors are still open in the sources and are not filed again. I edited no repo file and ran no
git command that changes anything.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item rests on an NVIDIA product claim (this round's changes are a rollback rule, a keep regex and counters). No catalog lookup, no install. |
| engineering:code-review | Grepped every `__file__` root in `src/objdel` (only `secrets.py:14`) and every live reference to `sample_vllm` / `hf_sample` in `src/`, `experiments/`, `tests/` (items 1, 2). Read `final_eval.py`, `policy.py:370-405, 460-500, 527-536, 604-605`. |
| engineering:debug | Walked the E2.5 commit: which `parents[n]` roots exist then, which test (U34) runs then, and what `make test` says (item 1). |
| engineering:system-design | Order walk B1.0 -> B1.6 -> E2.5 -> E8.1 for the repo-root rule; E4 step 3 -> V3 (post-E path) -> DOC6 `make headline` for the V4 label (items 1, 3). |
| engineering:architecture | Checked the layer order of E.3.1 and front 2.2 against the round-15 change: `data.splits.val_half` -> `core.provenance.repo_root()` points down (`core` is the bottom layer), so the change is legal. |
| engineering:documentation | Found stale sampler names in module docstrings that no planned check reads (item 2) and a label that has no carrier into the results docs (item 3). |
| huggingface-skills:huggingface-best | No model choice in scope. Nothing changed. |
| system-design-skills:system-design | Same order walks, framed as "which artifact does each step break, and which planned $0 check sees it". |
| system-design-skills:service-decomposition | One package with an exhaustive import-linter `layers` contract (a modular monolith) is still the right form. No split proposed. |
| mattpocock-skills:codebase-design | Seam check of `core.provenance.repo_root()` as the one repo-root interface: it has three planned adapters before E8.1 (`secrets.py`, `splits.val_half`, `final_eval.frozen_entry`), but E2.5 converts none of them (item 1). |
| WebSearch / WebFetch | Loaded. Not needed for a finding: every item rests on repo text and code. |

## Checked and found correct (no item)

- Round-15 arch-1 (V4-failure rule) is now one rule in five places: E4 step 1 ("fix forward, no restore"), E.7 row
  ("Fix forward in `objdel.policy.sampler` ... no restore of the old samplers"), E owner decision 7, front 6.3 ("E4
  step 1's V4-failure rule applies"), 99_back 3b. It agrees with E-G14 and U2 ("with no V4 exception").
- Round-15 data-2 keep regex: B1.6 `^docs/reference/(frozen_views|val_halves)\.json$`, D1 step 0 (same, plus
  `^docs/implementation2.0/kernels/K[123]/(?!build/|results/)`), CLEAN1 step 1 note and the CLEAN sketch `kept_docs`
  (looser on kernels, so it is a superset and passes). Today `code_files()` drops `^(data|docs)/` with no keep
  [code: experiments/phase2/vast/policy.py:379-380].
- `val_half` path: `src/objdel/data/splits.py` exists, so `parents[3]` is the repo root on the Mac and
  `/workspace/src` on the box, where `sync_up` rsyncs the tracked files and installs editable [code: policy.py:398-399].
  After E2.5, `core.provenance.repo_root()` finds `pyproject.toml` + `src/objdel` there. `data -> core` is a downward
  import in E.3.1.
- ADR numbers: front 4 table 0012-0022; B 0013/0014/0015, C 0022, D 0016, E 0017-0020, CLEAN 0021, DOC 0012; disk
  holds 0001-0011 [measured: ls docs/architecture/decisions]. No collision.
- V4 launch lines: `sync_up` defaults to the `lidc` bundle [code: policy.py:604]; `fetch_hf` takes a revision and
  writes `/workspace/runs/<run>/ckpt` [code: policy.py:527-536]; `job()` runs a non-`bash` argv as
  `/opt/venv/bin/python -I <argv>` [code: policy.py:486-489], so E.8's `("-m", "pytest", ...)` form is right;
  `experiments/cuda_checker/run_box.sh` exists.
- Round-15 cuda-1 counters and E4: D7 step 5's `final_eval` records (`reset_stats()` per view, `eval_scratch_peak`,
  `table_bytes`, `live_device_bytes`) are not named in E4's `evaluate` rewrite, but S12 (Mac, fakes, `final_eval.json
  ["engine"]` keys and one `reset_stats()` per view) fails at $0 if E4 drops them. No item.
- `make cov` gate 90 [code: Makefile:26; .github/workflows/tests.yml:44]; `testpaths = ["tests"]`
  [code: pyproject.toml:67], so `archive/phase2_samplers/tests/` is not collected.
- No scratch or home paths in the built doc or the log [measured: grep count 0 for both], so DOC1 step 1 (a) and DOC4
  check 14 stay green on the doc files.

## Items

### 1. minor: E2.5's AST test (U34) cannot pass on the E2.5 commit: two `__file__` roots that E2.5 does not convert still sit under `src/objdel/`.

Sources: E (E2 step 5, E8.1 "Repo root, one function", U34), B (B1.0).

- E.md "### E2" step 5: "Tests (U34): ... AST: 0 `__file__`-based roots (`parents[`, `.parent /`, `with_name`) under
  all of `src/objdel/`, except the search in `core/provenance.py`." The step converts only `CODE` and `code_hashes`
  ("A move commit brings `CODE` and `code_hashes`; a change commit adds `repo_root` and points `code_hashes` at it").
- E.md "### E8" step 1: "`secrets.py` sets `ROOT = Path(__file__).resolve().parents[2]` [code: secrets.py:14] ...
  Replace every `__file__`-based root in moved files ... the AST test of E2.5 (0 `__file__`-based roots under all of
  `src/objdel/`) stays green." But `src/objdel/secrets.py:14` is `ROOT = Path(__file__).resolve().parents[2]`
  [code], the only `__file__` use in `src/objdel` today [measured: grep], and it moves only in E8.1. So the test is
  red from E2.5 to E8.1, not "green".
- B.md B1.0 (round 15, data-2): `val_half` "reads `Path(__file__).resolve().parents[3] / "docs/reference/
  val_halves.json"` ... From E2.5 on it reads `core.provenance.repo_root() / ...`". E2 step 5 does not list this
  change, so the E2.5 engineer has no instruction to make it.
- Caught-before-spend: `make test` on the E2.5 commit fails and names both files. So minor.

Fix (E.md E2 step 5): add after "points `code_hashes` at it": "The same change commit replaces the two other
`__file__` roots under `src/objdel/` with `repo_root()`: `objdel/secrets.py:14` (`ROOT = repo_root()`; `secrets`
stays in `exhaustive_ignores` until E8.1 moves it) and B1.0's `data/splits.py` lookup of `docs/reference/
val_halves.json`. Then U34's AST check is green from this commit on." In E8.1, change "Replace every `__file__`-based
root in moved files (also `preflight.py:26`)" to "... (`preflight.py:26`, `policy.py:46-47`; `secrets.py` already uses
`repo_root()` since E2.5)".

### 2. minor: after the E4 step 3 archive, four module docstrings still name the archived samplers, and no planned check reads them.

Sources: E (E4 step 1, E4 step 3, U2, E-G15), DOC (DOC4 check 2).

- Live docstrings that name the old scripts [code]: `src/objdel/train/sample.py:18` ("the JSONL bridge to
  experiments/phase2/sample_vllm.py"), `src/objdel/train/sample.py:138` ("sample_vllm's completions"),
  `src/objdel/eval/metrics.py:62` ("train.sample / sample_vllm"), `experiments/phase2/select_ckpt.py:10,13`
  ("`sample_vllm.py` under /opt/vllm, or `hf_sample`" and "hf_sample keeps whole finished batches"). E moves these
  texts unchanged in pure move commits (R-1) to `policy/sampler.py`, `verify/metrics.py` and `sft/select.py`.
- E.md E4 step 1 lists only `tests/` hits ("[measured: grep of `tests/`, 2026-11-09]") and step 3 rewrites only
  `parity.py` and `r0_parity.py`. U2 checks that the files exist only under `archive/`; E-G15 resolves `objdel.*`
  strings and `tests/...py` paths; DOC4 check 2 reads `doc_files` only. None reads `src/` docstrings for
  `sample_vllm|hf_sample`. Owner rule: old code archived "so nobody is confused".
- No result changes. So minor.

Fix (E.md E4 step 3): add "The same commit rewrites the docstring lines that name the old scripts to name
`objdel.policy.sampler`: `policy/sampler.py` (from `sample.py:18,138`), `verify/metrics.py` (from `metrics.py:62`) and
`sft/select.py` (from `select_ckpt.py:10,13`)". Extend U2: "`grep -rn 'sample_vllm\|hf_sample' src experiments/phase2/*.py`
finds 0 lines outside `archive/` and frozen records".

### 3. minor: the V4-failure label "sampler not box-verified (E-G10 open)" has no carrier into the results docs.

Sources: E (E4 step 1), DOC (DOC6), 00_front (6.3 V4 paragraph).

- E.md E4 step 1 (round 15): "Until V4 passes, every post-E eval result carries the label "sampler not box-verified
  (E-G10 open)"." This sentence is the only place the label exists [measured: grep of all sources]. Front 6.3 and
  99_back 3b point to the rule but do not name the label.
- DOC.md DOC6: `headline.json` keys are `run, model_ref, ckpt, seed, K, ci_unit, test, val` (I03 asserts "keys
  exactly"), so `make headline` on a post-E V3 test eval writes README and HF-card numbers with no place for the
  label, and I03 forbids adding one.
- Gate B3 is safe (both arms use `VllmSampler`, B3.9 step 2). Only the published caveat is lost. So minor.

Fix: DOC.md DOC6 schema: add an optional key `"sampler": "trainer" | "policy.sampler" | "policy.sampler (not
box-verified, E-G10 open)"`, written by `headline.py` from a new `final_eval.json` key `sampler` that E4 step 3
writes (`"policy.sampler"`); I03's "keys exactly" list adds it; the README test renders the caveat when the value
ends in "open)". Front 6.3 V4 paragraph: name the label once and point to DOC6.

## Counts

blocking 0, major 0, minor 3 (items 1-3).
