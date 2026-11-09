# Final review round 13: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc reviewed: the built
`docs/implementation2.0/Implementation.md` (16,878 lines [measured: `wc -l`, 2026-11-09]). Findings are written against
the sources in `impl/phases/`. I read `review-log.md` (round 12 for every phase), `final_r12_arch.md`, and the
round-12 changes in B1.4a, B1.6, B1.7 (exclusion list), CLEAN1 and D1 as they meet E, DOC and CLEAN. No item repeats
an item that a log marks as handled. The round-11 and round-12 arch minors that are still open (for example r12
items 3-6: `RlStage` has no `--merge-last` form, stale "D3 puts the polygonizer in production" text, the DOC9 owner
row for B1.4a, `verify_target_hashes` missing from `LIVE_CLIS`) are not filed again. I edited no repo file and ran
no git command that changes anything (only `git ls-files`, `git check-ignore`). Probe output:
`reviews/r13arch_probe/riga_ignore.txt`.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item depends on an NVIDIA product claim. No catalog lookup, no install. |
| engineering:code-review | I read CLEAN5's move list and its owner `git add -A` line against the real `.gitignore` and the ignored files on disk (item 1). I read E2's `final_eval.py` name table against the names A3 and B1.6 add to that script (item 3). |
| engineering:debug | Reproduce-then-isolate on item 1: listed the ignored files in the trees that CLEAN5 moves, then ran `git check-ignore` at the old and the new paths (0 at the old, 1 at the new). |
| engineering:system-design | Walked the order DOC1 -> CLEAN wave 1 (in the `v1-b1` tag) -> E -> CLEAN waves 2-3 and asked which check sees a licensed file once it leaves `experiments/`. None does (item 1). |
| engineering:architecture | Seam check of the box code snapshot: `code_files()` drops `docs/` and (after CLEAN1) `archive/`, so the only RIGA leak test runs on a view of the tree that the moved files never enter (item 1). |
| engineering:documentation | Stale "after D3" attributions in E for the numpy polygonizer oracle (item 2). |
| huggingface-skills:huggingface-best | No model choice in scope. Nothing changed. |
| system-design-skills:system-design | Same order walk, framed as "where does each data file live after each wave, and who guards it". |
| system-design-skills:service-decomposition | One package with import-linter layers (a modular monolith) is still the right form; no split proposed. |
| mattpocock-skills:codebase-design | Deletion test on `final_eval.py` after E4: the names left in the script must be listed, or U10 cannot say where they go (item 3). |
| WebFetch | git documentation, gitignore: "If there is a separator at the beginning or middle (or both) of the pattern, then the pattern is relative to the directory level of the particular `.gitignore` file itself" [web: https://git-scm.com/docs/gitignore]. So `/experiments/phase0/results/riga/` cannot match the same folder after a move (item 1). |

## Items

### 1. blocking: CLEAN5 moves 31 git-ignored RIGA-derived files to paths that `.gitignore` does not cover, and CLEAN5 step 8 then stages them for the public repo.

Sources: CLEAN (CLEAN5, CLEAN12, CLEAN.9.3 CT-73 and CT-87).

- CLEAN.md "### CLEAN5": "Files to move (plain `mv`): ... `experiments/phase0/*.md` and `experiments/phase0/results/`
  -> `docs/research/records/phase0/` ... `experiments/phase1/lab/` (all remaining files) -> `archive/phase1_lab/`."
  Step 8: "Owner runs ...: `git add -A experiments archive docs/research/records tests pyproject.toml src`".
- The repo keeps RIGA-derived results out of git only by path. `.gitignore:37-49` [code]: "# RIGA-derived results
  stay local: no permission to share RIGA data (repo public, owner 2026-11-06)", then 12 patterns, all anchored to
  the old trees: `/experiments/phase0/results/medsam_riga.csv`, `/experiments/phase0/results/riga/`,
  `.../riga_outlines/`, `.../vlm_box/`, `.../vlm_polygon/`, `/experiments/phase1/lab/t_reward/results/h5_items.csv`,
  `/experiments/phase1/lab/t_data/riga/out_maxrg_check/`, `.../t_data/riga/maxrg_vs_prime.json`,
  `.../t_data/splits/out/riga_images.csv`, `.../t_reward/results/advantage_groups.json`,
  `.../t_reward/results/iter3_r5_raw.json`, `.../t_stack/results/t[24]_*.json`.
- On disk these patterns hide 31 files today (11 under `experiments/phase0/results/`, 20 under
  `experiments/phase1/lab/`) [measured: `git ls-files -o -i --exclude-standard experiments/phase0
  experiments/phase1/lab`, minus `__pycache__` and `*.log`; `r13arch_probe/riga_ignore.txt`]. A plain `mv` of the
  folder takes them along. At the new paths they are not ignored: `git check-ignore -q` gives 0 for
  `experiments/phase0/results/medsam_riga.csv` and 1 for `docs/research/records/phase0/results/medsam_riga.csv`; 1
  for `archive/phase1_lab/t_data/riga/maxrg_vs_prime.json`, `docs/research/records/phase0/results/vlm_box/replies.jsonl`
  and `archive/phase1_lab/t_stack/results/t2_qwen3-vl-2b_simple.json` [measured: same file]. A leading slash makes a
  pattern relative to the `.gitignore` folder, so no anchored line can follow a move [web: git-scm.com/docs/gitignore].
  Step 8's `git add -A archive docs/research/records` then stages them.
- The plan's own rule: CLEAN.md "## CLEAN.2" row "RIGA terms": "'copy, redistribution ... prohibited'; RIGA-derived
  tables already left the public MIT repo. So no tracked golden holds RIGA masks or per-case RIGA results". Front 6.4
  rule 7 holds the same rule for the box.
- No planned $0 check sees it:
  - `tests/test_p2_breakit_ops.py::test_code_snapshot_has_no_riga_derived_results` scans `P.code_files()` only
    [code: test_p2_breakit_ops.py:297-300]. `code_files()` drops `docs/` [code: policy.py:379], and CLEAN1 adds
    `archive/` to the drop. CLEAN5 step 4 itself rewrites line 299 to `docs/research/records/phase0/results/vlm_` and
    notes "(dropped by the `docs/` rule)", so that line can no longer fail.
  - CT-73 checks `git ls-files | grep -iE 'riga.*\.(png|tif|npy|npz)'`: no `.csv`, `.json` or `.jsonl`.
  - CT-87 checks counts and locations of the moved records, not their ignore state. CLEAN-G16 reads only
    `tests/goldens/*.npz`.
- So an engineer who follows CLEAN5 exactly commits RIGA-derived tables (per-case mask tables, VLM replies on RIGA
  images, per-model results) to the public repository. That breaks the owner rule written in `.gitignore` and the
  RIGA terms that CLEAN.2 quotes, and a push cannot be taken back.

Fix (CLEAN.md CLEAN5, CLEAN12, CT-73):
1. New CLEAN5 step 0, before any `mv`: record the ignored files of the trees that move,
   `git ls-files -z -o -i --exclude-standard experiments/phase0 experiments/phase1/lab experiments/phase2/lab`
   minus `__pycache__/` and `*.log`, into `runs/clean/ignored_before.txt` (31 paths today; read-only git).
2. In the same CLEAN5 commit, before the `mv`, rewrite `.gitignore:38-49` to the new paths (keep the comment line):
   `/docs/research/records/phase0/results/medsam_riga.csv`, `/docs/research/records/phase0/results/riga/`,
   `/docs/research/records/phase0/results/riga_outlines/`, `/docs/research/records/phase0/results/vlm_box/`,
   `/docs/research/records/phase0/results/vlm_polygon/`, `/archive/phase1_lab/t_reward/results/h5_items.csv`,
   `/archive/phase1_lab/t_data/riga/out_maxrg_check/`, `/archive/phase1_lab/t_data/riga/maxrg_vs_prime.json`,
   `/archive/phase1_lab/t_data/splits/out/riga_images.csv`, `/archive/phase1_lab/t_reward/results/advantage_groups.json`,
   `/archive/phase1_lab/t_reward/results/iter3_r5_raw.json`, `/archive/phase1_lab/t_stack/results/t[24]_*.json`.
   (Alternative with the same effect: move the 31 files into `data/local_records/` first, since `/data/` is ignored,
   and delete the 12 dead lines.)
3. CLEAN5 step 9 gate, run before step 8's `git add`: for every path in `ignored_before.txt`, mapped through the
   CLEAN5 move map, `git check-ignore -q <new path>` exits 0; `git ls-files -o --exclude-standard archive
   docs/research/records` prints none of the mapped paths. Add `.gitignore` to step 8's `git add` list.
4. CLEAN12: add `test_no_riga_derived_file_is_committable` to `tests/test_repo_hygiene.py`. Over
   `git ls-files -z -co --exclude-standard` of the whole repo (not `code_files()`), no path whose lower-case form
   contains `riga` ends in `.csv`, `.json` or `.jsonl`, and no path is under `docs/research/records/phase0/results/vlm_`.
   Allow by literal the two tracked files `archive/phase1_lab/t_data/riga/inventory.json` and `repo_parity.json`
   (today `experiments/phase1/lab/t_data/riga/{inventory,repo_parity}.json`, tracked [measured: `git ls-files`]). A
   negative case in a temp git repo: an untracked `docs/x/riga/a.csv` is reported.
5. CT-73: widen the grep to `riga.*\.(png|tif|npy|npz|csv|jsonl?)$` with the same two-file allow list.

### 2. minor: E still says the numpy polygonizer oracle moves at D3.

Sources: E (E3a table note, E.5, E3 seam paragraph).

- E.md E3a move-table row "(none)": "E-G16 lists them as \"moved by D3\" (both are gone from `src/` on the E0
  commit)". E.md "## E.5": "The numpy oracle is already `tests/oracles/polygonize_oracle.py` after D3." E3 text:
  "minus names the move table marks as archived, renamed or moved by D3".
- After round 12, B1.4a moves `checker/polygonize.py` (B.md B1.4a: "`git mv src/objdel/checker/polygonize.py
  tests/oracles/polygonize_oracle.py`"); D3 moves only `bitquad.py`. `old_names.json` is recorded on the E0 commit,
  when both files are already gone, so no test result changes.

Fix: E3a row: "E-G16 lists them as moved before E (B1.4a: `polygonize.py`; D3: `bitquad.py`)". E.5: "after B1.4a".
E3: "archived, renamed or moved before E0 (B1.4a, D3)".

### 3. minor: E2's name table for `final_eval.py` lists no home for the names A3 and B1.6 add.

Sources: E (E2 step 6), A (A3), B (B1.6).

- E.md "### E2" step 6, table row `final_eval.py`: moves `TERMS`, `per_item`, `summarize`; "Stays in the CLI: `main`".
  The rule above the tables: "Each top-level name of the script on the E0 commit has one home."
- On the E0 commit `final_eval.py` also defines A3's `view_content_sha256(view)` (A.md A3: "one function,
  `final_eval.view_content_sha256(view)`, that B1.6 reuses") and B1.6's `frozen_entry(view_name)`, which reads
  `Path(__file__).resolve().parents[2] / "docs/reference/frozen_views.json"` (B.md B1.6), plus the `--hash-view` and
  `--rescore` helpers.
- U10 (`every_name_has_a_home`) fails on the Mac at $0, so minor. But the engineer must choose a home. If
  `frozen_entry` moves into `src/objdel` as written, `parents[2]` points to `src/`; E2.5's AST test (0
  `__file__`-based roots under `src/objdel`) catches that, also at $0.

Fix (E.md E2 step 6 table): add `view_content_sha256`, `frozen_entry` and A3's `--hash-view` and `--rescore` helpers
to the "Stays in the CLI" cell of `final_eval.py`, with the note "`frozen_entry` keeps `parents[2]` only while it
lives in `experiments/phase2/`; if a `src/` caller ever needs it, it moves to `policy/evaluate.py` and uses
`core.provenance.repo_root()`".

## Counts

blocking 1 (item 1), major 0, minor 2 (items 2-3).
