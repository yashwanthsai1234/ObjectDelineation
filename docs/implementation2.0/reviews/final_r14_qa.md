# Final review round 14: QA and technical writing (Tests sections, gates, testability, cross-section consistency)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | I recommend no NVIDIA catalog skill: no NVIDIA product choice is in scope. The skill pointed me at the CUDA engine's test-only build variants and the compute-sanitizer step. I checked the sanitizer time budget against the K1 run log: the full group under memcheck took 47.8 s and under initcheck 55.1 s, below the 300 s cap [measured: runs/k1/cuda_checker/memcheck2.txt]. So there is no finding there. |
| engineering:code-review | I read `cuda/gpucheck.cpp:131-156` (`Block::get` grows and frees; `Lane`), the cited D2H lines 991, 1572, 1623, 1638, 1681 and 1799, `K1/k1poly.cpp:138-146`, `CMakeLists.txt`, the `Makefile`, `.github/workflows/tests.yml`, `pyproject.toml` (pytest options, extras, scikit-build `build-dir`), `.gitignore:14-49`, `experiments/phase2/vast/policy.py:498-507` (`cmd_pull` destination), `tests/test_cuda_checker.py` (imports and `-k masks` set), `tests/conftest.py:46-54` (synthetic RIGA items) and `tests/test_p2_breakit_ops.py:297-300`. |
| engineering:debug | I traced the round 13 variant mechanism step by step: build, then the `-I -S` child, then the load guard, then the pytest exit code. I traced the `alias` run through `Block::get` growth (item 5). I traced where the V1 A3 answers land on the Mac (`runs/<run>/`), against where A and C read them (`runs/a3/`) (item 2). |
| engineering:system-design | One owner per fact. The pinned-staging bound must be a constant that the engine knows at its first call (item 4). A path for V1 outputs must follow the rule that V2 and V3 already use, the `OBJDEL_V<n>_RUN` variable (item 2). |
| engineering:architecture | I checked the licence boundary of the public repo. LIDC-derived goldens are allowed (CC BY 3.0). RIGA-derived per-case data is not (`.gitignore:37`). The same boundary applies to D's always-run goldens as to CLEAN4's golden (item 1). |
| engineering:documentation | Carry-over check: I re-read the round 13 qa minors (items 3-9). All of them are still open in the sources (list at the end). |
| huggingface-skills:huggingface-best | No model choice is in scope. No change. |
| engineering:testing-strategy | For each new round 13 test (S16, I20, CT-89, CLEAN-G19), I asked three questions. Can it fail? Does a named gate run it? Does it still find its files after E8 moves the tests (item 3)? |
| mattpocock-skills:tdd | A planted defect must give a red run first, with a stated signal. For the `alias` variant, the signal "exit 1" is not certain (item 5). The expected values in D-R01 must come from data that the repo may publish (item 1). |

Web tools: loaded (WebSearch, WebFetch). One web check: the LIDC-IDRI licence on TCIA is CC BY 3.0
[web: https://wiki.cancerimagingarchive.net/x/hIRXAQ; https://huggingface.co/datasets/mrmrx/CADS-dataset/blob/68e55b5138c78cda4ece6bd98a2006b6df204b47/0016_lidc/README_0016_lidc.md].
So LIDC-derived goldens may be published. RIGA data may not: the owner rule is "RIGA-derived results stay local: no
permission to share RIGA data (repo public, owner 2026-11-06)" [code: .gitignore:37].

Measurements this round ($0, Mac, 2026-11-09):
- `docs/implementation2.0/Implementation.md` is 1,752,215 bytes and `review-log.md` is 375,248 bytes (`ls -l`). Both are
  under the DOC-G13 limit of 2 MB.
- No match of the check 14 pattern (`/[p]rivate/tmp`, `/[U]sers/`, `docs/implementation2.0`, `S[P]/[A-Z]`) in any `.md` file or kernel
  file under `docs/implementation2.0/` (`grep -rlE`).
- `git ls-files -co --exclude-standard | grep -i riga | grep -iE '\.(csv|json|jsonl)$'` gives only the two tables
  that CT-89 allows by literal (`t_data/riga/{inventory,repo_parity}.json`). So CT-89's name rule passes on today's
  tree.
- Sanitizer time: the K1 box run ran the full checker test group under memcheck in 47.8 s and under initcheck in
  55.1 s [measured: runs/k1/cuda_checker/memcheck2.txt]. That is well inside D8's `timeout 300` per tool.

## Items

### 1. blocking: D-R01's "always-run copy" commits 50 RIGA items (inputs and target texts) to the public repo, and no planned check catches it (D)

- Doc (D.md, D.9.4 "Regression, golden and data-backed", row D-R01): "CPU build on the Mac only (`data`; ... RIGA
  and most of these sets are not in the R-D2 payload). Always-run copy: 50 fixed LIDC and RIGA items in
  `tests/goldens/` equal the stored pre-D text".
- An "always-run" test carries no `data` marker. So it runs in `make test` and in CI, where `data/` does not exist.
  The golden therefore holds two things: the RIGA inputs (the rater outlines or masks of 50 real RIGA images) and
  the stored target text, which is a polygon text made from those raters. Both are per-case RIGA-derived data.
- Owner rule: "RIGA-derived results stay local: no permission to share RIGA data (repo public, owner 2026-11-06)"
  [code: .gitignore:37]. The plan states the same rule elsewhere, and keeps its other goldens LIDC-only:
  - D8: `make_d_goldens.py` "(CPU build, ..., LIDC only)".
  - A.9.6 G02b: `make_a_goldens.py` "LIDC only".
  - 6.4 rule 7: "RIGA never leaves".
  - CLEAN.R item 1 (blocking, accepted): "CLEAN4 golden would publish RIGA masks".
- No planned check catches a D-R01 golden:
  - CT-11 / CLEAN-G16 reads only `tests/goldens/*.npz` (CLEAN.md CLEAN-G16: "Arrays in `tests/goldens/*.npz`").
  - `test_code_snapshot_has_no_riga_derived_results` and the CT-89 name rule look only at the file name: `"riga" in
    f.lower()` [code: tests/test_p2_breakit_ops.py:298].
  - So a golden named, for example, `tests/goldens/d_r01_targets.json` passes every check. It reaches the public
    repo at the D3 commit.
- The repo's existing RIGA tests use synthetic items (`riga/MESSIDOR/image00<s>`, "3 raters, nested disc and cup"
  [code: tests/conftest.py:46-54]). So a compliant form already exists.
- Fix (D.md D-R01, and the D-R02 wording):
  - D-R01: "Always-run copy: 50 fixed LIDC items (LIDC-IDRI, CC BY 3.0) in `tests/goldens/d_r01_lidc.json` (input
    masks as run-length lists plus the pre-D text). RIGA items are rebuilt only in the `data` half on the Mac and are
    never written to `tests/`. The test asserts that no golden item id starts with `riga`."
  - D-R02: "RIGA disc/cup/rim" becomes "synthetic disc/cup/rim shapes (the `tests/conftest.py` RIGA generator), no
    real RIGA item".
  - Optional: extend CT-11 to `tests/goldens/*.json` item ids (no `riga` source id).

### 2. minor: A and C read A3's V1 answers from `runs/a3/eval_val/`, but V1's `pull.sh` writes `runs/<V1 run>/eval_val/`, and no step maps one to the other (A, C, 00_front)

- Doc (00_front 6.2, `a3_measure.yaml`): "`lidc-val-shifted` with saved answers, ... (item 5, `--out
  {rd}/eval_val`)". `{rd}` is `/workspace/runs/<run>`, and `pull.sh` copies it to the Mac as `runs/<run_id>/`
  [code: experiments/phase2/vast/policy.py:502-506, `dst = ... ROOT / "runs") / run_id`].
- The readers use a fixed `runs/a3/`:
  - A.md A2 CLI: `--answers runs/a3/eval_val/answers.jsonl`.
  - A.md A8 step 1: the same path.
  - C.md C2 "Mac G-EQ2 data": `(runs/a3/eval_val/answers.jsonl)`.
  - C.md C11 step 1: "A3's val answers (`runs/a3/eval_val/answers.jsonl` with `tokens[K]`)".
- A.9.6 G01 and G02a name no folder at all. A.9.6 also says "A test whose `runs/...` input is absent skips". So a
  wrong guess shows as a skip in `make test-data`, not as a fail.
- V2 and V3 already have a rule: `OBJDEL_V2_RUN=<run>` and `OBJDEL_V3_RUN=<run>` (front 6.3). V1 has none.
- Severity: the CLIs fail with "file not found" on the Mac at $0, so the engineer finds the gap before any spend.
- Fix:
  - Front 6.2, after the V1 launch lines: "After the V1 `pull.sh`, on the Mac: `OBJDEL_V1_RUN=<run> make
    test-data`."
  - A.9.6 G01/G02a: read `runs/$OBJDEL_V1_RUN/{eval_val,eval_val_edge}/`, and skip with "no V1 run pulled" when the
    variable is unset (the S13 pattern).
  - A2, A8, C2 and C11: write the path as `runs/$OBJDEL_V1_RUN/eval_val/answers.jsonl` (V2 item 1b:
    `runs/$OBJDEL_V2_RUN/...`), or add one Mac step that links `runs/a3` to the pulled run folder.

### 3. minor: E8's list of hard-coded test paths misses D's `tests/variants.py`, and `test_test_paths.py` does not scan `tests/` (D, E)

- Doc (D.md D2 step 7): `tests/variants.py` runs "`tests/test_cuda_checker.py -k masks`" in each variant child. D8
  adds a CI step `.venv/bin/python -I -m tests.variants`.
- Doc (E.md E8 step 3, "Hard-coded test paths"): the list is `pyproject.toml:58`, `tests.yml:54`, `preflight.py:29-34`,
  `policy.py:451-453`, `rl_full.yaml` and four configs. `tests/infra/test_test_paths.py` scans "`pyproject.toml`,
  `.github/`, `experiments/` (live), `specs/`, `src/` and `Makefile`", but not `tests/`.
- After E8 moves `tests/test_cuda_checker.py` into `tests/verify/`, the variant child gets pytest exit 4 (file not
  found). The pass rules then fail: `order1`, `order2` and `poison` need 0, and `alias` needs 1.
- Severity: `make test-order` and the CI step go red at $0, so this is caught before spend.
- Fix:
  - E8 step 3: add "`tests/variants.py` (D2 step 7: the `-k masks` file path)" to the hard-coded list.
  - `test_test_paths.py`: also scan the non-test helper modules `tests/*.py` that do not start with `test_` (today
    only `tests/variants.py`, `tests/conftest.py` and `tests/probe_models.py`).
  - Better: `tests/variants.py` takes the test file path from one constant that E8's move table rewrites.

### 4. minor: D5c's pinned bound uses "names per answer" and C_max from a run file, but the engine must set the bound inside its first call, whichever entry that is (D)

- Doc (D.md D5 step 5): "The engine sizes `in` and `out` of every lane to these bounds once, inside the first engine
  call of the process, whatever that entry is". The bounds are:
  - `in`: "256 answers x names per answer x `MAX_VERTICES` x `MAX_RINGS` x 8 bytes ... 20.5 MB per name".
  - `out`: "n_max = 256 x names ... C_max = 64 components per object [estimate], or twice the largest per-object
    count in the D1 replay set (`sig_start` steps in `runs/d1/bench.json`), whichever is larger".
- In RL, the first call is `add_tables_from_raters` (D5c "Why"). That call has no answers and no object names. So
  "names" has no value at that moment. LIDC answers have 1 name. RIGA answers have 2 names, `disc` and `cup`
  [code: src/objdel/train/dataset.py:75]. Also, `runs/` is git-ignored [code: .gitignore:23], so the engine cannot
  read C_max at run time.
- I20's control case ("a `score_polygons` call above the `out` bound adds exactly 1") needs the bound as a number,
  but the doc gives none.
- Severity: D-G15 is read only on the LIDC B3 run, where any names value >= 1 passes.
- Fix (D5c):
  - "`NAMES_MAX = 2` (RIGA disc and cup; the rim is derived) and `C_MAX` are compile-time constants in
    `gpucheck.cpp`. D1 sets `C_MAX = max(64, 2 x the replay maximum)` and writes the number into the ADR. A binding
    `gc.pinned_bounds() -> (in_bytes, out_bytes)` exposes both values."
  - I20's control case: "`score_polygons` with n = `out_bytes // 104 + 1` one-vertex answers (defects term 13 x 8 B
    per object)".

### 5. minor: the `alias` variant's red signal "exit 1" is not certain: the planted `scan` call is not named, and a growth inside `Block::get` frees layout A (D)

- Doc (D.md D2 step 7):
  - "`GC_POISON_SCAN` (the CPU `scan` calls `tmp.get(4 * (n + 1))` and fills that span with 0xFF ...)".
  - "`GC_ALIAS_SCAN_FOR_TEST` (one `scan` call in `chunk` is given `L.dev` in place of `L.scan` ...)".
  - "`alias` ...: the same run MUST exit 1 (tests failed, load guard passed). Exit 0 or 3 fails `make test-order`".
- Code: `Block::get` keeps the buffer only when `bytes <= cap`. Otherwise it releases the old buffer and allocates a
  new one [code: cuda/gpucheck.cpp:139-148].
- Two failure cases follow:
  - If the aliased call's `4 * (n + 1)` is larger than `L.dev`'s capacity, the poison frees layout A. The CPU build
    then reads freed memory. That is undefined behaviour, and a crash gives a negative return code, not 1.
  - If the aliased call's first `4 * (n + 1)` bytes of `L.dev` are dead at that point, the run stays green (exit 0).
- Both cases are red or visible at $0. But the doc gives no rule that picks a call that is certain to give exit 1.
- Fix (D2 step 7): "`GC_ALIAS_SCAN_FOR_TEST` aliases the first of the six scans (k1poly.cpp:824), where layout A's
  first array is live and larger than the scan span. The alias variant also asserts `4 * (n + 1) <= L.dev.cap`
  before the poison write. Pass rule: pytest exit 1 with at least one failed `-k masks` case named in the summary."

## Carry-overs still open (round 13 qa, minor; not re-argued)

- qa-3: B1.3 and DG1b still say "on the box" / "(box: the manifest's `excluded_ids`)" (B.md B.6.1 DG1b row).
- qa-4: C.9.4 still has no `limit_mm_per_prompt={"image": 0}` in the `vllm.LLM` call (C.md C.9.4).
- qa-5: C.9.4 still says "15.4 min < 30.9 min" (C.md:1968).
- qa-6: under decision 1 (a), it is still not stated that `tests/test_c_gpu.py` lands with B3.9 (B.md, C.md).
- qa-7: the C10 keep list still omits `--encode-check` and `--checks-from-step` (C.md:1050).
- qa-8: E.md:645 and E.md:747 still say "moved by D3", E.md:1760 still says "after D3", and the CLEAN-G14 fail
  column still says "(D3 not landed)" (CLEAN.md:1052).
- qa-9: 00_front.md line 4 still says "Revision: final round 3 plus the round-4 consistency pass".

## Checked, no finding

- The round 13 variant mechanism (D2 step 7, D8, S16):
  - `-I -S` with `sys.path` set inside the child keeps out the editable finder, which lives in a `.pth` file.
  - A regular `objdel` package in `build/variants/<name>/site` comes before the `site-packages/objdel/` folder of the
    editable install.
  - `build/` is git-ignored [code: .gitignore:16].
  - The CI lint + pytest job installs `-e ".[dev]"`. That pulls the base dependencies (opencv, scipy, numba), so
    I10's cv2 leg and the `-k masks` tests import in the child.
  - `tests/` has no `__init__.py`, so `python -I -m tests.variants` resolves as a namespace package through the
    repo-root `.pth` file.
  - `-C build-dir=...` overrides `[tool.scikit-build] build-dir = "build/{wheel_tag}"`.
- The round 13 pinned rule (D5c, D7b, I20, D-G15):
  - The cited D2H lines all use `L.out` today [code: gpucheck.cpp:991, 1572, 1623, 1638, 1681, 1799].
  - The 20.5 MB arithmetic is right: 256 x 1,000 x 10 x 8 = 20,480,000 B.
  - The defects arithmetic is right: 256 x 13 x 8 = 26,624 B.
  - I20, the D-G15 Mac half, GAP-24 and the coverage maps agree with each other.
- The round 13 CLEAN5 RIGA fix (step 0, step 9 gate, CT-89, CLEAN-G19):
  - The old 12 lines are anchored [code: .gitignore:38-49].
  - The new block maps the old paths one to one.
  - CT-89's name rule passes on today's tree (measurement above).
  - Step 8 runs only after the step 9 gate.
- The D test ids: every id that the D.9.8 coverage map cites is defined in D.9.1-D.9.7 (script check). The
  regression rows use `D-R0n` and the map uses `R0n`, which is cosmetic. The A test ids are all defined (script
  check).
- DOC4 check 3 (`make` mentions) skips the plan files (`PLAN_GLOBS`). So the plan's `make deploy` and `make nope`
  examples cannot fail DOC-G3.
