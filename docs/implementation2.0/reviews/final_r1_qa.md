# Final review round 1: QA lead and technical writer (ASD-STE100)

Target: REPO/docs/implementation2.0/Implementation.md (11,117 lines), with the kernel folders in
REPO/docs/implementation2.0/kernels/. Earlier final reviews: none (round 1). No REPO file was edited. No git write. No rental.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog lookup was needed: the doc already names the one match (`compute-sanitizer-diagnose-cuda-bugs`) as a recommendation only. I checked the toolchain and sanitizer claims against `run_box.sh` and the D8 sanitizer chain instead (items 1, 4). |
| engineering:code-review | Each claim was checked against the code. Results: `configs.py:73-77`, `gpu-run.md:129,137`, `run_box.sh:24`, the Makefile `cov` target, the pyproject markers, the R5 numbers in `final_eval.json`, `items.jsonl` line counts, `shared_rl.json` per_step 23.65. These became items 8, 11, 15 and 21. |
| engineering:debug | I traced each gate to its input files and asked: can it fail, can it pass, and can it run in the stated order? This found the gate that cannot pass in DOC1 (item 4) and the G1 rental that cannot run (item 1). |
| engineering:system-design | I put the phase order, the rentals and the dependencies on one time line. This found the order conflicts in items 1, 3 and 12. |
| engineering:architecture | I checked the ADR numbers against `docs/architecture/decisions/` (0011 is the last one). This found item 7. |
| engineering:documentation | I checked reader-first order, one source per number, and links to files that will vanish. This found items 5, 11, 14 and 17. |
| huggingface-skills:huggingface-best | I checked the model ids and sizes on the HF API: Qwen3-VL-4B-Instruct, MedMO-4B and MedMO-4B-Next are each 4,437,815,808 params. This found item 18 (two names for one model). |
| engineering:testing-strategy | I matched each Tests section to its gates: markers, baselines, data isolation and coverage maps. This found items 8, 9, 10 and 16. |
| mattpocock-skills:tdd | I looked for tests that pass however the code behaves, assertions that take either outcome, and seams nobody agreed. This found items 10 and 22. |
| Web (WebSearch/WebFetch) | Confirmed: import-linter has `exhaustive` on layers contracts (v2.8+). Zensical reads `mkdocs.yml`. arXiv 2601.13633 is EGM. |

Severity: blocking = wrong, unsafe or breaks an owner rule; major = missing or vague, so an engineer cannot implement it;
minor = clarity or consistency. Count: **4 blocking, 13 major, 11 minor**.

---

## Blocking

### 1. [blocking] [00_front, A, C, D] The G1 rental in section 6.2 cannot run as planned, and the phases plan other rentals
Problem.
- Section 6.2 puts A3, K1-K3, D `make gpu-check` (with G7/G8) and the C arms into one rental, capped at 2.9 h ($1.74).
- But item 5, D `make gpu-check`, needs the D2-D8 code and the pre-D2 engine (D-G7 builds `git archive <pre-D2 sha>`). In the section 1 order, D comes after the B3 run, and G2 waits for G1. So that code does not exist at G1 time.
- D itself plans two rentals: R-D1 (K1, K2, bench) and R-D2 ("run it first in the next Phase B/C rental"), with its own budget gate D-G13 ($1.40).
- A3 has its own launch block: `launch.sh OFFER a3 ... --cap-hours 1.5` (A3 header). The A "All" gate says "Rental spend <= $1.00 (cap-hours 1.5)".
- C plans its own session of $0.70-1.15 (C.8.2).
- Section 6.2 says "Item 7 also gives the A6 nsys baseline (same capture)". But A-G6 needs nsys on steps 4-5, compared with steps 3 and 6. The C0 arms capture steps 16-17.
- The A3 command block also leaves out `sweep_orphans.sh`, which section 6.4 rule 6 requires.

Exact fix.
- Make section 6 the one rental table. Split G1 into two rentals:
  - G1a, before B1: A3, the K1/K2/K3 `GPU_RUN.md` items, R-D1 `bench_checker.py --engine current --stages`, `c_speed_check` and the C0/fast/check arms.
  - G1b = R-D2: after D8 lands, run as the first item of a later B rental.
- Recompute the G1a minutes and cap without item 5. Move item 5 to G1b.
- In A3, replace the `launch.sh ... --cap-hours 1.5` header with "runs as G1a item 8; see section 6.2". Add `&& sweep_orphans.sh` after `teardown.sh`. Change the A "All" gate to "A3 spend is inside the G1a cap".
- In A-G6 and A6, choose one capture window:
  - (a) A-G6 takes steps 16-17 from the C0 baseline arm, with steps 15 and 18 as the overhead compare, or
  - (b) keep a separate 6-step A6 window inside G1a and add its minutes to the table.
- Change C.8.2 and D.8 to say "inside G1a/G1b (section 6)" and not their own rentals.
- Back matter decision 3: remove "share R-D2 inside G1". Replace it with "R-D2 rides the first B rental after D8".

### 2. [blocking] [00_front] The front states the wrong gate for "B uses the C code"
Problem.
- Section 1 says "Target 18 s per step or less, same rewards (Gate C)". The section 1 diagram says "Gate C: <= 18 s/step ... no -> run on old code".
- C.6 and C.8.3 say something else. B uses the C code when G-EQ1..6 and G-SP1 pass. G-SP1 is <= 20.0 s in mode A or <= 20.7 s in mode B.
- The 18 s value is G-C, and G-C counts only with C7 and C8 added inside Phase B.
- Other numbers disagree too: section 6.3 uses "20 s (Gate C passed)", and B.8.1 prices "Phase C's 17 s".
- An engineer who follows the front will throw away C code that passed at, for example, 19.5 s.

Exact fix.
- Section 1 bullet C: write "Gate C-SP1: <= 20.0 s (mode A) or <= 20.7 s (mode B) with G-EQ1..6 equal. G-C 18 s applies only after C7+C8 inside Phase B."
- Diagram: change the box text to "G-EQ1..6 + G-SP1 pass?".
- B.8.1: replace the "17 s" column with "20.0 s (G-SP1)", or label 17 s as "Direction estimate, not a gate".

### 3. [blocking] [00_front, C, D, E] The phases give opposite orders for the E moves and the C/D edits
Problem.
- Front order rule 2 says "E moves files after B1, C and D land".
- C.8.3 says "land C first, then E moves the new modules unchanged". D.7 R10 says "E's moves come after D3/D4b".
- E says the opposite:
  - E.7: "E5 before C coding starts, so C edits the split files".
  - E.8: "Phase C edits the split files after E5"; "B2 ... and B3 ... land after E3/E5"; Phase D "deletes `verify/pool.py` ... after E3".
- These cannot all be true. Each order changes which file paths the C, D and B3 tasks name.

Exact fix.
- Choose one order (the front rule matches C and D: B1 -> C1-C5 -> D2-D4b -> E0-E8 -> B2/B3 code in the E tree).
- Write it once as a numbered time line in section 1.
- Rewrite these to match: E.7 row 3, the E.8 bullets for B, C and D, and the D.8 bullet that names `verify/pool.py`.
- If E5 must come before C, then rewrite all the C file paths (`fast.py:...`) to the E paths instead.

### 4. [blocking] [DOC] The DOC1 commands cannot make gate DOC-G1 pass, and the obvious fix would publish binaries
Problem.
- The DOC1 `git add` list leaves out `docs/implementation2.0/`.
- The DOC1 step 3 and DOC-G1 check compares `git ls-files docs | grep -v ^docs/research/ | wc -l` with `find docs -type f -not -path 'docs/research/*' | wc -l`. That `find` counts every file under `docs/implementation2.0/`, so the check fails.
- If the owner adds that folder to fix the check, git takes 13.5 MB of build output:
  - `docs/implementation2.0/kernels/K1/build` (3.5 M) and `K2/build` (10 M) hold arm64 Mach-O binaries, `.dSYM` folders and docker logs.
  - `.gitignore` ignores only the root `/build/`.
- `K1/GPU_RUN.md`, `K2/GPU_RUN.md` and `K3/GPU_RUN.md` hold absolute paths with the owner's user name and the `/[p]rivate/tmp/...` scratch path. The repo is public.
- The DOC1 pre-publication scan does not cover these files.

Exact fix.
- Add `docs/implementation2.0/kernels/*/build/` to `.gitignore` and delete the three `build/` folders.
- In the three `GPU_RUN.md` files, replace the absolute paths with `$REPO` and `$K1/$K2/$K3 = $REPO/docs/implementation2.0/kernels/K*`.
- Add `docs/implementation2.0` to the DOC1 `git add` line.
- Add a scan step: `grep -rn "/[U]sers/\|/[p]rivate/tmp" docs/implementation2.0` must print nothing.

---

## Major

### 5. [major] [00_front, A, B, C, D, E, DOC, CLEAN, 99_back] The doc depends on scratch files that will vanish
Problem.
- The doc has 125 `S[P]/` references and 11 absolute `/[p]rivate/tmp` or `/[U]sers/` paths. SP is the session scratchpad under `/[p]rivate/tmp`. macOS clears it, so these links break.
- Load-bearing content lives only there:
  - the test-writer contracts `A.spec.md` to `E.spec.md` and `DOC.spec.md`. Examples: "Output keys: listed in A.spec.md §3.2", "interface in D.spec.md §3.1", DOC-G3 "all 13 checks of DOC.spec.md B5", Gate B3 "T1-T13 pass (B.spec.md)";
  - source code: B1.3 says "Start `quality.py` from `docs/implementation2.0/seeds/b1_scripts/*.py`";
  - measured evidence: `valsplit.json`, `CL1_baseline_pytest.txt`, the M/T/CL notes;
  - D1 step 4-5 copies results "back to docs/implementation2.0/kernels/K1/results/".
- Section 0 never defines SP or REPO. Each phase defines them again in a different way.
- This breaks the doc's own rule (CLEAN14: link SP notes only by their copy in `docs/research/notes/`).

Exact fix.
- Copy the 6 `*.spec.md` files to `docs/implementation2.0/specs/`. Copy the cited notes to `docs/research/notes/impl_2026-11-08/`. Copy `b1_scripts/` to `experiments/phase2/lab_b1/` (or inline the needed functions).
- Replace every `S[P]/...` path with its repo path. Add "Paths are relative to the repo root" to section 0.
- Change D1 step 4 and the D.5 table to `docs/implementation2.0/kernels/K*/results/` and `runs/d1/`.
- In the 99_back "Scratchpad notes" table, point at the promoted copies.

### 6. [major] [00_front, B, D, E, DOC, CLEAN] One id means different things in different phases
Problem.
- "G1" is the validation rental (section 6.2), a D gate, an E gate, a DOC gate and a CLEAN gate. "G5" is a D gate (sanitizers), an E gate (TRL), a DOC gate (stale strings) and a CLEAN gate.
- "D1-D16" are the B data gates. "D1-D9" are also the Phase D tasks. Section 1 says "Gate B1 data (D1-D16)" and "D gpu-check" in the same diagram.
- "S1/S2" are B cost sessions and also the A/B test-set origins.
- "R1-R14" are front risks, "R1-R13" are D risks, "R5" is a run and "R-D1" is a rental.
- A reader cannot tell which "G5" a cross-reference means. The back matter "G1 items 2-4" and A-G* gates already mix these.

Exact fix.
- Give every id its phase prefix: A-G1..8, B-DG1..16 (data gates), B-G1..3 (run gates), C-G-EQ1.., D-G1..16, E-G1..16, DOC-G1..13, CLEAN-G1..18.
- Rename the rentals V1 (validation) and V2 (B1 run). Rename the B sessions B-S1 and B-S2. Rename the front risks FR1..FR14.
- Apply the change with one search-and-replace pass per phase.
- Add the prefix rule to section 0 and to section 8.

### 7. [major] [B, D, E, DOC] Four phases claim the same ADR numbers
Problem.
- Today the last ADR is `0011-lora-weight-sync-without-merging.md` [measured: ls docs/architecture/decisions].
- ADR-0012 is claimed four times: B (numerics preset, "ADR-0012 rule"), D9 (`0012-cuda-polygonizer-in-production.md`), DOC (`0012-docs-site-with-zensical`) and E (`0012-layers-and-import-linter`, plus "ADR 0012" for `objdel.infra` in the wheel).
- ADR-0013 is claimed by B (reward v5) and E (`0013-typed-run-spec-pydantic`, plus "ADR 0013" vllm_native).
- ADR-0014 is claimed by B (polygonizer exception) and E (`0014-trl-adapter-and-startup-check`). E also claims 0015.
- The files would overwrite each other, and DOC-G11 ("0 index rows without a file") cannot be checked.

Exact fix.
- Add an ADR allocation table to section 4 (or DOC10) in time order of the work. For example:
  - 0012 numerics preset (B2)
  - 0013 reward v5 (B3)
  - 0014 polygonizer exception (B1)
  - 0015 docs site (DOC3)
  - 0016 CUDA polygonizer in production (D9)
  - 0017 layers and import-linter (E0)
  - 0018 typed run spec (E6)
  - 0019 TRL adapter (E5)
  - 0020 verifier API (E7)
  - 0021 library keep list (CLEAN14)
- Replace every phase reference with the allocated number.

### 8. [major] [00_front, A, B, C, D] The phases use different pytest markers, so some GPU tests never run
Problem.
- Section 4, A (`tests/test_a_gpu.py`, `@pytest.mark.gpu`), B (sketch `@pytest.mark.gpu`) and D8 use the marker `gpu`.
- C uses `cuda`: C.6 "no test needs a GPU unless marked `cuda`", C.9.0, C.9.4 "marker `cuda`" and C.9.5 "`cuda` + `slow`".
- `make gpu-check` runs `pytest -m gpu` (D8), so the C GPU tests never run in that chain.
- D8 registers the markers "with D2", but A and B use `gpu` earlier.
- C.9.0 says unregistered markers fail under `--strict-markers`. But `pyproject.toml` has no `addopts` and no `--strict-markers` [code: pyproject.toml:65-71], so they only warn.

Exact fix.
- Use one marker, `gpu`, everywhere. Replace `cuda` with `gpu` in C.6, C.9.0, C.9.4 and C.9.5.
- Move the marker registration (`gpu`, `sanitizer`) to the first task that uses it: A3/A6 `tests/test_a_gpu.py`.
- Add `addopts = "--strict-markers"` in that same task. Correct the C.9.0 sentence to match.

### 9. [major] [DOC, E, CLEAN] The test-count gates use one 2026-11-08 baseline, but each phase starts later
Problem.
- E-G1 and CLEAN-G2 say "1447 passed". DOC-G7 says "1,484 collected + new docs tests (no test removed)".
- A, B, C and D each add tests before E, DOC13 and CLEAN Wave 2-3 run. CLEAN and E also remove tests.
- So the literal numbers fail, and "no test removed" in DOC-G7 conflicts with CLEAN's archive lists.

Exact fix.
- Each phase records its own baseline when it starts: `pytest --collect-only -q | tail -1` and the `make test` summary line, in `runs/<phase>/baseline.txt`.
- Each gate then reads: "passed >= baseline at phase start - tests listed in this phase's archive READMEs + tests added".
- In DOC-G7, remove "no test removed".

### 10. [major] [A, B, C, E] Some behaviours are still open, and their tests are written to accept either outcome
Problem.
- A.9.8 lists 13 behaviours that are still undecided. Examples:
  - "`sel_vote` ... mixed presence raises or gives null";
  - "`summarize` with NaN raises or gives a NaN mean";
  - "clamp or reject" for window start < 2;
  - "skip with a count, or raise".
- B.9.12 has 9 more. Two of them change data gates:
  - #1: the D12 boundary differs by one token between spec 3.2 and section 6.1;
  - #6: may a window extend beyond the slice (padded windows)? This changes D4, D5 and D14.
- C.9.9 and E.9.10 add more.
- A test that takes "raises or returns null" pins nothing (tdd: an assertion that cannot fail). Under the data-first rule, B #6 is also a data-quality decision, not a test detail.

Exact fix.
- For each item, write one rule in the task text and make the test assert only that rule. Default: the first option listed, which is "raise" in most items.
- For B #1, use section 6.1 (`max + 1 < limit - 100`) and fix spec 3.2 to match.
- For B #6: "a scheduled window never extends beyond the slice; padded windows are dropped by the pool and counted in D14". Add a gate row: padded share = 0.
- List the decisions the owner must confirm in the 99_back owner list.

### 11. [major] [B] Data gate D3 (near duplicates) has no unit, and it lets known cross-split duplicates pass
Problem.
- D3 says "Near duplicates (cos >= 0.98) cross-split <= within-train baseline rate (4 / 7,341)", and today "1 test->train".
- The numerator and the denominator are not defined: pairs or items, and query items of which split. So the gate cannot be computed the same way twice.
- Under "dataset quality is the highest priority", a gate that passes while a test image has a train near-duplicate is too weak.
- The B1.3 planted tests cover D1, D2, D4, D5, D6 and D8, but not D3 (only U24 in B.9.10).

Exact fix.
- Define D3 as: "number of cross-split pairs (train item, val or test item) with cosine >= 0.98 on centred 32x32 HU crops. Every pair is listed in the report. Limit: 0 unreviewed pairs. A reviewed pair is resolved by dropping the train item from the pool, or by the owner's written OK."
- Report the within-train rate (pairs per train item) as information only.
- Add a planted test to B1.3: a train crop copied into val makes D3 fail.

### 12. [major] [00_front, A, B] The B gates do not name the K-aggregation, and the front puts Gate A after B3
Problem.
- Gate B1 and Gate B3 compare "val-select IoU" (R5 0.599). They do not say whether this is the mean over the K = 8 answers or the method that Gate A picks (medoid or vote).
- Gate A is decided in the A3 rental (G1, before G2). But the section 1 diagram draws "Gate A ... decides the answer-choice default" after the B3 run.
- If Gate A switches the default to medoid, the B1-vs-R5 pairs mix two aggregations.

Exact fix.
- In B.6.2 add: "All run gates use the per-item mean of IoU_tc over K = 8 answers (`sel_mean_of_k`), whatever Gate A picks. The Gate A method is reported as a second line."
- Move the "Gate A" box in the section 1 diagram to the Mac/G1 column, before the "owner approves the spend" box.

### 13. [major] [99_back] The owner decision list is missing open decisions, and its labels clash with task ids
Problem.
- These phase decisions are not in the 99_back list:
  - D.7 Q1 (may hole-ring order change; it changes hashes), Q2 (GPU CI runner), Q3 (`min_area` required in `build_target`) and Q4 (D ceiling $1.40);
  - C8 (beta 0 / KL weight 0), named as an owner decision in C.1 and C8;
  - A.8's "a failed legacy check (A-G2)" approval.
- Items 14-19 carry the labels "E1".."E6". These clash with tasks E1-E8. For example, back "E3" is the HTTP face, but task E3 is "lower layers".
- C.8.2 also says "Direction.md §8 shows $0.60 left under the $17.50 limit". This contradicts the binding $3 budget, and the list does not resolve it.

Exact fix.
- Add rows for D-Q1..Q4, C8 beta 0 and A-G2 acceptance, each with its default.
- Rename the E rows "E-Q1..E-Q6", matching the numbering in E.8.
- Delete the "$17.50" sentence in C.8.2, or add it as an explicit owner question: "which budget binds".

### 14. [major] [00_front, all phases] Section 0 says "one word for each thing", but the phases use synonyms
Problem.
- Section 0 says "One word for each thing. The phase sections use the same words." This is not true in the text:
  - "window/windows" 182 times and "frame", for the defined word "view";
  - "rater/raters" 92 times and "voters", for "doctor";
  - "exam", for "test" (B.1: "The exam is shifted");
  - "completion/completions" 32 times and "sample(s)", for "answer";
  - "crop", for "view" or for the 32x32 patch.
- Section 0 itself defines both "image" and "view" as "one window".
- "K" means 8 answers (final_eval), the number of window draws (B1.4 `--draws K`) and a cache count.

Exact fix.
- In section 0, define "image" as the pixels the model sees (one view), and keep "view" as the window.
- Add rows only where code names force a second word: "window = view (B1 code name)", "rater = doctor (code name)", "completion = answer (TRL name)".
- Replace "exam" with "test". Replace "K" in B1.4 with "D" (draws).

### 15. [major] [00_front, A, B, C, D, E] Price, stop rule and step time each have several values
Problem.
- Planning rate: $0.60/h (section 6), $0.641/h (B.8.1), $0.593/h (C.8.2) and $0.544-0.628/h (D.8). The costs are not comparable, and the $3 arithmetic changes with the rate.
- Stop rule: "no first sync in 50 min" (section 6.4 rule 3) and "35 min" (D.7 R13).
- R5 step time: 23.43 s (C), 23.65 s (B, section 6.3) and 23.6 s (section 3). `shared_rl.json` `timing_s.per_step` = 23.65 [measured]. 23.43 is the mean over steps 3-15.

Exact fix.
- Set one planning rate in section 6, for example $0.641/h, the highest measured. Make every phase table cite "section 6 rate".
- Set one stop rule. 35 min is safer, and the measured good launches took 18-31 min. Delete the other.
- Define in section 0: "R5 step = 23.65 s (all steps); 23.43 s = steps 3-15 mean, used only for the C speed gates".

### 16. [major] [A, B, D, CLEAN] The rule "run Python that reads data/ with -I" has no way to run in the tests
Problem.
- D.9.0 says "`data` = needs `data/` ... Python runs with `-I`". CLEAN L4 says "run with `-I`".
- The data-backed tests (D R01, R07, R08; B I-tests on the real store; A [box] "or on the Mac with `data/`") run under `make test`.
- `make test` is `$(PY) -m pytest -q -m "not slow"` with no `-I` [code: Makefile]. The project CLAUDE.md requires `-I` for any Python that reads `data/`.

Exact fix.
- Add a Makefile target: `test-data: $(PY) -I -m pytest -q -m data`. Change `test` to `-m "not slow and not data"`.
- State in each Tests section that `data` tests run only through `make test-data`.
- Add a repo-hygiene test: every test that opens `data/` carries the `data` marker.

### 17. [major] [all phases] Many sentences are too long for ASD-STE100
Problem.
- A check that joins wrapped lines finds 803 of about 4,965 prose sentences over 25 words. ASD-STE100 allows at most 20 words in procedures and 25 in descriptions.
- Front examples: line 262 (45 words), 265 ("If the gate fires ...", 60 words), 323 (46), 367 (57).
- Many phase table cells are paragraphs, for example A.R rows 2 and 5, the C.1 corrections and E R-1.
- An engineer must re-read these to find the action.

Exact fix.
- Split each sentence over 25 words into one action per sentence. Begin with the front (lines 240-380) and the task "Steps" lists.
- Move the reasons into a following "Why:" line.
- Run the same count after the edit and record it in DOC.6 as a non-blocking metric.

---

## Minor

### 18. [minor] [00_front, 99_back] One model has two names
Section 3 says `MBZUAI/MedMO-4B`. Back decision 13 says `MedMO-4B-Next`. Both exist on HF (4,437,815,808 params each).
Fix: use one id in both places. Prefer `MBZUAI/MedMO-4B`, the one section 3 checked.

### 19. [minor] [00_front] The section 6.1 step 3 commands do not run as written
- `./build_cpu.sh` needs the argument `<REPO/cuda>`. It uses `set -u`, so it exits without it [code: K1/build_cpu.sh:5].
- K2's script is `./run_cpu.sh`, not `k2test 20261108 1`.
- K3 `build/test_cpu` needs `make` first.

Fix: `(cd docs/implementation2.0/kernels/K1 && ./build_cpu.sh "$REPO/cuda")`, `(cd .../K2 && ./run_cpu.sh)`, `(cd .../K3 && make cpu && ./build/test_cpu)`. Use the target names from K3/Makefile.

### 20. [minor] [00_front, A, E] The NVTX module has two homes, and the tree lists two stores
- A.8 says "E moves `nvtx.py` to `infra/`". E.8 and section 2.2 say `core/marks.py`.
- The section 2.3 tree lists both `core/marks.py` and `infra/ ... nvtx (marks)`.
- It also lists both `verify/store.py` and `data/store.py`, and does not say what each holds.

Fix: use `core/marks.py` everywhere and delete the `infra` entry. Add a one-line purpose to `verify/store.py`, or remove it.

### 21. [minor] [A, E, CLEAN] Two phases each create the same test files
- A creates `tests/test_repo_hygiene.py` (S01-S11) as "new". CLEAN12 also calls it "new".
- A7 creates `tests/test_trl_adapter.py`, but E5 says "Add `tests/rl/test_trl_adapter.py`".
- E-G5 asks for "3 negative tests". A-G7 asks for "5 negative cases".

Fix: CLEAN12 "extend `tests/test_repo_hygiene.py`"; E5 "`git mv tests/test_trl_adapter.py tests/rl/` (owner runs)". Make E-G5 say ">= 5 negative cases (A-G7)".

### 22. [minor] [B] One B sketch can pass without testing anything, and its calls use two signatures
- `test_schedule_counts_and_rules` returns early on any `ValueError` that names a kind. If `schedule` always raises, all 40 examples pass.
- `test_unfillable_cell...` calls `schedule(ws, target, 4, 32)`. The other test calls it with 7 positional arguments.

Fix: add `hypothesis.event("raised")` and a separate non-hypothesis test that a feasible pool returns rows. Use keyword arguments `shares=, seed=, dw=` everywhere, matching the B1.2 signature.

### 23. [minor] [A] A-G1 needs bootstrap CIs equal to 16 digits, but numpy is not pinned
The bit-exact bootstrap CIs (seed 0) depend on the numpy `default_rng` stream and summation order. E.7 notes "numpy/scipy unpinned".
Fix: record `numpy.__version__` beside the golden. Make the exact check bind only on the recorded version and platform (as E-G3 does), else use `abs <= 1e-12`.

### 24. [minor] [A, C] A and C disagree on whether vLLM sampling repeats
A.7 says "vLLM sampling is not bit-reproducible run to run [estimate]". C G-EQ6 requires equal step-1 rewards across two arms. R5 measured "step 1 identical" with the same seed and host.
Fix: in A.7, cite the R5 measurement ("same seed and host: step 1 identical [measured: runs/r5]"). Keep the A-G3 ±0.01 rule for different hosts.

### 25. [minor] [00_front] The section 2.4 and 2.5 diagrams show optional work as the default
- Section 2.4 shows "GPU patchify (K3a)" as the default. Section 5 says C1 uses torch ops first, and K3a is only a fallback.
- Section 2.5 shows vLLM receiving "vision features once". C6b says feeding vLLM the trainer's features is an investigation, not in G-SP1, and today vLLM gets 256 PIL images.

Fix: label the box "GPU patchify (torch ops; K3a only if C1 misses its gate)". Label the vLLM arrow "views (uint8); feature handoff = C6b, not planned".

### 26. [minor] [A, B] A3 rebuilds the test views before B1.6 freezes them, and the doc gives no order
A3 step 1 rebuilds `data/eval_views` (test) and removes two roots. D15 / B1.6 later freezes the test shifted view hash. Nothing ties the A3 answers to the frozen hash.
Fix: in B1.6, "freeze the hash of the view root used by A3. `final_eval.json.view_manifest_sha256` from A3 must equal `frozen_views.json`." Add this to D15.

### 27. [minor] [00_front, D] The section 6.2 toolchain leaves out `nvtx`, which the D build needs
Section 6.2 installs `cuda-toolkit[nvcc,cccl]==13.0.3`. D1 makes the CUDA build fail with `FATAL_ERROR` when the `nvidia.nvtx` package is missing, and installs `[nvcc,cccl,nvtx]`.
Fix: write `cuda-toolkit[nvcc,cccl,nvtx]==13.0.3` in section 6.2.

### 28. [minor] [00_front, all phases] The doc uses tags that section 0 does not define
Section 0 defines four tags. The phases also use `[web: ...]`, `[unverified]`, `[measured now]`, `[measured: this session]`, `[measured: review]` and `[measured/estimate: ...]`.
Fix: add `[web: source]` and `[unverified]` to section 0. Change "measured now / this session / review" to `[measured: <file or command>, 2026-11-08]`. Split `[measured/estimate]` into two tags.
