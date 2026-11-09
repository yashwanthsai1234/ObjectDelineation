# Final review round 2: QA lead and technical writer (ASD-STE100)

Target: REPO/docs/implementation2.0/Implementation.md (13,065 lines, round 2), the kernel folders in
REPO/docs/implementation2.0/kernels/, and the code each claim names. Earlier review: `final_r1_qa.md` (28 items). Items
that round 2 fixed are not repeated. An item that comes back here was fixed in one section but not in another.
No REPO file was edited. No git write. No rental. Python that read `data/` ran with `-I`.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No new catalog lookup was needed. The front already names `compute-sanitizer-diagnose-cuda-bugs` and the 3D skills as recommendations only. Instead I checked the NVIDIA toolchain claims on PyPI and on the CUDA repo: the `cuda-toolkit 13.0.3` extras include `nvtx` and `sanitizer`; `nvidia-cuda-sanitizer-api 13.4.92` exists; the pinned nsys `.deb` URL answers HTTP 200 with 458,836,100 bytes, as A6 states. No item came from these checks. |
| engineering:code-review | I checked each claim against code: `Makefile` (`test`, `test-all`, `cov`), `pyproject.toml:65-71` markers, `policy.py:308-316` (launch checks one preflight config), `preflight.py:78-91`, `run_box.sh:24`, the K1/K2/K3 scripts. Results: items 7, 16, 18. |
| engineering:debug | For each gate I asked three questions: what input does it read, which rental or task makes that input, and can it fail? This found items 1, 2, 11 and 14. |
| engineering:system-design | I put the rentals V1-V3, the time line and the phase dependencies on one line. This found items 1, 2, 13 and 14. |
| engineering:architecture | I checked every ADR number against the front section 4 table. B and C still use the old numbers (item 4). |
| engineering:documentation | I checked one word per thing, one id per gate, and links to files that DOC1 step 0 copies. This found items 5, 12, 20, 21 and 22. |
| huggingface-skills:huggingface-best | The model ids and sizes in sections 3 and 99_back now agree (`MBZUAI/MedMO-4B`, 4.438 B). No item. |
| engineering:testing-strategy | I matched each Tests section to its gates: markers, Makefile targets, coverage gate, data isolation. This found items 7, 8, 11 and 15. |
| mattpocock-skills:tdd | I looked for tests that cannot fail (either-outcome assertions) and tests that are red from the start. This found items 9 and 15. |
| Web (WebFetch, WebSearch) | PyPI JSON of `cuda-toolkit 13.0.3` and `nvidia-cuda-sanitizer-api`; HEAD of the nsys `.deb` on `developer.download.nvidia.com/compute/cuda/repos/debian13`; a search for Nsight Systems 2026.5.x. |

Measured during this review (2026-11-08):
- R5 `items.jsonl`: val 1,084 rows, mean IoU_tc 0.61819, miss share 0.12638; test 2,200 rows, 0.64477, 0.10818.
  `items.jsonl` has no `group_id` column.
- `data/store/lidc/{val,test}.parquet`: 78 and 158 distinct `group_id` (matches A1 and A-G3).
- `runs/r5/shared_rl.json` `timing_s.per_step` = 23.653.
- Prose sentences over 25 words (code blocks and tables excluded; same rule as round 1): 993 of 5,203.

Severity: blocking = wrong, unsafe or breaks an owner rule; major = missing or vague, so an engineer cannot
implement it; minor = clarity or consistency.
Count: **3 blocking, 13 major, 10 minor**.

---

## Blocking

### 1. [blocking] [00_front, A, C, 99_back] Three sections give three different plans for the A3 slot and the A6 nsys capture

Problem.
- Front section 6.2, V1 item 4: A3 runs after the B1 run and holds "the A6 nsys window (6 steps, capture steps 4-5)".
  The 28 min of item 4 include that window. The front review log (line 577) says "A3 keeps its own 6-step window".
- Phase A says the opposite:
  - A3 is "the first item of the shared rental G1a, right after setup" (A.4, A3 header, A3 step 3, A.3.1 diagram).
  - `a3_measure.yaml` "has no nsys step". The A6 baseline is "the C0 baseline arm capture of the same rental".
  - A-G6 and G06 read steps 16-17 of the C0 baseline arm "of G1a".
- The rental "G1a" does not exist. Section 0 names only V1, V2, V3. The C0 baseline arm is V2 item 3, and V2 runs
  only with the money left (front 6.3).
- C0 step 3 still says "If the front picks QA option (a) for A-G6". Nobody has picked.
- 99_back open question 8 says "in the A3 nsys window of V1", which A no longer plans.
- Result: A-G6 has no box data in V1. The V1 minutes and cap include a capture that no plan file runs. A says A3
  goes first; the front says the B1 run goes first.

Exact fix.
- Choose option (a): A-G6 reads the C0 baseline arm in V2 item 3 (steps 16-17, overhead against 13-15).
- Front 6.2 item 4: delete "plus the A6 nsys window (6 steps, capture steps 4-5)". Set the minutes to
  17.5 sampling + 5 load = about 23 min [estimate]. Recompute the V1 total and cap.
- Front review log line 577: replace the text with "A6 capture = C0 baseline arm in V2 (option a)".
- Phase A: replace every "G1a" with "V1". Write "A3 is V1 item 4, after the B1 run and its val eval". Places:
  A.3.1 diagram, A.4 rental paragraph, A3 header and step 3, the `a3_measure.yaml` header comments, A.6 "All" row,
  A.7 rows, A.8 GPU list, the A.9.6 heading, G06 and the A.9.7 "All" row.
- A.3.1 diagram: move the "C0 baseline arm" box to a V2 column.
- A-G6: split it. The Mac half passes in Phase A. The box half is "open until V2 item 3".
- C0 step 3: delete "If the front picks QA option (a)". Write "the baseline arm's steps 16-17 are A6's capture".
- 99_back open question 8: change "the A3 nsys window of V1" to "the C0 baseline arm of V2 item 3".

### 2. [blocking] [00_front, A, B] V1 does not produce the inputs Gate B1 needs, and with them V1 is over its cap

Problem.
- Gate B1 has a hard guard row: "val IoU on the `edge` view: CI lower bound of (B1 - R5) > -0.02 (R5 value from A3)".
  B.8.4 lists "R5 evals on the val edge view" and "R5 val eval with a second sampling seed" as A3 outputs.
- A3 makes neither:
  - `a3_measure.yaml` evaluates only `lidc-val-shifted` and `lidc-test-shifted`.
  - A3 step 7 says that a second val sampling "is not planned". It uses a K = 4 half split instead.
  - B.6.2 "Noise floor" still says "Phase A3 evaluates R5 on val twice with two sampling seeds".
- V1 item 3 has one val eval of the B1 model (6.9 min). B-S1 in B.8.1 has two (shifted and edge) and B2.3b (15 min).
  B2.3b ("end of B-S1") has no slot in V1, but B2.4 step 4 reads its result.
- With the two edge evals (2 x 6.9 min) and B2.3b (15 min), V1 at 27.5 s is 202 + 28.8 = 231 min = $2.47
  [estimate: front 6 inputs]. That is above the V1 cap of $2.20, and it leaves $0.28 for V2.
- The A3 bundle ships only `lidc-test-shifted`. The edge view (B1.6) is not in any bundle.

Exact fix (the owner picks one; write the choice in front 6.2, A3, B.6.2, B.8.1 and 99_back decision 1):
- (a) Keep the guard row hard:
  - add a step to `a3_measure.yaml` for R5 on `{data}/eval_views/lidc-val-edge`;
  - add the B1 model on `lidc-val-edge` to V1 item 3 (2 x 6.9 min);
  - add B2.3b as V1 item 5 with the same skip rule as A3;
  - ship `lidc-val-edge` in the V1 bundle;
  - recompute the V1 total and cap, and the V2 money left.
- (b) Make the edge row "report only" in Gate B1, and move B2.3b to V3.
- In both cases:
  - B.6.2 "Noise floor": replace the two-seed sentence with "A3 step 7: K = 4 half split of the saved val answers".
  - B.8.1: make the B-S1 row list exactly the V1 items.
  - B.8.4: delete the A3 outputs that A3 does not make.

### 3. [blocking] [A, B, E] A and B use different decision intervals for the same gates, and no task builds the B interval

Problem.
- A1 step 5 says: "The percentile interval is the decision value on the whole val split (78 patients) and on test
  (158). A gate on fewer than 50 patients (for example a val half) uses a cluster t interval." A-G3 follows this rule.
- B.6.2 says: "The gate decision uses the BCa interval; the percentile interval is reported beside it." B.8.4 lists
  "BCa for the decision" as an A1 output.
- E-G10, E4 (line 8565) and GP5 also say BCa.
- A1 builds no BCa interval. Its API is `signal.ci`, `signal._ci` and `signal.cluster_t`.
- Gate B3 runs on val-select halves. A half has about 39 patients [estimate: 78 / 2]. So A says cluster t and
  B says BCa for the same gate.
- One gate with two decision rules can pass under one rule and fail under the other.

Exact fix.
- Use A1's rule everywhere: percentile interval at 50 patients or more; `signal.cluster_t` below 50.
- B.6.2: replace "The gate decision uses the BCa interval; the percentile interval is reported beside it" with
  "The decision interval is the A1 rule (percentile; cluster t under 50 patients)".
- Change B.8.4 row 1, B.R item 2, E-G10, E4 line 8565, GP5 and E review row "data 2" the same way.
- If the owner wants BCa: add `signal.bca(x, D)` to A1 with a unit test against `scipy.stats.bootstrap(method="BCa")`
  on unclustered data. Then change A1 step 5 and A-G3 to BCa as well.

---

## Major

### 4. [major] [B, C, DOC] B and C still use the old ADR numbers, so two phases write ADR 0012

Problem.
- Front section 4 sets: 0012 docs site (DOC10), 0013 numerics preset (B2), 0014 reward v5 (B3), 0015 scheduled rows
  and the numpy-polygonizer exception (B1).
- B still writes ADR-0012 = numerics, ADR-0013 = reward v5, ADR-0014 = polygonizer exception. There are about 50
  places, for example lines 2747, 2763, 2934, 2966, 2975, 3031, 3383, 3440-3464, 3494, 3508, 3601-3644, 3658, 3754,
  3761, 3824, 3847-3884, 3984, 4560, 4645-4691.
- C does the same at lines 5300, 5688 and 6573 ("ADR-0012 option F").
- DOC (line 9873, DOC-G11) writes `0012-docs-site-with-zensical.md`. B would write a second 0012 file. DOC-G11
  ("0 duplicate ADR numbers") then fails.

Exact fix.
- In B and C only, replace the numbers in this order, highest first, so that no number is changed twice:
  ADR-0014 -> ADR-0015, then ADR-0013 -> ADR-0014, then ADR-0012 -> ADR-0013.
- Do not touch DOC.
- Add a DOC4 check: every `ADR-00NN` in `docs/implementation2.0/Implementation.md` matches the title in the
  front section 4 table.

### 5. [major] [00_front, B, C] The front names the B data gates "B-DG1..16", but B has 21 gates without the prefix

Problem.
- The front says "B data gates B-DG1..16" at lines 98, 127, 251 and 455.
- B.6.1 has 21 gates: DG1..DG18 plus DG1b, DG9b and DG14b. U19 and E04 assert "21 gates".
- Three hard gates are outside "1..16":
  - DG1b (rows-root leakage);
  - DG17 (contour-to-slice check);
  - DG18 (padded windows).
- An engineer who follows front step 2 of section 6.1 can skip them.
- The prefix rule of section 0 ("every gate ... carries its phase prefix") is not used inside B (142 `DG` without
  `B-`) or inside C (92 `G-EQ`, 64 `G-SP` without `C-`).

Exact fix.
- Front lines 98, 127, 251, 455 and FR11: write "B data gates B-DG1..B-DG18 (21 gates, incl. 1b, 9b, 14b; B.6.1)".
- Section 0 prefix rule: add "Inside its own phase a gate may drop the prefix. A reference from another phase
  always uses it."
- Then fix the cross-phase references, for example A3 step 6 "C's G-EQ2 and G-EQ3" -> "C-G-EQ2 and C-G-EQ3".

### 6. [major] [B] B1.6 compares A3's manifest SHA with the frozen content hash, so DG15 fails by construction

Problem.
- B1.6 says "Freeze a content hash, not the manifest SHA", because the manifest holds Python and package versions.
- The same task then says "`final_eval.json` `view_manifest_sha256` from A3 must equal the frozen entry".
- A3 writes both keys. `view_manifest_sha256` is the SHA-256 of `manifest.json`, and that never equals a content hash.
- A3 step 1 and G01 already use `view_content_sha256`.

Exact fix.
- B1.6 "One hash from A3 to the end": replace `view_manifest_sha256` with `view_content_sha256`.
- DG15: add "(key `view_content_sha256`)".
- Add to I05: "a frozen entry compared with `view_manifest_sha256` is a test failure".

### 7. [major] [A, D] A changes `make test` for the `data` marker, but `make cov` and `make test-all` still run data tests without `-I`

Problem.
- A.9 changes `test` to `-m "not slow and not data"` and adds `test-data` with `-I`.
- `cov` is `$(PY) -m pytest -q -m "not slow" --cov=objdel ... --cov-fail-under=90` [code: Makefile:25-26].
  `test-all` is `$(PY) -m pytest -q` [code: Makefile:22-23].
- Both still collect the `data` tests (about 10 files today use the marker [measured: grep of tests/]).
  They read `data/` without `-I`, which breaks the repo rule.
- D-G12 gates coverage at 90 % with `make cov`. If `cov` deselects `data`, the coverage number changes, and nobody
  has measured the new value.
- S03 checks only the text of the `test` target.

Exact fix.
- A.9: also change `cov` to `-m "not slow and not data"` and `test-all` to `-m "not data"`.
- Add `test-all-data: $(PY) -I -m pytest -q -m data`.
- Before D-G12 binds, run `make cov` once with the new marker set. Record the percentage in A's report. If it is
  below 90 %, D-G12 uses that value as its floor.
- S03: check the `test`, `test-all` and `cov` lines.

### 8. [major] [E, 00_front] E uses a marker `box` that no task registers; with `--strict-markers` the suite stops at collection

Problem.
- E.9 conventions define "`box` = runs only at the next owner-approved rental". GP3, GP4 and GP5 carry it.
- A turns on `--strict-markers`. A registers only `gpu` and `sanitizer` (`data` and `slow` exist).
- An unregistered marker is then a collection error. `make test` fails from E on.
- The front Tests row lists only `gpu` and `sanitizer`. DOC adds `docs_build` and `post_commit`. No one place
  lists the full set.

Exact fix.
- E: mark GP3-GP5 `gpu` and `slow`, with `pytest.skip("box only: next owner-approved rental")` when
  `gc.build() != "cuda"`. Or, if `box` stays, E0 registers it in `pyproject.toml` with a `conftest.py` skip.
- Front section 4 Tests row: list every marker and its owner: `slow`, `data` (exist), `gpu`, `sanitizer` (A),
  `docs_build`, `post_commit` (DOC4), `box` (E0, if kept).

### 9. [major] [B, C] Four tests still accept either outcome

Problem. Round 2 fixed A.9.8 and wrote one rule per item in B.9.12 and C.9.9. Some test rows were not changed:
- B U12: "empty list raises or returns a pinned value". B.9.12 item 5 says `cell_shares([])` raises.
- B U32:
  - "group 4 with length 12 honoured or rejected consistently";
  - "empty arrays pinned (empty result or error)". B.9.12 item 5 says the result is an empty array;
  - "NaN in `iou_tc` raises or propagates as pinned".
- C U-CF2: "an existing non-empty directory raises or appends with distinct names".
- Such a test cannot fail (tdd: an assertion that takes both outcomes).

Exact fix.
- U12: "empty list raises `ValueError` (B.9.12 item 5)".
- U32:
  - "group 4 with length 12 raises `ValueError` when `spec.group == 8`";
  - "empty arrays return an empty float64 array";
  - "NaN in `iou_tc` raises `ValueError` naming the index".
- U-CF2: "an existing non-empty directory raises `FileExistsError`; nothing is written".
- Add the three new rules to B.9.12 item 5 and C.9.9.

### 10. [major] [C] G-SP1 has no outcome between the pass value and the fail value, and its threshold moves during the rental

Problem.
- G-SP1 passes at <= 20.0 s (mode A) or <= 20.7 s (mode B). It fails at > 21.2 s or > 21.9 s.
- A step time of 20.0-21.2 s (A) or 20.7-21.9 s (B) is neither a pass nor a fail. Gate C-adopt needs a pass. So
  nobody knows if B uses the C code at, for example, 20.6 s in mode A.
- G-SP1 is also "recomputed before the fast arm from the measured one-call push" (C5 step 5). This moves the
  gate after a measurement in the same rental. B.6.2 says "Rules are fixed before the run". The front fixes 20.0 s.

Exact fix.
- G-SP1: "pass <= 20.0 s (A) / <= 20.7 s (B). Anything above = not adopted; B keeps today's code. Above 21.2 s (A) /
  21.9 s (B) also means re-profile before more C work."
- Delete "Recomputed before the fast arm". Write "The recompute from C5 step 5 is reported beside the gate. It does
  not change the gate."
- Make the same change in the C.9.5 test text (line 6218).

### 11. [major] [A] G06 asks for an NVTX on/off step-time comparison on the box, but no arm runs with NVTX off

Problem.
- G06 says: "with `OBJDEL_NVTX=1` and no profiler `time/step_s` is within 3 % of unset (R5: 23.6 s/step)".
- The C0 arms all run with `OBJDEL_NVTX=1` (A.3.1 diagram, C0 step 3). No arm runs without it.
- R5 ran on another host and another day. Two hosts can differ by more than 3 %, so a compare against R5 can
  fail with no NVTX cost.
- The test has no data source in any rental.

Exact fix.
- Delete the 3 % box clause from G06 and A-G6.
- NVTX cost is covered by Q03 (Mac, 100k disabled ranges) and by the 16-17 against 13-15 overhead line in the
  C0 baseline arm. Keep those two.
- If the owner wants a box number: add a fourth C0 arm with `OBJDEL_NVTX=0` (17 steps, about 8 min = $0.09
  [estimate]) to V2 item 3, and write its minutes in the V2 table.

### 12. [major] [DOC, C, B] DOC1 step 0 does not copy all the evidence the plan cites, and C still cites the scratchpad

Problem.
- The DOC1 step 0 table has no row for `impl/c_rev/`. C cites these files as [measured]:
  `c_rev/render_bench.py`, `c_rev/pass2_text_bench.py` (C.1 facts, C2) and `c_rev/C.v1.md` (C header).
- The table skips all of `b1_out/` (12 MB). But DG3 "Today" cites `docs/implementation2.0/notes/b1_out/neardup.json`
  (529 bytes). B.1 cites `b1_out/geom.json` (9.1 MB).
- C still has 15 `S[P]/` paths (lines 4703-4708, 4755, 4796-4806, 5147, 5295, 5385, 5431, 5445, 5476, 5735-5736) and a
  header that says "`S[P]/...` names a session note. Those notes are not in the repo yet." This contradicts section 0.
- After step 0, DOC-G14 (check 14) still fails on C: no row maps `docs/implementation2.0/notes/c_rev/`.

Exact fix.
- DOC1 step 0 table:
  - add a row `impl/c_rev/*.py`, `impl/c_rev/C.v1.md` -> `docs/implementation2.0/notes/c_rev/`;
  - change the `b1_out/` skip to "copy `neardup.json` and `shacheck.json`; skip `geom.json` (9.1 MB) and
    `pix.json` (3.6 MB)".
- B.1 line 2791: cite geom.json as "not kept: 9.1 MB" (DOC1 README rule), or cite a small summary.
- C header: delete the `S[P]/` paragraph. Use the section 0 path rule.
- Replace each C `S[P]/` path with its `docs/implementation2.0/...` copy. Do this in the C phase edit, not only in DOC1.

### 13. [major] [D, 00_front] D.4 and D.8 do not follow the front rental order

Problem.
- D.8 says "R-D1 (D1) runs inside the first validation rental, before B1 (front section 6.2)". In the front,
  R-D1 is V2 item 4, after the B1 run. Decision 3 says "drop R-D1 from V2 if the money left is too small".
- D.4 says "Order: D1, then D2". D1 includes "one GPU run of K1-draft and K2-draft". Front time line step 5 lands
  D2-D8 code before V2, so before the D1 GPU run. If R-D1 is dropped, D1 never ends and D2 never starts.
- D.8 prices at $0.544-0.628/h. Section 6 says every phase cost table cites $0.641/h.
- D-G14 and D-G15 say "first B/C rental after D8". The front calls it V3.

Exact fix.
- D.4: "D1's Mac parts (engine counters, NVTX, `bench_checker.py`) gate D2. The R-D1 GPU numbers gate only speed
  claims (D-G7). D2-D8 use the CPU build as the oracle."
- D.8: R-D1 = V2 item 4 (optional, owner decision 3); R-D2 = V3 item 1. Price every row at $0.641/h. D-G13 becomes
  79 min x 0.641 / 60 = $0.84, the same number as today.
- D-G14, D-G15, D-G17, B04, B06: write "V3" for "first B/C rental after D8".

### 14. [major] [E, 00_front] The box parts of E-G10 and E-G13 need a rental after E, and the plan has none

Problem.
- E-G10 has a binding box part: prompt ids and SamplingParams equal, and TOST ±0.01 on 1,084 val items. E-G13 is a
  box smoke. E.8 says both "piggyback on the next owner-approved Phase B run".
- In the time line, E runs at step 7, after V3 (step 6). No rental comes after E.
- E.8 adds about 9 min ($0.10) to that run. No front rental table holds it. V1 leaves $0.04 under its cap.

Exact fix.
- Front section 6: add "V4 (after E; owner decision; outside the $3 unless raised): E-G10 box, E-G13, about
  15 min = $0.16 [estimate: E.8]".
- E.6: mark the box halves of E-G10 and E-G13 "open until V4". The E phase gate passes on the Mac halves.
- 99_back: add the V4 decision to "Money and GPU".

### 15. [major] [A, DOC, 00_front] Two A tests are red or have no run condition when A starts

Problem.
- S06 (no live file says "TRL 1.16") belongs to A7. It "passes only after DOC2" (A.9.4).
- The front time line step 1 runs only DOC1 before the phases. DOC.4 lets DOC2 run "before B1's preflight". Only one
  sentence inside DOC2 (line 9853) says "DOC1-DOC2 run before phase A edits code".
- If DOC2 is late, S06 makes `make test` red after A7. But A.4 says `make test` passes after each task.
- S08 "runs once DOC has written the files". It names no skip rule. If it is collected before DOC writes the
  pages, it fails.

Exact fix.
- Front time line step 1: "DOC1 and DOC2 run before any phase starts." Add the same order to DOC.4.
- S08: `pytest.skip("A9 pages not written yet")` when `docs/reference/metrics.md` does not hold the string
  "fair ceiling". Or move S08 to DOC's `tests/test_docs_content.py`.

### 16. [major] [00_front, A] The V1 launch names no config, and A3 preflights a config that V1 does not launch

Problem.
- `launch.sh` refuses a launch without a passing preflight stamp for the one config it is given
  [code: experiments/phase2/vast/policy.py:308-316; experiments/phase2/preflight.py:78-91].
- V1 runs two configs: the B1 run config (`b1_rows.yaml`, B1.5) and `a3_measure.yaml`.
- The A3 header comment preflights `r1_qwen3vl2b.yaml`.
- Front 6.2 says V1 "runs today's code". But V1 needs the A3 changes (`final_eval` write order, `CODE_SHA`), the A6
  NVTX code and the B1.5 `train_rl` changes (ordered rows, `rows_revision` check). The stamp hashes
  `objdel_and_script_sha256`, so all that code must be on the tree at preflight.

Exact fix.
- Front 6.1 step 4 and 6.2: "V1 launch config = `experiments/phase2/configs/b1_rows.yaml`. Run `make preflight` for
  it and for `a3_measure.yaml` on the same tree, after A3, A6 and B1.5 land and before DOC2/DOC5-type `src/` edits."
- Replace "today's code" with "the tagged pre-E tree with A3, A6 and B1.5; no C, D or kernel code".
- A3 header: replace `r1_qwen3vl2b.yaml` with "the V1 launch config (front 6.2)".

---

## Minor

### 17. [minor] [00_front, A] A3 reads the test split again, but section 0 says the test is read once

Section 0 says the test is "Read once, on the final candidate only." A3 re-evaluates R5 on all 2,200 test images
(A-G3 TOST). That is a repeat measurement of an already-read model, not a selection. But the text does not say so.
Fix: in section 0, add "Exception: A3 re-scores the R5 model on test to check equivalence (A-G3). No choice is made
from it."

### 18. [minor] [00_front] Section 6.1 step 3 builds K1 but does not run its tests, and `$REPO` is not defined

`K1/build_cpu.sh` only compiles and prints "build_cpu: OK" [code: K1/build_cpu.sh, last line]. The pass column asks
for "PASS". Section 0 says paths are relative, so `$REPO` is empty in the owner's shell. `set -u` then stops the
script with "unbound variable".
Fix: begin step 3 with `REPO=$(git rev-parse --show-toplevel)`. After the K1 build, add
`docs/implementation2.0/kernels/K1/build/cpu/test_k1_cpu 3000 41` and the `_o1` and `_o2` variants (K1 DESIGN.md:190).

### 19. [minor] [A, B, C] The longest train target is 1,182 tokens in A and C, and 1,181 + 1 in B

A (lines 709, 1440, A-G8, G05) and C (5384, U-RC1) say "max 1,182". B (DG12, B.8.4) says "max 1,181 + 1". The
B1 note counts "without end token" [measured: SP notes B1.md:90]. A1.md does not say which. If A8 counts the end
token, DG12's "+1" counts it twice.
Fix: in section 0, define "outline length = tokens without the end token". Write 1,181 in A, C and G05. Or keep
1,182 and label it "with end token" everywhere. Then DG12 drops the "+1".

### 20. [minor] [A, B] The test-source lists point at the repo `tests/` folder

A.9 cites `tests/A_tests_1.md` and `tests/A_tests_2.md`. B.9 cites `tests/B_tests_1.md` and `tests/B_tests_2.md`.
Under the section 0 path rule, these are files in the repo test folder. DOC1 step 0 copies them to
`docs/implementation2.0/tests/`.
Fix: write `docs/implementation2.0/tests/{A,B}_tests_{1,2}.md`, as C (after item 12), D, E and DOC do.

### 21. [minor] [A, D, E, CLEAN] "exam" is still used, against section 0

Section 0 says "Never called exam". The word remains at lines 791 ("exam path"), 1508 ("exam metrics"),
7127, 8033, 8807 ("exam metrics") and 10981 ("live exam harness").
Fix: replace "exam" with "test evaluation" (or "eval") at each place.

### 22. [minor] [C] C.9.0 still describes the marker and `-I` rules as proposals

C.9.0 says "QA review item 8 also proposes `--strict-markers`" and "QA review item 16 proposes a `make test-data`
target". A.9 adopts both. C.9.0 also says "an unregistered marker only warns", which stops being true when A lands.
Fix: "A registers `gpu` with `--strict-markers`; `data` tests run only through `make test-data` (A.9)."

### 23. [minor] [CLEAN, 99_back] CLEAN's open test questions are not in the owner list

CLEAN.9.1 lists Q-T1..Q-T4 as "not decided here". Back decision 32 lists the decided test rules of A.9.8, B.9.12,
C.9.9 and E.9.10 only. Q-T2 has a default ("add the filter"); the others do not.
Fix: write one rule each in CLEAN.9.1:
- Q-T1: a deeper `lab/` ships;
- Q-T2: add the `os.path.isfile` filter;
- Q-T3: exempt `data/`, `build/`, `.venv/`, `.claude/`;
- Q-T4: skip with a reason while `docs/` is untracked.
Add CLEAN.9.1 to decision 32.

### 24. [minor] [00_front, A] One R5 step time, three spellings

Section 0 defines the R5 step time as 23.65 s [measured: 23.653]. G06 says "R5: 23.6 s/step". The B-S2 option row
uses 23.6. C uses 23.43 only for its gates, as section 0 allows.
Fix: write 23.65 s in G06 and B.8.1.

### 25. [minor] [all phases] More sentences are over 25 words than in round 1

The round 1 rule (join wrapped lines, split on sentence ends, skip code and tables) now counts 993 of 5,203
prose sentences over 25 words [measured, 2026-11-08]. In round 1 it was 803 of 4,965. Most new long sentences are
in the review logs and the C/D gate cells. DOC.6 records this as a non-blocking metric, and the new text made it
worse.
Fix: split the long sentences in the task "Steps" lists and in the gate cells first. Example: G-SP4 is 96 words in
one cell. Use one rule per sentence, and put the reason on a "Why:" line.

### 26. [minor] [00_front, E] The B2/B3 path rule differs between the front and E

Front order rule 1 says "Every A-D task path is pre-E", and time line step 5 lands B2/B3 code before E. E.8 says "B2
and B3 land on pre-E paths if written before E; otherwise they land after E on post-E paths".
Fix: delete the "otherwise" clause in E.8, or add it to front rule 1 as a named exception.

---

Sources for the web checks: [cuda-toolkit 13.0.3 on PyPI](https://pypi.org/pypi/cuda-toolkit/13.0.3/json),
[nvidia-cuda-sanitizer-api on PyPI](https://pypi.org/pypi/nvidia-cuda-sanitizer-api/json),
[Nsight Systems 2026.5.1 announcement](https://forums.developer.nvidia.com/t/nsight-systems-2026-5-1-desktop-server-x86-64-arm-sbsa-arm64-native-and-x64-emulation-release-announcement/382869),
CUDA Debian 13 repo file `nsight-systems-2026.5.2_2026.5.2.323-1_amd64.deb` (HTTP 200, 458,836,100 bytes).
