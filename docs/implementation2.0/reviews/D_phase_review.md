# Review of Phase D plan (docs/implementation2.0/D.md and D.spec.md), 2026-11-08

Reviewer: fresh critical reviewer. REPO = $REPO.
Every finding below was checked against the code, a measurement file or a web page. Tags: [code: file:line],
[measured: file], [web: url], [estimate].

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | I checked the saved live catalog (docs/implementation2.0/notes/doc/nv_catalog.json). The only related skill is `compute-sanitizer-diagnose-cuda-bugs`. No catalog skill covers CCL, EDT or contour tracing. I installed nothing. This confirms D.md §0. I also found that racecheck and synccheck have never run in any repo chain (item 11). |
| engineering:code-review | I used it as the correctness lens for every line reference. Most line numbers in D.md are right. The review found these defects: the certificate legs are not independent (item 1), the rim check has no device path (item 2), and the duplicate rule is wrong (item 3). |
| engineering:debug | I traced the failure paths before any GPU spend: a `-n 0` flag with no xdist plugin (item 10), an NVTX include path (item 12), a missing nsys (item 13), and comparisons across two hosts (item 4). |
| engineering:system-design | I traced the data flow end to end. Result: `metrics()` runs inside the RL process (item 6). Expert caches are rebuilt on every load, so they are not offline (item 7). The pool-memory row of Architecture.md has no task (item 9). |
| engineering:architecture | I compared libraries with custom code. cuCIM EDT matches SciPy distances, so the stated reason to reject it is wrong (item 21). A disabled CI workflow is dead code (item 17). |
| engineering:documentation | I checked that each step is concrete. Vague inputs ("tests/goldens fixtures") do not exist (item 23). |
| huggingface-skills:huggingface-best | Loaded. No HF benchmark or model applies to Phase D kernels, so the review has no model finding. |
| engineering:testing-strategy | I found gaps. The chunk and scratch-split paths are not tested (item 18). The no-CPU-geometry gate is a grep that passes while cv2 still runs (item 8). Item 1 adds a mutation test for the certificate. |
| engineering:tech-debt | I found callers missing from the move lists (item 20) and the 3-D CPU path (item 25). |

Web checks done: compute-sanitizer docs (racecheck checks "only ... on-chip shared memory"
[web: docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html]), PyPI cuda-toolkit 13.0.3 extras (nvtx =
nvidia-nvtx 13.0.85, no Nsight Systems extra [web: pypi.org/project/cuda-toolkit/13.0.3]), cuCIM
`distance_transform_edt` (PBA+, `float64_distances=True` "to match SciPy", ties change only the indices
[web: docs.nvidia.com/cucim/latest/api/]), GitHub GPU runners (Tesla T4 only [web: docs.github.com/en/actions/reference/runners/larger-runners]),
TRUCO arXiv 2610.37297 (exists, v2 2026-11-01 [web: arxiv.org/abs/2610.37297]), and G-Contour (ICCAD 2025,
doi 10.1109/ICCAD66269.2025.11240932 [web]). The paper citations in D.md §2 that I checked are correct.

## Findings (most severe first)

1. **blocking**, D.md §4 D3 "certify.py ... Each leg is a separate kernel chain, so the legs stay independent", and
   D5 step 3 (D5b).
   Problem: the claim is false for SP-K1 option C. RingMark and RingInit enumerate the rings from the same run CCL
   (Merge, then Flatten roots), and the same CCL gives b0 and b1 [code: docs/implementation2.0/kernels/K1/k1poly.cpp:12-32, header
   "HoleKey ... RingMark ... RingInit ring g's start edge in closed form"]. So `(outer, holes) == (b0, b1)` is true by
   construction. The certificate (`certify.py:50-51`) loses one of its three independent checks. Today b0 and b1 come
   from cv2 Spaghetti and the rings come from numpy pointer jumping [code: checker/certify.py:17-18, 44-46]. D5b
   ("one CCL in the engine") would remove the last independent labeller.
   Fix: (a) take certify's b0/b1 leg from the scorer's pixel union-find chain (LabelInit, LabelMerge, LabelFlatten,
   Holes [code: gpucheck.cpp:1302-1392]). That chain is a different algorithm on bytes, not runs. Expose it through
   `certify_many`. (b) Drop D5b, or keep a second labeller. (c) Add a mutation test on the CPU build: compile Merge
   with foreground 4-connectivity under a test-only define, and require `certify(...).certified == False` on a
   diagonal-pair mask. (d) Rewrite the D3 sentence to say which legs share code.

2. **blocking**, D.md §4 D4 (`_build` change) and D.spec §3.2 `certify_polygons`.
   Problem: `_build` checks the RIGA rim on the HOST rasters: `kept("rim", ras["disc"] & ~ras["cup"])`
   [code: targets/build.py:239]. `certify_polygons` returns no pixels (spec C3), and it has no input that combines
   two shapes. So D4 as written cannot rebuild riga/train (3,282) or riga/val (348). G4 fails, or the plan keeps
   `to_mask` for RIGA, and then the claim "removes the per-rung to_mask host unpack" is false.
   Fix: add to `certify_polygons` an optional pair list `pairs: int32[P, 2]` and `pair_min_area: int64[P]`. For each
   pair (a, b), the device rasterizes a AND NOT b in words and returns the size-filtered signature of that mask. Do
   the same for the reference rim from `ref_bits`. Add a test: rim signatures equal `signature(disc & ~cup, ma["rim"])`
   on all RIGA ladder rungs of 50 items. Name the rim in D.spec B11.

3. **major**, D.md §4 D7 "Duplicates use mask-equality classes from the intersection matrix", and D.spec §3.2
   `eval_sets` `inter`.
   Problem: today the duplicate hash uses the RAW masks (`masks`, before the legal gate) plus the format and legal
   flags [code: eval/metrics.py:265-268]. `eval_sets` returns `inter` only for the DRAWN masks (an illegal sample
   becomes empty). Two illegal samples with different raw rasters are distinct today. In the new design they are
   duplicates, so the duplicate share changes and G9 fails. Also, the hash never leaves `metrics()` (it is used only
   at :335 and :348). The note about output-format.md is therefore not needed.
   Fix: add `inter_raw int64[n_sets, K, K]` (raw sample masks, per object name) to `eval_sets`. Define duplicates as
   equal (format_ok, legal) and, for every name, `inter_raw[i, j] == area_raw[i] == area_raw[j]`. Delete the
   output-format.md step from D7 and D9.

4. **major**, D.md §6 G7 and §8 (old bench in G1, new bench in G2).
   Problem: the old and new `check_masks` run on two different rentals, with different hosts, CPUs, drivers and
   clocks. A "new <= old on every workload" rule across hosts measures the host, not the code. The rule also has no
   noise margin.
   Fix: in the G2 rental, build the pre-D2 engine as a second module. Use `git archive <pre-D2 sha> cuda CMakeLists.txt`
   to get the old source; it is read-only, with no checkout. Build it with a module name override (`_gpucheck_old`).
   Run old and new interleaved: 10 repeats after 2 warm-ups, and report the median and IQR. Pass rule: new median
   <= 1.05 x old on every workload, and <= 0.9 x old on lesion512 b256 and on the 37,093 crops. Keep G1 for the
   SP-K1/SP-K2 standalone gates only.

5. **major**, D.md §6 G8 ("next RL run: mean time/rewards_s <= 0.0046 s").
   Problem: the R5 baseline is only 15 steps. Mean 0.004603 s, median 0.004509, range 0.003903-0.005383 s
   [measured: runs/r5/base/rl.json `log`]. The step-to-step spread is about ±15 %. The next run also trains on
   different rows (shifted and rescaled views), so the scoring work changes. The gate cannot pass or fail on the code.
   Fix: replace G8 with a replay benchmark in `bench_checker.py`. Score one fixed set of saved completions with their
   tables using the old module and the new module in the same rental (50 repeats each). Pass = new median <= 1.05 x
   old median, and the count of host waits per `score_polygons` stays 3 [code: gpucheck.cpp:1574, 1626, 1639; count
   them with an NVTX or CUPTI counter]. Keep the RL `time/rewards_s` value as information only.

6. **major**, D.md §1 item 7 ("Phase D is NOT an RL speed item") and §4 D7.
   Problem: `metrics()` runs INSIDE the RL process. The training callback calls `alarms_for` and then `metrics`
   [code: train/callbacks.py:63, 162]. After D7 that path calls `eval_sets`. `eval_sets` adds kernels, host waits and
   grow-only device arenas that are never freed [code: gpucheck.cpp:131-150 Block] next to vLLM and the trainer on one
   32 GB card. The plan gives no memory gate and no time gate for this.
   Fix: put `eval_sets` on the scorer lane, with a hard scratch cap (for example, chunk sets so that scratch stays
   <= 256 MiB [estimate; set the cap from the D1 measurement]). Add an NVTX range `eval_sets`. Add gate G14: in a
   6-step RL smoke run with the audit callback on, the audit time and the peak `torch.cuda.max_memory_reserved` plus
   the engine arena bytes are recorded, and the arena bytes stay <= the cap.

7. **major**, D.md §3.2 item 4 and §3.3 row "Expert-cache EDT ... (offline, once per item)".
   Problem: expert caches are not offline. `ExpertCache.__init__` runs `_band` (2 scipy EDTs), k cv2 signatures, a
   vote and `cv2.boundingRect` for every construction, and that includes `load` from npz (the npz stores only the
   raters) [code: reward/cache.py:222-238, 294]. D1.md C6 measured 1.0 / 3.6 / 7.5 ms per load at 128 / 256 / 400 px.
   With the D7 design, the device EDT copies a float32 map back to the host, and `add_tables` then copies it to the
   device again.
   Fix: add task D7b. A new binding `add_tables_from_raters(raters bits, k, H, W, min_area rule)` builds union, core,
   vote, the two EDTs, rater areas and rater signatures ON the device. It returns only the scalars that the host
   needs (perimeter, min_area, areas, signatures). Keep the host `ExpertCache` path as the oracle in `tests/`. Gate:
   the device tables are byte-equal to the host-built tables for 1,000 caches. Correct the "offline" text in §3.2 and
   §3.3.

8. **major**, D.md §3.3 table and §6 G11 ("No CPU geometry left").
   Problem: the checker's topology layer stays on cv2 in production: `instance_topology`, `size_filter` and
   `signature` (cv2 Spaghetti [code: checker/topology.py:77]). Production callers are `targets/build.py:151-155, 208`,
   `reward/cache.py:230`, `eval/probes.py:154, 308, 316, 571` and `experiments/phase2/build_views.py:57`. G11 greps only
   for the text `cv2.connectedComponents` in `targets/` and `metrics.py`. So G11 passes while cv2 still labels every
   target and every cache. This conflicts with the owner rule "CUDA is a must for the checker".
   Fix: (a) Add task D4b: `checker.signature`, `size_filter` and `instance_topology` call the device label chain
   (one batched binding, `signatures(masks, min_area)`). Merge the cv2 versions into the existing
   `tests/oracles/topology_oracle.py`. (b) Make G11 a runtime test, not a grep. Monkeypatch `cv2.connectedComponents*`,
   `scipy.ndimage.label` and `scipy.ndimage.distance_transform_edt` to raise. Then run `build_target` on one LIDC and
   one RIGA fixture, `ExpertCache.build`, and `metrics()`, and require no exception.

9. **major**, missing task (Architecture.md §5.2 row "Memory").
   Problem: Architecture.md §5.2 asks for pool memory and a batched table upload. D.md has no task for this.
   `add_tables` calls `cudaMalloc` and waits once per cache [code: gpucheck.cpp:1741-1748]. `drop_tables` calls
   `cudaFree` [code: :1774], which synchronises the whole device. `Block::get` frees and reallocates when it grows
   [code: :137-148]. These calls happen in the RL process.
   Fix: add task D5c. Use the stream-ordered allocator (`cudaMallocAsync` / `cudaFreeAsync` on the lane stream from a
   `cudaMemPool_t`, with a release threshold) for tables and arenas. Add a list form of `add_tables`: one H2D copy and
   one wait for a batch of caches. Gate: an nsys or CUPTI trace of 6 RL steps shows 0 `cudaFree` and 0 `cudaMalloc`
   calls after step 1.

10. **major**, D.md §4 D8 `Makefile` target `gpu-check` (`pytest -m gpu -n 0`) and D.spec §3.3.
    Problem: `-n` comes from pytest-xdist, and xdist is not a dependency (pyproject.toml:29 lists pytest, pytest-cov,
    hypothesis, ruff). On the box, pytest stops with "unrecognized arguments: -n". That wastes paid minutes.
    Fix: remove `-n 0`. Add `make gpu-check DRY=1`, which prints the chain and runs `pytest --collect-only -m gpu` on
    the Mac. Call it from `make preflight` with `experiments/phase2/configs/d_gpu_test.yaml`.

11. **major**, D.md §4 D8 sanitizer chain and §8 G2 time ("15-20 min [estimate]").
    Problem: racecheck and synccheck have never run in any chain [code: docs/testing.md:116]. The `sanitizer` marker
    set is not defined. racecheck is much slower than memcheck, so the 15-20 min figure has no basis. The SP-K1
    runbook caps each tool at 300 s with tiny inputs [code: K1/GPU_RUN.md:88-89]. D8 has no cap.
    Fix: list the `sanitizer`-marked tests by name. Use a small fixed corpus: <= 200 masks (including checkerboard
    and spiral), <= 500 polygons, 8 eval sets, 2 `certify_polygons` batches. Wrap each tool in `timeout 300`. Record
    each tool's time in `runs/d/sanitize.json`. Also run racecheck on the SP-K1 and SP-K2 binaries in G1 (their
    shared-memory ballots are where racecheck has coverage).

12. **major**, D.md §4 D1 Files (`#include <nvtx3/nvtx3.hpp>`) and common toolchain.
    Problem: the nvtx3 headers come from the `nvidia-nvtx` wheel, not from the nvcc wheel [web: PyPI cuda-toolkit
    13.0.3 extras]. D1 changes no include path in `CMakeLists.txt`. `cuda/docker_check.sh` (the only pre-rental
    compile check) does not install the extra either. The first failure would happen on a paid box.
    Fix: in `CMakeLists.txt`, add the wheel include directory, found from Python (`nvidia.nvtx` package path) with a
    clear error when it is missing. In `docker_check.sh` / `run_box.sh --build-only`, install `[nvcc,cccl,nvtx]`. Gate
    before G1: `cuda/docker_check.sh` passes with the NVTX ranges compiled. State how the C++ ranges relate to
    Phase A's `OBJDEL_NVTX` switch: they are always compiled in, and the Python side keeps its own switch.

13. **major**, D.md §4 D1 step 3 ("nsys ... if nsys exists") and D5 step 4 (path halving decided by the nsys profile).
    Problem: Nsight Systems is not in the cuda-toolkit pip extras [web: PyPI list], and the image is not shown to
    have it. If nsys is missing, the D5 step 4 decision has no input, and the "measure first" rule of D1 is half met.
    Fix: decide on per-stage CUDA-event timing, which needs no tool. Use SP-K1's `Timing` struct [code:
    K1/k1poly.h:54-56] and add the same event pairs around the scorer's label, count, filter and score stages in
    `bench_checker.py --stages`. Rewrite D5 step 4: "apply path halving only if LabelMerge is >= 20 % of score_words
    time on the replay benchmark". nsys stays optional.

14. **major**, D.md §4 D4 `_build` change ("per ladder level, one certify_polygons call for all names").
    Problem: `_build` advances the level of each name on its own, in an inner loop, and the rim check can raise all
    levels [code: targets/build.py:220-243]. "One call per level for all names" does not match this control flow. It
    also wastes the GPU batch.
    Fix: per name, build the parts of every remaining rung and certify them all in ONE call. The ladder has at most
    14 rungs: tau 3.0 down to 0 in 0.25 steps, plus lossless and exact rings [code: build.py:69-71]. Take the first
    rung that passes. This gives the same rung, so the text is unchanged (G4 checks it). The rim retry calls it again
    with the raised levels.

15. **major**, D.md §8 cost table and §6 G13.
    Problem: the budget ignores failed launches. The R5 first instance was destroyed after 15 min and cost $0.2385
    [measured: experiments/phase2/vast/ledger.csv, 2026-11-08T15:35:00Z]. Setup took 18, 20 and 31 min
    [measured: same file]. The $1.22 upper figure plus one failed launch passes the $1.25 cap.
    Fix: add a line for one failed-launch reserve ($0.25), and set the phase ceiling to $1.47, or ask the owner.
    Alternatively, keep $1.25 with a stop rule: abort the rental if the first sync takes more than 20 min. G13 counts
    every D ledger row, including destroyed ones. Option: run G1 on a plain CUDA image with no vLLM layers, because
    the G1 steps need only nvcc from pip.

16. **minor**, D.md §4 D8 `.github/workflows/gpu.yml` ("stays disabled until the owner picks a runner").
    Problem: a committed but disabled workflow is dead code (owner rule: remove dead code). A Vast lane driven from
    GitHub would need a Vast API key in repo secrets. It would also start a rental without the owner's go and
    without the preflight.
    Fix: do not create `gpu.yml` in Phase D. Keep `make gpu-check` as the manual Vast lane. Create the workflow only
    after the owner answers Q2.

17. **minor**, D.md §6 G7/G3 and D.spec B6: the chunk and scratch-split paths have no test.
    Problem: B6 says results do not depend on chunking. A chunk splits only above 2^25 pixels (more than 128 masks
    of 512x512), and the scratch halving only above 1 GiB [code: K1/k1poly.h:48-51]. No listed test reaches either
    path.
    Fix: expose `chunk_pixels` and `max_scratch_bytes` as test-only keyword arguments of the binding (default = the
    production values). On the CPU build, add a test that runs 40 mixed masks with `chunk_pixels=4096` and
    `max_scratch_bytes=65536` and requires output byte-equal to the one-chunk call.

18. **minor**, D.md §4 D6 step 1 (grid parked because "model output never reaches the grid").
    Problem: the reason covers only model output. Target-builder rungs go through the same `launch` validity (today
    through `to_mask`, and after D4 through `certify_polygons`). A rung can have up to 1,000 vertices per ring and
    10 rings [code: geometry/text.py:59; targets/build.py:225]. So shapes with >= 512 vertices do reach validity,
    though offline.
    Fix: in D1, log the vertex count per certified rung for 1,000 target rows (LIDC and RIGA). Report the EdgePairs
    time of the largest. Keep the grid parked only if those rungs are < 5 % of build time, and write that number in
    the ADR.

19. **minor**, D.md §4 D3 caller list and D4 "Update the 16 test calls".
    Problem: the lists are incomplete. `tests/test_breakit_geometry.py:17` imports `objdel.checker.bitquad` directly,
    so it breaks after the git mv. `experiments/phase2/lab/A6/cold_start.py:29` calls `build_target` with no
    `min_area`, so it breaks after D4. `experiments/phase1/bench/xarch_hashes.py:45` and
    `experiments/phase0/{lidc_census,riga_census,riga_outlines}.py` call `polygonize` or `certify` and must still run.
    Fix: add these files to the D3 and D4 lists, or archive the lab script under CLEAN. Add a check in
    `tests/test_experiment_scripts.py` that imports every non-lab script.

20. **minor**, D.md §2 "Cao et al., PBA+ (cuCIM distance_transform_edt). Rejected: ... indices can differ on ties".
    Problem: the reason does not apply. Phase D needs distances only. cuCIM says that ties change only the returned
    indices, and that `float64_distances=True` matches SciPy [web: docs.nvidia.com/cucim/latest/api/].
    Fix: give the real reasons. The Mac has no CUDA, so the one-source CPU build needs our own EDT as the oracle.
    cuCIM and CuPy are a large new dependency for one kernel. Optionally use cuCIM in a box-only test as a second
    oracle for the SP-K2 EDT.

21. **minor**, D.md §4 D7 step 3 and §6 G9 ("NSD and the soft-Dice mean agree within 1e-12").
    Problem: soft Dice is exact from integer counts. A numpy bool mean is the correctly rounded c / K. `_dice` uses
    integer counts. The mean of 5 floats keeps the same order [code: eval/metrics.py:155-163]. A tolerance would hide
    a real cut error.
    Fix: require exact equality for `soft_dice`. Keep the 1e-12 tolerance for `nsd` only, where float surfel lengths
    are summed in a different order.

22. **minor**, D.md §4 D1 step 2 and D7 step 4 ("tests/goldens fixtures").
    Problem: `tests/goldens/` holds only `polycore.npz`, with polygon keys, defects and even-odd words [measured: ls
    and np.load]. It has no completions, no tables and no EDT inputs. The inputs of the benchmark and of the EDT
    gate are therefore not defined.
    Fix: name the inputs. For example, use the saved R5 test answers (`runs/r5/eval/...`) for 256 completions on
    their items' caches. For the EDT gate, use `ExpertCache.build` on 1,000 LIDC items from `data/prepared/caches`.
    Write each input file's SHA-256 into `runs/d1/bench.json`.

23. **minor**, D.md §4 D5 step 2 (device bit-packing in `score_masks` and `add_tables`).
    Problem: sending bytes instead of words multiplies the H2D bytes and the pinned staging size by 8. The host loop
    it replaces is small (D1.md C12: "small").
    Fix: merge it only if the replay benchmark shows `add_tables` and `score_masks` time is not worse. Item 7's
    device-side table build makes most of this step unnecessary, so do item 7 first.

24. **minor**, SP-K1 scan `cub::DeviceScan::ExclusiveSum(..., int(n), st)` [code: K1/k1poly.cpp:140-142].
    Problem: the item count is cast to `int`. This is safe under 2^25-pixel chunks, but a single 4096 x 4096
    checkerboard has runs near 2^23 per class, and a future change of `chunk_pixels` could overflow it silently.
    Fix: when porting in D2, add `if (n > INT_MAX) throw std::logic_error(...)`, or pass an `i64` count if the
    installed CCCL 13.0.85 header declares `num_items` as a template type (check `cub/device/device_scan.cuh` on the
    box; not verified here).

25. **minor**, D.md §1 item 7 (3D direction) and §3.3 (no row for `checker/volume.py`).
    Problem: the plan uses the 3D direction as a reason for Phase D, but `betti3d` still uses `scipy.ndimage.label` on
    the CPU [code: checker/volume.py:56-57], and the §3.3 table does not list it. `reward/cache.py:127` (cv2.erode
    boundary) and `eval/probes.py:236-268` (ndimage label and EDT) are also missing from the table.
    Fix: add §3.3 rows with a reason for each one ("3D: deferred to the voxel phase; the run CCL extends to z-runs
    with 26/6 connectivity" and "probes: offline report"). Or move them to the device in a D7c extension.
