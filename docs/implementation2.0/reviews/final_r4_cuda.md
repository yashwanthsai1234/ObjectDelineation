# Final review, round 4: CUDA and GPU performance (whole Implementation.md)

Reviewer role: CUDA and GPU performance (NVIDIA CUDA guides, Nsight, Blackwell sm_120, CCCL). Target:
`REPO/docs/implementation2.0/Implementation.md` (14,787 lines, built from `docs/implementation2.0/*.md`). Read by section:
front 1, 2.4, 2.5, 5, 6 (all), A6, B1.2, B1.5, B2.2, B3.1, C header to C.8, D (all), E.5, CLEAN L5/L7 tests. Earlier
files read: `reviews/final_r3_cuda.md` and the round-3 and round-4-prep rows of `review-log.md`. Round-3 CUDA items
1-13 are handled (item 9 rejected with a code reason, which I accept), and this file does not repeat them.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | Read the saved live catalog (`(planning probe, not kept: nv_catalog_live.json)`). CUDA-related entries: `compute-sanitizer-diagnose-cuda-bugs` (already named by D), `nemo-mbridge-perf-cuda-graphs`, three `tilegym-*` kernel skills (E.5 already names one). None has an Nsight Systems or venv-toolchain recipe, so items 3, 5 and 6 rest on the NVIDIA docs and repo logs. Nothing installed. |
| engineering:code-review | Checked each GPU claim against code: `gpucheck.cpp` (`Block`, `Lane`, `Drain`, `alloc`, `ScoreRows`, line refs), `CMakeLists.txt`, `Makefile`, `run_box.sh`, `checker_gpu.yaml`, `rl_full.yaml`, `r5_shared.yaml`, `policy.py` (`cmd_sync_up`, `chain`, `cmd_launch`, `cmd_run`, `accrued`), `reward/cache.py`, `serve.fork_pool`, the box lock `vllm_native.txt`, torch's `CUDAFunctions.h`. Found items 2, 4, 5, 7. |
| engineering:debug | Reproduce-first lens on paid steps: which box step fails or passes vacuously with no bug (items 2, 4), and which Mac gate fails on a correct tree (item 1). |
| engineering:system-design | Rental timing from ledger rows instead of the planning constant (item 3); dependency order of D and the CLEAN waves (item 1). |
| engineering:architecture | One rule per fact: CLEAN's "cuda untouched" constraint against D's edits of the same files (item 1). |
| engineering:documentation | Stale or wrong facts: "the box image has no nvcc" (item 5), range names in 2.5 (item 7), rounding (item 8), `ObjOut` size after B3.1 (item 9). |
| huggingface-skills:huggingface-best | No model choice here. Used only to confirm the model under test (Qwen3-VL-2B, vLLM 0.32.0 path) when checking option F's base push (see web checks). |
| huggingface-skills:huggingface-papers | No CUDA claim in this round rests on a new arXiv paper; the D.2 CCL/EDT sources are journal papers. No item. |

Web and repo checks made for this review (2026-11-08):
- vLLM option F base push under level-2 sleep. vLLM issue #39934: with `enable_lora=True`, a reload after level-2
  sleep raised `KeyError` at `qwen2.py` `params_dict[name]`, because LoRA wrapping moves `qkv_proj.weight` to
  `qkv_proj.base_layer.weight` [web: github.com/vllm-project/vllm/issues/39934]. Fix PR #39935 merged 2026-08-06
  (`BaseLayerWithLoRA.load_weights` delegates to the base layer) [web: github.com/vllm-project/vllm/pull/39935].
  vLLM 0.32.0 was uploaded 2026-10-18 [measured: pypi.org JSON], and its `lora/layers/base.py` has the delegating
  `load_weights`, and `qwen2.py` `Qwen2Model.load_weights` uses `AutoWeightsLoader` [measured: raw files at tag v0.32.0,
  `(planning probe, not kept: r4cuda)/`]. So C5's F mode-B base push is expected to work on the pinned vLLM. No item; the G-EQ5 box check stays.
- `cuda-toolkit==13.0.3` has an `nvtx` extra (`nvidia-nvtx==13.0.85`) [measured: pypi.org JSON]. D1's install line is
  valid.
- Compute Sanitizer release notes: memcheck stream-ordered allocation race detection is on by default
  [web: docs.nvidia.com/compute-sanitizer/ReleaseNotes]. D5c's pool design (alloc, use and free on the same lane
  stream) fits it. No item.
- Sanitizer 13.4.92 runs on the box driver: `runs/rl/cuda_checker/memcheck.txt` ends `ERROR SUMMARY: 0 errors`; the
  box driver is 580.159.03, CUDA 13.0 [measured: runs/r5/env_check.json]. No item.
- `nsys status -e` on Linux prints a "CPU Profiling Environment Check", and "if the CPU Profiling Environment ... is
  not OK, you will still be able to run various trace operations" [web: docs.nvidia.com/nsight-systems/InstallationGuide].
  Its exit code in that case is not documented (item 6).
- Recomputed: V1 188.3 min = $2.233 and 169.1 min = $2.017; worst case $2.49; B-S2 177.9 min = $2.116; V2 required
  103 min = $1.267; items 1-5 151 min = $1.811; D share 79 x 0.641 / 60 = $0.844; pinned bound 256 x 1000 x 10 x 8 B
  = 20.48 MB; patches 32 x 1,024 x 1,536 x 4 B = 201.3 MB; images 32 x 512 x 512 x 3 B = 25.2 MB.

Counts: **blocking 2, major 2, minor 6.**

---

## Blocking

### 1. [blocking] CLEAN (L5 CT-55, L7 CT-72, CLEAN.spec C1) vs D (D1, D2-D8 file lists): the CLEAN wave gates require `cuda/` and `CMakeLists.txt` unchanged since a tag that D edits after

- Problem. CLEAN's owner gate script CT-72 (L7, "run at the end of each wave") checks:
  "`git diff --name-only pre-clean-2026-11-08 -- cuda CMakeLists.txt` is empty (CLEAN.spec constraint C1);
  `build/.gpucheck-stamp` not newer than the `cuda/` sources" (CLEAN, CLEAN.9.3 L7). CT-55 (L5, CLEAN16.1, gate
  CLEAN-G2) checks: "`git diff --stat pre-clean-2026-11-08 -- cuda/ tests/test_cuda_checker_gpu.py
  tests/test_gpu_scoring.py` is empty". The tag is made once, in CLEAN0 step 6 ("Owner runs: `git tag
  pre-clean-2026-11-08`"), in Wave 1, before V1. Steps 2-5 are re-run at each wave start; step 6 is not. Wave 2 is
  "after phase D3" and Wave 3 (CLEAN16) is "after phases D and E" (CLEAN.3.1 diagram; CLEAN16 heading).
- D edits exactly these paths before Waves 2 and 3: D1 "change `cuda/gpucheck.cpp`: add `#include
  <nvtx3/nvtx3.hpp>`" and "change `CMakeLists.txt`, inside `if(OBJDEL_CUDA)` ... `find_file(OBJDEL_NVTX3_HPP ...`";
  D2 deletes `check_mask_chunk` and friends from `cuda/gpucheck.cpp`; D8 "deletes the file-level
  `skipif(_gpus() == 0, …)` line [code: tests/test_cuda_checker_gpu.py:27]" (the line exists today: `pytestmark =
  pytest.mark.skipif(_gpus() == 0, ...)`, tests/test_cuda_checker_gpu.py:25); D adds GPU cases to
  `tests/test_gpu_scoring.py`. `make test` also refreshes `build/.gpucheck-stamp` after any `cuda/` edit
  [code: Makefile:14-16].
- So at the end of Wave 2, CT-72 fails on a correct tree, and at Wave 3, CT-55 and CT-72 both fail. The spec text is
  the cause: CLEAN.spec C1 says "No change to ... any file under `cuda/`" [SP/impl/phases/CLEAN.spec.md:167-168],
  which is right for CLEAN's own commits but wrong as a diff against a tag that predates D. Two parts that both claim
  to be the rule; the engineer cannot close the wave gate.
- The same applies to CT-72's "`uv pip freeze` hash before and after": DOC adds the docs-site packages and D/C may add
  tools between Wave 1 and Wave 3.
- Fix (CLEAN, one place each):
  - CT-72 and CT-55: replace the tag diff by a per-commit check of CLEAN's own commits. In CT-73 (run before each
    owner commit of a CLEAN task), add: "`git diff --cached --name-only` contains no path under `cuda/`, no
    `CMakeLists.txt`, no `tests/test_cuda_checker_gpu.py`, no `tests/test_gpu_scoring.py`". Drop the
    `pre-clean-2026-11-08` diff and the stamp rule from CT-55 and CT-72.
  - If a wave-level check is still wanted: the owner tags `pre-clean-w2` and `pre-clean-w3` at each wave start (CLEAN0
    step 6 repeated), and CT-72 diffs against the tag of the current wave; also state "no D or E commit lands while a
    CLEAN wave is open".
  - CT-72's `uv pip freeze` part: compare against the freeze recorded at the current wave start (`before_w2.json`,
    `before_w3.json`), not Wave 1.
  - CLEAN.spec C1 (and its copy in CLEAN.9): "No CLEAN commit changes any file under `cuda/`" instead of "No change to
    ... any file under `cuda/`".

### 2. [blocking] 00_front (6.2 plan files) and B (B1.5 `b1_rows.yaml`, G08): V1 never builds the CUDA checker, so the B1 run and every V1 eval score on the CPU build, and B's GPU smoke passes by skipping

(Found independently of the round-4 QA review, item 2; the venv measurement below is new.)
- Doc lines. 00_front 6.2 "Plan files": "`b1_rows.yaml` (B1.5 writes the commands): (1) `env_check.sh`; (2) B's GPU
  smoke, `pytest -m gpu tests/test_train_fast.py -k "G02 or G08"` (< 5 min, inside item 1); (3) `train_rl` ...;
  (4) one `final_eval` process". `a3_measure.yaml`: "one `final_eval` process on R5". B1.5 lists the same 4 steps.
  No step runs `run_box.sh --build-only` or `make cuda`.
- Code. `sync_up.sh` installs the package with `uv pip install --python /opt/venv/bin/python --no-deps
  --no-build-isolation -e .` [code: experiments/phase2/vast/policy.py:399], and the build defaults to the CPU:
  `option(OBJDEL_CUDA "build objdel._gpucheck with nvcc" OFF)` [code: CMakeLists.txt:15]. A `steps:` config runs its
  literal steps [code: policy.py:432-433]; the automatic `run_box.sh --build-only` plus `tests/test_gpu_scoring.py`
  is added only to non-`steps` configs with a scoring stage, and its docstring states the intent: "scoring runs on the
  GPU only" [code: policy.py:426-428, 450-452]. Measured: R5's `steps:` plan had no build step
  [code: experiments/phase2/configs/r5_shared.yaml], and its box reported `"nvcc": {"present": false}`
  [measured: runs/r5/env_check.json]; the full RL run did build it, as step 3 before `train_rl`
  [code: experiments/phase2/configs/rl_full.yaml:25, 39].
- Effect. The B1 run's rewards and items 3-5 run on the CPU build of `_gpucheck`. That breaks the owner rule "CUDA is
  a must for the checker, the scorer AND the polygonizer. GPU everywhere" (BRIEF section 3). The GPU smoke does not
  catch it: the root `tests/conftest.py` hook (A.9) skips a `gpu` test unless `build() == "cuda"` (D.9 conventions),
  so G02 and G08 skip and pytest exits 0. G08 also contains "G01 parity on 1,000 windows", and G01 checks `sum_x`,
  `sum_y` (B3.1), which is not on the `v1-b1` tag ("It holds no C, D, B2.2-B2.4, B3 or kernel code", front order
  rule 2).
- Fix.
  - 00_front 6.2 and B1.5: `b1_rows.yaml` step 0 = `[bash, -c, "bash experiments/cuda_checker/run_box.sh --build-only
    && /opt/venv/bin/python -I -c 'import objdel._gpucheck as g; assert g.build() == \"cuda\"' && /opt/venv/bin/python
    -m pytest -q -p no:cacheprovider tests/test_gpu_scoring.py"]`, the same pair `policy.chain` adds today.
    `a3_measure.yaml` runs later on the same box and the same editable install, so it needs only the `assert` line.
    Minutes: build about 3 min (D.8 estimate) inside item 1, or item 1 + 3 min.
  - B (G08): split it. G08a (V1) = G02 plus CPU-vs-CUDA parity of today's fields only; G08b (after B3.1) adds
    `sum_x`, `sum_y`. G08 asserts `build() == "cuda"` and fails, not skips, on a GPU box.
  - Side effect to state with the fix [measured: runs/k1/checker_gpu.log:12-22]: `run_box.sh` installs into the
    training venv with no constraints, because the `vllm_native` image has no `/opt/requirements.txt` (its Dockerfile
    copies `vllm_native.txt`). It swaps 4 wheels: `nvidia-cuda-cccl 13.3.4.3.1 -> 13.0.85`, `nvidia-cuda-crt
    13.4.92 -> 13.0.88`, `nvidia-cuda-nvcc 13.4.92 -> 13.0.88`, `nvidia-nvvm 13.4.92 -> 13.0.88`. The full RL run
    trained after the same swap, so it is known to work; write it down, or install the build toolchain into its own
    venv (item 5).
  - With a CUDA build, today's `final_eval` forks up to `POOL_CAP = 16` workers and each opens a CUDA context on its
    first checker call [code: serve.py:132, 158-176]. A2 step 2 already says "On a CUDA build it uses pool size 1";
    apply the same to `final_eval` in V1 (`--procs 1` or `pool_size` capped at 1 when `build() == "cuda"`) until D7
    step 6 lands.

---

## Major

### 3. [major] 00_front (6.2 item 1, notes, 6.4 rule 5): the 30-min overhead is below the measured launch-to-first-step time, so at 27.5 s the B1 run does not reach 300 steps

- Doc lines. 00_front 6.2 table: "1 | Launch, sync, checkpoint and rows-root download ..., `env_check`, B's GPU smoke
  ..., pull, adapter push, teardown; plus the fixed bandwidth b | 30 | 30 | 0.44 | session overhead + b". Notes: "At
  27.5 s the run ends at minute 167.5, before the stop. The stop trips only above about 28.5 s per step [estimate:
  (172.3 - 30) x 60 / 300]." Section 6 inputs: "Start-up per new training process | about 6 min ... not in the step
  minutes".
- Measurement. Launch to the first `run` row in the ledger: r5 15:35:45 -> 16:13:57 = 38.2 min; r4 13:42:58 ->
  14:07:06 = 24.1 min; rl 04:07:51 -> 04:25:30 = 17.7 min [measured: experiments/phase2/vast/ledger.csv]. V1 then
  runs `env_check` (about 1 min: R5's `run` row 16:13:57, vLLM init 16:14:55 [measured: runs/r5/r5_shared.log:39]),
  the GPU smoke (< 5 min, doc), and the `train_rl` start-up (about 6 min, doc). So step 1 starts at session minute
  about 25 to 50 (plus item 2's 3 min if the build is added). `BudgetStop` at $2.05 trips at session minute 172.3
  (doc). At 27.5 s per step and an r5-like launch: (172.3 - 50.2) x 60 / 27.5 = 266 steps. The stop trips for any
  step time above (172.3 - 50.2) x 60 / 300 = 24.4 s, not 28.5 s.
- Effect. The adapter is saved (money and the eval tail are safe), but Gate B1 judges a run of about 266 steps, not
  the planned 300; the B1.2 schedule is built for 300 (`--steps <planned_steps>`), so the last 1-2 blocks of rows never
  train. The engineer and the owner approve "300 steps" on a false premise.
- Fix (00_front 6.2, 6.4 rule 5; B1.5 follows):
  - Split item 1 into "before step 1" (launch to first `run` 17.7-38.2 min measured, + `env_check` 1 + smoke 5 +
    CUDA build 3 + start-up 6 = up to 53 min) and "after training" (pull, push, teardown, inside the `PULL_S` window).
  - Then pick one, as owner decision 2:
    (a) keep $2.40 / cap 3.35 h / `BudgetStop` $2.05 and plan `--steps 260` (a multiple of 20):
    53 + 260 x 27.5 / 60 = 172.2 min <= 172.3 [estimate]. This fits the $3 rule as written today.
    (b) keep 300 steps and raise the stop so it trips after minute 190.5 (53 + 137.5): `BudgetStop` >=
    0.10 + 190.5 x 0.680 / 60 = $2.26 [estimate], so $2.30 (minute 194.1). Then recompute the cap, the V1 amount, the
    worst case against HARD_STOP and the $3 total; this likely needs a reserve cut, so it is an owner decision.
  - Replace the "stop trips only above about 28.5 s" sentence with the measured-overhead version, and recompute the
    "money left after V1" lines of 6.3 from the chosen option.

### 4. [major] D (D1 `d_gpu_test.yaml`, D8 step 3, D.9.7 GPU02, D-G9, D-G10): R-D2's box payload is not defined, and `--gpu-check` fails on any skipped test

- Doc lines. D1 file list: "create `experiments/phase2/configs/d_gpu_test.yaml` for R-D2, based on `checker_gpu.yaml`
  (code-only sync, no vLLM model; the LIDC `data/checker_eval` set is copied as K1-draft GPU_RUN §0 does)". D8 step
  3: "`"$PY" -I -m pytest -q $D_GPU_TESTS -m gpu --junitxml runs/d/gpu.xml`; the step fails if `skipped` > 0". D.9.7
  GPU02: "GPU arrays equal CPU-build arrays (npz golden saved on the Mac) for `check_masks`, `topology`,
  `certify_polygons`, `score_polygons`, `eval_sets`, `add_tables_from_raters`". D-G10: "device EDT == scipy, bit for
  bit, on 1,000 LIDC caches ... | Mac + R-D2". D-G9: "`metrics()` old vs new on the saved val answers ... | Mac +
  R-D2". D.9 conventions: "`data` = needs `data/` (skip when missing)".
- Code. A `--code-only` sync ships no data; `checker_gpu.yaml` notes that its eval "stops on the missing checker_eval
  set; that set ships in the LIDC bundle" [code: experiments/phase2/configs/checker_gpu.yaml, header comment].
- Effect. The npz goldens (path not named), the 1,000 LIDC caches and the saved val answers are not in R-D2's payload.
  Any `gpu` test that reads them skips (or errors) on the box, step 3 fails by rule, and V3 item 1 ($0.28-0.50) is
  spent for nothing. `make gpu-check DRY=1` runs only `pytest --collect-only`, which does not evaluate run-time skips,
  so the gap is not seen at $0. Also, D.9 says `data` tests "run only through `make test-data`", while step 3's
  `-m gpu` selects a test that carries both markers.
- Fix (D):
  - D1 / D8: name the payload of `d_gpu_test.yaml` in one list: `data/checker_eval/` (exists), the GPU02 goldens at
    one path (for example `data/d_gpu/goldens/*.npz` with a SHA-256 manifest, written on the Mac by
    `experiments/phase2/perf/make_d_goldens.py`, LIDC only), and, only if D-G9/D-G10 keep a box part, the 1,000 LIDC
    caches and the saved val answers they read (LIDC only). RIGA never ships.
  - Or drop the box part of D-G9 and D-G10 (GPU02 already proves device == CPU build for `eval_sets` and
    `add_tables_from_raters`), so R-D2 needs only `checker_eval` and the goldens.
  - `--gpu-check --dry`: for each `gpu` test of `D_GPU_TESTS`, assert that every data path it declares (a module
    constant `BOX_INPUTS`) exists in the sync list; fail at $0 otherwise.
  - D.9 conventions: a test is `gpu` or `data`, never both; the box half of a data-backed check reads only the payload
    above.

---

## Minor

### 5. [minor] D (D8 "Why one shell entry"), 00_front (6.4 toolchain line): "the box image has no nvcc" is wrong, and the toolchain install changes the training venv before the B3 run

- Doc: D8 "Why one shell entry: the box image has no nvcc." 00_front 6.4: "The box toolchain for V2 and V3 is
  `cuda-toolkit[nvcc,cccl,nvtx]==13.0.3` and `nvidia-cuda-sanitizer-api==13.4.92`. This is what
  `experiments/cuda_checker/run_box.sh` installs".
- Code / measurement: the image lock already holds `nvidia-cuda-nvcc==13.4.92` (via humming-kernels),
  `nvidia-cuda-cccl==13.3.4.3.1`, `nvidia-nvtx==13.0.85` [code: experiments/phase2/docker/locks/vllm_native.txt:409-417,
  471-480]; nvcc is only off `PATH` (`"nvcc": {"present": false}`, runs/r5/env_check.json). The install swaps 4 wheels
  in `/opt/venv` [measured: runs/k1/checker_gpu.log:12-22]. In V3, item 1 (R-D2) runs this install, then item 4 (the B3
  run) uses the swapped venv, while V1 (unless item 2's fix adds the build) did not.
- Fix: D8 sentence: "nvcc is not on `PATH`; the image's wheel nvcc is 13.4.92, and the tested build uses 13.0.88".
  Prefer `uv venv /opt/cudabuild` for the toolkit and sanitizer, `CUDACXX=/opt/cudabuild/.../cu13/bin/nvcc`, and build
  the extension for `/opt/venv` from there, so `/opt/venv` stays byte-equal to the lock; or record the 4-wheel swap in
  6.4 and in `env_check.json` for every run that builds.

### 6. [minor] A (A6 step 8 `nsys_capture.sh`): the capture exits 3 on `nsys status -e`, and the profile command keeps CPU sampling on in a container

- Doc: `"$NSYS" --version > "$RD/nsys/version.txt" && "$NSYS" status -e >> "$RD/nsys/version.txt" || exit 3`, then
  `nsys profile -t cuda,nvtx,osrt --capture-range=cudaProfilerApi ...` (no `--sample` or `--cpuctxsw` flag). C0 step 6a
  already expects "a container often blocks `perf_event_open`".
- Web: NVIDIA says trace still works when the CPU profiling check is not OK [web: Nsight Systems Installation Guide];
  the exit code of `status -e` in that case is not documented. If it is non-zero, the whole nsys path falls back to
  torch.profiler with no reason.
- Fix: `"$NSYS" status -e >> "$RD/nsys/version.txt" || true` (record, do not gate); gate only on `--version`; add
  `--sample=none --cpuctxsw=none` to the profile line (CPU samples are not used by any gate), and S04 checks the text.

### 7. [minor] 00_front (2.5 diagram): NVTX range names differ from A6's fixed list

- Doc: 2.5 "NVTX ranges: generate, reward, old_logps, update, sync (Phase A6)". A6 step 6: "Range names, fixed ...:
  `step, loader, micro, batch, generate, sync, wake, sleep, rewards, checker, logps, forward, optimizer, log, ...`".
- Fix: 2.5 lists A6's names (`generate`, `rewards`, `checker`, `logps`, `forward`, `micro`, `sync`) and drops
  `old_logps` and `update`.

### 8. [minor] 00_front (6.2 total, 6.3 B-S2): two dollar values round down

- 169.1 min x $0.01133 + $0.10 = $2.016 -> $2.02 (doc: $2.01); 177.9 min -> $2.116 -> $2.12 (doc: $2.11). Money left
  at 23.65 s becomes $0.73. Fix the two cells and the sentences that reuse them (6.3 first line, totals for decision 1).

### 9. [minor] D (D5c out bound) and B (B3.1): `ObjOut` grows to 64 B after B3.1, but D5c states 48 B

- D5c: "`ObjOut` is 2 x i64 + 4 x double = 48 B [code: gpucheck.cpp:1199-1202]". B3.1: "`ObjTerms` sums the two
  partials over rows into `ObjOut` (new `i64 sum_x, sum_y`)". The bound still holds (signatures dominate at 256 B per
  object), but the stated term is stale once B3.1 lands.
- Fix: D5c writes "48 B today, 64 B after B3.1 (`sum_x`, `sum_y`)", and its bound test reads `sizeof(ObjOut)` from the
  binding instead of a literal.

### 10. [minor] C (C1 step 2, G-PX2): the per-row index of `per_image` False must be copied from pinned memory, or G-PX2 counts a sync

- Doc: C1 step 2 "assemble the 256-row output with one `torch.index_select` on the per-image rows (row -> image index
  built on the host)"; G-PX2 "`count_syncs()` ... around one pass-2 `ImageMemo.__call__` of 256 rows after `prefill`
  is 0".
- Code: a blocking host-to-device copy goes through `memcpy_and_sync`, which calls `warn_or_error_on_sync()` when the
  sync debug mode is on [code: .venv torch/include/c10/cuda/CUDAFunctions.h:79-88]. So `torch.tensor(idx,
  device="cuda")` or `.to(device)` from pageable memory reads as 1 sync.
- Fix: C1 step 2 adds "the index goes through a pinned `torch.empty(256, dtype=torch.int64, pin_memory=True)` and
  `.to(device, non_blocking=True)`, as C4 step 4 does", and U-PX tests assert the index tensor is built that way.
