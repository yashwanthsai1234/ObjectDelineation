# Final review round 8: CUDA and GPU performance (Implementation.md)

Reviewer role: CUDA / GPU performance (kernels, GPU placement, syncs, memory, rental plan). Sources read by heading:
00_front sections 2.4-2.5 and 6 (6.1-6.4); D.md D.0-D.8 (D1, D5 incl. D5c and the round-7 `label_pairs` rewrite,
D7, D7b, D8, D.6 gates, D.8 cost) and D.9.3/D.9.6/D.9.9 (I14, S10, GAP-5, GAP-21); C.md C1; B.md B3.9 / I22b (the
`b3_v5.yaml` content list). Review logs of rounds 1-7 (review-log.md, Phase D and Front) were read, so no fixed item is
repeated: r7 cuda-1 (dump-replay) and cuda-2 (`label_pairs` as a row count) are fixed; r7 minors are not re-filed.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No new catalog skill is needed; D.0 already cites `compute-sanitizer-diagnose-cuda-bugs`. Nothing installed. It pointed me at the sanitizer/driver pairing (item 6). |
| engineering:code-review | Read the cited code: `cuda/gpucheck.cpp` (alloc/release :78-84, `Block` :131-150, `Drain` :177-184, `Layout` :188-197, atomics :201-240, `ObjOut` :1199-1202, label chain :1301-1392, `score_words` :1549-1643, `add_tables`/`drop_tables` :1720-1776), `experiments/cuda_checker/run_box.sh`, `experiments/phase2/vast/policy.py` (`code_files` :374-383, `chain`, `cmd_run` :464-496, `cmd_pull` :499-508), `Makefile`, `.gitignore`, `reward/cache.py`, `configs/checker_gpu.yaml`, K2 `k2score.inc` warp helpers. This found items 1, 3, 4 and 7. |
| engineering:debug | Traced every R-D2 output file from the process that writes it to the Mac: `run_box.sh` cds to `/workspace/src`, `pull.sh` copies only `/workspace/runs/<run>` (item 1). Traced each V3 box gate to the code that would record it (item 2). |
| engineering:system-design | Checked the rental flow (sync -> run -> pull -> teardown) against the gate evidence each rental must bring home. Items 1 and 2 are both "data made on the box that never reaches the Mac". |
| engineering:architecture | The D5c pool design (own `cudaMemPoolCreate` pool, no `cudaDeviceSetMemPool`, `cudaFreeAsync` in lane order) and the two-labeller decision stay sound. Item 7 notes the API seam that D5c needs (`alloc` has no stream). |
| engineering:documentation | Each blocking or major item quotes the doc line (source file and heading) and the code line that proves it, and gives an exact fix. |
| huggingface-skills:huggingface-best | No model choice is in this scope; no change. |
| huggingface-skills:huggingface-papers | The D.2 papers are cited by venue; no paper fact decides an item; no change. |

Web checks (WebSearch, WebFetch, 2026-11-08): NVIDIA's minor-version-compatibility page lists three limits (features
that span toolkit and driver, PTX on older drivers, `-arch=sm_xx` needed) and says nothing about developer tools
[web: docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html]. CUDA 13.x minor compatibility starts
at R580 [web: docs.nvidia.com/cuda/cuda-toolkit-release-notes]. I found no NVIDIA page that states a driver floor for
compute-sanitizer 13.4. So item 6 stays minor.

Verified as correct (no item): the r7 `label_pairs` expected value (each object with area > 0 has bh + 2 crop rows
[code: gpucheck.cpp:1583]); the `i64` `atomic_add` works in both builds [code: :209-215]; the new `ObjOut` field rides
the existing D2H copy before the wait at :1626, so D-G8's 4 + 2 holds. The S10 import closure (`objdel/__init__.py`
imports nothing). The old venv runs the bench by path, not as `-m`, which is right because the tar holds no
`experiments/`. The box lock has `opencv-python-headless`, `scipy 1.18.1` and `numpy 2.5.4` for the oracle tests
[code: experiments/phase2/docker/locks/vllm_native.txt:376, 487, 762]. K2 `warp_kernel` returns whole warps together,
so the `0xFFFFFFFF` masks are valid [code: K2/k2score.inc:46-49, 66-93].

## Items

### 1. BLOCKING: every R-D2 result file is written inside `/workspace/src`, which `pull.sh` never copies, so D-G5, D-G7 and D-G8 evidence is lost at teardown (D; 00_front 6.3 V3 item 1)

- Doc, D.md "D8", `sanitize`: "It writes each tool's exit code and wall time to `runs/d/sanitize.json`."
  D8 `--gpu-check` step 3: "`"$PY" -I -m pytest -q $D_GPU_TESTS -m gpu --junitxml runs/d/gpu.xml`". Step 6: "`"$PY" -I
  -m experiments.phase2.perf.bench_checker --dump-replay runs/d/replay_inputs.npz` ... Then `bench_checker.py --engine
  old|current` (D-G7) and `bench_checker.py --engine old|current --replay runs/d/replay_inputs.npz` (D-G8)". No `--out`
  is given for the box benches. D1 step 2 defines the output as "`bench_checker.py --out runs/d1/bench.json` measures
  each workload: ... median and IQR". D-G5 row: "times recorded in `runs/d/sanitize.json`". D-G7 row: "Pass: new median
  <= 1.05 x old ... Report median and IQR". D-G17 row: "`runs/d/topology.json` records which method ran".
- Code: `run_box.sh` starts with `cd "$(dirname "$0")/../.."` [code: experiments/cuda_checker/run_box.sh:16], so on
  the box every relative `runs/d/...` path is `/workspace/src/runs/d/...`. `cmd_run` runs each job from
  `cd /workspace/src` and exports `OBJDEL_RUN_DIR=/workspace/runs/<run>` [code: experiments/phase2/vast/policy.py:489,
  492]. `cmd_pull` copies only `/workspace/runs/<run>/` to the Mac [code: policy.py:505-507], and teardown destroys
  the box after its final pull. `runs/` is git-ignored [code: .gitignore:23], so nothing under `/workspace/src/runs`
  exists in the code sync either. The repo's own GPU config shows the rule: `checker_gpu.yaml` copies its results into
  `{rd}` by hand before the job ends [code: experiments/phase2/configs/checker_gpu.yaml:24].
- Effect: an engineer who follows D8 exactly gets the pytest and sanitizer exit codes in the pulled job log, but the
  D-G7/D-G8 medians and IQRs, the sanitizer times and the JUnit file stay on the destroyed box. D-G7 and D-G8 have no
  numbers, so the bench minutes of V3 item 1 (about 6 min plus the dump, D.8) are paid for nothing, and D's speed
  claims stay ungated until another rental.
- Fix (D.md D8 `sanitize`, `--gpu-check` steps 3 and 6, D1 step 2, D-G5, D-G7, D-G8, D-G17, B03; 00_front 6.3 V3 item 1
  stays priced the same):
  1. `run_box.sh --gpu-check` sets `OUT="${OBJDEL_RUN_DIR:?run under run.sh}/d"` and `mkdir -p "$OUT"`. The `--dry`
     mode sets `OUT=$(mktemp -d)` instead.
  2. Every output uses `$OUT`: `--junitxml "$OUT/gpu.xml"`, `"$OUT/sanitize.json"`, `--dump-replay
     "$OUT/replay_inputs.npz"`, `--out "$OUT/bench_<engine>.json"` (D-G7) and `--out "$OUT/replay_<engine>.json"`
     (D-G8). `bench_checker.py` makes `--out` required.
  3. After the four bench calls, one `"$PY" -I -m experiments.phase2.perf.bench_checker --compare "$OUT"` prints the
     D-G7/D-G8 table (median, IQR, ratio, pass) to stdout, so it is also in the pulled job log, and exits 1 when a
     pass rule fails.
  4. D-G17's sampler writes `"$OBJDEL_RUN_DIR/d/topology.json"` (see item 2).
  5. New S-test in `tests/test_box_inputs.py` (no `data/`): `run_box.sh` has no output path that starts with `runs/`
     (grep), and `--gpu-check --dry` writes only under its temp `$OUT`.

### 2. MAJOR: the V3 box halves of D-G14, D-G15 and D-G17 have no code that records them, so the gates cannot be read after V3 (D; B; 00_front 6.3)

- Doc, D.md "D.6" D-G14: "V3, first 6 steps of the B3 run: the reward path's `arena_bytes` and
  `torch.cuda.max_memory_reserved()` go into the run JSON ... the B3 run's val eval ... records the `eval_sets` time,
  `arena_bytes` and peak memory". D-G15: "Box, V3 B3 run, steps 2-6: the `gc.stats()` delta of `sync_alloc` and of
  `sync_free` is 0 ... `pinned_grows == 0` over the whole run, read at the end". D-G17 (a): "record the row count of
  `nvidia-smi --query-compute-apps=...` ... sampled every 5 s ... `runs/d/topology.json` records which method ran";
  (c): "Peak `nvidia-smi --query-gpu=memory.used` <= 30 GB". D.8 row "D-G14/D-G15/D-G17": "0 extra". 00_front 6.3 V3
  notes: "D-G14, D-G15 and D-G17 read V3 ... They add no minutes." D7 step 5: "In V3, record the reward path's
  `gc.stats()["arena_bytes"]` ...".
- Proof: no task names the code that writes these numbers. In D.md, `train_rl` and `rl.json` appear only as evidence
  (D.1 item 9, D.3.3), never in a file list; `final_eval.py` appears only in D7 step 6 (pool removal) [grep of D.md,
  2026-11-08]. B.md I22b lists the full content of `b3_v5.yaml` ("same start checkpoint ... LocGuard and BudgetStop on,
  no test eval, `train_rl --out {rd}/b3-v5`") with no engine-stats key and no sampler step. Today's code has no
  `gc.stats()` at all (C.R3 item 27: "`gc.stats()` named by the review does not exist in today's engine"), and
  nothing samples `nvidia-smi` during `final_eval`. "Before scoring starts" (D-G17 a) is a point inside the
  `final_eval` process, after vLLM generation (D7 step 6), so a shell sampler outside the process cannot find it.
- Effect: V3 runs the B3 run and its val eval, but the engineer must invent where and when to read `gc.stats()`, how
  to take step deltas for steps 2-6, and how to start the D-G17 sampler at the scoring point. A guess made after V3
  starts cannot be fixed without another rental, so the gates are untestable as written.
- Fix (D.md D7 step 5 and D7 step 6 file lists, D-G14, D-G15, D-G17, B04, B06; B.md I22b; front 6.3 V3 note):
  1. `src/objdel/train/fast.py` (or the `train_rl` log callback): when `objdel._gpucheck.build() == "cuda"`, log
     `engine/<key>` for every `gc.stats()` key at the end of steps 1-6, and `engine/pinned_grows` at the end of the
     run, into `rl.json`; also `torch.cuda.max_memory_reserved()` for steps 1-6. D-G15 reads the steps 2-6 deltas
     from these keys.
  2. `experiments/phase2/final_eval.py`: around the `eval.device.item_counts` loop (D7 step 6), record the
     `eval_sets` wall time, `gc.stats()["arena_bytes"]` and `torch.cuda.max_memory_reserved()` into the
     `final_eval.json` `seconds`/`engine` blocks. Start one sampler thread just before the loop (`nvidia-smi
     --query-compute-apps` and `--query-gpu=memory.used` every 5 s, with the `/proc` fallback of D-G17 a); stop it
     after the loop; write `$OBJDEL_RUN_DIR/d/topology.json` with the method, the counts and the peak. The thread
     runs only when `OBJDEL_RUN_DIR` is set.
  3. Mac tests: a fake `gc.stats()` and a fake `nvidia-smi` on `PATH` give the expected keys and file. I22b asserts
     nothing new in `b3_v5.yaml` (the recording is in code, on by default on a CUDA build).

### 3. MAJOR: the D-G7 old-engine tar has no stated Mac path, no box extract step and an undefined `<dir>`; on most paths the code sync drops it (D)

- Doc, D.md "D.6" D-G7: "On the Mac (read-only, no checkout): `git archive -o old_engine.tar <pre-D2 sha> cuda src
  CMakeLists.txt pyproject.toml Makefile README.md`. The tar ships in the code sync. ... On the box, inside `run_box.sh
  --gpu-check` step 5 ...: `uv venv /opt/old && uv pip install --python /opt/old/bin/python numpy && make -C <dir> cuda
  PY=/opt/old/bin/python`. `make gpu-check DRY=1` builds the tar on the Mac CPU first". D8 step 5: "`uv venv /opt/old`,
  then `make -C <old dir> cuda PY=/opt/old/bin/python`". D8 `--dry`: "builds the old-engine tar (D-G7) on the CPU in a
  temp folder".
- Code: the code sync ships `git ls-files -co --exclude-standard` minus `^(data|docs)/`, lab and results folders
  [code: experiments/phase2/vast/policy.py:377-380]. `/runs/` and `/build/` are ignored [code: .gitignore:16, 23]. So
  the tar ships only if it sits at a non-ignored path outside `data/`, `docs/`, `runs/`, `build/` and the lab and
  results folders. The doc names none. The `--dry` mode builds its own tar in a temp folder, so it proves neither the
  shipped path nor the extract step. No step runs `tar -x`, and `<dir>` / `<old dir>` is a placeholder. `<pre-D2 sha>`
  has no tag or rule (the next commit after D1 merges, or the parent of the D2 merge?).
- Effect: the natural choices (`runs/d/old_engine.tar` like every other D output, or `build/`) are never synced, so
  step 5 fails in paid R-D2 time and, under `set -euo pipefail` [code: run_box.sh:15], step 6 (D-G7, D-G8) never runs.
- Fix (D.md D-G7, D8 steps 5 and `--dry`, D9 runbook):
  1. Old tree: the owner tags it (owner command, given to them: `git tag pre-d2 <parent of the D2 merge commit>`).
     D-G7 uses `pre-d2`, not "pre-D2 sha".
  2. Mac: `make old-engine-tar` writes `old_engine.tar` at the repo root (`git archive -o old_engine.tar pre-d2 cuda src
     CMakeLists.txt pyproject.toml Makefile README.md`) and `old_engine.tar.sha256`. Both ship because they are
     untracked, not ignored and outside the drop regex. Delete both after the rental.
  3. Box, step 5: `sha256sum -c old_engine.tar.sha256 && mkdir -p /opt/old_src && tar -xf old_engine.tar -C
     /opt/old_src && uv venv /opt/old && uv pip install --python /opt/old/bin/python numpy && make -C /opt/old_src cuda
     PY=/opt/old/bin/python`.
  4. `--dry` runs the same three lines on the repo-root tar (CPU build, temp copy of `/opt/old_src`), and S08 asserts
     that `policy.code_files()` lists `old_engine.tar` when it exists.

### 4. minor: the `relabelled` counter is said to be zeroed by the `changed` zero call, which is sized for `changed` only (D)

- D5 step 1: "each call runs `zero(changed, sizeof(i32) * n, st)` before `CropFill`" and "A per-object `i64
  relabelled[n]` sits in the `score` arena and is zeroed by the same `zero` call as `changed`". `Layout::put` puts
  each array at its own 256-byte aligned offset [code: gpucheck.cpp:188-197], so a call of `sizeof(i32) * n` bytes
  cannot reach `relabelled`. I14's second (stale-flag) check catches it at $0 on the Mac, so no money is at risk.
- Fix: "place `changed` (i32[n]) and `relabelled` (i64[n]) as one `Layout` slot `flags` of `12 n` bytes (relabelled
  first, then changed), and zero it with one `zero(flags, 12 * n, st)` before `CropFill`".

### 5. minor: `--dump-replay` cites the pre-D7b `device()` code that D7b removes (D)

- D1 file list: "per cache the 13 `add_tables` arguments ... exactly as `ExpertCache.device()` passes them today
  [code: reward/cache.py:247-253]". R-D2 runs after D7b, which changes `device()` to call `add_tables_from_raters` and
  moves the host `_band` path to `tests/oracles/reward_reference.py` (D7b). Both sources give the same bytes (D-G10),
  so no result changes.
- Fix: "from the host oracle `tests/oracles/reward_reference.py` (the pre-D7b `device()` arguments), imported inside
  the `--dump-replay` function".

### 6. minor: compute-sanitizer 13.4.92 has only run on driver 615.71.09, while the offer filter accepts any R580 driver (00_front, D)

- 00_front 6.1 step 5: "driver R580 or newer (CUDA 13)". D.4 toolchain: "`nvidia-cuda-sanitizer-api==13.4.92`". The
  only measured sanitizer run (k1, `ERROR SUMMARY: 0 errors` [measured: runs/k1/cuda_checker/memcheck.txt]) was on
  driver 615.71.09 [measured: runs/k1/box.txt]. NVIDIA's compatibility pages give no driver floor for developer tools
  (web checks above), so a 13.4 tool on a 580 driver is unproven. A tool failure stops V3 item 1 (D-G5).
- Fix: make `sanitize` start with a 10 s smoke (`compute-sanitizer --tool memcheck "$PY" -I -c "import objdel._gpucheck
  as g; g.check_masks(...)"` on one 8x8 mask). If it fails with a driver message, record "sanitizer unsupported on
  driver X" in `$OUT/sanitize.json` and fail the step with that reason. Or require driver >= 590 in the V3 offer pick.

### 7. minor: two D5 code seams are left implicit (D)

- D5c: "Device `alloc` calls `cudaMallocFromPoolAsync(&p, bytes, pool, lane_stream)`. `release` calls `cudaFreeAsync`
  on the lane stream." Today `alloc(bool, size_t)` and `release(bool, void*)` take no stream, and `Block<true>` does
  not know its lane [code: gpucheck.cpp:78-84, 131-150]. Fix: "`alloc`/`release` gain a `Stream` argument; `Block`
  stores its lane's stream at `Engine::start`; `drop_tables` passes `lanes[0].stream`."
- D5 step 1: "The second-pair `LabelInit{c, 1}` does `atomic_add(&relabelled[o], 1)`". `LabelInit` is the same functor
  for the first pass [code: gpucheck.cpp:1302-1312, 1612-1614]. Fix: "add `const i32* changed` and `i64* relabelled`
  fields to `LabelInit`, `LabelMerge` and `LabelFlatten`, null in the first pass; when non-null, a row whose object has
  `changed[o] == 0` returns at once".

### 8. minor: D-G8's cold arm says "1 wait per step", but the replay has no steps (D)

- D-G8: "Second replay with cold caches ... and 1 wait per step for table builds (D7b)". Fix: "1 wait per cold repeat
  for the table build (one `add_tables_from_raters` call for all caches of the replay)".

### 9. minor: D-G17 (b) checks pool workers that D7 step 6 removes from `final_eval` (D)

- D-G17 (b): "Each pool worker still alive before CLEAN16 logs `torch.cuda.is_initialized()` and
  `gc.stats()["contexts"]`". D7 step 6: "`final_eval`, `select_ckpt` and `g_probe` stop creating the pool for
  metrics". In the V3 val eval, (b) has no worker to check. Fix: "(b) applies only to a process that still forks
  (none in `final_eval` after D7 step 6); record `n_workers = 0`."

Counts: blocking 1, major 2, minor 6.

Sources (web): [CUDA minor version compatibility](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html),
[CUDA compatibility](https://docs.nvidia.com/deploy/cuda-compatibility/index.html),
[CUDA Toolkit release notes](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/).
