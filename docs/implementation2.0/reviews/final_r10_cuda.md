# Final review round 10: CUDA and GPU performance (Implementation.md)

Reviewer role: CUDA / GPU performance (kernels, GPU placement, syncs, memory, Blackwell sm_120, rental plan).
Read: 00_front (sections 2.4, 2.5, 5, 6.1-6.4), D (all tasks, gates, D.8), C (C.3, C0-C6, C.6), A (A2, A3, A6),
B (B2.1, B3.1), E (E3a, E6, E8, R-4, I29, E-G9), and the round 1-9 logs in review-log.md (items fixed there are
not repeated). Code opened: cuda/gpucheck.cpp, CMakeLists.txt, Makefile, experiments/cuda_checker/run_box.sh,
experiments/phase2/vast/policy.py, experiments/phase2/final_eval.py, select_ckpt.py, src/objdel/train/serve.py,
src/objdel/eval/metrics.py, experiments/phase2/docker/{Dockerfile.vllm_native,Dockerfile.phase2,locks/*},
docs/implementation2.0/kernels/K1/{k1poly.cpp,GPU_RUN.md,build_cpu.sh}, runs/k1/checker_gpu.log,
runs/r3/cuda_checker/memcheck.txt, .venv TRL vllm_generation.py.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog install (review only, no approval to install). It pointed me at the compute-sanitizer tool docs as the source for sanitizer claims; that check gave item 8 (memcheck `--track-stream-ordered-races`) and corrected my own first guess about synccheck and `__ballot_sync` (synccheck does not cover it, so item 7 is about CUB, not the ballot). |
| engineering:code-review | Line-by-line check of the D8 shell steps against run_box.sh:15-37 and the image Dockerfile; found that `CONSTRAINTS` is empty on the production image (item 2). |
| engineering:debug | Traced the post-E paths: `cmd_run` -> `chain()` (I29) for the V3 plan lines (item 1) and the D8 dry-run import of `experiments.phase2.vast.policy` after E8 (item 3). |
| engineering:system-design | Walked the V3 rental flow end to end (6.1 checklist, stamp, launch lines, R-D2 step 0) and found the missing `old_engine.tar` step (item 4). |
| engineering:architecture | Checked D5c's stream-ordered pool against the two-lane engine (which lane allocates and frees tables) for item 8. |
| engineering:documentation | Each item quotes the source file and heading, then the code line or measurement, then one exact fix. |
| huggingface-skills:huggingface-best | No model choice in scope this round; nothing changed. |
| huggingface-skills:huggingface-papers | No paper claim in my area needed a re-check this round (D.2 cites only kernel notes and K1/K2 designs); nothing changed. |

Web checks: compute-sanitizer documentation (docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html, read
2026-11-09): memcheck has `--track-stream-ordered-races all` for `cudaMallocAsync`/`cudaFreeAsync` use on another
stream; leaks are reported only when the context is destroyed (D8's reason stands); synccheck covers
`__syncthreads()`/`__syncwarp()`, not `__ballot_sync`.

Counts: blocking 1, major 3, minor 5.

---

## 1. BLOCKING. After E6, the V2/V3 launch lines name plan files that `chain()` refuses after the paid launch

Sources: 00_front, E (also A, D).

Problem.
- 00_front "6.3 V2, V3 and V4: only after a budget raise", the post-E6 rule: "After E6 (a raise that merges
  `after-raise` after E6, E R-4 "After E6"): `c_launch.yaml` is folded into `specs/c_speed.yaml` ... The V3 lines
  change the same way, to `specs/b3_v5.yaml` (also in `push_hf.sh`). Nothing else in this section changes."
  The V3 launch lines keep `run.sh <run> experiments/phase2/configs/d_gpu_test.yaml` (item 1) and
  `run.sh <run> experiments/phase2/configs/a3_test.yaml` (item 3); the V2 lines keep
  `run.sh <run> experiments/phase2/configs/a3_measure.yaml` (item 1b) and `.../r_d1.yaml` (item 4).
- E "R-4": "After E6 ... `chain()` reads only `specs/*.yaml` (I29)". E "I29": "`chain()` reads only
  `specs/*.yaml`: `chain()` on a path outside `specs/` raises before any `remote` call ... (`cmd_run` calls
  `chain()` only after launch, sync and fetch [code: experiments/phase2/vast/policy.py:469-470])". E-G9: "the box
  chain reads only `specs/*.yaml` (the 4 ported plans, the live A-D plans and V4's `specs/e_v4.yaml`; I29)".
- E "E6" step 2 lists the A-D plans that the owner marks live (ported) or finished (archived): "`a3_measure`,
  `b2_fp16`, `b3_v5`, `c_speed`, `d_gpu_test`". `a3_test.yaml` (A3, A.md "A3" file list: "create
  `experiments/phase2/configs/a3_test.yaml` (plan file): the test half, V3 item 3") is in no list, so no rule ports
  it. E R-4 (iv) ports only the branch plans `c_speed`, `b2_fp16`, `b3_v5`.
- Code today: `cmd_run` builds the chain locally, `stages, name = chain(ROOT / a.module, a.run_id)`
  [code: experiments/phase2/vast/policy.py:471], and it is called only after `launch.sh`, `sync_up.sh` and
  `fetch_hf.sh` in the 6.3 lines.

Failure: under the documented path "a raise after E6" (front "End state under (b')": "If E6 has also landed, the
merge is more than a path map"), an engineer who follows front 6.3 exactly launches V3, syncs, fetches, then
`run.sh <run> experiments/phase2/configs/d_gpu_test.yaml` raises (outside `specs/`), and so does item 3
(`a3_test.yaml`, never ported). The box bills while the engineer ports two plans by hand. Front 6.3 ("Nothing else
in this section changes") and E R-4 / I29 both claim to be the rule and disagree.

Fix.
- 00_front 6.3 "After E6" paragraph: replace "Nothing else in this section changes" by the full list: V2 item 1b
  `run.sh <run> specs/a3_measure.yaml`; V2 item 4 dropped (R-D1 needs a pre-D2 tree, so it never runs after E);
  V3 item 1 `run.sh <run> specs/d_gpu_test.yaml`; V3 item 2b `specs/b2_fp16.yaml`; V3 item 3
  `run.sh <run> specs/a3_test.yaml`.
- E "E6" step 2: add `a3_test` to the list of A-D plans that the owner marks live or finished on the E0 commit
  (live while V3 has not run).
- E I29 / E-G9: add one Mac test ($0) that scans the front 6.3 V2/V3/V4 code blocks for every `run.sh <run> PATH`
  and runs `chain(PATH, "x")` on each; on a post-E6 tree every path must be under `specs/` and load.
- D8 (`d_gpu_test.yaml` bullet and the preflight sentence) and D.8: add "after E6: `specs/d_gpu_test.yaml`".

---

## 2. MAJOR. D-G7's old venv does not get the main venv's numpy: `CONSTRAINTS` is empty on the box image

Sources: D.

Problem.
- D "D8", `--gpu-check` step 5: "`uv pip install --python /opt/old/bin/python "${CONSTRAINTS[@]}" numpy` (the
  image's lock [code: run_box.sh:22-23], so both venvs load the same numpy)". The D-G7 row repeats the same
  command.
- Code: `[ -f /opt/requirements.txt ] && CONSTRAINTS=(-c /opt/requirements.txt)` [code:
  experiments/cuda_checker/run_box.sh:23]. Only the old image writes that file: `COPY train.txt
  /opt/requirements.txt` [code: experiments/phase2/docker/Dockerfile.phase2:37]. The production image is
  `vllm_native` (`objdel-gpu:vllm_native-...` [code: docs/runbooks/gpu-run.md:17]; torch 2.13.0, the box venv
  that D1 names), and it copies `COPY vllm_native.txt vllm_native.override /opt/` [code:
  experiments/phase2/docker/Dockerfile.vllm_native:20]: no `/opt/requirements.txt`. The k1 box log confirms an
  unconstrained install on that image: the toolkit step moved `nvidia-cuda-nvcc` 13.4.92 -> 13.0.88 [measured:
  runs/k1/checker_gpu.log:14-22], which a `-c` lock pinning 13.4.92 [code: locks/vllm_native.txt:417] would refuse.
- So `CONSTRAINTS=()` and the old venv gets the newest PyPI numpy, while `/opt/venv` has numpy 2.5.4
  [code: locks/vllm_native.txt:376; Dockerfile.vllm_native:5 override].

Failure: D-G7 passes or fails on `new median <= 1.05 x old` per workload, timed in two processes that load
different numpy builds. The binding builds its output arrays through numpy at run time (nanobind numpy
ndarrays), so per-call overhead differs between the arms. On the small workloads (lesion512 batch 1) a few
microseconds is the whole 5 % margin. The doc states the opposite invariant, so nobody checks it.

Fix (D8 step 5, the D-G7 row and the `--dry` copy of step 5):
`NP=$("$PY" -c 'import numpy; print(numpy.__version__)')` then
`uv pip install --python /opt/old/bin/python "numpy==$NP"`. Delete "(the image's lock ..., so both venvs load the
same numpy)". `bench_checker.py` writes `numpy.__version__` into `bench_<engine>.json` and `replay_<engine>.json`,
and `--compare` exits 1 when the two arms differ.

---

## 3. MAJOR. After E8, `make gpu-check DRY=1` imports a module that E8 deletes

Sources: D, E.

Problem.
- D "D8", mode `--gpu-check --dry`: "`"$PY" -I -c` asserts that `experiments.phase2.vast.policy.code_files()`
  lists `old_engine.tar` and `old_engine.tar.sha256` (the shipped path)". D8 Steps make this dry run the Mac gate
  before every V3 launch ("On the Mac before the V3 launch ($0 ...): `make old-engine-tar`, then
  `make gpu-check DRY=1`").
- E "E8" step 1: `src/objdel/infra/vast/{offers,ledger,remote,chain,__main__}.py` is the "split of
  `experiments/phase2/vast/policy.py`"; the name table puts `code_files` in `infra/vast/remote.py`. E8 rewrites
  only the vast shell shims: "Shell shims stay in `experiments/phase2/vast/*.sh` ... after E8 they run
  `python -I -m objdel.infra.vast <cmd>`". `experiments/cuda_checker/run_box.sh` is not a shim and is not named.
  E's stale-string test U8 scans only literals matching `objdel(\.\w+)+` in Python files, so a
  `experiments.phase2.vast.policy` string inside a `.sh` file is not caught.
- 00_front "End state under (b')": a raise may come after E ("If E0 has started, the owner rebases the branch onto
  post-E paths"); D8's `run_box.sh` lives on `main`, so it is moved by E, not by the branch rebase.

Failure: on a post-E8 tree, `make gpu-check DRY=1` exits 1 with `ModuleNotFoundError:
experiments.phase2.vast` on the Mac. The engineer cannot finish the V3 pre-launch gate without guessing the new
module name. Cost is $0, but V3 item 1 is blocked until someone edits D's script.

Fix.
- E "E8" step 1: add `experiments/cuda_checker/run_box.sh` to the callers that E8 rewrites: the dry-run assert
  becomes `objdel.infra.vast.remote.code_files()`.
- E U8 (or a new U-test): also scan `*.sh`, `Makefile` and `experiments/**/*.yaml` for
  `experiments\.phase2\.vast\.(policy|bundle|preflight)` and fail on any hit after E8.
- D8: write the assert through a name that E keeps stable, for example
  `"$PY" -I -m experiments.phase2.vast.policy code-files` before E8 and `-m objdel.infra.vast code-files` after,
  or state "E8 rewrites this line (E8 step 1)".

---

## 4. MAJOR. The V3 checklist never runs `make old-engine-tar`; D8's preflight hook can never fire

Sources: D, 00_front.

Problem.
- D "D8", Makefile bullet: "`experiments/phase2/preflight.py` runs `make gpu-check DRY=1` when the config is
  `d_gpu_test.yaml`." D8 Steps: "`make old-engine-tar`, then `make gpu-check DRY=1` ... These are Mac steps before
  the launch, not rental items; no V3 line changes."
- 00_front "6.1 Before any rental (Mac, $0)", step 4: "before V3: `CONFIG=experiments/phase2/configs/b3_launch.yaml`
  ... The plan files `b1_rows.yaml`, `a3_measure.yaml`, `c_speed.yaml`, `r_d1.yaml`, `b3_v5.yaml` and the others
  never follow `preflight CONFIG=`." The V3 launch config holds `rows_plan: experiments/phase2/configs/b3_v5.yaml`
  (front 6.3 "V3 launch config"). So the stamped config is never `d_gpu_test.yaml`, and D8's hook is dead code.
  6.1 has no row for `make old-engine-tar` or `make gpu-check DRY=1` (grep of 00_front and 99_back: 0 hits).
- D "D8" `--gpu-check` step 0 on the box: "`[ -f old_engine.tar ] && [ -f old_engine.tar.sha256 ]` in the repo
  root, else exit 1". The tar reaches the box only through `sync_up.sh` [code: policy.py:396-398 rsync of
  `code_files()`].

Failure: an engineer who follows the front checklist (6.1 steps 0-6, then the V3 lines) never builds the tar. On
the box, `run.sh <run> .../d_gpu_test.yaml` stops at step 0, so D-G3, D-G5, D-G6, D-G7 and D-G8 are not read in
this rental. Recovery needs a mid-rental `make old-engine-tar`, `sync_up.sh <run> --code-only` and a second run,
billed at r_s.

Fix.
- 00_front 6.1 step 4, V3 part: add "when V3 item 1 is approved: `make old-engine-tar && make gpu-check DRY=1`
  (D8), both exit 0, before `sync_up.sh`". Pass column: "exit 0".
- D8: replace the dead hook by "preflight runs `make gpu-check DRY=1` when the stamped config is the V3 launch
  config (`b3_launch.yaml`, after E6 `specs/b3_v5.yaml`)", or delete the sentence and point to front 6.1 step 4.

---

## 5. MINOR. "Pool size 1" still forks one scoring worker, and that worker opens its own CUDA context

Sources: 00_front (6.2 "Why step 0"), A (A3 file list, A2 step 2).

Doc: 00_front 6.2: "On a CUDA build `final_eval` uses pool size 1, as A2 step 2 does, so no forked worker opens
its own CUDA context (A3 adds the cap)." A3: "caps the pool at 1 worker ..., so no forked worker opens its own
CUDA context."
Code: `pool = serve.fork_pool(serve.pool_size())` runs before the model load [code: final_eval.py:69, 72], and
`score_items` maps `_score` over the pool [code: select_ckpt.py:230], so with n = 1 one forked child runs
`metrics()` and the CUDA checker, and opens one context (legal, because the parent forked before CUDA init
[code: serve.py:164-169]). The cap limits extra contexts to one (not up to `POOL_CAP` = 16); it does not remove it.
Fix: write "at most one extra CUDA context (the one worker), about 0.5 GB beside vLLM and the trainer
[estimate]; D7 step 6 removes it". No code change.

## 6. MINOR. D5c growth test: "one call larger than the warm-up call" may not grow a power-of-two arena

Sources: D (D5 step 5, D-G15 Mac part).

Doc: "Then one call larger than the warm-up call, but inside the pinned bound, raises `pool_allocs` by at least 1
(one per grown arena)". Code: `Block::get` keeps the buffer when `bytes <= cap`, and grows to the next power of
two up to 64 MiB [code: gpucheck.cpp:135-148]. A call a little larger than the warm-up fits in the rounded
capacity, makes 0 `pool_allocs`, and the test fails with no bug.
Fix: "one call whose `score` arena need is more than twice the warm-up need (read `arena_bytes` after warm-up and
size the call above it)".

## 7. MINOR. "K1 uses no shared memory and no barriers" is not true of its CUB scan

Sources: D (D-G1 row, D8 sanitizer bullet).

Doc: D-G1: "racecheck/synccheck do not apply [K1-draft GPU_RUN.md:92]"; D8: "racecheck has coverage on the
engine's `raster_kernel` shared buffer ...; K1-draft has no shared memory." Code: K1 calls
`cub::DeviceScan::ExclusiveSum` [code: K1/k1poly.cpp:46, 140-142]; CUB's scan kernels use shared memory and
`__syncthreads`. After D2 the engine holds that scan (D2 file list), so D-G5's racecheck and synccheck also cover
it. The D-G1 conclusion (memcheck and initcheck are enough for the draft) still stands, because D-G5 runs all four
tools on the port.
Fix: D8: "racecheck covers `raster_kernel` and, after D2, the CUB scan kernels". D-G1: "racecheck/synccheck are
not run on the draft; D-G5 runs them on the port".

## 8. MINOR. D5c moves to stream-ordered pool memory but D8's memcheck does not check stream-ordered races

Sources: D (D5 step 5, D8 sanitizer bullet).

Doc: D5c: "Device `alloc` calls `cudaMallocFromPoolAsync(&p, bytes, pool, lane_stream)`. `release` calls
`cudaFreeAsync` on the lane stream ... `drop_tables` frees with `cudaFreeAsync` on the lane stream". The engine
has two lanes [code: gpucheck.cpp:152-166], and `check_polygons` uses both [code: gpucheck.cpp:1671-1692].
Web: memcheck "can detect stream-ordered allocations races using the `--track-stream-ordered-races all` option"
(an allocation freed with `cudaFreeAsync` on one stream must not be used on another without an event)
[web: compute-sanitizer docs, read 2026-11-09]. D8's memcheck line has no such flag, and D5c does not name the lane
for table allocation and free.
Fix: D5c: "tables are allocated, used (`score_words`, `eval_sets`) and freed on lane 0 only; arena `Block<true>`
growth frees on its own lane". D8: memcheck runs with `--track-stream-ordered-races all`.

## 9. MINOR. "The box image has no nvcc" is wrong for the production image (conclusion unchanged)

Sources: D (D8 "Why one shell entry"; D1 `run_box.sh` comment fix).

Doc: D8: "Why one shell entry: the box image has no nvcc." Code: the vllm_native lock installs
`nvidia-cuda-nvcc==13.4.92` and `nvidia-cuda-cccl==13.3.4.3.1` "via humming-kernels"
[code: locks/vllm_native.txt:411-418], and the k1 log shows `run_box.sh` replacing them with 13.0.88 / 13.0.85
[measured: runs/k1/checker_gpu.log:14-22]. The image has nvcc 13.4 inside the venv, not on `PATH`, and no
`CUDACXX`. `run_box.sh:5` says "(the image has none)".
Fix: D8 and the D1 comment fix: "the image's nvcc 13.4.92 sits in the venv, off `PATH`; only `run_box.sh` installs
13.0.88 and exports `CUDACXX`, and only in its own shell". Note the downgrade of a `humming-kernels` dependency in
front 6.2 "Why step 0" beside the "full RL run trained after the same swap" evidence.
