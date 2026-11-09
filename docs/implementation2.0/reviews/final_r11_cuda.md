# Final review round 11: CUDA and GPU performance (Implementation.md)

Reviewer role: CUDA / GPU performance (kernels, GPU placement, syncs, memory, Blackwell sm_120, rental plan).
Read: 00_front (2.4, 2.5, 6, 6.1-6.4), D (D.1-D.9.7, all tasks and gates), C (C.3, C0 step 0, C1), A (A6 step 8
`nsys_capture.sh`), and the round 1-10 logs in `docs/implementation2.0/review-log.md` plus
`reviews/final_r10_cuda.md` and its verdict (items handled there are not repeated).
Code opened: `cuda/gpucheck.cpp` (lanes, `alloc`/`release`, `Block`, `Drain`, `score_words`, label chain, `unite`,
`add_tables`, `drop_tables`, `score_polygons`, `score_masks`, `check_polygons`, `check_masks`, `kChunkSlots`),
`CMakeLists.txt`, `pyproject.toml`, `Makefile`, `experiments/cuda_checker/run_box.sh`,
`docs/implementation2.0/kernels/K1/k1poly.cpp` (scan, `pack_kernel`), `src/objdel/reward/cache.py` (`edt`, `_band`,
`device()`), `src/objdel/eval/metrics.py` (soft Dice, `_surface`, `_nsd`), `src/objdel/train/serve.py`,
`experiments/phase2/train_rl.py`, `experiments/phase2/configs/rl_full.yaml`,
`experiments/phase2/docker/{Dockerfile.vllm_native,build_images.sh,locks/vllm_native.txt}`,
`experiments/phase2/vast/policy.py:333`, `runs/{k1,r3}/cuda_checker/*`, `runs/*/env_check.json`.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Review only, no catalog install (no approval to install). It pointed me to NVIDIA's own docs for the tool-versus-driver question: I checked the CUDA compatibility guide and the compute-sanitizer pages for the driver needed by sanitizer 13.4.92 (item 1). Neither page gives a number, so item 1 stays minor. |
| engineering:code-review | Line-by-line check of D5's label skip against `score_words` and the label kernels (gpucheck.cpp:1597-1616, 1301-1345), and of D5 step 4's `atomicMin` text against `atomic_min_old` (gpucheck.cpp:233) (items 5, 6). |
| engineering:debug | Traced what happens when the sanitizer cannot start on the box: `--gpu-check` runs step 4 (sanitize) before steps 5-6 (old engine, benches), so a start failure would also lose D-G7/D-G8 (item 1). The K1 log shows this failure mode once already (`runs/k1/cuda_checker/memcheck.txt`). |
| engineering:system-design | Walked the V2 path from C0 step 0 (image layer) to the launch (`.vast-allowlist` image digest, `build_images.sh` tag), and found no rebuild/push step (item 2). |
| engineering:architecture | Checked D5c's stream-ordered pool against the two lanes and the `weakref.finalize` drop path (`cache.py:256`). No defect: tables are used and freed on lane 0 only, in stream order, under `e.mu`. |
| engineering:documentation | Each item quotes the source file and heading, then the code line or measurement, then one exact fix. |
| huggingface-skills:huggingface-best | No model choice in my area this round; nothing changed. |
| huggingface-skills:huggingface-papers | D.2 cites only kernel/CCL papers and NVIDIA docs; no arXiv claim in my area needed a re-check; nothing changed. |

Web checks (2026-11-09): CUDA compatibility guide (docs.nvidia.com/deploy/cuda-compatibility/latest): minor-version
compatibility "with some feature limitations", no word on developer tools and no driver number for the tools. A
search on compute-sanitizer driver requirements found only "Compute Sanitizer API calls will fail with
SANITIZER_ERROR_NOT_INITIALIZED if the CUDA driver version is not compatible with the Compute Sanitizer version" (sanitizer API
guide), with no version table. So item 1 is a risk with no proof of failure.

What I re-verified and found correct (no item): the D-G8 baseline (4 `wait()` calls at gpucheck.cpp:1802, 1574,
1626, 1639 plus 2 `Drain` syncs at 177-184); the `ObjOut` size (2 x i64 + 4 x double = 48 B, gpucheck.cpp:1199-1202);
the D5c pinned bound (256 x 1,000 x 10 x 8 B = 20.48 MB per name); the label skip semantics (fg and bg labels live in
separate arrays `o_fg`, `o_bg`; min-root union-find gives the same labels on a relabel, so keeping the first pass is
exact); the soft-Dice cut (numpy `mean` on bool sums in float64, one division) and the NSD cut (scipy returns
`sqrt(float64(d2))`, metrics.py:191-196); the K1 `pack_kernel` warp-uniform exit before `__ballot_sync`
(k1poly.cpp:378-391); `nvidia-nvtx==13.0.85` already in the box lock (locks/vllm_native.txt:479), so D1's `nvtx`
extra changes no wheel; every V1/V2/V3/V4 minute and dollar row in front 6.2-6.3 recomputed at r_s = $0.01133/min
(all match, including the cap 3.35 h, `BudgetStop` $2.05/$2.73 and V4 $0.71).

Counts: blocking 0, major 0, minor 7.

---

## 1. MINOR. Sanitizer 13.4.92 has run only on driver 615.71 hosts; the offer rule accepts R580, and a sanitizer start failure also loses D-G7/D-G8

Sources: 00_front (6.1 step 5, 6.4 box toolchain line), D (D8 `sanitize` and `--gpu-check` steps 4-6).

Doc:
- 00_front 6.1 step 5: "`offers.sh RTX_5090 --pick ID` | US/CA host, driver R580 or newer (CUDA 13)".
- 00_front 6.4: "The box toolchain for V2 and V3 is `cuda-toolkit[nvcc,cccl,nvtx]==13.0.3` and
  `nvidia-cuda-sanitizer-api==13.4.92`."
- D "D8", `--gpu-check` order: "4. `sanitize` (above); 5. the old-engine build for D-G7 ...; 6. the benches."

Measurements:
- The 13.4.92 sanitizer has run twice, both on driver 615.71.09: K1 (RTX 3090, `runs/k1/cuda_checker/RESULTS.md:87`,
  `memcheck2.txt` 0 errors) and R3 (RTX 5090, `runs/r3/cuda_checker/RESULTS.md:87`, `memcheck.txt` 0 errors).
- The two latest 5090 hosts had driver 580.159.03 (`runs/r4/env_check.json`, `runs/r5/env_check.json`). No sanitizer
  run exists on any 580 host.
- The K1 log shows how a sanitizer failure looks: "Target application terminated before first instrumented API
  call" (`runs/k1/cuda_checker/memcheck.txt`). That one came from the 13.0.85 wheel.

Risk: sanitizer 13.4 belongs to a newer toolkit than a 580 driver (CUDA 13.0). No NVIDIA page says it works there, and
no measurement shows it does. If it cannot start on the V3 host, D-G5 fails at step 4. Steps 5-6 come after it in the
same `set -euo pipefail` shell, so the D-G7/D-G8 benches never run. The rental minutes up to step 4 are lost.
No planned $0 check covers this, because the Mac has no GPU. It stays minor because there is no proof of failure.

Fix:
- D8 `--gpu-check`: after step 1 (toolkit install), add step 1b, a sanitizer smoke of about 10 s:
  `timeout 120 "$CUDA_HOME/bin/compute-sanitizer" --tool memcheck --error-exitcode 99 "$PY" -I -c "import numpy as np, objdel._gpucheck as g; print(g.device_count())"`.
  Run it after step 2's build so that a kernel loads. Write its exit code to `$OUT/sanitize.json`.
- D8: run `sanitize` as the last step (after steps 5-6). Record each tool's exit code and fail the mode at the end, so
  D-G7/D-G8 still produce numbers when the sanitizer cannot start.
- 00_front 6.1 step 5, V3 with item 1 approved: "prefer driver >= 615 (the only driver on which sanitizer 13.4.92 has
  run [measured: runs/k1, runs/r3 RESULTS.md]); on a 580 host, step 1b decides before any long step."

## 2. MINOR. C0 step 0 changes the box image, but no step rebuilds it, pushes it and updates the allowlist digest

Sources: C (C0 files and step 0), 00_front (6.1).

Doc: C "C0" step 0: "Box image: one nsys mechanism, A6's. In a new last layer of `Dockerfile.vllm_native` ... unpack
it with `dpkg-deb -x` to `/opt/nsys`. ... A6's run-time `nsys_capture.sh` then finds `/opt/nsys` and skips its
download." C.9 coverage row: "G-DK1 (`nsys --version`) | image build (CPU) and V2 item 1 | `env_check`".

Code: the box runs the image named in `.vast-allowlist` (`image = allow.get("image") or die(...)`
[code: experiments/phase2/vast/policy.py:333]). `build_images.sh` tags each image by a hash of the Dockerfile and
its inputs and says "record the pushed digest in .vast-allowlist" [code: experiments/phase2/docker/build_images.sh:14-15].
No phase source names `build_images.sh --push vllm_native` or an allowlist update (grep of `phases/*.md`: 0 hits for
`build_images`, `--push`, `allowlist` outside E's move tables).

Result: V2 boots the old image with no `/opt/nsys`. A6's script falls back to a 458 MB download inside the baseline
arm (A.md step 8 `get ...`). That costs minutes, not a result. Whether V2 item 1's G-DK1 `env_check` line then
fails is not said.

Fix: C0 step 0 add: "Mac, $0, before the V2 preflight: `build_images.sh --push vllm_native` (GHCR), then write the
pushed `<repo>:vllm_native-<hash>@sha256:...` into `.vast-allowlist` `image`; 6.1 step 4 (V2) checks that the
allowlist tag hash equals the hash of the committed `Dockerfile.vllm_native` and `locks/vllm_native.*`." Also say
that a first pull of the new 4.35 GB image counts inside V2 item 1's launch minutes.

## 3. MINOR. D-G13 uses the bare `dph_total`, while front section 6 prices every row at r_s; at r_s the D share can pass $0.85

Sources: D (D-G13, D.8), 00_front (6 planning inputs, 6.3 rows V2 item 4 and V3 item 1).

Doc:
- D "D-G13": "D's share = its work minutes ... x that rental's billed rate: <= $0.85 (§8)". D.8: "79 x 0.641 / 60 = $0.84".
- 00_front 6 "Session rate r_s": "$0.680/h = $0.01133 per minute ... Every $ in 6.2-6.3 uses it". The 6.3 rows price
  D's packages at r_s: V2 item 4 "33 | 0.37" and V3 item 1 "27-46 | 0.31-0.52".

Arithmetic: at r_s, D's 79 work minutes cost 79 x 0.01133 = $0.90 > $0.85. The front rows sum to 0.37 + 0.50 (44 of
the 46 min are D's) = $0.87. So the two sources give two answers for one gate. No money moves, because R-D1 is
dropped on the time-line order.

Fix: D-G13: "D's share = D's work minutes x r_s (front 6): <= $0.90", and D.8: "79 x 0.01133 = $0.90". Or keep $0.85
and say "bandwidth excluded, unlike front 6".

## 4. MINOR. D-G7's "about 6 min" of benches is not checked against its biggest workload on the old dense engine

Sources: D (D.8 R-D2 row, D-G7, D.9.7 B01).

Doc: D.8: "D-G7/D-G8 interleaved benches about 6 [estimate]". B01: "alternating old/new subprocess ..., 3 rounds x (2
warm-up + 10); workloads lesion512 b1/b256/b8192, 256 LIDC-like blobs, 37,093 crops, 4096x4096 sparse and dense,
100,000 8x8 masks, 1x4096".

Code and measurement: the old `check_masks` sweeps S = 4(W+1)(H+1) slots in 2L rounds per chunk
[code: gpucheck.cpp:931-937, 966-1057]. That is about 1.05 M slots and 44 rounds per 512x512 mask, so about
3.8e11 slot-rounds for b8192, per call [arithmetic]. The CPU build needs 15.3 ms per 512x512 mask
[measured: D.1 item 2], so one b8192 call is about 125 s on one core. B01 makes 36 old-engine calls of each workload. No
GPU number exists for the old engine on b8192. The bench steps have no time limit of their own; only the rental cap
stops them. An overrun shortens the B3 run: its `BudgetStop` is fixed from item 1's 46 min (front 6.3 "V3 amount
rule").

Fix: D-G7/B01: run b8192 for 1 round x (1 warm-up + 3) per engine, and record that in the report. Wrap step 6 in
`timeout 900`, and on a timeout write the partial table, so `--compare` reports "incomplete" and not a pass. D.8:
state the bench minutes as "<= 15 (timeout)", and front 6.3 V3 item 1 keeps 46 min as the cap of the row.

## 5. MINOR. D5 step 4's path-halving line is not valid CUDA and names no CPU-build path

Sources: D (D5 step 4).

Doc: "If so, use `atomicMin(par[x], par[par[x]])` (parents only decrease, so this stays a valid forest)."

Code: `atomicMin` takes an address. The engine already wraps it for both builds as
`atomic_min_old(i32* p, i32 v)` [code: gpucheck.cpp:233], and `unite` uses `atomic_min_old(par + a, b)` (:1278). The
`find` that D5 step 4 changes is `HD` (both builds) [code: gpucheck.cpp:1262-1265].

Fix: "inside `find`: `atomic_min_old(par + x, load(par + p))` before `x = p` (the existing `HD` helper, so the CPU
build compiles the same line); D-G6 (determinism) and I14 rerun after the change."

## 6. MINOR. D5 step 1 says "the kernel returns early", but the label pair is three kernels

Sources: D (D5 step 1).

Doc: "The final `label(1), label(0)` pair (:1614) runs per crop row only when `changed[o]` is set. ... The kernel
returns early otherwise, and the first-pass labels stay." Later: "The second-pair `LabelInit{c, 1}` does
`atomic_add(&relabelled[o], 1)`".

Code: `label(which)` launches `LabelInit`, `LabelMerge` and `LabelFlatten` [code: gpucheck.cpp:1597-1601]. If only
`LabelInit` returns early, `LabelMerge`/`LabelFlatten` still run on the kept labels. That is correct: the labels are
flattened min roots, so `unite` finds equal roots and returns (:1271-1283), and `find` returns `p`. But it saves only
one kernel's work. If all three return early, the result is the same, and the D-G8 "not worse" margin is easier to
meet. The text does not say which.

Fix: "All three kernels of the second pair (`LabelInit`, `LabelMerge`, `LabelFlatten`, both classes) return at row
entry when `changed[o] == 0`; only `LabelInit{c, 1}` counts the row into `relabelled[o]`. Test I14 checks the count."

## 7. MINOR. D-G7's "Outputs IDENTICAL" does not say which keys, and the new engine returns two keys the old one lacks

Sources: D (D-G7, D2 file list, D.9.7 B01).

Doc: D-G7: "... in alternating subprocesses, 3 rounds each of 2 warm-ups + 10 repeats. Outputs IDENTICAL." D2:
"`check_masks` ... return `ring_part` (int32 [R]) and `betti` (int64 [n, 2]). Keep the binding name and the old keys
(`quads`, `ring_start`, `pts_start`, `ring_area`, `pts`)."

Code: the `pre-d2` tree's `check_masks` returns only the old keys [code: cuda/gpucheck_py.cpp:152-156 builds the dict
from `quads` to `pts`; `grep ring_part\|betti cuda/gpucheck_py.cpp` = 0 hits, 2026-11-09]. A whole-dict equality in `--compare` fails with no defect. A compare over "the
common keys" passes, but then nothing checks that the new keys came back.

Fix: D-G7 and D8 step 6 `--compare`: "IDENTICAL = for each workload, the SHA-256 of `tobytes()` of `quads`,
`ring_start`, `pts_start`, `ring_area`, `pts` is equal across engines; the current engine's `ring_part` and `betti`
are checked by GPU02 (golden digests), not here."
