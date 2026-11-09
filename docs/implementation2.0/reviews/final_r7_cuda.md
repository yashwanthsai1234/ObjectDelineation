# Final review round 7: CUDA and GPU performance (Implementation.md)

Reviewer role: CUDA / GPU performance (kernels, GPU placement, syncs, memory, rental plan). Sources read by heading:
00_front sections 0, 1, 2.4-2.5, 4, 5, 6.1-6.4, 7; D.md D.1-D.8 (D1, D2, D5, D5c, D6, D7, D7b, D8, gates, cost);
C.md C1, C4, C5; A.md A3 (pool cap) and A6 (NVTX fork safety); E.md E.8 (V4). Review logs of rounds 1-6 were read
so that no fixed item is repeated (r6 cuda-1..7 are fixed or were filed as minor; not repeated).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill covers this review. The front already cites `compute-sanitizer-diagnose-cuda-bugs`. I did not install or recommend anything new. The sanitizer chain (D8) was checked against the compute-sanitizer facts it cites (racecheck covers shared memory only; leaks are reported only at context destroy). |
| engineering:code-review | Read `cuda/gpucheck.cpp` (alloc/release :78-84, `Block` :135-150, `Drain` :177-184, `raster_kernel` :418, `find`/`unite` :1262-1281, `score_words` label lambda :1597-1616, waits :1020/1027/1055/1574/1626/1639/1802), `CMakeLists.txt`, `Makefile` `cuda`, `run_box.sh`, `reward/cache.py` (imports, `device()` :241-259), `reward/score.py:104-120`, `train/serve.py:132-180`, `final_eval.py:46-90`, `train_rl.py:30-70`, `fast.py:440-475`. |
| engineering:debug | Traced each box gate to the process and venv that runs it. This found item 1: the old-engine venv has only numpy, but the D-G8 replay needs `objdel.reward.cache` (cv2, scipy). It also found item 2: the device-side label skip is invisible to a host launch counter. |
| engineering:system-design | Checked the rental order against the phase landing order under (b'). This found item 3: V4 is the first paid run of the D2-D8 device code, and it has no D box gate. |
| engineering:architecture | Checked one rule per fact. The D5c pool design (own `cudaMemPoolCreate` pool, no `cudaDeviceSetMemPool`, `cudaFreeAsync` in lane order) is sound. `Lane` priorities and the two-lane `Drain` match the code. |
| engineering:documentation | Each blocking or major item quotes the doc line (source file and heading) and the code line that proves it. |
| huggingface-skills:huggingface-best | No model choice is in this scope; no change. |
| huggingface-skills:huggingface-papers | The D.2 papers are cited by venue, not by new arXiv claims. No paper fact decides an item; no change. |

Web checks (WebSearch, WebFetch): CUDA 13.x driver floors (13.0 GA needs driver >= 580.65.06; 13.x minor-version
compatibility starts at R580), which agrees with front 6.1 step 5 ("driver R580 or newer"). Measured R580 boxes ran
the stack (runs/r2/box.txt 580.178.04, runs/r4/box.txt 580.159.03). I could not confirm a separate driver floor for
compute-sanitizer 13.4.92 on R580, so I file no item on it.

Verified as correct (no item): the V1 CUDA build step (`run_box.sh --build-only` -> `make cuda` -> `build() ==
'cuda'` assert [code: Makefile:41-45; run_box.sh:31]). The compiler-wheel swap evidence [measured:
runs/k1/checker_gpu.log:12-22]. `score_polygons` makes 4 `wait()` calls plus 2 `Drain` syncs [code: gpucheck.cpp:1574,
1626, 1639, 1802; 177-184]. `ImageMemo` keeps a reference to each image, so an `id()` is never reused [code:
fast.py:469-471]. The C1 divisor rule. The C5 parameter counts (1.41 B touched, 0.31 B embedding). The D5c bounds.
The V1 cap of 3.35 h and the $2.49 worst case.

## Items

### 1. MAJOR: the D-G8 replay cannot run in the old-engine venv, which holds only numpy (D)

- Doc, D.md "D.6. Acceptance gates", row D-G7: "`uv venv /opt/old && uv pip install --python /opt/old/bin/python
  numpy && make -C <dir> cuda PY=/opt/old/bin/python` ... `bench_checker.py` imports only numpy and
  `objdel._gpucheck`, so the old venv runs the current script by path. `bench_checker.py --engine old` runs with
  that venv's Python".
- Doc, row D-G8: "Same setup as D-G7, on the D1 replay set (256 completions + caches), 50 repeats per engine ...
  Second replay with cold caches (tables dropped before each repeat): new median <= 1.05 x the old `add_tables`
  path". D8 `--gpu-check` step 6: "the benches: `bench_checker.py --engine old|current` (D-G7) and the replay
  comparison (D-G8)". D1 step 2 defines the replay as "their target texts ... scored on their caches from
  `data/prepared/caches/lidc`".
- Code: `score_polygons` needs table handles from `cache.device()`, and that call computes `add_tables` inputs from
  an `ExpertCache` [code: src/objdel/reward/score.py:112-115; reward/cache.py:241-257]. `reward/cache.py` imports
  `cv2` and `scipy` at module level [code: src/objdel/reward/cache.py:84-86]. The old venv gets only numpy and a
  `--no-deps` editable install [code: Makefile:42]. So `--engine old` on the replay raises `ModuleNotFoundError: cv2`
  at the last step of `run_box.sh --gpu-check`, in paid R-D2 time (V3 item 1). The cold-cache arm ("old
  `add_tables` path") needs the same host table build. `make gpu-check DRY=1` builds only the tar and never runs the
  old arm, so the Mac check does not catch this.
- Fix (D.md D-G7, D-G8, D8 step 6, D1 file list for `bench_checker.py`): add `bench_checker.py --dump-replay OUT`.
  It runs in the main venv on the box, before the old arm. It writes one `replay_inputs.npz` with these arrays:
  the `pack_objects` output of the 256 target texts (the exact `score_polygons` arguments except the handles), and
  per cache the 13 `add_tables` argument arrays as `ExpertCache.device()` passes them today. It also writes the
  SHA-256 of the file. `--engine old|current --replay replay_inputs.npz` then imports only numpy and `_gpucheck`.
  It calls `add_tables` from the dumped arrays, then `score_polygons`. The cold arm times `add_tables` plus
  `score_polygons` on the dumped arrays for both engines. The new engine also times `add_tables_from_raters` on the
  dumped rater planes. The host table build is the same Python code for both engines, so it stays out of the
  comparison. Add one S-test to the dry run: `--engine current --replay` runs on the Mac CPU build in a venv with
  only numpy plus the package installed with `--no-deps`.

### 2. MAJOR: the `label_pairs` counter counts launches, so it cannot see the device-side label skip that I14 and GAP-5 claim it checks (D)

- Doc, D.md "D.4. Tasks", engine counters: "`label_pairs`: launches of the scorer's second label pair (D5 step 1,
  test I14)."
- Doc, D.md "D5", step 1: "The final `label(1), label(0)` pair (:1614) runs per crop row only when `changed[o]` is
  set ... The kernel returns early otherwise ... The flag lives in the arena. There is no host read, so no new sync."
- Doc, D.9.3 I14: "the second label pair runs only on rows with `changed[o]` (the `label_pairs` counter, D1)". GAP-5
  (D.9.9): "`label_pairs` ... is now in the D1 engine counters ... so I14 checks the skip directly".
- Code: `label(which)` always enqueues three kernels over all `crows` rows [code: cuda/gpucheck.cpp:1597-1601:
  `for_each(st, crows, LabelInit{...}); ... LabelMerge ...; ... LabelFlatten ...`], and line 1614 calls it
  unconditionally. Under D5 step 1 the skip is decided per row on the device, so the host launch count is the same
  whether the skip fires or not (2 pairs, 6 launches, every call). A launch counter cannot tell "skip taken" from
  "skip broken". I14's direct check and the D5 step 1 stale-flag test ("the same skip pattern (`label_pairs`) as
  fresh calls") are untestable as written, so an engineer must guess what to count. D-G8 also forbids any new wait.
- Fix (D.md D.4 counter list, D5 step 1, I14, GAP-5): redefine `label_pairs` as "crop rows relabelled by the second
  pair". The CPU build increments it inside the `LabelInit` functor when `changed[o]` is set. The CUDA build does
  `atomicAdd` into one `i64` in the `score` arena, zeroed with `changed`. That `i64` is placed in the `ObjOut` block,
  so the existing D2H copy before the `wait` at gpucheck.cpp:1626 brings it back. This adds no new wait, so D-G8's
  count of 4 still holds. I14 then asserts `label_pairs == sum of (bh + 2) over objects with changed[o]`, and 0 on
  the "size filter removes nothing" item.

### 3. MAJOR: V4 is the first paid run of the D2-D8 device code under (b'), but it has no D box gate; E-G10 is read from that unverified device path (00_front, E)

- Doc, 00_front section 1, "End state under (b')": "Lands on `main` under (b'): ... D1-D9 (D's box gates stay
  'open until a rental') ... then E". 00_front 6.3 "V3" item 1: "R-D2: `make gpu-check` ... only after D2-D8 merge".
  D.md D.8: "R-D2 (D2-D8 ...) is V3 item 1".
- Doc, 00_front 6.3 "V4": "Step 0 is the CUDA build ... Then come `env_check.sh` (E-G13), the E.9 box test list, and
  the sampler-arm `final_eval`". E.md "E.8", step 2 lists only `tests/verify/test_verifier_gpu.py
  tests/infra/test_box_smoke.py tests/policy/test_sampler_box.py`. E.8: "The trainer arm is already on disk
  (`runs/r5/eval/eval_val/items.jsonl`)".
- Proof: under (b') with a raise that funds V4 alone, V3 never runs, so R-D2 (GPU02 CPU-vs-GPU bit equality,
  sanitizers, determinism) never runs. V4's `final_eval` then scores with D7's `eval_sets`, D7b's
  `add_tables_from_raters` and D2's engine on a GPU for the first time. GPU02 is the only GPU equality proof for
  `eval_sets` (D.md D-G9: "No box half ... whose device output GPU02 proves equal to the CPU build"). V4's step 0
  runs only `tests/test_gpu_scoring.py` (the scorer, not `eval_sets` or the device tables). The TOST compares this
  arm with R5's on-disk items, which the old engine and the old host metrics scored. So a GPU-only defect in D's
  kernels moves the sampler arm, and E-G10 passes or fails for a reason outside E. The doc does not say whether to
  run `make gpu-check` first. That choice changes the outcome.
- Fix (00_front 6.3 "V4", E.8, D.8, 99_back decision 3b): add V4 step 0b, "`d_gpu_test.yaml` (R-D2, `run_box.sh
  --gpu-check`) when no earlier rental has passed R-D2; stop V4 on any failure (6.4 rule 4)". V4's `sync_up.sh <run>`
  already unpacks the LIDC bundle (default `--bundle lidc`), which is the R-D2 payload [code: policy.py:604,
  406-417]. Reprice V4: 54 + 25-44 = 79-98 min, which is $1.00-1.21 at r_s with b. Set the approved amount and the
  cap from that. Add the step to U51. D.8 then names "first rental after D8 merges (V3 item 1, else V4 step 0b)".

### 4. minor: BudgetStop minute arithmetic in V1 (00_front)

- 00_front 6.2 notes say "`BudgetStop` at $2.05 ends training at session minute (2.05 - 0.10) / 0.01133 = 172.3" and
  "260 steps end at minute 53 + 119.2 = 172.2, before the stop". The arithmetic gives 1.95 / 0.01133 = 172.1. So at
  27.5 s per step and a 53-min start, the stop can cut the last step (about step 259). The tail still fits
  (172.1 + 26 = 198.1 <= 201). Also "3.39 h" should read 3.38 h (2.30 / 0.680); the 3.35 h cap is correct.
- Fix: write 172.1, and say "the stop may end the run at step 259-260 at exactly 27.5 s; the adapter is saved". No
  money changes.

### 5. minor: stale B1 row counts (00_front, D)

- 00_front 2.4 diagram: "B1 rows root (9,600 rows)". D.md D7b: "The B1 run sees each prompt once (300 steps x 32
  prompts = 9,600 [arithmetic])". Section 0 "rows root" and B1.2 say 8,320 rows = 260 x 32.
- Fix: write 8,320 rows (260 steps x 32) in both places.

### 6. minor: D-G13 prices D's share without the bandwidth part of the session rate (D)

- D.md D.8: "79 min stays at or below $0.85 [estimate: 79 x 0.641 / 60 = $0.84]". Front section 6 inputs: "Session
  rate r_s ... $0.680/h ... Every $ in 6.2-6.3 uses it". At r_s, 79 min = $0.90.
- Fix: state D-G13 at r_s (bound $0.90), or say that D's share counts `dph_total` only and that bandwidth stays in the
  shared rental.

### 7. minor: "pool size 1 (as A3)" is stale for V4 (00_front, E)

- 00_front 6.3 "V4" ("the sampler-arm `final_eval` (pool size 1, as A3)") and E.8 step 3 ("on the CUDA build with
  pool size 1 (as A3)"). Under (b'), D7 step 6 has already landed before E. D7 step 6 removes the pool from
  `final_eval` ("`final_eval`, `select_ckpt` and `g_probe` stop creating the pool for metrics").
- Fix: write "scoring in-process (D7 step 6)".

### 8. minor: the `contexts` counter is undefined in a process where torch already owns the primary context (D)

- D.md D.4 counter list: "`contexts`: 1 after the first engine call that opens a CUDA context in this process". The
  engine uses the runtime API, so it attaches to the primary context that torch or vLLM already created; it does not
  open one. D-G17 (b) reads the counter only in pool workers, where 0 is right, so no gate result changes.
- Fix: define it as "1 after the first engine call that makes a CUDA runtime call in this process (`e.start()`),
  else 0".

Counts: blocking 0, major 3, minor 5.
