# Final review round 5: CUDA and GPU performance (Implementation.md, built from docs/implementation2.0)

Reviewer role: CUDA and GPU performance (kernels, GPU placement, syncs, memory, the rental plan).
Scope read: 00_front sections 1, 3, 5, 6 (all of 6.1-6.4) and 7; A3 and A6 (capture script); B1.5 and B3.9; C.3, C0,
C1, C5, C6, C.8; D.3 to D.8 (all tasks and gates); E4 and E.8 (V4); CLEAN2 and the CLEAN placement table. Earlier
items were checked against `docs/implementation2.0/review-log.md` (round 4 cuda-1..cuda-4 are fixed and not repeated).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Catalog routing for CUDA / Nsight / compute-sanitizer. I installed nothing (no approval to install). I checked the NVIDIA facts the doc relies on by hand: racecheck covers shared memory only (web check below), `cuda-toolkit==13.0.3` has an `nvtx` extra (PyPI JSON: `nvidia-nvtx==13.0.85; extra == "nvtx"`), and both pinned Nsight Systems `.deb` URLs return HTTP 200 (458,836,100 B and 232,877,734 B, the same as the doc). |
| engineering:code-review | I read each code line the doc cites: `policy.py:379` (code sync drops `data/` and `docs/`), `policy.py:424-492` (`chain`, `cmd_run`, `-I` prefix), `bundle.py:48-70` (bundle members; `lidc-train-*` views only when named), `gpucheck.cpp:66-71, 78-98, 131-158, 177-184, 1020-1055, 1574-1639, 1802` (lanes, alloc, grow-only `Block`, `Drain`, waits, the label pair at :1614), `final_eval.py:69-78`, `select_ckpt.py:172-231`, `serve.py:132-180`, `train_rl.py:35-61`, `Makefile` (`make cuda`), `CMakeLists.txt`, `run_box.sh`. This produced items 1, 2 and 3 and minors 5 and 8. |
| engineering:debug | I traced each failure path to the first step that breaks: the R-D1 config on a box with no `docs/` tree (item 1); `c_speed_check` with a bundle that has no train views (item 2); a grow-only arena on a step with a longer answer (item 3); a V4 eval on the CPU build (item 4). |
| engineering:system-design | I recomputed the V1 timeline and money: cap 3.35 h = 201 min; `BudgetStop` at minute (2.05 - 0.0975) / 0.011333 = 172.3; the 26-min eval tail; the pull window. Arithmetic checks pass (minors 2 and 3 are wording and margin). |
| engineering:architecture | I checked one-rule consistency: front 6.4 rule 7 (what leaves the Mac) against C1 ("train views in the V2 bundle"); time line step 5 against 6.3 ("D lands after the B runs") against D-G14/15/17 (read in the B3 run). Result: item 2 and minor 7. |
| engineering:documentation | Each blocking or major item quotes the source file and heading, and the code line or measurement that proves it. Each item ends with an exact fix. |
| huggingface-skills:huggingface-best | I checked the model facts of front section 3 against the HF `config.json` of `Qwen/Qwen3-VL-2B-Instruct` and `Qwen/Qwen3-VL-4B-Instruct`: vision depth 24, width 1024, DeepStack [5, 11, 17], merger out 2048 vs 2560, text 28 x 2048 (k/v 8 x 128 = 1024) vs 36 x 2560. All match the doc and the C5 LoRA arithmetic (196 tensors, 1.41 B params, 17,432,576 adapter params). No finding. |
| huggingface-skills:huggingface-papers | No paper claim is in my area. Model-choice citations (2601.13633, 2604.27720) were not re-read. No finding. |

Web check (required for review roles): NVIDIA Compute Sanitizer docs say racecheck covers on-chip shared memory only.
So D8's claim "Racecheck cannot see global-memory races", and the determinism test it adds, are correct.
Sources: [NVIDIA Compute Sanitizer](https://developer.nvidia.com/compute-sanitizer),
[HiRace, arXiv 2401.04701](https://www.arxiv.org/pdf/2401.04701).

Measured facts I used (repo): `runs/k1/checker_gpu.log` (the toolkit swap 13.4 to 13.0 works with the image lock;
`make cuda` takes 15.7 s); `runs/rl/cuda_checker/pytest.txt` (`tests/test_gpu_scoring.py` on the box: 12 passed in
6.88 s); `runs/r5/eval/eval_*/final_eval.json` (`total - load_and_sample` = 3.1 s for test and 0.7 s for val, so
pool scoring is a small part of the eval minutes).

---

## Items

### 1. BLOCKING: R-D1 (V2 item 4) and K3 (V2 item 5) cannot run: no step ships the kernel sources or the K1 masks to the box

Source files: 00_front (section 5 "Custom kernels index"; 6.3 V2 table items 4-5; 6.4 rule 7), D (D1 Files and step 0;
D1 step 4; D.8).

Problem.
- 00_front section 5: "The same step adds an R-D1 config that runs the K1 and K2 runbook steps through `run.sh`
  (heartbeat and ledger rows)". D.md D1 step 0: "Create `experiments/phase2/configs/r_d1.yaml` ... It runs K1-draft
  and K2-draft runbook steps 1-5 ... through `run.sh` ... No manual `scp`, no toolchain install by hand."
- The runbooks live under `docs/implementation2.0/kernels/K1`, `K2`, `K3`, and they need two inputs from the Mac.
  (a) The kernel sources. K1 `GPU_RUN.md:51`: "Upload from the Mac: scp -P <port> /tmp/k1_bundle.tgz root@<host>:/root/".
  K2 `GPU_RUN.md:33`: "rsync -av --exclude build .../kernels/K2/ BOX:~/K2/". K3 `GPU_RUN.md:16`:
  "rsync -az --exclude build "$K3"/ root@$HOST:/root/K3/". (b) K1's raw LIDC masks: `lidc/in.bits`, `in.hw` and
  `in.start`, made on the Mac by `cross_numpy.py`. Front section 5 puts them in "`$REPO/data/k1_bundle/`
  (git-ignored, LIDC only)".
- The code sync drops both trees: `drop = re.compile(r"^(data|docs)/|...")` [code:
  experiments/phase2/vast/policy.py:379]. The LIDC bundle holds only `prepared/{lidc,images/lidc,caches/lidc}`, the
  views, the eval views and `checker_eval/` [code: experiments/phase2/vast/bundle.py:52-70]. It holds no `k1_bundle`
  and no `docs/`. DOC1 puts `docs/` under git, but `code_files()` still drops it.
- So the first R-D1 step finds no `docs/implementation2.0/kernels/K1` on the box, and it fails. The doc forbids the manual
  scp, so the engineer cannot proceed. If the step runs anyway, V2 bills item 4 (33 min, $0.37) and gives no D-G1 or
  D-G2 result. K3 (item 5, $0.17) fails the same way. Front 6.4 rule 7 names "the K1 mask bundle" among the artifacts
  that leave the Mac, but no step says how it leaves.

Fix.
1. D1 step 0: move the three drafts to a path that the code sync keeps, for example
   `experiments/kernels/K1`, `K2`, `K3`. They are experiment code that runs on a box, not docs. Keep `DESIGN.md`
   beside the code, and leave a pointer in `docs/implementation2.0/kernels/README.md`. Or add a named allowlist for
   `docs/implementation2.0/kernels/**` (no `build/`) to `code_files()`, with a unit test.
2. `r_d1.yaml` step 1 builds K1's raw masks on the box from the bundle: `k1_dump_cpu` / `cross_numpy.py` on
   `{data}/checker_eval/masks.npz`. That file is LIDC-only and already in the bundle [code: bundle.py:66-70], so no
   `data/k1_bundle` must ship. Remove "the K1 mask bundle" from 6.4 rule 7, or name the bundle member that carries it.
3. Each runbook drops its scp/rsync lines and uses paths relative to the synced folder (`/workspace/src/...`).
4. Add a Mac test (front 6.1 step 4, before V2): every path that `r_d1.yaml` and the K3 runbook step read is in
   `policy.code_files()` or in `bundle.members("lidc", ...)`. Preflight runs it when the config is `r_d1.yaml`.

### 2. MAJOR: C1's "train views in the V2 bundle" are never put in a bundle; `c_speed_check` (V2 item 2) has no shifted or rescaled views on the box

Source files: C (C1 "`c_speed_check` batches"; C.8.3 "Phase B1 data"), 00_front (6.3 "Own batches"; 6.4 rule 7),
A (A3 header comments, the one bundle command).

Problem.
- C.md C1: "Draw 5 x 32 views, seeded (seed 0), from the train views in the V2 bundle: centred, `lidc-train-shifted`,
  `lidc-train-rescaled`." 00_front 6.3: "`c_speed_check` (item 2) makes its 5 batches on the box as its first step,
  from the train views in the V2 bundle (C1)." G-PX1 / G-EQ1 also need these view kinds (C.9.4).
- The bundle leaves out the train windows unless they are named: `if kind == "views" and not views: # RL train
  windows ... only when named` and `names = {n for n in names if not n.startswith(f"{ds}-train-")}` [code:
  experiments/phase2/vast/bundle.py:55-56]. `--view` also replaces the default list (bundle.py:54, 121). So the
  command must also name `lidc-bundle-v3`.
- The only bundle command in the doc is A.md's `.venv/bin/python -I -m experiments.phase2.vast.bundle lidc --view
  lidc-bundle-v3`. It ships no `lidc-train-*` view. No V2 step names a bundle command.
- Front 6.4 rule 7 lists what may leave the Mac: "the B1 rows root, the K1 mask bundle, the LIDC eval views named in
  the rental's table, and the LIDC val views ...". The train views are not in the list.
- Result: an engineer who uses the doc's bundle command starts V2 with no shifted or rescaled train views. Item 2
  then fails or draws centred views only, and G-EQ1 does not cover the view kinds that B1 trains on (C.8.3's stated
  reason). Item 2 costs $0.15, and the C-adopt gates it feeds are lost.

Fix. Pick one source and write it in C1 and in V2 item 1.
- (a) Preferred: draw the shifted and rescaled views from the B1 rows root. It is LIDC-only, on HF, already allowed by
  rule 7 and frozen by its manifest SHA-256. V2 item 1 adds `fetch_hf.sh <run> --rows <rows_revision> lidc-train-b1`,
  and `c_speed_check --rows {data}/rows/lidc-train-b1` draws by the rows' `view` kind.
- (b) Or: V2 item 1 names the bundle command `bundle lidc --view lidc-bundle-v3 --view lidc-train-shifted --view
  lidc-train-rescaled`, rule 7 adds "the LIDC train views named in the rental's table", and a Mac test asserts that
  the V2 bundle manifest lists both train views.

### 3. MAJOR: Gate D-G15's box half has no valid pass rule: it is vacuous after D5c, or it fails on correct code

Source files: D (D.6 D-G15; D5 step 5 "D5c"; D1 "Engine counters").

Problem.
- D-G15: "Box, V3 B3 run: `gc.stats()` over steps 2-6 shows 0 `cudaMalloc`/`cudaFree`/`cudaMallocHost`/`cudaFreeHost`
  calls. `pinned_grows == 0` over the whole run".
- D5c changes device memory: "Device `alloc` calls `cudaMallocFromPoolAsync(&p, bytes, pool, lane_stream)`. `release`
  calls `cudaFreeAsync`". After D5c no device path calls `cudaMalloc` or `cudaFree`. Read literally, the device half of
  D-G15 passes by construction and measures nothing.
- Read as the D1 counter ("device `alloc`/`release` calls"), it can fail with no bug. The device arenas are grow-only
  and sized by the call: `if (bytes <= cap) return p; ... while (c < bytes && c < (size_t(64) << 20)) c <<= 1; ...
  p = alloc(kDevice, c)` [code: cuda/gpucheck.cpp:131-150; `Block<true> dev, score` at :154]. D5c gives a first-call
  bound only for the pinned blocks ("So the engine sizes each pinned block at its first call from a stated bound",
  D5 step 5). It gives none for `Block<true>`. Answer lengths change from step to step (38 % of R5 steps hit
  1,536 tokens, BRIEF section 2). So a step 2-6 batch that is larger than every earlier batch grows `dev` or `score`
  by one async pool allocation. That is not a sync, and it is not a defect, but the gate fails.
- The engineer must guess which counter the gate reads. One choice makes the gate untestable; the other can fail a
  correct D merge in the paid B3 run.

Fix. Rewrite the D-G15 box rule so that it tests the intent ("no synchronizing allocation after warm-up"):
"Box, V3 B3 run, steps 2-6: `gc.stats()` shows 0 synchronizing allocator calls. That is 0 `cudaMalloc`, `cudaFree`,
`cudaMallocHost` and `cudaFreeHost` calls, counted by new counters `sync_alloc` / `sync_free` at the only call sites
left (pinned `Block<false>` growth and `add_tables` fallbacks). `pinned_grows == 0` over the whole run. Device pool
allocations (`pool_allocs`) are recorded, not gated." Add `sync_alloc`, `sync_free` and `pool_allocs` to the D1
counter list. Also add a Mac test: on the CPU build, a call larger than the warm-up call raises `pool_allocs` by 1
and leaves `sync_alloc` at 0.

### 4. MAJOR: V4 has no plan file and no CUDA build step, so its val eval scores on the CPU build

Source files: 00_front (6.3 "V4"), E (E4 gate "statistical (TOST)"; E.8 "V4"; E-G10).

Problem.
- 00_front 6.3: "**V4** ...: E-G10 box part (prompt ids and SamplingParams equal; TOST sampler arm on val), E4
  deterministic check and E-G13. Work: about 9 min [estimate: E.8]". E.8: "the box runs one eval only: 1,084 items x
  K 8 took 415 s".
- Round 4 (cuda-2, qa-2) added a CUDA build step, and an assert, to every V1-V3 plan. The fix list in
  `review-log.md` names `b1_rows.yaml`, `a3_measure.yaml`, `c_speed.yaml`, `a3_test.yaml` and `b3_v5.yaml`. V4 is not
  in it. No V4 config file is named anywhere.
- `sync_up` installs the CPU build [code: experiments/phase2/vast/policy.py:399 `uv pip install ... -e .`;
  CMakeLists.txt:15 `option(OBJDEL_CUDA ... OFF)`]. Only `run_box.sh --build-only` builds the CUDA one. So the V4
  sampler arm, written from the doc, scores 8,672 answers on the CPU checker. That breaks the owner rule "CUDA is a
  must for the checker, the scorer AND the polygonizer". The 9-min work estimate also leaves out the build.

Fix. In E.8 and front 6.3 V4, name the plan file `experiments/phase2/configs/e_v4.yaml`. Its step 0 is the CUDA build
of `b1_rows.yaml` step 0 (`run_box.sh --build-only`, `assert g.build() == 'cuda'`, `tests/test_gpu_scoring.py`). Then
come `env_check.sh` (E-G13), the E4 deterministic check, and the sampler-arm `final_eval`. Price V4 at 30 + 3 + 9 =
42 min. Extend the 6.1 step 4d grep test and B I22 to assert that the V4 step 0 equals `b1_rows.yaml` step 0.

---

### 5. minor: A3's pool cap wording is wrong: a pool of size 1 still forks one worker that opens a CUDA context

Source files: 00_front (6.2 "Why step 0"), A (A3 `final_eval` changes).
00_front: "On a CUDA build `final_eval` uses pool size 1 ... so no forked worker opens its own CUDA context". The code
forks the worker anyway: `pool = serve.fork_pool(serve.pool_size())` and `score_items(pool, ...)` maps `_score` over
the pool. `_score` makes "one CUDA checker call per item" [code: final_eval.py:69, 77; select_ckpt.py:224-230]. So
one child opens one context (about 0.5 GB beside vLLM at 0.35). The cap behaves as intended (one context, not 16),
but the sentence is false. Fix: "so at most one forked worker opens a CUDA context (one, not up to 16); D7 step 6 later
moves the device work into the main process."

### 6. minor: the V1 eval tail leaves out the post-training merge and save, and the one-worker scoring

Source files: 00_front (6.2 Notes, "Why items 3-5 always fit").
The tail is "items 3-5 (20.8 min) plus about 5.2 min for the untimed imports ...: 26 min", with 198.3 <= 201 (2.7 min
of margin). Not counted: (a) `train_rl`'s save of `final_adapter/` (new in B1.5), then `merge_and_unload()` and
`save_pretrained(final)` of a 2B model after the last step [code: train_rl.py:55-60]; (b) scoring on one pool worker
in place of N (R5: 3.1 s for 2,200 x 8 answers on the pool [measured: final_eval.json total - load_and_sample]; about
1 min serial for the 3 V1 views [estimate]). Both are small [estimate: 1-2 min in all], but they use most of the 2.7
min. Fix: add a row "post-train save and merge, about 1 min [estimate]". Record `time/final_save_s` in `rl.json` so
that V1 measures it. Say that items 3-5 assume pool size 1.

### 7. minor: three texts disagree on when D lands relative to the B3 run, and the B3 code tree is not pinned

Source files: 00_front (6.3 text under the V2 table; section 1 time line step 5), D (D-G14, D-G15, D-G17; D.8).
6.3: "R-D1 serves Phase D, which lands after the B runs." Time line step 5: "C1-C5, B2/B3 and D2-D8 code land on the
Mac", before V2 and V3. D-G14/15/17 are read in "V3, first 6 steps of the B3 run". R-D2 is "Skippable". Nothing says
which engine the B3 run uses if R-D2 is skipped or fails. `b3_v5.yaml` step 0 runs `tests/test_gpu_scoring.py`, which
guards the reward path on the GPU, so the risk is low. Fix: delete "which lands after the B runs". Add to B3.9: "the
B3 run uses the D engine only if R-D2 passed in the same V3; otherwise it runs on the pre-D2 tag, and D-G14/15/17 stay
open."

### 8. minor: CLEAN's placement table says the eval pool runs CPU numpy/scipy, but each worker also calls the CUDA checker

Source files: CLEAN (placement table, row "`serve.fork_pool` in `select_ckpt`, `g_probe`, `final_eval`").
"CPU (stays for now) | It runs the eval-only metrics (NSD, GED, HM-IoU, soft Dice) in numpy/scipy." `_score` ->
`metrics()` makes one CUDA checker call per item [code: select_ckpt.py:172-176, 224-225]. On a CUDA build every
worker opens a context. A3 caps only `final_eval`, and `select_ckpt` / `g_probe` keep `pool_size()` up to 16 [code:
serve.py:132]. No V1-V4 plan runs them on a box, so no money is at risk. Fix: correct the row, and add "on a CUDA build
`select_ckpt` and `g_probe` also cap the pool at 1 (A3's rule) until D7 step 6".

### 9. minor: the V1 worst case uses $2.38; for a cheaper offer it is up to $2.51

Source files: 00_front (6.2 Notes "Worst case").
"session spend at the cap is at most $2.40 ($2.38 at planning values) ... So the worst case is $2.49." The cap is
floored from the picked offer, so spend at the cap can come close to $2.40 for other offers. With up to $0.11 for
`PULL_S`, the worst case is $2.51. `$16.90 + $2.51 = $19.41 <= $19.90`, and with the reserve `$2.76 <= $3.00`, so no
rule changes. Fix: write "$2.51 at most ($2.49 at planning values)".

### 10. minor: stale numbers after the 300 to 260 step change and the V2 re-pricing

Source files: 00_front (section 7 FR1), D (D7b "Why"), C (C.8.2 "Payback").
FR1: "V2 items for C adoption (103 min = $1.26)" against the 6.3 table "90 | 1.12". D7b: "300 steps x 32 prompts =
9,600". C.8.2: "a 300-step B run ... 118 min". Fix: FR1 to 90 min = $1.12. D7b to "260 x 32 = 8,320". C.8.2 to 260
steps (102.5 min at 23.65 s, 86.7 min at 20.0 s, about 16 min saved per run).

### 11. minor: D-G7 does not say where the old-engine tar and venv go

Source files: D (D.6 D-G7; D8 `--gpu-check` step 5).
"`git archive -o old_engine.tar <pre-D2 sha> ...` The tar ships in the code sync" and "`make -C <dir> cuda
PY=/opt/old/bin/python`". The output path and `<dir>` are not named. A tar under `build/` (git-ignored, `/build/` in
`.gitignore`) would not ship, because `code_files` uses `git ls-files -co --exclude-standard` [code: policy.py:377].
`uv venv /opt/old` does not pin Python, and `CMakeLists.txt` needs `Python 3.12` with `Development.Module`. Fix:
"write `old_engine.tar` at the repo root (untracked, not ignored; S06 asserts `old_engine.tar in
policy.code_files()`); extract it to `/opt/old_src`; `uv venv --python 3.12 /opt/old`."

### 12. minor: D2's "IDENTICAL outputs" for D-G7 needs the compared keys

Source files: D (D-G7; D2 Files).
D2 adds `ring_part` and `betti` to `check_masks`. The pre-D2 engine returns neither. "Outputs IDENTICAL" should say
"on the keys both engines return (`quads`, `ring_start`, `pts_start`, `ring_area`, `pts`)".
