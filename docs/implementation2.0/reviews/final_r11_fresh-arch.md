# Final review round 11: fresh outside architect

Scope: the whole doc `docs/implementation2.0/Implementation.md` (round 11), read by sections. My focus is the end-to-end
architecture up to Phase E: the money path of front 6.2-6.4, the gate design of B.6.2, D3, E3.3 and E4, and the
owner decisions in 99_back. I checked every proof against the code, a measurement file or the web. I do not repeat
items that `review-log.md` (rounds 1-10) or the r10 verdicts already handle. In particular, the sealed V1 test eval,
the (b') sequencing, the D box-gate asymmetry and the doc size are not raised again.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item depends on an NVIDIA product choice. The front already lists the catalog skills for the 3D direction. I made no new catalog lookup and installed nothing. I did check one NVIDIA claim on the web (nvImageCodec has no GPU PNG decoder, front 2.1); it holds. |
| engineering:code-review | I read `experiments/phase2/train_rl.py:44-68` (no fallback after a failed `trainer.train`), `final_eval.py:46-86` (scoring is outside `seconds.load_and_sample`), `vast/policy.py:183-199, 308-335, 418-510, 560-620` (`chain`, `cmd_run` cap, `accrued`, sweep), `cuda/gpucheck_py.cpp:140-158` (`check_masks` returns rings: `pts`, `ring_start`, `ring_area`), `reward/score.py:100-124` and `targets/build.py:53, 93`. |
| engineering:debug | I traced two failure paths of V1 to the end: (a) the host dies mid-run; (b) `train_rl` exits non-zero at step N < 260. Neither path leaves a Gate B1 input or a resumable checkpoint on the Mac (item 2). |
| engineering:system-design | I checked whether each gate pairs its two arms under the same conditions. Gate B1 does (A3 re-evaluates R5 on the same box). Gate B3 does not: it pairs a V3 eval with a V1 eval, across code trees (D7, E4) and hosts (item 1). |
| engineering:architecture | I weighed "seam first, engine later" for the polygonizer. The existing engine already gives byte-identical rings, so the owner-rule exception of ADR-0015 is not needed (item 3). ADR-style options and costs are in each fix. |
| engineering:documentation | Every major item quotes the source file and heading, plus the code line or measurement that proves it. |
| huggingface-skills:huggingface-best | No item is about model choice. TRL 1.17.1 exists on PyPI (current release) [web: pypi.org/pypi/trl/json]. |
| system-design-skills:system-design | "Design for failure": the one paid rental has no story for a lost host (item 2). "Quantify before choosing": I priced each fix at r_s = $0.01133/min, b = $0.10. |
| engineering:tech-debt | Stale references left by the order of waves (items 5, 6, 7). |
| WebSearch / WebFetch / Exa | import-linter 2.16 uploaded 2026-11-06 (E.2 claim holds) [web: pypi.org/pypi/import-linter/json]. verl v0.7.1 lists "LoRA adapter only refit (TensorLoRARequest)" as a megatron PR #4632, so "on FSDP not verified" in front 4 is correct [web: github.com/verl-project/verl/releases/tag/v0.7.1]. nvImageCodec decodes PNG only through its CPU OpenCV extension [web: github.com/NVIDIA/nvImageCodec]. |

## Overall verdict (architecture)

The data-first order, the single source for rental money, the session-spend `BudgetStop`, the one-source CUDA engine
and the paired patient CIs are what a strong lab would do. Three things a top lab would do differently, each with
evidence below: (1) pair every gate's arms in one process on one box, not only Gate B1's; (2) keep the only paid
run recoverable, with checkpoints off the box during the run; (3) put the seam in place first with the engine that
already passes, instead of asking the owner to waive a binding rule until a new kernel lands.

## Items

### 1. major: Gate B3 pairs a V3 eval with a V1 eval, across two code trees and two hosts; on the planned post-E path the two arms even use different samplers

Sources: B, 00_front, E.

Problem.
- Gate B1 is built on same-box pairing. B "B.6.2 Run gates", Paired baseline: "Primary: A3's R5 val items (front
  6.2 item 5: ... the same code and the same box as the B1 eval). Then a view or host difference cannot enter the
  B1 - baseline difference."
- Gate B3 drops that rule. B "B.6.2", Gate B3: "after the run, on val-select (re-dealt halves, B1.0), paired against
  the B1 run". 00_front "6.3", V3 item 4: "then B3 on `lidc-val-shifted` and `lidc-val-edge` in one process
  (10.4 min)". No V3 item re-evaluates B1. So the B1 arm is V1 item 3 (`eval_b1_val_*`), made months of commits
  earlier on another host.
- Between V1 and V3 the eval code changes on `main`. D lands after V1 (00_front time line step 5) and D7 moves the
  eval pixel work into the engine: D "D.3.2 Data flow" item 3, "Eval (box, in the `final_eval` process): parse ->
  `score_many` (terms) plus `eval_sets` calls". Its GPU check (R-D2) is optional in V3 (00_front 6.3 item 1:
  "Optional: runs only if the owner names it").
- The doc plans a V3 after E. 00_front "End state under (b')": "If E0 has started, the owner rebases the branch onto
  post-E paths ... If E6 has also landed ... ports the branch plans and launch configs to `specs/`". 00_front 6.3
  "After E6": "The V3 lines change the same way, to `specs/b3_v5.yaml`". On that tree `final_eval` no longer uses the
  trainer path: E "E4" step 3, "`final_eval` calls `VllmSampler` + `evaluate`. ... It no longer calls `build_trainer`
  or `trainer._generate`." The check that the new sampler equals the trainer path is box-only and lives in V4: E
  "E4" Gate, "E4 is not done until E-G10's box part passes. No B run follows E in the front time line, so the box
  check needs its own owner-approved rental, V4". That sentence contradicts the post-E V3 path quoted above.
- No planned $0 check covers it. E4's Mac tests compare `build_prompt` ids with TRL's and `params_for` with
  `policy.models.sampling`; the comparison of `SamplingParams` with TRL's `VLLMGeneration` is the box part of E-G10
  (vLLM is not installed on the Mac).

Result. An engineer who follows the doc under a raise after E4 computes Gate B3 = (B3 with the new sampler, new
engine, host Y) minus (B1 with the trainer sampler, old engine, host X). Any sampler difference (penalties, stop ids,
prompt ids) is credited to reward v5, and the final-candidate choice "Best of B1 and B3 on val-select by the Gate B3
rules" (00_front 1.1) reads the same confounded difference. Before E the confound is smaller (host and D7 only) but
the design still breaks the rule that Gate B1 states.

Exact fix (B B.6.2 Gate B3, B3.9 and I22b; 00_front 6.3 V3 table, amount rule, launch lines, decision-1 totals;
E E4 Gate):
- Add V3 item 3b, "B1 re-eval": one `final_eval` process on `{rd}/b1/final` (already fetched by item 0) on
  `lidc-val-shifted` and `lidc-val-edge`, `--out {rd}/eval_b1v3_val_shifted` and `--out {rd}/eval_b1v3_val_edge`,
  as a step of `b3_v5.yaml` before `train_rl` (so the tail stays 26 min). Minutes: 205 s + 2 x 1,084 x 0.193 s =
  10.4 min + 2.6 min imports = 13.0 min = $0.15 [estimate: front 6 inputs].
- Gate B3 pairs B3 with item 3b, not with V1 item 3. V1 item 3 becomes a report line, as R5's `items.jsonl` is for
  Gate B1.
- V3 amount rule with items 0, 3b, 4: step 1 starts at minute 39 + 13 + 14 = 66, step 260 ends at minute 185.2, so
  `BudgetStop` = $0.10 + 185.2 x $0.01133 = $2.20 (rounded up) and the approved amount is $2.55. B-S2 = 195.6 min =
  $2.32 at 27.5 s. Decision 1 (a): $2.29 + $2.32 + $0.25 = $4.86 at 27.5 s.
- E4 Gate: replace "No B run follows E in the front time line" with "A B run follows E only on the post-E raise
  path; there Gate B3 pairs within V3 (item 3b), so the sampler change cannot enter it. V4 still closes E-G10."
- I22b: `b3_v5.yaml` holds the item 3b step before `train_rl`, with `--ckpt {rd}/b1/final`.

### 2. major: the only paid run has no off-box checkpoint and no partial-result path, so a lost host or a crash at step N < 260 spends about $2 for no Gate B1 input

Sources: 00_front, B.

Problem.
- 00_front "6.2", launch lines: one `pull.sh`, at the end only ("`pull.sh <run> && push_hf.sh ... && teardown.sh`").
  Nothing moves checkpoints off the box during the 2-3 h run.
- B "B1.5", Resume: "Save adapter, optimizer and step every 50 steps ... `train_rl --resume` passes the last
  checkpoint". The repo docstring states its purpose: "so --resume continues from the last one after a lost box"
  [code: experiments/phase2/train_rl.py:8-9]. After a lost box the checkpoints die with it, so `--resume` has nothing
  to read. The runbook recovery assumes pulled checkpoints: "Box lost during RL | ... sync the old `train_rl/`
  folder (the pulled checkpoints) up ... At most 50 steps are lost" [code: docs/runbooks/gpu-run.md:149].
- The tool exists and is not used: `com.objdel.vast.pull.plist` runs `pull.sh` every 1,200 s [code:
  experiments/phase2/vast/com.objdel.vast.pull.plist:3, 11]; `pull.sh` with no run id pulls every active run [code:
  policy.py `cmd_pull`, `for run_id in [a.run_id] if a.run_id else list(active())`]. It is not loaded on this Mac
  [measured: `launchctl list | grep objdel` empty, 2026-11-09]. The money is already planned: `UP_GB_PER_H = 6` is
  "checkpoint pulls of 3-5 GB per 20 min" [code: policy.py:59-61], and r_s includes it (front 6 inputs).
- A crash has no fallback. `b1_rows.yaml` runs its steps joined by `&&` [code: policy.py:487-494], and `train_rl`
  writes `final/` only after `trainer.train` returns [code: train_rl.py:55-60]. So a non-zero exit of step 3 skips
  step 4, the B1 eval, although `trainer/checkpoint-200` exists. 00_front 7 has no row for a host lost mid-run (FR12
  covers only "failed launch, slow image pull"), and B.7 has no row for a training crash.

Result. Under (b') nothing follows V1. A host loss or a crash at, say, step 210 (session minute about 150, $1.80
spent) leaves no B1 eval and, after a host loss, no checkpoint on the Mac. Gate B1 cannot be computed, and the
decision table 1.1 has no row for "no B1 output".

Exact fix (00_front 6.2 launch lines and notes, 7 risk register, 1.1; B B1.5 and B.7):
- 6.2: before `run.sh <run> .../b1_rows.yaml`, the main agent loads the pull job (`sed "s#__REPO__#$REPO#g"
  experiments/phase2/vast/com.objdel.vast.pull.plist > ~/Library/LaunchAgents/com.objdel.vast.pull.plist &&
  launchctl load ...`) and unloads it after `teardown.sh`. A note: the 30-min idle destroy after a pulled `EXIT`
  [code: policy.py:603] is why `a3_measure.yaml` starts within minutes of the `b1_rows` `EXIT`, as today.
- `b1_rows.yaml`: make step 3 soft (`["~", -m, experiments.phase2.train_rl, ...]`) and add step 3b,
  `train_rl --merge-last {rd}/b1-rows` (new flag, B1.5): if `final/` is missing, merge the last
  `trainer/checkpoint-N` adapter into `final/` and record `steps: N, partial: true` in `rl.json`. Step 4 then always
  has a model. U-test: a stub run that raises after step 2 still leaves `final/` and `partial: true`.
- 1.1: add a row "B1 stopped early (`partial: true`, N < 260 steps)": Gate B1 is computed and labelled "partial
  run, N steps"; a pass is not used for B-S2 until the owner reviews it.
- 7: add FR16 "Host lost mid-run": Medium / High; mitigation = the 20-min pull job (checkpoints on the Mac) plus
  `--resume` in a new rental after a raise; under (b') the money left ($0.46-0.65) does not buy a resume, and the
  owner is told that before V1.

### 3. major: ADR-0015 asks the owner to waive "CUDA for the polygonizer" for the B1 rows root "until K1 passes parity", but today's engine `check_masks` already gives byte-identical rings on every LIDC mask, on the CPU build and on a real RTX 5090

Sources: B, 00_front, D, 99_back.

Problem.
- B "B.3.3": "**Owner exception needed:** the targets use the numpy polygonizer, against the rule 'CUDA for the
  polygonizer'. The exception is recorded in ADR-0015 and ends when Phase D K1 passes parity". 99_back decision 9:
  "B1.4 targets use the numpy polygonizer until D K1 passes parity. Then B1.4 calls the CUDA polygonizer. [accept]".
  00_front 2.2 seam row repeats it. BRIEF section 3 makes the rule binding: "CUDA is a must for the checker, the
  scorer AND the polygonizer."
- The premise is wrong for the engine that exists now:
  - CPU build: "k1_dump: gpucheck::check_masks IDENTICAL / cross_numpy: 37093 masks, 38002 rings, 0 mismatches ->
    PASS" [measured: docs/implementation2.0/kernels/K1/results/cross_numpy_m3.txt]. `cross_numpy.py` compares "ring
    count, every ring's corners and signed area" against `objdel.checker.polygonize` [code: K1/cross_numpy.py:1-8].
    So `check_masks` equals the numpy polygonizer ring for ring on all 37,093 LIDC masks.
  - GPU: R3 on an RTX 5090, "Equivalence with the CPU checker's recorded results | 0 mismatches: ... 37,093 LIDC
    masks" [measured: docs/research/notes/issues_solved_2026-11-07.md:37].
  - Speed on the Mac CPU build: `check_masks` 4,833.9 ms for all 37,093 real masks [measured:
    K1/results/bench_cpu_m3.txt, "real (all)"], against "about 160 s with the numpy path" for 24,624 rows (D "D3"
    step 3).
  - The binding already returns rings: `check_masks` gives `pts`, `ring_start`, `pts_start`, `ring_area` [code:
    cuda/gpucheck_py.cpp:140-158]. Only tests and `experiments/cuda_checker/compare.py:76` call it today (D.1 item
    1).
- So the exception buys nothing: D3's seam (`polygonize(mask, offset)` over `_gpucheck.check_masks`) can exist before
  B1.4, with today's dense engine behind it, and D2 later swaps K1 in behind the same seam under D3 step 3's hash
  check. The doc instead gives the owner an option list without the compliant $0 option.
- The end condition is also undefined: "K1 passes parity" may mean the D2 CPU parity (Mac) or the GPU run (R-D1 or
  R-D2, which never run under (b')).

Result. If the owner accepts decision 9 as written, the one paid run trains on targets made against a binding owner
rule, on a premise that the repo's own measurements contradict. If the owner refuses it, the plan has no written
alternative and V1 is blocked.

Exact fix (B B.3.3, B.3.4 ADR-0015, B.8.3 item 2, B1.4; 00_front 2.2 seam row, order rule 2; D D3; 99_back
decision 9):
- New step B1.4a, before B1.4 (on the `v1-b1` tag; order rule 2 lists it as the one D-seam change): add
  `src/objdel/checker/rings.py` with `polygonize(mask, offset=(0, 0))` that packs the mask, calls
  `_gpucheck.check_masks` once and adds the offset on the host; `targets/build.py:93` calls it; the cv2 hole
  grouping in `_rings` stays as today. The numpy `polygonize.py` moves to `tests/oracles/` in the same commit.
- Gate: D3 step 3's `verify_target_hashes.py` (moved to B1.4a) rebuilds every recorded target set on the Mac CPU
  build; pass = every `target_set_sha256` byte-identical. Cost: seconds of Mac time, $0.
- D3 then reads: "D2 swaps the engine inside `check_masks` for K1; `rings.py` gains `trace` and `ring.part` (the cv2
  grouping goes); step 3 re-runs the same hash check." ADR-0015 drops "and the numpy-polygonizer exception";
  B.8.3 item 2 and 99_back decision 9 are deleted (keep the numbers of the other decisions).

### 4. minor (top-lab design): the GPU pixel path (C1) is tied to a speed gate that (b') never funds, so `main` keeps the CPU pixel path although the owner rule is about placement, not speed

Sources: 00_front, C.
00_front section 1: "B uses the C code only when C-G-EQ1..6 pass and C-G-SP1 passes", and "End state under (b')":
"today's CPU pixel path only, with no C flag". The owner rule (BRIEF 3) is "GPU everywhere; CPU only where
unavoidable". A top lab would adopt a placement change on exactness (C-G-EQ, item 2 of V2, 13 min = $0.15) plus
"no slower than baseline", and keep C-G-SP1 as the speed claim. The doc's choice is coherent and stated, so this
is minor. Fix (owner option, not a default): in 99_back decision 1, offer "(b'') C1 adoption on C-G-EQ1 alone in
the first raise", with its price ($0.58 launch + $0.15 = $0.73, front 6.3 item 1 and 2).

### 5. minor: three 99_back and diagram lines disagree with front 6.3 on the D rentals

Sources: 99_back, 00_front.
- 99_back Open questions row 12: "K1-K3 GPU results ... V2 item 4 (R-D1: K1, K2)". Front 6.3 item 4: "on the
  time-line order it is dropped". Fix: "R-D2 (V3 item 1, optional) covers K1 and K2 inside the ported engine; R-D1
  only on a tree before D2".
- 99_back decision 3: "R-D2 (`make gpu-check`) is the first box item of V3". Front 6.3 item 1: "Optional: runs only
  if the owner names it in the V3 approval". Fix: add "optional, only if named in the V3 approval".
- 00_front 1 diagram: "V3: R-D2 first, optional fp16 arm, B3 run". Fix: "V3: optional R-D2, optional fp16 arm, B3
  run".

### 6. minor: D3 names scripts that CLEAN5 archives in wave 1, before D starts

Sources: D, CLEAN.
D "D3" Files: "`experiments/phase2/lab/B3/crop_trace.py:43` and `experiments/phase0/{lidc_census,riga_census,
riga_outlines}.py` ... they must still import", and the `--help` check extends to "`experiments/phase0/*.py`".
CLEAN "CLEAN5" (Wave 1): "`experiments/phase0/*.py` ... -> `archive/phase0/`" and "`experiments/phase2/lab/` ->
`archive/phase2_lab/`". Wave 1 runs in time-line step 2; D runs in step 5. The glob then matches nothing (no
failure, a dead line). Fix: drop the phase0 and lab names from D3 and say "archived by CLEAN5 before D3".

### 7. minor: E.3.3 assumes C has landed before E, which is false under the default (b')

Sources: E.
E "E.3.3": the host-sync baseline is "recorded on the E0 commit (after C has removed syncs)", and the device table
says "Phase C removes the CPU work left in the step". Under (b') C stays on `after-raise` (00_front "End state under
(b')"). The gate itself still works (it records whatever count E0 has). Fix: "recorded on the E0 commit (with C's
syncs removed only if a raise merged `after-raise` before E0)"; same wording for the device row.

## Counts

blocking 0, major 3, minor 4.
