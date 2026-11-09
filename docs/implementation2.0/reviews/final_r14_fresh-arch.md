# Final review round 14: fresh outside architect

Reviewer role: a fresh architect with no prior opinion. Scope: the whole plan `docs/implementation2.0/Implementation.md`
(round 14), read by sections. The sources in `phases/` are the files to fix. Focus: is this the best end-to-end
architecture up to Phase E with GPU everywhere, and what would a top lab do differently. I opened the code that each
item cites. Before writing I read `docs/implementation2.0/review-log.md` (latest rounds of every phase) and the r13
fresh-arch review and verdict. I do not re-raise r13 fresh-arch-2 (test read under (b'), downgraded), -3 (pool size 1
wording), -4 (front header round), -5 (E5 row), or r12 fresh-arch-2/-3/-5/-6.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Loaded. The cached live catalog (`impl/nv_catalog_live.json`, read 2026-11-08) still lists `compute-sanitizer-diagnose-cuda-bugs`, `nv-segment-ct`, `dicom-series-to-volume`. No item here needs a new NVIDIA skill. Nothing installed. |
| engineering:code-review | Read `policy.py` (`chain` :422-461, `cmd_run` :464-495, `cmd_push_hf` :511-524, `CARDS` :48-55, `INET_CEIL` :61, `reasons` :122-124), `final_eval.py`, `train_rl.py:55-60`, `select_ckpt.py:54-96`, `sample_vllm.py` docstring, `run_box.sh`, `tests/test_gpu_scoring.py`. Found items 1, 3, 4, 7, 8. |
| engineering:debug | Traced the V1 failure path "crash before the first checkpoint" from `b1_rows.yaml` step 3 to step 3b to section 1.1 (item 2); traced the end-of-V1 push from `push_hf.sh` to `hf_push` (item 4). |
| engineering:system-design | Re-ran the V1 money arithmetic for the retry case of item 2 against `launch_refusal` and `HARD_STOP` $19.90. |
| engineering:architecture | Wrote the overall judgement as an ADR review (seams, layers, (b') end state, data-first order). |
| engineering:documentation | Every blocking and major item quotes the source file and heading plus the code line. Items 3 and 7 are wording or estimate gaps. |
| huggingface-skills:huggingface-best | No model choice is in question (2B now, 4B gate on val). Used for the HF side of item 6 (free private storage on the Hub). |
| system-design-skills:system-design | Ran the loop on "the default (b') end state after E": requirement (no old code beside new code), constraint (no rental after V1), failure mode (V4-gated archive never fires). Found item 1. |
| engineering:tech-debt | Scored what stays on `main` under (b') after E: two dead sampler scripts (item 1), a refactored trainer never run on a GPU (item 5). |

Web checks:
- Hugging Face free accounts: 100 GB private storage; PRO 1 TB [web: https://huggingface.co/docs/hub/storage-limits].
- A forum report of a 403 "Private repository storage limit reached" while the quota page showed space
  [web: https://discuss.huggingface.co/t/403-error-private-repository-storage-limit-reached-quota-shows-space-remaining/153121]
  (anecdotal; used only for minor item 6).

## Overall judgement

The architecture is sound. A small top-lab team with a $3 cap would build something close to it:
- Data first is the right answer to the measured diagnosis ("trained centred, tested shifted"). The scheduled rows
  root has 21 gates on the Mac and a pinned HF revision plus a manifest SHA.
- One CUDA engine source builds for the CPU and for nvcc. The CPU build is the oracle. `tests/test_gpu_scoring.py`
  runs as step 0 of every paid plan, so a bad GPU build stops in minutes. I checked that this test covers the whole
  reward path (`ExpertCache` -> device tables -> score, against `oracles/score_numpy.py`) [code:
  tests/test_gpu_scoring.py:1-5, 114-260]. So D code that lands on `main` without R-D2 still meets a GPU parity check
  before any B3 training.
- The three seams (`geom.polygonize`, `policy.sampler`, `rl.trl_adapter`) sit at the right places.
- The (b') end state is honest: the GPU pixel path stays off `main` until a raise.

What a top lab would do differently. None of these is a defect, except where an item says so.
1. Archive superseded code in the same commit that switches the callers, and keep rollback in a tag, as the plan
   already does for the trainer eval path (`pre-e4-eval`). The plan instead keeps two dead sampler scripts until
   a V4 that (b') never funds (item 1).
2. Save checkpoints at the data block size (20 steps) and decide the "no checkpoint" case before the one paid run
   (item 2). Early steps are where new code fails.
3. Run a 2-step GPU training smoke after a large refactor of the training path (E5 splits `fast.py`) before the next
   paid run trusts it (item 5).
4. Ship the V1 runbook as one short page (commands, gates, stop rules). Today it is spread over front 6.2, B1.5, A3
   and 1.1 inside a 17,050-line plan (no item; style).
5. Long term, move from materialised PNG windows to on-device random windows (ADR-0015 option (b)). At the 3D voxel
   step, materialising every scheduled sub-volume will not scale (no item; ADR-0015 already names the trigger).

## Items

### 1. BLOCKING. Under (b'), E keeps the old samplers `sample_vllm.py` and `hf_sample.py` beside `policy/sampler.py` forever. This breaks the owner archive rule and contradicts the front's "End state under (b')"

Sources: E (E4 step 1 "Archive after the box check"; E4 Gate; E.6 E-G10, E-G14; E.8), 00_front (time line "End state
under (b')"; 6.3 "V4"), CLEAN (CLEAN8 step 4 note), 99_back (decision 1 (b'), decision 3b).

Problem.
- 00_front, "End state under (b') (the one rule; ...)": "`main` then never holds a temporary flag, a second reward or
  a second pixel path, and no old code sits beside new code (owner archive rule)." Its list "Lands on `main` under
  (b')" ends with "then E, the DOC final sweep and CLEAN waves 2-3".
- 00_front 6.3: "Rule: under the default (b') ... V2, V3 and V4 each wait for a budget raise." 99_back decision 1
  (b'): "V2 (...), V3 (...) and V4 each wait for a budget raise." [(b')] is the default.
- E4 step 1: "**Archive after the box check** (owner rule: no second vLLM path): `experiments/phase2/sample_vllm.py`
  and `experiments/phase2/hf_sample.py` go to `archive/phase2_samplers/` ... This is a follow-up commit, made only
  after E-G10's box part passes (V4, E.8). Until then they stay".
- E.6 E-G14: "Until the V4 follow-up commit (E4), `sample_vllm.py` and `hf_sample.py` are the listed exceptions".
  00_front 6.3 "V4": "Until V4 passes, E4 keeps the old samplers (`sample_vllm.py`, `hf_sample.py`)." CLEAN8 note:
  "`hf_sample.py` is not in this list: E4 archives it with `sample_vllm.py` (E.8), after E-G10's box check."
- After E4 step 1 the two scripts are dead code. Today their only live caller is `select_ckpt.sampler_cmd`
  [code: experiments/phase2/select_ckpt.py:93-96: `str(ROOT / "experiments/phase2/sample_vllm.py")` and
  `-m experiments.phase2.hf_sample`]. E4 step 1 rewrites it: "`sft.select.sampler_cmd` now builds `[sys.executable,
  "-I", "-m", "objdel.policy.sampler", "--kind", kind, ...]`". E4 also moves their bodies into `policy/sampler.py`
  ("Carried over from `sample_vllm.py` into this module, with their tests"; `HfSampler` is "the body of
  hf_sample.py"). So after E4, `request_seed`, `check_memory`, `ENV` and the HF generate body exist twice. The old
  copies have no caller and no test.
- The plan already handles the same case well for the trainer eval path. E4 step 3 says: "the owner tags the last
  commit that still has the trainer path `pre-e4-eval`. The trainer path then lives in git history only, also before
  V4". The samplers get the opposite rule, and no reason is given for the difference. The fallback reason ("so a
  working vLLM path exists if the new sampler fails on the box") does not need the files in the working tree. They
  are in `pre-e4-eval`, and no caller uses them after E4 step 1 anyway.

Result. An engineer who follows E exactly under the default (b') ends with dead duplicate code on `main` next to
its replacement. That breaks the binding owner rules "Old code must be ARCHIVED (not left beside new code)" and
"Remove dead code" (BRIEF section 3). The front says "the one rule" holds no old code beside new code, and E's gate
table says the opposite. No planned check catches it: E-G14 names the two files as allowed exceptions, and E-G10's box
half shows "open until V4" for good.

Fix (E, 00_front, CLEAN, 99_back):
1. E4 step 1: replace the "Archive after the box check" paragraph with: "`experiments/phase2/sample_vllm.py` and
   `experiments/phase2/hf_sample.py` go to `archive/phase2_samplers/` (README: Was / Replaced by / Commit) in the
   E4 step 3 commit, the one that switches `final_eval` and `sampler_cmd`. Rollback is the tag `pre-e4-eval` (owner
   command: `git checkout pre-e4-eval -- experiments/phase2/sample_vllm.py experiments/phase2/hf_sample.py`). V4
   still closes E-G10's box half; if V4 fails, the owner restores them from the tag."
2. E-G14: delete "Until the V4 follow-up commit (E4), `sample_vllm.py` and `hf_sample.py` are the listed exceptions;
   after it they exist only under `archive/`". Write instead: "`sample_vllm.py` and `hf_sample.py` exist only under
   `archive/phase2_samplers/` from the E4 step 3 commit on". Extend U2/`test_old_packages_gone` to assert both paths
   are absent outside `archive/`.
3. E.8 and E Gate text: "E4 is done, and its follow-up commit archives ..., only after V4 passes" becomes "E4's
   code is done at the E4 step 3 commit; only E-G10's box half waits for V4".
4. 00_front 6.3 "V4": replace "Until V4 passes, E4 keeps the old samplers (`sample_vllm.py`, `hf_sample.py`)" with
   "The old samplers are archived at E4 step 3 and kept in the tag `pre-e4-eval`."
5. CLEAN8 note (CLEAN.md lines 767-768 and 1824): "E4 archives it with `sample_vllm.py` in the E4 step 3 commit".
6. Move any `parity.py` / `r0_parity.py` reliance on `sample_vllm.py --logprobs` into the same archive commit, or
   name their owner (CLEAN). Both read `sample_vllm` output only [code: experiments/phase2/parity.py:5, 71, 89].

### 2. MAJOR. A V1 crash before step 50 leaves no checkpoint, no B1 model and no row in the decision table; the plan does not say what follows

Sources: B (B1.5 `b1_rows.yaml` steps 3 and 3b), 00_front (1.1 decision table; 6.2 plan files; 7 FR16), 99_back
(decision 2b).

Problem.
- B1.5 step 3: "`... --steps 260 --save-steps 50 --shared-prefix --budget-usd 2.05`". Step 3b: "With no checkpoint
  (a crash before step 50) it exits non-zero, and the chain stops."
- 00_front 1.1: "The owner approves this table with the V1 item list ..., so the result cannot bend it." The only
  early-stop row is "B1 stopped early: `rl.json` has `partial: true` (a crash, then `train_rl --merge-last`, 6.2) and
  N < 260 steps | Gate B1 is computed on the merged checkpoint-N". A crash before step 50 writes no `partial: true`,
  because `--merge-last` has nothing to merge (I49: "on a folder with no checkpoint it exits non-zero").
- FR16 covers only a lost host or a crash with a checkpoint: "step 3b ... merges the last checkpoint into `final/`".
- Code: `train_rl` writes `final/` only after `trainer.train` returns [code: experiments/phase2/train_rl.py:55-60],
  and the `&&` join stops the chain [code: experiments/phase2/vast/policy.py:490].
- Steps 1-49 are the likely window for a new-code fault on the tag. Several paths run on a GPU for the first time
  there: the schedule-root loader and its asserts, `shuffle_dataset=False` with the sampler-shape check, the B2.1
  device stats, the session `BudgetStop` and the per-cell logs. The Mac fake-vLLM loop is the only earlier check.
- Money at that point: item 1 (53 min) plus up to 49 x 27.5 s = 22.5 min of steps, plus the `a3_measure.yaml`
  process (10.4 min + imports), comes to about $0.10 + 88 x $0.01133 = $1.10 [estimate: section 6 inputs]. A retry
  fits `launch_refusal` only if $16.90 + $1.10 + amount + $0.11 <= $19.90, so amount <= $1.79 [estimate]. That
  is below the approved $2.40 V1, so a retry needs a new owner approval and a shorter run.

Result. On the most likely failure, the engineer has no Gate B1 and no rule. They must guess whether to stop, to
re-launch with the remaining money (what amount, how many steps, which reserve), or to wait for a raise. Each choice
changes the outcome, and the owner approved the table believing it covered every result.

Fix (B1.5, 00_front 1.1 and 6.2, 99_back decision 2b):
1. B1.5 step 3 and I22a: `--save-steps 20` (one 20-step block; adapter plus optimizer for 17,432,576 parameters is
   about 0.2 GB, well inside the pull job's 6 GB/h). Step 3b's text: "With no checkpoint (a crash before step 20)".
   After that, only a crash in steps 1-19 (at most 9 min of steps) has no model.
2. 00_front 1.1: new row "B1 produced no checkpoint (crash before the first save; `--merge-last` exited non-zero)
   | No Gate B1. Mac diagnosis from `b1_rows.log` with the CPU build. A retry is a new rental that the owner
   approves with its own amount A, where $16.90 + session spend of V1 + A + $0.11 <= `HARD_STOP`. Its `--steps` is
   the largest multiple of 20 that fits (A - $0.35 - $0.10) / $0.01133 - 53 min at 27.5 s. Under (b') the reserve
   is spent first. The A3 answers from items 4-5 are kept | none | R5". Copy the summary into 99_back decision 2b.
3. 00_front 6.2 FR16 and B.7 crash row: add "a crash before the first checkpoint: 1.1 row 'no checkpoint'".

### 3. MINOR. I22a says `b1_rows.yaml` has "exactly the five steps of B1.5 (steps 0-4)", but B1.5 defines six (0, 1, 2, 3, 3b, 4)

Sources: B (B.9.3 I22a; B1.5).

I22a: "`b1_rows.yaml` has exactly the five steps of B1.5 (steps 0-4) in that order (step 0 ...; step 1 ...; step 2
...; step 3 ...; step 4 ...)". B1.5 adds "3b. `train_rl --merge-last ...` (a hard step)", and I47 parses "step 3b
(`train_rl --merge-last`) is hard". One of I22a or I47 fails on the Mac at $0 whatever the engineer writes, so it is
found during implementation (caught-before-spend rule). Fix: I22a "exactly the six steps of B1.5 (0, 1, 2, 3, 3b,
4)", plus "step 3 starts with `~`; step 3b is `train_rl --merge-last --model qwen3-vl-2b --ckpt {rd}/ckpt --out
{rd}/b1-rows` with no `~`".

### 4. MINOR. Pre-E, nothing says how `push_hf.sh` comes to call `hf_push --rl`, and no $0 test covers the call

Sources: B (B1.5 "B1 adapters to V3"; I41), 00_front (6.2 "The launch lines need five code changes").

B1.5: "B1.5 adds `hf_push --rl <train_rl out>` ... `push_hf.sh <run> b1_launch.yaml` runs it on the box before
teardown". Today `cmd_push_hf` builds `hf_push --run /workspace/runs/<run>/<output_dir> --repo ... --revision
<output_dir>` with no `--rl` [code: experiments/phase2/vast/policy.py:522-523]. The `--rl` rule only appears after E6
(E6 step 8b `push_target`). I41 tests `hf_push --rl` itself, not the command that `cmd_push_hf` builds. If the
engineer forgets the dispatch, the push fails at the end of V1, and `pull.sh && push_hf.sh && teardown.sh` skips the
teardown. No model is lost: `pull.sh` already copied `b1-rows/final` [code: policy.py:80 `PULL_EXCLUDE`], and the push
can be redone at $0. Fix: B1.5 "`cmd_push_hf` passes `--rl` when the config has `rows_plan` (a `train_rl` output)"; new
unit test I41b with `remote` mocked: `push_hf b1_launch.yaml` builds `hf_push --rl --run /workspace/runs/<run>/b1-rows
--revision b1-rows`, and `r1_qwen3vl2b.yaml` still builds today's command.

### 5. MINOR. After E, the refactored RL training path (E5 splits `fast.py`) has never run on a GPU, and V4 does not train

Sources: E (E.9.5 GP1-GP5; E.8 V4 steps), 00_front (6.3 "V4").

V4 runs `env_check.sh`, GP1-GP4 (scorer, sampler, imports) and one sampler eval. No V4 step runs `train_rl`. E's
move-invariance goldens run on the Mac with a tiny model and without vLLM. So the first GPU contact of the moved
shared-prefix forward, vLLM sync and caches is the next paid training run after a raise, after its 39-53 min fixed
start. Fix (optional, owner choice, priced in front 6.3): add a V4 step "`train_rl` 2 steps on the B1 rows root,
`--steps 20 --budget-usd <V4 stop>`, read only step 1 equality against a Mac-recorded step-1 golden (as C's
`--expect-step1`)". About 6 min start-up plus 1 min, so about $0.08 [estimate: section 6 inputs]. If the owner declines,
DOC11 Known issues states it.

### 6. MINOR. V1 pushes about 6 GB to a free HF account (100 GB private cap), and no preflight checks the quota

Sources: 00_front (4 "Data versioning"; 6.1 step 4; 6.2 item 2), B (B1.4 step 5; B1.5).

V1 uploads the rows root under `rows/` before the rental and pushes `final/` (a merged 2B model) plus
`final_adapter/` after item 5. The account is "storage only" [code: experiments/phase2/hf_push.py:8-9]. Free private
storage is 100 GB [web: huggingface.co/docs/hub/storage-limits]. I could not measure today's repo size without the
token, so this stays minor. Fix: a 6.1 step 4 check, $0: `HfApi().list_repo_refs` plus `get_paths_info` (or the
account's storage page) shows used + 7 GB < 100 GB; record the value in the stamp notes.

### 7. MINOR. The V1 eval-tail estimate leaves out `train_rl`'s merge and save of `final/`

Sources: 00_front (6.2 notes "Why items 3-5 always fit").

"The eval tail is items 3-5 (20.8 min) plus about 5.2 min for the untimed imports ...: 26 min". After `BudgetStop`,
`train_rl` still runs `merge_and_unload()` and `save_pretrained` of the 2B model and processor [code:
experiments/phase2/train_rl.py:58-60], and B1.5 adds the `final_adapter/` save before it. The margin is 201 - 198.3 =
2.7 min at planning values, or 30.9 - 29 = 1.9 min at the worst offer. If the save takes longer, the cap cuts item 5
first, and Gate B1 falls back to R5's `items.jsonl` (6.2 fallback). Fix: measure the save time on the Mac tiny-model
path or from the R5 log, and add it to the 26 min. If the margin goes below 0, set `--budget-usd` to $2.03.

### 8. MINOR. The plan's offer limit ($0.641/h, $0.0065/GB) is a manual check; the code admits $1.50/h and $0.03/GB

Sources: 00_front (6 "Rules"; 6.1 step 5).

00_front 6: "Offer limit: `dph_total` <= $0.641/h and `inet_up_cost`, `inet_down_cost` <= $0.0065/GB. A dearer offer
needs a new plan, because the `BudgetStop` values below assume these limits." The code lets through RTX 5090 offers up to
`"ceiling": 1.50` [code: experiments/phase2/vast/policy.py:51] and `INET_CEIL = 0.03` [code: policy.py:61, 123-124]. A
dearer offer that slips past the manual check of 6.1 step 5 makes `BudgetStop` trip before step 260. The run then
ends without `partial: true`, so no 1.1 row labels the shorter run. Fix: the 6.1 step 0 edit also adds
`policy.launch_refusal` reasons "dph_total > 0.641" and "inet_*_cost > 0.0065" when the config has `rows_plan` (V1, V3).
`rl.json` records `budget_stop_step` when the guard trips, and the 1.1 row "B1 stopped early" also matches
`budget_stop_step < 260`.

## Counts

blocking 1 (item 1), major 1 (item 2), minor 6 (items 3-8).
