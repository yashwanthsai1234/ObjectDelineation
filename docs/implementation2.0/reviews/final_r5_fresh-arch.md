# Final review round 5: fresh outside architect

Reviewer: a fresh outside architect with no earlier opinion. Doc: `docs/implementation2.0/Implementation.md` (round 5,
15,116 lines), built from `docs/implementation2.0/*.md`. I read the whole front matter (sections 0-8) and B.1 and B1.5. I
spot-read A2, A3, A.6, C.8.2-C.8.3, D7, E "Order of execution", E6 and the 99_back money decisions. I checked the
claims below against the repo code (`policy.py`, `final_eval.py`, `train_rl.py`, `select_ckpt.py`, `serve.py`,
`train_sft.py`, `run_box.sh`, `gpucheck_py.cpp`, the docker locks), the R5 run files, PyPI and the import-linter docs.
I did not repeat items that `review-log.md` lists as fixed, or the round-4 minors (eval without a trainer, async
GRPO, the V1 tail slack, the A3 test-view bundle).

## Skills used

| Skill | What it changed in this review |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | The front already records a live catalog read (2026-11-08). No item here needs a new NVIDIA tool, so I made no new catalog lookup and recommend no install. It made me check one GPU claim on the web: "no GPU PNG decoder" (section 2.1). Official nvImageCodec sources confirm GPU decode only for JPEG, JPEG 2000 and TIFF, so the claim holds. |
| `engineering:code-review` | Each item quotes a code line or a measurement: `policy.py:58-61, 189-193, 321-326, 424-490`, `final_eval.py:66-81`, `select_ckpt.py:172-231`, `serve.py:159-177`, `train_rl.py:47-61`, `docker/locks/vllm_native.txt`, `runs/r5/eval/eval_*/final_eval.json`. |
| `engineering:debug` | I traced the default money path end to end: V1 totals, then the money left, then the V2 items. V2 cannot finish its first item (item 1). I also traced the V1 eval scoring path (`final_eval` -> `fork_pool` -> `_score` -> checker) on a CUDA build (item 6). |
| `engineering:system-design` | I judged the system as a whole. The data path, the gates, the seams and the rental order are sound. What I found sits at the edges: a rental the budget cannot fund, and plan limits that live only in prose. |
| `engineering:architecture` | I checked the decision tables (sections 3 and 4) and the layer contract. The import-linter docs say `exhaustive` "is only supported for layers that define containers". E0 sets `containers = ["objdel"]`, so the target contract is valid. vastai 1.11.0 is on PyPI. |
| `engineering:documentation` | I looked for one rule stated in two places. Found: V4 "inside the $3" against (b') "V4 waits" (item 1); FR1 against 6.3 (item 2); the 2.4 diagram against section 0 (item 3). |
| `huggingface-skills:huggingface-best` | No Hub leaderboard covers polygon text for CT. The model ranking in section 3 matches the HF configs the doc cites. I found no reason to change "keep 2B; 4B only behind a val gate". No item here. |
| `system-design-skills:system-design` | I stress-tested the plan with the "what if the default decision is taken" question. Under (b'), V2 spends money and returns no gate result (item 1). |
| `engineering:tech-debt` | Two plan constraints exist only in prose: the offer $/GB limit and the Mac/box version match. A top lab encodes both in code or in a test (items 7, 8). |

Web: import-linter layers docs (`exhaustive`, `exhaustive_ignores`, containers rule); PyPI JSON for `vastai`
(1.11.0) and `import-linter` (2.16); nvImageCodec docs and search (GPU codecs: JPEG, JPEG 2000, TIFF).

## Verdict

This is the right end-to-end architecture for this budget and this stage. Data comes first. One paid experiment
(V1) carries every Gate B1 input, and the owner approves its decision table in advance. One CUDA engine source serves
both the CPU build and the GPU build. Three deep seams carry the design: `geom.polygonize`, `policy.sampler` and
`rl.trl_adapter`. TRL sits behind one adapter. A typed run spec and an exhaustive import-linter contract arrive in E.
A top lab would not restructure it.

A top lab would do these things differently:

1. Never plan a rental that the default budget cannot fund. Under the default (b'), V2 cannot finish its own launch
   item, and V4 is said to fit "inside the $3" (item 1). Under $3, a top lab spends every paid minute on the quality
   question. C and D GPU checks wait for a budget raise: they save about 18 min per run.
2. Score in one process from V1 on, not from D7 on. Today V1's evals still fork a CUDA worker (item 6).
3. Put plan limits in code: the offer $/GB limit (item 7), and the Mac test stack equal to the box stack (item 8).
4. Already logged in round 4 (not counted again): a vLLM-only eval path before V3, and the async-GRPO revisit
   trigger.

## Items

### 1. BLOCKING. 00_front (6.3 "Money left after V1", V2 table, "V4", "Totals for owner decision 1"), 99_back (decisions 1, 3, 3b), E (E.8): under the default decision (b') "keep $3", V2 runs "with the money left", but that money cannot pay even V2 item 1. V4 is also said to fit inside the $3.

Doc lines:
- 00_front, 6.3: "Money left after V1 (reserve $0.25): at 27.5 s: $3.00 - $2.29 - $0.25 = $0.46, which is 31 min of
  V2 after its fixed bandwidth b. At 23.65 s: $3.00 - $2.10 - $0.25 = $0.65, which is 48 min".
- 00_front, 6.3 V2 table, item 1: "Launch, sync, `env_check`, toolchain and CUDA build ... plus the fixed bandwidth b
  | 40 | 0.55". Item 2: "13 | 0.15". Item 4 (R-D1): "33 | 0.37".
- 00_front, 6.3 totals: "(b') keep $3: V1, then V2 with the money left ($0.46-0.65); C adoption does not fit; V3 and
  V4 wait."
- 00_front, 6.3 "V4": "(optional, after E, owner decision; inside the $3 unless the owner raises the cap) ... about
  39 min = $0.54".
- 99_back, decision 1 (default `[(b')]`): "(b') keep $3. ... V2 runs with the money left; C adoption needs a budget
  raise." Decision 3: "[as stated; drop R-D1 from V2 if the money left is too small]". Decision 3b: "It counts inside
  the $3 unless the owner raises the cap." E.md E.8: "V4 counts inside the $3 unless the owner raises the cap."
- A.md A.6, A-G6: "Box half (open until V2 item 3; Phase A closes without it)".

Arithmetic, using the doc's own inputs (r_s = $0.01133/min, b = $0.10, section 6 table):
- At 27.5 s, $0.46 buys (0.46 - 0.10) / 0.01133 = 31.8 min. V2 item 1 alone needs 40 min ($0.10 + 40 x 0.01133 =
  $0.553). Item 1 does not finish.
- At 23.65 s, $0.65 buys 48.5 min: item 1 (40) plus 8.5 min. Item 2 (13 min) does not fit. Item 3 (37 min) and
  R-D1 (33 min) do not fit either. Dropping R-D1 (decision 3) frees nothing, because item 1 is already the limit.
- Item 1's 40 min is also below V1's own basis for the same work. V1 uses launch to first `run` up to 38.2 min
  [measured: ledger.csv], plus the build (3 min) and `env_check` (1 min), which is 42.2 min = $0.58.
- V4 ($0.54) does not fit $0.46. It fits $0.65 only if V2 does not run. So "V2 runs with the money left" and "V4
  counts inside the $3" cannot both hold. The 6.3 totals line says V4 "waits" under (b'). Decision 3b and E.8 say it
  fits inside the $3.

Code: the cap formula ends every plan job at the cap, `cap = int(min(a.cap_hours * 3600, run["deadline"] -
time.time() - MARGIN_S))` [code: experiments/phase2/vast/policy.py:480]. With an approved amount of $0.46-0.65, the
cap kills V2 inside item 1 or item 2. No `c_speed_check` result, no C0 arm, no A-G6 box half and no R-D1 bench comes
out.

What goes wrong: the owner takes the default (b'). The engineer then follows 6.3 and rents V2 for $0.46-0.65. That
money is spent with no gate result: "spends GPU money for nothing". Decision 3b and E.8 let the engineer instead plan
V4 inside the same $3. That is the second claim to the same money.

Exact fix:
- 00_front 6.3 totals, (b'): "keep $3: V1 only. V2, V3 and V4 wait for a budget raise. The money left
  ($0.46-0.65) stays unspent as a second reserve. No V2 subset gives a gate result: item 1 alone costs $0.55-0.58,
  and items 1-3 cost $1.12." Change the 6.3 heading "V2: validation with the money left" to "V2: validation after a
  budget raise".
- 00_front 6.3 V2 item 1: use V1's basis, 42 min = $0.58 (38.2 + 3 + 1, plus b). Recompute "Required for C
  adoption" (92 min = $1.14) and "Items 1-5".
- 00_front 6.3 "V4", 99_back decision 3b and E.md E.8: replace "inside the $3 unless the owner raises the cap" with
  "needs a budget raise: after V1, $0.46-0.65 is left, and V4 costs $0.54".
- 99_back decision 1 (b'): "V2 runs with the money left" becomes "V2 waits for a raise, like V3 and V4". In decision
  3, delete "drop R-D1 from V2 if the money left is too small".
- A.md A.6 A-G6 and C.8.2 "Budget": "open until V2 item 3" becomes "open until V2 runs (needs a budget raise)".
- A Mac test in `tests/test_repo_hygiene.py` (DOC4 drift family): parse the 6.3 tables. Assert that under (b') no
  rental's first item costs more than the money left.

### 2. MINOR. 00_front (section 7, FR1): stale V2 price.

FR1: "the V2 items for C adoption (103 min = $1.26) do not fit". Section 6.3: "Required for C adoption (items 1, 2,
3) | 90 | 1.12". Section 6 says that it wins, so no engineer acts on FR1's number. Fix: FR1 points to 6.3 and holds
no number: "the V2 items for C adoption (front 6.3) do not fit".

### 3. MINOR. 00_front (2.4 diagram): stale rows-root size.

The 2.4 diagram says "-> B1 rows root (9,600 rows)". Section 0 says "rows root ... (B1: 8,320 rows = 260 steps x
32". B.md B1.2 says "the B1 rows root holds 8,320 rows. The 9,600-window sizes and times in B1 are the 300-step upper
bound." Fix: write "8,320 rows" in the diagram.

### 4. MINOR. 00_front (line 4): stale revision line.

"Revision: final round 3 plus the round-4 consistency pass." Round-4 review fixes are in the doc, and this is the
round-5 review. Fix: "Revision: final round 4". Better: let `assemble.py` write the round number, so the line cannot
drift.

### 5. MINOR. 00_front (6.2 notes "Why items 3-5 always fit", "Why 260 steps"), 99_back (open question 13): the BudgetStop minute is wrong by arithmetic.

The doc says "(2.05 - 0.10) / 0.01133 = 172.3". The true value is 1.95 / 0.01133 = 172.1. With the exact r_s
(0.680 / 60) it is 172.06. So "260 steps end at minute 53 + 119.2 = 172.2, before the stop" is false at the planning
values. No result changes, though. `BudgetStop` checks at step end. After step 259 the spend is $0.10 + (53 + 259 x
27.5 / 60) x 0.011333 = $2.046 < $2.05. After step 260 it is $2.051. So the stop trips exactly when step 260
ends. Fix: write 172.1 everywhere (front 6.2 twice, front 6.2 "Why 260 steps", 99_back question 13, the "28.5 s"
row of the old log is not live). Then write "260 steps end at about the stop minute; `BudgetStop` checks after each
step, so step 260 completes".

### 6. MINOR (top-lab difference). A (A2 step 2, A3 `final_eval` change), 00_front (6.2 "Why step 0"): "pool size 1, so no forked worker opens its own CUDA context" is wrong. One forked worker still opens one.

Doc: A.md A3: "on a CUDA build ... `final_eval` caps the pool at 1 worker ... so no forked worker opens its own CUDA
context". Front 6.2 says the same.

Code: `fork_pool(n)` builds `_Pool(n, mp_context=get_context("fork"))` and forks n workers [code:
src/objdel/train/serve.py:164-177]. `score_items` maps `_score` over the pool, "one CUDA checker call per item"
[code: experiments/phase2/select_ckpt.py:172-176, 223-231]. With n = 1, one forked child runs every checker call
and opens its own context, tables and scratch, beside vLLM and the trainer. D7 step 6 states the real fix: "`metrics()`
device work runs in ONE process: the main process".

Effect: one extra context, not sixteen. V1 still runs. The stated reason is false, though, and D-G17's "no new
context appears" would fail on this path. Fix (small, and in V1's interest): in A3, when `g.build() == "cuda"`,
`final_eval` calls `_score` serially in the main process after generation, with no pool. Scoring takes 0.7 s (val)
and 3.1 s (test) today [measured: runs/r5/eval/eval_{val,test}/final_eval.json, `total - load_and_sample`], so
serial scoring is cheap. Otherwise change the text to "at most one forked CUDA context until D7".

### 7. MINOR (top-lab difference). 00_front (section 6 rules, "Offer limit"): the plan's $/GB limit lives only in prose.

Doc: "Offer limit: `dph_total` <= $0.641/h and `inet_up_cost`, `inet_down_cost` <= $0.0065/GB. A dearer offer needs
a new plan, because the `BudgetStop` values below assume these limits." Code: `INET_CEIL ... = 0.03` [code:
experiments/phase2/vast/policy.py:61], which is 4.6x looser. At $0.03/GB, b = $0.45 and r_s rises by $0.18/h. Then
the "$0.35 buys at least 30.9 min" tail claim fails (0.35 / ((0.641 + 0.18) / 60) = 25.6 min < 29 min). 6.1 step 5
asks a person to check "inside the offer limit". Fix: the 6.1 step 0 edit also adds a `plan_limits` argument to
`policy.reasons()` (or a check in `launch_refusal`). It refuses `dph_total > 0.641` or a $/GB above 0.0065 for a
launch whose config has `rows_plan`. Add one unit test.

### 8. MINOR (top-lab difference). 00_front (6.1 step 4, order rule 2 pre-tag tests), A (A6): the pre-tag Mac gate runs on a different stack from the box.

The pre-tag tests (the tiny TRL loop `tests/test_p2_train_rl.py`, the A7 check, P13/R01) run on the Mac venv. The
Mac has torch 2.14.1, transformers 5.20.1, trl 1.17.1 and peft 0.21.2 [measured: `.venv` dist-info]. The box image
has torch 2.13.0 and transformers 5.19.0 [code: experiments/phase2/docker/locks/vllm_native.txt]. C.md notes "Two
stacks". A.md A6 describes its `get_batch_samples` NVTX wrap from transformers 5.20.1 ("In transformers 5.20.1 it
fetches all 16 ..."). A difference in transformers internals would surface only on the paid box, after about 53 min
of item 1 ($0.70). A top lab runs the go/no-go gate on the deploy lock. Fix: add 6.1 step 4e. Build a CPU test venv
with `torch==2.13.0 transformers==5.19.0 trl==1.17.1 peft==0.21.2 accelerate==1.15.0` (the box pins without vllm and
the CUDA wheels). Run the tiny TRL loop, the A7 check and the A6 NVTX tests there before the owner tags `v1-b1`.

## Counts

Blocking 1 (item 1). Major 0. Minor 7 (items 2-8).
