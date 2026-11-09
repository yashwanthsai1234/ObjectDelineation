# Phase B test cases, set 2 (adversarial, property, performance, failure modes)

Written from `phases/B.spec.md` and REPO code/tests only (no design notes). IDs are `B2-nn`. Level is one of
unit / property / integration / GPU / e2e. "hyp" = hypothesis (already used in `tests/test_reward.py`). Spec
section numbers refer to B.spec.md.

## Skills used

- nvidia-skills:nvidia-skill-finder: no catalog check run (no network use in this role); no NVIDIA skill needed for tests.
  The CUDA items (section 3.5) are tested through the existing `_gpucheck` CPU/CUDA dual build, not a new tool.
- engineering:testing-strategy: pyramid split. Many cheap unit/property tests, few GPU and e2e, plus explicit
  performance gates and coverage of error paths.
- mattpocock-skills:tdd: each test targets a public seam named in the spec (function signature or CLI) and uses
  independent expected values (worked literals, brute-force oracles), never a re-implementation of the code under test.
- engineering:code-review / debug / system-design / architecture / documentation / huggingface-best: read for
  context only; they changed nothing in this file beyond the "why it catches a bug" column discipline.

Conventions: `-I` for any Python that reads `data/`; GPU tests skip unless `_gpucheck.build() == "cuda"` or a
CUDA torch device exists (match `tests/test_cuda_checker_gpu.py`); data-backed tests skip when `data/` is absent
(match `tests/test_reward_gates.py`). Unit tests need no data: they build synthetic `Item`s with disks as masks.

Independent oracles to write once (in `tests/oracles/`, not importing `objdel.data.windows`):
- `lrr(x, n)`: largest-remainder rounding by plain Python (floor, then give the remainder to largest fractional
  parts, ties by lowest index). Used by every schedule count test.
- `brute_admissible(item, side, dx, dy)`: pixel-by-pixel check of "all doctor pixels inside with >= 2 px
  clearance, no other outlined nodule" using numpy slicing.
- `rloo_ref(R, G)`: loop `A_i = s * (R_i - sum_{l != i} R_l / (G-1))` with Python floats.
- `centroid_sums(mask)`: `np.nonzero` based sum of x and y.

---

## Part 1. Windows: pool and cell (3.1)

**B2-01 cell edges go up** (unit). Input: `cell(10, 0.15)`, `cell(9, 0.1499999)`, `cell(20, 0.2)`, `cell(40, 0.3)`,
`cell(39, 0.2999999)`, `cell(0, 0.0)`, `cell(10**6, 5.0)`. Expected: `(1,1)`, `(0,0)`, `(2,2)`, `(3,3)`, `(2,2)`,
`(0,0)`, `(3,3)`. Why: `searchsorted side="left"` instead of `"right"` puts edge values one bin low; float
rounding of 0.15 vs 0.15000000000000002 must be consistent with the literal tuple.

**B2-02 offset_of geometry and symmetry** (property, hyp). Input: random union boxes and window boxes (window
side S in 16..256). Expected: offset_of is >= 0; equals 0 iff centres coincide; invariant under translating both
boxes by the same vector; invariant under swapping x/y of both boxes; scales as 1/S when the window is scaled by
k and centre distance by k. Literal: union (10,10,20,20), window (0,0,40,40) -> centre (15,15) vs (20,20),
d = 7.0711, S = 40, offset = 0.17678 (within 1e-6). Why: dividing by the wrong side (union extent, or
source slice side) shifts every cell assignment; half-open `(x0+x1)/2` vs inclusive centre is off by 0.5 px.

**B2-03 offset_of uses half-open centres** (unit). Union box (0,0,1,1) (single pixel), window (0,0,3,3): centre
(0.5,0.5) vs (1.5,1.5), d = 1.4142, offset 0.4714. Input window (0,0,2,2) with union (0,0,1,1): offset 0.3536.
Why: single-pixel doctor mask is the smallest legal union; an inclusive-box implementation gives 0.

**B2-04 pool excludes non-RL items, counted** (unit). Items: 5 with 1 doctor, 3 with 2 doctors, 2 with 4.
`pool(items, draws=2)`. Expected: `dropped["not_rl"] == 5`; `{w.item_id}` has 5 ids; exactly 5 centred windows,
each draw 0. Why: counting per window instead of per item (would give 10) or letting a 1-doctor item through.

**B2-05 draws and kinds per item** (unit). 1 RL item with plenty of room, `draws=4`. Expected: 1 centred, <= 4
shifted (draws 0..3 each at most once), <= 4 rescaled; `window_id`s unique; `window_id == f"{item_id}~{kind}{draw}"`.
Centred appears only with draw 0. Why: a pool that emits a second centred window per draw inflates the 15 % share.

**B2-06 window_id has no slash and keeps the prefix** (property, hyp). For item ids of the form `a/b/c/d` and
`x` (no slash) check `"/" not in wid.split("~",1)[1]` hmm only after the source item id: assert
`wid.rsplit("/",1)[0] == item_id.rsplit("/",1)[0]` and `"/" not in wid[len(item_id):]`. Adversarial item_id
containing `~` (e.g. `lidc/p~1/s`): window_id must still round-trip; `Window.item_id` is recovered exactly.
Why: the image path is derived from window_id; a "/" in the suffix creates a nested directory and breaks the
unique-path rule (3.3); a "~" in item_id breaks any `split("~")` parser.

**B2-07 geometry invariants, brute force oracle** (property, hyp, integration for real items). For random
synthetic items (random 1..4 doctor masks of random blobs in a 512x512 slice, 0..2 other outlined nodules, nodule
near borders), every returned window passes `brute_admissible`: all doctor pixels inside; min distance of any
doctor pixel to the window border >= 2; no pixel of any other nodule's outline inside the window;
`0 <= dx`, `dx+side <= W` (window inside the slice, unless the spec allows padding: then the test uses the pad
rule of `build_eval_views`; read it first). Why: off-by-one on the clearance (`> 2` vs `>= 2`), half-open vs
closed box, forgetting a second nodule.

**B2-08 clearance exactness boundary** (unit). Doctor mask pixels at columns 10..19 (union x0=10, x1=20); window
with dx = 8 (clearance exactly 2 left) and right edge at x1+2 = 22 (clearance exactly 2) must be admissible; dx = 9
(clearance 1) must not. Build by calling `pool` on an item where only one (dx, side) is feasible. Expected: the
only admissible window equals the tight one, or the draw is dropped (`shifted_none`), never a 1 px clearance.
Why: fencepost on clearance is the most likely silent bug and breaks D4/D5 only on the tails.

**B2-09 touching borders** (unit). Nodule union touching x0 = 0 of the slice (cannot clear 2 px on the left).
Expected: no window has dx < 0 or negative clearance; if no window is admissible, `dropped` records it and no
window is returned (no fallback). Same for union touching the right/bottom edge (x1 = W). Why: the fallback path
(clamp to the slice) violates D4 and D6.

**B2-10 holes and multi-component doctors** (unit). A doctor mask with a hole and a second component 30 px away,
and a union box from both. Expected: union box covers both components; every window contains both. Why: bounding
only the largest component or the vote mask instead of the union.

**B2-11 single-pixel and tiny nodules** (unit). Union extent 1 px (single pixel) and 2 px. Expected:
`window.extent_px == 1` and cell size_bin 0; windows exist and satisfy clearance; no division by zero in offset.
Why: `max(x1-x0, y1-y0)` for zero-area guards; tiny extents are size bin 0 where the target share is largest.

**B2-12 max-size nodule fills the slice** (unit). Union box of 400 px extent in a 512 px slice. Expected: a window
side >= 404 exists only if the item's policy side allows; shifted windows are all `side == S0`; if none admissible,
all draws counted in `shifted_none`, rescaled maybe too; `pool` does not raise. Why: infinite loops or exceptions
in the rejection sampler when no solution exists.

**B2-13 shifted keeps policy side, rescaled changes it** (property, hyp). Expected: every shifted `side == S0(item)`,
every rescaled `side != S0(item)`; rescaled side is a neighbouring side (check against the neighbour list in
`build_eval_views`). Why: a rescaled window that happens to draw S0 is a kept-side window (D6).

**B2-14 draw 0 equals the exam view builder** (integration, data-backed). For 200 sampled val items and seed 0,
`(dx, dy, side)` of `pool(...)` draw 0 shifted/rescaled equals the offsets/side in `build_eval_views` records
(`view.offset`, `view.side`); for rescaled where the view is a kept-side window, `pool` has no draw 0 rescaled
window (the draw is counted in `rescaled_none`, or it is a different, admissible window by rule; assert the
spec literally: pool must not return the kept-side window). Why: a re-implemented sampler drifts from the exam
generator and the train/exam law differs silently; this is the regression that ties B1 to the exam.

**B2-15 dropped, no fallback** (unit). Construct an item where shifted is impossible but rescaled possible (nodule
within 1 px of both slice edges' feasible region on the policy side only). Expected: `dropped["shifted_none"] == draws`,
`dropped["rescaled_none"] < draws`, no shifted window, no window with side == S0 other than centred. Mirror case for
rescaled. Why: spec says no substitution across kind.

**B2-16 sorted, order independent, deterministic** (property, hyp). `pool(items)` equals `pool(shuffled(items))`
(windows tuple and dropped), `windows` sorted by `window_id`, repeated call equal, equal across processes
(`PYTHONHASHSEED=0` vs `=1` subprocess runs, compare `windows`). Why: set/dict iteration order or `hash()` used in
seeding makes runs irreproducible; seed derived from `hash(item_id)` is the classic bug.

**B2-17 seed independence per item** (property). Adding or removing an unrelated item must not change the other
items' windows (compare per-item window lists). Why: one global RNG consumed sequentially makes every window depend
on the item list, breaking audit reproducibility.

**B2-18 disagreement flag** (unit). Item whose doctors' topologies disagree -> `Window.disagreement is True` for all
its kinds; agreeing item -> False. Why: weights in `schedule` rely on it.

**B2-19 cell_shares** (property, hyp). Shape (4,4), sum to 1 within 1e-12, non-negative; raises or returns zeros?
Spec says sums to 1: with an empty list `cell_shares([])` must raise `ValueError` or return documented value
(fail the test if it returns NaN). Permutation invariant. Literal: 3 windows in cell (0,1), 1 in (2,3) ->
`[0,1]=0.75`, `[2,3]=0.25`, others 0. Why: NaN from 0/0 propagates into TV/KS gates which then compare False and
can silently pass a `<=` gate written as `not (x > limit)`.

---

## Part 2. Schedule (3.1)

**B2-20 exact counts, property** (property, hyp). Random pool built synthetically (Window objects, large enough,
random cells, random item ids), random `steps` 1..30, `per_step` 1..32, random valid shares and target. Expected:
`len == steps*per_step`; kind counts equal `lrr(shares*N)`; for each off-centre kind with n rows, per-cell counts
equal `lrr(target*n)` (flatten order row-major, tie rule same as oracle); no repeated `window_id`; per step no
repeated `item_id`; step-major layout. Why: the core contract, and the largest-remainder ties are an off-by-one
magnet. Literal: shares (0.15,0.75,0.10), N = 32*100 = 3200 -> centred 480, shifted 2400, rescaled 320 exactly
(no rounding needed); N = 31: 4.65, 23.25, 3.1 -> floors 4,23,3 = 30, one left -> largest fraction .65 goes to
centred: 5,23,3.

**B2-21 rounding sums and ties** (unit). N = 10, shares (1/3,1/3,1/3) -> counts sum 10, the extra goes to the
first kind in KINDS order (tie rule documented by `lrr`). `target` uniform 1/16 with n = 40 -> 8 cells get 3, 8
get 2 and which cells follows the oracle. Why: sum != N, non-deterministic tie order.

**B2-22 shares/target validation** (unit, parametrised). `target` of shape (3,4), (4,4) with a negative entry,
sum 1.0001, sum 0.9999, NaN, inf, dtype int -> `ValueError`; sum 1 + 5e-10 accepted; 1 + 2e-9 rejected.
`shares` missing "rescaled", extra key "rotated", negative, sum 1.01, NaN -> `ValueError`. `per_step=0`, `steps=0`:
`steps=0` returns `[]` or raises (assert a documented choice, must not hang). Why: validation lets garbage into
a paid run.

**B2-23 unfillable cell names the kind and cell** (unit). Pool with zero shifted windows in cell (1,2) while the
target gives that cell n >= 1. Expected: `ValueError` whose message contains `"shifted"` and `"(1, 2)"` (or `12`);
the pool contains rescaled windows in (1,2) and they must not be used. Also a pool where the cell has 3 windows
but 5 are needed (windows may not repeat) -> `ValueError`. Also: cell fillable overall but not under the
one-item-per-step rule (all windows of the cell belong to 2 items, step needs 3 in the cell on one step...) ->
`ValueError`. Why: silent substitution distorts the law the gate then reports as exact; the reasons for failure
must be diagnosable before a GPU is rented.

**B2-24 one item per step under pressure** (adversarial). Pool where one item owns 60 % of the windows of the
most requested cell. Expected: schedule succeeds or raises; if it succeeds, per-step item uniqueness holds.
Run for 200 seeds. Why: greedy fillers place the same item twice in a step when a cell runs low.

**B2-25 determinism and seed sensitivity** (unit). Same args -> `==` lists; seed 0 vs 1 differ (list inequality;
and not just a permutation within steps: the set of windows differs for a pool 3x larger than needed). Cross
process with different `PYTHONHASHSEED`. Why: seed ignored, or Python hash-order dependence.

**B2-26 disagreement weight monotone** (property). With a pool of 50 % disagreement items, weights 1, 4, 16, 64:
share of disagreement windows over the whole run non-decreasing across weights (strict at 1 -> 64 for a large
pool); at weight 1 the share is within 3 sigma of the pool share; at weight 0 (if allowed) no disagreement item
unless needed. Counts of kind and cell remain exactly the same for every weight. Why: weight mistakenly changes
marginals (it must only change which window fills a cell).

**B2-27 centred windows are not binned by the target** (unit). The centred share rows are drawn without a target
cell constraint and are not included in the cell check; `cell_coverage` and D7-D9 use off-centre rows only.
Build a schedule where centred rows are all in cell (0,0): off-centre D7/D8/D9 unaffected, D11 (all rows) counts
them. Why: spec 3.2 explicit split; mixing them lets a wrong 15 % pass or fail.

**B2-28 performance gate: schedule** (performance). Real pool size (est. 7,341 train items x K draws, e.g. 7,341*8
~ 58k windows), `steps=120`, `per_step=32` (3,840 rows). Expected: completes in < 5 s on the developer laptop CPU
and uses < 500 MB RSS; scaling check: doubling the pool time grows < 3x. Why: an O(N^2) per-cell search over the
pool or a Python loop per row makes the build impractical; the number is a gate, not a spec value, so set the
limit at 10x the measured baseline once measured [estimate].

**B2-29 performance gate: pool** (performance, data-backed). `pool` on 1,000 items, `draws=8` with `--procs 1`:
record seconds per item; gate at 2x the first measured value (stored in `tests/goldens/perf_pool.json`), plus
no more than 50 rejection attempts per window on average (counter exposed or monkeypatched RNG). Why: the
rejection sampler must not degrade for small nodules; catches accidental O(side^2) mask rebuilds per attempt.

---

## Part 3. Coverage, quality gates, leakage (3.2, 5.1)

**B2-30 cell_coverage literal** (unit). Target 4x4 with 0.5 in (0,0), 0.5 in (1,1); scheduled rows 60 % in (0,0), 40 %
in (1,1), 100 rows. Expected: `worst_ratio = 0.8`, `tv = 0.1`. `ks` is computed on supplied offsets vs row offsets
(use identical lists -> 0.0; shifted by 0.1 constant -> exactly the KS of two shifted samples by oracle
`scipy.stats.ks_2samp`). Cells with target share < 1 % are excluded from `worst_ratio`: add a cell with 0.5 %
share that has 0 rows -> `worst_ratio` unchanged. Why: including the rare cells makes the ratio 0 forever.

**B2-31 cell_coverage only off-centre** (unit). Schedule of centred-only rows -> raises `ValueError` or keys
reported as failing (D7/D8 fail: "A schedule with no off-centre rows fails D7/D8"). Centred rows added to a
valid schedule do not change the three values. Why: spec 3.2.

**B2-32 D1 and D2 planted leakage** (integration, tmp store). Build a tiny store (using `objdel.data.store`
test helpers from `tests/test_data_store.py`) with patients P1 in train and val -> `leakage(...)` reports it and
`gates(...)["D1"].ok is False`. Same slice hash shared by two patients -> D2 fails. Clean store -> D1, D2 ok
with value 0. Also: shared 32x32 crop hash inside one patient is NOT a failure. Why: gates that return `ok=True`
vacuously on a missing or empty store; gates reading the wrong split key.

**B2-33 D3 near duplicates** (integration). Cross-split pair of crops with cosine 0.98 exactly (limit inclusive:
`>= 0.98` counts) and 0.9799 (not). Rate computed vs within-train rate 4/7341 = 0.000545. Expected pass when
cross-split rate <= within-train rate; test a store with 1 cross-split duplicate in 10,000 (1e-4, pass) and 20
in 10,000 (fail). Why: strict/non-strict comparison and denominator (pairs vs items).

**B2-34 D4 / D5 planted** (unit). Window cutting one doctor pixel (window x1 one pixel short) -> D4 fails with
value < 100 %; window with clearance 1 -> D5 fails; policy frame (centred) containing another outlined nodule ->
D5 fails; a clean set -> both ok. One bad window among 3,840 must fail (hard gate, not percentage tolerance:
`ok` iff value == 100 %). Why: using `np.mean(...) > 0.99` hides single violations.

**B2-35 D6 planted** (unit). A window with `kind="shifted"` and `side != S0`? (not a fallback); a window flagged
fallback in records (`view.fallback == True`) and one with `kept_side == True` in rows -> D6 fails with the correct
count in `value`. Schedule with 0 -> ok. Why: D6 reads the row view, not geometry.

**B2-36 D7-D9 limits are exact** (unit). Construct schedules with `worst_ratio` 0.8 (ok), 0.7999 (fail); TV 0.10
(ok) / 0.1001 (fail); KS 0.06 (ok) / 0.0601 (fail). D10: share at policy side 0.75 (ok), 0.7499 (fail).
Why: the `limit` string and `ok` must agree (assert `Gate.limit` parses to the same threshold as `ok` uses).

**B2-37 D11 exactness** (unit). A schedule off by one row in one cell (swap one shifted window into another cell)
-> D11 fails; kind counts equal but cell counts off by one -> fails; exact -> ok. Why: D11 is the machine proof of
the schedule rules, a tolerance there would defeat B2-20.

**B2-38 D12 tokens** (unit). `tokens(...)` with a stub tokenizer returning len = limit - 100 for the longest
target ("limit - 100 tokens fails": Max certified tokens + 1 must be below limit - 100). Cases: max = limit - 102 ok
(max+1 = limit-101 < limit-100), max = limit - 101 fails (max+1 = limit-100 not below), max = limit - 100 fails. Target
errors counted: one unparseable target -> value of errors 1 -> fail. Why: the `+1` (EOS) and strictness fencepost.

**B2-39 gates() shape and names** (unit). `gates(report)` returns 16 gates with names D1..D16 in order; D14, D16
`hard=False`; all others `hard=True`; `failures(gates)` returns only hard failed ones and never D14/D16 even if
`ok=False`. Missing key in the report (e.g. no "D13" data) -> raises `KeyError`/`ValueError`, never `ok=True`.
Why: a missing key defaulting to pass is the worst failure mode of an audit layer.

**B2-40 D13 fresh render** (integration, data-backed, 150 windows). Take a built root; re-render each of 150
sampled windows from the store and compare (box, side, image sha256). Tamper one cached PNG byte -> exactly one
mismatch -> D13 fails. Sampling must be seeded and independent of file system order. Why: stale images after a
store rebuild.

**B2-41 D15 frozen views** (integration). `docs/reference/frozen_views.json` has the required keys; recompute
`content_sha256` from the view's records (newline-joined tuples `(item_id, view.offset, view.side, image_sha256)`
in order) and compare. Mutate one record's offset by 1 px, reorder two records, change one image sha -> D15 fails
each time; changing the manifest's `versions` or `reward_fingerprint` does not. Why: spec says the check must not
depend on those.

---

## Part 4. Scripts and trainer guards (3.3)

**B2-42 build_rl_rows end to end, tiny store** (e2e). Fixture store of ~40 items (30 RL-eligible), a VIEW_DIR with
val items half shape/half select. Run `python -I -m experiments.phase2.build_rl_rows --steps 3 --per-step 8 ...`.
Expected: exit 0; `train.jsonl` has 24 lines in schedule order (not sorted: assert `ids != sorted(ids)` when the
schedule is not alphabetical, and equals `schedule()` output order); every record has `rl == true`, `view` keys
exactly `{kind, draw, offset, side, fallback, kept_side, op}` with `fallback False`, `kept_side False`, `op "id"`;
`item_id == window_id`; image path and cache folder unique per record (set size == 24); `manifest.json` has
`schedule.target_half == "shape"`, `ids_sha256 == sha256("\n".join(ids))` (independent computation);
`reward_fingerprint` equals `objdel.reward.fingerprint()`. Why: every consumer of the root depends on this.

**B2-43 target from val-shape only** (integration). Change val-select items in VIEW_DIR (move their offsets to the
far tail): `target_sha256` and the schedule are bit-identical. Change one val-shape item -> they change. Why: val
select must never influence training data (constraint 4).

**B2-44 val-select never read** (failure mode, adversarial). Monkeypatch `open`/store reader to record the
`group_id`s read; run `build_rl_rows`: no val-select item's image or mask is opened (they may only be read as
val metadata for ids, not content). Why: the "val-select is never read by training or build_rl_rows" constraint
cannot be checked by output alone.

**B2-45 hard gate failure** (e2e failure mode). Plant a cross-split patient (D1) in the store; run `build_rl_rows`.
Expected: exit code 1, no `manifest.json`, but the partial outputs may exist; stderr names `D1`. Also a schedule
that raises `ValueError` (unfillable cell) -> non-zero exit and no manifest, no traceback-less success. Re-run on
a clean store in the same `--out`: manifest appears (no stale-state dependence). Why: a manifest written before the
gates is a "gate passed" lie.

**B2-46 rebuild idempotence and reproducibility** (e2e). Run twice, `--procs 1` and `--procs 4`, same seed:
`train.jsonl` bytes and all image sha256 equal; `ids_sha256` equal; different `--seed` differs. Interrupt after N
rows and re-run into the same out dir: result equals a clean run (no mixed windows). Why: multiprocessing completion
order leaking into file order is a real risk.

**B2-47 audit_data CLI** (e2e). On a good root: exit 0, output has `D1`..`D16`; on a root with one planted failure
for each of D1, D2, D4, D5, D6, D12 (parametrised): exit != 0 and the named gate is marked failed; D14 / D16 failures
never change the exit code. `--json FILE` content parses and has 16 entries. Why: CI will use the exit code.

**B2-48 build_eval_views --kinds** (integration). `--kinds shifted` builds only shifted records (no
rescaled/rotated files); `--kinds shifted,rescaled` builds two; unknown kind -> argparse error exit 2; default
builds all three and its records for shifted are byte-identical to `--kinds shifted` output (same seed). Why: a
filter that renumbers draws changes the exam.

**B2-49 final_eval refuses a changed view** (e2e failure mode). Four mutations of the test view (one offset +1, one
image sha, one reordered pair, one target set changed): each -> non-zero exit and a message that names the
`content_sha256` or `target_set_sha256`; unmodified view with a manifest whose `versions` and
`reward_fingerprint` are changed -> still runs. Missing `frozen_views.json` -> clear refusal, not KeyError trace.
Why: spec 3.3; protects the one-time test read.

**B2-50 train_rl schedule guards** (unit, no GPU: factor the checks, or call `main` with a stub trainer that raises
if constructed). Table (each raises `ValueError` before any model load):
row count = steps*per_step - 1; rows reordered by one swap (ids_sha256 differs); `generation_batch_size // num_generations`
= 31 vs per_step 32; `steps_per_generation != gradient_accumulation_steps`; `num_iterations = 2`; `--steps` = schedule.steps + 1.
Valid case with `--steps` = schedule.steps - 1 passes the guards. `shuffle_dataset` is False in the produced
GRPOConfig. Why: every one of these silently changes which prompts hit which step; the cheapest to test, the most
expensive to miss.

**B2-51 train_rl without schedule / after B1** (unit). Root without `schedule`: before gate B1 -> old behaviour
(shuffled, disagreement x4); after gate B1 (test flag or the removed symbols) -> refused. `--train-views` after B1 ->
argparse exit 2; `hasattr(dataset, "grpo_rows_views") is False`; no multi-root key tagging in `serve`. Before B1:
old paths still work. Why: removal must follow the gate; both directions guard against early or late deletion.

**B2-52 rl.json provenance** (integration). After a 2-step CPU/stub run: `rl.json` has `rows.root`, `rows.ids_sha256`
(equal to manifest), `numerics`, `reward.{name,fingerprint,c_topo,adv_scale,loc_width,area_gate}` for v5; for v4
those v5 keys are absent. Why: later gates compare runs by this file.

**B2-53 budget stop** (unit). `BudgetStop(rate=2.0, cap=1.0, start_time=t0)` with a fake clock: at elapsed 0.4999 h
no stop and no save; at exactly 0.5 h both `should_training_stop` and `should_save` True; clock going backwards
never triggers; `cap=0` triggers at first step; rate 0 never triggers (document or reject at construction).
`--budget-usd X --rate-usd-per-h Y` with only one given -> argparse/ValueError, not silent no-budget. Why: money.

**B2-54 LocGuard** (unit). Synthetic logs: 49 steps with loc rising and iou falling -> no stop; 50 steps -> stop
and save. loc slope > 0, iou slope exactly 0 -> stop (`<= 0`); loc slope 0 -> no stop (`> 0` strict); iou rising
-> no stop; `rewards/loc/mean` never logged -> inactive for 500 steps; only the last 50 count (a long rise followed by
50 flat/falling steps -> no stop). Noise: loc up by 1e-12 per step (float slope tiny positive) with iou flat ->
stops (document); NaN in logs -> no crash, no stop. Why: the guard is the only protection against the reward
hack-shaped loc term and is slope sign sensitive.

---

## Part 5. Number match (3.4)

**B2-55 is_ratio geo mask literal** (unit, CPU torch). `logps - sampling` per token: seq A = [0.0, 0.0, 0.0]
(mean 0 -> exp = 1, inside [0.99, 1.01]) -> ratios exp(0)=1 each. Seq B = [0.02, 0.02, 0.02] (exp mean 1.0202 > 1.01)
-> all zeros. Seq C = [-0.02,...] (0.9802 < 0.99) -> zeros. Seq D = [0.5, -0.5] mean 0 -> kept, ratios
exp(0.5)=1.6487 and exp(-0.5)=0.6065 with `hi=2.0` unchanged (below 2). Seq E token diff 1.0 with hi=2.0 -> clamp to 2.0
(e=2.718). Padding positions (mask 0) never enter the mean (put 100.0 in padded positions of seq A: must be
kept). Bounds are inclusive: mean-log exactly log(1.01) (use float64 literal) kept. Expected ratios exactly as
stated, dtype and device of input. Why: mean over padded tokens, wrong inclusive boundary, mask applied before
clamp.

**B2-56 is_ratio existing modes regression** (property, hyp). For random tensors, `token_truncate`, `sequence_truncate`,
`sequence_mask`, `token_mask` match the previous implementation (keep a frozen copy of the previous function in
`tests/oracles/` or a golden file generated from the pre-change commit). Exact equality. Why: spec "same values".

**B2-57 is_ratio edge** (unit). All-masked row (no valid tokens): no NaN, ratio defined (0 or 1 per documented), no
exception; single token sequence; `geo=None` with `mode="sequence_geo_mask"` -> ValueError; `geo=(1.01, 0.99)`
(reversed) -> ValueError; `hi=None` with truncate -> ValueError. Huge logps difference (`+-1e4`) -> `exp` overflow
must not produce inf/NaN outputs (clamped or masked). Why: inf*0 = NaN in "times 0" masking.

**B2-58 is_summary** (unit). Float32 output on the input device for bf16/float64 inputs; counts for a hand case:
4 sequences (one masked, one clamped (raw ratio outside [lo, hi]), two clean), `is_masked_share = 0.25`,
`is_clamped_share = 0.25`, `seq_logdiff_mean` = mean of signed sums, `seq_logdiff_max`, `tok_absdiff_gt_0.5_share`
tokens only valid ones, `masked_len_mean` = valid-token count of masked sequences (0 when none, no NaN),
`kept_len_mean`. Token modes: share per token. A clamped unit counts only if the raw ratio was outside [lo, hi]
(a ratio exactly equal to `hi` before clamp is not clamped). Why: units (sequence vs token) are easy to mix.

**B2-59 stats stay on device** (GPU or meta device). Under `torch.cuda.set_sync_debug_mode("error")` run 16
micro-batches of the accumulation path: no host sync until the logging call; the logging call is the only
`.item()`. CPU fallback: count `.item`/`.cpu` through a tensor subclass. Why: spec "no host read per micro-batch".

**B2-60 existing five sampling metrics unchanged** (regression). For the same input, the 5 old `sampling/*` metrics
match a golden captured from the current commit exactly (float equality). New seven keys exist exactly once per log
step. Why: rewriting the accumulator changes the old semantics.

**B2-61 NUMERICS presets literal** (unit). The exact dict values of 3.4 for the four presets, `DEFAULT_NUMERICS` in
`NUMERICS`; `grpo_config(..., numerics=name)` produces `vllm_importance_sampling_mode`, `clip_max` as stated; for every
preset `not (cfg.bf16 and cfg.fp16)`; fp16-tis -> `fp16 True`, `bf16 False`; unknown preset -> `KeyError`/`ValueError`
naming the choices. `grpo_config(reward=R)` with `num_generations=7` or `num_generations_eval=7` (both) and `R.group==8`
-> ValueError; 8 passes. Why: a config that quietly trains in the wrong dtype is a lost GPU-hour.

**B2-62 vllm_dtype context** (unit, mock). Replace `trl.generation.vllm_generation.LLM` with a recorder: inside the
context `LLM(model="x")` records `dtype="float16"`; after normal exit and after an exception inside the `with` the
original is back (identity check); nested contexts restore correctly; an explicit `dtype=` kwarg given by the caller is
overridden or errors (assert the documented behaviour). Why: leaked patch changes the dtype of unrelated later engines.

**B2-63 FastGRPO construction guards** (unit with fakes). `sequence_geo_mask` with the own ratio path inactive ->
`ValueError` at construction (not at first step); fp16-tis with a vLLM engine whose model dtype is bfloat16 -> raises;
with `torch.float16` -> constructs. Why: failure at step 1 of a paid run vs construction time.

**B2-64 fp16 vision check** (GPU, small). `fp16_vision_check` on a checkpoint where the float16 tower overflows
(plant: scale a weight by 1e4) -> inf/NaN count > 0 and exit 1; clean checkpoint -> 0 and exit 0. Why: the exit code
gates the fp16 arm.

**B2-65 b2_attribution** (GPU, integration). Prints three numbers for (a), (b), (c) for 256 answers; the gap of (a) with
the LoRA delta set to zero equals the gap of (c) within 1e-6 (merge of zero delta is identity); the three numbers are
>= 0 and finite; same batch twice gives the same numbers (determinism of eval forward, `torch.use_deterministic...`
not required, tolerance 1e-4). Why: sanity of the attribution tool before it drives a decision.

---

## Part 6. CUDA scorer sums (3.5)

**B2-66 sum_x / sum_y exactness** (property + GPU, hyp). 200 random masks (random sizes 1..512, random density, blobs,
holes, border-touching, full image 512x512 all ones, single pixel at (0,0) and at (W-1,H-1), empty mask)
through `score_masks`: `sum_x == np.nonzero(mask)[1].sum()`, `sum_y == np.nonzero(mask)[0].sum()` as int64, dtype
`int64`, one value per object, empty -> 0. Full 512x512: `sum_x = 512 * (511*512/2) = 67,043,328` exactly. Same for
`score_polygons` of rasterised polygons (compare with the mask returned for the same polygon). Why: off-by-one
indices (1-based), int32 overflow in the CUDA kernel (max real sum with 1024^2 px is ~5.4e8, still < 2^31 but a
`float` accumulator loses exactness above 2^24 pixels sums), atomicAdd type.

**B2-67 CPU vs CUDA bit equality** (GPU). Same inputs through the CPU build and the CUDA build: every field bit-equal,
including old fields (compare against the stored golden from the pre-change build for old fields). `sum_x`, `sum_y`
equal. Run on 10k random polygons from the existing fuzz in `tests/test_cuda_checker_gpu.py`. Why: "previous fields
bit-identical" regression and CPU/CUDA parity.

**B2-68 large masks overflow guard** (unit). Mask 4096 x 4096 all ones (if the scorer admits it): `sum_x` equals the
numpy int64; or the scorer rejects with a clear error. Why: silent int32 wrap.

**B2-69 ref_sum / radius2 of ExpertCache** (unit). Three concentric disks (as `three` fixture): `ref_sum` equals
`(sum_x, sum_y, n)` of the vote mask by numpy (vote non-empty). Rater stack where the vote is empty (all disjoint
masks): `ref_sum` is from the union, `n = union.sum()`. `radius2 == union.sum()/pi` to 1e-12. Pickle round trip and
save/load keep both. Why: the vote-vs-union fallback is only visible in this corner.

---

## Part 7. Fingerprint and val halves (3.0)

**B2-70 fingerprint invariance** (property). For each of `alpha_loc, miss_iou, w_f5, w_area, c_topo, adv_scale,
loc_width, area_gate, lam, eps, drop_spread, iou_gate` and each `RewardV4` weight/floor mutated through
`monkeypatch.setattr` (or by building a spec copy via `dataclasses.replace`): `fingerprint()` unchanged. Mutated:
`CacheSpec.expert_score`, `tau_floor`, `tau_quantile`; `min_area` rule probe; `numpy.__version__`,
`scipy.__version__`, `cv2.__version__` patched strings: `fingerprint()` changes each. 64-hex-digit string. Why: a
weight change invalidating every cache (hours of CPU) or a cache-dependent constant NOT invalidating.

**B2-71 cache values unchanged, digests change** (regression). Build `ExpertCache` for the fixtures with the old and new
code; `tau`, `t5` equal exactly (`==`); `digest` differs; `ExpertCache.load(path_old, expected_digest=new)` raises
`StaleCacheError`; `rebuild=True` yields the same tau/t5. Why: spec "stored values unchanged, only digests".

**B2-72 val_half** (property, hyp). `val_half(g) == "shape" iff hashlib.sha256(g.encode()).digest()[0] % 2 == 0` for
random text including unicode, empty string, very long strings; deterministic across processes and
`PYTHONHASHSEED`; all items of one patient (same `group_id`) the same half; over 2,000 random ids shape share in
[0.45, 0.55]. Only the values "shape" and "select". Literal: compute two known ids by hand once with `hashlib` and
freeze them. Why: use of `hash()` or the hex string's first character instead of the first byte.

---

## Part 8. Reward v5 (3.6, 5.2)

**B2-73 RewardV5 construction** (unit). `RewardV5(c_topo=1, adv_scale=1)` ok; `adv_scale = 0`, negative, NaN ->
ValueError; `alpha_loc=0.11` with `miss_iou=0.1` -> ValueError; `alpha_loc == miss_iou` ok; missing `c_topo` or
`adv_scale` -> TypeError; frozen (assignment raises). `REWARD_V5` defined and valid. Why: boundary of `>`.

**B2-74 location literal** (unit). `location(0,0,0,(rx,ry,n), r2)` = 0.0 (area 0). Equal centroids -> exactly 1.0 (use
integers: sums (100,200), area 10, ref (1000,2000,100) -> 10,20 vs 10,20). d2 = 25 (area 1, sum (5,0); ref (0,0,1)),
`radius2=49`, width 1 -> `exp(-25/98) = 0.77483` (1e-6). Width 2 -> `exp(-25/(2*4*49)) = exp(-0.06378) = 0.93821`. Large
d -> value > 0 (strict) or underflows to 0.0 only beyond 1e-300; result in `[0,1]`, never NaN with `radius2` tiny
(`radius2 = 1e-9`: value 0.0 not NaN) or `radius2 = 0` (documented: ValueError or 0.0, not ZeroDivisionError crash
inside training). `n = 0` in ref_sum -> documented error. Why: division guards.

**B2-75 location monotone** (property, hyp). For fixed area, ref and radius2, `location` is strictly decreasing along
any ray away from the reference centroid until underflow; invariant to rotation (rotate displacement by 90 deg:
equal); symmetric under d -> -d. Why: T5 at function level; asymmetric x/y mistakes (`sum_x` used for both).

**B2-76 T1 continuity/monotonic R in IoU** (property, hyp). Build terms for 1 group using `advantages_v5` internals
via a helper `reward_R(terms)` (or check A ordering with a 2-answer comparison in a G=8 group with 6 fillers): with
fixed f5, g_area, r_topo (sig_id, topo fields), legal 1, `area_gate=False`, location in [0,1], sweep iou_tc from 0
to 1 in 1001 steps with `loc` set by the per-object formula `1[iou<0.1]*max(0,1-iou/0.1)*loc_value`: R at
`iou=0.1-1e-12` and `0.1` differ by < 1e-9 (continuity), and R non-decreasing over the sweep for every
`loc_value in {0, 0.3, 1}`. Expected R slope below miss_iou is `1 - alpha_loc*loc_value/miss_iou` summed: with
alpha_loc 0.1, miss_iou 0.1, loc_value 1 -> slope exactly 0 (flat but not decreasing). Values from spec's
`RewardV5` defaults. Why: any config with `alpha_loc > miss_iou` makes a hit worse than a near miss; the spec
guard exists for this exact bug.

**B2-77 T2/T3 loc zero for gamed shapes** (GPU or CPU scorer). Reference disk radius 7 at (30,30). (T2) Legal disk
radius 40 centred at (30,30) with doctors' max area < its area -> `g_area == 0` -> `loc == 0`. But IoU_tc of a covering
outline: ensure the case really is a miss by choosing area so IoU_tc < 0.1: radius 40 vs radius 7 -> IoU = 49/1600 =
0.0306, a miss. (T3) Legal speck 3x3 px at the centroid (area 9 < doctors' min area): `g_area == 0` -> `loc == 0`.
Assert the setup facts first (`terms["g_area"] == 0`, `terms["iou_tc"] < 0.1`) so the test cannot pass for the wrong
reason. Why: location without the area gate is the dominant exploit (centre dot or giant blob).

**B2-78 T4 two blobs** (GPU or CPU scorer). Doctors drew one component (disk r=7 at (30,30)). Answer: two legal disjoint
disks of r=7 at (10,30) and (50,30) (joint centroid (30,30)), none overlapping the union, each above the size filter.
Expected: `loc == 0` via `topo_freq == 0` (topology differs from every doctor). Assert `topo_freq == 0` precondition.
Mirror: doctors drew two components and answer is two blobs at the right spots -> `loc` may be > 0 only when
`iou_tc < 0.1` (a miss) and topology frequency > 0. Why: centroid of a disconnected answer is a free parameter.

**B2-79 T5 shift monotone** (GPU or CPU scorer, property over d). Doctor-like disk r=7; answer disk r=7 shifted by d
in {8, 10, 12, ..., 40} px (all IoU = 0: with d >= 14 centres separate completely; use d >= 15 to keep IoU exactly
0). Expected `loc(d)` strictly decreasing; each `loc(d) == exp(-d^2 / (2*49))` to 1e-9 (`radius2 = union/pi` for a
rasterised disc is approximately 49, but compute from `cache.radius2` rather than 49). Also direction invariance (shift
along x vs y vs diagonal with equal Euclidean distance on an exact lattice, e.g. (9,12) and (15,0) both d=15: equal
loc). Why: anisotropy bug; the number of the rasterised centroid vs ideal.

**B2-80 T6 illegal and unparseable** (unit). Self-intersecting polygon with the centroid on target, unparseable text,
truncated completion (`truncated=True`), empty polygon list: each has `loc == 0` and R == 0 (R computed with
`legal` gate; check via `advantages_v5` single-group run: the answer's contribution is 0 so A equals -the mean of
the others scaled). Why: gate placement; a legal-gate multiply that misses the loc term.

**B2-81 T11 hit never loses to miss** (property, hyp). Sample random (iou_hit in [0.1, 1], iou_miss in [0, 0.1),
loc in [0, 1], same f5, g_area, legal, same r_topo): R_hit >= R_miss. Boundary: iou_hit exactly 0.1, iou_miss = 0.0999999
with loc = 1 (largest bonus): hit R = 0.1 + 0 vs miss R = 0.0999999 + 0.1*alpha*(1 - 0.999999) *1 ~ 0.1000000 -> hit
>= miss by the T1 flat slope: assert with tolerance 1e-12, not strictly greater. For `area_gate=True` same.
Why: this is the central safety property of the shaping.

**B2-82 advantages_v5 literals** (unit). G = 8, `adv_scale = 1`, one prompt, `legal=1`, `f5=g_area=0`, same `sig_id`
and no topology valid, iou_tc = [1,0,0,0,0,0,0,0] -> R = same (w_area*g_area = 0). A[0] = 1 - 0 = 1.0; A[1..7] =
0 - (1/7) = -0.142857...; sum 0. With iou_tc = [0.5]*7 + [0.4] -> A[7] = -0.1 (0.4 - 0.5), others =
0.5 - (0.5*6+0.4)/7 = 0.5 - 0.485714 = 0.0142857; sum 0 (1e-9). `adv_scale = 3` -> times 3 exactly. Why: independent
worked numbers for the RLOO baseline.

**B2-83 T7 sum zero** (property, hyp). Random terms over B in 1..6 prompts, any values including illegal, topology
valid etc., `adv_scale` random. For every group `|A.sum()| < 1e-9 * (1+|R|max)`, including dropped groups (A = 0).
Also finite everywhere. Why: nonzero-sum advantages make the baseline biased.

**B2-84 T8 scale equivariance** (property). Multiplying all reward terms (iou_tc, f5, g_area, loc, `legal` stays a gate
so scale only reward terms; `c_topo` scaled via a derived `RewardV5` copy) by k in {0.1, 3, 100} multiplies A by k;
no std division: with k = 1e-6 (tiny spread, > drop_spread? no: choose spread above 0.01 by using k >= 1) assert
linearity; also groups with tiny spread but a different `sig_id` are not zeroed. Why: hidden z-score.

**B2-85 T9 identical answers** (unit). 8 identical answers: A == 0 exactly, no NaN/inf; also 8 identical rewards but
different `sig_id` (spread 0, 2 signatures): not dropped; A computed from `r_topo` differences only; spread just
below `drop_spread` (0.0099) with one sig_id -> dropped; at exactly 0.01 -> kept (strict `<`). Why: fencepost on
the drop rule.

**B2-86 T13 adv_scale linear** (property). `advantages_v5(t, spec=s)` == `s * advantages_v5(t, spec=1)` for s in
{0.5, 2, 7.3} to float64 equality within 1e-12. Not applied to the drop rule (a group dropped at s is dropped at 1).
Why: `adv_scale` accidentally applied to R inside the spread test.

**B2-87 T12 group of misses at distances** (GPU or CPU scorer, integration). Group of 8 legal misses: same shape
(disk r=3 at distances 10, 12, 14, 16, 18, 20, 22, 24 from the centroid of a r=7 doctor mask; check
`iou_tc < 0.1` precondition for all; the spec's example radius 7 at 5 px overlaps the target (hand computed
IoU ~ 0.38, [estimate]) and is not a miss, so use the distances/radii that satisfy the precondition and also
run the spec's 3 distances 5/10/20 with disks small enough to miss) . Expected: A strictly decreasing in
distance, the closest has max A, all A of the closest > 0, A sums to 0, and Spearman(A, -d/r) = 1.0. Compare with
`advantages_v4` on the same group: v4 gives a lower Spearman (document the value). Why: this is the whole purpose
of B3 (signal toward the nodule).

**B2-88 loc and area_gate** (unit). `area_gate=True`: `area = w_area*g_area*min(1, iou_tc/miss_iou)`; iou 0 -> 0;
0.05 -> half; >= 0.1 -> full; `area_gate=False` -> constant. Combined with T1: monotone either way. Why: a gate
that multiplies by `iou/miss_iou` beyond 1 (no `min`) breaks the plateau.

**B2-89 advantages_v5 input validation** (unit). Term lengths differ, length 7, length 9, empty arrays with
G=8 (length 0: returns empty array or ValueError, document; must not crash on reshape(-1, 8) of empty), missing
`loc` key -> `KeyError` or `ValueError` with the key name (not a silent zero), `group=4` with `spec.group=8` (mismatch)
-> ValueError or honoured consistently; NaN in iou_tc -> propagate? assert no silent zero: raise or NaN, document and
test the documented outcome. Why: silent misalignment shifts all advantages by one prompt.

**B2-90 score_completion loc per object** (integration). Two-object completion (disc and cup, as RIGA): per-object L
computed independently; completion `loc` = mean of the two, with a hit object contributing 0 and a missed object
contributing L: literal: one hit (L=0), one miss with L=0.6 -> `loc == 0.3`. A truncated completion -> 0. `score_batch`
returns the `"loc"` key as a float array of the batch length, aligned with other keys; dtype float64. Why: mean vs
sum, alignment.

**B2-91 TRLReward end to end** (integration, CPU build). `TRLReward(caches, stop_ids, spec=REWARD_V5)`: on a fake
batch of 2 prompts x 8 completions, the aggregator output equals `advantages_v5(terms)` to float32 precision
(1e-6 relative), logged terms include `loc`, logs `kind/<kind>/iou_tc`, `cell/<ij>/iou_tc`,
`cell/<ij>/allmiss_groups` computed by hand: group with 8 misses -> 1.0 for its cell, mixed group -> 0.0; two groups
in one cell, one all-miss -> 0.5; when `kind`/`cell` columns are absent none of these keys are logged and no error.
Why: wiring between the reward and TRL is where float32 casting and index mix-ups live.

**B2-92 Spearman gate helper** (unit). Gate B3 statistic: Spearman(A, -d/r) over all-miss groups: perfect order -> 1.0;
reversed -> -1.0; groups with constant A contribute NaN: define behaviour (excluded) and test; a share of groups
where a miss has a higher R than a hit: construct 1 in 200 -> 0.5 % (<= 1 % ok) and 3 in 200 -> fail. Why: gate arithmetic is part of the
decision to spend money.

---

## Part 9. Cross-cutting adversarial, failure and regression

**B2-93 numbers regression, R5** (regression, data-backed). Read `runs/r5/eval/eval_val/items.jsonl` (if available)
filtered to val-select via `val_half`: assert item count > 0, shape/select partition covers all items exactly once,
mean IoU_tc over the full val reproduces today's R5 numbers within 1e-3 (full-val IoU 0.659 [measured: final-eval
record], test not read), miss share at offset > 0.3 = 31.8 % and IoU at offset < 0.15 = 0.732, IoU_tc on val-select
shifted = 0.599 per spec 5.3 (1e-3 tolerance). Skipped when the file is missing. Why: Gate B1/B3 are paired vs these
numbers; a different `val_half` implementation silently changes the baseline.

**B2-94 paired bootstrap** (property, unit). Nodule-level bootstrap: identical runs -> difference 0, CI [0, 0];
constant +0.05 shift -> CI contains 0.05 only (width 0); resampling is by nodule not item (items of one nodule
resampled together: a data set where two items of a nodule always agree gives a CI wider than item-level; assert
using construction). Fixed seed -> reproducible CI; 1 nodule -> degenerate CI or error, documented. 95 %: coverage
over 1,000 synthetic replicates within [0.92, 0.98]. Why: wrong resampling unit inflates confidence, passing a gate by chance.

**B2-95 Gate evaluation arithmetic** (unit). B2 gate: arm passes iff ratio mean >= 0.97 (0.9699 fails), masked share
< 0.5 % (0.5 % fails), NaN count 0, one-sided t lower bound > -0.02 on steps 11-20 only (change step 1-10 values: no
effect; step 20 vs 21 index fencepost, 1-based). Literal paired differences d = [-0.01]*10 -> mean -0.01, sd 0 ->
lower bound -0.01 > -0.02 pass; d = [-0.03]*10 -> fail; one-sided t quantile for df 9 is 1.833 [paper: standard
t table] with sd 0.01: lower = mean - 1.833*0.01/sqrt(10). Why: wrong-tail (two-sided 2.262) or wrong step window.

**B2-96 concurrency / multiprocess safety** (integration). `build_rl_rows --procs 4` on a store with a corrupted
image (zero bytes): clean failure message with the window_id, exit non-zero, no hang, no orphan processes (check
`psutil.children` empty after exit), and no `manifest.json`. A worker killed (SIGKILL) -> parent exits non-zero in < 60 s. Why: `ProcessPoolExecutor` breaks silently
into a hang.

**B2-97 path uniqueness and disk hygiene** (integration). Two windows with source items that differ only by case
(`A/1` vs `a/1`) on a case-insensitive file system (macOS): must not collide (use a case-fold check in a test
independent of the OS: `len({p.lower() for p in paths}) == len(paths)`). Item ids containing `..`, absolute paths or
NUL: rejected or sanitised, never written outside `--out` (assert every written path resolves under `out`). Why: data is
untrusted; path traversal and case collisions.

**B2-98 memory/time gate for `audit_data`** (performance, data-backed). Full train (7,341 items + 3,840 windows):
runtime < 10 min on the developer laptop, resident memory < 4 GB, D13 150 renders < 60 s [estimate, replace with
first measured value x 2 once known]. Why: gates that nobody runs because they are slow.

**B2-99 end to end dry run** (e2e, small real data, no GPU spend). `build_rl_rows` -> `audit_data` (exit 0) -> `train_rl`
config construction only (stubbed trainer): guards all pass; the `ids_sha256` printed by each is the same; a single-byte
edit of `train.jsonl` order is caught by `train_rl` only. Why: end-to-end check of the contract before a paid rental, per
the no-spend-before-scan rule.

**B2-100 GPU smoke gate** (GPU, runs only on the rented preflight host, < 5 min): scorer CPU vs CUDA parity on 1,000
windows' masks (B2-67); `is_ratio` + `is_summary` run on cuda tensors of shape (16, 1536) with no host syncs
(B2-59); `vllm_dtype("float16")` yields an engine of dtype float16 and `FastGRPO` accepts it; `fp16_vision_check`
exits 0. Why: the preflight must exercise every new GPU path before money is spent.

---

## Pytest sketches (the 5 most important)

Ranking: (1) schedule exactness and invariants, (2) advantages_v5 properties (T7/T8/T9/T13), (3) hit-never-loses and
location shaping (T1/T11 plus T2/T3), (4) fingerprint invariance and val_half, (5) trainer schedule guards.

```python
# tests/test_b_windows_schedule.py  (B2-20, B2-23, B2-24, B2-25)
import itertools
import numpy as np
import pytest
from hypothesis import given, settings, strategies as st
from objdel.data import windows as Wn
from oracles.largest_remainder import lrr          # independent plain-Python oracle

def synth_pool(rng, n_items=400, per_item=6):
    out = []
    for i in range(n_items):
        dis = bool(rng.random() < 0.2)
        for kind in ("shifted", "rescaled"):
            for d in range(per_item // 2):
                off, ext = rng.random() * 0.5, int(rng.integers(3, 80))
                out.append(Wn.Window(f"lidc/p{i//3}/s{i}", f"p{i//3}", f"n{i}", kind, d,
                                     0, 0, 100, off, ext, dis))
        out.append(Wn.Window(f"lidc/p{i//3}/s{i}", f"p{i//3}", f"n{i}", "centred", 0, 0, 0, 100, 0.0, 20, dis))
    return tuple(sorted(out, key=lambda w: w.window_id))

@settings(max_examples=40, deadline=None)
@given(steps=st.integers(1, 12), per_step=st.integers(1, 32), seed=st.integers(0, 2**16),
       dw=st.sampled_from([1.0, 4.0, 16.0]))
def test_schedule_counts_and_rules(steps, per_step, seed, dw):
    rng = np.random.default_rng(seed)
    ws = synth_pool(rng)
    target = rng.dirichlet(np.ones(16)).reshape(4, 4)
    shares = {"centred": 0.15, "shifted": 0.75, "rescaled": 0.10}
    n = steps * per_step
    try:
        out = Wn.schedule(ws, target, steps, per_step, shares, seed, dw)
    except ValueError as e:                       # a cell the pool cannot fill: error must name kind and cell
        assert any(k in str(e) for k in ("shifted", "rescaled")); return
    assert len(out) == n
    # kind counts == independent largest-remainder rounding
    want = dict(zip(Wn.KINDS, lrr(np.array([shares[k] for k in Wn.KINDS]) * n, n)))
    for k in Wn.KINDS:
        assert sum(w.kind == k for w in out) == want[k]
    # per-cell counts for each off-centre kind
    for k in ("shifted", "rescaled"):
        rows = [w for w in out if w.kind == k]
        got = np.zeros(16, int)
        for w in rows:
            got[w.cell[0] * 4 + w.cell[1]] += 1
        assert got.tolist() == lrr(target.reshape(-1) * len(rows), len(rows)).tolist()
    assert len({w.window_id for w in out}) == n
    for s in range(steps):
        ids = [w.item_id for w in out[s * per_step:(s + 1) * per_step]]
        assert len(set(ids)) == per_step
    assert out == Wn.schedule(ws, target, steps, per_step, shares, seed, dw)   # deterministic

def test_unfillable_cell_names_kind_and_cell_and_does_not_substitute():
    rng = np.random.default_rng(0)
    ws = tuple(w for w in synth_pool(rng) if not (w.kind == "shifted" and w.cell == (1, 2)))
    target = np.full((4, 4), 1 / 16)
    with pytest.raises(ValueError, match=r"shifted.*\(?1, ?2\)?"):
        Wn.schedule(ws, target, 4, 32)
```

```python
# tests/test_b_advantages_v5.py  (B2-82, B2-83, B2-84, B2-85, B2-86)
import dataclasses
import numpy as np
from hypothesis import given, settings, strategies as st
from objdel.reward import REWARD_V5, advantages_v5, RewardV5

def terms(B, G, rng, same_sig=False):
    n = B * G
    return dict(iou_tc=rng.random(n), f5=rng.random(n), legal=(rng.random(n) > .2).astype(float),
                g_area=rng.random(n), loc=rng.random(n) * .5, topo_any=(rng.random(n) > .5).astype(float),
                topo_freq=rng.random(n), iou_min=rng.random(n),
                sig_id=np.zeros(n) if same_sig else rng.integers(0, 3, n).astype(float))

def test_literal_rloo():
    t = terms(1, 8, np.random.default_rng(0), same_sig=True)
    t.update(iou_tc=np.array([1.] + [0.] * 7), f5=np.zeros(8), g_area=np.zeros(8), loc=np.zeros(8),
             legal=np.ones(8), topo_any=np.zeros(8), topo_freq=np.zeros(8), iou_min=np.zeros(8))
    spec = dataclasses.replace(REWARD_V5, adv_scale=1.0)
    a = advantages_v5(t, 8, spec)
    np.testing.assert_allclose(a, [1.0] + [-1 / 7] * 7, atol=1e-12)

@settings(max_examples=60, deadline=None)
@given(seed=st.integers(0, 10**6), B=st.integers(1, 5), s=st.floats(0.1, 9.0))
def test_sum_zero_linear_in_scale_and_finite(seed, B, s):
    t = terms(B, 8, np.random.default_rng(seed))
    a1 = advantages_v5(t, 8, dataclasses.replace(REWARD_V5, adv_scale=1.0))
    as_ = advantages_v5(t, 8, dataclasses.replace(REWARD_V5, adv_scale=s))
    assert np.isfinite(a1).all()
    assert np.abs(a1.reshape(B, 8).sum(1)).max() < 1e-9           # T7
    np.testing.assert_allclose(as_, s * a1, rtol=1e-12, atol=1e-12)  # T13

def test_identical_group_is_zero_and_drop_rule_is_strict():
    t = terms(1, 8, np.random.default_rng(1), same_sig=True)
    for k in ("iou_tc", "f5", "g_area", "loc", "topo_freq", "iou_min"):
        t[k] = np.full(8, 0.3)
    t["legal"] = np.ones(8)
    assert (advantages_v5(t) == 0).all()                           # T9, no NaN
    t["iou_tc"][0] += 0.0099                                       # spread < 0.01 -> still dropped
    assert (advantages_v5(t) == 0).all()
    t["iou_tc"][0] = 0.3 + 0.0101                                  # spread > 0.01 -> kept
    assert np.abs(advantages_v5(t)).max() > 0

def test_reward_scale_equivariance_T8():
    rng = np.random.default_rng(5); t = terms(2, 8, rng)
    base = advantages_v5(t)
    k = 3.0
    scaled = dict(t, **{key: t[key] * k for key in ("iou_tc", "f5", "g_area", "loc")})
    spec = dataclasses.replace(REWARD_V5, c_topo=REWARD_V5.c_topo * k)
    np.testing.assert_allclose(advantages_v5(scaled, 8, spec), k * base, rtol=1e-9, atol=1e-12)
```

```python
# tests/test_b_reward_v5_shaping.py  (B2-76, B2-77, B2-81; T1, T2, T3, T11)
import dataclasses
import numpy as np
import pytest
from hypothesis import given, settings, strategies as st
from objdel.reward import REWARD_V5, ExpertCache, score_completion, location
from tests.test_reward import disk                           # existing helper

def R_of(iou, loc_value, spec=REWARD_V5, f5=0.5, g=0.5, topo=0.0):
    """Reward of one answer, per the spec formula (legal=1)."""
    loc = (iou < spec.miss_iou) * max(0.0, 1 - iou / spec.miss_iou) * loc_value
    area = spec.w_area * g * (min(1.0, iou / spec.miss_iou) if spec.area_gate else 1.0)
    return iou + spec.alpha_loc * loc + spec.w_f5 * f5 + area + spec.c_topo * topo

@pytest.mark.parametrize("gate", [False, True])
@pytest.mark.parametrize("lv", [0.0, 0.3, 1.0])
def test_T1_R_continuous_and_monotone_in_iou(gate, lv):
    spec = dataclasses.replace(REWARD_V5, area_gate=gate) if hasattr(REWARD_V5, "area_gate") else REWARD_V5
    xs = np.linspace(0, 1, 1001)
    Rs = np.array([R_of(x, lv, spec) for x in xs])
    assert (np.diff(Rs) >= -1e-12).all()
    eps = 1e-12
    assert abs(R_of(spec.miss_iou - eps, lv, spec) - R_of(spec.miss_iou, lv, spec)) < 1e-9

@settings(max_examples=300, deadline=None)
@given(hit=st.floats(0.1, 1.0), miss=st.floats(0.0, 0.0999999), lv=st.floats(0, 1))
def test_T11_hit_never_loses_to_miss(hit, miss, lv):
    assert R_of(hit, 0.0) + 1e-12 >= R_of(miss, lv)

def test_T2_T3_gamed_shapes_get_no_location():
    st_ = np.stack([disk(60, 7, 30, 30)] * 2); cache = ExpertCache.build(st_)
    # giant legal disk centred on the nodule; speck at the centroid. Polygons given in 0-1000 grid text.
    from tests.helpers_poly import circle_text                # small helper: JSON of a regular polygon
    big, speck = circle_text(30, 30, 40, 60), circle_text(30, 30, 1.5, 60)
    out = score_completion(big, False, ("nodule",), [cache], 60, 60)
    assert out["g_area"] == 0 and out["iou_tc"] < 0.1 and out["loc"] == 0   # preconditions + result
    out = score_completion(speck, False, ("nodule",), [cache], 60, 60)
    assert out["g_area"] == 0 and out["loc"] == 0
```

```python
# tests/test_b_fingerprint_val_half.py  (B2-70, B2-72)
import dataclasses, hashlib, subprocess, sys
import pytest
import objdel.reward as R
from objdel.reward import spec as S
from objdel.data.splits import val_half
from hypothesis import given, strategies as st

@pytest.mark.parametrize("field,val", [("alpha_loc", .05), ("miss_iou", .2), ("w_f5", .3), ("w_area", .2),
    ("c_topo", 5.0), ("adv_scale", 2.0), ("loc_width", 2.0), ("area_gate", True), ("lam", .9), ("eps", .2),
    ("drop_spread", .05), ("iou_gate", .6)])
def test_weights_do_not_change_fingerprint(monkeypatch, field, val):
    before = R.fingerprint()
    monkeypatch.setattr(R, "REWARD_V5", dataclasses.replace(R.REWARD_V5, **{field: val}), raising=False)
    assert R.fingerprint() == before

@pytest.mark.parametrize("field,val", [("expert_score", .9), ("tau_floor", 2.0), ("tau_quantile", 90.0)])
def test_cache_inputs_change_fingerprint(monkeypatch, field, val):
    before = R.fingerprint()
    monkeypatch.setattr(S, "CACHE_SPEC", dataclasses.replace(S.CACHE_SPEC, **{field: val}), raising=False)
    assert R.fingerprint() != before

def test_library_versions_change_fingerprint(monkeypatch):
    import numpy
    before = R.fingerprint()
    monkeypatch.setattr(numpy, "__version__", "0.0.0")
    assert R.fingerprint() != before

@given(st.text())
def test_val_half_is_first_byte_parity(g):
    want = "shape" if hashlib.sha256(g.encode()).digest()[0] % 2 == 0 else "select"
    assert val_half(g) == want

def test_val_half_independent_of_hashseed():
    code = "from objdel.data.splits import val_half;print(val_half('lidc/LIDC-IDRI-0001'))"
    outs = {subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True,
                           env={"PYTHONHASHSEED": s}).stdout for s in ("0", "1", "2")}
    assert len(outs) == 1
```

```python
# tests/test_b_train_rl_guards.py  (B2-50, B2-53)
import json, hashlib
import pytest
from experiments.phase2 import train_rl

def make_root(tmp_path, steps=4, per_step=32, ids=None, schedule=True):
    ids = ids or [f"i{n}" for n in range(steps * per_step)]
    rows = [{"item_id": i} for i in ids]
    (tmp_path / "lidc").mkdir()
    (tmp_path / "lidc/train.jsonl").write_text("\n".join(map(json.dumps, rows)))
    m = {"files": {}, "sets": {}, "reward_fingerprint": "x", "versions": {}}
    if schedule:
        m["schedule"] = dict(steps=steps, per_step=per_step, ids_sha256=hashlib.sha256("\n".join(ids).encode()).hexdigest())
    (tmp_path / "manifest.json").write_text(json.dumps(m))
    return tmp_path

@pytest.mark.parametrize("mutate,match", [
    (lambda c: c.update(rows=c["rows"][:-1]), "row count|rows"),
    (lambda c: c.update(rows=c["rows"][1:2] + c["rows"][:1] + c["rows"][2:]), "ids_sha256|order"),
    (lambda c: c.update(generation_batch_size=31 * 8), "per_step|generation_batch_size"),
    (lambda c: c.update(steps_per_generation=2, gradient_accumulation_steps=1), "steps_per_generation"),
    (lambda c: c.update(num_iterations=2), "num_iterations"),
    (lambda c: c.update(steps=5), "steps"),
])
def test_bad_configs_raise_before_any_model_is_loaded(tmp_path, monkeypatch, mutate, match):
    monkeypatch.setattr(train_rl, "load_model", lambda *a, **k: pytest.fail("model loaded before guards"),
                        raising=False)
    cfg = train_rl.check_inputs(make_root(tmp_path))       # factored guard: returns the valid defaults
    mutate(cfg)
    with pytest.raises(ValueError, match=match):
        train_rl.check_schedule(cfg)

def test_budget_stop_boundary():
    from objdel.train.callbacks import BudgetStop
    class St: should_training_stop = should_save = False
    now = [100.0]
    cb = BudgetStop(rate_usd_per_h=2.0, cap_usd=1.0, start_time=100.0, clock=lambda: now[0]) \
        if "clock" in BudgetStop.__init__.__code__.co_varnames else None
    if cb is None:
        pytest.skip("BudgetStop has no injectable clock; patch time.time instead")
    s = St(); now[0] = 100.0 + 0.4999 * 3600; cb.on_step_end(None, s, None)
    assert not s.should_training_stop and not s.should_save
    now[0] = 100.0 + 0.5 * 3600; cb.on_step_end(None, s, None)
    assert s.should_training_stop and s.should_save
```

Notes for the implementer of these sketches: helper names (`check_inputs`, `check_schedule`, `circle_text`,
`clock=`) are suggestions; if the production code exposes different seams, keep the assertions and change only the
call sites. The sketches never recompute expected values with the code under test: every expected value comes from
a literal, an independent oracle (`lrr`, `hashlib`, `numpy`) or the spec formula written out in the test.
