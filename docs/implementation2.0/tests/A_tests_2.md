# Phase A tests, set 2 (adversarial, property, performance, failure mode)

Source: `impl/phases/A.spec.md` (rev 2) and REPO code/tests (style of `tests/test_p2_signal.py`, `tests/test_final_eval.py`,
`tests/test_profile.py`). Written without the design doc. Spec section numbers are in brackets.

## Skills used
Skill tool calls were not made in this subagent run (no skill content was loaded). The structure follows the
testing-strategy and TDD habits the task named: one failing, specific assertion per behaviour; property tests with
seeded RNG; exact (`==`) comparisons wherever the spec says "bit" or "exactly"; no mocks of our own code, only of
`torch.cuda.nvtx`, subprocess and the clock. Fixtures: `tests/goldens/` (r5 golden), `tmp_path`, tiny synthetic masks.

Levels: U unit, P property, I integration, G GPU (`@pytest.mark.gpu`, skipped on the Mac), E end-to-end.
All Mac tests need no GPU, network, vLLM, `runs/` or `data/`. Every random test uses `np.random.default_rng(<fixed>)`.

---------------------------------------------------------------------------------------------------------------
## 1. Cluster bootstrap `objdel.eval.signal.draws` [3.1]

**T01 golden CI numbers (regression of today).** U/I. Setup: load `tests/goldens/r5_test_items.jsonl.gz`.
Input: ids, `iou_tc`; `final_eval.summarize` items built from them (only `iou_tc` needed; use a stub-item path or a
helper that fills other keys with 0). Expect: `ci95_item == [0.6334010210209712, 0.6561989221793691]`,
`ci95_nodule == [0.6279017988891661, 0.6620867640136667]`, `miss_rate == 0.10818181818181818`, all with `==`.
Also assert the golden has 2200 items and 401 distinct clusters (guards a changed fixture). Catches: any change of
RNG order, `seed`, quantile method (`np.percentile` vs `nanquantile` differ in the last digit), or cluster rule.

**T02 legacy bit-identity.** P. Setup: a verbatim copy of the pre-change `draws` (frozen in the test file as
`_legacy_draws`). Inputs: 50 random id lists (n 1..120; mix of `lidc/s/nXXX/zYYY` with 1..6 slices per nodule, and
non-lidc ids; shuffled order, so clusters interleave), `boot` in {1, 7, 100}, `seed` in {0, 1, 12345}.
Expect `len` equal and every array `np.array_equal` and equal dtype. Catches a refactor that changes first-appearance
order (e.g. `sorted(set())`) or draws the group index with a different call (`choice` instead of `integers`).

**T03 one item per cluster equals the plain item bootstrap.** P. Input: `ids = [f"x{i}" for i in range(n)]`, n in
{1, 2, 17, 401}; boot 50, seed 3. Expect `np.stack(draws) == default_rng(3).integers(0, n, (50, n))` exactly. Also
`lidc/` ids with unique nodules (one slice each) give the same. Catches a draw order mismatch between the item and
cluster paths; this is the spec's own invariant and the reason `ci95_nodule == ci95_item` in gate 1.

**T04 length and membership invariants.** P. Input: random clustered ids (cluster sizes 1..9, 2..60 clusters).
For each draw: every index in [0, n); the array is the concatenation of whole clusters (split it by cluster sizes in
the drawn order: each block equals `np.flatnonzero(cluster == c)` exactly, in ascending order inside the block); the
length equals the sum of drawn cluster sizes. Catches slices of one nodule being resampled separately (leakage of
independence), partial clusters, or unsorted blocks.

**T05 custom `key` (patient).** U. Input: ids `lidc/s1/n000/z1, lidc/s1/n000/z2, lidc/s1/n001/z5, lidc/s2/n000/z9`;
`groups = {id: "P1" for first three, "P2" for the 4th}`; `key=groups.__getitem__`. Expect C = 2 clusters; every draw
is a concatenation of blocks `[0,1,2]` and `[3]` only. Edge: `key` that maps all ids to one value gives arrays equal
to `arange(n)` for every draw (C = 1). `key` that is injective equals the default on non-lidc ids.
Catches the key being ignored (nodule fallback) or applied to the cluster string instead of the id.

**T06 patient interval is wider or equal on average.** P. Synthetic: 40 patients x 3 nodules x 4 slices; per-patient
random offset added to IoU so the data correlate inside patients. Expect width(`ci95_patient`) >
width(`ci95_nodule`) > width(`ci95_item`) for the seeded fixture (assert each strict). Catches the three CIs sharing
one code path or `groups` not used. Choose the seed once, then freeze it; if it flips with a seed it is the test
that is wrong, so also assert the weaker "not equal" form for 5 seeds.

**T07 edge sizes.** U. `draws([], boot=3)` : either returns `[]` arrays of length 0 or raises `ValueError`; fix one
behaviour in the spec and assert it (the current code raises `ValueError` from `integers(0, 0)`; assert that the
failure is a clear `ValueError`, never a silent NaN). `draws(["a"], boot=4)` gives four arrays `[0]`.
`boot=0` gives `[]`. Catches a divide by zero or an empty-array NaN CI in `summarize` of 1 item (assert
`summarize([one_item])` returns CI `[v, v]`, no NaN, no exception).

**T08 determinism and independence of `seed`.** P. Same `(ids, boot, seed)` twice gives equal arrays; `seed` 0 vs 1
differ; calling `draws` does not touch the global `np.random` state (`np.random.get_state()` before and after equal).
Catches use of the legacy global RNG.

**T09 performance gate.** Perf. 2,200 items, 401 clusters, boot 2000: `draws` under 2.0 s on the Mac (current is a
Python loop; record the time of the legacy copy in the test and assert new <= 1.5 x legacy). `summarize` of the golden
with all 12 METRICS plus 5 selection keys, nodule + patient + item: under 30 s and memory (tracemalloc peak) under
600 MB (2000 index arrays of ~2200 int64 = 35 MB per arm; assert draws are not cached per metric: count calls to
`signal.draws` with a spy, expect exactly 2 per `summarize` call, one nodule one patient). Catches an O(metrics x boot)
re-draw.

---------------------------------------------------------------------------------------------------------------
## 2. Boundary IoU and `metrics()` additions [3.4]

**T10 band shapes by hand.** U. 5x5 square in 9x9 frame, d = 1: `band(m,1).sum() == 16`; d = 2: 24 (outer ring of the
5x5 plus the second ring: 25 - 1 = 24); d = 3: 25 (nothing survives, whole mask is band). 1x1 mask any d: band is the
pixel (1). 1xN line: band is the whole line (erode of a 1-high line is empty). Empty mask: band empty, dtype bool,
same shape. Full-ones frame (touching all four borders), d = 1: band is the 1-pixel frame ring, `== 4*(n-1)`;
this is the property that the 1-pixel zero pad treats the frame edge as background. Catches an implementation that
lets `cv2.erode` use its default border (replicate or infinity), which would give an empty band for a full frame.

**T11 d iterations versus 1-pixel pad (adversarial).** P. For d in {2, 9, 16} and masks that touch every edge, compare
`band` with a reference that pads by `d` pixels, erodes d times, and crops. They must be equal. Reason: the spec pads
only one pixel; after the first erosion the pad pixel is still 0 and 3x3 erosion shrinks from it, so the results are
equal; the test proves the spec claim for large d. If it fails, the spec is wrong, and that is a real finding.
Input: 200 random blobs (union of random rectangles and discs) 64x64, some clipped at the frame.

**T12 reference equality on random masks.** P. Reference (`full-frame` in the test): `cv2.copyMakeBorder(m,1,1,1,1,
BORDER_CONSTANT,0)`, `cv2.erode(kernel ones(3,3), iterations=d)`, crop, `m & ~e`. Inputs: 500 random masks with
density in {0.02, 0.5, 0.98}, shapes in {(1,1), (1,7), (7,1), (16,16), (64,48)}, d in {1,2,9,16}, all touching the
edge by construction. Expect `boundary_iou(a,b,d) ==` the reference formula exactly (`==`, not approx), computed
with integer counts then one division. Catches float accumulation (`np.mean`) differences, wrong dtype (uint8 overflow
in `sum`), and a cropped-box implementation (see T13).

**T13 crop-box shortcut must equal full frame.** P. The implementation may work in a joint bounding box for speed
(the repo does for `_pairwise_iou`). Input: a mask at the exact frame corner (rows 0..4, cols 0..4), another
shifted by 1; box crop without padding would drop the zero pad and under-count the band. Expect equality with the
full-frame reference, d = 2 and 9. Also a 512x512 frame with a 3-pixel nodule at (0,0). Catches the one real bug class
the spec names ("masks that touch the frame edge").

**T14 algebraic properties.** P. Symmetry `bi(a,b,d) == bi(b,a,d)`; `bi(a,a,d) == 1.0` including empty and including
masks thinner than 2d (band = whole mask); range [0,1]; disjoint bands give 0.0; `bi(empty, nonempty) == 0.0`
(exactly one empty band); `bi(empty, empty) == 1.0`. Adversarial: a non-empty mask whose band is non-empty vs a mask
that is non-empty but whose band is empty cannot happen (the band of a non-empty mask is never empty, since the erosion
is a strict subset) — assert `band(m,d).any()` for every non-empty random m. That last property would fail if the
pad were missing and a full-frame mask were passed.

**T15 `metrics()` keys and shapes.** I. Use the existing synthetic item builders from `tests/test_eval.py` (disc
rasters rendered to polygon text), K = 4, k = 3, W = H = 128, names = ("nodule",). Expect `boundary_iou` and
`boundary_iou_vote` lists of length K, `sample_iou` a K x K nested list, diagonal 1.0 for legal answers, symmetric,
every entry in [0,1]. Existing keys equal a frozen copy of the pre-change output (store the pre-change dict for this
fixture in `tests/goldens/` as JSON; compare with `==` after removing the three new keys). Catches a regression of
today's outputs (spec: "Existing keys and values are unchanged").

**T16 illegal and failed answers.** I. K = 4: answers = [good, "not json", good-but-self-intersecting (illegal),
"" (empty)]. Expect `boundary_iou[1] == boundary_iou[2] == boundary_iou[3] == 0.0` and same for `_vote`;
`sample_iou[1][1]` == 1.0 only if its mask counts as empty (two empty masks IoU 1.0); `sample_iou[1][2] == 1.0`
(both empty), `sample_iou[0][1] == 0.0`. Catches an illegal answer being scored by its raster (the checker returns
masks for illegal answers too; `metrics` already zeroes via `drawn = masks & legal`), and the 0/0 division for two empty
masks giving NaN (assert `not np.isnan(np.array(sample_iou)).any()`).

**T17 max over raters.** I. Item with three raters where the answer equals rater 2 only: `boundary_iou == 1.0` for
that answer although it differs from raters 0 and 1. Permute the rater order (`np.random.permutation`): the value does
not change. Catches argmax over the wrong axis or "first rater" scoring.

**T18 vote mask strictness.** U/P. k = 2 raters whose masks overlap in half: a pixel in only one rater has
2*1 > 2 false, so is not in the vote; k = 4 with exactly 2 of 4 true: not in the vote (strict). Property: for random
rater stacks, `ExpertCache.vote == (raters.sum(0)*2 > k)` and `select.vote_mask(raters, [True]*k)` equals it.
Catches `>=` instead of `>`.

**T19 `sample_iou` with several objects (disc, cup).** U. In `per_item`, `sample_iou` is averaged over objects first:
disc IoU matrix A and cup IoU matrix B give `(A+B)/2` elementwise, exactly. Catches using only the first object.

---------------------------------------------------------------------------------------------------------------
## 3. Selection `objdel.eval.select` [3.5]

**T20 medoid hand cases.** U. (a) Tie: 3 answers all pairwise IoU 0.5 -> index 0. Tie at the top between 1 and 2
with 0 lower -> 1. (b) No legal answer -> 0 even though the matrix is not trivial. (c) One legal answer at index 5
-> 5. (d) Three legal empty answers (IoU 1.0 with each other) and two legal non-empty distinct answers -> index of the
first empty one (a majority of legal empty answers wins). (e) K = 2 both legal -> 0 (tie, equal means). (f) Illegal
answer with IoU 1.0 against everything must be ignored: it must never be returned, and must not raise another's score.
Expected values are hand-computed in the test, written as literals.

**T21 eligible-only mean (adversarial).** U. 4 answers, answer 3 illegal. Matrix such that the mean over ALL j != i
picks index 0 but the mean over eligible j != i picks index 2. Expect 2. Catches including illegal columns, or dividing by K - 1 instead of
(#eligible - 1), or including the diagonal (diagonal 1.0 would shift every score equally only when the eligible set is
the same size; it differs here, so the bug shows).

**T22 permutation property.** P. For random symmetric `S` (K = 8, entries in [0,1], diagonal 1) and random legal:
permuting answers by pi gives `medoid(S[pi][:,pi], legal[pi])` whose position maps back to a score-maximal answer:
`score[pi[m']] == score[m]` (the index may differ on ties; compare the score, not the index). With random float
matrices ties are improbable, so also assert the equal index mapping with `pi[m'] == m`. Catches
position bias beyond the "lowest index wins ties" rule.

**T23 no doctor input (leak guard).** U. `inspect.signature(medoid).parameters` is exactly `(sample_iou, legal)`;
`per_item` computes `sel_medoid` with the doctor-mask caches replaced by `None`/garbage in the part that feeds
`medoid` (monkeypatch `ExpertCache.raters` to raise on access inside the medoid branch only; the rest of `per_item`
may read them). Hard check: build two `result` dicts with identical `sample_iou` and `samples.legal` but different `iou_tc` and rater
data; `medoid` index is equal; `sel_medoid` differs only by the `iou_tc` at that index. Catches a label leak into
selection (spec constraint "Selection never reads doctor masks").

**T24 permuting rater masks leaves the medoid unchanged (gate 4).** I. Run `metrics()` twice, once with raters in
order (0,1,2) and once (2,0,1); `sample_iou` is equal exactly; `medoid` equal; `sel_vote` equal; `sel_oracle` equal.

**T25 `vote_mask`.** U/P. No legal mask -> all-False array of shape (H, W), dtype bool. One legal mask -> equals it.
2 legal masks, strict: pixel true only if both true. 3 legal of 8, with illegal masks all-True: the illegal ones are
excluded (result differs from the all-8 vote). Property: `vote_mask(m, legal)` equals `m[legal].sum(0)*2 > legal.sum()`
for 300 random cases, K in 1..8, H, W in 1..9. Catches counting illegal masks in the threshold.

**T26 per_item selection keys.** U. Hand result with K = 4: iou_tc = [0.2, 0.9, 0.0, 0.6], legal = [1,1,0,1].
`sel_mean_of_k == 0.425`, `sel_legal_mean == (0.2+0.9+0.6)/3` (compare `pytest.approx`, rel 1e-12), `sel_oracle == 0.9`,
and with the matrix favouring answer 3: `sel_medoid == 0.6`, `sel_medoid_legal == 1.0`. No legal answers: `sel_legal_mean == 0.0`,
`sel_medoid_legal == 0.0` (medoid index 0 is illegal), no NaN. `sel_oracle >= sel_medoid >= 0` and
`sel_oracle >= sel_mean_of_k` always (property over 200 random items). Missing `sel_vote` data: the key is absent,
and `summarize` of a mix of items with and without `sel_vote` either omits `vote_minus_mean` as null or raises a
clear error; assert one, fix it in the spec (spec says null).

---------------------------------------------------------------------------------------------------------------
## 4. Per-item and summary schema [3.2]

**T27 METRICS contract.** U. `set(FE.METRICS) == {legal, iou_tc, f5, precision_vote, recall_vote, dice_vote, nsd,
boundary_iou, boundary_iou_vote, soft_dice, hm_iou, ged}`, exactly 12; only `ged` has `False`; each source is one of
`samples`, `objects`, `item`. Catches an added or dropped metric silently changing the table.

**T28 per_item averages.** U. Extend the existing `result()` helper: K = 3 lists and item scalars; `miss` is 1.0 at
`iou_tc == 0.0999999` and 0.0 at exactly `0.1` (threshold is strict `< 0.1`). `iou_tc` exactly 0.1 in float:
build from list `[0.1, 0.1, 0.1]` whose mean is `0.10000000000000002` or `0.1`; assert the code compares the
mean, and pin the value with a literal. Catches `<=` and a per-sample miss rule.

**T29 summary schema.** U. `summarize(items, groups)`: `schema == 2`, `boot == draws`, `ci_level == 0.95`, `nodules`
and `patients` counts, no key `ci95` anywhere (recursive walk of the dict), every metric has exactly
{mean, ci95_nodule, ci95_patient, ci95_item, higher_is_better} (iou_tc plus `median`). Without `groups`:
`patients is None` and every `ci95_patient is None`. JSON round trip (`json.dumps` without a default) works, so there
is no NaN or numpy scalar left (`allow_nan=False`).

**T30 CI contains the mean (property).** P. For random item sets, `lo <= mean <= hi` for the item bootstrap; for
nodule and patient bootstrap too (percentile CI contains the point estimate for n >= 30 in practice; test with
n >= 60 and fixed seeds). Constant metric: `[v, v]` exactly. `ged` has `higher_is_better False` in the output.

**T31 `paired`.** U. Items with `sel_medoid = sel_mean_of_k + 0.05` exactly: `medoid_minus_mean == [d, lo, hi]`
with d == 0.05 (approx 1e-12), lo and hi both approx 0.05. A vote-less set gives `vote_minus_mean is None`. Items
where the difference is symmetric noise around 0: lo < 0 < hi. Catches a sign flip and the CI taken on the two arms
separately instead of the per-item difference (the unpaired interval would be wider for correlated arms: assert
paired width < unpaired width on correlated data).

**T32 `iou_found`.** U. 10 items, 3 misses with iou_tc 0.0, 0.05, 0.09; found items mean computed by hand; CI from
nodule draws restricted to the found subset (ids subset passed to `draws`). All items missed: `iou_found` is null or
a clear value, never NaN-in-JSON (pick one in spec; assert JSON-safe). `miss_rate` with all found is 0.0.

**T33 groups missing an item.** F. `groups` lacking one item id -> `KeyError`/`ValueError` naming the item id, not a
silent skip (a silent skip makes patient CI cover fewer items than the mean). Extra groups keys are ignored.

---------------------------------------------------------------------------------------------------------------
## 5. Ceiling [3.6]

**T34 loo_caches shapes.** U. k = 4 all voting: 4 entries, indices 0..3 increasing, each cache has k == 3, `single`
False; entry r built from `np.delete(stack, r, 0)` (compare `cache.digest` against `ExpertCache.build(that).digest`).
k = 3: caches have k == 2. k = 2: two entries, each `single` True, cache k == 2 with both rows equal to the other
rater (`np.array_equal(cache.raters[0], cache.raters[1])`). Non-voting rater (signature rejected by `voters`): no
entry for it, but it still takes part in the other entries' stacks. Catches indexing by position instead of by
`voters` entry, and the leave-out of non-voters.

**T35 k = 2 equals plain IoU (property).** P. 100 random mask pairs (including empty vs non-empty, identical, thin):
`score_mask(a, loo[0].cache).iou_tc` equals the numpy IoU of `a` and `b` (`(a&b).sum()/(a|b).sum()`), within 1e-12;
exact `==` when both are non-empty with integer areas. Note the exception: both empty (IoU 1.0 by the spec for samples,
undefined for raters) — assert the chosen behaviour is not NaN.

**T36 doctor equal to another scores 1.0.** U. Three identical rater masks (k = 3): `doctor == 1.0`; `legacy_dmax ==
[1.0, 1.0, 1.0]`. Four identical masks with the same signature: doctor IoU_tc equals the legacy max IoU (gate 3).

**T37 model answer equals a doctor.** I. Model answer text = rater r's own polygon (`geometry` rasterization equals
the rater mask, see T41). Then the entry r contributes `score(answer, loo[r].cache)`; assert `diff` for that answer
against `loo[r]` equals the doctor term of r, so the per-entry difference is exactly 0. Run for every r. The
mean `diff` of an answer equal to the best doctor can be non-zero for the others; assert only the entry-r equality.

**T38 illegal and failed answers.** F. `texts = ["garbage", "{}", legal_text]`, `truncated = [False, True, False]`.
Model contribution of the first two is 0.0 exactly; `model == legal_part / 3` (hand). `truncated=None` accepted.
`texts` shorter than `truncated` raises `ValueError`. `texts=None`: `model, model_all_k, diff` are `None` and `doctor`
still present.

**T39 model_all_k versus model.** U. When every answer is the vote-like consensus that scores IoU_tc 1.0 against
the full cache but < 1 against each leave-one-out cache, `model_all_k > model`. Hand-built 3-rater case. Documents the
optimism gap the phase exists to remove. Property: `model_all_k >= model - 1e-9` is NOT always true; do not assert it.

**T40 summarize strata.** U. Rows with k in {2,2,3,4,4,4}: `n(all) == n(k>=3) + n(k=2)` (6 == 4 + 2). An empty
stratum (all k = 4): `k=2` has `n == 0` and `[None]*3` or `[nan]*3` serialised as null; `all == k>=3` numerically
(exact). Same k everywhere: `legacy_loo_dmax_pooled == legacy_loo_dmax_item mean` (exact up to 1e-12). With
different k: pooled differs from the item mean (weighted by voters); assert they differ on a 2-item example by hand:
item A dmax [0.5,0.5], item B dmax [1,1,1,1] -> item mean 0.75, pooled (0.5*2 + 4)/6 = 0.8333333. `diff` CI equals the
CI of `model - doctor` per item, not a difference of two CIs.

**T41 target text rasterizes to the rater mask.** P. 200 random simple polygons rasterized by the repo's own routine,
turned into target text by the repo's target writer, `geometry.parse` + `to_masks` gives a mask equal to the original
on a 512 frame (`np.array_equal`). Include a polygon touching all frame borders, a 1-pixel-wide sliver, and a hole
case (a hole must either round trip or be listed as unsupported in the test; assert the repo behaviour). Catches an
off-by-one in the 0-1000 grid scaling (the pre-existing root of many mismatches).

**T42 legacy ceiling regression on the box data [box].** G/E. val `legacy_loo_dmax_item` within +-0.005 of 0.734 and
test 0.732 (gate 3), and `n(all) == 1084`, `n(k=2) == 281` on val. Skipped without `data/`.

---------------------------------------------------------------------------------------------------------------
## 6. Final eval CLI, rescore, failure modes [3.2]

**T43 write order and partial file.** I. Monkeypatch the sampler to return 3 items x K = 8 fixed texts, and monkeypatch
the scorer to raise `RuntimeError("boom")`. Run `main([...])`. Expect `SystemExit` or exception propagates; and
`answers.jsonl` has 3 lines with the required fields (`item_id, group_id, W, H, samples (8), truncated (8), tokens (8)`),
`final_eval.json["partial"] is True` with `ckpt, view, split, K, seed, view_manifest_sha256, code_sha, checker_build`
and `seconds.load_and_sample`; no `metrics.jsonl` and no `items.jsonl`. The same run without the raise ends with all four
files and `"partial": false` plus `schema == 2`. Catches scoring-before-saving (the loss of a paid GPU sample run).

**T44 atomic write.** F. Kill (raise inside) between writing `answers.jsonl` and the end of `final_eval.json` writing
(monkeypatch `Path.write_text` to raise on the 2nd call): no truncated half-JSON `final_eval.json` is left (the file is
either absent or fully valid). Use write to temp then `os.replace`. Catches a crash that leaves invalid JSON, which
`--rescore` would then fail to parse.

**T45 rescore equals direct (bit equality).** I. Direct run on a 5-item synthetic fixture, then `--rescore DIR`: for
each item, `samples.legal, iou_tc, f5, g_area` lists compare with `==` (no approx). Item lines of `items.jsonl` equal
modulo the `seconds`. `metrics.jsonl` lines `json.dumps(sort_keys)` equal. Rescore loads no model: monkeypatch
`transformers.AutoModel*`/`AutoProcessor` and `torch.load` to raise; the test passes only when never called.

**T46 rescore with a partial `final_eval.json`.** I. DIR has `"partial": true`. The rescore finishes and writes
`partial: false` into `O`. Also a full one. A `DIR/final_eval.json` missing: exit with non-zero and a message (not a
traceback with `KeyError`); spec silent, fix it.

**T47 manifest SHA mismatch exits 2.** F. Change one byte in `V/manifest.json` after the sampling run. `main` returns /
`SystemExit` code 2; `O` does not exist or is empty (`not any(O.iterdir())`). A manifest with trailing newline changed
also fails (hash is over bytes). Same manifest: passes. Catches hashing parsed JSON (reformatting would pass) and
writing partial outputs before checking.

**T48 `code_sha` precedence.** U. (a) `CODE_SHA` file present with `"abc123-dirty-0123456789ab\n"` in a tmp repo root that
also has `.git`: result is the stripped content. (b) no file, a `.git` dir: `git rev-parse HEAD` output (monkeypatch
`subprocess.run`; assert the command list). (c) neither: `None`, serialised as JSON null. (d) `CODE_SHA` present but empty:
define it (treat as absent, then fall back); assert. (e) `git` binary missing (`FileNotFoundError`): `None`, no crash.
Catches a crash on the GPU box where `.git` is absent.

**T49 sync_up writes CODE_SHA text.** U. Pure function that builds the string: clean tree -> `<40 hex>`; dirty ->
`<40 hex>-dirty-<12 hex>` where the 12 hex equals `sha256(diff).hexdigest()[:12]` of the supplied diff text; the same
diff gives the same suffix; a diff differing by one byte gives a different suffix; empty diff means clean (no
`-dirty` suffix). Test with a `git` fake (a tmp git repo made with `git init` is also fine). Untracked-only changes: spec
says "differs from HEAD" through `git diff HEAD`, which ignores untracked files; assert that and flag it in the review as a
known blind spot (an untracked new module will not mark dirty).

**T50 answers.jsonl contents.** I. `tokens` equals the list returned by the stubbed `trainer._generate` ids lengths, K ints;
`truncated` K bools (JSON `true/false`, not 0/1); `samples` are str. One line per item in dataset order; reading the
file twice gives equal output. `W` and `H` are ints.

**T51 K != 8 and k variations.** E. `--k 2` and `--k 8`: metrics shapes follow K. K = 1 raises (metrics requires >= 2).
Item with k = 2 raters works for every metric (hm_iou, ged finite). Item with a non-voting rater works.

---------------------------------------------------------------------------------------------------------------
## 7. Sampling `objdel.train.sample.sample` [3.3]

**T52 tokens are lengths of the completion id sequences.** U. Fake trainer whose `_generate(prompts)` returns
`(prompt_ids, completion_ids_list, ...)` with known lengths [5, 1, 1536, 0, ...]; `sample` returns the same ints.
Ordering is by answer index for each item, `items` chunking (items=16, 40 rows -> 3 chunks) does not reorder rows. An
EOS-padded answer counts its real length. Truncated = (length == max_completion_length and no EOS) as the existing code
defines; assert the existing behaviour is unchanged by comparing with a frozen copy of the old output on the fake.
Catches `len(text)` instead of `len(ids)` (chars vs tokens).

---------------------------------------------------------------------------------------------------------------
## 8. NVTX [3.7]

**T53 import is inert.** U. In a subprocess (`python -I -c`): set `OBJDEL_NVTX=1`, install a `torch.cuda` stub that raises
on any attribute access (`sys.modules` patch before import), import `objdel.train.nvtx`: no exception. `enabled()` then
raises or returns False? Spec: `enabled()` queries `is_available()` at first call; so only the import is checked. Also
with `OBJDEL_NVTX` unset, `"0"`, `"true"`, `""`: `enabled()` is False and `is_available` is never called (counter).

**T54 enabled caches the CUDA query.** U. Monkeypatch `torch.cuda.is_available` to a counter returning True; call
`enabled()` 1000 times: counter == 1. `OBJDEL_NVTX=1` set after the first call does not change the result (cached). Provide
a `reset` through `importlib.reload` in the test.

**T55 range balance on exception.** U. With a fake `nvtx` recording `push/pop` events and enabled True:
`with range("a"): raise ValueError` leaves events `[push a, pop]`. Nested `range("a")`, `range("b")`: `[push a, push b,
pop, pop]`. Disabled: zero calls to the fake (spy raises on call). `push`/`pop` mirrored. Catches pop skipped on error
(leaves the NVTX stack unbalanced and corrupts every later nsys range).

**T56 NvtxSteps push/pop balance (property).** P. Drive callbacks with fake `state.global_step`:
`on_train_begin`, then n in 0..30 x (`on_pre_optimizer_step`, `on_optimizer_step`, `on_step_end`), `on_train_end`.
Assert the event list: pushes == pops, and the stack depth never goes below 0 or above 2 (`step` + `optimizer`), `step`
is always the outermost, and the number of `step` ranges == n + 1... exactly n instances close with a step end plus the
final one opened at the last `on_step_end` is closed at `on_train_end`. n = 0: a train with no steps is balanced. `on_train_end`
without `on_train_begin` (an exception earlier): no pop on an empty stack (assert no unbalanced pop, or a defined no-op).
`on_optimizer_step` without a pre-step (gradient accumulation edge): defined, no crash.

**T57 profiler window.** U. window (4, 6), CUDA faked available and `torch.cuda.profiler` replaced by a recorder; run
steps 1..8 with `global_step` set as the HF Trainer does. `start` is called exactly once, in the `on_step_end` with
`global_step == 3`, after that call's `pop` and before its `push` (compare positions in one shared event list);
`stop` exactly once at `global_step == 6`, after the pop (and nothing pushed after stop in the same callback? spec says
after the pop; the next push still follows: assert stop precedes the next push). No CUDA: no profiler call but ranges still
balanced when `enabled()` is True... (no CUDA means `enabled()` False so no range calls at all; assert the callback does
not crash). `window=None`: no profiler call ever. Window (1, 1): `s - 1 == 0` never equals a post-increment
global_step >= 1; document the behaviour (start never fires) and either reject s < 2 in the constructor with `ValueError`
or assert no start; choose and pin. Resume from checkpoint with `global_step` already > s: no start, no stop. A trainer
that stops before `e`: `on_train_end` must call `stop()` if started and not yet stopped (otherwise nsys capture never
ends; the spec lists no such rule, so this test is expected to FAIL against the spec as written: report it).

**T58 make_trainer wiring.** U. `a.profile_window=None`: no `NvtxSteps` in callbacks; set: exactly one, constructed with
a tuple of ints. `a.trace` set: `TraceStep` present with `wait == window[0] - 2` (window (4,6) -> 2); window (2,5) gives
wait 0 (`>= 0`; a window starting at 1 gives wait -1: assert a clear error or clamp to 0, pin it); no window: wait 1.
`a.trace` unset: no TraceStep even with a window. Both flags absent: callback list equals the old list exactly.

**T59 FastGRPO range wiring.** U. Patch `nvtx.range` with a recorder; run `FastGRPO.timed("rewards")` (the tiny offline TRL
model from conftest): events `[enter rewards, exit rewards]` and `time/rewards_s` appended with a positive float. Run
`training_step`, `get_batch_samples`, `log`: ranges `micro`, `loader`, `log` entered once each, return values identical
to the same call with NVTX disabled (compare tensors with `torch.equal`). Range names used in the code are all in the fixed
list (grep the source with `re` for `range("` / `timed("` literals and assert subset of the fixed 17 names). Catches a typo
name (`rewad`) that nsys analysis would silently miss.

**T60 NVTX overhead gate.** Perf. Disabled `range()` costs under 2 microseconds per call (timeit 100k calls; assert
< 5 us on a loaded CI Mac); `enabled()` after caching under 1 us. It must add no `torch.cuda.synchronize`: patch
`torch.cuda.synchronize` and `torch.Tensor.item`/`.cpu` to count and run `log()` with the counters on a tiny model: 0 new
host reads versus the baseline (count baseline first, assert `new == baseline`).

**T61 ratio diagnostics.** U/P. `raw_ratio` + `bounded` equal `is_ratio` for all four modes on 300 random tensors
(shapes (B, T), random masks including an all-zero mask row and an all-ones row, logps differences in [-20, 20] to include
`exp` overflow to inf and underflow to 0; `torch.equal` where both finite, and `isinf` positions equal). Counts: on a toy
tensor with known values `zero_share`, `hi_share`, `lo_share` equal numpy counts: raw = [0, 0.05, 0.5, 1, 3.5, 5, inf],
clip_min 0.1, clip_max 3: lo = {0, 0.05} count 2 (0 counted once), hi = {3.5, 5, inf} 3, zero = {0} 1 in `token_mask`; in
`sequence_*` modes the base is sequence rows and the denominator is rows. In `*_truncate` the final ratio is clipped
but never zero except raw 0 -> `zero_share` by final ratio. In `*_mask` modes the final ratio is 0 above clip_max, so
`zero_share = hi_share + (raw == 0 share)` minus overlap; assert this identity. A raw of exactly `clip_min` or `clip_max`
is neither lo nor hi (strict). Empty mask (no counted entries): share is 0.0 or NaN? Pin it, no ZeroDivisionError. The
shares add no host read: run on a CUDA tensor in the G test; the CPU variant asserts the return values are tensors, not
Python floats.

**T62 `time/step_s`.** U. Fake monotonic clock: log calls at t = 10, 40, 75; the first has no `time/step_s` (or 0/None; pin),
the second 30.0, the third 35.0. It must not call `torch.cuda.synchronize`. Two logs with evaluation logs in between
(HF logs `eval_*` through the same method) are not mixed: the interval is between consecutive calls as the spec says;
a flagged limitation, assert the literal spec.

---------------------------------------------------------------------------------------------------------------
## 9. Trace analysis [3.7a]

**T63 `profile.covered`.** U/P. Hand cases: disjoint [(0,1),(2,3)] on [0,10] -> 2.0; overlapping [(0,2),(1,3)] -> 3.0;
nested [(0,10),(2,3)] -> 10; clipped [(−5,2),(8,20)] on [0,10] -> 4; empty iterable -> 0.0; zero-length interval 0;
interval entirely outside -> 0; t1 == t0 -> 0.0; unsorted input equals sorted input. Property: for 500 random interval
sets, `covered` equals the area of the rasterised union on a fine grid (1e-3 step, tolerance 2e-3 x n), is <= t1 - t0,
and `covered(A + B) <= covered(A) + covered(B)` and >= max of the two. Generator input (an iterator consumed once) works.
`busy` output on an existing fixture equals the pre-change output, `==`.

**T64 nsys CSV fixture, hand values.** I. Build in `tmp_path` the four CSVs with a 2 step ranges (`step` 0..1000 us and
1000..2000 us, same thread), inside the first: `loader` 0..300, `micro` 300..900 (coverage 0.9), and inside `micro` a
`forward` 400..800 (must not count: only the ranges directly inside `step`). Kernels: 100..200, 150..400 (overlap),
700..1100 (crosses the step boundary: clipped to the first step at 1000 and also counted in the second: 1000..1100).
Memcpy and memset intervals included. API trace rows: two `cudaStreamSynchronize_v3020`, one `cudaMemcpyAsync` (not
counted), one `cudaMemcpy` (counted, name suffix removed), one `cudaEventSynchronize_v3020` inside step 1, one
`cudaDeviceSynchronize` outside any step. Expect per instance: `coverage == 0.9`, `gpu_busy == union/1000` computed by
hand, `gpu_idle == 1 - gpu_busy`, `host_syncs == {"cudaStreamSynchronize": 2, "cudaMemcpy": 1, "cudaEventSynchronize": 1}`
(no `_v3020`), `per_range_s` totals in seconds (us / 1e6; pin the unit), `gpu_proj_s` from the projection CSV.
Overlapping kernels on different streams must be counted once (union, not sum). Catches double counting, sums over 1.0,
unit errors, and counting `Async` copies.

**T65 nesting by thread.** U. Two threads each with `step`: a `loader` range on thread B inside the time of thread A's `step`
is not counted for A (same-thread containment only). Ranges that straddle the step end are clipped or excluded (pin:
spec says "nesting by time containment", so a straddler is excluded; assert coverage excludes it). Back-to-back identical
timestamps (child starts exactly at the parent start, ends at parent end) are contained (closed interval).

**T66 overhead.** U. `--log` JSON with `time/step_s` per step: steps 1..10 values 10 except window steps 4, 5 at 12
(20 percent slower): `overhead_affected is True`; at 10.9 (9 percent): False; at exactly 11.0 (10 percent): False (the
spec: more than 10 percent, strict; float 11.0/10 = 1.1 exactly? compute `1.1000000000000001`: pin by integer
arithmetic: `11.0 > 10.0 * 1.1` evaluates False; assert the literal behaviour). The reference is the steps just before and
after (3 and 6), not the whole run. Missing `time/step_s` in the log: clear error, not KeyError traceback.

**T67 malformed or empty CSVs.** F. Missing file -> exit non-zero naming the file. Header only -> `[]` step entries and
valid JSON. A `step` instance with zero kernels: `gpu_busy == 0.0`, `gpu_idle == 1.0`, no divide by zero. A `step` with
zero duration: skipped or NaN-free. Column names with `nsys` version differences (`Start (ns)` vs `Start(ns)`): fail with a
clear message rather than reading a wrong column; assert one real header from the nsys 2025 docs.

---------------------------------------------------------------------------------------------------------------
## 10. TRL check [3.8]

**T68 passes on the pinned TRL.** I. `check()` returns None on trl 1.17.1 with no GPU and with `vllm` made unimportable
(`sys.modules["vllm"] = None`); the call completes in under 1.0 s (it uses `inspect` on source, no model load). Run twice:
identical. Skip with a clear reason if installed trl != 1.17.1 (CI pins it; assert in a separate test that the dev
environment has the pin so the skip cannot hide a silent drift).

**T69 each failure names its item (adversarial, parametrised).** F. Using `monkeypatch` on the real classes (patch,
never edit the library): (1) delete `GRPOTrainer._prepare_inputs` -> message contains `_prepare_inputs`; (2) wrap
`_compute_loss` with a function of signature `(self, model, inputs, extra)` -> contains `_compute_loss`; (3) swap two
parameter names in `_get_per_token_logps_and_entropies` (order matters) -> names it; (4) delete `shuffle_sequence_dict` ->
names it; (5) make `inspect.getsource` of `_generate_and_score_completions` return text without `"tool_mask"` -> names
`tool_mask`; same for each of the 4 other quoted keys and for the substring `self.use_vllm and
self.vllm_importance_sampling_correction`; (6) source of `__init__` without `self._metrics`; `VLLMGeneration` without
`self._dist`, `def sync_weights`, `self.llm`; (7) each of the 10 `self.<name> =` attributes removed in turn (loop). Each
raises `TRLContractError` (a `RuntimeError`) and the message contains the failing name exactly. Order test: break items 1 and 5 at once,
the message names item 1 (first failure). Catches a check that compares only the version, and one that reports a
generic message.

**T70 version rule.** U. `trl.__version__ = "1.18.0"` with items 1-5 passing: raises, message contains `1.18.0` and the
pin; with `OBJDEL_TRL_UNPINNED=1`: returns None and exactly one `warnings.warn` (use `pytest.warns` count == 1) that names `1.18.0`; with
`OBJDEL_TRL_UNPINNED=0` or `"true"`: raises (only `"1"` counts). Version changed and item 3 broken together with
UNPINNED=1: still raises (items 1-5 are never waived). A dev version string `1.17.1.dev0` or `1.17.1+cu128`: does not equal the PIN,
so it raises; pin that (strict equality), and note the Mac wheel sometimes has local suffixes.

**T71 TRL's ratio block equals `fast.is_ratio` (the key coupling).** P. Extract the block from the TRL source by the two
anchor lines (the spec's start and end anchors; assert both anchors are found exactly once, else fail with the anchor
text, so a TRL change fails loudly), `textwrap.dedent` and `exec` it in a namespace with tensors
`old_per_token_logps, sampling_per_token_logps, mask (completion_mask), self` stub (`vllm_importance_sampling_mode`,
`..._clip_min`, `..._clip_max`). For 4 modes x 200 random inputs (include rows with all-zero mask, NaN-free logps,
differences up to +-15) compare to `fast.is_ratio` with `torch.equal` (same ops, same order, so bit equality is
expected; allow `torch.allclose(rtol=0, atol=0)` and flag any 1-ulp drift in the review as a real finding). Catches
a TRL release that re-orders the ops (e.g. `exp(a-b)` vs `exp(a)/exp(b)`), where the diagnostic share counts would drift.

**T72 `FastGRPO.__init__` calls `check()` first.** U. Patch `trl_adapter.check` to raise `TRLContractError`; patch every
other init dependency (`GRPOTrainer.__init__`) to record a call. Constructing `FastGRPO` raises `TRLContractError` and the
recorder is empty (nothing ran before). Patch `check` to record the call order relative to `GRPOTrainer.__init__`: check is first.

**T73 preflight entry.** U. `preflight.CHECKS` contains `"CODE TRL private names"`, it is the first entry whose name starts with
`CODE`; `fn(a, cfg)` returns `[]` normally, and with a broken `check` returns `[msg]` where msg is the string of the error
(list of length exactly 1). A non-`TRLContractError` exception (e.g. `ImportError` because trl is missing) also becomes a
one-element list, not a crash (catch `Exception`? the spec says "when it raises": assert that behaviour).

---------------------------------------------------------------------------------------------------------------
## 11. Outline lengths [3.9]

**T74 `count` with a char-stub.** U. `tokens_of = lambda s: len(s) + 1`... the CLI adds the +1, so the library `count` must NOT
add it: `tokens_of = len` gives exactly `len(text)`. Records: 3 records, targets `["abcd", None, "ab"]`, `["a"*300]`, `[None]`;
`n_records == 3`, `n_targets == 3` (None skipped), lengths [4, 2, 300]; `mean == 102.0`, `max == 300`, `above[256] == 1`,
`above[512] == 0`; `above` has exactly the 6 keys; strict `>`: a length exactly 256 is NOT above 256. Voting callable returns
`{0}` for record 1: `n_voting_targets` counts only index-0 targets that are non-None; `above_voting` consistent;
a voting set containing an index whose target is None counts nothing. `longest` sorted decreasing, ties stable, capped at 20
with 25 targets, entries `{item_id, rater, tokens}` with rater the target index. Percentile definition pinned with `np.percentile`
default linear method on a 5-value list (p50 3.0; p90 4.6).

**T75 empty and degenerate inputs.** F. No records: `n_records 0`, mean etc. None or 0.0 and not a crash on `np.percentile([])`
(pin: spec silent; require no exception, JSON-safe). All-None targets: `n_targets 0`. `record` without a `targets` key:
clear `KeyError`/`ValueError` naming the key. `voting=None`: `n_voting_targets == n_targets` or null (pin).
`tokens_of` returns numpy ints: JSON serialisable (`int()` cast) -> `json.dumps` passes.

**T76 `answers` join.** U. 4 answers with tokens [100, 300, 600, 1600], metric rows keyed by `item_id` with per-sample lists:
the answer rows hold K entries per item; join by `item_id` and then by index inside the list. Hand fixture: K = 2, 2 items.
For L = 256: share_above 0.5; `legal_share`, `mean_iou_tc` over those answers; `n_legal_good` counts legal and IoU_tc >= 0.5
(0.5 inclusive: an answer with exactly 0.5 counts). For L = 1536 with one answer of 1600 tokens: values for that answer only.
No answer above L: `share_above == 0.0`, and `legal_share`, `mean_iou_tc` are `None` (no division by zero), `n_legal_good == 0`.
A metric row missing for an answer item: `KeyError` naming the id or an explicit skip count; a silent drop is a bug (pin).
Different K between `answers.jsonl` tokens and metrics lists: `ValueError`.

**T77 real data counts [box data].** G/E (skipped without `data/`). `data/prepared` LIDC train: 7,341 records, 23,469 targets, max 1,182,
`above[1200] == 0`. Views rows carry `manifest_sha256` equal to sha256 of the file bytes. Time: the CLI over train under 120 s CPU.
Hard assertion that the count of targets equals the independent sum over `len(rec["targets"])` minus None entries.

**T78 `target_tokens` and `voting` replace the privates.** U. `objdel.train.dataset` has `target_tokens` and `voting`; neither
`_tokens` nor `_voting` exists (`hasattr` False; no alias); `grep` of the repo (via `pathlib` walk over `src`, `experiments`, `tests`) finds no
reference to the old names. Behaviour equal to the old functions on a frozen set of 20 records (store expected outputs
as literals before the rename).

---------------------------------------------------------------------------------------------------------------
## 12. CPU versus GPU exactness, performance, system tests

**T79 [GPU] CPU/GPU bit equality of the checker terms.** G. 500 random items (golden masks from `tests/oracles`, plus
adversarial polygons: self-touching, hole-containing, frame-touching, single pixel, 1000 x 1000 grid corners) scored
through `objdel._gpucheck` GPU path and CPU build: `legal, iou_tc, f5, g_area` equal with `==` per sample, and the
`sample_iou` matrix (if computed on device) equal. The Mac rescore of box `answers.jsonl` matches the box metrics (gate 6).
The pixels-stay-on-device rule: while scoring, patch `torch.Tensor.cpu` and `.numpy` to record sizes: no transfer larger
than `K * 8` bytes of scores per item for the checker terms (boundary IoU and `sample_iou` use host masks and are exempt:
list them explicitly in the test so an extra exemption is a visible diff).

**T80 [GPU] boundary_iou on device equals host.** G. If boundary IoU is implemented in the CUDA build, 300 random masks
(T12 generator): `==` with the host reference. Edge masks and d = 16 on a 512 frame included.

**T81 performance gates, Mac.** Perf. (a) `boundary_iou` 512 x 512, d = 9, 8 answers x 4 raters x 1 object: under 20 ms per item
(the metrics path must use a bounding box crop plus the 1-pixel pad: assert equality with the full frame T13 and speed
> 3x the full frame on a 40 x 40 nodule). (b) `metrics()` for K = 8, k = 4: new keys add under 25 percent to the old time
(measure old keys via a flag or the frozen copy; assert ratio < 1.25). (c) Rescore of 300 synthetic items < 60 s. (d)
`sample_iou` K = 8 costs one `_pairwise_iou` already computed (spy: `_pairwise_iou` called once per object per item, not twice). (e)
Memory: `answers.jsonl` of 2,200 items x 8 answers x ~1,000 chars is ~18 MB: reading it line by line keeps peak RSS under 200 MB
(`tracemalloc` over the rescore reader).

**T82 [box] performance and cost gates.** G. nsys run: `step` ranges exist for steps 4 and 5; `coverage >= 0.90` each; the NVTX
overhead gate: `time/step_s` with `OBJDEL_NVTX=1` and no profiler within 3 percent of `OBJDEL_NVTX` unset (R5 baseline 23.6 s/step
[measured: memory R5]); the nsys capture ends within 900 s + 120 s (the wrapper uses `timeout`; test the wrapper command line
contains `timeout 1020` or two stages); total spend assertion in the preflight output <= $1.00 and North America host.

**T83 E2E mini pipeline on the Mac.** E. Synthetic dataset of 12 items (3 patients, 2 nodules each, 2 slices), a fake
sampler that returns for each item: [truth text, truth text shifted by 3 px, garbage, empty polygon, truth, truth, small blob, truth]. Run final
eval full path (`main` with the sampler patched), then `ceiling` CLI with `--answers`, then `outline_lengths` with a stub tokenizer. Invariants across
files: items count equal everywhere; `final_eval.json["items"] == 12`, `nodules == 6`, `patients == 3`; `legal` mean == 7/8
(hand: garbage illegal; the empty polygon is illegal or legal per checker rule, pin with the checker); `sel_oracle >= sel_mean_of_k` for
every item; the medoid of the example picks a "truth" copy (4 or 5 copies agree) so `sel_medoid == iou of truth` and `paired.medoid_minus_mean > 0`;
the ceiling's `diff` with answers equal to doctors is >= 0 on the entry-r term. Everything JSON round-trips with `allow_nan=False`.

**T84 repeatability.** E. The whole Mac pipeline twice in a fresh subprocess with `PYTHONHASHSEED` different: `final_eval.json`
files equal after removing `seconds` and `code_sha`. Catches dict or set iteration order in clusters, thread count in the CPU
checker (set `OMP_NUM_THREADS` 1 vs 4 and compare), and use of the global RNG.

---------------------------------------------------------------------------------------------------------------
## Pytest sketches: the 5 most important tests

Ranked: (1) golden CI regression plus legacy bit identity (T01, T02), (2) write order and rescore bit equality (T43, T45, T47),
(3) boundary IoU full-frame equality at the edges (T12, T13), (4) TRL ratio block equals `fast.is_ratio` and per-item failure names (T71, T69),
(5) NvtxSteps balance and profiler window (T56, T57). Also strongly worth it: T21 (medoid eligible-only mean) and T64 (nsys hand values).

```python
# ---- 1. tests/test_a_signal_golden.py -------------------------------------------------------------------
import gzip, json
from pathlib import Path
import numpy as np
from experiments.phase2 import final_eval as FE
from objdel.eval import signal as S

GOLD = Path(__file__).parent / "goldens" / "r5_test_items.jsonl.gz"

def _legacy_draws(ids, boot=2000, seed=0):            # frozen copy of the pre-change implementation
    nod = [S.cluster(i) for i in ids]
    groups = [np.array([j for j, x in enumerate(nod) if x == u]) for u in dict.fromkeys(nod)]
    rng = np.random.default_rng(seed)
    return [np.concatenate([groups[g] for g in rng.integers(0, len(groups), len(groups))]) for _ in range(boot)]

def _gold():
    rows = [json.loads(l) for l in gzip.open(GOLD, "rt")]
    return [r["item_id"] for r in rows], np.array([r["iou_tc"] for r in rows])

def test_golden_intervals_are_exact():
    ids, v = _gold()
    assert len(ids) == 2200 and len({S.cluster(i) for i in ids}) == 401
    items = [{"item_id": i, "iou_tc": x, "miss": float(x < 0.1)} for i, x in zip(ids, v)]
    s = FE.summarize_golden_subset(items)             # helper: runs summarize with the other metrics zero-filled
    assert s["iou_tc"]["ci95_item"] == [0.6334010210209712, 0.6561989221793691]
    assert s["iou_tc"]["ci95_nodule"] == [0.6279017988891661, 0.6620867640136667]
    assert s["miss_rate"]["mean"] == 0.10818181818181818

def test_draws_bit_identical_to_legacy_for_random_id_lists():
    rng = np.random.default_rng(0)
    for t in range(50):
        n_nod = int(rng.integers(1, 30))
        ids = [f"lidc/s{u}/n{u % 3:03d}/z{z}" for u in range(n_nod) for z in range(int(rng.integers(1, 6)))]
        ids += [f"other{j}" for j in range(int(rng.integers(0, 4)))]
        rng.shuffle(ids)
        for boot, seed in [(1, 0), (7, 1), (100, 12345)]:
            new, old = S.draws(ids, boot, seed), _legacy_draws(ids, boot, seed)
            assert len(new) == len(old) and all(np.array_equal(a, b) and a.dtype == b.dtype for a, b in zip(new, old))

# ---- 2. tests/test_a_final_eval_io.py ---------------------------------------------------------------------
import hashlib, json, pytest
from experiments.phase2 import final_eval as FE

def _run(tmp_path, monkeypatch, scorer=None):
    view = tmp_path / "view"; view.mkdir(); (view / "manifest.json").write_bytes(b'{"n": 3}\n')
    monkeypatch.setattr(FE, "sample_items", lambda *a, **k: FAKE_ANSWERS)       # 3 items x K=8, fixed texts
    if scorer: monkeypatch.setattr(FE, "score_items", scorer)
    out = tmp_path / "out"
    FE.main(["--ckpt", "c", "--view", str(view), "--split", "val", "--train", "t", "--out", str(out)])
    return view, out

def test_answers_are_on_disk_before_scoring_fails(tmp_path, monkeypatch):
    def boom(*a, **k): raise RuntimeError("boom")
    with pytest.raises(RuntimeError):
        _run(tmp_path, monkeypatch, scorer=boom)
    out = tmp_path / "out"
    lines = [json.loads(l) for l in (out / "answers.jsonl").read_text().splitlines()]
    assert len(lines) == 3 and all(len(r[k]) == 8 for r in lines for k in ("samples", "truncated", "tokens"))
    fe = json.loads((out / "final_eval.json").read_text())
    assert fe["partial"] is True and fe["view_manifest_sha256"] == hashlib.sha256(b'{"n": 3}\n').hexdigest()
    assert not (out / "metrics.jsonl").exists() and not (out / "items.jsonl").exists()

def test_rescore_is_bit_identical_and_refuses_a_changed_manifest(tmp_path, monkeypatch):
    view, out = _run(tmp_path, monkeypatch)
    out2 = tmp_path / "out2"
    assert FE.main(["--rescore", str(out), "--view", str(view), "--split", "val", "--out", str(out2)]) in (None, 0)
    for a, b in zip(open(out / "metrics.jsonl"), open(out2 / "metrics.jsonl")):
        sa, sb = json.loads(a)["samples"], json.loads(b)["samples"]
        assert all(sa[k] == sb[k] for k in ("legal", "iou_tc", "f5", "g_area"))      # exact
    (view / "manifest.json").write_bytes(b'{"n": 4}\n')
    out3 = tmp_path / "out3"
    with pytest.raises(SystemExit) as e:
        FE.main(["--rescore", str(out), "--view", str(view), "--split", "val", "--out", str(out3)])
    assert e.value.code == 2 and (not out3.exists() or not any(out3.iterdir()))

# ---- 3. tests/test_a_boundary_iou.py ------------------------------------------------------------------------
import cv2, numpy as np
from objdel.eval import metrics as M

def _ref_band(m, d):
    p = cv2.copyMakeBorder(m.astype(np.uint8), 1, 1, 1, 1, cv2.BORDER_CONSTANT, value=0)
    e = cv2.erode(p, np.ones((3, 3), np.uint8), iterations=d)[1:-1, 1:-1].astype(bool)
    return m & ~e

def _ref_biou(a, b, d):
    A, B = _ref_band(a, d), _ref_band(b, d)
    u = int((A | B).sum())
    return 1.0 if u == 0 else int((A & B).sum()) / u

def test_band_hand_values():
    m = np.zeros((9, 9), bool); m[2:7, 2:7] = True
    assert M.band(m, 1).sum() == 16 and M.band(m, 2).sum() == 24 and M.band(m, 3).sum() == 25
    full = np.ones((8, 8), bool)
    assert M.band(full, 1).sum() == 4 * 7                    # frame edge is background, so the ring survives

def test_boundary_iou_equals_full_frame_reference_on_edge_touching_masks():
    rng = np.random.default_rng(0)
    for _ in range(500):
        h, w = [(1, 1), (1, 7), (7, 1), (16, 16), (64, 48)][rng.integers(5)]
        dens = [0.02, 0.5, 0.98][rng.integers(3)]
        a, b = rng.random((h, w)) < dens, rng.random((h, w)) < dens
        a[0, :] = True; b[:, -1] = True; a[-1, 0] = True           # touch the frame
        d = [1, 2, 9, 16][rng.integers(4)]
        assert M.boundary_iou(a, b, d) == _ref_biou(a, b, d)
        assert M.boundary_iou(a, b, d) == M.boundary_iou(b, a, d)
    z = np.zeros((5, 5), bool)
    assert M.boundary_iou(z, z, 2) == 1.0 and M.boundary_iou(z, ~z, 2) == 0.0

# ---- 4. tests/test_a_trl_adapter.py -------------------------------------------------------------------------
import inspect, textwrap, pytest, torch, trl
from objdel.train import fast, trl_adapter as T
from trl import GRPOTrainer

def _block():
    src = inspect.getsource(GRPOTrainer._generate_and_score_completions)
    a = "per_token_logps_diff = (old_per_token_logps - sampling_per_token_logps) * mask"
    b = "Unknown vLLM importance sampling level"
    assert src.count(a) == 1 and src.count(b) == 1, "TRL ratio block anchors moved"
    lines = src.splitlines(); i0 = next(i for i, l in enumerate(lines) if a in l)
    i1 = next(i for i, l in enumerate(lines) if b in l)
    while "else:" not in lines[i1 - 1]: i1 -= 1
    return textwrap.dedent("\n".join(lines[i0:i1 - 1]))

@pytest.mark.parametrize("mode", ["sequence_mask", "sequence_truncate", "token_mask", "token_truncate"])
def test_trl_ratio_block_equals_fast_is_ratio(mode):
    code = _block(); g = torch.Generator().manual_seed(0)
    for _ in range(100):
        old = torch.randn(4, 12, generator=g) * 5; samp = old + torch.randn(4, 12, generator=g) * 5
        mask = (torch.rand(4, 12, generator=g) < 0.8).float(); mask[0] = 0.0
        ns = dict(old_per_token_logps=old, sampling_per_token_logps=samp, mask=mask, torch=torch,
                  self=type("S", (), dict(vllm_importance_sampling_mode=mode,
                       vllm_importance_sampling_clip_min=0.1, vllm_importance_sampling_clip_max=3.0))())
        exec(code, ns)
        ref = fast.is_ratio(old, samp, mask, mode, 0.1, 3.0)
        assert torch.equal(ns["importance_sampling_ratio"], ref)

@pytest.mark.parametrize("breaker, needle", [
    (lambda m: m.delattr(GRPOTrainer, "_prepare_inputs"), "_prepare_inputs"),
    (lambda m: m.setattr(trl, "__version__", "9.9.9"), "9.9.9"),
])
def test_check_names_the_failing_item(monkeypatch, breaker, needle):
    monkeypatch.delenv("OBJDEL_TRL_UNPINNED", raising=False)
    breaker(monkeypatch)
    with pytest.raises(T.TRLContractError, match=needle):
        T.check()

def test_unpinned_env_only_waives_the_version(monkeypatch):
    monkeypatch.setenv("OBJDEL_TRL_UNPINNED", "1"); monkeypatch.setattr(trl, "__version__", "9.9.9")
    with pytest.warns(UserWarning, match="9.9.9"):
        T.check()
    monkeypatch.delattr(GRPOTrainer, "_prepare_inputs")
    with pytest.raises(T.TRLContractError):
        T.check()

# ---- 5. tests/test_a_nvtx_steps.py ----------------------------------------------------------------------------
import types
from objdel.train import nvtx
from objdel.train.callbacks import NvtxSteps

class Rec:
    def __init__(self): self.ev = []
    def range_push(self, n): self.ev.append(("push", n))
    def range_pop(self): self.ev.append(("pop",))
    def start(self): self.ev.append(("start",))
    def stop(self): self.ev.append(("stop",))

@pytest.fixture
def rec(monkeypatch):
    r = Rec()
    monkeypatch.setattr(nvtx, "enabled", lambda: True)
    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setattr(torch.cuda.nvtx, "range_push", r.range_push, raising=False)
    monkeypatch.setattr(torch.cuda.nvtx, "range_pop", r.range_pop, raising=False)
    monkeypatch.setattr(torch.cuda.profiler, "start", r.start, raising=False)
    monkeypatch.setattr(torch.cuda.profiler, "stop", r.stop, raising=False)
    return r

def _drive(cb, n):
    st = types.SimpleNamespace(global_step=0)
    cb.on_train_begin(None, st, None)
    for i in range(1, n + 1):
        cb.on_pre_optimizer_step(None, st, None); cb.on_optimizer_step(None, st, None)
        st.global_step = i
        cb.on_step_end(None, st, None)
    cb.on_train_end(None, st, None)

@pytest.mark.parametrize("n", [0, 1, 2, 10, 30])
def test_pushes_and_pops_balance_and_nest(rec, n):
    _drive(NvtxSteps(None), n)
    depth = 0
    for e in rec.ev:
        depth += {"push": 1, "pop": -1}.get(e[0], 0)
        assert 0 <= depth <= 2
    assert depth == 0 and not any(e[0] in ("start", "stop") for e in rec.ev)

def test_profiler_window_starts_after_pop_before_push_and_stops_once(rec):
    _drive(NvtxSteps((4, 6)), 8)
    names = [e[0] for e in rec.ev]
    assert names.count("start") == 1 and names.count("stop") == 1
    i = names.index("start")
    assert names[i - 1] == "pop" and names[i + 1] == "push"        # after the pop, before the next push
    j = names.index("stop"); assert names[j - 1] == "pop" and j > i
    steps_before_start = [e for e in rec.ev[:i] if e == ("push", "step")]
    assert len(steps_before_start) == 4                              # started at the end of global_step 3 (3 closes + initial)
```

Notes for the implementer of the tests (known risks, flagged because they may be spec gaps):
1. T57: a run that ends before step `e` leaves the profiler started; the spec has no `on_train_end` stop rule.
2. T49: `git diff HEAD` ignores untracked files, so a new untracked module does not mark `-dirty`.
3. T58: window start 1 gives `wait = -1`; the spec does not say clamp or reject.
4. T11: the 1-pixel pad for d > 1 is claimed equal to a d-pixel pad; the test proves or disproves it.
5. T66: "more than 10 %" at exactly 10 % depends on float rounding; pin it by test.
6. Sketches use helper names (`sample_items`, `score_items`, `FAKE_ANSWERS`, `summarize_golden_subset`) that do not exist yet; they mark the seam the implementer must provide or the test must adapt to.
