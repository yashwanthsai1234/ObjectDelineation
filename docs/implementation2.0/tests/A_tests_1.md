# Phase A tests, set 1 (written from A.spec.md only)

## Skills used
- nvidia-skills:nvidia-skill-finder: checked for a CUDA-testing skill; no change (tests use the CPU build `objdel._gpucheck` as the spec says).
- engineering:testing-strategy: levels (unit/property/integration/GPU/e2e), gates, and the "Mac vs [box]" split.
- mattpocock-skills:tdd: tests assert behaviour at public interfaces only (spec section 3), never internals; each test is red before the feature exists.
- engineering:code-review / debug / system-design / architecture / documentation / huggingface-best: used as checklists (failure modes, hidden couplings, file-order contracts, doc clarity). No other change.

Conventions. Marks: `[mac]` runs on the Mac with no GPU/data; `[box]` needs GPU or `data/`. Level in brackets.
Each test: name, level, setup, input, expected, why it can catch a bug. Fixture words: "mask" = bool (H,W) array.
Test files proposed: `tests/test_a_signal.py`, `test_a_metrics.py`, `test_a_select.py`, `test_a_ceiling.py`, `test_a_final_eval.py`, `test_a_nvtx.py`, `test_a_trl_adapter.py`, `test_a_outline_lengths.py`, `test_a_profile.py`.

---------------------------------------------------------------------------------------------------------------
## 1. Cluster bootstrap `signal.draws` (spec 3.1, gate 1, 2)

**T1.1 default-key bit-identity to the legacy implementation** [unit, mac]
- Setup: copy the OLD `draws` body into the test as `legacy_draws` (verbatim from today's signal.py: `nod=[cluster(i)...]`, `groups`, `rng.integers(0, len(groups), len(groups))` per draw).
- Input: ids = the 2,200 ids of `tests/goldens/r5_test_items.jsonl.gz`; also a 7-id toy with 3 clusters; boot=50, 2000; seed 0 and 3.
- Expected: `len(new)==len(old)==boot`; for every b, `np.array_equal(new[b], old[b])` AND same dtype (int64).
- Why: the refactor to a `key` parameter can silently change RNG call order (e.g. drawing all rows at once instead of per draw changes the stream). Gate 1's exact CIs depend on it.

**T1.2 one item per cluster equals plain `integers` rows** [unit, mac]
- Input: ids = `["x0".."x9"]` (non-`lidc/`, so each is its own cluster), boot=100, seed=5.
- Expected: `np.stack(draws(...))` equals `np.random.default_rng(5).integers(0, 10, (100, 10))` exactly.
- Why: the spec states it; it also proves per-draw `integers(0,C,C)` equals one 2-D call (true for numpy Generator; the test pins it so a numpy upgrade or a changed call pattern is caught). Gate 1 sentence "ci95_nodule equals ci95_item exactly" follows from it.

**T1.3 custom key, first-appearance order, variable cluster sizes** [unit, mac]
- Input: ids `["a1","a2","b1","c1","c2","c3"]`, `key=lambda s: s[0]`, boot=1, seed=0. Let g = rng.integers(0,3,3) from a fresh rng(0).
- Expected: the array is `np.concatenate([[0,1],[2],[3,4,5]][i] for i in g)`. Clusters ordered a, b, c (first appearance), not sorted by key.
- Why: using `set()` or `sorted()` of keys would reorder clusters and change results; hand-check catches it.

**T1.4 interleaved ids keep within-cluster index order** [unit, mac]
- Input: ids `["a","b","a","b"]`, key identity. Expected: any draw containing cluster a contains `[0,2]` in that order, b `[1,3]`.
- Why: groups must list indices ascending even when a cluster is not contiguous.

**T1.5 draw length is variable, not n** [unit, mac]
- Input: clusters of sizes 1 and 10 (11 ids). Expected: array lengths range over {2,11,20} only, i.e. `len == sum(size of drawn cluster)`; all indices in [0,n).
- Why: a bug that pads or truncates to n would bias CIs.

**T1.6 empty and single-cluster edge** [unit, mac]
- Input: `ids=[]` ; `ids=["lidc/s/n000/z1","lidc/s/n000/z2"]` (one cluster).
- Expected: `[]` -> each draw an empty array without exception (or a documented ValueError, but NOT an infinite loop / IndexError from `integers(0,0)`); one cluster -> every draw is `[0,1]`.
- Why: `rng.integers(0,0,0)` raises in numpy ("high <= 0"); the spec is silent, so this pins the behaviour. Mark as "pin current behaviour, then keep".

**T1.7 `cluster()` id parsing** [unit, mac]
- Input: `lidc/1.3.6.1/n000/z185` -> `lidc/1.3.6.1/n000`; `riga/x/img1` -> itself; `lidc/` -> `lidc`?? (pin: `"lidc/".rsplit("/",1)[0]=="lidc"`).
- Why: patient key via `groups` must not be confused with `cluster`.

**T1.8 CI coverage property (adversarial: correlated slices)** [property, mac]
- Setup: simulate 200 nodules x 6 slices; per-nodule effect N(0,1) shared by its slices, tiny slice noise. 300 Monte-Carlo replications, boot=400.
- Expected: empirical coverage of the true mean by `ci95_nodule` in [0.90,0.99]; `ci95_item` coverage < 0.80; and `width(nodule) > 1.8 * width(item)`.
- Why: verifies the cluster CI is the wider, honest interval and that nobody swapped the keys `ci95_item` / `ci95_nodule`.

**T1.9 gate 1 exact regression** [integration, mac]
- Input: golden gz, `summarize` on items built from (`item_id`, `iou_tc`) with `groups=None`, seed 0, draws 2000.
- Expected exactly (`==`, not approx): `ci95_item == [0.6334010210209712, 0.6561989221793691]`, `ci95_nodule == [0.6279017988891661, 0.6620867640136667]`, `miss_rate.mean == 0.10818181818181818`, `nodules == 401`, `items == 2200`.
- Why: today's numbers; any change in quantile method, RNG order or float summation order breaks it. Float summation: compute means with `v[idx].mean(1)` exactly as the spec writes it.

---------------------------------------------------------------------------------------------------------------
## 2. final_eval `METRICS`, `per_item`, `summarize` (spec 3.2)

**T2.1 METRICS table is exact** [unit, mac]
- Expected: `set(METRICS) == {legal, iou_tc, f5, precision_vote, recall_vote, dice_vote, nsd, boundary_iou, boundary_iou_vote, soft_dice, hm_iou, ged}` (12 names); `{k for k,(s,h) in METRICS.items() if not h} == {"ged"}`.
- Why: a thirteenth or missing name changes every downstream table.

**T2.2 per_item averages per-sample lists, copies scalars** [unit, mac]
- Input: K=4 hand result: `samples.legal=[1,1,0,1]`, `iou_tc=[.8,.6,0,.4]`, `f5=[.9,.7,0,.5]`; `objects.nodule`: `precision_vote=[.5,.5,0,1]`, `recall_vote=[1,.5,0,.5]`, `dice_vote=[.6,.5,0,.7]`, `nsd=[.2,.4,0,.2]`, `boundary_iou=[.5,.5,0,.5]`, `boundary_iou_vote=[.1,.2,0,.3]`; scalars `soft_dice=.33`, `hm_iou=.44`, `ged=.55`.
- Expected: legal .75, iou_tc .45, f5 .525, precision_vote .5, recall_vote .5, dice_vote .45, nsd .2, boundary_iou .375, boundary_iou_vote .15; soft_dice/hm_iou/ged exactly .33/.44/.55; `miss == 0.0` (0.45 >= 0.1). All values python `float`.
- Why: wrong source (samples vs objects), or averaging the scalar again.

**T2.3 miss threshold is strict at 0.1** [unit, mac]
- Input: K-mean iou_tc = 0.0999999, 0.1 exactly, 0.1000001 (build samples so the mean is exact: `[0.1]*K` for exactly 0.1; note `np.mean([0.1]*8)` may be 0.1 or 0.10000000000000002 - assert the implementation uses `iou_tc < 0.1` on the same float it reports).
- Expected: miss = 1.0, 0.0, 0.0. `miss == float(item["iou_tc"] < 0.1)` for all.
- Why: `<=` vs `<`; computing miss from a different rounding than the reported iou_tc.

**T2.4 per_item uses the `obj` argument** [unit, mac]
- Input: result with objects `disc` and `cup`; call `per_item(r, "cup")`.
- Expected: reads `objects["cup"]`. Default is "nodule" and raises KeyError (not silently 0) for a result with no nodule.
- Why: RIGA path.

**T2.5 K=1 and all-illegal items** [unit, mac]
- Input: K=1; and K=8 all illegal with all zeros.
- Expected: no division errors; `legal==0.0`, `sel_legal_mean==0.0`, `sel_medoid==` IoU_tc of answer 0, `sel_medoid_legal==0.0`, `sel_oracle==0.0`, miss 1.0.
- Why: empty-eligible branch of selection; K=1 medoid has no pairs.

**T2.6 per_item selection keys (hand case)** [unit, mac]
- Input: K=4, `iou_tc=[.9,.8,.1,0]`, `legal=[1,1,1,0]`, `sample_iou` = [[1,.8,.1,0],[.8,1,.2,0],[.1,.2,1,0],[0,0,0,1]].
- Expected: `sel_mean_of_k = 0.45`; `sel_legal_mean = (.9+.8+.1)/3 = 0.6`; medoid scores: i0=(.8+.1)/2=.45, i1=(.8+.2)/2=.5, i2=.15 -> medoid 1 -> `sel_medoid=.8`, `sel_medoid_legal=1.0`; `sel_oracle=.9`. `sel_vote` present only when the vote key is available (assert key absent, not None, when `objects` lack what vote needs, if spec "when available" means absent; pin whichever, document).
- Why: medoid must ignore the illegal answer 3 even though its IoU row with others is 0.

**T2.7 selection key invariants (property)** [property, mac]
- Random per-item results (50 items, K=8). Expected for each: `sel_oracle >= sel_medoid`, `sel_oracle >= sel_legal_mean`, `sel_oracle >= sel_mean_of_k`; `sel_legal_mean >= sel_mean_of_k` (illegal score 0 so removing zeros cannot lower the mean); `sel_medoid_legal in {0.0,1.0}`; every selection value in [0,1].
- Why: cheap algebraic checks that catch mix-ups of keys.

**T2.8 summarize schema** [unit, mac]
- Input: 20 items with groups for 8 patients. Expected: `schema == 2`, `boot == draws`, `ci_level == 0.95`, `items==20`, `nodules` and `patients` counts right; for every METRICS key and selection key: dict keys exactly `{"mean","ci95_nodule","ci95_patient","ci95_item","higher_is_better"}`; `"ci95" not in` any entry; `iou_tc["median"]` equals `np.median`; `miss_rate` and `iou_found` present; `paired` has `medoid_minus_mean`, `vote_minus_mean`.
- Without `groups`: `patients is None` and every `ci95_patient is None`.
- Why: stale `ci95` consumers must fail loudly; null contract.

**T2.9 CI definitions, recomputed independently** [unit, mac]
- Input: 30 items, random iou_tc, ids in 10 nodules, groups mapping to 4 patients; seed 0, draws 300.
- Expected (exact): `ci95_item == np.percentile(v[idx].mean(1),[2.5,97.5])` with `idx = rng(0).integers(0,n,(300,n))`; `ci95_nodule` from `signal.draws(ids,300,0)` with percentile 2.5/97.5; `ci95_patient` from `draws(ids,300,0,key=groups.__getitem__)`. Note: spec says "quantiles"; assert equality to `np.quantile(...,[.025,.975])` (default linear = percentile). `mean == v.mean()`.
- Why: the three CIs share a seed but use different draws; mixing them is the likeliest bug.

**T2.10 `ci95_nodule` collapses to `ci95_item` with unique clusters** [unit, mac]
- Input: ids non-lidc (each its own cluster). Expected exact equality of the two CI lists, for every metric.

**T2.11 patient CI wider than nodule CI when patients have several nodules (sanity)** [property, mac]
- Input: synthetic: 40 patients x 3 nodules, patient effect large. Expected `width(ci95_patient) > width(ci95_nodule)` for iou_tc.

**T2.12 miss_rate and iou_found** [unit, mac]
- Input: iou_tc `[0.05,0.5,0.7,0.0]` (clusters distinct). Expected `miss_rate.mean == 0.5`; `iou_found.mean == 0.6`; its CI is computed on the found subset only with draws over the found items' ids (assert by recomputing). All found/none found: no found items -> `iou_found.mean is None` (or documented NaN); must not raise ZeroDivision/warn-crash.
- Why: subset CI built with indices from the full set (off-by-subset) is a classic bug.

**T2.13 `paired` differences** [unit, mac]
- Input: items where `sel_medoid - sel_mean_of_k` = `[.1,.1,-.05,0]`; no `sel_vote` anywhere.
- Expected: `medoid_minus_mean[0] == 0.0125`; lo<=d<=hi; `vote_minus_mean is None`. With `sel_vote` present: list of 3 floats. Mixed presence (some items missing `sel_vote`): raise or null - pin, never average over a subset silently.

**T2.14 NaN / non-finite input** [adversarial, mac]
- Input: one item with `ged=nan`. Expected: summarize raises ValueError naming the item, or yields mean nan - NOT a silently finite CI (np.percentile of nan gives nan; assert not a finite number). Pin the decision.
- Why: a NaN from an empty rater set must not hide.

**T2.15 `higher_is_better`** [unit, mac] `summarize(...)["ged"]["higher_is_better"] is False`, all others True; selection keys True.

---------------------------------------------------------------------------------------------------------------
## 3. final_eval CLI, answers, rescore (spec 3.2, gate 6, 7)

**T3.1 write order: answers before scoring** [integration, mac, monkeypatched sampler]
- Setup: monkeypatch the model/sampler layer (the part marked [box]) to return fixed answers for a tiny fixture view (3 items); monkeypatch the scoring function to raise `RuntimeError("boom")`.
- Expected: run raises; `answers.jsonl` exists with 3 lines, schema `{item_id, group_id, W, H, samples(K), truncated(K), tokens(K)}`; `final_eval.json` exists with `partial == true`, has `ckpt, view, split, K, seed, view_manifest_sha256, code_sha, checker_build`, `seconds.load_and_sample`; `metrics.jsonl`, `items.jsonl` absent.
- Why: the whole point of the order: GPU minutes of sampling must survive a scoring crash.

**T3.2 full run replaces partial** [integration, mac]
- Expected after success: 4 files; `final_eval.json["partial"] is False`, `schema==2`, still carries ckpt/view/split/K/seed/seconds/sha fields; `metrics.jsonl` line count == items; `items.jsonl` lines equal `per_item` of each metrics line (recompute).

**T3.3 answers.jsonl lengths** [unit, mac] each line: `len(samples)==len(truncated)==len(tokens)==K`; `tokens` are ints; `truncated` bools; `W,H` ints.

**T3.4 rescore bit-equality** [integration, mac]
- Setup: fixture view (3-6 synthetic items from conftest `_disc`), answers = valid polygon texts, one illegal text, one unparseable "garbage", one empty `"[]"`.
- Direct-score vs `--rescore DIR --view V --split S --out O`.
- Expected: for every item, per-sample `legal, iou_tc, f5, g_area` compared with `==` on the python floats (and `np.array_equal`), not approx; `items.jsonl` byte-identical modulo `seconds`.
- Why: a rescore path that re-derives caches in a different order or with a different tau would drift in the last bit.

**T3.5 rescore loads no model** [unit, mac] monkeypatch `sys.modules["torch"]`/`vllm`/`transformers` import to raise; rescore must still succeed. Also assert `g_probe`/`serve` are not imported (`"experiments.phase2.g_probe" not in sys.modules` after the call in a subprocess).

**T3.6 manifest SHA guard** [integration, mac]
- Input: same answers, view whose `manifest.json` differs by one byte (a trailing space).
- Expected: exit code exactly 2 (`pytest.raises(SystemExit) as e; e.value.code == 2`, or subprocess returncode 2), and the out dir has NO new file (empty or not created). With matching SHA exit 0.
- Also: accepts a partial and a full `final_eval.json` in DIR; DIR without `final_eval.json`, or one lacking `view_manifest_sha256` -> nonzero exit, clear message (pin: not exit 2 silently producing output).
- Why: scoring against the wrong view gives plausible but wrong numbers.

**T3.7 `view_manifest_sha256` value** [unit, mac] equals `hashlib.sha256(manifest.read_bytes()).hexdigest()`; manifest with CRLF or unicode: hash over raw bytes, not parsed JSON.

**T3.8 `code_sha` precedence** [unit, mac]
- (a) `<root>/CODE_SHA` containing `"abc123-dirty-0123456789ab\n"` plus a `.git` dir -> `"abc123-dirty-0123456789ab"` (stripped, CODE_SHA wins).
- (b) no CODE_SHA, repo with `.git` (tmp git repo created by the test with `git init`, one commit) -> `git rev-parse HEAD` value.
- (c) neither -> `None` (JSON null). Empty CODE_SHA file: pin (None or empty? spec says "stripped content ... when exists": empty string).
- Why: on the box there is no `.git`; wrong precedence stamps wrong code identity.

**T3.9 `checker_build`** [unit, mac] on the Mac equals `"cpu"`; assert it is one of {"cuda","cpu"} and matches `objdel._gpucheck`'s own flag.

**T3.10 truncated answers score as given** [integration, mac] a sample with `truncated=True` and an otherwise legal text: spec says illegal/unparseable score 0; pin whether truncation forces illegal (see sample's path). Assert rescore and direct agree either way.

**T3.11 sampling determinism** [box] two runs with `--seed 11` on 16 items: `answers.jsonl` equal? (vLLM may not be bit-deterministic; assert only schema and that scoring of identical answers is identical.) Mark xfail-tolerant; the hard gate is T3.4.

**T3.12 box counts (gate 7)** [box, e2e] `answers.jsonl` has 1,084 (val) / 2,200 (test) lines, every list length 8; paired with `runs/r5/eval/eval_test/items.jsonl` by `item_id`: set equality of ids; mean per-item `iou_tc` difference within +-0.01 and the nodule CI of the difference contains 0.

---------------------------------------------------------------------------------------------------------------
## 4. Sampling `sample.sample` (spec 3.3)

**T4.1 tokens are completion id lengths** [unit, mac, fake trainer]
- Setup: fake `trainer._generate(prompts)` returning `(prompt_ids, completion_ids, ..., )` in the real signature order (read from `objdel.train.sample` for the tuple shape); completion ids lists of lengths [5, 1536, 7].
- Expected: `tokens == [5,1536,7]`; `len(tokens)==K==len(samples)==len(truncated)` per item; item order and `item_id` preserved; `items=16` chunking: 40 rows, K=8, items=16 -> 3 generate calls, outputs in input order.
- Why: tokens counted from decoded text (wrong) vs id sequence (right); padding counted; off-by-one with EOS. Use ids containing a pad id and assert the pad is not counted if `_generate` already strips (assert whatever the fake returns is the length).

**T4.2 truncated flag consistent with tokens** [property, mac] fake generate with `max_completion_length=1536`: `truncated[i] == (tokens[i] >= 1536 and no EOS)`; pin the rule from the code; never `tokens>1536`.

**T4.3 empty rows** [unit] `sample(trainer, [], 8)` returns `[]` without calling generate.

---------------------------------------------------------------------------------------------------------------
## 5. Boundary IoU (spec 3.4, gate 5)

**T5.1 5x5 square band** [unit, mac] `m` = 5x5 True block in a 20x20 frame at (7,7); `band(m,1).sum() == 16`; band is exactly the ring (all True except the center 3x3). `band(m,2).sum() == 24` (d=2: erosion twice leaves the center 1 pixel; 25-1). `band(m,3).sum()==25` (erode of 1 pixel -> empty). Why: erosion iterations semantics.

**T5.2 identical, disjoint, empty** [unit, mac] identical -> 1.0 exactly; disjoint squares -> 0.0; both empty -> 1.0; exactly one empty (a nonempty mask vs all-zeros) -> 0.0 (band of a nonempty mask is non-empty, so union >0, inter 0). Also two empty bands from NON-empty masks: none exist (a nonempty mask always has a nonempty band); assert `band(m,d).any()` for every nonempty m, d (property, random masks).

**T5.3 equals full-frame reference formula exactly on random masks** [property, mac]
- Reference in test: `ref_band(m,d)`: `p=cv2.copyMakeBorder(m.astype(u8),1,1,1,1,cv2.BORDER_CONSTANT,value=0)`; `e=cv2.erode(p, np.ones((3,3),u8), iterations=d)[1:-1,1:-1]`; band = `m & ~e.astype(bool)`; IoU from that.
- Input: 500 random masks 64x64 from (a) random blobs via `cv2.circle/ellipse/polylines`, (b) 30 % salt-and-pepper, (c) masks with holes, (d) masks touching/covering borders: full-frame True, top row only, single corner pixel, 1-pixel-wide full column; d in {1,2,9,16}, also d larger than the mask.
- Expected: `==` exactly for every pair.
- Why: the key bug is cv2.erode default border (BORDER_CONSTANT with `+inf`, i.e. erosion does NOT erode at the frame edge), which makes edge-touching masks lose their boundary band at the frame edge. The padding step must give the edge a band. Assert explicitly: full-frame True mask 8x8, d=1 -> band is the 28-pixel outer ring, not empty.

**T5.4 single pixel** [unit] one True pixel: band == itself for d>=1 (erosion removes it); boundary_iou(single, same) = 1.0; vs neighbour pixel = 0.0.

**T5.5 symmetry and range (property)** [property] `boundary_iou(a,b,d)==boundary_iou(b,a,d)`, in [0,1]; monotone: `band(m,d)` grows with d (superset), `band(m,d1) <= band(m,d2)` for d1<d2.

**T5.6 dtype/shape/input mutation** [unit] output bool, same shape; input not modified (`m.copy()` equal after); non-contiguous input (a transpose) and uint8 0/1 input behave as bool (or raise TypeError - pin).

**T5.7 frame invariance (property)** [property] embedding the same mask in a larger zero frame (shifted by (3,5)) gives the same boundary_iou. Why: pad logic is not frame-size dependent; a mask touching the edge within a small frame equals its centered version in a bigger frame.

**T5.8 metrics() additions** [integration, mac]
- Input: synthetic item (conftest disc), K=4: one exact copy of a rater, one shifted by 3 px, one illegal (self-intersecting polygon text), one unparseable string; 2 raters.
- Expected: `objects["nodule"]["boundary_iou"]` length 4; copy -> 1.0; illegal & unparseable -> exactly 0.0 for BOTH `boundary_iou` and `boundary_iou_vote`; shifted value equals max over raters of `boundary_iou(sample_mask, rater, 2)`, recomputed in the test with the reference; `_vote` is against `vote` mask of strict majority (`2*count > k`): with k=2 vote = intersection of both raters (assert!). Disc uses d=9, cup 16 (RIGA fixture or monkeypatched).
- `sample_iou` is 4x4 nested list, diagonal 1.0 for legal answers; symmetric; illegal/unparseable rows IoU 0 with legal-nonempty, and 1.0 with each other (both empty) - including diagonal entries of illegal = 1.0.
- Existing keys unchanged: compare against a stored golden of the old output (key set `old <= new`, values equal) - a regression golden is dumped from the current code BEFORE the change.

**T5.9 legal-but-empty answer** [unit] a legal answer rasterizing to an empty mask (e.g. `"[]"` if legal): boundary_iou 0.0 or 1.0? Spec: "illegal ... get 0.0"; a legal empty vs nonempty rater is "exactly one empty band" -> 0.0; assert 0.0 and `sample_iou` with another empty = 1.0.

**T5.10 speed** [perf, mac] `boundary_iou` on 512x512 masks, d=16: < 20 ms per call (generous); metrics() for K=8 items does not exceed old time by more than 25 %.

---------------------------------------------------------------------------------------------------------------
## 6. Selection `medoid`, `vote_mask` (spec 3.5, gate 4)

**T6.1 hand cases for medoid** [unit, mac]
- K=3, iou = [[1,.5,.5],[.5,1,.9],[.5,.9,1]], all legal -> scores [.5,.7,.7] tie -> lowest index 1.
- K=3, `legal=[0,0,0]` -> 0 regardless of matrix (even if matrix favours 2).
- K=3, `legal=[0,1,0]` -> 1. `legal=[0,0,1]` -> 2 (one eligible returns it).
- K=4, legal=[1,1,1,1], answers 0,1 empty (iou(0,1)=1) and 2,3 distinct nonempty with zero IoU to everything else: scores: 0:(1+0+0)/3=1/3, 1: 1/3, 2: 0, 3: 0 -> returns 0 (tie, lowest). K=5 with three legal empty answers (IoU 1 among them) and two nonempty -> medoid in the empty set, lowest of them. Spec: "a majority of legal empty answers wins".
- Diagonal ignored: set diagonal to 0 or 5; result unchanged. (Why: j != i.)
- Why: tie rule and diagonal.

**T6.2 argmax/tie property** [property] random symmetric matrices with entries on a coarse grid (multiples of 0.25 to force ties) - result equals a brute-force Python implementation (loop with `max(..., key=(score, -i))`) for 2000 cases, K in 1..8, random legal masks. Include float-tie hazard: scores computed as mean of floats like .1+.2 vs .3: the spec says ties go to lowest index; assert the implementation's tie decision equals the brute-force using the SAME arithmetic (`np.mean`); document exact-vs-approx. Adversarial: a case where two scores differ by 1e-17: the spec is silent, so pin current numpy behaviour.

**T6.3 permutation of answers equivariance** [property] permuting answers (and `legal`, and the matrix rows+cols) maps the medoid index through the permutation when there is no tie.

**T6.4 independence of rater masks (gate 4)** [property/integration] build an item with raters R; compute `per_item` selection keys; permute the raters, and replace rater masks with other random masks (same K answers): `medoid` index identical (it takes only `sample_iou`, `legal`). Assert `inspect.signature(medoid)` has exactly (`sample_iou`, `legal`) - no rater parameter. Also monkeypatch `ExpertCache`/rater access to raise inside the selection module call: medoid must not touch it.

**T6.5 input validation** [unit] non-square matrix, `len(legal) != K` -> ValueError (not silent truncation or broadcasting). K=0 -> pin (return 0 or raise).

**T6.6 vote_mask hand cases** [unit] masks (3,H,W): three legal -> pixel true iff >=2 true; four legal -> >=3 (2*2=4 > 4 false, so a 2-2 tie is False); one legal masks of the 4 illegal-flagged ones with arbitrary content ignored: legal=[1,1,0,0] with masks 0,1 equal -> equals mask 0; two legal masks differing -> intersection only (2*1 > 2 false); no legal -> all False, shape (H,W), dtype bool. K=1 legal -> that mask. Why: strict majority off-by-one (>= vs >).

**T6.7 vote_mask invariants** [property] vote of `[m, m, m]` is m; vote of `[m, ~m]` empty; result is subset of union and superset of intersection of legal masks; invariant under permutation; does not mutate input; accepts uint8? pin.

**T6.8 vote vs ExpertCache.vote consistency** [property] feeding the k rater stack as `masks` with all legal gives `ExpertCache.vote` exactly (`2*sum > k`). Why: two implementations of one rule.

**T6.9 `sel_vote` IoU_tc uses the item's cache** [integration] score a `vote_mask` through `score_mask(vote, cache).iou_tc`: when vote == one rater mask and raters identical, equals 1.0. If vote is empty (no legal) -> `sel_vote == 0.0`.

---------------------------------------------------------------------------------------------------------------
## 7. Fair ceiling `loo_caches`, `item`, `summarize` (spec 3.6, gate 3)

Setup helper: `cache_from(masks)` = `ExpertCache.build` of a stack with conftest `_disc` shapes (same call the repo uses; read tests/test_reward*.py for tau/t5 defaults).

**T7.1 loo structure k=4** [unit, mac] four distinct discs: `loo_caches` returns 4 `LeaveOut`, `index == [0,1,2,3]` (increasing), each `single is False`, each `cache.k == 3`, and `cache.raters` equals the stack with row r removed (array_equal). `isinstance(x, tuple)` NamedTuple fields (`index, cache, single`).

**T7.2 loo k=3 and k=2** [unit, mac] k=3 -> caches of 2 raters, `single False`; k=2 -> `single True`, `cache.k == 2`, raters `[other, other]`. For k=2 scoring any mask m against `loo.cache`: `score_mask(m, loo.cache).iou_tc == plain_iou(m, other)` with `np` IoU `(m&o).sum()/(m|o).sum()` - tested on 200 random legal masks (random discs shifted/scaled), to `==`/`1e-12`. Also `m == other` -> 1.0, disjoint -> 0.0. Why: the [other, other] trick relies on a band equal to the mask; if band widening applies, it breaks.

**T7.3 only voting raters** [unit] item with 4 raters where rater 3 is a non-voter (tiny or empty mask, per `cache.voters` rule): `loo_caches` has an entry for each of `cache.voters` only, but the stack removed `r` still contains the non-voter row (spec: "stack without r"). Assert the entry count == `len(voters)`, indices are voters.

**T7.4 doctor == duplicate scores 1.0 (gate 3)** [unit] k=3 where rater 0 equals rater 1 exactly: `doctor` contribution for r=0 vs `loo(0)` (raters 1,2): `score_mask(raters[0], loo0.cache).iou_tc` = 1.0 only if... careful: rater0 == rater1 and loo0 contains rater1, but also rater2. Per spec "a doctor equal to another doctor scores 1.0": construct k=2 identical masks: `item(...)["doctor"] == 1.0`; and k=4 all identical: `doctor == 1.0`.

**T7.5 all same signatures k=4, doctor equals legacy max IoU** [unit] four copies of the same disc at tiny jitter that keeps one topology signature: `doctor` (mean over r of IoU_tc vs LOO) equals `legacy_loo_dmax_item`? Gate wording: "doctor IoU_tc equals the legacy max IoU". Use exactly identical masks: both 1.0. And a hand case: k=2, masks A, B: `doctor == plain_iou(A,B)` (both entries equal) and `legacy_dmax == [plain_iou, plain_iou]`; `legacy_loo_dmax_item == plain_iou`.

**T7.6 legacy_dmax definition** [unit] k=3 hand-drawn rectangles with known IoUs: A-B .5, A-C .2, B-C .8 -> legacy_dmax over r=A,B,C = [.5,.8,.8]; `legacy_loo_dmax_item == 0.7`. Why: "max plain IoU with any other rater".

**T7.7 model = doctor-equal answer gives diff 0 (gate 3)** [integration] k=3, texts = [polygon text rasterizing to rater r] for K=1 repeated: `model` = mean over entries of IoU_tc(text r vs loo_e): equals `doctor` iff the text equals every rater, i.e. only for identical raters. Use identical raters (k=3 all same): `model == doctor`, `diff == 0.0` exactly. For distinct raters: with the answer == rater r only, `diff` is 0 for that entry's term (assert via recompute using entry-level scoring, not the mean).

**T7.8 illegal / failed answer scores 0.0** [unit] texts = [valid, self-intersecting, "not json"]; `model` = mean of (x, 0, 0) per entry. `truncated=[False, False, True]` with a text that is legal: spec only says "failed"; pin whether truncated forces 0 (follow the code of `score_completion`). `model_all_k` uses full `cache` with the same zeros.

**T7.9 without texts** [unit] `texts=None, truncated=None`: `model`, `model_all_k`, `diff` are `None` (key present); `doctor`, `legacy_*`, `k`, `single` present. `texts=[]`: pin (None or error).

**T7.10 `single` flag and k** [unit] k=2 -> `single True`, k=3,4 -> False. A cache with k=2 raters but one non-voter: voters has 1 -> pin: `loo_caches` returns 1 entry; k==2 so `single` still True? Spec: `single` = (k==2). Assert via `item()`.

**T7.11 k=1** [edge] single rater: `voters` has 1 voter; loo with no other rater. Spec covers 2..4; assert a clear ValueError (not an IndexError or `ExpertCache.build` of an empty stack crash with an obscure message). LIDC has k>=2 so this is a guard test.

**T7.12 summarize strata** [unit, mac] 10 rows: 6 with k=3/4, 4 with k=2 plus hand values; `n all = 10 = 6 + 4`; each metric `[mean, lo, hi]` with lo<=mean<=hi, `mean == np.mean`; `diff` mean equals mean of per-item `diff` (not `model - doctor` of means, which is equal for means but the CI differs: test by recomputing with the same draws); stratum `k=2` excludes k>=3 rows; stratum with zero rows: pin (n 0, entries None) - must not raise on `k>=3` empty (e.g. only k=2 items on a tiny fixture).
- `legacy_loo_dmax_pooled` = mean of all legacy_dmax values pooled (weights by number of voters): hand case rows with legacy_dmax `[.5,.5]` and `[.9,.9,.9,.9]` -> pooled `(1+3.6)/6 = .7667` while `legacy_loo_dmax_item` mean = `(.5+.9)/2 = .7`. When every row has the same k, pooled == mean of item values (1e-12).
- Rows with `model=None` (no texts): `model`, `model_all_k`, `diff` strata entries None, rather than crash.
- With `groups`: clusters for CIs by nodule via `cluster()`; assert draws equals `signal.draws(ids)` (rows must carry `item_id`? the spec's row has no `item_id` - rows come with it from the CLI; assert KeyError-free behaviour with the field the code uses and document).

**T7.13 target text rasterizes to its rater mask** [integration, mac] for 30 synthetic rater masks (discs, ellipses, a C shape with a hole, a thin diagonal): build text through the repo's target writer, `parse` + `to_masks` -> `np.array_equal` with the rater mask. Include touching-frame masks (border pixel) and a mask with a hole (if the format supports holes) - for those where the format loses information pin the known loss, don't skip silently.

**T7.14 adversarial: LOO is label-fair** [property, mac] A random answer set scored vs LOO caches never exceeds IoU_tc 1.0 and a doctor stack is scored on exactly the same caches as the model: monkeypatch `score_mask` to record `id(cache)`; for each entry r, `doctor` and `model` calls both use `loo.cache` of that r (set equality of ids). Why: the spec's entire point is the same k - 1 doctors.

**T7.15 box data gate (gate 3)** [box] on val: `n(all)=1084`, `n(k=2)=281`; `legacy_loo_dmax_item` or `_pooled` within +-0.005 of 0.734 (val) / 0.732 (test): assert `min(abs(item-0.734), abs(pooled-0.734)) <= 0.005`; fail otherwise (spec says either).

**T7.16 CPU vs GPU bit-equality of ceiling** [GPU, box] `item()` on 50 val items with both builds: `doctor`, `model` equal to `==` (the checker is exact, one source).

**T7.17 ceiling CLI** [integration, mac] `--answers FILE` optional: without it `model*`, `diff` null for all rows; with it uses answers by `item_id`; missing item in answers -> skip with count or raise (pin); outputs `ceiling.json` (has `seconds`) and `ceiling_items.jsonl`.

---------------------------------------------------------------------------------------------------------------
## 8. NVTX (spec 3.7)

Use a fake `torch.cuda` via `monkeypatch`; record calls in a list.

**T8.1 import is inert** [unit, mac] in a subprocess with `OBJDEL_NVTX=1`: replace `torch.cuda.is_available`/`nvtx.range_push` with functions that raise before `import objdel.train.nvtx`; import succeeds. (Spec: reads env, calls no torch.cuda.)

**T8.2 `enabled()` truth table and caching** [unit] env unset/"0"/"true"/"1" with cuda True/False: only ("1", True) -> True. `torch.cuda.is_available` call count after 5 `enabled()` calls == 1 (<=1 per process; use `importlib.reload` or a reset hook; with env != "1" the CUDA query need not run at all - assert <=1). Why: `is_available` is slow and can init CUDA; env "1" with cuda False -> False.

**T8.3 range pushes/pops, also on exception** [unit] enabled: `with range("x"): pass` -> calls `[("push","x"),("pop",)]`; `with pytest.raises(ValueError): with range("x"): raise ValueError` -> pop still recorded. Nested ranges -> LIFO order. Disabled: zero calls to any nvtx function (`range`, `push`, `pop`); `push/pop` mirrored.

**T8.4 `timed` still logs and opens range** [unit, mac] with a stub FastGRPO-like object (use the real `FastGRPO.timed` unbound with a fake self; read `fast.py` to build it): enabled -> range push "generate"; always -> `_metrics`/log gets `time/generate_s` appended once, float >= 0; exception inside the phase: time still recorded (pin current behaviour) and range popped.

**T8.5 range names** [unit] every name in the fixed list is used somewhere: grep-style test: `{step, loader, micro, batch, generate, sync, wake, sleep, rewards, checker, logps, forward, optimizer, log, vision_cache, rope_cache, image_memo}` all appear as string literals in a `nvtx.range(...)`/`timed(...)`/`push(...)` call under `src/objdel/train` (AST scan); and no call uses a name outside this list. Why: typos make nsys analysis silently miss ranges.

**T8.6 `training_step`, `get_batch_samples`, `log` ranges, results unchanged** [integration, mac, tiny offline TRL model from conftest] enabled with fake nvtx: a call to each records `micro`, `loader`, `log` ranges, balanced; return values equal to those with NVTX disabled (same seed, tensors `torch.equal`).

**T8.7 `NvtxSteps` balance** [unit] drive callbacks manually with a fake `state`: `on_train_begin`, N x `on_step_end` (global_step 1..N), `on_train_end`: calls = push step, then (pop, push) x N, final pop; `count(push)==count(pop)==N+1`; stack depth never negative or >1 at the `step` level. `on_train_end` without `on_train_begin` (crash recovery): no unbalanced pop (pin: no pop recorded / no exception).

**T8.8 optimizer range** [unit] `on_pre_optimizer_step` pushes `optimizer`, `on_optimizer_step` pops it; order relative to `step` range nested (push step, push optimizer, pop optimizer, pop step). Missing `on_optimizer_step` (e.g. grad overflow skipping) - pin.

**T8.9 profiler window** [unit] window=(4,6), CUDA available (faked), steps 1..8: `profiler.start` once, fired at `on_step_end` with `global_step == 3` AFTER the pop and BEFORE the next push (assert the order of recorded events: `pop, start, push`); `profiler.stop` once at `global_step == 6` after the pop (`pop, stop, push` - the push after the stop still occurs since training continues). Never start twice if `on_step_end` is called twice with the same global_step (idempotence - gradient accumulation can repeat; pin). `window=None`: no profiler calls. `window=(1,2)`: start at global_step 0, i.e. `on_step_end` with step 0 never happens -> start never fires; assert the spec: `s-1 == 0` - decide: start also in `on_train_begin`? Spec says only on_step_end, so test documents "window start < 2 never starts" and `make_trainer` guard (wait = s-2 < 0?) is checked.
- No CUDA: window set, zero profiler calls, push/pop of `step` also zero (`nvtx.enabled()` False).

**T8.10 `make_trainer` wiring** [unit, mac, with monkeypatched trainer class] `a.profile_window=None` -> callbacks contain no `NvtxSteps`, no `TraceStep` unless `a.trace`; `profile_window=[4,6]` -> exactly one `NvtxSteps` with window `(4,6)` (tuple); `trace="DIR"` with window -> `TraceStep(wait=2)`; with `trace` and no window -> `wait=1`; `trace=None` with window -> no TraceStep. Window `[1,3]` with trace -> `wait=-1`: pin (clamp to 0 or raise).

**T8.11 train_rl flags** [unit] `train_rl` argparse: `--profile-window 4 6` -> `[4,6]`; without flags the namespace attributes are None; mutually bad input (`--profile-window 6 4`) -> error exit (pin: must reject s>=e).

**T8.12 `raw_ratio` + `bounded` == `is_ratio` all four modes** [property, mac] torch random `logps, sampling` (B=6, T=40) with `mask` having random ragged lengths including an all-masked row and a fully-valid row; lo/hi in {None, (0.1, 3.0)}; for each mode in `sequence_mask, sequence_truncate, token_mask, token_truncate`: `torch.equal(bounded(raw_ratio(...), mode, lo, hi), is_ratio(...))` (exact, same ops). Adversarial: values giving ratios exactly at `lo`/`hi`; `+-inf` log-differences (`sampling = -inf`), NaN absent. Why: guards the refactor splitting one function into two.

**T8.13 share counters equal numpy counts** [unit, mac] toy `raw` ratios per sequence `[0.0, 0.05, 0.5, 1.0, 2.9, 3.0, 3.1, 100.0]`, `clip_min=0.1`, `clip_max=3.0`, sequence_mask mode (final ratio 0 where outside bounds): expected: `zero_share` = count(final==0)/8 ; `hi_share` = count(raw > 3.0)/8 = 2/8 (3.0 NOT above, strict); `lo_share` = count(raw < 0.1 or raw == 0)/8 = 2/8 (0.0 and 0.05; the 0.0 counted once); assert numpy recomputation equality and exact values `zero_share = 4/8` in sequence_mask when bounds mask-out both sides (0.0, 0.05, 3.1, 100.0), `sequence_truncate`: zero_share = 1/8 (only raw 0.0), hi_share 2/8, lo_share 2/8. Token modes: counted entries = masked tokens only (padding excluded): row with mask zeros must not enter the denominator - test with a toy where padding has ratio 0 (would inflate zero_share).
- `clip_min=None/clip_max=None`: `hi_share`/`lo_share` pin (0 or absent), lo_share still counts raw==0? Pin.
- No host read: run the computation under `torch.cuda.set_sync_debug_mode("error")` [GPU] or assert on CPU that returned values are tensors (not python floats) until logged; on Mac: monkeypatch `torch.Tensor.item`/`__float__`/`tolist` to raise during the counter computation, then allow in logging.
- Why: strictness of the inequalities, denominators.

**T8.14 `time/step_s`** [unit] call `log` three times with a fake clock (monkeypatch `time.perf_counter` or whatever the code uses) at t=0, 5, 12 -> first call has no `time/step_s` (or 0? pin), second 5.0, third 7.0; no `torch.cuda.synchronize` call (monkeypatch to raise).

**T8.15 NVTX overhead when disabled** [perf, mac] 100k `with range("x")` disabled < 0.3 s total (<3 us each).

**T8.16 CUDA, nsys smoke** [box GPU] run a 3-step toy with `OBJDEL_NVTX=1`: the exported ranges contain every name that fires in a normal step; `step` instances == steps.

---------------------------------------------------------------------------------------------------------------
## 9. Trace analysis `profile.covered`, `nsys_baseline` (spec 3.7a, gate 9)

**T9.1 `covered` hand cases** [unit, mac] intervals `[(0,2),(1,3),(5,6)]`, t0=0,t1=10 -> 4.0; clipped `t0=1.5,t1=5.5` -> (1.5..3)=1.5 + (5..5.5)=0.5 = 2.0; empty list -> 0.0; interval fully outside -> 0; touching `(0,1),(1,2)` -> 2.0; contained duplicates; reversed interval `(3,1)` -> pin (ignore or error); zero-length -> 0; unsorted input; generator input (iterator, consumed once); t1 <= t0 -> 0.0 (not negative).

**T9.2 `covered` property** [property] random intervals: `covered(...)` equals a brute-force count on an integer grid (intervals on integers, measure = number of unit cells), `0 <= covered <= t1 - t0`, `covered(A+B) <= covered(A)+covered(B)`, monotone in t1, idempotent under duplication. Never exceeds sum of lengths.

**T9.3 `busy` unchanged** [regression] stored cases: `profile.busy` on the existing tests' inputs unchanged (the existing tests/test_profile.py must pass untouched); add a golden of 5 hand inputs.

**T9.4 nsys_baseline CSV fixture, hand values** [integration, mac] create tiny CSVs in tmp: `stats_nvtx_pushpop_trace.csv` (columns as nsys exports; read the real column names from `experiments/phase2/profile.py` and `nsys_baseline` loader; if absent, fixture header = nsys `nvtxpushpop_trace`: `Start (ns),End (ns),Duration (ns),Name,Level?,TID...` - the test author copies the real header from a real nsys export checked in `tests/goldens/` first).
- Fixture: thread 1: `step` [0,100] ; inside: `generate` [0,40], `sync` [40,55], `rewards` [60,90]; a second `step` [100,200] with `forward` [100,150] only. A range on thread 2 `generate` [0,100] must NOT count for thread 1's step (same-thread containment). A nested range `batch` inside `generate` must not count (only directly inside).
- Expected step 0: `coverage = (40+15+30)/100 = 0.85`, step 1: `0.5`.
- GPU trace: kernels [10,30],[20,50] (overlap), memcpy [60,70], memset [95,120] (straddles step end): `gpu_busy` step 0: union [10,50] + [60,70] + [95,100] = 40+10+5 = 55 -> 0.55; `gpu_idle = 0.45` exactly 1 - busy. Step 1: [100,120] = 0.2.
- API trace: calls inside step 0: `cudaStreamSynchronize_v3020` x2, `cudaDeviceSynchronize` x1, `cudaMemcpyAsync` (excluded), `cudaMemcpy_v3020` (included, name `cudaMemcpy`), `cudaMemcpyToSymbol` (non-Async memcpy-ish: included by "non-Async cudaMemcpy*" rule), `cudaMemcpy2DAsync` (excluded), `cudaEventSynchronize` x1, `cudaLaunchKernel` (excluded), a call straddling step boundary (decide: counted by start time inside the step; pin) -> `host_syncs == {"cudaStreamSynchronize":2,"cudaDeviceSynchronize":1,"cudaEventSynchronize":1,"cudaMemcpy":1,"cudaMemcpyToSymbol":1}`; suffix removal regex `_v\d+$` only at the end (a name like `cudaMemcpy_v2_foo`: not touched).
- `per_range_s`: total duration by range name in seconds (ns / 1e9), including nested ranges; `gpu_proj_s` from `stats_nvtx_gpu_proj_trace.csv`: by range name.
- `overhead`: log JSON with `time/step_s` per step: steps 1..8 `[10,10,10,10,12,12,10,10]` window (5,6) -> window mean 12, neighbours mean 10 -> `overhead_affected True` (20 % > 10 %); window equal to 10.9 -> False; exactly 11.0 (10.0 % slower): False (">" strict; mind float: 11.0/10.0-1 = 0.10000000000000009 > 0.1 True! pin the arithmetic - compare `win > 1.1*base` -> 11.0 > 11.000000000000002 False; the test pins whichever implementation, flagged as a float hazard).
- Why: nesting by time containment is easy to get wrong (siblings with identical start/end, ranges ending exactly at the step end).

**T9.5 coverage edge cases** [unit] step with no children -> 0.0; a child ending exactly at step end and one starting exactly at step start count fully; a child that overlaps the step boundary is clipped or excluded (pin); two children overlapping each other: sum of durations > step; coverage can exceed 1? Spec says summed duration / step duration, so >1 possible; assert no clamp (pin) or flagged. Step of zero duration -> no ZeroDivisionError (skip or 0.0).

**T9.6 missing/empty CSVs** [failure] missing file -> nonzero exit with the file name in the message; header-only `stats_cuda_gpu_trace.csv` -> `gpu_busy 0.0`, `gpu_idle 1.0`, no crash; no `step` ranges -> `[]` entries and clear message; units: nsys CSV has ns; assert conversion (a step of 1e9 ns is 1.0 s).

**T9.7 box (gate 9)** [box] real export: steps 4 and 5 each `coverage >= 0.90`.

---------------------------------------------------------------------------------------------------------------
## 10. TRL contract `trl_adapter.check` (spec 3.8, gate 10)

**T10.1 passes on pinned trl** [unit, mac] `check()` returns None on trl 1.17.1 with no GPU/vLLM (assert `"vllm" not in sys.modules` afterwards in a subprocess).

**T10.2 each failure names the item** [unit, mac] monkeypatch (restore after): one case per row, each `pytest.raises(TRLContractError, match=<name>)`:
- delete `GRPOTrainer._prepare_inputs` -> match `_prepare_inputs`;
- replace `_compute_loss` with `def _compute_loss(self, model, inputs, extra)` -> match `_compute_loss`; also a renamed parameter (same arity) and reordered parameters (exact order matters);
- replace `_get_per_token_logps_and_entropies` dropping `sampling_mask` -> match it;
- remove `trl.trainer.grpo_trainer.shuffle_sequence_dict` -> match `shuffle_sequence_dict`; change its param to `d`;
- remove `VLLMGeneration._iter_named_params`;
- `GRPOTrainer.__init__` source without `self._metrics` (patch `inspect.getsource` to return altered text) -> `_metrics`; `VLLMGeneration` source without `self._dist`;
- for each of the ten `self.<name> =` names, source with that name removed -> match the name;
- `VLLMGeneration` without `def sync_weights` / `self.llm` / `self.model`;
- `_generate_and_score_completions` source without `"tool_mask"` -> match `tool_mask`; without `self.use_vllm and self.vllm_importance_sampling_correction`; without each of the 3 other keys; `_compute_loss` source without `"importance_sampling_ratio"`.
- Why: `check` is the only guard against silent TRL drift; each substring must be live, not vacuous (a typo in the checked substring would make the check pass on any source: the passing test alone cannot catch that, the negative test does).

**T10.3 order of checks** [unit] break both item 1 (missing `_prepare_inputs`) and item 5 (missing `"tool_mask"`): the message names `_prepare_inputs` (first failure). Break item 2 and 4: names item 2's.

**T10.4 version mismatch** [unit] `trl.__version__ = "1.18.0"` with items 1-5 intact: raises (msg contains `1.18.0` and `1.17.1`) when `OBJDEL_TRL_UNPINNED` unset, `"0"`, `"true"`; with `"1"`: returns None and `pytest.warns(UserWarning)` exactly ONE warning that names `1.18.0` (`len(w)==1`). Mismatch AND item failure with UNPINNED=1: still raises (item error wins). Version check last: with a changed version and a missing method, the message names the method.

**T10.5 ratio block equals `fast.is_ratio`** [property, mac] extract TRL's block from `inspect.getsource(GRPOTrainer._generate_and_score_completions)` between the marker line `per_token_logps_diff = (old_per_token_logps - sampling_per_token_logps) * mask` (inclusive of that line) and the `else:` before `Unknown vLLM importance sampling level` (exclusive); `exec` it in a namespace with torch tensors and a stub `self` (attrs: `vllm_importance_sampling_mode`, `_clip_min/_clip_max` names per source), 4 modes x 20 random seeds incl. all-masked rows; assert `torch.equal` with `fast.is_ratio(old_logps, sampling_logps, mask, mode, lo, hi)`. Test fails (not skips) if the markers are not found. (Indentation: `textwrap.dedent` after slicing.)
- Why: our patched ratio could differ from TRL on edge rows (all-masked sequence: mean of empty).

**T10.6 `FastGRPO.__init__` calls check first** [unit] monkeypatch `trl_adapter.check` to record and raise `TRLContractError`, patch `GRPOTrainer.__init__` to record: constructing raises the contract error and `GRPOTrainer.__init__` was never called (and no `torch` cuda touch).

**T10.7 preflight entry** [unit, mac] `preflight.CHECKS` contains `("CODE TRL private names", fn)` (the structure per existing tests/test_preflight.py); it is the first entry whose name starts with `CODE` (index of first "CODE" prefixed name); `fn(a, cfg) == []` when passing; when check raises `TRLContractError("boom: _compute_loss")`: returns `["boom: _compute_loss"]` (a one-element list of `str`). A non-TRLContractError exception (e.g. ImportError) propagates or is wrapped? pin.

---------------------------------------------------------------------------------------------------------------
## 11. Outline lengths (spec 3.9, gate 11)

**T11.1 count with a char stub** [unit, mac] `tokens_of=len`; records: A targets `["abcd", None, "ab"]`, B targets `["abcdefghij"]`; Expected: wait: the +1 end token is applied in the CLI's `tokens_of`, not in `count`. `count` takes `tokens_of` as is: lengths [4,2,10]; `n_records==2`, `n_targets==3` (None skipped), mean 5.333.., `max 10`, `above[256]==0`. Test with `tokens_of=lambda s: len(s)+1` for the +1: [5,3,11].
- Percentiles: with lengths 1..1000 (a record each) p50 == `np.percentile` (linear) 500.5, p90 900.1, p99 990.01, p99_9 999.001 exactly equal numpy; pin method (linear).
- `above = {256:..., 512:..., 800:..., 1000:..., 1200:..., 1536:...}` keys ints; strict `>`: a length of exactly 1200 is not above 1200; exactly 1201 is. Test with lengths [1200,1201,1536,1537].

**T11.2 voting subset** [unit] `voting=lambda r: {0,2}` over records with 4 targets each: `n_voting_targets` counts only indices {0,2} and non-None; `above_voting` same; a voting index whose target is None: skipped; `voting=None` -> `n_voting_targets == n_targets`? spec: `voting` optional - pin (equal or absent).

**T11.3 longest** [unit] 30 targets -> `longest` has exactly 20 entries, strictly sorted by tokens descending, ties stable (pin), each `{item_id, rater, tokens}` with `rater` = target index; fewer than 20 targets -> all.

**T11.4 empty input** [edge] no records -> zeros/None without crash (`np.percentile([])` raises): pin `n_records == 0` and statistics null; all-None targets same.

**T11.5 `answers` hand fixture** [unit] answer rows: item i1 tokens `[100, 300, 600, 1600]`, item i2 `[50, 1537]`; metric rows with per-sample `legal`, `iou_tc` (find structure `metrics` output: `samples.legal/iou_tc`); join by `item_id`; for each L: `share_above[L]` = share with tokens > L over ALL answers (6): L=256 -> 4/6, L=512 -> 3/6, L=800 -> 2/6, L=1536 -> 2/6, L=1536 > 1536: 1600 and 1537 -> 2/6; `legal_share`, `mean_iou_tc` over answers above L; `n_legal_good` = legal and iou_tc >= 0.5 (0.5 exactly counts). L with no answers above: `legal_share` None, `mean_iou_tc` None (not ZeroDivision), `n_legal_good 0`. An answer row with no matching metric row -> raise or skip (pin); duplicate item_id -> pin. Different K per item ok.

**T11.6 refactor names** [unit] `objdel.train.dataset.target_tokens` and `.voting` exist; `_tokens` and `_voting` do NOT (no aliases): `not hasattr`. Search repo for the old names (`grep -rn "dataset._tokens\|dataset._voting"`) in `src`, `experiments`, `tests`: no match. Behaviour of `target_tokens` equals the old `_tokens` on a fixture (golden captured before rename).

**T11.7 CLI tokenizer formula** [unit] stub tokenizer `lambda text, add_special_tokens: {"input_ids": list(text)}` -> `tokens_of(text) == len(text)+1`; `add_special_tokens=False` is passed (spy).

**T11.8 box data (gate 11)** [box] `data/prepared` LIDC train: `n_records == 7341`, `n_targets == 23469`, `max == 1182`, `above[1200] == 0`; outline_lengths.json includes the rows for `lidc-train-shifted` and `lidc-train-rescaled` with `manifest_sha256` matching the files' SHA-256; keyed by root, dataset, split; `model_answers` section per `--answers` file.

---------------------------------------------------------------------------------------------------------------
## 12. Code identity on the box (spec 3.2a)

**T12.1 CODE_SHA format** [unit, mac, tmp git repo; monkeypatch the ssh/rsync part] clean tree: file content `<40 hex>` with no `-dirty`; modified tracked file: `<HEAD>-dirty-<12 hex>` where the 12 hex are `sha256(git diff HEAD output).hexdigest()[:12]`, recomputed in the test; a modification then revert -> clean again; staged-only change counts as dirty (`git diff HEAD` includes staged); an untracked file: `git diff HEAD` ignores untracked - pin as clean (documented gap: untracked files are not in the stamp; flag as a design risk to the owner). Written for `--code-only` too, path `/workspace/src/CODE_SHA` (assert the remote command string/rsync target). Never ends with extra bytes beyond newline (consumer strips).

**T12.2 no git** [failure] `git` failing -> sync_up raises or writes nothing, never writes an empty/garbage CODE_SHA.

---------------------------------------------------------------------------------------------------------------
## 13. Cross-cutting, adversarial, and performance gates

**X1 no pixels leave the device** [GPU, box] with `torch.cuda.set_sync_debug_mode` unavailable for custom code: assert the python binding's returns are scalars (`objdel._gpucheck` result arrays have size O(items), not O(pixels)); a memory-transfer spy (nsys `cuda_gpu_mem_time_sum` / cudaMemcpy DtoH total bytes during `score_mask` for 64 items < 64 * 1 KiB). The new ceiling code (`loo_caches` scoring via `loo.cache.device()`) must use device tables: check `cache.device()` is called on loo caches and that DtoH bytes per scoring call stay below threshold.

**X2 CPU vs GPU bit-equality of all new metrics** [GPU, box] `final_eval --rescore` of the same answers with cuda build vs cpu build: per-sample `legal, iou_tc, f5, g_area` `==` for 1,084 val items x 8 (gate 6 box line). `boundary_iou` is CPU-side (cv2): equal by construction; assert computed value independent of build.

**X3 deterministic reruns** [integration, mac] run `summarize` twice, and with items shuffled: `mean` equal to 1e-12 (float order), CIs differ after shuffle because draws depend on order - assert and document (shuffle changes ids order -> different clusters order); sorting items by id before summarize is advisable; test that `final_eval` sorts or keeps answer order and that `rescore` preserves `answers.jsonl` order, so `items.jsonl` order equals answers order.

**X4 large input** [perf, mac] `summarize` on 2,200 items, 12 metrics + 6 selection keys, draws 2000: < 20 s on the Mac (3 CI kinds); memory: `idx` of shape (2000, 2200) int64 = 35 MB - OK; `signal.draws` for the 2,200 golden: < 3 s. Gate: no per-metric recomputation of draws (compute once; spy `signal.draws` call count == 2 (nodule, patient), not 2 x 18).

**X5 max-size frame** [perf, mac] boundary_iou/band for 512x512 masks and `sample_iou` K=8: metrics() for one item with K=8, 4 raters < 2 s on CPU build (pin from today's timing + 25 %).

**X6 JSON-native outputs** [unit] `json.dumps(summarize(...))` and `json.dumps(per_item(...))`, `ceiling.summarize` succeed with no `numpy` scalars (np.float64 is a float subclass and dumps OK; np.int64 and np.bool_ fail): assert `json.dumps` passes and `type(x) in (float,int,bool,str,list,dict,NoneType)` recursively.

**X7 no new dependency** [unit] AST-scan new modules' imports `⊆ {stdlib, numpy, scipy, cv2, torch, trl, transformers, objdel.*, experiments.*, PIL}`.

**X8 reads untrusted data safely** [unit] CLI modules run via `python -I -m ...`: `subprocess.run([sys.executable, "-I", "-m", mod, "--help"])` returns 0 for each new CLI (`final_eval`, `ceiling`, `nsys_baseline`, `outline_lengths`) with a planted `json.py` in cwd that raises on import (must not run). 

**X9 make lint/test** [e2e] `make test` and `make lint` pass; new tests count > 0; no test reads `runs/` or `data/` unless marked [box] (fixture-scan: grep tests for `"runs/"`/`"data/"` outside marked tests).

**X10 golden regression of today's numbers** [regression, mac] tests/goldens `r5_test_items.jsonl.gz` has 2200 rows with `item_id`, `iou_tc`; assert 2,200 rows, 401 clusters (`len({cluster(i)})`), unique ids, `0<=iou<=1`; whole Gate 1 values (T1.9). Also the old `final_eval.summarize` output on a fixture is a subset of the new (`mean` unchanged, old `ci95` == new `ci95_item`).

---------------------------------------------------------------------------------------------------------------
## Pytest sketches for the 5 most important tests

Ranking by bug-catching power on the gates: (1) golden CI regression + exactness of draws, (2) boundary IoU vs full-frame reference incl. edge, (3) rescore bit-equality + SHA guard, (4) medoid hand cases + rater independence, (5) fair ceiling k=2 plain IoU and LOO fairness.

```python
# tests/test_a_signal.py  (T1.1 + T1.9)
import gzip, json
from pathlib import Path
import numpy as np
from objdel.eval import signal as S
from experiments.phase2 import final_eval as FE

GOLD = Path(__file__).parent / "goldens" / "r5_test_items.jsonl.gz"

def load():
    with gzip.open(GOLD, "rt") as f:
        return [json.loads(l) for l in f]

def legacy_draws(ids, boot, seed):                      # verbatim copy of today's body
    nod = [S.cluster(i) for i in ids]
    groups = [np.array([j for j, x in enumerate(nod) if x == u]) for u in dict.fromkeys(nod)]
    rng = np.random.default_rng(seed)
    return [np.concatenate([groups[g] for g in rng.integers(0, len(groups), len(groups))]) for _ in range(boot)]

def test_draws_default_key_is_bit_identical_to_legacy():
    ids = [r["item_id"] for r in load()]
    for boot, seed in [(50, 0), (200, 3)]:
        new, old = S.draws(ids, boot, seed), legacy_draws(ids, boot, seed)
        assert len(new) == len(old) == boot
        assert all(np.array_equal(a, b) and a.dtype == b.dtype for a, b in zip(new, old))

def test_one_item_per_cluster_equals_integers_rows():
    ids = [f"x{i}" for i in range(10)]
    got = np.stack(S.draws(ids, 100, 5))
    assert np.array_equal(got, np.random.default_rng(5).integers(0, 10, (100, 10)))

def test_gate1_exact_numbers_on_r5_test_items():
    rows = load()
    items = [{"item_id": r["item_id"], "iou_tc": r["iou_tc"]} for r in rows]   # plus the other keys summarize needs, set to iou_tc
    for it in items:
        for k in FE.METRICS: it.setdefault(k, it["iou_tc"])
        it["miss"] = 1.0 if it["iou_tc"] < 0.1 else 0.0
    s = FE.summarize(items, draws=2000, seed=0)
    assert s["iou_tc"]["ci95_item"] == [0.6334010210209712, 0.6561989221793691]
    assert s["iou_tc"]["ci95_nodule"] == [0.6279017988891661, 0.6620867640136667]
    assert s["miss_rate"]["mean"] == 0.10818181818181818
    assert (s["items"], s["nodules"]) == (2200, 401)
```

```python
# tests/test_a_metrics.py  (T5.1 + T5.3)
import cv2, numpy as np, pytest
from objdel.eval import metrics as M

def ref_band(m, d):
    p = cv2.copyMakeBorder(m.astype(np.uint8), 1, 1, 1, 1, cv2.BORDER_CONSTANT, value=0)
    e = cv2.erode(p, np.ones((3, 3), np.uint8), iterations=d)[1:-1, 1:-1].astype(bool)
    return m & ~e

def ref_biou(a, b, d):
    A, B = ref_band(a, d), ref_band(b, d)
    u = (A | B).sum()
    return 1.0 if u == 0 else (A & B).sum() / u

def test_five_by_five_square_has_sixteen_band_pixels():
    m = np.zeros((20, 20), bool); m[7:12, 7:12] = True
    assert M.band(m, 1).sum() == 16 and not M.band(m, 1)[8:11, 8:11].any()

def test_full_frame_mask_keeps_its_outer_ring():           # cv2's default border would erase it
    m = np.ones((8, 8), bool)
    assert M.band(m, 1).sum() == 28

def rand_mask(rng, kind, n=64):
    m = np.zeros((n, n), np.uint8)
    if kind == 0:   cv2.circle(m, tuple(rng.integers(0, n, 2)), int(rng.integers(1, 30)), 1, -1)
    elif kind == 1: m = (rng.random((n, n)) < 0.3).astype(np.uint8)
    elif kind == 2: m[: rng.integers(1, n), rng.integers(0, n)] = 1                      # touches the edge
    else:           m[:] = 1; cv2.circle(m, (n // 2, n // 2), int(rng.integers(1, 20)), 0, -1)  # hole
    return m.astype(bool)

@pytest.mark.parametrize("d", [1, 2, 9, 16])
def test_boundary_iou_equals_full_frame_reference_exactly(d):
    rng = np.random.default_rng(d)
    for i in range(125):
        a, b = rand_mask(rng, i % 4), rand_mask(rng, (i + 1) % 4)
        assert M.boundary_iou(a, b, d) == ref_biou(a, b, d)
        assert np.array_equal(M.band(a, d), ref_band(a, d))

def test_degenerate_pairs():
    z = np.zeros((10, 10), bool); s = z.copy(); s[2:6, 2:6] = True; t = z.copy(); t[7:9, 7:9] = True
    assert M.boundary_iou(z, z, 2) == 1.0 and M.boundary_iou(z, s, 2) == 0.0
    assert M.boundary_iou(s, s, 2) == 1.0 and M.boundary_iou(s, t, 1) == 0.0
```

```python
# tests/test_a_final_eval.py  (T3.1 + T3.4 + T3.6)  -- fixture helper `make_view(tmp_path)` builds 4 synthetic items
import hashlib, json, subprocess, sys
import pytest
from experiments.phase2 import final_eval as FE

def run_rescore(src, view, out):
    return subprocess.run([sys.executable, "-I", "-m", "experiments.phase2.final_eval", "--rescore", str(src),
                           "--view", str(view), "--split", "val", "--out", str(out)], capture_output=True, text=True)

def test_scoring_failure_leaves_answers_and_partial_summary(tmp_path, monkeypatch, fake_sampler, make_view):
    view = make_view(tmp_path)                                  # fake_sampler patches the [box] sampling layer
    monkeypatch.setattr(FE, "score_items", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")))
    with pytest.raises(RuntimeError, match="boom"):
        FE.main(["--ckpt", "c", "--view", str(view), "--split", "val", "--train", "t", "--out", str(tmp_path / "o")])
    out = tmp_path / "o"
    lines = [json.loads(l) for l in (out / "answers.jsonl").read_text().splitlines()]
    assert all(len(r["samples"]) == len(r["truncated"]) == len(r["tokens"]) == 8 for r in lines)
    doc = json.loads((out / "final_eval.json").read_text())
    assert doc["partial"] is True and doc["view_manifest_sha256"] == hashlib.sha256((view / "manifest.json").read_bytes()).hexdigest()
    assert not (out / "metrics.jsonl").exists() and not (out / "items.jsonl").exists()

def test_rescore_reproduces_direct_scoring_to_the_bit(tmp_path, direct_run, make_view):
    view, src = direct_run                                       # a finished run: answers + metrics.jsonl
    r = run_rescore(src, view, tmp_path / "re")
    assert r.returncode == 0, r.stderr
    old = [json.loads(l) for l in (src / "metrics.jsonl").read_text().splitlines()]
    new = [json.loads(l) for l in (tmp_path / "re" / "metrics.jsonl").read_text().splitlines()]
    for a, b in zip(old, new, strict=True):
        for key in ("legal", "iou_tc", "f5", "g_area"):
            assert a["samples"][key] == b["samples"][key]        # == on floats: bit equality

def test_wrong_manifest_exits_2_and_writes_nothing(tmp_path, direct_run):
    view, src = direct_run
    (view / "manifest.json").write_bytes((view / "manifest.json").read_bytes() + b" ")
    r = run_rescore(src, view, tmp_path / "re")
    assert r.returncode == 2
    assert not (tmp_path / "re").exists() or not any((tmp_path / "re").iterdir())
```

```python
# tests/test_a_select.py  (T6.1 + T6.4)
import inspect
import numpy as np
import pytest
from objdel.eval import select as SEL

M3 = np.array([[1, .5, .5], [.5, 1, .9], [.5, .9, 1]])

def test_medoid_tie_goes_to_lowest_index():
    assert SEL.medoid(M3, [True, True, True]) == 1               # scores .5, .7, .7

@pytest.mark.parametrize("legal,want", [([0, 0, 0], 0), ([0, 1, 0], 1), ([0, 0, 1], 2), ([1, 0, 1], 0)])
def test_medoid_eligibility(legal, want):
    assert SEL.medoid(M3, legal) == want                         # none -> 0, one -> itself, two -> tie -> lowest

def test_majority_of_legal_empty_answers_wins():
    K = 5
    iou = np.zeros((K, K)); iou[:3, :3] = 1.0; np.fill_diagonal(iou, 1.0)   # answers 0-2 empty, 3-4 distinct nonempty
    assert SEL.medoid(iou, [True] * K) == 0

def test_diagonal_is_ignored():
    for diag in (0.0, 1.0, 7.0):
        m = M3.copy(); np.fill_diagonal(m, diag)
        assert SEL.medoid(m, [True] * 3) == 1

def test_medoid_takes_no_rater_input_and_matches_brute_force():
    assert list(inspect.signature(SEL.medoid).parameters) == ["sample_iou", "legal"]
    rng = np.random.default_rng(0)
    for _ in range(2000):
        K = int(rng.integers(1, 9))
        A = rng.integers(0, 5, (K, K)) / 4; A = (A + A.T) / 2
        legal = rng.random(K) < 0.7
        el = [i for i in range(K) if legal[i]]
        if not el: want = 0
        elif len(el) == 1: want = el[0]
        else:
            sc = {i: np.mean([A[i, j] for j in el if j != i]) for i in el}
            want = max(el, key=lambda i: (sc[i], -i))
        assert SEL.medoid(A, legal) == want
```

```python
# tests/test_a_ceiling.py  (T7.2 + T7.12 + T7.14)  -- `disc_stack(...)` from conftest._disc; `cache_of` = ExpertCache.build
import numpy as np
from objdel.eval import ceiling as C
from objdel.reward.cache import ExpertCache
from objdel.reward import rewards   # or the module that has score_mask (read tests/test_reward*.py)

def plain_iou(a, b):
    return (a & b).sum() / (a | b).sum()

def test_k2_loo_scores_equal_plain_iou(two_rater_cache, rng_masks):
    cache, (A, B) = two_rater_cache
    loos = C.loo_caches(cache)
    assert [l.index for l in loos] == [0, 1] and all(l.single and l.cache.k == 2 for l in loos)
    for m in rng_masks(200):                                     # legal random discs
        assert abs(score_mask(m, loos[0].cache).iou_tc - plain_iou(m, B)) < 1e-12   # rater 0 left out -> other = B
        assert abs(score_mask(m, loos[1].cache).iou_tc - plain_iou(m, A)) < 1e-12

def test_k4_loo_removes_one_rater_and_keeps_order(four_rater_cache):
    cache, stack = four_rater_cache
    for l in C.loo_caches(cache):
        assert not l.single and l.cache.k == 3
        assert np.array_equal(l.cache.raters, np.delete(stack, l.index, axis=0))
    assert [l.index for l in C.loo_caches(cache)] == sorted(cache.voters)

def test_doctor_and_model_use_the_same_k_minus_1_doctors(monkeypatch, four_rater_cache, texts_k8):
    cache, _ = four_rater_cache
    seen = []
    real = score_mask
    monkeypatch.setattr(C, "score_mask", lambda m, c, *a, **k: (seen.append(id(c)), real(m, c, *a, **k))[1])
    C.item(cache, texts_k8, [False] * 8, 512, 512)
    loo_ids = {id(l.cache) for l in C.loo_caches(cache)}         # NOTE: loo_caches rebuilds; compare by digest instead
    assert {c for c in seen} <= loo_ids | {id(cache)}            # replace by `.digest` comparison in the real test

def test_identical_raters_give_doctor_one_and_diff_zero(identical_cache, identical_text):
    r = C.item(identical_cache, [identical_text] * 8, [False] * 8, 512, 512)
    assert r["doctor"] == 1.0 and r["model"] == 1.0 and r["diff"] == 0.0

def test_summarize_strata_add_up_and_pooled_differs_from_item_mean():
    rows = [{"k": 3, "single": False, "doctor": .8, "model": .7, "model_all_k": .7, "diff": -.1,
             "legacy_dmax": [.9] * 4, "legacy_loo_dmax_item": .9, "item_id": f"x{i}"} for i in range(6)] + \
           [{"k": 2, "single": True, "doctor": .6, "model": .5, "model_all_k": .5, "diff": -.1,
             "legacy_dmax": [.5, .5], "legacy_loo_dmax_item": .5, "item_id": f"y{i}"} for i in range(4)]
    s = C.summarize(rows, boot=200)
    assert s["all"]["n"] == s["k>=3"]["n"] + s["k=2"]["n"] == 10
    assert abs(s["all"]["legacy_loo_dmax_pooled"] - (6 * 4 * .9 + 4 * 2 * .5) / (6 * 4 + 4 * 2)) < 1e-12
    assert abs(s["all"]["legacy_loo_dmax_item"][0] - (6 * .9 + 4 * .5) / 10) < 1e-12
