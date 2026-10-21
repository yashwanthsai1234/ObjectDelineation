# BK1 break-it findings

Scope: the data path store → `objdel.targets.frame` → `augment` → `build_views` / `build_eval_views` →
`dataset.view_rows` → TRL vision collator → `vast/bundle.py`. Tests: `tests/test_p2_breakit_data.py` (28 tests:
21 pass as regression guards, 7 strict xfail = the defects below; 2 of them `data`; 14 s on the Mac). Scratch
evidence in `<scratchpad>/bk1/` (`s11.py`, `s11fix.py`, `find200.py`, `a0real.py`, `coll.py`).

## Findings

No high-severity finding: no wrong training data at scale and no crash in a planned Phase 2 run.

| id | sev | component | repro (test) | expected vs actual | suggested fix | owner |
|---|---|---|---|---|---|---|
| BK1-1 | low | M1 targets (S1-1 root cause) | `test_s200_lidc_targets_are_pixel_exact` (xfail); fix evidence `test_s200_targets_are_exact_with_tau_below_half_a_pixel` | D9 / §6.0 say LIDC targets are pixel-exact. Actual: certified (IoU ≥ 0.99) but not exact. All 11 rows are at **S = 200**, which is 11 of 55 S = 200 rows (20 %). S = 100 and 125 are 100 % exact (6,576 rows). Cause: `BuildConfig.tau = 3` is in grid units, and one pixel is 1000 / S units. So tau is 0.15 / 0.30 / 0.375 px at S = 50 / 100 / 125, but 0.6 px at S = 200. Every ladder rung down to 2.5 is still ≥ 0.5 px, and an L-inf shortcut can cross a pixel centre once tau ≥ ½ px. The ladder steps down only on IoU, signature or legality, never on exactness. Side effect: D4 variants of one mask get different rasters (z240: `id` exact at tau 2.25, `r270` not at 2.5). Synthetic: 26 / 40 rotated ellipses fail at S = 200, 0 / 40 at S = 125 | In `build_views._item`, use `BuildConfig(tau=min(3.0, 0.45 * 1000 / S))` for LIDC (2.25 at S = 200). Real data: 55 / 55 S = 200 rows become exact, +2.5 % characters on those rows. Then assert 100 % exact for LIDC in the `data` test | V |
| BK1-2 | medium | M7 bundle × M2 `cache_loader` | `test_bundled_prepared_root_loads_with_cache_loader` (xfail) | Expected: the bundled `prepared/` root loads the same as on the Mac. Actual: `bundle-lidc.tar` ships `prepared/manifest.json`, which lists `riga/{train,val}.jsonl`, but not the riga files themselves. `dataset.cache_loader` reads every listed file, so on the box it raises `FileNotFoundError: prepared/riga/train.jsonl`. Hit by `eval_run` on the centred view and by `smoke_grpo`. `serve.preload(root, ["lidc"])` works (guard `test_bundled_prepared_root_preloads_on_the_box`). `select_ckpt` and `g_probe` use preload, so no planned Stage-2 script crashes today | Either `bundle` writes a dataset-filtered `prepared/manifest.json` (only the `files` / `sets` of `ds`), or `cache_loader(root, datasets=None)` skips manifest entries of other datasets, as preload does | O (bundle) or T |
| BK1-3 | medium | M1 / M2 view integrity | `test_view_rows_refuse_an_edited_context` (xfail) | Expected: as for prepared PNGs (B3-5, `image_sha256`), a context that is not the one the rows were built on is refused. Actual: `view_hash` covers only `rows.jsonl`. `manifest.json` has no context hashes, so an edited, swapped or rebuilt `context/<item>.npz` with the same origin trains silently. A missing one crashes inside a DataLoader worker at the first render, possibly hours into a run. On the box, `SHA256SUMS` catches transport corruption only, not a context rebuilt from a changed store | `build_views`: add `context_sha256: {item_id: sha}` to the manifest, or fold the sorted context hashes into `view_hash`. `view_rows`: check that every context of the rows exists and matches (7,341 files, 210 MB: well under a second of SHA-256) | V (manifest), T (check) |
| BK1-4 | low | M2 `view_rows(a0=True)` | `test_a0_island_inside_a_hole_stays_legal` (xfail) | Expected: A0 targets are legal. Actual: a part inside another part's hole (a ring with an island) becomes two nested outer rings when holes are emptied. `to_mask` then reports `legal=False, nested_outers=1`, and that illegal text would be taught. Real LIDC train: 0 of 1,656 non-simple rows (`test_real_a0_targets_stay_legal` passes), so the defect is latent | `_no_holes`: drop parts whose outer ring lies inside another part's outer ring (fill = union), or rebuild A0 from the filled mask with `build_target` (`binary_fill_holes`) | T |
| BK1-5 | low | M2 `view_rows` | `test_max_tokens_without_a_tokenizer_is_refused` (xfail) | Expected: `max_tokens` is a token budget, so without a tokenizer the call raises. `sft_rows` raises `TypeError` in that case. Actual: it silently filters on character count (`length` = characters). Characters outnumber tokens, so it drops rows that fit the token budget | Raise `ValueError` when `max_tokens` is set and `tokenizer is None` | T |
| BK1-6 | low | M1 `frame` / `hu_window` | `test_degenerate_hu_window_is_refused` (xfail) | Expected: `ValueError` for a window with `hi <= lo`. Actual: `(0, 0)` gives a black image with a numpy divide-by-zero warning. An inverted window gives garbage values. `jitter` never draws these, so the risk comes only from hand-written specs or configs | `if hi <= lo: raise ValueError` in `hu_window` | V |
| BK1-7 | low | M7 bundle | `test_bundle_carries_no_wall_clock_file` (xfail) | Expected: two builds of the same view from the same store give the same bundle hash. Actual: `views/<v>/build_report.json`, which holds wall-clock seconds, is bundled, so the bundle SHA changes on every view rebuild | Skip `build_report.json` in `bundle.walk`, or have `build_views` write the report outside the view dir | O |

## Guards that pass (edge cases tried, no defect)

- Border crops: unions on or next to every image edge and corner. Moved windows pad, and `render` equals the offline
  reference for all 8 ops (`test_border_windows_pad_and_render_like_the_reference`).
- Clearance limit: offsets at exactly 2 px on two sides are accepted, and one px more is refused on each of 4 sides,
  at S = 50, 100 and 125.
- Side across all LADDER divisors: masks are whole and the box is S. Non-divisors and too-small sides are refused.
- D4 × (multi-part + hole, ring + island) at S = 50, 100 and 125: target → raster → inverse op equals the original mask
  pixel for pixel, with the same signature. Inverse round trip also holds on non-square arrays. RIGA cup ⊂ disc
  under `fh` keeps the rim signature.
- HU extremes (int16 −32768 / 32767), σ = 10 and the window jitter's four corner values: render equals the reference.
  2,000 `jitter` draws stay in bounds.
- Contexts cover every admissible window: 60 random draws per item, plus the four extreme corners at the largest
  drawable S, on 12 items, including crowded and edge items.
- `build_views` with 12 items (more than the `imap` chunksize of 8):
  - epoch 0 of a 3-epoch / 1-process build equals a 1-epoch / 2-process build, and contexts are byte-identical;
  - every epoch has fresh draws;
  - row count = Σ voting raters × (1 | 4) × epochs.
- End to end per row (S = 50, 100, 125, 200):
  - `box = frame(offset, side).box`;
  - render equals the reference;
  - the rendered 512-px image, sampled at frame pixel centres, is bright exactly on the op-transformed mask, so
    image and target agree in orientation and position;
  - the target is certified, and exact for S ≤ 125.
- `max_tokens` boundary: `n + 1 ≤ budget` is kept and `n + 1 > budget` is dropped. A0 on multi-part + hole keeps 2
  parts, fills the hole and stays legal. A reordered or extended `rows.jsonl` is refused.
- RIGA view through `view_rows`: `id` / `fh` only, no window or noise, the RIGA prompt, an (S, S, 3) image equal to the
  reference.
- TRL 1.16 `DataCollatorForVisionLanguageModeling` (cached SmolVLM-256M processor, offline) on
  `ds.__getitems__([0, 1, 2])`: `pixel_values` batch 3, labels present.
- Real data (`data`):
  - the first 25 store items' rows (box, render vs the store HU reference, legality and signature);
  - `view_hash` and the row count;
  - A0 legal on all 1,656 non-simple rows.

Not re-reported (known): S1-2, S1-7, S1-8. Note for S1-8: `train_sft` raises when the view is too short, so a
1-epoch view fails fast and does not train wrongly.

## Simplicity audit (data path, ~800 lines: frame 182, augment 159, build_views 174, build_eval_views 107,
dataset view part ~60, bundle 117)

No 200 → 50 case. The code is already tight. These cuts save about 50 lines and remove two duplicated rules:

1. `build_views.py:46-60, 96-99` `recheck`: `build_target` already certifies legal, IoU, signature and rim, and
   raises otherwise. So `certified` is always True (24,701 / 24,701), and the whole re-parse can become one exactness
   line: `exact = all(np.array_equal(to_mask(t.objs[n], S, S).mask, src[n]) for n in names)`. About −14 lines.
2. `build_eval_views.py:69-86` repeats `dataset.prepare`'s per-set summary (`dataset.py:143-152`) and calls the
   private `_prepare` (S1-2). A shared `dataset.summarize_set(recs, skipped)` would remove about 12 lines and the
   coupling.
3. `augment.py:142-159` `render` re-implements `frame.py:171-181`: crop, HU window, black padding, channel repeat,
   upscale. A shared `frame._pixels(image, x0, y0, S, window)` would serve both and keep them from drifting. About
   −8 lines.
4. `augment.py:152` builds a full-image `np.ones((1, H, W))` mask on every render just to find the padding. The
   in-image slice `ix0:ix1, iy0:iy1` is already computed on line 147, so 2 lines can replace it. No functional
   change; saves an H × W allocation per row.
5. `build_views.py:69` repeats `reward.rater_answers`' voting rule (`score.py:109`, the `... or range(k)` fallback). A
   `reward.voting(sigs_by_name)` helper used by both would keep one rule.
6. `build_views.py:105-112` `_jobs` passes back a mutable `skipped` list and does its own `limit` counting.
   `itertools.islice` plus counting `image is None` in the parent is about −4 lines.
7. `augment.py:71-79` `_margins` and `_side` are always called together. One `_max_side(box, dataset)` helper for
   `context_box`, with `draw` inlining the margin, saves about 4 lines.
8. `bundle.py:44-45, 60` checks `is_raw` twice: the per-file check makes the root check redundant, except for its
   message.
