# Break-it B1: geometry (parser, polycore, topology, targets and their wiring)

Code under attack: `cpp/polycore/polycore.hpp`, `src/objdel/_polycore.cpp`, `src/objdel/checker/{polycore,topology,
bitquad,certify,polygonize}.py`, `src/objdel/geometry/`, `src/objdel/targets/`, and the paths parse → to_mask →
topology/signature and mask → build_target → parse → to_mask (→ reward).

Tests: `tests/test_breakit_geometry.py`. Each confirmed breakdown has a strict xfail (`B1-<n>` in the reason); the
unmarked tests are attacks that held, kept as regression tests (37 tests, ~5 s in total).
Scratch scripts: `scratchpad/breakit/B1/p*.py`.

## Findings

| id | severity | module | one line |
|---|---|---|---|
| B1-1 | low | geometry.parse | Non-ASCII and control "whitespace" around the document is stripped and the text accepted |
| B1-2 | low | checker.polycore.pack | Rings of shape (k, 3), (k, 1) or 3-D are silently reshaped to (-1, 2): wrong vertices, no error |
| B1-3 | low | targets.build_target | The default `min_area` (own mask) is not "stricter" than a larger union rule, as the contract claims |
| B1-4 | low | targets.build_target | `MAX_RINGS` is checked on the raw mask before the size rule, so sub-filter pinholes make a certifiable mask uncertifiable |
| B1-5 | low (deepens I1) | targets, reward.cache | I1 reaches `build_target` (refuses the mask) and `ExpertCache` (rater signatures lose their holes); any True byte ≠ 1, not only 255 |
| B1-6 | low | targets.frame | A 4-channel image becomes a 0-channel "RGB" frame image; a 2-channel one stays 2-channel |

None of them is reached by the production pipeline today (store masks are canonical bools from `np.unpackbits`,
callers pass the union `min_area`, RIGA images are RGB, no store item raised in M5's build); all six are contract
gaps that a new caller or data source would fall into silently.

### B1-1 (low) geometry.parse: Unicode / control whitespace is stripped

- **Reproduction:** `test_b1_1_only_json_whitespace_is_stripped` (7 cases).
  `parse("　" + doc + "　", ("nodule",)).ok` is `True`; so are U+2028, U+00A0, U+0085, `\x1c`–`\x1f`, `\x0b`.
- **Expected:** only JSON whitespace (space, tab, LF, CR) is "surrounding whitespace"; any other non-ASCII character
  gives `not_ascii` (as `"﻿" + doc` already does), any other control character `syntax`.
- **Actual:** accepted. The accept set is therefore larger than "grammar + JSON whitespace", and differs from the
  constrained-decoding regex in more than the documented ASCII-whitespace strip.
- **Root cause:** `text.strip()` with no argument strips every character with `str.isspace()`, which in CPython
  includes the information separators `\x1c`–`\x1f` and Unicode `Zs`/`B`/`S` characters.
- **Suggested fix:** `s = text.strip(" \t\n\r")`. Sources: RFC 8259 §2 (insignificant whitespace is exactly
  `%x20 / %x09 / %x0A / %x0D`) https://datatracker.ietf.org/doc/html/rfc8259#section-2 ; CPython `strip` uses
  `_PyUnicode_IsWhitespace`, which includes U+001C..U+001F beyond the Unicode White_Space property
  https://bugs.python.org/msg147603 , https://github.com/python/cpython/issues/83599 .

### B1-2 (low) checker.polycore.pack: any integer array is reshaped to (-1, 2)

- **Reproduction:** `test_b1_2_pack_rejects_rings_that_are_not_k_by_2`.
  `polycore.rasterize([np.array([[0, 0, 1000], [0, 1000, 1000]])], 4, 4)["area_evenodd"] == 10` — the 2×3 array is
  read as the triangle (0,0), (1000,0), (1000,1000).
- **Expected:** a ring is (k, 2); anything else raises `ValueError`.
- **Actual:** `_as_i32(a, 2)` calls `a.reshape(-1, 2)` for every shape with an even size, so (k, 3), (k, 1), (2, 2, 2)
  inputs silently become different polygons. (A related ambiguity: a `tuple` of two ndarray rings is taken as a packed
  `(xy, off)` pair; it usually ends in a confusing offsets `ValueError`.)
- **Suggested fix:** in `_as_i32(..., ndim=2)` accept only `a.ndim == 2 and a.shape[1] == 2` (and, if flat lists must
  stay supported, `a.ndim == 1` with even length); raise `ValueError` otherwise.

### B1-3 (low) targets.build_target: the default min_area is not stricter

- **Reproduction:** `test_b1_3_default_min_area_target_keeps_signature_under_a_larger_union_rule`. A 500 px crop, a
  disc of radius 150 with one elliptic hole (371 px). With the default (`min_area` of the mask = 352) the target is
  certified, but its τ = 3 polygon rasterizes the hole as 369 px. Under a union rule of 371 (a union ≈ 6 % larger than
  this rater's mask) the rater keeps the hole, signature (1,), and the target loses it, (0,).
- **Expected (contract, EXECUTION.md §3 and M5 note):** "default = the object's own mask, never larger than the
  union rule, hence stricter".
- **Actual:** signature equality at threshold a does not imply equality at b > a: any hole or component whose size
  lies between the mask's and the polygon's pixel counts flips at some b. Random search: 15 / 240 random 500 px blobs
  with holes break at some union rule in [a, a + 60] (`p11.py`).
- **Impact:** production callers (`build_targets.py`, `train/dataset.py`) pass the union `min_area`, so the shipped
  targets are unaffected; a caller who relies on the documented default gets targets the reward can score as
  topologically wrong.
- **Suggested fix:** make `min_area` required (or certify at every threshold from the default up to the largest one
  in use, which is what "stricter" would need), and correct the docstring and EXECUTION.md §3.

### B1-4 (low) targets.build_target: ring cap before the size rule

- **Reproduction:** `test_b1_4_sub_filter_pinholes_do_not_make_a_mask_uncertifiable`. A disc with ten 1-px pinholes
  has 11 raw rings, so `build_target` raises `11 rings > MAX_RINGS = 10`. The filtered disc is a legal polygon with
  IoU 0.998 against the raw mask and the same size-filtered signature (0,), which is everything certification asks.
- **Expected:** the mask is certifiable (certification uses only the size-filtered signature, M5 decision 2).
- **Actual:** refused before any rung is tried.
- **Root cause:** `len(_rings(src[n])) > MAX_RINGS` is checked on the raw mask; sub-filter specks and pinholes count
  toward the cap although the size rule ignores them.
- **Suggested fix:** when the raw ring count exceeds the cap, build from `size_filter(mask, size[n])` (still certified
  against the raw `ref[n]` via signature and IoU); raise only if the filtered mask still has > MAX_RINGS rings.
  Source: the size rule as an area opening/closing (Vincent 1994, cited in `checker/topology.py`).

### B1-5 (low, deepens I1) non-canonical bool masks in targets and the reward cache

- **Reproduction:** `test_b1_5_targets_normalise_bool_bytes`, `test_b1_5_expert_cache_normalises_bool_bytes` (bytes
  2 and 255). A 50 px square with one hole whose True bytes are 2: `np.array_equal(nb, m)` is True, but
  `build_target` raises `no certified target` and `ExpertCache.build(np.stack([nb, nb])).rater_sigs == ((0,), (0,))`
  instead of `((1,), (1,))`. `size_filter` also leaves a 1-px hole unfilled at `min_area` 2 (`p8.py`).
- **Expected:** a bool array is a bool array; results equal those of the canonical mask.
- **Actual:** I1 (topology reads raw bytes) propagates: `np.ascontiguousarray(x, bool)` is a no-op on a bool input, so
  the 2/255 bytes reach `.view(np.uint8)` and `1 - mp` in `checker.topology` (1 − 2 = 255 marks the foreground as
  background in the hole labelling). It is not specific to 255.
- **Suggested fix:** fix I1 once at the entry of `_crop_pad`/`size_filter` (`m.view(np.uint8) != 0`, or
  `np.asarray(mask) != 0` for non-bool input); add a byte-2 case to the I1 regression test. Source: numpy gives no
  guarantee for bool values other than 0/1 (comparable note in JAX `Array.view`:
  https://docs.jax.dev/en/latest/_autosummary/jax.Array.view.html ).

### B1-6 (low) targets.frame: channel handling

- **Reproduction:** `test_b1_6_frame_image_is_rgb_or_refused`. A (80, 80, 4) uint8 image gives a (50, 50, 0) frame
  image; (80, 80, 2) gives (50, 50, 2).
- **Expected:** uint8 RGB (contract), or a `ValueError` for unsupported layouts.
- **Actual:** `np.repeat(image, 3 // len(image), 0)` repeats 0 times for 4 channels and once for 2.
- **Suggested fix:** accept 1 or 3 channels (drop alpha explicitly if wanted) and raise `ValueError` otherwise.

## Attacks that held

Parser
- Independent recursive-descent reference vs `parse` on 210,000 mutated serializations (insert / delete / replace /
  cut; alphabet includes Unicode digits `١`, full-width `０`, NUL, BOM, backslash, `-0`, `00`, `1000`): 0 accept/reject
  or value mismatches (`p3.py`; 3,000 kept as a test).
- Timing, 4 doublings each (2k–16k units, up to 192 kB): long ring with bad end, many holes + junk, many parts + wrong
  key, over-cap ring, 100-digit runs, `[`×n`]`×n, `{"disc":`×n, deep JSON objects, 4-digit tail, valid-JSON floats,
  far leading zero, deep valid lists: all linear (≤ 5.6 ms at 192 kB). Error codes vary with depth (`keys` vs
  `syntax` past the recursion limit) but parse never raises.
- Lone surrogate, BOM, NUL, CRLF, `-0`, escapes in keys, `1e2`, `1.0`, `true`, `null`, `0000`, `1001`, non-`str`
  inputs (None, bytes, bytearray, int, list): all rejected with a code, no exception.
- Duplicate keys / reordered keys / JSON-valid-but-grammar-invalid texts: rejected (covered by the reference fuzz).

polycore
- Exactness vs an independent `Fraction` pixel-centre oracle on 4,500 boundary-biased shapes (vertices at 0, 1000,
  pixel corners and centre lines, slightly out of range; W, H ∈ {1, 2, 3, 5, 7, 8, 13, 27, 40, 64}): 0 wrong pixels,
  both fill rules (`p6.py`; 150 kept as a test).
- W/H at 1 and 4096 (1×1, 4096×1, 4096²), 0, −1, 4097, 2³¹, floats, numpy ints, G at 0 / 2³⁰ / 2³⁰+1 / 2⁶³: correct
  values or `ValueError`/`TypeError`.
- Float / NaN / bool / object / huge (uint64, int64 > int32, Python 2⁷⁰) rings: `TypeError`/`OverflowError`.
  Fortran, negative-stride and strided rings and packed arrays: correct.
- Malformed packed offsets (negative, beyond N, empty, decreasing, 11 rings, float), 10 vs 11 rings, 1,000 vs 1,001
  vertices, empty rings, all-equal rings, roles (wrong length, self, chain, float, bool, 2⁴⁰, −2): all as specified.
- `score_batch`: empty batch, length mismatch, wrong gt shape, 3-D gt, NaN gt, `None` gt, threads −5 / 10,000 / 2³¹,
  numpy / generator W lists.
- 12 Python threads × 3 passes of rasterize + validate + `score_batch(threads=0..28)` on 300 random shapes: 0
  mismatches (thread-local workspace safe). 100k rasterize + validate calls: max RSS unchanged (180.6 MB).
- Validity hand cases: hole touching outer at a vertex, hole vertex on outer edge, hole outside, nested holes, island
  in hole (legal), island touching hole, identical outers, CCW outer, spike, consecutive/closing duplicates (legal),
  collinear 3, all-equal, zero-area hole, straight-through vertex (legal), coordinates −1, hole = outer, outer in
  outer, parts sharing an edge or a corner, figure-8 pinch: all classified as specified.
- Cost at the caps: 10 rings × 1,000 vertices, every edge pair crossing (5 M crossings): `to_mask` 0.64 s at 500²
  (validate 0.54 s). Bounded, and the text is 95 k characters, far beyond the 1,536-token completion budget.

Topology
- Every dtype (uint8/16, int8/16/32/64, uint64, float16/32/64, complex64, bool, object, uint8 0/255, ×256, negative,
  0.5, NaN), lists, Fortran, negative strides: consistent with `!= 0` semantics (NaN counts as foreground).
- 0×0, 0×5, 1×1, 1-D / 3-D / 0-d (ValueError; the 0-d message reports shape (1,)), 40,000×40, 3×70,000.
- Checkerboard 50×50 (one component, 1,152 holes, χ8 certificate holds), diagonal-only diamond, notch at the border
  (no hole), background diagonal to the border (a hole under 4-connectivity), island in hole in island.
- 450 masks up to 700×700 from 5 families (i.i.d., smooth, checker, grid lines, noisy) at OpenCV 1, 4 and 8 threads
  vs the scipy oracle: topology, hole areas, χ8, `signature(m, a)` and `size_filter(m, a)` for random a — 0 mismatches.
- Size rule exactly at the boundary (a and a − 1 px holes and components, a ∈ {5, 6, 17}).

Targets and wiring
- 4,200 + 2,100 `build_target` calls on 14 adversarial families (blobs, 1-px horizontal/vertical/diagonal lines,
  width-1 rings, combs, spirals, many 1-px holes, objects in image corners, full masks, checker patches, isolated
  pixels, noisy blobs) at sides 1–500 incl. non-divisors (3, 7, 13, 73, 333, 499) and rectangles, with the default or
  a larger `min_area`: no exception other than `ValueError` (> MAX_RINGS, or the vertex cap on a 250-tooth comb), and
  every returned text re-certified from the text alone (parse, legal, IoU ≥ 0.99, size-filtered signature).
- mask → target → `score_completion` on its own 2-rater `ExpertCache` (493 targets): format, legal, topo-any and
  IoU_tc ≥ 0.99 every time; the dataset's `min_area` equals the cache's.
- Determinism: 60 targets hashed identically across `PYTHONHASHSEED` 1 / 987 and `NUMBA_NUM_THREADS` 1.

Counts: **27 attack groups tried, 21 held (the bullets above), 6 broke (B1-1 … B1-6)**; about 160 individual cases
in the held groups (fuzz runs counted once).
