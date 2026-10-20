# M5: SFT targets and the model frame (`objdel.targets`)

Plan: D11, D12 (LIDC 512 px), §4.4, gates K5, K6, K7, K10, R1 (target part). Promotes `lab/t_targets`
(`builder.py` 626 + `simplify.py` 447 + `fastgeom.py` 307 + `crops.py` 61 lines).

## What was built

| File | Lines (non-blank, incl. docstrings) | Contents |
|---|---|---|
| `src/objdel/targets/simplify.py` | 187 | `lattice_path`, `to_grid`, `simplify` (numba kernels `_graph`, `_best_cycle`) |
| `src/objdel/targets/build.py` | 155 | `BuildConfig`, `Target`, `build_target`, `target_set_hash` |
| `src/objdel/targets/frame.py` | 93 | `FrameConfig`, `FRAMES`, `Frame`, `frame` |
| `src/objdel/targets/__init__.py` | 8 | re-exports |
| `tests/test_targets.py`, `tests/test_targets_frame.py` | | 67 tests (see below) |
| `experiments/phase1/build_targets.py` | | full train/val build, K5–K7 / R1 summary, manifest with hash |

Dropped from the lab: `fastgeom.py` (certification is polycore through `objdel.geometry.to_mask`; polycore is
already checked against the exact oracles in M1, so a third geometry implementation adds nothing), the lab's
serializer/parser/canonical form (M2 owns them), Douglas–Peucker, Visvalingam–Whyatt, classic potrace, pin/split
repairs, `max_rounds`, `max_level`, `floor`, the `schema`/`derived` parameters and the non-divisor crop policies.

## API and semantics (contract, EXECUTION.md §3)

```python
frame(item, cfg=None) -> Frame(image, masks, W, H, box, scale)      # cfg defaults to FRAMES[item.dataset]
build_target(masks, names, cfg=None, min_area=None) -> Target(text, objs, info)   # cfg None = BuildConfig()
target_set_hash(rows: Iterable[(item_id, rater, text)]) -> str     # SHA-256, order-independent
```

**Frame.** Square crop, W = H = S.
- `masks[name]`: bool (k, H, W), one layer per rater: the ground truth for SFT and reward alike.
- `W`, `H`: the frame side S in pixels; S divides 1000, 50 ≤ S ≤ 500. The reward rasterizes at W × H.
  One frame pixel is exactly 1000 / S grid units, so every pixel corner is a grid point.
- `box = (x0, y0, x1, y1)`: the crop window in item pixels, half-open, `x1 - x0 = scale * W`; it may extend past
  the image (padding: masks False, image black). Python ints.
- `scale`: integer item pixels per frame pixel (1 unless the object needs a downscale).
- grid vertex g ↔ frame pixel coordinate g·W/1000 ↔ item pixel coordinate x0 + scale·g·W/1000 (same for y).
- `image`: uint8 RGB model input. RIGA: the crop at frame resolution (block mean if scale > 1). LIDC: Hounsfield
  units through the lung window [−1000, 400], then bicubically upscaled to 512 × 512 (image only); `None` when the
  item has no pixels.

Crop policy: side s = larger side of the union box of all raters and structures; need = s + 2·max(8 px, ⌈25 % s⌉);
f = ⌈need / 500⌉; S = smallest divisor of 1000 ≥ max(⌈need / f⌉, 50). The window of f·S item pixels is centred on
the union box and shifted inside the image when it fits. f > 1: f × f block majority on masks (ties to the
foreground), block mean on the image.

**Target.** `text` is the canonical compact text (M2 grammar) of all `names`; `objs` is that text parsed back;
`info[name] = {"tau": ladder rung (None = exact rings), "vertices", "iou"}`, `info["time"]` (CPU s).
`min_area` maps object names (and `"rim"` for RIGA) to the reward's size filter (min_area of the raters' union, as
M6's `ExpertCache` computes it); the target must keep the size-filtered signature. Default: `min_area` of the object's
own mask, never larger than the union rule, hence stricter. `ValueError` when no target can be certified: masks larger than
500 px, more than `MAX_RINGS` = 10 rings in one object, a ring that needs more than `MAX_VERTICES` = 1,000 vertices
on every rung (the text would not parse), or (never observed) a ladder that runs out. Empty mask → `{"nodule":[]}`, the legal empty object of M2.

## Research (online)

- **Potrace polygon stage** (Selinger 2003, §2.2, read from https://potrace.sourceforge.net/potrace.pdf):
  straight path = approximated within max-distance ½ by a segment and not all four directions used; all straight
  subpaths in O(n²) with a two-inequality constraint updated in O(1); optimal polygon lexicographic (#segments,
  penalty) where the penalty is the segment length times the RMS distance, computed from prefix sums; possible
  segments are clipped by one vertex at each end; potrace's `bestpolygon` fixes vertex 0 (non-cyclic).
  Kept: the constraint cone, min-#, penalty tie-break via prefix sums. Changed (as in the lab): tolerance τ in grid
  units instead of ½ px, no clipping, no vertex adjustment (vertices stay on pixel corners = exact grid points),
  cyclic optimum via an anchor (every cycle has a vertex at the anchor or a segment over it).
- **Digital straight segments**: Debled-Rennesson & Reveillès 1995, arithmetical DSS recognition, incremental and
  linear (https://doi.org/10.1142/S0218001495000249; DGtal implementation
  https://www.dgtal.org/doc/1.3/moduleArithDSSReco.html); maximal segments / tangential cover, Lachaud et al.
  (https://arxiv.org/abs/0906.2716). Our lossless test is the half-plane form of a DSS: the segment's line must
  separate the inner and outer pixel centres of every step, maintained as an exact integer cone (O(1) per step),
  which replaced the lab's float atan2 filter plus an O(l) integer re-check.
- **Min-# with an error bound**: Imai–Iri shortest path in the shortcut graph; survey and O(n²) min-# (Chan & Chin)
  in https://arxiv.org/abs/cs/0209007 and Kolesnikov & Fränti
  (https://cs.uef.fi/pub/franti/papers/Kolesnikov/Kolesnikov_Paper4.pdf). Simplifying a subdivision without
  crossings is hard to approximate (Estkowski & Mitchell), hence per-ring optimal + exact verification + repair.
- **Polygon outputs of VLMs, 2024–2026**: SimpleSeg (https://arxiv.org/abs/2601.19228, contour tracing + optional
  tolerance sparsification, unnamed simplifier, SFT→RL with IoU reward), Groundbench
  (https://arxiv.org/pdf/2609.27821, Visvalingam–Whyatt to exactly N vertices), procedural-geometry VLM grounding
  (https://arxiv.org/pdf/2604.02893, Douglas–Peucker, 0–255 tokens), VectorLLM (https://arxiv.org/abs/2507.04664,
  corner-by-corner regression of regular contours), PolyFormer (https://arxiv.org/abs/2302.07387), GeoGround
  (https://arxiv.org/abs/2411.11904, Text-Mask: an N×N grid with run-length encoding, not polygons). All of them
  use DP/VW or fixed-N resampling, which the lab already measured as worse than min-# on staircase boundaries (more
  vertices at lower IoU). No candidate better than minseg_hybrid was found, so nothing new was sandboxed against it
  except the repair strategy below.
- **Mask downscaling**: bilinear/bicubic-and-threshold, pooling and adaptive-crossing-number methods can all change
  topology; Chen & Peng (ECCV 2024, https://arxiv.org/abs/2407.17786) guarantee Betti numbers via integer
  programming. Not adopted: the downscaled mask *is* the ground truth for SFT and reward (plan), block majority is
  monotone (cup ⊆ disc survives), and no RIGA/LIDC item needs f > 1 (largest union 309 px), so the path is tested
  but unused.

## Sandbox (scratch `impl/M5/`)

**Repair strategy.** Lab: defect-localized pins for topology/legality failures, per-ring tolerance halving for IoU
failures, then also a pure halving run, keep the fewer vertices ("best"). Alternative tried: one tolerance ladder
per object, rebuild the object one rung lower on any failure (no localization). 3,169 fixture targets (pad_div
crops, all raters, `proto/compare.py`):

| method | RIGA disc vertices med / p95 / max (sum) | RIGA cup med / p95 / max (sum) | LIDC | certified | CPU / RIGA target |
|---|---|---|---|---|---|
| lab hyb3 best | 21 / 30 / 80 (26,147) | 19 / 40 / 89 (26,399) | 25,023 vertices | 100 % | 64.5 ms |
| ladder 3, 1.5, 0.75, 0.375, 0 (lab's rungs) | 21 / 30 / 144 (1.004×) | 19 / 41 / 87 (1.011×) | identical texts | 100 % | 22.3 ms |
| ladder 3, 2.5, …, 0.5, 0 | 18 / 21 / 33 (0.831×) | 17 / 30 / 82 (0.875×) | identical texts | 100 % | 28.2 ms |
| ladder 3, 2.75, …, 0.25, 0 (prototype) | 18 / 19 / 29 (0.811×) | 16 / 23 / 82 (0.787×) | identical texts | 100 % | 35.3 ms |
| **production `build_target`** (same ladder, crack-edge bevel, size rule) | **18 / 19 / 21 (0.809×)** | **16 / 22 / 37 (0.772×)** | 99.85 % identical, same vertex total | 100 % | 22 ms |

With the lab's own rungs the ladder reproduces the lab text on 98.1 % of RIGA targets and 100 % of LIDC targets;
the 0.25-step ladder needs 19–21 % fewer vertices than the lab at the same certification (IoU ≥ 0.99, signature,
legality; IoU medians 0.9914 / 0.9917). Decision: per-object ladder, step 0.25. This replaces the plan's
"defect-localized repair" with evidence (plan change, reported to the coordinator).

**numba vs alternatives** (`bench_kernel.py`, 904-step RIGA-like ring, S = 500): numba `simplify` 3.7 ms
(τ = 3), 1.1 ms (τ = 1), 0.8 ms (lossless); the same code without JIT (pure Python) 1,115 / 280 / 87 ms (≈ 300×);
a numpy-vectorized DP from a single start takes 8.8 ms, already 4× the numba cycle search over all anchor starts
(2.1 ms). C++ would need a second extension (M1 owns the build). Decision: numba (plan D14), `cache=True`.

**AB5 preview (lossless vs hybrid)** on fixtures (every 4th RIGA, every 3rd LIDC item, `ab5.py`): RIGA lossless
needs 3.66× the vertices (median 59 vs 17 per object); on LIDC both are identical (crops of 40–100 px: τ = 3 grid
units is ≤ 0.3 px, so every admissible shortcut is lossless anyway).

**Two defects found on the store build, each fixed test-first.**
1. *Bevel on slanted edges (lab bug).* The lab split a pinch vertex one grid unit along the simplified edges; on a
   45° edge that rounds to (1, 1), which at S = 500 is a pixel centre, so the tie rule flips it. Every rung then
   failed and 7 RIGA objects and 5 LIDC nodules fell back to exact rings (answers up to 9,260 tokens). Now the
   split is one unit back along the pinch's incoming crack step and one unit along the outgoing one: the cut stays
   in the background pixel's corner. Regression test: two ellipses touching at a corner (found by random search;
   the lossless rung must reproduce the mask exactly). LIDC fallbacks went from 5 to 0.
2. *Raw vs size-filtered signature.* The lab certified the raw (b0, b1). On jagged pen-stroke masks (RIGA BinRushed,
   rebuilt by M4's F1/F4 stages) the τ = 3 polygon rasterizes with 1–3 px pinholes or specks; the filtered signature
   is unchanged but the raw one is not, so 1–2 % of RIGA objects went down to τ = 0.5 (150–550 vertices, answers
   up to 6,036 tokens, beyond `max_new_tokens` = 1,536). Plan §4.2 applies the size rule to targets, reward and eval,
   and K5 is defined on it, so certification now uses the size-filtered signature only. Regression test: a seeded
   jagged disc whose τ = 3 raster has two sub-filter pinholes is certified at τ = 3 with < 150 vertices (it took
   τ = 1.25 and 350 vertices before). Raw signature is still reported: 99.2 % (disc), 99.5 % (cup), 100 % (LIDC).

## Tests (written first, run red, then green)

`tests/test_targets.py` (57) and `tests/test_targets_frame.py` (10), ~15 s:
- every target is re-certified from its text alone (`check`): parses, canonical, legal (polycore), size-filtered
  signature equal, IoU ≥ 0.99, RIGA rim signature equal (R1 target part, K5, K6);
- regression tests for the two defects above; a ring beyond `MAX_VERTICES` is refused;
- edge cases × crop sides (own frame, 50, 500): 1-px object, 1-px horizontal / vertical / diagonal lines, both
  pinches, checkerboard, holes touching diagonally, a hole touching the outside, island in a hole (two parts),
  staircase; empty masks (`{"nodule":[]}`); objects touching every border with a 1-px pinhole (exact text);
  exact pixel-corner vertices of a rectangle; RIGA-like disc + cup at 500 px (< 60 vertices each); cup touching the
  disc (rim); a wiggly boundary that needs τ < 3;
- Hypothesis: 150 random lab-family masks (rects, holes, pixels, diagonals, checkers, lines) placed in crops of
  own / 50 / 125 / 333 / 500 px with random min_area; 40 smooth-noise blobs (components, holes, necks) at 100–500 px;
  masks beyond `MAX_RINGS` must raise;
- determinism: same text twice, input untouched, same SHA-256 in two processes with different `PYTHONHASHSEED`;
- invalid inputs; `target_set_hash` order-independence and sensitivity;
- `@pytest.mark.data` differential vs the lab builder on fixtures (every 20th RIGA, every 6th LIDC item, all
  raters): LIDC texts identical, RIGA total vertices ≤ the lab's;
- frame: exact slice at scale 1 (masks and image), integer pixel ↔ grid mapping of a target, padding at the image
  border, LIDC window + 512 px RGB image with masks at frame size, `None` image, block-majority downscale (vs a
  direct block vote) with nesting kept and block-mean image, even-factor ties to foreground, reproducibility, empty
  union raises, Hypothesis policy properties (side divides 1000, 50 ≤ S ≤ 500, union + margin inside the box).

## Gates on the full train/val store

Command: `.venv/bin/python -I -m experiments.phase1.build_targets` (store manifest SHA-256 `56c1d339…`, 3 worker
processes; 140 s wall under load from other agents, 76–97 s otherwise). Outputs (git-ignored): `data/targets/<dataset>/<split>.jsonl`
(item_id, rater, W, H, box, scale, text), `data/targets/manifest.json` (hash per file and overall, config,
versions), `data/targets/summary.json`. No item raised; every text parses.

| split | items | targets | frame side S (items) | K5 (size-filtered signature) | legal | IoU min / median | vertices med / p95 / max |
|---|---|---|---|---|---|---|---|
| RIGA train | 591 | 3,546 | 200: 2, 250: 162, 500: 427 | disc, cup, rim 100 % | 100 % | disc 0.9900 / 0.9922, cup 0.9900 / 0.9915 | disc 19 / 24 / 186, cup 17 / 24 / 198 |
| RIGA val | 64 | 384 | 250: 22, 500: 42 | 100 % | 100 % | disc 0.9901 / 0.9924, cup 0.9900 / 0.9916 | disc 19 / 25 / 40, cup 17 / 25 / 97 |
| LIDC train | 7,565 | 24,301 | 50: 6,550, 100: 1,009, 125: 6 | 100 % | 100 % | 1.0 / 1.0 | 11 / 31 / 145 |
| LIDC val | 961 | 3,126 | 50: 844, 100: 117 | 100 % | 100 % | 1.0 / 1.0 | 12 / 33 / 155 |

**Superseded by** the rebuild on the final store (`fixes/F4.md`, store manifest `37fe7665…`): `data/targets` RIGA
train 592 items / 3,552 targets, val 63 / 378; LIDC train 7,341 / 23,469, val 1,084 / 3,444 (30,843 targets);
target_set_sha256 `1af819f4…` (`data/targets/manifest.json`); errors 0, parse 100 %.

No item needs a downscale (scale 1 everywhere). Tolerance reached (RIGA train): disc τ = 3 on 78 %, ≥ 2.75 on 94 %;
cup spread over τ = 3 … 0.75 (median 1.75–2): the cup's IoU ≥ 0.99 binds. LIDC: τ = 3 on 100 % (crops ≤ 125 px make
every admissible shortcut lossless, IoU 1.0).

**K7: answer tokens** (whole target text, no special tokens), median / p95 / p99 / max:

| tokenizer | RIGA train | RIGA val | LIDC train | LIDC val |
|---|---|---|---|---|
| SmolVLM-256M (`HuggingFaceTB/SmolVLM-256M-Instruct`) | 316 / 388 / 480 / 2,268 | 316 / 396 / 667 / 1,124 | 103 / 264 / 431 / 1,186 | 111 / 279 / 416 / 1,259 |
| LFM2.5-VL-450M (`LiquidAI/LFM2.5-VL-450M`) | 207 / 252 / 310 / 1,427 | 207 / 257 / 426 / 712 | 70 / 171 / 275 / 750 | 75 / 180 / 266 / 794 |
| Qwen3-VL-2B (`Qwen/Qwen3-VL-2B-Instruct`) | 309 / 381 / 473 / 2,261 | 309 / 389 / 660 / 1,117 | 100 / 260 / 428 / 1,181 | 108 / 276 / 413 / 1,259 |

Per object (K7 gate ≤ 250), median Qwen / SmolVLM / LFM tokens: disc 163 / 166 / 109, cup 147 / 150 / 99,
nodule 100 / 103 / 70 (train; val within 8). Lab: disc 157, cup 197, nodule 92. Answers above 1,536 Qwen tokens:
10 / 3,546 RIGA train, 0 elsewhere; above 1,024: 24 RIGA train, 1 RIGA val, 3 LIDC.

**Build time** (CPU per target): RIGA median 29 ms, p95 47 ms, max 0.37 s; LIDC median 1.2 ms, max 6 ms.
Total 145 CPU s for 31,357 targets.

**K10: target-set SHA-256** `d9d4d932b764f433570c5924b65b46e6a834881689aa8ee9e236ee5ab1fa25af` (31,357 targets,
train + val; superseded: `1af819f4…`, 30,843 targets, final store), identical in two full runs with 2 and 3 worker processes. Per file: in `data/targets/manifest.json`.
Pins: Python 3.12.11, numpy 2.5.4, numba 0.68.0, opencv-python-headless 5.0.0.93, objdel 0.1.0, transformers 5.20.1
(tokens only). The hash is over (item_id, rater, text), sorted, so it is independent of order and worker count.
Penalties are exact integers divided once (IEEE), so the min-# tie-break should not depend on FMA or the platform;
not yet verified on x86 (M9).

**R1 (target part, plus M6's terms)** in the sandbox (`r1_full.py`, M6's `ExpertCache` / `score_completion` as of
this run; M6 is still in progress): every train/val target (31,357) parsed and scored against its item's rater
caches: legal 100 %, size-filtered signature = the rater's 100 %, IoU_tc ≥ 0.99 100 % (min 0.99005), F5 ≥ 0.95
100 % (min 0.988), topo_any 100 %, every term identical after rotating ring starts and reordering holes 100 %.

`make lint`: clean. `make test`: 363 passed, 4 failed, all in `tests/test_reward_gates.py` (M6, in progress).

## Open issues

- **numba is not a declared dependency** (installed in `.venv`, pinned in EXECUTION.md §4, used by plan D14). CI
  installs `.[dev]` only, so `tests/test_targets*.py` fail to import there until `numba==0.68.0` is added to
  `pyproject.toml` (coordinator). *(Resolved: `numba==0.68.0` is a core dependency in `pyproject.toml`.)*
- **Long RIGA answers:** 10 RIGA train targets exceed 1,536 Qwen tokens (max 2,261); the longest are BinRushed1
  (jagged reconstructed masks; their cup IoU needs τ ≤ 1). M7 should drop or cap them, or raise `max_new_tokens`.
  *(Resolved: M7 drops them from SFT (I5); since B3-1 BinRushed1 has no model image and is not prepared, so 1 /
  3,282 prepared RIGA train targets is over budget for every model.)*
- **Rim certification** is kept (lab K5 counted rims). M6's reward has no rim term (joint disc/cup signatures only),
  so the rim check is stricter than the reward needs; how many targets it lengthens was not measured. Drop it if M6
  keeps no rim term.
- **Frame jitter:** the crop is deterministic and nearly centred on the union (shifted only at image borders); SFT
  may want random crop jitter (VectorLLM 1.1–1.5× box) to avoid a centre prior. Not built: it is an M7 augmentation
  choice and would change the ground-truth frame per sample.
- **Margin and LIDC side:** margin max(8 px, 25 %) and min side 50 px were chosen, not tuned; LIDC frames are
  50–125 px (87 % at 50 px) and the model image is 512 px. The Phase 2 precision probe should confirm.
- **Lossy targets are not idempotent** (rebuilding from a target's raster can give a different text), as in the lab.
- **x86 determinism** of the hash is untested (M9). *(Checked in `M9_x86.md`: target texts and `target_set_hash` identical on x86-64 in Docker, a sample corpus; LIDC frame images differ.)*
