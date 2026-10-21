# M4b implementation note (agent P): probes and baselines

## API as built

`objdel.eval.probes` (imports without torch / datasets):
- Build: `lidc_stats(store) -> {"diams", "prior", "holes", "multis"}` (train voting raters; hole / multi-part records of
  the non-simple targets after the size rule, the same order as lab/B3 `census` + `hole_census`);
  `patches(store, split, half=None, n, ...) -> [(group_id, HU patch)]`; `precision_half(group_id) -> "train" | "val"`
  (seeded patient hash, SHA-256 of `probe-precision:<gid>`, low bit);
  `build_precision(out, patches, diams, n=1000, seed=2026, split="val")` (ValueError if a patch comes from the other
  half); `build_topology(out, patches, holes, multis, diams, seed=2027, per_bin=120, n_two=160, n_simple=160)`;
  `fidelity(holes, multis, seed)` (P5.2 KS + bin shares).
- Layout (C-PROBE): `data/eval_views/probe-precision/{images,caches,lidc/val.jsonl[,lidc/train.jsonl],manifest.json}`,
  `probe-topology/` alike. Record = the `dataset.prepare` record (`item_id` `lidc/probe-<name>/<split>/<k>`, `group_id`
  = patient of the background patch, `W = H = S`, `raters ["gt"]`, `targets [gt text]`, `rl false`, `caches {}`,
  `image_sha256`) + `"view": {"offset": [0, 0], "side": S}` + `"probe": {family, diam[, bin, hole_area]}`. Manifest
  `files[rel].sha256` as `prepare` writes it, so `serve.preload` accepts the directory (tested: 0 RL items; single-rater
  items cannot have an `ExpertCache`, k ≥ 2, so there are no caches, exactly as `prepare` does for single-rater
  items). Probe-train (`--n-train 20000`) is a plain `lidc/train.jsonl`, so `sft_rows(root, "lidc", "train")` serves
  the capacity run unchanged.
- Score: `score_precision(rows, root=..., split="val")` → `{n_items, n_samples, legal, assd, iou, centroid, hd95,
  win_share, checks, pass}` on 6–16 px items; every sample counts for legal, medians over legal samples, the paired win
  uses the per-item mean IoU (illegal = 0) vs the equal-area centroid disc. `score_topology(rows, root=...)` →
  `{"bins": [[per-item mean hit] × 3], "hole_recall": {bins, n_bins, pooled [pt, lo, hi], any_of_k, vs_report_line},
  part_recall, controls_correct, false_nonsimple, legal, gate, pass}` (`pass` = rule (a) only, read by
  `gates.phase2_entry`; `hole_recall` report-only). `score_precision(..., capacity=rows)` adds `capacity: {pass, ...}`. Rows of one item (seeds) are pooled; items are in
  record order so per-seed `bins` align (G's `c2_analyze.analyze` input; tested). `false_nonsimple` = mean over the
  160 controls of the per-item share of samples that are legal with a signature ≠ (0,); `controls_correct` = share
  legal with (0,) (gate (a) ≥ 95 %); `legal` = share of legal samples over all items.
  `stratified_ci(per_bin, weights, rng, B)` = `power_c2.pooled` (same RNG consumption).
- Synthetic answer models (tests, sanity): `model_exact`, `model_noise(σ)`, `model_shift`, `model_centroid_disc`,
  `model_centred_prior(prior)`, `model_a0`.

`experiments/phase2/probes.py`: `build [--n-test 1000] [--n-train 0] [--only ...]`, `score --probe ... --completions`.
`experiments/phase2/baselines.py`: `fit_prior`, `score_prior`, `win_share` (sign test via `scipy.stats.binomtest`),
`summary`, `_Sam` (MedSAM official min-max, SAM 2.1-tiny via `transformers.Sam2Model`; lazy torch), jittered-box
MedSAM row; per view (centred = `data/prepared`, shifted, rescaled eval views) and dataset; centred per-item win share
over the prior for MedSAM and SAM 2.1-tiny (G-prior (a) reference). The A4 jitter script is subsumed: the prior on
the shifted view is the prior scored on V's `lidc-val-shifted` directory.

## Design choices and evidence

- **Generator v3e only** (simplicity pass; v3a–v3d dropped). Equivalence to the lab is exact: the real-data test
  rebuilds the hole records from the store and reproduces `topo_probe_v3e_seed{7,8}.json` KS p-values and bin shares to
  the 4 recorded decimals (seed 7 fails, hole area p 0.0037; seed 8 passes). One bug found by this test: a
  `next(x for x in DIVISORS if x >= f(rng.uniform()))` drew the margin once per divisor; fixed (one draw).
- Ground truth stored as the certified target text, not a raster: B3 P2.1 measured 100 % pixel-exact texts; the
  writer asserts it per item.
- Prior fit from the prepared train records' certified targets (no store pass) vs A4's raw masks every 10th item:
  chosen for one data path through the view layout; differs only by the size rule.
- Background patches larger than the 125-px patch (S ≥ 200, rare) are reflect-padded (B3's renderer would crash).
- Precision family mix kept at B3's (fourier 60 / lobulated 30 / hole 5 / two 5) because the calibration was measured
  on it; topology is scored only on the topology probe.

## Tests (`tests/test_p2_probes.py`, 29: 27 CI + 2 `data`)

Includes `phase2_entry` on a minimal report: σ 0.5 passes G-precision, the disc fails it, disc + exact capacity
passes with the capacity warning; G-topology-a reads `pass`, G-topology-b reports `hole_recall`.

Layout/manifest; exactness (ASSD = HD95 = 0, IoU 1 on GT); surface = brute force (1e-9); calibration (disc IoU within
0.806 ± 0.03, synthetic log-normal diameters); discrimination ×8 (exact / σ 0.25 / σ 0.5 pass; σ 1 / σ 2 / disc /
centred prior / 1-px shift fail; disc fails only via the win rule, win share 0); 1-px shift ASSD > 0, centroid 1;
topology composition (bins, size rule); sanity ×3 (exact 1.0, A0 0 in every bin, disc 0; controls 1.0); σ 0.5 in (0, 1);
stratified CI = `power_c2.pooled`; patient halves disjoint + build refuses the wrong half; `c2_analyze.analyze` on 3
seeds; `serve.preload` accepts the directory; prior / win share on a prepared-layout fixture; `--help` ×2; data: v3e
fidelity seeds 7 / 8, built directories patient-disjoint. Runtime ~10 s CI part, ~70 s with data.

## Real build (`.venv/bin/python -I -m experiments.phase2.probes build`)

117 s wall, 85 s CPU, peak 0.6 GB, 1 process. `lidc_stats` 50 s: 23,459 train diameters, 285 holes, 182 multi-part
targets (= B3 A1; bin shares 0.4386 / 0.3614 / 0.2000). probe-precision `lidc/val.jsonl`: 1,000 items (597 at
6–16 px), 71 patches, 155 MB; probe-train not built (`--n-train 20000` on demand, ~3 GB est.). probe-topology: 680
items (360 hole = 120 / 120 / 120, 160 two, 160 simple), 19 val patients, 115 MB; build-seed (2027) fidelity
passes (lowest KS p 0.018, hole area; shares 39 / 39 / 22). Sanity on the real directories (B3 values in brackets):
centroid-disc median IoU 0.823 [0.822]; precision pass for exact and σ 0.5 (legal 96.6 %, IoU 0.891, win 78 %
[78.5 %]), fail for σ 1 (win 28 %), disc (win 0, every other check passes), 1-px shift; topology exact 1 / 1 / 1, A0
0 / 0 / 0, disc 0 / 0 / 0, σ 0.5 mean over 8 0.64 / 0.83 / 0.85, pooled 0.749 [0.723, 0.774] [0.76 [0.73, 0.78]],
0 false non-simple, controls correct 97.3 %.

## Line counts

Lab sources promoted: probe_lib 293 + precision_probe 313 + topo_probe_v3 274 + probe_rule 54 + power_c2 119 +
hole_census 133 + A4 probe_proto 141 = 1,327; A4 null_prior 65 + null_prior_jitter 51 + medsam_check 121 = 237.
Shipped: `src/objdel/eval/probes.py` 664 (≈ 520 before `ruff format` split long calls), `experiments/phase2/probes.py`
76, `baselines.py` 213, tests 289. Above the ~200 target: the module carries two generators, metrics, writer and two
scorers that the plan puts in one file.

## Open issues

- Precision probe-test patches: only 71 lung-parenchyma patches pass B3's background criteria on the "val" half of
  the train patients at ≤ 3 per patient (B3 used 200 val patches); items reuse them at random.
- MedSAM / SAM 2.1-tiny paths are untested here (no downloads by rule); `Sam2Processor.post_process_masks` signature
  should be checked on first run (transformers 5.2).
- Prompts are not written into the probe manifests (`PROMPTS` lives in `objdel.train.dataset`, which needs
  `datasets`); R's export uses `PROMPTS["lidc"]` as for real items.
- `score_precision` ignores `truncated`; a truncated sample is simply illegal/unparsable.

## Fix round (FX-EGS)

- `stratified_ci` is the one pooled-bootstrap implementation (`c2_analyze.pooled_ci` calls it); it rounds point / lo /
  hi to 12 decimals so an exact 288 / 960 = 0.30 meets the line (BK2-4). power_c2 equality tests kept (atol 1e-12).
- Lines: `probes.py` 666 → 669.
