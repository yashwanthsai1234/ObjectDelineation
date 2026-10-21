# BK2 break-it findings

Scope: the scoring and decision path, JSONL bridge (`sample.export_prompts` / `import_completions`, `sample_vllm.py`)
→ `serve` (preload, CacheRef, fork_pool) → `eval.metrics` / `summarize` → `eval.probes` scorers → `eval.signal`
(panel, `select_rb`) → `prereg.py`, `c2_analyze.py` → `gates.phase2_entry` → `select_ckpt.py`, `phase1/eval_run.py`.

Tests: `tests/test_p2_breakit_scoring.py`, 38 tests: **24 passing guards, 14 strict xfails** (12 defects, BK2-1 to
BK2-12). Command: `.venv/bin/python -I -m pytest tests/test_p2_breakit_scoring.py` (7 s; the 5 `data` tests read
`data/eval_views`, skip without it). ruff-clean. No `src/` file or existing test or script was edited.

Owners: T = train bridge (sample), S = serve, E = eval (metrics, gates, eval_run), P = probes, G = signal / prereg /
c2_analyze / select_ckpt (M8 / M11).

## Findings

| id | sev | component | repro (test) | expected vs actual | suggested fix | owner |
|---|---|---|---|---|---|---|
| BK2-1 | **high** | select_ckpt `--probes` → `dataset.grpo_rows` | `test_select_ckpt_probe_rows_are_the_probe_items` | `sample_view(a, best, a.probes / "probe-…")` (`select_ckpt.py:210`) builds prompts with `grpo_rows`, which keeps only `rl` records (`dataset.py:220`). Every probe record has `rl: false` (M4b: single-rater, `caches {}`), so 0 of 680 topology / 0 of 1,000 precision items are sampled. `score_precision([])` and `score_topology([])` then **return `pass: False` without an error** (legal / controls = NaN). Result: G-precision and G-topology-a FAIL on every R1 branch whatever the model does; C1 cannot pass | a prompt-row builder that ignores `rl` (e.g. `grpo_rows(..., rl_only=False)` or `prompt_rows(root, ds, split)`) for probe dirs; make both scorers raise on 0 scored items or an empty bin / control set | G (+ T for the builder, P for the empty-input guard) |
| BK2-2 | **high** (HF / Qwen3.5 path) | select_ckpt `--items` on probes | `test_topology_probe_on_select_ckpt_items_subset_can_pass` | `--items N` (`select_ckpt.py:119-120`, HF path documented as `--items 256`) also cuts the probe views to their first N records. Topology probe record order is bin 0 / bin 1 / bin 2 / two / simple: the first 256 hold 0 simple controls and 16 bin-2 holes, so a **perfect model** gets `controls_correct = NaN` → G-topology-a FAIL; precision keeps 152 of the 1,000 items | never subset probe views (the gate is defined on the full probe), or subset by a seeded stratified draw; same root cause as S1-12's "first N, not a seeded subset" but here it flips a gate | G |
| BK2-3 | medium | `sample.import_completions` | `test_bridge_short_item_is_refused`, `test_bridge_duplicate_id_k_is_refused` | Only checks that k is 0..m-1 per id. A completions file cut at a line boundary (an item with 7 of 8 samples) is accepted: metrics then use one group of 7 (`group = K`), G-signal / eligibility silently change. A duplicated `(id, k)` line silently overwrites the first. Items missing entirely are dropped (eval_run scores fewer items; select_ckpt `score_view` KeyErrors at `select_ckpt.py:138`) | `import_completions(path, stop_ids, n=None, ids=None)`: ValueError on a duplicate `(id, k)`, on a count ≠ n, and on ids ≠ the exported set (select_ckpt has both from `export_prompts`) | T |
| BK2-4 | medium | `c2_analyze.pooled_ci` / `prereg.c2_outcome` | `test_c2_pooled_exactly_on_the_floor_is_not_a_fail`, `test_c2_pooled_exactly_on_the_pass_line_counts_as_050` | Bin totals are integers over 960 samples, so an exact pooled 0.30 / 0.50 is reachable: bins 349 / 349 / 44 hits give exactly 288/960 = 0.30 but the float sum is 0.29999999999999993 → **outcome "fail" (R-ctrl) instead of "partial" (R-x8)**; bins 597 / 597 / 12 give 0.49999999999999994 → clause (1) and the seed clause miss the pass line. 8 (0.30) and 31 (0.50) such bin triples exist. Same class as S1-15 | compute the point estimates from integer hit counts with `fractions.Fraction` weights, or compare against the lines with a 1e-12 tolerance (`round(x, 12)`) in `c2_clauses` / `c2_seed_clause` / `c2_outcome`; same for the per-bin clause (3) | G |
| BK2-5 | medium | select_ckpt eligibility | `test_select_ckpt_eligibility_needs_a_vision_result` | `eligible(..., res["vision"].get("pass", True))` (`select_ckpt.py:189`) counts an unmeasured G-vision as passed; `--vision` defaults to `{}` = off. prereg `RB["eligible"]["g_vision"] = True` and `gates.phase2_entry` fail a missing vision, so R-B can select a branch that fails G-vision, and selection ≠ the gate table | `vision.get("pass") is True`; make `--vision` required (or refuse to select without it) | G |
| BK2-6 | low | select_ckpt eligibility | `test_select_ckpt_eligibility_needs_groups_of_8` | `eligible` reads `dup_rate` / `signal_share` without checking `summary["group"] == 8`; gates refuse G-signal unless groups are 8 (`gates.py:95`). `--k 12` → eligible branch, failed gate | check `group == REWARD_V4.group`, or compute eligibility from `gates.phase2_entry(report, ["lidc"])["passed"]` (G-format, G-signal, G-vision): one rule, no copy (also fixes BK2-5) | G |
| BK2-7 | low | `signal._pair` R-B veto margin | `test_rb_veto_fires_at_exactly_15_percent` | sd_geom 0.30 → 0.345 (exactly +15 %) gives `rel = 0.1499999…` → veto does not fire; 38 of the 99 two-decimal bases fail the same way (`signal.py:161-162`) | `diff[0] >= VETO_REL * base - 1e-12` (or round rel to 12 digits) | G |
| BK2-8 | low | `signal.select_rb` | `test_rb_refuses_duplicate_depths` | Two eligible branches of one depth are accepted; their order (and so the pick) is the input order. select_ckpt derives depth from the directory name (`float(b.name[1:])`), so `e1` and `e1.0` would collide | ValueError on duplicate depth or name | G |
| BK2-9 | low | `metrics._ged` | `test_ged_is_never_negative` | 8 samples = each of 4 raters twice (GED exactly 0) gives −1.1e-16 (`metrics.py:226`); report-only, but the same no-clamp pattern as S1-15 (`_nsd`, `metrics.py:197`) | `max(0.0, …)` in `_ged`; clamp `_nsd` to [0, 1] (S1-15) | E |
| BK2-10 | medium | `gates.phase2_entry` | `test_empty_report_passes_no_c1_gate` | With no judged dataset (datasets=None and no shifted summary: e.g. `eval_run --view centred / rescaled`, or a report missing its shifted view) the overall `pass` is False but **every per-gate entry of `passed` is True** (gates start `pass: True`, `gates.py:129`, and nothing ANDs into them). `passed` is `prereg.c1_stop`'s input, so with G-vision passing `c1_stop` accepts the rung; eval_run's report.md prints PASS for G-format … G-signal | start each dataset gate at `bool(judged)` (or set them False when `rows` is empty) | E |
| BK2-11 | low | `serve.preload` after `fork_pool` | `test_preload_after_fork_is_refused` (guard `test_pool_worker_after_a_later_preload_fails_loudly`) | A preload after the fork is accepted; the parent's `CacheRef` resolves, the worker raises KeyError (loud for a new root), but a **re-preload of the same root** with rebuilt caches is silently scored on the old ones in the workers | module flag set by `fork_pool`; `preload` raises RuntimeError while a pool forked from this STORE is alive | S |
| BK2-12 | medium | prereg C1 ↔ gates (G-prior (a) slot) | `test_prior_a_threshold_is_applied_when_set` | Setting `G_PRIOR_A_THRESHOLD` (the plan's owner option before V1) adds "G-prior-a" to `c1_gates`, but `phase2_entry` never writes `passed["G-prior-a"]` and no code compares `win_share` with the threshold: `c1_stop` then fails **every** rung (`gates.get("G-prior-a", False)`) | `phase2_entry(report, datasets, prior_a_threshold=None)` adds a G-prior-a gate (every judged structure's centred win share ≥ threshold) when set; prereg passes its constant | E / G |

## Guards that hold (passing tests)

Bridge: unicode / empty / `finish_reason "length"` round trip (truncation = TRL's id rule, empty completion
truncated, tokens = id count); a gap in k and a cut last line are loud; `sample_vllm.outputs` refuses a missing vLLM
result (strict zip). Serve: CacheRef for a key / root not preloaded and for a key in two roots raise KeyError;
`loader(root)` disambiguates; preload refuses an edited jsonl, a stray jsonl, and another item's cache; pool scores
equal serial scores byte for byte on 6 real shifted-view items (select_ckpt `_score`, including a truncated garbage
sample). Metrics: an all-failed group is finite and `allow_nan=False`-serialisable through summarize; an identical
group has dup 7/8, no σ, zero advantages; K = 12 is one group and gates then refuse G-signal; mixed token presence →
`tokens: None`; TV / JS / minimum-mode retention stay in [0, 1] over 3,000 random histograms. Probes: an exact model
passes the full topology probe (controls 1.0, bins 120 / 120 / 120). R-B: single eligible branch, no eligible branch,
veto at +16 % with an IoU CI touching 0 (fires: not a significant loss), CI lower bound exactly 0 (does not fire).
Prereg / C2: seed clause, fourth-seed clause and outcome table at 0.30 / 0.50 exactly; `c2_analyze` reads
`score_topology`'s keys. Gates: gate names equal `prereg.C1_GATES` (and `passed` = C1_GATES + G-probe, `reported` ⊇
REPORT_ONLY); prereg's eligibility constants equal the gate constants; G-prior (a) at win share 0 does not block a
pass; NaN values fail their gates; a model dropped from LIDC is judged on RIGA alone and no LIDC value is read.

## Simplicity audit (this path)

1. `select_ckpt.py:96-106` `eligible` re-implements G-format + G-signal from `gates._checks` with its own copy of
   the thresholds (via `prereg.RB`); call `gates.phase2_entry({"views": {"shifted": …}, "vision": …}, ["lidc"])` and
   read `passed` for G-format / G-signal / G-vision (−10 lines, fixes BK2-5 / BK2-6).
2. `prereg.py:51` `RB["eligible"]` repeats `gates.py:50-52` literals (0.98 / 0.90 / 0.10 / 0.50); reference
   `gates.FORMAT` etc. (the BK2 guard test currently pins their equality).
3. `c2_analyze.py:25-32` `pooled_ci` is `probes.stratified_ci` (`probes.py:569-574`) line for line (both "B3
   power_c2.pooled"); keep one (−8 lines). `prereg.py:98` `C2["bin_weights"]` duplicates `probes.py:58`
   `BIN_WEIGHTS`; import one from the other. `probes.TOPOLOGY_RULE["recall_report"]` (0.30) duplicates
   `C2["floor"]`.
4. `select_ckpt.py:58-83` `sampler_cmd`: 26 lines of one-token-per-line lists (black style, unlike the rest of the
   repo) → one shared tail `["--prompts", str(prompts), "--model", str(model), "--out", str(out)]` and two heads
   (~6 lines).
5. `metrics.py:160-163` and `metrics.py:208-212` both compute the joint bounding-box crop; one `_crop(masks)` helper.
6. `serve.py:101-103` `ref_loader(key)` is `CacheRef(key)`; `ref_loader = CacheRef` (or document `CacheRef` as the
   loader) removes a wrapper.
7. `probes.score_precision` / `score_topology` (`probes.py:530`, `577`) silently produce NaN on empty input and
   `pass: False`; a 2-line `if not rows: raise ValueError` is simpler than callers reasoning about NaN, and would have
   surfaced BK2-1 at once.
8. `sample.import_completions` (`sample.py:137-153`) builds `{id: {k: rec}}`, then re-sorts; reading `(id, k)` into
   one dict with a duplicate check and comparing counts against the exported `n` is the same size and closes BK2-3.

## Notes (no test, informational)

- `select_ckpt.hf_sample` seeds every prompt with the same `torch.manual_seed(seed)` (`select_ckpt.py:276`): sample
  streams are correlated across items (vLLM gets one seed per request too, so parity holds; fine for selection).
- `score_precision` counts an item twice in `win_share` if rows repeat an item (e.g. two seeds concatenated); its
  docstring does not say rows must be one per item.
