# Phase 2 execution: Stage 1 (build on the Mac)

Stage 1 of `PLAN.md` §7: milestones M1–M5, M7, M8, M9, M11 (M6 removed; M10 after V1). No GPU, no Vast, no git
commits. Owner time box: 1–2 hours, so every agent works in parallel against the contracts below instead of
waiting for upstream milestones. Plan sections referenced: §3 (architecture), §5 (gates), §6.0 (bundle), §7 (tests).

## 1. Rules for every agent

1. **Think, test, implement.** For every design choice the plan does not fix, try at least two options in your
   sandbox (`<scratchpad>/phase2/impl/<agent>/`), measure them (speed, line count, correctness on real data where
   possible), pick one, and record the evidence in your impl note. No proposal goes into `src/` without a test
   that shows it works.
2. **TDD.** Write the tests first (from PLAN §7's "Tests" column), run them and see them fail, then implement.
3. **Simplicity.** After the tests pass, do one simplification pass: could 200 lines be 50? Delete dead options,
   merge helpers, prefer the standard library / numpy / TRL / transformers over custom code. Report line counts
   before and after in your impl note. Target: a new module ≤ ~200 lines unless the plan forces more.
4. **Own your files only** (§3 table). Need a change in a file you do not own? Ask the coordinator in your final
   message, or code against the contract and say what the owner must provide. Never edit `tests/conftest.py`,
   `pyproject.toml`, `Makefile`, `src/objdel/*/__init__.py` (except the owner listed), `experiments/phase1/`
   (except the M9 files listed), `PLAN.md`, `ROADMAP.md`.
5. **Resources.** Mac: 8 cores, 16 GB shared by ~8 agents. At most 2 processes and ~2 GB per agent. **Run only your
   own test files** (`.venv/bin/python -m pytest -q tests/test_p2_<you>*.py`); the coordinator runs the full
   suite. No MPS, no model downloads (use the tiny offline model in `tests/conftest.py` or a tiny random config).
   Wait loops must not match their own command line (wait on a PID, or `pgrep -f "[p]attern"`).
6. **CI.** CI runs `pytest -m "not slow"` with only `[dev]` installed (no torch, transformers, trl, datasets,
   pydicom) and without `data/`. Tests that need them: `pytest.importorskip(...)` at the top of the test (or the
   module), and the `data` marker plus a skip when `data/` is absent. Library modules under `src/objdel/eval`,
   `targets`, `reward`, `checker`, `geometry` must import without torch.
7. **Data** in `data/` is untrusted: Python that reads it runs with `.venv/bin/python -I`. Real-data checks go in
   tests marked `data` (skipped in CI) or in your impl note's evidence.
8. **Lint**: `.venv/bin/python -m ruff check src tests experiments` and `ruff format` on your files.
9. **Report**: edit your impl note `experiments/phase2/impl/<M>.md` (already created; use Edit, not Write):
   design choices with evidence, API as built, tests, line counts, open issues. Your final message: ≤ 25 lines.

## 2. Dependency graph and task graph

```
 wave 1 (parallel, contracts §4 break the waits)                          wave 2              wave 3      wave 4
 ┌──────────────────────────────────────────────┐
 │ V  M1 views: frame(offset, side, window),    │──┐ C-VIEW
 │    augment, build_views, build_eval_views,   │  │
 │    M9 fetch_lidc_pixels fixes                │  │
 ├──────────────────────────────────────────────┤  ▼
 │ T  M2 + M5 + M9(dataset): view_rows, configs │──┐ C-ROWS, C-BRIDGE
 │    (qwen3.5, LM-only lora, sampling, WSD),   │  │
 │    sample export/import, sample_vllm.py      │  │
 ├──────────────────────────────────────────────┤  │
 │ S  M3 serve.py + k3_harness                  │──┤ C-SERVE          ┌─────────────┐  ┌─────────┐  ┌────────┐
 ├──────────────────────────────────────────────┤  │                  │ 4 break-it  │  │ fix     │  │ final  │
 │ E  M4a metrics (HM-IoU, NSD, stats), vision, │──┤ C-EVAL     ────► │ agents:     │─►│ fleet   │─►│ verify │
 │    gates (revised phase2_entry), M9 eval_run │  │                  │ wiring +    │  │ (owners │  │ make   │
 ├──────────────────────────────────────────────┤  │                  │ edge cases  │  │ of the  │  │ test / │
 │ P  M4b probes + experiments probes/baselines │──┤ C-PROBE          │ → ISSUES.md │  │ issues) │  │ cov    │
 ├──────────────────────────────────────────────┤  │                  └─────────────┘  └─────────┘  └────────┘
 │ G  M11 signal + prereg.py + c2_analyze       │──┤ C-SIGNAL
 ├──────────────────────────────────────────────┤  │
 │ O  M7 Vast ops: scripts, policy.py, Docker   │  │ (independent)
 ├──────────────────────────────────────────────┤  ▼
 │ R  M8 runners: train_sft, select_ckpt,       │◄─┘ integrates every contract; CPU smoke train_sft → select_ckpt
 │    g_probe, run configs, CPU smoke           │    → g_probe on the tiny model (stubs until owners land)
 └──────────────────────────────────────────────┘
```

Plan order (M2 after M1, M5 after M2, M11 after M4, M8 after all) is kept **at the interface level**: downstream
agents code against §4 now and switch from stubs to the real functions when the owner lands them. R runs the
integration smoke last.

## 3. File ownership

| Agent | Milestones | Owns (create or edit) | Impl note |
|---|---|---|---|
| V | M1, M9 (LIDC fetch) | `src/objdel/targets/{frame,build,augment}.py`, `src/objdel/targets/__init__.py`, `experiments/phase2/{build_views,build_eval_views}.py`, `experiments/phase1/fetch_lidc_pixels.py`, `tests/test_p2_views.py` | `impl/M1.md` |
| T | M2, M5, M9 (dataset) | `src/objdel/train/{configs,dataset,sample}.py`, `experiments/phase2/sample_vllm.py`, `tests/test_p2_train.py`, `tests/test_p2_bridge.py` (existing `tests/test_train_*.py` only where your change requires) | `impl/M2.md` |
| S | M3 | `src/objdel/train/serve.py`, `src/objdel/train/rewards.py` (pool integration only), `experiments/phase2/k3_harness.py`, `tests/test_p2_serve.py` | `impl/M3.md` |
| E | M4a, M9 (eval) | `src/objdel/eval/{metrics,vision,gates}.py`, `src/objdel/eval/__init__.py`, `experiments/phase1/eval_run.py`, `tests/test_p2_eval.py` (existing `tests/test_eval*.py` where the gate change requires) | `impl/M4.md` |
| P | M4b | `src/objdel/eval/probes.py`, `experiments/phase2/{probes,baselines}.py`, `tests/test_p2_probes.py` | `impl/M4b.md` |
| G | M11, M8 (stats) | `src/objdel/eval/signal.py`, `experiments/phase2/{prereg,c2_analyze}.py`, `tests/test_p2_signal.py` | `impl/M11.md` |
| O | M7 | `experiments/phase2/vast/**`, `experiments/phase2/docker/{Dockerfile.phase2,Dockerfile.phase3,build_images.sh,image_check.sh}`, `tests/test_p2_vast.py` | `impl/M7.md` |
| R | M8 (runners) | `experiments/phase2/{train_sft,select_ckpt,g_probe}.py`, `experiments/phase2/configs/*.yaml`, `tests/test_p2_scripts.py` | `impl/M8.md` |

`experiments/phase2/__init__.py` exists (coordinator). Each script has a `--help` test in its owner's test file.

## 4. Contracts (fixed; an owner may extend, never break; changes go through the coordinator)

**C-VIEW (V → T, R, P).**
- `objdel.targets.frame(item, cfg=None, *, offset=(0, 0), side=None, window=None) -> Frame`. Defaults reproduce
  today's frame byte for byte. `offset=(dx, dy)` moves the crop window in item pixels; `side=S` forces the frame
  side (a divisor of 1000, ≤ 500; ValueError if the rater union does not fit with ≥ 2 px clearance);
  `window=(lo, hi)` is the HU window for LIDC (default `LUNG_WINDOW = (-1000, 400)`), integers.
- `objdel.targets.augment`: `OPS = ("id", "r90", "r180", "r270", "fh", "fv", "tr", "at")` (D4; RIGA uses
  `("id", "fh")`); `apply(a, op)` on the last two axes; `render(spec: dict, context) -> np.ndarray` uint8
  (side_img, side_img, 3), equal to `frame(...)` with the spec's offset / side / window, the op applied at frame
  resolution, then the Pillow upscale; `load_context(view_dir, item_id)`.
- Training view directory `data/views/<view>/`: `manifest.json` (view name, config, source manifests' SHA-256,
  counts, `view_hash` = SHA-256 of `rows.jsonl`), `rows.jsonl` sorted by (epoch, item_id, rater_index, rep), one
  line per training row: `{"row_id", "epoch", "item_id", "dataset", "rater", "rater_index", "op", "offset": [dx,
  dy], "side", "window": [lo, hi] | null, "noise_sigma", "noise_seed", "simple": bool, "rep", "target"}` (only
  voting raters; non-simple rows repeated `nonsimple_repeat` times, each repeat with its own draw), and
  `context/<item_id>.npz` (`pixels`: int16 HU for LIDC, uint8 RGB for RIGA; `x0`, `y0` in item pixels).
- Eval / probe view directories `data/eval_views/<name>/`: **the `data/prepared` layout** (`images/`, `caches/`,
  `<dataset>/<split>.jsonl`, `manifest.json`), so `records`, `cache_loader`, `grpo_rows` and `eval_run` work
  unchanged; each record adds `"view": {"offset": [dx, dy], "side": S}`. Names: `lidc-val-shifted`,
  `lidc-val-rescaled` (RIGA alike when permitted).

**C-ROWS (T → R).** `objdel.train.dataset.view_rows(view_dir, tokenizer=None, max_tokens=None, a0=False) ->
Dataset`: columns `messages` (as `sft_rows`), `image` (rendered lazily via `augment.render`, PIL), `key`,
`rater`, `row_id`, `length` (tokens of the target). `a0=True` empties every `"holes"` list and keeps parts.
`configs`: `MODELS["qwen3.5-0.8b"]`, `MAX_COMPLETION` (1,024 for LFM and Qwen3.5), `sampling(model) -> dict`
(T 1, top_p 1, top_k 0, min_p 0, repetition 1, presence 0), `lora(model)` language-model only,
`param_groups(model, base_lr)` (vision tower + projector at 0.1×), `sft_config(model, out, **o)` with WSD.

**C-BRIDGE (T → R).** `objdel.train.sample.export_prompts(rows, path, model, dataset, n)` → JSONL
`{"id", "prompt", "image", "processor_kwargs", "sampling", "n"}`; `experiments/phase2/sample_vllm.py --prompts P
--model DIR --out C` (runs in `/opt/vllm`, imports nothing from objdel; refuses to start when CUDA free memory
< 90 %) → JSONL `{"id", "k", "text", "token_ids", "finish_reason"}`; `import_completions(path, stop_ids)` → eval
JSONL rows `{"item_id", "samples", "truncated"}` (truncation by `rewards.truncated`).

**C-SERVE (S → R).** `objdel.train.serve`: `preload(root, datasets) -> int` (fills the module store; LIDC-only
loads no RIGA cache), `CacheRef(key, root)` (Mapping, pickles as key and root), `loader(root)`, `fork_pool(n)` (fork
context, warmed, `gc.freeze()` first; asserts `torch.cuda.is_initialized()` is False when torch is imported).
`TRLReward(caches=serve.loader(root), ..., pool=...)` is the scoring path.

**C-EVAL (E → R, G).** `metrics(...)` gains HM-IoU16 and NSD@τ; `summarize` gains truncation and token stats.
`objdel.eval.vision.check(model, processor, images, items) -> {"ratio", "tv", "pass", ...}` (G-vision, §5.1;
ratio vs the base reference). `gates.phase2_entry(report, datasets)` revised per §5.1; E documents the `report`
schema in its module docstring (views centred / shifted / rescaled, prior, probes, vision, signal, g_probe).

**C-PROBE (P → R).** `objdel.eval.probes`: build the precision and topology probe items **as eval-view
directories** (C-VIEW eval layout, `data/eval_views/probe-precision`, `probe-topology`) and score them:
`score_precision(rows) -> dict`, `score_topology(rows) -> dict` (pass rules of §5.1 / D20).

**C-SIGNAL (G → R).** `objdel.eval.signal`: `panel(groups) -> dict` (per-group metrics, the variance split,
sd_geom), `select_rb(branches) -> (name, reasons)` (D12 rule R-B), `spearman_ci(x, y)`. `experiments/phase2/
prereg.py`: every pre-registered rule as data + pure functions (C1 stop rule, R-B, G-probe, C2 outcomes, §6.3).

## 5. Status

| Agent | State | Notes |
|---|---|---|
| V, T, S, E, P, G, O, R | wave 1 done 2026-10-20 | notes `impl/M*.md`; full fast suite 994 passed, 3 failed (S1-14, S1-15) |
| BK1–BK4 | wave 2 (break-it) done: 55 findings, 11 high | data path, scoring path, runners, ops + CI simulation; findings `impl/BK*.md` → `impl/ISSUES.md` |
| fix fleet (6) | wave 3 done | all findings fixed test-first; open items in `impl/ISSUES.md` (status section) |
