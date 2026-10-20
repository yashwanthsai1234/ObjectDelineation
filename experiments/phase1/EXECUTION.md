# Phase 1 execution: task graph, dependencies, contracts and agent protocol

This file coordinates the agents that implement [PLAN.md](PLAN.md). Read it fully before touching code.
Milestone numbers (M1–M9) and gates (K*, R*) are the plan's.

## 1. Agent protocol (every task, no exceptions)

Follows the Karpathy coding guidelines (github.com/forrestchang/andrej-karpathy-skills): think before coding,
simplicity first, surgical changes, goal-driven execution with verifiable success criteria.

1. **Think.** Read this file, the plan section for your milestone, the lab code and the lab report it promotes.
   Write down assumptions and the success criteria (the gate) before writing code.
2. **Research online (mandatory).** Before each non-trivial choice, search the web (load `WebSearch`/`WebFetch`
   with ToolSearch `select:WebSearch,WebFetch`). Look for a better algorithm, library or API than the lab used,
   and for known pitfalls (version changes, edge cases). Record every source you rely on.
3. **Sandbox-test the alternatives.** Prototype in your scratch dir (§6), never in the repo. A proposal enters the
   repo only with evidence: a test or measurement that shows it is correct and at least as good as the lab version.
4. **Test first.** In the repo, write the failing tests (unit, property with Hypothesis, oracle/differential,
   data-backed with `@pytest.mark.data`), run them red, then write the minimum code that makes them green.
5. **Simplify.** After green, re-read your code: could 200 lines be 50? Remove dead code, unused parameters,
   speculative options and copies of lab code. One production implementation per operation; lab prototypes may
   survive only as slow oracles under `tests/oracles/`.
6. **Cite.** Every module docstring has a `References:` list of URLs (papers, docs, issues) for the algorithm or
   API it implements, and non-obvious lines cite the specific source in a comment. The same sources go in the
   implementation note.
7. **Document.** Write `experiments/phase1/impl/<task>.md`: what was researched (with links), alternatives tried
   in the sandbox and their numbers, the decision, the tests, the gate result, open issues.
8. **Verify.** Your own tests, then `make lint` and the full `make test`. A failure in a file you do not own is
   reported in your final message, not fixed.

**No git commands that change state** (no commit, add, stash, checkout, reset). Read-only `git diff`/`status` is fine.

## 2. Task graph

```
 wave 1 (parallel)          wave 2                 wave 3              wave 4      wave 5–7
 ┌────────────────┐
 │ M1 polycore +  │──┬───────────────┐
 │ build + CI     │  │               ▼
 └────────────────┘  │        ┌──────────────┐
 ┌────────────────┐  │        │ M5 targets   │──────┐
 │ M2 geometry    │◄─┘(to_mask)└──────────────┘      │
 │ spec, parser   │──────────┐        ▲             ▼
 └────────────────┘          │        │      ┌──────────────┐
 ┌────────────────┐          ▼        │      │ M7 TRL       │──┐
 │ M3 topology +  │───┬──►┌──────────────┐   │ integration  │  │   ┌─────────┐   ┌──────────┐
 │ size rule      │   │   │ M6 reward v4 │──►│ + smoke run  │  ├──►│ M9 x86  │──►│ break-it │
 └────────────────┘   │   └──────────────┘   └──────────────┘  │   │ Docker +│   │ → fix →  │
 ┌────────────────┐   │          │                  ▲          │   │ RESULTS │   │ simplify │
 │ M4 data: RIGA, │───┼──────────┼──────────────────┘          │   └─────────┘   │ → verify │
 │ LIDC, splits,  │   │          ▼                             │                 │ → A/B    │
 │ store          │   │   ┌──────────────┐                     │                 │ → docs   │
 └────────────────┘   └──►│ M8 eval +    │─────────────────────┘                 └──────────┘
                          │ hack audit   │
                          └──────────────┘
```

| Task | Depends on | Starts when | Owner paths (only the owner edits these) |
|---|---|---|---|
| M1 | — | now | `cpp/polycore/`, `src/objdel/_polycore*`, `src/objdel/checker/polycore.py`, `CMakeLists.txt` (root), `pyproject.toml` `[build-system]`/`[tool.scikit-build]`, `Makefile`, `.github/workflows/`, `tests/test_polycore*.py`, `tests/oracles/raster_oracle.py`, `tests/oracles/validity_oracle.py` |
| M2 | M1 (only for `to_mask`) | now | `src/objdel/geometry/`, `tests/test_geometry*.py`; deletes `src/objdel/checker/polygon_io.py`, `tests/test_polygon_io.py` |
| M3 | — | now | `src/objdel/checker/topology.py`, `src/objdel/checker/__init__.py`, `tests/test_topology*.py`, `experiments/phase1/census_sizerule.py` |
| M4 | — | now | `src/objdel/data/`, `tests/test_data*.py`, `experiments/phase1/build_data.py`, `data/store/` (output) |
| M5 | M1, M3 | M1 and M3 done | `src/objdel/targets/`, `tests/test_targets*.py`, `experiments/phase1/build_targets.py` |
| M6 | M2, M3 (M4 for store-backed runs; fixtures suffice for tests) | M2 and M3 done | `src/objdel/reward/`, `tests/test_reward*.py`, `tests/oracles/reward_*.py`, `experiments/phase1/reward_gates.py`; deletes `reward/annotators.py` |
| M7 | M2, M4, M5, M6 | all four done | `src/objdel/train/`, `tests/test_train*.py`, `experiments/phase1/smoke_*.py` |
| M8 | M6 | M6 done | `src/objdel/eval/`, `tests/test_eval*.py`, `experiments/phase1/eval_*.py` |
| M9 | M1–M8 | all done | `experiments/phase1/bench/`, `experiments/phase1/RESULTS.md` (x86-64 tests in Docker; no Vast rental in Phase 1 — x86 timing is a Phase 2 entry gate) |

Shared files (`README.md`, `CLAUDE.md`, `ROADMAP.md`, `src/objdel/__init__.py`, `tests/conftest.py`,
`pyproject.toml` dependencies and pytest/ruff sections) belong to the coordinator. If you need a dependency or
marker, say so in your final message; the pinned stack (below) is already installed in `.venv`.

**Change to the plan:** TRL chat-row export moves from M4 to M7, because SFT rows need M5's target text and
GRPO rows need M6's caches. M4 delivers items; M5 delivers the model frame (crop/scale) and targets.

## 3. Interface contracts

Names and types below are binding between tasks. An owner may improve a contract only with evidence, and must
update this section and tell the coordinator in the final message.

```python
# M1  objdel.checker.polycore   (wraps the nanobind extension objdel._polycore; see impl/M1_polycore.md)
DEFECTS: tuple[str, ...]                       # 13 defect classes, fixed order
MAX_RINGS = 10; MAX_RING_VERTICES = 1000; MAX_DIM = 4096   # beyond the caps: ValueError
pack(rings: list[array (k,2) int]) -> tuple[xy int32 (N,2), off int32 (R+1,)]   # shape = rings list or (xy, off)
    # a ring that is not (k, 2) (or empty) raises ValueError; non-integers TypeError (F1, B1-2)
rasterize(shape, W, H, G=1000) -> {"evenodd": bool (H,W), "nonzero": bool (H,W), "ring_wmin"/"ring_wmax": int32 (R,),
                                   "area_evenodd", "area_nonzero", "px_multi", "rings_mixed_sign": int}
validate(shape, roles=None, G=1000) -> {"total": int64 (13,), "per_ring": int32 (R,13), "area2", "n_clean",
                                        "pair_cross", "pair_touch"}   # roles: -1 outer, p = hole of outer p
signed_area2(ring) -> int                      # exact; > 0 = clockwise on screen (outer)
score_batch(shapes, Ws, Hs, roles=None, gts=None, G=1000, threads=0) -> dict of (n,) arrays
    # area_evenodd, area_nonzero, px_multi, rings_mixed_sign, gt_area, inter_evenodd, inter_nonzero; defects (n,13)
    # an exception in a worker thread (std::bad_alloc -> MemoryError) is raised in Python after all join (F1, I8)

# M2  objdel.geometry
GRID = 1000
SCHEMAS = {"riga": ("disc", "cup"), "lidc": ("nodule",)}
Part = tuple[np.ndarray, list[np.ndarray]]     # (outer (k,2) int64, holes)
Objects = dict[str, list[Part]]
serialize(objs: Objects, names) -> str         # compact grammar of lab/t_stack/prod.py
parse(text, names, truncated=False) -> Parsed  # Parsed(ok: bool, objs: Objects | None, error: str); never raises
canonicalize(objs: Objects) -> Objects         # CW outer on screen, CCW holes, top-most-then-left-most start, parts sorted
to_mask(parts: list[Part], W, H) -> Geometry   # Geometry(mask bool (H,W), legal bool, defects dict[str,int]); polycore only
MAX_VERTICES = 1000 ; MAX_RINGS = 10                # per ring / per object; beyond = parse error cap_vertices / cap_rings
MAX_CHARS = 16_384                                  # longer text = parse error too_long (F1, I4); only JSON whitespace
                                                    # (space, tab, LF, CR) is stripped around the document (B1-1)
pattern(names) -> str                          # regex of the grammar (caps excluded); portable to vLLM/xgrammar structured outputs

# M3  objdel.checker.topology
min_area(union_area: int) -> int               # max(5, ceil(0.005 * union_area))
instance_topology(mask) -> InstanceTopology     # per component: area, holes (8-conn fg / 4-conn bg, Spaghetti)
    # InstanceTopology(areas, holes, hole_areas: int64 arrays); canonical order (area, holes) desc, hole areas desc
signature(mask, min_area=0) -> tuple[int, ...] # filtered signature, same definition as lab filtered_signature
size_filter(mask, min_area) -> mask             # area opening (8-conn) then area closing (4-conn hole pixel count)
betti(mask) -> (b0, b1)                         # kept for certify; betti3d, FG8, BG4 removed
as_mask(mask, ndim=2) -> bool array            # the one mask normalisation: any nonzero value/byte is foreground,
    # bytes 0/1, C-contiguous, new array; other ndim raises ValueError. Every function above calls it (F1, I1)

# M4  objdel.data
Item(item_id: str, dataset: "riga"|"lidc", group_id: str, split: "train"|"val"|"test",
     image: np.ndarray | None, masks: dict[str, np.ndarray bool (k,H,W)], raters: list[str], meta: dict)
    # image: RIGA uint8 (H,W,3), LIDC int16 (H,W) HU, None = no model image (BinRushed1 reconstructions, B3-1;
    # LIDC without CT); checked on write and load. LIDC meta adds copied_raters {dropped: kept} (identical layers
    # count once, B3-3) and other_nodules [[x0,y0,x1,y1]] (items whose frame overlaps one are not stored, B3-2)
store.write(items, root) -> manifest ; store.verify(root) -> manifest ; store.load(dataset, split, root) -> Iterator[Item]
    # manifest: SHA-256 + rows per file; write: item_id unique across splits (B3-6); load streams one row at a time (I20)
lidc.sop_z(xml_root) -> {series: {sop_uid: z}} (CT fetch) ; lidc.crowded(slice, box) -> bool (frame overlaps another
    nodule, B3-2) ; riga.keypoints(ref) -> (xy, ORB descriptors | None) ; splits.keypoint_duplicates(features,
    threshold=22, pairs=None) -> [(i, j, inliers)] (rotation/scale-tolerant duplicate groups, B3-8; next to
    splits.near_duplicates NCC)

# M5  objdel.targets
frame(item, cfg=None) -> Frame(image, masks: dict[str, (k,H,W)], W, H, box, scale)   # cfg default FRAMES[dataset]
    # Frame.image None when Item.image is None (masks and box still built; prepare skips the item)
    # W == H == S, S divides 1000, 50 ≤ S ≤ 500: one frame pixel = 1000/S grid units; reward rasterizes at W×H.
    # box = (x0, y0, x1, y1) half-open in item pixels, x1-x0 = scale·W (may extend past the image: zero-padded).
    # scale = int item px per frame px; grid g ↦ item px x0 + scale·g·W/1000. All plain ints.
    # image: RIGA uint8 RGB at S; LIDC lung window [−1000, 400] then bicubic to 512 px (image only).
    # image input gray (H,W[,1]) / RGB / RGBA (alpha dropped), else ValueError; masks canonical bools (F1, B1-6, I1).
build_target(masks: dict[str, bool (H,W)], names, cfg=None, min_area=None) -> Target(text, objs, info)
    # min_area: per-object reward size filter (+ "rim"), as ExpertCache computes it: the target keeps the filtered
    # signature at that rule. A name without an entry is certified at every rule a union on the crop can have
    # (own min_area .. min_area(H*W)): safe for any union, possibly more vertices (F1, B1-3). MAX_RINGS counts rings
    # after the size rule (B1-4). Raises ValueError when uncertifiable (side > 500, > 10 rings after the size rule,
    # > 1,000 vertices per ring).
    # Over-long answers (10 / 3,546 RIGA train > 1,536 Qwen tokens; after B3-1 dropped BinRushed1: 1 / 3,282 prepared
    # RIGA train targets for every model, 0 elsewhere): M7 drops them from SFT; RL keeps the items.
target_set_hash(rows (item_id, rater, text)) -> sha256   # order-independent; data/targets manifest + prepare manifest
min_areas(stacks: dict[str, bool (k,H,W)]) -> dict[str, int]   # build_target's min_area for every rater of an item:
    # min_area of each structure's rater union, plus "rim" (disc & ~cup) when disc and cup are present (S1)

# M6  objdel.reward
REWARD_V4                                      # frozen constants (plan §4.3)
ExpertCache.build(stack bool (k,H,W)) ; ExpertCache.save(path) -> digest
ExpertCache.load(path, expected_digest=None, rebuild=False)   # ValueError: corrupted / "digests differ" from expected;
    # StaleCacheError(ValueError): built under another fingerprint (REWARD_V4, size rule, numpy/scipy/OpenCV), unless
    # rebuild=True (tau, T5 rebuilt from the verified masks). cache.voters: raters with a topology vote (F2, B2-4)
rater_answers(caches, names) -> {joint signature tuple: [rater indices]}   # raters voting in every object (all
    # raters if none does); ValueError on unequal k. voters(sigs): a rater whose size-filtered signature is empty
    # abstains (every rater when all do) (B2-4)
score_completion(text, truncated, names, caches: dict[str, ExpertCache], W, H) -> dict
    # keys: format_ok, legal, iou_tc, iou_min (min per-object IoU_tc: the validity gate), f5, g_area, topo_any,
    # topo_freq (share of voting raters), sig_id, error, defects
advantages_v4(terms: dict[str, array (B*G,)], group=8) -> array (B*G,)   # iou_min optional (single object: = iou_tc)
# also exported: score_objects (score_completion + masks + ObjectTerms, for eval), topology_valid (the topology gate
# shared by advantages_v4 and eval), legal_gate (the legality gate shared by advantages_v4, topology_valid and eval;
# R2), score_batch(jobs, pool), rl_eligible (k >= 2), fingerprint() (REWARD_V4 + size rule + library versions;
# cache staleness)

# M7  objdel.train   (import submodules directly; see impl/M7_train.md)
dataset.prepare(store_root, out_root, splits=("train", "val"), procs=3) -> manifest   # data/prepared: images/<id>.png,
    # caches/<id>/<object>.npz (ExpertCache.save), <ds>/<split>.jsonl records, manifest.json; items without pixels skipped
dataset.records(root, dataset, split) -> list[dict]   # item_id, group_id, W, H, box, scale, raters, targets (text | None
    # per rater), rl, disagreement, image_sha256, caches {object: digest}, errors
dataset.sft_rows(root, dataset, split, tokenizer=None, max_tokens=None) -> datasets.Dataset
    # messages [user: {"type":"image"} + PROMPTS[dataset], assistant: target], image (PNG path, decoded lazily), key, rater;
    # one row per *voting* rater (B2-4, reward.rater_answers: 10 / 23,469 LIDC train targets from sub-filter raters
    # dropped; RL keeps the items); answers with tokens + 1 > max_tokens dropped (I5)
dataset.grpo_rows(root, dataset, split, oversample=1) -> datasets.Dataset   # prompt, image, key, W, H, names; rl items;
    # disagreement items repeated `oversample` times (DISAGREEMENT_OVERSAMPLE = 4)
    # both row builders decode every PNG and check it against image_sha256 (ValueError; B3-5)
dataset.cache_loader(root, maxsize=64) -> (key -> {object: ExpertCache})    # LRU; ExpertCache.load(path, record digest):
    # ValueError on a corrupted / swapped file, StaleCacheError after a reward fingerprint change (B2-8)
rewards.stop_ids(model, processing_class) -> set[int]    # TRL's EOS + pad rule for truncation
rewards.truncated(ids, stop) -> bool    # the one truncation rule (reward and sample.sample)
rewards.TRLReward(caches, stop_ids, pool=None).funcs() -> (funcs, weights)
    # funcs named LOGGED = (format_ok, legal, iou_tc, f5, g_area, topo_any, r_topo) with weight 0, then
    # (GATED = iou_tc, f5, g_area are logged times reward.legal_gate: 0 for an illegal rollout, R3)
    # AGGREGATOR = "advantage_v4" (advantages_v4) with weight 1; needs dataset columns key, W, H, names;
    # raises ValueError unless every call holds whole, in-order groups of 8; logs error/<code>, defect/<class>;
    # completions table columns key, truncated, group, error (B4-3); .recorded = [] (set by Controller) collects
    # (key, text, truncated, W, H, names) per scored completion
callbacks.Controller(trainer, metric="eval_rewards/iou_tc/mean", alarm=audit_alarm, baseline=None, audit=True)
    # default metric = mean of the legality-gated iou_tc logger: illegal held-out answers count 0 (R3)
    # TrainerCallback; refuses num_generations / num_generations_eval != 8 (B2-5)
    # on_evaluate: eval.alarms_for on TRL's held-out completions (gathered over processes) -> metrics
    #   eval_audit/<structure|dataset>/<stat>, eval_audit/fired = checks firing now and not on `baseline`
    #   (default: the step-0 evaluation, GRPOConfig eval_on_start=True) (B4-4, I21); then alarm(metrics) -> reason
    #   => .rollback(reason) BEFORE best is updated (B2-1); best trainable params -> <output_dir>/best/
    #   audit statistics are written to controller.jsonl (`audit` event, `stats`) and logged once at the eval
    #   step via Trainer.log (on_evaluate runs after transformers logs eval metrics; review M3 / fixes/R3.md)
    # initial weights are the first best (rollback before any evaluation restores them, B2-2)
    # .switch(reward_weights={name: w}, beta=b) at the next step; <output_dir>/control.json read by rank 0 and
    #   broadcast; unparseable -> retried next step, invalid -> control_error event (B2-3)
configs.MODELS, MAX_COMPLETION (1536 / 1024 LFM), lora(model), processor_kwargs(model, dataset),
configs.sft_config(model, output_dir, **overrides) -> SFTConfig, grpo_config(model, output_dir, **overrides) -> GRPOConfig
    # grpo_config raises ValueError for num_generations / num_generations_eval other than 8 (B2-5)
sample.sample(trainer, rows, k) -> [{"item_id", "samples", "truncated"}]   # k (multiple of 8) completions per grpo_rows
    # row through GRPOTrainer._generate in eval mode (TRL's own generation path)
sample.export(parquet_paths, keys=None) -> rows    # TRL completions tables -> eval_run rows; sample.write(rows, path)
experiments/phase1/sample_eval.py: --checkpoint (K samples per item) | --completions (export) -> eval_run JSONL
    # bf16 only when torch.cuda.is_available() and torch.cuda.is_bf16_supported() (S4)

# M8  objdel.eval   (see impl/M8_eval.md)
metrics(samples: Sequence[str], raters: dict[str, bool (k>=2,H,W)], W, H, names=None, truncated=None, caches=None)
    -> dict    # one item, K >= 2 completions; JSON-native. "samples": per-sample lists format_ok, legal, iou_tc, f5,
               # g_area, topo_any, passed (legal & topo-any & IoU_tc >= 0.5 every object), advantage (v4, groups of 8
               # or one group of K), error; "objects"[name]: per-sample iou_tc, f5, g_area, iou_vote, dice_vote, d_max,
               # area, a_amax, a_amin, m_e (= A/|U|), components, extra_components + soft_dice, ged; item: k, K, group,
               # modes, tv, js, min_mode, coverage, pass_at_k {"1","2",...,"64"}, groups {dup, signal}
               # illegal samples (reward.legal_gate) score 0 on every geometric score and are empty masks in soft
               # Dice / GED; areas and components keep their raster for the audit (R2, REVIEW #2)
summarize(items: list[{"item_id", "dataset", **metrics}]) -> {"datasets": {ds: rates, medians, per-image means}}
hack_audit(items, baseline: hack_audit-result | None = None) -> {"structures", "datasets", "alarms", "flags", "fired"}
    # alarms (plan §9): area_median, area_tail (with a baseline: > baseline share + 10 points), components (needs
    # baseline), duplicates (identical rasters + flags, not text);
    # flags: inflation, deflation, fragmentation, newly_rewarded_failures, topology_collapse (K >= 64);
    # each {"fired": bool | None (not measured), "on": [structures / datasets], "values"}
alarms_for(records [(key, text, truncated, W, H, names)], caches: key -> {object: ExpertCache}, baseline=None)
    -> hack_audit result   # training evaluations (controller alarm hook)
gates.phase2_entry({"summary": summarize(...), ...}, datasets=None) -> {"pass", "datasets": {ds: pass}, "gates": ...}
    # per dataset; default = datasets in the report; ("riga",) for a model dropped from LIDC (plan §7)
    # format, geometry_median, geometry_pass, topology, grpo_duplicates, grpo_sigma
experiments/phase1/eval_run.py: JSONL {"item_id", "samples", "truncated"} + store -> report.json / report.md / items.jsonl
```

## 4. Pinned stack (installed in `.venv`)

Python 3.12; numpy 2.5.4; opencv-python-headless 5.0.0; torch 2.14.1; transformers 5.20.1; trl 1.16.0;
datasets 5.1.0; peft 0.21.2; accelerate 1.15.0; numba 0.68.0; nanobind 3.1.0; scikit-build-core 1.1.1;
hypothesis 6.170; pyarrow 26; pylidc 0.2.3; pydicom 3.0.2; scipy 1.18.1 (part of the reward fingerprint).

## 5. Data

- Fixtures (git-ignored): `data/fixtures/{riga_experts.npz, lidc_experts.npz, lidc_full.npz, lidc_full_meta.csv}`.
  Tests reach them through the `fixtures_dir` fixture in `tests/conftest.py`; mark such tests `@pytest.mark.data`.
- Raw: `data/raw/riga`, `data/raw/lidc_xml`, `data/lidc_slices`; Phase 0 RIGA masks in `data/riga_true`.
- Built (final): `data/store` (build_data: LIDC 7,341 / 1,084 / 2,200 items, patients 709 / 100 / 201; RIGA
  592 / 63 / 95; `build_report.json`), `data/targets` (build_targets: target_set_sha256 `1af819f4…`),
  `data/prepared` (prepare_rl: LIDC train 7,341 items / 23,469 targets / 289 disagreement, val 1,084 / 3,444 / 44;
  RIGA train 547 / 3,282, val 58 / 348, BinRushed1 skipped for lack of pixels). Counts and hashes: `manifest.json`.
- Untrusted: any script that reads `data/` runs as `.venv/bin/python -I ...`. Never commit or copy data out.
- RIGA licence forbids redistribution: the RIGA store stays local.

## 6. Sandbox and resources

- Scratch: `/private/tmp/claude-501/-Users-sreeyashwanthsaivenkatesh-Workspace-projects-CUDA-ObjectDelineation/a7b21ead-acef-47ee-86e8-9b88a7bc677b/scratchpad/impl/<task>/`.
  Docker (Desktop, linux arm64 and amd64 via emulation) is available for build/sanitizer isolation.
- The Mac has 8 cores and 16 GB shared by all agents: use ≤ 3 worker processes and ≤ 3 GB RAM per task.
- If `import objdel` fails briefly while M1 swaps the build backend, wait a minute and retry; do not reinstall.

## 7. Later waves (coordinator)

- **Break-it (wave 5):** 3–4 agents probe wiring and try to break the system with edge cases; every failure goes
  in `experiments/phase1/impl/ISSUES.md` with a reproducing test.
- **Fix (wave 6):** the fleet fixes every issue test-first.
- **Simplify and verify (wave 7):** dead/stale code sweep, line-count review, `/verify`, coverage (K13),
  A/B tests (§8), `experiments/phase1/RESULTS.md`.

### Wave 7 task graph

```
 data/store (F4, manifest 37fe7665…) ──► prepare_rl ──► build_targets ──► goldens frozen (coordinator)
        │                                                                     │
        ▼                                                                     ▼
 S1 geometry ─┐                                                     M9 final x86 run (Docker amd64)
 S2 reward    ├─ parallel, disjoint paths; behaviour-preserving ──► full make test / lint ──► /verify
 S3 data      │  (goldens identical before and after)                   ──► A/B table + ISSUES merge
 S4 train    ─┘                                                         ──► RESULTS.md, docs
 M9 Docker setup ── parallel (first run on pre-sweep code, final run after the sweep)
```

| Task | Owner paths |
|---|---|
| S1 | `cpp/`, `src/objdel/_polycore.cpp`, `src/objdel/checker/`, `src/objdel/geometry/`, `src/objdel/targets/`, `experiments/phase1/{build_targets,census_sizerule}.py`, `tests/test_{polycore,geometry,topology,targets,checker}*.py`, `tests/test_breakit_geometry.py`, `tests/oracles/{raster,validity}_oracle.py` |
| S2 | `src/objdel/reward/`, `src/objdel/eval/`, `experiments/phase1/{reward_gates,eval_run,eval_demo}.py`, `tests/test_{reward,eval}*.py`, `tests/test_breakit_reward.py`, `tests/oracles/` (other files) |
| S3 | `src/objdel/data/`, `experiments/phase1/{build_data,fetch_lidc_pixels}.py`, `tests/test_data*.py`, `tests/test_breakit_data.py` |
| S4 | `src/objdel/train/`, `experiments/phase1/{prepare_rl,smoke_grpo,sample_eval}.py`, `tests/test_train*.py`, `tests/test_breakit_wiring.py`, `tests/test_experiment_scripts.py`, `tests/test_readme.py` |
| M9 | `experiments/phase1/bench/`, `experiments/phase1/impl/M9_x86.md` |
| D1 | `experiments/phase1/impl/{ISSUES,AB}.md` (issue ledger, A/B table AB1–AB7) |
| D2 | `experiments/phase1/RESULTS.md`, `README.md` status, `ROADMAP.md` Phase 1 row |
| D3 | `experiments/phase1/PLAN.md`, this file §3–§6, `impl/M*.md`, `impl/{breakit,fixes}/*.md` (stale references only) |

Sweep rules: no behaviour change (every golden output byte-identical, all tests green); public names in §3 keep
their meaning, and a name may be removed only when nothing in the repo uses it; `experiments/phase1/lab/` is
committed evidence and is out of scope. Report each file's line count before and after in `impl/S<n>.md`.

## 8. A/B tests

| ID | A | B | Measure |
|---|---|---|---|
| AB1 | Phase 0 PIL rasterizer | polycore | wrong pixels vs oracle, time per batch |
| AB2 | Phase 0 reward (`AnnotatorSet`) / reward v2 | reward v4 | adversary suite, neutrality R4, K9 |
| AB3 | reward v4 | v4 + hole-matching term (M6 spike) | stamped-holes win rate, R4 |
| AB4 | g_area zero tolerance | τ·P tolerance | leave-one-rater-out g_area, R3, R4 |
| AB5 | lossless DSS targets | `minseg_hybrid` | tokens per tokenizer, IoU, topology |
| AB6 | scipy labeling | OpenCV Spaghetti | equality, time |
| AB7 | our advantages | TRL advantages (aggregator path) | max abs difference (R2) |
