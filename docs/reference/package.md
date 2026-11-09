# Package reference: `src/objdel`

This page lists every module of the `objdel` package: its purpose, its public API, its inputs and outputs, and its dependencies.
Signatures come from the source. Behaviour comes from the module docstrings and the tests named in the last column of each section.
Names with a leading underscore are private and are not listed.

Package version: `0.1.0` (`objdel.__version__`). Python 3.12 or newer.

## Layers and dependencies

```text
  secrets (no dependencies)

  checker  <-- numpy, scipy, cv2          pure mask geometry, no model, no GPU
     ^
  data     <-- checker, pyarrow, PIL      items and the store
     ^
  geometry <-- objdel._gpucheck           text grammar, canonical form, polygon -> mask
     ^
  targets  <-- data, checker, geometry, numba      frames and certified targets
     ^
  reward   <-- checker, geometry, objdel._gpucheck, scipy, cv2   expert tables, terms, advantages
     ^
  eval     <-- reward, geometry, targets, checker, scipy
     ^
  train    <-- data, geometry, targets, reward, eval, torch, transformers, trl, peft, datasets
```

`objdel._gpucheck` is the compiled extension of [cuda/gpucheck.cpp](cuda-engine.md).
Importing `objdel.train` submodules is lazy on purpose: `objdel/train/__init__.py` holds only a docstring, so using `objdel.train.configs` does not import `objdel.train.fast`.
`objdel.eval.vision` and `objdel.eval.coverage` import without torch.

---

## `objdel.checker`

Exact geometry of binary masks and voxel volumes. No polygon logic: integer polygons live in the CUDA checker.
Exports (`checker/__init__.py`): `InstanceTopology`, `as_mask`, `betti`, `betti3d`, `bitquads`, `certify`, `instance_topology`, `min_area`, `polygonize`, `signature`, `size_filter`.

### `checker.bitquad`

Gray's 1971 bit-quad counts. Depends on numpy.

| Name | Signature | Returns |
|---|---|---|
| `bitquads` | `(mask)` | `BitQuads(n1, n2, nD, n3, n4)` over the zero-padded mask |
| `BitQuads.chi8` | property | `(n1 - n3 - 2 nD) / 4`. Euler number, 8-connected foreground |
| `BitQuads.area` | property | `(n1 + 2 n2 + 3 n3 + 4 n4 + 2 nD) / 4` |
| `BitQuads.perimeter` | property | `n1 + n2 + n3 + 2 nD` unit edges |
| `BitQuads.corners` | property | `n1 + n3 + 2 nD`, the rectilinear ring vertices |

### `checker.topology`

Per-instance topology with the (8, 4) pairing, and the size rule. Depends on cv2 (Spaghetti labelling) and numpy.

| Name | Signature | Returns |
|---|---|---|
| `as_mask` | `(mask, ndim=2) -> ndarray` | Canonical bool array. Any nonzero value is foreground, NaN included |
| `bbox` | `(m) -> (x, y, w, h)` | Foreground bounding box |
| `instance_topology` | `(mask) -> InstanceTopology` | `areas`, `holes` (per 8-connected part, sorted by area then holes, descending) and `hole_areas` |
| `betti` | `(mask) -> (b0, b1)` | Parts and holes |
| `min_area` | `(union_area) -> int` | `max(5, ceil(union_area / 200))` |
| `size_filter` | `(mask, min_area) -> ndarray` | Area opening, then area closing. Idempotent. Never writes the caller's array |
| `signature` | `(mask, min_area=0) -> tuple[int, ...]` | Holes per part, descending, after the size filter |

Callers in worker pools call `cv2.setNumThreads(0)`. On the macOS wheels `setNumThreads(1)` is a no-op (issue I14, in the module docstring).

### `checker.polygonize`

Mask to rectilinear rings with array operations only (edge linking, pointer jumping, list ranking). It mirrors the GPU linker. Depends on numpy.

| Name | Signature | Returns |
|---|---|---|
| `polygonize` | `(mask, offset=(0, 0)) -> list[Ring]` | Rings. Outer rings have positive shoelace area. Holes are negative |
| `Ring` | dataclass | `pts` (k, 2) int64 corners, `area`, property `hole` |
| `boundary_edges`, `link`, `ring_labels`, `list_rank` | | Building blocks. Exposed for tests |

### `checker.certify`

| Name | Signature | Returns |
|---|---|---|
| `certify` | `(mask, min_size=0) -> Report` | `Report(b0, b1, area, perimeter, corners, chi8, certified, failures)`. `certified` is true when bit quads, labelled topology and traced rings agree on chi8, ring counts, area and corners |
| `crop` | `(mask) -> (crop, (x, y))` | Crop to the foreground box |

### `checker.volume`

3-D topology of stacked slices. Depends on numpy and scipy.

| Name | Signature | Returns |
|---|---|---|
| `stack` | `(slices: Mapping[int, ndarray]) -> (z0, volume)` | Slice masks placed by index. Missing slices are empty. Shapes must match |
| `euler` | `(vol) -> int` | `V - E + F - C` of the cube complex |
| `betti3d` | `(vol) -> (b0, b1, b2)` | Parts (26-connected), tunnels (from chi), cavities (6-connected background not touching the border) |
| `topology_any` | `(pred, readers) -> bool` | True when some reader's volume has the predicted `(b0, b1, b2)` |

---

## `objdel.data`

Exact multi-rater data. Exports only `Item`.

### `data.store`

Parquet item store. See [Data](data.md) for layout.

| Name | Signature | Notes |
|---|---|---|
| `Item` | dataclass | `item_id`, `dataset`, `group_id`, `split`, `image` (uint8 RGB for RIGA, int16 HU for LIDC, or `None`), `masks` (`{structure: bool (k, H, W)}`), `raters`, `meta` |
| `write` | `(items, root) -> manifest` | Items ascending by id within each dataset and split. Raises when an id is in two splits or out of order |
| `load` | `(dataset, split, root) -> Iterator[Item]` | Verifies the file's SHA-256 first. Streams one row at a time |
| `verify` | `(root) -> manifest` | Checks every file against the manifest |
| `encode_masks`, `decode_masks` | `(stack)`, `(bits, width)` | `np.packbits` along W, as NPY bytes |

Constants: `DATASETS = ("riga", "lidc")`, `SPLITS = ("train", "val", "test")`, `ROW_GROUP_BYTES = 128 MiB`, zstd level 6.

### `data.lidc`

LIDC-IDRI XML to per-slice multi-reader masks. Rules P1 to P6 are in the module docstring. Depends on PIL and numpy.

| Name | Signature | Notes |
|---|---|---|
| `parse` | `(path)` | One reading-session XML file to sessions of `Roi` |
| `rasterize` | `(rois, zvals) -> {slice index: bool mask}` | Inclusions add, exclusions remove |
| `clusters` | `(sessions, zvals)` | Groups reader annotations into nodules by shared voxels |
| `slices` | `(sessions, zvals) -> list[Slice]` | One `Slice` per multi-reader slice. Identical masks count once (P5). Other nodules are listed (P6) |
| `crowded` | `(s, box) -> bool` | A window that shows another outlined nodule |
| `scans` | `(db=None) -> {series uid: Scan}` | Reads the pylidc SQLite database read-only |
| `ct_slices` | `(npz_path, zvals) -> {slice index: int16 HU}` | |
| `item` | `(series, scan, s, split, image, sessions, xml) -> Item` | Item id `lidc/<series>/n<cluster>/z<slice>`, group id = patient |

Constants: `GRID = 512`, `MARGIN = 8`.

### `data.riga`

Six-rater pen drawings to disc and cup masks, in stages F1 to F7. Depends on `objdel.checker`, scipy, PIL.

| Name | Signature | Notes |
|---|---|---|
| `discover` | `(root) -> list[Photo]` | |
| `extract` | `(photo)` | Staged extractor. Stage A reproduces Phase 0 bit for bit |
| `certified` | `(disc, cup) -> bool` | Each one part, no hole, cup inside disc, rim certified |
| `descriptor`, `keypoints` | `(ref)` | Near-duplicate features for splits |
| `item` | `(photo, ref, kind, masks, split, group, duplicates) -> Item` | Item id `riga/<subset>/image<NNN>` |

### `data.splits`

Seeded, group-aware splits. Depends on numpy.

| Name | Signature | Notes |
|---|---|---|
| `lidc_split` | `(patients, seed=0) -> {patient: split}` | 70 / 10 / 20 |
| `riga_split` | `(subsets, groups, seed=0) -> {image: split}` | Test = Magrabia. 10 % validation carve-out |
| `near_duplicates`, `keypoint_duplicates` | | NCC >= 0.45, or >= 22 RANSAC inliers |
| `group_ids`, `leaks`, `crossing` | | Union-find groups and leak checks |

---

## `objdel.geometry`

The one spec of model output text. See [Output format](output-format.md).
Exports: `GRID`, `MAX_CHARS`, `MAX_RINGS`, `MAX_VERTICES`, `SCHEMAS`, `Geometry`, `Objects`, `Parsed`, `Part`, `canonicalize`, `parse`, `pattern`, `serialize`, `to_mask`, `to_masks`.

| Name | Signature | Returns |
|---|---|---|
| `parse` | `(text, names, truncated=False) -> Parsed` | `Parsed(ok, objs, error)`. Never raises |
| `serialize` | `(objs, names) -> str` | Compact JSON in schema order |
| `pattern` | `(names) -> str` | The grammar as a regex |
| `canonicalize` | `(objs) -> Objects` | Clockwise outers, counter-clockwise holes, top-left starts, sorted parts. Copies. Never edits input |
| `to_mask` | `(parts, W, H) -> Geometry` | `Geometry(mask, legal, defects)` for one object |
| `to_masks` | `(objs: [(parts, W, H)]) -> list[Geometry]` | One `check_polygons` call for all objects |
| `pack_objects`, `pack_rings` | | Pack rings into the checker's arrays. Refuse bad shapes, non-integers and values outside int32 |
| `unpack_rows` | `(words, H, W) -> bool (H, W)` | Unpack the checker's bit-packed mask |
| `signed_area2` | `(ring) -> int` | Twice the signed shoelace area, exact |

---

## `objdel.targets`

Frames and certified targets. Exports: `FRAMES`, `LUNG_WINDOW`, `BuildConfig`, `Frame`, `FrameConfig`, `Target`, `augment`, `build_target`, `frame`, `min_areas`, `target_set_hash`.

### `targets.frame`

| Name | Signature | Notes |
|---|---|---|
| `frame` | `(item, cfg=None, *, offset=(0, 0), side=None, window=None) -> Frame` | `Frame(image, masks, W, H, box, scale)`. Defaults reproduce the policy frame. A moved or resized window must hold the union with 2 pixels of clearance and meet no other outlined nodule, else `ValueError` |
| `policy_side` | `(box, cfg) -> (f, S)` | Reduction factor and side. S divides 1000 |
| `hu_window` | `(a, window=LUNG_WINDOW)` | Hounsfield units to 0..255 (int64, rounded half up) with the lung window (-1000, 400). `ValueError` for an inverted window |
| `FrameConfig` | dataclass | `max_side=500`, `min_side=50`, `margin_px=8`, `margin_pct=25`, `image_side` |
| `FRAMES` | dict | `riga`: defaults. `lidc`: `image_side=512` |

### `targets.build`

| Name | Signature | Notes |
|---|---|---|
| `build_target` | `(masks, names, cfg=None, min_area=None) -> Target` | One rater's masks on one frame to certified canonical text. `Target(text, objs, info)`. Raises `ValueError` when no target certifies |
| `min_areas` | `(stacks) -> {name: int}` | The reward's size rule per structure, plus `rim` for RIGA |
| `target_set_hash` | `(rows) -> sha256` | Order-independent |
| `BuildConfig` | dataclass | `tau=3.0`, `step=0.25`, `iou_min=0.99`, `exact=("nodule",)` |

The tolerance ladder runs tau = 3, 2.75, ..., 0.25, then lossless, then exact rings.

### `targets.augment`

D4 operations, window draws, CT jitter, context crops and the renderer. `OPS = ("id", "r90", "r180", "r270", "fh", "fv", "tr", "at")`. LIDC uses all eight, RIGA uses `id` and `fh`.

| Name | Signature | Notes |
|---|---|---|
| `draw` | `(item, rng, side=None, tries=20)` | One window: side from `sides(item)`, offset uniform with clearance |
| `neighbour_side` | `(item, side, rng)` | The rescaled eval view's side |
| `row_rng` | `(*keys) -> Generator` | Seeded from SHA-256 of the keys, so a draw does not depend on order |
| `render` | `(spec, context) -> uint8 image` | Image of a view row |
| `jitter` | `(rng) -> dict` | CT level -300 +- 60, width 1400 +- 80, noise sigma <= 10 HU |

### `targets.simplify`

Numba kernels. `simplify(path, dirs, W, H, tau, G=1000)` returns a minimum-vertex ring within tolerance. `lattice_path`, `to_grid` are helpers.

---

## `objdel.reward`

Reward v4. See [Reward](reward.md). Exports: `NUMERIC`, `REWARD_V4`, `TERMS`, `ExpertCache`, `ObjectTerms`, `StaleCacheError`, `advantages_v4`, `fingerprint`, `legal_gate`, `loo_bias`, `rater_answers`, `rl_eligible`, `score_batch`, `score_completion`, `score_many`, `score_mask`, `sig_id`, `topology_reward`, `topology_valid`, `voters`.

| Name | Signature | Returns |
|---|---|---|
| `ExpertCache.build` | `(stack) -> ExpertCache` | From a bool stack `(k >= 2, H, W)`. Computes tau and T5 |
| `ExpertCache.device` | `() -> int` | Handle of the resident tables. Uploads on first call in a process |
| `ExpertCache.save` / `load` | `(path)` / `(path, expected_digest=None, rebuild=False)` | `.npz` with SHA-256 and the reward fingerprint. `load` raises `ValueError` on corruption and `StaleCacheError` on a fingerprint change |
| `score_completion` | `(text, truncated, names, caches, W, H) -> dict` | Terms of one answer |
| `score_many` | `(jobs) -> [(terms, {object: ObjectTerms})]` | One checker call for all parsed objects |
| `score_batch` | `(jobs) -> dict` | float64 arrays for `NUMERIC` keys plus `error` and `defects` lists |
| `score_mask` | `(mask, cache) -> ObjectTerms` | Terms of a mask that is not from an answer |
| `advantages_v4` | `(terms, group=8) -> ndarray (B*G,)` | Group advantages, prompt by prompt |
| `topology_reward` | `(terms, group=8) -> ndarray` | `r_topo` |
| `legal_gate`, `topology_valid` | `(T) -> bool array` | The two gates |
| `fingerprint` | `() -> sha256` | Of `REWARD_V4`, the size rule and library versions |

A job is the tuple `(text, truncated, names, caches, W, H)`.
`NUMERIC = ("format_ok", "legal", "iou_tc", "iou_min", "f5", "g_area", "topo_any", "topo_freq", "sig_id")`.

---

## `objdel.eval`

Exam harness. Exports: `alarms_for`, `gates`, `hack_audit`, `metrics`, `summarize`, `vision`. See [Evaluation](../guides/evaluation.md).

| Module | Public API | Purpose |
|---|---|---|
| `metrics` | `metrics(samples, raters, W, H, names=None, truncated=None, caches=None, tokens=None, nsd_tau=...)`, `summarize(items)`, `pass_at_k`, `tv_js`, `min_mode_retention`, `soft_dice`, `hm_iou`, `nsd`, `ged` | Per-item metrics of K sampled answers, and their summary |
| `audit` | `hack_audit(items, baseline=None)`, `alarms_for(records, caches, baseline=None)` | Hedging, inflation, deflation, fragmentation, duplicates |
| `gates` | `phase2_entry(report, datasets=None, prior_a_threshold=None)`, `eligible(summary, vision, dataset)` | Phase 2 entry gates |
| `probes` | `build_precision`, `build_topology`, `score_precision`, `score_topology`, `lidc_stats`, `patches`, `make_item`, `make_hole`, `make_two` | Synthetic precision and topology probes |
| `signal` | `group`, `panel`, `per_item`, `paired`, `select_rb`, `spearman_ci`, `draws`, `ci` | GRPO-signal panel and checkpoint rule R-B |
| `vision` | `decide`, `fixed_set`, `projector_distance`, `next_token_probs`, `image_use`, `check` | G-vision |
| `coverage` | `props`, `distribution`, `gaps` | Train-versus-exam coverage |

---

## `objdel.train`

TRL integration. Import submodules directly. Requires the `train` extra.

| Module | Public API | Purpose |
|---|---|---|
| `configs` | `MODELS`, `MAX_COMPLETION`, `sft_config(model, output_dir, steps=None, branch_at=None, **overrides)`, `grpo_config(model, output_dir, **overrides)`, `lora(model)`, `param_groups(model, base_lr)`, `sampling(model, api="vllm")`, `as_names`, `processor_kwargs(model, dataset)`, `dtype(lora)` | SFT and GRPO configs of the three screened models |
| `dataset` | `prepare(store_root, out_root, splits, procs)`, `records`, `sft_rows`, `prompts`, `prompt_rows`, `grpo_rows`, `grpo_rows_views`, `view_rows`, `cache_loader`, `OpenOnce`, `PROMPTS`, `DISAGREEMENT_OVERSAMPLE` | Prepared set and TRL rows |
| `rewards` | `TRLReward(caches, stop_ids)`, `.funcs()`, `.prefetch(...)`, `stop_ids(model, processing_class)`, `truncated(ids, stop)`, `LOGGED`, `GATED`, `AGGREGATOR` | Reward v4 as TRL reward functions |
| `fast` | `FastGRPO`, `FastForward`, `VisionCache`, `RopeCache`, `ImageMemo`, `lora_merged_params`, `compile_text_layers`, `window`, `by_length`, `by_group_length`, `one_pad_per_image`, `is_ratio`, `packing_wins`, `no_lora_input_cast` | Speed patches of TRL |
| `callbacks` | `Controller(trainer, metric="eval_rewards/iou_tc/mean", alarm=audit_alarm, baseline=None, audit=True)`, `.switch`, `.rollback`, `TraceStep`, `audit_alarm` | Hack audit, best-checkpoint tracking, rollback, guard switches |
| `serve` | `preload(root, datasets, threads=4)`, `loader(*roots)`, `CacheRef`, `fork_pool(n, modules=())`, `pool_size()`, `cpus()`, `quota` | Preloaded expert tables and the warmed fork pool |
| `sample` | `sample(trainer, rows, k, items=16)`, `export`, `write`, `export_prompts`, `import_completions` | Checkpoint to K answers per item, and the JSONL bridge to vLLM |

Key behaviours:

- `grpo_config` raises `ValueError` for `num_generations` other than 8. Reward v4 scores groups of exactly 8.
- `TRLReward` raises `ValueError` when a reward call does not hold whole groups of 8, in order.
- `FastGRPO.get_train_dataloader` raises `ValueError` unless `dataloader_num_workers` is 0.
- `Controller` refuses a `num_generations` other than 8 and refuses to raise `beta` from 0 without a reference model.
- `serve.fork_pool` is refused once CUDA is initialised. Call it before the model loads.

---

## `objdel.secrets`

| Name | Signature | Notes |
|---|---|---|
| `get` | `(name) -> str or None` | Looks in `<repo>/.secrets` first, then the macOS keychain. Raises `PermissionError` when `.secrets` is readable by others (needs `chmod 600`). Cached |
| `path` | `() -> Path` | `$OBJDEL_SECRETS` or `<repo>/.secrets` |

Values are never printed. `.secrets` is git-ignored and never sent to a box.

---

## `objdel._gpucheck` (extension)

Nine functions and one constant, built from `cuda/`. See [CUDA engine](cuda-engine.md).
