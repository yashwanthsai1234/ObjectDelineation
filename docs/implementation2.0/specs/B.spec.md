# Phase B specification: quality (B1 data, B2 number match, B3 reward v5)

Revision 2 (2026-11-08), kept consistent with phases/B.md revision 2. This file is the contract for test writers. It
lists the goal, required behaviours, interfaces, constraints and acceptance gates. It gives no design reasons and no
implementation detail.

---

## 1. Goal

1. B1: RL rows match the shifted exam's (object size x object offset) distribution, pass data-quality gates, and
   reach the trainer in a fixed, audited order.
2. B2: the trainer-vLLM mismatch is measured and attributed, and a numerics preset is chosen that keeps almost every
   answer's signal.
3. B3: reward version 5 gives answers that miss the nodule a learning signal toward it, uses a leave-one-out (RLOO)
   baseline, and cannot be gamed by the listed shapes.

---

## 2. Definitions

- **Union box**: bounding box of the union of all doctors' masks of an item, half-open pixel box (x0, y0, x1, y1).
- **Extent** (px): `max(x1 - x0, y1 - y0)` of the union box.
- **Offset**: Euclidean distance between the union-box centre `((x0+x1)/2, (y0+y1)/2)` and the window centre,
  divided by the window side S (pixels of the source slice).
- **Cell**: `(size_bin, offset_bin)`, with `size_bin = searchsorted((10, 20, 40), extent, side="right")` and
  `offset_bin = searchsorted((0.15, 0.2, 0.3), offset, side="right")`; 4 x 4 cells.
- **Kinds**: `centred` (the policy frame), `shifted` (policy side S0, moved), `rescaled` (neighbouring side, moved).
- **RL-eligible item**: an item with 2 or more doctors (`rl_eligible(raters)`).
- **Val-shape / val-select**: the two patient halves of the val split. A val item is in val-shape when the first
  byte of SHA-256(group_id) is even, else in val-select.
- **Miss**: an answer with IoU_tc < 0.1. **Found**: IoU_tc >= 0.1.
- **Step**: 32 prompts x 8 answers = 256 answers.
- **Paired comparison**: per-item differences of two runs on the same val-select items, with a nodule-level
  bootstrap (resample nodules, 95 % CI).

---

## 3. Interfaces and behaviours

### 3.0 Cache fingerprint and val halves

```python
@dataclass(frozen=True)
class CacheSpec:                       # objdel.reward.spec
    expert_score: float = 0.95
    tau_floor: float = 1.0
    tau_quantile: float = 95.0

def fingerprint() -> str               # objdel.reward
def val_half(group_id: str) -> str     # objdel.data.splits; "shape" | "select"
```

- `fingerprint()` changes when `expert_score`, `tau_floor`, `tau_quantile`, the size rule, or the numpy / scipy /
  OpenCV version changes; it does NOT change when any reward weight changes (`alpha_loc`, `miss_iou`, `w_f5`,
  `w_area`, `c_topo`, `adv_scale`, `loc_width`, `area_gate`, `lam`, `eps`, `drop_spread`, `iou_gate`, or any
  `RewardV4` weight or floor).
- Stored cache values (tau, T5) for the same masks are unchanged by the new fingerprint; only digests change.
- `val_half` is deterministic, returns `"shape"` iff `hashlib.sha256(group_id.encode()).digest()[0] % 2 == 0`, and
  keeps all items of one patient in one half.

### 3.1 `objdel.data.windows`

```python
SIZE_EDGES: tuple[int, int, int] = (10, 20, 40)
OFFSET_EDGES: tuple[float, float, float] = (0.15, 0.2, 0.3)
KINDS: tuple[str, str, str] = ("centred", "shifted", "rescaled")

@dataclass(frozen=True)
class Window:
    item_id: str; group_id: str; nodule: str; kind: str; draw: int
    dx: int; dy: int; side: int
    offset: float; extent_px: int
    disagreement: bool
    @property
    def window_id(self) -> str: ...          # f"{item_id}~{kind}{draw}"
    @property
    def cell(self) -> tuple[int, int]: ...

@dataclass(frozen=True)
class WindowPool:
    windows: tuple[Window, ...]
    dropped: Mapping[str, int]               # keys "not_rl", "shifted_none", "rescaled_none"

def offset_of(union_box: tuple[int, int, int, int], window_box: tuple[int, int, int, int]) -> float
def cell(extent_px: int, offset: float) -> tuple[int, int]
def pool(items: Iterable[Item], draws: int, seed: int = 0) -> WindowPool
def cell_shares(windows: Sequence[Window]) -> np.ndarray            # shape (4, 4), sums to 1
def schedule(windows: Sequence[Window], target: np.ndarray, steps: int, per_step: int = 32,
             shares: Mapping[str, float] = {"centred": 0.15, "shifted": 0.75, "rescaled": 0.10},
             seed: int = 0, disagreement_weight: float = 4.0) -> list[Window]
```

Behaviours:

- `cell`: values on a bin edge go to the upper bin (offset 0.15 -> bin 1; extent 10 -> bin 1).
- `window_id`: the part after the source `item_id` contains no "/"; `window_id.rsplit("/", 1)[0]` equals
  `item_id.rsplit("/", 1)[0]`.
- `pool`:
  - items that are not RL-eligible give no window and are counted in `dropped["not_rl"]`;
  - exactly one `centred` window (draw 0) per RL-eligible item;
  - for each `d in range(draws)`, at most one `shifted` and at most one `rescaled` window per item;
  - every `shifted` window has `side` equal to the item's policy side S0; every `rescaled` window has a side other
    than S0;
  - draw 0 has the same (dx, dy, side) as `experiments.phase2.build_eval_views` builds for that item and kind with
    the same seed (for rescaled: where that view is not a kept-side window);
  - a draw with no admissible window is left out and counted in `dropped["shifted_none"]` or
    `dropped["rescaled_none"]`; no fallback or kept-side window is returned;
  - every window keeps every doctor pixel inside, with >= 2 px clearance, and contains no other outlined nodule;
  - `windows` is sorted by `window_id`; the result does not depend on the order of `items`;
  - same inputs give an identical result.
- `schedule`:
  - returns exactly `steps * per_step` windows, step-major (rows `per_step*s .. per_step*(s+1)-1` are step s);
  - kind counts over the run equal the largest-remainder rounding of `shares * steps * per_step`;
  - for each off-centre kind with n rows, the per-cell counts equal the largest-remainder rounding of `target * n`;
  - no two windows in one step share `item_id`;
  - no `window_id` appears twice in the run;
  - `target` must have shape (4, 4), be non-negative and sum to 1 (within 1e-9), else `ValueError`;
  - `shares` must have the three kinds, be non-negative and sum to 1, else `ValueError`;
  - raises `ValueError` naming the kind and the cell when a (kind, cell) cannot be filled without breaking the
    rules above; no window of another kind is substituted;
  - same arguments give an identical list; a different seed gives a different list (for a pool larger than needed);
  - a larger `disagreement_weight` gives a larger (or equal) share of disagreement items, all else equal.

### 3.2 `objdel.data.coverage` and `objdel.data.quality`

`objdel.eval.coverage` moves to `objdel.data.coverage` (same public functions, same behaviour); `objdel.eval.coverage`
no longer exists. New:

```python
def cell_coverage(scheduled: Sequence[Window], target: np.ndarray, target_offsets: Sequence[float]) -> dict
    # keys "worst_ratio", "tv", "ks" over off-centre rows, cells with target share >= 1 %
```

```python
@dataclass(frozen=True)
class Gate:                                        # objdel.data.quality
    name: str; value: float; limit: str; hard: bool; ok: bool

def leakage(store_root) -> dict
def window_checks(windows: Sequence[Window], store_root) -> dict
def schedule_report(scheduled: Sequence[Window], pool: WindowPool, target: np.ndarray, train_law: np.ndarray) -> dict
def tokens(rows_root, tokenizer, limit: int) -> dict
def gates(report: dict) -> list[Gate]          # report = merged dicts of the four functions above
def failures(gates: Sequence[Gate]) -> list[Gate]   # gates with hard=True and ok=False
```

Behaviours:

- `gates` returns one `Gate` per row D1-D16 of section 5.1 (names `"D1"` .. `"D16"`); D14 and D16 have `hard=False`.
- A planted cross-split patient fails D1. A planted shared pixel hash across patients fails D2. A window that cuts a
  doctor pixel fails D4. A policy frame that contains another outlined nodule fails D5. A planted fallback or
  kept-side window fails D6. A schedule with no off-centre rows fails D7/D8. A target of more than `limit - 100`
  tokens fails D12.
- `cell_coverage` and D7-D9 use the off-centre rows only (shifted + rescaled); D11 uses all rows.
- The `objdel.data.quality` module name is not `audit`; `objdel.eval.audit` is unchanged.

### 3.3 Scripts (run as modules with `python -I -m`)

- `experiments.phase2.build_rl_rows --store DIR --target VIEW_DIR --out DIR --steps N --per-step 32 --draws K
  --shares C S R --seed N --procs N`
  - the target shares come from the val-shape items of `VIEW_DIR` only;
  - writes a prepared-layout root: `images/`, `caches/`, `lidc/train.jsonl`, `manifest.json`;
  - `lidc/train.jsonl` lines are in schedule order (not sorted);
  - each record: `item_id` = the window's `window_id`; `source_item_id`; `cell` (string `"ij"`); `group_id`;
    `view = {kind, draw, offset: [dx, dy], side, fallback: false, kept_side: false, op: "id"}`; plus every field
    today's prepared records have (targets, caches, image_sha256, rl, disagreement, raters, box, W, H, scale);
  - every record has `rl = true`;
  - image and cache paths are unique per window (no two records share an image file or a cache folder);
  - `manifest.json` has the usual `files`, `sets`, `reward_fingerprint`, `versions`, plus
    `schedule = {steps, per_step, seed, draws, shares, target_view, target_half, target_sha256, ids_sha256,
    dropped, gates}` where `target_half = "shape"` and `ids_sha256` = SHA-256 of the newline-joined ordered
    `item_id`s;
  - if any hard gate fails: exit code 1 and no `manifest.json`.
- `experiments.phase2.audit_data --rows DIR --target VIEW_DIR --store DIR [--json FILE]`: prints D1-D16; exit 0 iff
  no hard failure.
- `experiments.phase2.build_eval_views ... --kinds shifted[,rescaled,rotated]`: builds only the named kinds (default:
  all three, as today).
- `docs/reference/frozen_views.json` holds, for `lidc-test-shifted`: `view`, `items`, `content_sha256`,
  `target_set_sha256`, `seed`, `kinds`, `store_manifest_sha256`. `content_sha256` = SHA-256 of the newline-joined
  ordered tuples `(item_id, view.offset, view.side, image_sha256)` of the view's records.
- `experiments.phase2.final_eval`: for a test view, recomputes `content_sha256` and `target_set_sha256` and refuses
  (non-zero exit, clear message) when either differs from `frozen_views.json`. The check does not depend on the
  manifest's `versions` or `reward_fingerprint`.
- `experiments.phase2.check_series_z`: lists, for the 4 named LIDC series, z values shared by different SOP UIDs in
  the LIDC XML; reads no DICOM unless a DICOM folder is given.
- `experiments.phase2.fp16_vision_check`: runs the checkpoint's vision tower in float16 on 32 train images and
  prints the inf/NaN count; exit 1 if the count is > 0.
- `experiments.phase2.b2_attribution`: prints the mean |trainer - vLLM| logprob gap for one batch of 256 answers under
  (a) the FastGRPO default forward, (b) the plain TRL forward, (c) the forward with the LoRA delta merged into bf16
  weights.
- `experiments.phase2.train_rl`:
  - `--train DIR` whose manifest has `schedule`:
    - rows are used in file order; the GRPOConfig has `shuffle_dataset=False`;
    - `ValueError` before training when the row count differs from `schedule.steps * schedule.per_step`, or the
      SHA-256 of the ordered row keys differs from `schedule.ids_sha256`;
    - `ValueError` before training when `generation_batch_size // num_generations != schedule.per_step`, or
      `steps_per_generation != gradient_accumulation_steps`, or `num_iterations != 1`, or `--steps` >
      `schedule.steps`;
  - `--train DIR` without `schedule`: behaviour as before (shuffled rows, disagreement x4), until gate B1 passes;
  - after gate B1 passes: `--train-views` no longer exists (argparse error), a root without `schedule` is refused,
    `dataset.grpo_rows_views` and `serve` multi-root key tagging no longer exist;
  - `--numerics {bf16-seqmask,bf16-tis,bf16-geomask,fp16-tis}` (default `configs.DEFAULT_NUMERICS`; after B2.4 only
    the chosen preset remains);
  - `--reward {v4,v5}` (until v4 is archived);
  - `--budget-usd X --rate-usd-per-h Y`: training stops and saves once elapsed hours x Y >= X;
  - `rl.json` adds `rows = {root, ids_sha256}`, `numerics`, `reward = {name, fingerprint, c_topo, adv_scale,
    loc_width, area_gate}` (v5 fields only for v5).
- Prompt rows built from a root whose records have `view` carry columns `kind` and `cell`.
- `TRLReward` logs, per reward call, `kind/<kind>/iou_tc`, `cell/<ij>/iou_tc` and `cell/<ij>/allmiss_groups` (share
  of groups whose 8 answers are all misses), when `kind`/`cell` columns are present.

### 3.4 Number match: `objdel.train.fast` and `objdel.train.configs`

```python
def is_ratio(logps, sampling, mask, mode: str, lo: float | None, hi: float | None,
             geo: tuple[float, float] | None = None) -> Tensor
def is_summary(logps: Tensor, sampling: Tensor, mask: Tensor, ratio: Tensor, mode: str,
               lo: float | None, hi: float | None) -> Tensor   # float32, on the input device
NUMERICS: dict[str, dict]       # keys "bf16-seqmask", "bf16-tis", "bf16-geomask", "fp16-tis"
DEFAULT_NUMERICS: str
def grpo_config(model: str, output_dir: str, numerics: str = DEFAULT_NUMERICS, reward=REWARD_V5,
                **overrides) -> GRPOConfig
def vllm_dtype(dtype: str) -> ContextManager   # in objdel.train.fast
```

Behaviours:

- `NUMERICS["bf16-seqmask"]` sets `vllm_importance_sampling_mode="sequence_mask"`, `clip_max=3.0`, `clip_min=None`.
- `NUMERICS["bf16-tis"]` sets `vllm_importance_sampling_mode="token_truncate"`, `clip_max=2.0`, `clip_min=None`.
- `NUMERICS["bf16-geomask"]` sets mode `"sequence_geo_mask"`, `clip_max=2.0`, geo bounds (0.99, 1.01).
- `NUMERICS["fp16-tis"]` = bf16-tis plus `fp16=True`, `bf16=False`; the base model loads in float16; vLLM runs in
  float16. A GRPOConfig never has both `bf16` and `fp16` true.
- `is_ratio(..., mode="sequence_geo_mask", hi=C, geo=(a, b))`: per token `min(exp(logps - sampling), C)` on valid
  tokens, times 0 for every sequence whose `exp(mean over valid tokens of (logps - sampling))` is outside `[a, b]`.
  Existing modes give the same values as before.
- `FastGRPO` with mode `"sequence_geo_mask"` raises `ValueError` at construction when its own ratio path is not
  active.
- `vllm_dtype("float16")`: inside the context, `trl.generation.vllm_generation.LLM(...)` receives `dtype="float16"`;
  after the context, the original name is restored.
- With fp16-tis, `FastGRPO` raises at construction unless the vLLM engine's model dtype is `torch.float16`.
- `FastGRPO` logs, once per logging step, in addition to the existing 5 `sampling/*` metrics:
  `sampling/is_masked_share`, `sampling/is_clamped_share`, `sampling/seq_logdiff_mean`, `sampling/seq_logdiff_max`,
  `sampling/tok_absdiff_gt_0.5_share`, `sampling/masked_len_mean`, `sampling/kept_len_mean`.
  - masked = ratio units set to 0 by a `*_mask` mode (incl. `sequence_geo_mask`, counted per sequence); clamped =
    ratio units equal to `lo` or `hi` after a `*_truncate` clamp where the raw ratio was outside `[lo, hi]`; ratio
    unit = sequence in `sequence_*` modes, valid token in `token_*` modes;
  - seq_logdiff = signed sum over a completion's valid tokens of `(logps - sampling)`;
  - tok_absdiff_gt_0.5_share = share of valid tokens with `|logps - sampling| > 0.5`;
  - masked_len_mean / kept_len_mean = mean valid-token count of masked / kept sequences (0 when none).
- The accumulated statistics stay on the device until the logging call (no host read per micro-batch).
- The existing 5 `sampling/*` metrics keep their values for the same inputs.
- `grpo_config(..., reward=R)` raises `ValueError` when `num_generations` or `num_generations_eval` differs from
  `R.group`.

### 3.5 CUDA scorer (`objdel._gpucheck`)

- `score_polygons(...)` and `score_masks(...)` results contain two new int64 arrays, one value per object:
  - `sum_x` = sum of column indices x of the object's mask pixels;
  - `sum_y` = sum of row indices y of the object's mask pixels.
- An empty mask gives `sum_x = sum_y = 0`.
- Every previously returned field is bit-identical to the previous build for the same input.
- The CPU build and the CUDA build return identical `sum_x`, `sum_y`.
- This section applies only when the scorer change is made.

### 3.6 Reward v5: `objdel.reward`

```python
@dataclass(frozen=True)
class RewardV5:
    c_topo: float                    # required
    adv_scale: float                 # required, > 0
    name: str = "reward_v5"
    alpha_loc: float = 0.1
    miss_iou: float = 0.1
    w_f5: float = 0.1
    w_area: float = 0.1
    loc_width: float = 1.0
    area_gate: bool = False
    lam: float = 0.5
    eps: float = 0.1
    group: int = 8
    drop_spread: float = 0.01
    iou_gate: float = 0.5

REWARD_V5: RewardV5
def location(sum_x: int, sum_y: int, area: int, ref_sum: tuple[int, int, int], radius2: float,
             width: float = 1.0) -> float
def advantages_v5(terms: dict, group: int = 8, spec: RewardV5 = REWARD_V5) -> np.ndarray
```

`ExpertCache` gains:

- `ref_sum -> tuple[int, int, int]`: `(sum of x, sum of y, pixel count)` of the vote mask if it has pixels, else of
  the union mask.
- `radius2 -> float`: union area / pi.

Behaviours:

- `RewardV5(...)` with `alpha_loc > miss_iou`, or `adv_scale <= 0`, raises `ValueError`.
- `location`: returns 0.0 when `area == 0`; else a value in (0, 1] that is 1.0 when the centroids are equal and
  strictly decreasing in `d2 = (sum_x/area - rx/n)^2 + (sum_y/area - ry/n)^2`, `(rx, ry, n) = ref_sum`; with the
  Gaussian kernel it equals `exp(-d2 / (2 * width^2 * radius2))`.
- `score_completion` / `score_batch` output gains numeric key `"loc"`:
  - per object: `L = 1[iou_tc < miss_iou] * max(0, 1 - iou_tc / miss_iou) * g_area * 1[topo_freq > 0] * location(...)`;
  - completion `loc` = mean of the objects' L;
  - unparseable or truncated completion: `loc = 0`.
- `advantages_v5` (input: flat terms, prompt by prompt, G contiguous answers per prompt):
  - `area = w_area * g_area`, multiplied by `min(1, iou_tc / miss_iou)` when `area_gate` is true;
  - `R = legal * (iou_tc + alpha_loc*loc + w_f5*f5 + area) + c_topo * r_topo`, where `legal` is the hard gate
    (`legal >= 0.5`) and `r_topo` is `topology_reward(terms, group)`;
  - `A_i = adv_scale * (R_i - mean_{l != i} R_l)` within each group;
  - groups with `max(R) - min(R) < drop_spread` and one `sig_id` get `A = 0`;
  - raises `ValueError` when term lengths differ or are not a multiple of `group`.
- `TRLReward(caches, stop_ids, spec=REWARD_V5)`: logged terms include `loc`; the aggregator returns
  `advantages_v5`; TRL's resulting advantage equals `advantages_v5` to float32 precision.

### 3.7 Run guards (`objdel.train.callbacks`)

- `BudgetStop(rate_usd_per_h, cap_usd, start_time)`: sets `should_training_stop` and `should_save` at the first step
  where elapsed hours x rate >= cap; never before.
- `LocGuard(window=50)`: sets `should_training_stop` and `should_save` when, over the last `window` logged steps, the
  least-squares slope of `rewards/loc/mean` is > 0 and that of `rewards/iou_tc/mean` is <= 0; never with fewer than
  `window` logged steps; inactive when `rewards/loc/mean` is not logged.

---

## 4. Constraints

- Library code under `src/objdel`, tests under `tests/`, scripts under `experiments/phase2/` run as modules.
  `make test` and `make lint` pass.
- Python that reads `data/` runs with `-I`.
- No new third-party dependency. TRL stays 1.17.1. Any new replacement of a private TRL name is listed in the
  start-up private-name check.
- Masks never leave the GPU during training; only numbers return.
- Old paths are removed only after the gate that replaces them passes: `dataset.grpo_rows_views`, `serve` multi-root
  key tagging, `train_rl --train-views` and the no-schedule branch after gate B1; reward v4 moves to
  `archive/reward_v4/` after gate B3. `experiments/phase2/configs/rl_shifted.yaml` moves to `archive/`.
- Val-select items are never read by training, by `build_rl_rows` or by the B3.7 lab round.
- The test view is scored once in Phase B, for the final candidate.
- No paid GPU run without a passing preflight, owner approval and a North America host. The Phase B total must
  stay within the owner's approved budget (today $3 for GPU scripts).
- Exception pending owner sign-off: B1 targets use the numpy polygonizer until the CUDA polygonizer passes parity.
- No commit message or PR text with AI attribution.

---

## 5. Acceptance gates

### 5.1 Data gates (`quality.gates`; hard unless "report")

| Name | Check | Pass |
|---|---|---|
| D1 | Patient, series, item overlap between train, val, test | 0 each |
| D2 | Exact slice or 32x32 nodule-crop SHA-256 shared across patients or splits | 0 |
| D3 | Near duplicates (cosine >= 0.98 of centred 32x32 HU crops) cross-split | rate <= within-train rate (4 / 7,341) |
| D4 | Every doctor pixel inside each scheduled window | 100 % |
| D5 | Clearance >= 2 px and no other outlined nodule, in every scheduled window incl. centred ones | 100 % |
| D6 | Fallback windows and kept-side windows in the schedule | 0 and 0 |
| D7 | Implementation check. Off-centre rows: min over cells with val-shape target share >= 1 % of (row share / target share) | >= 0.8 |
| D8 | Implementation check. Off-centre rows: total variation distance to val-shape target cell shares | <= 0.10 |
| D9 | Implementation check. Off-centre rows: KS distance of offset to the val-shape offsets | <= 0.06 |
| D10 | Share of off-centre rows at the item's policy side | >= 0.75 |
| D11 | Kind counts and per-kind cell counts equal the largest-remainder counts | exact |
| D12 | Max certified target tokens + 1 below `MAX_COMPLETION - 100`; target errors | true; 0 |
| D13 | Built window equals a fresh render (box, side, image SHA-256) on 150 sampled windows | 0 mismatches |
| D14 (report) | Disagreement share, top-10 % nodule share, low-contrast share (< 5 grey levels), padded share, pool dropped counts | reported |
| D15 | Test shifted view `content_sha256` and `target_set_sha256` equal `docs/reference/frozen_views.json` | equal |
| D16 (report) | Off-centre coverage (worst ratio, TV, KS) against the train-split shifted law | reported |

### 5.2 Reward v5 anti-gaming tests (before any paid B3 run)

| Id | Test | Pass |
|---|---|---|
| T1 | For fixed f5, g_area, r_topo, legal and `area_gate = False`: R is continuous at IoU_tc = miss_iou and non-decreasing in IoU_tc over [0, 1] for any location value in [0, 1] | true |
| T2 | A legal outline that covers the union and is centred on the reference centroid, with area above the doctors' max area so that g_area = 0 | loc = 0 |
| T3 | A legal speck (area below the doctors' min area, g_area = 0) at the reference centroid | loc = 0 |
| T4 | Two legal disjoint blobs, each above the size filter, none overlapping the union, whose joint centroid equals the reference centroid, for an item where every doctor drew one component | loc = 0 |
| T5 | One legal shape of doctor-like area shifted by d px away from the nodule, all positions with IoU_tc = 0 | loc strictly decreasing in d |
| T6 | Illegal or unparseable answer with a perfect centroid | loc = 0 and R = 0 |
| T7 | Any group: sum of A over the group | 0 (abs < 1e-9) |
| T8 | Group rewards multiplied by k > 0 (all terms scaled) | A multiplied by k (no std division) |
| T9 | Group of 8 identical answers | A = 0, no NaN |
| T10 | `sum_x`, `sum_y` vs numpy on 200 random masks and all fixture masks; CPU build vs CUDA build | equal |
| T11 | Any answer with IoU_tc >= miss_iou vs any miss, equal f5 and g_area, all other terms equal | hit R >= miss R |
| T12 | Group of 8 misses at different distances, as in Direction.md s5 (radius 7 px, distances 5, 10, 20 px) | the closest answer has the largest A, and A > 0 |
| T13 | The same group under `adv_scale = s` and `adv_scale = 1` | A(s) = s x A(1) |

### 5.3 Run gates

All run comparisons are paired comparisons on val-select (shifted view, K = 8). R5 per-item values come from
`runs/r5/eval/eval_val/items.jsonl` filtered to val-select. A gate passes when both its point rule and its CI rule
hold.

- **Gate B2** (probe, 20 steps per arm, same start checkpoint, same 640 rows in the same order): an arm passes when
  `importance_sampling_ratio/mean >= 0.97`; `is_masked_share < 0.5 %` (mask modes) or `is_clamped_share < 0.5 %`
  (truncate modes); NaN/inf count 0; for fp16-tis the dtype check passed; the one-sided 95 % t lower bound of the
  mean paired per-step difference (arm - control) of `rewards/iou_tc/mean` over steps 11-20 is > -0.02.
- **Gate B1** (after the B1 run, against R5): CI lower bound of (B1 - R5) IoU_tc > 0 (R5 0.599); miss share at
  offset > 0.3 <= 23.9 % (R5 31.8 %) and CI upper bound of (B1 - R5) < 0; CI lower bound of (B1 - R5) IoU_tc at
  offset < 0.15 > -0.02 (R5 0.732); legal share >= 0.99. Full-val IoU is reported against the target 0.66 (not a
  gate). The test view is not read.
- **Gate B3** (before the run): T1-T13 pass; on val-shape answers, Spearman(A, -d / r) over all-miss groups > 0.3
  under v5 and greater than under v4; share of groups where a miss has a higher R than a hit <= 1 %.
  (after the run, against the B1 run): miss share falls by >= 25 % relative and CI upper bound of (B3 - B1) miss
  share < 0; on items found by both runs, CI lower bound of (B3 - B1) IoU_tc > -0.01; legal share >= 0.99.
