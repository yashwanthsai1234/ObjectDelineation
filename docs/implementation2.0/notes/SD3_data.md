# SD3: Data system design, data-centric lens (raw LIDC/RIGA to versioned, gated, GPU-fed rows), up to Phase E

Author role: fresh system designer, no prior opinion. Date 2026-11-08. Diagram:
`docs/implementation2.0/notes/SD3_data_diagram.html` (self-contained HTML + SVG).
Fact tags: [code: file:line], [measured: source], [paper/web: link], [estimate]. SP = the session scratchpad.

---

## 0. Skills used (and what each changed)

| Skill | What it changed in this design |
|---|---|
| system-design-skills:system-design + :design | Gave the method: clarify, estimate, design, trade-offs, failure, iterate, score. I ran the loop inline and did not start the orchestrator agent. I am the fresh designer, so a second agent would only repeat the work. Section order follows its design-doc template. |
| system-design-skills:back-of-the-envelope | Made me size the work from measured counts before I chose any part. Result: the per-step data need is about 1.4 images/s, so throughput never drives this design. Latency on the CPU, VRAM bytes and transfer bytes do (section 3). |
| system-design-skills:data-storage | "Start from the queries, write the key first." So every table below has a primary key and its one access path. Only a few small tables and a single writer exist, so it rejects a database server (YAGNI) and keeps Parquet files plus manifests. |
| system-design-skills:blob-store | "Keep the bytes in a blob and a pointer in the table. Too many tiny objects make the index the ceiling." So 9,635 cache `.npz` files and 9,030 PNGs (the training path) become one flat context blob plus an index. Exams stay as frozen, immutable objects. |
| system-design-skills:caching | Stops device tables from being cached for the whole run. Today they are uploaded once and stay resident, which grows to about 1.2 GB with 3 views. Instead, tables are built per step and freed. Keys carry versions: the digest includes the reward fingerprint. |
| system-design-skills:task-scheduling | Its YAGNI check says no scheduler and no Airflow. The build is a one-shot DAG run by `make`, and every task in it is idempotent and content-addressed. |
| system-design-skills:architecture-diagram | Rendered the settled design to `SD3_data_diagram.html`. |
| huggingface-skills:huggingface-datasets | Chose transport: LIDC-derived Parquet goes to a private HF dataset repo, pinned by git revision. The Dataset Viewer API does not apply here (private repo, gated data). RIGA is never uploaded. |
| nvidia-skills:nvidia-skill-finder | Catalog check (`npx skills add nvidia/skills --list` printed nothing here, so I used the catalog dump `(planning probe, not kept: nv.json)`). It has `dicom-series-to-volume`, `dicom-series-preflight`, `dicom-metadata-extract`, `dali-dynamic-mode`. None of them is needed up to Phase E, because decoding happens once. `dicom-series-to-volume` is the candidate for the 3D-voxel phase. I installed nothing. |
| engineering:system-design / :architecture | Requirements and constraints framing. ADR-style trade-off rows (section 7). |
| engineering:code-review | I read the builders as a reviewer. Findings that shaped the design: static one-draw train windows; device-table growth; the coverage gap measured in section 5. |
| engineering:debug | Reproduce, then isolate. I re-ran the geometry audit (`docs/implementation2.0/seeds/b1_scripts/b1_analyze.py`) and did not trust the summary numbers. |
| engineering:documentation | Architecture-doc structure: context, design, decisions, data flow. |
| huggingface-skills:huggingface-best | Checked for scope only. Model choice is outside this data design, so it changed nothing. |

Web research (section 9) covered: lakeFS taking over DVC; HF Xet; Lance; Croissant-RAI; leakage studies; nvImageCodec on small images; cleanlab segmentation.

---

## 1. Requirements

### Functional (core, in scope)
1. **Certify.** Turn raw LIDC XML + CT slices and RIGA photos into the item store, then into a "doctor record" per item. The record holds the rater masks (exact), the size rule, the voters, the topology answers and a certified target for the policy frame.
2. **Views.** Turn each item into window specs: centred, shifted, rescaled (train and exams), plus rotated (LIDC exam only). Each window carries its own doctor masks.
3. **Row sets.** RL prompt rows (spec id + mixture weight), SFT rows (spec + certified target), and frozen eval splits (val/test: shifted, rescaled, rotated).
4. **Quality gates and audits.** Coverage against the exam distribution, leakage, doctor masks inside windows, duplicates, CPU/CUDA parity, outline lengths. A failed gate blocks the version.
5. **GPU-fed training store.** On the rented box, images are rendered and doctor tables are built per step on the GPU from compact resident data. No PNG decode, no per-window cache files.
6. **Versions.** One dataset version id per build. Every run and every result JSON names it, and a run refuses to start on a mismatch.

### Non-functional
- **Deterministic.** The same inputs give byte-identical outputs and the same version id. This already holds for the store, views and bundles [code: src/objdel/data/store.py:3-5; experiments/phase2/build_views.py:18-22].
- **Integrity.** Every file has a SHA-256 in a manifest. Data is checked at build and again at load.
- **Cost.** Offline builds run on the Mac at $0. Rented GPU time goes only to parity preflight (cents) and training. Total budget is $3 [BRIEF §3].
- **Build time.** A full rebuild takes 10 min or less on the Mac. One view takes about 1 min [measured: data/views/build_reports/lidc-bundle-v3.json wall_s 49.4 with 6 procs].
- **VRAM.** The data store's resident share is 0.1 GB or less on the 32 GB RTX 5090 (the model and vLLM own the rest).
- **Host work per step.** No PNG decode, no fp32 pixel build on the CPU, and the data path adds no host syncs.
- **Safety.** Raw data is untrusted: readers run under `python -I`, and NPY loads use `allow_pickle=False` [code: store.py:80,85]. RIGA never goes to a rented host [docs/research/research2.0/Architecture.md §7.1].
- **One source.** Target polygonize, the checker and the scorer use `cuda/gpucheck.cpp` (built for CPU on the Mac, nvcc on the box) [code: CMakeLists.txt:1-16].

### Out of scope (explicit)
- New data sources or relabelling.
- Model, reward and trainer internals beyond the row/table interface.
- The 3D-voxel data model. Only forward-compatible keys are kept (section 8).
- CAD (parked).
- RIGA RL on the box (RIGA stays local).

### Assumptions (load-bearing, stated, not asked)
- A1. The exam that matters is the shifted LIDC test view (2,200 items) [measured: runs/r5/eval/eval_test].
- A2. All windows are scale 1 ("None of the RIGA or LIDC items needs f > 1") [code: src/objdel/targets/frame.py:14-16]. **This assumption carries the most weight:** it lets every window's doctor tables be derived from one per-item record.
- A3. RL views use op `id` and no CT jitter or noise (build_eval_views draws no jitter) [code: experiments/phase2/build_eval_views.py:52-66]. SFT bundles use D4 ops, jitter and noise [code: build_views.py:70-73].

---

## 2. Current state (read, then measured)

```
data/raw (26 GB RIGA, 177 MB LIDC XML) + data/lidc_slices (5.0 GB, 999 series npz)      [measured: du]
   | readers lidc.py / riga.py (P1-P6, F1-F7), splits.py (patient / near-dup groups)
   v
data/store  6 Parquet files, 5.08 GB, SHA-256 manifest                                    [measured: manifest.json]
   |-- prepare()        -> data/prepared: 9,030 PNG (LIDC 1.3 GB) + 9,635 cache npz (LIDC 33 MB)
   |-- build_views      -> SFT bundle v3: 98,804 rows, 4 epochs, 482.5 MB int16 contexts, 49 s
   |-- build_eval_views -> one fixed window per item and kind: PNG + caches per view
   |                        train-shifted 1.2 GB, train-rescaled 1.5 GB, val views 0.18-0.23 GB
   v
train_rl --train-views: grpo_rows_views = 8,208 centred + 8,208 shifted + 8,208 rescaled = 24,624
   (8,208 = 7,341 items + 289 disagreement items x 3 extra copies)  [measured: prepared manifest; code: dataset.py:240-253]
```

Strong points to keep: patient-level splits with no overlap (patients, series, items all 0 across splits) [measured: b1_analyze over docs/implementation2.0/notes/b1_out/geom.json]. Certified targets: 0 errors, exact share 1.0 on 98,804 rows [measured: lidc-bundle-v3 manifest]. Windows keep all rater pixels: inside = 1.0000 on every split and kind; crowded = 0 % [measured: geom.json]. Content-hash chains: store sha → prepared sha → view_hash / target_set_sha256 / reward_fingerprint.

Weak points found (data-centric):
- **W1. The coverage gap is still there after the fix.** The 3-view mix under-covers the failure cells. See section 5: small (<20 px) objects at offset 0.15-0.3 get 0.37-0.48 of their test share, and centred rows get 2.5-3.2x their share.
- **W2. Train windows are static.** Each item has one shifted and one rescaled window for every epoch, drawn with the exam's draw keys [code: build_eval_views.py:57]. This does not matter for 300 steps (9,600 prompts = 0.39 epoch of 24,624 rows). It does matter for multi-epoch runs.
- **W3. Each new view costs a full materialisation:** 1.2-1.5 GB of PNG + caches on disk, a bundle sync, and caches preloaded on the host (`serve.preload`).
- **W4. Device tables grow.** `ExpertCache.device()` uploads tables at first use and keeps them until garbage collection [code: src/objdel/reward/cache.py:241-258]. For all 3 train views that is about 1.23 GB [measured: S^2*(k+3)+4*S^2 summed over geom.json windows].
- **W5. The production targets use the numpy polygonizer** while the CUDA `check_masks` exists [BRIEF §2; code: cuda/gpucheck.cpp:1697]. The owner wants CUDA.
- **W6. True outline lengths above 1,200 tokens are uncounted** [BRIEF §2].
- **W7. No single dataset version.** A run lists roots and hashes are scattered over several manifests. `rl_shifted.yaml` already warns "check that lidc-test-shifted is listed in [the bundle]" [code: experiments/phase2/configs/rl_shifted.yaml:16], so a missing artifact is a real past risk.

---

## 3. Estimates (numbers decide the design)

| Quantity | Value | Source |
|---|---|---|
| LIDC items train / val / test | 7,341 / 1,084 / 2,200 (patients 557 / 78 / 158; nodules 1,247 / 200 / 401) | [measured: store manifest; geom.json] |
| RIGA items train / val / test | 592 / 63 / 95 (547 / 58 with a model image) | [measured: store + prepared manifests] |
| RL rows today | 24,624 = 3 x 8,208 | [measured: prepared manifest + code] |
| Prompts per step / step time | 32 / 23.6 s (R5, shared prefix) → **1.36 images/s** | [BRIEF §2; rl_shifted.yaml] |
| Prompts per 300-step run | 9,600 = **0.39 epoch** of 24,624 | [estimate from the above] |
| Context crops (int16), train | **476 MB raw**, median 33,488 px (about 183^2), max 111,456 px | [measured: npz headers of lidc-bundle-v3/context] |
| Context crops val + test | about 213 MB (scaled by item count) | [estimate] |
| Doctor masks, union box + 2 px, packed bits, train | **2.0 MB** | [measured: geom.json ext_px, k] |
| Full device tables, 3 train views, materialised | **1.23 GB** (median 26.9 KB, max 391 KB per window) | [measured: geom.json] |
| Per-step render output | 32 x 512 x 512 x 3 B = 25 MB uint8 | [estimate] |
| Per-step context transfer | 32 crops of 125^2 or less x 2 B ≤ 1 MB → about 40 µs at about 25 GB/s PCIe | [estimate] |
| Pillow bicubic upscale (S → 512, RGB) | 0.97 / 1.08 / 1.14 ms at S = 50 / 100 / 125 → about 35 ms per step on one core | [measured: Mac M3, synthetic arrays] |
| PNG decode per image (today's check) | about 1 ms with 3 threads | [code: dataset.py:32-33 docstring] |
| Spec table, RL (24,624 rows) | about 2 MB Parquet | [estimate: about 80 B/row] |
| SFT targets, 98,804 rows | 60.2 MB jsonl | [measured: build report rows_mb] |
| Upload to box today vs proposed | prepared 1.5 + views 2.7 GB vs contexts 0.69 GB + about 0.1 GB tables | [measured du; estimate] |
| Offline build cost | Mac CPU, minutes, $0. View build 49 s; train windows 2 min 11 s | [measured: build reports; rl_plateau note] |

**What the numbers force:**
1. Throughput is irrelevant (1.4 images/s). The data system must spend near zero **host time**, near zero **VRAM** and near zero **host syncs** per step. It must also make a new window kind **free** (no new GB on disk). The **gates** carry the quality.
2. With 476 MB of contexts against 2 MB of doctor bits, put the doctor data on the GPU and keep contexts in pinned host RAM. Copy about 1 MB per step. A resident context blob would buy nothing measurable.
3. Data is decoded once into the store, so GPU decode (nvImageCodec/DALI) does not pay here. NVIDIA's own docs say GPU decode can be slower than CPU decode on small images [web: docs.nvidia.com/cuda/nvimagecodec/samples/DICOM-pydicom.html].
4. Offline builds take minutes on the Mac at $0. Running them on a rented GPU costs money and gains nothing. GPU work pays **on the box, per step** (render + tables), and as the single CUDA source for polygonize and certify (parity).

---

## 4. High-level design

### 4.1 Layers (each content-addressed, each with a manifest)

```
 L0 RAW (local only, untrusted, immutable)            data/raw, data/lidc_slices
     | readers (python -I, no pickle)                   [exists]
 L1 ITEM STORE  items/{lidc,riga}/{split}.parquet       [exists] + column nodule_id (forward key)
     | certify (gpucheck CPU build on Mac == nvcc build on box)
 L2 DOCTOR RECORDS  doctors/{ds}/{split}.parquet        per item: union-box bit stack, invariants,
     |                                                   policy frame, certified policy target
     | draw (augment.draw / neighbour_side, CPU, seeded by SHA-256 keys)
 L3 VIEW SPECS  specs/{ds}/{split}/{kind}.parquet        windows only: (dx, dy, S, op, window, noise_seed)
     |                                                   + cell (size x offset) + weight
     | targets (SFT only): CUDA ring trace -> simplify -> certify
 L4 ROW SETS   rowsets/{name}.parquet                     RL: spec ids + mixture weights; SFT: spec + target
               exams/{ds}-{split}-{kind}/                 FROZEN: PNG + caches + records (as today)
     | audit (gates G1-G15) -> datacard.json
 L5 LOCK       data.lock.json                            version id = SHA-256 of all parent hashes,
                                                          code fingerprints, gate results
     | transport: private HF dataset repo (LIDC only), pinned revision; rsync fallback
 L6 GPU-FED STORE (on the box)                           contexts: pinned host blob (476 MB)
                                                          doctors: device-resident (about 2 MB + scalars)
                                                          per step: render + tables on GPU, freed after
```

### 4.2 Data flow (ASCII)

```
                   MAC (CPU, $0)                                       |        VAST BOX (RTX 5090)
                                                                       |
 raw LIDC/RIGA -> readers -> ITEM STORE ---------+                     |
                                |                |                     |
                                v                v                     |
                       DOCTOR RECORDS      CONTEXT BLOB (int16 HU,     |
                       (bits + invariants)  one flat file + index)     |
                                |                |                     |
                  draw specs <--+                |                     |
                      |                          |                     |
                      v                          |                     |
                 VIEW SPECS --(SFT)--> TARGETS   |                     |
                      |       (CUDA rings,       |                     |
                      |        CPU build)        |                     |
                      v                          |                     |
          ROW SETS + FROZEN EXAMS                |                     |
                      |                          |                     |
                      v                          v                     |
               AUDIT G1-G15  ----> DATACARD + data.lock.json ----------+--> download, verify SHA-256
                      |  (fail = no version)                           |       |
                                                                       |       v
                                                                       |  preflight parity G8/G9 (cents)
                                                                       |       |
                                                                       |       v
                                                                       |  sampler (weights by cell)
                                                                       |    | 32 spec ids
                                                                       |    v
                                                                       |  pinned ctx -> H2D (1 MB) -> RENDER
                                                                       |  kernel: window, op, bicubic -> uint8
                                                                       |    |                       |
                                                                       |    v                       v
                                                                       |  GPU pixel build      DOCTOR TABLES
                                                                       |  (Phase C2)           kernel: place bits
                                                                       |    |                  at offset, op,
                                                                       |    v                  union/core/vote,
                                                                       |  vLLM -> answers ---> EDT -> SCORER
                                                                       |                       (freed after step)
                                                                       |  results.json { data_version, ... }
```

### 4.3 Sequence (one version, one run)

```
Owner       Mac build (make data)        Audit          HF repo (private)      Box preflight        Trainer step loop
  |  make data VERSION=next  |              |                  |                     |                     |
  |------------------------->| L1->L2->L3->L4                  |                     |                     |
  |                          |------------->| G1-G15           |                     |                     |
  |                          |              |--fail--> stop, print failing gate      |                     |
  |                          |<--datacard + data.lock.json-----|                     |                     |
  |  review datacard (diff vs previous version)                |                     |                     |
  |  hf upload (LIDC only) --+------------------------------->| revision r           |                     |
  |  launch (owner-approved, NA host, preflight) -------------------------------->   |                     |
  |                          |              |                  |<-- hf download r ---|                     |
  |                          |              |                  |  verify every SHA vs lock (G1)            |
  |                          |              |                  |  parity: nvcc tables/targets/render       |
  |                          |              |                  |  == CPU-built references (G8, G9)         |
  |                          |              |                  |------ pass --------->| load contexts      |
  |                          |              |                  |                      | (pinned) + doctors |
  |                          |              |                  |                      | (device)           |
  |                          |              |                  |        each step: sample 32 spec ids ->   |
  |                          |              |                  |        render + tables -> vLLM -> score   |
  |                          |              |                  |        log spec ids + data_version        |
  |                          |              |                  |  final_eval on FROZEN exam (exam version)  |
  |<------------- results.json {data_version, exam_version, code, reward} -----------------------------|
```

### 4.4 Interfaces (concrete)

**Spec row (L3), primary key `spec_id`.** `spec_id = sha256(seed|kind|epoch|item_id|draw)[:16]`. Read by spec_id (sampler) and scanned by (split, kind) (audit).
```
spec_id str | dataset str | split str | kind {centred,shifted,rescaled,rotated} | epoch i16 | draw i16
item_id str | dx i32 | dy i32 | S i16 | op {id,r90,...} | win_lo i16? | win_hi i16? | noise_sigma f32 | noise_seed i64
off_euc f32 | ext_px i16 | cell str (size_bin x offset_bin) | weight f32 | fallback bool | kept_side bool
```
**Doctor record (L2), primary key `item_id`.** Read by item_id (table build) and scanned by split (audit).
```
item_id | group_id (patient) | nodule_id (series/cluster) | z | k | raters[k] | ub (x0,y0,x1,y1) | bits (k x h x w packed)
tau f64 | t5 f64 | min_area i32 | perimeter i32 | sigs[k] | voters[] | disagreement bool | copied{} | other_nodules[]
policy_S i16 | policy_box | policy_target str | fingerprint str (reward)
```
**Context index.** `item_id → (offset, x0, y0, h, w)` into `contexts.i16` (one flat little-endian file). The shape and dtype are in the manifest.

**Box API (Python, one module `objdel.data.gpustore`):**
```
store = GpuStore.open(lock_path)                 # verifies every SHA vs lock; refuses fingerprint mismatch
imgs  = store.render(spec_ids) -> uint8 (B,512,512,3) on device   # no host copy back
tabs  = store.tables(spec_ids) -> device handles for the scorer   # built now, freed by store.release(step)
meta  = store.rows(rowset) -> (spec_ids, weights)                 # for the sampler; tiny
```
RL dataset rows carry `key = spec_id` (no more `<i>::item_id` tags [code: dataset.py:240-253]).

**data.lock.json** (version identity):
```
{ "version": "<sha256>", "parents": {"items": sha, "doctors": sha, "contexts": sha, "specs": {...}, "targets": sha,
  "rowsets": {...}}, "exam": {"lidc-test-shifted": sha, ...}, "code": {"gpucheck": build_hash, "augment": module_sha,
  "frame": module_sha, "reward_fingerprint": fp}, "gates": {"G1": "pass", ..., "G10": {"min_ratio": 0.83, ...}},
  "built_on": "mac-cpu", "created": "<iso>" }
```
`version` = SHA-256 of the canonical JSON without `version` and `created`. Exams have their own `exam_version`, so a training-data change never moves the exam.

### 4.5 Key mechanisms

1. **Specs, not pixels, are the dataset.** A window is 10 numbers, so any number of fresh windows costs KB, not GB (this fixes W2 and W3). Pixels and tables are pure functions of (spec, doctor record, context) and are derived on the box.
2. **One doctor record per item gives every window's tables (A2).** Window tables (union, core, vote, distance map) are a pure function of the union-box bits placed at the window offset. The invariants (tau, t5, min_area, perimeter, signatures, voters) do not depend on the window while the union stays inside with 2 px clearance. Reasons: tau and t5 crop to the union first [code: cache.py:176-200 `_crop`]; min_area is taken over the union area; the distance maps are exact EDT to sets fully inside the window [code: cache.py:130-188]. **This is a hypothesis until gate G8 proves digest and score equality** on every existing train and exam window.
3. **Coverage-matched sampling.** The audit bins every spec by (object size, offset) and sets the mixture weights. Then the train-mix share in each cell reaches at least 0.8 of the **val** shifted share, never the test share (test statistics must not tune training). A larger spec pool (for example 4 shifted draws per item per epoch) makes the weights feasible (section 5).
4. **Frozen exams.** Exam PNGs, caches and records are rendered once and kept immutable. They are the golden reference for render and table parity, and they keep every past result comparable (test IoU 0.645 was measured on them).
5. **Single CUDA source.** Targets come from `check_masks` rings ("in polygonize's order" [code: cuda/gpucheck.h:106]), then the simplifier, then certify. On the Mac this uses the CPU build of the same file, on the box the nvcc build, and the gate compares `target_set_sha256`. After parity, the numpy polygonizer is archived (owner rule).
6. **Render parity, in stages.** Stage 1: render from pinned contexts with Pillow on the CPU, about 35 ms per step [measured above], which is byte-identical by construction. Stage 2: an integer CUDA port of Pillow's 8-bit fixed-point bicubic. Pillow uses 22-bit integer coefficients [code: frame.py:20-22 docstring], so an exact port is possible. It goes into gpucheck.cpp (one source), and G9 requires byte equality on every frozen exam image before Stage 1 is deleted.

---

## 5. Coverage audit (measured now: the gate this design adds would FAIL today's B1 mix)

Cells = object extent (px) x window offset (fraction of the side). Shares are % of items [measured: geom.json, seed-0 draws that equal build_eval_views].

```
TEST shifted (the exam), %            TRAIN 3-view mix (24,624 rows' windows), %
size\off  <.15 .15-.2 .2-.3  >.3      size\off  <.15 .15-.2 .2-.3  >.3
<10        2.1   2.0   4.7   9.3      <10        6.0   0.8   2.2   6.0
10-20      8.0   6.6  17.6  23.0      10-20     19.6   3.0   8.5  15.9
20-40      5.9   3.8   8.0   3.9      20-40     15.4   3.3   6.3   5.6
>=40       1.3   0.7   2.5   0.6      >=40       4.0   0.8   1.7   1.0

ratio train-mix / test                 ratio train-shifted-only / test
<10        2.79  0.37  0.47  0.65      <10        0.84  0.57  0.86  0.85
10-20      2.45  0.46  0.48  0.69      10-20      0.91  0.83  0.86  0.83
```

- Total variation distance (TV) to test: train mix 0.299, train shifted alone 0.114, **val shifted 0.074**. Val is a good stand-in for the exam.
- KS distance: offset (train shifted vs test) 0.059; object size (train vs test) 0.118. Train objects are larger: median extent 16 vs 14 px, median area 84 vs 64 px.
- Mixture weights centred:shifted:rescaled tried against val (min cell ratio over cells with share ≥ 1 %, TV): 1:1:1 → 0.48, 0.285; 1:2:1 → 0.58, 0.222; 0.5:2:1 → 0.66, 0.147; 0:2:1 → 0.77, 0.073. **No weighting of the three fixed windows reaches 0.8.** Size reweighting (val/train size marginal) does not reach it either (0.47). The cause: the one-draw windows themselves are sparse in the small x 0.15-0.2 cells.
- Design consequence: the spec pool needs more draws per item (cheap now) and cell-stratified sampling. The centred share becomes an explicit knob (for example about 20 % to keep the saturated centred skill; the model gate "centred ≥ 0.73" decides).
- The miss pattern this targets: test IoU 0.545 and 23 % misses at offset > 0.3. Off-centre objects of 0-40 px score 0.326 with 47 % misses [measured: docs/research/notes/rl_plateau_diagnosis_2026-11-08.md].
- Redundancy: 4 slices per nodule (median), p90 11, max 52. The top 10 % of nodules give 30.8 % of train items and 28.3 % of test items [measured: geom.json]. This is the same structure in both splits, so slice weighting stays. It is reported, and evaluation uses nodule-level bootstrap (Phase A1).

---

## 6. Quality-gate table

"Box" = rented RTX 5090 preflight (parity, cents). "Mac" = offline build ($0). A gate that fails produces no version.

| # | Metric | Check | Threshold | Where | Today |
|---|---|---|---|---|---|
| G1 | Artifact integrity | SHA-256 + row count of every file vs manifest and lock | 0 mismatches | Mac at build; Box at load, before any training | partial (per manifest) [code: store.py:171-182, dataset.py:290-301] |
| G2 | Determinism | Rebuild twice (sample of 500 items on every build, full on release) → same hashes | 0 diffs | Mac | holds for views [code: build_views.py:21] |
| G3 | Group leakage | patient, series, nodule_id, item_id overlap between splits | 0 | Mac | 0 / 0 / - / 0 [measured] |
| G4 | Duplicate leakage | exact hash of int16 slice pixels across splits; RIGA near-dup groups (NCC ≥ 0.45 or ORB inliers ≥ 22) crossing splits | 0 cross-split; within-split count reported | Mac | RIGA groups exist [code: splits.py:45-50]; LIDC pixel hash new |
| G5 | Selection leakage | coverage targets and checkpoint selection read **val** only; `final_eval` on test logs one read per model | 0 test reads outside final_eval | Mac (audit), Box (run log) | not enforced |
| G6 | Doctor masks inside window | per spec: every rater's pixel count inside = full count; clearance ≥ 2 px; no other outlined nodule | 100 % of specs | Mac at draw; Box re-check (device area sums) | 1.0000, crowded 0 % [measured] |
| G7 | Target certification | legal, ≤ MAX_VERTICES, signature kept, IoU ≥ 0.99, LIDC pixel-exact | 0 errors; exact share = 1.0 (LIDC) | Mac (gpucheck CPU build) | 0 errors, 1.0 [measured: bundle v3] |
| G8 | CPU↔CUDA parity | (a) target_set_sha256 CPU == nvcc; (b) GPU-built doctor table digest == CPU ExpertCache digest for every train + exam window; (c) identical scores on a fixed set of 2,000 saved answers | 0 diffs | Box preflight | new |
| G9 | Render parity | GPU-rendered bytes == Pillow reference (image_sha256) on all frozen exam images + 1,000 train specs | 0 diffs (Stage 2 only) | Box preflight | new |
| G10 | Coverage vs val | per (size x offset) cell with val share ≥ 1 %: train-mix share / val share; TV; KS on offset | min ratio ≥ 0.8; TV ≤ 0.10; KS ≤ 0.06 | Mac | **fails**: 0.48 / 0.285 [measured] |
| G11 | Side / kind coverage | each exam side S present at ≥ 0.5x its exam share; fallback windows; kept_side share | fallback = 0; kept_side ≤ 30 % (plan §5.0) | Mac | fallback 0, kept 0 [measured] |
| G12 | Outline length | tokens per certified target (per tokenizer): count > 1,200, count > max_completion_length, p99 | 0 RL reference targets over the limit; counts recorded in the lock | Mac (tokenizer only) | 0 of 23,469 > 1,536; > 1,200 unknown [BRIEF §2] |
| G13 | Label quality | per item: rater leave-one-out IoU; copied raters dropped (P5); abstaining raters (B2-4); disagreement items | report; block if outlier-rater share (LOO IoU < 0.2) moves > 1 pt vs the previous version | Mac | P5 / B2-4 exist [code: lidc.py:19-25; dataset.py:20-25] |
| G14 | Redundancy | slices per nodule; top-10 % nodule share; padded-window share | report; padded ≤ 5 % | Mac | 30.8 % / 28.3 %; padded 3.4 % train rescaled [measured] |
| G15 | Version stamp | run start: loaded hashes == lock; every result JSON has data_version + exam_version + reward fingerprint | missing or mismatch = run refuses | Box start-up | new |

---

## 7. Trade-offs (solves / worsens / change when)

| Decision | Solves | Worsens | Change when |
|---|---|---|---|
| **Specs are the dataset; pixels and tables derived on the box** | New windows cost KB; per-epoch redraws; cell-matched pools; upload 4.2 → about 0.8 GB | Needs table and render parity (G8, G9); a derivation bug hits every row | If parity cannot be reached: ship materialised tables for the RL rows only (today's path, about 1.2 GB) |
| **Doctor tables built per step on GPU from union-box bits (about 2 MB)** | Removes device-table growth (W4) and 9,635 npz files; any window is scorable | New GPU EDT kernel (brute force over boundary pixels, S ≤ 200) and parity work | A2 breaks (a window with f > 1) → cache those windows explicitly |
| **Contexts in pinned host RAM, not VRAM** | Saves 0.48-0.69 GB VRAM for vLLM and KV; costs about 40 µs per step of H2D | One async copy per step to manage on a side stream | VRAM shows free headroom after Phase C → make resident (one flag) |
| **Render Stage 1 Pillow CPU, then Stage 2 CUDA Pillow-exact** | Byte-identical from day one; all old hashes and results stay valid | Stage 1 keeps about 35 ms per step on the CPU; Stage 2 is custom code (no library is Pillow-exact) | If Stage 2 parity fails: re-freeze exams with a library GPU renderer (torchvision/DALI) as a new exam_version and re-baseline (about $0.25) |
| **Coverage target = val distribution, sampled by cell** | Fixes W1 without touching test; measurable gate | Uses val labels (object boxes) to shape training; val is then less clean for selection | If val becomes the selection signal for coverage weights: split val into val-select / val-shape by patient |
| **Frozen exams, separate exam_version** | Results stay comparable across training-data versions | Two version ids to track | Never for the headline exam; add new exam kinds as new frozen sets |
| **Version = content-addressed lock + private HF dataset repo (Xet dedup)** | Reuses the existing manifests; one library already in use; chunk dedup across versions [web: huggingface.co/docs/hub/xet/deduplication]; revision pinning | HF is a third party; LIDC terms must allow it (CC BY 3.0 per TCIA, owner to confirm); RIGA excluded | Team > 1, data > 100 GB, or branching data → DVC (now stewarded by lakeFS, Nov 2025) or lakeFS [web: lakefs.io/blog/lakefs-acquires-dvc] |
| **Parquet + one flat blob, not Lance/WebDataset** | Tables are KB-MB; pyarrow already pinned; byte-deterministic writers | No built-in row versioning (the lock does it) | 3D volumes or random access over > 10 GB → Lance (versioned columnar, fast random access) [web: lancedb.com/docs/overview/lance] |
| **Offline build on Mac CPU (gpucheck CPU build)** | $0; minutes; same source as CUDA | CPU and CUDA could drift → G8 needed | Build > 30 min (3D) → run builds on the box inside paid time |
| **No GPU decode (nvImageCodec/DALI)** | No extra dependency; decode happens once | None up to Phase E | 3D phase with compressed DICOM series → `dicom-series-to-volume` / nvImageCodec batch decode |
| **No scheduler (make DAG)** | Nothing to operate | No run history UI | Recurring rebuilds or several data owners |

---

## 8. Failure modes (what breaks, how it degrades, recovery)

| # | Failure | Detection | Degradation (stale beats error) | Recovery |
|---|---|---|---|---|
| F1 | GPU tables differ from CPU caches (EDT float, op handling) | G8 at preflight | Fall back to CPU-built caches for the rows in use (today's path) | Fix kernel; G8 again; no training on unproven tables |
| F2 | Render bytes differ | G9 | Stage 1 Pillow path (about 35 ms per step) | Fix kernel or re-freeze exams as a new exam_version |
| F3 | Missing or partial artifact on the box (the test view gap seen before) | G1 against the lock at start-up; Mac preflight lists the lock | Run refuses before training (minutes of paid time, not hours) | Re-download the pinned revision |
| F4 | Reward fingerprint changed, tables stale | Fingerprint in lock and digest; StaleCacheError [code: dataset.py:36-38] | Refuse | Rebuild L2 (minutes) |
| F5 | Coverage drift after a view-code change | G10 diff vs previous version | No version | Fix the draw or the weights |
| F6 | Test leakage through tuning | G5 | n/a | Re-derive weights from val |
| F7 | VRAM pressure from the data store | Device allocation log; data budget ≤ 0.1 GB | Contexts already in pinned RAM | Lower in-flight specs |
| F8 | Device-table growth over a long run (today's W4) | VRAM trend per step | n/a | Per-step build and free (design) |
| F9 | Sampler over-weights long nodules | G14 report | n/a | Optional per-nodule cap; nodule-level CIs |
| F10 | Lost box mid-run | Checkpoint every 50 steps; spec ids logged per step | Resume from the last checkpoint | Specs are deterministic, so the same rows replay |
| F11 | Corrupt or tampered download | SHA-256 vs lock (G1) | Refuse | Re-download |
| F12 | Untrusted raw data | Readers parse XML/NPY only, `allow_pickle=False`, `python -I` | n/a | n/a |

No single point of failure costs data: the Mac build is reproducible from L0 plus code, and the lock names every byte. The single point of failure for time is the rented box, handled by resume.

---

## 9. Web research (late 2026) that shaped decisions

- lakeFS took over the DVC open-source project from Iterative on 2025-11-18. DVC stays for smaller datasets; lakeFS targets enterprise [web: https://lakefs.io/blog/lakefs-acquires-dvc/]. For one owner with about 6 GB of derived data, both add a second version-control system without a gain over the existing manifests.
- HF Xet: chunk-level dedup into 64 MB content-addressed blocks; only changed chunks upload. One example cut two dataset versions by 53 % against Git LFS [web: https://huggingface.co/docs/hub/xet/deduplication]. This fits frequent spec and target versions over unchanged contexts.
- Lance: every write makes an immutable version; vendor-reported fast random access [web: https://lancedb.com/docs/overview/lance]. Kept for the 3D phase.
- Croissant + RAI 1.0: machine-readable dataset card (data life cycle, labelling, limitations). HF embeds Croissant JSON-LD [web: https://docs.mlcommons.org/croissant/docs/croissant-rai-spec.html; https://github.com/mlcommons/croissant]. `datacard.json` uses its field names (no tooling dependency).
- Leakage: slice-level splits inflate accuracy by 25-55 % on brain MRI [web: https://www.nature.com/articles/s41598-021-01681-w]. Subject-wise splits report 66-90 % against 95-99 % with leakage in a 2025 review [web: https://www.mdpi.com/2075-4418/15/18/2348]. This supports G3 at patient level (kept) and nodule-level CIs.
- Data-centric checklists: DC-Check (split rationale, subgroup evaluation, shift) [web: https://arxiv.org/pdf/2211.05764]. A 2026 audit checklist: datasheet, dedup + split overlap, distribution vs deployment, versioned audits stored next to evals [web: https://safeguard.sh/resources/blog/ai-data-quality-audit-checklist]. These map to G3/G4/G10 and the datacard.
- nvImageCodec: batch GPU decode helps large or compressed images. Small images can be slower than CPU [web: https://docs.nvidia.com/cuda/nvimagecodec/samples/DICOM-pydicom.html]. Hence no GPU decode here.
- cleanlab segmentation label-quality scores need model pixel probabilities [web: https://docs.cleanlab.ai/stable/tutorials/segmentation.html]. Deferred: G13 uses rater leave-one-out IoU, which needs no model. The pixel vote over the 8 policy answers could later act as probabilities.

---

## 10. Scale evolution (what breaks next)

- **Now (2D, 10^4 items):** the bottleneck is coverage quality (W1), not bytes. This design fixes it with specs + G10.
- **Next (3D voxels after Phase E):** per-nodule int16 sub-volumes, for example 64 x 128 x 128 x 2 B = 2 MB each. For 1,848 nodules that is about 3.7 GB [estimate]. Contexts no longer fit beside the model, so streaming is needed (Lance or Parquet shards, prefetch). DICOM series to volume (NVIDIA `dicom-series-to-volume` skill; nvImageCodec batch decode) starts to pay. The doctor record becomes a voxel bit stack per rater. Forward keys are already in place: `nodule_id`, `z` in L1/L2; `kind` can gain `volume`.
- **10x data (new datasets):** the Mac build stops being minutes, so move L2-L4 builds to the box in paid time, and the lock still proves equality.

---

## 11. Score and diagnosis (quality bar, 0-5)

| Dimension | Score | To reach 5 |
|---|---|---|
| Requirements | 5 | - |
| Numbers | 4 | Measure the box CPU's Pillow time and the H2D time (estimates now) |
| Trade-offs | 5 | - |
| Building blocks under stress | 4 | Prove A2 (window tables = f(record, offset)) on all windows. It is the base of the design and still a hypothesis |
| Failure | 4 | Rehearse F1 and F3 once on a box (preflight script) |
| Interfaces | 4 | Write the GpuStore ↔ Phase C2 pixel-build contract (dtype, layout, stream) with the Phase C owner |
| Collaboration | 3 | Owner decisions needed: centred share; val-shaping split; HF hosting of LIDC derivatives |

**Weakest dimension:** an unproven base (A2 + G8). If GPU-built tables do not match CPU caches bit for bit, the derivation idea falls back to materialised caches. The rest still stands: specs, gates, lock, coverage sampling. **Current bottleneck:** coverage (G10 fails at 0.48 against 0.8). **At the next order of magnitude (3D),** contexts leave RAM and I/O becomes the bottleneck.

---

## 12. Ten checks any implementation plan must pass

1. **Coverage gate before spend.** The plan runs G10 on the exact RL row set it will train. It shows min cell ratio ≥ 0.8 against **val** (not test), TV ≤ 0.10 and KS ≤ 0.06, plus the before/after cell table, before any GPU rental.
2. **Version stamp.** Every training and eval command takes one `data.lock.json`. The run refuses on any hash or fingerprint mismatch. Every result JSON carries `data_version`, `exam_version`, the reward fingerprint and the gpucheck build hash.
3. **CPU↔CUDA parity.** One preflight step on the box proves `target_set_sha256` equality, doctor-table digest equality on every RL and exam window, and identical scores on a fixed answer set. It runs before training steps, and it costs cents.
4. **Leakage zero.** The plan shows patient, series, nodule and item overlap = 0. It shows cross-split exact pixel duplicates = 0 (LIDC) and no crossing RIGA near-duplicate group. Test is read only by `final_eval`, and coverage weights come from val.
5. **Doctor masks inside.** 100 % of specs keep every rater's full pixel count with ≥ 2 px clearance and no other outlined nodule. This is checked at draw time and re-checked on the device.
6. **Exams frozen.** The headline exam (lidc-test-shifted) bytes and caches do not change. Any renderer or table change proves byte or score equality on it, or declares a new exam_version and re-baselines.
7. **Outline-length count.** The plan reports per-tokenizer counts of certified targets above 1,200 and above the completion limit. Any limit change cites them.
8. **No silent materialisation growth.** The RL path has no per-window PNG or cache files. VRAM for data stays ≤ 0.1 GB and does not grow over steps (log it). Old prepared/views training paths are archived, not left beside the new ones.
9. **One CUDA source for geometry.** Production targets use `check_masks` rings (CPU build on the Mac, nvcc on the box) with a parity test against the numpy polygonizer on all 23,469 policy targets and the bundle. Then the numpy polygonizer is archived. Tests live in `tests/`, and `make test` / `make lint` pass.
10. **Determinism and cost.** Two rebuilds give the same version id. The full Mac rebuild takes ≤ 10 min. The plan states its GPU dollars (preflight + run) within the $3 budget, North-America hosts only, and preflight before every rental.
