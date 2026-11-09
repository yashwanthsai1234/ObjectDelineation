# Architecture: the new system

Date: 2026-11-08. Language: ASD-STE100 Simplified Technical English.

This report gives the design of the new system. The design is ready for production use, and it can grow from one GPU
to many, and from 2D outlines to 3D voxel volumes and robot assets. `Direction.md` gives the order of the work. `Research2.0.md` gives
the facts.

---

## 0. Words in this report

We use the same words as the other reports: model, image, outline, doctor, checker, reward, answer, group, step,
miss, offset, ceiling, GPU, verifier. We add these words:

| Word | Meaning |
|---|---|
| layer | One part of the code with one job. A layer uses only the layers below it. |
| engine | The CUDA geometry code. It is the lowest layer. |
| shape kind | One kind of shape that the verifier knows: 2D outline, 3D voxel volume, closed mesh. |
| doctor store | The doctor masks and the tables made from them. They stay in GPU memory. |
| generation | The writing of answers by vLLM. |
| warp | A team of 32 GPU threads that run the same instruction together. |
| ballot, popcount | GPU instructions: ballot packs 32 yes/no values into one 32-bit word; popcount counts the 1 bits. |
| Euler number | Parts minus holes. The engine uses it to count holes. |
| medoid | Of 8 answers, the answer with the highest mean IoU against the other 7. |
| GED, HM-IoU, HD95 | Measures from papers: two compare a set of answers with a set of doctors; HD95 measures the worst edge distance. |
| CISPO, GSPO | Two RL update methods that stay stable when generation and training run at the same time. |
| FSDP2, Megatron | Tools that split one model across many GPUs. |
| LoRA | A small set of extra weights that RL changes. The main model weights stay fixed. |
| trainer | The program that calculates the loss and updates the adapter. |
| run specification | One typed file that tells a run what to do. |

---

## 1. What the system must do

### 1.1 Functions

1. Train a model to draw shapes as text.
2. Score each answer exactly against a set of doctors.
3. Test a trained model with fair and standard metrics.
4. Offer the verifier to other programs: other trainers, other models, data filters.
5. Grow to 3D voxel volumes and robot assets with the same verifier.

### 1.2 Qualities

| Quality | Target | How we measure it |
|---|---|---|
| Exact | CPU and GPU give the same bits for every score. | Golden files and oracle tests. |
| GPU first | The GPU is busy for 80 % or more of each step. | NVTX time marks and Nsight Systems. |
| Fast | 16-18 s for each step on one RTX 5090 (12-14 s later). | Step logs. |
| Reproducible | Each result names its code version, reward version and data version. | Run records. |
| Framework-free | The verifier and the data work with TRL, verl or NeMo-RL. | Two trainers pass the same reward test. |
| Safe with money | No run starts without a check and the owner's approval. | Preflight and the host policy. |

---

## 2. The layers

Each layer uses only the layers below it: upper layers call lower layers, never the reverse. A tool checks this rule
in `make lint`.

```
+---------------------------------------------------------------------------+
|  experiments/        thin command-line tools: read arguments, call library |
+---------------------------------------------------------------------------+
|  infra               rent hosts, send data, run, pull, stop; cost guard    |
+--------------------------------+------------------------------------------+
|  rl                            |  sft                                     |
|  build trainer, TRL adapter,   |  supervised training state, resume       |
|  pure forward functions,       |                                          |
|  reward adapter                |                                          |
+--------------------------------+------------------------------------------+
|  policy       model specifications, prompt, sampler (vLLM and HF),        |
|               model caches (vision features, positions)                    |
+---------------------------------------------------------------------------+
|  data         store, LIDC and RIGA readers, splits, views, targets        |
+---------------------------------------------------------------------------+
|  verify       the verifier: reward version, doctor store, score, metrics, |
|               statistics, service.  NO torch, NO TRL.                      |
+---------------------------------------------------------------------------+
|  geom         the engine: CUDA checker, polygonizer, labels, topology,    |
|               distance; Python binding; CPU build for tests                |
+---------------------------------------------------------------------------+
```

### 2.1 What each layer owns

| Layer | Owns | Does not own |
|---|---|---|
| geom | All exact geometry, on CPU and GPU, from one source file. | Any reward weight or any model. |
| verify | The reward version (one source for weights and limits), the doctor store, batch scoring, evaluation metrics, confidence intervals. | Torch, TRL, vLLM. |
| data | Raw data readers, the Parquet store, group-aware splits, shifted and rescaled views, certified targets. | Training. |
| policy | What a model is: name, processor settings, prompt, LoRA targets, sampler. | The loss. |
| rl | The RL algorithm and the adapter to one framework. | Reward maths. |
| sft | Supervised training. | Reward maths. |
| infra | Cloud hosts and money. | Library code. |

### 2.2 Why this split

- The verifier is the part with other users. It must not need torch or TRL.
- The trainer adapter is the part that breaks when TRL changes. It must be small and in one place.
- Evaluation must not need a trainer. Today it builds a full trainer only to generate answers.

---

## 3. One RL step in the new system

### 3.1 Data flow

```
                          GPU (RTX 5090, 32 GB)
 +------------------------------------------------------------------------------+
 |                                                                              |
 |  32 images (uint8) --> GPU pixel build --> vision encoder (once per image)   |
 |        ^                                         |                           |
 |        |                                         v                           |
 |  [CPU: read 32 PNG]                      vision features (cached)            |
 |                                                  |                           |
 |  prompt tokens + positions (calculated once per run, kept on the GPU)        |
 |                                                  |                           |
 |                                                  v                           |
 |  vLLM generation: 32 requests x 8 answers, CUDA graphs, answer limit ~1,200  |
 |                                                  |                           |
 |              +-----------------------------------+---------------+           |
 |              v                                                   v           |
 |  verifier: parse, raster, score          trainer: packed forward, one        |
 |  against doctor store (side stream)      prompt per group, loss, backward    |
 |              |                                                   ^           |
 |              +-------------- rewards (256 numbers) --------------+           |
 |                                                                  |           |
 |  adapter update --> in-place weight push into vLLM (LoRA layers only)        |
 |                                                                              |
 +------------------------------------------------------------------------------+
   Only these cross to the CPU: 32 image files in, 256 answers as text, log numbers.
```

### 3.2 Where each part runs

| Part | Today | New |
|---|---|---|
| Image file read | CPU | CPU (no GPU decoder exists for PNG) |
| Pixel build and resize | CPU, 1.61 GB fp32, twice | GPU, from 32 uint8 images |
| Prompt tokens and positions | CPU, 256 rows, twice | Once for each run, on the GPU |
| Generation | GPU (vLLM) | GPU (vLLM), 8 answers in one request |
| Parse of the answer text | CPU, milliseconds, hidden | CPU (no gain from a move) |
| Raster and score | GPU (checker) | GPU (verifier, same engine) |
| Reference pass | GPU | GPU, or removed if the KL weight is 0 |
| Update | GPU | GPU, micro-batches of 32 x 8 |
| Weight push | about 700 calls, all weights | In place, LoRA layers only, no vision weights |
| Polygonizer for targets | CPU (numpy) | GPU (CUDA) |
| Evaluation metrics | CPU pool | GPU |

### 3.3 Who talks to whom in one step

```
Data loader   Policy cache    vLLM         Verifier       Trainer       vLLM weights
    |              |            |              |              |              |
    |--32 images-->|            |              |              |              |
    |              |--prompts,--|              |              |              |
    |              |  features->|              |              |              |
    |              |            |--256 answers>|              |              |
    |              |            |--256 ids, log-probs-------->|              |
    |              |            |              |--256 rewards>|              |
    |              |            |              |              |--forward,    |
    |              |            |              |              |  backward    |
    |              |            |              |              |--LoRA push-->|
    |              |            |<-------------------------------ready-------|
```

### 3.4 Memory on one GPU

**Worked example 1: does it fit?** The full RL run measured about 11.2 GB for vLLM and a trainer peak of about 17 GB.
Today vLLM sleeps during the update and gives back most of its memory. If vLLM stays awake, both peaks occur at the
same time: 11.2 + 17 = about 28 GB of 32 GB. That leaves about 4 GB. Micro-batches of 32 x 8 need more memory than
16 x 16. Thus we test "vLLM stays awake" and "larger micro-batches" one at a time, and we measure the peak each
time.

---

## 4. The verifier

### 4.1 Parts

```
+-------------------+     +---------------------+     +----------------------+
| client            |     | batch joiner        |     | engine (CUDA)        |
| TRL / verl /      | --> | joins requests,     | --> | parse, raster,       |
| NeMo-RL / filter  |     | one GPU call each   |     | validity, labels,    |
+-------------------+     | batch               |     | topology, score      |
         ^                +---------------------+     +----------+-----------+
         |                                                       |
         |                +---------------------+                |
         +--- results --- | result packer       | <--------------+
                          | terms, reward,      |      doctor store (GPU):
                          | versions            |      doctor bit planes,
                          +---------------------+      union, core, vote,
                                                       distance tables
```

### 4.2 Contract

Inside one process, the trainer calls a Python function. Across processes, the same call goes over HTTP. Both give the
same bits.

```
POST /v1/score
{ "spec": "reward_v5",
  "task": "outline2d",
  "items": [ { "id": "a_000123",
               "key": "lidc/train/0042::shift",
               "W": 512, "H": 512,
               "names": ["nodule"],
               "text": "{...outline...}",
               "truncated": false } ] }

200
{ "spec": "reward_v5", "checker": "<code version>", "cache": "<data version>",
  "results": [ { "id": "a_000123", "reward": 0.73,
                 "terms": { "format_ok": 1, "legal": 1, "iou_tc": 0.71,
                            "f5": 0.90, "location": 0.0, "topo_any": 1 },
                 "error": "" } ] }

Errors: unknown_key, spec_mismatch, batch_too_large, device_error.
Other calls: GET /v1/spec, GET /healthz, GET /metrics.
```

Rules:

- Scoring is a pure function of the reward version, the data and the text. A repeated request gives the same result.
- A new reward is a new name, for example `reward_v5`. One server can hold two versions for a comparison.
- Each result names the code version and the data version.
- Doctor masks never leave the GPU. Only numbers come back.

**Worked example 2: the load on the verifier.** One step sends 256 answers of about 1 KB each. One step takes about
18 s. Thus the load is 256 KB / 18 s = about 14 KB each second. At 13 s for each step, one GPU makes about
256 / 13 = 20 answers each second. With 64 GPUs in async mode at this rate, the load is about 1,300 answers each second, or about 1.3 MB each second. Plain HTTP with
JSON carries this. We do not need a stream protocol.

### 4.3 Shape kinds

| Shape kind | Input | Exact outputs | State |
|---|---|---|---|
| outline2d | Rings and holes on a grid | Legality (13 defect types), IoU against the set, boundary F, area, b0, b1 | Exists |
| volume3d | A stack of slice outlines, slice thickness | 3D labels (26 and 6 neighbours), b0, b1, b2, IoU against the set | CPU only. GPU in Phase G. |
| mesh3d | Closed triangle mesh | Watertight, manifold, self-intersection, genus, volume IoU interval | Phase H and I |

### 4.4 Why the results are exact on any GPU

Every output is an integer, or a correctly rounded operation on integers. The engine uses no floating-point atomic
sums. Thus the order of the threads does not change the result. Any correct parallel algorithm gives the same bits.

---

## 5. The engine

### 5.1 Inside the engine

```
 outlines --> [raster: 1 warp per row, ballot] --> bit words --+--> [validity: grid pairs,
                                                               |     warp parity] --> defects
 masks ------------------------------------------> bit words --+
                                                               +--> [labels: runs in bit words,
                                                               |     8 and 4 neighbours] --> parts
                                                               |        |
                                                               |        +--> [Euler number per part] --> holes
                                                               |        +--> [polygonizer: true edges only,
                                                               |              ring start from labels,
                                                               |              corner order by ranking] --> rings
                                                               +--> [score: popcount planes,
                                                                     distance tables] --> terms
```

### 5.2 Engine changes

| Part | Today | New | Gain |
|---|---|---|---|
| Polygonizer | Visits 4(W+1)(H+1) edge places. About 40 jump rounds. 3 host waits. Not used in production. | One block for each mask. Emits only true edges. Ring start in closed form. One host wait. Used by the target builder. | About 500 times less work on a 256 x 256 mask. Obeys the CUDA rule. |
| Labels | One thread for each row, byte pixels, always 4 passes. | Runs in bit words. Skip passes when nothing is small. Holes from the Euler number of each part. | 5 to 20 times faster (estimate) |
| Validity totals | One thread for each shape. | One warp for each shape. Grid for candidate pairs. | No single slow shape stops a batch. |
| Memory | Grows with `cudaMalloc` and `cudaFree`, which stop the GPU. | Pool memory. Upload expert tables in one batch. | No device-wide stop |
| Stream priority | Same as the trainer. | Truly lowest. | Trainer is not delayed. |

### 5.3 The polygonizer in detail

**Worked example 3: a 3 x 3 square with a 1-pixel hole.** The square has 8 object pixels around 1 hole pixel.

```
  x=0   1   2   3
 y=0 +---+---+---+      # = object pixel
     | # | # | # |      . = hole pixel
   1 +---+---+---+
     | # | . | # |      outer ring: (0,0) -> east -> (3,0) -> (3,3) -> (0,3) -> (0,0)
   2 +---+---+---+      hole ring:  (1,1) -> south -> (1,2) -> (2,2) -> (2,1) -> (1,1)
     | # | # | # |
   3 +---+---+---+
```

- Edges: the outer ring has 12 unit edges. The hole ring has 4 unit edges. The total is 16 edges.
- The labels give 1 object part and 2 background parts (outside and the hole).
- The outer ring starts at the top-left corner of the first object pixel and goes east.
- The hole ring starts at the top-left corner of the hole pixel and goes south.
- Each edge finds its next edge from the 2 x 2 pixel pattern at its end. At a diagonal touch, the edge turns left.
- A ranking pass puts the corners in order along each ring. The outer ring has 4 corners. The hole ring has 4
  corners.
- The signed area of the outer ring is +9. The signed area of the hole is -1. The object area is 8. This matches the
  8 object pixels.

The old method finds each ring by repeated pointer jumps over all edge places. The new method takes the ring start
from the labels, which the scorer already calculates.

---

## 6. The training system

### 6.1 Parts

```
+--------------------+      +---------------------+      +--------------------+
| run specification  | ---> | rl.build_trainer    | ---> | framework adapter  |
| (one typed file)   |      | model, data, reward |      | TRL now, verl later|
+--------------------+      +---------------------+      +---------+----------+
                                                                    |
            +-------------------------+-----------------------------+
            v                         v                             v
+--------------------+   +-----------------------+   +--------------------------+
| policy.sampler     |   | verify (reward)       |   | rl.forward (pure)        |
| vLLM / HF          |   | in process or HTTP    |   | packing, shared prompt,  |
+--------------------+   +-----------------------+   | ratio, advantage         |
                                                     +--------------------------+
```

### 6.2 Rules for the TRL adapter

- One file holds all replacements of private TRL names. Today there are 9.
- A start-up test checks that each name still exists. If one name is missing, the run stops with a clear message.
- TRL stays pinned at 1.17.1. An update is a planned job with the full equivalence tests.

### 6.3 The update rule

```
reward     = legal x ( IoU + 0.1 x location ) + topology
advantage  = reward - mean of the other 7 rewards in the group
loss       = Dr. GRPO token loss (no length bias)
correction = cut each token's ratio at 2
             (other choice: FP16 numbers in the trainer and in vLLM)
async (2+ GPUs) = CISPO or GSPO, plus the correction
```

### 6.4 How the system grows

```
 Stage 1: 1 GPU, in turns                 Stage 2: 2 GPUs, at the same time
 +-----------------------------+          +--------------+    LoRA     +--------------+
 | vLLM | trainer | verifier   |          | GPU 0: vLLM  | <---------- | GPU 1:       |
 | (in turns on one card)      |          | generation   | ----------> | trainer +    |
 +-----------------------------+          +--------------+   answers   | verifier     |
  16-18 s per step (estimate)              13-14 s per step (estimate) +--------------+

 Stage 3: 8+ GPUs, many generation GPUs
 +---------+ +---------+ +---------+        +-----------------+     +-----------------+
 | vLLM    | | vLLM    | | vLLM    |  ...   | trainer (FSDP2  | --> | verifier        |
 | GPU     | | GPU     | | GPU     | -----> | or Megatron)    |     | service (HTTP)  |
 +---------+ +---------+ +---------+        +-----------------+     +-----------------+
  framework: verl (first choice) or NeMo-RL; recipe: ScaleRL (CISPO, batch-level scale)
```

| Stage | GPUs | Framework | Weight push | Gate to the next stage |
|---|---|---|---|---|
| 1 | 1 | TRL + our adapter | In place, LoRA only | Phase C done; more than 30 % of the step is still generation |
| 2 | 2 | verl async | LoRA only, CUDA IPC or NCCL | Same reward curve as stage 1; time for each batch 60 % or less |
| 3 | 8+ | verl or NeMo-RL | NCCL, sparse delta | A larger model or a 3D shape kind needs it |

---

## 7. Data and evaluation

### 7.1 Data

```
raw LIDC / RIGA --> readers --> Parquet store --> group-aware splits
                                     |
                                     +--> views: centred, shifted, rescaled
                                     |        (each view has its own doctor masks)
                                     +--> targets: certified outlines (GPU polygonizer)
                                     +--> doctor store (GPU tables)
```

Rules:

- Downloaded data is untrusted. Python that reads it runs in isolated mode.
- RIGA data never goes to a rented host.
- RL uses all three views. The test uses the shifted view.

### 7.2 Evaluation

| Metric | Why | Source of the answers |
|---|---|---|
| IoU against the best matching doctor | Our main metric | 8 answers |
| Dice against the doctors' majority | Other papers use it | 8 answers |
| GED, HM-IoU | Standard for many doctors | 8 answers as a distribution |
| Boundary F, boundary IoU, HD95 | Edge quality | 8 answers |
| Legal answers | Safety | 8 answers |
| Fair ceiling | When a nodule has k doctors, each doctor is compared with the other k - 1. The model must also be compared with k - 1 doctors, one left out each time. | Doctor masks |
| Confidence interval | Resample whole nodules, not single images | All of the above |

The final answer for each image comes from a label-free choice: keep legal answers, then take the medoid of the 8.

---

## 8. Measurement and logs

### 8.1 Time marks in each step

```
step
 |-- load
 |-- generate
 |     |-- push weights
 |     |-- decode
 |     |-- tail (1-2 answers alone)
 |-- prepare (must be near 0 s on the CPU)
 |-- reference
 |-- update
 |     |-- micro-batch 1 .. n
 |-- score (side stream)
```

### 8.2 Numbers in each step log

| Number | Why | Alarm |
|---|---|---|
| Time of each phase above | Find idle time | prepare above 0.5 s |
| Mean, lowest and share-at-0 of the correction ratio | Find lost signal | share at 0 above 0.5 % |
| Share of groups with no spread | Find dead groups | above 25 % |
| Entropy for each view | Find collapse on hard images | below 0.03 on shifted images |
| IoU and misses for each offset bucket | Find the shift problem | far misses above 20 % |
| Longest answer, share at the limit | Find the long tail | above 5 % |
| Reward version, code version, data version | Reproduce the run | missing value stops the run |

---

## 9. Tests and quality

| Test | What it proves |
|---|---|
| Oracle tests (numpy, exact integer) | The engine is correct on the CPU. |
| Golden files (195,809 recorded cases) | The engine has not changed. |
| GPU job (build with nvcc, run all checker tests) | The GPU raster, a separate program, is correct. |
| Sanitizer (memcheck, initcheck) | No memory errors. 0 errors on 40 GPU tests today. |
| TRL contract test | All 9 replaced TRL names exist. |
| Equivalence tests for the fast trainer | Each speed change keeps the same loss and gradient. |
| Import rule check in `make lint` | No layer imports a layer above it. |
| Reward A/B test through two trainers | The verifier gives the same bits to TRL and verl. |

---

## 10. Operations

```
Owner           Preflight          Host policy        Rented GPU         Ledger
  |                 |                   |                  |                |
  |--plan, cost---->|                   |                  |                |
  |                 |--data, checker,   |                  |                |
  |                 |  code tests OK--->|                  |                |
  |<--GPU hours, $--|                   |                  |                |
  |--approve----------------------------->|                |                |
  |                 |                   |--US/CA host only>|                |
  |                 |                   |--hard stop $---->|                |
  |                 |                   |                  |--run, pull---->|
  |                 |                   |                  |--teardown----->|
  |<--results, cost-----------------------------------------------------cost|
```

---

## 11. Decision records (short)

| # | Decision | Status | Reason |
|---|---|---|---|
| D1 | Keep TRL plus our adapter on one GPU. Move to verl at 2 GPUs. | Proposed | No framework is faster on one GPU. TRL async cannot use images. |
| D2 | The verifier is a separate layer with no torch, and later a service. | Proposed | Other programs need it. Only numbers cross the boundary. |
| D3 | The production polygonizer is the CUDA polygonizer, rebuilt to visit only true edges. | Proposed | Owner rule. About 500 times less work. |
| D4 | Ring start in closed form from component labels. | Proposed | No repeated jumps: about 10 times less work for a ring of 1,000 edges. Must pass the golden files first. |
| D5 | Reward version 5: legality x (IoU + 0.1 x location) + topology; no division by the spread. | Proposed | Signal for misses. No change of objective in the middle of a run. |
| D6 | Fix the number mismatch: token-level truncation, or FP16. | Proposed | About half of each step's signal is lost today. |
| D7 | Next model: Qwen3-VL-4B, not Qwen3.5. | Proposed | +3.4 RefCOCO, same structure. Qwen3.5 breaks the trainer's shortcuts. |
| D8 | One modular program, not many services. | Proposed | Only the verifier has a second user. |
| D9 | 3D by boundary faces and run labels, not marching cubes. | Proposed | Exact, and no ambiguous cases. |
| D10 | CAD solids as a shape kind. | Parked by the owner (2026-11-08) | The main direction is 3D voxels. |

---

## Summary

The new system has one exact GPU engine at the bottom and a verifier above it that needs no training framework. The
model, the data and the trainer sit above the verifier. Data stays on the GPU, and only short texts and numbers cross
to the CPU. The same verifier serves 2D now and 3D voxel volumes and robot assets later, and the trainer can move from one GPU
to many without a change to the reward.
