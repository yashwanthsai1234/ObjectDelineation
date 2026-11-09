# Previous against new

Date: 2026-11-08. Language: ASD-STE100 Simplified Technical English.

This report compares the previous system and plan with the new system and plan. "Previous" means the repo on
2026-11-08, before Research 2.0. "New" means the design in `Architecture.md` and the plan in `Direction.md`. Values
in the "New" column are targets or estimates until a run measures them.

---

## 0. Words in this report

We use the same words as the other reports: model, image, outline, doctor, checker, reward, answer, group, step,
miss, offset, ceiling, GPU, verifier, layer, engine, task.

---

## 1. The story in short

The project planned a long study with many RL arms on large GPUs. The real project took a different road. Two small
models failed the entry tests. The owner then chose one model, Qwen3-VL-2B, and one RL run on one RTX 5090. The team
built a fast CUDA checker and a fast trainer. The RL run lifted IoU on 128 held-out images from 0.511 to 0.659 for about $4.56.

Research 2.0 then found three causes for the limit: centred training images, a number mismatch that drops about half
of the learning signal, and a reward that gives no signal to misses. It also found about 40 % idle GPU time in each
step, and one CPU polygonizer that breaks the CUDA rule.

```
 PREVIOUS                                     NEW
 +------------------------------+             +------------------------------+
 | plan: 9 phases, ~150 GPU-h,  |             | plan: 9 phases with gates,   |
 | many RL arms, A100/H100      |             | each run approved, 1 GPU     |
 +--------------+---------------+             +--------------+---------------+
                |                                            |
                v                                            v
 +------------------------------+             +------------------------------+
 | what ran: 1 model, 1 RL run, |   Research  | measure fairly, fix 3 causes,|
 | centred images, $16.90 total | ----2.0---> | cut idle time, CUDA          |
 | test IoU 0.645               |             | everywhere, then 3D voxels   |
 +------------------------------+             +------------------------------+
```

---

## 2. Plan against reality against new plan

| Item | First plan (9 phases) | What happened | New plan |
|---|---|---|---|
| Model search | SmolVLM, then LFM, then Qwen; stop at the first pass | LFM2.5-VL-450M failed (IoU 0.495 shifted). Qwen3-VL-2B also failed the gate (median 0.583 against 0.60). The owner approved RL anyway. | Keep Qwen3-VL-2B. Add Qwen3-VL-4B as a scaling point. |
| RL arms | 9 arms (A1-A9), best-of-8, DPO, 3 seeds, about 50 GPU-hours | 1 run with the full reward, 770 steps | 3 small runs with gates (B1, B2, B3), then baselines |
| GPUs | A100 or H100, about 150 GPU-hours | RTX 4090 for SFT, RTX 5090 for RL, $16.90 total | RTX 5090; a second GPU only after Phase C |
| Checker | CPU first, CUDA "later" | CUDA built at once, bit-exact, all scoring on the GPU | Same engine, work-optimal, CUDA polygonizer in production |
| 3D | Phase 7 | CPU volume code only | GPU 3D verifier in Phase G |
| 3D voxels, Physical AI | 3D in the first plan | Not started | Phases G and H (3D voxels), Phase I (robot assets). CAD: studied, then parked by the owner |

---

## 3. Quality

| Measurement | Previous | New (target, estimate) | What changes it |
|---|---|---|---|
| Test IoU | 0.645 | 0.70-0.73 | Shifted images, number fix, reward version 5 |
| Test boundary F | 0.816 | 0.86-0.89 | Fewer misses. Found nodules already score 0.896 (validation split). |
| Legal answers | 99.4 % | about 100 % | Keep only legal answers at test time |
| Misses at offset above 0.3 | 23 % | below 10 % | Shifted training images |
| Ceiling | Favours the model (k against k - 1) | Fair (k - 1 against k - 1) | Fair ceiling in Phase A |
| Confidence intervals | Resample single images | Resample whole nodules | Phase A fix |
| Metrics reported | IoU, F, precision, recall, Dice, legal | Plus GED, HM-IoU, boundary IoU, HD95 | 8 answers used as a distribution |
| Baselines on our splits | MedSAM only (Phase 0) | nnU-Net 2D, MedSAM2, SAM 3, SimpleSeg | Phase F2 |

**Worked example 1: what the shifted images change.** On test, 825 images have offset above 0.3. Today 23 % of them
are misses: about 190 images. The target is below 10 %: fewer than 83 images. Thus about 107 more nodules are found.
A found nodule scores about 0.70 and a miss about 0.05. The gain is 107 x 0.65 / 2,200 = about +0.032 test IoU from
the far images alone. The B1 target (0.69-0.71) also counts gains at other offsets.

---

## 4. The reward and the update

| Part | Previous | New |
|---|---|---|
| Reward terms | Format gate, legality, IoU against the best doctor, boundary F, area, topology | Legality x (IoU + 0.1 x location for misses) + topology. Boundary F and area are limits. |
| Scaling | Each term divided by its group spread, with floors | Subtract the mean of the other 7 answers. No division. |
| Signal for a group of 8 misses | About 0 | The answer nearest to the nodule gets a positive advantage |
| Objective during the run | Changes when the floors start to control | Fixed for the whole run |
| Correction ratio | Mean 0.70-0.85, 0 in 98.1 % of steps | Near 1 (token truncation or FP16) |
| Answers that teach in a step | About 40-55 % of 256 (estimate) | About 100 % (target) |
| Hard images | Same rate as easy images | Shown more often, with 16 answers |
| Diversity option | None | GED reward (strictly proper), for a model that gives a range of outlines |

**Worked example 2: a group of 8 misses.**

| Answer | Distance to nodule centre | Previous reward | New location term (radius 7 px) |
|---|---|---|---|
| 1 | 5 px | about 0.02 | 0.78 |
| 2 | 10 px | about 0.02 | 0.36 |
| 3 to 8 | 20 px | about 0.02 | 0.02 |

Previous: all 8 rewards are almost equal, so the advantages are about 0. New: answer 1 gets an advantage of about
+0.1 x (0.78 - 0.07) = +0.07 (0.07 is the mean of the other 7 location terms), and the model learns to move toward the nodule.

---

## 5. Speed and GPU use

| Part of the step | Previous (measured) | New (estimate) | Change |
|---|---|---|---|
| Weight push, wake, sleep | 1.7 s, about 700 calls, all weights | about 0.05 s, in place, LoRA layers only | C1 |
| CPU work before the reference pass | 3.0-3.2 s | about 1-2 s | C2, C4 |
| Long tail of generation | 4-5 s (R5 run) | 1-2 s with the limit at 1,200 | C6 (owner) |
| Update | 7.7 s | about 5.8 s | C3 |
| Reference pass | 1.4-1.7 s | 0 s if the KL weight is 0 | C7 (owner) |
| Total step | 23.4 s | 16-18 s; 13-14 s with C6 and C7 | Phase C |
| GPU idle | about 40 % | below 20 % (target) | Phase C |
| Pixels moved per step | 1.61 GB fp32, built twice | 32 small uint8 images | C2 |
| Small host waits per step | about 1,250 | a few tens (target) | C4 |

**Worked example 3: hours for one pass over the shifted data.** The shifted data has 24,624 rows (a row is one image
with its doctor masks; some images repeat). One step uses 32 rows, so one pass is 770 steps. Previous: 770 x 23.4 s
= 5.0 hours. New: 770 x 17 s = 3.6 hours. The R5 host cost $0.59 for each hour. Thus the cost falls from about
$2.95 to about $2.15 for each pass. (The full RL run cost $4.56 for 770 steps, because its steps took 27.5 s and the
host price was higher.)

---

## 6. The checker and the engine

| Part | Previous | New |
|---|---|---|
| Polygonizer used for targets | CPU numpy (breaks the CUDA rule) | CUDA |
| CUDA polygonizer work | All 4(W+1)(H+1) edge places, about 40 jump rounds, 3 host waits | Only true edges, ring start from labels, 1 host wait |
| Labels in the scorer | Always 4 passes, one thread for each row | Skip when nothing is small; runs in bit words |
| Validity totals | One thread for each shape | One warp for each shape |
| Memory growth | `cudaMalloc` and `cudaFree` stop the device | Pool and batched upload |
| Evaluation metrics | CPU pool, masks come back to the CPU | GPU, only numbers come back |
| GPU tests | Only on rented boxes | A GPU test job for every change |
| 3D | CPU volume code, skipped slices split a nodule, no slice thickness | GPU 3D verifier, skipped slices kept, slice thickness stored |
| Correctness | Exact, bit-identical, 0 sanitizer errors | Same, kept by the same oracles and golden files |

**Worked example 4: polygonizer work on a 256 x 256 mask.** Previous: 4 x 257 x 257 = 264,196 edge places, each
visited in about 40 rounds, about 10.6 million visits. New: about 1,000 true edges, about 20 rounds, about 20,000
visits. The ratio is about 500.

---

## 7. Code structure

| Part | Previous | New |
|---|---|---|
| Trainer builder | In an experiment script, with a loose argument object | In the library, with typed inputs |
| Evaluation | Builds a full trainer and needs the training set | Uses a sampler. No trainer. |
| Changes to TRL internals | 9 private names in several places, no check | One adapter file, a start-up check, TRL pinned |
| Reward weights | Copied by hand in 3 places | One source |
| Learning rate | Two sources (5e-5 and 1e-5) | One source |
| Run settings | 16 YAML files with settings inside command strings | One typed run specification with inheritance |
| Tests that depend on experiment scripts | 35 of 71 test files | 0 |
| Layers | Not enforced | Enforced by an import rule in `make lint` |
| Documentation | README, plans, notes | Plus a full manual: getting started, architecture, decision records, reference, guides, runbooks, testing, glossary |

---

## 8. RL infrastructure

| Part | Previous | New |
|---|---|---|
| Framework | TRL 1.17.1 with our fast trainer | Same on 1 GPU; verl async at 2 GPUs; verl or NeMo-RL at 8+ |
| Generation and training | In turns | In turns on 1 GPU; at the same time on 2+ GPUs |
| Requests to vLLM | 256 single requests | 32 requests of 8 answers |
| Reward call | Inside the trainer process | Same by default; HTTP service with versions when other programs need it |
| Update method | Dr. GRPO, KL 0.001, clip-higher (never active) | Dr. GRPO; CISPO or GSPO with async; KL 0 (owner) |
| Run record | Logs and model | Plus reward version, code version and data version in each result |

---

## 9. Research direction

| Area | Previous | New |
|---|---|---|
| Main claim | Exact GPU checker as an RL reward for outlines | Certified Geometry: one exact verifier as reward, evaluation and data filter, from 2D outlines to 3D voxel volumes and robot assets |
| Novelty check | Not done against 2026 papers | Done. SimpleSeg is the nearest work. Four open claims (Direction section 1.3). |
| Model choice | Smallest that passes | Same, plus Qwen3-VL-4B as a scaling point; Qwen3.5 rejected for now |
| Test-time method | Mean of 8 answers | Legal filter, then the medoid of 8 |
| How the model writes a shape | JSON outline on a 0-1000 grid | Same (costs no IoU); test a centre plus 16 radii for 60 % fewer tokens |
| 3D | Volume Betti numbers on the CPU | GPU 3D verifier and a 3D set benchmark |
| 3D voxels | Volume code with two bugs (empty filled slices, no slice thickness) | Exact GPU 3D verifier (b0, b1, b2), then a whole-nodule model trained with a 3D reward |
| Physical AI | None | Certify generated robot assets: watertight, manifold, collision shape keeps holes |

---

## 10. One step, previous against new: who talks to whom

```
PREVIOUS                                         NEW
CPU        vLLM       checker    trainer         CPU        vLLM       checker    trainer
 |           |           |          |             |           |           |          |
 |--pixels x2, prompts x2 (CPU)---->|             |--32 images (uint8)------------------>|
 |--256 single requests->|          |             |--32 requests x 8------>|           |
 |           |--256 answers-------->|             |           |--256 answers-------->|
 |           |           |--rewards>|             |           |           |--rewards>|
 |<--1,250 small waits--------------|             |<--a few tens of waits------------|
 |           |<--700 weight calls---|             |           |<--LoRA layers only---|
 |           |  (wake, sleep: 1.7 s)|             |           |  (awake, ~0.05 s)    |
```

## 11. What stays the same

Some parts were right from the start. We keep them.

- The exact integer geometry and the one-source CPU and GPU engine.
- The set-valued reward idea: an outline is correct if any doctor drew that shape.
- The rule that doctor masks never leave the GPU during scoring.
- The fast trainer's shared-prompt forward pass, image-feature cache and position cache.
- The preflight check, the North America rule and the money limit in the code.
- The owner's decision for each run.

---

## Summary

The previous system proved the idea: an exact CUDA checker can train a small model to draw medical outlines. The
new system keeps that core, fixes the three measured causes of the IoU limit, removes about half of the idle GPU
time, and moves the last CPU geometry to CUDA. It then grows the same verifier into 3D voxel volumes and robot assets.
