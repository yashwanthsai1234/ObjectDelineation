# Direction: what we do next

Date: 2026-11-08. Language: ASD-STE100 Simplified Technical English.

This report tells what to do next, in what order, and why. The facts come from `Research2.0.md`. All gains are
estimates until a GPU run measures them.

---

## 0. Words in this report

We use the same words as `Research2.0.md`: model, image, outline, doctor, checker, reward, answer, group, step,
miss, offset, ceiling, GPU. We add these words:

| Word | Meaning |
|---|---|
| track | One line of work with one goal. |
| gate | A test with a number. Work goes to the next phase only when the gate passes. |
| run | One rented GPU job. |
| verifier | The checker when other programs use it as a service: as a reward, as an evaluation, or as a data filter. |
| nodule | One lung nodule, with all its images. |
| row | One training item: one image with its doctor masks. Some images appear in more than one row. |
| medoid | Of 8 answers, the answer with the highest mean IoU against the other 7. It needs no doctor masks. |
| pixel vote | A mask where each pixel is on if most of the 8 answers cover it. |
| GED, HM-IoU | Measures that compare a set of answers with a set of doctors. Papers on many doctors use them. |
| NVTX | Time marks in the code that the NVIDIA profiler shows. |
| reward version 5 | The new reward in section 5. |

---

## 1. The goal in one sentence

**Build one exact GPU verifier for geometry, and use it to train, test and filter models that draw shapes: first 2D
outlines, then 3D voxel volumes, then assets for robot simulation.**

The owner decided on 3D voxels as the main direction (2026-11-08). CAD was a study only. It is not in the plan
(see "Parked: CAD" in section 3).

We call this goal **Certified Geometry**.

### 1.1 Why this goal

- The verifier already exists for 2D. It is exact, bit-identical on CPU and GPU, and fast (7.74 million outlines
  each second in a speed test).
- No published RL work uses an exact GPU verifier with topology.
- No published RL work uses a set of doctors as the reward. Real shapes often have more than one correct answer.
- The same verifier can grow from 2D rings to 3D voxels to closed meshes. Each step reuses the same maths.

### 1.2 What we must not claim

- We must not say "first RL for models that write outlines". SimpleSeg (January 2026) did this on natural images.
- We must not say "exact verifier for world-model video". Video tests lift geometry with learned models first.
- We must not compare our IoU with published LIDC numbers. The test methods differ.

### 1.3 What we can claim when the work is done

1. The first RL reward that scores an answer against a set of doctors, with exact topology.
2. The first medical model that writes outlines directly as text, without a mask tool.
3. An exact GPU verifier as a reward service for 2D outlines and 3D voxel volumes.
4. A test of the question: "Does an exact reward teach a model better than an approximate reward?"

---

## 2. Rules that do not change

These are the owner's rules. Every phase obeys them.

| Rule | How we obey it |
|---|---|
| The checker, the scorer and the polygonizer use CUDA. | Phase D moves the last CPU polygonizer calls to CUDA. |
| The GPU does the work when it can. | Phase C removes CPU work from the step. Phase D moves evaluation to the GPU. |
| No money before a check. | Each run starts after the preflight check passes and the owner approves the cost. |
| Rent hosts in the USA or Canada only. | The rule stays in the code. |
| Use the smallest model that passes. | We keep the 2B model as the main result. The 4B model is a scaling point. |
| The owner decides each run. | Each phase ends with a gate and an owner decision. |
| Simple code. | Each change removes the code it replaces. |

---

## 3. The plan in one picture

```
  Phase A          Phase B           Phase C          Phase D
  Measure     -->  Quality      -->  Speed       -->  CUDA core
  (Mac, $0.25)     (GPU, ~$4.5)      (with B runs)    (Mac + 1 short GPU test)
      |                |                 |                 |
      v                v                 v                 v
  honest numbers   IoU 0.70-0.73     16-18 s/step     polygonizer on GPU
                                                           |
      +----------------------------------------------------+
      |
      v
  Phase E          Phase F           Phase G          Phase H          Phase I
  Structure   -->  Scale       -->   3D voxel    -->  3D voxel    -->  Physical AI
  (Mac)            (4B, baselines,   verifier         model + RL       (asset
                    2 GPUs)           (GPU)            (whole nodule)   certifier)
```

Phases A to E fix and finish the 2D system. Phases F to I extend it. Phase E can run beside B, C and D because it is
Mac work.

---

## 4. The phases

### Phase A: measure first (Mac, one short run, about $0.25)

**Goal.** Make every number fair before we spend money on training.

| # | Task | Why | Cost |
|---|---|---|---|
| A1 | Resample whole nodules, not single images, in the final test. | The present confidence intervals are too narrow. | Mac |
| A2 | Make the ceiling fair: compare the model with k - 1 doctors, as each doctor is. | The present ceiling favours the model. | Mac |
| A3 | Save every answer's text in the final test. Run the test again on the RL model. | We need single answers for best-of-8 and selection. | about $0.25 |
| A4 | Select an answer with no doctor masks: keep legal answers, then take the medoid. | Legality then goes to about 100 %. IoU estimate: +0.01 to +0.04. | Mac, after A3 |
| A5 | Report GED, Dice against the majority, HM-IoU and boundary IoU beside our IoU. | Other papers use these values. | Mac, after A3 |
| A6 | Add time marks (NVTX) to each phase of the step. | We must measure before we change speed. | Mac |
| A7 | Add a start-up test that checks the 9 private TRL names that our trainer replaces. | A TRL update can stop our trainer without a message. | Mac |

**Gate A.** If the medoid or the pixel vote raises validation IoU by 0.015 or more, with a confidence interval above
0, it becomes the default test method.

### Phase B: quality (3 runs, $4.4-5.0)

**Goal.** Find the far nodules, and use the full learning signal of each step.

| # | Run | What changes | Estimate | Cost |
|---|---|---|---|---|
| B1 | Shifted images | RL uses centred, shifted and rescaled training images (24,624 rows). The code exists. Log the share of answers with ratio 0 and the share of groups with no spread. | Test IoU 0.645 to 0.69-0.71. Far misses from 23 % to below 10 %. | about $1.6 |
| B2 | Number match | Fix the trainer-vLLM mismatch. Try two fixes against the present method, 150 steps each: (a) cut each token's ratio at 2, in place of a drop of the whole answer; (b) use FP16 numbers in both programs. | Correction ratio near 1. Signal from about 50 % to about 100 % of answers. +0.01 to +0.02 IoU. | about $1.2-1.8 |
| B3 | Reward version 5 | New reward (section 5). Show hard images more often. Use 16 answers for hard images. | Misses down by 30-50 % more than B1. +0.02 to +0.04 IoU. | about $1.6 |

**Gate B1.** IoU on the whole validation split is 0.66 or more (today 0.618). Validation misses at offset above 0.3
fall below 20 % (today 31.8 %; the same value on the test split is 23 %). Centred images stay at 0.73 or more.

**Gate B2.** The mean ratio is 0.97 or more. Fewer than 0.5 % of answers are masked. Validation IoU is not lower.

**Gate B3.** Misses fall by 25 % or more against B1. IoU on found nodules does not fall by more than 0.01. Before the
run, test that a very large outline does not get a high location term.

### Phase C: speed (inside the Phase B runs)

**Goal.** Cut the idle GPU time from about 40 % of each step.

| # | Change | Saving per step (estimate) | Risk |
|---|---|---|---|
| C1 | Push only the changed adapter weights into vLLM, in place. Keep vLLM awake. Skip the vision weights. | about 1.5 s | Low. The result is bit-identical. |
| C2 | Calculate the prompt tokens and positions once. All LIDC images are 512 x 512. Build pixels on the GPU from the 32 images. | 0.6-1.5 s | Low |
| C3 | Use micro-batches of 32 x 8 in place of 16 x 16. The update is the same. | 1.3-2.5 s | Memory test needed |
| C4 | Remove small host waits: keep loss values on the GPU, use index lists in place of masks. | 0.3-0.8 s | Low |
| C5 | Send the 8 answers of one image as one request to vLLM. | 0.2-1 s | Low |
| C6 | Lower the answer limit from 1,536 to about 1,200 tokens. A repo note says that the longest true LIDC outline needs 1,182; we did not check it again. | 0.3-4 s | Owner decision |
| C7 | Set the KL weight to 0. This removes the reference pass. | 1.4-1.7 s (the reference pass in R5) | Owner decision |
| C8 | Keep some layers without gradient checkpointing. | 1-2.5 s | Competes with C1 for memory |

**Gate C.** NVTX timings show 18 s or less for each step with C1 to C5, and the same rewards on the first step.

**Worked example 1: the new step time.** Today one step takes 23.4 s. C1 saves 1.5 s. C2 saves 1.0 s (middle of
0.6-1.5). C3 saves 1.9 s (middle of 1.3-2.5). C4 saves 0.5 s. C5 saves 0.6 s. The sum is 5.5 s. The new step is
23.4 - 5.5 = 17.9 s. C7 then saves 1.5 s. C6 cuts the long tail from 4-5 s to 1-2 s, a saving of about 3 s. The
estimate is 17.9 - 1.5 - 3 = 13.4 s, thus about 13-14 s.

### Phase D: CUDA core (Mac, one short GPU test)

**Goal.** Obey the CUDA rule everywhere, and make the checker work-optimal.

| # | Task | Why |
|---|---|---|
| D1 | Connect the target builder and the certificate code to the CUDA polygonizer. Move the numpy polygonizer to the test oracles. | The rule says the polygonizer uses CUDA. Today production uses the CPU. |
| D2 | Make the CUDA polygonizer visit only true edges. Find each ring start in closed form from the component labels. One GPU block for each mask. One host wait for each call. | The present version does about 500 times the necessary work on a 256 x 256 mask. |
| D3 | Skip labelling passes when the size filter changes nothing. Label runs of pixels in bit words. | The scorer does 4 passes where 2 are usually enough. |
| D4 | Use one warp for each shape in the validity totals. Use a grid to find candidate edge pairs. | One bad shape can stop a batch for about 0.1 s. |
| D5 | Move evaluation metrics to the GPU: exact integer distance transform, GED, HM-IoU. | Evaluation brings masks back to the CPU today. |
| D6 | Add a GPU test job that builds the CUDA path and runs its tests. | The GPU raster is a separate program from the CPU raster. Only rented boxes test it now. |

**Gate D.** The target hashes stay byte-identical. All checker tests pass on the GPU. The sanitizer finds 0 errors.

### Phase E: structure (Mac, beside B, C and D)

**Goal.** Make the code a clean product with clear parts. `Architecture.md` gives the design.

| # | Task |
|---|---|
| E1 | Fix small debts first: one learning-rate source, no unused worker pool, one source for reward weights. |
| E2 | Move the trainer builder from an experiment script into the library. |
| E3 | Make one verifier package. It has no torch and no TRL. |
| E4 | Make one sampler interface with a vLLM version and an HF version. Evaluation then needs no trainer. |
| E5 | Split the fast trainer file into pure functions, model caches and TRL glue. |
| E6 | Replace the 16 command-line YAML files with one typed run specification. |
| E7 | Offer the verifier as a batch service with a version number. |

**Gate E.** `make test` and `make lint` pass after each step. No result changes.

### Phase F: scale (more money, owner decision)

**Goal.** Show that the method scales with model size, beats strong baselines, and works on a second data set.

| # | Task | Cost (estimate) |
|---|---|---|
| F1 | Train Qwen3-VL-4B with the same code. It is the scaling point: 2B against 4B. | $8-14 |
| F2 | Run baselines on our splits: nnU-Net 2D, MedSAM2 with a box, SAM 3, SimpleSeg. | $2-5 |
| F3 | Run the same recipe on RIGA (optic disc and cup). Compare with MRNet (cup Dice 87.2). | $2-4 |
| F4 | Test a shorter shape code: a centre and 16 radii. It is always legal and uses about 60 % fewer tokens. Use it only for shapes without holes. | $2-3 |
| F5 | On 2 GPUs: test verl in async mode against our trainer. Pass if the reward curve is the same and the time for each batch is 60 % or less. | a few GPU-hours |

### Phase G: 3D

**Goal.** Score whole nodules in 3D on the GPU, with the same exactness as in 2D.

| # | Task |
|---|---|
| G1 | Fix the volume code: do not fill a skipped slice with empty pixels. Store the slice thickness. |
| G2 | Build the 3D checker in CUDA: run labels, exact Euler number, boundary faces. Check it against the CPU volume code on all LIDC volumes. |
| G3 | Make a 3D set benchmark: doctor volumes, best-match IoU, topology against any doctor, presence F1. |
| G4 | Add a whole-nodule reward. Measure false tunnels and false cavities before and after RL. |

### Phase H: 3D voxel model and RL

**Goal.** Train the model to draw a whole nodule in 3D, and reward it with the exact 3D voxel verifier.

| # | Task |
|---|---|
| H1 | The model writes a whole nodule as one outline on each slice. The 3D verifier fills the outlines into one voxel volume and checks it: one part (b0), no false tunnels (b1), no false cavities (b2). |
| H2 | RL with the whole-volume reward: volume IoU against the set of doctor volumes, plus 3D topology. Compare it with the sum of slice rewards. Measure false tunnels and false cavities. |
| H3 | Test a shorter 3D code (for example, a centre and radii on each slice) for fewer tokens. Keep the exact verifier as the judge. |
| H4 | Put the 3D evaluation measures on the GPU: volume IoU against a set, surface distance, presence F1. |

Cost: estimate after Phase G, from the measured speed of the 3D verifier.

### Parked: CAD

CAD was a study only. Section 8.4 of `Research2.0.md` keeps the findings: CAD RL papers score with approximate CPU
geometry, so an exact verifier could help there. The owner chose 3D voxels. CAD is not in the plan.

### Phase I: Physical AI

**Goal.** Use the 3D voxel verifier to certify assets for robot simulation.

| # | Task |
|---|---|
| I1 | Run the 3D verifier on public generated assets. Count assets that are not watertight, not manifold, self-intersecting, or with a collision shape that fills a hole. |
| I2 | Show one robot task where certified assets give a better result than uncertified assets. |

---

## 5. The new reward (version 5)

Agent H2 derived this reward. It matches the metric that we report.

```
reward = legal x ( IoU + 0.1 x location ) + topology

location = only for misses (IoU < 0.1):
           exp( - (distance between centres)^2 / (2 x radius^2) )

advantage = reward - mean reward of the other 7 answers in the group
            (no division by the group spread)

boundary F and area: limits, not goals
            (small fixed weights, or a gate at 0.9)
```

**Worked example 2: a group of 8 misses.** Today all 8 answers get IoU near 0.02, so all 8 get almost the same
reward. The advantage is about 0 for each answer. Nothing is learned. With version 5 (radius 7 px):

| Answer | Distance to the centre | Location term |
|---|---|---|
| 1 | 5 px | exp(-25/98) = 0.78 |
| 2 | 10 px | exp(-100/98) = 0.36 |
| 3 to 8 | 20 px | exp(-400/98) = 0.02 |

The mean location term of the other 7 answers is (0.36 + 6 x 0.02) / 7 = 0.07. Answer 1 gets an advantage of
0.1 x (0.78 - 0.07) = +0.07. The model learns to move toward the nodule.

**Worked example 3: why no division by the spread.** Group 1 has IoU 0.61, 0.61, 0.62, 0.62, 0.62, 0.62, 0.63,
0.63. Its spread (standard deviation) is 0.007. Group 2 has IoU 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9. Its spread is
0.229. With division by the spread, both groups give advantages of the same size. Group 1 teaches noise as strongly
as group 2 teaches a real difference. Without division, the advantages of group 2 are 0.229 / 0.007 = about 32 times
larger.

---

## 6. Decisions for the owner

| # | Decision | Options | Our recommendation |
|---|---|---|---|
| 1 | Budget for Phase A and Phase B | $0.25, then $4.4-5.0 | Approve Phase A now. Approve B1 after Gate A. |
| 2 | KL weight | Keep 0.001, or 0 | Try 0 in B2. It saves 1.4-1.7 s. The reward is exact, so noise cannot trick it. |
| 3 | Answer limit | 1,536 or about 1,200 tokens | About 1,200. A repo note says that the longest true LIDC outline needs 1,182; we did not check it again. |
| 4 | Second GPU | Stay on 1, or rent 2 | Stay on 1 until Phase C is done. Then test 2 in F5. |
| 5 | Model size | 2B, or 2B plus 4B | Keep 2B as the main model. Add 4B in F1. |
| 6 | NVIDIA skills | Install, or not | Install only when a phase needs them: the NeMo-RL skills at F5, the NVIDIA CT segmentation skill at F2. |
| 7 | First direction after 2D | 3D voxels | Decided by the owner: 3D voxels (Phases G and H). CAD is parked. |

---

## 7. How the work flows each time

```
Owner        Main / agent      Mac tests        Preflight        GPU host         Checker
  |               |                |                |                |                |
  |--approve----->|                |                |                |                |
  |               |--change code-->|                |                |                |
  |               |<--tests pass---|                |                |                |
  |               |--run preflight---------------->|                |                |
  |               |<--------------cost + data OK----|                |                |
  |<--cost--------|                |                |                |                |
  |--go---------->|                |                |                |                |
  |               |--launch run------------------------------------>|                |
  |               |                |                |                |--score-------->|
  |               |                |                |                |<--rewards------|
  |               |<--logs, model, test results-----------------------|                |
  |               |--check gate--->|                |                |                |
  |<--result------|                |                |                |                |
```

---

## 8. Expected results (estimates)

The Phase B target is near the present ceiling. The present ceiling favours the model (section 4.8 of
`Research2.0.md`). Phase A gives the fair ceiling. Then the owner sets the final target below the fair ceiling.

| Measurement | Today | After Phase A | After Phase B | After Phase C |
|---|---|---|---|---|
| Test IoU | 0.645 | 0.655-0.685 | 0.70-0.73 | same as B |
| Ceiling (present method) | 0.732 | fair value from A2 | - | - |
| Boundary F | 0.816 | about 0.83 | 0.86-0.89 | same as B |
| Legal answers | 99.4 % | about 100 % | about 100 % | about 100 % |
| Far misses (offset above 0.3) | 23 % | 23 % | below 10 % | same as B |
| Time for each step | 23.4 s | 23.4 s | 23.4 s | 16-18 s (13-14 s with C6, C7) |
| Money | - | about $0.25 | $4.4-5.0 | $0 extra |

**Worked example 4: the money for Phases A to C.** Phase A costs $0.25. Phase B costs $1.6 + ($1.2 to $1.8) + $1.6 =
$4.4 to $5.0. Phase C
uses the Phase B runs, so it costs no extra money. The total is about $5. The spend ledger shows $16.90 spent and $0.60 left under the
$17.50 limit. Thus Phase A fits, but the owner must raise the limit before Phase B.

---

## 9. Papers and public results

| # | Result | Where | Needs |
|---|---|---|---|
| 1 | Set-valued, topology-exact RL for a model that writes outlines, on LIDC and RIGA | MICCAI or MIDL | Phases A, B, F2, F3 |
| 2 | Exact 3D topology reward for whole-nodule voxel volumes (false tunnels and cavities) | Workshop or MICCAI | Phases G, H |
| 3 | Does an exact reward teach better than an approximate reward? 2D and 3D voxels | Main conference | Phases G, H |
| 4 | An open verifier library and benchmark for geometry rewards | Open source | Phases D, E7 |

---

## Summary

First measure fairly, then fix the three known causes: centred training images, the number mismatch, and a reward
with no signal for misses. Then cut the idle GPU time and move the last CPU geometry to CUDA. After that, grow the
same exact verifier from 2D outlines to 3D voxel volumes, then to robot assets.
