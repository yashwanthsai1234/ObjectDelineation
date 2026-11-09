# Research 2.0: what we studied, how, and what we found

Date: 2026-11-08. Language: ASD-STE100 Simplified Technical English.

This report tells the story of the Research 2.0 study. It says what we asked, how we found the answers, and what
the answers are. The next reports use these facts:

- `Direction.md` tells what to do next.
- `Architecture.md` tells how to build the new system.
- `previous vs new.md` compares the old system with the new system.

---

## 0. Words in this report

We use one word for each thing. We do not change the word.

| Word | Meaning |
|---|---|
| model | The vision-language model that we train. Today it is Qwen3-VL-2B with a LoRA adapter. |
| image | The square part of a CT slice or of an eye (fundus) image that the model sees. |
| outline | The polygon that the model writes as text. It has rings and holes. |
| doctor | A person who drew a true outline. LIDC has up to 4 doctors for each nodule. |
| checker | The CUDA program that checks an outline and gives it a score. |
| reward | The number that the checker gives to one answer during RL. |
| answer | One outline that the model writes for one image. |
| group | The 8 answers that the model writes for one image in one step. |
| step | One RL update. One step uses 32 images and 256 answers. |
| miss | An answer with an IoU below 0.1. The model did not find the nodule. |
| offset | The distance from the image centre to the nodule centre, as a part of the image side. |
| ceiling | The mean IoU that a doctor gets against the best of the other doctors. It is the best score we can expect. |
| nodule | One lung nodule. It has many images (slices) and up to 4 doctor outlines on each image. |
| IoU | Overlap of two shapes: shared area divided by total area. 1 is a perfect match. |
| Dice, GED | Two other overlap measures that papers use. GED compares a set of answers with a set of doctors. |
| b0, b1, b2 | Betti numbers: the count of parts, holes (or tunnels in 3D) and closed cavities. |
| correction ratio | The probability of an answer in the trainer, divided by its probability in vLLM. The ideal value is 1. |
| entropy | A measure of how much the answers of the model vary. Low entropy means almost the same answer each time. |
| LoRA | A small set of extra weights that RL changes. The main model weights stay fixed. |
| KL weight | The strength of a term that keeps the model near its start point. |
| GPU | The graphics card. Today it is one rented RTX 5090 with 32 GB. |

---

## 1. The question

The owner asked five questions:

1. Why are IoU, legality and boundary F not better after RL?
2. Where does the system lose GPU time? Where does data move between CPU and GPU without need?
3. Is our RL system the best? What do the best teams and papers use?
4. What is the best mathematics and the best algorithm for each part?
5. Where can the project go next: 3D voxels, CAD, and Physical AI?

The owner also gave two fixed rules:

- The checker, the scorer and the polygonizer must use CUDA.
- The GPU must do the work when it can. The CPU does only the work that the GPU cannot do.

---

## 2. Where we started

The model writes an outline. The checker scores the outline. RL makes the model better.

```
+-----------+     +-----------+     +-----------+     +-----------+
|  image    | --> |  model    | --> |  outline  | --> |  checker  |
| (CT image)|     | Qwen3-VL  |     |  (text)   |     |  (CUDA)   |
+-----------+     +-----------+     +-----------+     +-----+-----+
                        ^                                   |
                        |          reward (one number)      |
                        +-----------------------------------+
```

These are the results before this study. We read each number from the run files.

| Measurement | Value |
|---|---|
| Test IoU (2,200 images, 8 answers each) | 0.645 (median 0.756) |
| Test boundary F | 0.816 |
| Test Dice against the doctors' majority | 0.701 |
| Test legal answers | 99.4 % |
| Ceiling on the test split | about 0.732 |
| Time for one step | 23.4 s |
| Money spent on GPUs so far | $16.90 of a $17.50 limit ($0.60 left) |

---

## 3. How we did the study

### 3.1 The team

We used 13 agents. Each agent had one duty. No agent read the history of older sessions. Each agent used the code,
the run files, the papers and the web. Each agent wrote one findings file.

| # | Agent | Duty | Main method |
|---|---|---|---|
| A | History | What we planned, what we did, and why it changed | Read plans, results and 87 commits |
| B | GPU audit | Find each CPU-GPU transfer and each idle second | Read the code line by line, and the step logs |
| C | RL quality | Find why the IoU stops at 0.645 | Read the logs; measure each image again |
| D | RL systems | Compare 17 RL frameworks and 2025-2026 methods | Web, GitHub, papers |
| E | Vision models | Find the best methods and models for outlines | Web, papers, Hugging Face leaderboards |
| F | 3D, CAD, Physical AI | Find where an exact checker can go next | Web, papers, NVIDIA catalog |
| G | GPU methods | Find the best GPU tools for this loop | Web, NVIDIA skills, installed code |
| H1 | Geometry maths | Find the best algorithm for each geometry part | Proofs and papers |
| H2 | RL maths | Put the reward and the update on firm maths | Proofs and papers |
| H3 | Shape coding | Find the best way to write a shape as text | Measurements on 3,444 doctor masks |
| H4 | Code quality | Find weak structure and slow algorithms | Read all code |
| Docs | Documentation | Write the repo manual | Read all code |
| Main | Coordinator | Give duties, check facts, write these reports | All of the above |
| R1-R3 | Paper readers | Read the 2025-2026 papers that the agents cited, plus newer 2026 papers | Full text and tables; section 12 |

After the main study, three reading agents read the papers again, one agent for each topic. One reader split its work
into four smaller batches. Section 12 lists every paper that they read.

### 3.2 The order of work

```
  owner question
        |
        v
+-------------------+
| Main: read repo   |
| give 13 duties    |
+---------+---------+
          |
          +--------------------+--------------------+--------------------+
          v                    v                    v                    v
+-----------------+  +-----------------+  +-----------------+  +-----------------+
| Repo audits     |  | Web research    |  | Maths and       |  | Documentation   |
| A, B, C, H4     |  | D, E, F, G      |  | algorithms      |  | Docs agent      |
|                 |  |                 |  | H1, H2, H3      |  |                 |
+--------+--------+  +--------+--------+  +--------+--------+  +--------+--------+
         |                    |                    |                    |
         +--------------------+---------+----------+--------------------+
                                        v
                             +---------------------+
                             | Main: check facts   |
                             | against code + logs |
                             +----------+----------+
                                        v
                             +---------------------+
                             | 4 reports + review  |
                             +---------------------+
```

### 3.3 Who talked to whom

```
Owner        Main           Agents A-H4         Web / papers       Repo + run logs
  |            |                 |                    |                   |
  |--question->|                 |                    |                   |
  |            |--13 duties----->|                    |                   |
  |--rules---->|--rules--------->|                    |                   |
  |            |                 |--search----------->|                   |
  |            |                 |<-----papers, docs--|                   |
  |            |                 |--read code, logs-------------------->|
  |            |                 |<--------------------------numbers----|
  |            |<--findings------|                    |                   |
  |            |--check numbers---------------------------------------->|
  |            |<-------------------------------------------same values--|
  |<-reports---|                 |                    |                   |
```

### 3.4 How we kept the facts true

- Each agent marked each number as measured, read from code, or estimated.
- The main agent checked the key numbers again in the run files. Section 9 lists them.
- When two agents gave different numbers, we used the measured value. Section 9.2 shows one correction.
- An estimate is always called an estimate in these reports.

### 3.5 The skills we used

The main session and every agent loaded these skills: NVIDIA skill finder, engineering code review, debug, system
design, architecture, documentation, and Hugging Face best model. Agents also used Hugging Face papers, TRL
training, testing strategy, tech debt and codebase design where they fit.

Two skills gave live data:

- **Hugging Face best model.** It read the official leaderboards and the model cards. It found the new Qwen3.5
  family (0.8B, 2B, 4B, 9B, February 2026).
- **NVIDIA skill finder.** It read the live NVIDIA catalog. It found skills for CT segmentation, DICOM volumes,
  CAD to SimReady, Cosmos and NeMo-RL. We did not install them. The owner must approve each installation.

---

## 4. Finding 1: why the IoU stops at 0.645

### 4.1 Most of the gap is misses, not bad outlines

The model has two failure types. It can find the nodule and draw a poor outline. Or it can miss the nodule. We
measured both on the validation split (1,084 images).

| Group of images (validation split) | Share of images | IoU | Ceiling |
|---|---|---|---|
| All images | 100 % | 0.618 | 0.734 |
| Images where the model finds the nodule | 87.4 % | 0.706 | 0.734 |
| Misses | 12.6 % | below 0.1 | - |

The misses are 12.6 % of the images, but they cause 76 % of the gap to the ceiling. When the model finds the
nodule, its boundary F is 0.896.

**Worked example 1: what one miss costs.** Take 100 images. 87 images score 0.706. 13 images score about 0.05.
The mean is (87 x 0.706 + 13 x 0.05) / 100 = 0.621. If the model finds all 100 nodules at 0.706, the mean is 0.706.
The misses alone take away 0.085 IoU.

### 4.2 The misses are far from the centre (test split)

| Offset | Images (test) | Model IoU | Ceiling | Misses |
|---|---|---|---|---|
| below 0.15 | 382 | 0.76-0.77 | 0.71-0.74 | 0-0.5 % |
| 0.15 to 0.2 | 280 | 0.710 | 0.729 | 1.8 % |
| 0.2 to 0.3 | 713 | 0.670 | 0.742 | 5.9 % |
| above 0.3 | 825 | 0.545 | 0.724 | 23.0 % |

Near the centre, the model is as good as the doctors. Far from the centre, the model misses one nodule in four.

### 4.3 RL caused the problem

RL used only centred images. The tests use shifted images.

| Data | Median offset | Share above 0.2 |
|---|---|---|
| RL training images | 0.012 | 0.0 % |
| Validation images | 0.230 | 60.5 % |

Before RL, the model scored about the same at every offset (0.49 to 0.54). After RL, the model gained 0.24 IoU on
centred images, but only 0.05 on far images. Thus RL taught the model a habit: "the nodule is in the centre".

The fix exists in the code but has not run on a GPU. It adds shifted and rescaled training images, each with its own
doctor masks. It needs about $1.6.

### 4.4 Half of each step's learning signal is lost

This is a new finding. vLLM writes the answers. The trainer then calculates the probability of each answer again.
The two probabilities must be equal, because the weights do not change between them. Then the correction ratio is
1. The logs show a different result.

| Steps | Mean ratio | Lowest ratio | Highest ratio |
|---|---|---|---|
| 1-10 | 0.994 | 0.24 | 2.34 |
| 100-199 | 0.704 | 0.0 | 2.89 |
| 700-770 | 0.847 | 0.0 | 2.84 |

The lowest ratio is 0 in 98.1 % of steps. An answer with ratio 0 gives no learning signal. The cause is a small
difference in the numbers of the two programs. Both use 16-bit numbers (bf16), and vLLM adds the LoRA weights into
the main weights, but the trainer keeps them apart. The estimate is that only 40 % to 55 % of the 256 answers in a
step teach the model.

**Worked example 2: the lost signal.** A step has 256 answers. The trainer sees 1.0 as the ideal ratio. Suppose that
120 answers get a ratio near 0 or a cut ratio. Then only 136 answers carry full signal. That is 53 %. The step costs
23.4 s, but it does the work of about 12.4 s.

### 4.5 The reward stops giving signal

On centred training images, the reward is almost at its maximum. The table shows the logs.

| Steps | Training IoU | Boundary F | Entropy | Reward spread |
|---|---|---|---|---|
| 1-10 | 0.523 | 0.717 | 0.36 | 1.44 |
| 100-199 | 0.738 | 0.920 | 0.085 | 0.52 |
| 700-770 | 0.771 | 0.943 | 0.067 | 0.34 |

When all 8 answers in a group miss, all 8 get the same low reward. The group then gives no signal. The topology part
of the reward gives points only when IoU is 0.5 or more. Thus a group of misses learns nothing.

### 4.6 The output format is not the limit

We measured this on 3,444 doctor masks.

- Every true outline turns back into its mask exactly. The 0-1000 grid is 10 to 20 times finer than one pixel.
- Only 0.12 % of answers reach the length limit.
- Only 2.5 % of the model's points fall between pixel corners. Moving them to the corners changes nothing.

### 4.7 Small nodules are hard for everyone

| Nodule size | Mean IoU between two doctors | IoU when the true mask moves 1 pixel |
|---|---|---|
| below 5 px | 0.363 | 0.455 |
| 5 to 9 px | 0.644 | 0.665 |
| all sizes | 0.683 | - |

The value 0.683 is a mean over pairs of doctors. The ceiling (0.732 on test, 0.734 on validation) is higher, because
each doctor gets the best of the other doctors.

For a round nodule with diameter d, a shift of s pixels gives IoU of about 1 - 8s/(pi x d).

**Worked example 3: an 8-pixel nodule.** For IoU 0.8, we need 1 - 8s/(pi x 8) = 0.8. Then s = 0.2 x pi = 0.63 px.
With the shape error included, the outline must be within about 0.4 px. One visual token of the model covers about
3.1 px. Thus the model must place the outline to one eighth of a token.

### 4.8 Two metric problems

- **The ceiling favours the model.** The model gets the best IoU of all k doctors. Each doctor gets the best of only
  k - 1 other doctors. The true gap is larger than 0.732 - 0.645.
- **The confidence intervals are too narrow.** The final test resamples single images. All other tests in the repo
  resample whole nodules, with all their images together.

---

## 5. Finding 2: where the GPU waits

### 5.1 One step today

We measured one step in the R5 run (steps 3-15). The timed step is 23.43 s. Work outside the timed step adds 0.22 s.

```
time (s)  0          5          10         15         20      23.4
          |----------|----------|----------|----------|---------|
generate  [=====================]                                  11.1 s
            sync 1.3 | decode 9 incl. long tail | sleep 0.3
batch                            [=========]                        4.7 s
                                  CPU 3.0-3.2 | reference 1.4-1.7
update                                      [==============]        7.7 s
outside the timed step                                      [.]     0.2 s
```

About 9 to 10 s of the 23.4 s have the GPU idle or almost idle. That is about 40 %.

| Idle cause | Seconds per step | Type |
|---|---|---|
| Long tail of generation: 1 or 2 answers decode alone | 4-5 | estimate |
| CPU work between generation and the reference pass | 3.0-3.2 | estimate from 2 runs |
| Weight push to vLLM, wake, sleep | 1.7 | measured |
| Prompt work and other CPU work | 0.5-1.0 | estimate |

### 5.2 The long tail

In the R5 run, generation time = 5.61 s + 4.3 ms x the longest answer. The mean answer has 131 tokens. The longest
answer in a step has 1,271 tokens on average. In 38 % of R5 steps, one answer reaches the limit of 1,536 tokens. The trainer then removes
that answer from the loss. Thus the GPU decodes up to 1,400 tokens for nothing.

**Worked example 4: one runaway answer (full run).** In the full run, generation time = 6.06 s + 4.1 ms x the
longest answer. All finished answers are below 600 tokens. One answer runs to 1,536. The
extra time is 4.1 ms x (1,536 - 600) = 3.8 s. In the full run this happened in 19 % of steps. The mean loss is about
0.7 s for each step.

### 5.3 The data that moves

| What moves | Size or count per step | Problem |
|---|---|---|
| Image pixels (fp32) | 1.61 GB | Built twice on the CPU. Copied in a blocking way. |
| Weights to vLLM | about 700 calls, 4.26 GB | Many small calls. Vision weights do not change, but they move. |
| Small host waits | about 1,250 | Each wait stops the GPU queue. |
| Position numbers of each token | 256 answers | Calculated on the CPU in a Python loop. |

All 7,341 LIDC training images are 512 x 512. Thus every prompt has the same tokens and the same positions. The
system can calculate them once.

### 5.4 The checker is not the problem

In each step, the trainer waits about 5 ms for the checker to score 256 answers. The checker runs on its own streams,
beside other work. In a speed test with large batches (1,024 outlines of 16 corners), the checker checks 7.74
million outlines each second on the RTX 5090.

### 5.5 The polygonizer rule is broken

The CUDA polygonizer exists in the checker. It matches the CPU polygonizer on 195,809 of 195,809 cases. But the
target builder and the certificate code still call the CPU version. Also, the CUDA version does too much work. It
visits every possible edge place in the mask, not only the true edges. For a 256 x 256 mask, it does about 500 times
the work that is necessary.

---

## 6. Finding 3: the RL system against the world

### 6.1 On one GPU, our system is near the best

We compared 16 RL frameworks and modes: TRL, verl, NeMo-RL, slime, Miles, OpenRLHF, AReaL, ROLL, SkyRL, prime-rl, Unsloth,
EasyR1, LlamaFactory, Oat, torchforge, and TRL's new async trainer.

- On one GPU, all frameworks use the same loop: generate, push weights, train.
- Our trainer already has the advanced parts: the shared-prompt forward pass, the image-feature cache, the
  position cache, no extra old-probability pass, and the probability correction.
- A move to another framework on one GPU gives about the same step time (estimate: 20-25 s).

### 6.2 What we do not have: generation and training at the same time

The best systems generate the next batch while they train on the current batch. All of them need two or more GPUs
for this. With 2 GPUs, the estimate is 13-14 s for each step. TRL's async trainer cannot use images yet. Thus for
async RL with images, the choice is verl or NeMo-RL.

### 6.3 The best update methods for a smooth reward

| Method | Use for this project |
|---|---|
| Dr. GRPO (in use) | Good. No length bias. |
| Scale by group spread | Bad for IoU. A group with IoU 0.61-0.63 gets the same weight as a group with 0.2-0.9. |
| CISPO, GSPO | Best in a study of 400,000 GPU-hours (ScaleRL). Use them with async RL. |
| GDPO | Scales each reward part alone. Try it against the current sum. |
| KL-Cov, Clip-Cov | Use if entropy falls too low on hard images. |

### 6.4 The maths of the reward

Agent H2 proved these points:

1. **Any reward that scores one answer at a time teaches the model to write the same outline every time.** The fall
   of entropy is the expected result, not an error. The KL weight is 0.001 and IoU has a weight of about 6. Thus a
   gain of 0.01 IoU multiplies the odds of an outline by e^(6 x 0.01 / 0.001) = e^60.
2. **The objective changed in the middle of the run.** The reward divides each part by its spread in the group.
   The spread has a fixed lower limit. Early in the run, the spreads are large. Late in the run, most spreads are
   below the limit. Then each part is multiplied by a fixed number. Agent H2 calculated a sum of about
   6 x IoU + 6 x boundary F + 6 x area.
3. **For the current metric, use a simpler reward.** Use the legality gate times (IoU plus a small location term for
   misses), plus topology. Subtract the mean of the other 7 answers. Do not divide by the spread.
4. **For a model that gives a range of plausible outlines, use a GED reward.** GED with d = 1 - IoU is strictly
   proper. The diversity weight must be 1, not 0.5.
5. **The shift problem has one correct fix.** Training images and test images overlap too little to re-weight them
   (0.1 % against 70 %). The correct fix is to train on shifted images. The code already has this fix.

---

## 7. Finding 4: the best algorithms

Agent H1 checked each geometry part. It found no error in the maths. All results are integers or correctly rounded
operations on integers. Thus any correct parallel order gives the same bits.

| Part | Today | Best design |
|---|---|---|
| Polygonizer | CPU in production. CUDA version visits all possible edge places. | One GPU block for each mask. Emit only true edges. Find the ring start in closed form from the component labels. |
| Connected components | One thread for each row. Always 4 passes. | Skip passes when the size filter changes nothing. Label runs of pixels in bit words. |
| Shape validity | One thread for each shape. Can take about 0.1 s for a bad shape. | A grid to find candidate edge pairs. One warp for each shape. |
| Distance transform | CPU in evaluation | Exact separable integer transform on the GPU |
| 3D topology | CPU only | Boundary faces from the same 2 x 2 table. Exact b1 from the Euler number. |

**The key idea for the polygonizer.** Each ring is the border between one object part and one background part. Its
start is the top-left corner of the first pixel of that part. The checker already finds these parts. Thus the ring
search needs no repeated jumps. For a ring with 1,000 edges, the work falls from about 10,000 operations to about
1,000.

---

## 8. Finding 5: models, novelty, 3D, CAD and Physical AI

### 8.1 Which model

We read the Hugging Face leaderboards and the model cards.

| Model | Size | RefCOCO mean | Cost to change |
|---|---|---|---|
| Qwen3-VL-2B (today) | 2.13 B | 84.8 | none |
| Qwen3-VL-4B | about 4 B | 88.2 | Configuration only. Same structure. About 2 times the step time. |
| Qwen3.5-2B | 2.27 B | 84.8 | High. New linear-attention layers break our trainer's shortcuts. |
| Qwen3.5-4B | 4.66 B | 88.1 | High, same reason. |

The best next model is Qwen3-VL-4B. Qwen3.5-2B gives no gain in grounding.

### 8.2 What is new in this project

The nearest published work is SimpleSeg (January 2026). It also writes outlines as text and uses an IoU reward. Thus
we cannot say "first RL for models that write outlines". We found no published work for these four points:

1. A set of doctors as the RL reward. Other work uses one mask or a supervised loss.
2. Exact Betti numbers as a reward for a model that draws outlines.
3. A medical model that writes outlines directly. Other medical work uses SAM as a mask tool.
4. An exact, bit-identical GPU checker as a reward service.

Also: medical models cannot find objects well. On one test, MedGemma-4B gets box IoU 0.05. A fixed centre box gets
0.10.

### 8.3 Our numbers need fair comparison

Published LIDC Dice runs from 0.49 to 0.92. The value changes with the test method. Thus we must run the main
baselines on our own splits: nnU-Net 2D, MedSAM2, SAM 3, and SimpleSeg. We must also report GED, Dice against the
majority and HM-IoU. Our 8 answers for each image give these values at no extra cost.

### 8.4 3D, CAD and Physical AI

| Area | What we found | Real or hype |
|---|---|---|
| 3D CT | The checker can grow to voxels. Today 3D topology runs only on the CPU. Two errors: a skipped slice splits a nodule into two parts, and no slice thickness is stored. | Real. Fix first. |
| 3D cavities | Only 29 of 2,654 nodules have a real cavity after the size filter. A b2 reward mostly removes false cavities. | Small effect |
| CAD | All CAD RL papers score with approximate CPU geometry. One paper estimates IoU from 8,192 random points. Its noise is about the size of its whole RL gain. | Real gap |
| CAD to 3D | Solved tools exist. The open problem is the reverse: image to CAD. | Aim at image to CAD |

**Owner decision (2026-11-08).** The main direction is 3D voxels. CAD was a study only, and it is parked. The CAD
findings stay here as a record.
| Physical AI video | Every video test first lifts geometry with a learned model. Exactness is lost. | Hype today |
| Sim-ready assets | NVIDIA's own list asks for watertight meshes but has no checker for it. A better collision shape raised one robot task from 49 % to 80 % success. | Real gap |

---

## 9. Facts that the main agent checked again

### 9.1 Checked in the run files and code

| Fact | Value | Where |
|---|---|---|
| Test IoU, median, boundary F, Dice, legal | 0.645, 0.756, 0.816, 0.701, 99.4 % | final test results file |
| Lowest correction ratio is 0 | 98.1 % of 770 steps | full RL log |
| Mean correction ratio | 0.994, 0.704, 0.847 | full RL log (steps 1-10, 100-199, 700-770) |
| Entropy | 0.36 to 0.067 | full RL log |
| Training IoU | 0.523 to 0.771 | full RL log |
| Mean time for each step in the full run | 27.47 s | full RL log |
| GPU money spent | $16.90 in 27 rental records; limit $17.50 | spend ledger, money rule in the host policy code |
| KL weight, length limit, learning-rate rule | 0.001, 1,536 tokens, constant | training configuration code |
| CUDA polygonizer exists | yes, in the checker | checker source |
| Production polygonizer callers | target builder and certificate code, both CPU | target and certificate source |

### 9.2 A correction during the study

The main agent first read 88.2 as the RefCOCO score of Qwen3.5-2B. Agent E checked the columns of the model card.
The 88.2 is the score of Qwen3-VL-4B. Qwen3.5-2B scores 84.8. We corrected the shared file before any report used
it.

---

## 10. Limits of this study

- We did not run anything on a GPU. All speed gains are estimates.
- The test view is not on the Mac. Some tables use the validation split.
- The final RL model's single answers were not saved. Thus we cannot measure best-of-8 yet.
- We did not install any NVIDIA skill or any new package.

---

## 11. Source files

All findings files are in the session scratchpad. The main facts come from these repo files:

- Final test results: `runs/r5/eval/eval_test/final_eval.json`
- Full RL log: `runs/rl/train_rl/rl.json`
- Step timing: `runs/r5/shared_rl.json`
- Configuration: `src/objdel/train/configs.py`
- Checker and CUDA polygonizer: `cuda/gpucheck.cpp`
- Production polygonizer callers: `src/objdel/targets/build.py`, `src/objdel/checker/certify.py`

---

## 12. Appendix: papers read

Reading agents read these papers after the main study, to check the findings against the newest work. Each row
gives the paper, how the agent read it, the main measured result, and its use for this project.

- Papers: 104 unique papers, all from 2025 or 2026 (51 from 2026).
- "full" means that the agent fetched the full text and read the abstract, method and results passages. It did not
  always read every page or every appendix. Where a table lost its numbers in the text conversion, the agent used
  only the numbers in the prose. Read in this way: 74 papers.
- "abstract + tables" means that the agent read the abstract and the result tables or result sentences only. Read
  in this way: 30 papers.
- The agents wrote no number that is not in the paper. "not stated" means that the paper text did not give it.

### 12.A RL algorithms and RL systems (42 papers)

| # | Paper | Read | Key result | Use for us |
|---|---|---|---|---|
| 1 | [2503.14476](https://arxiv.org/abs/2503.14476): DAPO: An Open-Source LLM Reinforcement Learning System at Scale (March 2025) | full | Qwen2.5-32B base, AIME 2024: plain GRPO reaches 30. DAPO reaches 50, with 50 % of the training steps of DeepSeek-R1-Zero-Qwen-32B. | Dynamic sampling drops groups with zero gradient. It saves compute on all-miss groups but gives no signal. Clip-Higher is the cheap first test against entropy collapse. |
| 2 | [2503.20783](https://arxiv.org/abs/2503.20783): Understanding R1-Zero-Like Training: A Critical Perspective (March 2025) | abstract + tables | 7B base, AIME 2024: 43.3 % with the minimalist recipe. Dr.GRPO improves token efficiency with the same reasoning accuracy. | Std normalization blows up tiny reward differences. In an all-miss group with near-zero IoU, Dr.GRPO keeps the advantage small instead of large. Test it as a one-line change. |
| 3 | [2505.22617](https://arxiv.org/abs/2505.22617): The Entropy Mechanism of Reinforcement Learning for Reasoning Language Models (May 2025) | full | Qwen2.5, DAPO-MATH training, average score (Table 2): 7B GRPO 38.6, Clip-Cov 40.4, KL-Cov 40.6. 32B GRPO 45.8, Clip-Cov 50.3, KL-Cov 52.2. Entropy stays higher. | Direct fix for our entropy collapse. We can log the covariance between log-probability and advantage and apply KL-Cov to the top 0.2 % tokens (paper setting for 32B: k = 2e-4). |
| 4 | [2506.01939](https://arxiv.org/abs/2506.01939): High-Entropy Minority Tokens Drive Effective Reinforcement Learning for LLM Reasoning (June 2025) | abstract + tables | Only the 20 % tokens: Qwen3-32B +11.04 on AIME'25 and +7.71 on AIME'24 over full-gradient updates. Qwen3-14B +4.79 and +5.21. Qwen3-8B stays comparable. | In polygon text, the fork tokens are probably the vertex coordinates. We can restrict the gradient to high-entropy tokens and see if entropy stays alive. |
| 5 | [2505.24864](https://arxiv.org/abs/2505.24864): ProRL: Prolonged Reinforcement Learning Expands Reasoning Boundaries in Large Language Models (May 2025) | full | Nemotron-Research-Reasoning-Qwen-1.5B vs DeepSeek-R1-1.5B: average pass@1 +14.7 % (math), +13.9 % (coding). Gains continue after 2,000 steps. | A small 1.5B model can keep improving. The KL-reset trick gives a plan if our policy drifts from the reference after many steps. |
| 6 | [2508.10751](https://arxiv.org/abs/2508.10751): Pass@k Training for Adaptively Balancing Exploration and Exploitation of Large Reasoning Models (August 2025) | abstract + tables | Gain is shown on Enigmata with Qwen2.5-7B-Instruct. I did not extract a single headline number. | Our reward is a best-of-group style check on a set of doctors. A pass@k-shaped advantage can reward a group where one answer is a hit. Low priority until entropy is fixed. |
| 7 | [2507.18071](https://arxiv.org/abs/2507.18071): Group Sequence Policy Optimization (July 2025) | full | More stable than GRPO and stabilizes MoE training. No single number extracted. The paper says sequence likelihoods tolerate train-inference precision differences better, so it can skip recomputing old-policy likelihoods. | Direct help for the mismatch problem. SimpleSeg (2601.19228) trains a polygon-output model with GSPO. |
| 8 | [2511.20347](https://arxiv.org/abs/2511.20347): Soft Adaptive Policy Optimization (November 2025) | full | Better Pass@1 than GSPO and GRPO under the same budget on math. It also improves Qwen3-VL-30B-A3B training (AIME25, LiveCodeBench v6, ZebraLogic, MathVision). No numbers extracted. | A hard clip kills our whole sequence when a few tokens mismatch. SAPO keeps the good tokens. It was tested on a VLM family. |
| 9 | [2507.20673](https://arxiv.org/abs/2507.20673): Geometric-Mean Policy Optimization (July 2025) | abstract + tables | GMPO-7B improves average Pass@1 over GRPO by up to 4.1 % on math. | One-line way to damp the extreme ratios we see. Compare with GSPO and SAPO. |
| 10 | [2602.05261](https://arxiv.org/abs/2602.05261): Revealing and Controlling Response Length Variation in RLVR (February 2026) | full | Qwen2.5-VL-7B-Instruct: LUSPO beats GSPO by 5.1 % (WeMath) and 6.0 % (LogicVista). | Polygon length matters for us. If we adopt GSPO, use the length-unbiased form so the model does not shrink to short polygons. |
| 11 | [2501.03262](https://arxiv.org/abs/2501.03262): 2501.03262 — Stabilizing Critic-Free Policy Optimization with Global Normalization (January 2025; REINFORCE++) () | full | Chat-Arena-Hard: single-sample variant 46.7 vs GRPO 46.8, with shorter answers (832 vs 860 tokens). On a small dataset GRPO with local normalization overfits completely. | Batch-level normalization keeps a small but real signal from all-miss groups. ScaleRL (2510.13786) also uses batch-level advantage normalization. |
| 12 | [2601.05242](https://arxiv.org/abs/2601.05242): GDPO: Group reward-Decoupled Normalization Policy Optimization for Multi-reward RL Optimization (January 2026) | full | Qwen2.5-Instruct-1.5B tool calling: about +2.7 % average accuracy and over 4 % in format ratio over GRPO. Up to 80 % less length-exceeding on AIME. | Our reward has IoU, boundary F and Betti topology. Normalize each term per group before the sum so a rare topology term is not drowned. |
| 13 | [2503.06639](https://arxiv.org/abs/2503.06639): GRPO's Effective Loss, Dynamics, and Success Amplification (March 2025) | abstract + tables | Proof only. No benchmark number. | Explains why GRPO needs some successes in the group. All-miss groups do not move. Background for the reward-signal fix. |
| 14 | [2504.05118](https://arxiv.org/abs/2504.05118): VAPO: Efficient and Reliable Reinforcement Learning for Advanced Reasoning Tasks (April 2025) | abstract + tables | Qwen2.5-32B, AIME 2024: 60.4. More than 10 points above DAPO and R1-Zero-Qwen-32B. 5,000 steps, no crashes across runs. | A value model costs memory on one RTX 5090. Does not apply now. Keep as a fallback if GRPO stays flat. |
| 15 | [2501.12599](https://arxiv.org/abs/2501.12599): 2501.12599 — Scaling Reinforcement Learning with LLMs (January 2025; Kimi k1.5) () | full | AIME 77.5, MATH-500 96.2, MathVista 74.9. Prioritized sampling: pick problem i with weight 1 - s_i (s_i = its success rate). | Prioritized sampling sends compute to hard crops. Partial rollouts address the long-tail idle time. |
| 16 | [2512.23075](https://arxiv.org/abs/2512.23075): Trust Region Masking for Long-Horizon LLM Reinforcement Learning (December 2025) | abstract + tables | Theory with a stability demonstration. No number extracted. | Gives a principled rule for the mismatch: mask the whole sequence if any token diverges too much, and log the acceptance rate as a diagnostic. |
| 17 | [2602.10693](https://arxiv.org/abs/2602.10693): Variational Sequence-Level Soft Policy Optimization for Stable Off-Policy LLM Training (February 2026) | abstract + tables | Stable up to 64x staleness and fully asynchronous runs, on dense and MoE models. | Another way to tame ratios without clipping. Useful if we later go async. |
| 18 | [2510.26788](https://arxiv.org/abs/2510.26788): Defeating the Training-Inference Mismatch via FP16 (October 2025) | full | FP16 mismatch is about 24x smaller than BF16. Plain importance-weighted policy gradient in FP16 beats all BF16 baselines (GRPO, GSPO, TIS, MIS) in the paper's sanity test. LoRA and larger models are in the tested set. | Cheapest test for the "IS ratio min 0 in 98 % of steps" problem. Run the trainer and vLLM in FP16 and compare the ratio logs. |
| 19 | [2609.20807](https://arxiv.org/abs/2609.20807): Score Centering Stabilizes Off-policy Reinforcement Learning (September 2026) | abstract + tables | 0.6B to 30B models. It matches or beats importance sampling under quantization. Gap grows with mismatch (Qwen3-30B-A3B-Base, INTELLECT-2 math data). | Works when mismatch cannot be removed. It composes with importance sampling. Try after the FP16 test. |
| 20 | [2601.18150](https://arxiv.org/abs/2601.18150): FP8-RL: A Practical and Stable Low-Precision Stack for LLM Reinforcement Learning (January 2026) | full | Up to 44 % rollout throughput gain. FP8 without TIS loses accuracy against the BF16 baseline. FP8 with TIS matches it. | Confirms that a rollout correction is required when the sampler is quantized. Speed gain is for long contexts. Our sequences are short, so expect less. |
| 21 | [2510.18855](https://arxiv.org/abs/2510.18855): 2510.18855 — Every Step Evolves: Scaling Reinforcement Learning for Trillion-Scale Thinking Model (October 2025; Ring-1T) () | full | Ring-1T: AIME-2025 93.4, HMMT-2025 86.72, CodeForces 2088, ARC-AGI-1 55.94. | IcePop is a token-level mask for the mismatch. It is the closest fix for tokens with ratio 0. |
| 22 | [2506.13585](https://arxiv.org/abs/2506.13585): MiniMax-M1: Scaling Test-Time Compute Efficiently with Lightning Attention (June 2025) | full | On Qwen2.5-32B, CISPO gives a 2x speedup over DAPO. | Two cheap items: CISPO as the loss, and an FP32 LM head in both trainer and vLLM. |
| 23 | [2605.14220](https://arxiv.org/abs/2605.14220): 2605.14220 — Diagnosing Training Inference Mismatch in LLM Reinforcement Learning (May 2026) [new] () | full | Experiments on Qwen3-1.7B GRPO and a R1D-1.5B sanity test. Plain recompute or bypass fail while the exact engine works. I did not extract a single headline number. | Tells us to treat our 98 % min-0 finding as a first-order bug, not noise. Gives a method: build a zero-mismatch reference and compare. |
| 24 | [2602.19362](https://arxiv.org/abs/2602.19362): 2602.19362 — LLMs Can Learn to Reason Via Off-Policy RL (February 2026) [new] () | abstract + tables | Beats GRPO with importance sampling on competition math. Matches DeepCoder on LiveCodeBench with 3x fewer generations. Works with lags above 400 gradient steps. | If the ratio is unusable (min 0), this removes the ratio altogether. |
| 25 | [2510.13786](https://arxiv.org/abs/2510.13786): The Art of Scaling Reinforcement Learning Compute for LLMs (October 2025) | full | More than 400,000 GPU-hours study. ScaleRL asymptotic reward A = 0.61. One run scaled to 100,000 GPU-hours and the curve predicted it. CISPO and GSPO give a higher asymptote than DAPO. FP32 at the LM head raises the asymptote. | Best single recipe checklist. Most items are small code changes in TRL. Also gives a way to compare our runs by curve fit, not by one step. |
| 26 | [2505.24298](https://arxiv.org/abs/2505.24298): AReaL: A Large-Scale Asynchronous Reinforcement Learning System for Language Reasoning (May 2025) | full | Up to 2.77x training speedup over synchronous systems with the same GPU count, with equal or better final score. The paper says the gain is smaller at short context length. | Our idle GPU is the sync pattern. But we have one GPU, so full decoupling does not fit. The short-context caveat applies to us. |
| 27 | [2509.19128](https://arxiv.org/abs/2509.19128): 2509.19128 — Faster On-policy Reinforcement Learning for Long Sequence Generation (September 2025; PipelineRL) () | full | About 2x faster learning than a conventional baseline. The abstract says 128 H100 GPUs. The introduction says 4 DGX-H100 nodes. Data stays highly on-policy. | Pattern for removing idle time while staying near on-policy. Needs more than one GPU or a clean split of the one GPU. |
| 28 | [2506.05433](https://arxiv.org/abs/2506.05433): Prefix Grouper: Efficient GRPO Training through Shared-Prefix Forward (June 2025) | full | Same outputs and gradients as GRPO. Lower GPU memory at every tested group size. No speed number extracted. | Already in our run as the shared-prefix forward. This paper is the reference for the equivalence proof. |
| 29 | [2509.23232](https://arxiv.org/abs/2509.23232): SPEC-RL: Accelerating On-Policy Reinforcement Learning via Speculative Rollouts (September 2025) | full | Rollout time cut by 2 to 3x with no loss in score (Qwen-3-8B-Base, PPO, GRPO, DAPO). | We revisit the same crops across epochs. Our outputs are short, so check the gain before building. |
| 30 | [2509.21792](https://arxiv.org/abs/2509.21792): FastGRPO: Accelerating Policy Optimization via Concurrency-aware Speculative Decoding and Online Draft Learning (September 2025) | abstract + tables | End-to-end speedup of 2.35x to 2.72x on math datasets. | Needs a draft model. Poor fit for a 2B model with short outputs. Low priority. |
| 31 | [2503.06520](https://arxiv.org/abs/2503.06520): Seg-Zero: Reasoning-Chain Guided Segmentation via Cognitive Reinforcement (March 2025) | full | Seg-Zero-7B reaches 57.5 on ReasonSeg zero-shot, 18 % above LISA-7B. Soft reward gives a small gain on ReasonSeg and a clear loss on RefCOCOg. A strict format reward helps out-of-domain. | Supports a hard-threshold reward next to our soft IoU. Test both. Strict format reward is good for polygon syntax. |
| 32 | [2504.07954](https://arxiv.org/abs/2504.07954): Perception-R1: Pioneering Perception Policy with Reinforcement Learning (April 2025) | full | Qwen2-VL-2B-Instruct: +4.2 % RefCOCO+, +17.9 % PixMo-Count, +4.2 % PageOCR, 31.9 AP on COCO2017 val. | Same model size as ours. Matching is the right pattern for scoring against a set of doctors. |
| 33 | [2505.12081](https://arxiv.org/abs/2505.12081): VisionReasoner: Unified Reasoning-Integrated Visual Perception via Reinforcement Learning (May 2025) | abstract + tables | Relative to Qwen2.5-VL: +29.1 % COCO, +22.1 % ReasonSeg, +15.3 % CountBench. | Confirms that rewards on box and point accuracy transfer. Little new for polygons. |
| 34 | [2505.22596](https://arxiv.org/abs/2505.22596): SAM-R1: Leveraging SAM for Reward Feedback in Multimodal Segmentation via Reinforcement Learning (May 2025) | abstract + tables | Strong results with only 3k training samples. Exact numbers not extracted. | Shows RL needs few samples. Our exact CUDA checker plays the role of SAM here and is cheaper. |
| 35 | [2506.22624](https://arxiv.org/abs/2506.22624): Seg-R1: Segmentation Can Be Surprisingly Simple with Reinforcement Learning (June 2025) | full | 0.873 S-measure on COD10K. 71.4 cIoU on RefCOCOg test and 56.7 gIoU on ReasonSeg test, zero-shot. S-measure alone led to reward hacking. The fix: 0.7 IoU + 0.3 S-measure. | Direct warning for our multi-term reward. Watch for hacking in the boundary F and Betti terms. Keep IoU as the main weight. |
| 36 | [2509.06321](https://arxiv.org/abs/2509.06321): Text4Seg++: Advancing Image Segmentation via Generative Language Modeling (September 2025) | full | Row-wise RLE cuts descriptor length by 74 % and speeds inference 3x with no loss. | Gives a text-output precedent. Shows output length cuts cost. Our polygons are already short. |
| 37 | [2512.24023](https://arxiv.org/abs/2512.24023): RSAgent: Learning to Reason and Act for Text-Guided Segmentation via Multi-Turn Tool Invocations (December 2025) | abstract + tables | 66.5 % gIoU on ReasonSeg test (+9 % over Seg-Zero-7B). 81.5 % cIoU on RefCOCOg. | Multi-turn refinement is costly. Read together with 2606.13156 before trying. |
| 38 | [2606.13156](https://arxiv.org/abs/2606.13156): Iterative Visual Thinking and the Self-Correction Mirage in VLM Grounding (June 2026) | full | Qwen3-VL-4B: stopping at step 0 matches the base and beats every shippable stopping rule. The gain size is not legible in my extract. | Do not add a self-correction loop to fix shifted crops. Score only with a fixed rule that does not use the label. |
| 39 | [2608.13549](https://arxiv.org/abs/2608.13549): Exponential Convex Calibration Dimension for the Multi-Label Jaccard Measure (August 2026) | abstract + tables | Proofs only. No benchmark number. | Supports using the exact IoU as a reward instead of a smooth surrogate loss. No code change. |
| 40 | [2601.19228](https://arxiv.org/abs/2601.19228): 2601.19228 — Towards Pixel-Level VLM Perception via Simple Points Prediction (January 2026; SimpleSeg) () | full | Reward table (RefCOCO / RefCOCO+ / RefCOCOg gIoU): IoU 76.9 / 71.9 / 72.9. IoU + distance 77.1 / 72.2 / 73.1. Adding a hard length penalty drops to 66.7 / 62.4 / 62.0. Polygon length sweet spot: 221 tokens best, 35.6 cIoU at 78 tokens, 72.5 cIoU at 859 tokens. Pre-training raises gIoU by 4.6 (SFT) and 13.0 (SFT+RL) on RefCOCO. | Closest match to our task (polygon text from a VLM with IoU reward). Warns against a hard length penalty. Shows a point-count sweep we can copy. |
| 41 | [2605.15951](https://arxiv.org/abs/2605.15951): From Failure to Feedback: Group Revision Unlocks Hard Cases in Object-Level Grounding (May 2026) | full | Gains over GRPO-based models on referring and reasoning segmentation, REC and counting. No single number extracted. | Aimed at our all-miss problem. A revision gain gives signal when every raw IoU is near zero. Check cost: it adds rollouts. |
| 42 | [2601.07695](https://arxiv.org/abs/2601.07695): Smooth Operator: Smooth Verifiable Reward Activates Spatial Reasoning Ability of Vision-Language Models (January 2026) | abstract + tables | AP-GRPO matches large supervised methods with higher data efficiency on 3D tasks (Numerical3D-50k). No number extracted. | Names our failure: near-zero IoU groups lose signal after normalization. A dense reward plus an absolute scale term is a small change. |

### 12.B Vision models that write shapes, medical and multi-doctor segmentation (33 papers)

| # | Paper | Read | Key result | Use for us |
|---|---|---|---|---|
| 1 | [2604.02593](https://arxiv.org/abs/2604.02593): 2604.02593 — Moondream Segmentation: From Words to Masks (Apr 2026) [new addition] () | full | cIoU 80.2 % on RefCOCO val and 62.6 % mIoU on LVIS val. The RL reward is piecewise: box IoU first, then mask IoU, with a boundary-IoU band (width 5 % of image diagonal). It penalizes false negatives more than false positives. Invalid traces are filtered from the group. | test a false-negative-heavy IoU term plus boundary IoU in our reward, and filter invalid polygons out of the GRPO group instead of giving them a score. |
| 2 | [2507.04664](https://arxiv.org/abs/2507.04664): VectorLLM: Human-like Extraction of Structured Building Contours via Multimodal LLMs (Jul 2025) | full | +5.6 AP on WHU, +7.1 AP on WHU-Mix, +13.6 AP on CrowdAI over the prior best. With oracle boxes, WHU AP is 90.2. DPO adds 1.8 AP and 1.4 AP75. DPO pairs come from SFT outputs with IoU below 0.8, plus labels with random point deletion or insertion. | a cheap DPO round on our model's own bad polygons (IoU below 0.8, or corrupted ground truth) is a fallback if GRPO stalls. |
| 3 | [2608.03471](https://arxiv.org/abs/2608.03471): Hi-Token: Hierarchical Coordinate Tokenization for Generative Visual Grounding (Aug 2026) | full | gains on three backbones in matched comparisons. The paper reports that Hi-Token beats an extra-tuned flat baseline on RefCOCO before RL. On small boxes Hi-R1 gets mIoU 58.0 and P@0.9 14.7 (375 small test boxes). A one-bin coordinate shift hurts small boxes more than large ones. I did not recover the exact backbone gain table. | small nodules are sensitive to one grid step on the 0-1000 grid. Test finer grid digits or a coarse-then-fine coordinate order, and a validity gate in the reward. |
| 4 | [2603.28069](https://arxiv.org/abs/2603.28069): MolmoPoint: Better Pointing for VLMs with Grounding Tokens (Mar 2026) | abstract + tables | 70.7 % on PointBench, 61.1 % on ScreenSpot-Pro, +6.3 % on Molmo2Track. The paper also claims higher sample efficiency (no number captured). | a "no more vertices" token and a fixed vertex order are cheap to copy. The patch-selection head needs architecture changes, which our text-only LoRA plan avoids. |
| 5 | [2510.03230](https://arxiv.org/abs/2510.03230): Improving GUI Grounding with Explicit Position-to-Coordinate Mapping (Oct 2025) | abstract + tables | consistent gains on ScreenSpot, V2 and Pro, largest on high-resolution screens. No headline number captured. The paper says it does not reach state of the art. | a grid-hint idea for the 0-1000 grid (print a ruler row in the prompt). Low priority, needs a test on our crops. |
| 6 | [2601.13633](https://arxiv.org/abs/2601.13633): EGM: Efficient Visual Grounding Language Models (Jan 2026) | full | Qwen3-VL-2B-Thinking goes from 83.6 to 85.6 (SFT) to 89.6 (RL) on the RefCOCO family average; RL adds 4.0 for 2B versus 2.6 for 8B. EGM-8B reaches 91.4 IoU at 737 ms versus 90.5 IoU at 4,320 ms for Qwen3-VL-235B. Removing reasoning costs only 0.9 for 2B (88.7 vs 89.6). | same 2B backbone as ours. RL gains are largest at 2B, which supports our RL path. Reasoning text helps little at 2B, so skip it to save tokens. |
| 7 | [2602.09701](https://arxiv.org/abs/2602.09701): RL-Driven Vision–Language Grounding for Fine-Grained Referring Segmentation (Feb 2026) | full | 0.7127 cIoU on RefCOCOg val (+15.3 over the Qwen3-VL Instruct baseline, +3.3 over Seg-Zero-7B). On GRefCOCO, 82.40 % accuracy on no-target prompts with 0.5 % false-negative rate. Qwen3-VL uses a [0,1000] coordinate grid. | shows the no-target output ("empty polygon") can be trained by RL. Our nodule crops always contain a nodule, but the RIGA cup can be absent, so test an explicit empty answer. |
| 8 | [2605.12953](https://arxiv.org/abs/2605.12953): Seg-Agent: Test-Time Multimodal Reasoning for Training-Free Language-Guided Segmentation (May 2026) | full | +5.1 gIoU over the single-step baseline on reasoning segmentation. It matches the trained method Seg-Zero (Qwen2.5-VL-7B + SAM2-L). I did not capture the absolute numbers. | a render-and-reselect pass is a test-time option. Draw the model's polygon on the crop and ask a second pass to accept or fix it. |
| 9 | [2606.09303](https://arxiv.org/abs/2606.09303): Reason Twice: Segmentation via Candidate Discovery and Comparative Reasoning (Jun 2026) | full | the paper reports gains on its new ReasonSeg-SGDR benchmark and on ReasonSeg. I did not recover the numbers. | a test-time selector over K sampled polygons. Our checker can score each candidate exactly, so use it as the selector on the dev split to find the upper bound of best-of-K. |
| 10 | [2603.14382](https://arxiv.org/abs/2603.14382): StAR: Segment Anything Reasoner (Mar 2026) | full | StAR-7B with mask-level majority voting is "overall comparable" to the SAM 3 Agent (72B) on the ReasonSeg-X test set (paper text). Trained with 5k samples. ReasonSeg-X has 1,169 samples (240/156/773 train/val/test). | our best test-time candidate. Sample K polygons, cluster by IoU, return the medoid of the largest cluster. No new training. |
| 11 | [2503.01785](https://arxiv.org/abs/2503.01785): Visual-RFT: Visual Reinforcement Fine-Tuning (Mar 2025) | abstract + tables | +24.3 % one-shot fine-grained classification (about 100 samples); +21.9 on COCO two-shot; +15.4 on LVIS. | early proof that an IoU reward works with few samples. Our reward is already stronger. |
| 12 | [2504.07615](https://arxiv.org/abs/2504.07615): VLM-R1: A Stable and Generalizable R1-style Large Vision-Language Model (Apr 2025) | full | Qwen2.5-VL-3B RL scores 63.16 on out-of-domain LISA-Grounding (SFT 54.82). On OVDEval it reaches 31.01 nms-AP (SFT 26.50). | RL transfers better than SFT to new data, which fits our off-centre shift. Watch for reward hacking (for example many tiny polygons); log failure cases. |
| 13 | [2508.14153](https://arxiv.org/abs/2508.14153): LENS: Learning to Segment Anything with Unified Reinforced Reasoning (Aug 2025) | abstract + tables | 81.2 % average cIoU on RefCOCO/+/g with Qwen2.5-VL-3B, up to +5.6 over GLaMM. | shows a 3B VLM can match larger models with multi-level rewards. We already have box-free reward; no new action. |
| 14 | [2501.04001](https://arxiv.org/abs/2501.04001): Sa2VA: Marrying SAM2 with MLLM for Dense Grounded Understanding of Images and Videos (Jan 2025) | abstract + tables | Sa2VA-8B gets cIoU 81.9 / 76.5 / 78.9 on RefCOCO / + / g, +2.4, +3.9, +4.7 over GLaMM-7B. | the decoder route (token to SAM) is the alternative to our text polygons. Not needed unless polygons fail. |
| 15 | [2511.16719](https://arxiv.org/abs/2511.16719): SAM 3: Segment Anything with Concepts (Nov 2025) | abstract + tables | the paper claims it doubles prior accuracy on concept segmentation. On SA-Co/VEval SA-V it scores cgF1 30.3 against a human 53.1. | the verifier idea matches test-time selection. A generic SAM does not know nodule rater style. |
| 16 | [2511.21631](https://arxiv.org/abs/2511.21631): Qwen3-VL Technical Report (Nov 2025) | abstract + tables | no single grounding number used here. The report states the [0,1000] grid improves robustness to image size and aspect ratio. | confirms our 0-1000 grid is the native format of Qwen3-VL-2B. |
| 17 | [2504.09975](https://arxiv.org/abs/2504.09975): OctGPT: Octree-based Multiscale Autoregressive Models for 3D Shape Generation (Apr 2025) | abstract + tables | training 13x faster and generation 69x faster than prior autoregressive 3D models. | off topic. Coarse-to-fine tokens echo Hi-Token, nothing to run. |
| 18 | [2505.22914](https://arxiv.org/abs/2505.22914): cadrille: Multi-modal CAD Reconstruction with Online Reinforcement Learning (May 2025) | full | SFT-only CAD-Recode has CC3D IoU about 60 % and 10 % invalid outputs. The cadrille reward is 10 x IoU minus 10 for an invalid program. Training uses only hard examples (mean SFT reward below 7.5 over 3 samples). Online GRPO beats offline DPO. | a model for our reward: scale IoU, add a big negative for invalid polygons, and mine hard crops (small off-centre nodules) for RL. |
| 19 | [2508.08177](https://arxiv.org/abs/2508.08177): MedReasoner: Reinforcement Learning Drives Reasoning Grounding from Clinical Thought to Pixel-Level Precision (Aug 2025) | full | overall IoU 32.42 with Lingshu-7B and the soft reward (IoU and pDice only); the strict reward scores 31.69, and a hard-threshold baseline scores 15.85. With reasoning, 32.42 vs 30.29 without. Seg-Zero-7B gets IoU 16.14. All models are weak on histology (11.66). | shows that soft (continuous) rewards beat strict ones in medical grounding. Our checker is already continuous. |
| 20 | [2603.05911](https://arxiv.org/abs/2603.05911): CORE-Seg: Reasoning-Driven Segmentation for Complex Lesions via Reinforcement Learning (Mar 2026) | full | 37.06 % mean Dice and 18.42 % failure rate on ComLesion-14K, +14.89 points over baselines. Random rationales drop Dice to 33.85 % and raise failures to 23.64 %. | a connectivity penalty already sits in a medical RL reward. Our Betti reward goes further. The failure rate (no valid mask) is a metric worth logging. |
| 21 | [2602.03320](https://arxiv.org/abs/2602.03320): MedSAM-Agent: Empowering Interactive Medical Image Segmentation with Multi-turn Agentic Reinforcement Learning (Feb 2026) | full | 6 modalities, 21 datasets. In the ablation with MedSAM2, Dice is 0.794 (0.772 without the improvement reward) and about 2.11 turns. Direct RL without SFT fails on CT and MRI. | RL alone struggles on low-contrast CT. Our SFT-then-RL order is right. A second "fix" turn is a later option. |
| 22 | [2603.28120](https://arxiv.org/abs/2603.28120): 2603.28120 — MedLoc-R1: Performance-Aware Curriculum Reward Scheduling for GRPO-Based Medical Visual Grounding (Mar 2026) [new addition] () | full | improves A@0.5 and A@0.8 and training stability over GRPO baselines on three medical grounding benchmarks, over 3 seeds. I did not capture the numbers. | direct fix for off-centre small nodules. Start with a loose reward and tighten as the group IoU rises. |
| 23 | [2605.26621](https://arxiv.org/abs/2605.26621): 2605.26621 — MedVol-R1: Reward-Driven Evidence Grounding for Volumetric Reasoning Segmentation (May 2026) [new addition] () | abstract + tables | GRPO adds +4.34 DSC (AbdomenCT-1K), +2.20 (CT-ORG), +9.25 (KiTS23) over SFT. Ablation: SFT-only 85.52 to full model 89.86 DSC. | the largest RL gain is on the hardest set (KiTS23, tumours), a hint that RL helps small lesions most. |
| 24 | [2604.27720](https://arxiv.org/abs/2604.27720): Why Does Grounding Hurt Medical VQA? Benchmarking, Diagnosis, and Fine-Tuning of Vision-Language Models (Apr 2026) | full | all off-the-shelf systems have mean IoU 0.05 to 0.24 on SLAKE grounding (center-box baseline 0.10). Answer-only SFT of Qwen-2.5-VL-7B gives 0/418 parseable boxes. Adding a 15 % mix of real box supervision gives IoU 0.390 (0.360 on unseen images). | confirms that off-the-shelf VLMs cannot do our task and that localization supervision must stay in the mix. Cropping hides context, which links to our centred-crop training gap. |
| 25 | [2504.05977](https://arxiv.org/abs/2504.05977): Diffusion Based Ambiguous Image Segmentation (Apr 2025) | full | GED (16 samples) 0.233 on 128x128 centred crops; on random 64x64 crops, GED 0.312 with 16 predictions (the best prior was 0.320) and dice 0.442 with 16 predictions (0.459 for 4 predictions with post-processing). The paper says centred LIDC is "too easy", because the lesion location is known. Typical failure: predictions that are too small. A relative-area post-process (threshold 0.5) helps. | this is the same effect as our off-centre problem. Report our metrics on this random-crop protocol too. Test a relative-area filter on our polygons. |
| 26 | [2512.00748](https://arxiv.org/abs/2512.00748): Disentangling Preference and Ambiguity for Controllable Medical Lesion Segmentation (Dec 2025) | full | LIDC-IDRI GED 0.1152, soft Dice 91.53 %, mean personalised Dice 90.25 % (per-rater 89.49 to 91.01). NPC GED 0.2272. The four LIDC "raters" are virtual: 12 real readers are assigned to 4 slots by area rank. | LIDC rater identity is not stable, so a rater-conditioned output is weak. Scoring against the SET of raters (our checker) is the safer choice. Use ProSeg's GED and Dice-max as extra metrics. |
| 27 | [2603.04024](https://arxiv.org/abs/2603.04024): Volumetric Directional Diffusion: Anchoring Uncertainty Quantification in Anatomical Consensus for Ambiguous Medical Image Segmentation (Mar 2026) | full | 3D LIDC-IDRI MaxDice 0.8031, GED 0.2081 (3D nnUNet MaxDice 0.7622). Without the anchor, MaxDice (16 samples) drops to 0.6670 and GED rises to 0.5809. Results hold from 2 to 4 raters (MaxDice about 0.803). | a consensus anchor is the right frame for set-of-raters scoring. Our model could predict a consensus polygon first, then variants. |
| 28 | [2605.08210](https://arxiv.org/abs/2605.08210): Harmonized Feature Conditioning and Frequency-Prompt Personalization for Multi-Rater Medical Segmentation (May 2026) | full | LIDC GED 0.105 (D-Persona 0.136); NPC-170 GED 0.175 (0.197). Personalised LIDC Dice-max 92.65 %, Dice-match 90.0 %, +1.61 mean Dice over D-Persona. Metrics given for K = 10 to 50 samples. | the SOTA table for LIDC multi-rater. Compare our Dice-max and GED against these numbers, with the caveat that they use 128x128 centred crops. |
| 29 | [2605.02437](https://arxiv.org/abs/2605.02437): Multi-Rater Calibrated Segmentation Models (May 2026) | full | lowest multi-rater ECE on all four sets (CoCaHis, REFUGE disc, REFUGE cup, LongCIU). REFUGE cup 6.67 vs 12.01 for random sampling; disc 1.09 vs 3.41 (values x 10^-2). AUC is unchanged (cup 92.57 vs 92.06). | uses REFUGE, the same family as RIGA. If we output a confidence per polygon, the ordinal-agreement target is the calibration baseline. |
| 30 | [2603.27197](https://arxiv.org/abs/2603.27197): 2603.27197 — K-alpha-LOS finds Consensus: A Meta-Algorithm for Evaluating Inter-Annotator Agreement in Complex Vision Tasks (Mar 2026) [new addition] () | full | on LIDC-IDRI (nodules at least 3 mm, four reading sessions treated as unnamed raters) alpha is stable at about 0.37 for localization thresholds 0.2 to 0.5, then drops to about 0.33 at 0.6. Raters agree where the nodule is but not how far it extends. | gives the evidence that location is consistent and extent is not. This supports a miss (wrong location) being a model error, not rater noise. |
| 31 | [2603.13054](https://arxiv.org/abs/2603.13054): Topo-R1: Detecting Topological Anomalies via Vision-Language Models (Mar 2026) | full | leading VLMs score near random. Topo-R1 (Qwen3-VL-4B) reaches F1@0.5 of 45.2 on the test set. With only 25 % of the SFT data it keeps over 91 % of the F1. The reward weights are 0.10 format, 0.85 accuracy, 0.05 topology. | the closest recipe for a Betti-based reward in GRPO. Verifiable Betti labels are the right target. Note the small topology weight (0.05). |
| 32 | [2603.19169](https://arxiv.org/abs/2603.19169): 2603.19169 — ARIADNE: A Perception-Reasoning Synergy Framework for Trustworthy Coronary Angiography Analysis (Mar 2026) [new addition] () | abstract + tables | centerline Dice 0.838 on 1,400 angiograms, 41 % fewer false positives than geometric baselines. External checks on ARCADE and XCAD. | Betti numbers work as preference labels (chosen vs rejected). This is a cheaper alternative to a Betti reward in GRPO. |
| 33 | [2504.03600](https://arxiv.org/abs/2504.03600): MedSAM2: Segment Anything in 3D Medical Images and Videos (Apr 2025) | abstract + tables | the human-in-the-loop pipeline cuts manual annotation cost by more than 85 % (5,000 CT lesions, 3,984 liver MRI lesions, 251,550 echo frames). | MedSAM2 is the segmenter used by MedReasoner and MedVol-R1. We could use it as an optional box-to-mask baseline on the same crops. |

### 12.C 3D, CAD, Physical AI and GPU methods (29 papers)

| # | Paper | Read | Key result | Use for us |
|---|---|---|---|---|
| 1 | [2501.19054](https://arxiv.org/abs/2501.19054): Text-to-CAD Generation Through Infusing Visual Feedback in Large Language Models (January 2025) | full | LLaMA-3-8B with LoRA on DeepCAD text pairs (4 A6000 GPUs): LVM score 8.96, invalidity ratio 6.20%, versus GPT-4o 5.13 and 37.93%. | It shows that a non-differentiable renderer can drive training, and that alternating reward steps with supervised steps stops format collapse; our exact checker should give a stronger reward than their VLM judge. |
| 2 | [2505.06507](https://arxiv.org/abs/2505.06507): Text-to-CadQuery: A New Paradigm for CAD Generation with Scalable Large Model Capabilities (May 2025) | full | Qwen2.5-3B reaches 69.3% Gemini-judged match (10.5% above the Text2CAD baseline) and volumetric IoU median 0.9868. The abstract states a 48.6% Chamfer Distance drop. Training is plain cross-entropy, no RL. | It gives a ready CadQuery dataset recipe for the SFT start of an image-to-CAD run, and shows that a small model suffices; its voxel IoU (0.02 resolution) is the metric our GPU verifier should replace exactly. |
| 3 | [2505.19713](https://arxiv.org/abs/2505.19713): CAD-Coder: Text-to-CAD Generation with Chain-of-Thought and Geometric Reward (May 2025) | full | Qwen2.5-7B-Instruct, 8 A800 GPUs: mean CD 6.54 (x10^3) and median CD 0.17, versus 74.55 and 0.33 for SFT only. GRPO took 146 hours; SFT took 7 hours. IR values are not stated in the text. | It is the closest recipe to our setup (code output, executable geometry reward, GRPO, k=8). Its CD-only reward failed after about 200 steps when code became invalid, so we must keep a validity term; its 146 GPU-hours show reward execution cost matters. |
| 4 | [2510.17157](https://arxiv.org/abs/2510.17157): GACO-CAD: Geometry-Augmented and Conciseness-Optimized CAD Model Generation from Single Image (October 2025) | full | Qwen2-VL-2B, one H800 for SFT and four H800 for RL: generated code is about 6.1% shorter on DeepCAD and 5.7% shorter on Fusion360. IoU, CD and IR values are not stated in the extracted text. | It shows a 2B VLM can do image-to-CAD with an IoU reward (I_high 0.8, I_low 0.4), and the scheduled reward weights and length reward are cheap tricks to copy into our reward. |
| 5 | [2605.10865](https://arxiv.org/abs/2605.10865): BenchCAD: A Comprehensive, Industry-Standard Benchmark for Programmatic CAD (May 2026) | full | Best frontier model (Gemini 3.1 Pro thinking) scores 0.318 total on image-to-code. Code QA reaches about 0.838 versus 0.587 for image QA. A Qwen3-VL-2B trained on it scores 0.7631 on in-distribution code generation, with a clear gap on held-out families. | It is a ready test set and training set for a Qwen3-VL-2B CAD run, with a reward of 0.2 x essential-op recall + 0.8 x IoU, 16 rollouts and GRPO. It also supplies a rotation-invariant voxel IoU (24 rotations) for the GPU verifier. |
| 6 | [2609.34494](https://arxiv.org/abs/2609.34494): ConCAD: Constraint-Aware Image-to-CAD Generation with Dual-Granularity Rewards (September 2026) | full | Qwen3-VL-2B, one A800 for training, DeepCAD test: IoU 0.734, CD 0.009, IR 0%, G-CSR 0.913. The SFT start gave 0.717, 0.011, 1% and 0.903. On Zero2CAD, IoU rose from 0.443 to 0.507. | It gives a 2B-model, single-GPU RL recipe that fits our case. It warns that a pure execution reward collapses to trivial cubes. Its IoU uses 8,192 Monte Carlo points, so our exact GPU IoU is a direct upgrade. |
| 7 | [2604.26779](https://arxiv.org/abs/2604.26779): Accelerating RL Post-Training Rollouts via System-Integrated Speculative Decoding (April 2026) | full | 8B reasoning RL, synchronous, 8 GB200 nodes (4 GPUs each). Generation time falls from 100.0 s to 56.6 s (1.8x) on RL-Zero, and from 133.6 s to 87.0 s (1.5x) on RL-Think. Step speedup is 1.41x and 1.35x. Generation is 65-72% of step time. Draft length 3 is best; longer drafts are slower. N-gram drafting is slower than plain decoding. Async at policy lag 1: step time 75.0 s to 60.5 s (1.24x). The 2.5x at 235B is a simulator projection, not a measurement. | Our rollouts are short polygon outputs from a 2B VLM, so decode share and acceptance length are likely small; use its step-time formula to check the gain before any work. |
| 8 | [2604.07853](https://arxiv.org/abs/2604.07853): QaRL: Rollout-Aligned Quantization-Aware RL for Fast and Stable Training under Training–Inference Mismatch (April 2026) | full | Qwen3-30B-A3B-Base, math, OpenR1-Math-46K: average math score 51.2 for QaRL, 52.1 for BF16, 45.7 for quantized rollout alone. Training is 1.3x faster than BF16 (1.4x for quantized rollout alone). Qwen3-8B-Base: 48.9 vs 51.8 (BF16) vs 43.9 (quantized rollout). GPU type is not stated in the sections read. | It names the sampler-learner mismatch that can hurt our RL; the cheap part to copy is the mismatch weight and sequence-level clipping, not the W4 rollout. Their approach tolerates mismatch; ours needs bit-exact scoring only for the reward, which is separate. |
| 9 | [2606.18967](https://arxiv.org/abs/2606.18967): EfficientRollout: System-Aware Self-Speculative Decoding for RL Rollouts (June 2026) | full | 8x A100-80GB, veRL + vLLM, GRPO, batch 128, group 8, max 8k tokens, 100 steps. Qwen2.5-7B: rollout 82.4 s to 66.3 s (-19.6%), step 132.6 s to 115.7 s (-12.7%). Qwen2.5-14B: step -10.8%. Llama3.1-8B: step -7.9%. A learned EAGLE3 drafter made Llama3.1-8B step 25.6% slower. Reward curves match the baseline. Rollout is about 70% of step time. | It gives a concrete rule to decide on one RTX 5090 whether decode is long-tailed enough for speculation; with short outputs it likely is not, and the roofline toggle is a cheap test. |
| 10 | [2608.14498](https://arxiv.org/abs/2608.14498): Rollplex: Cross-Phase GPU Spatial Sharing for Vision Language Model Post-Training (August 2026) | full | Qwen2.5-VL-32B, Video-R1 setup, 32 H800 GPUs: 1.23x-1.30x faster than serial colocation and 1.57x-2.24x faster than disaggregation. Decode uses under 19% of active SMs. Default MPS is within 2.7% of the best policy; 70/30 limits cost 1.30x-3.44x. Rewards track the baseline. Prompt tokens are 79%-98% of video samples (11%-19% for text RL). Gain is capped by the shorter of decode time and prefix time. | Our 2B VLM sees image-heavy prompts and short outputs, so prefix work is large; overlapping image encoding and prefix passes with generation on one 5090 under default MPS is the closest match to our GPU-use goal. |
| 11 | [2511.11450](https://arxiv.org/abs/2511.11450): VoxTell: Free-Text Promptable Universal 3D Medical Image Segmentation (November 2025) | full | Trained on 62K+ CT, MRI and PET volumes (158 public datasets, 1K+ classes), 64 A100 GPUs, about six days. Zero-shot mean Dice 70.85 vs 51.23 for SAT (next best) over 11 unseen datasets. Ablation: late fusion 51.68-55.11, multi-stage with deep supervision 62.55, plus batch 128 gives 69.43. Real report prompts: 50.2 Dice. ReXGroundingCT after fine-tuning: 28.2 Dice vs 13.1 for SAT. Limit: fails on knee MRI, not seen in training. No Betti or topology measure is reported. | Gives a Dice baseline and a prompt-robustness test idea for a future 3D phase; it is a segmentation network, not a polygon writer or RL reward, so it is a comparison point, not a method to copy. |
| 12 | [2605.01536](https://arxiv.org/abs/2605.01536): The Antipodal Method: Fast, Accurate, and Robust 3D Generalized Winding Numbers (May 2026) | abstract + tables | Speedup factors and queries-per-second are blank in the text (stripped), so not stated here. It reports "up to three orders of magnitude" faster than the exact hierarchical method on CPU. Setup: AMD Ryzen 9 7950X3D CPU, NVIDIA RTX 4090 GPU (OpenGL). It claims 120 FPS for 4K slices on moderate meshes. | It gives an insidedness test for open or non-manifold meshes, so a GPU voxel verifier could fill voxels from CAD meshes that are not watertight. |
| 13 | [2609.28731](https://arxiv.org/abs/2609.28731): CuACD: A Fully GPU-Resident Approximate Convex Decomposition (September 2026) | abstract + tables | "More than an order of magnitude" faster than prior search-based ACD, at matched or lower mean concavity. Exact factors are blank in the text. Setup: single RTX 4090, V-HACD, PartNet-Mobility and Objaverse subset. Zero failures on all inputs. | It shows how to keep irregular geometry kernels resident on one GPU (warp per task, device heap, no host loop), which we can copy for our CUDA checker. |
| 14 | [2607.07459](https://arxiv.org/abs/2607.07459): EmbodiedGen V2: An Agentic, Simulation-Ready 3D World Engine for Embodied AI (July 2026) | full | Single RTX 4090, 200 held-out assets: 96.5% human acceptance and 98.6% collision success at 2.6 min per asset. Without mesh fixing, runtime rises from 2.6 to 21.3 min. 83.3% of 150 task worlds usable. Cited RL study: simulation success 9.7% to 79.8%, real robot 21.7% to 75.0%. | It is a template for our Physical AI asset certification, with explicit gates and a single-attempt pass rate per gate (76.2%, 75.9%, 91.0%) that an exact checker could replace. |
| 15 | [2608.23869](https://arxiv.org/abs/2608.23869): Gen2Physics: Grounding Generated 3D Meshes in Physics via Multi-View Material Decomposition (August 2026) | full | Material mIoU 48.27% on PartNet-Material (100 objects, 1100 masks), against 15.56% for PUGS and 14.54% for NeRF2Physics. ViT only: 22.58%. Mass estimation on ABO-500 is within 0.2% ADE of the best method. Hardware not stated in the paper. | Its per-face majority vote across views and VLM-critique step is a pattern for our image-to-outline work, and its watertight sub-meshes are a target for an exact watertightness check. |
| 16 | [2610.40333](https://arxiv.org/abs/2610.40333): RFW-A: Agentic Procedural Generation of 400K Articulated Assets over 500 Categories (October 2026) | full | 402,410 assets, 500+ categories. Mean generator cost $10.62; $0.053 per asset over 200 assets, against $1.237 per asset for per-asset agentic synthesis (23-fold). Strict URDF validity 1.0000. Stable 1-second rollout 0.995 in four simulators (Core-200 subset). Real-robot success 6.2/10 to 8.0/10 with edited assets. | Its validity checks (watertight-like collision coverage, kinematic sweeps, cross-simulator rollout) and amortized-generator idea suit CAD code and asset certification with our exact checker as the gate. |
| 17 | [2604.02580](https://arxiv.org/abs/2604.02580): VoxelCodeBench: Benchmarking 3D World Modeling Through Code Generation (April 2026) | full | GPT-5 reaches 87.9% shape correctness, Claude Sonnet 4.5 and GPT-5 Mini 80.4%, Gemini 3 Pro 19.5%. Geometric construction is hardest (GPT-5 66.7%). Hardware not stated in the paper. | It shows that human rating of renders is the weak point of code-to-3D tasks, and an exact voxel verifier (b0/b1/b2, solid IoU) could replace it as a reward. |
| 18 | [2503.06800](https://arxiv.org/abs/2503.06800): VideoPhy-2: A Challenging Action-Centric Physical Commonsense Evaluation in Video Generation (March 2025) | abstract + tables | Best model Wan2.1-14B scores 32.6% joint score on the full set and 21.9% on the hard subset (human rating, 590 test prompts per model). Conservation of momentum and mass are violated about 40% of the time. The auto-rater beats Gemini-2.0-Flash-Exp by 81% (meaning match) and 236% (physical sense) in relative correlation. GPU and hardware cost: not stated in the paper. | It shows that learned video judges are weak and human labels are costly, so an exact checker as reward is a strong choice for our CAD and sim-asset work. |
| 19 | [2504.00983](https://arxiv.org/abs/2504.00983): WorldScore: A Unified Evaluation Benchmark for World Generation (April 2025) | full | 3,000 test examples, 20 models. WorldScore-Static: WonderWorld 72.69, LucidDreamer 70.40, best video model CogVideoX-I2V 62.15. Best video camera control is 40.22 (CogVideoX-T2V). Hardware cost: not stated in the paper. | It gives a template for a geometric 3D consistency check (reprojection error from SLAM depth), which is a soft cousin of our exact voxel and Betti verifier for sim-ready assets. |
| 20 | [2505.00337](https://arxiv.org/abs/2505.00337): T2VPhysBench: A First‑Principles Benchmark for Physical Consistency in Text‑to‑Video Generation (May 2025) | full | Best average score is 0.42 (Wan 2.1), worst 0.19 (SD Video), across 10 models. Conservation laws score lowest (best 0.29). Hints rarely help: Newton's third law rises from 0.25 to 0.88, but Newton's first law drops from 0.12 to 0.00. The numeric value of each rating level was stripped from the text. | It shows that prompt hints cannot fix physics errors and that humans are the only judge here, which supports training with an exact verifier reward instead. |
| 21 | [2507.15824](https://arxiv.org/abs/2507.15824): PhysVidBench: Language-Grounded Evaluation of Physical Commonsense in Text-to-Video Models (July 2025) | full | Best average accuracy across 12 models is 36.2%. Pearson correlation with humans on 120 videos is 0.45 to 0.69 for this method and 0.21 to 0.47 for the best direct VLM judge. Two video VLMs answer "yes" to nearly all questions. Open models ran on a single NVIDIA A40. The refinement gains (+11.1 and +11.9 points) are printed, but the starting values were stripped. | Direct VLM judges collapse or score poorly on physics, which backs our choice of an exact checker over a VLM judge for the reward. |
| 22 | [2601.21282](https://arxiv.org/abs/2601.21282): WorldBench: Disambiguating Physics for Diagnostic Evaluation of World Models (January 2026) | full | Foreground mIoU averages 0.4508 for Cosmos-1 (33 frames, 7B), 0.2573 for Cosmos-1 (121 frames), and 0.3488 for Cosmos-2.5 (2B). Models give realistic paths but wrong constants. Free-fall gravity estimates run from 10.831 (Cosmos-1 AR) to 0.939 (Hunyuan) against 9.81. Dataset: 469 intuitive-physics videos, 279 parameter videos. Hardware: not stated in the paper. | It scores generated video with mask IoU against ground truth and fits exact physical constants; the same per-frame IoU could run as a GPU kernel inside our checker for Physical AI tests. |
| 23 | [2609.35641](https://arxiv.org/abs/2609.35641): Verifiable Visual Rewards Transfer from Synthetic Scenes to Natural Prompts (September 2026) | full | SD3.5-M accuracy on VVRBench rises from 2.81% to 28.27% (trained on VVR-Easy) and to 46.60% (VVR-Matched). GenEval is 0.729 vs 0.688 for the GenEval2-trained model. Human win rate is 71.6% on 160 natural prompts. GPT-Image-2.5-Sunburst solves 21.39% of VVRBench-Challenge. Training hardware: not stated in the text I read. | It is the closest match to our design: an exact program verifier with partial credit and a penalty factor against reward hacking, plus a generator with set difficulty for curricula in CAD and voxel tasks. |
| 24 | [2610.10959](https://arxiv.org/abs/2610.10959): 2610.10959 — GPU-Accelerated Computation of Persistent Homology for Topological Analysis of Image Data (October 2026) [NEW] () | full | average end-to-end speedup 53.24x, maximum 198.01x, over Cubical Ripser on every benchmark. In a PARSE pulmonary artery network (128^3 volumes, 480 epochs), epoch time falls from 309.76 s to 48.56 s. | It gives a ready GPU path for b0, b1, b2 of voxel masks, but it is a persistence tool, not an exact integer Betti counter. |
| 25 | [2610.37297](https://arxiv.org/abs/2610.37297): 2610.37297 — TRUCO: A Scalable Lock-Free Algorithm for Parallel Contour Extraction (October 2026) [NEW] () | full | 2.3 to 17.8 times faster than the fastest sequential method with the same output, on multi-core CPUs. Output matched Suzuki-Abe in 5.2 million executions. | It shows how to prove a parallel contour tracer bit-identical to a reference; the ownership rule can guide a CUDA hole and outline checker. It targets CPUs, not GPUs. |
| 26 | [2607.08891](https://arxiv.org/abs/2607.08891): 2607.08891 — Ortho2CAD: 3D CAD generation from orthographic drawings using vision language models (July 2026) [NEW] () | full | RL model reaches 100 % valid code on four test sets (100 examples each). On F360Recon the IoU is 0.5601, against 0.3697 for the SFT model. | It is the closest recipe to our CAD plan: execute code, compute exact solid IoU, give 0 on failure. Its IoU needs a CAD kernel on CPU. |
| 27 | [2606.13368](https://arxiv.org/abs/2606.13368): 2606.13368 — IterCAD: An Iterative Multimodal Agent for Visually-Grounded CAD Generation and Editing (June 2026) [NEW] () | full | IterCAD-Draw benchmark, invalid ratio 64.60 % (Qwen3.5-4B base) to 0.30 % (full), mean Chamfer 10.48 to 5.09. Average 2.97 turns. | The prefix mask fixes credit assignment when an early valid geometry is followed by errors; reuse it for multi-turn CAD repair. |
| 28 | [2610.03453](https://arxiv.org/abs/2610.03453): 2610.03453 — I2CD: Direct Image-to-Convex Decomposition for Simulation-Ready Collision Geometry (October 2026) [NEW] () | full | On 227 held-out objects it has the highest volumetric IoU of eight reconstruct-then-decompose pipelines, 6 to 37 times faster. MuJoCo, PyBullet and PhysX silently replace a raw mesh by its convex hull, the wrong collision shape for about 60 % of objects; I2CD geometry is preserved on at least 98.4 % of objects in four engines. | It shows that "it loads" hides silent collision-shape swaps; a sim-ready certificate must test the shape the engine really uses. |
| 29 | [2607.09207](https://arxiv.org/abs/2607.09207): 2607.09207 — Bidirectional Resource Scheduling for Disaggregated and Asynchronous RL Post-Training (July 2026) [NEW] () | full | Up to 1.94x throughput over veRL, AReaL and ROLL on two 32-GPU testbeds (A6000, H100), with unchanged convergence. | It confirms that idle time is the main loss, but it needs many GPUs; on one RTX 5090 only the idea of filling bubbles transfers. |

---

## Summary

The model is as good as the doctors near the image centre. It misses far nodules because RL used only centred
images. About half of each step's learning signal is lost to a number mismatch, and about 40 % of each step leaves
the GPU idle. The fixes are known, cheap, and mostly built.
