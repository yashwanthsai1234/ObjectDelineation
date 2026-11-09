# The implementation architecture, explained

Date: 2026-11-09. Status: plan only. Nobody changed code or rented a GPU for this plan.

This text explains the full plan in simple English. It uses diagrams made of text. The full plan has the exact
tasks, files and tests. This text tells you how the parts fit together and why.

## Words in this text

We use one word for each thing. The table gives the word once. The text does not use a second word for it.

| Word | Meaning |
|---|---|
| CPU | The normal chip of a computer. |
| GPU | The fast chip that does most of the work in this project. |
| CUDA | The language for GPU programs from NVIDIA. |
| model | The program that we train. It looks at an image and writes an outline. Today it has 2 billion numbers. |
| CT slice | One flat cross-section of the body from a CT scan. |
| image | The square piece of a CT slice that the model sees: 512 x 512 pixels. |
| centred image | An image with the nodule at the middle. |
| shifted image | An image with the nodule away from the middle. |
| resized image | An image cut larger or smaller, so the nodule looks smaller or larger. |
| candidate image | A possible image. The plan stores only its position and size, not its pixels. |
| image numbers | The numbers that the model builds from an image before it can read it. |
| nodule | The small round lump in the lung that the model must outline. |
| mask | A black-and-white map of the image. The white pixels are the nodule. |
| outline | A closed line around the nodule. |
| answer | One outline from the model, written as text. |
| legal answer | An answer with none of the 13 kinds of defect. |
| cut answer | An answer that reaches the outline limit and stops before it is complete. |
| runaway answer | An answer that does not stop by itself. It runs until the outline limit. |
| doctor | One of up to 4 doctors who drew an outline on the same image. |
| checker | The exact GPU program that tests each answer and gives it a score. |
| score | How well an answer matches the doctors. The main part is overlap: 0 is none, 1 is a full match. |
| ceiling | How well the doctors match each other: 0.732. It is the best fair score. |
| miss | An image where the mean overlap of the 8 answers is below 0.1. The model looked in the wrong place. |
| offset | How far the nodule is from the middle of the image, as a part of the image side. 0.5 is at the edge. |
| bin | A set of images with one range of nodule size and one range of offset. |
| step | One training round: 32 images, 8 answers for each image, 256 answers in total. |
| group | The 8 answers for one image in one step. |
| reward | The number that the teacher learns from. The plan builds it from the score. |
| lesson size | The reward of one answer minus the mean reward of its group. |
| writer | The GPU part that writes the answers. Its name is vLLM. |
| teacher | The GPU part that changes the model after the checker gives the scores. Its name is TRL. |
| adapter | The small trained part of the model. Training changes only the adapter. |
| trial training | Training where the model writes answers, gets scores and learns from them. |
| copy training | Training where the model learns to copy the doctor outlines. |
| token | One small piece of text. A mean answer has 131 tokens. |
| outline limit | The largest number of tokens that one answer can have. Today it is 1,536. |
| split | One of three sets of images: training, validation and test. No patient is in two splits. |
| validation images | Images that are not used in training. We use them to choose between models. |
| test images | Images that are not used in training. We use them once, for the final score. |
| fingerprint | A short code made from a file. If the file changes, the code changes. |
| store | A private folder on Hugging Face that holds files. |
| Mac | The computer of the owner. Work on the Mac costs nothing. |
| box | A GPU computer that we rent from Vast AI: one RTX 5090 with 32 GB of memory. It costs money. |
| rental | One session on one box, from start to stop. |
| gate | A test with a number. Work goes on only when the gate passes. |
| hard gate | A gate that stops the work if it fails. |
| layer | One part of the code. A layer can use only the layers below it. |
| kernel | A small program that we write for the GPU. |
| archive | A folder for old code. Nothing uses it. |

## 1. The short story

The model reaches 0.645 overlap on 2,200 test images. The ceiling is 0.732. On the images where the model finds
the nodule, it is almost as good as a doctor. It scores 0.722. The doctors score 0.740 on these images.

The gap comes from misses. 10.8 % of the test images are misses. Each miss scores almost 0.

The misses have one main cause. In trial training, every image was a centred image. In the test, the nodule is
often far from the middle. The model never practised that case.

The plan fixes this in seven phases. It fixes the training data first. Then it fixes how the teacher learns. Then
it removes the GPU waits. Then it moves all outline and mask work to the GPU. Last, it puts the code in order.

```
+-------------+   +-------------+   +-------------+   +-------------+   +-------------+   +--------+   +---------+
| A: measure  |-->| B: quality  |-->| C: speed    |-->| D: CUDA     |-->| E: code     |-->| DOC:   |-->| CLEAN:  |
| fairly      |   | DATA FIRST, |   | remove GPU  |   | outline and |   | layers and  |   | docs   |   | archive |
|             |   | writer and  |   | waits       |   | mask work   |   | one run     |   | site   |   | old and |
|             |   | teacher     |   |             |   | on the GPU  |   | file        |   |        |   | unused  |
|             |   | agree, new  |   |             |   |             |   |             |   |        |   | code    |
|             |   | reward      |   |             |   |             |   |             |   |        |   |         |
+-------------+   +-------------+   +-------------+   +-------------+   +-------------+   +--------+   +---------+
       Mac           Mac, then        Mac, then         Mac, then          Mac             Mac          Mac
                     rental 1         a later rental    a later rental
```

A to E change the work. DOC writes the docs. CLEAN moves old code to the archive and removes unused code.

After phase E, the work goes to 3D voxels. This plan does not cover 3D. It keeps the parts that 3D needs.

## 2. The whole system

This diagram shows all the parts and where they run.

```
 MAC (no cost)                                  STORE (private Hugging Face)        BOX (rented RTX 5090)
+------------------------------------+         +---------------------------+       +-------------------------------+
| CT slices and doctor outlines      |         | training images and their |       | read the images; compare the  |
|   |                                |         | list, with one fingerprint|------>| fingerprint; stop if it is    |
|   v                                |         | (one fixed version)       |       | different                     |
| make candidate images              |         +---------------------------+       |   |                           |
| (centred, shifted, resized)        |                      ^                      |   v                           |
|   |                                |                      |                      | writer: 8 answers per image   |
|   v                                |                      |                      |   |                           |
| choose 8,320 images by a plan      |                      |                      |   v                           |
|   |                                |                      |                      | checker: legal? overlap,      |
|   v                                |                      |                      | edge match, area, holes       |
| quality gates (leakage, coverage,  |----------------------+                      |   |                           |
| outline inside image, tokens)      |   upload only if all hard gates pass        |   v                           |
+------------------------------------+                                             | teacher: change the adapter   |
                                                                                   |   |                           |
+------------------------------------+         +---------------------------+       |   v                           |
| read the results; decide by the    |<--------| adapter, scores, logs,    |<------| send the new adapter to the   |
| gates; the owner approves money    |         | saved answers             |       | writer; next step             |
+------------------------------------+         +---------------------------+       +-------------------------------+
```

The Mac does all the data work. The box does only training and scoring. Only the box costs money.

The fingerprint protects the data. If one training image changes after the checks, the fingerprint changes. Then
the box stops before it trains.

## 3. The layers of the code

Today the training code sits inside one experiment script. Many tests depend on experiment scripts. The plan puts
the code in eight layers. A layer can use only the layers below it. A tool (import-linter) makes this rule a test.

```
+-----------------------------------------------------------------------------------+
| experiment scripts (short)  +  run files (one checked file for each run)          |
+-----------------------------------------------------------------------------------+
| infra: box rental rules, money record, upload, start-up checks                    |
+----------------------------------------+------------------------------------------+
| rl: trial training (teacher set-up,    | sft: copy training                       |
|     reward call, adapter push)         |                                          |
|     rl and sft never use each other    |                                          |
+----------------------------------------+------------------------------------------+
| policy: the model, the prompt, the image numbers on the GPU, time marks           |
+-----------------------------------------------------------------------------------+
| data: readers, splits, candidate images, coverage, quality gates                  |
+-----------------------------------------------------------------------------------+
| verify: reward, scores, statistics, gates, ceiling                                |
|         uses NO training library                                                  |
+-----------------------------------------------------------------------------------+
| geom: the CUDA checker, the answer text reader, the outline builder               |
+-----------------------------------------------------------------------------------+
| core: the form of a run file, and the fingerprints of the code for each run       |
+-----------------------------------------------------------------------------------+
        upper layers use lower layers; a lower layer never uses an upper layer
```

"Checked file" means that a tool (pydantic) finds a wrong or missing setting before the run starts.

Why this order:

1. The checker and the verify layer are at the bottom. They need no training library. So the same checker can
   score 2D outlines today and 3D voxels later.
2. Our code uses some internal names of the teacher library. One file holds all of them. A start-up check tests
   them. A new library version then fails on the Mac in seconds. It does not fail on the box after money is spent.
3. Old code goes to the archive. Nobody sees two versions of the same thing.

The three most important joins in the code:

| Join | What goes in and out | What can change behind it |
|---|---|---|
| outline builder | a mask in, closed outlines out | the GPU builder (kernel K1) replaces the old CPU builder |
| answer maker (for evaluation only) | images in, answers out | the writer, or a plain call to the model |
| teacher join | one file for every internal name of the teacher | the version of the teacher library |

## 4. The data path (phase B1, the highest priority)

### 4.1 What is wrong with the data today

Today trial training uses 8,208 centred images. The code also has a first fix: 24,624 images. It adds one shifted
image and one resized image for each centred image. The plan measured this first fix. It is better, but it still
does not match the validation images.

```
                                   median offset    share with offset above 0.2
 trial training today                  0.018             almost 0 %
 the first fix (24,624 images)         0.181             47.2 %
 validation images                     0.257             68.4 %
```

The cause: the first fix cuts only one shifted image and one resized image for each nodule. Some bins stay too
small. These are mainly small nodules far from the middle.

### 4.2 The new data path

The new path makes many candidate images. Then it chooses images by a plan. Only the chosen images get pixels.

```
+--------------------+   +-----------------------+   +------------------------+   +----------------------+
| store of CT slices |-->| candidate pool: many  |-->| schedule: 260 steps x  |-->| cut only the chosen  |
| and doctor         |   | shifted and resized   |   | 32 = 8,320 images, in  |   | images; build their  |
| outlines (training |   | candidate images for  |   | blocks of 20 steps;    |   | doctor outlines      |
| patients only)     |   | each nodule; position |   | each block fills each  |   |                      |
+--------------------+   | and size only         |   | bin by its target share|   +----------+-----------+
                         +-----------------------+   +------------------------+              |
                                                                                             v
                                                     +------------------------+   +----------------------+
                                                     | upload to the private  |<--| quality gates:       |
                                                     | store; the run file    |   | leakage, coverage,   |
                                                     | holds the version and  |   | outline inside the   |
                                                     | the fingerprint        |   | image, token count   |
                                                     +------------------------+   | (21 gates)           |
                                                                                  +----------------------+
```

The target shares come from the training patients only. No validation image or test image shapes the training
data. A block of 20 steps already matches the target shares. So a run that stops early still trained on the right
mix.

### 4.3 The quality gates

There are 21 gates of 6 kinds. A failed hard gate writes no image list. Then the teacher cannot load the images.
So bad data cannot reach a box.

| Gate kind | What it checks | Number to pass |
|---|---|---|
| leakage | no patient, CT series or pixel copy is in two splits | 0 overlaps |
| inside | every doctor outline is fully inside its image | 100 % |
| coverage | each bin has enough images | the smallest bin has at least 0.8 of its target share |
| coverage | the whole mix of bins is near the target | it differs by at most 0.10 (a measure called TV) |
| coverage | the spread of offsets is near the target | it differs by at most 0.06 (a measure called KS) |
| tokens | no doctor outline is longer than the outline limit | 0 above the outline limit |

The plan already checked today's data. The patient, series and pixel overlap across splits is 0. Every doctor
outline is inside every image. Two near-copy pairs cross splits. A gate lists them for review.

## 5. One training step after all phases

### 5.1 Flow diagram (what happens in order)

```
+------------------+   +------------------+   +------------------+   +------------------+   +------------------+
| 1. 32 images     |-->| 2. the teacher   |-->| 3. writer: 8     |-->| 4. checker:      |-->| 5. reward and    |
| read on the CPU; |   | builds the image |   | answers for each |   | legal? overlap,  |   | lesson size for  |
| 25.2 MB copied   |   | numbers once, on |   | image; each      |   | edge match, area,|   | each answer      |
| to the GPU       |   | the GPU (not     |   | answer has its   |   | holes; only the  |   |                  |
|                  |   | 1.61 GB on the   |   | own outline limit|   | scores leave the |   |                  |
|                  |   | CPU, twice)      |   |                  |   | GPU              |   |                  |
+------------------+   +------------------+   +------------------+   +------------------+   +--------+---------+
                                                                                                    |
+------------------+   +------------------+   +------------------+                                  |
| 8. next step     |<--| 7. send the new  |<--| 6. teacher:      |<---------------------------------+
|                  |   | adapter to the   |   | change the       |
|                  |   | writer in one    |   | adapter; one     |
|                  |   | call             |   | read from the GPU|
+------------------+   +------------------+   +------------------+
```

### 5.2 Interaction diagram (who talks to whom, in order)

```
 CPU                    writer (GPU)          checker (GPU)          teacher (GPU)
  |                          |                      |                      |
  |--32 images (25.2 MB)------------------------------------------------->|  builds the image numbers once
  |--32 prompts with their images-->|               |                      |
  |                          |  writes 8 answers for each image            |
  |<--256 answers (text)-----|                      |                      |
  |--reads the text into outlines------------------>|                      |
  |                          |                      |--256 scores--------->|
  |                          |                      |                      |--changes the adapter
  |                          |<------------- new adapter, one call --------|
  |<-------------------------------------------- one log read per step ----|
```

### 5.3 What runs where

The rule is that all work runs on the GPU. The CPU keeps only the work that it must do.

| Work | Today | After the plan | Why the CPU keeps it |
|---|---|---|---|
| image file decode | CPU | CPU | no GPU decoder exists for this file type (PNG); 32 files for each step |
| image numbers | CPU, twice, 1.61 GB | GPU, once, from 25.2 MB | - |
| position numbers (where each part of the image is) | CPU loop | GPU, one formula | - |
| image part of the model | GPU, many times | GPU, once for each image | - |
| answer text to outline | CPU | CPU | the writer gives text to the CPU; a GPU reader saves about 0 s |
| legal check and scores | GPU | GPU | - |
| outline builder for the doctor outlines | CPU | GPU (kernel K1) | - |
| statistics for evaluation | CPU | GPU (kernel K2) | - |
| adapter push | about 700 calls to the writer | 1 call | - |
| small reads from the GPU | about 1,250 in each step | one read in each step, for the logs | - |

## 6. How the teacher learns (phases B2 and B3)

### 6.1 The writer and the teacher must agree

The writer and the teacher each give every token a chance value. These values must be equal. They are not.

The plan found one measured cause. After each step, the teacher adds the adapter change into the writer's copy of
the model. The writer keeps its numbers in a short 16-bit format, called bf16. Most small changes get lost in the
rounding.

Example: a stored number can be 1.000 or 1.008, but nothing between them. A change of 0.002 rounds back to 1.000.
The change is lost.

```
 today                                                   option F (expected)
+------------------+                                    +------------------+
| teacher: model + |                                    | teacher: model + |
| adapter, kept    |                                    | adapter, kept    |
| apart            |                                    | apart            |
+--------+---------+                                    +--------+---------+
         | adds the change into the model,                       | sends the adapter only,
         | then rounds to 16 bits                                | kept apart
         v                                                       v
+------------------+                                    +------------------+
| writer: 46-69 %  |                                    | writer: the full |
| of the small     |                                    | change, as the   |
| changes are lost |                                    | teacher has it   |
+------------------+                                    +------------------+
```

At step 100, the rounding error is 0.53 of the size of the adapter change in bf16. In the other 16-bit format
(FP16), it is 0.10. Option F removes the rounding. In option F, the writer keeps the adapter apart, as the teacher
does. A rule in phase B picks the option from a measured number. With today's numbers, the rule picks option F.

### 6.2 The new reward

Today, if all 8 answers of a group score near 0, the group teaches nothing about where the nodule is. The new
reward adds a small score for "near the nodule". So such a group can still push the model toward the nodule.

```
 reward      = legal x (overlap + 0.1 x nearness) + a small term for the right number of holes
 lesson size = reward - mean reward of the group        (we do not divide by the spread)
```

The new reward also removes a push toward cut answers. Today a cut answer counts in the group mean, but it gives
no lesson. This makes cut answers more likely over time. In the full earlier run, the share of steps with an answer
at the outline limit rose from 13.8 % to 31.5 %.

## 7. The checker and the three kernel drafts

```
+---------------------------------------------------------------------------------+
| checker (one source file for the CPU and the GPU)                               |
|                                                                                 |
|  +--------------------+  +--------------------+  +--------------------------+   |
|  | legal check:       |  | scores: overlap    |  | outline builder (K1):    |   |
|  | 13 defect kinds,   |  | with each doctor,  |  | mask in, outlines out;   |   |
|  | one group of GPU   |  | edge match, area,  |  | looks only at pixels on  |   |
|  | threads for each   |  | holes              |  | the edge of the nodule   |   |
|  | outline (K2)       |  |                    |  | (today: about 46 million |   |
|  |                    |  |                    |  | pixel visits per mask)   |   |
|  +--------------------+  +--------------------+  +--------------------------+   |
|                                                                                 |
|  +--------------------------------------------+  +--------------------------+   |
|  | evaluation statistics (K2): distance map,  |  | image numbers and        |   |
|  | scores against the set of all doctors      |  | position numbers (K3),   |   |
|  |                                            |  | used only if the library |   |
|  |                                            |  | way fails on the box     |   |
|  +--------------------------------------------+  +--------------------------+   |
+---------------------------------------------------------------------------------+
```

| Kernel | What it does | What is proved | What is not proved |
|---|---|---|---|
| K1 | builds outlines from masks on the GPU | its CPU build gives the same outlines as the old builder on all 37,093 real LIDC masks | no GPU run yet |
| K2 (both boxes) | legal totals, distance map, scores against all doctors | the same results on the Mac and on Linux; it compiles for the RTX 5090 | no GPU run yet |
| K3 | image numbers, position numbers | 0 of 64.5 million values differ from the library, on the CPU | no GPU run yet |

A kernel stays only if it is faster than the library in a measured test. If it is not, we drop it.

## 8. The model choice

The plan keeps the 2B model for the run with the new training images. It changes one thing at a time.

Reasons:

1. The cause is data. On the images where the model finds the nodule, it is near doctor level.
2. In a 2026 paper (EGM, arXiv 2601.13633), training cut the gap between the 2B and 4B models. On a standard
   test of pointing at objects (RefCOCO), the gap went from 3.6 to 1.4 points.
3. A 4B run costs about $2.9 by estimate. It does not fit the $3 budget beside the 2B run.

The 4B model (Qwen3-VL-4B) is the next model that we try only if two things are true on the validation images:

- the overlap is still below about 0.68;
- the misses are still above about 7 %.

The owner sets these numbers before any money is spent. The 4B model has 4.44 billion numbers. The owner wrote "up
to 4B", so the owner must also decide if 4.44 billion is allowed.

## 9. Tools: what we add and what we do not add

| Add now | Why |
|---|---|
| checked run files (pydantic) | a wrong setting fails on the Mac, not on the paid box |
| a layer test (import-linter) | the layer rule becomes a test |
| time marks in each part of the step (NVTX) | we measure each wait before we fix it |

We do not add these now:

- Redis (a fast shared memory) or a work queue;
- Ray (a tool to spread work over many computers);
- DVC (a tool to track data versions);
- a database or a tracking server.

Each is a system that someone must run, and no measured need exists. Each has a written condition. For example, we
add a work queue only when more than one box writes answers.

## 10. The order of work and the money

```
 Mac ($0)                                                        box ($)
+--------------------------------------------+
| docs under version control                 |
| A: measurement code; true outline lengths  |
| B1: build and check the training images    |
| CLEAN, part 1                              |
+---------------------+----------------------+
                      |
      all 21 data gates pass? --no--> fix on the Mac; no money spent
                      | yes
                      v
          owner approves rental 1 -----------------------------> rental 1: training with the
                      |                                         new images, then evaluation
                      |<--------------------------------------- scores, saved answers, adapter
                      v
          Mac: the gate result decides the next step
                      |
                      v
          Mac: code for C, B2, B3 and D -----------------------> rental 2: speed checks and
                      |                                         kernel checks (needs a budget raise)
                      |<---------------------------------------
                      v
          owner budget decision -------------------------------> rental 3: new reward run
                      |                                         (needs a budget raise)
                      v
          E: move the code into layers -> DOC -> CLEAN, parts 2 and 3
```

### 10.1 The money for rental 1

Two money words:

- **code stop**: the spending limit in the code. Today it is $17.50. The code refuses a rental that can pass it.
- **money guard**: a rule inside rental 1. It stops the training when the rental has cost $2.05.

The $3 budget pays for one rental: rental 1. The money record already shows $16.90 of earlier spend. So today the
code stop refuses any rental. The owner must approve three things in one message:

1. Raise the code stop to $16.90 + $3.00 = $19.90.
2. Approve $2.40 for rental 1, with its item list.
3. Approve the table that decides the next step from the result of rental 1.

| # | Item of rental 1 | Minutes | Cost |
|---|---|---|---|
| 1 | start the box, copy the code and the images, build the CUDA checker, run the start-up checks (includes a fixed $0.10 data-transfer fee) | 53 | $0.70 |
| 2 | train for 260 steps on the new training images | 119 | $1.35 |
| 3 | score the new model on two sets of validation images: the shifted set and the edge set (nodules near the image edge) | 10 | $0.12 |
| 4 | score the old model on the edge set (for comparison) | 7 | $0.08 |
| 5 | score the old model on the shifted set, and save all its answers | 4 | $0.04 |
| | total | 193 | $2.29 |

These numbers use 27.5 s for each step. That is the mean step time of the full 3-pass run of 2026-11-08.

How the three money numbers fit together:

- $2.29 is the planned cost.
- $2.40 is the time limit of the box, as money. The box stops at this limit. The money guard at $2.05 makes sure
  that items 3 to 5 finish before it.
- $2.49 is the worst case. It is $2.40 plus $0.11 for the last 10 minutes, when the box sends back the results
  and shuts down.

With a reserve of $0.25, the worst case is $2.49 + $0.25 = $2.74. That is below $3.00.

### 10.2 After rental 1

After rental 1, $3.00 - $2.29 - $0.25 = $0.46 is left. If the steps are faster (23.65 s), the cost is $2.10, and
$0.65 is left. That money does not pay for any useful part of rental 2. The speed checks of rental 2 alone cost
$1.21. So rentals 2 and 3 each wait for the owner to raise the budget. Until then, the money that is left stays
unspent.

Rental 1 runs today's code with the new training images. It has no new speed code and no new kernel code. So its
result shows the effect of the new data only.

## 11. Worked examples

**Worked example 1: why the data comes first.** The plan measured the misses by offset on the 2,200 test images.

```
 offset          images   overlap   misses
 below 0.15         375     0.751      -
 0.15 to 0.3      1,018       -        -
 above 0.3          807     0.552   22.2 %  = 179 images
```

So 807 of 2,200 images (36.7 %) have the nodule far from the middle. On those images, 22.2 % are misses: 179
images. In trial training, the median offset was 0.018. That is 0.018 x 512 = 9 pixels from the middle. In the
test, the median is 0.262. That is 134 pixels from the middle. The model did not practise this case. More GPU
speed or a larger model does not fix it. Better training images do.

**Worked example 2: the image numbers.** The model cannot read an image file. It needs image numbers. One image
becomes 1,024 rows of 1,536 numbers, at 4 bytes each: 1,024 x 1,536 x 4 = 6.3 MB. (This 1,536 is a size inside
the model. It is not the outline limit.) The teacher needs a copy for each of the 256 answers: 256 x 6.3 MB =
1.61 GB. Today the CPU builds this twice in each step. After phase C, the CPU sends only the 32 plain images:
32 x 512 x 512 x 3 bytes = 25.2 MB. The GPU builds the numbers once.

**Worked example 3: the outline limit.** You asked if a higher outline limit of 1,800 removes the wait. It makes
the wait longer, not shorter. The writer must wait for the longest answer in the step. The measured rule is:
writer time = 5.61 s + 4.3 ms for each token of the longest answer. A runaway answer runs to the outline limit.

```
 all other answers end before 600 tokens; one runaway answer runs to the outline limit

 outline limit 1,536:   4.3 ms x (1,536 - 600) = 4.0 s of extra wait
 outline limit 1,800:   4.3 ms x (1,800 - 600) = 5.2 s of extra wait     (longer, not shorter)
 own outline limit:     the median doctor outline has 101 tokens. That image gets
                        1.5 x 101 + 64 = 216 tokens. A runaway answer on it stops at 216.
                        That is below 600, so it adds 0 s of extra wait.
```

A higher outline limit lets a runaway answer run longer. So the wait grows. The plan does the opposite. It counted
the true outline lengths. The longest doctor outline in the training images has 1,182 tokens. Only 13 of 23,469
doctor outlines have 800 tokens or more. Phase C can give each image its own outline limit: 1.5 times its longest
doctor outline, plus 64 tokens, but never more than 1,536. For the image with the 1,182-token outline, that gives
1,837, so its outline limit stays at 1,536. A runaway answer on a normal image then stops early. This is an owner
decision. A gate on the validation answers must pass first.

## Summary

The model misses far nodules because it practised only on centred images. So the plan first builds and checks
better training images. Then it makes the writer and the teacher agree, and moves all image and outline work to the
GPU in a layered code base. Only the box costs money, and every rental waits for the gates and for the owner.
