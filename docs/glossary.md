# Glossary

Each term has one definition and one word.
The other documents never use a synonym for a term below.
Terms are in groups. Inside a group they follow the order in which you meet them.
Identifiers in `code font` are the names in the source.

## Data and views

| Term | Definition |
|---|---|
| **image** | The picture the model sees. A CT slice (LIDC-IDRI) or a fundus photograph (RIGA). In the store, LIDC images are int16 Hounsfield units and RIGA images are uint8 RGB. |
| **series** | One LIDC CT scan: a stack of slices with one series uid. A patient has one series, or two for 8 patients. |
| **rater** | One expert who drew an outline of the object. LIDC has up to four per nodule. RIGA has six per photograph. The tables of an item keep one mask per rater. |
| **item** | One image together with the exact mask of every rater (`objdel.data.Item`). One LIDC item is one slice of one nodule that at least two raters outlined. |
| **structure** | A named thing the model outlines. `nodule` for LIDC. `disc` and `cup` for RIGA (`SCHEMAS`). |
| **object** | The outline of one structure inside one answer. Also called by its name, for example the `nodule` object. |
| **frame** | The square region of the image that the model sees and that the scorer uses. Its side S divides 1000 and is at most 500 pixels (`objdel.targets.frame`). |
| **view** | One window of an item. `centred` (the policy frame), `shifted` (the window moves), `rescaled` (the window side changes), `rotated` (a quarter turn, LIDC exam only). |
| **window** | The frame position and side, `(offset, side)`. For CT, "window" also names the Hounsfield range, written "CT window". |
| **grid** | The 0 to 1000 integer coordinate system of every outline. Grid point (x, y) is pixel point (x W / 1000, y H / 1000). |
| **view directory** | A directory that holds the rows of one training or exam set (`data/views/<name>`, `data/eval_views/<name>`). The word "view" alone always means one window kind. |
| **split** | `train`, `val` or `test`. Splits are group-aware: one patient, or one set of near-duplicate photographs, stays in one split. |
| **store** | The Parquet files in `data/store` that hold exact images and masks (`objdel.data.store`). |
| **prepared set** | The files in `data/prepared`: frame images as PNG, certified targets and expert tables, per item. |
| **disagreement item** | An item whose raters do not give one topology answer. RL repeats these items 4 times (`DISAGREEMENT_OVERSAMPLE`). |
| **RL-eligible** | An item with at least two raters. Single-rater items do not enter RL (`rl_eligible`). |

## Geometry

| Term | Definition |
|---|---|
| **outline** | The polygon description of an object: parts, each with one outer ring and zero or more holes. The model writes outlines as JSON text. |
| **polygon** | The checker's word for an outline in packed form: rings of integer points. It appears in function names (`check_polygons`) and in "polygon JSON". |
| **part** | One outer ring with its holes (`Part`). An object has zero or more parts. |
| **ring** | A closed list of at least three integer points. An outer ring runs clockwise on screen. A hole runs counter-clockwise. |
| **hole** | A ring inside an outer ring that cuts out area. The checker tests each hole against its outer ring. |
| **legal** | An answer is legal when it parses and the checker finds no defect in any ring of any object. |
| **defect** | One of 13 named faults that make a ring illegal (`DEFECTS`): for example `self_crossing`, `wrong_winding`, `out_of_range`. |
| **even-odd mask** | The pixel mask of an outline. A pixel is inside when its centre is inside under the even-odd rule with a half-open edge rule. |
| **b0, b1, b2** | Betti numbers. b0 counts parts. b1 counts holes in 2-D, or tunnels in 3-D. b2 counts cavities in a 3-D volume. |
| **chi8** | The Euler number of a mask with 8-connected foreground and 4-connected background. chi8 = b0 - b1. |
| **bit quad** | A 2x2 pixel pattern. Counts of the 16 patterns give area, perimeter, corners and chi8 (Gray 1971). |
| **certificate** | The check that three independent measures of one mask agree: bit quads, labelled topology and traced rings (`checker.certify`). |
| **size filter** | The rule that removes parts and fills holes smaller than `min_area` = max(5, ceil(0.5 % of the union area)) pixels. |
| **signature** | The holes per part after the size filter, sorted high to low. `(0,)` is one part without a hole. `(0, 0)` is two parts. `(1,)` is one part with one hole. `()` is empty. |
| **topology answer** | The signature of an answer (for several objects, the tuple of signatures). The reward compares it with the signatures of the raters. |

## Reward and scoring

| Term | Definition |
|---|---|
| **crop (engine)** | Inside the checker only: the bounding box of one mask plus a one-pixel border, used to label parts and holes. |
| **checker** | The whole CUDA engine `gpucheck`. It validates outlines, draws masks, traces rings and runs the scorer. One source file builds for the CPU and for nvcc. |
| **scorer** | The part of the checker that turns an outline into the reward terms, against the expert table of the item. |
| **expert table** | The per-object data of the raters: masks, union, core, vote, distance map, tau, T5, signatures. Built on the CPU once (`ExpertCache`), then kept on the device. |
| **union, core, vote** | The pixels drawn by at least one rater (union), by all raters (core), and by more than half of them (vote). |
| **band** | The region between core and union. Pixels inside the band are neither right nor wrong. Also called the egg-yolk band. |
| **IoU_tc** | The overlap term. The highest IoU against the raters who drew the same signature as the answer. If no rater did, the IoU against the vote. |
| **F5** | The boundary term, `1 / (1 + (2 D / P) / T5^2)`. D is the depth of error beyond the band. P is the union contour length. A new rater scores 0.95 on average. |
| **g_area** | The area guard, `min(1, A_max / A) * min(1, A / A_min)`. A is the answer area. A_min and A_max are the smallest and largest rater areas. |
| **T5, tau** | Tolerances in pixels. tau is the 95th percentile of rater-pair contour distances, at least 1 pixel. T5 is the self-calibrated F5 tolerance. |
| **voter, abstain** | A rater votes on topology when the size-filtered mask is not empty. A rater with an empty signature abstains. |
| **topo_any** | 1 when some voting rater drew the topology answer of the answer, else 0. |
| **topo_freq** | The share of voting raters who drew that topology answer. |
| **r_topo** | The topology reward. It matches the rater frequency of the answer with a leave-one-out estimate of the model's own probability. |
| **term** | One number of the reward: `iou_tc`, `f5`, `legal`, `g_area`, `topo_any`, `topo_freq`, `r_topo`. |
| **gate** | A rule that sets terms to zero. The legality gate zeroes every term of an illegal answer. The validity gate decides if `r_topo` counts. |
| **reward v4** | The frozen reward definition (`REWARD_V4`). A change of any value is a new reward. |
| **precision, recall, Dice** | Pixel scores of the answer against the strict-majority outline of the raters. Used in evaluation only. |

## Training

| Term | Definition |
|---|---|
| **answer** | The text the model writes for one image. In code, a completion. |
| **group** | The G = 8 answers for one image in one step. Advantages are computed inside a group. |
| **step** | One optimizer step of RL: 32 images, 8 answers each, 256 answers. |
| **advantage** | The number that scales the gradient of one answer. In reward v4 it is the weighted sum of group z-scores plus the centred topology reward. |
| **SFT** | Supervised fine-tuning. The model learns to write certified targets. |
| **target** | The certified, canonical outline text that SFT teaches (`build_target`). One per voting rater. |
| **GRPO** | Group relative policy optimisation. The RL method. TRL implements it. |
| **Dr. GRPO** | The loss variant with a constant length normaliser and no standard-deviation scaling (`loss_type="dr_grpo"`). |
| **LoRA** | Low-rank adapters. Qwen3-VL-2B trains rank-16 adapters on the language model only. The vision tower stays frozen. |
| **FastGRPO** | The subclass of TRL's `GRPOTrainer` in `objdel.train.fast`. It trims forwards, caches vision features and packs shared prompts. |
| **colocated vLLM** | vLLM runs in the same process and on the same GPU as the trainer. Sleep mode frees its memory during the update. |
| **shared prefix** | The packed forward that runs the prompt of a group once for its 8 answers (`--shared-prefix`). |
| **bucket** | A rounded width of a forward chunk. FastGRPO rounds to 64 columns unless the text layers are compiled. |
| **WSD** | The learning-rate schedule warmup, stable, decay. SFT keeps decay branches at 0.5, 1, 2 and 3 epochs. |
| **branch** | A decay continuation of the stable SFT run, saved as a model folder (`branches/e<epochs>`). |
| **G-probe** | The 20-step RL test that asks whether RL improves a checkpoint. |
| **G-vision** | The check that the model still looks at the image (`objdel.eval.vision`). |
| **hack audit** | The statistics that detect hedging, inflation, fragmentation and duplicate answers (`objdel.eval.audit`). |

## Operations

| Term | Definition |
|---|---|
| **box** | A rented GPU machine on Vast AI. |
| **offer** | A box listing that passes the filter in `policy.terms`. |
| **run** | One launch of a chain on a box. It has a run id and a folder `runs/<id>`. |
| **chain** | The list of steps in a config YAML that `run.sh` executes in order. |
| **bundle** | One deterministic tar with the prepared set, views and eval views. It is the only data that goes to a box. |
| **preflight** | The set of checks on the Mac that must pass before a box is rented. It writes a stamp. |
| **stamp** | The file `runs/preflight/<config>.json`. Launch refuses without a fresh passing stamp. |
| **cap** | The maximum hours of a run. The watchdog destroys the box when the cap passes. |
| **ledger** | `experiments/phase2/vast/ledger.csv`. It records every offer, launch, sync and invoice. |
| **oracle** | A slow, simple, exact reference implementation used only in tests (`tests/oracles`). |
| **golden** | A recorded result of a retired implementation, kept as a test reference (`tests/goldens/polycore.npz`). |
| **polycore** | The retired CPU checker. It is archived in `archive/cpu_checker`. |
