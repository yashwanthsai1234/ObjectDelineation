# Reward exact topology against annotator disagreement

The strongest open mix-and-match for RL fine-tuning is **ambiguity-calibrated exact-topology RLVR**. A model may emit polygon text, mask tokens or a mask. Whatever it emits is first canonicalized into exact pixel-edge rings with holes. The model is then rewarded on three things: per-instance β0/β1 matching the *set* of topologies that independent annotators drew, graded multi-ring legality, and a boundary term that is neutral inside the annotators' disagreement band (between their intersection and union) and penalized outside it. Each term is normalized separately within the GRPO group, and the topology terms count only once the overlap term clears a threshold. Each half of this already has neighbours in print, so neither half alone is a safe "first":

- Topo-R1 put a topology term into VLM GRPO, but it computes clDice on the *input* mask and weights it at 0.05 ([arXiv 2603.13054](https://arxiv.org/html/2603.13054v1)).
- FloorplanVLM already rewards binary "watertight closed polygon" validity ([arXiv 2602.06507](https://arxiv.org/html/2602.06507)).
- TraversRL already rewards a fixed-width buffered IoU ([arXiv 2607.17479](https://arxiv.org/pdf/2607.17479)).

Nobody found in this research rewards **exact β0/β1/χ of the policy's own output**. Nobody uses a **spatially varying tolerance derived from multiple raters** as an RL reward. Nobody **scores topology against a rater set rather than a single ground truth**. The combination matters because the two halves guard each other. The band term blocks the cheapest topology hack (dilating until components merge). The topology term catches the small-area errors (spurious holes, specks) that any overlap term, banded or not, barely sees.

The GPU engine is the enabler rather than the headline. Reward compute is a documented bottleneck only in lithography ([LithoGRPO](https://arxiv.org/html/2606.00228)). In VLM RL the engine's value is exact fill-rule semantics, ring hierarchy and determinism. Residual risk is real:

- VectorLLM++ and Learning to Trace remain unreadable at the reward level. My subjective estimate is a 20–30% chance that an unread paper already uses β-exact output rewards.
- Binary validity gates already push invalidity to about 0% in CAD and SVG RL. "Validity" by itself is therefore nearly solved; *topology correctness under ambiguity* is not.

## Prior art already holds four of the ten claims you might make

The 2025–2026 reward vocabulary for perception RL is overlap and distance terms plus format gates. Across roughly 30 papers, rewards are IoU, GIoU, L1, Dice, S-measure, mAP or Hungarian-matched F1. Mask pipelines nearly always hand geometry to a frozen SAM, so the policy never controls the topology it is scored on. Two examples:

- SAM-R1 tiers mask IoU at 4/3/2/0 ([arXiv 2505.22596](https://arxiv.org/html/2505.22596)).
- VisionReasoner says outright that it scores "bounding boxes and points … rather than binary masks" ([arXiv 2505.12081](https://arxiv.org/html/2505.12081v5)).

The important exception is **SimpleSeg**, which emits polygons as text and uses GSPO. Its rewards are rasterized IoU with a threshold, negative centroid MSE, and a format gate. Without clockwise ordering "no segmentation mask can be derived" ([arXiv 2601.19228](https://arxiv.org/pdf/2601.19228)). It is the natural host for the proposal.

Geometric and validity rewards in *generative geometry* are now crowded:

| Domain | Paper | Reward recipe |
|---|---|---|
| Lithography | LithoGRPO | Shot count inside flow-matching GRPO ([arXiv 2606.00228](https://arxiv.org/html/2606.00228)) |
| CAD | cadrille (ICLR 2026 oral) | R = 10·IoU − 10 for invalid code; invalidity falls from 2.1% to 0.0% on DeepCAD ([arXiv 2505.22914](https://arxiv.org/html/2505.22914v2)) |
| SVG | CTRL-S | Multiplies everything by a render-success gate and reaches a 99.85% success rate ([arXiv 2603.16189](https://arxiv.org/html/2603.16189)) |
| SVG diagrams | GeoSVG-RL | Browser-backed verifier for containment, canvas overflow and anchor accuracy, at about 18 SVGs/s on 8 CPU workers ([arXiv 2605.25447](https://arxiv.org/html/2605.25447)) |
| Floor-plan generation | Mila's floor-plan RLVR | Overlapping rooms get reward 0 ([arXiv 2605.14117](https://arxiv.org/html/2605.14117v1)) |
| 3D meshes | Mesh-RFT | Hole-sensitive Boundary Edge Ratio to build masked-DPO pairs ([arXiv 2505.16761](https://arxiv.org/html/2505.16761v1)) |

A loophole shows up in the floor-plan RLVR code: it repairs self-intersecting rooms with `buffer(0)` instead of penalizing them, silently skips rooms that throw exceptions, and never checks for gaps, holes or a single connected outline ([repo](https://github.com/ludolara/floor-plan-rlvr)).

The table below is the verdict on each ingredient you listed. "Confirmed novel" means targeted searches found nothing. It does not mean an exhaustive proof of absence. The last column gives the reason the verdict could flip.

| Ingredient | Verdict | Closest prior art, and what it actually does | What remains open | Why it could flip |
|---|---|---|---|---|
| Exact β0/β1/χ of the **policy's own output** as reward | **Confirmed novel** (moderate–high confidence) | Topo-R1: clDice on the corrupted *input* mask, w=0.05; Betti numbers used only to verify synthetic training errors ([arXiv 2603.13054](https://arxiv.org/html/2603.13054v1)). Mesh-RFT: BER in 3D, DPO. | Per-instance β0/β1 match, χ-error, per-ring hole matching on outputs | VectorLLM++, Learning to Trace, unindexed MICCAI/ISPRS papers |
| Polygon **validity** reward | **Already done (binary)** | FloorplanVLM: binary watertight gate, validity 90.2%→96.1%. Floor-plan RLVR overlap gate. cadrille, CTRL-S, GeoSVG exec/render gates. | *Graded*, hole-aware, multi-ring legality; fill-rule disagreement as a validity signal | FloorplanVLM's checker is unspecified and could count rings |
| **Annotator/ensemble band**-tolerant reward | **Confirmed novel as an RL reward** | NSD sets τ from the inter-observer 95th percentile, but as a *metric* ([arXiv 1809.04430](https://arxiv.org/abs/1809.04430)). TraversRL uses a fixed buffer. SAM-R1 and Topo-R1 use tiers. | Spatially varying band from raters or an independent ensemble, inside GRPO | A workshop paper could exist; the search was web-only |
| **Temporal** topology-consistency reward | **Confirmed novel** | Veason-R1's "consistency" reward is GT IoU over SAM2-propagated frames ([arXiv 2508.11538](https://arxiv.org/html/2508.11538)). VGGRPO uses geometric consistency for video *generation* ([arXiv 2603.26599](https://arxiv.org/html/2603.26599v2)). | Label-free Δβ flicker, ring-hash stability, flicker-area penalty | Low; but see the gaming risk below |
| **Canonical ring hash** for dedup/caching | **Novel but low value** | No RL paper hashes rollouts | Exact-match verifier; duplicate-rate diagnostic; caching for expensive rewards | Trivial to replicate; not a paper on its own |
| **Hole-aware SFT polygon data** | **Confirmed novel** (moderate) | SimpleSeg traces one clockwise outer ring and lists holes as a failure. Groundbench keeps only the largest part for 148 of 1,500 referents, covering 71.76% of the object ([arXiv 2609.27821](https://arxiv.org/pdf/2609.27821)). | Ring- and part-delimited canonical tokens | PolyFormer, VistaLLM and VectorLLM++ multi-polygon formats were not verified |
| **Batched GPU rasterize + polygonize** | **Novel as infrastructure** | NPP gives per-label pixel chains with no hole or parent field, and NVIDIA admits "a few contour pixels can escape" ([sample](https://github.com/NVIDIA/CUDALibrarySamples/blob/master/NPP%2B/findContour/findContourNPPPlus.cpp)). nvdiffrast and Kaolin are triangle-only ([nvdiffrast](https://nvlabs.github.io/nvdiffrast/)). | Fill-rule polygon rasterizer plus ring-hierarchy polygonizer on device | An engineering claim, not an ML claim |
| **Cross-domain verifier suite** | **Novel as a shared library** | Every paper ships bespoke verifiers: shapely, CadQuery subprocess, CairoSVG/Chromium, Numba plus CBC ILP. | One API: validity, β0, β1, χ, area, perimeter, vertex count, jogs, shot proxy, hash | Demand is inferred, not stated |
| Non-differentiable geometric RL rewards in general | **Already done** | LithoGRPO shot count; GeoSVG containment; OP-HRG part-in-object containment ([arXiv 2607.15374](https://arxiv.org/html/2607.15374)) | — | — |
| Topology rewards for **mask-generating diffusion/flow** | **Open** | Only LithoGRPO among geometric-output flow RL; Flow-GRPO targets T2I ([arXiv 2505.05470](https://arxiv.org/pdf/2505.05470)) | Flow-GRPO with β rewards on ambiguous segmenters | A medical compositional-reward paper uses GDPO with unknown segmentation sub-rewards ([arXiv 2609.05028](https://arxiv.org/html/2609.05028)) |

Your verifier-suite list breaks down unevenly by domain.

- **VLM grounding** is open for topology and legality rewards.
- **Floorplans** already have binary validity, but exact gap, hole and single-outline checks are open.
- **CAD** has execution gates. Sketch-profile topology, slice-wise β and volumetric χ are open. I found no CAD-RL reward with genus or watertightness terms beyond "executes and is non-degenerate".
- **SVG** is render-and-compare in pixel or embedding space. Path-level topology (component and hole counts of the rendered fill) is unchecked.
- **Lithography** has LithoGRPO's CPU shot-count proxy. NVIDIA's ILT-RL rewards only EPE and enforces curvilinear mask-rule compliance by morphological opening, not by reward ([arXiv 2602.19027](https://arxiv.org/html/2602.19027)).
- **Spatial omics and microscopy** have live RL and bandit acquisition loops, but their rewards are classification confidence or acquisition cost ([arXiv 2512.14930](https://arxiv.org/pdf/2512.14930)).
- **Rock fragmentation** has the clearest metric alignment: merge and split errors directly change fragment counts and the size distribution ([RIT](https://www.rit.edu/dirs/research/structure-aware-segmentation-post-blast-rock-fragmentation-analysis)). I found no RL loop there at all.

### The uncertain ones stayed uncertain

Several gaps were closed in this round, but these remain open:

- **VectorLLM++** advertises "reinforcement learning with verifiable rewards that directly optimize geometric quality" for polygons, masks and boxes. Only the Research Square landing page was reachable, so its reward terms are unknown ([Research Square](https://www.researchsquare.com/article/rs-10949232/v1)). Its predecessor, VectorLLM, used DPO with rejected answers made by deleting or inserting polygon points ([arXiv 2507.04664](https://arxiv.org/html/2507.04664v1)).
- **Learning to Trace** (NeurIPS 2026) cuts DRIVE vessel components from 29.1 to 4.5 with a procedural RL tracer. OpenReview blocked the fetch, so whether its reward counts components is unknown ([NeurIPS](https://neurips.cc/virtual/2026/161763)).
- **FloorplanVLM** does not name the algorithm behind its watertightness check.
- **Topo-R1's supplement** was not read.
- **GDPO** details come from snippets only.

Because of the first two, every novelty sentence in a paper should carry a dated "to our knowledge" and a footnote naming these two works.

```
                      where topology/validity is computed
                  on input or curation         on the policy's own output
                ┌──────────────────────────┬──────────────────────────────────┐
 exact counts   │ Topo-R1: Betti numbers   │ ★ PROPOSAL: per-instance β0/β1/χ,│
 (β, χ, rings)  │ verify synthetic errors  │   hole matching, graded legality,│
                │ (data curation only)     │   scored vs a RATER SET          │
                ├──────────────────────────┼──────────────────────────────────┤
 soft / proxy   │ Topo-R1 R_topo = clDice  │ Mesh-RFT BER (3D, DPO)           │
                │ on input mask, w = 0.05  │ Learning to Trace (reward ?)     │
                ├──────────────────────────┼──────────────────────────────────┤
 binary gate    │            —             │ FloorplanVLM watertight;         │
                │                          │ floor-plan RLVR overlap = 0;     │
                │                          │ cadrille / CTRL-S / GeoSVG gates │
                └──────────────────────────┴──────────────────────────────────┘
 tolerance axis:  single GT, fixed τ (SAM-R1 tiers, TraversRL buffer, NSD metric)
                  ──────────────►  ★ rater-derived, spatially varying band (open)
```

## One idea survives every filter: score topology against what raters disagree about

### The precise statement

Take a policy that emits a region: polygon text, mask tokens, or mask logits from a segmenter or generator. Canonicalize its output into exact pixel-edge rings with outer/hole labels and parent links. Compute a vector of verifiable terms against a reference set of R ≥ 2 independent annotations, or a frozen independent ensemble. From that set, define:

- the **yolk** Y = ∩A_r, the region every rater includes;
- the **egg** E = ∪A_r, the region at least one rater includes;
- the **topology set** T = {(β0(A_r), β1(A_r))}.

The reward has five parts:

1. a graded legality term;
2. a topology term that is maximal when the output's per-instance (β0, β1) lies in T and decays with the distance to its nearest member;
3. a hole-matching term over individual rings;
4. a band term that ignores pixels in E∖Y and penalizes mass outside E or missing from Y;
5. per-term group normalization (GDPO-style), with the topology and legality terms gated on a minimum band overlap.

### Why combine rather than ship either half

Topology targets are themselves ambiguous on exactly the data where topology matters. A single-ground-truth β reward punishes a policy for agreeing with rater 2 rather than rater 1. A band-tolerant overlap term alone tolerates fragmented or perforated masks inside the band, because a few spurious holes cost almost no area.

The pairing closes both gaps:

- The **egg bound** caps the "dilate until components merge" hack. That hack is the mask-level analogue of the oversized-box exploit Topo-R1 had to guard against with LocPen.
- The **topology term** catches the specks and pinholes that Ground-R1's analysis predicts area-dominated rewards will ignore: "training rewards are dominated by large visual regions" ([arXiv 2505.20272](https://arxiv.org/html/2505.20272)).

The egg-yolk construction has a thirty-year pedigree. Cohn and Gotts defined the egg as the maximal extent and the yolk as the certain core ([ResearchGate](https://www.researchgate.net/publication/2827200_The_Egg-Yolk'_Representation_Of_Regions_with_Indeterminate_Boundaries)). That helps reviewers accept the band definition.

The rater data to build it exists:

| Dataset | Raters | Size |
|---|---|---|
| LIDC-IDRI | 4 | 1,575 nodules |
| RIGA | 6 | 750 fundus images |
| QUBIQ brain growth | 7 | 34 training subjects |

Source: [multi-rater literature summary](https://www.researchgate.net/publication/355882110_Learning_Calibrated_Medical_Image_Segmentation_via_Multi-rater_Agreement_Modeling). The disagreement is large enough to matter: pairwise Dice between skin-lesion annotators is only 0.753–0.791 ([IMA++](https://arxiv.org/html/2508.09381v1)). For natural images, Boundary IoU found that two experts' contours "seldom" differ by more than 1% of the image diagonal ([arXiv 2103.16562](https://arxiv.org/abs/2103.16562)). That gives a citable default band where only one annotation exists.

### Safe claims and claims to avoid

| Claim | Status | Required citation or caveat |
|---|---|---|
| "Exact topological invariants (per-instance β0, β1, χ) and hole-aware ring legality **of the model's own output**, used as verifiable GRPO reward terms" | **Safe** | Cite Topo-R1 (input-side, soft) and Mesh-RFT (3D, DPO) as nearest; footnote VectorLLM++ and Learning to Trace |
| "Rewards whose boundary tolerance is **spatially varying and derived from annotator or independent-ensemble bands**, generalizing NSD's fixed τ to RL" | **Safe** | Cite NSD, Boundary IoU, TraversRL's fixed buffer |
| "Topology scored against the **set** of rater topologies rather than a single ground truth" | **Safe** (nothing close found) | Show first that raters actually disagree on topology (Phase 0 census) |
| "Graded topology rewards outperform binary validity (FloorplanVLM-style) and soft skeleton proxies (Topo-R1-style clDice)" | **Defensible only as an experimental result** | Needs those exact baselines |
| "Canonical ring hashing for polygon-emitting VLMs" | Safe but minor | Present as tooling |
| "First topology-aware reward in VLM RL" | **Do not claim** | Topo-R1 holds it |
| "First polygon-validity reward" | **Do not claim** | FloorplanVLM, floor-plan RLVR, cadrille |
| "First RL for connectivity / vessel topology" | **Do not claim** | Learning to Trace |
| "First automated topology score for RL fine-tuning of geometry generators" | **Do not claim** | Mesh-RFT (3D) |
| "First geometric verifier for RL" | **Do not claim** | GeoSVG-RL, LithoGRPO |
| "GPU speed makes VLM RL feasible" | **Do not claim** | Reward cost is negligible beside decoding outside lithography |

## The reward has five terms, one gate, and no self-derived bands

### Terms, normalization and gating

The design follows from five pieces of evidence:

1. **Summed multi-term rewards collapse under GRPO.** Distinct reward combinations end up with identical advantages. GDPO normalizes each term per group, sums, then batch-normalizes ([arXiv 2601.05242](https://arxiv.org/pdf/2601.05242); snippet-level).
2. **Binary terms go silent once most rollouts pass them.** FloorplanVLM ends at 96.1% validity, where a binary gate is constant across most groups.
3. **Near-miss samples treated as failures cause "advantage collapse"**, which sigmoid-smoothed numeric rewards fix ([arXiv 2601.07695](https://arxiv.org/abs/2601.07695)).
4. **Gating overlap terms by a threshold works in practice.** FloorplanVLM shrinks its interior-wall weight to 0.1 when exterior IoU < 0.3.
5. **Normalizing each instance by its own size counters the large-object bias** that Ground-R1 documents.

From these follows the pipeline below. Weights are equal after GDPO normalization unless an ablation says otherwise.

```
per prompt x: rollouts O_1..O_G ; references A_1..A_R (raters, or a FROZEN independent ensemble)
precompute once per prompt (cache by GT hash):
    Y = ∩ A_r   (yolk)      E = ∪ A_r   (egg)      T = {(β0(A_r), β1(A_r))}   rater topology set
    H* = consensus hole rings (holes present in ≥ ⌈R/2⌉ raters), via ring IoU

per rollout i:
 0  parse          fmt_i ∈ {0,1}; if 0 → every r_k,i := 0 (bottom of every term)
 1  canonicalize   M_i  = fill(O_i, even-odd) ; M'_i = fill(O_i, nonzero)   (polygon outputs)
                   M_i  = threshold(O_i)                                    (mask outputs)
                   engine(M_i) → rings, outer/hole, parents, β0, β1, χ, area, perim, V, hash
 2  terms (all in [0,1], graded)
    r_band  = 1 − ( |M_i \ E| + |Y \ M_i| ) / |E|        pixels in E\Y are neutral; per-instance
    r_topo  = exp( −min_{t∈T} ( |β0_i − t0| + |β1_i − t1| ) )        set membership, graded
    r_hole  = mean ring-IoU of Hungarian-matched holes vs H*  (1 if H* = ∅ and no holes)
    r_legal = 1 − ( n_illegal_rings + 1[ area(M_i ⊕ M'_i) > 0 ] ) / ( n_rings + 1 )
    r_bnd   = fraction of predicted boundary pixels inside (E\Y) ⊕ τ        (optional, NSD-style)
 3  gate           if r_band < κ (start κ = 0.3): r_topo, r_hole, r_legal := 0
 4  GDPO           Â_k,i = (r_k,i − mean_G r_k) / (std_G r_k + ε) ;  A_i = Σ_k w_k Â_k,i ;  batch-normalize A
```

**Graded versus binary.** Every structural term is graded, because binary terms carry no advantage once saturated. The only hard gate is parse/format, which is a precondition for having geometry at all. The even-odd versus nonzero disagreement in r_legal is my own design proposal, not a published technique. A self-overlapping or wrongly wound polygon fills differently under the two rules, so a nonzero XOR area is a cheap, exact flag that needs no segment-intersection search.

**Where the band comes from.** For datasets with one annotation (refCOCO, SA-1B), build Y and E by eroding and dilating the ground truth by d. Set d to 2% of the image diagonal from the Boundary IoU expert-agreement result, or 0.5% for Cityscapes-quality labels. T then becomes the β of the eroded, original and dilated masks. That is a morphologically robust topology set, which is also my own proposal.

**Never derive the band from the policy itself.** Conformal bands are calibrated on the model's *own* errors. Used as a reward tolerance, they create a feedback loop: a worse model gets a wider band and a more lenient reward. Bands must come from annotators, an independent frozen ensemble, or a generator such as SSB++. SSB++ reaches a GED of 0.208 on LIDC at 50 function evaluations ([arXiv 2509.17187](https://arxiv.org/html/2509.17187)).

### Hacking risks and guards

The reward-hacking record for IoU-family rewards is specific:

- VLM-R1's mAP reward made the model "predict all categories" ([arXiv 2504.07615](https://arxiv.org/html/2504.07615)).
- Seg-R1 found "S-measure alone leads to reward hacking" ([arXiv 2506.22624](https://arxiv.org/html/2506.22624)).
- Topo-R1 added LocPen against "oversized bounding boxes that trivially capture topological changes".

The table maps each foreseeable exploit to the term that invites it and the guard.

| Exploit | Term that invites it | Guard |
|---|---|---|
| Dilate or bridge until components merge, hitting β0* = 1 | r_topo | r_band penalizes mass outside E; the gate zeroes topology when the band fails |
| Emit specks to inflate β0, or pinholes to inflate β1 | r_topo | Hungarian matching of components and holes by IoU (unmatched ones count as errors); r_hole scores individual rings, not counts |
| Fill every hole to avoid hole errors | r_band, r_legal | r_topo and r_hole penalize β1 = 0 when every rater drew a hole |
| Self-intersecting polygons whose fill happens to score well | r_band | Score under both fill rules; disagreement lowers r_legal; never call `make_valid` |
| Trivial legal shape (a square) to bank r_legal | r_legal | Gate on r_band ≥ κ |
| Wide-band leniency | r_band | Bands only from raters or a frozen independent reference; log band width against reward |
| Hug the yolk (minimum-area safe answer) | r_band | Neutral band means no bonus for shrinking; add r_bnd so the boundary must sit inside E∖Y |
| Overlong outputs or vertex padding | all | Dr. GRPO-style length-bias correction ([arXiv 2503.20783](https://arxiv.org/abs/2503.20783)); avoid hard vertex caps, which SimpleSeg found degrade performance |

### Controls

On Qwen-family models, random rewards gave +21.4 points on MATH-500 against +29.1 for ground truth. The same effect "often fail[s]" on Llama or OLMo ([arXiv 2506.10947](https://arxiv.org/abs/2506.10947)). Every arm therefore needs:

- a **random-reward control**;
- a **format-only control**;
- ideally **one non-Qwen backbone**;
- **pass@k curves against the base model**, because RLVR mostly sharpens what the base can already sample ([arXiv 2504.13837](https://arxiv.org/abs/2504.13837v4)).

## The system needs one disaggregated GPU reward server, not a faster findContours

### The two input paths

The reward server has to handle two output types.

**Polygon-emitting policies** (SimpleSeg, FloorplanVLM, VectorLLM-style) need polygon→mask rasterization with explicit fill rules. They should never go through CPU repair. In this research's own benchmark (Apple M3, Shapely 2.2.0 / GEOS 3.14.1, single thread):

- `is_valid` costs about 150 ns per vertex on valid rectilinear rings.
- Default `make_valid` on random self-intersecting polygons takes **9.3 ms at 32 vertices, 71 ms at 64 and 687 ms at 128**.
- A GEOS IoU reward on jittered 64-vertex polygons costs 0.88 ms per pair.

These are in-house measurements, not published figures. The conclusion is a design rule, not a speed pitch: score invalidity directly and rasterize with a fill rule. Off-the-shelf GPU options do not fit:

- nvdiffrast, Kaolin and PyTorch3D are triangle rasterizers. They need triangulation, which is ill-defined for self-intersecting output.
- torch_contour is differentiable but dense and has no hole argument ([GitHub](https://github.com/antoinehabis/torch_contour)).
- Siemens reports up to 290× on H100 for Manhattan masks with a closed-source CUDA polygon rasterizer ([Siemens](https://blogs.sw.siemens.com/calibre/2026/04/13/why-a-gpu-rasterizer-matters-for-computational-lithography-in-both-performance-and-precision-at-scale/)). That shows a custom tile kernel is the right shape.

**Mask-emitting policies** (mask tokens, SAM-in-loop, diffusion or flow generators) need the engine's native direction: bit-quads for χ, area, perimeter and vertex count, then ring linking for β0, β1, parents and the canonical start. No open GPU library emits rings with hole hierarchy. GitHub code search for CUDA `findContours` and `polygonize` found only application files and marching cubes.

```
┌──────────────────────── rollout tier (vLLM/SGLang workers) ────────────────────────┐
│ π_θ → G outputs per prompt:  polygon text │ mask tokens │ mask logits (diffusion)  │
└──────────────┬──────────────────────────────────────────────┬──────────────────────┘
               │ strings (host)                               │ tensors (DLPack, device)
               ▼                                              ▼
┌──────────────────── reward server: 1 GPU, disaggregated (ROLLART-style) ───────────┐
│ P  parse + format gate → vertex buffers [B·G, rings, K]                            │
│ R  fill-rule rasterizer: one block per polygon bbox tile, edges in shared mem,     │
│    even-odd AND nonzero in one pass  → M, M'                                       │
│ Q  bit-quad classify + block histogram → V, area, perimeter, χ  (per instance)     │
│ L  compact corners → stable x-sort → pair → pointer-jump ring IDs →                │
│    outer/hole by min corner → CCL parent → β0, β1, ring areas, canonical hash      │
│ M  match vs reference set: Hungarian on component IoU; hole-ring IoU vs H*         │
│ B  band ops vs cached (Y, E, T) for the prompt: |M\E|, |Y\M|, boundary-in-band     │
│ T  term tensor r[B, G, K] → per-term group normalize (GDPO) → advantages           │
│ C  hash cache: ring hash → expensive terms (litho sim, CAD exec); dup-rate metric  │
└────────────────────────────────────────────────────────────────────────────────────┘
               │ advantages [B·G] + per-term logs (β error, legality, band width)
               ▼
                                   learner (GRPO / GSPO / Flow-GRPO)
```

### How much compute it saves

The honest throughput argument is modest for VLMs. At SimpleSeg scale (about 6.4M reward evaluations), CPU polygon rasterization plus IoU comes to roughly 1.8 CPU-hours. That is my earlier estimate, not a measurement. Against it, LENS reports about 1.5 s of MLLM latency per sample ([arXiv 2508.14153](https://arxiv.org/html/2508.14153)).

What does carry weight is keeping the reward tier on its own GPU. In agentic LLM RL, ROLLART moved reward computation off the main pipeline. Reward-GPU utilization rose from 6% to 88%, and per-step rollout time halved ([arXiv 2512.22560](https://arxiv.org/pdf/2512.22560)).

The hash cache matters only where a single reward evaluation is expensive. Lithography and CAD are the examples. LithoGRPO's released shot counter is still a Numba plus PuLP/CBC ILP with a 300 s time limit ([repo](https://github.com/laiyao1/LithoGRPO)). cadrille runs every CadQuery exec in a separate process with a 3 s timeout "to overcome CadQuery memory leaks" ([repo](https://github.com/col14m/cadrille)).

A CPU reference (OpenCV or skimage fill plus the existing CPU ring linker) is fine for Phases 0–1 below. The GPU server earns its place in multi-instance dense prompts, in diffusion-scale RL, and in lithography.

## A four-phase program decides the question for about 8 GPUs

The program is built to fail fast. Phase 0 can kill the idea before any RL run if raters almost never disagree on topology, or if hole and multi-part targets are rare. In either case the "rater-set" half collapses into single-GT topology, which is still novel but narrower.

| Phase | Model(s) | Data | Arms | Primary metrics | Kill / continue criterion |
|---|---|---|---|---|---|
| **0. Census and engine validation** (≈2–3 wk) | none | SA-1B subset; LIDC, QUBIQ, RIGA rater sets; Groundbench's 148 multi-part referents | — | Share of masks with β1 ≥ 1 or β0 ≥ 2; rater topology-disagreement rate per dataset; engine certificate pass rate (#outer = #fg components, #outer − #holes = χ) | Continue if hole or multi-part targets are ≥ a few % and raters disagree on β in a non-trivial share of cases |
| **1. Polygon VLM, single-GT** (≈4–6 wk) | Qwen2.5-VL-3B/7B or Qwen3-VL-4B, plus one non-Qwen backbone | refCOCO/+/g with a hole and multi-part enriched split; SA-1B topology-stratified subset | A0 SFT outer-ring (SimpleSeg-style) vs A0′ SFT hole-aware tokens; then RL: A1 random, A2 format-only, A3 IoU-only (SimpleSeg recipe), A4 IoU + binary validity (FloorplanVLM-style), A5 IoU + clDice on output (Topo-R1-style), A6 IoU + fixed buffer (TraversRL/NSD τ), A7 full summed, A8 full + GDPO | cIoU; Groundbench legality rate; β0/β1 error; hole recall; pass@k vs base | A8 must cut β error at equal cIoU (within 1 point) vs A3 and A4, and beat A1/A2 by more than seed variance |
| **2. Multi-rater ambiguity** (≈4–6 wk) | Same VLM (polygon) and a mask-token policy | LIDC-IDRI, QUBIQ, RIGA (disc minus cup is a natural annulus, so β1 = 1: my observation), IMA++ | A3, A6, A8 with single-GT topology, A8 with rater-set topology and band, A9 = A8 with a **self-derived conformal band** (negative control) | NSD at rater τ; rater-set topology hit rate; β error vs each rater; Dice vs the rater average (QUBIQ-style) | Rater-set A8 must beat single-GT A8 on hit rate without losing NSD; A9 should show leniency drift |
| **3. Generative and connectivity** (≈6–8 wk) | Flow-GRPO on an ambiguous-segmentation generator (SSB or PHiSeg class); a mask-emitting tubular segmenter | LIDC (GED, D_max); DRIVE (Learning to Trace's benchmark) | Base vs IoU-reward vs full | GED, D_max, CI; components per image vs 29.1→4.5 reference | Must improve GED or Betti error without a GED regression |

Setup notes:

- **Group size and hardware.** Topo-R1 ran GRPO with G=4 on Qwen3-VL-4B. SAM-R1 trained on 3,000 samples with 8 rollouts on 8×A100. So G=8 on 8 A100/H100-class GPUs is a realistic minimum for Phases 1–2.
- **Seeds.** Report at least three seeds per arm.
- **Diagnostics to log.** Per-term reward curves, the distribution of |M|/|E|, the component-count histogram and the in-group duplicate-hash rate. These expose hacking before the headline metrics do.

The phases are deliberately ordered from cheapest to most expensive. Phase 0 is also a publishable contribution on its own. No published number exists for the share of SA-1B's 1.1B masks that have holes or multiple parts ([Meta SA-1B](https://ai.meta.com/datasets/segment-anything/)).

## Eight secondary combinations, ranked by evidence over effort

| Rank | Idea | Evidence for | Main risk | Effort |
|---|---|---|---|---|
| 1 | **GPU rank-faithful mask-cost proxy for lithography RL**: vertex, jog, hole and island counts plus a rectangle-partition bound per rollout, replacing LithoGRPO's CPU ILP and adding the manufacturability terms NVIDIA's ILT-RL leaves out | Only *documented* reward bottleneck: 34–148 s per mask cut to 0.16–0.30 s ([arXiv 2606.00228](https://arxiv.org/html/2606.00228)); NVIDIA rewards EPE only | The n/2 + h − g − 1 partition formula is unverified; overlapping-shot semantics differ; must hit Spearman ρ ≈ 1 against adabox | Medium |
| 2 | **Hole-aware canonical polygon tokens plus a SA-1B topology census** as the SFT foundation of the main idea | Holes are SimpleSeg's named failure; Groundbench truncation to 71.76% coverage | Token-length growth; SimpleSeg found hard length limits hurt | Low–medium |
| 3 | **Cross-domain "geometry unit tests" library**: replace `buffer(0)` in floor-plan RLVR, add path-topology checks to SVG RL, sketch-profile β to CAD RL | No shared verifier exists; documented loopholes | Demand inferred; each domain wants bespoke semantics | Medium–high |
| 4 | **Label-free temporal topology reward** for VOS RL: Δβ flicker, flicker-XOR area, ring-hash stability for static objects, added to Veason-R1 or VideoSEG-O3 pipelines | No self-consistency reward exists; the mTC metric dates from 2020 ([CVPRW](https://openaccess.thecvf.com/content_CVPRW_2020/papers/w20/Varghese_Unsupervised_Temporal_Consistency_Metric_for_Video_Segmentation_in_Highly-Automated_Driving_CVPRW_2020_paper.pdf)) | Static or over-smoothed masks score perfectly; must pair with GT IoU and a flow warp | Medium |
| 5 | **Topology-rewarded Flow-GRPO for mask generators** (ambiguous segmenters, synthetic-data generators, label-map diffusion) | Open slot besides LithoGRPO; SPO shows early denoising steps set layout ([arXiv 2406.04314](https://arxiv.org/html/2406.04314v1)) | Credit assignment over steps; compute | High |
| 6 | **Fragment-count/PSD reward for rock fragmentation segmenters** | Merges and splits change P80 and fines directly; YOLO-seg runs at about 15 FPS ([arXiv 2507.20126](https://arxiv.org/html/2507.20126v1)) | No RL loop or multi-rater data exists; small datasets | Medium |
| 7 | **Topology events as state or reward in smart-microscopy RL** (division = 1→2 components; vacuole = a new hole) | Live RL and bandit loops with 100–200 ms budgets ([arXiv 2504.00047](https://arxiv.org/html/2504.00047v2)) | Current rewards are acquisition-centric; speculative fit | Medium |
| 8 | **Hash-based rollout dedup and caching** | No paper does it | Value only where rewards are expensive; otherwise a diagnostic | Low |

The straight-through polygon layer (exact forward pass, soft bit-quad backward pass) stays off this list. It remains unbuilt anywhere and far from the current code. Topograph-style topology losses already cost 96–603 ms per call ([arXiv 2411.03228](https://arxiv.org/html/2411.03228v2)), so an exact GPU β gate that skips the loss when topology already matches is a cheaper first step.

## Counter-evidence: what argues against the main idea

**Speed only matters in lithography.** LithoGRPO is the only paper that documents reward cost as a blocker. GeoSVG-RL's 18 SVGs/s CPU verifier is the next candidate, but its batch semantics make the cost unconfirmed. In VLM segmentation RL the decoder dominates. Any pitch that leads with GPU throughput for VLMs will be rejected by reviewers who do this arithmetic.

**IoU may already capture most of the signal.** Several findings support this:

- SimpleSeg argues that contours are many-to-one with respect to masks and that RL handles that flexibility on its own.
- Perception-R1 finds RL comparable to SFT on grounding ([arXiv 2504.07954](https://arxiv.org/html/2504.07954)).
- Topology may improve as a side effect of IoU optimization; no paper measures this.

The strongest rebuttal is Groundbench. At 64 vertices Gemini's filled IoU is 37.7, but legal-polygon coverage is 7.5. Overlap and legality therefore diverge sharply for polygon emitters. Even so, the A3 versus A8 comparison is the experiment that decides the question, and it could come out null.

**Validity is close to solved by binary gates.** cadrille drives invalidity to 0.0–0.2%, CTRL-S reaches 99.85% render success, and FloorplanVLM reaches 96.1% validity. A graded-validity term alone is a marginal contribution. That is why the core idea rests on *topology under ambiguity*, with legality as a supporting term.

**Raters may rarely disagree on topology.** Annotator disagreement is documented in Dice and contour distance, not in Betti numbers. If LIDC or QUBIQ raters almost always agree on β, the rater-set term reduces to single-GT topology. A further caution: KαLOS argues that raw IoU "reads around 0.95 even for annotators who never looked at the image" on typical detection corpora ([arXiv 2603.27197](https://arxiv.org/html/2603.27197); snippet-level). That is a reminder that agreement statistics need chance correction before they set a band.

**The temporal reward is gameable by static masks.** A model that freezes its output across frames maximizes any self-consistency term. A flow-warped mask-level mTC term may already capture most of what ring-level stability adds. The polygon contribution must be shown empirically, for example through boundary-jitter reduction that mTC misses.

**Multi-rater data is small and medical.** QUBIQ brain growth has 34 training subjects. Remote sensing has no multi-annotator building or field dataset, and 2026 building-extraction papers say explicitly that no inter-annotator analysis was available ([Remote Sens.](https://doi.org/10.3390/rs18172971)). Natural-image work falls back to fixed-d bands, which are TraversRL-adjacent.

**SAM-in-the-loop pipelines dilute the idea.** When the policy emits only boxes and points and SAM draws the mask, a topology reward mostly grades SAM. The idea bites hardest on policies that emit geometry directly: SimpleSeg-style text polygons, mask tokens, and generators.

**Residual prior-art risk is concentrated.** It sits in two unreadable papers (VectorLLM++, Learning to Trace), one unspecified checker (FloorplanVLM), and unindexed venues (MICCAI 2026, ISPRS, TGRS). The Exa search tool was unavailable for the whole session, and OpenReview and CVF blocked fetches.

## Conclusion

The research moved the defensible novelty from "topology rewards" to **"what counts as correct topology when experts disagree."** That shift makes the idea stronger, not weaker. The exact-counts engine supplies something soft proxies such as clDice and Fast-χ cannot: which ring is the spurious hole, and whether a prediction's (β0, β1) matches *any* rater. The rater set supplies something no exact verifier has: a principled tolerance that does not punish legitimate alternatives.

The practical implication is to sequence the work by information value, not engineering appeal:

1. Run the Phase 0 census (rater topology disagreement plus SA-1B holes and multi-part share). It costs weeks and either opens or closes the main bet.
2. Run the lithography proxy in parallel. It is the one place where the GPU engine's speed is the product.

If both land, the engine's identity becomes "the verifier layer for spatial RLVR". In that role, determinism, fill-rule semantics and ring hierarchy are the specification, and throughput is a property it happens to have.
