# State of the art, November 2026: VLMs that output geometry, RL with verifiable visual rewards, and medical (multi-rater) segmentation

## Skills used

| Skill | What it contributed |
|---|---|
| `huggingface-skills:huggingface-best` | Pulled the official HF benchmark list (58 datasets). The only grounding leaderboard is `likaixin/ScreenSpot-Pro`; I fetched its ranking. I searched the Hub for every current small-VLM family and got parameter counts and licenses from the HF API (`safetensors.total`, `cardData.license`). This found the Qwen3.5 (Feb 2026), Gemma 4, LFM3-VL-2B (Oct 2026), MiniCPM-V 4.7, Molmo2 and Moondream 3.1 releases. I parsed the RefCOCO, CountBench and ScreenSpot rows out of the HTML tables on the Qwen3.5-2B, Qwen3.5-4B and LFM3-VL-2B model cards (Section 2). |
| `huggingface-skills:huggingface-papers` | Used `/api/papers/2601.19228` for SimpleSeg metadata (published 2026-01-27, 19 upvotes, GitHub, project page). `/api/models?filter=arxiv:2601.19228` found the released checkpoints **`sthui/SimpleSeg-Qwen2.5-VL`** and **`sthui/SimpleSeg-Kimi-VL`**: a ready zero-shot and fine-tune baseline. Two novelty searches on `/api/papers/search` ("multi-rater … reward reinforcement … segmentation", "polygon VLM medical segmentation GRPO") returned no matching prior work. |
| `huggingface-skills:huggingface-vision-trainer` | SAM/SAM2 fine-tuning knowledge for the specialist baselines: `facebook/sam2.1-hiera-small` is the suggested start; by default only the mask decoder trains; box prompts in xyxy; DiceCE loss. This defines the "SAM2 fine-tuned with box prompts on our LIDC/RIGA crops" baseline in Section 6 (local run on the 5090, not HF Jobs). |
| `nvidia-skills:nvidia-skill-finder` | Live catalog check (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`). Relevant entries, recommended only and not installed: `nv-segment-ct`, `nv-segment-ct-finetune`, `nv-segment-ctmr` (VISTA-3D-style CT segmentation, a possible extra LIDC specialist baseline); `dicom-series-to-volume`, `dicom-series-preflight`, `dicom-metadata-extract`; `nv-reason-ct`, `nv-reason-cxr`; `tao-generate-image-grounding`, `tao-generate-referring-expressions`. No catalog skill covers geometry verifiers or VLM RL. |
| `engineering:debug` | Structured reproduce/isolate/diagnose pass on "our 0.645 against literature numbers". Root cause: IoU_tc is a band-tolerant metric (overlap with the doctors' band, from the code notes `docs/research/notes/gpu_scoring.md` and `final_eval_2026-11-08.md`). Literature LIDC Dice ranges from 0.49 to 0.92 by protocol. Fix: report Dice against majority (0.701), GED, D_max and HM-IoU on a stated protocol (Section 3.2, "Debug note"). |
| `engineering:code-review` | Checked the survey's claims about the project against the repo, not against memory: metric definitions (IoU_tc, f5, precision/recall against the strict majority), test-set numbers and caveats (the SFT baseline was not measured on the full test split). I also checked the main thread's lookup: it mislabelled the Qwen3.5-2B RefCOCO value. **88.2 is the Qwen3-VL-4B column; Qwen3.5-2B is 84.8 / 84.3**, the same as Qwen3-VL-2B (84.8). |
| `engineering:architecture` | ADR-style decision record for the backbone (Section 2.1, "ADR: backbone"). |
| `engineering:system-design` | Requirements → components → trade-offs sketch for "the exact geometry verifier as a reusable RL reward service" (Section 6.1). |
| `engineering:documentation` | Report structure: bottom line first, numbers in tables, sources at the end, unverified items labelled. |

Survey for ObjectDelineation. The project: Qwen3-VL-2B with LoRA writes polygon JSON on a 0-1000 grid for LIDC-IDRI nodules
and RIGA disc/cup. It is trained with SFT and then GRPO against an exact CUDA checker that scores legality, IoU against a
set of expert masks, boundary F, and Betti-number topology.
Repo numbers (from `docs/research/notes/final_eval_2026-11-08.md`, shifted crops, 8 samples per item, 2,200 test items):
IoU_tc 0.645 (median item 0.756), f5 0.816, Dice against the majority outline 0.701, precision 0.695, recall 0.748,
legal 99.4 %. The coordinator gives the doctor ceiling as about 0.73 IoU.

Sourcing: every number below comes from the cited page as fetched on 2026-11-08, unless it is marked
"(from memory, not re-checked)". Vendor model-card numbers are flagged as vendor-reported.

---

## 0. Bottom line

1. **Closest prior art is SimpleSeg** (arXiv 2601.19228, Jan 2026, ICLR 2026 submission). It is a decoder-free MLLM that
   writes polygon points as text, trained SFT then RL (GSPO) with an IoU reward. It lists holes as a limitation.
   ObjectDelineation must cite it and set itself apart on four points: (a) a set-valued multi-rater reward, (b) exact
   topology with holes (b0/b1), (c) a certified, bit-exact GPU verifier, (d) medical data with real inter-rater ambiguity.
2. **Nobody found uses a multi-rater set as an RL reward.** Multi-rater segmentation work uses supervised losses (MRNet,
   D-Persona, ProSeg, RPDL) or generative samplers (Prob U-Net, PHiSeg, diffusion). Medical RL segmentation papers
   (MedReasoner, CORE-Seg, IBISAgent, MedSAM-Agent, MedVol-R1, PathChat-SegR1) all score against one mask, and all
   route geometry through a SAM-family decoder. None writes polygons.
3. **Nobody found uses exact Betti numbers as a reward for a delineating VLM.** Topo-R1 (2603.13054) is the nearest: GRPO
   with a clDice term at weight 0.05, for *detecting* topology errors in masks, not for drawing outlines. CORE-Seg
   penalises disconnected regions with a heuristic.
4. **Medical VLMs cannot localise out of the box.** On a SLAKE grounding split, MedGemma-4B gets box IoU 0.05 and
   Lingshu-7B 0.07. Both are below a fixed centre box (0.10). GPT-5.5 gets 0.24 (2604.27720). A 2B model that reaches
   Dice 0.70 on LIDC with polygons is a strong story, provided it is benchmarked against the right specialists.
5. **The comparison that decides "top of field" is missing.** We need nnU-Net 2D, MedSAM/MedSAM2 with box prompts, and
   SAM 3 on the *same* LIDC/RIGA splits. We also need the standard multi-rater metrics (GED, HM-IoU, D_max). IoU_tc is a
   band-tolerant metric of this project. It cannot be compared with published Dice or IoU.
6. **Model choice.** The safest upgrade is **Qwen3-VL-4B**: same architecture as the current model, and RefCOCO-avg 88.2
   against 84.8 for Qwen3-VL-2B on Qwen's own card. **Qwen3.5-2B scores 84.8/84.3, a tie with the current model.**
The "88.2" in the main-thread lookup is the Qwen3-VL-4B column. **Qwen3.5-4B** (Feb 2026, Apache-2.0, RefCOCO-avg 88.1) is a hybrid
   Gated-DeltaNet model. It would break several trainer optimisations and needs vLLM ≥ 0.17 and flash-linear-attention.
   LFM3-VL-2B (Oct 2026) posts the best vendor-reported small-model grounding, but its license is not Apache.
7. **Test-time scaling.** Rendered-overlay self-correction is mostly a mirage for grounding (2606.13156v2: the model
   ignores the overlay). Use parallel sampling plus selection instead: IoU-cluster voting (StAR), a learned IoU critic
   (Reason-Twice, SegAgent PRM), or checker gates that need no ground truth (legality, topology prior). The 8 samples per
   image can also serve as a *distribution* and be scored with GED against the raters. This is a unique angle.

---

## 1. VLM → geometry output: methods, numbers, and trade-offs

### 1.1 Taxonomy

| Output form | Representative work | How geometry leaves the model | Strengths | Weaknesses |
|---|---|---|---|---|
| Polygon / contour as text coordinates | PolyFormer, VisionLLM, VistaLLM, Florence-2, **SimpleSeg**, **this project** | Vertex sequence as digits or loc tokens | No decoder; human-readable; verifiable exactly; RL-friendly (sequence reward); token cost grows with boundary length, not resolution | Long sequences fail (exposure bias); vertex order matters; holes and multi-part shapes need grammar; VisionLLM v1 was weak, and v2 went back to decoders |
| Text-as-mask (patch labels) | Text4Seg (ICLR 2025) | 16×16 or 32×32 grid of labels, R-RLE compressed | No decoder; easy SFT | Coarse (needs SAM/CRF refiner to compete); no topology guarantees |
| VQ mask tokens | PaliGemma / PaliGemma 2; Gemini 2.5 (likely similar) | 4 loc tokens + 16 VQ-VAE codes → 64×64 mask inside the box | Compact (20 tokens) | 64×64 inside the box is coarse; codes are opaque; the decoder is frozen |
| [SEG] token + mask decoder | LISA, LISA++, GLaMM, PixelLM, OMG-LLaVA, Sa2VA, UniPixel | Hidden state prompts SAM/SAM2 decoder | Best RefCOCO numbers; pixel-precise | Extra module; catastrophic forgetting risk (LISA-7B MME 1/1 per Sa2VA); geometry is not verifiable token by token |
| MLLM → box/points → frozen SAM, trained by RL | Seg-Zero, Seg-R1, SAM-R1, VisionReasoner, StAR, GenSeg-R1, LENS, all medical RL papers | Text box + points; SAM makes the mask | Strong zero-shot; RL only on the MLLM | Upper-bounded by SAM; SAM in the reward loop is expensive and approximate; error cascades from the box |
| Concept-prompted segmenter | SAM 3 (ICLR 2026) | Noun phrase → all instance masks | SOTA open-vocabulary segmentation | No reasoning or free-form language beyond noun phrases (SAM 3 Agent adds an MLLM) |

### 1.2 Key papers and numbers

**Polygon / point text output (our family)**
- **PolyFormer** (CVPR 2023, arXiv 2302.07387, https://github.com/amazon-science/polygon-transformer). Seq2seq with a
  *regression* decoder: floating-point vertices, no quantisation. Multi-polygon via separator tokens. RefCOCO val
  mIoU 75.96 (B) / 76.94 (L); oIoU 74.82 / 75.96. It is not an LLM, but it was the first to show polygons beating mask
  methods on all three RefCOCO sets.
- **VistaLLM** (CVPR 2024, https://shramanpramanick.github.io/VistaLLM/). LLM writes contour points. *Gradient-aware
  adaptive contour sampling* beats uniform sampling by 3-4 mIoU. Relevant to how our targets are polygonised.
- **Florence-2** (CVPR 2024, arXiv 2311.06242; from memory, not re-checked). Outputs referring-segmentation polygons with
  1,000 quantised location tokens.
- **SimpleSeg** (arXiv 2601.19228, 27 Jan 2026; https://simpleseg.github.io, https://github.com/songtianhui/SimpleSeg).
  Points as normalised text `[[x,y],...]`, clockwise, Suzuki-Abe tracing plus ε-simplification. SFT on about 800k
  RefCOCO-family samples, then **GSPO** with mask-IoU reward (zero below a threshold), centroid-distance reward, and a
  format gate. Results: RL adds **+9.7 / +9.8 / +10.5 gIoU** over SFT on the three refCOCO val sets (up to +13.0 with
  pre-training). Kimi-VL backbone: average cIoU **74.8**, above LISA, PixelLM, GSVA and decoder-free Text4Seg/UFO. Token
  sweet spot is about **221 tokens**: 35.6 cIoU at about 78 tokens, falling to 72.5 at about 859. A hard length penalty
  *hurts* (66.7 against 76.9 gIoU). Unordered points fail to parse. Stated limitations: highly curved shapes, sharp
  corners, thin parts, **holes**. → This is our most direct comparator.

**Text-as-mask / VQ tokens**
- **Text4Seg** (ICLR 2025, arXiv 2410.09855). 16×16 semantic descriptors, R-RLE cuts length by 74 % and speeds inference
  3×. RefCOCO/+/g average cIoU 75.4 with InternVL2-8B, *with a SAM refiner*. Without the refiner: 67.5 → 70.1 with CRF,
  71.4 with 32×32 descriptors.
- **PaliGemma / PaliGemma 2** (arXiv 2412.03555). 1,024 `<locXXXX>` tokens plus 128 `<segXXX>` VQ-VAE codewords. A
  segmentation answer is 4 loc tokens and 16 seg tokens, decoded to a 64×64 mask inside the box
  (big_vision evaluator code).
- **Gemini 2.5** (Apr 2025) returns `box_2d` on a 0-1000 grid plus a base64 PNG mask. Simon Willison measured
  123-353 output tokens, so a compact internal mask code (PaliGemma-like) is likely (inference, not confirmed). On
  SA-Co/Gold, Gemini 2.5 gets cgF1 13.0 against 54.1 for SAM 3.

**[SEG]-token families** (RefCOCO val / RefCOCO+ val / RefCOCOg val cIoU, from the UniPixel NeurIPS 2025 tables)
- LISA-7B (CVPR 2024) 74.9 / 65.1 / 67.9; ReasonSeg val gIoU 61.3 (fine-tuned).
- GLaMM-7B (CVPR 2024) 79.5 / 72.6 / 74.2. OMG-LLaVA-7B (NeurIPS 2024) 78.0 / 69.1 / 72.9. PixelLM-7B 73.0 / 66.3 / 69.3.
- Sa2VA-8B (arXiv 2501.04001) 81.9 / 76.5 / 78.9; Sa2VA-26B 82.5 / 78.8 / 79.7.
- **UniPixel** (NeurIPS 2025, https://github.com/PolyU-ChenLab/UniPixel): 3B 80.5 / 74.3 / 76.3 co-trained, 81.9 / 75.3 / 77.2
  fine-tuned. 7B fine-tuned 83.0 / 77.8 / 78.7. ReasonSeg val gIoU 64.0 (3B).

**MLLM-prompts-SAM with RL** (ReasonSeg test gIoU, zero-shot)
- **Seg-Zero** (arXiv 2503.06520): Qwen2.5-VL + SAM2, pure GRPO, 9k RefCOCOg samples. 3B 56.1, 7B **57.5**, against
  LISA-7B 48.7. SFT *hurt* out-of-domain results: ReasonSeg 47.6 → 44.9. RL+CoT gave 53.8. Qwen2-VL-2B as reasoner:
  only 37.2. RefCOCO testA cIoU 80.3 (7B).
- **Seg-R1** (arXiv 2506.22624): GRPO trained on COD/SOD only. Zero-shot RefCOCOg test 71.4 cIoU, ReasonSeg 56.7 gIoU.
  **Reward-hacking lesson:** S-measure alone led to all-black masks. The fix was 0.7·IoU + 0.3·S.
- **SAM-R1** (arXiv 2505.22596): SAM as the reward provider, tiered IoU reward (4 if >0.8, 3, 2, else 0), 3k samples.
- **StAR** (arXiv 2603.14382, Jul 2026): adds a tier (5 if IoU >0.9, 1 if 0.3-0.5) and IoU-clustered majority-vote
  test-time scaling.
- **LENS** (2508.14153), **GenSeg-R1** (2602.09701; SAM 2 IoU reward at weight 5, negative-point validity at weight 10),
  **RSAgent** (2512.24023; multi-turn tools, 66.5 gIoU ReasonSeg zero-shot).

**Grounding RL (boxes)**
- **Visual-RFT** (ICCV 2025, arXiv 2503.01785): IoU reward. LISA reasoning grounding: Qwen2-VL-2B bbox mIoU 26.9 →
  SFT 28.3 → RFT 37.6, from 239 images.
- **VLM-R1** (arXiv 2504.07615): 3B RL scores 63.16 on LISA-Grounding OOD against 54.82 for SFT. OVDEval nms-AP 31.01
  against 26.50.
- **Perception-R1** (arXiv 2504.07954, StepFun): Qwen2-VL-2B. RefCOCO 89.1, COCO AP 31.9. **Key ablations:**
  "thinking" *hurts* perception (RefCOCO 75.1 with thinking against 89.1 without). RL beats SFT mainly on
  high-perplexity tasks (counting, detection). For low-perplexity grounding, RL ≈ SFT. Polygon delineation with
  ambiguous boundaries is high-perplexity, which supports our RL choice.
- **Qwen3-VL** (tech report arXiv 2511.21631, Nov 2025): boxes and points, relative 0-1000 coordinates. RefCOCO-avg
  (8 splits, from the EGM paper's table): 2B-Instruct 84.2, 4B 87.3, 8B 88.3, 235B 90.5. Thinking variants do *not* beat
  Instruct on grounding (−0.6 / −0.1 / −0.5).

**SAM 3** (arXiv 2511.16719, Meta, Nov 2025; ICLR 2026; https://github.com/facebookresearch/sam3). Promptable concept
segmentation. SA-Co/Gold instance cgF1 **54.1** (human 72.8, OWLv2 24.6, Gemini 2.5 13.0). LVIS mask AP 48.5. COCO box
AP 53.6. CountBench 93.8 %. Note: **SA-Co/Gold has 3 annotations per example and scores against the best-matching one
(oracle)**. This is the closest mainstream precedent for set-valued *evaluation*, and a good citation for our set-valued
reward. SAM 3.1 (Mar 2026) is faster.

### 1.3 Coordinate tokenisation research

- **GETok** (CVPR 2026, Ren et al.): same Qwen2.5-VL-7B and data. RefCOCO testA: text coordinates 92.9, bin tokens
  92.3, 2-D grid tokens 93.0. *Text digits are as good as bin tokens*, and grid tokens win only slightly.
- **Hi-Token** (arXiv 2608.03471): axis-specific hundreds/tens/ones digit tokens (60 types) against flat
  `<0>…<999>`. **+14.1 mIoU / +14.9 P@0.5** over flat location tokens under matched training. Supervision per token type
  is about 50× denser. Their GRPO reward (Hi-GAR) has multi-scale IoU tiers and a validity gate.
- **MolmoPoint** (arXiv 2603.28069): grounding tokens select an image patch, then a sub-patch, then a 3×3 location,
  coarse to fine. PointBench 70.7 (SOTA), ScreenSpot-Pro 61.1. 3 tokens per point instead of 8.
- **Ruler tokens** (arXiv 2510.03230): auxiliary tokens share RoPE positions with image patches. ScreenSpot-Pro
  31.1 → 37.2.
- **Implication for us:** a 0-1000 text grid is a sound default. Plain text digits are not worse than bin tokens.
  The 0-1000 grid matches Qwen's pretraining. Hi-Token-style digit structure already exists in our format: "537" tokenises
  to digits in Qwen. A cheap gain may be per-vertex grounding to patches (MolmoPoint-style), but that changes the
  architecture. Keep it as future work.

### 1.4 "Polygon as text" against "mask decoder head" for this project

The polygon-as-text route is the right choice here because the reward is **exact and verifiable on the token output
itself**. Legality, Betti numbers and ring validity are properties of the polygon, not of a decoder's soft mask. Decoder
heads would make topology an emergent property of SAM. Costs to manage, from SimpleSeg: sequence length (keep about
100-250 tokens per object), vertex order (canonical start and clockwise order are essential), and thin or concave parts.

---

## 2. Best small VLMs to RL-train on one RTX 5090 (32 GB)

`huggingface-best` run: the only official HF leaderboard tied to grounding is **ScreenSpot-Pro** (GUI point grounding).
There is no official RefCOCO or segmentation leaderboard on the Hub, so model-card numbers fill the gap. Parameter counts
come from HF `safetensors.total`.

| Model (HF) | Params | Released | License | RefCOCO-avg (source) | Point / GUI grounding | Notes for our RL stack |
|---|---|---|---|---|---|---|
| [Qwen/Qwen3-VL-2B-Instruct](https://huggingface.co/Qwen/Qwen3-VL-2B-Instruct) (current) | 2.13B | Oct 2025 | Apache-2.0 | 84.8 (Qwen card); 84.2 (EGM) | ScreenSpot-Pro 48.5 (Qwen card) | Everything is built for it |
| ⭐ [Qwen/Qwen3-VL-4B-Instruct](https://huggingface.co/Qwen/Qwen3-VL-4B-Instruct) | ~4.4B | Oct 2025 | Apache-2.0 | **88.2** (Qwen3.5-2B card); 87.3 (EGM) | ScreenSpot-Pro 59.5 | **Drop-in**: same arch, same 0-1000 grid, same tokenizer. Fits 32 GB with LoRA + colocated vLLM |
| [Qwen/Qwen3.5-2B](https://huggingface.co/Qwen/Qwen3.5-2B) | 2.27B | Feb 2026 | Apache-2.0 | 84.8 think / 84.3 (Qwen card); 78.3 (LFM card) | ScreenSpot-Pro 54.5 | Hybrid: 3 Gated-DeltaNet + 1 gated-attention per block |
| [Qwen/Qwen3.5-4B](https://huggingface.co/Qwen/Qwen3.5-4B) | 4.66B | Feb 2026 | Apache-2.0 | 88.1 (Qwen card); 87.4 (LFM card) | ScreenSpot-Pro 60.3 | Unsloth: bf16 LoRA about 10 GB; needs vLLM ≥ 0.17, flash-linear-attention, causal-conv1d; known LoRA packing (QKVZ) and text-only vLLM issues (vllm#36275, ms-swift#9800) |
| [Qwen/Qwen3.5-9B](https://huggingface.co/Qwen/Qwen3.5-9B) | 9.65B | Feb 2026 | Apache-2.0 | 89.7 | ScreenSpot-Pro 65.2 | bf16 LoRA about 22 GB: tight with colocated vLLM |
| [Qwen/Qwen3.5-0.8B](https://huggingface.co/Qwen/Qwen3.5-0.8B) | 0.87B | Feb 2026 | Apache-2.0 | 79.3 / 77.8 | ScreenSpot-Pro 46.5 | Fast ablations |
| [LiquidAI/LFM3-VL-2B](https://huggingface.co/LiquidAI/LFM3-VL-2B) | 1.92B | Oct 2026 | "other" (LFM license) | **90.5** think / 87.4 (vendor) | Point-Bench 73.8, ScreenSpot-Pro 55.5 (vendor) | Best vendor-reported small grounding. Conv+attention hybrid LM and new vision encoder: new stack, vendor-only numbers |
| [openbmb/MiniCPM-V-4.7-1B](https://huggingface.co/openbmb/MiniCPM-V-4.7-1B) | 1.55B | Oct 2026 | Apache-2.0 | 86.3 (LFM card) | ScreenSpot-Pro 4.9 | Good REC, weak GUI pointing |
| [OpenGVLab/InternVL3_5-4B](https://huggingface.co/OpenGVLab/InternVL3_5-4B) | ~4.7B | Aug 2025 | Apache-2.0 | 88.8 (LFM card) | ScreenSpot-Pro 21.5 | Mature; different coordinate conventions |
| [allenai/Molmo2-4B](https://huggingface.co/allenai/Molmo2-4B) / 8B | 4.85B | Dec 2025 | Apache-2.0 | n/a | PixMo/Point SOTA lineage; MolmoPoint 70.7 PointBench | Points, not polygons; fully open data |
| [google/gemma-4-E2B-it](https://huggingface.co/google/gemma-4-E2B-it) / E4B | ~5B / 8.0B (with PLE) | Mar 2026 | Apache-2.0 | 68.0 (E2B, LFM card); 12B: 82.2 | Point-Bench 14.3 (E2B) | Weak grounding |
| [moondream/moondream3.1-9B-A2B](https://huggingface.co/moondream/moondream3.1-9B-A2B) | 9.3B MoE (2B active) | Jun 2026 | "other" | n/a | native point/detect skills | MoE + custom code: poor fit for TRL/vLLM LoRA |
| Kimi-VL-A3B (16B MoE), GLM-4.6V-Flash (9B), SmolVLM2-2.2B | – | 2025 | mixed | – | – | Kimi-VL is SimpleSeg's backbone but too large for our loop; SmolVLM2 has no grounding focus |

Top of the official ScreenSpot-Pro leaderboard among small models (HF API, 2026-11-08): UI-Venus-1.5-8B 68.4,
UI-Venus-1.5-2B 57.7, UI-TARS-1.5-7B 61.6, GTA1-7B 55.5, Qwen2.5-VL-3B 16.1. These are GUI-specialised fine-tunes.
They show that RL-tuned 2B models can point precisely, but they are not polygon models.

**Verification of the Qwen3.5 numbers (owner request).** I parsed the vision table on the Qwen3.5-2B model card. Its
columns are `Qwen3-VL-4B | Qwen3-VL-2B | Qwen3.5-2B | Qwen3.5-0.8B`. The RefCOCO(avg) row reads
`88.2 | 84.8 | 84.8/84.3 | 79.3/77.8`. So **88.2 belongs to Qwen3-VL-4B, not Qwen3.5-2B**. Qwen3.5-2B (84.8 thinking /
84.3 non-thinking) **ties the current Qwen3-VL-2B (84.8)** on referring grounding. Other rows, same column order:
- CountBench 89.4 | 84.1 | 91.4/86.8 | 77.0/68.6
- ODinW13 39.4 | 36.0 | 35.9/40.5 | 31.6/33.2
- ScreenSpot-Pro 59.5 | 48.5 | 54.5 | 46.5
- MMMU 70.8 | 61.4 | 64.2/64.2 | 49/47.4

Qwen3.5-2B gains over Qwen3-VL-2B on GUI pointing (+6.0 ScreenSpot-Pro), counting and reasoning, but **not on box
grounding**. On LFM's card (third-party eval), Qwen3.5-2B scores only 78.3 RefCOCO-avg. For polygon delineation, box/point
grounding precision is the most relevant proxy, so Qwen3.5-2B is **not a clear upgrade** over Qwen3-VL-2B. Its
architecture change also has a real engineering cost.

### 2.1 ADR: backbone for the next RL run

**Status:** Proposed. **Context:** one RTX 5090 (32 GB). The custom GRPO trainer is tuned for Qwen3-VL: per-chunk
trimmed forwards, RoPE cache, shared-prefix forward, vLLM colocate with sleep. Polygon outputs run 100-250 tokens per
object. The goal is the best IoU/topology at the smallest size that passes.

| Option | Grounding proxy | Engineering cost | Risk |
|---|---|---|---|
| A. Stay Qwen3-VL-2B | RefCOCO 84.8 | none | ceiling |
| **B. Qwen3-VL-4B** | **88.2 (+3.4)** | config change only; about 2× step time | memory with colocated vLLM (should fit in bf16 LoRA) |
| C. Qwen3.5-2B | 84.8 / 84.3 (= A) | high: Gated-DeltaNet breaks RoPE-cache, chunk-trim and packed-prefix code; needs vLLM ≥ 0.17, flash-linear-attention; QKVZ LoRA packing issues | no grounding gain to pay for it |
| D. Qwen3.5-4B | 88.1 | high (as C) | as C |
| E. LFM3-VL-2B | 90.5 (vendor only) | high: new LM + vision encoder; non-Apache license | vendor numbers inconsistent with Qwen's |

**Decision:** B as the scaling point. A stays the headline small model. C/D only as an "architecture generality"
ablation after a 20-step smoke test on sm_120. **Consequences:** B doubles the per-step cost but gives a clean
2B→4B scaling plot with identical code. C/D would need trainer forks.

**Recommendation.**
1. Keep Qwen3-VL-2B as the headline "2B" result for continuity.
2. Run **Qwen3-VL-4B** as the scaling point. It needs no code change and the base gains about +3.4 RefCOCO-avg.
3. Treat **Qwen3.5-4B** as a "newest architecture" ablation only after vLLM/TRL support is checked on sm_120. The
   Gated-DeltaNet layers invalidate chunk-trimmed forwards, the RoPE cache and the packed-prefix tricks that assume
   standard attention.
4. Do not move to LFM3-VL without an independent RefCOCO check. Its card shows Qwen3.5-2B at 78.3, while Qwen's own
   card says 84.8. Vendor numbers do not match.

---

## 3. Medical segmentation: what to beat

### 3.1 Foundation and specialist models
- **MedSAM** (Nature Communications 2024; from memory, not re-checked: 1.57M image-mask pairs, box-prompted).
- **MedSAM2** (arXiv 2504.03600, Apr 2025; https://github.com/bowang-lab/MedSAM2). SAM2.1-Tiny fine-tuned on 455k 3D
  pairs. Box on the middle slice, then propagation. Median DSC: CT lesions **86.68 %** (IQR 74.32-91.14), CT organs
  88.84, MRI lesions 88.37, PET lesions 87.22. EfficientMedSAM-Top1 got 77.95 on CT lesions. Its annotation pipeline cut
  CT lesion time from 525.9 s to 74.3 s. Note: box-prompted, so it is an "oracle box" upper reference, not a fair
  text-only comparison.
- **BiomedParse** (Nature Methods 2025, arXiv 2405.12971). Text-prompted, 9 modalities. It beats MedSAM even when MedSAM
  gets oracle boxes. **v2** (BoltzFormer, end-to-end 3D, built-in object-existence check):
  https://github.com/microsoft/BiomedParse/tree/v2. Pipelines that chain a medical VLM to BiomedParse reach only Dice
  30.3 (MedGemma), 44.8 (Lingshu) and 52.7 (HuluMed) on grounded disease recognition, against 66.7 for UniBiomed
  (Nature Communications 2026 supplement).
- **SAM-Med3D** (arXiv 2310.15161; from memory, not re-checked). 3D promptable.
- **nnU-Net** (Nature Methods 2021). LIDC numbers depend heavily on protocol:
  - VDD paper (arXiv 2603.04024), patient-level 8:1:1: nnU-Net **3D Dice 0.767**, **2D Dice 0.638**. Prob U-Net 3D 0.722.
  - SOTA2 leaderboard (2023): nnU-Net 74.86 Dice; nnU-Net + multi-scale + SaTTCA 77.4.
  - A modified 3D nnU-Net with focal loss: 83.4 Dice, 5-fold (ICECET 2023).
  - Our 2D single-slice setting is closest to "nnU-Net 2D ≈ 0.64 Dice". Our **Dice 0.701 against the majority**
    looks competitive, but only a same-split run settles it. **Action: train nnU-Net 2D on our exact crops/splits,
    against the majority and against each rater. This is the ceiling-to-beat number reviewers will ask for.**

### 3.2 Multi-rater / ambiguous segmentation on LIDC (128×128 crops, 4 raters)
GED uses d = 1 − IoU and is lower-is-better. Protocols differ, so compare within a paper.

| Method | GED | Other | Source |
|---|---|---|---|
| Prob U-Net (NeurIPS 2018) | 0.252 (16 samples, from memory, not re-checked); 0.31 (SSB table); 0.2168 (D-Persona protocol, 50 samples) | | Kohl et al.; MICCAI 2025 SSB; 2605.08210 |
| PHiSeg (MICCAI 2019) | 0.270 | | SSB table |
| CIMD (2023) | 0.318-0.321 | | SSB table |
| CCDM | 0.264 | | SSB table |
| SSB / SSB++ (MICCAI 2025, Schrödinger bridge) | 0.245 / **0.208** | CI 0.958 | papers.miccai.org/miccai-2025/paper/4859 |
| D-Persona (CVPR 2024) | 0.1358 (50 samples) | Dice_max 91.37 | 2605.08210 table |
| ProSeg (2512.00748) | 0.1077 (prior) / 0.1152 | D_mean 90.25 % personalised | arXiv |
| Harmonizer (2605.08210, May 2026) | **0.1048** | Dice_soft 91.81 | arXiv |
| Diffusion ambiguous seg (2504.05977) | 0.233 (16 preds) | dice 0.486 | arXiv (this protocol counts empty masks, so Dice is low) |

The table shows that "Dice on LIDC" ranges from 0.49 to 0.92 by protocol: empty-mask handling, crop, 2D or 3D, and
majority against best-rater. **We must state our protocol and report GED/D_max/HM-IoU computed from our K = 8 samples
against the 4 raters.** GPRO already draws these 8 samples, so the cost is zero.

**Debug note: why our 0.645 cannot be compared with published numbers (engineering:debug).**
- *Expected:* compare with the literature. *Actual:* published LIDC Dice runs from 0.486 to 0.918.
- *Isolate:* our IoU_tc counts pixels in the doctors' disagreement band as neither right nor wrong
  (`gpu_scoring.md`: union / core / vote planes, band depth). Literature uses Dice/IoU against one rater, the majority,
  or the best rater (D_max), on 128×128 crops, with or without empty masks, in 2D or 3D.
- *Root cause:* different metric and different protocol.
- *Fix:* lead with Dice against the strict majority (ours: 0.701), plus per-rater D_max, GED and HM-IoU from the K = 8
  samples. Name the protocol explicitly: Kohl-style 128×128 central crops against our shifted crops. Run baselines on
  our split.
- *Prevention:* one eval script that emits both our metrics and the literature's.

### 3.3 RIGA (6 ophthalmologists)
- MRNet (CVPR 2021) on Magrabia test, soft Dice: disc 97.55, **cup 87.20**. Inter-rater evidence: a U-Net trained on one
  rater scores cup Dice 0.79-0.87 against its own rater and lower against others.
- Best cup on the SOTA2 leaderboard (mean voting): TAB 88.22, AVAP 87.6. Self-calibrated MrPrism: 88.24 (MV).
- MPNN (2311.18496): IoU_cup 78.11 / Dice_cup 87.25 against majority vote.
- Our README measures the RIGA cup ambiguous band at a median of 56 % of the union. Under this ambiguity, a set-valued
  reward is principled. A single-label Dice of about 87 % is the ceiling for single-rater-supervised models.

### 3.4 Medical VLMs and medical RL segmentation (2025-2026)
- **Off-the-shelf grounding is near random**. SLAKE grounding split, box IoU (2604.27720, Jul 2026): GPT-5.5 0.24,
  Gemini-3-Flash 0.20, Lingshu-7B 0.07, **MedGemma-4B 0.05**, Grounding DINO 0.08, centre box 0.10. Fine-tuning with a
  grounding mix raises it to 0.39.
- RL medical segmentation, all through SAM-family decoders and single-mask rewards:
  - **MedReasoner** (2508.08177): GRPO box + 2 points → SAM; U-MRG-14K.
  - **CORE-Seg** (2603.05911): SFT→GRPO, a "dual-granularity mask reward" (GIoU when Dice <0.05, else 1 + 2·Dice), and a
    penalty for disconnected regions. Mean Dice 37.06 on ComLesion-14K.
  - **IBISAgent** (CVPR 2026): multi-turn clicks → MedSAM2, step-wise IoU-improvement rewards.
  - **MedSAM-Agent** (2602.03320, Stanford): multi-turn agentic RL.
  - **MedVol-R1** (MICCAI 2026): key slice + boxes → MedSAM2.
  - **PathChat-SegR1** (ICLR 2026): SO-GRPO.
  - **MedMO** (CVPR 2026): box rewards with Hungarian matching.
- **None outputs polygons directly, none uses multi-rater sets, none certifies topology.**

---

## 4. Topology- and boundary-aware learning, and how to use it as RL reward

| Method | Venue | Idea | As an RL reward |
|---|---|---|---|
| TopoLoss (Hu et al.) | NeurIPS 2019 | Persistent-homology matching of critical points | Not needed: we can compute exact Betti numbers |
| Clough et al. | TPAMI 2022 (PMC9721526) | PH loss with a prior Betti number, no GT needed | Prior-only topology: usable *at test time* as a verifier (e.g. nodule b0 = 1) |
| clDice (Shit et al.) | CVPR 2021 | Skeleton-overlap Dice for tubular shapes | Topo-R1 uses it (weight 0.05) |
| Betti matching (Stucki et al.) | ICML 2023 | Induced matching of barcodes, spatially correct | **Best upgrade**: it penalises "right count, wrong place" |
| Topograph | ICLR 2025 | Graph-based, strict guarantees, fast | Alternative to Betti matching |
| Boundary IoU (Cheng et al.) | CVPR 2021 | IoU within a distance d of the contours | Complements our f5 |
| HD loss (Karimi & Salcudean) | IEEE TMI 2020 | Hausdorff surrogate | Use capped HD95 as a penalty: it catches stray vertices |
| Boundary loss (Kervadec) | MIDL 2019 | Distance-map integral | Our band depth D is similar |

**Using these as RL rewards: lessons from the literature**
- With RL the reward needs no differentiability, so exact discrete metrics can replace their soft surrogates directly.
  This is our structural advantage over the loss-function literature.
- Reward hacking is real. S-measure alone gave black masks (Seg-R1). Box-proxy rewards drift from mask quality (StAR on
  VisionReasoner). A hard length penalty hurt by about 10 gIoU (SimpleSeg). Rule: **gate topology and boundary terms on
  an overlap threshold**, as Topo-R1 does with its type/size gating and Hi-GAR with its validity gate. Otherwise "a
  legal blob with b0 = 1" earns reward far from the target.
- Tiered or graded IoU (SAM-R1, StAR) gives denser signal than a binary threshold. Our continuous IoU_tc is already
  graded.
- Proposed reward stack:
  `gate(format, legal) × [IoU_set + λ_b·BoundaryIoU_set + λ_t·1{(b0,b1) ∈ rater topology set} − λ_h·min(HD95,cap)]`.
  The topology term should become **Betti matching against the nearest rater** once IoU passes about 0.3.

---

## 5. Test-time scaling with verifiers, and visual chain of thought

- **Parallel sampling + selection works.**
  - **EGM** (NVlabs, arXiv 2601.13633): scaling generated tokens in small models. Qwen3-VL-2B RefCOCO-avg 83.6 → 89.6.
    EGM-8B gets 91.4 at 737 ms, against 90.5 at 4,320 ms for the 235B model.
  - **StAR** (2603.14382): IoU-cluster majority voting over masks (threshold 0.85).
  - **Reason Twice** (2606.09303): an MLLM mask evaluator predicts IoU for each candidate and reranks.
  - **SegAgent** (CVPR 2025): a PRM that predicts the current mIoU for each step.
- **Sequential visual self-correction mostly fails for grounding.** **IVT** (arXiv 2606.13156, v2 Jul 2026, "Self-Correction
  Mirage"): predict, render a red box, refine. It gains +2.4 pp Acc@0.5. But painting the *ground-truth* box on the
  overlay barely helps (IoU 0.67 → 0.68), and painting a wrong box degrades the answer by 5.6-8.1 pp. The model does not
  read the rendered evidence. Training-free SoM loops (Seg-Agent 2605.12953) help selection among candidates more than
  refinement.
- **"Thinking" hurts plain perception** (Perception-R1: RefCOCO 89.1 → 75.1 with thinking; Qwen3-VL Thinking ≤ Instruct).
  It helps only when reasoning is trained for the task (EGM, Seg-Zero CoT on OOD: 51.3 → 53.8).
- **What fits this project.**
  1. Self-consistency medoid: from the K = 8 samples, pick the polygon with the highest mean IoU to the others. No ground
     truth needed; one GPU checker call (pairwise IoU).
  2. Checker gates that need no ground truth: illegal rings, wrong prior topology (nodule b0 ≠ 1 unless the image
     supports it), area outliers.
  3. A learned critic: distil checker scores (IoU_tc, f5) into a small head or the VLM itself ("predict your IoU").
     Then best-of-N. The exact checker supplies unlimited labels for the critic.
  4. **Distributional output**: report the 8 samples as a predictive distribution and score GED / D_max against the
     raters. This turns sampling temperature into a calibrated ambiguity model, directly comparable with Prob U-Net,
     PHiSeg and D-Persona. A VLM evaluated as an ambiguous-segmentation model would be new.

---

## 6. Novelty gap: what the project can claim (searched, no prior found)

| Claim | Nearest prior art | Why ours is different |
|---|---|---|
| **Set-valued multi-rater RL reward** (IoU against the band between core and union; topology counts as correct if any rater drew it) | RPDL loss (reward-penalty Dice map from multiple experts, supervised); MRNet / D-Persona / ProSeg (supervised); SA-Co/Gold oracle-best-of-3 *evaluation* | First use of an expert *set* as the verifiable reward in policy optimisation; no consensus collapse |
| **Exact topology (b0, b1, holes) as a reward for a delineating VLM** | Topo-R1 (clDice, anomaly *detection*); CORE-Seg (disconnection heuristic); TopoLoss/Betti-matching (supervised losses) | Exact integer Betti numbers from certified polygons; holes supported (SimpleSeg cannot); measured 6.1 % LIDC topology disagreement |
| **Decoder-free polygon VLM on medical images** | All medical RL work routes through SAM/MedSAM2; SimpleSeg works on natural images | Verifiable token-level geometry for clinical audit; 2B model |
| **Certified, bit-exact GPU geometry verifier as a reward service** | SAM-in-the-loop rewards (SAM-R1, GenSeg-R1, StAR): approximate and expensive; Shapely-style CPU checks | One source for CPU/GPU, bit-identical, masks never leave the device. It can be packaged as an RLVR verifier library for geometry, like math verifiers for math RL |
| **VLM evaluated as an ambiguity model (GED from K samples)** | Prob U-Net / PHiSeg / diffusion samplers | Language-model sampling as the ambiguity source, with a reward that aims at the expert set |

### 6.1 Verifier as a reusable reward service: design sketch (engineering:system-design)

- **Requirements.** Functional: polygons (rings, holes, multiple objects) plus an expert-set handle in; legality
  defects, IoU against the set or band, boundary F / Boundary IoU, HD95, Betti signature and Betti matching out. Non-functional:
  bit-exact between CPU and GPU, no host round-trip of masks, at least 256 completions per step under 1 s on one GPU,
  deterministic.
- **Components.**
  1. Expert-set store: device-resident rater bit planes, union / core / vote planes, distance maps. Built once per item.
  2. Batch scorer: packed polygons → warp-per-row raster → reductions.
  3. Thin Python API: `score_many`.
  4. Optional gRPC/HTTP front for external RL frameworks (verl, TRL, OpenRLHF), batching requests per step.
- **Trade-offs.**
  - In-process (current) is lowest latency, but tied to one trainer. A service is reusable across frameworks and
    models, but adds about 1 ms of IPC per batch. That is negligible against 23-36 s steps.
  - Exact integer geometry costs more code than a SAM-in-the-loop reward, but has no model error and is about 1000×
    cheaper than a SAM2 forward per rollout.
- **What to revisit.** 3D (voxel / mesh) topology. Multi-class shared boundaries. A public benchmark harness
  ("GeoVerify-Med") so other polygon-VLM papers can report set-valued scores.

**Must-do experiments for a top-venue paper (MICCAI / MIDL / CVPR workshop to main track):**
1. Same-split baselines:
   - nnU-Net 2D (majority and per-rater)
   - MedSAM2 with oracle box; SAM2.1-small with the decoder fine-tuned on box prompts (vision-trainer recipe)
   - SAM 3 ("lung nodule", "optic cup")
   - NVIDIA `nv-segment-ct` (VISTA-3D-style), optional
   - Qwen3-VL zero-shot
   - **released SimpleSeg checkpoints (`sthui/SimpleSeg-Qwen2.5-VL`), zero-shot and fine-tuned with its own
     IoU+centroid reward**, as the closest-method baseline.
2. Reward ablations: single-rater against majority against **set-valued**. Report per-rater D_max, GED, and topology
   error rate (b0/b1 mismatch against every rater). Topology on/off, especially on the 65 hole-bearing nodules.
3. Standard metrics: Dice/IoU against majority, Boundary IoU, HD95, GED, HM-IoU, plus our IoU_tc and f5 as secondary.
4. Scaling: 0.8B → 2B → 4B (Qwen3-VL-4B), with and without RL. Test-time selection (medoid, critic) against K.
5. RIGA: report cup Dice against MV in the MRNet protocol (Magrabia held out) for direct comparison with 87.20 / 88.22.

**Risks.** SimpleSeg (and its ICLR 2026 outcome) sets the "first RL for decoder-free polygon segmentation" priority. Do
not claim that. Claim the set-valued, topology-certified, medical, and verifier contributions. Literature Dice on LIDC
is protocol-chaotic, so lead with GED and same-split baselines, not cross-paper Dice.

---

## Sources (primary URLs)

- SimpleSeg: https://arxiv.org/abs/2601.19228, https://simpleseg.github.io/, https://openreview.net/forum?id=CcQFKeUK5K
- PolyFormer: https://arxiv.org/abs/2302.07387, https://github.com/amazon-science/polygon-transformer
- VistaLLM: https://shramanpramanick.github.io/VistaLLM/
- Text4Seg: https://arxiv.org/abs/2410.09855
- PaliGemma 2: https://arxiv.org/abs/2412.03555; evaluator: https://github.com/google-research/big_vision/blob/main/big_vision/evaluators/proj/paligemma/transfers/segmentation.py
- Gemini 2.5 segmentation: https://simonwillison.net/2025/Apr/18/gemini-image-segmentation/
- UniPixel: https://github.com/PolyU-ChenLab/UniPixel (NeurIPS 2025 paper PDF on proceedings.neurips.cc)
- Sa2VA: https://arxiv.org/abs/2501.04001
- Seg-Zero: https://arxiv.org/abs/2503.06520; Seg-R1: https://arxiv.org/abs/2506.22624; SAM-R1: https://arxiv.org/abs/2505.22596; StAR: https://arxiv.org/abs/2603.14382; GenSeg-R1: https://arxiv.org/abs/2602.09701; LENS: https://arxiv.org/abs/2508.14153
- Visual-RFT: https://arxiv.org/abs/2503.01785; VLM-R1: https://arxiv.org/abs/2504.07615; Perception-R1: https://arxiv.org/abs/2504.07954
- Qwen3-VL report: https://arxiv.org/abs/2511.21631; EGM: https://arxiv.org/abs/2601.13633, https://github.com/NVlabs/EGM
- SAM 3: https://arxiv.org/abs/2511.16719, https://github.com/facebookresearch/sam3
- Hi-Token: https://arxiv.org/abs/2608.03471; GETok (CVPR 2026): https://openaccess.thecvf.com/content/CVPR2026/papers/Ren_Grounding_Everything_in_Tokens_for_Multimodal_Large_Language_Models_CVPR_2026_paper.pdf; MolmoPoint: https://arxiv.org/abs/2603.28069; Ruler: https://arxiv.org/abs/2510.03230
- Model cards: https://huggingface.co/Qwen/Qwen3.5-2B, https://huggingface.co/Qwen/Qwen3.5-4B, https://huggingface.co/LiquidAI/LFM3-VL-2B; HF ScreenSpot-Pro leaderboard: https://huggingface.co/datasets/likaixin/ScreenSpot-Pro
- Qwen3.5 RL support: https://www.unsloth.ai/docs/models/qwen3.5/fine-tune, https://github.com/vllm-project/vllm/issues/36275, https://github.com/modelscope/ms-swift/pull/9800
- MedSAM2: https://arxiv.org/abs/2504.03600; BiomedParse: https://arxiv.org/abs/2405.12971, v2 https://github.com/microsoft/BiomedParse/blob/v2/README.md
- LIDC: VDD https://arxiv.org/abs/2603.04024; SSB https://papers.miccai.org/miccai-2025/paper/4859_paper.pdf; Harmonizer https://arxiv.org/abs/2605.08210; ProSeg https://arxiv.org/abs/2512.00748; Gen. Prob U-Net https://arxiv.org/abs/2207.12872; leaderboards https://www.sota2.com/research/sota/medical-image-segmentation-on-lidc-idri-test, https://www.sota2.com/research/sota/lung-nodule-segmentation-on-lidc-idri-test
- RIGA: MRNet https://openaccess.thecvf.com/content/CVPR2021/papers/Ji_Learning_Calibrated_Medical_Image_Segmentation_via_Multi-Rater_Agreement_Modeling_CVPR_2021_paper.pdf; MrPrism https://arxiv.org/abs/2206.05092; MPNN https://arxiv.org/abs/2311.18496; leaderboard https://www.sota2.com/research/sota/optic-disc-and-optic-cup-segmentation-on-riga
- Medical grounding: https://arxiv.org/abs/2604.27720; MedReasoner https://arxiv.org/abs/2508.08177; CORE-Seg https://arxiv.org/abs/2603.05911; IBISAgent (CVPR 2026); MedSAM-Agent https://arxiv.org/abs/2602.03320; MedVol-R1 https://papers.miccai.org/miccai-2026/paper/2593_paper.pdf; PathChat-SegR1 (ICLR 2026)
- Topology: Topo-R1 https://arxiv.org/abs/2603.13054; Betti matching https://proceedings.mlr.press/v202/stucki23a.html; TopoLoss https://proceedings.neurips.cc/paper_files/paper/2019/file/2d95666e2649fcfc6e3af75e09f5adb9-Paper.pdf; Clough https://pmc.ncbi.nlm.nih.gov/articles/PMC9721526/; Topograph (ICLR 2025)
- Multi-rater losses: RPDL https://lacuna.tiptreesystems.com/work/learning-non-unique-segmentation-with-reward-penalty-dice; ordinal multi-rater calibration https://arxiv.org/abs/2605.02437
- TTS / self-correction: IVT https://arxiv.org/abs/2606.13156; Seg-Agent https://arxiv.org/abs/2605.12953; Reason Twice https://arxiv.org/abs/2606.09303
