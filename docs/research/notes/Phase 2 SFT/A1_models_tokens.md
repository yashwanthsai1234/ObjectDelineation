# A1: small VLMs and polygon output representation for Phase 2 SFT

Task A1 of `experiments/phase2/RESEARCH.md`. Written 2026-10-19.

The question: which small VLMs, and which polygon output representation, give the best SFT for legal, precise polygons? The polygons have outer rings, holes and several parts. The data are RIGA disc and cup crops and LIDC nodule crops. A second question is what the newest systems do.

This note builds on Phase 1 rather than repeating it:
- `docs/research/notes/Phase 1 stack and algorithms/output_format.md` (D10 evidence: SimpleSeg, Groundbench, PolyFormer, VistaLLM, VectorLLM, FloorplanVLM, Rex-Omni, Hi-Token);
- `training_stack.md` (first model table).

Evidence labels:
- **[local]** measured in this repo today. Scripts and JSON are in `experiments/phase2/lab/A1/`. Every script that reads `data/` runs with `.venv/bin/python -I` from the repo root.
- **[paper] / [card] / [code]** read in the primary source this session.
- **[agent]** gathered by a research sub-agent from primary pages. I spot-checked the load-bearing ones and say which.
- **[unverified]** not confirmed.

---

## 0. Bottom line

1. **Keep the D10 text format** (compact JSON, digits 0–1000, clockwise outers and counter-clockwise holes, top-most-then-left-most start) for every model. The field has converged on it:
   - Qwen3-VL, Qwen3.5, InternVL3/3.5, Gemma 4 and LFM2.5-VL-3B all emit 0–1000 integers natively.
   - In SFT, text coordinates match special tokens (GETok Table 9; Shikra).
   - A canonical start vertex is worth about 40 IoU points (LPM Table 3).
2. **Change the model ladder.** Add **Qwen3.5-0.8B** (Apache-2.0, native 0–1000 grounding, the only small family in TRL's GRPO VLM CI test). Demote **SmolVLM-256M** to an optional control: no prior, excluded from TRL CI over bf16, and 1 visual token per 64 px tile patch. The ladder to screen becomes LFM2.5-VL-450M → Qwen3.5-0.8B → Qwen3-VL-2B (or Qwen3.5-2B), with InternVL3.5-1B as an optional fourth.
3. **The biggest lever on GRPO cost is output tokens per vertex, and it differs 1.5× between tokenizers.** Measured on all 30,843 targets:
   - LFM: 5.7 / 6.1 tokens per vertex (RIGA / LIDC), from 1–3-digit number tokens;
   - Qwen, SmolVLM, InternVL and Gemma: 8.6 / 8.9, from single-digit tokens.
   - Added coordinate tokens would give 2.5–2.8 (3.1–3.4× fewer). This is the only representation change worth a test (H-A1-4).
4. **Reject delta encoding and pixel-unit coordinates** with numbers:
   - Delta costs 25 % fewer tokens, but one wrong digit spreads to every later vertex. A ±1 slip in delta-pixel units costs 13 IoU points on LIDC.
   - Pixel units save 23 % on LIDC, but turn a last-digit slip into a 1 px error (−2.8 IoU points, against 0 for 0–1000). They also cannot express the target builder's bevel vertices (0.4 % of targets).
5. **A0 vs A0′ is a LIDC-only question, with little data and no val power:**
   - RIGA is identical under both arms, because the rim is derived from disc and cup.
   - LIDC train: 2.2 % of targets are not simple; 5.0 % of items have at least one non-simple rater; only 23 items are unanimously non-simple.
   - LIDC val: 48 items with a non-simple rater and **0 unanimous**, so the entry gate "topology-any on unanimous items" cannot separate A0 from A0′.
   - We need a topology probe (H-A1-7) and oversampling (H-A1-6).
   - No VLM paper trains holes in text polygons. Even explicit-hole specialists reach only Hole-F1 ≤ 0.53 (PolyTopoBench).

---

## 1. Local measurements [local]

### 1.1 Tokens per vertex

The table covers all 30,843 RIGA and LIDC targets (train + val) and seven candidate encodings, for nine tokenizers. Script: `measure_enc.py`; results: `enc_tokens.json`.

Each cell is tokens per vertex, RIGA / LIDC, including all JSON skeleton tokens. The encodings:
- **json1k**: today's D10 format.
- **json_px**: integers in frame pixels (0..S).
- **flat1k**: rings as `[x,y,x,y,...]`.
- **space1k**: rings as `[x y x y]`.
- **delta**: first vertex absolute, then signed differences.
- **loc1k**: 2,002 added tokens `<x0..1000><y0..1000>`, written `[<x500><y440>...]`.

| tokenizer (models) | json1k | json_px | flat1k | delta1k | delta_px | loc1k |
|---|---|---|---|---|---|---|
| Qwen3-VL / Qwen3.5 / InternVL3.5 / FastVLM (single digits) | 8.55 / 8.85 | 8.37 / 6.85 | 8.45 / 8.71 | 6.35 / 6.73 | 5.95 / 5.01 | 2.50 / 2.78 |
| SmolVLM (SmolLM2, single digits) | 8.74 / 9.06 | 8.56 / 7.06 | 8.74 / 9.06 | 7.01 / 7.34 | 6.62 / 5.62 | 2.79 / 3.13 |
| Gemma 3n / MedGemma (single digits) | 8.55 / 8.85 | 8.37 / 6.85 | 8.45 / 8.71 | 6.83 / 7.12 | 6.43 / 5.41 | 2.50 / 2.78 |
| MiniCPM-V-4 | 8.66 / 8.99 | 8.48 / 6.99 | 8.55 / 8.85 | 7.36 / 7.59 | 6.97 / 5.87 | 2.66 / 2.99 |
| **LFM2.5-VL 450M / 1.6B** (`\p{N}{1,3}`) | **5.71 / 6.06** | 5.71 / 6.06 | **4.76 / 5.13** | 5.71 / 6.06 | 5.71 / 6.06 | 2.82 / 3.21 |
| Florence-2 / Moondream2 (GPT-2 BPE, irregular merges) | 4.69 / 4.92 | 4.55 / 4.92 | 4.59 / 4.78 | 5.03 / 5.31 | 5.03 / 5.31 | 2.50 / 2.85 |
| Phi-4-MM (o200k) | 4.61 / 4.92 | 4.61 / 4.92 | 4.50 / 4.78 | 5.09 / 5.31 | 5.09 / 5.31 | 2.55 / 2.85 |

Median answer tokens, RIGA / LIDC (p95 in brackets):

| | json1k | loc1k |
|---|---|---|
| Qwen | 309 / 100 (381 / 268) | 91 / 33 |
| SmolVLM | 316 / 103 | 102 / 38 |
| LFM | 207 / 70 (252 / 175) | 103 / 39 |

What the table shows:
- **Decoding is lossless everywhere.** `decode(encode(t)) == t` for 100 % of targets, for every tokenizer and encoding.
- **Why "flat" helps only LFM:**
  - Qwen and SmolVLM already merge `],[` into one token, so flat rings save 1 % there.
  - LFM splits `],` and `[` into two tokens, so flat rings save 17 % on LFM.
  - Space-separated rings save nothing on the digit tokenizers.
- **LFM numeric vocabulary.** It has 1,110 numeric tokens: all 1,000 three-digit strings, 100 two-digit and 10 one-digit. Of these, 110 start with `0` (e.g. `042`), a leading-zero failure mode only LFM can produce. `1000` becomes `100|0`. This never matters for targets (next item), but the grammar allows it.
- **Coordinate range.** Every target coordinate lies in **[168, 832]**: no 0, no 1000, always 3 digits. The 25 % frame margin causes this.
  - LIDC frames are 50 px (6,267 train items), 100 px (1,068) or 125 px (6). Their coordinates are multiples of 20 or 10, so the last digit is always `0`. On digit tokenizers that is about one LIDC coordinate token in three, almost free to predict.
  - RIGA frames are 500 px (424 items), 250 px (167) or 200 px (1), with even coordinates.
- **Bevel vertices.** 0.34 % (RIGA) and 0.47 % (LIDC train) of targets contain vertices off the pixel corners: the builder's pinch bevel, one grid unit off a corner. These need the full 0–1000 resolution. A pixel-unit encoding (or a grid of G=500) cannot express them and makes about 0.13 % of LIDC objects illegal (`coarse_grid.py`).

### 1.2 One-digit sensitivity of each encoding

Setup: 2,000 random val targets; one random coordinate of one ring is corrupted, then decoded and rasterized with polycore. Each cell is IoU mean (p5) against the uncorrupted mask, then the share that became illegal. Script: `perturb.py`.

| error | dataset | abs 0–1000 (D10) | abs pixels | delta 0–1000 | delta pixels |
|---|---|---|---|---|---|
| ±1 in the last unit | LIDC | **1.000 (1.000), 0 %** | 0.972 (0.915), 0.2 % | 0.9998, 0 % | 0.871 (0.706), 0.3 % |
| ±1 in the last unit | RIGA | 0.9995, 0 % | 0.9987, 0 % | 0.9955, 0 % | 0.9887, 0 % |
| random digit replaced | LIDC | 0.894 (0.493), 16.9 % | 0.785 (0.288), 32.2 % | 0.872 (0.383), 1.4 % | 0.642 (0.001), 10.9 % |
| random digit replaced | RIGA | 0.949 (0.760), 13.7 % | 0.929 (0.707), 17.9 % | 0.926 (0.706), 1.1 % | 0.888 (0.505), 3.2 % |

Reading:
- On the 0–1000 grid the last digit of an LIDC coordinate is a sub-pixel digit (1 unit = 0.05 px on a 50 px frame), so near-miss slips are free.
- Delta encodings rarely go illegal: a shifted tail stays simple. But they lose more IoU, because the error moves every later vertex. This is the "error accumulation" that STAMP and VectorArk cite against long coordinate chains.

### 1.3 Visual tokens per crop

Counted with each model's processor on blank images; includes the image start/end markers. Script: `vis_tokens.py`.

| processor | 250 px | 500 / 512 px | image px per visual token |
|---|---|---|---|
| Qwen3-VL-2B, Qwen3.5-0.8B/2B, LFM2.5-VL-450M/1.6B | 66 | 258 | 32 × 32 |
| InternVL3.5-1B-HF | 259 | 259 | 28 × 28 on a 448 tile |
| Gemma 3n E2B, MedGemma 4B | 260 | 260 | fixed |
| SmolVLM-256M/500M, default processor | 1,127 | 1,127 | upsampled to 2048 and tiled |
| SmolVLM, `longest_edge` 512 / 1024 / 1536 | 76 / 342 / 673 | same | 64 × 64 (pixel shuffle 4) |

Consequences:
- **250 px RIGA frames get only 64 visual tokens** on Qwen and LFM, i.e. 28 % of RIGA train items. The LIDC upscale to 512 px (`FrameConfig(image_side=512)`) gives a 50 px nodule frame 256 tokens, one per 3.1 frame px.
- RIGA is not upscaled, so a 500 px frame gets one token per 32 frame px and a 250 px frame one per 31 frame px. See H-A1-5.
- Phase 1 already sets SmolVLM's `longest_edge` (`src/objdel/train/configs.py:79`). The default would cost 1,127 tokens.

### 1.4 A0 vs A0′ ceiling on real targets

Script: `a0_ceiling.py`.

| | LIDC train | LIDC val |
|---|---|---|
| targets | 23,469 | 3,444 |
| multi-part | 200 | 6 |
| with holes | 317 | 47 |
| both | 14 | 1 |
| items | 7,341 | 1,084 |
| items with ≥ 1 non-simple rater | 366 (5.0 %) | 48 (4.4 %) |
| items with **all** raters non-simple | 23 | **0** |
| IoU of A0-single (largest outer ring only) on non-simple targets | mean 0.838, p5 0.527 | 0.943 |
| IoU of A0-multi (all outer rings, no holes; SimpleSeg-style) on non-simple targets | mean 0.964, median 0.993 | 0.968 |

- RIGA has 0 holes and 0 multi-part objects, so A0 = A0′ there.
- Most LIDC holes are tiny; multiple parts carry most of the IoU difference.

### 1.5 Coarser grid and parametric ellipse

- **Coarser grids** (`coarse_grid.py`; IoU against the target; legality). Re-quantizing to G=100 is lossless for LIDC, because 50 and 100 px frames divide it, but loses about 2.7 points on RIGA. G=250 costs RIGA 1 point and 1.5 % legality.

| grid | RIGA IoU mean (p5) | RIGA legal | LIDC IoU mean | LIDC legal |
|---|---|---|---|---|
| G = 100 | 0.973 (0.949) | 97.8 % | 1.000 | 99.87 % |
| G = 250 | 0.990 | 98.5 % | 0.9992 | — |
| G = 500 | 0.9999 | 100 % | 1.000 | 99.87 % (bevels) |

- **Best-fit ellipse** for disc and cup (`ellipse.py`, RIGA val). Disc IoU median 0.974 (p5 0.956); cup 0.966 (p5 0.942). The certified targets reach 0.99 or more. An ellipse would cost about 15 tokens, but it cannot carry topology or irregular nodules, so it is a baseline, not a format.

---

## 2. Small VLMs (≤ ~4B), October 2026

Sources: an [agent] survey of model cards, configs, the vLLM and TRL docs, and the TRL test source.

What I verified myself:
- the TRL `test_grpo_trainer.py` model lists;
- the Qwen3.5-0.8B and LFM2.5-VL-450M cards;
- Gemma 4's `box_2d`;
- the configs and tokenizers of Qwen3.5-0.8B/2B, InternVL3.5-1B-HF and LFM2.5-VL-1.6B, which I downloaded.

Visual tokens in §1.3 are [local]. "Digits" is [local] from tokenizer.json.

| model | params (vision + LM) | licence | native coordinates | grounding prior | TRL GRPO VLM CI | vLLM (LoRA) | digits |
|---|---|---|---|---|---|---|---|
| **LFM2.5-VL-450M** | SigLIP2 NaFlex 86M + LFM2.5-350M (conv + attention hybrid) | LFM Open 1.0 (commercial use under $10M revenue) | JSON `bbox`, floats in [0,1] | RefCOCO-M 81.28 [card] | no (SFT tests only; "LFM2.5-VL support" in TRL v1.11) | yes | 1–3-digit chunks |
| LFM2.5-VL-1.6B / 3B | SigLIP2-so400m + 1.2B / 2.6B | same | 3B: **integers 0–1000** | 3B RefCOCO 87.9 vs 57.1 for LFM2-VL-3B [agent] | no | yes | 1–3 |
| **Qwen3.5-0.8B / 2B** | ViT 12×768 / 24×1024, patch 16, merge 2 + Gated-DeltaNet/attention hybrid, 24 layers, vocab 248k (tied embeddings: about 254M of the 0.8B) [local config] | Apache-2.0 | JSON `bbox_2d`, 0–1000 | RefCOCO avg 79.3 / 77.8 (0.8B, thinking / non-thinking), 2B 84.8 / 84.3 [card] | **yes** (`tiny-Qwen3_5…-NoThink`) [code] | yes | single |
| **Qwen3-VL-2B** / 4B | SigLIP2 + 1.7B / 4B | Apache-2.0 | `bbox_2d`, `point_2d`, 0–1000 | RefCOCO avg 85.6 (2B Instruct) [agent]; our zero-shot box IoU 0.76 (Phase 1 D12) | no (SFT only; but ran our K12 smoke) | yes | single |
| InternVL3.5-1B / 2B -HF | InternViT-300M + Qwen3 0.8B / 2B | Apache-2.0 | `<box>[[x1,y1,x2,y2]]</box>`, 0–1000 | RefCOCO val 85.4 / 88.7 [agent: arXiv 2508.18265] | yes (`tiny-InternVL…`) [code] | yes | single |
| Gemma 4 E2B | ~150M ViT + 2.3B effective (5.1B with embeddings) | Apache-2.0 | `box_2d`, 1000×1000 grid [blog] | no published numbers | yes (`tiny-Gemma4…`, transformers ≥ 5.5) [code] | yes | single [unverified for Gemma 4] |
| SmolVLM-256M / 500M | SigLIP-B/16 93M + SmolLM2 135M / 360M | Apache-2.0 | none | none | **excluded** ("seems not to support bf16 properly") [code] | yes, but TRL+vLLM placeholder bug (PR #6718 open) | single |
| Florence-2 B / L | DaViT + BART, 0.23B / 0.77B, encoder-decoder | MIT | `<loc_0..999>` **polygons**, clockwise | RES val 78.0 / 80.5 fine-tuned [agent] | no ("only causal LMs") | plugin only | irregular BPE |
| Ovis2.5-2B | SigLIP2-so400m + Qwen3-1.7B | Apache-2.0 | `<box>` floats in [0,1) | RefCOCO val 92.8 (best ≤ 2B) [agent] | no (remote code) | no LoRA | single |
| PaliGemma 2 3B | So400m + Gemma 2 2B | Gemma ToU | `<loc####>` + `<seg###>` VQ masks | RefCOCO seg mIoU 76.3 at 448 | no | yes | loc tokens |
| Rex-Omni 3B | Qwen2.5-VL-3B | IDEA research licence | `<0>…<999>` tokens; boxes, points, OCR polygons | — | own pipeline | own | special |
| Isaac-0.2-1B | ~0.6B-class LM | CC-BY-NC | `<polygon>` on 0–1000 in the SDK; whether the open weights learned it is [unverified] | 70–77 (0.1) | no | yes | single |
| Moondream 3 segment | 9B MoE (2B active) | BSL | SVG path, refined by a head | RefCOCO val 80.2 cIoU ([arXiv 2604.02593](https://arxiv.org/html/2604.02593)) | no | no LoRA | — |
| MedGemma 4B / 1.5 | SigLIP 417M + Gemma 3 4B | HAI-DEF | 1.5 adds chest X-ray boxes (IoU 38.0) | no polygons | Gemma3 tested | yes | single |

Facts that matter for us:
- **No model of 4B or less emits text polygons with published referring-segmentation numbers, in a TRL- and vLLM-friendly form.** Florence-2 (polygons) is encoder-decoder. Rex-Omni and Isaac have restrictive licences. Polygons will be taught by SFT in every case.
- **Medical variants of 4B or less** (MedGemma, MediX-R1-2B, Hulu-Med-4B, MedVLM-R1, Fundus-R1-3B) do VQA or boxes, never polygons. MedMO (Qwen3-VL-based, 4B/8B) outputs absolute-pixel boxes only, after 25 days on 64 GPUs ([arXiv 2602.06965](https://arxiv.org/html/2602.06965)).
- **A small-VLM medical precedent.** A 2026 GI-endoscopy study fine-tuned Florence-2 (0.2B), Qwen3.5-0.8B and InternVL3.5-1B with LoRA on one RTX 3090 Ti. Florence-2's native location-token polygons gave the best grounding (out-of-distribution pointing-game 0.645) ([arXiv 2607.27122](https://arxiv.org/html/2607.27122)). It is evidence that a native polygon prior helps at under 1B, and that Qwen3.5-0.8B and InternVL3.5-1B are what others pick at this size.
- **SmolVLM-256M can learn a coordinate format, but slowly.** CropVLM: the base "is not originally capable of generating valid bounding boxes". It needed about 3 A100-hours of LoRA SFT (r=128) on about 62k synthetic boxes in percent units, then 24 GPU-hours of GRPO ([arXiv 2511.19820](https://arxiv.org/html/2511.19820)). That is 4 numbers per box, not 10–25 vertices.
- **Qwen3.5 risks:**
  - The Gated-DeltaNet layers need the `fla` / `causal_conv1d` kernels, or they fall back to slow PyTorch.
  - vLLM GDN kernel bug [vllm#34948](https://github.com/vllm-project/vllm/issues/34948).
  - `AsyncGRPOTrainer` rejects linear attention ([trl#7761](https://github.com/huggingface/trl/issues/7761)).
  - Weight sync needs vLLM ≥ 0.27 ([trl#5269](https://github.com/huggingface/trl/issues/5269)).
  - The 248k vocabulary makes the logits tensor about 1.6× Qwen3-VL's 152k (A2: fused or chunked CE).
  - The 2B card warns of thinking loops: use non-thinking mode ([card](https://huggingface.co/Qwen/Qwen3.5-2B)).
- **Open issue for every VLM:** vLLM's encoder cache is not reset after a weight sync, so freeze the ViT or reset the cache ([trl#7198](https://github.com/huggingface/trl/issues/7198)). This is for A5.

---

## 3. Output representations: what the newest systems do

The Phase 1 survey (output_format.md §1–3) still stands. New or corrected findings follow.

**Coordinate tokenization**
- **Hi-Token is weaker evidence than Phase 1 assumed** [paper, verified] ([arXiv 2608.03471](https://arxiv.org/html/2608.03471)):
  - Fixed recipe: flat tokens 58.3 mIoU vs Hi-Token 72.4.
  - Against a flat baseline with tuned recipe: 68.5 vs 72.4 (P@0.95 26.3 vs 31.7). The authors call the tuned comparison "a conservative check rather than a causal comparison".
  - There is **no native-digit-text baseline**.
  - Qwen3-VL-2B P@0.95: 18.8 → 27.8.
  - Hi-GAR GRPO adds +12 mIoU over SFT (72.4 → 84.3). RL matters more than the token scheme.
- **GETok (CVPR 2026)** [paper, verified] ([arXiv 2512.10554](https://arxiv.org/html/2512.10554v1)):
  - SFT with text coordinates vs grid tokens: RefCOCO testA 92.9 vs 93.0, RefCOCOg 87.4 vs 87.6. **No SFT gain from special tokens.**
  - A coarse grid plus offset tokens beats a finer grid (32² + offset 9.2 tokens per mask, gIoU 68.2; 64² 14.6 tokens, 67.1).
  - Qwen2.5-VL-7B; masks come from SAM.
- **Shikra:** plain-text numbers ≥ added vocabulary (RefCOCO val 81.47 vs 81.03) ([ar5iv 2306.15195](https://ar5iv.labs.arxiv.org/html/2306.15195)).
- **New-token initialization dominates.** HiVG (2,384 coordinate tokens, Qwen2.5-VL-3B): SSIM 0.440 with noise init vs 0.720 with a hierarchical mean plus numeric-order projection ([arXiv 2604.05072](https://arxiv.org/html/2604.05072)) [agent].
- **VectorArk** (InternVL2-1B, SVG): reusing the vocabulary "leads to more stable and efficient convergence compared to introducing new tokens"; no numbers ([arXiv 2605.24398](https://arxiv.org/pdf/2605.24398)) [agent].
- **OmniSVG:** one token per (x, y) pair, x·w + y on a 200² canvas (about 40k tokens). Needs 250B training tokens; useless at our scale ([arXiv 2504.06263](https://arxiv.org/html/2504.06263v2)) [agent].
- **Qwen3-VL's official reason for relative 0–1000:** it "improves robustness to variations in image resolution and aspect ratio … simplifying post-processing" ([arXiv 2511.21631](https://arxiv.org/html/2511.21631)) [agent]. **LFM2.5-VL-3B moved from [0,1] floats to 0–1000 integers** ([card](https://huggingface.co/LiquidAI/LFM2.5-VL-3B)) [agent], and Gemma 4 uses a 1000 grid ([blog](https://huggingface.co/blog/gemma4)) [verified]. **The industry has converged on our grid.**

**Distance-aware losses on coordinate tokens** (for A2)
- **Number Token Loss** (ICML 2025): adds MSE or Wasserstein on digit-token values to CE ([arXiv 2411.02083](https://arxiv.org/html/2411.02083v3)) [agent].
  - Math interpolation: 64 → 75 %.
  - Regression RMSE 2.01 → 1.81, the same as a regression head.
  - It is designed for single-digit tokenizers (Qwen, SmolVLM) and helps less with multi-digit ones (LFM).
  - No VLM-coordinate follow-up found.
- **LPM (Amazon, CVPRW 2026)**, DistilGPT2 with 224 index tokens ([PDF](https://cdn.amazon.science/aa/3c/d4ea78084f199eaf08743fc57e50/scipub-approval152134-45892548-rethinking-language-models-for-building-outline-extraction-from-remote-sensing-imagery.pdf)) [agent].
  - Gaussian label smoothing on coordinate tokens leaves IoU flat (83.77 → 83.67) but improves boundary error, PoLiS 1.403 → 0.955 on SpaceNet2.
  - It hurts on WHU at σ=0.5. The right σ depends on label noise.
- **R-VLM's** IoU-weighted CE gave +8.3 on ScreenSpot. Sources disagree on which component gave that gain [agent, conflicting].

**Ordering**
- **LPM Table 3** [paper, verified in the PDF text]. First vertex / building order, IoU:

| first vertex / order | SpaceNet2 | WHU | INRIA |
|---|---|---|---|
| random / random | 42.19 | 46.71 | 70.72 |
| top-left / random | 83.64 | 88.10 | 72.86 |
| top-left / top-left | 83.80 | 88.69 | 75.33 |

  - A canonical start vertex is worth about 40 IoU points; object order is worth 0–2.5.
  - This is the strongest ordering result to date, and it confirms D10.
  - VectorLLM's random start helps only in pretraining.

**Precision per token**
- **Text4Seg++** box-wise 64² "semantic bricks": 150 tokens per mask, and it beats VistaLLM polygons at 7B (79.3 vs 71.2 cIoU) ([arXiv 2509.06321](https://arxiv.org/html/2509.06321)) [agent]. Our polygons cost 70–207 tokens per answer on LFM and 100–309 on Qwen (medians, §1.1), the same range. Only shown at 7B and above.
- **STAMP** (Qwen2-VL-2B): non-autoregressive patch-classification mask tokens, RefCOCO 79.1 cIoU with SAM, 1.3 s per image. It argues sparse polygons are "unreliable due to error accumulation, where a single inaccurate point can corrupt the entire mask" ([arXiv 2512.00395](https://arxiv.org/html/2512.00395)) [paper]. §1.2 quantifies that error for us: a random digit costs 0.894 mean IoU on LIDC with absolute coordinates.
- **Groundbench:** frontier models' exact-N polygons collapse at N=64 (legal coverage 7.5 %) ([arXiv 2609.27821](https://arxiv.org/abs/2609.27821)). Our targets' p95 is 24 vertices (RIGA) and 32 (LIDC), below that regime.
- **RLRF** (SVG GRPO, Qwen2.5-VL-3B/7B): length collapse and tiny-viewBox hacks under rendering rewards ([arXiv 2505.20793](https://arxiv.org/html/2505.20793v2)) [agent]. Our area guard and validity gate already target these.
- No VLM or LLM-segmentation ablation of delta, chain-code or turtle encodings was found [agent; search budget ran out]. §1.2 is our own evidence against delta.

**Topology** (holes and multiple parts)
- **No VLM paper trains or reports holes in text-polygon output.**
  - SimpleSeg lists holes as a failure case.
  - VectorLLM "cannot handle hollow buildings … nor … multiple parts" ([arXiv 2507.04664](https://arxiv.org/html/2507.04664)).
  - Groundbench keeps only the largest part.
- **PolyTopoBench** (11 polygon extractors, specialist models) ([arXiv 2609.32856](https://arxiv.org/html/2609.32856)) [agent]:
  - Building Hole-F1 at most 0.527 (ACPV-Net); Pix2Poly 0.047.
  - Even at IoU ≥ 0.9, 84–100 % of predictions have the wrong topology.
  - Truncating the vertex budget silently drops holes.
- **Keyhole or bridge hole encodings** create self-intersections and break under cropping ([arXiv 2505.15379](https://arxiv.org/html/2505.15379), [arXiv 2603.14764](https://arxiv.org/abs/2603.14764)) [agent].
- Our explicit `holes` list plus legality gate is therefore the right design and a novelty. The expected base rate of hole recovery is low, so A0′ needs oversampling and a probe (§1.4).

---

## 4. Patterns, and what fits ≤ 2B models at about 1 GPU-hour each

1. **Convergence on 0–1000 integer text.** Qwen3-VL/3.5, InternVL, Gemma 4, LFM2.5-VL-3B, Isaac, Florence-2's bins and Groundbench all use it. Models that use it zero-shot are the ones whose priors transfer. D10 matches, so keep it.
2. **Special coordinate tokens win at pretraining scale, not at SFT scale.** The wins come from Rex-Omni (22M samples), VectorLLM (pretraining), OmniSVG (250B tokens) and Hi-Token (80k samples, flat baseline only). At SFT scale text ties (GETok, Shikra), and new tokens are sensitive to initialization (HiVG). With about 30k targets and 1 GPU-hour, adding tokens is a risk. The upside is a **3.1–3.4× cheaper GRPO rollout** (§1.1), which only matters if rollouts dominate step time. Test it (H-A1-4); do not adopt it by default.
3. **Tokenizer choice is a free 1.5× on output length.** LFM's 3-digit chunks cost 5.7–6.1 tokens per vertex against 8.6–8.9. Among the candidates, LFM2.5-VL-450M is the smallest and the cheapest per rollout.
4. **Order and canonical start matter more than token format** (about 40 IoU points in LPM; SimpleSeg fails without them). RL then adds +10–12 points over SFT (SimpleSeg, Hi-GAR). SFT's job is format, legality and a good prior; polish comes from GRPO.
5. **Every compression that changes the meaning of a token raises error propagation** (delta, chain codes) or loses sub-pixel precision (pixel units, coarse grids). The 0–1000 grid's spare digit is cheap insurance on LIDC.
6. **Topology is the frontier.** Nobody has shown a VLM emitting holes. Our A0′ arm is the experiment, but our val split cannot measure it, so a probe is needed.
7. **Visual resolution, not output resolution, limits small-object precision.** One Qwen or LFM visual token covers 32 image px. LIDC is already upscaled to 512; RIGA 250 px frames get 64 tokens.

---

## 5. Recommendations (ranked)

### 5.1 Models to screen in Phase 2

Order of screening, smallest first; stop at the first that passes `phase2_entry`.

1. **LFM2.5-VL-450M (keep).**
   - Cheapest rollouts: 207 / 70 median answer tokens.
   - Real grounding prior (RefCOCO-M 81.28).
   - Phase 1 SFT and K12 already work.
   - Risks: not in TRL's GRPO CI; licence (research use OK; check before release); [0,1] float prior vs our integers (its 3B sibling switched to 0–1000, so expect fast adaptation).
2. **Qwen3.5-0.8B (add).**
   - Apache-2.0; native 0–1000 `bbox_2d`.
   - In TRL GRPO CI; vLLM with LoRA.
   - Gives a Qwen-family point between 450M and 2B.
   - Risks: GDN kernels and vLLM issues (§2); 248k vocabulary.
   - First run H-A1-1 (zero-shot prior) and H-A1-2 (Mac pilot).
3. **Qwen3-VL-2B (keep as the 2B rung and the "non-smallest" control).**
   - Proven in our K12 smoke; zero-shot box IoU 0.76.
   - **Qwen3.5-2B** is the alternative if Qwen3.5-0.8B proves the stack works (RefCOCO 84.3 vs 85.6, same tokenizer family).
4. **InternVL3.5-1B-HF (optional).** RefCOCO 85.4 at 1B with 0–1000 boxes and in TRL CI. Its cost is a fixed 256 visual tokens even for 250 px crops. Screen it only if 450M and 0.8B both fail.
5. **SmolVLM-256M (demote to optional control).**
   - No prior.
   - Excluded from TRL CI over bf16.
   - Open vLLM bug.
   - One visual token per 64 px at `longest_edge` 512.
   - CropVLM needed about 3 GPU-hours to learn 4-number boxes.
   - The "non-Qwen control" role (ROADMAP) is already covered by LFM, whose LM is not Qwen.
   - Keep it only if a 15-minute screen on the 4090 is essentially free.
6. **Watch list:** Gemma 4 E2B (1000 grid, Apache, TRL CI, but about 5B in memory and no grounding numbers); MiniCPM-V 4.7-1B (Qwen3.5-0.8B backbone).
7. **Baseline only:** Florence-2-large fine-tuned with SFT on its native `<loc>` polygons. It is the strongest polygon prior under 1B, but it is encoder-decoder, so it has no TRL GRPO path. It belongs to A4's baseline set.

### 5.2 Output representation options

1. **D10 as is** (compact JSON, digits 0–1000, canonical rings) for every model. Default.
2. **D10 with added coordinate tokens** `<x0..1000><y0..1000>`, initialized from the mean of each number's digit-token embeddings. Only if H-A1-4 passes; payoff is 3.1–3.4× fewer rollout tokens.
3. **LFM-only flat rings** `[x,y,x,y,...]`: −17 % tokens on LFM, about 0 on others. Not worth a second grammar unless rollout time dominates and H-A1-4 fails.
4. **A0′ = D10 grammar.** A0 = same grammar with `holes` emptied and parts kept (SimpleSeg-faithful, "A0-multi"). Report A0-single (largest outer ring) only as a ceiling.
5. **A0′ SFT oversamples non-simple LIDC targets** (×4–×8, H-A1-6), mirroring Phase 1's ×4 disagreement oversampling in GRPO rows.
6. **Distance-aware digit loss** (NTL-WAS with place-value weights, or Gaussian smoothing on the last digit): hand to A2 as an SFT loss option. The expected gain is boundary error, not IoU.

---

## 6. Rejected alternatives

| alternative | reason |
|---|---|
| Delta or relative vertex encoding | −25 % tokens (Qwen), but a ±1 slip in delta pixels costs 13 IoU points on LIDC and a random digit 36 (§1.2); errors accumulate (STAMP, VectorArk) |
| Pixel-unit coordinates (grid S in the prompt) | −23 % tokens on LIDC only. Loses the 0–1000 prior of Qwen, InternVL, Gemma 4 and LFM-3B. Last-digit slip = 1 px (−2.8 IoU). Cannot express bevel vertices (0.4 % of targets) without changing `targets/` |
| Coarser grid (G = 100 / 250) | Lossless on LIDC but −1 to −2.7 IoU and −1.5 to −2.2 % legality on RIGA (§1.5); there is no shared grid |
| Coarse grid + offset tokens (GETok) | Built for 7B with SAM decoding; our grid already matches the frame pixels |
| Joint xy tokens (OmniSVG) | Vocabulary of w·h (10⁶ at 1000²); needs pretraining-scale data |
| Hierarchical digit tokens (Hi-Token) | Only beats atomic tokens, never native digits; the gap shrinks to +3.9 mIoU against a tuned flat baseline; native single-digit tokenizers are already hierarchical digits |
| Mask tokens or decoders (LISA [SEG], Sa2VA, UniPixel, PaliGemma `<seg>`, STAMP, Text4Seg) | No polygon, so no exact legality certificate and no checker reward on the model's output; most need SAM and none has a TRL GRPO path |
| Ellipse or rounded primitives for disc and cup | Ceiling 0.974 / 0.966 IoU against ≥ 0.99 for the targets; cannot carry topology. Kept as a baseline |
| Keyhole or bridge encoding of holes | Self-intersections, fragile under crops and flips (P³, 2603.14764); our explicit holes are legal by grammar |
| Random start vertex in SFT targets | About −40 IoU points (LPM); random start only as a pretraining augmentation |
| Fixed vertex count | Collapses (Groundbench); our adaptive targets have p95 24 / 32 vertices |
| Spaced or float JSON | +22 % / 2× tokens (D10); space-separated numbers save nothing on digit tokenizers (§1.1) |
| Ovis2.5-2B | Best prior (RefCOCO 92.8), but remote code, no TRL path, no vLLM LoRA |
| Rex-Omni-3B, Isaac, Moondream 3, FastVLM, Eagle | Licence (research or non-commercial), own pipelines, or size |
| Phi-4-MM, Granite Vision, Kimi-VL | 1.3k–3.7k visual tokens or more than 5B total |
| Gemma 3n, Gemma 3 4B, MedGemma 4B | No grounding prior (MedGemma 1.5: chest X-ray boxes only); no vLLM LoRA for 3n; fixed 896 input |
| PaliGemma 2 | Gemma ToU; `<loc>`/`<seg>` vocabulary unrelated to polygons; 3B |

---

## 7. Hypotheses to test

All tests are sandbox tests on this Mac (CPU/MPS). Each uses ≤ 3 GB unless it says MPS training (≤ 7 GB, one at a time, scheduled by the coordinator). Data is read with `.venv/bin/python -I`. Lab code goes to `experiments/phase2/lab/A1/`. Priority: P1 decides the plan, P2 sharpens it, P3 is optional.

**H-A1-1 (P1): Qwen3.5-0.8B has a usable zero-shot coordinate prior on our crops, close to Qwen3-VL-2B.**
- Test: inference only, no training.
- Models: Qwen3.5-0.8B (about 1.8 GB bf16 download) and LFM2.5-VL-450M on MPS. Re-use Phase 1's Qwen3-VL-2B number (box IoU 0.76) or re-run it.
- Data: 100 val items per dataset from `data/prepared`, with their own prompt conventions (Qwen `bbox_2d` 0–1000; LFM `bbox` [0,1]).
- Ask for the box of the optic disc / the nodule. Score box IoU against the union-of-raters box, and the parse rate.
- InternVL3.5-1B-HF (about 2.1 GB) needs agreement to download.
- **Pass:** Qwen3.5-0.8B parse ≥ 90 % and median box IoU ≥ 0.6 on RIGA disc → keep it on the ladder. Below 0.4 on both datasets → it has no prior advantage over LFM; screen it only after LFM.

**H-A1-2 (P1): at equal steps, format learnability ranks LFM2.5-VL-450M ≥ Qwen3.5-0.8B > SmolVLM-256M.**
- Test: MPS LoRA SFT pilot, same harness as K12 (`objdel.train`).
- Training: r=64 on the LM (plus the projector), 300 steps, micro-batch 2 × grad-accum 4, on a fixed 2,000-target train subset (1,000 RIGA, 1,000 LIDC; one random rater per item).
- Every 100 steps: greedy decode of 64 val items. Record parse, legal, per-object IoU_tc, duplicate-vertex runs, and mean answer tokens. Also record seconds per step.
- **Pass:** a model reaching parse ≥ 95 % and legal ≥ 90 % within 300 steps is "learns the format". The ranking by median IoU_tc at step 300 orders the Vast screen.
- If SmolVLM is not within 0.05 IoU of the best, drop it from the Vast screen.

**H-A1-3 (P2): the 0–1000 digit grid is not worse than pixel-unit coordinates for LIDC, despite 23 % more tokens.**
- Test: same pilot as H-A1-2 on one model (the H-A1-2 winner) and the LIDC subset only. Arm (a) json1k; arm (b) json_px with "grid S" in the prompt, built from targets without bevels.
- Measure IoU_tc and boundary D_max (`objdel.eval`).
- **Pass (keep D10):** arm (a) D_max ≤ arm (b) D_max + 0.25 px and IoU within 0.01.
- **Fail:** arm (b) better by > 0.02 IoU → reopen D10 for LIDC.
- Priors (§1.2) predict pass.

**H-A1-4 (P1 if rollouts dominate GRPO step time, else P3): added coordinate tokens learn as fast as digits and cut completion tokens about 3×.**
- Test: pilot as in H-A1-2 on LFM2.5-VL-450M and Qwen3.5-0.8B.
- Arm (a): D10 digits.
- Arm (b): 2,002 tokens `<x0..1000><y0..1000>` (or 1,001 shared `<c…>`), each initialized to the mean of its digit-token embeddings plus a small sinusoidal numeric-order term (HiVG-style). Train the embedding and LM-head rows.
- Equal steps.
- **Pass:** legal and parse within 2 points, and median IoU_tc within 0.01 of digits at step 300. Then escalate to one Vast A/B run.
- **Fail:** keep digits. Report completion tokens per answer (expected: Qwen 91 / 33 vs 309 / 100, §1.1).

**H-A1-5 (P2): upscaling RIGA 250 px frames to 512 px improves cup and disc precision on Qwen and LFM.**
- 250 px frames get 64 visual tokens today; 512 px gets 256 (§1.3).
- Test: pilot on the H-A1-2 winner with the RIGA subset only. Arm (a) current frames; arm (b) model image upscaled to 512 with Pillow bicubic, as for LIDC. Masks, targets and W stay at the frame resolution, so this is only an image change: `FrameConfig(image_side=512)` for RIGA.
- Compare IoU_tc on val items with 250 px frames (17) and 500 px frames (45).
- **Pass:** median IoU_tc on 250 px items +0.01 or more, with no loss on 500 px items → adopt `image_side=512` for RIGA (a Phase 2 config decision).
- Also report the step-time cost.

**H-A1-6 (P1): A0′ learns holes and parts only when non-simple targets are oversampled.**
- Test: pilot on the H-A1-2 winner with the LIDC subset enriched with all 517 non-simple train targets.
- Arms: A0 (holes emptied, parts kept), A0′ ×1, A0′ ×4, A0′ ×8 (oversampling of non-simple targets). Same steps.
- Evaluate on the topology probe (H-A1-7) and on 200 simple val items.
- **Pass:** A0′ ×k recovers the probe's (β0, β1) on ≥ 50 % of probe items, with legal ≥ 90 %, and simple-item IoU_tc within 0.01 of A0. The smallest such k is the SFT setting.
- **Fail** (no k reaches 30 %) → holes are left to GRPO and the plan says so.

**H-A1-7 (P1, no training): the val split cannot separate A0 from A0′; a probe of about 300 items can.**
- Test, part 1 (computation): binomial power on the 48 val items with a non-simple rater; 0 are unanimous (§1.4).
- Test, part 2 (data): build a probe from **held-out LIDC train items with any non-simple rater** (366 items, kept out of SFT for this purpose). Optionally add synthetic carve-outs on real CT crops: annulus or two-blob masks, certified with `targets.build`, with geometry chosen inside the nodule size range.
- **Pass:** probe size gives a 95 % CI half-width ≤ 6 points on "rater-topology hit rate", and the probe is disjoint by patient from the SFT rows.
- Expected: the 48-item val slice gives about ±14 points at p = 0.5, which fails.

**H-A1-8 (P2): SFT checkpoints sample legal polygons at temperature 1 often enough for GRPO without constrained decoding.**
- Test: on the H-A1-2 pilot checkpoints, K = 8 samples at T = 1 on 64 val items. Measure parse, legal, duplicate rate and in-group σ of IoU_tc.
- Repeat with the LFM-specific failure counted separately: zero-led number chunks such as `042` give `leading_zero` (§1.1).
- **Pass:** legal ≥ 80 % at T = 1 and duplicates ≤ 10 % → no constrained decoding needed at GRPO start.
- **Fail** → adopt the `pattern()` regex in vLLM (D16 / A5).

**H-A1-9 (P2, no training): `pattern(names)` compiles with xgrammar for every candidate tokenizer and accepts every target's canonical token sequence.**
- Test: scratch venv with `xgrammar`. For LFM2.5-VL, Qwen3.5, Qwen3-VL and SmolVLM: compile `objdel.geometry.pattern` per schema, then walk the token IDs of all 30,843 targets through the matcher.
- Also check that LFM's non-canonical chunkings (`5|00`) stay allowed and zero-led chunks are masked.
- Report compile time and mask time per token.
- **Pass:** 100 % acceptance and < 50 µs per token for the mask. Coordinate with A5.

**H-A1-10 (P3, cross-reference A2): a place-value-weighted number-token loss lowers boundary error at equal IoU.**
- Test: pilot on Qwen3.5-0.8B (single-digit tokenizer): CE vs CE + NTL-WAS on coordinate digits, with weights 100 / 10 / 1 by place, λ ∈ {0.1, 0.3}.
- **Pass:** median D_max −10 % or more with IoU_tc within 0.01.

Hypotheses deliberately not proposed:
- Canonical start vs random start (LPM: about 40 points; settled).
- Fixed-N vs adaptive vertex count (Groundbench, VistaLLM; settled in D10/D11).
- Delta encoding (§1.2 is enough to reject it).

---

## 8. Sources

Primary sources read this session (verified):
- Hi-Token, arXiv 2608.03471 (tables 2, 4, 5; limitations): https://arxiv.org/html/2608.03471
- GETok, arXiv 2512.10554 (table 3 grid ablation; table 9 text vs tokens): https://arxiv.org/html/2512.10554v1
- LPM (Amazon, CVPRW 2026), table 3 ordering: https://cdn.amazon.science/aa/3c/d4ea78084f199eaf08743fc57e50/scipub-approval152134-45892548-rethinking-language-models-for-building-outline-extraction-from-remote-sensing-imagery.pdf
- STAMP, arXiv 2512.00395: https://arxiv.org/html/2512.00395
- CropVLM (SmolVLM-256M boxes), arXiv 2511.19820: https://arxiv.org/html/2511.19820
- GI endoscopy small VLMs (Florence-2, Qwen3.5-0.8B, InternVL3.5-1B), arXiv 2607.27122: https://arxiv.org/html/2607.27122
- MedMO, arXiv 2602.06965: https://arxiv.org/html/2602.06965
- STAMPlus abstract, arXiv 2608.02791: https://arxiv.org/abs/2608.02791
- Qwen3.5-0.8B card: https://huggingface.co/Qwen/Qwen3.5-0.8B
- LFM2.5-VL-450M card: https://huggingface.co/LiquidAI/LFM2.5-VL-450M
- Gemma 4 blog: https://huggingface.co/blog/gemma4
- Moondream segmentation update (benchmarks only): https://moondream.ai/blog/segmenting-update-2026-03-10
- TRL GRPO tests (`test_train_vlm` model list; SmolVLM excluded over bf16): https://github.com/huggingface/trl/blob/main/tests/test_grpo_trainer.py
- Configs and tokenizers (local download, tokenizer files only): Qwen/Qwen3.5-0.8B, Qwen/Qwen3.5-2B, OpenGVLab/InternVL3_5-1B-HF, LiquidAI/LFM2.5-VL-1.6B, unsloth/gemma-3n-E2B-it, unsloth/medgemma-4b-it, microsoft/Florence-2-base, vikhyatk/moondream2, microsoft/Phi-4-multimodal-instruct, openbmb/MiniCPM-V-4, apple/FastVLM-0.5B, HuggingFaceTB/SmolVLM-500M-Instruct

From sub-agents, primary pages, not re-read by me:
- Number Token Loss: https://arxiv.org/html/2411.02083v3
- R-VLM: https://arxiv.org/html/2507.05673v1
- Shikra: https://ar5iv.labs.arxiv.org/html/2306.15195
- HiVG: https://arxiv.org/html/2604.05072
- VectorArk: https://arxiv.org/pdf/2605.24398
- OmniSVG: https://arxiv.org/html/2504.06263v2
- LLM4SVG: https://arxiv.org/html/2412.11102
- Text4Seg++: https://arxiv.org/html/2509.06321
- RLRF: https://arxiv.org/html/2505.20793v2
- Floorplan vectorization: https://arxiv.org/html/2608.25608
- PolyTopoBench: https://arxiv.org/html/2609.32856
- P³: https://arxiv.org/html/2505.15379
- Bridge-vertex clipping: https://arxiv.org/abs/2603.14764
- Click2Poly (Florence-2 polygons): https://arxiv.org/html/2608.11424
- Qwen3-VL report: https://arxiv.org/html/2511.21631
- Qwen3-VL coordinates issue: https://github.com/QwenLM/Qwen3-VL/issues/1616
- Qwen3.5-2B card: https://huggingface.co/Qwen/Qwen3.5-2B
- Qwen3.5 in transformers: https://huggingface.co/docs/transformers/main/en/model_doc/qwen3_5
- LFM2.5-VL-3B card: https://huggingface.co/LiquidAI/LFM2.5-VL-3B
- LFM licence: https://huggingface.co/LiquidAI/LFM2-VL-450M/blob/main/LICENSE
- InternVL3.5: https://arxiv.org/html/2508.18265
- InternVL3: https://arxiv.org/html/2504.10479
- InternVL data format: https://internvl.readthedocs.io/en/latest/get_started/chat_data_format.html
- Ovis2.5: https://huggingface.co/AIDC-AI/Ovis2.5-2B ; https://arxiv.org/html/2508.11737
- Gemma 4 docs: https://huggingface.co/docs/transformers/main/en/model_doc/gemma4
- Gemma 4 report: https://arxiv.org/html/2607.02770
- MedGemma 1.5: https://huggingface.co/google/medgemma-1.5-4b-it
- PaliGemma 2: https://arxiv.org/html/2412.03555
- Florence-2: https://arxiv.org/abs/2311.06242 ; https://huggingface.co/microsoft/Florence-2-large/blob/main/processing_florence2.py
- Rex-Omni: https://arxiv.org/html/2510.12798
- Isaac SDK: https://github.com/perceptron-ai-inc/perceptron
- Moondream segmentation paper: https://arxiv.org/html/2604.02593
- Sa2VA-Qwen3-VL-4B: https://huggingface.co/ByteDance/Sa2VA-Qwen3-VL-4B
- X-SAM: https://arxiv.org/html/2508.04655
- UniPixel: https://github.com/PolyU-ChenLab/UniPixel
- DenseMLLM: https://arxiv.org/html/2602.14134
- MediX-R1-2B: https://huggingface.co/MBZUAI/MediX-R1-2B
- vLLM supported models: https://docs.vllm.ai/en/latest/models/supported_models.html
- TRL vLLM integration: https://huggingface.co/docs/trl/main/en/vllm_integration
- TRL issues and PRs: #6718 https://github.com/huggingface/trl/pull/6718 ; #7198 https://github.com/huggingface/trl/issues/7198 ; #7761 https://github.com/huggingface/trl/issues/7761 ; #5269 https://github.com/huggingface/trl/issues/5269 ; #4061 https://github.com/huggingface/trl/issues/4061
- vLLM GDN kernel bug: https://github.com/vllm-project/vllm/issues/34948
- SmolVLM image processor (longest_edge 2048 default): https://raw.githubusercontent.com/huggingface/transformers/main/src/transformers/models/smolvlm/image_processing_smolvlm.py
- LFM2-VL image processor: https://raw.githubusercontent.com/huggingface/transformers/main/src/transformers/models/lfm2_vl/image_processing_lfm2_vl.py

Phase 1 evidence re-used: `docs/research/notes/Phase 1 stack and algorithms/output_format.md`, `training_stack.md`, `experiments/phase1/PLAN.md` §2 (D10, D11, D12), `experiments/phase1/RESULTS.md` §8, `data/targets/summary.json` (target_set_sha256 `1af819f4…`).
