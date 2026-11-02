# Phase 2 results: SFT on Vast AI

Plan: [PLAN.md](PLAN.md). Build log: [V0_V1.md](V0_V1.md). Issue ledger: [impl/ISSUES.md](impl/ISSUES.md).

**Status: in progress.** R1 rung 1 (LFM2.5-VL-450M) is done and frozen below: trained, selected, judged; it fails the
§5.1 entry gates. Later rungs are added as new sections; nothing in a frozen section changes. Hashes are SHA-256 (first
16 hex digits); `runs/` and `data/` are git-ignored, so the hashes pin the local files.

## 1. R1, rung 1: LFM2.5-VL-450M, full fine-tune (2026-11-02, frozen)

**Run.** Box 58836372 (RTX 4090, Quebec, $0.401/h; billed $0.727 for R0 rerun + R1). Config
`configs/r1_lfm.yaml`: bundle view, A0′, seed 1, WSD 2 epochs + branches 0.5 / 1 / 2, effective batch 32, HF sampler
(D13 fallback; R0 parity failed). Training code at `fa86f5b`, select and probes at `024089f` (CPU G-vision fix).

| Phase | Steps | s / step | Peak CUDA memory | G-vision (step: ratio) |
|---|---|---|---|---|
| stable | 1,542 | 0.69 | 20.8 GiB | 386: 1.00, 771: 0.96, 1542: 0.84 (all ok) |
| e0.5 / e1 / e2 decays | to 444 / 887 / 1,773 | — | 17.8 / 20.8 / 20.9 GiB | — |

**Selection (D12 R-B, 256 shifted-val items, K = 8).** e2 ineligible: legal 87.0 % < 90 % (its pass rate 0.69 was
the best: e0.5 0.40, e1 0.48). Eligible e0.5 and e1; the signal veto did not fire (e0.5's sd_geom 9 % above e1's, under the 15 % needed).
**Selected: e1.**

**Entry gates on e1 (§5.1; full views, 1,084 items, K = 8): FAIL.**

| Gate | Result | Value (requirement) |
|---|---|---|
| G-format | pass | parse 98.9 %, legal 91.4 % (≥ 98 %, ≥ 90 %) |
| G-geometry | **fail** | median IoU_tc shifted 0.495, rescaled 0.092 (≥ 0.60); rescaled legal 71 % |
| G-prior-b | pass | paired gain over the image-blind prior: shifted +0.37 (CI 0.34–0.40), rescaled +0.18 (0.16–0.21) |
| G-precision | **fail** | probe IoU 0.61, legal 97 %, centroid 1.55 px, HD95 2.8 px; checks assd / iou / centroid / win fail |
| G-topology-a | **fail** | simple controls correct 93.7 % (≥ 95 %); hole recall 0.5 % (report line 30 %), any-of-8 4 %; two-part recall 0.2 % |
| G-signal | pass | duplicates 0.07 %, σ-share 99.6 % |
| G-vision | pass | projector ratio 0.92, image-use TV 0.10 |

G-probe (R2) not run: §5.1 runs it only after the gates pass (`g_probe.json`: skipped).

**Diagnosis (Mac, no GPU).**
- *Crop size, not the rescaled view itself.* By crop side W: shifted W=50 median IoU 0.52 (legal 94 %), shifted
  W=100 0.27 (64 %); rescaled W=100 0.13 (71 %). Rescaling moves 978 nodules from W=50 to W=100 crops.
- *Small in frame.* Nodules spanning < 10 % of the crop (99 % of them W=100 crops): median IoU 0.00; 15–30 %: 0.44–0.50.
- *Detection, not drawing.* MedSAM given the expert box (D17 baseline, every 2nd item, 542 per view): shifted 0.77,
  rescaled 0.80, jittered box 0.75 / 0.75; image-blind prior 0.00. MedSAM does not drop with the crop size because the
  box removes the search. On the same items e1 reaches MedSAM on 3 % (shifted) and 0 % (rescaled W=100).
- *Gallery* (CT crop, every expert mask, the model's samples, drawn from objdel.geometry's own masks): overlays agree
  with the scores (no frame or axis error). W=50: the right nodule, outline 1–2 px off. W=100: often the wrong bright
  spot (vessel) or a sprawling outline; on the right spot IoU 0.42–0.45. Holes: never drawn, even 111–200 px holes;
  the outline fills them.

**Artifacts.**

| What | Where | SHA-256 |
|---|---|---|
| e1 weights | `runs/v1/r1-lfm/branches/e1/model.safetensors`; HF `svenk029/objdel-sft` (private) branch `r1-lfm`, commit `f5f7cf039a60` | `6962aa203734edf7` |
| run manifest | `runs/v1/r1-lfm/manifest.json` (code `c16c409072ce1957`, view `f240358ccc8ea2d8`) | `388e652fc5e280f7` |
| selection | `runs/v1/r1-lfm/select/selected.json` | `c1db0ad234d80019` |
| report / gates | `runs/v1/r1-lfm/select/report.json` / `gates.json` | `868776fb4d31490b` / `adff14d53b319a50` |
| G-probe | `runs/v1/r1-lfm/g_probe/g_probe.json` | `fd5ac056beb008a0` |
| HF sampling rate (R0) | `runs/v1/hf-rate-lfm/hf_rate.json`: 0.132 s / sample on base weights | `150846b312388b96` |
| G-vision reference (R0) | `runs/v1/vision_ref/lfm2.5-vl-450m.json` | `7551e33d173ea1fe` |
| baselines | `data/eval/baselines_lidc_v1.json` (prior, MedSAM, MedSAM jitter; shifted + rescaled, every 2nd item) | `486ef19699d775bc` |
| gallery | `runs/v1/r1-lfm/gallery/*.png` | — |

**Open (decided outside this section).** The ladder's next rung (Qwen3.5-0.8B, LoRA) or a data change for small
nodules in large crops; holes are C2 / Phase 3 work by design (§5.2).
