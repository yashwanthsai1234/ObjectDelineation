"""Bottom-up GPU-hour and $ estimate for the Phase 2 SFT screening (lab evidence for A3_vast_gpus.md).

    python3 -I experiments/phase2/lab/A3/estimate_cost.py

Model: training FLOPs = 6 * (N_lm * lm_tokens + N_vit * vit_patches) per row (forward + backward, dense; attention
and the loss head ignored), rows = every voting rater target (M7: "every rater is a row"), divided by
peak dense BF16 (FP32-accumulate) tensor throughput * MFU. MFU is a range, not a measurement: small models at short
sequences run far below peak (launch overhead, data loading, no packing), so 0.10-0.30 is assumed; the first Vast
session replaces it with measured s/step (RESULTS.md section 8, "Cost" gate).

Inputs and their sources (see the note for URLs):
  rows: data/prepared/manifest.json (LIDC train 23,469 targets, RIGA train 3,282)
  image tokens: M7_train.md (SmolVLM 64 tokens at 512 px, 320 at 1024 = 5 tiles); LFM2.5-VL and Qwen3-VL ~256 at
      512 px (16 px patches, 2x2 merge) -- UNVERIFIED here, A1 owns the exact counts
  answer tokens: M5/M7 (nodule ~100, RIGA p50 310 Qwen tokens; LFM 5.2 vs 8.2 tokens per vertex, D10)
  parameter counts: model cards (SmolVLM-256M = SmolLM2-135M + SigLIP-B/16 93M; LFM2.5-VL-450M = ~350M LM + ~86M
      SigLIP2-B; Qwen3-VL-2B = ~1.7B LM + ~0.4B ViT, the ViT size UNVERIFIED)
  peak dense BF16 TFLOPS: NVIDIA whitepapers / product pages (GeForce cards at FP32 accumulate);
      RTX PRO 6000 is derived (4000 sparse FP4 TOPS / 8), UNVERIFIED
  prices: Vast public offer search, 2026-10-19 ~24:00 UTC, verified, >= 8 effective cores, reliability > 0.98,
      CUDA >= 12.8, 1 GPU, 100 GB storage included (vast_snapshot.json)
"""

ROWS = {"lidc": 23469, "riga": 3282}
PROMPT = 70

# per model: LM params, ViT params, (image tokens, ViT patches, answer tokens) per dataset, full-FT memory note
MODELS = {
    "smolvlm-256m": dict(n_lm=135e6, n_vit=93e6, lidc=(64, 1024, 100), riga=(320, 5120, 330)),
    "lfm2.5-vl-450m": dict(n_lm=350e6, n_vit=86e6, lidc=(256, 1024, 63), riga=(256, 1024, 200)),
    "qwen3-vl-2b": dict(n_lm=1.7e9, n_vit=0.4e9, lidc=(256, 1024, 100), riga=(240, 960, 330)),
}

# peak dense BF16 tensor TFLOPS (FP32 accumulate), GPU memory GB, on-demand $/h (min, median) from the snapshot
GPUS = {
    "RTX 3090": dict(tflops=71, mem=24, price=(0.24, 0.33)),
    "RTX 4090": dict(tflops=165.2, mem=24, price=(0.48, 0.66)),
    "RTX 5090": dict(tflops=209.5, mem=32, price=(0.56, 0.96)),
    "RTX PRO 6000": dict(tflops=500, mem=96, price=(1.36, 1.57)),
    "A100 SXM4 80GB": dict(tflops=312, mem=80, price=(1.09, 1.19)),
    "H100 SXM": dict(tflops=989, mem=80, price=(2.20, 5.49)),
}

# full fine-tuning footprint ~16 B/param (bf16 weights + grads, fp32 Adam moments + master) + activations
FULL_FT_GB = {"smolvlm-256m": 6, "lfm2.5-vl-450m": 10, "qwen3-vl-2b": 38}


def flops_per_epoch(m):
    p = MODELS[m]
    total = 0.0
    for ds, rows in ROWS.items():
        img_tok, patches, ans = p[ds]
        lm_tokens = img_tok + PROMPT + ans
        total += rows * 6 * (p["n_lm"] * lm_tokens + p["n_vit"] * patches)
    return total


def main(epochs=2, mfu=(0.10, 0.30), eval_h=0.25):
    print(f"epochs={epochs}, MFU {mfu[0]:.2f}-{mfu[1]:.2f}, + {eval_h} h eval sampling per run\n")
    for m in MODELS:
        f = flops_per_epoch(m)
        print(f"{m}: {f:.2e} FLOPs/epoch, full FT ~{FULL_FT_GB[m]} GB")
        for g, s in GPUS.items():
            fits = FULL_FT_GB[m] <= s["mem"] * 0.85
            hi = epochs * f / (s["tflops"] * 1e12 * mfu[0]) / 3600 + eval_h
            lo = epochs * f / (s["tflops"] * 1e12 * mfu[1]) / 3600 + eval_h
            usd = (lo * s["price"][0], hi * s["price"][1])
            note = "" if fits else "  (full FT does not fit: LoRA / 8-bit optimizer only)"
            print(f"  {g:15s} {lo:5.2f}-{hi:5.2f} h  ${usd[0]:5.2f}-{usd[1]:6.2f} per run{note}")
        print()
    # Phase 2 package: 3 screening runs (one per model, A0') + A0 vs A0' on the chosen model (2 more runs)
    print("Phase 2 package = screening (3 models x 1 run) + A0 and A0' on the chosen model (2 runs):")
    for chosen in ("smolvlm-256m", "qwen3-vl-2b"):
        for g, s in GPUS.items():
            h_lo = h_hi = 0.0
            for m in list(MODELS) + [chosen, chosen]:
                f = flops_per_epoch(m)
                h_hi += epochs * f / (s["tflops"] * 1e12 * mfu[0]) / 3600 + eval_h
                h_lo += epochs * f / (s["tflops"] * 1e12 * mfu[1]) / 3600 + eval_h
            print(f"  chosen={chosen:14s} {g:15s} {h_lo:5.1f}-{h_hi:5.1f} GPU-h  ${h_lo * s['price'][0]:6.2f}-{h_hi * s['price'][1]:6.2f}")


if __name__ == "__main__":
    main()
