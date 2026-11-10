"""R7 feasibility: image tokens, effective input resolution and crop-pixels per token / per encoder patch.

    HF_HOME=./hf HF_HUB_OFFLINE=1 .venv/bin/python -I r7_resolution.py
Processors only (no weights). For each model x crop side S (RIGA disc crops <= 500 px, LIDC 64-128 px), optionally
pre-upscaled to U px, report: LLM image tokens, the resolution the encoder actually sees (content, not padding),
scale = encoder px per crop px, crop px per encoder patch (16 px at encoder scale) and per LLM token.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

from common import HERE, model_path  # noqa: E402

CONFIGS = {
    "smolvlm-256m@512": ("smolvlm-256m", {"size": {"longest_edge": 512}}),
    "smolvlm-256m@1024": ("smolvlm-256m", {"size": {"longest_edge": 1024}}),
    "lfm2.5-vl-450m": ("lfm2.5-vl-450m", {}),
    "lfm2.5-vl-450m@max1024": ("lfm2.5-vl-450m", {"max_image_tokens": 1024}),
    "qwen3-vl-2b": ("qwen3-vl-2b", {}),
}
CASES = [(500, None), (250, None), (448, None), (128, None), (100, None), (64, None),
         (128, 512), (64, 512), (100, 500), (64, 256), (250, 500)]


def probe(proc, name, img):
    tok = proc.tokenizer
    msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": "x"}]}]
    text = proc.apply_chat_template(msgs, add_generation_prompt=True)
    enc = proc(text=[text], images=[[img]], return_tensors="pt")
    tid = getattr(proc, "image_token_id", None) or tok.convert_tokens_to_ids(proc.image_token)
    n_tok = int((enc["input_ids"][0] == tid).sum())
    S = img.size[0]
    if name.startswith("smolvlm"):
        pam = enc["pixel_attention_mask"][0]  # (tiles, H, W) valid pixels
        # the last tile is the global image when split; with one tile it is the image itself
        tiles = pam.shape[0]
        g = pam[-1].numpy()
        content = int(g.any(1).sum()), int(g.any(0).sum())
        if tiles > 1:  # local tiles cover the image at (tiles-1) x 512 resolution
            rows = int(np.sqrt(tiles - 1))
            content = (rows * 512, rows * 512)
        eff = content[0]
        patch_enc = 16
    elif name.startswith("lfm"):
        ss = enc["spatial_shapes"].numpy()  # (tiles, 2) patch grid per tile (incl. thumbnail)
        if len(ss) > 1:  # tiles at native 512 + a thumbnail
            eff = int(round(np.sqrt(sum(h * w for h, w in ss[:-1]))) * 16)
        else:
            eff = int(ss[0][0] * 16)
        content = tuple(int(x) * 16 for x in ss[-1])
        patch_enc = 16
    else:
        thw = enc["image_grid_thw"][0].numpy()
        eff = int(thw[1] * 16)
        content = (int(thw[1] * 16), int(thw[2] * 16))
        patch_enc = 16
    scale = eff / S
    return {"tokens": n_tok, "encoder_side_px": eff, "scale_enc_px_per_crop_px": round(scale, 3),
            "crop_px_per_encoder_patch": round(patch_enc / scale, 2),
            "crop_px_per_llm_token_side": round(S / np.sqrt(n_tok), 2), "content_or_thumb": content}


def main():
    from transformers import AutoProcessor
    rng = np.random.default_rng(0)
    out = {}
    for cname, (m, kw) in CONFIGS.items():
        try:
            proc = AutoProcessor.from_pretrained(model_path(m), **kw)
        except Exception as e:  # noqa: BLE001
            out[cname] = {"error": repr(e)[:200]}
            continue
        out[cname] = {}
        for S, U in CASES:
            img = Image.fromarray(rng.integers(0, 255, (S, S, 3), dtype=np.uint8))
            if U:
                img = img.resize((U, U), Image.BICUBIC)
            r = probe(proc, cname, img)
            if U:  # express everything per ORIGINAL crop pixel
                f = U / S
                r["scale_enc_px_per_crop_px"] = round(r["scale_enc_px_per_crop_px"] * f, 3)
                r["crop_px_per_encoder_patch"] = round(r["crop_px_per_encoder_patch"] / f, 2)
                r["crop_px_per_llm_token_side"] = round(r["crop_px_per_llm_token_side"] / f, 2)
            key = f"{S}" + (f"->up{U}" if U else "")
            out[cname][key] = r
            print(cname, key, r, flush=True)
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / "r7_resolution.json").write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
