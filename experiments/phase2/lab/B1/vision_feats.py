"""POST-HOC diagnostic (CPU): are the visual features that reach the language model still image-dependent after LoRA
SFT? For 4 eval items (centred view) and a black image, the merged checkpoint's image features (model.get_image_features,
the projector output the LM sees) are compared: mean pairwise cosine distance between images, feature norm, NaN count.
Base model for reference. Reads untrusted data: run with python -I.

    .venv/bin/python -I -m experiments.phase2.lab.B1.vision_feats --arm qwen3.5-0.8b --run qwen3.5-0.8b_s0_n150
"""
import argparse
import json

import torch
from peft import PeftModel
from PIL import Image
from transformers import AutoModelForImageTextToText

from objdel.train import dataset as D

from . import common as C


def feats(model, proc, image):
    msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": D.PROMPTS["lidc"]}]}]
    prompt = proc.apply_chat_template(msgs, tokenize=False, add_generation_prompt=True)
    b = proc(text=[prompt], images=[[image]], return_tensors="pt")
    kw = {k: v for k, v in b.items() if k not in ("input_ids", "attention_mask", "mm_token_type_ids")}
    with torch.no_grad():
        f = model.model.get_image_features(**kw)
    f = getattr(f, "pooler_output", f)
    f = torch.cat([x.reshape(-1, x.shape[-1]) for x in f]) if isinstance(f, (list, tuple)) else f.reshape(-1, f.shape[-1])
    return f.float()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=list(C.ARMS))
    ap.add_argument("--run", required=True, help="run directory name, or base")
    a = ap.parse_args()
    torch.set_num_threads(2)
    proc = C.processor(a.arm, "lidc")
    model = AutoModelForImageTextToText.from_pretrained(C.ARMS[a.arm]["hf"], dtype=torch.float32)
    if a.run != "base":
        model = PeftModel.from_pretrained(model, str(C.SCRATCH / "runs" / a.run / "adapter")).merge_and_unload()
    model.eval()
    V = C.views("lidc", C.eval_records("lidc", 64)[:4])
    ims = [x["image"] for x in V["centred"]] + [Image.new("RGB", V["centred"][0]["image"].size)]
    F = [feats(model, proc, im) for im in ims]
    d = [float(1 - torch.nn.functional.cosine_similarity(F[i].flatten(), F[j].flatten(), 0))
         for i in range(len(F)) for j in range(i + 1, len(F))]
    print(json.dumps({"arm": a.arm, "run": a.run, "shape": list(F[0].shape),
                      "mean_cosine_distance_between_images": round(sum(d) / len(d), 5),
                      "feature_norm_mean": round(float(torch.stack([f.norm(dim=-1).mean() for f in F]).mean()), 3),
                      "nan": int(sum(torch.isnan(f).sum() for f in F))}))


if __name__ == "__main__":
    main()
