"""Where one MPS training step goes: forward + backward time of one SFT row per arm, split into the vision tower and
the decoder layers by module (linear-attention / GDN vs full attention vs conv), with module-level synchronised timers.

    HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B1.prof_step --arm qwen3.5-0.8b

Reads untrusted data: run with python -I.
"""
import argparse
import json
import time
from collections import defaultdict

import torch
from peft import get_peft_model
from transformers import AutoModelForImageTextToText

from objdel.train import dataset as D

from . import common as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=list(C.ARMS))
    ap.add_argument("--dtype", default="float32")
    ap.add_argument("--reps", type=int, default=3)
    a = ap.parse_args()
    C.cap_mps()
    proc = C.processor(a.arm, "lidc")
    rec = D.records(C.PREPARED, "lidc", "train")[0]
    text = next(t for t in rec["targets"] if t)
    from PIL import Image
    image = Image.open(D._image(C.PREPARED, rec["item_id"])).convert("RGB")
    msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": D.PROMPTS["lidc"]}]},
            {"role": "assistant", "content": [{"type": "text", "text": text}]}]
    s = proc.apply_chat_template(msgs, tokenize=False)
    batch = proc(text=[s], images=[[image]], return_tensors="pt")
    batch = {k: v.to("mps") for k, v in batch.items() if hasattr(v, "to")}
    model = AutoModelForImageTextToText.from_pretrained(C.ARMS[a.arm]["hf"], dtype=getattr(torch, a.dtype)).to("mps")
    model = get_peft_model(model, C.LORA)
    model.train()
    times = defaultdict(float)

    def wrap(name, mod):
        def pre(*_):
            torch.mps.synchronize()
            mod._t = time.time()

        def post(*_):
            torch.mps.synchronize()
            times[name] += time.time() - mod._t
        mod.register_forward_pre_hook(pre)
        mod.register_forward_hook(post)

    base = model.base_model.model
    vis = getattr(base.model, "visual", None) or getattr(base.model, "vision_tower", None) \
        or getattr(base.model, "vision_model", None)
    if vis is not None:
        wrap("vision", vis)
    lm = getattr(base.model, "language_model", None) or getattr(base.model, "text_model", None)
    for layer in lm.layers:
        kind = type(getattr(layer, "linear_attn", None) or getattr(layer, "self_attn", None)
                    or getattr(layer, "conv", None) or layer).__name__
        wrap(f"layer:{kind}", layer)
    out = {"arm": a.arm, "dtype": a.dtype, "tokens": int(batch["input_ids"].shape[1])}
    for r in range(a.reps + 1):
        times.clear()
        torch.mps.synchronize()
        t = time.time()
        loss = model(**batch, logits_to_keep=128).logits.logsumexp(-1).mean()  # answer-sized logits only, as TRL
        torch.mps.synchronize()
        fwd = time.time() - t
        loss.backward()
        torch.mps.synchronize()
        if r:  # first rep is warm-up
            out[f"rep{r}"] = {"forward_s": round(fwd, 3), "fwd_bwd_s": round(time.time() - t, 3),
                              **{k: round(v, 3) for k, v in times.items()}}
        model.zero_grad()
    out["peak_mps_gb"] = round(C.mps_gb(), 2)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
