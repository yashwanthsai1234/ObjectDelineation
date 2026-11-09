import numpy as np
from PIL import Image
from transformers import AutoProcessor

R = {
    "smolvlm-256m": "HuggingFaceTB/SmolVLM-256M-Instruct",
    "smolvlm-500m": "HuggingFaceTB/SmolVLM-500M-Instruct",
    "lfm2.5-vl-450m": "LiquidAI/LFM2.5-VL-450M",
    "lfm2.5-vl-1.6b": "LiquidAI/LFM2.5-VL-1.6B",
    "qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct",
    "qwen3.5-0.8b": "Qwen/Qwen3.5-0.8B",
    "internvl3.5-1b": "OpenGVLab/InternVL3_5-1B-HF",
    "gemma3n-e2b": "unsloth/gemma-3n-E2B-it",
    "medgemma-4b": "unsloth/medgemma-4b-it",
}
for k, r in R.items():
    try:
        p = AutoProcessor.from_pretrained(r)
    except Exception as e:
        print(k, "FAIL", type(e).__name__, str(e)[:120])
        continue
    row = []
    for side in (250, 500, 512):
        img = Image.fromarray(np.zeros((side, side, 3), np.uint8))
        msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": "x"}]}]
        try:
            prompt = p.apply_chat_template(msgs, add_generation_prompt=True)
            n_all = len(p(text=prompt, images=[img], return_tensors="pt")["input_ids"][0])
            n_txt = len(
                p.tokenizer(
                    p.apply_chat_template(
                        [{"role": "user", "content": [{"type": "text", "text": "x"}]}], add_generation_prompt=True
                    ),
                    add_special_tokens=False,
                )["input_ids"]
            )
            row.append(f"{side}px:{n_all - n_txt}")
        except Exception as e:
            row.append(f"{side}px:ERR {type(e).__name__} {str(e)[:80]}")
    print(k, " ".join(row))
