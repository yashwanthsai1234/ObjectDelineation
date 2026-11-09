import time, sys
from PIL import Image
from transformers import AutoProcessor
from trl.data_utils import apply_chat_template
proc = AutoProcessor.from_pretrained("Qwen/Qwen3-VL-2B-Instruct", local_files_only=True)
img = Image.new("RGB", (512, 512))
text = "Outline the lung nodule as a JSON polygon on a 0-1000 grid. " * 6
prompt = [{"role": "user", "content": [{"type": "image", "image": img}, {"type": "text", "text": text}]}]
for _ in range(3):
    apply_chat_template({"prompt": prompt}, proc)
t = time.perf_counter()
texts = [apply_chat_template({"prompt": prompt}, proc)["prompt"] for _ in range(256)]
t1 = time.perf_counter()
out = proc(images=[img] * 256, text=texts, padding=True, return_tensors="pt")
t2 = time.perf_counter()
proc.tokenizer(texts and [proc.image_token * 256 + text] * 256, padding=True)
t3 = time.perf_counter()
print(f"256 renders {t1-t:.3f} s; processor call (pixels+text) {t2-t1:.3f} s; tokenizer-only 256 rows {t3-t2:.3f} s; ids shape {tuple(out['input_ids'].shape)}")
