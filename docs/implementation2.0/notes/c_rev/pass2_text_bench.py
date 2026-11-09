import time
from PIL import Image
from transformers import AutoProcessor
from trl.data_utils import apply_chat_template
proc = AutoProcessor.from_pretrained("Qwen/Qwen3-VL-2B-Instruct", local_files_only=True)
img = Image.new("RGB", (512, 512))
text = "Outline the lung nodule as a JSON polygon on a 0-1000 grid. " * 6
prompt = [{"role": "user", "content": [{"type": "image", "image": img}, {"type": "text", "text": text}]}]
texts = [apply_chat_template({"prompt": prompt}, proc)["prompt"] for _ in range(256)]
inner = proc.image_processor
one = inner([img], return_tensors="pt")
import torch
from transformers.feature_extraction_utils import BatchFeature
class Cached:  # stands in for ImageMemo with every image a hit
    def __init__(s, inner): s.inner = inner
    def __getattr__(s, k): return getattr(s.inner, k)
    def __call__(s, images, **kw):
        n = len(images)
        return BatchFeature({"pixel_values": one["pixel_values"][:1], "image_grid_thw": one["image_grid_thw"].repeat(n, 1)}, tensor_type="pt")
proc.image_processor = Cached(inner)
for _ in range(2): proc(images=[img] * 256, text=texts, padding=True, return_tensors="pt")
ts = []
for _ in range(5):
    t = time.perf_counter(); proc(images=[img] * 256, text=texts, padding=True, return_tensors="pt"); ts.append(time.perf_counter() - t)
print("pass-2 processor call with all-hit pixels, 5 runs (s):", [round(x, 3) for x in ts])
