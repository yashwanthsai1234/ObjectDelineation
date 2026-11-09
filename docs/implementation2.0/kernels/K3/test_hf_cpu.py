"""K3 (a) and (c) against Hugging Face, on the CPU (no GPU needed; runs on macOS).

(a) The HF Qwen3-VL fast image processor's pixel_values (fp32, and their bf16 rounding) vs the K3 CPU reference and
    the host emulation of the GPU kernel (k3_cpu_cli preprocess): byte equality on 32 random + 9 adversarial images
    (3 of them 'digit' images that prove the layout at every position).
(c) HF Qwen3VLModel.get_rope_index on a real Qwen3-VL prompt (chat template + one 512 x 512 image) vs the K3 template +
    extrapolation (host emulation of the GPU kernel) and the C++ port of get_rope_index: equality on 16 rows with
    random left pads, answer lengths 0..1536, a hole in a mask, plus adversarial rows that must raise the error flag.

Run (the model's processor files must be in the HF cache; nothing is downloaded):
    HF_HUB_OFFLINE=1 REPO/.venv/bin/python -I test_hf_cpu.py --cli build/k3_cpu_cli --tmp /some/scratch/dir
"""

import argparse
import os
import subprocess
import sys
import types as pytypes

import numpy as np
import torch
from PIL import Image
from transformers import AutoProcessor
from transformers.models.qwen3_vl.modeling_qwen3_vl import Qwen3VLModel

MODEL = "Qwen/Qwen3-VL-2B-Instruct"
PROMPT = ('Outline the lung nodule in this CT crop. Answer with JSON only: '
          '{"nodule":[{"outer":[[x,y],...],"holes":[[[x,y],...]]}]}, '
          'integers 0-1000 relative to the crop, outer rings clockwise, holes counter-clockwise.')  # dataset.PROMPTS


def images(rng):
    out = [rng.integers(0, 256, (512, 512, 3), dtype=np.uint8) for _ in range(32)]
    out.append(np.zeros((512, 512, 3), np.uint8))
    out.append(np.full((512, 512, 3), 255, np.uint8))
    i = np.arange(512 * 512 * 3, dtype=np.int64).reshape(512, 512, 3)
    out.append(((i * 7 + i // 3) & 255).astype(np.uint8))
    # digit images: byte k of each source element's flat index. The value map is injective, so equal outputs on
    # these 3 images prove the index map (patch layout + temporal copy) equal at every output position
    out += [((i >> (8 * k)) & 255).astype(np.uint8) for k in range(3)]  # every value in every channel, scattered
    yy, xx = np.mgrid[:512, :512]
    out.append(np.repeat((((yy // 16 + xx // 16) & 1) * 255).astype(np.uint8)[..., None], 3, 2))  # patch checkerboard
    out.append(np.stack([xx & 255, yy & 255, (xx + yy) & 255], -1).astype(np.uint8))
    g = rng.integers(0, 256, (512, 512), dtype=np.uint8)
    out.append(np.repeat(g[..., None], 3, 2))  # grayscale CT stored as RGB, as the dataset's PNGs
    return out


def test_preprocess(proc, cli, tmp, rng):
    ip = proc.image_processor
    imgs = images(rng)
    hf = ip([Image.fromarray(a) for a in imgs], return_tensors="pt")
    pv, grid = hf["pixel_values"], hf["image_grid_thw"]
    assert pv.dtype == torch.float32 and tuple(pv.shape) == (len(imgs) * 1024, 1536), pv.shape
    assert (grid == torch.tensor([1, 32, 32])).all()
    mean, std, _ = ip._fuse_mean_std_and_rescale_factor(True, tuple(ip.image_mean), tuple(ip.image_std), True,
                                                        ip.rescale_factor, torch.device("cpu"))
    mean, std = [float(v) for v in mean], [float(v) for v in std]
    src = os.path.join(tmp, "imgs.u8")
    np.stack(imgs).tofile(src)
    f32, b16 = os.path.join(tmp, "out.f32"), os.path.join(tmp, "out.bf16")
    subprocess.run([cli, "preprocess", src, str(len(imgs)), "512", "512", "1", *[f"{v:.9g}" for v in mean + std],
                    f32, b16], check=True)
    ours = torch.from_numpy(np.fromfile(f32, dtype=np.float32)).view(pv.shape)
    ours16 = torch.from_numpy(np.fromfile(b16, dtype=np.int16)).view(pv.shape)
    eq = torch.equal(ours.view(torch.int32), pv.view(torch.int32))
    eq16 = torch.equal(ours16, pv.to(torch.bfloat16).view(torch.int16))
    diff = int((ours.view(torch.int32) != pv.view(torch.int32)).sum())
    print(f"(a) preprocess: {len(imgs)} images, mean' {mean} std' {std}: fp32 bit-equal {eq} ({diff} differ), "
          f"bf16 bit-equal {eq16}")
    return eq and eq16


def hf_rope(ids, mask, types, grid):
    fake = pytypes.SimpleNamespace(
        config=pytypes.SimpleNamespace(vision_config=pytypes.SimpleNamespace(spatial_merge_size=2)),
        get_vision_position_ids=lambda *a, **k: Qwen3VLModel.get_vision_position_ids(None, *a, **k))
    pos, _ = Qwen3VLModel.get_rope_index(fake, ids, mm_token_type_ids=types, image_grid_thw=grid,
                                         attention_mask=mask)
    return pos


def run_rope(cli, tmp, tmpl_ids, tmpl_types, grid1, ids, mask, types):
    B, W = ids.shape
    Lp = tmpl_ids.numel()
    head = torch.tensor([B, W, Lp, 2, *grid1.tolist()], dtype=torch.int64)
    blob = torch.cat([head, tmpl_ids, tmpl_types, ids.flatten(), mask.flatten(), types.flatten()]).numpy()
    fin, fout = os.path.join(tmp, "rope.in"), os.path.join(tmp, "rope.out")
    blob.tofile(fin)
    subprocess.run([cli, "rope", fin, fout], check=True)
    out = torch.from_numpy(np.fromfile(fout, dtype=np.int64))
    n = 3 * B * W
    return int(out[0]), out[1:1 + n].view(3, B, W), out[1 + n:].view(3, B, W)


def test_rope(proc, cli, tmp, rng):
    msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": PROMPT}]}]
    text = proc.apply_chat_template(msgs, add_generation_prompt=True, tokenize=False)
    img = Image.fromarray(rng.integers(0, 256, (512, 512, 3), dtype=np.uint8))
    x = proc(text=[text], images=[img], return_tensors="pt")
    tids, ttypes, grid1 = x["input_ids"][0], x["mm_token_type_ids"][0], x["image_grid_thw"][0]
    Lp, pad = tids.numel(), proc.tokenizer.pad_token_id
    B, maxpad, T = 16, 7, 1536
    P = Lp + maxpad
    W = P + T
    ids = torch.full((B, W), pad, dtype=torch.int64)
    mask = torch.zeros((B, W), dtype=torch.int64)
    types = torch.zeros((B, W), dtype=torch.int64)
    for b in range(B):
        lead = maxpad if b < 2 else int(rng.integers(0, maxpad + 1))  # TRL left-pads prompts to the batch's longest
        L = 0 if b == 0 else (T if b == 1 else int(rng.integers(1, T + 1)))
        ids[b, lead:lead + Lp], types[b, lead:lead + Lp], mask[b, lead:lead + Lp] = tids, ttypes, 1
        ids[b, P:P + L] = torch.from_numpy(rng.integers(0, 151643, L))
        mask[b, P:P + L] = 1
    mask[2, P + 3] = 0  # a hole inside an answer
    hf = hf_rope(ids, mask, types, grid1.repeat(B, 1))
    err, emul, ref = run_rope(cli, tmp, tids, ttypes, grid1, ids, mask, types)
    ok = err == 0 and torch.equal(emul, hf) and torch.equal(ref, hf)
    print(f"(c) rope: Lp {Lp}, grid {grid1.tolist()}, {B} rows x {W}: err {err}, template+extrapolation == HF "
          f"{torch.equal(emul, hf)}, C++ port == HF {torch.equal(ref, hf)}")
    # adversarial rows must be flagged
    flags = []
    ids2 = ids.clone(); ids2[5, P - 3] += 1  # a prompt token differs (last prompt tokens sit before P)
    flags.append(run_rope(cli, tmp, tids, ttypes, grid1, ids2, mask, types)[0] & 1)
    ty2 = types.clone(); ty2[1, W - 1] = 1  # an image token in an answer
    flags.append(run_rope(cli, tmp, tids, ttypes, grid1, ids, mask, ty2)[0] & 2)
    m2 = mask.clone(); m2[4] = 0  # an all-padding row
    flags.append(run_rope(cli, tmp, tids, ttypes, grid1, ids, m2, types)[0] & 4)
    print(f"(c) rope adversarial flags (mismatch, image in answer, short row): {[bool(f) for f in flags]}")
    return ok and all(flags)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cli", required=True)
    ap.add_argument("--tmp", required=True)
    a = ap.parse_args()
    os.makedirs(a.tmp, exist_ok=True)
    proc = AutoProcessor.from_pretrained(MODEL)
    rng = np.random.default_rng(0)
    ok = test_preprocess(proc, a.cli, a.tmp, rng)
    ok = test_rope(proc, a.cli, a.tmp, rng) and ok
    print("test_hf_cpu:", "all passed" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
