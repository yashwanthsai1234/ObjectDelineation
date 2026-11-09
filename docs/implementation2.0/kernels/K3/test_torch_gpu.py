"""K3 on the GPU against the production references (Hugging Face, TRL, objdel.train.fast), plus timings.

(a) preprocess (GPU) vs the HF Qwen3-VL fast image processor (CPU, as ImageMemo calls it): bit equality, fp32 and bf16,
    32 random + 6 adversarial images, and the 256-row generation-batch layout.
(c) rope (GPU) vs HF Qwen3VLModel.get_rope_index on a real prompt: equality on 256 rows; error flags.
(b) lse (GPU) vs fp64 torch; tensor-core and exact modes; TRL's ChunkedLogProbFunction measured alongside.
(b) is_ratio (GPU, libdevice exp) vs objdel.train.fast.is_ratio on CUDA: bit equality for the token modes; the
    sequence modes sum in another order (max relative difference reported).
Timings use CUDA events (median). Results go to --out (JSON).

Run on the GPU box from this folder (REPO's venv; the model's processor files in the HF cache):
    TORCH_CUDA_ARCH_LIST=12.0 HF_HUB_OFFLINE=1 REPO/.venv/bin/python -I test_torch_gpu.py --out build/torch_gpu.json
"""

import argparse
import json
import math
import os
import sys
import time
import types as pytypes

import numpy as np
import torch
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = "Qwen/Qwen3-VL-2B-Instruct"
PROMPT = ('Outline the lung nodule in this CT crop. Answer with JSON only: '
          '{"nodule":[{"outer":[[x,y],...],"holes":[[[x,y],...]]}]}, '
          'integers 0-1000 relative to the crop, outer rings clockwise, holes counter-clockwise.')  # dataset.PROMPTS
RESULTS = {}


def build():
    from torch.utils.cpp_extension import load

    t = time.perf_counter()
    bdir = os.path.join(HERE, "build", "ext")
    os.makedirs(bdir, exist_ok=True)
    ext = load(name="k3_ext", sources=[os.path.join(HERE, "k3_torch.cpp"), os.path.join(HERE, "k3_kernels.cu")],
               extra_cflags=["-O3", "-ffp-contract=off"],
               extra_cuda_cflags=["-O3", "--fmad=false", "-lineinfo", "-U__CUDA_NO_BFLOAT16_CONVERSIONS__",
                                  "-U__CUDA_NO_HALF_CONVERSIONS__", "-U__CUDA_NO_HALF_OPERATORS__"], build_directory=bdir,
               verbose=False)
    RESULTS["build_s"] = round(time.perf_counter() - t, 1)
    return ext


def cuda_ms(fn, reps=10):
    for _ in range(2):
        fn()
    torch.cuda.synchronize()
    ts = []
    for _ in range(reps):
        a, b = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
        a.record()
        fn()
        b.record()
        b.synchronize()
        ts.append(a.elapsed_time(b))
    return float(np.median(ts))


def wall_ms(fn, reps=3):
    fn()
    ts = []
    for _ in range(reps):
        torch.cuda.synchronize()
        t = time.perf_counter()
        fn()
        torch.cuda.synchronize()
        ts.append((time.perf_counter() - t) * 1e3)
    return float(np.median(ts))


def images(rng):
    out = [rng.integers(0, 256, (512, 512, 3), dtype=np.uint8) for _ in range(32)]
    out += [np.zeros((512, 512, 3), np.uint8), np.full((512, 512, 3), 255, np.uint8)]
    i = np.arange(512 * 512 * 3, dtype=np.int64).reshape(512, 512, 3)
    out.append(((i * 7 + i // 3) & 255).astype(np.uint8))
    # digit images: byte k of each source element's flat index. The value map is injective, so equal outputs on
    # these 3 images prove the index map (patch layout + temporal copy) equal at every output position
    out += [((i >> (8 * k)) & 255).astype(np.uint8) for k in range(3)]
    yy, xx = np.mgrid[:512, :512]
    out.append(np.repeat((((yy // 16 + xx // 16) & 1) * 255).astype(np.uint8)[..., None], 3, 2))
    out.append(np.stack([xx & 255, yy & 255, (xx + yy) & 255], -1).astype(np.uint8))
    g = rng.integers(0, 256, (512, 512), dtype=np.uint8)
    out.append(np.repeat(g[..., None], 3, 2))
    return out


def test_preprocess(ext, proc, rng):
    ip = proc.image_processor
    imgs = images(rng)
    pil = [Image.fromarray(a) for a in imgs]
    hf = ip(pil, return_tensors="pt")["pixel_values"]
    mean, std, _ = ip._fuse_mean_std_and_rescale_factor(True, tuple(ip.image_mean), tuple(ip.image_std), True,
                                                        ip.rescale_factor, torch.device("cpu"))
    mean, std = [float(v) for v in mean], [float(v) for v in std]
    src = torch.from_numpy(np.stack(imgs)).cuda()
    ours = ext.preprocess(src, None, mean, std, True, False)
    ours16 = ext.preprocess(src, None, mean, std, True, True)
    eq32 = torch.equal(ours.cpu().view(torch.int32), hf.view(torch.int32))
    eq16 = torch.equal(ours16.cpu().view(torch.int16), hf.to(torch.bfloat16).view(torch.int16))
    # generation batch: 256 rows, 8 per image (ImageMemo's assembly: per-image slices concatenated in row order)
    idx = torch.arange(256, dtype=torch.int32) // 8
    rows = ext.preprocess(src[:32].contiguous(), idx.cuda(), mean, std, True, False)
    memo = torch.cat([hf[i * 1024:(i + 1) * 1024] for i in idx.tolist()])
    eqrows = torch.equal(rows.cpu().view(torch.int32), memo.view(torch.int32))
    print(f"(a) preprocess vs HF: fp32 {eq32}, bf16 {eq16}, 256-row layout {eqrows}")
    # the library alternative: the same HF processor run on the GPU (device="cuda")
    hf_gpu = ip(pil, return_tensors="pt", device="cuda")["pixel_values"]
    eq_hfgpu = torch.equal(hf_gpu.cpu().view(torch.int32), hf.view(torch.int32))
    t_hfgpu = wall_ms(lambda: ip(pil[:32], return_tensors="pt", device="cuda"))
    print(f"    HF processor with device='cuda': bit-equal to HF CPU {eq_hfgpu}, {t_hfgpu:.1f} ms for 32 images")
    # timings: HF on CPU (32 images, as ImageMemo's miss path) + the CPU row assembly, vs H2D + kernel
    t_hf = wall_ms(lambda: ip(pil[:32], return_tensors="pt"))
    t_cat = wall_ms(lambda: torch.cat([hf[i * 1024:(i + 1) * 1024] for i in idx.tolist()]))
    pinned = torch.from_numpy(np.stack(imgs[:32])).pin_memory()
    dst = torch.empty_like(pinned, device="cuda")
    t_h2d = cuda_ms(lambda: dst.copy_(pinned, non_blocking=True))
    t_k32 = cuda_ms(lambda: ext.preprocess(dst, None, mean, std, True, False))
    didx = idx.cuda()
    t_rows16 = cuda_ms(lambda: ext.preprocess(dst, didx, mean, std, True, True))
    t_rows32 = cuda_ms(lambda: ext.preprocess(dst, didx, mean, std, True, False))
    RESULTS["preprocess"] = dict(eq_fp32=eq32, eq_bf16=eq16, eq_rows=eqrows, hf_cpu_32img_ms=t_hf,
                                 hf_gpu_eq_hf_cpu=eq_hfgpu, hf_gpu_32img_wall_ms=t_hfgpu,
                                 cpu_row_cat_ms=t_cat, h2d_32img_ms=t_h2d, gpu_32img_fp32_ms=t_k32,
                                 gpu_256rows_bf16_ms=t_rows16, gpu_256rows_fp32_ms=t_rows32)
    print(f"    HF CPU 32 images {t_hf:.1f} ms + CPU 256-row cat {t_cat:.1f} ms | GPU: H2D {t_h2d:.2f} ms, "
          f"32 images fp32 {t_k32:.3f} ms, 256 rows bf16 {t_rows16:.3f} ms / fp32 {t_rows32:.3f} ms")
    return eq32 and eq16 and eqrows


def hf_rope(ids, mask, types, grid):
    from transformers.models.qwen3_vl.modeling_qwen3_vl import Qwen3VLModel

    fake = pytypes.SimpleNamespace(
        config=pytypes.SimpleNamespace(vision_config=pytypes.SimpleNamespace(spatial_merge_size=2)),
        get_vision_position_ids=lambda *a, **k: Qwen3VLModel.get_vision_position_ids(None, *a, **k))
    return Qwen3VLModel.get_rope_index(fake, ids, mm_token_type_ids=types, image_grid_thw=grid,
                                       attention_mask=mask)[0]


def rope_torch(ids, mask, tpos, tend):
    """The same template + extrapolation with plain torch ops on the GPU (no custom kernel, no validation)."""
    Lp = tpos.size(1)
    r = mask.cumsum(1) - 1
    head = tpos[:, r.clamp(0, Lp - 1)]  # (3, B, W)
    pos = torch.where(r < Lp, head, (tend + r - Lp).expand(3, -1, -1))
    return torch.where(mask.bool(), pos, 0)


def test_rope(ext, proc, rng):
    msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": PROMPT}]}]
    text = proc.apply_chat_template(msgs, add_generation_prompt=True, tokenize=False)
    x = proc(text=[text], images=[Image.fromarray(rng.integers(0, 256, (512, 512, 3), dtype=np.uint8))],
             return_tensors="pt")
    tids, ttypes, grid1 = x["input_ids"][0], x["mm_token_type_ids"][0], x["image_grid_thw"][0]
    # the template, once: HF on the bare prompt
    tpos = hf_rope(tids[None], torch.ones_like(tids)[None], ttypes[None], grid1[None])[:, 0]
    tend = int(tpos.max()) + 1  # current_pos after the prompt (positions are increasing groups)
    Lp, pad, B, maxpad, T = tids.numel(), proc.tokenizer.pad_token_id, 256, 7, 1536
    P, W = Lp + maxpad, Lp + maxpad + T
    ids = torch.full((B, W), pad, dtype=torch.int64)
    mask = torch.zeros((B, W), dtype=torch.int64)
    types = torch.zeros((B, W), dtype=torch.int64)
    for b in range(B):
        lead = maxpad if b < 2 else int(rng.integers(0, maxpad + 1))
        L = 0 if b == 0 else (T if b == 1 else int(rng.integers(1, T + 1)))
        ids[b, lead:lead + Lp], types[b, lead:lead + Lp], mask[b, lead:lead + Lp] = tids, ttypes, 1
        ids[b, P:P + L] = torch.from_numpy(rng.integers(0, 151643, L))
        mask[b, P:P + L] = 1
    mask[2, P + 3] = 0
    t = time.perf_counter()
    ref = hf_rope(ids, mask, types, grid1.repeat(B, 1))
    t_hf = (time.perf_counter() - t) * 1e3
    g = [v.cuda() for v in (ids, mask, types, tids, ttypes, tpos.contiguous())]
    pos, err = ext.rope(*g, tend)
    eq = torch.equal(pos.cpu(), ref) and int(err) == 0
    flags = []
    i2 = g[0].clone(); i2[5, P - 3] += 1
    flags.append(int(ext.rope(i2, *g[1:], tend)[1]) & 1)
    t2 = g[2].clone(); t2[1, W - 1] = 1
    flags.append(int(ext.rope(g[0], g[1], t2, *g[3:], tend)[1]) & 2)
    m2 = g[1].clone(); m2[4] = 0
    flags.append(int(ext.rope(g[0], m2, *g[2:], tend)[1]) & 4)
    t_gpu = cuda_ms(lambda: ext.rope(*g, tend), reps=20)
    eq_torch = torch.equal(rope_torch(g[0], g[1], g[5], tend).cpu(), ref)
    t_torch = cuda_ms(lambda: rope_torch(g[0], g[1], g[5], tend), reps=20)
    RESULTS["rope"] = dict(eq=eq, flags=[bool(f) for f in flags], hf_cpu_256rows_ms=t_hf, gpu_256rows_ms=t_gpu,
                           torch_ops_eq=eq_torch, torch_ops_ms=t_torch,
                           Lp=Lp, W=W)
    print(f"(c) rope vs HF get_rope_index ({B} rows x {W}, Lp {Lp}): {eq}; flags {[bool(f) for f in flags]}; "
          f"HF CPU {t_hf:.0f} ms vs GPU {t_gpu:.3f} ms; plain torch ops {eq_torch}, {t_torch:.3f} ms")
    return eq and all(flags)


def lse_reference(H, W, tgt, inv_t=1.0):
    """fp64 log-probs and entropies (chunks of 256 tokens)."""
    Wd = W.double()
    lps, ents = [], []
    for s in range(0, H.size(0), 256):
        z = (H[s:s + 256].double() @ Wd.T) * inv_t
        lsm = torch.log_softmax(z, -1)
        t = tgt[s:s + 256].long()
        ok = (t >= 0) & (t < W.size(0))
        lps.append(torch.where(ok, lsm.gather(1, t.clamp(0, W.size(0) - 1)[:, None])[:, 0], 0.0))
        ents.append(-(lsm.exp() * lsm).sum(-1))
    return torch.cat(lps), torch.cat(ents)


def test_lse(ext, rng):
    torch.manual_seed(0)
    V, D = 151936, 2048
    W = (torch.randn(V, D, device="cuda") * 0.05).bfloat16()
    N = 1024
    H = torch.randn(N, D, device="cuda").bfloat16()
    tgt = torch.randint(0, V, (N,), device="cuda", dtype=torch.int32)
    tgt[0], tgt[1], tgt[2] = 0, V - 1, -100
    rl, re = lse_reference(H, W, tgt)
    out = {}
    lp, en = ext.lse(H, W, tgt, 1.0, 0, False, True)
    out["tc"] = ((lp.double() - rl).abs().max().item(), (en.double() - re).abs().max().item())
    lpx, enx = ext.lse(H[:256].contiguous(), W, tgt[:256].contiguous(), 1.0, 0, True, True)
    out["exact"] = ((lpx.double() - rl[:256]).abs().max().item(), (enx.double() - re[:256]).abs().max().item())
    try:
        from trl.kernels import ChunkedLogProbFunction

        def trl_fn(h, t):
            return ChunkedLogProbFunction.apply(h, W, None, t.long().clamp(min=0), 1.0, 32768, None, 1.0,
                                                ("entropy",))[:2]
        with torch.no_grad():
            tl, te = trl_fn(H, tgt)
        tl = torch.where(tgt >= 0, tl, 0.0)
        out["trl"] = ((tl.double() - rl).abs().max().item(), (te.double() - re).abs().max().item())
    except Exception as e:  # noqa: BLE001 - measured alongside, not required
        trl_fn = None
        out["trl_error"] = repr(e)
    naive = torch.log_softmax(torch.nn.functional.linear(H[:256], W).float(), -1)  # bf16 logits (the naive head)
    out["naive_bf16_logits"] = ((naive.gather(1, tgt[:256].long().clamp(min=0)[:, None])[:, 0].double()
                                 - torch.where(tgt[:256] >= 0, rl[:256], 0)).abs()[3:].max().item(),)
    print("(b) lse max |error| vs fp64 (logp, entropy): " + ", ".join(f"{k} {v}" for k, v in out.items()))
    ok = out["tc"][0] < 2e-3 and out["tc"][1] < 2e-3 and out["exact"][0] < 1e-4
    times = {}
    for n in (2048, 8192, 24576):
        Hn = torch.randn(n, D, device="cuda").bfloat16()
        tn = torch.randint(0, V, (n,), device="cuda", dtype=torch.int32)
        times[f"k3_tc_N{n}_ms"] = cuda_ms(lambda h=Hn, t=tn: ext.lse(h, W, t, 1.0, 0, False, True), reps=5)
        if trl_fn is not None:
            with torch.no_grad():
                times[f"trl_N{n}_ms"] = cuda_ms(lambda h=Hn, t=tn: trl_fn(h, t), reps=5)
        torch.cuda.reset_peak_memory_stats()
        base = torch.cuda.memory_allocated()
        ext.lse(Hn, W, tn, 1.0, 0, False, True)
        torch.cuda.synchronize()
        times[f"k3_tc_N{n}_extra_MB"] = (torch.cuda.max_memory_allocated() - base) / 1e6
        if trl_fn is not None:
            torch.cuda.reset_peak_memory_stats()
            with torch.no_grad():
                trl_fn(Hn, tn)
            torch.cuda.synchronize()
            times[f"trl_N{n}_extra_MB"] = (torch.cuda.max_memory_allocated() - base) / 1e6
        del Hn, tn
    print("    timings / extra memory: " + ", ".join(f"{k} {v:.2f}" for k, v in times.items()))
    RESULTS["lse"] = dict(err=out, ok=ok, **times)
    return ok


def ref_is_ratio():
    try:
        from objdel.train.fast import is_ratio

        return is_ratio, "objdel.train.fast.is_ratio"
    except Exception:  # noqa: BLE001 - verbatim copy of REPO/src/objdel/train/fast.py:351-362
        def is_ratio(logps, sampling, mask, mode, lo, hi):
            diff = torch.nan_to_num((logps - sampling) * mask, nan=0.0)
            if mode.startswith("sequence"):
                diff = diff.sum(-1, keepdim=True)
            r = torch.exp(diff)
            if mode.endswith("truncate"):
                return torch.clamp(r, min=lo, max=hi)
            lo, hi = -math.inf if lo is None else lo, math.inf if hi is None else hi
            return r.masked_fill((r < lo) | (r > hi), 0.0)

        return is_ratio, "verbatim copy of fast.is_ratio"


def test_is(ext, rng):
    fn, src = ref_is_ratio()
    torch.manual_seed(1)
    B, T = 256, 1536
    lp = -(torch.randn(B, T, device="cuda") * 0.05).abs() * 20
    sp = lp + torch.randn(B, T, device="cuda") * 0.05
    mask = torch.rand(B, T, device="cuda") > 0.1
    sp[0, 3] = float("nan"); mask[0, 4] = False; sp[0, 4] = float("nan"); sp[0, 5] = -float("inf")
    lp[0, 6], sp[0, 6] = -80.0, 0.0
    mask[7] = False
    modes = {"token_truncate": 0, "token_mask": 1, "sequence_truncate": 2, "sequence_mask": 3}
    res, ok = {}, True
    for name, (lo, hi) in [("token_truncate", (None, 3.0)), ("token_mask", (0.0, 3.0)),
                           ("token_truncate", (0.5, 2.0)), ("sequence_truncate", (None, 3.0)),
                           ("sequence_mask", (None, 3.0)), ("sequence_mask", (0.0, 3.0))]:
        ref = fn(lp, sp, mask, name, lo, hi)
        ours, stats = ext.is_ratio(lp, sp, mask, modes[name], lo, hi, True)
        # the 7 logged stats as FastGRPO.set_is_ratio computes them (sums in another order: relative tolerance)
        ok_ = mask & ~torch.isnan(sp)
        d = (lp - sp).abs()
        rk = torch.ones_like(ref, dtype=torch.bool) if name.startswith("sequence") else mask
        want = torch.stack([d.masked_fill(~ok_, 0).sum(), ok_.sum(), d.masked_fill(~ok_, -math.inf).amax(),
                            ref.masked_fill(~rk, 0).sum(), rk.sum(), ref.masked_fill(~rk, math.inf).amin(),
                            ref.masked_fill(~rk, -math.inf).amax()]).float()
        srel = ((stats - want).abs() / want.abs().clamp(min=1e-30)).nan_to_num(0.0)  # inf - inf
        stats_ok = bool(((stats == want) | (srel < 1e-4)).all())
        ok = ok and stats_ok
        if name.startswith("token"):
            eq = torch.equal(ours, ref)
            res[f"{name}[{lo},{hi}]"] = ("bit-equal" if eq else f"{int((ours != ref).sum())} differ") + f", stats ok {stats_ok}"
            ok = ok and eq
        else:
            rel = ((ours - ref).abs() / ref.abs().clamp(min=1e-30)).max().item()
            res[f"{name}[{lo},{hi}]"] = f"max rel diff {rel:.2e}, stats ok {stats_ok}"
            ok = ok and rel < 1e-4
    t_ref = cuda_ms(lambda: fn(lp, sp, mask, "token_truncate", None, 3.0), reps=50)
    t_ours = cuda_ms(lambda: ext.is_ratio(lp, sp, mask, 0, None, 3.0, True), reps=50)
    RESULTS["is_ratio"] = dict(reference=src, results=res, ok=ok, torch_ms=t_ref, k3_ms_with_stats=t_ours)
    print(f"(b) is_ratio vs {src}: {res}; torch {t_ref:.3f} ms vs K3 (ratio + 7 stats) {t_ours:.3f} ms")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "build", "torch_gpu.json"))
    a = ap.parse_args()
    from transformers import AutoProcessor

    RESULTS["device"] = torch.cuda.get_device_name()
    RESULTS["torch"] = torch.__version__
    ext = build()
    proc = AutoProcessor.from_pretrained(MODEL)
    rng = np.random.default_rng(0)
    ok = test_preprocess(ext, proc, rng)
    ok = test_rope(ext, proc, rng) and ok
    ok = test_is(ext, rng) and ok
    ok = test_lse(ext, rng) and ok
    RESULTS["all_passed"] = ok
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(RESULTS, f, indent=1)
    print("test_torch_gpu:", "all passed" if ok else "FAILED", "->", a.out)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
