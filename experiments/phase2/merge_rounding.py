"""B2.0 merge-rounding table ($0, CPU): how much of each LoRA delta survives in the weights vLLM receives.

    python -I -m experiments.phase2.merge_rounding --run runs/rl/train_rl --out RD/merge_rounding.json

For every checkpoint step, over all text layers and every LoRA module of `<run>/trainer/checkpoint-N`, the pushed
weight is `W + d.to(dtype)` in `dtype`, as fast.lora_merged_params computes it, with d = (B @ A) x alpha / r in fp32.
W proxy: `<run>/final/model.safetensors` (the merged final model; |d| << |W|, so W's exponent is the same).
- erased_share: non-zero delta elements whose pushed weight equals W (exact rounding: ties to even and the half
  spacing below a power of two are handled; 0 when no element is non-zero);
- rel_err: ||pushed - W - d|| / ||d|| (0 when d = 0);
- median_d_over_ulp: median of |d| / ulp(W) over non-zero d (ulp = the spacing above |W|; at most 50,000 evenly
  spaced elements per module).
B2.4 rule 0 reads the step-100 bf16 rel_err.
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import torch
from safetensors import safe_open

from objdel.core.io import read_json

SAMPLE = 50_000  # elements per module for the median
DTYPES = {"bf16": torch.bfloat16, "fp16": torch.float16}


def table_row(pairs, dtype: torch.dtype) -> dict:
    """pairs: (W, d) with W in any float dtype and d in fp32; the three table numbers for pushes in `dtype`."""
    erased = nonzero = 0
    num = den = 0.0
    ratios = []
    for w, d in pairs:
        w = w.to(dtype)
        pushed = (w + d.to(dtype)).double()
        wd, dd = w.double(), d.double()
        nz = dd != 0
        erased += int((nz & (pushed == wd)).sum())
        nonzero += int(nz.sum())
        num += float((pushed - wd - dd).pow(2).sum())
        den += float(dd.pow(2).sum())
        a = w.abs()
        ulp = (torch.nextafter(a, torch.full_like(a, torch.inf)) - a).double()
        r = (dd.abs() / ulp)[nz]
        ratios.append(r[:: max(1, r.numel() // SAMPLE)].numpy())
    r = np.concatenate(ratios) if ratios else np.zeros(0)
    return {"erased_share": erased / nonzero if nonzero else 0.0, "rel_err": (num / den) ** 0.5 if den else 0.0,
            "median_d_over_ulp": float(np.median(r)) if r.size else 0.0}  # fmt: skip


def pairs(model, adapter, scale: float):
    """(W, d) for every LoRA module of the adapter, in the adapter file's key order."""
    for k in adapter.keys():
        if k.endswith(".lora_A.weight"):
            p = k.removesuffix(".lora_A.weight")
            d = (adapter.get_tensor(p + ".lora_B.weight").float() @ adapter.get_tensor(k).float()) * scale
            yield model.get_tensor(p.removeprefix("base_model.model.") + ".weight"), d


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", type=Path, required=True, help="train_rl output: trainer/checkpoint-N and final/")
    ap.add_argument("--steps", type=int, nargs="+", default=[50, 100, 300, 500, 750])
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    model = safe_open(a.run / "final/model.safetensors", "pt")
    rows, names = [], set()
    for step in a.steps:
        ck = a.run / f"trainer/checkpoint-{step}"
        cfg = read_json(ck / "adapter_config.json")
        adapter = safe_open(ck / "adapter_model.safetensors", "pt")
        names |= {k.removesuffix(".lora_A.weight") for k in adapter.keys() if k.endswith(".lora_A.weight")}
        scale = cfg["lora_alpha"] / cfg["r"]
        rows.append({"step": step} | {k: table_row(pairs(model, adapter, scale), dt) for k, dt in DTYPES.items()})
    layers = len({re.search(r"layers\.(\d+)\.", n)[1] for n in names})
    modules = sorted({n.rsplit(".", 1)[1] for n in names})
    doc = {"run": str(a.run), "layers": layers, "modules": modules, "rows": rows}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(doc, indent=1) + "\n")
    print(f"{a.run}: {layers} layers x {len(modules)} modules ({', '.join(modules)})")
    print("step  erased_bf16  rel_err_bf16  rel_err_fp16  median|d|/ulp_bf16")
    for r in rows:
        b, f = r["bf16"], r["fp16"]
        print(f"{r['step']:4d}  {b['erased_share']:11.3f}  {b['rel_err']:12.3f}  {f['rel_err']:12.3f}  "
              f"{b['median_d_over_ulp']:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
