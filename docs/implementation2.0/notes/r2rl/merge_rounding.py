"""Share of LoRA delta lost when W + d is rounded to bf16 (what lora_merged_params pushes to vLLM).
W proxy: the final merged model (W_base + d_770; |d| << |W| so the bf16 exponent of W is the same for almost all
elements). Usage: python -I merge_rounding.py REPO"""
import sys, torch
from safetensors import safe_open
repo = sys.argv[1]
W = safe_open(f"{repo}/runs/rl/train_rl/final/model.safetensors", "pt")
scale = 32 / 16
layers = [0, 9, 18, 27]
mods = ["self_attn.q_proj", "self_attn.v_proj", "mlp.down_proj", "mlp.up_proj"]
print("step  lost_share(|d|<ulp/2)  rel_err ||round(W+d)-W-d||/||d||   median|d|/ulp")
for step in [50, 100, 300, 500, 770]:
    A = safe_open(f"{repo}/runs/rl/train_rl/trainer/checkpoint-{step}/adapter_model.safetensors", "pt")
    lost = tot = 0; num = den = 0.0; ratios = []
    for l in layers:
        for m in mods:
            p = f"base_model.model.model.language_model.layers.{l}.{m}"
            a = A.get_tensor(p + ".lora_A.weight").float(); b = A.get_tensor(p + ".lora_B.weight").float()
            d = (b @ a) * scale
            w = W.get_tensor(f"model.language_model.layers.{l}.{m}.weight").float()
            wb = w.to(torch.bfloat16).float()
            ulp = (torch.nextafter(wb.to(torch.bfloat16), torch.tensor(float('inf'), dtype=torch.bfloat16)).float() - wb).abs()
            merged = (wb + d).to(torch.bfloat16).float()
            err = merged - wb - d
            lost += (d.abs() < ulp / 2).sum().item(); tot += d.numel()
            num += err.pow(2).sum().item(); den += d.pow(2).sum().item()
            ratios.append((d.abs() / ulp).flatten()[::97])
    r = torch.cat(ratios).median().item()
    print(f"{step:4d}  {lost/tot:.3f}                 {(num/den)**0.5:.3f}                                {r:.3f}")
