#!/usr/bin/env bash
# PLAN §4.7 step 3 (A3 §8.3, D15): first-session image check, run ON THE BOX after sync_up (code in /workspace/src).
#   image_check.sh [OUT_JSON]   (default /workspace/runs/image_check.json; exit 1 unless every check passes)
# Records: CUDA available; bf16 8192^3 matmul TFLOPS vs the card's dense BF16 peak (pass >= 0.6; A3 §3.1); the
# Qwen3.5 Gated-DeltaNet route per function (hub `kernels` / package fla, causal_conv1d / reference) as transformers
# binds it at import (a reference route is recorded, not a failure: §4.7 step 8 prices it); objdel._polycore and the
# /opt/vllm sampler imports; the size of /dev/shm (BK4-19). The use_kernels=True hub route resolves only at model load:
# R0 records that one.
set -euo pipefail
OUT=${1:-/workspace/runs/image_check.json}
mkdir -p "$(dirname "$OUT")"
/opt/venv/bin/python -I - "$OUT" <<'PY'
import importlib, importlib.util, json, re, shutil, sys, time
import torch
r = {"torch": torch.__version__, "cuda_build": torch.version.cuda, "arch": torch._C._cuda_getArchFlags(),
     "cuda": torch.cuda.is_available()}
PEAK = {"3090": 71, "4090": 165.2, "5090": 209.5, "PRO 6000": 500, "A100": 312, "L40S": 362, "H100": 989}
if r["cuda"]:
    name = torch.cuda.get_device_name(0)
    a = torch.randn(8192, 8192, device="cuda", dtype=torch.bfloat16); b = torch.randn_like(a)
    for _ in range(3):
        a @ b
    torch.cuda.synchronize(); t = time.time()
    for _ in range(20):
        a @ b
    torch.cuda.synchronize()
    tf = 2 * 8192**3 / ((time.time() - t) / 20) / 1e12
    peak = next((v for k, v in PEAK.items() if k in name), None)
    r.update(gpu=name, bf16_tflops=round(tf, 1), peak=peak, matmul_ok=bool(peak and tf >= 0.6 * peak))
FUNCS = ["chunk_gated_delta_rule", "fused_recurrent_gated_delta_rule", "causal_conv1d_fn", "causal_conv1d_update"]
r["packages"] = {m: importlib.util.find_spec(m) is not None for m in ("fla", "causal_conv1d", "kernels")}
routes = {}
for modname in ("transformers.models.qwen3_5.modeling_qwen3_5", "transformers.models.qwen3_next.modeling_qwen3_next"):
    try:
        mod = importlib.import_module(modname)
    except ImportError as e:
        routes[modname] = f"import failed: {e}"
        continue
    for f in FUNCS:
        fn = getattr(mod, f, None)
        src = getattr(fn, "__module__", "") or ""
        route = ("reference" if fn is None or src.startswith("transformers") else
                 "hub" if re.search(r"kernels|hub", src) else "package" if re.match(r"(fla|causal_conv1d)", src)
                 else "other")
        routes[f"{modname.rsplit('.', 1)[-1]}.{f}"] = {"route": route, "module": src}
    break
r["gdn_routes"] = routes
importlib.import_module("objdel._polycore")
r["polycore"] = True
r["shm_bytes"] = shutil.disk_usage("/dev/shm").total  # RIGA decrypts into /dev/shm (~280 MB; Docker default 64 MB)
r["pass"] = bool(r["cuda"] and r.get("matmul_ok"))
json.dump(r, open(sys.argv[1], "w"), indent=1)
print(json.dumps(r, indent=1))
PY
/opt/vllm/bin/python -I - "$OUT" <<'PY'
import json, sys
from vllm import LLM, SamplingParams  # noqa: F401
import vllm
r = json.load(open(sys.argv[1]))
r.update(vllm=vllm.__version__, vllm_import=True)
json.dump(r, open(sys.argv[1], "w"), indent=1)
print("vllm", vllm.__version__, "pass", r["pass"])
sys.exit(0 if r["pass"] else 1)
PY
