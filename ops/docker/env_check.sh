#!/usr/bin/env bash
# perf/implementation.md C4: start-of-rental environment check, run ON THE BOX (RTX 5090 first; other cards record only).
#   env_check.sh MODEL_DIR [OUT_JSON]   (default /workspace/runs/env_check.json; exit 1 unless every check passes)
# Checks: nvidia_smi (name, driver, PCIe gen/width, power limit), capability ((12, 0) and sm_120 on a 5090),
# bf16_matmul (>= 0.6 of the dense BF16 peak, as image_check.sh), triton_compile, torch_compile (small MLP fwd+bwd),
# vllm_engine (Qwen3-VL-2B at MODEL_DIR, CUDA graphs on, VLLM_USE_FLASHINFER_SAMPLER=0, attention backend logged),
# nvcc (--list-gpu-arch, when nvcc exists), trl_private_names (objdel.rl.trl_adapter.check(), E5).
set -euo pipefail
MODEL=${1:?usage: env_check.sh MODEL_DIR [OUT_JSON]}
OUT=${2:-/workspace/runs/env_check.json}
mkdir -p "$(dirname "$OUT")"
LOG="${OUT%.json}_vllm.log"
rm -f "$OUT.vllm"
CHECKS=$(mktemp -d)/checks.py  # a file, not stdin: Triton's @jit reads the kernel's source
cat > "$CHECKS" <<'PY'
import json, shutil, subprocess, sys, time
r, bad = {}, []

def check(key):
    def wrap(fn):
        try:
            r[key] = fn()
        except Exception as e:  # a check that crashes is a failed check
            r[key] = {"error": f"{type(e).__name__}: {e}"}
            bad.append(key)
        return fn
    return wrap

@check("nvidia_smi")
def _():
    q = "name,driver_version,pcie.link.gen.current,pcie.link.width.current,power.limit"
    row = subprocess.check_output(["nvidia-smi", f"--query-gpu={q}", "--format=csv,noheader"], text=True).split(",")
    return dict(zip(("name", "driver_version", "pcie_gen", "pcie_width", "power_limit"), (x.strip() for x in row)))

import torch
is5090 = "5090" in r["nvidia_smi"].get("name", "")

@check("capability")
def _():
    cap, archs = torch.cuda.get_device_capability(0), torch.cuda.get_arch_list()
    ok = (cap == (12, 0) and "sm_120" in archs) or not is5090
    assert ok, f"5090 needs capability (12, 0) and sm_120 in the torch build; got {cap}, {archs}"
    return {"capability": list(cap), "arch_list": archs}

@check("bf16_matmul")
def _():
    peak = {"3090": 71, "4090": 165.2, "5090": 209.5}.get(next((k for k in ("3090", "4090", "5090")
                                                               if k in r["nvidia_smi"]["name"]), ""))
    a = torch.randn(8192, 8192, device="cuda", dtype=torch.bfloat16); b = torch.randn_like(a)
    for _ in range(3):
        a @ b
    torch.cuda.synchronize(); t = time.time()
    for _ in range(20):
        a @ b
    torch.cuda.synchronize()
    tf = 2 * 8192**3 / ((time.time() - t) / 20) / 1e12
    assert not peak or tf >= 0.6 * peak, f"{tf:.1f} TFLOPS is under 0.6 of the {peak} peak"
    return {"tflops": round(tf, 1), "peak": peak}

@check("triton_compile")
def _():
    import triton, triton.language as tl

    @triton.jit
    def add(x, y, o, n, B: tl.constexpr):
        i = tl.program_id(0) * B + tl.arange(0, B)
        tl.store(o + i, tl.load(x + i, mask=i < n) + tl.load(y + i, mask=i < n), mask=i < n)
    x = torch.ones(1000, device="cuda"); o = torch.empty_like(x)
    add[(1,)](x, x, o, 1000, B=1024)
    assert bool((o == 2).all())
    return {"triton": triton.__version__}

@check("torch_compile")
def _():
    m = torch.compile(torch.nn.Sequential(torch.nn.Linear(256, 256), torch.nn.GELU(), torch.nn.Linear(256, 8)).cuda())
    m(torch.randn(32, 256, device="cuda")).sum().backward()
    torch.cuda.synchronize()
    return {"torch": torch.__version__}

@check("nvcc")
def _():
    if not shutil.which("nvcc"):
        return {"present": False}  # the toolkit is installed at run time by run_box.sh
    out = subprocess.check_output(["nvcc", "--list-gpu-arch"], text=True).split()
    assert not is5090 or "compute_120" in out, f"nvcc has no compute_120: {out}"
    return {"present": True, "gpu_arch": out}

@check("trl_private_names")
def _():
    import trl
    from objdel.rl import trl_adapter
    trl_adapter.check()  # a TRL change that breaks the adapter fails the rental here, before any training (E-G13)
    return {"trl": trl.__version__, "items": len(trl_adapter.ALL_ITEMS)}

json.dump({"bad": bad, **r}, open(sys.argv[1], "w"), indent=1)
print(json.dumps(r, indent=1))
PY
/opt/venv/bin/python -I "$CHECKS" "$OUT"
if [ -x /opt/vllm/bin/python ]; then VLLM_PY=/opt/vllm/bin/python; else VLLM_PY=/opt/venv/bin/python; fi
VLLM_USE_FLASHINFER_SAMPLER=0 VLLM_LOGGING_LEVEL=INFO "$VLLM_PY" -I - "$MODEL" "$OUT.vllm" <<'PY' >"$LOG" 2>&1 || true
import json, sys
from vllm import LLM, SamplingParams
llm = LLM(sys.argv[1], max_model_len=2048, gpu_memory_utilization=0.5, enforce_eager=False)  # CUDA graphs on
out = llm.generate(["Hello"], SamplingParams(max_tokens=8))[0].outputs[0].text
json.dump({"generated": out}, open(sys.argv[2], "w"))
PY
/opt/venv/bin/python -I - "$OUT" "$LOG" <<'PY'
import json, os, sys
r = json.load(open(sys.argv[1]))
log = open(sys.argv[2], errors="replace").read().splitlines()
backend = next((ln.strip() for ln in log if "backend" in ln.lower() and "attention" in ln.lower()), "not logged")
v = sys.argv[1] + ".vllm"
r["vllm_engine"] = {"started": os.path.exists(v), "attention_backend": backend, "log": sys.argv[2],
                    "tail": log[-5:] if not os.path.exists(v) else []}
if not os.path.exists(v):
    r["bad"].append("vllm_engine")
r["pass"] = not r["bad"]
json.dump(r, open(sys.argv[1], "w"), indent=1)
print("vllm_engine", r["vllm_engine"]["started"], "| attention:", backend, "| bad:", r["bad"] or "none")
sys.exit(0 if r["pass"] else 1)
PY
