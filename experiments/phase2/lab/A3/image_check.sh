#!/usr/bin/env bash
# Runs INSIDE the candidate Vast image (Dockerfile.vast) with /io mounted (src.tar.gz in, logs out). Lab evidence for
# A3_vast_gpus.md, hypothesis H-IMG. On the Mac this runs under linux/amd64 emulation without a GPU: it checks that the
# stack imports, which CUDA build and GPU architectures torch carries, that objdel builds from the source snapshot
# and that the CPU reward path hashes like the validated x86 CPU image (experiments/phase1/bench). On a Vast box the
# same script is the first-session smoke test (it then also reports the GPU, driver and a bf16 matmul).
set -uo pipefail
IO=${IO:-/io}
mkdir -p /work && cd /work && tar -xzf "$IO/src.tar.gz" -C /work
{
    echo "uname: $(uname -m)  nproc: $(nproc)  cgroup cpu.max: $(cat /sys/fs/cgroup/cpu.max 2>/dev/null)  mem.max: $(cat /sys/fs/cgroup/memory.max 2>/dev/null)"
    grep -m1 'model name' /proc/cpuinfo
    echo "image size of /opt/venv: $(du -sh /opt/venv | cut -f1)"
} > "$IO/env.txt" 2>&1
t=$SECONDS
uv pip install --python /opt/venv/bin/python --offline --no-deps --no-build-isolation -e . > "$IO/build.log" 2>&1
echo "editable build rc=$? in $((SECONDS - t)) s" >> "$IO/env.txt"
python - > "$IO/stack.txt" 2>&1 <<'EOF'
import importlib, time
t0 = time.time()
import torch
print("torch", torch.__version__, "cuda build", torch.version.cuda, "cudnn", torch.backends.cudnn.version() if torch.backends.cudnn.is_available() else None)
print("arch flags", torch._C._cuda_getArchFlags())
print("cuda available", torch.cuda.is_available())
if torch.cuda.is_available():
    p = torch.cuda.get_device_properties(0)
    print("gpu", p.name, f"{p.total_memory/2**30:.1f} GiB", "cc", f"{p.major}.{p.minor}", "bf16", torch.cuda.is_bf16_supported())
    a = torch.randn(8192, 8192, device="cuda", dtype=torch.bfloat16); b = torch.randn_like(a)
    for _ in range(3): a @ b
    torch.cuda.synchronize(); t = time.time(); n = 20
    for _ in range(n): a @ b
    torch.cuda.synchronize(); dt = (time.time() - t) / n
    print(f"bf16 8192^3 matmul: {2*8192**3/dt/1e12:.1f} TFLOPS")
    print("sdpa flash available", torch.backends.cuda.flash_sdp_enabled(), "cudnn sdpa", torch.backends.cuda.cudnn_sdp_enabled())
for m in ["transformers", "trl", "peft", "accelerate", "datasets", "numpy", "scipy", "cv2", "numba", "objdel", "objdel._polycore"]:
    mod = importlib.import_module(m)
    print(m, getattr(mod, "__version__", ""), getattr(mod, "__file__", ""))
print(f"imports {time.time()-t0:.1f} s")
EOF
echo "stack rc=$?" >> "$IO/env.txt"
t=$SECONDS
python -I -m experiments.phase1.bench.xarch_hashes --out "$IO/hashes-vastimg.json" > "$IO/hashes.log" 2>&1
echo "hashes rc=$? in $((SECONDS - t)) s" >> "$IO/env.txt"
uv pip freeze --python /opt/venv/bin/python > "$IO/freeze.txt" 2>&1
cat "$IO/env.txt" "$IO/stack.txt"
