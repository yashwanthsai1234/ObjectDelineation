#!/usr/bin/env bash
# B5: runs INSIDE a candidate image (Dockerfile.phase2 / Dockerfile.phase3), linux/amd64, no GPU, --network none.
# /io holds src.tar.gz (git ls-files snapshot, no data/) and receives the logs. LABEL names the outputs.
#   1. stack: per venv (/opt/venv, /opt/vllm if present) torch / CUDA build / arch list and the imports Phase 2/3 need
#      (vLLM, TRL's vLLM paths, objdel.train); Python warnings are printed (TRL's "supports vLLM 0.22.0 to 0.31.0")
#   2. objdel editable build into /opt/venv (offline, no deps, no build isolation, as x86_inside.sh)
#   3. the fast suite exactly as experiments/phase1/bench/x86_inside.sh runs it (same marker and deselect)
#   4. experiments.phase1.bench.xarch_hashes -> /io/hashes-$LABEL.json (compared on the host, --strict)
set -uo pipefail
IO=/io; LABEL=${LABEL:?}
mkdir -p /work && cd /work && tar -xzf "$IO/src.tar.gz" -C /work
step() {
    local name=$1; shift; local t0=$SECONDS
    "$@" > "$IO/$LABEL-$name.log" 2>&1; local rc=$?
    printf '%s\t%s\t%s\t%ss\n' "$LABEL" "$name" "$rc" "$((SECONDS - t0))" | tee -a "$IO/status.tsv"
}
for v in /opt/venv /opt/vllm; do
    [ -x "$v/bin/python" ] || continue
    step "stack-$(basename $v)" "$v/bin/python" -W always - <<'PY'
import importlib, importlib.util, sys, time
t0 = time.time()
import torch
print(sys.prefix, "torch", torch.__version__, "cuda build", torch.version.cuda, "arch", torch._C._cuda_getArchFlags(),
      "cuda available", torch.cuda.is_available())
mods = ["transformers", "numpy", "numba", "PIL"]
if importlib.util.find_spec("trl"):
    mods += ["trl", "peft", "accelerate", "datasets", "trl.trainer.grpo_trainer", "trl.trainer.sft_trainer"]
if importlib.util.find_spec("vllm"):
    mods += ["vllm", "xgrammar", "llguidance"]
    if importlib.util.find_spec("trl"):
        mods += ["trl.generation.vllm_generation", "trl.generation.vllm_client"]
for m in mods:
    mod = importlib.import_module(m)
    print(m, getattr(mod, "__version__", ""))
if importlib.util.find_spec("vllm"):
    from vllm import LLM, SamplingParams  # noqa: F401
    from vllm.sampling_params import StructuredOutputsParams  # noqa: F401
    print("vllm.LLM, SamplingParams, StructuredOutputsParams import OK")
    if importlib.util.find_spec("trl"):
        import trl.generation.vllm_client as c
        print("TRL client: stateful trainer engine (vLLM >= 0.28 path)", c._HAS_STATEFUL_TRAINER_ENGINE)
print(f"imports {time.time() - t0:.1f} s")
PY
done
step build uv pip install --python /opt/venv/bin/python --offline --no-deps --no-build-isolation -e .
step objdel /opt/venv/bin/python -c "import objdel, objdel._polycore, objdel.train.configs, objdel.train.dataset, objdel.train.rewards, objdel.train.sample, objdel.train.callbacks; print('objdel.train imports OK')"
DESELECT=(); grep -q VirtualApple /proc/cpuinfo && DESELECT=(--deselect tests/test_polycore.py::test_score_batch_out_of_memory_raises_in_python)
step fast /opt/venv/bin/python -m pytest -q -m "not slow and not data" -p no:cacheprovider --junitxml="$IO/junit-$LABEL.xml" "${DESELECT[@]}"
step hashes /opt/venv/bin/python -I -m experiments.phase1.bench.xarch_hashes --out "$IO/hashes-$LABEL.json"
uv pip freeze --python /opt/venv/bin/python > "$IO/$LABEL-freeze-venv.txt" 2>&1
[ -x /opt/vllm/bin/python ] && uv pip freeze --python /opt/vllm/bin/python > "$IO/$LABEL-freeze-vllm.txt" 2>&1
du -sh /opt/venv /opt/vllm 2>/dev/null | tee "$IO/$LABEL-du.txt"
