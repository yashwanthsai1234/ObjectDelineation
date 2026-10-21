"""B5 / H-A5-7: which training + serving stacks resolve together for Linux x86_64 (Vast) at Python 3.12?

Runs `uv pip compile --python-platform x86_64-manylinux_2_40 --python-version 3.12` (the Vast image is Debian trixie,
glibc 2.41) on candidate requirement sets, optionally with uv overrides (`--override`),
(no install, no network beyond the package index) and records, per candidate, whether it resolves, the resolved
versions of the packages that matter, and uv's conflict explanation when it does not. Linux x86_64 wheels of torch on
PyPI are the CUDA builds (torch 2.13.0 / 2.14.1: CUDA 13.0), so the resolution is the one a Vast image would get.

Usage: python experiments/phase2/lab/B5/resolve_matrix.py <out_dir>
Writes <out_dir>/<candidate>.in, <candidate>.txt (lock or uv's error), <candidate>.override and summary.json.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

# objdel core + dev pins (pyproject.toml) as the Mac's .venv resolves them
CORE = ["numpy==2.5.4", "scipy==1.18.1", "pillow==13.0.0", "opencv-python-headless==5.0.0.93", "pyarrow==26.0.0",
        "numba==0.68.0", "pytest>=8.0", "hypothesis>=6.100", "scikit-build-core==1.1.1", "nanobind==3.1.0"]
TRAIN_PINNED = ["torch==2.14.1", "transformers==5.20.1", "trl==1.16.0", "datasets==5.1.0", "peft==0.21.2",
                "accelerate==1.15.0"]
TRAIN_REPIN = ["torch==2.13.0", "transformers==5.19.0", "trl==1.16.0", "datasets==5.1.0", "peft==0.21.2",
               "accelerate==1.15.0"]

CORE_NO_NUMBA = [r for r in CORE if not r.startswith(("numba", "numpy"))]
NUMBA_068 = ["numba==0.68.0"]  # override of vLLM's numba==0.65.0 (vLLM uses numba only for n-gram speculation)
# ... plus numpy: vLLM (and SGLang) require mistral-common[image], which caps numpy < 2.4 on Python < 3.13
NUMBA_NUMPY = NUMBA_068 + ["numpy==2.5.4"]

CANDIDATES = {
    # (b) the two environments separately
    "b_train_as_pinned": CORE + TRAIN_PINNED,
    "b_sampler_vllm032": CORE + ["vllm==0.32.0", "trl==1.16.0", "peft==0.21.2", "datasets==5.1.0"],
    "b_sampler_vllm032_bare": ["vllm==0.32.0"],
    # (a) one environment
    "a_repin_vllm032_tf519": CORE + TRAIN_REPIN + ["vllm==0.32.0"],
    "a_trl_extra_vllm_tf518plus": CORE + ["trl[vllm]==1.16.0", "transformers>=5.18", "datasets==5.1.0",
                                          "peft==0.21.2", "accelerate==1.15.0"],
    "a_trl_extra_vllm_any_tf": CORE + ["trl[vllm]==1.16.0", "datasets==5.1.0", "peft==0.21.2",
                                       "accelerate==1.15.0"],
    "a_vllm032_torch2141": CORE + ["torch==2.14.1", "vllm==0.32.0", "trl==1.16.0"],
    "a_vllm032_tf5201": CORE + ["transformers==5.20.1", "vllm==0.32.0", "trl==1.16.0"],
    "a_repin_vllm032_numba_override": (CORE + TRAIN_REPIN + ["vllm==0.32.0"], NUMBA_068),
    "a_repin_vllm032_numba065": CORE_NO_NUMBA + TRAIN_REPIN + ["vllm==0.32.0"],
    "a_repin_vllm032_numba_numpy_override": (CORE + TRAIN_REPIN + ["vllm==0.32.0"], NUMBA_NUMPY),
    "a_trl_extra_vllm_tf518plus_numba_numpy_override": (CORE + ["trl[vllm]==1.16.0", "transformers>=5.18"],
                                                        NUMBA_NUMPY),
    "a_trl_extra_vllm_numba_numpy_override": (CORE + ["trl[vllm]==1.16.0", "datasets==5.1.0", "peft==0.21.2",
                                                      "accelerate==1.15.0"], NUMBA_NUMPY),
    "a_trl_extra_vllm_numba_override": (CORE + ["trl[vllm]==1.16.0", "datasets==5.1.0", "peft==0.21.2",
                                                "accelerate==1.15.0"], NUMBA_068),
    "a_trl_extra_vllm_tf518plus_numba_override": (CORE + ["trl[vllm]==1.16.0", "transformers>=5.18"], NUMBA_068),
    # (c) SGLang next to the training pins
    "c_sglang0522_train": CORE + ["torch==2.14.1", "transformers==5.19.0", "trl==1.16.0", "datasets==5.1.0",
                                  "peft==0.21.2", "accelerate==1.15.0", "sglang==0.5.22"],
    "c_sglang0522_train_numba_override": (CORE + ["torch==2.14.1", "transformers==5.19.0", "trl==1.16.0",
                                                  "datasets==5.1.0", "peft==0.21.2", "accelerate==1.15.0",
                                                  "sglang==0.5.22"], NUMBA_068),
    "c_sglang0522_train_numba_numpy_override": (CORE + ["torch==2.14.1", "transformers==5.19.0", "trl==1.16.0",
                                                        "datasets==5.1.0", "peft==0.21.2", "accelerate==1.15.0",
                                                        "sglang==0.5.22"], NUMBA_NUMPY),
}
WATCH = ("torch", "torchvision", "transformers", "trl", "vllm", "sglang", "xgrammar", "llguidance", "flashinfer-python",
         "numpy", "numba", "datasets", "peft", "accelerate", "tokenizers", "huggingface-hub", "triton", "outlines",
         "cuda-toolkit", "nvidia-cudnn-cu13", "opencv-python-headless", "scipy", "pillow")


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    summary = {}
    for name, spec in CANDIDATES.items():
        reqs, overrides = spec if isinstance(spec, tuple) else (spec, [])
        inp = out / f"{name}.in"
        inp.write_text("\n".join(reqs) + "\n")
        cmd = ["uv", "pip", "compile", str(inp), "--python-platform", "x86_64-manylinux_2_40", "--python-version",
               "3.12", "--no-header", "--no-annotate", "--quiet"]
        if overrides:
            (out / f"{name}.override").write_text("\n".join(overrides) + "\n")
            cmd += ["--override", str(out / f"{name}.override")]
        p = subprocess.run(cmd, capture_output=True, text=True)
        (out / f"{name}.txt").write_text(p.stdout if p.returncode == 0 else p.stderr)
        rec = {"resolves": p.returncode == 0, "overrides": overrides}
        if p.returncode == 0:
            pins = dict(m.groups() for m in re.finditer(r"^([A-Za-z0-9_.\-\[\]]+)==(\S+)", p.stdout, re.M))
            rec["n_packages"] = len(pins)
            rec["versions"] = {k: pins[k] for k in WATCH if k in pins}
        else:
            rec["error"] = " ".join(p.stderr.split())[:1500]
        summary[name] = rec
        print(name, "OK" if rec["resolves"] else "UNSAT", flush=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
