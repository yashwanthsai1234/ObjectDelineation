"""Limit causal-conv1d 1.7.0's CUDA targets to our cards: sm_80 (its code runs on the 3090 sm_86 and the 4090 sm_89),
sm_120 (5090, RTX PRO 6000) and compute_120 PTX. setup.py hard-codes nine targets, which multiplies the compile time
under the Mac's amd64 emulation. Also nvcc --threads 4 -> 1: compiling three targets in parallel OOM-killed
causal_conv1d_bwd.cu in the 4 GB builder (V0 rebuild, 2026-11-01).

    python patch_arch.py setup.py
"""

import sys

START = '        cc_flag.append("-gencode")\n        cc_flag.append("arch=compute_75'
END = "    # HACK: The compiler flag -D_GLIBCXX_USE_CXX11_ABI"
THREADS = '["--threads", "4"]'
OURS = (
    '        cc_flag += ["-gencode", "arch=compute_80,code=sm_80", "-gencode", "arch=compute_120,code=sm_120",\n'
    '                    "-gencode", "arch=compute_120,code=compute_120"]\n\n'
)


def patch(text: str) -> str:
    if START not in text or END not in text or THREADS not in text:
        raise ValueError("causal-conv1d setup.py: gencode block or nvcc --threads not found (version changed?)")
    i, j = text.index(START), text.index(END)
    return (text[:i] + OURS + text[j:]).replace(THREADS, '["--threads", "1"]')


if __name__ == "__main__":
    p = sys.argv[1]
    with open(p) as f:
        out = patch(f.read())
    with open(p, "w") as f:
        f.write(out)
    print("patch_arch: sm_80, sm_120, compute_120; nvcc --threads 1")
