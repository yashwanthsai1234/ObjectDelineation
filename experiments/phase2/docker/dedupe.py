"""Hardlink identical files across the image's two venvs (V0: /opt/vllm and /opt/venv pin the same cuBLAS, cuFFT,
cuSOLVER, cuSPARSE, ... wheels, and the first push stored them twice). A layer stores a hardlinked file once, so this
runs in the RUN that creates both venvs; a hardlink across layers is a copy. Files differing in mode stay apart.

    python dedupe.py /opt/vllm /opt/venv
"""

import hashlib
import os
import sys
from collections import defaultdict
from pathlib import Path


def digest(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def dedupe(roots, min_size: int = 1 << 20) -> int:
    """Hardlink every regular file of >= min_size bytes to the first identical one (size, mode, SHA-256); returns the
    bytes saved. Only files that share a size are hashed."""
    by_size = defaultdict(list)
    for root in roots:
        for p in sorted(Path(root).rglob("*")):
            if not p.is_symlink() and p.is_file() and (st := p.stat()).st_size >= min_size:
                by_size[(st.st_size, st.st_mode)].append(p)
    saved = 0
    for (size, _), paths in by_size.items():
        first = {}
        for p in paths if len(paths) > 1 else ():
            keep = first.setdefault(digest(p), p)
            if keep != p and not os.path.samefile(keep, p):
                tmp = p.with_name(p.name + ".dedupe-tmp")
                os.link(keep, tmp)
                os.replace(tmp, p)
                saved += size
    return saved


if __name__ == "__main__":
    print(f"dedupe: {dedupe(sys.argv[1:]) / 1e9:.2f} GB hardlinked")
