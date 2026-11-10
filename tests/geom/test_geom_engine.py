"""geom.engine is the one importer of the CUDA extension (Implementation.md E3a): its __all__ covers every
`_gpucheck.<name>` that src, tests and experiments use, and each name exists in the built extension."""

import re
from pathlib import Path

from objdel.geom import engine
from objdel.geom import engine as ext

ROOT = Path(__file__).resolve().parents[2]
USE = re.compile(r"\b_gpucheck\.([A-Za-z]\w*)\b")


def test_engine_covers_every_extension_name_the_repo_uses():
    files = [p for d in ("src", "tests", "experiments") for p in (ROOT / d).rglob("*.py")]
    used = {n for p in files for n in USE.findall(p.read_text())} - {"py"}  # "_gpucheck.py": a file name in prose
    assert used - set(engine.__all__) == set()
    assert all(hasattr(ext, n) and getattr(engine, n) is getattr(ext, n) for n in engine.__all__)
