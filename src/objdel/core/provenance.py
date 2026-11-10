"""Where the source tree is: `repo_root()`, the one root every layer uses (no `__file__`-based roots elsewhere)."""

import os
from pathlib import Path


def repo_root() -> Path:
    """OBJDEL_REPO if set; else the first parent of this file that holds pyproject.toml and src/objdel (the Mac
    checkout, or the box's editable install, which has no .git)."""
    if env := os.environ.get("OBJDEL_REPO"):
        return Path(env)
    for p in Path(__file__).resolve().parents:
        if (p / "pyproject.toml").is_file() and (p / "src" / "objdel").is_dir():
            return p
    raise RuntimeError(f"no repo root above {__file__}: set OBJDEL_REPO to the source tree")
