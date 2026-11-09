"""JSONL reading and atomic file writing, shared (no heavy imports)."""

import json
import os
import subprocess
from pathlib import Path


def read_jsonl(path) -> list[dict]:
    """One dict per line, read line by line."""
    with open(path) as f:
        return [json.loads(line) for line in f]


def write_atomic(path, text: str) -> None:
    """Write through a temp file and os.replace: a crash leaves the old file or none, never half a file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def git(*args: str, root, required: bool = True) -> str | None:
    """git's stdout in root. A missing .git (the box gets the code by rsync) or a git failure raises when required,
    else gives None."""
    try:
        if not (Path(root) / ".git").exists():
            raise FileNotFoundError(f"{root}: no .git")
        return subprocess.run(["git", *args], cwd=root, capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        if required:
            raise
        return None
