"""Local secrets for every phase: <repo>/.secrets first, then the macOS keychain (owner 2026-11-05).

.secrets holds one `name=value` per line (# comments); the names are the keychain service names
(objdel-hf-write, objdel-vast-scoped, ...). It is git-ignored, never sent to a rented box (sync_up ships tracked
files only) and refused unless only the owner can read it (chmod 600). OBJDEL_SECRETS overrides its path (tests).
Values are never printed.
"""

import functools
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def path() -> Path:
    """The secrets file path: `$OBJDEL_SECRETS` if set, else `.secrets` at the repo root."""
    return Path(os.environ.get("OBJDEL_SECRETS", ROOT / ".secrets"))


@functools.cache
def get(name: str) -> str | None:
    """The value of `name` from .secrets, else from the keychain; None when neither has it. PermissionError when
    .secrets is readable by others."""
    p = path()
    if p.exists():
        if p.stat().st_mode & 0o077:
            raise PermissionError(f"{p}: readable by others; run chmod 600 {p}")
        lines = (x.strip() for x in p.read_text().splitlines())
        pairs = (x.split("=", 1) for x in lines if x and not x.startswith("#") and "=" in x)
        if value := {k.strip(): v.strip() for k, v in pairs}.get(name):
            return value
    try:
        r = subprocess.run(["security", "find-generic-password", "-s", name, "-w"], capture_output=True, text=True,
                           stdin=subprocess.DEVNULL)  # fmt: skip
    except FileNotFoundError:  # not macOS
        return None
    return r.stdout.strip() or None if r.returncode == 0 else None
