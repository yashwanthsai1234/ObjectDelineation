"""JSON and JSONL reading, atomic file and folder writing and the git code identity, shared (no heavy imports)."""

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path


def read_jsonl(path) -> list[dict]:
    """One dict per line, read line by line."""
    with open(path) as f:
        return [json.loads(line) for line in f]


def write_jsonl(path, rows) -> str:
    """rows as sorted-key JSON lines to path; the SHA-256 hex digest of the text written."""
    text = "".join(json.dumps(r, sort_keys=True) + "\n" for r in rows)
    Path(path).write_text(text)
    return hashlib.sha256(text.encode()).hexdigest()


def write_atomic(path, text: str) -> None:
    """Write through a temp file and os.replace: a crash leaves the old file or none, never half a file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def replace_dir(write, dest, tmp=None) -> None:
    """write(tmp) then rename tmp over dest: a crash never leaves a half-written kept folder (resume trusts what
    exists). tmp defaults to .<dest name>.tmp beside dest; a stale tmp is removed first."""
    dest = Path(dest)
    tmp = dest.with_name(f".{dest.name}.tmp") if tmp is None else Path(tmp)
    shutil.rmtree(tmp, ignore_errors=True)
    write(tmp)
    shutil.rmtree(dest, ignore_errors=True)
    os.replace(tmp, dest)


def read_json(path):
    """A whole JSON file."""
    return json.loads(Path(path).read_text())


def plain(x):
    """A numpy scalar or array as a plain JSON value (json.dumps default=plain)."""
    return x.tolist() if hasattr(x, "tolist") else x.item()


def write_json(path, doc) -> None:
    """doc as JSON with indent 1, sorted keys and a trailing newline (the manifests' bytes), through write_atomic."""
    write_atomic(path, json.dumps(doc, indent=1, sort_keys=True) + "\n")


def _git(args, root) -> bytes:
    if not (Path(root) / ".git").exists():
        raise FileNotFoundError(f"{root}: no .git")
    return subprocess.run(["git", *args], cwd=root, capture_output=True, check=True).stdout


def _or_none(get, required: bool):
    try:
        return get()
    except (OSError, subprocess.CalledProcessError):
        if required:
            raise
        return None


def git(*args: str, root, required: bool = True) -> str | None:
    """git's stdout in root, stripped. A missing .git (the box gets the code by rsync) or a git failure raises when
    required, else gives None."""
    return _or_none(lambda: _git(args, root).decode().strip(), required)


def git_files(root, *flags: str) -> list[str]:
    """`git ls-files -z FLAGS` in root as paths (NUL split, no strip: a name with spaces stays whole). A missing .git or
    a git failure raises."""
    return [f.decode() for f in _git(("ls-files", "-z", *flags), root).split(b"\0")[:-1]]


def git_sha(root, required: bool = True) -> str | None:
    """The code identity of the tree at root (vast sync_up writes it to the box as CODE_SHA): git rev-parse HEAD, plus
    -dirty-<first 12 hex of sha256(D)> when git status --porcelain is not empty; D = git diff HEAD, then the sorted
    untracked, not-ignored paths, one per line (a new module marks the tree dirty). No git: as git()."""
    def sha() -> str:
        head = _git(("rev-parse", "HEAD"), root).decode().strip()
        if not _git(("status", "--porcelain"), root):
            return head
        new = sorted(_git(("ls-files", "-z", "--others", "--exclude-standard"), root).split(b"\0")[:-1])
        d = _git(("diff", "HEAD"), root) + b"".join(p + b"\n" for p in new)
        return f"{head}-dirty-{hashlib.sha256(d).hexdigest()[:12]}"

    return _or_none(sha, required)


def code_sha(root, required: bool = True) -> str | None:
    """The code identity a run records: <root>/CODE_SHA stripped (sync_up writes git_sha there on the box, which has no
    .git; empty = absent), else git_sha(root, required)."""
    f = Path(root) / "CODE_SHA"
    return (f.read_text().strip() if f.is_file() else "") or git_sha(root, required)
