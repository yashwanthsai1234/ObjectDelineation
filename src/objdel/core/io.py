"""JSON and JSONL reading, atomic file writing and the git code identity, shared (no heavy imports)."""

import hashlib
import json
import os
import subprocess
from pathlib import Path

from objdel.core.hashing import bytes_sha256


def read_jsonl(path, *, skip_blank: bool = False) -> list[dict]:
    """One dict per line, read line by line; a blank line raises unless skip_blank."""
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip() or not skip_blank]


def read_jsonl_cut(path, need=()) -> list[tuple[str, dict]]:
    """(line, record) of every line that is a JSON object holding all the `need` keys, in file order; a missing file
    gives []. A line cut by a kill (bad JSON, not an object, a key missing) is skipped, anywhere in the file. A
    resumable file reads this way; a file that must be whole reads with read_jsonl (any bad line raises)."""
    out = []
    for x in Path(path).read_text().splitlines() if Path(path).exists() else []:
        try:
            c = json.loads(x)
        except ValueError:
            continue
        if isinstance(c, dict) and all(k in c for k in need):
            out.append((x, c))
    return out


def files_under(root, skip=()) -> list[tuple[str, Path]]:
    """(relative posix path, path) of every file under root except the `skip` paths, sorted by the posix string (not
    by Path: `a-b/x` sorts before `a/x` by string, after it by Path; rows_tar and files_digest hashes rest on this)."""
    return sorted((rel, p) for p in Path(root).rglob("*") if p.is_file()
                  if (rel := p.relative_to(root).as_posix()) not in skip)


def read_split(root, dataset: str, split: str) -> list[dict]:
    """The rows of `<root>/<dataset>/<split>.jsonl`."""
    return read_jsonl(Path(root, dataset, f"{split}.jsonl"))


def read_jsonl_glob(root, pattern: str) -> list[dict]:
    """The rows of every file under root that matches the glob pattern, file by file in name order."""
    return [r for p in sorted(Path(root).glob(pattern)) for r in read_jsonl(p)]


def write_jsonl(path, rows, *, sort_keys: bool = True, atomic: bool = False, **dumps) -> str:
    """rows as JSON lines (sorted keys unless sort_keys=False; `dumps`: json.dumps options) to path, line by line, its
    folder made first; through a temp file and os.replace when atomic. The SHA-256 hex digest of the text written."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp") if atomic else path
    h = hashlib.sha256()
    with open(tmp, "w", encoding="utf-8") as f:
        for r in rows:
            line = json.dumps(r, sort_keys=sort_keys, **dumps) + "\n"
            f.write(line)
            h.update(line.encode())
    if atomic:
        os.replace(tmp, path)
    return h.hexdigest()


def write_atomic(path, text: str) -> None:
    """Write through a temp file and os.replace: a crash leaves the old file or none, never half a file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


_RAISE = object()


def read_json(path, default=_RAISE):
    """A whole JSON file; `default` instead when the path is not a file (without a default: the OSError)."""
    if default is not _RAISE and not Path(path).is_file():
        return default
    return json.loads(Path(path).read_text())


def plain(x):
    """A numpy scalar or array as a plain JSON value (json.dumps default=plain)."""
    return x.tolist() if hasattr(x, "tolist") else x.item()


def json_text(doc, *, sort_keys: bool = True, strict: bool = False, default=None, newline: bool = True) -> str:
    """doc as JSON with indent 1 (the manifests' text): sorted keys, a trailing newline; strict raises on NaN, inf."""
    return json.dumps(doc, indent=1, sort_keys=sort_keys, allow_nan=not strict, default=default) + "\n" * newline


def write_json(path, doc, **options) -> None:
    """doc as json_text(doc, **options), through write_atomic (the defaults are the manifests' bytes)."""
    write_atomic(path, json_text(doc, **options))


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
        return f"{head}-dirty-{bytes_sha256(d)[:12]}"

    return _or_none(sha, required)


def code_sha(root, required: bool = True) -> str | None:
    """The code identity a run records: <root>/CODE_SHA stripped (sync_up writes git_sha there on the box, which has no
    .git; empty = absent), else git_sha(root, required)."""
    f = Path(root) / "CODE_SHA"
    return (f.read_text().strip() if f.is_file() else "") or git_sha(root, required)
