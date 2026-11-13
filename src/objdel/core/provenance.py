"""Where the source tree is and what code a run used (stdlib only, plus the two stdlib-only helpers of objdel).

`repo_root()` is $OBJDEL_REPO when set, else the first parent of this file that holds both pyproject.toml and
src/objdel: a search, not a fixed parents[n], so it holds on the Mac and on a box (the editable install keeps the
tree; the box has no .git). `code_hashes()` hashes every file of `CODE` under it: preflight's stamp and the run
manifests, so the launch gate sees a code edit.
"""

import hashlib
import os
from pathlib import Path

from objdel.core.hashing import file_sha256
from objdel.core.io import git_sha

CODE = ("src/objdel/**/*.py", "src/objdel/*.cpp", "cpp/**/*.hpp", "cuda/*.h", "cuda/*.cpp", "tools/**/*.py",
        "CMakeLists.txt", "pyproject.toml")  # everything that runs on the box


def repo_root() -> Path:
    """$OBJDEL_REPO (it must hold pyproject.toml), else the nearest parent with pyproject.toml and src/objdel."""
    if env := os.environ.get("OBJDEL_REPO"):
        if not (Path(env) / "pyproject.toml").is_file():
            raise FileNotFoundError(f"OBJDEL_REPO={env} holds no pyproject.toml")
        return Path(env)
    here = Path(__file__).resolve()
    for d in here.parents:
        if (d / "pyproject.toml").is_file() and (d / "src" / "objdel").is_dir():
            return d
    raise FileNotFoundError(f"no parent of {here} holds pyproject.toml and src/objdel: set OBJDEL_REPO")


STORE, EVAL_VIEWS, PREPARED = "data/store", "data/eval_views", "data/prepared"  # default CLI folders, under the repo
SERIES_Z = "runs/b1/series_z.json"


def reference_dir() -> Path:
    """The repo's reference/ folder (the pinned JSON records)."""
    return repo_root() / "reference"


def data_root(default: str = "data") -> Path:
    """The Mac value of {data}: $OBJDEL_DATA_ROOT (empty counts as unset), else `default`."""
    return Path(os.environ.get("OBJDEL_DATA_ROOT") or default)


def run_dir(default: str = "runs") -> Path:
    """The runs folder: $OBJDEL_RUN_DIR (empty counts as unset), else `default`."""
    return Path(os.environ.get("OBJDEL_RUN_DIR") or default)


def code_hashes() -> dict:
    """SHA-256 over (path, file SHA-256) of every CODE file, the git code identity (objdel.core.io.git_sha; None without
    .git) and the file count."""
    root = repo_root()
    files = sorted({p for g in CODE for p in root.glob(g)})
    h = hashlib.sha256(b"".join(p.relative_to(root).as_posix().encode() + bytes.fromhex(file_sha256(p)) for p in files))
    return {"objdel_and_script_sha256": h.hexdigest(), "git_head": git_sha(root, required=False),
            "files": len(files)}
