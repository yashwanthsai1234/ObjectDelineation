"""Throwaway git repos in tmp_path for the tests of code that reads `git ls-files`."""
import subprocess


def init_repo(root, files, *, content="x", add=False):
    """Write each of `files` (a path relative to `root`) holding `content`, then `git init -q` in `root` (and
    `git add -A` when `add`). Only ever call it with a tmp_path: it runs git there, never in the real repo."""
    for f in files:
        (root / f).parent.mkdir(parents=True, exist_ok=True)
        (root / f).write_text(content)
    subprocess.run(["git", "init", "-q"], cwd=root, check=True)
    if add:
        subprocess.run(["git", "add", "-A"], cwd=root, check=True)
    return root
