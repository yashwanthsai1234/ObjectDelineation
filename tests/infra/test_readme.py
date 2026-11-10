"""README.md stays true: its Python quick start runs and prints what its comment says, its make targets exist (B4-6)."""
import contextlib
import io
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
README = (ROOT / "README.md").read_text()


def test_readme_python_quick_start_runs():
    blocks = re.findall(r"```python\n(.*?)```", README, re.S)
    assert blocks
    for code in blocks:
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            exec(compile(code, "README.md", "exec"), {})
        expected = re.findall(r"print\(.*#\s*(.+)$", code, re.M)
        assert out.getvalue().split() == " ".join(expected).split()


def test_readme_make_targets_exist():
    targets = set(re.findall(r"^([a-z][a-z-]*):", (ROOT / "Makefile").read_text(), re.M))
    used = set(re.findall(r"^make ([a-z-]+)", README, re.M)) | set(re.findall(r"make ([a-z-]+) adds", README))
    assert used and used <= targets
