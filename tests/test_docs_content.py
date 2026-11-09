"""Content tests for the documentation tasks DOC1 and DOC2 (plan section DOC.9: I15, I18, X1, U19, I02).

Bad strings are built at run time so that this file never trips the scans it implements.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# Frozen paths (DOC.spec "Definitions"); experiments/phase1 is live after DOC4, except its lab folder.
FROZEN_PARTS = ("/lab/", "/impl/")
FROZEN_PREFIXES = ("docs/research/", "experiments/phase0/", "experiments/cuda_checker/", "docs/implementation2.0/")
BANNER = "> Frozen record (written "


def _read(rel: str) -> str:
    return (REPO / rel).read_text(encoding="utf-8")


def _git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, check=True).stdout


def _tracked_or_new(*globs: str) -> list[str]:
    out = _git("ls-files", "-z", "--cached", "--others", "--exclude-standard", "--", *globs)
    return [p for p in out.split("\0") if p and (REPO / p).is_file()]


def _is_frozen(rel: str) -> bool:
    return rel.startswith(FROZEN_PREFIXES) or any(part in "/" + rel for part in FROZEN_PARTS)


# ---------------------------------------------------------------- DOC1


def test_plan_moved_out_of_frozen_folder():
    names = sorted(p.name for p in (REPO / "docs/plan").iterdir())
    assert names == ["Architecture.md", "Direction.md", "Research2.0.md", "previous-vs-new.md"]
    assert not (REPO / "docs/research" / ("research" + "2.0")).exists()


def test_no_reference_to_the_old_plan_path():
    old = ("research" + "2.0", "previous " + "vs new", "previous%20" + "vs%20new")
    targets = ["README.md", "ROADMAP.md", "Makefile", "pyproject.toml", "docs", "src", "tests", "experiments"]
    targets.append(".github")
    hits = []
    for rel in _tracked_or_new(*targets):
        if rel.startswith(("docs/research/", "docs/implementation2.0/")):
            continue  # frozen, and the implementation plan names the old path as history
        try:
            text = _read(rel)
        except UnicodeDecodeError:
            continue
        hits += [f"{rel}: {o}" for o in old if o in text]
    assert hits == []


def test_research_page_corrected_ceiling():
    row = next(ln for ln in _read("docs/plan/Research2.0.md").splitlines() if ln.startswith("| All images"))
    assert "0.723" in row and "0.734" not in row


def test_gitignore_has_docs_build_rules():
    rules = set(_read(".gitignore").splitlines())
    wanted = ("/site/", "/.cache/", ".lycheecache", "docs/research/*.pdf", "docs/implementation2.0/kernels/*/build/")
    for rule in wanted:
        assert rule in rules, rule
    ignored = subprocess.run(
        ["git", "-C", str(REPO), "check-ignore", "-q", "docs/research/2212.07020v1.pdf"], capture_output=True
    )
    assert ignored.returncode == 0


def test_no_scratch_or_home_paths_in_implementation_docs():
    patterns = [
        "/" + "private/tmp",
        "/" + "Users/",
        "S" + "P/impl",
        "S" + "P/[A-Z]",
    ]
    rx = re.compile("|".join(patterns))
    hits = []
    for rel in _tracked_or_new("docs/implementation2.0"):
        if rel.startswith("docs/implementation2.0/notes/research/"):
            continue  # written by the research wave, scanned by its own owner
        try:
            lines = _read(rel).splitlines()
        except UnicodeDecodeError:
            continue
        hits += [f"{rel}:{i}" for i, ln in enumerate(lines, 1) if rx.search(ln)]
    assert hits == []


def test_kernel_build_output_is_not_tracked_or_cited():
    assert [p for p in _git("ls-files", "docs/implementation2.0/kernels").splitlines() if "/build/" in p] == []
    for k in ("K2", "K3"):
        assert not re.search(r"\bbuild/", _read(f"docs/implementation2.0/kernels/{k}/DESIGN.md"))
    results = REPO / "docs/implementation2.0/kernels"
    assert not list(results.rglob("*.log"))  # *.log is git-ignored: evidence is kept as .log.txt
