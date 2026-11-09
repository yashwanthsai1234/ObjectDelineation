"""Pure helpers of tests/test_docs.py (DOC4): what each docs check reads. No I/O outside the repo root it is given."""

from __future__ import annotations

import re
import urllib.parse
from pathlib import Path

import yaml

# Frozen records are exempt by path, never by banner (a banner cannot go into .py data).
FROZEN_GLOBS = [
    "docs/research/**", "experiments/phase0/**", "experiments/**/lab/**", "experiments/**/impl/**",
    "experiments/cuda_checker/results/**", "experiments/*/PLAN*.md", "experiments/*/JUDGE_*.md",
    "experiments/*/V0_V1.md", "experiments/*/RESULTS*.md", "experiments/*/EXECUTION.md", "experiments/*/RESEARCH.md",
]  # fmt: skip
# Local-only files (git-ignored: .gitignore) are not public docs; the checks read tracked docs only.
LOCAL_ONLY_GLOBS = ["ROADMAP.md", "archive/**", "docs/plan/**", "docs/implementation2.0/**"]
POLICY_MODULE = "experiments.phase2.vast.policy"


def glob_re(glob: str) -> re.Pattern:
    """`**` matches any path (slashes included); `*` and `?` match within one path part (not fnmatch)."""
    out, i = "", 0
    while i < len(glob):
        if glob.startswith("**", i):
            out, i = out + ".*", i + 2
        else:
            out += {"*": "[^/]*", "?": "[^/]"}.get(glob[i], re.escape(glob[i]))
            i += 1
    return re.compile(out + r"\Z")


def _matches(rel: str, globs: list[str]) -> bool:
    return any(glob_re(g).match(rel) for g in globs)


def is_frozen(rel: str) -> bool:
    return _matches(rel, FROZEN_GLOBS)


def is_local_only(rel: str) -> bool:
    return _matches(rel, LOCAL_ONLY_GLOBS)


def doc_files(root: Path) -> list[Path]:
    """README.md, docs/**/*.md and experiments/README.md, minus frozen and local-only paths."""
    found = [root / "README.md", root / "experiments/README.md", *sorted((root / "docs").rglob("*.md"))]
    return [p for p in found if p.is_file()
            and not is_frozen(rel := p.relative_to(root).as_posix()) and not is_local_only(rel)]  # fmt: skip


def outside_docs_dir(files: list[Path], root: Path) -> list[Path]:
    """Members not under docs/ (Zensical checks the links under docs/)."""
    return [p for p in files if p.relative_to(root).parts[0] != "docs"]


def relative_links(text: str) -> list[str]:
    """Markdown link targets that are not http(s):, mailto: or #...; handles <...> targets."""
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    out = []
    for t in re.findall(r"\]\((<[^>]*>|[^)\s]*)", text):
        t = t.strip("<>")
        if t and not t.startswith(("http:", "https:", "mailto:", "#")):
            out.append(t)
    return out


def broken_links(md: Path, root: Path) -> list[str]:
    bad = []
    for t in relative_links(md.read_text(encoding="utf-8")):
        target = urllib.parse.unquote(t.split("#")[0])
        if target and not (root / target.lstrip("/") if target.startswith("/") else md.parent / target).exists():
            bad.append(t)
    return bad


def backtick_paths(text: str) -> list[str]:
    spans = re.findall(r"`([^`\n]+)`", text)
    return [s for s in spans if s.startswith(("src/", "experiments/", "tests/", "cuda/", "docs/", "archive/"))
            and not re.search(r"[\s*<{]", s)]  # fmt: skip


def make_targets(makefile: Path) -> set[str]:
    return set(re.findall(r"^([a-zA-Z][\w-]*):", makefile.read_text(), re.M))


def make_mentions(text: str) -> set[str]:
    """`make X` in code spans and fences."""
    code = re.findall(r"`([^`\n]+)`", text) + re.findall(r"```.*?```", text, flags=re.S)
    return {m for c in code for m in re.findall(r"\bmake ([a-z][\w-]*)", c)}


def adr_files(d: Path) -> list[Path]:
    return sorted(p for p in d.glob("[0-9][0-9][0-9][0-9]-*.md"))


def adr_index(readme: Path) -> set[str]:
    return set(re.findall(r"\]\(([0-9]{4}-[^)]+\.md)\)", readme.read_text()))


def adr_headings(md: Path) -> set[str]:
    text = md.read_text()
    return {m.strip() for m in re.findall(r"^#{2,3} (.+)$", text, re.M)} | {
        m for m in ("Status", "Date", "Deciders") if re.search(rf"^\*\*{m}:\*\*", text, re.M)}


def config_files(root: Path) -> list[Path]:
    """specs/*.yaml and experiments/**/configs/*.yaml; raises when there are none."""
    found = sorted((root / "specs").glob("*.yaml")) + sorted((root / "experiments").glob("**/configs/*.yaml"))
    if not found:
        raise ValueError("no run specs or experiment configs found")
    return found


def adr_table(plan: Path) -> dict[str, tuple[str, str]]:
    """Rows `| NNNN | title | task |` after the line '**ADR numbers (one table'."""
    lines = plan.read_text().splitlines()
    start = next((i for i, ln in enumerate(lines) if ln.startswith("**ADR numbers (one table")), None)
    if start is None:
        raise ValueError("ADR table not found")
    rows: dict[str, tuple[str, str]] = {}
    for ln in lines[start + 1:]:
        m = re.match(r"\|\s*(\d{4})\s*\|\s*(.+?)\s*\|\s*(.+?)\s*\|\s*$", ln)
        if m:
            if m[1] in rows:
                raise ValueError(f"ADR number {m[1]} repeated")
            rows[m[1]] = (m[2], m[3])
        elif rows and not ln.startswith("|"):
            break
    if not rows:
        raise ValueError("ADR table has no rows")
    return rows


def adr_mentions(text: str) -> list[str]:
    nums = re.findall(r"ADR[- ]?(\d{4})", text)
    for span in re.findall(r"`([^`\n]+)`", text):
        m = re.match(r"(\d{4})-[a-z0-9]+(?:-[a-z0-9]+)+\.md$", span.rsplit("/", 1)[-1])
        if m:
            nums.append(m[1])
    return nums


def allowed_versions(root: Path, name: str) -> set[str]:
    """Every exact `name==X` pin in pyproject.toml and the lock files."""
    texts = [(root / "pyproject.toml").read_text()]
    texts += [p.read_text() for p in sorted((root / "experiments/phase2/docker/locks").glob("*.txt"))]
    return {v for t in texts for v in re.findall(rf"(?im)^\W*{re.escape(name)}==([\w.]+)", t)}


def version_mentions(text: str, name: str) -> list[str]:
    pat = rf"(?i)\b{name}(?:\[\w+\])?(?:==|\s+)v?(\d+\.\d+(?:\.\d+)?)\b(?!\s*(?:s|ms|GB|MB|%|x)\b)"
    out = re.findall(pat, text)
    if name in ("trl", "transformers"):
        out += re.findall(rf"github\.com/huggingface/{name}/blob/v(\d+\.\d+\.\d+)/", text)
    return out


def version_legal(mention: str, allowed: set[str]) -> bool:
    return any(a == mention or a.startswith(mention + ".") for a in allowed)


def cost_mentions(text: str) -> list[tuple[str, float]]:
    """Dollar values on `Alarm` / `Hard stop` table rows and in 'hard stop is $X' sentences."""
    out = []
    for ln in text.splitlines():
        cell = ln.split("|")[1].strip().lower() if ln.startswith("|") and ln.count("|") > 2 else ""
        if cell in ("alarm", "hard stop"):
            out += [(cell, float(v)) for v in re.findall(r"\$(\d+(?:\.\d+)?)", ln)]
    for kind, v in re.findall(r"(?i)\b(hard stop|alarm)\s+(?:is|=|of|at)\s+\$(\d+(?:\.\d+)?)", text):
        out.append((kind.lower(), float(v)))
    return out


def headline_strings(h: dict) -> list[str]:
    """The test numbers as README.md writes them (rounded only here), plus '<unit>-level' when not item-level."""
    m = {k: v["mean"] for k, v in h["test"]["metrics"].items()}
    out = [f"IoU {m['iou_tc']:.3f}", f"boundary F-score {m['f5']:.3f}", f"{m['legal'] * 100:.1f} %"]
    return out + ([f"{h['ci_unit']}-level"] if h["ci_unit"] != "item" else [])


def nav_pages(mkdocs_yml: Path) -> set[str]:
    cfg = yaml.load(mkdocs_yml.read_text(), Loader=yaml.BaseLoader)  # BaseLoader: the python/name tags stay text

    def walk(n):
        if isinstance(n, str):
            yield n
        elif isinstance(n, list):
            for x in n:
                yield from walk(x)
        elif isinstance(n, dict):
            for x in n.values():
                yield from walk(x)

    return {p for p in walk(cfg.get("nav", [])) if p.endswith(".md")}


def excluded_docs(mkdocs_yml: Path) -> list[str]:
    cfg = yaml.load(mkdocs_yml.read_text(), Loader=yaml.BaseLoader)
    return [ln.strip() for ln in cfg.get("exclude_docs", "").splitlines() if ln.strip()]


def tables(text: str) -> list[tuple[int, list[str]]]:
    """(start line, header cells) of each Markdown table."""
    lines, out = text.splitlines(), []
    for i, ln in enumerate(lines[:-1]):
        header = re.match(r"\|[\s:|-]+\|\s*$", lines[i + 1])
        if ln.startswith("|") and header and (i == 0 or not lines[i - 1].startswith("|")):
            out.append((i + 1, [c.strip() for c in ln.strip().strip("|").split("|")]))
    return out

