"""DOC4: the docs stay true. One test per check; each has a negative test on tmp_path that must report its fault.

Fixture strings that would trip a scan are built at run time, so this file never matches itself.
"""

from __future__ import annotations

import importlib
import json
import re
import subprocess
from pathlib import Path

import pytest
import yaml
from docs_check import (
    POLICY_MODULE,
    adr_files,
    adr_headings,
    adr_index,
    adr_mentions,
    adr_table,
    allowed_versions,
    backtick_paths,
    broken_links,
    config_files,
    cost_mentions,
    doc_files,
    excluded_docs,
    glob_re,
    headline_strings,
    is_frozen,
    make_mentions,
    make_targets,
    nav_pages,
    outside_docs_dir,
    project_version,
    tables,
    version_legal,
    version_mentions,
)

REPO = Path(__file__).resolve().parents[2]
PLAN = REPO / "docs/implementation2.0/Implementation.md"
NEEDS_PLAN = pytest.mark.skipif(not PLAN.exists(), reason="the plan is local-only (git-ignored)")
ADR_HEADINGS = {"Context", "Decision", "Options considered", "Trade-off analysis", "Consequences", "Action items",
                "Status", "Date", "Deciders"}  # fmt: skip


def put(root: Path, rel: str, text: str = "") -> Path:
    p = root / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


# ---------------------------------------------------------------- the glob matcher
def test_frozen_globs_match_by_path_part():
    assert is_frozen("docs/research/reports/x.md") and is_frozen("experiments/phase2/impl/B5/x.py")
    assert is_frozen("experiments/phase2/PLAN.md") and not is_frozen("experiments/phase2/sub/PLAN.md")
    assert not is_frozen("experiments/phase1/prepare_rl.py") and not is_frozen("experiments/cuda_checker/run_box.sh")
    assert glob_re("a/*.md").match("a/b.md") and not glob_re("a/*.md").match("a/b/c.md")


# ---------------------------------------------------------------- 1 links
def test_relative_links_resolve():
    files = doc_files(REPO)
    assert files
    bad = {str(p.relative_to(REPO)): b for p in outside_docs_dir(files, REPO) if (b := broken_links(p, REPO))}
    assert not bad


def test_broken_link_is_reported(tmp_path):
    put(tmp_path, "docs/a.md", "")
    md = put(tmp_path, "README.md", "[ok](docs/a.md) [no](docs/b.md) [web](https://x.y) [t](<docs/a%2Eb.md#h>)")
    assert broken_links(md, tmp_path) == ["docs/b.md", "docs/a%2Eb.md#h"]


# ---------------------------------------------------------------- 2 backtick paths
def local_only(path: str) -> bool:
    """A path in a git-ignored top folder that this checkout lacks (the owner keeps archive/ on his disk only,
    a255988): a worktree cannot check it, the owner's checkout does (the folder is there, so the path must be)."""
    top = path.split("/", 1)[0]
    if (REPO / top).exists():
        return False
    return subprocess.run(["git", "check-ignore", "-q", path], cwd=REPO).returncode == 0


def test_backtick_paths_exist():
    paths = {(p, s.split(":")[0]) for p in doc_files(REPO) for s in backtick_paths(p.read_text())}
    missing = {f"{p.relative_to(REPO)}: {s}" for p, s in paths
               if not (REPO / s.rstrip("/")).exists() and not local_only(s)}  # "archive/": a folder, for git
    assert not missing


def test_only_an_ignored_absent_folder_is_local_only():
    assert not local_only("src/objdel/no_such.py")  # tracked folder: a missing file is a fault
    assert not local_only("no_such_folder/x.md")  # not ignored


def test_backtick_paths_rule():
    assert backtick_paths("`src/x/y.py` `src/a b` `docs/*.md` `tests/<x>` `foo/bar`") == ["src/x/y.py"]


# ---------------------------------------------------------------- 3 make targets
def test_make_mentions_are_targets():
    targets = make_targets(REPO / "Makefile")
    bad = {f"{p.relative_to(REPO)}: make {m}" for p in doc_files(REPO) for m in make_mentions(p.read_text())
           if m not in targets}  # fmt: skip
    assert not bad


def test_unknown_make_target_is_reported(tmp_path):
    mk = put(tmp_path, "Makefile", "test:  ## x\n\techo\n")
    assert make_targets(mk) == {"test"}
    assert make_mentions("run `make test` and\n```\nmake nope\n```") == {"test", "nope"}


# ---------------------------------------------------------------- 4 and 5 ADRs
ADR_DIR = REPO / "docs/architecture/decisions"


def test_adr_index_complete():
    files = adr_files(ADR_DIR)
    numbers = [f.name[:4] for f in files]
    assert len(numbers) == len(set(numbers))
    assert {f.name for f in files} == adr_index(ADR_DIR / "README.md")


def test_adr_index_reports_missing_and_duplicate(tmp_path):
    put(tmp_path, "0012-a.md")
    put(tmp_path, "0012-b.md")
    put(tmp_path, "_template.md")
    idx = put(tmp_path, "README.md", "| [0012](0012-a.md) |")
    files = adr_files(tmp_path)
    assert [f.name for f in files] == ["0012-a.md", "0012-b.md"]
    assert {f.name for f in files} != adr_index(idx)
    assert len({f.name[:4] for f in files}) < len(files)


@pytest.mark.parametrize("adr", adr_files(ADR_DIR), ids=lambda p: p.name[:4])
def test_adr_headings(adr):
    """An Accepted record has every heading; a Proposed one (still being written) has the three header lines."""
    need = ADR_HEADINGS if "**Status:** Accepted" in adr.read_text() else {"Status", "Date", "Deciders"}
    assert need <= adr_headings(adr)


# ---------------------------------------------------------------- 6 configs
def test_configs_documented():
    text = (REPO / "docs/reference/configs.md").read_text()
    missing = [str(c.relative_to(REPO)) for c in config_files(REPO) if c.name not in text]
    assert not missing


def test_configs_raise_on_nothing(tmp_path):
    with pytest.raises(ValueError):
        config_files(tmp_path)


# ---------------------------------------------------------------- 7 version strings
PACKAGES = ("trl", "transformers", "torch", "vllm")


def _version_sources():
    own = {"tests/infra/test_docs.py", "tests/docs_check.py"}
    files = [p for g in ("src/**/*.py", "tests/**/*.py", "experiments/**/*.py", "experiments/**/Dockerfile*",
                         "experiments/**/configs/*.yaml", "specs/*.yaml") for p in REPO.glob(g)]  # fmt: skip
    files = [p for p in files if p.is_file() and not is_frozen(rel := p.relative_to(REPO).as_posix())
             and rel not in own]
    return files + doc_files(REPO)


def test_version_strings_match_pins():
    allowed = {n: allowed_versions(REPO, n) for n in PACKAGES}
    assert all(allowed.values())
    bad = sorted({f"{p.relative_to(REPO)}: {n} {v}" for p in _version_sources()
                  for n in PACKAGES for v in version_mentions(p.read_text(errors="ignore"), n)
                  if not version_legal(v, allowed[n])})  # fmt: skip
    assert not bad


def test_version_rule():
    old, torch, vllm_s = "TR" + "L 1.16", "torch 2.13.0", "vL" + "LM 0.09 s"
    assert version_mentions(old, "trl") == ["1.16"] and version_mentions(vllm_s, "vllm") == []
    assert version_mentions("trl==" + "1.16.0", "trl") == ["1.16.0"]
    assert version_legal("2.13.0", {"2.13.0"}) and version_legal("1.17", {"1.17.1"})
    assert not version_legal("1.16", {"1.17.1"})
    assert version_mentions(torch, "torch") == ["2.13.0"]


def test_allowed_versions_read_pins(tmp_path):
    put(tmp_path, "pyproject.toml", '    "trl==1.17.1",\n')
    put(tmp_path, "experiments/phase2/docker/locks/a.txt", "trl==1.18.0\n")
    assert allowed_versions(tmp_path, "trl") == {"1.17.1", "1.18.0"}


# ---------------------------------------------------------------- 8 cost constants
def test_cost_constants_in_docs():
    policy = importlib.import_module(POLICY_MODULE)
    const = {"alarm": policy.ALARM, "hard stop": policy.HARD_STOP}
    seen = {k: 0 for k in const}
    for p in doc_files(REPO):
        for kind, v in cost_mentions(p.read_text()):
            assert v == const[kind], f"{p.relative_to(REPO)}: {kind} ${v}"
            seen[kind] += p.name == "gpu-run.md"
    assert all(seen.values())


def test_cost_mentions_rule():
    assert cost_mentions("the hard stop is $15") == [("hard stop", 15.0)]
    assert cost_mentions("| Hard stop | refused past $15 | x |") == [("hard stop", 15.0)]
    assert cost_mentions("rose from $15 to $16.60") == []


# ---------------------------------------------------------------- 9 headline
def test_readme_headline_matches():
    h = json.loads((REPO / "docs/results/headline.json").read_text())
    readme = (REPO / "README.md").read_text()
    assert [s for s in headline_strings(h) if s not in readme] == []


# ---------------------------------------------------------------- 10 citation, 12 attribution
def test_citation_version():
    cff = yaml.safe_load((REPO / "CITATION.cff").read_text())
    assert str(cff["version"]) == project_version(REPO)


def test_no_ai_attribution():
    bad = "Co-Authored-" + "By: Claude"
    gen = "Generated with " + "Claude"
    for rel in (".github/pull_request_template.md", "docs/contributing.md", "CITATION.cff"):
        text = (REPO / rel).read_text()
        assert bad not in text and gen not in text, rel


# ---------------------------------------------------------------- 11 archive (git-ignored: checked only when tracked)
ARCHIVE_PARTS = ("What it was", "What replaced it", "Retiring commit", "How to restore")


def archive_problems(root: Path) -> list[str]:
    out = []
    for d in sorted(p for p in (root / "archive").glob("*") if p.is_dir()):
        readme = d / "README.md"
        heads = set(re.findall(r"^#{1,3} (.+)$", readme.read_text(), re.M)) if readme.exists() else set()
        out += [f"{d.name}: no heading {h!r}" for h in ARCHIVE_PARTS if h not in heads]
    return out


def test_archive_readmes():
    tracked = subprocess.run(["git", "-C", str(REPO), "ls-files", "archive"], capture_output=True, text=True).stdout
    if not tracked.strip():
        pytest.skip("archive/ is local-only (git-ignored)")
    assert archive_problems(REPO) == []


def test_archive_folder_without_readme_is_reported(tmp_path):
    (tmp_path / "archive/x").mkdir(parents=True)
    assert len(archive_problems(tmp_path)) == 4


# ---------------------------------------------------------------- 13 nav
def test_nav_complete():
    skip = [glob_re(g) for g in excluded_docs(REPO / "mkdocs.yml")]
    nav = nav_pages(REPO / "mkdocs.yml")
    pages = {p.relative_to(REPO / "docs").as_posix() for p in (REPO / "docs").rglob("*.md")}
    missing = sorted(p for p in pages if p not in nav and not any(g.match(p) for g in skip)
                     and not p.startswith("reference/objdel/"))
    assert not missing
    assert not sorted(n for n in nav if n not in pages)


def test_nav_pages_reads_nested_nav(tmp_path):
    y = put(tmp_path, "mkdocs.yml", "nav:\n  - A:\n      - a.md\n      - B: b/b.md\n  - c.md\n"
                                    "exclude_docs: |\n  x/**\n")
    assert nav_pages(y) == {"a.md", "b/b.md", "c.md"} and excluded_docs(y) == ["x/**"]


# ---------------------------------------------------------------- 14 scratch and personal paths
def test_no_scratch_or_local_paths():
    pats = [re.compile(p) for p in ("/[p]rivate/tmp", "/[U]sers/", "S[P]/impl", "S[P]/[A-Z]")]
    files = doc_files(REPO) + [p for p in (REPO / "docs/implementation2.0").rglob("*.md") if p.is_file()
                               and not is_frozen(p.relative_to(REPO).as_posix())]  # fmt: skip
    bad = [f"{p.relative_to(REPO)}:{i}" for p in files for i, ln in enumerate(p.read_text().splitlines(), 1)
           if any(r.search(ln) for r in pats)]  # fmt: skip
    assert not bad


# ---------------------------------------------------------------- 15 and 16 plan checks (local-only plan)
@NEEDS_PLAN
def test_adr_numbers_match_front_table():
    table = adr_table(PLAN)
    on_disk = {f.name[:4]: f for f in adr_files(ADR_DIR)}
    text = PLAN.read_text()
    assert [n for n in adr_mentions(text) if n not in on_disk and n not in table] == []
    names: dict[str, set[str]] = {}
    for span in re.findall(r"`([^`\n]+)`", text):
        m = re.match(r"(\d{4})-[a-z0-9]+(?:-[a-z0-9]+)+\.md$", span.rsplit("/", 1)[-1])
        if m:
            names.setdefault(m[1], set()).add(span.rsplit("/", 1)[-1])
    assert {n: s for n, s in names.items() if len(s) > 1} == {}
    index = (ADR_DIR / "README.md").read_text()
    for n, (title, _task) in table.items():
        if n in on_disk and n in names:
            assert names[n] == {on_disk[n].name}
        if n in on_disk:
            row = next((ln for ln in index.splitlines() if f"[{n}]" in ln), "")
            assert title in row, (n, title)


def test_adr_table_rejects_repeat_and_empty(tmp_path):
    head = "**ADR numbers (one table, x)**\n\n"
    with pytest.raises(ValueError):
        adr_table(put(tmp_path, "a.md", head + "| 0012 | A | t |\n| 0012 | B | t |\n"))
    with pytest.raises(ValueError):
        adr_table(put(tmp_path, "b.md", head + "no rows\n"))
    assert adr_table(put(tmp_path, "c.md", head + "| 0012 | A | t |\n| 0013 | B | u |\n")) == {
        "0012": ("A", "t"), "0013": ("B", "u")}  # fmt: skip


def rental_tables_outside_front(text: str) -> list[int]:
    lines = text.splitlines()
    start = next((i for i, ln in enumerate(lines) if ln.startswith("## 6. GPU rentals")), None)
    if start is None or not any(ln.startswith("### 6.2 V1") for ln in lines):
        raise ValueError("front section 6 or its 6.2 V1 heading not found")
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("## ")), len(lines))
    return [n for n, head in tables(text) if any("$" in c for c in head) and not start < n - 1 < end]


@NEEDS_PLAN
def test_rental_tables_only_in_front():
    assert rental_tables_outside_front(PLAN.read_text()) == []


def test_rental_tables_rule():
    front = "## 6. GPU rentals\n\n### 6.2 V1\n\n| item | $ |\n|---|---|\n| a | 1 |\n\n## 7. Next\n"
    assert rental_tables_outside_front(front) == []
    assert rental_tables_outside_front(front + "\n### B.8.1\n\n| $ at 27.5 s |\n|---|\n| 1 |\n") == [13]
    assert rental_tables_outside_front(front + "\n| Cost |\n|---|\n| 1 |\n") == []
    with pytest.raises(ValueError):
        rental_tables_outside_front("## 6. GPU rentals\n")


# ---------------------------------------------------------------- 14 negative
def test_local_path_scan_pattern_is_built_at_run_time():
    assert not re.search("/[U]sers/", Path(__file__).read_text().replace("/[U]sers/", ""))


# ---------------------------------------------------------------- the site build (needs the docs extra)
@pytest.mark.docs_build
def test_site_builds_strict(tmp_path):
    r = subprocess.run([str(REPO / ".venv/bin/zensical"), "build", "--strict", "-f", str(REPO / "mkdocs.yml")],
                       capture_output=True, text=True, cwd=REPO)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
