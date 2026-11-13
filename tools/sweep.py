"""Final sweep of the tree: CLEAN16 (pool names, dead code, before/after record), CLEAN-G5/G7/G10/G11 (retired trees,
models, Controller), DOC13 (docs links) and DOC5 (docstring gate), in one report.

    .venv/bin/python -I -m tools.sweep [--only a,b] [--skip a] [--json] [--limit 15]
                                                                 [--before runs/clean/before.json] [--write-after] [--wave 3]

`--wave 1|2` makes the Wave 3 checks (LFM names, pool leftovers) advisory, so the sweep can run on a tree before D7/E.

Each check is hard (a finding fails the sweep, exit 1) or advisory (listed, never fails). Reads git, source and docs
only: no `data/`, no network. "Public" files are `git ls-files -co --exclude-standard`: tracked or new and not ignored.
Internal plans, reviews, notes, ROADMAP and archive/ are git-ignored and local-only, so public docs must not link them.
"""

import argparse
import ast
import json
import re
import subprocess
import sys
import tomllib
import urllib.parse
from collections.abc import Callable
from pathlib import Path

from objdel.core.io import git_files, read_json
from objdel.core.provenance import repo_root

ROOT = repo_root()
TEXT = (".py", ".md", ".yaml", ".yml", ".toml", ".txt", ".sh", ".cfg", ".json")
SELF = {"tools/sweep.py", "tests/experiments/test_final_sweep.py", "tests/experiments/test_repo_hygiene.py"}
# Frozen records and history keep old names by design (DOC4 FROZEN_GLOBS, in short).
FROZEN = re.compile(r"^(docs/research/|docs/changelog\.md|results/|experiments/phase[0]/|experiments/[^/]+/(RESULTS|PLAN|EXECUTION|"
                    r"RESEARCH|JUDGE_|V0_V1)[^/]*\.md|experiments/.*/(lab|impl)/|experiments/phase2/perf/|archive/)")
RETIRED_TREES = tuple("experiments/" + d for d in ("phase0", "phase1/lab", "phase2/lab"))  # split: CT-15 scans this file
SMOLVLM = re.compile(r"(?i)smolvlm")
LFM = re.compile(r"(?i)\blfm|liquidai")
POOL_NAMES = re.compile(r"\b(_Pool|fork_pool|pool_size)\b")
POOL_IDENTS = {"cpus", "quota"}  # CLEAN16 step 1: AST identifiers, not words in prose
PROCS_FILES = ("src/objdel/cli/train/select.py", "src/objdel/cli/train/probe.py", "src/objdel/cli/eval/run.py",
               "src/objdel/infra/preflight.py", "src/objdel/cli/eval/ceiling.py")  # build_eval_views keeps its own --procs
CONTROLLER = re.compile(r"\b(Controller|alarms_for)\b")
LINK = re.compile(r"(?<!!)\[[^\]]*\]\((<[^>]+>|[^)\s]+)(?:\s+\"[^\"]*\")?\)|!\[[^\]]*\]\((<[^>]+>|[^)\s]+)\)")
TICK = re.compile(r"`((?:src|tools|ops|results|experiments|tests|cuda|docs|archive)/[^`\s*<{]+|ROADMAP\.md)`")
DOTTED = re.compile(r"\b(?:objdel|tools)(?:\.[A-Za-z_]\w*)+")
PYPATH = re.compile(r"(?:src/objdel|tools)/[\w/]+\.py")
D1 = ["D100", "D101", "D102", "D103", "D104", "D106"]

Finding = str
Check = tuple[str, bool, Callable[..., list[Finding]]]  # name, hard, fn(ctx) -> findings


class Ctx:
    """The tree under test: its public files and a cache of their text."""

    def __init__(self, root: Path, files: list[str] | None = None, before: Path | None = None, write_after=False):
        self.root, self.before, self.write_after = root, before, write_after
        self.files = sorted(files if files is not None else public_files(root))
        self._text: dict[str, str] = {}

    def text(self, rel: str) -> str:
        if rel not in self._text:
            try:
                self._text[rel] = (self.root / rel).read_text(errors="ignore")
            except OSError:
                self._text[rel] = ""
        return self._text[rel]

    def texts(self, pred: Callable[[str], bool]) -> list[str]:
        return [f for f in self.files if f.endswith(TEXT + ("Makefile",)) or "Dockerfile" in f if pred(f)]


def public_files(root: Path) -> list[str]:
    return [f for f in git_files(root, "-co", "--exclude-standard") if f and (root / f).is_file()]


def live(f: str) -> bool:
    return f not in SELF and not FROZEN.match(f)


def grep_rule(ctx: Ctx, pat: re.Pattern, scope: Callable[[str], bool], label: str) -> list[Finding]:
    out = []
    for f in ctx.texts(lambda f: live(f) and scope(f)):
        for n, line in enumerate(ctx.text(f).splitlines(), 1):
            if pat.search(line):
                out.append(f"{f}:{n}: {label}: {line.strip()[:100]}")
    return out


def under(*prefixes: str) -> Callable[[str], bool]:
    return lambda f: f.startswith(prefixes)


def retired_trees(ctx: Ctx) -> list[Finding]:
    """A live file names a retired tree (CLEAN-G5, CT-15), or one is still tracked."""
    pat = re.compile("|".join(map(re.escape, RETIRED_TREES)))
    out = grep_rule(ctx, pat, lambda f: True, "retired tree")
    return out + [f"{f}: file in a retired tree" for f in ctx.files if f.startswith(RETIRED_TREES[1:]) or
                  f.startswith(RETIRED_TREES[0] + "/")]


def _model_scope(f: str) -> bool:
    return f.startswith(("src/", "tests/", "tools/"))


def retired_256m(ctx: Ctx) -> list[Finding]:
    """SmolVLM names in src, tools or tests (CLEAN-G10, Wave 1)."""
    return grep_rule(ctx, SMOLVLM, _model_scope, "retired model SmolVLM")


def retired_lfm(ctx: Ctx) -> list[Finding]:
    """LFM / LiquidAI names in src, tools or tests (CLEAN-G10, Wave 3, after CLEAN8)."""
    return grep_rule(ctx, LFM, _model_scope, "retired model LFM")


def controller(ctx: Ctx) -> list[Finding]:
    """`Controller` or `alarms_for` in src or tools (CLEAN-G11)."""
    return grep_rule(ctx, CONTROLLER, under("src/", "tools/"), "retired Controller")


def _idents(tree: ast.AST):
    for n in ast.walk(tree):
        names = [n.id] if isinstance(n, ast.Name) else [n.attr] if isinstance(n, ast.Attribute) else \
            [n.arg] if isinstance(n, ast.arg) else [n.name] if isinstance(n, (ast.FunctionDef, ast.ClassDef)) else \
            [n.arg] if isinstance(n, ast.keyword) and n.arg else []
        for name in names:
            yield name, getattr(n, "lineno", getattr(getattr(n, "value", None), "lineno", 0))


def pool_leftovers(ctx: Ctx) -> list[Finding]:
    """CLEAN16 step 1: pool names, `--procs` on the CLIs that lost it, `procs` on the stage classes."""
    out = grep_rule(ctx, POOL_NAMES, lambda f: f.endswith(".py"), "pool name")
    out += grep_rule(ctx, re.compile(r"--procs\b"), lambda f: f in PROCS_FILES, "--procs flag")
    for f in ctx.texts(lambda f: live(f) and f.endswith(".py") and f.startswith(("src/", "tools/"))):
        try:
            tree = ast.parse(ctx.text(f))
        except SyntaxError:
            out.append(f"{f}: does not parse")
            continue
        out += [f"{f}:{ln}: pool identifier `{name}`" for name, ln in _idents(tree) if name in POOL_IDENTS]
        for c in ast.walk(tree):
            if isinstance(c, ast.ClassDef) and c.name in ("SelectStage", "ProbeStage"):
                out += [f"{f}:{s.lineno}: {c.name}.procs" for s in c.body
                        if isinstance(s, ast.AnnAssign) and getattr(s.target, "id", "") == "procs"]
    return out


def _resolve(ctx: Ctx, base: str, target: str) -> str | None:
    """Repo-relative path of a link target, or None when the target is not a local path."""
    target = target.strip("<>")
    if re.match(r"^(\w+:|#|/)", target):
        return None
    path = urllib.parse.unquote(target.split("#", 1)[0].split("?", 1)[0])
    if not path:
        return None
    return Path(*[p for p in (Path(base).parent / path).parts]).as_posix() if ".." not in path else \
        _norm(Path(base).parent / path)


def _norm(p: Path) -> str:
    parts: list[str] = []
    for part in p.parts:
        if part == "..":
            if parts:
                parts.pop()
        elif part != ".":
            parts.append(part)
    return "/".join(parts)


def _present(rel: str, public: set[str]) -> bool:
    rel = rel.rstrip("/")
    return rel in public or any(p.startswith(rel + "/") for p in public)


def _ignored(ctx: Ctx, rels: list[str]) -> set[str]:
    """Which of these paths git-ignores (a fresh worktree lacks the ignored files, so existence cannot tell)."""
    if not rels:
        return set()
    r = subprocess.run(["git", "-C", str(ctx.root), "check-ignore", "-z", "--stdin"], input="\0".join(rels),
                       capture_output=True, text=True)
    return set(filter(None, r.stdout.split("\0")))


def _doc_refs(ctx: Ctx, pattern: re.Pattern, tick: bool) -> list[Finding]:
    """Refs in public docs to files that are not public: `local-only` (git-ignored, or present but untracked) or `missing`."""
    public, bad = set(ctx.files), []
    for f in ctx.texts(lambda f: f.endswith(".md") and live(f)):
        for n, line in enumerate(ctx.text(f).splitlines(), 1):
            for m in pattern.finditer(line):
                raw = m.group(1) if tick else (m.group(1) or m.group(2))
                if tick:
                    raw = re.sub(r"(::\w+|:\d+(-\d+)?|[.,;)]+)$", "", raw)
                    if "..." in raw or "<" in raw:
                        continue
                rel = raw if tick else _resolve(ctx, f, raw)
                if rel and not _present(rel, public):
                    bad.append((f, n, raw, rel.rstrip("/")))
    ign = _ignored(ctx, sorted({b[3] for b in bad}))
    return [f"{f}:{n}: {'local-only' if rel in ign or (ctx.root / rel).exists() else 'missing'}: {raw}"
            for f, n, raw, rel in bad]


def doc_links(ctx: Ctx) -> list[Finding]:
    """Markdown links in public docs that point to a git-ignored, untracked or missing file (DOC13 step 1)."""
    return _doc_refs(ctx, LINK, tick=False)


def doc_mentions(ctx: Ctx) -> list[Finding]:
    """Backticked repo paths in public docs that are local-only or missing (DOC4 check 2; advisory: some are planned)."""
    return _doc_refs(ctx, TICK, tick=True)


def _modules(ctx: Ctx) -> dict[str, str]:
    mods = {}
    for f in ctx.files:
        if f.endswith(".py") and f.startswith(("src/objdel/", "tools/")) and not FROZEN.match(f) and f not in SELF:
            name = f.removeprefix("src/").removesuffix(".py").replace("/", ".")
            if not name.endswith(("__init__", "__main__")):
                mods[name] = f
    return mods


def _imported(tree: ast.AST, pkg: str) -> set[str]:
    out = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            out |= {a.name for a in n.names}
        elif isinstance(n, ast.ImportFrom):
            base = n.module or ""
            if n.level:
                head = pkg.split(".")[: len(pkg.split(".")) - n.level + 1]
                base = ".".join(head + ([n.module] if n.module else []))
            out |= {base} | {f"{base}.{a.name}" for a in n.names}
    return out


def importers(ctx: Ctx) -> tuple[dict[str, set[str]], dict[str, str]]:
    """Module name -> files that import it or name it (dotted name or `.py` path in any text file); tests included."""
    mods = _modules(ctx)
    by_path = {f: m for m, f in mods.items()}
    refs: dict[str, set[str]] = {}
    for f in ctx.texts(lambda f: f not in SELF):
        text = ctx.text(f)
        seen: set[str] = set()
        if f.endswith(".py"):
            try:
                pkg = f.removeprefix("src/").removesuffix(".py").replace("/", ".")
                seen |= _imported(ast.parse(text), pkg.rsplit(".", 1)[0] if "." in pkg else pkg)
            except SyntaxError:
                pass
        seen |= set(DOTTED.findall(text))
        seen |= {by_path[p] for p in PYPATH.findall(text) if p in by_path}
        for m in seen & mods.keys():
            if mods[m] != f:
                refs.setdefault(m, set()).add(f)
    return refs, mods


def unimported(ctx: Ctx) -> list[Finding]:
    """Modules of src/objdel and tools that no file imports or names, and those only tests use (CLEAN16 step 2).

    A test-only module that says `Test surface:` in its docstring has its decision recorded and is not reported."""
    refs, mods = importers(ctx)
    out = []
    for m, f in sorted(mods.items()):
        r = refs.get(m, set())
        if not r:
            out.append(f"{f}: nothing imports or names it")
        elif all(x.startswith("tests/") for x in r) and "Test surface:" not in ctx.text(f)[:3000]:
            out.append(f"{f}: only tests use it ({len(r)} file(s)): keep as test surface with a docstring, or archive")
    return out


def docstrings(ctx: Ctx) -> list[Finding]:
    """DOC5: D100-D106 are selected in pyproject.toml and `src/objdel` has none missing."""
    out = []
    sel = tomllib.loads((ctx.root / "pyproject.toml").read_text())["tool"]["ruff"]["lint"].get("select", [])
    if missing := [d for d in D1 if d not in sel and "D" not in sel]:
        out.append(f"pyproject.toml: ruff select lacks {','.join(missing)} (the DOC5 gate is not on)")
    r = subprocess.run([sys.executable, "-m", "ruff", "check", "src/objdel", "--select", ",".join(D1), "--no-cache",
                        "--output-format", "json"], cwd=ctx.root, capture_output=True, text=True)
    try:
        items = json.loads(r.stdout or "[]")
    except json.JSONDecodeError:
        return out + [f"ruff gave no JSON: {r.stderr.strip()[:200]}"]
    for i in items:
        rel = Path(i["filename"]).resolve().relative_to(ctx.root.resolve()).as_posix()
        out.append(f"{rel}:{i['location']['row']}: {i['code']} {i['message']}")
    return out


def _vulture(ctx: Ctx, *paths: str) -> tuple[int, list[Finding]]:
    r = subprocess.run([sys.executable, "-m", "vulture", *paths], cwd=ctx.root, capture_output=True, text=True)
    lines = [ln for ln in r.stdout.splitlines() if ln.strip()]
    if r.returncode not in (0, 3):
        lines.append(f"vulture failed (exit {r.returncode}): {r.stderr.strip()[:200]}")
    return r.returncode, lines


def vulture_gate(ctx: Ctx) -> list[Finding]:
    """`make dead`: vulture with the pyproject config over src, tools and tests prints nothing (CLEAN-G9)."""
    return _vulture(ctx)[1]


def vulture_src(ctx: Ctx) -> list[Finding]:
    """vulture over src and tools without tests: code that only tests use needs a decision (CLEAN16 step 2)."""
    return _vulture(ctx, "src", "tools")[1]


def record(root: Path) -> dict:
    """The CLEAN0 step 2 record of the tree now (fingerprint and counts); a failed import leaves the key out."""
    rec: dict = {}
    tracked = git_files(root)
    rec.update(tracked=len(tracked), lab=sum("/lab/" in x for x in tracked),
               phase0=sum(x.startswith(RETIRED_TREES[0] + "/") for x in tracked))
    try:
        import objdel.verify.reward.spec as spec

        rec["fingerprint"] = spec.fingerprint()
    except Exception as e:  # noqa: BLE001 - the record is a report, not a gate on import health
        rec["fingerprint_error"] = repr(e)[:120]
    try:
        from objdel.infra.vast import remote

        files = remote.code_files()
        rec.update(code_files=len(files), claude=sum(x.startswith(".claude/") for x in files),
                   archive=sum(x.startswith("archive/") for x in files))
    except Exception as e:  # noqa: BLE001
        rec["code_files_error"] = repr(e)[:120]
    return rec


def compare(before: dict, after: dict) -> list[Finding]:
    """Hard rules of the end-of-wave comparison: G3 fingerprint equal, G5 no lab or phase0 file, G7 bundle clean."""
    out = []
    if "fingerprint" not in after:
        out.append(f"no fingerprint now: {after.get('fingerprint_error')}")
    elif before.get("fingerprint") and before["fingerprint"] != after["fingerprint"]:
        out.append(f"fingerprint changed: {before['fingerprint'][:12]} -> {after['fingerprint'][:12]} "
                   "(expected only where B1.0/B3 changed the reward; compare with the record of this wave's start)")
    out += [f"{k} = {after[k]} (must be 0)" for k in ("lab", "phase0", "claude", "archive") if after.get(k)]
    out += [f"{k} not measurable: {after[k + '_error']}" for k in ("code_files",) if k + "_error" in after]
    return out


def clean_record(ctx: Ctx) -> list[Finding]:
    """CLEAN16 step 3: re-run CLEAN0's record into runs/clean/after.json (with --write-after) and compare."""
    after = record(ctx.root)
    path = ctx.before
    if path is None:
        found = [ctx.root / "runs/clean" / n for n in ("before_w3.json", "before_w2.json", "before.json")]
        path = next((p for p in found if p.is_file()), None)
    if ctx.write_after:
        dest = ctx.root / "runs/clean/after.json"
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(after, indent=1))
    out = compare(read_json(path), after) if path else ["no runs/clean/before*.json: nothing to compare"]
    if path:
        b = read_json(path)
        out.append(f"info: vs {path.name}: tracked {b.get('tracked')} -> {after['tracked']}, "
                   f"code_files {b.get('code_files')} -> {after.get('code_files')}")
    return out


WAVE3 = {"retired_lfm", "pool_leftovers"}  # CLEAN-G10 (lfm, Wave 3) and CLEAN16 step 1 (D7, E3b): late by design

CHECKS: list[Check] = [
    ("retired_trees", True, retired_trees),
    ("retired_256m", True, retired_256m),
    ("retired_lfm", True, retired_lfm),
    ("controller", True, controller),
    ("pool_leftovers", True, pool_leftovers),
    ("doc_links", True, doc_links),
    ("docstrings", True, docstrings),
    ("vulture_gate", True, vulture_gate),
    ("clean_record", True, clean_record),
    ("doc_mentions", False, doc_mentions),
    ("unimported", False, unimported),
    ("vulture_src", False, vulture_src),
]


def run(ctx: Ctx, only: set[str] | None = None, skip: set[str] = frozenset(), wave: int = 3) -> list[dict]:
    """Run the selected checks; each result has name, hard, findings. An `info:` finding never fails a hard check.

    Before wave 3 the WAVE3 checks are advisory.
    """
    out = []
    for name, hard, fn in CHECKS:
        if (only and name not in only) or name in skip:
            continue
        hard = hard and (wave >= 3 or name not in WAVE3)
        found = fn(ctx)
        out.append({"name": name, "hard": hard, "findings": found,
                    "failed": hard and any(not x.startswith("info:") for x in found)})
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--only", default="", help="comma list of check names")
    ap.add_argument("--skip", default="", help="comma list of check names")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--limit", type=int, default=15, help="findings printed per check (0 = all)")
    ap.add_argument("--before", type=Path, help="CLEAN0 record to compare with (default: newest runs/clean/before*.json)")
    ap.add_argument("--wave", type=int, choices=(1, 2, 3), default=3, help="before wave 3 the lfm and pool checks are advisory")
    ap.add_argument("--write-after", action="store_true", help="write runs/clean/after.json")
    a = ap.parse_args(argv)
    names = {n for n, _, _ in CHECKS}
    only, skip = set(filter(None, a.only.split(","))), set(filter(None, a.skip.split(",")))
    if bad := (only | skip) - names:
        ap.error(f"unknown check(s) {sorted(bad)}; known: {sorted(names)}")
    results = run(Ctx(ROOT, before=a.before, write_after=a.write_after), only or None, skip, a.wave)
    if a.json:
        print(json.dumps(results, indent=1))
    else:
        for r in results:
            tag = "FAIL" if r["failed"] else ("PASS" if r["hard"] else "note")
            real = [x for x in r["findings"] if not x.startswith("info:")]
            print(f"[{tag}] {r['name']}: {len(real)} finding(s)" + ("" if r["hard"] else " (advisory)"))
            shown = r["findings"] if not a.limit else r["findings"][: a.limit]
            for x in shown:
                print("    " + x)
            if len(r["findings"]) > len(shown):
                print(f"    ... {len(r['findings']) - len(shown)} more (--limit 0 prints all)")
        print("SWEEP " + ("FAILED: " + ", ".join(r["name"] for r in results if r["failed"])
                          if any(r["failed"] for r in results) else "PASSED"))
    return 1 if any(r["failed"] for r in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
