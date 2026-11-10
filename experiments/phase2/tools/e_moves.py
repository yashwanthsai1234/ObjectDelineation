"""Phase E file moves as one idempotent command (Implementation.md E3a, E3b, E3c, E5; rule R-1, no shims R-2).

Run from the repo root:

    .venv/bin/python -I -m experiments.phase2.tools.e_moves --dry-run   # print the plan, change nothing
    .venv/bin/python -I -m experiments.phase2.tools.e_moves             # apply it
    .venv/bin/python -I -m experiments.phase2.tools.e_moves --stage E3a # only the moves of one sub-move

Stages (one commit each, in this order): E3a geom, E3b verify (and tests/oracles), E3c data, E5 the rest (policy, rl,
sft). A stage applies every table row whose new module sits in that layer or an earlier one; a split module whose
names are only partly due keeps the rest, and imports from the new modules the names its rest still uses. The old
package rows (`drop_reexports`, `delete_packages`) carry their stage in `stages`.

The move table is `e_moves.json` next to this file. Steps, in order:

1. Resolve. Build the old -> new map of modules: whole-file moves, split modules (one old module -> several new
   modules by an explicit name map) and the re-exports of each old package `__init__` (a caller of a re-exported
   name is sent to the defining module, R-1). Every top-level name of a split module must have a row, or the script
   stops and names the missing ones (a module that another lane added after the table was written).
2. Rewrite imports in every `.py` under src, tests and experiments. The stdlib `ast` finds each import statement and
   its exact span; only that span is replaced (no regex over code). A module alias that pointed at a split module or
   at an old package (`from objdel.train import dataset as D`, `D.prompts`) is rewritten per attribute to the new
   module (`policy_rows.prompts`); `monkeypatch.setattr(D, "name", ...)` follows the same rule.
3. Rewrite strings: dotted names (`"objdel.train.fast.FastGRPO"`, `"objdel.eval.metrics"`) and repo paths
   (`src/objdel/eval/gates.py`) in code and in the text files of `text_extra` (docs, Makefile, pyproject).
4. Split: cut each split module by its top-level statements (with decorators and the comment lines above them),
   write each piece to its new module with only the imports it uses plus imports of moved sibling names.
5. `git mv` the whole-file moves, remove the re-exports of the old package `__init__` files, delete emptied old
   packages, create the new package `__init__` files and `geom/engine.py` (the only importer of `objdel._gpucheck`).
6. `ruff check --fix --select I001,F401` on the touched files, wrap the comment and docstring lines the longer
   names pushed past 120 columns; then report every old name and every long code line still found
   (the part that needs a hand edit).

Idempotent: a second run finds no old module, no old import and nothing to move, and prints "nothing to do".
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

TABLE = Path(__file__).with_name("e_moves.json")
STAGES = ("E3a", "E3b", "E3c", "E5")
LAYER_STAGE = (("objdel.geom", "E3a"), ("objdel.verify", "E3b"), ("tests.oracles", "E3b"), ("objdel.data", "E3c"))
SKIP_DIRS = {".venv", "build", "__pycache__", ".git", "data", "runs", "archive", "node_modules"}


def mod_of(path: str) -> str:
    """Module name of a repo path: src/objdel/a/b.py -> objdel.a.b; a tests/ path maps the same way."""
    p = path[4:] if path.startswith("src/") else path
    p = p[: -len(".py")] if p.endswith(".py") else p
    p = p[: -len("/__init__")] if p.endswith("/__init__") else p
    return p.replace("/", ".")


def path_of(mod: str) -> str:
    """Repo path of a module name (file form), the inverse of mod_of for src/objdel modules."""
    return ("src/" if mod.startswith("objdel") else "") + mod.replace(".", "/") + ".py"


def top_names(node: ast.stmt) -> list[str]:
    """Names a top-level statement defines (def, class, assignment); [] for anything else."""
    if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
        return [node.name]
    if isinstance(node, ast.Assign):
        out = []
        for t in node.targets:
            out += [n.id for n in ast.walk(t) if isinstance(n, ast.Name)]
        return out
    if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
        return [node.target.id]
    return []


def stage_of(target: str) -> str:
    """Sub-move of a new module or path: the layer it lands in (geom E3a, verify E3b, data E3c, else E5)."""
    m = mod_of(target) if "/" in target else target
    return next((st for pre, st in LAYER_STAGE if m == pre or m.startswith(pre + ".")), "E5")


def stage_table(table: dict, stage: str) -> dict:
    """The rows of `table` due at `stage` (that stage and the earlier ones). Split modules keep every name in
    `split_names` (the missing-row check reads the whole table)."""
    due = set(STAGES[: STAGES.index(stage) + 1])
    stages = table.get("stages", {})
    out = dict(table)
    out["files"] = {o: n for o, n in table["files"].items() if stage_of(n) in due}
    out["splits"] = {s: {t: ns for t, ns in tg.items() if stage_of(t) in due} for s, tg in table["splits"].items()}
    out["split_names"] = {s: [n for ns in tg.values() for n in ns] for s, tg in table["splits"].items()}
    out["drop_reexports"] = [i for i in table["drop_reexports"] if stages.get(i, "E5") in due]
    out["delete_packages"] = [d for d in table["delete_packages"] if stages.get(d, "E5") in due]
    out["packages"] = {p: c for p, c in table["packages"].items() if stage_of(p) in due}
    return out


@dataclass
class Plan:
    """The resolved old -> new map, computed from the table and the live tree."""

    root: Path
    table: dict
    modmap: dict[str, str] = field(default_factory=dict)  # old module -> new module (whole moves, packages)
    symmap: dict[str, dict[str, str]] = field(default_factory=dict)  # old module -> {name: new defining module}
    split_mods: set[str] = field(default_factory=set)
    partial: dict[str, list[str]] = field(default_factory=dict)  # split module path -> names that stay this stage
    old_pkgs: set[str] = field(default_factory=set)  # packages that disappear (eval, train, checker, ...)
    moves: list[tuple[str, str]] = field(default_factory=list)  # (old path, new path) still to do
    splits: dict[str, dict[str, list[str]]] = field(default_factory=dict)  # still to do
    errors: list[str] = field(default_factory=list)
    errors_soft: list[str] = field(default_factory=list)


def build_plan(root: Path, table: dict) -> Plan:
    plan = Plan(root, table)
    for old, new in table["files"].items():
        plan.modmap[mod_of(old)] = mod_of(new)
        if (root / old).exists() and not (root / new).exists():
            plan.moves.append((old, new))
        elif (root / old).exists() and (root / new).exists():
            plan.errors.append(f"both {old} and {new} exist: merge them by hand, then re-run")
    for pkg in table["delete_packages"]:
        plan.old_pkgs.add(mod_of(pkg))
    # packages moved whole keep a package name (reward -> verify.reward, targets -> data.targets)
    for old, new in table["files"].items():
        if old.endswith("/__init__.py"):
            plan.modmap[mod_of(old)] = mod_of(new)
    for src, targets in table["splits"].items():
        m = mod_of(src)
        f = root / src
        if not f.exists():
            continue
        name_to = {n: t for t, names in targets.items() for n in names}
        tree = ast.parse(f.read_text())
        defined = [n for node in tree.body for n in top_names(node)]
        rows = set(table.get("split_names", {}).get(src, name_to))
        missing = sorted({n for n in defined if n not in rows and n != "__all__"})
        if missing:
            plan.errors.append(f"{src}: no row in e_moves.json for {missing}; add each name to a new module")
        due = {n for n in defined if n in name_to}
        if not due:
            continue  # nothing of this module moves at this stage
        plan.split_mods.add(m)
        plan.symmap[m] = name_to
        plan.splits[src] = targets
        stay = sorted({n for n in defined if n not in name_to and n != "__all__"})
        if stay:
            plan.partial[src] = stay
            name_to.update(dict.fromkeys(stay, m))  # names that stay resolve to the old module
            continue
        # names the old module imported (patched by tests as `fast.time`): home = the one new module that uses them
        users: dict[str, set[str]] = {}
        for node in tree.body:
            for n in top_names(node):
                if n in name_to:
                    for u in used_names([node]):
                        users.setdefault(u, set()).add(name_to[n])
        for node in tree.body:
            for b in bound_names(node):
                homes = users.get(b, set())
                if b not in name_to and homes:
                    name_to[b] = sorted(homes)[0]
                    if len(homes) > 1:
                        plan.errors_soft.append(f"{src}: imported name `{b}` is used by {sorted(homes)}; "
                                                f"`<alias>.{b}` goes to {name_to[b]}")
    # re-exports of old package __init__ files: name -> defining module (then mapped to its new name)
    for init in table["drop_reexports"]:
        f = root / init
        if not f.exists():
            continue
        pkg = mod_of(init)
        tree = ast.parse(f.read_text())
        exports: dict[str, str] = {}
        for node in tree.body:
            if isinstance(node, ast.ImportFrom):
                base = resolve_from(pkg, True, node)
                for a in node.names:
                    sub = f"{base}.{a.name}"
                    if module_exists(root, sub, plan):
                        continue  # `from . import gates`: a submodule, resolved as a module
                    exports[a.asname or a.name] = final_symbol_module(plan, base, a.name)
        if exports:  # an __init__ with no re-exports is not an old module (keeps a second run a no-op)
            plan.symmap.setdefault(pkg, {}).update(exports)
    return plan


def module_exists(root: Path, mod: str, plan: Plan) -> bool:
    """True if `mod` is a module or package in the old tree or a known old module of the table."""
    if mod in plan.modmap or mod in plan.split_mods or mod in plan.old_pkgs:
        return True
    p = root / path_of(mod)
    return p.exists() or (p.with_suffix("") / "__init__.py").exists()


def resolve_from(cur_mod: str, is_pkg: bool, node: ast.ImportFrom) -> str:
    """Absolute module of an ImportFrom (relative levels resolved against the importing module)."""
    if not node.level:
        return node.module or ""
    parts = cur_mod.split(".")
    base = parts if is_pkg else parts[:-1]
    base = base[: len(base) - (node.level - 1)]
    return ".".join(base + ([node.module] if node.module else []))


def final_symbol_module(plan: Plan, mod: str, name: str) -> str:
    """New defining module of `mod.name` (follows splits, re-exports and moves)."""
    if mod in plan.symmap and name in plan.symmap[mod]:
        return plan.symmap[mod][name]
    return plan.modmap.get(mod, mod)


def map_module(plan: Plan, mod: str) -> str:
    """New name of a module (whole move or a module inside a moved package)."""
    if mod in plan.modmap:
        return plan.modmap[mod]
    parts = mod.split(".")
    for i in range(len(parts) - 1, 0, -1):
        pre = ".".join(parts[:i])
        if pre in plan.modmap and pre + ".__init__" not in plan.modmap:
            return plan.modmap[pre] + "." + ".".join(parts[i:])
    return mod


def is_old(plan: Plan, mod: str) -> bool:
    return mod in plan.modmap or mod in plan.symmap or mod in plan.old_pkgs


# ---------------------------------------------------------------- import rewrite (ast spans, no regex on code)


def line_offsets(data: bytes) -> list[int]:
    offs, pos = [0], 0
    for line in data.splitlines(keepends=True):
        pos += len(line)
        offs.append(pos)
    return offs


def span(offs: list[int], node: ast.AST) -> tuple[int, int]:
    return offs[node.lineno - 1] + node.col_offset, offs[node.end_lineno - 1] + node.end_col_offset


def alias_for(mod: str) -> str:
    """Generated alias for a new module: objdel.policy.rows -> policy_rows."""
    parts = mod.split(".")
    return "_".join(parts[-2:]) if len(parts) > 2 else parts[-1]


def from_line(mod: str, names: list[tuple[str, str | None]]) -> str:
    items = ", ".join(n if not a or a == n else f"{n} as {a}" for n, a in names)
    return f"from {mod} import {items}"


def module_line(mod: str, alias: str) -> str:
    parent, _, leaf = mod.rpartition(".")
    if parent and leaf == alias:
        return f"from {parent} import {leaf}"
    if parent:
        return f"from {parent} import {leaf} as {alias}"
    return f"import {mod} as {alias}"


def rewrite_python(plan: Plan, rel: str, text: str, keep_gpucheck: bool) -> tuple[str, list[str]]:
    """Rewrite the import statements and module-alias attributes of one file. Returns (text, notes)."""
    try:
        tree = ast.parse(text)
    except SyntaxError as e:
        return text, [f"{rel}: not parsed ({e.msg})"]
    data = text.encode()
    offs = line_offsets(data)
    cur = mod_of(rel)
    is_pkg = rel.endswith("__init__.py")
    moving = rel in dict(plan.moves) or rel in plan.splits
    edits: list[tuple[int, int, str]] = []
    notes: list[str] = []
    split_aliases: dict[str, str] = {}  # local alias -> old split module / old package
    pending: list = []  # (start, end, node, lines, split aliases bound there): rendered once aliases are known
    parents: dict[int, list] = {}  # id(statement) -> the body list that holds it
    for blk in ast.walk(tree):
        for fld in ("body", "orelse", "finalbody", "handlers"):
            for child in (getattr(blk, fld, None) if isinstance(getattr(blk, fld, None), list) else []):
                if isinstance(child, ast.stmt):
                    parents[id(child)] = getattr(blk, fld)

    def indent_join(node: ast.stmt, lines: list[str]) -> str:
        line = data[offs[node.lineno - 1] : offs[node.lineno]].decode()
        lead = line[: len(line) - len(line.lstrip())]
        own_line = line.strip().startswith(("from", "import"))
        if not own_line:
            return "; ".join(lines)
        tail = data[offs[node.end_lineno - 1] + node.end_col_offset : offs[node.end_lineno]].decode().strip()
        noqa = re.search(r"#\s*noqa(:[\sA-Z0-9,]+)?", tail)
        if noqa and len(lines) > 1:  # every line of a split statement keeps the noqa comment of the old line
            lines = [ln + "  " + noqa.group(0).strip() for ln in lines[:-1]] + lines[-1:]
        return ("\n" + lead).join(lines)

    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            base = resolve_from(cur, is_pkg, node)
            if not keep_gpucheck and base == "objdel._gpucheck":
                edits.append((*span(offs, node), from_line("objdel.geom.engine", [(a.name, a.asname) for a in node.names])))
                continue
            if not keep_gpucheck and base == "objdel" and any(a.name == "_gpucheck" for a in node.names):
                lines = []
                rest = [(a.name, a.asname) for a in node.names if a.name != "_gpucheck"]
                for a in node.names:
                    if a.name == "_gpucheck":
                        lines.append(module_line("objdel.geom.engine", a.asname or "_gpucheck"))
                if rest:
                    lines.insert(0, from_line("objdel", rest))
                edits.append((*span(offs, node), indent_join(node, lines)))
                continue
            touched = is_old(plan, base) or (node.level and moving) or any(is_old(plan, f"{base}.{a.name}") for a in node.names)
            if not touched:
                continue
            groups: dict[str, list[tuple[str, str | None]]] = {}
            mods: list[str] = []
            node_aliases: list[str] = []
            for a in node.names:
                sub = f"{base}.{a.name}"
                local = a.asname or a.name
                reexport = base in plan.symmap and a.name in plan.symmap[base]
                if not reexport and module_exists(plan.root, sub, plan):
                    if sub in plan.split_mods or (sub in plan.old_pkgs and sub in plan.symmap):
                        split_aliases[local] = sub
                        node_aliases.append(local)
                        continue
                    mods.append(module_line(map_module(plan, sub), local))
                    continue
                new = final_symbol_module(plan, base, a.name)
                new = map_module(plan, new)
                groups.setdefault(new, []).append((a.name, a.asname))
            lines = [from_line(m, ns) for m, ns in groups.items()] + mods
            pending.append((*span(offs, node), node, lines, node_aliases))
        elif isinstance(node, ast.Import):
            lines, changed, node_aliases = [], False, []
            for a in node.names:
                if not keep_gpucheck and a.name == "objdel._gpucheck" and a.asname:
                    lines.append(module_line("objdel.geom.engine", a.asname))
                    changed = True
                elif a.name in plan.modmap:  # a package or module moved whole: `import objdel.reward`
                    new = plan.modmap[a.name]
                    lines.append(f"import {new} as {a.asname}" if a.asname else f"import {new}")
                    changed = True
                elif a.name in plan.split_mods or (a.name in plan.old_pkgs and a.name in plan.symmap):
                    if a.asname:
                        split_aliases[a.asname] = a.name
                        node_aliases.append(a.asname)
                    else:  # resolved after the string pass: one `import <new>` per new module the file names
                        lines.append(f"import __E_PLAIN__{a.name}")
                    changed = True
                elif is_old(plan, a.name) or map_module(plan, a.name) != a.name:
                    new = map_module(plan, a.name)
                    lines.append(f"import {new} as {a.asname}" if a.asname else f"import {new}")
                    changed = True
                else:
                    lines.append(f"import {a.name} as {a.asname}" if a.asname else f"import {a.name}")
            if changed:
                pending.append((*span(offs, node), node, lines, node_aliases))

    # module aliases of split modules / old packages: rewrite each `alias.attr` and setattr(alias, "attr")
    needed: dict[str, dict[str, str]] = {}  # local alias -> {new module: alias used}
    for local, old in split_aliases.items():
        uses: list[tuple[ast.AST, str]] = []
        bare = 0
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == local:
                uses.append((node, node.attr))
            elif (isinstance(node, ast.Call) and len(node.args) >= 2 and isinstance(node.args[0], ast.Name)
                  and node.args[0].id == local and isinstance(node.args[1], ast.Constant)
                  and isinstance(node.args[1].value, str)):
                uses.append((node.args[0], node.args[1].value))
        names_attr = {id(n) for n, _ in uses}
        for node in ast.walk(tree):
            if isinstance(node, ast.Name) and node.id == local and id(node) not in names_attr:
                parents_ok = any(isinstance(n, ast.Attribute) and n.value is node for n, _ in uses)
                if not parents_ok:
                    bare += 1
        if bare:
            notes.append(f"{rel}: alias `{local}` of {old} is used without an attribute {bare}x: fix by hand")
        targets: dict[str, str] = {}
        resolved: list[tuple[ast.AST, str, str, bool]] = []
        for n, attr in uses:
            sub = f"{old}.{attr}"
            if old in plan.symmap and attr in plan.symmap[old]:
                resolved.append((n, attr, map_module(plan, plan.symmap[old][attr]), False))
            elif module_exists(plan.root, sub, plan) and sub not in plan.split_mods:
                resolved.append((n, attr, map_module(plan, sub), True))
            else:
                notes.append(f"{rel}:{n.lineno}: `{local}.{attr}`: {old}.{attr} has no row: fix by hand")
        mods = sorted({m for _, _, m, is_mod in resolved if not is_mod})
        for m in mods:
            targets[m] = local if len(mods) == 1 or m == old else alias_for(m)  # the rest of a partial split
        needed[local] = targets
        for n, attr, m, is_mod in resolved:
            if isinstance(n, ast.Attribute):
                if is_mod:
                    ali = alias_for(m)
                    needed[local][m] = ali
                    edits.append((*span(offs, n), ali + "__MODALIAS__"))
                else:
                    edits.append((*span(offs, n.value), targets[m]))
            else:
                if is_mod:
                    notes.append(f"{rel}:{n.lineno}: setattr on submodule `{local}.{attr}`: fix by hand")
                else:
                    edits.append((*span(offs, n), targets[m]))
    # render each rewritten import: its plain lines plus the module imports its split aliases now need
    for st, en, node, lines, aliases in pending:
        extra = [module_line(m, ali) for local in aliases for m, ali in needed.get(local, {}).items()]
        allines = list(dict.fromkeys(lines + extra))
        if allines:
            edits.append((st, en, indent_join(node, allines)))
        elif len(parents.get(id(node), [])) > 1 and data[offs[node.lineno - 1] : st].strip() == b"":
            edits.append((offs[node.lineno - 1], offs[node.end_lineno], ""))  # drop the whole line
        else:
            edits.append((st, en, "pass"))
    out = bytearray(data)
    for s_, e_, new in sorted(edits, key=lambda x: x[0], reverse=True):
        out[s_:e_] = new.replace("__MODALIAS__", "").encode()
    result = out.decode()
    return result, notes


def resolve_plain(plan: Plan, text: str) -> str:
    """Replace each `import __E_PLAIN__<old>` marker by `import <new>` for every new module the file now names."""
    def rep(m: re.Match) -> str:
        lead, old = m.group(1), m.group(2)
        cands = {map_module(plan, v) for v in plan.symmap.get(old, {}).values()}
        cands |= {v for k, v in plan.modmap.items() if k.startswith(old + ".")}
        found = sorted(c for c in cands if re.search(rf"(?<![\w.]){re.escape(c)}\b", text))
        return "\n".join(f"{lead}import {c}" for c in found) or f"{lead}pass"
    return re.sub(r"^([ \t]*)import __E_PLAIN__([\w.]+)", rep, text, flags=re.M)


def add_imports(text: str, lines: list[str]) -> str:
    """Insert import lines after the last top-level import (or after the docstring)."""
    tree = ast.parse(text)
    at = 0
    for node in tree.body:
        if isinstance(node, ast.Import | ast.ImportFrom):
            at = node.end_lineno
        elif at == 0 and isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
            at = node.end_lineno
    src = text.splitlines(keepends=True)
    new = [ln + "\n" for ln in lines if ln + "\n" not in src]
    return "".join(src[:at] + new + src[at:])


# ---------------------------------------------------------------- string rewrite (dotted names and repo paths)

DOTTED = re.compile(r"(?<![\w.])objdel(?:\.[A-Za-z_]\w*)+")
PATHS = re.compile(r"(?<![\w/])(?:src/)?objdel(?:/[A-Za-z_]\w*)+(?:\.py)?")


def map_dotted(plan: Plan, name: str, keep_gpucheck: bool) -> str | None:
    """New spelling of a dotted name, or None when it cannot be resolved (bare split module)."""
    if name.startswith("objdel._gpucheck."):
        attr = name.split(".")[2]
        if keep_gpucheck or attr not in plan.table["engine"]["names"]:
            return name
        return "objdel.geom.engine." + name[len("objdel._gpucheck."):]
    parts = name.split(".")
    for i in range(len(parts), 1, -1):
        pre = ".".join(parts[:i])
        rest = parts[i:]
        if pre in plan.split_mods or (pre in plan.old_pkgs and pre not in plan.modmap):
            if not rest:
                return name if path_of(pre) in plan.partial else None
            if pre in plan.symmap and rest[0] in plan.symmap[pre]:
                return ".".join([map_module(plan, plan.symmap[pre][rest[0]]), *rest])
            if module_exists(plan.root, f"{pre}.{rest[0]}", plan):
                continue
            return None
        if pre in plan.modmap or map_module(plan, pre) != pre:
            if pre in plan.symmap and rest and rest[0] in plan.symmap[pre]:
                return ".".join([map_module(plan, plan.symmap[pre][rest[0]]), *rest])
            return ".".join([map_module(plan, pre), *rest])
    return name


def rewrite_strings(plan: Plan, rel: str, text: str, keep_gpucheck: bool) -> tuple[str, list[str]]:
    notes: list[str] = []

    def dotted(m: re.Match) -> str:
        new = map_dotted(plan, m.group(0), keep_gpucheck)
        if new is None:
            notes.append(f"{rel}: `{m.group(0)}` names a split module or a removed package: fix by hand")
            return m.group(0)
        return new

    def pathed(m: re.Match) -> str:
        s = m.group(0)
        pre = "src/" if s.startswith("src/") else ""
        body = s[len(pre):]
        is_py = body.endswith(".py")
        mod = body[: -3].replace("/", ".") if is_py else body.replace("/", ".")
        new = map_dotted(plan, mod, keep_gpucheck)
        if new is None or (not is_py and new == mod):
            if new is None:
                notes.append(f"{rel}: path `{s}` names a split module or a removed package: fix by hand")
            return s
        if new == mod:
            return s
        return pre + new.replace(".", "/") + (".py" if is_py else "")

    text = DOTTED.sub(dotted, text)
    text = PATHS.sub(pathed, text)
    return text, notes


# ---------------------------------------------------------------- split one module by its name map


def chunk_start(lines: list[str], node: ast.stmt) -> int:
    """0-based first line of a statement, including decorators and the comment lines right above it."""
    start = min([node.lineno] + [d.lineno for d in getattr(node, "decorator_list", [])]) - 1
    while start > 0 and lines[start - 1].lstrip().startswith("#"):
        start -= 1
    return start


def bound_names(node: ast.stmt) -> list[str]:
    if isinstance(node, ast.Import):
        return [a.asname or a.name.split(".")[0] for a in node.names]
    if isinstance(node, ast.ImportFrom):
        return [a.asname or a.name for a in node.names]
    return []


def used_names(nodes: list[ast.stmt]) -> set[str]:
    out = set()
    for n in nodes:
        for x in ast.walk(n):
            if isinstance(x, ast.Name):
                out.add(x.id)
            elif isinstance(x, ast.Constant) and isinstance(x.value, str) and x.value.isidentifier():
                out.add(x.value)  # string annotations
    return out


def trimmed_import(node: ast.stmt, keep: set[str]) -> str | None:
    if isinstance(node, ast.ImportFrom):
        names = [(a.name, a.asname) for a in node.names if (a.asname or a.name) in keep]
        if not names:
            return None
        mod = "." * node.level + (node.module or "")
        return from_line(mod, names)
    names = [a for a in node.names if (a.asname or a.name.split(".")[0]) in keep]
    if not names:
        return None
    return "\n".join(f"import {a.name} as {a.asname}" if a.asname else f"import {a.name}" for a in names)


def is_type_checking(node: ast.stmt) -> bool:
    """A top-level `if TYPE_CHECKING:` block that holds only imports (an import group of its own)."""
    test = getattr(node, "test", None)
    name = test.id if isinstance(test, ast.Name) else getattr(test, "attr", "")
    return (isinstance(node, ast.If) and name == "TYPE_CHECKING" and not node.orelse
            and all(isinstance(i, ast.Import | ast.ImportFrom) for i in node.body))


def split_module(plan: Plan, rel: str, targets: dict[str, list[str]], apply: bool) -> list[str]:
    """Cut one old module into its new modules. Returns the touched paths."""
    root = plan.root
    text = (root / rel).read_text()
    lines = text.splitlines(keepends=True)
    tree = ast.parse(text)
    old_mod = mod_of(rel)
    stay = plan.partial.get(rel, [])
    name_to = {n: t for t, names in targets.items() for n in names} | dict.fromkeys(stay, old_mod)
    imports = [n for n in tree.body if isinstance(n, ast.Import | ast.ImportFrom)]
    future = [n for n in imports if isinstance(n, ast.ImportFrom) and n.module == "__future__"]
    guarded = [i for n in tree.body if is_type_checking(n) for i in n.body]  # `if TYPE_CHECKING:` imports
    chunks: dict[str, list[tuple[ast.stmt, str]]] = {}
    body = [n for n in tree.body if not isinstance(n, ast.Import | ast.ImportFrom) and not is_type_checking(n)]
    doc = ast.get_docstring(tree)
    for node in body:
        names = top_names(node)
        if not names:
            if node is tree.body[0] and doc is not None:
                continue
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
                continue
            calls = {name_to[n] for n in used_names([node]) if n in name_to}  # `if __name__ == "__main__": main()`
            if isinstance(node, ast.If) and len(calls) == 1:
                names = [n for n in used_names([node]) if n in name_to]
        if not names:
            raise SystemExit(f"{rel}:{node.lineno}: top-level statement with no name ({type(node).__name__}): "
                             "split it by hand")
        homes = {name_to[n] for n in names if n != "__all__"}
        if not homes or homes == {old_mod}:
            continue
        if len(homes) > 1:
            raise SystemExit(f"{rel}:{node.lineno}: one statement binds names of two new modules {sorted(homes)}")
        start = chunk_start(lines, node)
        chunks.setdefault(homes.pop(), []).append((node, "".join(lines[start : node.end_lineno])))
    touched = []
    for target, items in chunks.items():
        nodes = [n for n, _ in items]
        used = used_names(nodes)
        own = {n for node in nodes for n in top_names(node)}
        header = []
        self_alias = None
        typed = [ln for imp in guarded if (ln := trimmed_import(imp, used - own))]
        if typed:
            used.add("TYPE_CHECKING")
        for imp in imports:
            if imp in future:
                continue
            if isinstance(imp, ast.ImportFrom):
                base = resolve_from(old_mod, False, imp)
                if base == target.rpartition(".")[0] and any(a.name == target.rpartition(".")[2] for a in imp.names):
                    self_alias = next(a.asname or a.name for a in imp.names if a.name == target.rpartition(".")[2])
                    if self_alias in used:
                        used.discard(self_alias)
                if base == target:
                    continue
            line = trimmed_import(imp, used - own)
            if line and isinstance(imp, ast.ImportFrom) and imp.level:
                line = from_line(resolve_from(old_mod, False, imp),
                                 [(a.name, a.asname) for a in imp.names if (a.asname or a.name) in used])
            if line:
                header.append(line)
        sib: dict[str, list[tuple[str, None]]] = {}
        for n in sorted(used - own):
            if n in name_to and name_to[n] != target:
                sib.setdefault(name_to[n], []).append((n, None))
        header += [from_line(m, ns) for m, ns in sorted(sib.items())]
        if typed:
            header.append("\nif TYPE_CHECKING:\n" + "\n".join("    " + ln for t in typed for ln in t.splitlines()))
        code = "\n\n".join(c.rstrip("\n") for _, c in items) + "\n"
        if self_alias:
            code = re.sub(rf"\b{re.escape(self_alias)}\.", "", code)
        dest = root / path_of(target)
        touched.append(path_of(target))
        if not apply:
            print(f"  split {rel} -> {path_of(target)}: {', '.join(n for node in nodes for n in top_names(node))}")
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            existing = dest.read_text()
            existing = add_imports(existing, header) if header else existing
            dest.write_text(existing.rstrip("\n") + "\n\n\n" + code)
        else:
            fut = "from __future__ import annotations\n\n" if future else ""
            if plan.table.get("docstring_to", {}).get(rel) == target and doc is not None and not stay:
                d0 = tree.body[0]
                dest.write_text("".join(lines[: d0.end_lineno]).rstrip("\n") + "\n\n" + fut + "\n".join(header)
                                + ("\n\n\n" if header else "\n") + code)
                if apply:
                    touched.append(path_of(target))
                continue
            names = ", ".join(n for node in nodes for n in top_names(node))
            first = (doc or "").strip().splitlines()[0] if doc else ""
            old = rel.removeprefix("src/objdel/")  # a path, not a dotted name: the old module no longer resolves
            docline = f'"""{target.rpartition(".")[2].capitalize()}: split from the old {old} by phase E ({names}).'
            if len(docline) > 116:
                docline = f'"""Split from the old {old} by phase E (E5 move table: experiments/phase2/tools/e_moves.json).'
            del first
            dest.write_text(docline + '"""\n\n' + fut + "\n".join(header) + ("\n\n\n" if header else "\n") + code)
    if stay:  # the old module keeps the rest and imports the moved names the rest still uses
        cut = {i for items in chunks.values() for node, _ in items
               for i in range(chunk_start(lines, node), node.end_lineno)}
        rest = "".join(ln for i, ln in enumerate(lines) if i not in cut)
        body_used = used_names([n for n in ast.parse(rest).body if not isinstance(n, ast.Import | ast.ImportFrom)])
        moved: dict[str, list[tuple[str, None]]] = {}
        for n in sorted(body_used):
            if n in name_to and name_to[n] != old_mod:
                moved.setdefault(name_to[n], []).append((n, None))
        rest = add_imports(rest, [from_line(m, ns) for m, ns in sorted(moved.items())])
        if apply:
            (root / rel).write_text(re.sub(r"\n{4,}", "\n\n\n", rest))
            touched.append(rel)
        else:
            print(f"  keep in {rel}: {', '.join(stay)}")
    elif apply:
        git(root, "rm", "-q", "-f", rel)
    else:
        print(f"  remove {rel}")
    return touched


# ---------------------------------------------------------------- driver


def git(root: Path, *args: str) -> None:
    subprocess.run(["git", *args], cwd=root, check=True)


def iter_files(root: Path, roots: list[str], globs: list[str]) -> list[Path]:
    out = []
    for r in roots:
        base = root / r
        if base.is_file():
            out.append(base)
            continue
        if not base.exists():
            continue
        for g in globs:
            for p in base.rglob(g):
                parts = p.relative_to(root).parts
                if parts[0] not in SKIP_DIRS and not any(x in ("__pycache__", ".venv", ".git") for x in parts):
                    out.append(p)
    return sorted(set(out))


def strip_init(text: str) -> str:
    """Keep only the docstring of an old package __init__ (re-exports and __all__ go)."""
    tree = ast.parse(text)
    doc = ast.get_docstring(tree, clean=False)
    return f'"""{doc}"""\n' if doc is not None else ""


def wrap_prose(text: str, limit: int = 120, keep: frozenset[str] | set[str] = frozenset()) -> tuple[str, list[int]]:
    """Wrap lines that a longer module name pushed past `limit`. Comment and docstring lines: the tail after the last
    space that fits goes to a new line with the same indent (and `# ` for a comment). Code lines: a trailing comment
    moves to its own line above; else the line breaks after its last comma inside brackets. A line in `keep` (the
    tree's lines before the move) was long already: it stays as it is. Returns the text and the numbers of the long
    lines it left (the hand edit)."""
    import io
    import tokenize

    lines = text.splitlines(keepends=True)
    prose: set[int] = set()  # 0-based line numbers inside a docstring or holding only a comment
    comma: dict[int, int] = {}  # 0-based line -> column after the last comma inside brackets that keeps it short
    trailing: dict[int, int] = {}  # 0-based line -> column of a trailing comment
    depth = 0
    try:
        for tok in tokenize.generate_tokens(io.StringIO(text).readline):
            if tok.type == tokenize.OP and tok.string in "([{":
                depth += 1
            elif tok.type == tokenize.OP and tok.string in ")]}":
                depth -= 1
            elif tok.type == tokenize.OP and tok.string == "," and depth > 0 and tok.end[1] < limit - 1:
                comma[tok.start[0] - 1] = tok.end[1]
            if tok.type == tokenize.COMMENT and lines[tok.start[0] - 1][: tok.start[1]].strip():
                trailing[tok.start[0] - 1] = tok.start[1]
            if tok.type == tokenize.COMMENT and not lines[tok.start[0] - 1][: tok.start[1]].strip():
                prose.add(tok.start[0] - 1)
            elif tok.type == tokenize.STRING and tok.string.lstrip("rRbBuU").startswith(('"""', "'''")):
                prose.update(range(tok.start[0] - 1, tok.end[0]))
    except (tokenize.TokenError, IndentationError):
        return text, []
    out: list[str] = []
    left: list[int] = []
    for i, line in enumerate(lines):
        body = line.rstrip("\n")
        if len(body) <= limit or body in keep:
            out.append(line)
            continue
        cut = body.rfind(" ", 0, limit + 1)
        indent = len(body) - len(body.lstrip())
        head, tail = body[:cut].rstrip(), body[cut + 1:]
        nl = line[len(body):] or "\n"
        if i not in prose and i in trailing and len(body[: trailing[i]].rstrip()) <= limit:
            c = trailing[i]  # the trailing comment goes on its own line above the code
            out += [" " * indent + body[c:] + "\n", body[:c].rstrip() + nl]
            continue
        if i not in prose and i not in trailing and i in comma and len(body) - comma[i] + indent + 4 <= limit:
            c = comma[i]  # break after the last comma inside brackets; the rest is indented one level deeper
            out += [body[:c] + "\n", " " * (indent + 4) + body[c:].lstrip() + nl]
            continue
        if i not in prose or cut <= indent or not head.strip() or tail.strip() in ('"""', "'''"):
            out.append(line)
            left.append(i + 1)
            continue
        is_comment = body.lstrip().startswith("#")
        out += [head + "\n", " " * indent + ("# " if is_comment else "") + tail + nl]
    return "".join(out), left


def run(root: Path, apply: bool, stage: str = STAGES[-1]) -> int:
    table = stage_table(json.loads(TABLE.read_text()), stage)
    plan = build_plan(root, table)
    if plan.errors:
        for e in plan.errors:
            print("STOP:", e)
        return 2
    keep = set(table["gpucheck_keep"])
    roots = table["scan_roots"]
    py_files = iter_files(root, roots, ["*.py"])
    text_files = [p for p in iter_files(root, roots + table["text_extra"], table["text_globs"]) if p.suffix != ".py"]
    text_files += [root / x for x in table["text_extra"] if (root / x).is_file()]
    notes: list[str] = []
    changed: list[str] = []
    before: set[str] = set()  # every line of the scanned Python files before the move (wrap_prose keeps them)
    for p in py_files:
        rel = str(p.relative_to(root))
        if rel in ("experiments/phase2/tools/e_moves.py", "tests/test_e_moves.py"):  # their old names are planted
            continue
        text = p.read_text()
        before.update(text.splitlines())
        new, n1 = rewrite_python(plan, rel, text, rel in keep)
        new, n2 = rewrite_strings(plan, rel, new, rel in keep)
        new = resolve_plain(plan, new)
        notes += n1 + n2
        if new != text:
            changed.append(rel)
            if apply:
                p.write_text(new)
    for p in sorted(set(text_files)):
        rel = str(p.relative_to(root))
        # phase_e goldens record old names on purpose (old_names.json is the E-G16 baseline; module_strings_stale is
        # re-recorded after the moves)
        if rel.startswith(("docs/implementation2.0", "docs/research", "tests/goldens/phase_e/")) or "tools/e_moves" in rel:
            continue
        text = p.read_text()
        new, n2 = rewrite_strings(plan, rel, text, rel in keep)
        notes += n2
        if new != text:
            changed.append(rel)
            if apply:
                p.write_text(new)
    print(f"rewrite: {len(changed)} files" + (":" if not apply else ""))
    if not apply:
        for c in changed:
            print("  edit", c)
    touched = []
    for old, new in plan.moves:
        if not apply:
            print(f"  git mv {old} {new}")
            continue
        (root / new).parent.mkdir(parents=True, exist_ok=True)
        git(root, "mv", old, new)
        touched.append(new)
    for rel, targets in plan.splits.items():
        touched += split_module(plan, rel, targets, apply)
    for init in table["drop_reexports"] + [n for o, n in table["files"].items() if o.endswith("__init__.py")]:
        f = root / init
        if f.exists() and apply:
            s = strip_init(f.read_text())
            if s != f.read_text():
                f.write_text(s)
                touched.append(init)
    for pkg in table["delete_packages"]:
        d = root / pkg
        if not d.exists():
            continue
        left = [p for p in d.rglob("*.py") if p.name != "__init__.py" and "__pycache__" not in p.parts]
        if left:
            print(f"  keep {pkg}: still holds {[str(p.relative_to(root)) for p in left]} (no row in the table)")
            continue
        if not apply:
            print(f"  delete package {pkg}")
            continue
        if (d / "__init__.py").exists():
            git(root, "rm", "-q", "-f", str((d / "__init__.py").relative_to(root)))
        shutil.rmtree(d, ignore_errors=True)
    for path, content in table["packages"].items():
        f = root / path
        if f.exists():
            continue
        if not apply:
            print(f"  create {path}")
            continue
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(content)
        touched.append(path)
    eng = root / table["engine"]["path"]
    if not eng.exists():
        names = table["engine"]["names"]
        if not apply:
            print(f"  create {table['engine']['path']} ({len(names)} names of objdel._gpucheck)")
        else:
            eng.parent.mkdir(parents=True, exist_ok=True)
            body = ",\n".join(f"    {n}" for n in names)
            alls = ",\n".join(f'    "{n}"' for n in names)
            eng.write_text('"""The CUDA geometry engine: the only importer of the `objdel._gpucheck` extension (E3a).\n\n'
                           "Callers import it as `from objdel.geom import engine`. A test checks that `__all__` covers "
                           "every\n`_gpucheck.<name>` used in the repo.\n\"\"\"\n\n"
                           f"from objdel._gpucheck import (  # noqa: TID251\n{body},\n)\n\n__all__ = [\n{alls},\n]\n")
            touched.append(table["engine"]["path"])
    if apply:
        fix = [f for f in set(changed + touched) if f.endswith(".py") and (root / f).exists()]
        if fix:
            subprocess.run([sys.executable, "-m", "ruff", "check", "--fix", "--select", "I001,F401", "-q", *fix],
                           cwd=root, check=False)
    for f in sorted(set(changed + touched)):
        if f.endswith(".py") and (root / f).exists():
            text = (root / f).read_text()
            new, left = wrap_prose(text, keep=before)
            if apply and new != text:
                (root / f).write_text(new)
            notes += [f"{f}:{n}: line longer than 120 after the rename: wrap it by hand" for n in left]
    for n in sorted(set(notes)):
        print("MANUAL:", n)
    for n in plan.errors_soft:
        print("CHECK:", n)
    if not changed and not plan.moves and not plan.splits:
        print("nothing to do: the tree already has the phase E layout")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true", help="print the plan and change nothing")
    ap.add_argument("--root", type=Path, default=Path.cwd(), help="repo root (default: the current folder)")
    ap.add_argument("--stage", choices=STAGES, default=STAGES[-1], help="apply the moves up to this sub-move")
    a = ap.parse_args(argv)
    return run(a.root.resolve(), apply=not a.dry_run, stage=a.stage)


if __name__ == "__main__":
    sys.exit(main())
