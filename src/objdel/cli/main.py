"""The `objdel` command: `objdel <group> <command> [ARGS]` runs objdel.cli.<group>.<command> (a dash in a name stands
for an underscore) as `python -I -m <module> ARGS`. A console script cannot pass -I, so the command always starts a
fresh isolated interpreter: the working directory never lands on sys.path (the CLIs read untrusted data), and the
module runs exactly as its documented `python -I -m` form. `objdel box` also runs the Vast commands of
objdel.infra.vast.

    objdel                                            the groups
    objdel train                                      the commands of a group, one line each
    objdel train rl --spec specs/b1_rows.yaml --stage train
    objdel box offers RTX_5090 --pick top             = python -I -m objdel.infra.vast offers RTX_5090 --pick top
"""

import ast
import importlib
import importlib.util
import os
import pkgutil
import sys

from objdel.core.spec import CLI

VAST = "objdel.infra.vast"  # the box group's commands beside its own modules
HELP = ("-h", "--help")


class UsageError(Exception):
    """No command named: the message is the help to print, `code` the exit status (0 when help was asked for)."""

    def __init__(self, text: str, code: int = 2):
        super().__init__(text)
        self.code = code


def summary(doc: str | None, width: int = 100) -> str:
    """The first paragraph of a docstring on one line, cut at a word to `width` characters."""
    text = " ".join((doc or "").split("\n\n")[0].split())
    return text if len(text) <= width else text[:width - 4].rsplit(" ", 1)[0] + " ..."


def module_doc(module: str) -> str | None:
    """A module's docstring, read from its source without importing it."""
    with open(importlib.util.find_spec(module).origin) as f:
        return ast.get_docstring(ast.parse(f.read()))


def groups() -> dict[str, str]:
    """Group name -> its summary: the packages under objdel.cli."""
    pkg = importlib.import_module(CLI)
    return {m.name: summary(module_doc(f"{CLI}.{m.name}")) for m in pkgutil.iter_modules(pkg.__path__) if m.ispkg}


def commands(group: str) -> dict[str, str]:
    """Command name -> its summary: the modules of the group; the box group adds the Vast commands."""
    pkg = importlib.import_module(f"{CLI}.{group}")
    out = {m.name.replace("_", "-"): summary(module_doc(f"{CLI}.{group}.{m.name}"))
           for m in pkgutil.iter_modules(pkg.__path__)}
    if group == "box":
        vast = importlib.import_module(f"{VAST}.__main__")
        out |= {c.replace("_", "-"): summary(getattr(vast, f"cmd_{c}").__doc__) for c in vast.CLI}
    return dict(sorted(out.items()))


def usage(group: str | None = None) -> str:
    """The help text: the groups, or the commands of one group."""
    head, rows = ("objdel <group> <command> [ARGS]", groups()) if group is None else (
        f"objdel {group} <command> [ARGS]   (ARGS --help: the command's own flags)", commands(group))
    return "\n".join([f"usage: {head}", "", *(f"  {name:<14} {text}" for name, text in rows.items())])


def resolve(argv: list[str]) -> list[str]:
    """The arguments after `python -I` for `objdel ARGV`: ["-m", module, *ARGS]; UsageError when ARGV names no
    command."""
    if not argv or argv[0] in HELP:
        raise UsageError(usage(), 0 if argv else 2)
    group, rest = argv[0], argv[1:]
    if group not in groups():
        raise UsageError(f"objdel: no group {group!r}\n\n{usage()}")
    if not rest or rest[0] in HELP:
        raise UsageError(usage(group), 0 if rest else 2)
    if rest[0] not in commands(group):
        raise UsageError(f"objdel: no command {group} {rest[0]!r}\n\n{usage(group)}")
    name, args = rest[0].replace("-", "_"), rest[1:]
    module = f"{CLI}.{group}.{name}"
    if group == "box" and importlib.util.find_spec(module) is None:
        return ["-m", VAST, name, *args]
    return ["-m", module, *args]


def main(argv: list[str] | None = None) -> int:
    """Resolve the command and replace this process with `python -I -m <module> ARGS`; print help otherwise."""
    try:
        cmd = resolve(sys.argv[1:] if argv is None else list(argv))
    except UsageError as e:
        print(e, file=sys.stderr if e.code else sys.stdout)
        return e.code
    os.execv(sys.executable, [sys.executable, "-I", *cmd])


if __name__ == "__main__":
    sys.exit(main())
