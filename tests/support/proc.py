"""Run Python in a fresh isolated interpreter (`python -I`): the data and the cwd are untrusted (repo rule)."""
import os
import subprocess
import sys


def run_isolated(*args, code=None, cwd=None, env=None, timeout=None, check=False) -> subprocess.CompletedProcess:
    """`python -I <args>` (or `python -I -c <code>` when `code` is given), stdout and stderr captured as text.
    `env` is merged over os.environ, never a replacement."""
    argv = [sys.executable, "-I", *(["-c", code] if code is not None else []), *map(str, args)]
    return subprocess.run(argv, cwd=cwd, env=None if env is None else os.environ | env, capture_output=True, text=True,
                          timeout=timeout, check=check)


def run_module(*args, env=None, timeout=900) -> subprocess.CompletedProcess:
    """`python -I -m <args>` from the repo root, output captured."""
    from tests.support.paths import ROOT

    return run_isolated("-m", *args, cwd=ROOT, env=env, timeout=timeout)


class InlinePool:
    """A multiprocessing.Pool stand-in that runs every job in this process, in order: coverage and monkeypatches
    reach the worker function (a spawned child sees neither)."""

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def imap(self, fn, jobs, **_):  # chunksize: one process, no chunks
        return map(fn, jobs)

    def map(self, fn, jobs, **_):
        return list(map(fn, jobs))


class _InlineContext:
    Pool = InlinePool


def inline_context(*_):  # the start method ("spawn") does not apply
    """Drop-in for multiprocessing.get_context: monkeypatch.setattr(module, "get_context", inline_context)."""
    return _InlineContext()
