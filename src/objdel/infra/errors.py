"""The one exception family of objdel.infra: library code raises it, the CLI entry points turn it into a message and
an exit status."""

import functools
import sys


class InfraError(Exception):
    """A refusal or a failed external call (vastai, ssh, rsync, the Hub, a bundle input); `str(e)` is the message."""


def exits(prefix: str):
    """A decorator: fn with an InfraError turned into `sys.exit("<prefix>: <message>")` (the CLI edge)."""
    def deco(fn):
        @functools.wraps(fn)
        def run(*args, **kwargs):
            try:
                return fn(*args, **kwargs)
            except InfraError as e:
                sys.exit(f"{prefix}: {e}")

        return run

    return deco
