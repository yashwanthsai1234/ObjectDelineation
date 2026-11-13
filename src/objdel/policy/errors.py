"""The refusal the policy and SFT libraries raise; the CLI edge turns it into a message and an exit code."""

import functools
import sys


class Refusal(Exception):
    """A run refused for a bad flag, view, checkpoint or GPU state. `code` is the process exit code."""

    def __init__(self, message: str, code: int = 1):
        super().__init__(message)
        self.code = code

    @classmethod
    def usage(cls, prog: str, msg: str) -> "Refusal":
        """A usage error as argparse words it ("<prog>: error: <msg>"), exit code 2."""
        return cls(f"{prog}: error: {msg}", 2)

    def system_exit(self) -> SystemExit:
        """What a CLI raises: code 1 carries the message (Python prints it); another code prints it to stderr first."""
        if self.code == 1:
            return SystemExit(str(self))
        print(self, file=sys.stderr)
        return SystemExit(self.code)



def exits_on_refusal(fn):
    """fn with a Refusal turned into the process exit (the CLI edge: Refusal.system_exit, no traceback)."""
    @functools.wraps(fn)
    def run(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except Refusal as e:
            raise e.system_exit() from None

    return run
