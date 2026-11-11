"""The refusal the policy and SFT libraries raise; the CLI edge turns it into a message and an exit code."""

import sys


class Refusal(Exception):
    """A run refused for a bad flag, view, checkpoint or GPU state. `code` is the process exit code."""

    def __init__(self, message: str, code: int = 1):
        super().__init__(message)
        self.code = code

    def system_exit(self) -> SystemExit:
        """What a CLI raises: code 1 carries the message (Python prints it); another code prints it to stderr first."""
        if self.code == 1:
            return SystemExit(str(self))
        print(self, file=sys.stderr)
        return SystemExit(self.code)

