"""The one exception family of objdel.infra: library code raises it, the CLI entry points turn it into a message and
an exit status."""

import re


class InfraError(Exception):
    """A refusal or a failed external call (vastai, ssh, rsync, the Hub, a bundle input); `str(e)` is the message."""


_TOKEN_RE = re.compile(r"hf_[A-Za-z0-9]{20,}|gh[pousr]_[A-Za-z0-9]{20,}|github_pat_\w{20,}"
                       r"|(?<=-p )[A-Za-z0-9+/=_-]{40,}|(?i:(?<=token[ =:]))[A-Za-z0-9+/=_-]{40,}")


def redact(text: str) -> str:
    """text with token-shaped strings (Hugging Face and GitHub tokens, a 40+ character run after `-p ` or `token`)
    replaced by <redacted>: for error and ledger text built from a command line."""
    return _TOKEN_RE.sub("<redacted>", text)
