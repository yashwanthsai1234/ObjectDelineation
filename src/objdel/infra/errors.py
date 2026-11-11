"""The one exception family of objdel.infra: library code raises it, the CLI entry points turn it into a message and
an exit status."""


class InfraError(Exception):
    """A refusal or a failed external call (vastai, ssh, rsync, the Hub, a bundle input); `str(e)` is the message."""
