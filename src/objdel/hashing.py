"""The two SHA-256 helpers every module shares (no heavy imports, so any module may use them)."""

import hashlib


def file_sha256(path) -> str:
    """SHA-256 hex digest of the file's bytes."""
    with open(path, "rb") as fh:
        return hashlib.file_digest(fh, "sha256").hexdigest()


def pixels_sha256(image) -> str:
    """The image hash recorded in every row and audit: SHA-256 of repr(shape) plus the raw bytes."""
    return hashlib.sha256(repr(image.shape).encode() + image.tobytes()).hexdigest()
