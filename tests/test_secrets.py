"""objdel.secrets: <repo>/.secrets first (name=value, mode 600 only), then the macOS keychain."""
from types import SimpleNamespace

import pytest

from objdel import secrets

ROOT = secrets.ROOT


def test_secrets_file_first_then_the_keychain(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(secrets.subprocess, "run", lambda cmd, **k: calls.append(cmd) or
                        SimpleNamespace(returncode=0, stdout="from-keychain\n"))  # fmt: skip
    f = tmp_path / ".secrets"
    monkeypatch.setenv("OBJDEL_SECRETS", str(f))
    f.write_text("# HF\nobjdel-hf-write = hf_abc \nobjdel-ghcr-read=\n")
    f.chmod(0o600)
    assert secrets.get("objdel-hf-write") == "hf_abc" and not calls
    assert secrets.get("objdel-ghcr-read") == "from-keychain" and calls[-1][:2] == ["security", "find-generic-password"]


def test_a_secrets_file_readable_by_others_is_refused(tmp_path, monkeypatch):
    f = tmp_path / ".secrets"
    monkeypatch.setenv("OBJDEL_SECRETS", str(f))
    f.write_text("objdel-hf-write=hf_abc\n")
    f.chmod(0o644)
    with pytest.raises(PermissionError, match="chmod 600"):
        secrets.get("objdel-hf-write")


def test_missing_everywhere_is_none(monkeypatch):
    monkeypatch.setattr(secrets.subprocess, "run", lambda cmd, **k: SimpleNamespace(returncode=44, stdout=""))
    assert secrets.get("objdel-nothing") is None


def test_the_secrets_file_is_git_ignored_at_the_repo_root():
    assert secrets.path().name == ".secrets" and "\n/.secrets\n" in (ROOT / ".gitignore").read_text()
