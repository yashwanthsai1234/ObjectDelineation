"""hf_pull (box) and fetch_hf.sh (Mac): a model revision from the private HF repo onto the box for RL; the token
travels on the ssh session's stdin only, as for push_hf (never on argv, in the environment or on disk)."""

import io
from pathlib import Path
from types import SimpleNamespace

import pytest

from experiments.phase2 import hf_pull as L
from experiments.phase2 import hf_push as H
from experiments.phase2.vast import policy as P
from tests.test_p2_hf_push import REPO, token


def test_hf_pull_downloads_the_revision_with_the_stdin_token(tmp_path, monkeypatch):
    got = {}
    monkeypatch.setattr(L, "snapshot_download", lambda **kw: got.update(kw) or str(tmp_path / "m"))
    monkeypatch.setattr("sys.stdin", io.StringIO("hf_secret\n"))
    assert L.main(["--repo", REPO, "--revision", "r1-qwen3vl2b", "--out", str(tmp_path / "m")]) == 0
    assert got == {"repo_id": REPO, "revision": "r1-qwen3vl2b", "local_dir": str(tmp_path / "m"), "token": "hf_secret"}


def test_hf_pull_refuses_without_a_token(monkeypatch):
    monkeypatch.setattr("sys.stdin", io.StringIO(""))
    with pytest.raises(SystemExit, match="no token on stdin"):
        L.main(["--repo", REPO, "--revision", "r", "--out", "o"])


def test_fetch_hf_streams_the_token_over_stdin_only(monkeypatch):
    sent = []
    monkeypatch.setattr(P, "get_run", lambda rid: {"instance": 7, "ip": "1.2.3.4", "port": 22})
    monkeypatch.setattr(P, "secret", lambda s: {"objdel-hf-write": "hf_secret"}[s])
    monkeypatch.setattr(H, "whoami", lambda tok: token())
    monkeypatch.setattr(P, "remote", lambda run, cmd, stdin=None, what="": sent.append((cmd, stdin, what)))
    P.cmd_fetch_hf(SimpleNamespace(run_id="p1", revision="r1-qwen3vl2b"))
    ((cmd, stdin, what),) = sent
    assert "hf_secret" not in cmd and stdin == "hf_secret\n" and "objdel.secrets" in what
    assert f"-m experiments.phase2.hf_pull --repo {REPO} --revision r1-qwen3vl2b --out /workspace/runs/p1/ckpt" in cmd


def test_fetch_hf_has_its_shim():
    assert (Path(P.__file__).parent / "fetch_hf.sh").resolve().name == "vast.sh"
