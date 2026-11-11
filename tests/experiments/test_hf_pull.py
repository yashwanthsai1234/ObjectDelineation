"""hf_pull (box) and fetch_hf.sh (Mac): a model revision from the private HF repo onto the box for RL; the token
travels on the ssh session's stdin only, as for push_hf (never on argv, in the environment or on disk)."""

import hashlib
import io
import json
import shlex
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

hub = pytest.importorskip("huggingface_hub")
from objdel.cli.hub import pull as L  # noqa: E402
from objdel.infra import hf as H  # noqa: E402
from objdel.infra.errors import InfraError  # noqa: E402
from objdel.infra.vast import __main__ as P  # noqa: E402
from tests.experiments.test_p2_hf_push import REPO, token  # noqa: E402


def test_hf_pull_downloads_the_revision_with_the_stdin_token(tmp_path, monkeypatch):
    got = {}
    monkeypatch.setattr(hub, "snapshot_download", lambda **kw: got.update(kw) or str(tmp_path / "m"))
    monkeypatch.setattr("sys.stdin", io.StringIO("hf_secret\n"))
    assert L.main(["--repo", REPO, "--revision", "r1-qwen3vl2b", "--out", str(tmp_path / "m")]) == 0
    assert got == {"repo_id": REPO, "revision": "r1-qwen3vl2b", "local_dir": str(tmp_path / "m"), "token": "hf_secret",
                   "allow_patterns": None}
    assert json.loads((tmp_path / "m/hf_revision.json").read_text()) == {"repo": REPO, "revision": "r1-qwen3vl2b"}


def test_hf_pull_writes_no_revision_marker_for_a_rows_pull(tmp_path, monkeypatch):
    monkeypatch.setattr(hub, "snapshot_download", lambda **kw: str(tmp_path))
    monkeypatch.setattr("sys.stdin", io.StringIO("hf_secret\n"))
    assert L.main(["--repo", REPO, "--revision", "a" * 40, "--out", str(tmp_path), "--allow", "rows/x/**"]) == 0
    assert not (tmp_path / "hf_revision.json").exists()


def test_hf_pull_refuses_without_a_token(monkeypatch):
    monkeypatch.setattr("sys.stdin", io.StringIO(""))
    with pytest.raises(SystemExit, match="no token on stdin"):
        L.main(["--repo", REPO, "--revision", "r", "--out", "o"])


def box(monkeypatch) -> list:
    sent = []
    monkeypatch.setattr(P, "get_run", lambda rid: {"instance": 7, "ip": "1.2.3.4", "port": 22})
    monkeypatch.setattr(P, "secret", lambda s: {"objdel-hf-write": "hf_secret"}[s])
    monkeypatch.setattr(H, "whoami", lambda tok: token())
    monkeypatch.setattr(P, "remote", lambda run, cmd, stdin=None, what="": sent.append((cmd, stdin, what)))
    return sent


def fetch(revision=None, rows=None, name="ckpt"):
    P.cmd_fetch_hf(SimpleNamespace(run_id="p1", revision=revision, rows=rows, name=name))


def test_fetch_hf_streams_the_token_over_stdin_only(monkeypatch):
    sent = box(monkeypatch)
    fetch("r1-qwen3vl2b")
    ((cmd, stdin, what),) = sent
    assert "hf_secret" not in cmd and stdin == "hf_secret\n" and "objdel.infra.secrets" in what
    assert f"-m objdel.cli.hub.pull --repo {REPO} --revision r1-qwen3vl2b --out /workspace/runs/p1/ckpt" in cmd


def test_fetch_hf_has_its_shim():
    assert (P.repo_root() / "experiments/phase2/vast/fetch_hf.sh").resolve().name == "vast.sh"


# --- B1.5: rows root and checkpoint folder names (I40, I40b) ------------------------------------------------------

ROWS_REV = "d2c8edb8e15534363ebbced8f839cecace49c495"  # B1.4 step 5 upload (rows/lidc-train-b1) and its manifest pin
ROWS_SHA = "a0f2c1d212b8c685e2b85aaac1aa9fa3de9c004ab310af13412d5aa630d4a6b9"


def test_hf_pull_allow_keeps_repo_paths_under_local_dir(tmp_path, monkeypatch):
    """I40 / I40b (a): --allow passes allow_patterns; snapshot_download keeps repo paths under local_dir, so the repo
    folder rows/lidc-train-b1/ lands at {data}/rows/lidc-train-b1, the root train_rl --train reads."""
    def fake(repo_id, revision, local_dir, token, allow_patterns):
        assert allow_patterns == "rows/lidc-train-b1/**" and revision == ROWS_REV
        d = Path(local_dir) / "rows/lidc-train-b1"
        (d / "lidc").mkdir(parents=True)
        (d / "manifest.json").write_text("{}")
        (d / "lidc/train.jsonl").write_text("")
        return local_dir

    monkeypatch.setattr(hub, "snapshot_download", fake)
    monkeypatch.setattr("sys.stdin", io.StringIO("hf_secret\n"))
    data = tmp_path / "data/lidc"
    assert L.main(["--repo", REPO, "--revision", ROWS_REV, "--out", str(data), "--allow", "rows/lidc-train-b1/**"]) == 0
    root = data / "rows/lidc-train-b1"
    assert (root / "manifest.json").is_file() and (root / "lidc/train.jsonl").is_file()


def test_fetch_rows_downloads_the_root_and_checks_it_on_the_box(monkeypatch):
    sent = box(monkeypatch)
    fetch(rows=[ROWS_REV, "lidc-train-b1"])
    ((cmd, stdin, _),) = sent
    assert "hf_secret" not in cmd and stdin == "hf_secret\n"
    assert (f"-m objdel.cli.hub.pull --repo {REPO} --revision {ROWS_REV} --out /workspace/data/lidc "
            "--allow 'rows/lidc-train-b1/**' && ") in cmd
    untar, check = (shlex.split(c) for c in cmd.split(" && ")[-2:])
    assert untar[:4] == ["/opt/venv/bin/python", "-I", "-c", P.UNTAR]
    assert untar[4:] == ["/workspace/data/lidc/rows/lidc-train-b1/root.tar", "/workspace/data/lidc/rows/lidc-train-b1"]
    assert check[:3] == ["/opt/venv/bin/python", "-I", "-c"] and check[3] == P.ROWS_CHECK
    assert check[4:] == ["/workspace/data/lidc/rows/lidc-train-b1", ROWS_SHA]  # the SHA the plan pins for that revision


def rows_root(tmp: Path, files_count=2) -> tuple[Path, str]:
    d = tmp / "rows/lidc-train-b1"
    (d / "lidc").mkdir(parents=True)
    (d / "lidc/train.jsonl").write_text("{}")
    (d / "images.png").write_bytes(b"x")
    (d / "manifest.json").write_text(json.dumps({"schedule": {"files_count": files_count}}))
    return d, hashlib.sha256((d / "manifest.json").read_bytes()).hexdigest()


@pytest.mark.parametrize("case", ["ok", "sha", "count", "missing"])
def test_box_check_exits_non_zero_on_a_wrong_root(tmp_path, case):
    d, sha = rows_root(tmp_path, files_count=3 if case == "count" else 2)
    if case == "missing":
        (d / "manifest.json").unlink()
    r = subprocess.run([sys.executable, "-I", "-c", P.ROWS_CHECK, str(d), "0" * 64 if case == "sha" else sha],
                       capture_output=True, text=True)
    assert (r.returncode == 0) == (case == "ok"), r.stderr


def test_box_unpacks_the_pushed_tar_and_the_check_passes(tmp_path):
    """hf.rows_tar on the Mac, UNTAR then ROWS_CHECK on the box: the same files, the tar removed, exit 0."""
    d, sha = rows_root(tmp_path)
    box = tmp_path / "box"
    box.mkdir()
    with (box / H.ROWS_TAR).open("wb") as fh:
        H.rows_tar(d, fh)
    (box / "manifest.json").write_bytes((d / "manifest.json").read_bytes())
    for code, args in ((P.UNTAR, [str(box / H.ROWS_TAR), str(box)]), (P.ROWS_CHECK, [str(box), sha])):
        r = subprocess.run([sys.executable, "-I", "-c", code, *args], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
    files = lambda root: sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file())  # noqa: E731
    assert files(box) == files(d)


def test_fetch_rows_refuses_a_revision_no_plan_pins(monkeypatch):
    sent = box(monkeypatch)
    with pytest.raises(InfraError, match="manifest SHA-256"):
        fetch(rows=["f" * 40, "lidc-train-b1"])
    assert sent == []


def test_fetch_as_names_the_checkpoint_folder(monkeypatch):
    """I40b (b): --as b1 puts revision b1-rows at {rd}/b1, the default stays {rd}/ckpt; / or .. is refused."""
    sent = box(monkeypatch)
    fetch("b1-rows", name="b1")
    fetch("rl-full-3ep")
    assert sent[0][0].endswith("--revision b1-rows --out /workspace/runs/p1/b1")
    assert sent[1][0].endswith("--revision rl-full-3ep --out /workspace/runs/p1/ckpt")
    for bad in ("a/b", "..", "../x"):
        with pytest.raises(InfraError):
            fetch("b1-rows", name=bad)
        with pytest.raises(InfraError):
            fetch(rows=[ROWS_REV, bad])
    assert len(sent) == 2


def test_fetch_hf_takes_a_revision_or_rows(monkeypatch):
    sent = box(monkeypatch)
    with pytest.raises(InfraError):
        fetch()  # neither
    with pytest.raises(InfraError):
        fetch("x", rows=[ROWS_REV, "lidc-train-b1"])  # both
    assert sent == []
    got = []
    monkeypatch.setattr(P, "cmd_fetch_hf", got.append)
    P.main(["fetch_hf", "p1", "--rows", ROWS_REV, "lidc-train-b1"])
    P.main(["fetch_hf", "p1", "b1-rows", "--as", "b1"])
    assert (got[0].rows, got[0].revision, got[1].revision, got[1].name) == ([ROWS_REV, "lidc-train-b1"], None,
                                                                             "b1-rows", "b1")
