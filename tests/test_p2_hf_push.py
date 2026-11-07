"""hf_push (box) and push_hf.sh (Mac): the selected branch to the private HF repo, nothing else, never RIGA."""

import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from experiments.phase2 import hf_push as H
from experiments.phase2.vast import policy as P

REPO = "svenk029/objdel-sft"


def token(role="fineGrained", scoped=None, global_=()):
    scoped = [{"entity": {"type": "model", "name": REPO}, "permissions": ["repo.content.read", "repo.write"]}] \
        if scoped is None else scoped  # fmt: skip
    return {"auth": {"type": "access_token", "accessToken": {
        "role": role, "fineGrained": {"global": list(global_), "scoped": scoped}}}}  # fmt: skip


def test_a_token_scoped_to_the_one_repo_passes():
    assert H.token_problems(token(), REPO) == []


@pytest.mark.parametrize("who, why", [
    (token(role="write"), "role write"),
    (token(scoped=[{"entity": {"type": "user", "name": "svenk029"},
                    "permissions": ["repo.write", "inference.serverless.write", "job.write"]}]), "user svenk029"),
    (token(scoped=[{"entity": {"type": "model", "name": REPO},
                    "permissions": ["repo.write", "inference.endpoints.write"]}]), "inference.endpoints.write"),
    (token(scoped=[{"entity": {"type": "model", "name": "svenk029/other"}, "permissions": ["repo.write"]}]),
     "svenk029/other"),
    (token(global_=["discussion.write"]), "discussion.write"),
    (token(scoped=[]), "no write access"),
])  # fmt: skip
def test_an_over_broad_token_is_refused(who, why):
    """The keychain token reaches a rented box: it may write this one repo and nothing else (no inference, jobs or
    billing: the owner's HF account has no card, and a leaked broad token could act on every repo)."""
    assert any(why in p for p in H.token_problems(who, REPO))


class FakeApi:
    def __init__(self, private=True):
        self.private, self.calls = private, []

    def repo_info(self, repo, repo_type):
        self.calls.append(("repo_info", repo, repo_type))
        return SimpleNamespace(private=self.private)

    def create_branch(self, repo, branch, repo_type, exist_ok):
        self.calls.append(("create_branch", repo, branch, exist_ok))

    def create_commit(self, repo, operations, commit_message, revision, repo_type):
        files = {op.path_in_repo: Path(op.path_or_fileobj).read_bytes() if isinstance(op.path_or_fileobj, (str, Path))
                 else op.path_or_fileobj for op in operations}  # fmt: skip
        self.calls.append(("create_commit", repo, revision, files, commit_message))


def run_dir(tmp: Path, dataset="lidc", selected="e1") -> Path:
    view = tmp / "views" / f"{dataset}-bundle"
    view.mkdir(parents=True)
    (view / "manifest.json").write_text(json.dumps({"dataset": dataset, "view_hash": "vh"}))
    run = tmp / "runs" / "v1" / "r1-lfm"
    for b in ("e0.5", "e1"):
        (run / "branches" / b).mkdir(parents=True)
        (run / "branches" / b / "model.safetensors").write_bytes(b"w-" + b.encode())
        (run / "branches" / b / "config.json").write_text("{}")
    (run / "select").mkdir()
    (run / "select" / "selected.json").write_text(json.dumps({
        "selected": selected, "gates_branch": selected or "e1", "gates_pass": True, "gates_failed": [],
        "sampler": "hf", "k": 8, "items": 256}))  # fmt: skip
    (run / "manifest.json").write_text(json.dumps({
        "config": {"model": "lfm2.5-vl-450m", "seed": 1, "lora": False, "branches": [0.5, 1.0, 2.0],
                   "view": str(view)},
        "view": {"view_hash": "vh"}, "code": {"objdel_and_script_sha256": "c0de"}, "status": "done"}))  # fmt: skip
    return run


@pytest.fixture
def fake(monkeypatch):
    pytest.importorskip("trl")  # hf_push.main reads objdel.train.configs (the training stack)
    api = FakeApi()
    monkeypatch.setattr(H, "api", lambda tok: api.calls.append(("token", tok)) or api)
    monkeypatch.setattr(H, "base_card", lambda base: (b"LFM OPEN LICENSE", {"license": "other",
                                                                             "license_name": "lfm1.0"}))
    return api


def push(run, monkeypatch, ok=True):
    monkeypatch.setattr("sys.stdin", io.StringIO("hf_secret\n"))
    argv = ["--run", str(run), "--repo", REPO, "--revision", "r1-lfm"]
    if ok:
        return H.main(argv)
    with pytest.raises(SystemExit) as e:
        H.main(argv)
    return str(e.value)


def test_push_uploads_the_selected_branch_with_licence_and_card_in_one_commit(tmp_path, monkeypatch, fake):
    assert push(run_dir(tmp_path), monkeypatch) == 0
    (commit,) = [c for c in fake.calls if c[0] == "create_commit"]
    _, repo, rev, files, _ = commit
    assert repo == REPO and rev == "r1-lfm" and ("token", "hf_secret") in fake.calls
    assert ("create_branch", REPO, "r1-lfm", True) in fake.calls
    assert files["model.safetensors"] == b"w-e1" and set(files) == {"model.safetensors", "config.json", "README.md",
                                                                    "LICENSE"}  # fmt: skip
    assert files["LICENSE"] == b"LFM OPEN LICENSE"
    card = files["README.md"].decode()
    for want in ("LiquidAI/LFM2.5-VL-450M", "license_name: lfm1.0", "base_model: LiquidAI/LFM2.5-VL-450M", "e1",
                 "full fine-tune", "LIDC-IDRI", "10.7937/K9/TCIA.2015.LO9QL9SX", "c0de"):
        assert want in card, want


@pytest.mark.parametrize("case, want", [
    ({"dataset": "riga"}, "RIGA"),
    ({"selected": None}, "no branch selected"),
])  # fmt: skip
def test_push_refuses_riga_and_an_unselected_run(tmp_path, monkeypatch, fake, case, want):
    assert want in push(run_dir(tmp_path, **case), monkeypatch, ok=False)
    assert not [c for c in fake.calls if c[0] in ("create_branch", "create_commit")]


def test_push_refuses_a_public_repo(tmp_path, monkeypatch, fake):
    fake.private = False
    assert "not private" in push(run_dir(tmp_path), monkeypatch, ok=False)
    assert not [c for c in fake.calls if c[0] in ("create_branch", "create_commit")]


# --- Mac side: push_hf.sh RUN_ID CONFIG ---------------------------------------------------------------------------


def _mac(monkeypatch, who):
    sent = []
    monkeypatch.setattr(P, "get_run", lambda rid: {"instance": 7, "ip": "1.2.3.4", "port": 22})
    monkeypatch.setattr(P, "secret", lambda s: {"objdel-hf-write": "hf_secret"}[s])
    monkeypatch.setattr(H, "whoami", lambda tok: who)
    monkeypatch.setattr(P, "remote", lambda run, cmd, stdin=None, what="": sent.append((cmd, stdin, what)))
    return sent


def test_push_hf_streams_the_token_over_stdin_only(monkeypatch):
    sent = _mac(monkeypatch, token())
    P.cmd_push_hf(SimpleNamespace(run_id="v1", config="experiments/phase2/configs/r1_lfm.yaml"))
    ((cmd, stdin, what),) = sent
    assert "hf_secret" not in cmd and stdin == "hf_secret\n" and "objdel.secrets" in what
    assert "-m experiments.phase2.hf_push --run /workspace/runs/v1/r1-lfm" in cmd
    assert f"--repo {REPO} --revision r1-lfm" in cmd


def test_push_hf_refuses_an_over_broad_token_before_touching_the_box(monkeypatch):
    sent = _mac(monkeypatch, token(role="write"))
    with pytest.raises(SystemExit, match="objdel-hf-write.*role write"):
        P.cmd_push_hf(SimpleNamespace(run_id="v1", config="experiments/phase2/configs/r1_lfm.yaml"))
    assert sent == []


def test_push_hf_has_its_shim():
    assert (Path(P.__file__).parent / "push_hf.sh").resolve().name == "vast.sh"


def test_push_from_the_mac_takes_a_local_copy_of_the_view_with_the_same_hash(tmp_path, monkeypatch, fake):
    """The run manifest names the box's view path; on the Mac --view points at the local copy, accepted only when its
    view_hash is the one the run trained on (and its dataset is LIDC)."""
    run = run_dir(tmp_path)
    man = json.loads((run / "manifest.json").read_text())
    local = Path(man["config"]["view"])
    man["config"]["view"] = "/workspace/data/lidc/views/lidc-bundle"  # the box path
    (run / "manifest.json").write_text(json.dumps(man))
    monkeypatch.setattr("sys.stdin", io.StringIO("hf_secret\n"))
    assert H.main(["--run", str(run), "--repo", REPO, "--revision", "r1-lfm", "--view", str(local)]) == 0
    (local / "manifest.json").write_text(json.dumps({"dataset": "lidc", "view_hash": "other"}))
    monkeypatch.setattr("sys.stdin", io.StringIO("hf_secret\n"))
    with pytest.raises(SystemExit, match="view_hash"):
        H.main(["--run", str(run), "--repo", REPO, "--revision", "r1-lfm", "--view", str(local)])


@pytest.mark.parametrize("license", ["apache-2.0", "other"])
def test_base_card_without_a_license_file_uses_apache_only_when_the_card_declares_it(monkeypatch, license):
    """Qwen/Qwen3-VL-2B-Instruct (2026-11-05) has no LICENSE file; its card says apache-2.0."""
    huggingface_hub = pytest.importorskip("huggingface_hub")
    from huggingface_hub.errors import EntryNotFoundError

    def missing(*a, **k):
        raise EntryNotFoundError("no LICENSE")

    monkeypatch.setattr(huggingface_hub, "hf_hub_download", missing)
    card = SimpleNamespace(data=SimpleNamespace(to_dict=lambda: {"license": license}))
    monkeypatch.setattr(huggingface_hub.ModelCard, "load", lambda base, token=None: card)
    if license == "other":
        with pytest.raises(EntryNotFoundError):
            H.base_card("x/y")
        return
    text, meta = H.base_card("Qwen/Qwen3-VL-2B-Instruct")
    assert meta == {"license": "apache-2.0"} and b"Apache License" in text and b"Version 2.0, January 2004" in text
    assert hashlib.sha256(text).hexdigest().startswith("cfc7749b")
