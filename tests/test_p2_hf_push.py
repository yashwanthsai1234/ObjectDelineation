"""hf_push (box) and push_hf.sh (Mac): the selected branch to the private HF repo, nothing else, never RIGA."""

import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from experiments.phase2 import hf_push as HP
from experiments.phase2.vast import policy as P
from objdel.infra import hf as H

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
        return HP.main(argv)
    with pytest.raises(SystemExit) as e:
        HP.main(argv)
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
    assert "hf_secret" not in cmd and stdin == "hf_secret\n" and "objdel.infra.secrets" in what
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
    assert HP.main(["--run", str(run), "--repo", REPO, "--revision", "r1-lfm", "--view", str(local)]) == 0
    (local / "manifest.json").write_text(json.dumps({"dataset": "lidc", "view_hash": "other"}))
    monkeypatch.setattr("sys.stdin", io.StringIO("hf_secret\n"))
    with pytest.raises(SystemExit, match="view_hash"):
        HP.main(["--run", str(run), "--repo", REPO, "--revision", "r1-lfm", "--view", str(local)])


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


# --- hf_push --rl: a train_rl output (B1.5, I41) ------------------------------------------------------------------


def rl_dir(tmp: Path, final=True, adapter=True) -> Path:
    out = tmp / "runs" / "v1" / "b1-rows"
    out.mkdir(parents=True)
    if final:
        (out / "final").mkdir()
        (out / "final" / "model.safetensors").write_bytes(b"merged")
        (out / "final" / "config.json").write_text("{}")
    if adapter:
        (out / "final_adapter").mkdir()
        (out / "final_adapter" / "adapter_model.safetensors").write_bytes(b"lora")
        (out / "final_adapter" / "adapter_config.json").write_text("{}")
    (out / "trainer" / "checkpoint-50").mkdir(parents=True)
    (out / "trainer" / "checkpoint-50" / "optimizer.pt").write_bytes(b"opt")
    (out / "rl.json").write_text(json.dumps({"model": "qwen3-vl-2b", "steps": 260, "ckpt": "/workspace/runs/v1/ckpt",
                                             "rows": {"revision": "a" * 40, "ids_sha256": "1" * 64}}))
    return out


def push_rl(out, monkeypatch, ok=True):
    monkeypatch.setattr("sys.stdin", io.StringIO("hf_secret\n"))
    argv = ["--rl", str(out), "--repo", REPO, "--revision", "b1-rows"]
    if ok:
        return HP.main(argv)
    with pytest.raises(SystemExit) as e:
        HP.main(argv)
    return str(e.value)


def test_push_rl_uploads_the_merged_model_and_the_adapter(tmp_path, monkeypatch, fake):
    """final/ at the revision's root (the folder fetch_hf.sh and final_eval --ckpt read), final_adapter/ beside it,
    the base licence and an LIDC card; no checkpoint or optimizer state."""
    assert push_rl(rl_dir(tmp_path), monkeypatch) == 0
    ((_, repo, rev, files, _),) = [c for c in fake.calls if c[0] == "create_commit"]
    assert repo == REPO and rev == "b1-rows" and ("create_branch", REPO, "b1-rows", True) in fake.calls
    assert set(files) == {"model.safetensors", "config.json", "final_adapter/adapter_model.safetensors",
                          "final_adapter/adapter_config.json", "LICENSE", "README.md"}
    assert files["model.safetensors"] == b"merged" and files["final_adapter/adapter_model.safetensors"] == b"lora"
    card = files["README.md"].decode()
    for want in ("Qwen/Qwen3-VL-2B-Instruct", "base_model: Qwen/Qwen3-VL-2B-Instruct", "GRPO", "260 steps",
                 "final_adapter", "a" * 40, "LIDC-IDRI", "10.7937/K9/TCIA.2015.LO9QL9SX"):
        assert want in card, want


@pytest.mark.parametrize("case, want", [({"final": False}, "final"), ({"adapter": False}, "final_adapter")])
def test_push_rl_refuses_an_incomplete_run(tmp_path, monkeypatch, fake, case, want):
    assert want in push_rl(rl_dir(tmp_path, **case), monkeypatch, ok=False)
    assert not [c for c in fake.calls if c[0] in ("create_branch", "create_commit")]


def test_push_rl_refuses_a_public_repo(tmp_path, monkeypatch, fake):
    fake.private = False
    assert "not private" in push_rl(rl_dir(tmp_path), monkeypatch, ok=False)
    assert not [c for c in fake.calls if c[0] in ("create_branch", "create_commit")]


def test_push_hf_of_the_b1_spec_pushes_the_train_rl_output(monkeypatch):
    """push_hf.sh <run> specs/b1_rows.yaml: the last rl stage's out is a train_rl output (revision b1-rows)."""
    sent = _mac(monkeypatch, token())
    P.cmd_push_hf(SimpleNamespace(run_id="v1", config="specs/b1_rows.yaml"))
    ((cmd, stdin, _),) = sent
    assert "-m experiments.phase2.hf_push --rl /workspace/runs/v1/b1-rows" in cmd and "hf_secret" not in cmd
    assert f"--repo {REPO} --revision b1-rows" in cmd and stdin == "hf_secret\n"


# ---------------------------------------------------------------- DOC11: --card-only
HEADLINE = {"run": "r5", "model_ref": "svenk029/objdel-sft@rl-full-3ep", "seed": 11, "K": 8, "ci_unit": "item",
            "test": {"source": "runs/r5/eval/eval_test/final_eval.json", "view": "lidc-test-shifted", "items": 2200,
                     "metrics": {k: {"mean": 0.5 + i / 100, "ci95": [0.4, 0.6]}
                                 for i, k in enumerate(("iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote",
                                                        "legal"))}},  # fmt: skip
            "val": {"source": "runs/r5/eval/eval_val/final_eval.json", "view": "lidc-val-shifted", "items": 1084,
                    "metrics": {k: {"mean": 0.3, "ci95": [0.2, 0.4]}
                                for k in ("iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote",
                                          "legal")}}}  # fmt: skip
OLD_CARD = "---\nlicense: apache-2.0\n---\n\n# objdel RL r5\n\nBase: Qwen3-VL-2B.\n\n## Data\n\n" + H.LIDC + "\n"


def test_card_update_appends_sections_and_keeps_every_other_line():
    new = H.card_update(OLD_CARD, HEADLINE)
    assert set(OLD_CARD.splitlines()) <= set(new.splitlines())
    assert "## Evaluation" in new and "## Limitations" in new and "0.500" in new and "2200" in new
    assert "shifted views" in new and "single-run" in new and "RIGA" in new


def test_card_update_is_idempotent_and_replaces_an_old_evaluation():
    once = H.card_update(OLD_CARD, HEADLINE)
    assert H.card_update(once, HEADLINE) == once
    other = HEADLINE | {"seed": 12}
    again = H.card_update(once, other)
    assert again.count("## Evaluation") == 1 and "seed 12" in again and "seed 11" not in again


def card_only(monkeypatch, tmp_path, who=None, private=True):
    api = FakeApi(private)
    api.upload_file = lambda **kw: api.calls.append(("upload_file", kw["path_in_repo"], kw["revision"],
                                                      kw["path_or_fileobj"]))  # fmt: skip
    (tmp_path / "h.json").write_text(json.dumps(HEADLINE))
    (tmp_path / "README.md").write_text(OLD_CARD)
    monkeypatch.setattr(H, "api", lambda tok: api)
    monkeypatch.setattr(H, "whoami", lambda tok: who or token())
    monkeypatch.setattr("huggingface_hub.hf_hub_download", lambda *a, **k: str(tmp_path / "README.md"))
    monkeypatch.setattr("sys.stdin", io.StringIO("hf_secret\n"))
    argv = ["--card-only", "--repo", REPO, "--revision", "rl-full-3ep", "--headline", str(tmp_path / "h.json")]
    return api, argv


def test_card_only_uploads_exactly_the_readme(tmp_path, monkeypatch):
    pytest.importorskip("trl")
    api, argv = card_only(monkeypatch, tmp_path)
    assert HP.main(argv) == 0
    (up,) = [c for c in api.calls if c[0] == "upload_file"]
    assert up[1:3] == ("README.md", "rl-full-3ep") and b"## Evaluation" in up[3]


def test_card_only_with_a_bad_token_or_public_repo_uploads_nothing(tmp_path, monkeypatch):
    pytest.importorskip("trl")
    api, argv = card_only(monkeypatch, tmp_path, who=token(role="write"))
    with pytest.raises(SystemExit):
        HP.main(argv)
    api, argv = card_only(monkeypatch, tmp_path, private=False)
    with pytest.raises(SystemExit):
        HP.main(argv)
    assert not [c for c in api.calls if c[0] == "upload_file"]
