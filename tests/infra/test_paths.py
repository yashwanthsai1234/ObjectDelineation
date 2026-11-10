"""U34 / E-G15 (E8.1): one repo root for every moved file; the Vast state stays in experiments/phase2/vast."""



from objdel.core.provenance import repo_root
from objdel.infra import secrets
from objdel.infra.vast import ledger, remote


def test_repo_root_holds_the_source_tree_from_any_cwd(tmp_path, monkeypatch):
    monkeypatch.delenv("OBJDEL_REPO", raising=False)
    monkeypatch.chdir(tmp_path)
    root = repo_root()
    assert (root / "pyproject.toml").is_file() and (root / "src/objdel").is_dir()


def test_objdel_repo_wins(tmp_path, monkeypatch):
    (tmp_path / "pyproject.toml").touch()  # E2: an OBJDEL_REPO without one is refused (tests/core/test_provenance.py)
    monkeypatch.setenv("OBJDEL_REPO", str(tmp_path))
    assert repo_root() == tmp_path


def test_secrets_default_to_the_repo_root(monkeypatch):
    monkeypatch.delenv("OBJDEL_SECRETS", raising=False)
    assert secrets.path() == repo_root() / ".secrets"


def test_vast_state_stays_in_experiments_phase2_vast(monkeypatch, tmp_path):
    """After a plain move home() would be the package folder: spend() would read $0 and active() {} (unsafe)."""
    monkeypatch.delenv("VAST_HOME", raising=False)
    assert ledger.home() == repo_root() / "experiments/phase2/vast"
    assert "image" in ledger.load(".vast-allowlist")
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    assert ledger.home() == tmp_path


def test_onstart_reads_the_watchdog_from_the_code_tree_not_the_vast_home(monkeypatch, tmp_path):
    """d050dcc: VAST_HOME holds state only; an empty VAST_HOME still launches the code tree's watchdog.sh."""
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    assert remote.watchdog_script() == repo_root() / "experiments/phase2/vast/watchdog.sh"
    assert "watchdog.sh" in remote.onstart_cmd("r1", 600)
