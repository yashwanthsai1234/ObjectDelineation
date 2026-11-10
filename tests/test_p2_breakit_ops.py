"""Break-it BK4 (Phase 2 Stage 1): Vast operations and the CI `[dev]`-only job. Real defects are strict xfails
(`BK4-<n>` in experiments/phase2/impl/BK4.md); the rest are guards that pass today. No real Vast call: the fake
executables of tests/test_vast.py, extended with two knobs (a transient vastai failure, a running reply without
ports)."""

import ast
import json
import os
import re
import subprocess
import sys
import time
from fnmatch import fnmatchcase
from pathlib import Path, PurePosixPath

import numpy as np
import pytest

from objdel.infra.vast import ledger as LG
from objdel.infra.vast import offers as OF
from objdel.infra.vast import remote as RM
from objdel.sft import train as sft_train
from tests import test_vast as V
from tests.test_p2_scripts import stage

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = sorted(c for c in (ROOT / "experiments" / "phase2" / "configs").glob("*.yaml")
                 if "\nsteps:" not in c.read_text())  # r0.yaml is a step list for run.sh, not a train_sft config
ANCHOR = '    save = lambda: (d / "state.json").write_text(json.dumps(st))\n'
HOOK = """    k = " ".join(argv[:2])
    if st.get("fail", {}).get(k, 0) > 0:
        st["fail"][k] -= 1
        save()
        sys.exit("fake vastai: transient HTTP 502")
    if k == "create instance" and st.get("create_lost"):  # the box is created, the reply is lost
        iid = st["next_id"]; st["next_id"] += 1
        st["instances"][str(iid)] = {"id": iid, "label": argv[argv.index("--label") + 1], "dph_total": 0.6,
                                     "start_date": time.time()}
        save(); print(st["create_lost"]["out"]); sys.exit(st["create_lost"]["code"])
    if k == "show instance" and st.get("no_ports"):
        print(json.dumps(dict(st["instances"][argv[2]], actual_status="running")))
        sys.exit(0)
"""


@pytest.fixture
def E(tmp_path):
    assert ANCHOR in V.FAKE
    e = V.Env(tmp_path)
    fake = V.FAKE.replace(ANCHOR, ANCHOR + HOOK)
    (tmp_path / "bin" / "_fake.py").write_text(f"#!{sys.executable}\nSECRETS = {V.SECRETS!r}\n" + fake)
    return e


def knob(E, **kw):
    f = E.fake / "state.json"
    f.write_text(json.dumps(dict(json.loads(f.read_text()), **kw)))


# --- CI: the [dev]-only job (no torch / transformers / trl / datasets / pydicom, no data/) ---------------------------

BLOCKED = ("torch", "transformers", "trl", "datasets", "peft", "accelerate", "pydicom", "pylidc")
CI_RUN = (
    "import sys\n"
    f"for m in {BLOCKED!r}:\n"
    "    sys.modules[m] = None  # import raises ModuleNotFoundError, find_spec returns None: the [dev]-only env\n"
    "import pytest\n"
    "sys.exit(pytest.main(sys.argv[1:]))\n"
)


@pytest.mark.parametrize("node", ["tests/test_p2_bridge.py", "tests/test_p2_views.py::test_build_eval_views_help"])
def test_ci_dev_only_job_skips_instead_of_failing(node):
    env = dict(os.environ, OBJDEL_FIXTURES="/nonexistent", OBJDEL_STORE="/nonexistent")
    cmd = [sys.executable, "-c", CI_RUN, "-q", "-p", "no:cacheprovider", "-m", "not slow", node]
    p = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, timeout=600)
    assert p.returncode == 0, p.stdout[-2000:]


# --- launch: every failure after `create instance` must destroy (G8) -------------------------------------------------


def test_launch_destroys_when_polling_fails(E):
    E.pick()
    knob(E, fail={"show instance": 3})
    E.sh("launch.sh", "111", "r1", ok=False)
    assert E.instances() == {}


def test_launch_survives_one_status_poll_hiccup(E):
    """One failed `show instance` destroyed a good box: the poll retries (tries=3)."""
    E.pick()
    knob(E, fail={"show instance": 1})
    E.sh("launch.sh", "111", "r1")
    assert list(E.instances()) == ["4242"] and E.active()["r1"]["port"] == 40022


@pytest.mark.parametrize("out,code", [('{"success": true}', 0), ("", 1), ("not json", 0)])
def test_launch_destroys_a_box_whose_create_reply_is_lost(E, out, code):
    """An error, an interrupt or a reply without `new_contract` from `create` came before the cleanup: the box was
    created and nothing destroyed it. The cleanup finds it by its label."""
    E.state(instances={"9": {"id": 9, "label": "other", "dph_total": 0.3, "start_date": time.time()}})
    E.pick()
    knob(E, create_lost={"out": out, "code": code})
    E.sh("launch.sh", "111", "r1", "--allow-parallel", ok=False)
    assert list(E.instances()) == ["9"] and E.vast("destroy instance") == [["destroy", "instance", "4242", "-y"]]
    assert E.ledger()[-1]["event"] == "destroy" and E.active() == {}


def test_launch_destroys_when_running_reply_lacks_ports(E):
    E.pick()
    knob(E, no_ports=True)
    E.sh("launch.sh", "111", "r1", ok=None)
    assert E.instances() == {}


def test_launch_bad_status_and_timeout_destroy_even_without_status(E):
    E.state(statuses=[None])
    E.pick()
    p = E.sh("launch.sh", "111", "r1", ok=False, VAST_LAUNCH_TIMEOUT_S="0.2")
    assert "timeout" in p.stderr and E.instances() == {} and E.active() == {}


# --- run layout on the box: what run.sh writes is what pull.sh pulls -------------------------------------------------


def config_value(path: Path, key: str) -> str | None:
    m = re.search(rf"^{key}:\s*(\S+)", path.read_text(), re.M)
    return m.group(1) if m else None


def test_config_outputs_land_where_pull_reads(monkeypatch):
    # run.sh: `cd /workspace/src && ... python -m <module>`; pull.sh: rsync /workspace/runs/<RUN_ID>/ only.
    # train_sft resolves output_dir under run_dir = $OBJDEL_RUN_DIR (run.sh's /workspace/runs/<RUN_ID>), else runs/.

    assert all(config_value(c, "output_dir") and config_value(c, "run_dir") is None for c in CONFIGS)
    monkeypatch.setenv("OBJDEL_RUN_DIR", "/workspace/runs/v1")
    for c in CONFIGS:
        out = sft_train.resolve(stage("m", "v", config_value(c, "output_dir"), vision=False)).output_dir
        assert PurePosixPath("/workspace/runs/v1") in PurePosixPath("/workspace/src", out).parents, c.name
    monkeypatch.delenv("OBJDEL_RUN_DIR")
    assert sft_train.resolve(stage("m", "v", "r1-x", vision=False)).output_dir == str(Path("runs/r1-x"))


def test_config_data_root_is_the_box_bundle_dir(tmp_path, monkeypatch):
    # sync_up extracts the bundle to /workspace/data/lidc: the default data_root there, data/ on the Mac.

    assert not [c.name for c in CONFIGS if config_value(c, "data_root")]
    monkeypatch.delenv("OBJDEL_DATA_ROOT", raising=False)
    monkeypatch.setattr(sft_train, "BOX_DATA", tmp_path / "missing")
    assert sft_train.resolve(stage("m", "v", "o")).data_root == "data"
    monkeypatch.setattr(sft_train, "BOX_DATA", tmp_path)
    assert sft_train.resolve(stage("m", "v", "o")).data_root == str(tmp_path)
    monkeypatch.setenv("OBJDEL_DATA_ROOT", "/elsewhere")
    assert sft_train.resolve(stage("m", "v", "o")).data_root == "/elsewhere"


# --- cost guard (G9, G10) --------------------------------------------------------------------------------------------


def test_sweep_reconciles_a_run_whose_instance_vanished(E):
    E.launched(cap="4")
    act = E.active()
    act["r1"]["started"] -= 30 * 3600  # watchdog self-destroyed it at its cap long ago (G5); the Mac never heard
    act["r1"]["deadline"] -= 30 * 3600
    (E.home / "active_runs.json").write_text(json.dumps(act))
    E.state(instances={})
    E.sh("sweep_orphans.sh")
    assert E.active() == {} and any(r["run_id"] == "r1" and r["cost_usd"] for r in E.ledger())


def test_parallel_launches_cannot_jointly_pass_the_hard_stop(E):
    E.pick()  # $0.55/h
    for r in ("r1", "r2", "r3"):  # 3 x 0.55 x 12 h = $19.80 worst case
        E.sh("launch.sh", "111", r, "--cap-hours", "12", "--allow-parallel", ok=None)
    assert len(E.vast("create instance")) <= 2


def test_launch_projected_past_the_alarm_needs_approval(E):
    V.spent(E, 9.0)
    E.pick()
    E.sh("launch.sh", "111", "r1", "--cap-hours", "4", ok=None)  # 9 + 0.55 x 4 = 11.2 > 10
    assert not E.vast("create instance")


def test_launch_refuses_an_active_run_id(E):
    E.launched()
    E.sh("launch.sh", "111", "r1", "--allow-parallel", ok=None)
    assert len(E.vast("create instance")) == 1


@pytest.mark.parametrize("run_id", ["r1;touch${IFS}/tmp/x", "../../etc", "r1.2", "r" * 41])
def test_launch_validates_run_id(E, run_id):
    E.pick()
    E.sh("launch.sh", "111", run_id, ok=None)
    assert not E.vast("create instance")


def test_spend_counts_closed_rows_only_once_and_ignores_blank_costs(E, monkeypatch):
    monkeypatch.setenv("VAST_HOME", str(E.home))
    LG.ledger(event="offer", offer_id=1)
    LG.ledger(event="teardown", run_id="a", cost_usd=1.25)
    LG.ledger(event="orphan", run_id="b", cost_usd=0.5)
    assert LG.spend() == pytest.approx(1.75)


# --- watchdog (G5) ---------------------------------------------------------------------------------------------------


def test_watchdog_keeps_box_when_outputs_are_newer_than_the_pull(E):
    run = E.tmp / "box" / "runs" / "r1"
    (run / "checkpoint-100").mkdir(parents=True)
    ok, ckpt, now = run / "PULLED_OK", run / "checkpoint-100" / "model.safetensors", time.time()
    ok.touch()
    ckpt.write_bytes(b"w")
    (run / "HEARTBEAT").write_text("1\n")  # a job ran (without one the never-started rule applies)
    os.utime(ok, (now + 60,) * 2)  # the Mac pulled during the idle period ...
    os.utime(ckpt, (now + 120,) * 2)  # ... and a CPU stage wrote a checkpoint after it
    p, curls = V.watchdog(E, ok, max_ticks=3, WATCHDOG_IDLE_S="0")
    assert curls == [] and p.returncode == 3


# --- H-GOV (D16), offers (§4.2) --------------------------------------------------------------------------------------


@pytest.mark.parametrize("granted", [True, False, "true"])
@pytest.mark.parametrize("evidence", ["", "written permission 2026-11-01"])
@pytest.mark.parametrize("hosting_type", [0, 1, None])
def test_governance_matrix(granted, evidence, hosting_type):
    o = V.offer(hosting_type=hosting_type)
    allow = {"riga_permission": {"granted": granted, "evidence": evidence}}
    ok = not OF.governance("riga", o, allow)
    assert ok == (granted is True and bool(evidence) and hosting_type == 1)
    assert OF.governance("lidc", o, {}) == [] and OF.governance("other", o, allow)


@pytest.mark.parametrize("gpu", sorted(OF.CARDS))
def test_card_boundaries(gpu):
    c = OF.CARDS[gpu]
    lo, hi = c["ram"]
    base = V.offer(gpu=gpu.replace("_", " "), gpu_ram=hi * 1024, dph_total=c["ceiling"], hosting_type=1)
    assert OF.reasons(base, gpu, "lidc") == [] and OF.reasons(base, gpu, "riga") == []
    assert OF.reasons(dict(base, dph_total=c["ceiling"] + 0.001), gpu, "lidc") == ["price fails the §4.2 filter"]
    assert OF.reasons(dict(base, gpu_ram=hi * 1024 + 1), gpu, "lidc") == ["gpu_ram fails the §4.2 filter"]
    if lo:
        assert OF.reasons(dict(base, gpu_ram=lo * 1024 - 1), gpu, "lidc") == ["gpu_ram fails the §4.2 filter"]
    need = max(32, c.get("cpu_ram", 0))
    assert OF.reasons(dict(base, cpu_ram=need * 1000 - 1), gpu, "lidc") == ["cpu_ram fails the §4.2 filter"]  # MB, V0
    assert OF.reasons(dict(base, hosting_type=0), gpu, "riga") == ["datacenter fails the §4.2 filter"]
    q = OF.query(gpu, "lidc")
    assert f"dph<={c['ceiling']}" in q and f"gpu_ram<={hi}" in q and "datacenter" not in q
    assert "datacenter=True" in OF.query(gpu, "riga") and "cpu_ram>=48" in OF.query(gpu, "riga")


@pytest.mark.parametrize("bad", [{"reliability2": 0.5}, {"rentable": False}, {"inet_down": 10},
                                 {"direct_port_count": 0}, {"disk_space": 10}])  # fmt: skip
def test_offer_recheck_covers_the_whole_query(bad):
    assert OF.reasons(V.offer(**bad), "RTX_4090", "lidc")


# --- secrets (D16, G9), destroy-never-stop (G4) ----------------------------------------------------------------------


def test_riga_identity_only_on_the_decrypt_stdin_for_the_whole_session(E):
    V.riga_ready(E)
    outs = [
        E.sh("sync_up.sh", "r1", "--bundle", "riga", "--bundle-dir", E.bundles, DRY_RUN="0"),
        E.sh("run.sh", "r1", "experiments.phase2.k3_harness", "--out", "/workspace/runs/r1/k3.json", DRY_RUN="0"),
        E.sh("pull.sh", DRY_RUN="0"),
        E.sh("teardown.sh", "r1", DRY_RUN="0"),
        E.sh("sweep_orphans.sh", DRY_RUN="0"),
    ]
    carriers = [c for c in E.calls() if V.IDENTITY in json.dumps(c)]
    (dec,) = carriers  # the fake `security` prints the identity but records only its argv
    assert dec["tool"] == "ssh" and "age -d -i /dev/stdin" in dec["argv"][-1]
    assert V.IDENTITY not in json.dumps(dec["argv"]) and dec["stdin"].strip() == V.IDENTITY
    assert not any(c["env_has_identity"] for c in E.calls())
    assert V.IDENTITY not in "".join(o.stdout + o.stderr for o in outs)
    state = V.files_text(E.home, E.runs, E.bundles)  # not bin/ (the fakes)
    assert state and V.IDENTITY not in state


def test_destroy_never_stop(E):
    E.launched()
    E.sh("sync_up.sh", "r1", "--bundle-dir", E.bundles)
    E.sh("teardown.sh", "r1")
    E.launched()
    knob(E, instances=dict(E.instances(), **{"9": {"id": 9, "label": "x", "dph_total": 0.3, "start_date": 0}}))
    V.spent(E, 17.5)
    E.sh("sweep_orphans.sh")
    verbs = {tuple(a[:2]) for a in (c["argv"] for c in E.calls("vastai"))}
    assert not [v for v in verbs if "stop" in v[0] or "stop" in v[1]] and ("destroy", "instance") in verbs
    assert E.instances() == {}
    src = (V.VAST / "policy.py").read_text() + (V.VAST / "watchdog.sh").read_text()
    assert not re.search(r'vast\(\s*"stop"|/stop|"stop"', src) and "-X DELETE" in src


# --- transfer: data/raw, RIGA-derived files, symlinks ----------------------------------------------------------------


def test_code_snapshot_has_no_riga_derived_results():
    leaks = [f for f in RM.code_files() if "riga" in f.lower() and f.endswith((".csv", ".json", ".jsonl"))]
    leaks += [f for f in RM.code_files() if f.startswith("docs/research/records/phase0/results/vlm_")]
    assert not leaks, leaks[:5]


RIGA_BLOCK = """/docs/research/records/phase0/results/medsam_riga.csv
/docs/research/records/phase0/results/riga/
/docs/research/records/phase0/results/riga_outlines/
/docs/research/records/phase0/results/vlm_box/
/docs/research/records/phase0/results/vlm_polygon/
/archive/phase1_lab/t_reward/results/h5_items.csv
/archive/phase1_lab/t_data/riga/out_maxrg_check/
/archive/phase1_lab/t_data/riga/maxrg_vs_prime.json
/archive/phase1_lab/t_data/splits/out/riga_images.csv
/archive/phase1_lab/t_reward/results/advantage_groups.json
/archive/phase1_lab/t_reward/results/iter3_r5_raw.json
/archive/phase1_lab/t_stack/results/t[24]_*.json""".splitlines()
RIGA_TABLES = {"archive/phase1_lab/t_data/riga/inventory.json", "archive/phase1_lab/t_data/riga/repo_parity.json"}


def riga_block_hits(paths, block=RIGA_BLOCK):
    """Paths that a line of the RIGA block of .gitignore covers (file line, folder line, or the t[24] glob)."""
    return [p for p in paths if any(p.startswith(b[1:]) if b.endswith("/") else fnmatchcase(p, b[1:]) for b in block)]


def riga_named(paths):
    return [p for p in paths if "riga" in p.lower() and p.lower().endswith((".csv", ".json", ".jsonl"))
            and p not in RIGA_TABLES]  # fmt: skip


def committable(root):
    out = subprocess.run(["git", "ls-files", "-z", "-co", "--exclude-standard"], cwd=root, capture_output=True,
                         text=True, check=True).stdout  # what `git add -A` stages
    return [f for f in out.split("\0") if f]


def test_no_riga_derived_file_is_committable(tmp_path):
    """CLEAN5 step 0: the RIGA block of .gitignore follows the moved trees (RIGA terms: results stay local)."""
    assert (ROOT / ".gitignore").read_text().count("\n".join(RIGA_BLOCK)) == 1  # the constant is the real block
    files = committable(ROOT)
    assert files and riga_block_hits(files) == [] and riga_named(files) == []
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    new = ["docs/research/records/phase0/results/riga/masks.csv", "archive/phase1_lab/t_stack/results/t2_x.json",
           "docs/x/riga/a.csv"]  # fmt: skip
    for n in new:
        (tmp_path / n).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / n).write_text("x")
    gi = tmp_path / ".gitignore"
    gi.write_text((ROOT / ".gitignore").read_text())
    got = committable(tmp_path)
    assert riga_named(got) == ["docs/x/riga/a.csv"] and riga_block_hits(got) == []
    gi.write_text("")
    got = committable(tmp_path)
    assert set(riga_block_hits(got)) == set(new[:2])


def golden_leaks(folder):
    """Arrays in a golden that could hold a mask: more than one dimension, or bool (polycore.npz is the exception)."""
    leaks = []
    for f in sorted(folder.glob("*.npz")):
        if f.name != "polycore.npz":
            z = np.load(f, allow_pickle=False)
            leaks += [f"{f.name}:{k}" for k in z.files if z[k].ndim > 1 or z[k].dtype == bool]
    return leaks


def test_goldens_hold_no_masks(tmp_path):
    """CLEAN4 licence rule: RIGA terms forbid copying masks, and tests/ ships in every bundle."""
    assert (ROOT / "tests" / "goldens" / "polycore.npz").exists() and golden_leaks(ROOT / "tests" / "goldens") == []
    np.savez(tmp_path / "a.npz", m=np.zeros((2, 2), np.uint8), v=np.zeros(3, bool), ok=np.zeros(3, np.uint8))
    assert golden_leaks(tmp_path) == ["a.npz:m", "a.npz:v"]
    np.savez(tmp_path / "polycore.npz", m=np.zeros((2, 2)))
    assert golden_leaks(tmp_path) == ["a.npz:m", "a.npz:v"]  # polycore.npz alone is exempt


def test_bundle_refuses_a_symlinked_view_directory(tmp_path):
    from objdel.infra import bundle as B

    d = V.data_tree(tmp_path)
    (d / "views" / "lidc-extra").symlink_to(d / "views" / "riga-train", target_is_directory=True)
    with pytest.raises(SystemExit):
        B.build("lidc", d, tmp_path / "out")


# --- image (D15) -----------------------------------------------------------------------------------------------------


def test_image_check_peaks_and_threshold():
    text = (V.DOCKER / "image_check.sh").read_text()
    peak = ast.literal_eval(re.search(r"^PEAK = (\{.*\})$", text, re.M).group(1))
    lookup = lambda name: next((v for k, v in peak.items() if k in name), None)  # noqa: E731 (as image_check.sh)
    assert lookup("NVIDIA GeForce RTX 4090") == 165.2 and lookup("NVIDIA GeForce RTX 3090") == 71
    assert lookup("NVIDIA GeForce RTX 5090") == 209.5 and lookup("NVIDIA A100-SXM4-80GB") == 312
    assert lookup("NVIDIA RTX PRO 6000 Blackwell Workstation Edition") == 500 and lookup("Tesla T4") is None
    assert "tf >= 0.6 * peak" in text and "8192" in text
