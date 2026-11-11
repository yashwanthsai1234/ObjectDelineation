"""M7: Vast operations (PLAN §4, D16, G1-G10). No real Vast call: every script runs against fake `vastai`, `ssh`,
`rsync`, `security`, `curl` and `nvidia-smi` executables on PATH that record their argv, stdin and environment."""

import base64
import csv
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tarfile
import time
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest

from objdel.core import datasets
from objdel.core.io import git_sha
from objdel.infra import bundle as B
from objdel.infra import preflight as PF
from objdel.infra.errors import InfraError
from objdel.infra.vast import __main__ as P
from objdel.infra.vast import chain as CH
from objdel.infra.vast import ledger as LG
from objdel.infra.vast import offers as OF
from objdel.infra.vast import remote as RM

ROOT = Path(__file__).resolve().parents[2]
VAST = ROOT / "experiments" / "phase2" / "vast"
DOCKER = ROOT / "experiments" / "phase2" / "docker"
SNAP = ROOT / "tests" / "fixtures" / "vast_snapshot.json"
IDENTITY = "AGE-SECRET-KEY-1FAKEIDENTITYFORTESTSQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ"
API_KEY = "fakevastscopedkey0123456789abcdef"
PAT = "ghp_fakeGHCRreadPAT0123456789"
SECRETS = {"objdel-vast-scoped": API_KEY, "objdel-ghcr-read": PAT, "objdel-riga-age": IDENTITY}
SPENDING_SCRIPTS = ("launch.sh", "teardown.sh", "sweep_orphans.sh")
SCRIPTS = ["offers.sh", "launch.sh", "sync_up.sh", "run.sh", "pull.sh", "teardown.sh", "sweep_orphans.sh"]

FAKE = r"""
import json, os, sys, time
from pathlib import Path
d = Path(os.environ["FAKE_DIR"]); tool = Path(sys.argv[0]).name; argv = sys.argv[1:]
stdin = "" if sys.stdin is None or sys.stdin.isatty() else sys.stdin.read()
rec = {"tool": tool, "argv": argv, "stdin": stdin, "env_has_identity": SECRETS["objdel-riga-age"] in str(os.environ),
       "api_key": os.environ.get("VAST_API_KEY")}
with open(d / "calls.jsonl", "a") as f:
    f.write(json.dumps(rec) + "\n")
if tool == "security":
    print(SECRETS[argv[argv.index("-s") + 1]])
elif tool == "nvidia-smi":
    print(os.environ.get("FAKE_GPU_UTIL", "0"))
    if os.environ.get("FAKE_TOUCH"):
        Path(os.environ["FAKE_TOUCH"]).touch()
elif tool == "curl":
    print('{"success": true}')
elif tool == "vastai":
    st = json.loads((d / "state.json").read_text())
    save = lambda: (d / "state.json").write_text(json.dumps(st))
    if argv[:2] == ["search", "offers"]:
        print(json.dumps(st["offers"]))
    elif argv[:2] == ["create", "instance"]:
        iid = st["next_id"]; st["next_id"] += 1
        st["instances"][str(iid)] = {"id": iid, "label": argv[argv.index("--label") + 1], "dph_total": 0.6,
                                     "start_date": time.time() - 3600}
        save(); print(json.dumps({"success": True, "new_contract": iid}))
    elif argv[:2] == ["show", "instance"]:
        s = st["statuses"].pop(0) if len(st["statuses"]) > 1 else st["statuses"][0]; save()
        print(json.dumps(dict(st["instances"][argv[2]], actual_status=s, public_ipaddr="203.0.113.7",
                              ports={"22/tcp": [{"HostPort": "40022"}]})))
    elif argv[:2] == ["show", "instances"]:
        print(json.dumps(list(st["instances"].values())))
    elif argv[:2] == ["destroy", "instance"]:
        st["instances"].pop(argv[2], None); save(); print("destroying instance", argv[2])
    elif argv[:2] == ["show", "invoices-v1"]:
        print(json.dumps(st.get("charges", {"results": [], "count": 0, "next_token": None})))
    else:
        sys.exit(f"fake vastai: unknown {argv}")
"""


def offer(oid=111, gpu="RTX 4090", **kw):
    """An offer as Vast's search returns it, admissible for RTX_4090 / LIDC (V0: `verification` is a string)."""
    return dict(id=oid, ask_contract_id=oid, gpu_name=gpu, num_gpus=1, gpu_ram=24564, cpu_ram=64000,
                cpu_cores_effective=16, cuda_max_good=13.2, dph_total=0.55, reliability2=0.995, disk_space=200,
                inet_down=900, direct_port_count=10, verification="verified", rentable=True, hosting_type=0, host_id=7,
                machine_id=70, inet_up_cost=0.001, inet_down_cost=0.001, geolocation="Texas, US") | kw  # fmt: skip


def read_json(f: Path) -> dict:
    return json.loads(f.read_text()) if f.exists() else {}


def files_text(*dirs: Path) -> str:
    return "".join(p.read_text(errors="replace") for d in dirs for p in d.rglob("*") if p.is_file())


def after(argv: list[str], flag: str) -> str:
    return argv[argv.index(flag) + 1]


def PF_config(tmp: Path) -> Path:
    """A run spec over a tiny data root (tmp/data: a bundle view manifest and the exam views' manifests): enough to
    hash."""
    data = tmp / "data"
    for d in ("views/v", "checker_eval", *(f"eval_views/{v}" for v in datasets.LIDC.exams)):
        (data / d).mkdir(parents=True)
        (data / d / "manifest.json").write_text(json.dumps({"view_hash": "h", "dir": d}))
    cfg = tmp / "pf_run.yaml"
    cfg.write_text("schema_version: 1\nname: pf_run\ncompute: {gpu: RTX_4090, cap_hours: 24.0}\nbundle_view: views/v\n"
                   "stages:\n- {kind: cmd, name: box, argv: [bash, -c, 'true']}\n")
    return cfg


class Env:
    def __init__(self, tmp: Path):
        self.tmp, self.fake, self.home, self.runs = tmp, tmp / "fake", tmp / "home", tmp / "runs"
        self.bundles = tmp / "bundles"
        bin_ = tmp / "bin"
        for p in (bin_, self.fake, self.home, self.runs, self.bundles):
            p.mkdir()
        src = bin_ / "_fake.py"
        src.write_text(f"#!{sys.executable}\nSECRETS = {SECRETS!r}\n" + FAKE)
        src.chmod(0o755)
        for t in ("vastai", "ssh", "rsync", "security", "curl", "nvidia-smi", "osascript"):
            (bin_ / t).symlink_to(src)
        self.state(offers=[offer(), offer(112, dph_total=0.70), offer(113, hosting_type=1, dph_total=0.9)])
        self.allow(image="ghcr.io/owner/objdel-vast:p2-abc@sha256:" + "0" * 64, ghcr_user="owner",
                   riga_permission={"granted": False, "evidence": ""})
        for f in ("bundle-lidc.tar", "bundle-lidc.tar.sha256"):
            (self.bundles / f).write_text("x")
        env = {k: v for k, v in os.environ.items() if not k.startswith(("VAST", "DRY_RUN"))}
        self.env = dict(env, PATH=f"{bin_}{os.pathsep}{env['PATH']}", PYTHON=sys.executable, FAKE_DIR=str(self.fake),
                        VAST_HOME=str(self.home), VAST_RUNS=str(self.runs), DRY_RUN="1", VAST_POLL_S="0.01",
                        VAST_LAUNCH_TIMEOUT_S="3", OBJDEL_RUN_DIR=str(tmp / "pf"),
                        OBJDEL_DATA_ROOT=str(tmp / "data"))
        (self.home / "ledger.csv").write_text(",".join(LG.FIELDS) + "\n")  # launch refuses without a ledger
        (self.home / "watchdog.sh").symlink_to(VAST / "watchdog.sh")  # onstart reads it from the Vast home
        self.config = PF_config(tmp)
        self.stamp()

    def stamp(self, ok=True, **hashes):
        """The preflight stamp launch.sh checks (the gate): passing and current, or as given."""
        path = PF.repo_root() / self.env["OBJDEL_RUN_DIR"] / "preflight" / f"{self.config.stem}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        with mock.patch.dict(os.environ, OBJDEL_DATA_ROOT=self.env["OBJDEL_DATA_ROOT"]):
            now = PF.hashes(self.config)
        stamp = {"config": str(self.config.resolve()), "pass": ok, "checks": {"x": {"pass": ok}},
                 "hashes": now | hashes}
        path.write_text(json.dumps(stamp))

    def state(self, statuses=("loading", "running"), instances=None, offers=None, charges=None):
        old = read_json(self.fake / "state.json")
        st = {"next_id": 4242, "offers": old.get("offers", []) if offers is None else offers,
              "statuses": list(statuses), "instances": old.get("instances", {}) if instances is None else instances}
        if c := charges or old.get("charges"):
            st["charges"] = c
        (self.fake / "state.json").write_text(json.dumps(st))

    def allow(self, **kw):
        """Update .vast-allowlist (None removes a key)."""
        new = read_json(self.home / ".vast-allowlist") | kw
        (self.home / ".vast-allowlist").write_text(json.dumps({k: v for k, v in new.items() if v is not None}))

    def sh(self, script, *args, ok=True, **env):
        if script == "launch.sh" and "--config" not in args:
            args = (*args, "--config", self.config)
        if script in SPENDING_SCRIPTS:  # DRY_RUN=1 never rents or destroys: these run against the fake vastai
            env = {"DRY_RUN": "0", **env}
        p = subprocess.run(["bash", str(VAST / script), *map(str, args)], env=dict(self.env, **env),
                           capture_output=True, text=True, timeout=60)  # fmt: skip
        assert ok is None or (p.returncode == 0) == ok, p.stdout + p.stderr
        return p

    def calls(self, tool=None):
        f = self.fake / "calls.jsonl"
        rows = [json.loads(x) for x in f.read_text().splitlines()] if f.exists() else []
        return [r for r in rows if tool in (None, r["tool"])]

    def vast(self, verb):
        return [c["argv"] for c in self.calls("vastai") if c["argv"][:2] == verb.split()]

    def instances(self):
        return read_json(self.fake / "state.json")["instances"]

    def active(self):
        return json.loads((self.home / "active_runs.json").read_text())

    def ledger(self):
        f = self.home / "ledger.csv"
        return list(csv.DictReader(f.open())) if f.exists() else []

    def pick(self, oid="111", *flags, ok=True):
        return self.sh("offers.sh", "RTX_4090", *flags, "--pick", oid, ok=ok)

    def launched(self, cap="2"):
        self.pick()
        self.sh("launch.sh", "111", "r1", "--cap-hours", cap)

    def tmux(self, *run_args):
        """The tmux line `run.sh r1 ARGS` prints under DRY_RUN."""
        return next(ln for ln in self.sh("run.sh", "r1", *run_args).stdout.splitlines() if "tmux" in ln)


@pytest.fixture
def E(tmp_path):
    return Env(tmp_path)


def all_text(E, *outs):
    """Everything a secret could leak into: state files, ledger, recorded argv/env, script output."""
    # the GHCR PAT is sent to Vast only as `create instance --login` (V10); everything else must be clean
    argv = [[a if i == 0 or c["argv"][i - 1] != "--login" else "<login>" for i, a in enumerate(c["argv"])]
            for c in E.calls()]  # fmt: skip
    env = [c["env_has_identity"] for c in E.calls()]
    return files_text(E.home) + json.dumps(argv + env) + "".join(o.stdout + o.stderr for o in outs)


# --- offer filter and ceilings (§4.2, §4.8) -------------------------------------------------------------------------


def test_query_per_card():
    q = OF.query("RTX_4090")
    terms = ("num_gpus=1 gpu_name=RTX_4090 gpu_ram<=26 cpu_cores_effective>=12 cpu_ram>=32 cuda_vers>=13.0"
             " verified=true reliability>0.98 direct_port_count>=2 dph<=1.0")  # fmt: skip
    assert set(terms.split()) <= set(q.split()) and "datacenter=True" not in q
    assert {"cpu_ram>=32", "dph<=1.5"} <= set(OF.query("RTX_5090").split())
    assert "gpu_ram>=79" in OF.query("A100_SXM4").split()
    assert "cpu_ram>=48" in OF.query("RTX_PRO_6000_WS").split()
    assert "disk_space>=150" in OF.query("RTX_PRO_6000_WS", disk=150).split()
    with pytest.raises(ValueError):
        OF.query("RTX_2080")


def test_ceilings_are_1p5x_snapshot_median():
    snap = json.loads(SNAP.read_text())["offers"]
    for card in ("RTX_3090", "RTX_4090", "RTX_5090", "RTX_PRO_6000_WS", "RTX_PRO_6000_S"):
        med = snap[f"{card.replace('_', ' ')}|on-demand|any"]["median"]
        assert abs(OF.CARDS[card]["ceiling"] - 1.5 * med) <= 0.08, (card, med)
    want = {"RTX_4090": 1.00, "RTX_5090": 1.50, "A100_SXM4": 1.80, "RTX_3090": 0.50}
    assert {c: OF.CARDS[c]["ceiling"] for c in want} == want


@pytest.mark.parametrize("bad,why", [
    ({"dph_total": 1.01}, "price"), ({"gpu_ram": 49140}, "gpu_ram"), ({"cuda_max_good": 12.8}, "cuda"),
    ({"cpu_cores_effective": 8}, "cores"), ({"cpu_ram": 16000}, "cpu_ram"), ({"num_gpus": 2}, "num_gpus"),
    ({"gpu_name": "RTX 3090"}, "gpu_name"), ({"verification": "unverified"}, "verified"),
    ({"verification": None}, "verified"),
])  # fmt: skip
def test_admissible_rejects(bad, why):
    assert OF.reasons(offer(), "RTX_4090") == []
    assert any(why in r for r in OF.reasons(offer(**bad), "RTX_4090"))


def test_vast_raw_error_with_exit_0_stops_with_the_message(monkeypatch):
    """V0: on an API error vastai prints {"error": ...} to stderr, nothing to stdout, and exits 0 (cli/main.py
    run_command). The first launch crashed on json.loads("") and lost the message."""
    err = '{"error": true, "status_code": 400, "msg": "offer no longer available"}\n'
    monkeypatch.setenv("VAST_API_KEY", "k")
    monkeypatch.setattr(RM.subprocess, "run", lambda cmd, **kw: subprocess.CompletedProcess(cmd, 0, "", err))
    with pytest.raises(InfraError, match="offer no longer available"):
        RM.vast("create", "instance", "1", "--raw")


def test_vast_call_has_a_timeout_and_a_timeout_is_a_failed_try(monkeypatch):
    seen = []

    def run(cmd, **kw):
        seen.append(kw["timeout"])
        raise subprocess.TimeoutExpired(cmd, kw["timeout"])

    monkeypatch.setenv("VAST_API_KEY", "k")
    monkeypatch.setenv("VAST_POLL_S", "0")
    monkeypatch.setattr(RM.subprocess, "run", run)
    with pytest.raises(InfraError, match="no answer in 120 s"):
        RM.vast("show", "instances", "--raw", tries=3)
    assert seen == [RM.VAST_TIMEOUT_S] * 3


def test_vast_retries_back_off_exponentially_with_jitter_and_a_cap(monkeypatch):
    monkeypatch.setenv("VAST_POLL_S", "10")
    monkeypatch.setattr(RM.random, "uniform", lambda lo, hi: 1.0)
    assert [RM.backoff_s(i) for i in range(5)] == [10, 20, 40, 60, 60]
    monkeypatch.setattr(RM.random, "uniform", lambda lo, hi: hi)
    assert RM.backoff_s(0) == 15 and RM.backoff_s(9) == 90
    sleeps, tries = [], []
    monkeypatch.setenv("VAST_API_KEY", "k")
    monkeypatch.setattr(RM.time, "sleep", sleeps.append)
    monkeypatch.setattr(RM.subprocess, "run", lambda cmd, **kw: tries.append(1) or subprocess.CompletedProcess(
        cmd, 1, "", "429"))
    with pytest.raises(InfraError, match="429"):
        RM.vast("destroy", "instance", "1", tries=3)
    assert len(tries) == 3 and sleeps == [15, 30]


def test_destroy_tries_six_times(monkeypatch, tmp_path):
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    calls = []
    monkeypatch.setattr(RM, "vast", lambda *a, **k: calls.append((a[:2], k.get("tries", 1))) or [])
    monkeypatch.setattr(RM, "active", lambda *a: {})
    monkeypatch.setattr(RM.time, "sleep", lambda s: None)
    RM.destroy(5, "r", "destroy", {"id": 5})
    assert (("destroy", "instance"), 6) in calls


def test_ssh_and_rsync_carry_timeouts_and_sh_has_one(monkeypatch, tmp_path):
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    argv = RM.ssh({"port": 1})
    assert {"ConnectTimeout=30", "ServerAliveInterval=30", "ServerAliveCountMax=4"} <= set(argv)
    seen = []
    monkeypatch.setattr(RM.subprocess, "run", lambda cmd, **kw: seen.append((cmd, kw["timeout"])) or SimpleNamespace(
        returncode=0))
    RM.rsync({"port": 1}, "a", "b")
    assert {"--timeout=300", "--partial"} <= set(seen[0][0]) and seen[0][1] == RM.SH_TIMEOUT_S

    def hang(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, kw["timeout"])

    monkeypatch.setattr(RM.subprocess, "run", hang)
    with pytest.raises(InfraError, match="rsync timed out"):
        RM.rsync({"port": 1}, "a", "b")


def test_watchdog_curl_has_a_max_time():
    assert "curl -fsS --max-time 30" in (VAST / "watchdog.sh").read_text()


def test_the_lock_is_reentrant_exclusive_and_waits_a_bounded_time(tmp_path, monkeypatch):
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    with LG.locked(), LG.locked():  # reentrant: destroy() holds it and calls ledger()
        LG.ledger(event="x")
        probe = subprocess.run([sys.executable, "-c", (
            "import fcntl, sys; f = open(sys.argv[1], 'a')\n"
            "try: fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)\nexcept BlockingIOError: sys.exit(7)"),
            str(tmp_path / ".lock")])
        assert probe.returncode == 7  # another process cannot take it meanwhile
        monkeypatch.setattr(LG, "LOCK_WAIT_S", 0)
    assert subprocess.run([sys.executable, "-c", (
        "import fcntl, sys; f = open(sys.argv[1], 'a'); fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)"),
        str(tmp_path / ".lock")]).returncode == 0  # released after the block
    holder = open(tmp_path / ".lock", "a")
    import fcntl

    fcntl.flock(holder, fcntl.LOCK_EX)
    with pytest.raises(InfraError, match="held by another"), LG.locked():
        pass
    holder.close()


def test_active_runs_json_is_written_atomically(tmp_path, monkeypatch):
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    LG.update_active(lambda runs: {**runs, "r1": {"instance": 1}})
    LG.update_active(lambda runs: {**runs, "r2": {"instance": 2}})
    assert sorted(LG.active()) == ["r1", "r2"] and not list(tmp_path.glob("*.tmp"))


def test_real_offer_from_v0_is_admissible():
    """A 3090 offer as Vast's search returned it at V0 (2026-11-01): `verification` is a string, and cpu_ram is in
    MB with 1 GB = 1000 MB (Vast returned cpu_ram 32160 for the query cpu_ram>=32). All offers failed before."""
    o = (json.loads((Path(__file__).parents[1] / "fixtures" / "vast_offer_3090.json").read_text())
         | {"geolocation": "Texas, US"})
    assert OF.reasons(o, "RTX_3090") == []
    assert any("cpu_ram" in r for r in OF.reasons(dict(o, cpu_ram=31999), "RTX_3090"))


def test_policy_reads_only_offer_fields_the_pinned_cli_returns():
    """CT-86. Key names only, no offer data: the V0 offer recording (2026-11-01) plus `geolocation`, which the test
    above adds by hand. Re-record from `uvx vastai==1.11.0 search offers ... --raw` when a Vast key is at hand."""
    keys = set(json.loads((Path(__file__).parents[1] / "fixtures" / "vast_offer_keys_1_11_0.json").read_text()))
    read = {t[2] for g in OF.CARDS for t in OF.terms(g) if t[2]}
    read |= {"id", "geolocation", *OF.offer_cols({})}
    assert len(read) > 12 and read <= keys, sorted(read - keys)


# --- ledger and secrets (G9, D16) -----------------------------------------------------------------------------------


def test_ledger_has_no_secret_fields(tmp_path, monkeypatch):
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    assert not [f for f in LG.FIELDS if re.search(r"key|token|secret|pass|identity|login|pat\b", f, re.I)]
    LG.ledger(event="offer", run_id="r", dph_total=0.5)
    with pytest.raises(KeyError):
        LG.ledger(event="x", api_key="k")
    for s in (IDENTITY, PAT, "github_pat_11ABC", "-----BEGIN OPENSSH PRIVATE KEY-----"):
        with pytest.raises(ValueError):
            LG.ledger(event="x", detail=f"oops {s}")
    assert (tmp_path / "ledger.csv").read_text().splitlines()[0] == ",".join(LG.FIELDS)


def test_onstart_length_is_checked_and_run_ids_are_short():
    """Vast refuses an onstart over 4048 characters; the run_id appears in it twice and had no length limit."""
    assert RM.RUN_ID_RE.fullmatch("a" * 40) and not RM.RUN_ID_RE.fullmatch("a" * 41)
    assert len(RM.onstart_cmd("a" * 40, 10**7, idle_s=10**6)) <= 4048
    with pytest.raises(InfraError, match="4048"):
        RM.onstart_cmd("a" * 1000, 7200)


def test_onstart_carries_no_secret_and_is_the_watchdog():
    cmd = RM.onstart_cmd("r1", 7200)
    assert len(cmd) <= 4000 and not [s for s in SECRETS.values() if s in cmd]
    b64 = re.search(r"printf %s (\S+) \| base64 -d \| gunzip >", cmd).group(1)
    assert gzip.decompress(base64.b64decode(b64)).decode() == (VAST / "watchdog.sh").read_text()
    assert "WATCHDOG_CAP_S=7200" in cmd and "CONTAINER_API_KEY" not in cmd


# --- end to end in DRY_RUN=1 (create -> poll -> sync -> run -> pull -> destroy -> verify absent) --------------------


def test_end_to_end_dry_run(E):
    o1 = E.sh("offers.sh", "RTX_4090")
    assert "111" in o1.stdout and not (E.home / "offers").exists()
    (q,) = E.vast("search offers")
    assert "dph<=1.0" in q[2] and "--storage" in q
    E.pick()
    dry = E.sh("launch.sh", "111", "r1", "--cap-hours", "2", DRY_RUN="1")  # printed, never rented
    assert "DRY_RUN: vastai create instance 111" in dry.stdout and PAT not in dry.stdout + dry.stderr
    assert not E.vast("create instance") and E.instances() == {} and not (E.home / "active_runs.json").exists()
    o2 = E.sh("launch.sh", "111", "r1", "--cap-hours", "2")
    (c,) = E.vast("create instance")
    assert {"--ssh", "--direct", "--cancel-unavail", "--onstart-cmd", "--disk", "--image", "--login"} <= set(c)
    assert after(c, "--label") == "r1" and c[2] == "111"
    assert PAT in after(c, "--login") and PAT not in after(c, "--onstart-cmd")
    assert len(E.vast("show instance")) == 2  # polled loading -> running
    assert all(x["api_key"] == API_KEY for x in E.calls("vastai"))  # scoped key via env, never argv
    o3 = E.sh("sync_up.sh", "r1", "--bundle-dir", E.bundles)
    o4 = E.sh("run.sh", "r1", "objdel.cli.train.sft", "--spec", "specs/r1_qwen3vl2b.yaml", "--stage", "train")
    assert "timeout" in o4.stdout and "tmux new-session" in o4.stdout and "HEARTBEAT" in o4.stdout
    o5 = E.sh("pull.sh", "r1")
    assert "PULLED_OK" in o5.stdout
    assert not E.calls("ssh") and not E.calls("rsync")  # dry run: printed, never executed
    o6 = E.sh("teardown.sh", "r1")  # a spending script: DRY_RUN=0 against the fake tools
    assert E.vast("destroy instance") == [["destroy", "instance", "4242", "-y"]]
    assert E.instances() == {} and E.vast("show instances")  # verified absent after destroy
    rsyncs = [ln for o in (o3, o5) for ln in o.stdout.splitlines() if ln.startswith("DRY_RUN: rsync")]
    assert len(rsyncs) >= 3 and all("--exclude=data/raw/" in ln for ln in rsyncs)
    assert [r["event"] for r in E.ledger()] == ["offer", "launch", "sync", "run", "teardown", "invoice"]
    assert float(E.ledger()[-2]["cost_usd"]) == pytest.approx(0.6, abs=0.05)
    assert E.active() == {}
    text = all_text(E, o1, o2, o3, o4, o5, o6)
    assert not [s for s in SECRETS.values() if s in text]


def test_dry_run_never_destroys(E):
    E.launched()
    p = E.sh("teardown.sh", "r1", DRY_RUN="1")
    assert "DRY_RUN: vastai destroy instance 4242 -y" in p.stdout
    assert not E.vast("destroy instance") and list(E.instances()) == ["4242"]  # still rented
    assert "r1" in E.active()  # and still recorded


def test_the_cli_turns_an_infra_error_into_its_message_and_exit_status(E):
    p = E.sh("pull.sh", "nope", ok=False)
    assert p.returncode == 1 and p.stderr.strip() == "vast: no active run nope"


def test_make_without_a_target_lists_the_targets_and_smoke_names_the_preflight_test():
    out = subprocess.run(["make"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    assert all(t in out for t in ("preflight", "smoke", "lint", "test "))
    dry = subprocess.run(["make", "-n", "smoke"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
    assert PF.DRY_RUN in dry


def test_launch_destroys_on_sigterm_while_polling(E):
    # BK4-3: an interrupt (Ctrl-C, SIGTERM, SIGHUP) after `create instance` destroys the instance
    E.state(statuses=["loading"])
    E.pick()
    env = dict(E.env, VAST_LAUNCH_TIMEOUT_S="60", VAST_POLL_S="0.05", DRY_RUN="0")
    argv = ["bash", str(VAST / "launch.sh"), "111", "r1", "--config", str(E.config)]
    p = subprocess.Popen(argv, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    t0 = time.time()
    while len(E.vast("show instance")) < 2 and time.time() - t0 < 30:
        time.sleep(0.05)
    p.send_signal(signal.SIGTERM)
    p.communicate(timeout=60)
    assert p.returncode != 0 and E.instances() == {} and E.ledger()[-1]["event"] == "destroy"


def test_bandwidth_in_filter_cost_guard_and_ledger(E):
    # BK4-18: per-GB bandwidth is in the filter (ceiling), the G10 projection (estimate) and the ledger
    q = OF.query("RTX_4090").split()
    assert f"inet_up_cost<={LG.INET_CEIL}" in q and f"inet_down_cost<={LG.INET_CEIL}" in q
    assert OF.reasons(offer(inet_up_cost=LG.INET_CEIL * 2), "RTX_4090") == ["bandwidth fails the §4.2 filter"]
    assert OF.reasons({k: v for k, v in offer().items() if k != "inet_down_cost"}, "RTX_4090")
    E.state(offers=[offer(inet_up_cost=LG.INET_CEIL, inet_down_cost=LG.INET_CEIL)])
    E.pick()
    spent(E, LG.HARD_STOP - 2.8)  # + 0.55 x 4 = HARD_STOP - 0.6 by price alone; the bandwidth estimate takes it past
    p = E.sh("launch.sh", "111", "r1", "--owner-approved", ok=False)
    assert "hard stop" in p.stderr and not E.vast("create instance")
    E.state(offers=[offer()])
    E.pick()
    (E.home / "ledger.csv").write_text(",".join(LG.FIELDS) + "\n")
    E.sh("launch.sh", "111", "r1")
    row = E.ledger()[-1]
    assert row["event"] == "launch" and float(row["inet_up_cost"]) == 0.001 and float(row["inet_down_cost"]) == 0.001
    E.sh("teardown.sh", "r1")
    inv = E.ledger()[-1]
    assert inv["event"] == "invoice" and "billed=pending" in inv["detail"] and not inv["cost_usd"]
    assert (E.home / "invoices" / "r1.json").exists()


CHARGES_V0 = json.loads((ROOT / "tests" / "fixtures" / "vast_charges_v0.json").read_text())


def test_billed_reads_the_instance_charges_vast_returned_in_v0():
    """V0: `invoices-v1 --charges --raw` lists one `instance` row per instance, `source` = "instance-<id>", amount in
    USD with gpu (from `running`), disk (from create) and bandwidth items; there is no instance_id key."""
    assert LG.billed(CHARGES_V0, 58715600) == pytest.approx(0.029)
    assert LG.billed(CHARGES_V0, 1) is None and LG.billed({"results": []}, 58715600) is None


def test_teardown_records_the_billed_amount_beside_the_estimate(E):
    E.launched()
    E.state(charges={"results": [dict(CHARGES_V0["results"][0], source="instance-4242")], "next_token": None})
    E.sh("teardown.sh", "r1")
    (call,) = E.vast("show invoices-v1")
    assert "--charges" in call and "--raw" in call
    inv = E.ledger()[-1]
    assert inv["event"] == "invoice" and "billed=$0.029" in inv["detail"]
    assert not inv["cost_usd"]  # spend() counts the teardown row once


def test_launch_cap_adds_the_image_pull(E):
    """V0: the watchdog counts from launch, and the 5.88 GB pull took 5.7 min of a 10-min cap. GPU billing starts at
    `running`, so the pull costs disk only; --cap-hours is job time and PULL_S is added to the box's cap."""
    E.launched(cap="2")
    (c,) = E.vast("create instance")
    assert f"WATCHDOG_CAP_S={7200 + RM.PULL_S}" in after(c, "--onstart-cmd") and RM.PULL_S >= 360
    run = E.active()["r1"]
    assert run["deadline"] - run["started"] == pytest.approx(7200 + RM.PULL_S)


def test_launch_timeout_destroys(E):
    E.state(statuses=["loading"])
    E.pick()
    p = E.sh("launch.sh", "111", "r1", ok=False, VAST_LAUNCH_TIMEOUT_S="0.3")
    assert "timeout" in p.stderr and E.vast("destroy instance") and E.instances() == {}
    assert E.ledger()[-1]["event"] == "destroy" and E.active() == {}


@pytest.mark.parametrize("status", ["exited", "unknown", "offline"])
def test_launch_bad_status_destroys(E, status):
    E.state(statuses=["loading", status])
    E.pick()
    E.sh("launch.sh", "111", "r1", ok=False)
    assert E.vast("destroy instance") and E.instances() == {}


def test_one_instance_rule(E):
    E.state(instances={"9": {"id": 9, "label": "other", "dph_total": 0.3, "start_date": time.time()}})
    E.pick()
    p = E.sh("launch.sh", "111", "r1", ok=False)
    assert "--allow-parallel" in p.stderr and not E.vast("create instance")
    E.sh("launch.sh", "111", "r1", "--allow-parallel")
    assert E.vast("create instance")


def test_pick_top_takes_the_cheapest_admissible_offer_of_the_same_search(E):
    """V0: Vast re-lists offers between searches, so `offers.sh` then `--pick <id>` raced and lost; `--pick top`
    picks from the one search."""
    E.state(offers=[offer(222, dph_total=0.4), offer(111, dph_total=0.5)])
    p = E.pick("top")
    assert (E.home / "offers" / "222.json").exists() and "picked 222" in p.stdout


def test_pick_id_takes_any_admissible_offer_not_only_the_five_shown(E):
    """V1 (2026-11-02): the five cheapest admissible 4090s were all in Asia, and a Thailand host had not pulled the
    5.9 GB image in 15 min; a US host beyond the cheapest five could not be picked. Any admissible offer can be."""
    cheap = [offer(100 + i, dph_total=0.3 + i / 100) for i in range(6)]
    E.state(offers=[*cheap, offer(999, verification="unverified")])
    p = E.pick("105")
    assert "105" not in p.stdout.split("picked")[0] and "picked 105" in p.stdout  # not listed, still pickable
    assert "not among the admissible" in E.pick("999", ok=False).stderr  # the filter still decides


@pytest.mark.parametrize("geo,ok", [("Texas, US", True), ("Quebec, CA", True), ("US", True),
                                    ("Thailand, TH", False), ("Mexico, MX", False), ("", False), (None, False)])
def test_only_north_american_hosts_are_admissible(geo, ok):
    """Owner (2026-11-02): US and CA hosts only. V1: a Thailand host had not pulled the 5.9 GB image in 15 min. A
    fixed rule of the filter, so launch's re-check refuses a saved offer from elsewhere; no location fails."""
    assert (OF.reasons(offer(geolocation=geo), "RTX_4090") == []) == ok


def test_a_bad_machine_is_never_admissible():
    """2026-11-03 v2: machine 28360's 4090 was busy and power-capped (7.4 of 165 bf16 TFLOPS)."""
    assert OF.reasons(offer(machine_id=70), "RTX_4090") == []
    assert any("machine 28360" in r for r in OF.reasons(offer(machine_id=28360), "RTX_4090"))


def test_sync_up_checks_the_gpu_before_the_bundle(E):
    """2026-11-03 v2: the GPU check ran only in R0, after a 12 min bundle upload; image_check.sh now gates it."""
    E.launched()
    out = E.sh("sync_up.sh", "r1").stdout.splitlines()
    check = next(i for i, ln in enumerate(out) if "image_check.sh /workspace/runs/r1/image_check.json" in ln)
    bundle = next(i for i, ln in enumerate(out) if ln.startswith("DRY_RUN: rsync") and "bundle-lidc.tar" in ln)
    assert check < bundle
    assert "image_check.sh" not in E.sh("sync_up.sh", "r1", "--code-only").stdout  # mid-session the GPU is busy


def test_sync_up_writes_code_sha_also_with_code_only(E):
    """A3 U49: after the code rsync, the box gets /workspace/src/CODE_SHA (final_eval.code_sha reads it)."""
    E.launched()
    for args in ((), ("--code-only",)):
        out = E.sh("sync_up.sh", "r1", *args).stdout.splitlines()
        code = next(i for i, ln in enumerate(out) if ln.startswith("DRY_RUN: rsync") and "--files-from" in ln)
        sha = next(i for i, ln in enumerate(out) if "/workspace/src/CODE_SHA" in ln)
        assert code < sha and git_sha(ROOT) in out[sha]


def test_offers_search_and_launch_keep_to_north_america(E):
    E.state(offers=[offer(1, dph_total=0.3, geolocation="Thailand, TH"), offer(3, dph_total=0.35, geolocation=None),
                    offer(2, dph_total=0.4, geolocation="Quebec, CA")])  # fmt: skip
    p = E.pick("top")
    assert "picked 2" in p.stdout and " 1 " not in p.stdout.split("picked")[0]
    (q,) = E.vast("search offers")
    assert "geolocation in [US,CA]" in q[q.index("offers") + 1]
    assert "--region" in E.pick("top", "--region", "TH", ok=False).stderr  # no flag widens the rule
    saved = dict(offer(7, geolocation="Thailand, TH"), _gpu="RTX_4090", _disk=80)
    (E.home / "offers" / "7.json").write_text(json.dumps(saved))  # saved before the rule, as offers/49545596.json
    p = E.sh("launch.sh", "7", "r1", ok=False)
    assert "region" in p.stderr and not E.vast("create instance")


def test_launch_needs_a_picked_admissible_offer(E):
    E.sh("launch.sh", "111", "r1", ok=False)  # never picked
    E.state(offers=[offer(dph_total=1.2)])
    E.pick(ok=False)  # above the ceiling: not offered, not saved
    assert not E.vast("create instance")


# --- transfer rules in the scripts -----------------------------------------------------------------------


def test_code_snapshot_ships_the_gpu_scoring_test(E):
    """E8: the chain runs tests/verify/test_gpu_scoring.py on the box; the install's import check names library modules
    and the box stack only (no tests.* and no experiments.*)."""
    files = RM.code_files()
    assert "tests/verify/test_gpu_scoring.py" in files and "tests/oracles/reward_reference.py" in files
    assert "src/objdel/__init__.py" in files and not [f for f in files if f.startswith(("data/", ".venv/", "build/"))]
    E.launched()
    o = E.sh("sync_up.sh", "r1", "--bundle-dir", E.bundles)
    assert "tests/oracles/reward_reference.py" in (E.home / "code_files.txt").read_text().split()
    (inst,) = [ln for ln in o.stdout.splitlines() if "uv pip install" in ln]
    assert "-I -c" in inst and RM.IMPORT_CHECK in inst
    assert "tests." not in RM.IMPORT_CHECK and "experiments." not in RM.IMPORT_CHECK


def test_rsync_never_sends_data_raw(E):
    E.launched()
    p = E.sh("sync_up.sh", "r1", "--bundle-dir", ROOT / "data" / "raw" / "bundles", ok=False)
    assert "data/raw" in p.stderr and "rsync" not in p.stdout
    assert all(not f.startswith("data/") for f in RM.code_files())


SHIPPED = ("reference/frozen_views.json", "reference/val_halves.json")


def _child(cwd: Path, code: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-I", "-c", code], cwd=cwd, capture_output=True, text=True, timeout=120,
                          env=dict(os.environ, OBJDEL_REPO=str(cwd)))


def test_the_money_guard_reads_its_limits_from_reference():
    """ALARM and HARD_STOP come from reference/spend_limits.json (an owner edit, not a code edit); the file ships,
    since the box imports the ledger (run_spec)."""
    lim = json.loads((ROOT / "reference/spend_limits.json").read_text())
    assert (LG.ALARM, LG.HARD_STOP) == (lim["alarm_usd"], lim["hard_stop_usd"]) == (10.0, 31.45)
    assert "reference/spend_limits.json" in RM.code_files()


def test_frozen_views_ship(tmp_path, monkeypatch):
    """I50: the code snapshot keeps the two reference/ records the box reads (final_eval's frozen check, the
    held-out val_half) and drops docs/; on a copy of the snapshot both lookups find their file."""
    repo = tmp_path / "repo"
    for f in (*SHIPPED, "docs/a.md", "docs/reference/data.md", "tests/verify/test_gpu_scoring.py", "pyproject.toml"):
        (repo / f).parent.mkdir(parents=True, exist_ok=True)
        (repo / f).write_text("{}")
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    with monkeypatch.context() as m:
        m.setenv("OBJDEL_REPO", str(repo))
        assert sorted(RM.code_files()) == sorted([*SHIPPED, "tests/verify/test_gpu_scoring.py", "pyproject.toml"])
    files = RM.code_files()  # (b) the real repo
    assert set(SHIPPED) <= set(files)
    assert not [f for f in files if f.startswith("docs/")]

    box = tmp_path / "box"  # (c) the box's /workspace/src
    for f in files:
        if (ROOT / f).is_file():
            (box / f).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / f, box / f)
    view = box / "views" / "lidc-val-edge"
    (view / "lidc").mkdir(parents=True)
    (view / "lidc" / "val.jsonl").write_text("")
    (view / "manifest.json").write_text(json.dumps({"sets": {"lidc/val": {"target_set_sha256": "0"}}}))
    fe = ("import json; from objdel.policy import final_eval as FE; "  # OBJDEL_REPO=box: the box's repo_root()
          f"assert str(FE.FROZEN_VIEWS).startswith({str(box)!r}), FE.FROZEN_VIEWS; "
          "print(json.dumps(FE.frozen_entry('lidc-val-edge')), flush=True); "
          f"FE.check_view(FE.Path({str(view)!r}), 'lidc', 'val')")
    p = _child(box, fe)
    assert p.returncode != 0 and "content_sha256" in json.loads(p.stdout.splitlines()[0])  # the entry, then E05
    assert "content_sha256" in p.stderr
    (box / SHIPPED[0]).unlink()
    p = _child(box, fe)
    assert p.stdout.splitlines()[0] == "null" and p.returncode != 0 and "frozen_views.json" in p.stderr

    # (d) objdel from the copy (the box installs /workspace/src; here the worktree's editable finder steps aside)
    vh = (f"import sys; sys.path.insert(0, {str(box / 'src')!r}); "
          "sys.meta_path[:] = [f for f in sys.meta_path if 'ScikitBuild' not in type(f).__name__]; "
          f"from objdel.data import splits; assert splits.__file__.startswith({str(box)!r}), splits.__file__; "
          "print(splits.val_half('LIDC-IDRI-0026'))")
    p = _child(box, vh)
    assert p.returncode == 0, p.stderr
    assert p.stdout.strip() == json.loads((ROOT / SHIPPED[1]).read_text())["halves"]["LIDC-IDRI-0026"]
    (box / SHIPPED[1]).unlink()
    assert _child(box, vh).returncode != 0


# --- cost guard (G10) and orphans (G6) -----------------------------------------------------------------------------


def spent(E, usd):
    with mock.patch.dict(os.environ, VAST_HOME=str(E.home)):
        LG.ledger(event="teardown", run_id="old", cost_usd=usd)


def test_cost_guard_hard_stop(E):
    """Owner 2026-11-07: a 10 h cap on the one host that had run the image before (three others failed to pull it)."""
    spent(E, LG.HARD_STOP)
    E.pick()
    p = E.sh("launch.sh", "111", "r1", "--owner-approved", ok=False)
    assert "hard stop" in p.stderr and not E.vast("create instance")


def test_cost_guard_allows_the_owners_extra_7_dollars(E):
    spent(E, LG.HARD_STOP - 1.0)  # + 0.55 x (1 h + the 10 min pull) < HARD_STOP
    E.pick()
    E.sh("launch.sh", "111", "r1", "--cap-hours", "1", "--owner-approved")
    assert E.vast("create instance")


def test_cost_guard_projection_and_alarm(E):
    spent(E, LG.HARD_STOP - 12.0)
    E.pick()
    p = E.sh("launch.sh", "111", "r1", "--cap-hours", "23", ok=False)  # (HARD_STOP - 12) + 0.55 * 23 > HARD_STOP
    assert "hard stop" in p.stderr
    spent(E, 1.5)
    p = E.sh("launch.sh", "111", "r1", "--cap-hours", "1", ok=False)  # 13.5 >= alarm
    assert "alarm" in p.stderr and not E.vast("create instance")
    E.sh("launch.sh", "111", "r1", "--cap-hours", "1", "--owner-approved")
    assert E.vast("create instance")


def test_sweep_destroys_orphans_and_overdue(E):
    E.launched()
    now = time.time()
    E.state(instances=E.instances() | {"9": {"id": 9, "label": "r1", "dph_total": 0.3, "start_date": now - 7200}})
    E.sh("sweep_orphans.sh")
    assert E.vast("destroy instance") == [["destroy", "instance", "9", "-y"]] and "4242" in E.instances()
    act = E.active()
    assert list(act) == ["r1"]  # the orphan's label does not remove the live run
    act["r1"]["deadline"] = now - 1
    (E.home / "active_runs.json").write_text(json.dumps(act))
    E.sh("sweep_orphans.sh")
    assert E.instances() == {} and [r["event"] for r in E.ledger()][-2:] == ["orphan", "overdue"]


def test_a_held_pull_lock_never_blocks_the_sweep(E):
    import fcntl

    E.launched()
    now = time.time()
    E.state(instances=E.instances() | {str(i): {"id": i, "label": "x", "dph_total": 0.3, "start_date": now - 7200}
                                       for i in (8, 9)})
    with (E.home / ".pull-r1.lock").open("a") as held:  # a long pull of r1
        fcntl.flock(held, fcntl.LOCK_EX)
        E.sh("sweep_orphans.sh")
    assert sorted(E.instances()) == ["4242"] and [r["event"] for r in E.ledger()][-2:] == ["orphan", "orphan"]


def test_two_pulls_of_the_same_run_serialize(tmp_path, monkeypatch):
    import threading

    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    order, started = [], threading.Event()

    def pull(tag):
        with LG.pull_lock("r1"):
            order.append(f"{tag}+")
            started.set()
            time.sleep(0.3)
            order.append(f"{tag}-")

    def child():
        import fcntl

        with (tmp_path / ".pull-r1.lock").open("a") as fh:  # another process: the lock is per file, not per thread
            started.wait()
            t0 = time.monotonic()
            fcntl.flock(fh, fcntl.LOCK_EX)
            order.append(f"waited {time.monotonic() - t0 > 0.1}")

    th = threading.Thread(target=child)
    th.start()
    pull("a")
    th.join()
    assert order[:2] == ["a+", "a-"] and order[2] == "waited True"
    with LG.pull_lock("r1"), LG.locked():  # other names do not interfere
        pass


def test_sweep_hard_stop_destroys_everything(E):
    E.launched()
    spent(E, LG.HARD_STOP)
    p = E.sh("sweep_orphans.sh")
    assert E.instances() == {} and "hard stop" in p.stderr


# --- onstart watchdog ----------------------------------------------------------------------------------------------


def watchdog(E, ok_file, max_ticks=40, **env):
    e = dict(E.env, WATCHDOG_CAP_S="3600", WATCHDOG_IDLE_S="1", WATCHDOG_TICK_S="0.05", CONTAINER_ID="4242",
             WATCHDOG_MAX_TICKS=str(max_ticks), WATCHDOG_PULLED_OK=str(ok_file), CONTAINER_API_KEY="inst-key",
             RUN_ID="r1", WATCHDOG_ENV_FILE=str(E.tmp / "environment"), WATCHDOG_STATE=str(E.tmp / "watchdog.start"),
             WATCHDOG_SKIP=re.search(r"WATCHDOG_SKIP=(\S+)", RM.onstart_cmd("r1", 1)).group(1))
    e.update(env)
    p = subprocess.run(["bash", str(VAST / "watchdog.sh")], env=e, capture_output=True, text=True,
                       stdin=subprocess.DEVNULL, timeout=60)  # fmt: skip
    return p, [c["argv"] for c in E.calls("curl")]


def test_watchdog_destroys_on_cap(E):
    p, curls = watchdog(E, E.tmp / "PULLED_OK", WATCHDOG_CAP_S="0", FAKE_GPU_UTIL="97")
    assert p.returncode == 0 and len(curls) == 1
    assert "DELETE" in curls[0] and curls[0][-1].endswith("/instances/4242/")
    env = (E.tmp / "environment").read_text()
    assert "RUN_ID=r1" in env and "inst-key" not in env


def test_watchdog_destroys_when_idle_after_pulled_ok(E):
    run = E.tmp / "box" / "runs" / "r1"
    run.mkdir(parents=True)
    ok = run / "PULLED_OK"
    (run / "HEARTBEAT").write_text("1\n")
    p, curls = watchdog(E, ok, max_ticks=5, WATCHDOG_IDLE_S="0", FAKE_TOUCH=str(ok))  # no EXIT: job running
    assert p.returncode == 3 and curls == []
    (run / "EXIT").write_text("0\n")
    p, curls = watchdog(E, ok, FAKE_TOUCH=str(ok))  # the Mac pulls while the GPU is idle, after the job exited
    assert p.returncode == 0 and len(curls) == 1


def test_watchdog_cap_counts_from_the_first_start(E):
    # BK4-17: Vast re-runs the onstart command on every container restart; the cap must not reset
    ok, two_h_ago = E.tmp / "PULLED_OK", int(time.time()) - 7200
    (E.tmp / "watchdog.start").write_text(f"{two_h_ago}\n")
    p, curls = watchdog(E, ok, FAKE_GPU_UTIL="97")
    assert p.returncode == 0 and len(curls) == 1 and "cap" in p.stdout
    (E.tmp / "watchdog.start").unlink()
    p, curls = watchdog(E, ok, FAKE_GPU_UTIL="97", WATCHDOG_CREATED=str(two_h_ago))
    assert p.returncode == 0 and len(curls) == 2  # the Mac's create time is earlier than this start: it counts
    p, curls = watchdog(E, ok, max_ticks=2, FAKE_GPU_UTIL="97", WATCHDOG_STATE=str(E.tmp / "s2"))
    assert p.returncode == 3 and len(curls) == 2  # a fresh start under the cap keeps the box
    env = (E.tmp / "environment").read_text().splitlines()
    assert len(env) == len(set(env))  # restarts do not duplicate /etc/environment lines
    assert "WATCHDOG_CREATED=" in RM.onstart_cmd("r1", 7200)


def test_watchdog_keeps_box_when_not_pulled_or_busy(E):
    ok = E.tmp / "PULLED_OK"
    (E.tmp / "HEARTBEAT").write_text("1\n")
    p, curls = watchdog(E, ok, max_ticks=5, WATCHDOG_IDLE_S="0")  # idle, never pulled
    assert p.returncode == 3 and curls == []
    ok.touch()
    os.utime(ok, (time.time() - 86400,) * 2)  # pulled before the GPU went idle: outputs since then not on the Mac
    p, curls = watchdog(E, ok, max_ticks=5, WATCHDOG_IDLE_S="0")
    assert p.returncode == 3 and curls == []
    p, curls = watchdog(E, ok, max_ticks=5, WATCHDOG_IDLE_S="0", FAKE_GPU_UTIL="88", FAKE_TOUCH=str(ok))  # busy
    assert p.returncode == 3 and curls == []


def _box_run(E):
    run = E.tmp / "box" / "runs" / "r1"
    run.mkdir(parents=True)
    for f in ("HEARTBEAT", "EXIT"):
        (run / f).write_text("0\n")
    return run


def test_watchdog_destroys_a_box_whose_job_never_started(E):
    """No HEARTBEAT ever (run.sh never ran) and the GPU idle 2 x IDLE_S: nothing to lose, so no wait for the cap."""
    ok = E.tmp / "box" / "runs" / "r1" / "PULLED_OK"
    ok.parent.mkdir(parents=True)
    p, curls = watchdog(E, ok, max_ticks=5, WATCHDOG_IDLE_S="2")  # idle < 4 s
    assert p.returncode == 3 and curls == []
    p, curls = watchdog(E, ok, max_ticks=5, WATCHDOG_IDLE_S="0", FAKE_GPU_UTIL="97")  # busy
    assert p.returncode == 3 and curls == []
    p, curls = watchdog(E, ok, max_ticks=5, WATCHDOG_IDLE_S="0")
    assert p.returncode == 0 and len(curls) == 1 and "no job ever started" in p.stdout
    (ok.parent / "HEARTBEAT").write_text("1\n")  # a job ran: only the EXIT + PULLED_OK rule may destroy
    p, curls = watchdog(E, ok, max_ticks=5, WATCHDOG_IDLE_S="0")
    assert p.returncode == 3 and len(curls) == 1


def test_watchdog_ignores_the_run_dir_mtime_and_skipped_dirs(E):
    """pull.sh's `mv PULLED_OK.next PULLED_OK` makes the run dir newer than PULLED_OK, and rolling/ decay/ stable/
    are never pulled: neither may block the idle destroy forever."""
    run, later = _box_run(E), time.time() + 60
    ok = run / "PULLED_OK"
    (run / "r1-x" / "rolling" / "checkpoint-9").mkdir(parents=True)
    (run / "r1-x" / "rolling" / "checkpoint-9" / "optimizer.pt").write_bytes(b"o")
    ok.touch()
    os.utime(ok, (later,) * 2)
    for p in (run / "r1-x/rolling/checkpoint-9/optimizer.pt", run / "r1-x/rolling", run / "r1-x", run):
        os.utime(p, (later + 60,) * 2)
    p, curls = watchdog(E, ok, max_ticks=5, WATCHDOG_IDLE_S="0")
    assert p.returncode == 0 and len(curls) == 1


def test_watchdog_keeps_files_written_during_the_pull(E):
    """PULLED_OK carries the pull's start: a file written while rsync ran is newer, so the box waits for the next
    pull (stamped at the end, it looked pulled and the box could go with final results unpulled)."""
    run, later = _box_run(E), time.time() + 60
    ok = run / "PULLED_OK"
    ok.touch()
    os.utime(ok, (later,) * 2)  # touched as PULLED_OK.next when the pull started, then renamed
    (run / "branches").mkdir()
    (run / "branches" / "report.json").write_text("{}")
    os.utime(run / "branches" / "report.json", (later + 30,) * 2)  # written during the pull
    p, curls = watchdog(E, ok, max_ticks=5, WATCHDOG_IDLE_S="0")
    assert p.returncode == 3 and curls == []


# --- scripts, launchd, docker ---------------------------------------------------------------------------------------


@pytest.mark.parametrize("script", SCRIPTS)
def test_script_help_and_strict_mode(E, script):
    text = (VAST / script).read_text()
    assert "set -euo pipefail" in text and len(text.splitlines()) <= 8  # thin shim
    assert "usage" in E.sh(script, "--help").stdout.lower()


def test_shell_syntax():
    files = [*VAST.glob("*.sh"), *DOCKER.glob("*.sh")]
    assert len(files) >= 10
    for f in files:
        assert subprocess.run(["bash", "-n", str(f)]).returncode == 0, f
        assert "set -euo pipefail" in f.read_text(), f


def test_launchd_plists():
    for name, interval in [("pull", 1200), ("sweep", 3600)]:
        plist = VAST / f"com.objdel.vast.{name}.plist"
        assert f"<integer>{interval}</integer>" in plist.read_text()
        assert "__REPO__/experiments/phase2/vast/" in plist.read_text()
        assert not shutil.which("plutil") or subprocess.run(["plutil", "-lint", str(plist)]).returncode == 0


def test_image_check_python_compiles():
    text = (DOCKER / "image_check.sh").read_text()
    blocks = re.findall(r"<<'PY'\n(.*?)\nPY\n", text, re.S)
    assert blocks
    for b in blocks:
        compile(b, "image_check", "exec")
    words = ("0.6", "objdel._gpucheck", "/opt/vllm/bin/python")
    assert [w for w in words if w not in text] == []


BASE_DIGEST = "python:3.12-slim-trixie@sha256:a6e34c598f2467ed0e9a8d349809fcd8b5c603269512df273a0bb1784edc11b1"


def test_dockerfile_phase2():
    text = (DOCKER / "Dockerfile.phase2").read_text()
    words = ("/opt/venv", "/opt/vllm", "train.txt", "vllm.txt", " age ", "openssh-server", "rsync", "tmux", BASE_DIGEST)
    assert [w for w in words if w not in text] == []
    assert "COPY src" not in text and "COPY . " not in text  # code is never baked in
    assert BASE_DIGEST in (DOCKER / "Dockerfile.vllm_native").read_text()


# --- bundle builder (§4.4) ------------------------------------------------------------------------------------------

LIDC_EVAL = datasets.LIDC.exams  # the bundle requires every exam view of the entry


def data_tree(root: Path) -> Path:
    d = root / "data"
    files = ["prepared/manifest.json", "raw/LIDC/secret.dcm", "eval_views/riga-val-shifted/manifest.json"]
    files += [f"prepared/{s}" for s in ("lidc/val.jsonl", "images/lidc/1.2.3/a.png", "caches/lidc/c.npz",
                                        "riga/val.jsonl", "images/riga/r.png", "caches/riga/c.npz")]  # fmt: skip
    files += [f"views/{s}" for s in ("lidc-train/rows.jsonl", "lidc-train/context/1.2.3.npz", "riga-train/rows.jsonl")]
    files += [f"eval_views/{n}/{f}" for n in LIDC_EVAL for f in ("manifest.json", "lidc/val.jsonl", "images/x.png")]
    files += ["views/lidc-old/rows.jsonl", "checker_eval/manifest.json", "checker_eval/polygons.npz"]
    for f in files:
        (d / f).parent.mkdir(parents=True, exist_ok=True)
        (d / f).write_bytes(f.encode() * 3)
    sets = {f"{ds}/val": {"items": 1} for ds in ("lidc", "riga")}
    man = {"files": {f"{ds}/val.jsonl": {"sha256": "0"} for ds in ("lidc", "riga")}, "sets": sets, "prompts": {}}
    (d / "prepared/manifest.json").write_text(json.dumps(man))
    return d


def test_bundle_lidc_members_manifest_and_determinism(tmp_path):
    d = data_tree(tmp_path)
    a = B.build("lidc", d, tmp_path / "out1")
    for p in d.rglob("*"):
        os.utime(p, (1e9, 1e9))  # mtimes must not matter
    b = B.build("lidc", d, tmp_path / "out2")
    assert a.name == "bundle-lidc.tar" and a.read_bytes() == b.read_bytes()  # byte-identical builds
    side = (tmp_path / "out1" / "bundle-lidc.tar.sha256").read_text()
    assert side == f"{hashlib.sha256(a.read_bytes()).hexdigest()}  bundle-lidc.tar\n"
    with tarfile.open(a) as t:
        names = t.getnames()
        sums = dict(ln.split("  ")[::-1] for ln in t.extractfile("SHA256SUMS").read().decode().splitlines())
        for m in t.getmembers():
            assert m.mtime == 0 and m.uid == m.gid == 0 and m.isfile()
            assert m.name == "SHA256SUMS" or sums[m.name] == hashlib.sha256(t.extractfile(m).read()).hexdigest()
        man = json.loads(t.extractfile("prepared/manifest.json").read())
    assert names[0] == "SHA256SUMS" and names[1:] == sorted(names[1:]) and set(sums) == set(names[1:])
    assert "prepared/images/lidc/1.2.3/a.png" in names and "views/lidc-train/context/1.2.3.npz" in names
    assert all(f"eval_views/{n}/manifest.json" in names for n in LIDC_EVAL)
    assert "checker_eval/polygons.npz" in names and "views/lidc-old/rows.jsonl" in names  # default: every view
    assert not [n for n in names if "riga" in n or "raw" in n]
    assert list(man["files"]) == ["lidc/val.jsonl"] and list(man["sets"]) == ["lidc/val"] and "prompts" in man


def test_bundle_refuses_raw_and_incomplete(tmp_path):
    d = data_tree(tmp_path)
    (d / "views/lidc-train/leak.dcm").symlink_to(d / "raw/LIDC/secret.dcm")
    with pytest.raises(InfraError, match="data/raw|symlink"):
        B.build("lidc", d, tmp_path / "out")
    (d / "views/lidc-train/leak.dcm").unlink()
    with pytest.raises(InfraError, match="data/raw"):
        B.build("lidc", d / "raw", tmp_path / "out")
    shutil.rmtree(d / "eval_views/probe-topology")
    with pytest.raises(InfraError, match="probe-topology"):
        B.build("lidc", d, tmp_path / "out")
    assert not (tmp_path / "out" / "bundle-lidc.tar").exists()


def test_bundle_ships_only_the_named_view(tmp_path):
    d = data_tree(tmp_path)
    names = [n for n, _ in B.members("lidc", d, ["lidc-train"])]
    assert "views/lidc-train/rows.jsonl" in names and not [n for n in names if n.startswith("views/lidc-old/")]
    with pytest.raises(InfraError):
        B.members("lidc", d, ["lidc-nope"])


def test_bundle_ships_rl_train_windows_only_when_named(tmp_path):
    """build_eval_views --split train windows (views/<ds>-train-shifted, -rescaled: 2.7 GB for LIDC) go into a bundle
    only when a run names them; the default bundle stays the SFT views (sync time)."""
    d = data_tree(tmp_path)
    (d / "views/lidc-train-shifted/lidc").mkdir(parents=True)
    (d / "views/lidc-train-shifted/lidc/train.jsonl").write_text("{}")
    default = [n for n, _ in B.members("lidc", d)]
    assert "views/lidc-train/rows.jsonl" in default and not [n for n in default if "lidc-train-shifted" in n]
    named = [n for n, _ in B.members("lidc", d, ["lidc-train", "lidc-train-shifted"])]
    assert "views/lidc-train-shifted/lidc/train.jsonl" in named


def test_bundle_lidc_needs_the_checker_eval_set(tmp_path):
    d = data_tree(tmp_path)
    shutil.rmtree(d / "checker_eval")
    with pytest.raises(InfraError):
        B.members("lidc", d)


@pytest.mark.data
def test_bundle_real_lidc_selection():
    d = ROOT / "data"
    if not all((d / p).is_dir() for p in ("prepared", "views", *(f"eval_views/{n}" for n in LIDC_EVAL))):
        pytest.skip("real prepared / view directories absent")
    files = B.members("lidc", d, ["lidc-bundle-v3"])  # the launch view; no names ships every lidc-* view
    size = sum(len(p) if isinstance(p, bytes) else p.stat().st_size for _, p in files)
    assert 0.5e9 < size < 4e9, size  # §4.4: images 1.3 GB + views + eval views
    assert not [n for n, _ in files if "riga" in n] and len({n for n, _ in files}) == len(files)


def test_bundle_help():
    p = subprocess.run([sys.executable, "-I", "-m", "objdel.infra.bundle", "--help"], capture_output=True,
                       text=True)  # fmt: skip
    assert p.returncode == 0 and "usage" in p.stdout


# --- run chain: train_sft -> select_ckpt -> g_probe from one config (D12) ------------------------------------------


def test_run_chain_for_an_r1_config(E):
    cfg, box = ROOT / "specs/r1_qwen3vl2b.yaml", "/workspace/data/lidc"
    gpu, train, select, probe = CH.chain(cfg, "r1")
    assert gpu[:2] == ["bash", "-c"] and "run_box.sh --build-only" in gpu[2]  # scoring runs on the CUDA checker
    assert train[:2] == ["-m", "objdel.cli.train.sft"] and after(train, "--output-dir") == "r1-qwen3vl2b"
    assert select[:2] == ["-m", "objdel.cli.train.select"] and "--no-vision" not in select
    vision = json.loads(after(select, "--vision"))
    assert vision == {"prepared": f"{box}/prepared", "reference": "/workspace/runs/r1/vision_ref/qwen3-vl-2b.json"}
    assert after(select, "--run") == "/workspace/runs/r1/r1-qwen3vl2b"
    assert probe[:2] == ["-m", "objdel.cli.train.probe"] and "--ckpt" not in probe
    assert after(probe, "--select") == "/workspace/runs/r1/r1-qwen3vl2b/select"
    assert after(probe, "--model") == "qwen3-vl-2b"
    E.launched()
    line = E.tmux(cfg.relative_to(ROOT))
    words = ("objdel.cli.train.sft", "objdel.cli.train.select", "--vision", "objdel.cli.train.probe",
             "/workspace/runs/r1/r1-qwen3vl2b/select", "OBJDEL_RUN_DIR=/workspace/runs/r1")  # fmt: skip
    assert [w for w in words if w not in line] == []


def test_r1_refuses_on_the_box_without_its_g_vision_reference(E):
    E.launched()
    line, ref = E.tmux("specs/r1_qwen3vl2b.yaml"), "/workspace/runs/r1/vision_ref/qwen3-vl-2b.json"
    assert f"test -f {ref}" in line and "--vision-reference-only first" in line  # refused on the box, with the remedy
    assert line.index(f"test -f {ref}") < line.index("objdel.cli.train.sft")
    assert CH.vision_refs(ROOT / "specs/r1_qwen3vl2b.yaml", "r1") == [ref]


def test_each_run_sh_job_has_its_own_tmux_session_and_log(E):
    """Audit: one session and one run.log per run_id made a second job (K3 under R1's load, §4.7 step 6) fail and
    overwrote logs."""
    E.launched()
    out = E.sh("run.sh", "r1", "objdel.cli.eval.run", "--out", "/workspace/runs/r1/k3_load.json").stdout
    assert "tmux new-session -d -s r1-eval_run" in out and "/workspace/runs/r1/eval_run.log" in out
    out = E.sh("run.sh", "r1", "specs/r1_qwen3vl2b.yaml").stdout
    assert "-s r1-r1_qwen3vl2b" in out and "/workspace/runs/r1/r1_qwen3vl2b.log" in out
    E.sh("run.sh", "r1", "experiments.phase1.bench.xarch_hashes", "--out", "x.json")
    E.sh("run.sh", "r1", "objdel.cli.data.store", "--help")  # every command module runs
    E.sh("run.sh", "r1", "objdel.infra.bundle", ok=False)  # a library module is no command


def _exe(path: Path, text: str) -> None:
    path.write_text(text)
    path.chmod(0o755)


def _on_box(E, tmp, *run_args):
    """Run the script run.sh hands to tmux with bash, box paths under tmp (fake python and timeout): the calls made
    as (cwd, OBJDEL_RUN_DIR, argv) lines, the EXIT code and the run dir."""
    inner = shlex.split(shlex.split(shlex.split(E.tmux(*run_args).split("DRY_RUN: ", 1)[1])[-1])[-1])[2]
    (tmp / "src").mkdir(exist_ok=True)
    (tmp / "bin").mkdir(exist_ok=True)
    py, calls, rd = tmp / "bin/python", tmp / "calls", tmp / "runs/r1"
    _exe(py, f'#!/bin/sh\necho "$PWD|$OBJDEL_RUN_DIR|$*" >> {calls}\n')
    _exe(tmp / "bin/timeout", '#!/bin/sh\nshift 3\nexec "$@"\n')
    for box, local in (("/workspace/src", tmp / "src"), ("/workspace/runs", tmp / "runs"), ("/opt/venv/bin/python", py),
                       ("sleep 60", "sleep 1")):  # fmt: skip
        inner = inner.replace(box, str(local))
    env = dict(os.environ, PATH=f"{tmp / 'bin'}:{os.environ['PATH']}")
    subprocess.run(["bash", "-c", inner], env=env, timeout=30, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    lines = calls.read_text().splitlines() if calls.exists() else []
    return [ln.split("|") for ln in lines], (rd / "EXIT").read_text().strip(), rd


def test_run_sh_job_runs_in_the_code_dir_with_the_run_dir_exported(E, tmp_path):
    """Audit: `cd ... && export ... && (heartbeat) &` backgrounded the whole && list, so the job ran in ssh's home dir
    without OBJDEL_RUN_DIR (relative --config paths and train_sft's /workspace/runs check fail on the box)."""
    E.launched()
    calls, code, rd = _on_box(E, tmp_path, "objdel.cli.eval.run", "--out", "x.json")
    assert calls == [[str(tmp_path / "src"), str(rd), "-I -m objdel.cli.eval.run --out x.json"]]
    assert code == "0" and (rd / "eval_run.log").exists() and (rd / "HEARTBEAT").exists()


def test_soft_steps_continue_and_hard_steps_stop_the_chain(E, tmp_path):
    """A soft stage's failure is logged and the chain goes on; a hard stage's failure stops it. chain() reads only a
    run spec in specs/ (E6.3), so the spec sits in a throwaway repo that OBJDEL_REPO names."""
    repo = tmp_path / "repo"
    (repo / "specs").mkdir(parents=True)
    (repo / "pyproject.toml").write_text("")
    (repo / "reference").mkdir()  # the ledger reads its money limits from the repo it runs in
    (repo / "reference/spend_limits.json").write_bytes((ROOT / "reference/spend_limits.json").read_bytes())
    argvs = [(True, ["bash", "-c", "exit 3  # it's"]), (True, ["-m", "experiments.phase2.a", "it's"]),
             (False, ["-m", "experiments.phase2.b"]), (False, ["bash", "-c", "exit 4"]),
             (False, ["-m", "experiments.phase2.c"])]
    (repo / "specs" / "soft.yaml").write_text(json.dumps({
        "schema_version": 1, "name": "soft", "compute": {"gpu": "RTX_5090", "cap_hours": 1.0},
        "stages": [{"kind": "cmd", "name": f"s{i}", "soft": soft, "argv": argv}
                   for i, (soft, argv) in enumerate(argvs)]}))  # fmt: skip
    E.launched()  # against the real repo (its preflight stamp); only the run reads the throwaway one
    E.env["OBJDEL_REPO"] = str(repo)
    calls, code, rd = _on_box(E, tmp_path, "specs/soft.yaml")
    guard, *stages = [c[2] for c in calls]  # E6.5: the spec record runs before any stage
    assert guard == f"-I -m objdel.infra.vast run_spec specs/soft.yaml {rd}/run_spec.json"
    assert stages == ["-I -m experiments.phase2.a it's", "-I -m experiments.phase2.b"]
    log = (rd / "soft.log").read_text()
    assert code == "4" and log.count("step failed, continuing: bash -c 'exit 3") == 1


def test_chain_passes_a_numeric_g_vision_reference_as_is():
    """R0's base distance may be a number; select_ckpt reads a str reference as a file path."""
    from objdel.core import spec as S

    stage = S.SelectStage(name="select", run="{rd}/o", prepared="{data}/prepared", probes="{data}/eval_views",
                          shifted="{data}/eval_views/lidc-val-shifted",
                          vision=S.VisionRef(prepared="{data}/prepared", reference=0.25))
    select = S.to_argv(stage, "/workspace/runs/v1", CH.BOX_DATA.format(dataset="lidc"))
    assert json.loads(after(select, "--vision")) == {"prepared": "/workspace/data/lidc/prepared", "reference": 0.25}


def test_sync_up_code_only_skips_the_bundle(E):
    E.launched()
    out = E.sh("sync_up.sh", "r1", "--code-only").stdout
    rsyncs = [ln for ln in out.splitlines() if ln.startswith("DRY_RUN: rsync")]
    assert len(rsyncs) == 1 and "--files-from" in rsyncs[0] and "bundle-lidc" not in out
    assert "pip install" in out and E.ledger()[-1]["event"] == "sync" and not E.ledger()[-1]["bundle"]


def test_code_snapshot_ships_no_local_vast_state():
    """Audit: `git ls-files -co --exclude-standard` also lists untracked files, so the ledger, offers, invoices and run
    records would reach the rented host. They are ignored; the allowlist (image pin, no secret) stays tracked."""
    state = ("ledger.csv", "offers/", "invoices/", "active_runs.json", "code_files.txt", "known_hosts", ".lock")
    assert not {f for f in RM.code_files() for s in state if f.startswith(f"experiments/phase2/vast/{s}")}


def test_launch_idle_minutes_reach_the_watchdog(E):
    """V1 has decision gaps (stop rule, gates, C2) between jobs; the watchdog's idle destroy (default 30 min after a
    pulled EXIT) can be widened per session."""
    E.pick()
    E.sh("launch.sh", "111", "r1", "--idle-minutes", "90")
    (c,) = E.vast("create instance")
    assert "WATCHDOG_IDLE_S=5400" in after(c, "--onstart-cmd") and len(after(c, "--onstart-cmd")) <= 4048
    assert "WATCHDOG_IDLE_S=1800" in RM.onstart_cmd("r2", 600)  # the default, as watchdog.sh's own


def test_onstart_repairs_the_ssh_key_file_vast_writes_with_bad_modes():
    """V1 (2026-11-02, Quebec host): sshd refused our registered key, "Authentication refused: bad ownership or modes
    for file /root/.ssh/authorized_keys" (vastai logs); attach ssh said "already associated". Known Vast host issue
    (vast-cli #336): onstart runs after Vast writes the key, so it repairs owner and modes; a no-op on good hosts."""
    cmd = RM.onstart_cmd("v1", 600)
    fix = cmd[: cmd.index("watchdog.sh")]
    for part in ("chown root:root /root/.ssh/authorized_keys", "chmod 600 /root/.ssh/authorized_keys",
                 "chmod 700 /root/.ssh", "chmod go-w /root"):  # fmt: skip
        assert part in fix, part
    assert ") &" in fix  # in the background: the watchdog starts at once
    assert len(cmd) <= 4048 and subprocess.run(["bash", "-n", "-c", cmd]).returncode == 0


def test_onstart_turns_off_vast_auto_tmux():
    """V1 (2026-11-02, Texas host): Vast's login script starts tmux on every ssh login, so `ssh host CMD` (remote,
    rsync, sftp) failed with "open terminal failed: not a terminal". ~/.no_auto_tmux is Vast's documented opt-out."""
    cmd = RM.onstart_cmd("v1", 600)
    assert cmd.startswith("touch /root/.no_auto_tmux; ")


def test_teardown_destroys_even_when_the_final_pull_fails(monkeypatch, tmp_path):
    """V1 (2026-11-02): the final pull needs ssh; with ssh refused, teardown stopped before destroy and the box kept
    billing. G8: destroy first-class, then report the failed pull."""
    calls = []
    monkeypatch.setattr(P, "home", lambda: tmp_path)
    monkeypatch.setattr(P, "get_run", lambda rid: {"instance": 7, "port": 40022})
    monkeypatch.setattr(P, "cmd_pull", lambda a: (_ for _ in ()).throw(InfraError("rsync failed: ...")))
    monkeypatch.setattr(P, "destroy", lambda iid, rid, ev: calls.append((iid, ev)))
    monkeypatch.setattr(P, "vast", lambda *a, **k: {"results": []})
    monkeypatch.setattr(P, "ledger", lambda **kw: calls.append(kw["event"]))
    with pytest.raises(InfraError, match="final pull failed.*destroyed"):
        P.cmd_teardown(SimpleNamespace(run_id="v1"))
    assert calls[0] == (7, "teardown") and "invoice" in calls


def _teardown(monkeypatch, tmp_path, pulls):
    """cmd_teardown with each pull attempt taking the next of `pulls` (None: success, else raised); the calls made."""
    calls, pulls = [], list(pulls)

    def pull(a):
        calls.append("pull")
        if (e := pulls.pop(0)) is not None:
            raise e

    monkeypatch.setattr(P, "home", lambda: tmp_path)
    monkeypatch.setattr(P, "get_run", lambda rid: {"instance": 7, "port": 40022})
    monkeypatch.setattr(P, "cmd_pull", pull)
    monkeypatch.setattr(P, "destroy", lambda iid, rid, ev: calls.append((iid, ev)))
    monkeypatch.setattr(P, "vast", lambda *a, **k: {"results": []})
    monkeypatch.setattr(P, "ledger", lambda **kw: calls.append(kw["event"]))
    return calls


def test_teardown_retries_the_final_pull_twice(monkeypatch, tmp_path):
    calls = _teardown(monkeypatch, tmp_path, [InfraError("ssh"), InfraError("ssh"), None])
    P.cmd_teardown(SimpleNamespace(run_id="v1"))
    assert calls[:4] == ["pull", "pull", "pull", (7, "teardown")]
    calls = _teardown(monkeypatch, tmp_path, [InfraError("ssh")] * 3)
    with pytest.raises(InfraError, match="final pull failed"):
        P.cmd_teardown(SimpleNamespace(run_id="v1"))
    assert calls[:4] == ["pull", "pull", "pull", (7, "teardown")]


def test_teardown_destroys_on_ctrl_c_during_the_final_pull(monkeypatch, tmp_path):
    """Only InfraError was caught: a Ctrl-C (KeyboardInterrupt) in a long final pull skipped the destroy (G8)."""
    calls = _teardown(monkeypatch, tmp_path, [KeyboardInterrupt()])
    with pytest.raises(KeyboardInterrupt):
        P.cmd_teardown(SimpleNamespace(run_id="v1"))
    assert calls == ["pull", (7, "teardown")]


def test_run_cap_kills_a_job_that_ignores_sigterm(E):
    """A process stuck in CUDA can ignore timeout's SIGTERM; -k sends SIGKILL 120 s later."""
    E.launched()
    assert re.search(r"timeout -k 120 \d+s bash -c", E.tmux("objdel.cli.eval.run"))


def test_pull_stamps_its_start_and_skips_full_checkpoints(E):
    """PULLED_OK carries the pull's start (files written during a long pull are newer than it), and rolling/,
    decay/ and stable/ (full checkpoints with optimizer state, 29 GB in R1) stay on the box."""
    E.launched()
    out = E.sh("pull.sh", "r1").stdout.splitlines()
    rd = "/workspace/runs/r1"
    (touch,) = [i for i, ln in enumerate(out) if f"touch {rd}/PULLED_OK.next" in ln]
    (pull,) = [i for i, ln in enumerate(out) if ln.startswith("DRY_RUN: rsync")]
    (mv,) = [i for i, ln in enumerate(out) if f"mv {rd}/PULLED_OK.next {rd}/PULLED_OK" in ln]
    assert touch < pull < mv and not [ln for ln in out if f"touch {rd}/PULLED_OK'" in ln]
    argv = shlex.split(out[pull].removeprefix("DRY_RUN: "))
    assert {f"--exclude={x}" for x in ("rolling/", "decay/", "stable/", "PULLED_OK", "PULLED_OK.next")} <= set(argv)
    assert "--exclude=branches/" not in argv  # model-only branch checkpoints: select_ckpt and the Mac need them
    assert f"WATCHDOG_SKIP={','.join(RM.PULL_EXCLUDE)}" in RM.onstart_cmd("r1", 600)  # one list for pull and watchdog


def test_code_files_keep_paths_with_spaces_whole(tmp_path, monkeypatch):
    """V1 (2026-11-02): `git ls-files` output was split on any whitespace, so a tracked path with spaces (docs/research
    notes) became words that rsync --files-from could not stat, and sync_up failed before sending a byte."""
    for p in ("tests/verify/test_gpu_scoring.py", "experiments/my dir/a b.py", "docs/notes with space/x.md",
              "pyproject.toml"):
        (tmp_path / p).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / p).write_text("x")
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "add", "-A"], cwd=tmp_path, check=True)
    monkeypatch.setenv("OBJDEL_REPO", str(tmp_path))
    assert sorted(RM.code_files()) == ["experiments/my dir/a b.py", "pyproject.toml",
                                       "tests/verify/test_gpu_scoring.py"]


def _git_repo(root, files):
    for f in files:
        (root / f).parent.mkdir(parents=True, exist_ok=True)
        (root / f).write_text("x")
    subprocess.run(["git", "init", "-q"], cwd=root, check=True)


def test_code_snapshot_drops_agent_and_archive_folders(tmp_path, monkeypatch):
    """CLEAN1: .claude/ and archive/ never ship, even untracked; tracked-but-moved files and look-alikes are handled."""
    keep = ["tests/verify/test_gpu_scoring.py", "src/a.py", "src/archive_util.py", "tests/data_helpers.py",
            "experiments/p/riga.py", "pyproject.toml"]
    drop = [".claude/x", ".claude/worktrees/w/a.py", "archive/x/y.py", "data/d.npz", "docs/a.md",
            "experiments/p/lab/l.py", "experiments/p/results/r.json", "experiments/p/RIGA_x.CSV"]  # fmt: skip
    _git_repo(tmp_path, keep + drop + ["src/gone.py"])
    subprocess.run(["git", "add", "src/gone.py"], cwd=tmp_path, check=True)
    (tmp_path / "src/gone.py").unlink()  # tracked, deleted on disk: rsync --files-from would abort on it
    monkeypatch.setenv("OBJDEL_REPO", str(tmp_path))
    assert sorted(RM.code_files()) == sorted(keep)
    monkeypatch.undo()
    assert not [f for f in RM.code_files() if f.startswith((".claude/", "archive/"))]  # the real repo


def test_agent_worktrees_are_git_ignored():
    assert subprocess.run(["git", "check-ignore", "-q", ".claude/worktrees/x/y"], cwd=ROOT).returncode == 0


def test_bundle_imports_standalone(tmp_path):
    """CLEAN1: the files code_files() returns, plus the built extension, import by themselves ($0 IMPORT_CHECK)."""
    pytest.importorskip("cv2")
    pytest.importorskip("scipy")
    import importlib.util
    import site

    for f in RM.code_files():
        (tmp_path / f).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / f, tmp_path / f)
    ext = Path(importlib.util.find_spec("objdel._gpucheck").origin)
    shutil.copy2(ext, tmp_path / "src" / "objdel" / ext.name)
    check = [m.strip() for m in RM.IMPORT_CHECK.removeprefix("import ").split(",")]
    mods = [m for m in check if m != "vllm"] + ["objdel", *RM.LIVE_CLIS]  # vllm: box only (the image has it)
    # -S: without it the editable install's .pth files serve objdel and tests from the repo, whatever sys.path says
    code = (
        "import sys, json, importlib\n"
        "root, sp, mods = sys.argv[1], json.loads(sys.argv[2]), json.loads(sys.argv[3])\n"
        "sys.path[:0] = [root + '/src', root, *sp]\n"
        "[importlib.import_module(m) for m in mods]\n"
        "out = []\n"
        "for n, m in sorted(sys.modules.items()):\n"
        "    if n.split('.')[0] in ('objdel', 'experiments', 'tests'):\n"
        "        out += [m.__file__] if getattr(m, '__file__', None) else list(m.__path__)\n"
        "print(json.dumps(out))\n"
    )
    argv = [sys.executable, "-I", "-S", "-c", code, str(tmp_path), json.dumps(site.getsitepackages()), json.dumps(mods)]
    p = subprocess.run(argv, cwd=tmp_path, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-2000:]
    files = json.loads(p.stdout)
    assert len(files) > 10
    assert [f for f in files if not Path(f).resolve().is_relative_to(tmp_path.resolve())] == []


def test_sync_up_sends_both_bundle_files_to_the_box_data_dir(E):
    """V1 (2026-11-02): rsync(run, src, dst, *extra) was called with two sources, so the .sha256 became the destination
    and the remote path an option ("local file in list of remote sources"). Sources first, the remote dir last."""
    E.launched()
    out = E.sh("sync_up.sh", "r1").stdout
    (line,) = [ln for ln in out.splitlines() if ln.startswith("DRY_RUN: rsync") and "bundle-lidc.tar" in ln]
    argv = shlex.split(line.removeprefix("DRY_RUN: "))
    assert argv[-1].endswith(":/workspace/data/") and [Path(a).name for a in argv[-3:-1]] == [
        "bundle-lidc.tar", "bundle-lidc.tar.sha256"]
    assert not any(a.startswith("root@") for a in argv[:-1])


def test_run_refuses_extra_arguments_after_a_config(E):
    """V1 (2026-11-02): `run.sh v1 configs/r0.yaml --cap-hours 1.5` put the flag into the module's arguments, which a
    config run ignores, so R0 ran under the instance cap instead of 1.5 h. Flags go before the config."""
    E.launched()
    cfg = "specs/r1_qwen3vl2b.yaml"
    r = E.sh("run.sh", "r1", cfg, "--cap-hours", "1.5", ok=False)
    assert "before the config" in r.stderr and "--cap-hours" in r.stderr
    assert "cap_s=5400" in E.sh("run.sh", "r1", "--cap-hours", "1.5", cfg).stdout + \
        E.ledger()[-1]["detail"]


def test_launch_gate_refuses_without_a_passing_current_preflight(E):
    """Missing, failed and stale (code, data or config changed) stamps each stop the rental before any Vast call."""
    E.pick()
    for stamp, why in (("missing", "no preflight stamp"), ("failed", "FAILED preflight"),
                       ("data", "stale (data"), ("code", "stale (code"), ("config", "stale (config")):  # fmt: skip
        if stamp == "missing":
            (Path(E.env["OBJDEL_RUN_DIR"]) / "preflight" / f"{E.config.stem}.json").unlink()
        elif stamp == "failed":
            E.stamp(ok=False)
        else:
            E.stamp(**{stamp: "0" * 64})
        p = E.sh("launch.sh", "111", "r1", ok=False)
        assert why in p.stdout + p.stderr and not E.vast("create instance")
    E.stamp()
    E.sh("launch.sh", "111", "r1")  # a passing, current stamp launches
    assert E.vast("create instance")


def test_launch_is_bounded_by_the_spec_compute(E):
    """E6.8b: --cap-hours above compute.cap_hours, or an offer whose card is not compute.gpu, dies before any Vast
    call."""
    E.pick()
    p = E.sh("launch.sh", "111", "r1", "--cap-hours", "24.5", ok=False)
    assert "compute.cap_hours 24.0" in p.stderr and not E.vast("create instance")
    E.config.write_text(E.config.read_text().replace("RTX_4090", "RTX_5090"))
    E.stamp()
    p = E.sh("launch.sh", "111", "r1", ok=False)
    assert "RTX_4090 is not compute.gpu RTX_5090" in p.stderr and not E.vast("create instance")


# --- vllm_native image: one venv (perf/implementation.md B3, B4, C4) ------------------------------------------------


def test_allowlist_image_is_the_vllm_native_pin():
    image = json.loads((VAST / ".vast-allowlist").read_text())["image"]
    assert re.fullmatch(r"ghcr\.io/yashwanthsai1234/objdel-gpu:vllm_native-[0-9a-f]{12}@sha256:[0-9a-f]{64}", image)


def test_image_check_handles_one_and_two_venvs():
    text = (DOCKER / "image_check.sh").read_text()
    assert "/opt/vllm/bin/python" in text and "/opt/venv/bin/python" in text
    assert re.search(r"if \[ -x /opt/vllm/bin/python \]", text)  # two venvs: vLLM in the second
    assert "else" in text  # one venv: vLLM in /opt/venv


ENV_CHECKS = ("nvidia_smi", "capability", "bf16_matmul", "triton_compile", "torch_compile", "vllm_engine", "nvcc",
              "trl_private_names")


def test_env_check_lists_every_check():
    text = (DOCKER / "env_check.sh").read_text()
    assert "set -euo pipefail" in text and "env_check.json" in text
    for key in ENV_CHECKS:
        assert f'"{key}"' in text or f"'{key}'" in text, key
    for word in ("(12, 0)", "sm_120", "VLLM_USE_FLASHINFER_SAMPLER=0", "--list-gpu-arch", "attention_backend",
                 "pcie", "power.limit", "driver_version"):
        assert word in text, word
    for b in re.findall(r"<<'PY'\n(.*?)\nPY\n", text, re.S):
        compile(b, "env_check", "exec")
    # Triton's @jit reads its kernel's source: the checks run from a file, not stdin (box p1, 2026-11-06)
    assert 'cat > "$CHECKS" <<\'PY\'' in text and '/opt/venv/bin/python -I "$CHECKS" "$OUT"' in text


def test_teardown_of_a_box_that_never_got_ssh_destroys_it_without_a_pull(monkeypatch, tmp_path):
    """2026-11-07: a box torn down while its image was still loading has no ssh port; the pull raised KeyError
    'port' (destroy still ran, from the finally) and the invoice was never recorded. Nothing is on such a box."""
    calls = _teardown(monkeypatch, tmp_path, [])
    monkeypatch.setattr(P, "get_run", lambda rid: {"instance": 7})
    P.cmd_teardown(SimpleNamespace(run_id="v1"))
    assert calls == [(7, "teardown"), "invoice"]


def test_launch_refusal_lets_the_v1_offer_through_and_stops_at_the_hard_stop(tmp_path, monkeypatch):
    """6.1 step 0: with $16.90 closed spend, a $0.90/h offer and a 2.4 h cap (plus PULL_S) is within HARD_STOP."""
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    (tmp_path / "ledger.csv").write_text(",".join(LG.FIELDS) + "\n" + "," * 13 + "16.90,\n")
    assert LG.HARD_STOP == 31.45 and LG.ALARM == 10.0
    offer = {"dph_total": 0.90, "inet_up_cost": 0.01, "inet_down_cost": 0.01}
    assert P.launch_refusal(offer, 2.4) == ""
    assert "hard stop" in P.launch_refusal(dict(offer, dph_total=7.0), 2.4)


def test_launch_refuses_without_a_ledger(E):
    """E8.1: a missing ledger.csv is not $0 spend (a moved VAST home would read none): launch refuses before create."""
    E.pick()
    (E.home / "ledger.csv").unlink()
    p = E.sh("launch.sh", "111", "r1", ok=False)
    assert "no ledger.csv" in p.stderr and not E.vast("create instance")


def test_code_files_ship_reference_but_nothing_from_docs(tmp_path, monkeypatch):
    """docs/ is local only (owner, 2026-11-11): no docs/ path ships, the kernel drafts included (they shipped for
    R-D1 until then); the two reference/ tables do."""
    keep = ["reference/frozen_views.json", "reference/val_halves.json", "tests/verify/test_gpu_scoring.py",
            "pyproject.toml"]
    drop = ["docs/implementation2.0/kernels/K1/k1poly.cpp", "docs/implementation2.0/kernels/K3/GPU_RUN.md",
            "docs/implementation2.0/kernels/K1/build/x", "docs/implementation2.0/kernels/K1/results/x.txt",
            "docs/implementation2.0/kernels/README.md", "docs/implementation2.0/kernels/K4/a.cpp",
            "docs/reference/other.json", "docs/implementation2.0/Implementation.md"]
    for p in keep + drop:
        (tmp_path / p).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / p).write_text("x")
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "add", "-A"], cwd=tmp_path, check=True)
    monkeypatch.setenv("OBJDEL_REPO", str(tmp_path))
    assert sorted(RM.code_files()) == sorted(keep)


# --- E6.5 (I16): runs/<id>/run_spec.json ------------------------------------------------------------------------------


def test_a_spec_run_writes_run_spec_json_first(E):
    """run.sh on a spec: the first (hard) job on the box writes /workspace/runs/<id>/run_spec.json beside the log."""
    E.launched()
    line = E.tmux("specs/b1_rows.yaml")
    assert "-m objdel.infra.vast run_spec specs/b1_rows.yaml /workspace/runs/r1/run_spec.json &&" in line
    assert line.index("run_spec.json") < line.index("run_box.sh --build-only")


def test_run_spec_json_holds_the_spec_code_and_versions_and_never_overwrites_another_spec(tmp_path, monkeypatch):
    from importlib.metadata import version

    from objdel.core.spec import RunSpec, load_spec
    from objdel.verify.reward.spec import fingerprint

    monkeypatch.setenv("OBJDEL_REPO", str(tmp_path))
    (tmp_path / "pyproject.toml").touch()  # repo_root refuses an OBJDEL_REPO without one (E2)
    (tmp_path / "specs").mkdir()
    spec = tmp_path / "specs/b1_rows.yaml"
    spec.write_text((ROOT / "specs/b1_rows.yaml").read_text())
    (tmp_path / "CODE_SHA").write_text("a" * 40 + "\n")
    out = tmp_path / "rd/run_spec.json"
    P.cmd_run_spec(SimpleNamespace(spec="specs/b1_rows.yaml", out=str(out)))
    rec = json.loads(out.read_text(), parse_constant=lambda c: pytest.fail(f"not strict JSON: {c}"))
    assert RunSpec.model_validate(rec["spec"]) == load_spec(spec)
    assert rec["git_sha"] == "a" * 40 and rec["reward_fingerprint"] == fingerprint()
    libs = ("trl", "torch", "transformers")
    assert {k: rec[k] for k in libs} == {k: version(k) for k in libs}
    P.cmd_run_spec(SimpleNamespace(spec="specs/b1_rows.yaml", out=str(out)))  # the same spec again: fine
    (tmp_path / "CODE_SHA").unlink()
    spec.write_text(spec.read_text().replace("cap_hours: 9.5", "cap_hours: 2.0"))
    with pytest.raises(InfraError, match="different spec"):
        P.cmd_run_spec(SimpleNamespace(spec="specs/b1_rows.yaml", out=str(out)))
    out.unlink()
    P.cmd_run_spec(SimpleNamespace(spec="specs/b1_rows.yaml", out=str(out)))
    assert json.loads(out.read_text())["git_sha"] == "unknown"


def test_launch_without_an_image_in_the_allowlist_refuses(monkeypatch, tmp_path):
    """The launch dies on a missing `image` after every money check passed (empty allowlist, empty ledger)."""
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    (tmp_path / "ledger.csv").write_text(",".join(LG.FIELDS) + "\n")
    (tmp_path / "offers").mkdir()
    offer = {"_gpu": "RTX_4090", "_disk": 80, "dph_total": 0.5, "inet_up_cost": 0.0,
             "inet_down_cost": 0.0}
    (tmp_path / "offers/1.json").write_text(json.dumps(offer))
    compute = SimpleNamespace(cap_hours=4, gpu="RTX_4090")
    monkeypatch.setattr(P, "load_spec", lambda p: SimpleNamespace(compute=compute))
    monkeypatch.setattr(P, "reasons", lambda *a: [])
    monkeypatch.setattr(P, "vast", lambda *a, **k: [])
    monkeypatch.setattr("objdel.infra.preflight.refusal", lambda cfg: "")
    a = SimpleNamespace(run_id="r1", config="c.yaml", offer=1, cap_hours=2.0, owner_approved=False,
                        allow_parallel=False, idle_minutes=30)
    with pytest.raises(InfraError, match="no image in .vast-allowlist"):
        P.cmd_launch(a)


def test_onstart_ships_the_code_trees_watchdog_not_vast_homes(tmp_path, monkeypatch):
    """VAST_HOME holds state only: a launch from another tree must not ship that folder's older watchdog.sh."""
    (tmp_path / "watchdog.sh").write_text("#!/usr/bin/env bash\necho stale\n")
    monkeypatch.setenv("VAST_HOME", str(tmp_path))
    b64 = re.search(r"printf %s (\S+) \| base64 -d \| gunzip >", RM.onstart_cmd("r1", 7200)).group(1)
    assert gzip.decompress(base64.b64decode(b64)).decode() == (VAST / "watchdog.sh").read_text()
