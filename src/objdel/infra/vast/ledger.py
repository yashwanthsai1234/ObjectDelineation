"""Vast state under VAST_HOME (default <repo>/experiments/phase2/vast): the ledger, active runs, offers, the allowlist,
and the money: what runs accrue, what was spent, what Vast billed."""

import contextlib
import csv
import fcntl
import json
import os
import re
import time
from pathlib import Path

from objdel.core.io import read_json, write_atomic
from objdel.core.provenance import repo_root
from objdel.infra.errors import InfraError

# The money guard: alarm at ALARM USD, refuse past HARD_STOP USD; the owner edits them in reference/spend_limits.json
ALARM, HARD_STOP = (float(read_json(repo_root() / "reference/spend_limits.json")[k]) for k in ("alarm_usd",
                                                                                              "hard_stop_usd"))
# $/GB ceiling, up and down (over 2x the median of the five cards, $4-13.65/TB). Estimate: image 7.5 GB + bundle 2 GB
# + model downloads once; checkpoint pulls of 3-5 GB per 20 min.
INET_CEIL, DOWN_GB, UP_GB_PER_H = 0.03, 15, 6
FIELDS = tuple(
    "ts event run_id instance_id offer_id gpu dph_total inet_up_cost inet_down_cost hosting_type host_id machine_id "
    "bundle cost_usd detail".split()
)
LOCK_WAIT_S = 60  # a read-modify-write takes milliseconds
_held: dict[str, list] = {}  # lock name -> [open file, depth]
SECRET_RE = re.compile(r"AGE-SECRET-KEY-|ghp_|github_pat_|PRIVATE KEY|Bearer ")


def home() -> Path:
    """The Vast state folder: VAST_HOME, else <repo>/experiments/phase2/vast (ledger, active runs, allowlist, offers,
    known_hosts; watchdog.sh ships from the code tree: remote.watchdog_script), never this package's folder (spend()
    would read $0 there)."""
    return Path(os.environ.get("VAST_HOME", repo_root() / "experiments/phase2/vast"))


@contextlib.contextmanager
def locked(name: str = ".lock", wait_s: float | None = None):
    """fcntl.flock on <home>/<name>, reentrant per name in one process; InfraError after `wait_s`. The default lock
    guards the quick read-modify-writes of active_runs.json and ledger.csv; `pull_lock` is the per-run one."""
    held, wait_s = _held.setdefault(name, [None, 0]), LOCK_WAIT_S if wait_s is None else wait_s
    if held[1] == 0:
        fh = (home() / name).open("a")
        deadline = time.monotonic() + wait_s
        while True:
            try:
                fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.monotonic() > deadline:
                    fh.close()
                    raise InfraError(f"{fh.name} is held by another vast command for over {wait_s:g} s") from None
                time.sleep(0.2)
        held[0] = fh
    held[1] += 1
    try:
        yield
    finally:
        held[1] -= 1
        if held[1] == 0:
            held[0].close()  # closing releases the flock
            held[0] = None


def pull_lock(run_id: str):
    """The lock of one run's pulls (a pull and a teardown's final pull never share PULLED_OK); waits up to 1 h, the
    longest a pull may run."""
    return locked(f".pull-{run_id}.lock", 3600)


def ledger(**row) -> None:
    """Append one row to ledger.csv: known FIELDS only, never a secret-looking value."""
    if bad := set(row) - set(FIELDS):
        raise KeyError(f"not a ledger field: {sorted(bad)}")
    if any(SECRET_RE.search(str(v)) for v in row.values()):
        raise ValueError("refusing to write a secret-looking value to the ledger")
    f = home() / "ledger.csv"
    with locked(), f.open("a", newline="") as fh:
        w = csv.DictWriter(fh, FIELDS)
        if fh.tell() == 0:
            w.writeheader()
        w.writerow({"ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), **row})


def load(name: str) -> dict:
    """A JSON state file under VAST_HOME ({} when absent): active_runs.json, .vast-allowlist, offers/<id>.json."""
    f = home() / name
    return read_json(f) if f.exists() else {}


def active(update: dict | None = None) -> dict:
    """The live runs (active_runs.json); with `update`, write that dict first (atomically). A caller that reads, then
    writes, holds `locked()` around both."""
    if update is not None:
        write_atomic(home() / "active_runs.json", json.dumps(update, indent=1))
    return load("active_runs.json")


def update_active(change) -> None:
    """active_runs.json := change(its current dict), read and written under the lock."""
    with locked():
        active(change(active()))


def accrued(r: dict, t: float) -> float:
    """A run's cost from its start to t: dph x hours + the bandwidth estimate."""
    h = max(0.0, t - r["started"]) / 3600
    return r["dph"] * h + r.get("inet_down_cost", 0) * DOWN_GB + r.get("inet_up_cost", 0) * UP_GB_PER_H * h


def spend(committed: bool = False) -> float:
    """Closed sessions from the ledger plus what live runs have accrued so far (committed: to their deadlines)."""
    f, now = home() / "ledger.csv", time.time()
    closed = sum(float(r["cost_usd"] or 0) for r in csv.DictReader(f.open())) if f.exists() else 0.0
    return closed + sum(accrued(r, max(now, r["deadline"]) if committed else now) for r in active().values())


def billed(charges: dict, iid: int) -> float | None:
    """USD Vast charged instance `iid` (`invoices-v1 --charges` rows carry it only as `source` "instance-<id>")."""
    rows = [r["amount"] for r in charges.get("results", []) if r.get("source") == f"instance-{iid}"]
    return round(sum(rows), 4) if rows else None
