"""Vast state under VAST_HOME (default <repo>/experiments/phase2/vast): the ledger (G9), active runs, offers, the
allowlist, and the money: what runs accrue, what was spent, what Vast billed (G10, BK4-18)."""

import csv
import json
import os
import re
import time
from pathlib import Path

from objdel.core.io import read_json
from objdel.core.provenance import repo_root

# G10, USD spent in Phase 2: $16.90 at plan start + the approved $3.00 + $1.50 (2026-11-09)
ALARM, HARD_STOP = 10.0, 31.45
# BK4-18: $/GB ceiling, up and down (> 2x the A3 snapshot medians of the five cards, $4-13.65/TB; unit UNVERIFIED until
# V0). Estimate: image 7.5 GB + bundle 2 GB + model downloads once; checkpoint pulls of 3-5 GB per 20 min (§4.4).
INET_CEIL, DOWN_GB, UP_GB_PER_H = 0.03, 15, 6
FIELDS = tuple(
    "ts event run_id instance_id offer_id gpu dph_total inet_up_cost inet_down_cost hosting_type host_id machine_id "
    "bundle cost_usd detail".split()
)
SECRET_RE = re.compile(r"AGE-SECRET-KEY-|ghp_|github_pat_|PRIVATE KEY|Bearer ")


def home() -> Path:
    """The Vast state folder: VAST_HOME, else <repo>/experiments/phase2/vast (ledger, active runs, allowlist, offers,
    known_hosts; watchdog.sh ships from the code tree: remote.watchdog_script), never this package's folder (spend()
    would read $0 there)."""
    return Path(os.environ.get("VAST_HOME", repo_root() / "experiments/phase2/vast"))


def ledger(**row) -> None:
    """Append one row to ledger.csv (G9): known FIELDS only, never a secret-looking value."""
    if bad := set(row) - set(FIELDS):
        raise KeyError(f"not a ledger field: {sorted(bad)}")
    if any(SECRET_RE.search(str(v)) for v in row.values()):
        raise ValueError("refusing to write a secret-looking value to the ledger")
    f = home() / "ledger.csv"
    new = not f.exists()
    with f.open("a", newline="") as fh:
        w = csv.DictWriter(fh, FIELDS)
        if new:
            w.writeheader()
        w.writerow({"ts": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), **row})


def load(name: str) -> dict:
    """A JSON state file under VAST_HOME ({} when absent): active_runs.json, .vast-allowlist, offers/<id>.json."""
    f = home() / name
    return read_json(f) if f.exists() else {}


def active(update: dict | None = None) -> dict:
    """The live runs (active_runs.json); with `update`, write that dict first."""
    if update is not None:
        (home() / "active_runs.json").write_text(json.dumps(update, indent=1))
    return load("active_runs.json")


def accrued(r: dict, t: float) -> float:
    """A run's cost from its start to t: dph x hours + the bandwidth estimate (BK4-18)."""
    h = max(0.0, t - r["started"]) / 3600
    return r["dph"] * h + r.get("inet_down_cost", 0) * DOWN_GB + r.get("inet_up_cost", 0) * UP_GB_PER_H * h


def spend(committed: bool = False) -> float:
    """Closed sessions from the ledger plus what live runs have accrued so far (committed: to their deadlines)."""
    f, now = home() / "ledger.csv", time.time()
    closed = sum(float(r["cost_usd"] or 0) for r in csv.DictReader(f.open())) if f.exists() else 0.0
    return closed + sum(accrued(r, max(now, r["deadline"]) if committed else now) for r in active().values())


def billed(charges: dict, iid: int) -> float | None:
    """USD Vast charged instance `iid` (V0: `invoices-v1 --charges` rows carry it only as `source` "instance-<id>")."""
    rows = [r["amount"] for r in charges.get("results", []) if r.get("source") == f"instance-{iid}"]
    return round(sum(rows), 4) if rows else None
