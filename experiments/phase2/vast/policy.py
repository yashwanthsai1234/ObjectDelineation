"""Vast AI operations for Phase 2 (PLAN §4, D16, G1-G10). Every decision and every command is built here; the
`*.sh` files next to this module are symlinks to the one shim `vast.sh` (the link's name is the command), and
`watchdog.sh` is the only logic that runs on the box.

    offers.sh GPU [--riga] [--disk GB] [--pick OFFER_ID|top]   search (§4.2 filter, US/CA hosts), top 5; --pick saves it
    launch.sh OFFER_ID RUN_ID --config CONFIG.yaml [--cap-hours H] [--idle-minutes M] [--allow-parallel] [--owner-approved]
    sync_up.sh RUN_ID [--bundle lidc|riga] [--bundle-dir DIR] [--code-only]   code (git ls-files) + bundle, SHA-256
    run.sh RUN_ID MODULE|CONFIG.yaml [ARGS]                    timeout + tmux + heartbeat; a config runs its chain
    pull.sh [RUN_ID]                                           rsync runs/<id>/ to the Mac (not PULL_EXCLUDE), then
                                                               PULLED_OK = the pull's start time (no id: all)
    push_hf.sh RUN_ID CONFIG.yaml                              the run's selected branch to the private HF repo
    fetch_hf.sh RUN_ID REVISION                                a revision of that repo onto the box (runs/RUN_ID/ckpt)
                                                               (hf_push; after pull.sh, before teardown.sh)
    teardown.sh RUN_ID                                         final pull, destroy, verify absent, invoices, ledger
    sweep_orphans.sh                                           destroy unknown / overdue instances, close vanished runs;
                                                               alarm and hard stop: see ALARM and HARD_STOP below

G10: a launch is refused when the money committed (closed sessions + every active run to its deadline) plus its own
worst case passes HARD_STOP, and needs --owner-approved (recorded in the ledger) when that projection reaches ALARM.
Costs are dph_total x time plus a bandwidth estimate (DOWN_GB once, UP_GB_PER_H per hour, at the offer's $/GB); teardown saves
the raw charges and records what Vast billed beside that estimate (`billed`, BK4-18).

Environment: VAST_HOME (ledger.csv, active_runs.json, offers/, .vast-allowlist; default: this directory), VAST_RUNS
(pull target, default <repo>/runs), DRY_RUN=1 (ssh / rsync are printed, not run; `vastai` is still called, so tests
put a fake one on PATH), VAST_POLL_S, VAST_LAUNCH_TIMEOUT_S. Secrets live in the macOS keychain and are read with
`security`: objdel-vast-scoped (scoped API key, passed to vastai as VAST_API_KEY), objdel-ghcr-read (read-only GHCR
PAT, Vast --login), objdel-riga-age (age identity: only ever the stdin of the ssh decrypt command; D16),
objdel-hf-write (HF token: checked by hf_push.token_problems, then only ever the stdin of the ssh upload command).
"""

import argparse
import base64
import csv
import gzip
import hashlib
import json
import math
import operator
import os
import re
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path

from objdel.core.spec import RunSpec, fill, load_spec, to_argv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = "/workspace/data/lidc"  # {data} on the box
# §4.2 / §4.8: gpu_ram bounds in GB (stock + 2, or >= 79 for the 80 GB A100), price ceilings = 1.5 x snapshot median.
CARDS = {
    "RTX_3090": {"ram": (None, 26), "ceiling": 0.50},
    "RTX_4090": {"ram": (None, 26), "ceiling": 1.00},
    "RTX_5090": {"ram": (None, 34), "ceiling": 1.50},
    "RTX_PRO_6000_WS": {"ram": (None, 98), "ceiling": 2.40, "cpu_ram": 48},
    "RTX_PRO_6000_S": {"ram": (None, 98), "ceiling": 2.40, "cpu_ram": 48},
    "A100_SXM4": {"ram": (79, 82), "ceiling": 1.80},
}
CPU_RAM = {"lidc": 32, "riga": 48}  # GB; RIGA caches ~4.5 GB RSS per scoring process (B4 §1)
ALARM, HARD_STOP = 10.0, 19.9  # G10, USD spent in Phase 2: $16.90 at plan start + the approved $3.00 (2026-11-09)
# BK4-18: $/GB ceiling, up and down (> 2x the A3 snapshot medians of the five cards, $4-13.65/TB; unit UNVERIFIED until
# V0). Estimate: image 7.5 GB + bundle 2 GB + model downloads once; checkpoint pulls of 3-5 GB per 20 min (§4.4).
INET_CEIL, DOWN_GB, UP_GB_PER_H = 0.03, 15, 6
FIELDS = tuple(
    "ts event run_id instance_id offer_id gpu dph_total inet_up_cost inet_down_cost hosting_type host_id machine_id "
    "bundle cost_usd detail".split()
)
SECRET_RE = re.compile(r"AGE-SECRET-KEY-|ghp_|github_pat_|PRIVATE KEY|Bearer ")
RUN_ID_RE = re.compile(r"[a-z0-9][a-z0-9_-]{0,39}")  # Vast label, shell, tmux name, local path (BK4-16); onstart size
BAD_STATUS = ("exited", "unknown", "offline")  # never reach running (V19)
IMPORT_CHECK = "import objdel._gpucheck, experiments.phase1.reward_gates, tests.oracles.reward_harness"
LIVE_CLIS = ("experiments.phase2.train_rl", "experiments.phase2.g_probe", "experiments.phase2.select_ckpt",
             "experiments.phase2.final_eval")  # the CLIs the box runs; test_bundle_imports_standalone imports them
# MARGIN_S: run.sh's timeout ends this long before the watchdog cap, leaving time for the last pull. GRACE_S: sweep
# never closes a run younger than this as vanished (it may not be listed yet). PULL_S: added to --cap-hours, as the
# watchdog counts from launch; V0's 5.88 GB image pull took 5.7 min.
MARGIN_S = GRACE_S = PULL_S = 600
OPS = {"=": operator.eq, ">": operator.gt, ">=": operator.ge, "<=": operator.le}
REGION = ("US", "CA")  # owner (2026-11-02): North American hosts only; V1: a Thailand host could not pull the image
# Machines that failed the GPU check on a paid box; never offered again
BAD_MACHINES = {28360: "2026-11-03 v2: GPU busy before our job, power limit 150 W, 7.4 of 165 bf16 TFLOPS"}
# Not pulled, and ignored by the watchdog's "nothing newer than PULLED_OK" (it gets this list in the onstart env):
# full checkpoints with optimizer state (R1: 29 GB; branches/ is model-only and is pulled) and the pull's own markers.
PULL_EXCLUDE = ("rolling/", "decay/", "stable/", "PULLED_OK", "PULLED_OK.next")
ONSTART_MAX = 4048  # Vast API spec


def home() -> Path:
    return Path(os.environ.get("VAST_HOME", HERE))


def die(msg: str):
    sys.exit(f"vast: {msg}")


def is_raw(p: Path) -> bool:
    return "/data/raw/" in f"{Path(p).resolve().as_posix()}/"


# --- offers: the §4.2 filter as data; `query` sends it, `reasons` re-checks it on the returned JSON -------------------


def terms(gpu: str, dataset: str, disk: int = 80) -> list[tuple]:
    """(reason, query term, offer JSON key, op, limit, scale) per term. Vast reports gpu_ram (MiB) and cpu_ram (MB) and
    the query takes GB; V0 showed cpu_ram>=32 returns 32160, so cpu_ram counts 1 GB = 1000 MB. `verification` is the
    string "verified" in the JSON (the query term is `verified=true`). `duration` is query-only: its unit is UNVERIFIED."""
    if gpu not in CARDS:
        raise ValueError(f"unknown card {gpu}; one of {sorted(CARDS)}")
    c, mb = CARDS[gpu], 1 / 1024
    (lo, hi), cpu = c["ram"], max(CPU_RAM[dataset], c.get("cpu_ram", 0))
    rows = [
        ("num_gpus", "num_gpus=1", "num_gpus", "=", 1, 1),
        ("gpu_name", f"gpu_name={gpu}", "gpu_name", "=", gpu, 1),
        ("gpu_ram", f"gpu_ram>={lo}", "gpu_ram", ">=", lo, mb) if lo else None,
        ("gpu_ram", f"gpu_ram<={hi}", "gpu_ram", "<=", hi, mb),
        ("verified", "verified=true", "verification", "=", "verified", 1),
        ("rentable", "rentable=true", "rentable", "=", True, 1),
        ("reliability", "reliability>0.98", "reliability2", ">", 0.98, 1),
        ("cores", "cpu_cores_effective>=12", "cpu_cores_effective", ">=", 12, 1),
        ("cpu_ram", f"cpu_ram>={cpu}", "cpu_ram", ">=", cpu, 1 / 1000),
        ("cuda", "cuda_vers>=13.0", "cuda_max_good", ">=", 13.0, 1),
        ("disk", f"disk_space>={disk}", "disk_space", ">=", disk, 1),
        ("inet_down", "inet_down>=500", "inet_down", ">=", 500, 1),
        ("ports", "direct_port_count>=2", "direct_port_count", ">=", 2, 1),
        ("duration", "duration>=1", None, None, None, None),
        ("price", f"dph<={c['ceiling']}", "dph_total", "<=", c["ceiling"], 1),
        ("bandwidth", f"inet_up_cost<={INET_CEIL}", "inet_up_cost", "<=", INET_CEIL, 1),
        ("bandwidth", f"inet_down_cost<={INET_CEIL}", "inet_down_cost", "<=", INET_CEIL, 1),
        ("datacenter", "datacenter=True", "hosting_type", "=", 1, 1) if dataset == "riga" else None,
        ("region", f"geolocation in [{','.join(REGION)}]", None, None, None, None),  # checked in reasons()
    ]
    return [r for r in rows if r]


def query(gpu: str, dataset: str, disk: int = 80) -> str:
    return " ".join(t[1] for t in terms(gpu, dataset, disk))


def reasons(o: dict, gpu: str, dataset: str, disk: int = 80) -> list[str]:
    """Why an offer JSON fails the §4.2 filter ([] = admissible); a missing or mistyped field fails."""
    bad = []
    for why, _, key, op, limit, scale in terms(gpu, dataset, disk):
        v = o.get(key)
        if key and (v is None or isinstance(v, str) != isinstance(limit, str) or not OPS[op](
                v.replace(" ", "_") if isinstance(v, str) else v * scale, limit)):  # fmt: skip
            bad.append(why)
    if (o.get("geolocation") or "").rsplit(", ", 1)[-1] not in REGION:  # "Texas, US"; none or unknown fails
        bad.append("region")
    if o.get("machine_id") in BAD_MACHINES:
        bad.append(f"machine {o['machine_id']} ({BAD_MACHINES[o['machine_id']]})")
    return [f"{r} fails the §4.2 filter" for r in dict.fromkeys(bad)]


def governance(dataset: str, o: dict, allow: dict) -> list[str]:
    """H-GOV (D16): LIDC anywhere verified; RIGA only with recorded written permission AND a datacenter offer."""
    if dataset != "riga":
        return [] if dataset == "lidc" else [f"unknown bundle {dataset}"]
    perm = allow.get("riga_permission") or {}
    why = [] if perm.get("granted") is True and perm.get("evidence") else ["no recorded RIGA permission"]
    if o.get("hosting_type") != 1 or o.get("datacenter", True) is not True:
        why.append(f"offer {o.get('id')} is not Secure Cloud (hosting_type={o.get('hosting_type')})")
    return why


# --- state: ledger (G9), active runs, money, secrets -----------------------------------------------------------------


def ledger(**row) -> None:
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
    return json.loads(f.read_text()) if f.exists() else {}


def active(update: dict | None = None) -> dict:
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


def secret(service: str) -> str:
    """objdel.secrets: <repo>/.secrets, else the macOS keychain; dies when neither has it."""
    from objdel import secrets

    try:
        value = secrets.get(service)
    except PermissionError as e:
        die(str(e))
    return value or die(f"secret {service} missing: not in {secrets.path()}, not in the keychain")


# --- commands ------------------------------------------------------------------------------------------------------


def vast(*args: str, tries: int = 1):
    env = dict(os.environ, VAST_API_KEY=os.environ.get("VAST_API_KEY") or secret("objdel-vast-scoped"))
    for i in range(tries):
        p = subprocess.run(["vastai", *args], capture_output=True, text=True, env=env, stdin=subprocess.DEVNULL)
        if not p.returncode and (p.stdout.strip() or "--raw" not in args):  # an API error exits 0 with empty stdout (V0)
            return json.loads(p.stdout) if "--raw" in args else p.stdout
        if i < tries - 1:
            time.sleep(float(os.environ.get("VAST_POLL_S", 10)))
    die(f"vastai {args[:2]} failed: {p.stderr.strip()[-300:]}")


def sh(cmd: list[str], stdin: str | None = None, what: str = "") -> None:
    """ssh / rsync: printed under DRY_RUN=1 (stdin shown only as a label), otherwise run; a failure stops."""
    if os.environ.get("DRY_RUN") == "1":
        print("DRY_RUN:", shlex.join(cmd), f"<stdin: {what}>" if stdin is not None else "")
    elif subprocess.run(cmd, input=stdin, text=True, stdin=None if stdin else subprocess.DEVNULL).returncode:
        die(f"{cmd[0]} failed: {shlex.join(cmd[-2:])[:160]}")


def ssh(run: dict) -> list[str]:
    known = f"UserKnownHostsFile={home() / 'known_hosts'}"
    return ["ssh", "-p", str(run["port"]), "-o", "StrictHostKeyChecking=accept-new", "-o", known]


def remote(run: dict, command: str, stdin: str | None = None, what: str = "") -> None:
    sh([*ssh(run), f"root@{run['ip']}", command], stdin, what)


def rsync(run: dict, *paths: str, opts: tuple = ()) -> None:
    """paths: sources, then the destination (V1: two sources once made the .sha256 the destination)."""
    sh(["rsync", "-az", "--exclude=data/raw/", *opts, "-e", shlex.join(ssh(run)), *paths])


def destroy(iid: int, run_id: str, event: str, inst: dict | None = None) -> None:
    """Destroy (never stop, G4), verify the instance is no longer listed, record the session's cost (Vast's dph and
    start when it lists the instance, else the run record's)."""
    inst = inst or next((i for i in vast("show", "instances", "--raw", tries=3) if i["id"] == iid), {})
    run = next((r for r in active().values() if r["instance"] == iid), {})
    r = dict(run, **({"dph": inst.get("dph_total", 0), "started": inst["start_date"]} if "start_date" in inst else {}))
    cost = accrued(r, time.time()) if r else 0.0
    vast("destroy", "instance", str(iid), "-y", tries=3)
    for _ in range(10):
        if all(i["id"] != iid for i in vast("show", "instances", "--raw", tries=3)):
            break
        time.sleep(float(os.environ.get("VAST_POLL_S", 10)))
    else:
        die(f"instance {iid} still listed after destroy: destroy it in the console NOW")
    active({k: x for k, x in active().items() if x["instance"] != iid})
    ledger(event=event, run_id=run_id, instance_id=iid, cost_usd=round(cost, 4))


def onstart_cmd(run_id: str, cap_s: int, created: float | None = None, idle_s: int = 1800) -> str:
    """--onstart-cmd: unpack the watchdog and start it; carries no secret (CONTAINER_API_KEY is injected by Vast).
    Vast's API spec limits onstart to ONSTART_MAX characters."""
    b64 = base64.b64encode(gzip.compress((HERE / "watchdog.sh").read_bytes(), mtime=0)).decode()  # gzip: Debian essential
    created = int(created or time.time())
    env = (f"RUN_ID={run_id} WATCHDOG_CAP_S={cap_s} WATCHDOG_IDLE_S={idle_s} WATCHDOG_CREATED={created} "
           f"WATCHDOG_PULLED_OK=/workspace/runs/{run_id}/PULLED_OK WATCHDOG_SKIP={','.join(PULL_EXCLUDE)}")  # fmt: skip
    k = "/root/.ssh/authorized_keys"  # V1: some hosts write it with a foreign owner and sshd refuses it (vast-cli #336)
    fix = (f"(for i in $(seq 60); do [ -f {k} ] && chown root:root {k} && chmod 600 {k} && chmod 700 /root/.ssh "
           f"&& chmod go-w /root; sleep 2; done) & ")
    tmux = "touch /root/.no_auto_tmux; "  # V1: Vast's login tmux broke every `ssh host CMD` ("not a terminal")
    cmd = tmux + fix + f"printf %s {b64} | base64 -d | gunzip > /root/watchdog.sh && {env} nohup bash /root/watchdog.sh >> /root/watchdog.log 2>&1 &"
    if len(cmd) > ONSTART_MAX:
        die(f"onstart is {len(cmd)} characters, over Vast's {ONSTART_MAX} (shorten watchdog.sh or the run_id)")
    return cmd


def offer_cols(o: dict) -> dict:
    return {k: o.get(k) for k in ("dph_total", "inet_up_cost", "inet_down_cost", "hosting_type", "host_id", "machine_id")}


def cmd_offers(a) -> None:
    ds = "riga" if a.riga else "lidc"
    got = vast("search", "offers", query(a.gpu, ds, a.disk), "-o", "dph_total", "--storage", str(a.disk), "--raw")
    ok = [o for o in got if not reasons(o, a.gpu, ds, a.disk)]  # --pick ID: any of them, not only the 5 shown
    for o in ok[:5]:
        print(f"{o['id']:>10}  ${o['dph_total']:.3f}/h  {o['gpu_ram'] / 1024:.0f} GB  {o['cpu_cores_effective']} cores"
              f"  {o['cpu_ram'] / 1000:.0f} GB RAM  cuda {o['cuda_max_good']}  hosting_type {o.get('hosting_type')}"
              f"  ${o['inet_up_cost']}/GB up")  # fmt: skip
    if a.pick is None:
        return
    o = ok[0] if a.pick == "top" and ok else next((o for o in ok if str(o["id"]) == a.pick), None)  # top: no race
    if o is None:
        die(f"offer {a.pick} is not among the admissible offers")
    (home() / "offers").mkdir(exist_ok=True)
    (home() / "offers" / f"{o['id']}.json").write_text(json.dumps(dict(o, _gpu=a.gpu, _dataset=ds, _disk=a.disk)))
    ledger(event="offer", offer_id=o["id"], gpu=a.gpu, **offer_cols(o), detail=query(a.gpu, ds, a.disk))
    print(f"picked {o['id']}")


def launch_refusal(offer: dict, cap_hours: float) -> str:
    """Why this offer may not launch (empty when it may): committed spend + the worst case to the deadline > HARD_STOP."""
    money = [offer["dph_total"], offer["inet_up_cost"], offer["inet_down_cost"]]
    if not all(math.isfinite(x) and x >= 0 for x in money) or not (math.isfinite(cap_hours) and cap_hours > 0):
        return f"refused: price {money} and --cap-hours {cap_hours} must be finite, non-negative and the cap positive"
    now = time.time()
    bw = {k: offer[k] for k in ("inet_up_cost", "inet_down_cost")}
    run = {"dph": offer["dph_total"], **bw, "started": now, "deadline": now + int(cap_hours * 3600) + PULL_S}
    committed, worst = spend(committed=True), accrued(run, run["deadline"])
    if not committed + worst <= HARD_STOP:  # not "> ": a NaN price must refuse too
        return f"hard stop: ${committed:.2f} committed + ${worst:.2f} worst case > ${HARD_STOP} (G10)"
    return ""


def cmd_launch(a) -> None:
    if not RUN_ID_RE.fullmatch(a.run_id):
        die(f"run_id must match {RUN_ID_RE.pattern}")
    if a.run_id in active():
        die(f"run {a.run_id} is still active (teardown.sh first)")
    from experiments.phase2 import preflight

    if why := preflight.refusal(Path(a.config)):  # the mandatory Mac preflight (make preflight): no pass, no rental
        die(f"preflight gate: {why}")
    if not (o := load(f"offers/{a.offer}.json")):
        die(f"offer {a.offer} was not picked with offers.sh --pick")
    if why := reasons(o, o["_gpu"], o["_dataset"], o["_disk"]):
        die("; ".join(why))
    if why := launch_refusal(o, a.cap_hours):
        die(why)
    cap_s, now = int(a.cap_hours * 3600) + PULL_S, time.time()
    bw = {k: o[k] for k in ("inet_up_cost", "inet_down_cost")}
    run = {"instance": None, "offer": a.offer, "dph": o["dph_total"], **bw, "started": now, "deadline": now + cap_s}
    committed, worst = spend(committed=True), accrued(run, run["deadline"])
    if committed + worst >= ALARM and not a.owner_approved:
        die(f"alarm: ${committed:.2f} committed + ${worst:.2f} worst case >= ${ALARM}; pass --owner-approved only "
            "with the owner's written approval (G10, §4.6)")  # fmt: skip
    if (live := vast("show", "instances", "--raw")) and not a.allow_parallel:
        die(f"{len(live)} instance(s) already exist; one instance at a time unless --allow-parallel")
    allow = load(".vast-allowlist")
    image = allow.get("image") or die("no image in .vast-allowlist")
    login = f"-u {allow.get('ghcr_user', '')} -p {secret('objdel-ghcr-read')} ghcr.io"
    flags = ["--ssh", "--direct", "--cancel-unavail", "--label", a.run_id, "--raw"]
    sigs = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
    onstart = onstart_cmd(a.run_id, cap_s, now, int(a.idle_minutes * 60))
    for s in sigs:  # an interrupt from `create` on must reach the destroy below (BK4-3)
        signal.signal(s, signal.default_int_handler)
    iid = None
    try:  # G8: any failure or interrupt from here on destroys the instance
        out = vast("create", "instance", str(a.offer), "--image", image, "--login", login, "--disk", str(o["_disk"]),
                   "--onstart-cmd", onstart, *flags)  # fmt: skip
        run["instance"] = iid = int(out["new_contract"])
        active({**active(), a.run_id: run})
        ledger(event="launch", run_id=a.run_id, instance_id=iid, offer_id=a.offer, gpu=o["_gpu"], **offer_cols(o),
               detail=f"cap_s={cap_s} image={image}" + " owner_approved=yes" * a.owner_approved)  # fmt: skip
        t0, timeout = time.time(), float(os.environ.get("VAST_LAUNCH_TIMEOUT_S", 900))
        while True:
            inst = vast("show", "instance", str(iid), "--raw", tries=3)
            status, port = inst.get("actual_status"), ((inst.get("ports") or {}).get("22/tcp") or [{}])[0]
            if status == "running" and inst.get("public_ipaddr") and port.get("HostPort"):  # BK4-4: ports come late
                run.update(ip=inst["public_ipaddr"], port=int(port["HostPort"]))
                active({**active(), a.run_id: run})
                print(f"running: ssh -p {run['port']} root@{run['ip']}")
                return
            if status in BAD_STATUS or time.time() - t0 > timeout:
                die(f"instance {iid} {status if status in BAD_STATUS else 'timeout'} before running; destroyed")
            time.sleep(float(os.environ.get("VAST_POLL_S", 10)))
    except BaseException:
        for s in sigs:
            signal.signal(s, signal.SIG_IGN)
        # create failed or its reply was lost: a box may exist all the same, found by its label
        ids = [iid] if iid else [i["id"] for i in vast("show", "instances", "--raw", tries=3) if i.get("label") == a.run_id]
        for i in ids:
            destroy(i, a.run_id, "destroy")
        raise


def get_run(run_id: str) -> dict:
    return active().get(run_id) or die(f"no active run {run_id}")


def code_files() -> list[str]:
    """The code snapshot (tracked and new files), without the roots data/, docs/, archive/ and .claude/, lab or results
    dirs, or RIGA-derived tables (D16, BK4-13); only files on disk (rsync --files-from aborts on a missing one).
    tests/ ships: k3_harness -> phase1.reward_gates -> tests.oracles. Kept from docs/: the two records the box reads
    (final_eval's frozen view hashes, val_half's table; LIDC ids, counts and hashes, no pixels, no RIGA; B1.6) and the
    kernel drafts K1-K3 without their build/ and results/ folders (the R-D1 runbooks run them on the box)."""
    git = ["git", "ls-files", "-z", "-co", "--exclude-standard"]  # -z: paths with spaces stay whole (V1)
    files = subprocess.run(git, cwd=ROOT, capture_output=True, text=True, check=True).stdout.split("\0")[:-1]
    riga = re.compile(r"riga[^/]*\.(csv|jsonl?)$", re.I)
    drop = re.compile(r"^(data|docs|archive|\.claude)/|^experiments/[^/]+/(lab|results)/", re.I)
    keep = re.compile(r"^docs/reference/(frozen_views|val_halves)\.json$"
                      r"|^docs/implementation2\.0/kernels/K[123]/(?!build/|results/)")
    files = [f for f in dict.fromkeys(files) if not riga.search(f) and (keep.match(f) or not drop.search(f)) and (ROOT / f).is_file()]
    if not any(f.startswith("tests/oracles/") for f in files):
        die("tests/oracles missing from the code snapshot")
    return files


def code_sha_text(root: Path = ROOT) -> str:
    """The box's CODE_SHA (final_eval.code_sha): git rev-parse HEAD, plus -dirty-<first 12 hex of sha256(D)> when git
    status --porcelain is not empty; D = git diff HEAD, then the sorted untracked, not-ignored paths, one per line
    (a new module marks the tree dirty). A git failure raises: never an empty or wrong CODE_SHA."""
    def git(*args: str) -> bytes:
        return subprocess.run(["git", *args], cwd=root, capture_output=True, check=True).stdout

    head = git("rev-parse", "HEAD").decode().strip()
    if not git("status", "--porcelain"):
        return head
    new = sorted(git("ls-files", "-z", "--others", "--exclude-standard").split(b"\0")[:-1])
    d = git("diff", "HEAD") + b"".join(p + b"\n" for p in new)
    return f"{head}-dirty-{hashlib.sha256(d).hexdigest()[:12]}"


def cmd_sync_up(a) -> None:
    run = get_run(a.run_id)
    o, allow = json.loads((home() / "offers" / f"{run['offer']}.json").read_text()), load(".vast-allowlist")
    if why := governance(a.bundle, o, allow):
        die("H-GOV refuses: " + "; ".join(why))
    name = "bundle-riga.tar.age" if a.bundle == "riga" else f"bundle-{a.bundle}.tar"
    files = [a.bundle_dir / name, a.bundle_dir / f"{name}.sha256"]
    if any(map(is_raw, [a.bundle_dir, *files])):
        die("refusing to send anything under data/raw")
    lst, sha = home() / "code_files.txt", code_sha_text()
    lst.write_text("\n".join(code_files()) + "\n")
    host = f"root@{run['ip']}:"
    rsync(run, f"{ROOT}/", f"{host}/workspace/src/", opts=(f"--files-from={lst}",))
    remote(run, f"printf '%s\\n' {shlex.quote(sha)} > /workspace/src/CODE_SHA")  # the box has no .git (A3)
    remote(run, "cd /workspace/src && uv pip install --python /opt/venv/bin/python --no-deps --no-build-isolation -e ."
                f" && cd / && /opt/venv/bin/python -I -c {shlex.quote(IMPORT_CHECK)}")  # fmt: skip
    if a.code_only:  # a code fix mid-session: the bundle is already on the box
        return ledger(event="sync", run_id=a.run_id, instance_id=run["instance"], detail="code only")
    # GPU check before the bundle (2026-11-03 v2: a busy, power-capped card was found only after a 12 min upload)
    check = f"bash experiments/phase2/docker/image_check.sh /workspace/runs/{a.run_id}/image_check.json"
    remote(run, f"cd /workspace/src && {check}")
    rsync(run, *map(str, files), f"{host}/workspace/data/")
    verify = f"cd /workspace/data && sha256sum -c {name}.sha256"
    unpack = "mkdir -p {d} && {x} -C {d} && cd {d} && sha256sum -c --quiet SHA256SUMS"
    if a.bundle == "riga":  # D16: identity streamed into age over this session's stdin; plaintext only in /dev/shm
        # BK4-19: Docker's /dev/shm defaults to 64 MB and V0 has not shown Vast's; refuse before decrypting
        need = -(-(a.bundle_dir / name).stat().st_size * 11 // 10240)  # KB, ciphertext + 10 %
        shm = f"avail=$(df -Pk /dev/shm | awk 'NR==2 {{print $4}}') && [ \"$avail\" -ge {need} ]"
        unpack = unpack.format(d="/dev/shm/riga", x=f"age -d -i /dev/stdin {name} | tar -x")
        remote(run, f"{verify} && {shm} && {unpack}", secret("objdel-riga-age") + "\n", "age identity from the keychain")
        detail = f"recipient={allow.get('riga_recipient', '')} shm=UNVERIFIED until V0 (BK4-19)"
    else:
        remote(run, f"{verify} && " + unpack.format(d=a.bundle, x=f"tar -xf {name}"))
        detail = ""
    ledger(event="sync", run_id=a.run_id, instance_id=run["instance"], bundle=name, detail=detail)


def chain(cfg: Path, run_id: str) -> list[list[str]]:
    """A config's stages on the box, as argv lists (`-m module ...` runs in /opt/venv, `bash script ...` as is).
    `steps:` (R0) lists them literally with {rd} / {data} filled in; a leading "~" lets a step fail without stopping
    the chain. Otherwise `stages:` (default train, select, probe; C2: train, topology) picks from train_sft, then on
    the saved branches select_ckpt (G-vision block, rotated view, `sampler:`) or its topology-only form, and g_probe;
    a config without branches (timing arms) only trains. A chain that scores (select, topology or probe) first builds
    the CUDA checker (run_box.sh --build-only) and checks the scorers on it against the host scorer (tests/test_gpu_scoring.py),
    so a failed build or a difference stops the chain before training: scoring runs on the GPU only. Box
    paths: data /workspace/data/lidc, runs /workspace/runs/<id>."""
    import yaml

    c, data, rd = yaml.safe_load(Path(cfg).read_text()), DATA, f"/workspace/runs/{run_id}"
    if "schema_version" in c:  # a run spec (specs/*.yaml)
        return [to_argv(s, rd, data) for s in load_spec(Path(cfg)).stages]
    if "steps" in c:
        return [[str(x).format(rd=rd, data=data) for x in s] for s in c["steps"]]
    run, ev, prep = f"{rd}/{c['output_dir']}", f"{data}/eval_views", f"{data}/prepared"
    v, samp = c.get("vision"), ["--sampler", c["sampler"]] if c.get("sampler") else []
    items = ["--items", "256"] if c.get("sampler") == "hf" else []  # D13: HF path selects on 256 shifted-val items
    ref = v and (f"{rd}/{v['reference']}" if isinstance(v["reference"], str) else v["reference"])  # a number as is
    vis = ["--vision", json.dumps({"prepared": f"{data}/{v['prepared']}", "reference": ref})] if v else ["--no-vision"]
    select = ["-m", "experiments.phase2.select_ckpt", "--run", run, "--prepared", prep, "--probes", ev, *samp]
    table = {
        "train": ["-m", "experiments.phase2.train_sft", "--config", Path(cfg).resolve().relative_to(ROOT).as_posix()],
        "select": [*select, "--shifted", f"{ev}/lidc-val-shifted", "--rescaled", f"{ev}/lidc-val-rescaled",
                   "--rotated", f"{ev}/lidc-val-rotated", *items, *vis],
        "topology": [*select, "--only", "topology"],
        "probe": ["-m", "experiments.phase2.g_probe", "--model", c["model"], "--select", f"{run}/select", "--train",
                  prep, "--probe", f"{ev}/lidc-val-shifted", "--out", f"{run}/g_probe"],
    }  # fmt: skip
    stages = c.get("stages", ["train", "select", "probe"]) if c.get("branches") or c.get("branches_from") else ["train"]
    gpu = [["bash", "-c", "bash experiments/cuda_checker/run_box.sh --build-only && /opt/venv/bin/python -m pytest -q "
                          "-p no:cacheprovider tests/test_gpu_scoring.py"]] if set(stages) - {"train"} else []
    return gpu + [table[s] for s in stages]


def vision_refs(cfg: Path, run_id: str) -> list[str]:
    """The G-vision reference JSON a config reads on the box ([] for a number or `vision: off`)."""
    import yaml

    c = yaml.safe_load(Path(cfg).read_text())
    if "schema_version" in c:  # a run spec: the reference of each sft and select stage, in stage order, once
        rd, refs = f"/workspace/runs/{run_id}", []
        for s in load_spec(Path(cfg)).stages:
            if s.kind == "sft" and s.vision and isinstance(s.vision["reference"], str):
                refs.append(f"{rd}/{s.vision['reference']}")
            elif s.kind == "select" and s.vision and isinstance(s.vision.reference, str):
                refs.append(fill(s.vision.reference, rd, DATA))
        return list(dict.fromkeys(refs))
    ref = (c.get("vision") or {}).get("reference")
    return [f"/workspace/runs/{run_id}/{ref}"] if isinstance(ref, str) else []


def cmd_run(a) -> None:
    """MODULE [ARGS] runs one module (experiments.phase2.*, or R0's xarch_hashes); a configs/*.yaml path runs its
    `chain` under one cap, refused on the box while a G-vision reference it reads is missing (R0 writes them). Each job
    has its own tmux session and log (<run_id>-<name>, <name>.log), so a second job can run beside the first."""
    run, guard = get_run(a.run_id), []
    if a.module.endswith(".yaml") and a.args:  # V1: a trailing --cap-hours was silently ignored
        die(f"flags go before the config (run.sh RUN_ID [--cap-hours H] CONFIG); got {' '.join(a.args)}")
    if a.module.endswith(".yaml"):
        stages, name = chain(ROOT / a.module, a.run_id), Path(a.module).stem
        guard = [f"{{ test -f {r} || {{ echo 'missing G-vision reference {r}: run.sh {a.run_id} "
                 f"experiments/phase2/configs/r0.yaml first' >&2; false; }}; }}"
                 for r in vision_refs(ROOT / a.module, a.run_id)]  # fmt: skip
    elif re.fullmatch(r"experiments\.phase2\.\w+|experiments\.phase1\.bench\.xarch_hashes", a.module):
        stages, name = [["-m", a.module, *a.args]], a.module.rsplit(".", 1)[1]
    else:
        die("module must be experiments.phase2.<script> or experiments/phase2/configs/<name>.yaml")
    cap = int(min(a.cap_hours * 3600, run["deadline"] - time.time() - MARGIN_S))
    if cap <= 0:
        die("no time left before the instance cap")

    def job(s: list[str]) -> str:  # "~": a failure is logged and the chain goes on
        soft = s[0] == "~"
        cmd = shlex.join(s[soft:]) if s[soft] == "bash" else f"/opt/venv/bin/python -I {shlex.join(s[soft:])}"
        return f"{{ {cmd} || echo {shlex.quote('step failed, continuing: ' + cmd)}; }}" if soft else cmd

    rd = f"/workspace/runs/{a.run_id}"
    jobs = " && ".join([*guard, *map(job, stages)])
    # the session clock of train_rl --budget-usd (callbacks.BudgetStop): accrued(run, t) = fixed + rate x hours
    rate, fixed = run["dph"] + run.get("inet_up_cost", 0) * UP_GB_PER_H, run.get("inet_down_cost", 0) * DOWN_GB
    clock = f"OBJDEL_SESSION_START={run['started']!r} OBJDEL_RATE_USD_PER_H={rate!r} OBJDEL_FIXED_USD={fixed!r}"
    # EXIT is removed first: the watchdog's idle destroy needs it (BK4-11); `;` keeps cd/export in this shell
    inner = (f"cd /workspace/src && export OBJDEL_RUN_DIR={rd} {clock} && mkdir -p {rd} && rm -f {rd}/EXIT "
             "|| exit 1; "
             f"(while :; do date -u +%s > {rd}/HEARTBEAT; sleep 60; done) & "
             f"timeout -k 120 {cap}s bash -c {shlex.quote(jobs)} > {rd}/{name}.log 2>&1; echo $? > {rd}/EXIT; kill %1")  # fmt: skip
    remote(run, f"tmux new-session -d -s {a.run_id}-{name} {shlex.quote('bash -c ' + shlex.quote(inner))}")
    ledger(event="run", run_id=a.run_id, instance_id=run["instance"], detail=f"{a.module} cap_s={cap}")


def cmd_pull(a) -> None:
    """rsync verifies every transferred file with a whole-file checksum; only then does PULLED_OK appear, carrying
    the pull's start time (PULLED_OK.next, touched before rsync), so a file written during the pull is newer."""
    for run_id in [a.run_id] if a.run_id else list(active()):
        run, dst = get_run(run_id), Path(os.environ.get("VAST_RUNS", ROOT / "runs")) / run_id
        dst.mkdir(parents=True, exist_ok=True)
        rd = f"/workspace/runs/{run_id}"
        remote(run, f"touch {rd}/PULLED_OK.next")
        rsync(run, f"root@{run['ip']}:{rd}/", f"{dst}/", opts=tuple(f"--exclude={x}" for x in PULL_EXCLUDE))
        remote(run, f"mv {rd}/PULLED_OK.next {rd}/PULLED_OK")


def push_target(spec: RunSpec, run_id: str) -> tuple[str, str, bool]:
    """(folder on the box, repo revision, is an rl output) of a spec's last rl or sft stage. For rl the revision is the
    last part of its out folder; for sft the output_dir."""
    last = next((s for s in reversed(spec.stages) if s.kind in ("rl", "sft")), None)
    if last is None:
        die("nothing to push: the spec has no rl or sft stage")
    rd = f"/workspace/runs/{run_id}"
    if last.kind == "rl":
        path = fill(last.out, rd, DATA)
        return path, path.rsplit("/", 1)[1], True
    return f"{rd}/{last.output_dir}", last.output_dir, False


def cmd_push_hf(a) -> None:
    """hf_push on the box for CONFIG's run (CONFIG is a run spec, or an SFT launch config: repo branch = its
    output_dir; a config with rows_plan names a train_rl output: --rl); the token is checked here first."""
    import yaml

    from experiments.phase2 import hf_push

    run, c = get_run(a.run_id), yaml.safe_load((ROOT / a.config).read_text())
    if "schema_version" in c:
        path, out, rl = push_target(load_spec(ROOT / a.config), a.run_id)
        kind = "--rl" if rl else "--run"
    else:
        out, kind = c["output_dir"], "--rl" if c.get("rows_plan") else "--run"
        path = f"/workspace/runs/{a.run_id}/{out}"
    token = secret("objdel-hf-write")
    if bad := hf_push.token_problems(hf_push.whoami(token), hf_push.REPO):
        die(f"secret objdel-hf-write (.secrets or keychain) is too broad for a rented box: {'; '.join(bad)}. Recreate it at "
            f"huggingface.co/settings/tokens as fine-grained with only {hf_push.REPO}: write (no other boxes ticked)")
    cmd = (f"cd /workspace/src && /opt/venv/bin/python -I -m experiments.phase2.hf_push {kind} {path} "
           f"--repo {hf_push.REPO} --revision {out}")  # fmt: skip
    remote(run, cmd, token + "\n", "HF token from objdel.secrets")


ROWS_CHECK = ("import hashlib, json, pathlib, sys; d, want = pathlib.Path(sys.argv[1]), sys.argv[2]; "
              "m = d / 'manifest.json'; ok = m.is_file() and hashlib.sha256(m.read_bytes()).hexdigest() == want "
              "and sum(p.is_file() for p in d.rglob('*')) - 1 == json.loads(m.read_text())['schedule']['files_count']; "
              "sys.exit(0 if ok else f'rows root {d}: no manifest.json, not SHA-256 {want}, or not files_count files')")


def rows_pin(revision: str) -> str:
    """The --rows-manifest-sha256 that the plans in experiments/phase2/configs and specs/ pin beside REVISION."""
    import yaml

    plans = sorted((ROOT / "experiments/phase2/configs").glob("*.yaml"))
    pins = {str(s[s.index("--rows-manifest-sha256") + 1]) for p in plans
            for s in (yaml.safe_load(p.read_text()) or {}).get("steps", [])
            if "--rows-revision" in s and str(s[s.index("--rows-revision") + 1]) == revision}
    pins |= {s.rows_manifest_sha256 for p in sorted((ROOT / "specs").glob("*.yaml")) for s in load_spec(p).stages
             if s.kind == "rl" and s.rows_revision == revision}
    if len(pins) != 1:
        die(f"--rows {revision}: the plans pin {len(pins)} manifest SHA-256 values for it (need exactly one)")
    return pins.pop()


def cmd_fetch_hf(a) -> None:
    """hf_pull on the box: REVISION of the private HF repo to /workspace/runs/RUN_ID/NAME (--as; default ckpt, RL's
    start checkpoint); or --rows REVISION ROOT: the repo folder rows/ROOT/ to /workspace/data/lidc/rows/ROOT, then
    the same ssh call checks it there (manifest.json with the SHA-256 the plans pin for REVISION, files_count files)."""
    from experiments.phase2 import hf_push

    if (a.revision is None) == (a.rows is None):
        die("fetch_hf: give REVISION or --rows REVISION ROOT")
    rev, name = a.rows or (a.revision, a.name)
    if "/" in name or ".." in name or name in ("", "."):
        die(f"fetch_hf: folder name {name!r} may not hold / or ..")
    pull = (f"cd /workspace/src && /opt/venv/bin/python -I -m experiments.phase2.hf_pull --repo {hf_push.REPO} "
            f"--revision {shlex.quote(rev)}")
    if a.rows:
        d = f"/workspace/data/lidc/rows/{name}"
        cmd = (f"{pull} --out /workspace/data/lidc --allow {shlex.quote(f'rows/{name}/**')} && /opt/venv/bin/python "
               f"-I -c {shlex.quote(ROWS_CHECK)} {shlex.quote(d)} {rows_pin(rev)}")
    else:
        cmd = f"{pull} --out /workspace/runs/{a.run_id}/{name}"
    run, token = get_run(a.run_id), secret("objdel-hf-write")
    if bad := hf_push.token_problems(hf_push.whoami(token), hf_push.REPO):
        die(f"secret objdel-hf-write is too broad for a rented box: {'; '.join(bad)}")
    remote(run, cmd, token + "\n", "HF token from objdel.secrets")


def billed(charges: dict, iid: int) -> float | None:
    """USD Vast charged instance `iid` (V0: `invoices-v1 --charges` rows carry it only as `source` "instance-<id>")."""
    rows = [r["amount"] for r in charges.get("results", []) if r.get("source") == f"instance-{iid}"]
    return round(sum(rows), 4) if rows else None


def cmd_teardown(a) -> None:
    run, lost = get_run(a.run_id), None
    try:
        for _ in range(3 if "port" in run else 0):  # no port: the box never got ssh, so nothing is on it to pull
            try:
                cmd_pull(a)
                lost = None
                break
            except SystemExit as e:  # V1: ssh refused, so the pull failed; the box must still go (G8)
                lost = e
    finally:  # Ctrl-C included
        destroy(run["instance"], a.run_id, "teardown")
    inv = vast("show", "invoices-v1", "--charges", "--raw", "--limit", "100", tries=3)
    (home() / "invoices").mkdir(exist_ok=True)
    (home() / "invoices" / f"{a.run_id}.json").write_text(json.dumps(inv, indent=1))
    usd = billed(inv, run["instance"])  # in detail, not cost_usd: spend() already counts the teardown row
    usd = "pending" if usd is None else f"${usd:.3f}"
    ledger(event="invoice", run_id=a.run_id, instance_id=run["instance"],
           detail=f"billed={usd} (Vast charges accrue late; invoices/{a.run_id}.json)")
    if lost:
        die(f"final pull failed ({lost}); instance {run['instance']} destroyed")


def alert(msg: str) -> None:
    print(f"vast ALERT: {msg}", file=sys.stderr)
    if sys.platform == "darwin" and os.environ.get("DRY_RUN") != "1":
        subprocess.run(["osascript", "-e", f'display notification {json.dumps(msg)} with title "objdel vast"'])


def cmd_sweep_orphans(a) -> None:
    listed, now = vast("show", "instances", "--raw", tries=3), time.time()
    for k, r in active().items():  # BK4-7: gone without our destroy (watchdog G5, host loss): close at its last state
        if r["instance"] not in {i["id"] for i in listed} and now - r["started"] > GRACE_S:
            alert(f"run {k}: instance {r['instance']} no longer listed; closing it")
            active({j: x for j, x in active().items() if j != k})
            ledger(event="vanished", run_id=k, instance_id=r["instance"], detail="cost to min(now, deadline)",
                   cost_usd=round(accrued(r, min(now, r["deadline"])), 4))  # fmt: skip
    runs = {r["instance"]: (k, r) for k, r in active().items()}
    for inst in listed:
        k, r = runs.get(inst["id"], (inst.get("label") or "?", None))
        if r is None or now > r["deadline"]:
            alert(f"destroying {'orphan' if r is None else 'overdue'} instance {inst['id']} ({k})")
            destroy(inst["id"], k, "orphan" if r is None else "overdue", inst)
    usd = spend()
    if usd >= HARD_STOP:
        alert(f"hard stop: ${usd:.2f} spent >= ${HARD_STOP}; destroying every instance")
        for inst in vast("show", "instances", "--raw", tries=3):
            destroy(inst["id"], runs.get(inst["id"], (inst.get("label", "?"),))[0], "hard_stop", inst)
    elif usd >= ALARM:
        alert(f"alarm: ${usd:.2f} spent >= ${ALARM} (G10)")


FLAG = {"action": "store_true"}
CLI = {  # subcommand -> its arguments, in order
    "offers": {"gpu": {}, "--riga": FLAG, "--pick": {}, "--disk": {"type": int, "default": 80}},
    "launch": {"offer": {"type": int}, "run_id": {}, "--config": {"required": True, "help": "the run's config.yaml"},
               "--cap-hours": {"type": float, "default": 4},
               "--allow-parallel": FLAG, "--owner-approved": FLAG,
               "--idle-minutes": {"type": float, "default": 30, "help": "watchdog idle destroy after a pulled EXIT"}},
    "sync_up": {"run_id": {}, "--bundle": {"default": "lidc", "choices": ["lidc", "riga"]},
                "--bundle-dir": {"type": Path, "default": ROOT / "build" / "bundles"}, "--code-only": FLAG},
    "run": {"run_id": {}, "module": {}, "--cap-hours": {"type": float, "default": 24},
            "args": {"nargs": argparse.REMAINDER}},
    "pull": {"run_id": {"nargs": "?"}},
    "push_hf": {"run_id": {}, "config": {}},
    "fetch_hf": {"run_id": {}, "revision": {"nargs": "?"}, "--rows": {"nargs": 2, "metavar": ("REVISION", "ROOT")},
                 "--as": {"dest": "name", "default": "ckpt", "help": "folder under the run dir (default ckpt)"}},
    "teardown": {"run_id": {}},
    "sweep_orphans": {},
}  # fmt: skip


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="vast", description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for cmd, args in CLI.items():
        p = sub.add_parser(cmd)
        for name, kw in args.items():
            p.add_argument(name, **kw)
    a = ap.parse_args(argv)
    globals()[f"cmd_{a.cmd}"](a)


if __name__ == "__main__":
    main()
