"""Vast AI operations for Phase 2 (PLAN §4, D16, G1-G10). Every decision and every command is built here; the
`*.sh` files next to this module are symlinks to the one shim `vast.sh` (the link's name is the command), and
`watchdog.sh` is the only logic that runs on the box.

    offers.sh GPU [--riga] [--disk GB] [--pick OFFER_ID]       search (§4.2 filter), top 5; --pick saves the offer
    launch.sh OFFER_ID RUN_ID [--cap-hours H] [--allow-parallel] [--owner-approved]
    sync_up.sh RUN_ID [--bundle lidc|riga] [--bundle-dir DIR]  code (git ls-files) + bundle, SHA-256 verified on the box
    run.sh RUN_ID MODULE|CONFIG.yaml|vision-ref CONFIG... [ARGS]  timeout + tmux + heartbeat; a config runs its chain
    pull.sh [RUN_ID]                                           rsync runs/<id>/ to the Mac, touch PULLED_OK (no id: all)
    teardown.sh RUN_ID                                         final pull, destroy, verify absent, invoices, ledger
    sweep_orphans.sh                                           destroy unknown / overdue instances, close vanished runs;
                                                               $10 alarm, $15 stop

G10: a launch is refused when the money committed (closed sessions + every active run to its deadline) plus its own
worst case passes $15, and needs --owner-approved (recorded in the ledger) when that projection reaches $10. Costs are
dph_total x time plus a bandwidth estimate (DOWN_GB once, UP_GB_PER_H per hour, at the offer's $/GB); invoices are
saved raw and flagged UNPARSED until V0 shows their format (BK4-18).

Environment: VAST_HOME (ledger.csv, active_runs.json, offers/, .vast-allowlist; default: this directory), VAST_RUNS
(pull target, default <repo>/runs), DRY_RUN=1 (ssh / rsync are printed, not run; `vastai` is still called, so tests
put a fake one on PATH), VAST_POLL_S, VAST_LAUNCH_TIMEOUT_S. Secrets live in the macOS keychain and are read with
`security`: objdel-vast-scoped (scoped API key, passed to vastai as VAST_API_KEY), objdel-ghcr-read (read-only GHCR
PAT, Vast --login), objdel-riga-age (age identity: only ever the stdin of the ssh decrypt command; D16).
"""

import argparse
import base64
import csv
import functools
import json
import operator
import os
import re
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
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
ALARM, HARD_STOP = 10.0, 15.0  # G10, USD spent in Phase 2
# BK4-18: $/GB ceiling, up and down (> 2x the A3 snapshot medians of the five cards, $4-13.65/TB; unit UNVERIFIED until
# V0). Estimate: image 7.5 GB + bundle 2 GB + model downloads once; checkpoint pulls of 3-5 GB per 20 min (§4.4).
INET_CEIL, DOWN_GB, UP_GB_PER_H = 0.03, 15, 6
FIELDS = tuple(
    "ts event run_id instance_id offer_id gpu dph_total inet_up_cost inet_down_cost hosting_type host_id machine_id "
    "bundle cost_usd detail".split()
)
SECRET_RE = re.compile(r"AGE-SECRET-KEY-|ghp_|github_pat_|PRIVATE KEY|Bearer ")
RUN_ID_RE = re.compile(r"[a-z0-9][a-z0-9_-]*")  # Vast label, shell text and tmux name on the box, local path (BK4-16)
BAD_STATUS = ("exited", "unknown", "offline")  # never reach running (V19)
IMPORT_CHECK = "import objdel._polycore, experiments.phase1.reward_gates, tests.oracles.reward_harness"
MARGIN_S = 600  # run.sh's timeout ends this long before the watchdog cap, leaving time for the last pull
GRACE_S = 600  # sweep: a run younger than this may not be listed yet, so it is never closed as vanished
OPS = {"=": operator.eq, ">": operator.gt, ">=": operator.ge, "<=": operator.le}


def home() -> Path:
    return Path(os.environ.get("VAST_HOME", HERE))


def dry() -> bool:
    return os.environ.get("DRY_RUN") == "1"


def die(msg: str):
    sys.exit(f"vast: {msg}")


def is_raw(p: Path) -> bool:
    parts = Path(p).resolve().parts
    return any(a == "data" and b == "raw" for a, b in zip(parts, parts[1:], strict=False))


# --- offers: the §4.2 filter as data; `query` sends it, `reasons` re-checks it on the returned JSON -------------------


def terms(gpu: str, dataset: str, disk: int = 80) -> list[tuple]:
    """(reason, query term, offer JSON key, op, limit, scale) per term. Vast reports gpu_ram / cpu_ram in MB (the query
    takes GB). `duration` is query-only: its JSON unit is UNVERIFIED."""
    if gpu not in CARDS:
        raise ValueError(f"unknown card {gpu}; one of {sorted(CARDS)}")
    c, mb = CARDS[gpu], 1 / 1024
    (lo, hi), cpu = c["ram"], max(CPU_RAM[dataset], c.get("cpu_ram", 0))
    rows = [
        ("num_gpus", "num_gpus=1", "num_gpus", "=", 1, 1),
        ("gpu_name", f"gpu_name={gpu}", "gpu_name", "=", gpu, 1),
        ("gpu_ram", f"gpu_ram>={lo}", "gpu_ram", ">=", lo, mb) if lo else None,
        ("gpu_ram", f"gpu_ram<={hi}", "gpu_ram", "<=", hi, mb),
        ("verified", "verified=true", "verified", "=", True, 1),
        ("rentable", "rentable=true", "rentable", "=", True, 1),
        ("reliability", "reliability>0.98", "reliability2", ">", 0.98, 1),
        ("cores", "cpu_cores_effective>=12", "cpu_cores_effective", ">=", 12, 1),
        ("cpu_ram", f"cpu_ram>={cpu}", "cpu_ram", ">=", cpu, mb),
        ("cuda", "cuda_vers>=13.0", "cuda_max_good", ">=", 13.0, 1),
        ("disk", f"disk_space>={disk}", "disk_space", ">=", disk, 1),
        ("inet_down", "inet_down>=500", "inet_down", ">=", 500, 1),
        ("ports", "direct_port_count>=2", "direct_port_count", ">=", 2, 1),
        ("duration", "duration>=1", None, None, None, None),
        ("price", f"dph<={c['ceiling']}", "dph_total", "<=", c["ceiling"], 1),
        ("bandwidth", f"inet_up_cost<={INET_CEIL}", "inet_up_cost", "<=", INET_CEIL, 1),
        ("bandwidth", f"inet_down_cost<={INET_CEIL}", "inet_down_cost", "<=", INET_CEIL, 1),
        ("datacenter", "datacenter=True", "hosting_type", "=", 1, 1) if dataset == "riga" else None,
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
    return [f"{r} fails the §4.2 filter" for r in dict.fromkeys(bad)]


def governance(dataset: str, o: dict, allow: dict) -> list[str]:
    """H-GOV (D16): LIDC anywhere verified; RIGA only with recorded written permission AND a datacenter offer."""
    if dataset == "lidc":
        return []
    if dataset != "riga":
        return [f"unknown bundle {dataset}"]
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


def active(update: dict | None = None) -> dict:
    f = home() / "active_runs.json"
    if update is not None:
        f.write_text(json.dumps(update, indent=1))
    return json.loads(f.read_text()) if f.exists() else {}


def accrued(r: dict, t: float) -> float:
    """A run's cost from its start to t: dph x hours + the bandwidth estimate (BK4-18)."""
    h = max(0.0, t - r["started"]) / 3600
    return r["dph"] * h + r.get("inet_down_cost", 0) * DOWN_GB + r.get("inet_up_cost", 0) * UP_GB_PER_H * h


def spend(committed: bool = False) -> float:
    """Closed sessions from the ledger plus what live runs have accrued so far (committed: to their deadlines)."""
    f, now = home() / "ledger.csv", time.time()
    closed = sum(float(r["cost_usd"] or 0) for r in csv.DictReader(f.open())) if f.exists() else 0.0
    return closed + sum(accrued(r, max(now, r["deadline"]) if committed else now) for r in active().values())


def allowlist() -> dict:
    f = home() / ".vast-allowlist"
    return json.loads(f.read_text()) if f.exists() else {}


@functools.cache
def secret(service: str) -> str:
    cmd = ["security", "find-generic-password", "-s", service, "-w"]
    p = subprocess.run(cmd, capture_output=True, text=True, stdin=subprocess.DEVNULL)
    if p.returncode:
        die(f"keychain item {service} missing (security find-generic-password)")
    return p.stdout.strip()


# --- commands ------------------------------------------------------------------------------------------------------


def vast(*args: str, tries: int = 1):
    env = dict(os.environ)
    if "VAST_API_KEY" not in env:
        env["VAST_API_KEY"] = secret("objdel-vast-scoped")
    for i in range(tries):
        p = subprocess.run(["vastai", *args], capture_output=True, text=True, env=env, stdin=subprocess.DEVNULL)
        if not p.returncode:
            return json.loads(p.stdout) if "--raw" in args else p.stdout
        if i < tries - 1:
            time.sleep(float(os.environ.get("VAST_POLL_S", 10)))
    die(f"vastai {args[:2]} failed: {p.stderr.strip()[-300:]}")


def sh(cmd: list[str], stdin: str | None = None, what: str = "") -> int:
    """ssh / rsync: printed under DRY_RUN=1 (stdin shown only as a label), otherwise run."""
    if dry():
        print("DRY_RUN:", shlex.join(cmd), f"<stdin: {what}>" if stdin is not None else "")
        return 0
    return subprocess.run(cmd, input=stdin, text=True, stdin=None if stdin else subprocess.DEVNULL).returncode


def ssh_cmd(run: dict) -> list[str]:
    opts = f"-o StrictHostKeyChecking=accept-new -o UserKnownHostsFile={home() / 'known_hosts'}"
    return ["ssh", "-p", str(run["port"]), *opts.split(), f"root@{run['ip']}"]


def remote(run: dict, command: str, stdin: str | None = None, what: str = "") -> None:
    if sh([*ssh_cmd(run), command], stdin, what):
        die(f"remote command failed: {command[:120]}")


def rsync(run: dict, src: str, dst: str, *extra: str) -> None:
    e = shlex.join(ssh_cmd(run)[:-1])
    if sh(["rsync", "-az", "--exclude=data/raw/", *extra, "-e", e, src, dst]):
        die(f"rsync {src} -> {dst} failed")


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


def onstart_cmd(run_id: str, cap_s: int, created: float | None = None) -> str:
    """--onstart-cmd: unpack the watchdog and start it; carries no secret (CONTAINER_API_KEY is injected by Vast)."""
    b64 = base64.b64encode((HERE / "watchdog.sh").read_bytes()).decode()
    env = (f"RUN_ID={run_id} WATCHDOG_CAP_S={cap_s} WATCHDOG_CREATED={int(created or time.time())} "
           f"WATCHDOG_PULLED_OK=/workspace/runs/{run_id}/PULLED_OK")  # fmt: skip
    return f"printf %s {b64} | base64 -d > /root/watchdog.sh && {env} nohup bash /root/watchdog.sh >> /root/watchdog.log 2>&1 &"


def offer_cols(o: dict) -> dict:
    return {k: o.get(k) for k in ("dph_total", "inet_up_cost", "inet_down_cost", "hosting_type", "host_id", "machine_id")}


def cmd_offers(a) -> None:
    ds = "riga" if a.riga else "lidc"
    got = vast("search", "offers", query(a.gpu, ds, a.disk), "-o", "dph_total", "--storage", str(a.disk), "--raw")
    ok = [o for o in got if not reasons(o, a.gpu, ds, a.disk)][:5]
    for o in ok:
        print(
            f"{o['id']:>10}  ${o['dph_total']:.3f}/h  {o['gpu_ram'] / 1024:.0f} GB  {o['cpu_cores_effective']} cores"
            f"  {o['cpu_ram'] / 1024:.0f} GB RAM  cuda {o['cuda_max_good']}  hosting_type {o.get('hosting_type')}"
            f"  ${o['inet_up_cost']}/GB up"
        )
    if a.pick is None:
        return
    o = next((o for o in ok if o["id"] == a.pick), None)
    if o is None:
        die(f"offer {a.pick} is not among the admissible offers above")
    (home() / "offers").mkdir(exist_ok=True)
    (home() / "offers" / f"{o['id']}.json").write_text(json.dumps(dict(o, _gpu=a.gpu, _dataset=ds, _disk=a.disk)))
    ledger(event="offer", offer_id=o["id"], gpu=a.gpu, **offer_cols(o), detail=query(a.gpu, ds, a.disk))


def cmd_launch(a) -> None:
    f = home() / "offers" / f"{a.offer}.json"
    if not RUN_ID_RE.fullmatch(a.run_id):
        die(f"run_id must match {RUN_ID_RE.pattern}")
    if a.run_id in active():
        die(f"run {a.run_id} is still active (teardown.sh first)")
    if not f.exists():
        die(f"offer {a.offer} was not picked with offers.sh --pick")
    o = json.loads(f.read_text())
    if why := reasons(o, o["_gpu"], o["_dataset"], o["_disk"]):
        die("; ".join(why))
    cap_s, now = int(a.cap_hours * 3600), time.time()
    bw = {k: o[k] for k in ("inet_up_cost", "inet_down_cost")}
    run = {"instance": None, "offer": a.offer, "dph": o["dph_total"], **bw, "started": now, "deadline": now + cap_s}
    committed, worst = spend(committed=True), accrued(run, run["deadline"])
    if committed + worst > HARD_STOP:
        die(f"hard stop: ${committed:.2f} committed + ${worst:.2f} worst case > ${HARD_STOP} (G10)")
    if committed + worst >= ALARM and not a.owner_approved:
        die(f"alarm: ${committed:.2f} committed + ${worst:.2f} worst case >= ${ALARM}; pass --owner-approved only "
            "with the owner's written approval (G10, §4.6)")  # fmt: skip
    if (live := vast("show", "instances", "--raw")) and not a.allow_parallel:
        die(f"{len(live)} instance(s) already exist; one instance at a time unless --allow-parallel")
    allow = allowlist()
    image = allow.get("image") or die("no image in .vast-allowlist")
    login = f"-u {allow.get('ghcr_user', '')} -p {secret('objdel-ghcr-read')} ghcr.io"
    flags = ["--ssh", "--direct", "--cancel-unavail", "--label", a.run_id, "--raw"]
    sigs = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
    for s in sigs:  # an interrupt after `create` must reach the destroy below (BK4-3)
        signal.signal(s, signal.default_int_handler)
    out = vast("create", "instance", str(a.offer), "--image", image, "--login", login, "--disk", str(o["_disk"]),
               "--onstart-cmd", onstart_cmd(a.run_id, cap_s, now), *flags)  # fmt: skip
    run["instance"] = iid = int(out["new_contract"])
    try:  # G8: any failure or interrupt from here on destroys the instance
        active({**active(), a.run_id: run})
        ledger(event="launch", run_id=a.run_id, instance_id=iid, offer_id=a.offer, gpu=o["_gpu"], **offer_cols(o),
               detail=f"cap_s={cap_s} image={image}" + " owner_approved=yes" * a.owner_approved)  # fmt: skip
        t0, timeout = time.time(), float(os.environ.get("VAST_LAUNCH_TIMEOUT_S", 900))
        while True:
            inst = vast("show", "instance", str(iid), "--raw")
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
        destroy(iid, a.run_id, "destroy")
        raise


def get_run(run_id: str) -> dict:
    return active().get(run_id) or die(f"no active run {run_id}")


def code_files() -> list[str]:
    """The tracked code snapshot, without data, docs, lab or results dirs, or RIGA-derived tables (D16, BK4-13).
    tests/ ships: k3_harness -> phase1.reward_gates -> tests.oracles."""
    git = ["git", "ls-files", "-co", "--exclude-standard"]
    files = subprocess.run(git, cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
    drop = re.compile(r"^(data|docs)/|^experiments/[^/]+/(lab|results)/|riga[^/]*\.(csv|jsonl?)$", re.I)
    files = [f for f in files if not drop.search(f)]
    if not any(f.startswith("tests/oracles/") for f in files):
        die("tests/oracles missing from the code snapshot")
    return files


def cmd_sync_up(a) -> None:
    run = get_run(a.run_id)
    o = json.loads((home() / "offers" / f"{run['offer']}.json").read_text())
    allow = allowlist()
    if why := governance(a.bundle, o, allow):
        die("H-GOV refuses: " + "; ".join(why))
    name = "bundle-riga.tar.age" if a.bundle == "riga" else f"bundle-{a.bundle}.tar"
    files = [a.bundle_dir / name, a.bundle_dir / f"{name}.sha256"]
    if is_raw(a.bundle_dir) or any(is_raw(f) for f in files):
        die("refusing to send anything under data/raw")
    lst = home() / "code_files.txt"
    lst.write_text("\n".join(code_files()) + "\n")
    host = f"root@{run['ip']}:"
    rsync(run, f"{ROOT}/", f"{host}/workspace/src/", f"--files-from={lst}")
    remote(run, "cd /workspace/src && uv pip install --python /opt/venv/bin/python --no-deps --no-build-isolation -e ."
                f" && cd / && /opt/venv/bin/python -I -c {shlex.quote(IMPORT_CHECK)}")  # fmt: skip
    rsync(run, *map(str, files), f"{host}/workspace/data/")
    verify, detail = f"cd /workspace/data && sha256sum -c {name}.sha256", ""
    if a.bundle == "riga":  # D16: identity streamed into age over this session's stdin; plaintext only in /dev/shm
        # BK4-19: Docker's /dev/shm defaults to 64 MB and V0 has not shown Vast's; refuse before decrypting
        need = -(-(a.bundle_dir / name).stat().st_size * 11 // 10240)  # KB, ciphertext + 10 %
        shm = f"avail=$(df -Pk /dev/shm | awk 'NR==2 {{print $4}}') && [ \"$avail\" -ge {need} ]"
        cmd = (f"{verify} && {shm} && mkdir -p /dev/shm/riga && age -d -i /dev/stdin {name} | tar -x -C /dev/shm/riga"
               " && cd /dev/shm/riga && sha256sum -c --quiet SHA256SUMS")  # fmt: skip
        remote(run, cmd, secret("objdel-riga-age") + "\n", "age identity from the keychain")
        detail = f"recipient={allow.get('riga_recipient', '')} shm=UNVERIFIED until V0 (BK4-19)"
    else:
        remote(run, f"{verify} && mkdir -p {a.bundle} && tar -xf {name} -C {a.bundle}"
                    f" && cd {a.bundle} && sha256sum -c --quiet SHA256SUMS")  # fmt: skip
    ledger(event="sync", run_id=a.run_id, instance_id=run["instance"], bundle=name, detail=detail)


def chain(cfg: Path, run_id: str) -> list[list[str]]:
    """A config's stages on the box: train_sft, then (branches saved) select_ckpt with the config's G-vision block
    (`--no-vision` only for `vision: off` timing arms) and g_probe on <run>/select (selected branch, then the D12(c)
    repair). Paths are the box's: data under /workspace/data/lidc, runs under /workspace/runs/<id>."""
    import yaml

    c, data, rd = yaml.safe_load(Path(cfg).read_text()), "/workspace/data/lidc", f"/workspace/runs/{run_id}"
    run, ev = f"{rd}/{c['output_dir']}", f"{data}/eval_views"
    stages = [["-m", "experiments.phase2.train_sft", "--config", Path(cfg).resolve().relative_to(ROOT).as_posix()]]
    if not c.get("branches"):
        return stages
    v = c.get("vision")
    vis = ["--vision", json.dumps({"prepared": f"{data}/{v['prepared']}", "reference": f"{rd}/{v['reference']}"})
           ] if v else ["--no-vision"]  # fmt: skip
    return [*stages,
            ["-m", "experiments.phase2.select_ckpt", "--run", run, "--shifted", f"{ev}/lidc-val-shifted", "--rescaled",
             f"{ev}/lidc-val-rescaled", "--prepared", f"{data}/prepared", "--probes", ev, *vis],
            ["-m", "experiments.phase2.g_probe", "--model", c["model"], "--select", f"{run}/select", "--train",
             f"{data}/prepared", "--probe", f"{ev}/lidc-val-shifted", "--out", f"{run}/g_probe"]]  # fmt: skip


def vision_refs(cfg: Path, run_id: str) -> list[str]:
    """The G-vision reference JSON a config reads on the box ([] for a number or `vision: off`)."""
    import yaml

    ref = (yaml.safe_load(Path(cfg).read_text()).get("vision") or {}).get("reference")
    return [f"/workspace/runs/{run_id}/{ref}"] if isinstance(ref, str) else []


def cmd_run(a) -> None:
    """MODULE [ARGS] runs one module; a configs/*.yaml path runs its whole chain (`chain`) under one cap, refused on
    the box while a G-vision reference it reads is missing; `vision-ref CONFIG...` writes those references (R0, same
    run dir: train_sft --vision_reference_only true)."""
    run, guard = get_run(a.run_id), []
    if a.module == "vision-ref":
        if not a.args or not all(vision_refs(ROOT / c, a.run_id) for c in a.args):
            die("vision-ref needs configs whose vision.reference is a .json path")
        stages = [["-m", "experiments.phase2.train_sft", "--config", c, "--vision_reference_only", "true"] for c in a.args]
    elif a.module.endswith(".yaml"):
        stages = chain(ROOT / a.module, a.run_id)
        guard = [f"{{ test -f {r} || {{ echo 'missing G-vision reference {r}: run.sh {a.run_id} vision-ref "
                 f"{a.module} first' >&2; false; }}; }}" for r in vision_refs(ROOT / a.module, a.run_id)]  # fmt: skip
    elif re.fullmatch(r"experiments\.phase2\.\w+", a.module):
        stages = [["-m", a.module, *a.args]]
    else:
        die("module must be experiments.phase2.<script> or experiments/phase2/configs/<name>.yaml")
    cap = int(min(a.cap_hours * 3600, run["deadline"] - time.time() - MARGIN_S))
    if cap <= 0:
        die("no time left before the instance cap")
    rd = f"/workspace/runs/{a.run_id}"
    jobs = " && ".join([*guard, *(f"/opt/venv/bin/python -I {shlex.join(s)}" for s in stages)])
    inner = (  # EXIT is removed first: the watchdog's idle destroy needs it (BK4-11)
        f"cd /workspace/src && export OBJDEL_RUN_DIR={rd} && mkdir -p {rd} && rm -f {rd}/EXIT && "
        f"(while :; do date -u +%s > {rd}/HEARTBEAT; sleep 60; done) & "
        f"timeout {cap}s bash -c {shlex.quote(jobs)} > {rd}/run.log 2>&1; echo $? > {rd}/EXIT; kill %1"
    )
    remote(run, f"tmux new-session -d -s {a.run_id} {shlex.quote('bash -c ' + shlex.quote(inner))}")
    ledger(event="run", run_id=a.run_id, instance_id=run["instance"], detail=f"{a.module} cap_s={cap}")


def cmd_pull(a) -> None:
    """rsync verifies every transferred file with a whole-file checksum; only then is PULLED_OK touched."""
    for run_id in [a.run_id] if a.run_id else list(active()):
        run, dst = get_run(run_id), Path(os.environ.get("VAST_RUNS", ROOT / "runs")) / run_id
        dst.mkdir(parents=True, exist_ok=True)
        rsync(run, f"root@{run['ip']}:/workspace/runs/{run_id}/", f"{dst}/", "--exclude=PULLED_OK")
        remote(run, f"touch /workspace/runs/{run_id}/PULLED_OK")


def cmd_teardown(a) -> None:
    run = get_run(a.run_id)
    cmd_pull(a)
    destroy(run["instance"], a.run_id, "teardown")
    inv = vast("show", "invoices-v1", "--raw", tries=3)
    (home() / "invoices").mkdir(exist_ok=True)
    (home() / "invoices" / f"{a.run_id}.json").write_text(json.dumps(inv, indent=1))
    ledger(event="invoice", run_id=a.run_id, instance_id=run["instance"],
           detail=f"UNPARSED: raw invoices/{a.run_id}.json (format UNVERIFIED until V0); reconcile by hand")  # fmt: skip


def alert(msg: str) -> None:
    print(f"vast ALERT: {msg}", file=sys.stderr)
    if sys.platform == "darwin" and not dry():
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


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="vast", description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("offers")
    p.add_argument("gpu"), p.add_argument("--riga", action="store_true"), p.add_argument("--pick", type=int)
    p.add_argument("--disk", type=int, default=80)
    p = sub.add_parser("launch")
    p.add_argument("offer", type=int), p.add_argument("run_id"), p.add_argument("--cap-hours", type=float, default=4)
    p.add_argument("--allow-parallel", action="store_true"), p.add_argument("--owner-approved", action="store_true")
    p = sub.add_parser("sync_up")
    p.add_argument("run_id"), p.add_argument("--bundle", default="lidc", choices=["lidc", "riga"])
    p.add_argument("--bundle-dir", type=Path, default=ROOT / "build" / "bundles")
    p = sub.add_parser("run")
    p.add_argument("run_id"), p.add_argument("module"), p.add_argument("--cap-hours", type=float, default=24)
    p.add_argument("args", nargs=argparse.REMAINDER)
    sub.add_parser("pull").add_argument("run_id", nargs="?")
    sub.add_parser("teardown").add_argument("run_id")
    sub.add_parser("sweep_orphans")
    a = ap.parse_args(argv)
    globals()[f"cmd_{a.cmd}"](a)


if __name__ == "__main__":
    main()
