"""Vast AI operations: the CLI, `python -I -m objdel.infra.vast <cmd>`. Every
decision and every command is built here or in objdel.infra.vast (offers, ledger, remote, chain); the `*.sh` files in
ops/vast/ are symlinks to the one shim `vast.sh` (the link's name is the command), and `watchdog.sh`
there is the only logic that runs on the box.

    offers.sh GPU [--disk GB] [--pick OFFER_ID|top]            search (US/CA hosts), top 5; --pick saves it
    launch.sh OFFER_ID RUN_ID --config CONFIG.yaml [--cap-hours H] [--idle-minutes M] [--allow-parallel]
              [--owner-approved]
    sync_up.sh RUN_ID [--bundle lidc] [--bundle-dir DIR] [--code-only]   code (git ls-files) + bundle, SHA-256
    run.sh RUN_ID MODULE|CONFIG.yaml [ARGS]                    timeout + tmux + heartbeat; a config runs its chain
    pull.sh [RUN_ID]                                           rsync runs/<id>/ to the Mac (not PULL_EXCLUDE), then
                                                               PULLED_OK = the pull's start time (no id: all)
    push_hf.sh RUN_ID CONFIG.yaml                              the run's selected branch to the private HF repo
    fetch_hf.sh RUN_ID REVISION                                a revision of that repo onto the box (runs/RUN_ID/ckpt)
                                                               (hf_push; after pull.sh, before teardown.sh)
    teardown.sh RUN_ID                                         final pull, destroy, verify absent, invoices, ledger
    sweep_orphans.sh                                           destroy unknown / overdue instances, close vanished runs;
                                                               alarm and hard stop: see ALARM and HARD_STOP below

Money: a launch is refused when the money committed (closed sessions + every active run to its deadline) plus its own
worst case passes HARD_STOP, and needs --owner-approved (recorded in the ledger) when that projection reaches ALARM
(both in reference/spend_limits.json: the owner edits the file, not the code).
Costs are dph_total x time plus a bandwidth estimate (DOWN_GB once, UP_GB_PER_H per hour, at the offer's $/GB);
teardown saves the raw charges and records what Vast billed beside that estimate (`billed`).

Environment: VAST_HOME (ledger.csv, active_runs.json, offers/, .vast-allowlist; default:
<repo>/ops/vast), VAST_RUNS (pull target, default <repo>/runs), DRY_RUN=1 (ssh, rsync and `vastai
create` / `destroy` are printed, never run; `vastai` calls that only read still run, so tests put a fake one on
PATH), VAST_POLL_S, VAST_LAUNCH_TIMEOUT_S. Secrets live in the macOS keychain and are read with `security`:
objdel-vast-scoped (scoped API key, passed to vastai as VAST_API_KEY), objdel-ghcr-read (read-only GHCR PAT, Vast
--login), objdel-hf-write (HF token: checked by objdel.infra.hf.token_problems, then only ever the stdin of the ssh
upload command).
"""

import argparse
import json
import math
import os
import re
import shlex
import signal
import sys
import time
from pathlib import Path

from objdel.core import datasets
from objdel.core.io import code_sha, git_sha, read_json
from objdel.core.provenance import repo_root
from objdel.core.spec import CLI as CLI_PACKAGE
from objdel.core.spec import load_spec
from objdel.infra import hf
from objdel.infra.errors import InfraError
from objdel.infra.vast.chain import BOX_DATA, chain, push_target, vision_refs
from objdel.infra.vast.ledger import (
    ALARM,
    DOWN_GB,
    HARD_STOP,
    UP_GB_PER_H,
    accrued,
    active,
    billed,
    home,
    ledger,
    load,
    locked,
    pull_lock,
    spend,
    update_active,
)
from objdel.infra.vast.offers import is_raw, offer_cols, query, reasons
from objdel.infra.vast.remote import (
    BAD_STATUS,
    GRACE_S,
    IMPORT_CHECK,
    MARGIN_S,
    PULL_EXCLUDE,
    PULL_S,
    RUN_ID_RE,
    alert,
    code_files,
    destroy,
    dry_run,
    get_run,
    onstart_cmd,
    remote,
    rsync,
    secret,
    vast,
)


def cmd_offers(a) -> None:
    """offers.sh: list the five cheapest admissible offers; --pick ID saves one for launch."""
    got = vast("search", "offers", query(a.gpu, a.disk), "-o", "dph_total", "--storage", str(a.disk), "--raw")
    ok = [o for o in got if not reasons(o, a.gpu, a.disk)]  # --pick ID: any of them, not only the 5 shown
    for o in ok[:5]:
        print(f"{o['id']:>10}  ${o['dph_total']:.3f}/h  {o['gpu_ram'] / 1024:.0f} GB  {o['cpu_cores_effective']} cores"
              f"  {o['cpu_ram'] / 1000:.0f} GB RAM  cuda {o['cuda_max_good']}  hosting_type {o.get('hosting_type')}"
              f"  ${o['inet_up_cost']}/GB up")  # fmt: skip
    if a.pick is None:
        return
    o = ok[0] if a.pick == "top" and ok else next((o for o in ok if str(o["id"]) == a.pick), None)  # top: no race
    if o is None:
        raise InfraError(f"offer {a.pick} is not among the admissible offers")
    (home() / "offers").mkdir(exist_ok=True)
    (home() / "offers" / f"{o['id']}.json").write_text(json.dumps(dict(o, _gpu=a.gpu, _disk=a.disk)))
    ledger(event="offer", offer_id=o["id"], gpu=a.gpu, **offer_cols(o), detail=query(a.gpu, a.disk))
    print(f"picked {o['id']}")


def launch_refusal(offer: dict, cap_hours: float) -> str:
    """Why this offer may not launch (empty when it may): no ledger.csv, or committed spend + the worst case to the
    deadline > HARD_STOP."""
    money = [offer["dph_total"], offer["inet_up_cost"], offer["inet_down_cost"]]
    if not all(math.isfinite(x) and x >= 0 for x in money) or not (math.isfinite(cap_hours) and cap_hours > 0):
        return f"refused: price {money} and --cap-hours {cap_hours} must be finite, non-negative and the cap positive"
    if not (home() / "ledger.csv").is_file():  # a missing ledger is unknown spend, not $0
        return f"no ledger.csv under {home()}: the spend guard cannot count what was spent (G10)"
    now = time.time()
    bw = {k: offer[k] for k in ("inet_up_cost", "inet_down_cost")}
    run = {"dph": offer["dph_total"], **bw, "started": now, "deadline": now + int(cap_hours * 3600) + PULL_S}
    committed, worst = spend(committed=True), accrued(run, run["deadline"])
    if not committed + worst <= HARD_STOP:  # not "> ": a NaN price must refuse too
        return f"hard stop: ${committed:.2f} committed + ${worst:.2f} worst case > ${HARD_STOP} (G10)"
    return ""


def cmd_launch(a) -> None:
    """launch.sh: every refusal first (run id, preflight stamp, offer, spec caps, money, one instance), then create
    the instance, wait for ssh and record it; any failure destroys it."""
    if not RUN_ID_RE.fullmatch(a.run_id):
        raise InfraError(f"run_id must match {RUN_ID_RE.pattern}")
    if a.run_id in active():
        raise InfraError(f"run {a.run_id} is still active (teardown.sh first)")
    from objdel.infra import preflight

    if why := preflight.refusal(Path(a.config)):  # the mandatory Mac preflight (make preflight): no pass, no rental
        raise InfraError(f"preflight gate: {why}")
    if not (o := load(f"offers/{a.offer}.json")):
        raise InfraError(f"offer {a.offer} was not picked with offers.sh --pick")
    if why := reasons(o, o["_gpu"], o["_disk"]):
        raise InfraError("; ".join(why))
    compute = load_spec(Path(a.config)).compute  # the spec bounds the flag: one cap source per launch
    if a.cap_hours > compute.cap_hours:
        raise InfraError(f"--cap-hours {a.cap_hours} is above the spec's compute.cap_hours {compute.cap_hours}")
    if o["_gpu"] != compute.gpu:
        raise InfraError(f"offer {a.offer} is a {o['_gpu']}: {o['_gpu']} is not compute.gpu {compute.gpu}")
    if why := launch_refusal(o, a.cap_hours):
        raise InfraError(why)
    cap_s, now = int(a.cap_hours * 3600) + PULL_S, time.time()
    bw = {k: o[k] for k in ("inet_up_cost", "inet_down_cost")}
    run = {"instance": None, "offer": a.offer, "dph": o["dph_total"], **bw, "started": now, "deadline": now + cap_s}
    committed, worst = spend(committed=True), accrued(run, run["deadline"])
    if committed + worst >= ALARM and not a.owner_approved:
        raise InfraError(f"alarm: ${committed:.2f} committed + ${worst:.2f} worst case >= ${ALARM}; pass "
                         "--owner-approved only with the owner's written approval (G10, §4.6)")
    if (live := vast("show", "instances", "--raw")) and not a.allow_parallel:
        raise InfraError(f"{len(live)} instance(s) already exist; one instance at a time unless --allow-parallel")
    allow = load(".vast-allowlist")
    image = allow.get("image")
    if not image:
        raise InfraError("no image in .vast-allowlist")
    login = f"-u {allow.get('ghcr_user', '')} -p {secret('objdel-ghcr-read')} ghcr.io"
    flags = ["--ssh", "--direct", "--cancel-unavail", "--label", a.run_id, "--raw"]
    sigs = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
    onstart = onstart_cmd(a.run_id, cap_s, now, int(a.idle_minutes * 60))
    for s in sigs:  # an interrupt from `create` on must reach the destroy below
        signal.signal(s, signal.default_int_handler)
    iid = None
    try:  # any failure or interrupt from here on destroys the instance
        out = vast("create", "instance", str(a.offer), "--image", image, "--login", login, "--disk", str(o["_disk"]),
                   "--onstart-cmd", onstart, *flags)  # fmt: skip
        if dry_run():  # printed, not rented: nothing to record, poll or destroy
            return
        run["instance"] = iid = int(out["new_contract"])
        update_active(lambda runs: {**runs, a.run_id: run})
        ledger(event="launch", run_id=a.run_id, instance_id=iid, offer_id=a.offer, gpu=o["_gpu"], **offer_cols(o),
               detail=f"cap_s={cap_s} image={image}" + " owner_approved=yes" * a.owner_approved)  # fmt: skip
        t0, timeout = time.time(), float(os.environ.get("VAST_LAUNCH_TIMEOUT_S", 900))
        while True:
            inst = vast("show", "instance", str(iid), "--raw", tries=3)
            status, port = inst.get("actual_status"), ((inst.get("ports") or {}).get("22/tcp") or [{}])[0]
            if status == "running" and inst.get("public_ipaddr") and port.get("HostPort"):  # ports come late
                run.update(ip=inst["public_ipaddr"], port=int(port["HostPort"]))
                update_active(lambda runs: {**runs, a.run_id: run})
                print(f"running: ssh -p {run['port']} root@{run['ip']}")
                return
            if status in BAD_STATUS or time.time() - t0 > timeout:
                why = status if status in BAD_STATUS else "timeout"
                raise InfraError(f"instance {iid} {why} before running; destroyed")
            time.sleep(float(os.environ.get("VAST_POLL_S", 10)))
    except (Exception, KeyboardInterrupt):  # an interrupt is a KeyboardInterrupt (signal handlers above)
        for s in sigs:
            signal.signal(s, signal.SIG_IGN)
        # create failed or its reply was lost: a box may exist all the same, found by its label
        ids = [iid] if iid else [i["id"] for i in vast("show", "instances", "--raw", tries=3)
                                 if i.get("label") == a.run_id]
        for i in ids:
            destroy(i, a.run_id, "destroy")
        raise


def cmd_sync_up(a) -> None:
    """sync_up.sh: send and install the code snapshot, check the GPU, then send and unpack the bundle."""
    run = get_run(a.run_id)
    name = f"bundle-{a.bundle}.tar"
    files = [a.bundle_dir / name, a.bundle_dir / f"{name}.sha256"]
    if any(map(is_raw, [a.bundle_dir, *files])):
        raise InfraError("refusing to send anything under data/raw")
    lst, sha = home() / "code_files.txt", git_sha(repo_root())
    lst.write_text("\n".join(code_files()) + "\n")
    host = f"root@{run['ip']}:"
    rsync(run, f"{repo_root()}/", f"{host}/workspace/src/", opts=(f"--files-from={lst}",))
    remote(run, f"printf '%s\\n' {shlex.quote(sha)} > /workspace/src/CODE_SHA")  # the box has no .git
    remote(run, "cd /workspace/src && uv pip install --python /opt/venv/bin/python --no-deps --no-build-isolation -e ."
                f" && cd / && /opt/venv/bin/python -I -c {shlex.quote(IMPORT_CHECK)}")  # fmt: skip
    if a.code_only:  # a code fix mid-session: the bundle is already on the box
        return ledger(event="sync", run_id=a.run_id, instance_id=run["instance"], detail="code only")
    # GPU check before the bundle (a busy, power-capped card was once found only after a 12 min upload)
    check = f"bash ops/docker/image_check.sh /workspace/runs/{a.run_id}/image_check.json"
    remote(run, f"cd /workspace/src && {check}")
    rsync(run, *map(str, files), f"{host}/workspace/data/")
    remote(run, f"cd /workspace/data && sha256sum -c {name}.sha256 && mkdir -p {a.bundle} && tar -xf {name} -C "
                f"{a.bundle} && cd {a.bundle} && sha256sum -c --quiet SHA256SUMS")  # fmt: skip
    ledger(event="sync", run_id=a.run_id, instance_id=run["instance"], bundle=name)


# The modules run.sh runs by name: a command module (objdel.cli.<group>.<command>) or a tool script (tools.*)
RUNNABLE = re.compile(rf"{re.escape(CLI_PACKAGE)}\.\w+\.\w+|tools(\.\w+)+")


def cmd_run(a) -> None:
    """MODULE [ARGS] runs one module (RUNNABLE); a run spec specs/<name>.yaml runs its `chain` under one cap,
    refused on the box while a G-vision reference it reads is missing.
    Each job has its own tmux session and log (<run_id>-<name>, <name>.log; a module's name is its last two parts,
    e.g. train_rl), so a second job can run beside the first."""
    run, guard = get_run(a.run_id), []
    if a.module.endswith(".yaml") and a.args:  # a trailing --cap-hours would be silently ignored
        raise InfraError(f"flags go before the config (run.sh RUN_ID [--cap-hours H] CONFIG); got {' '.join(a.args)}")
    if a.module.endswith(".yaml"):
        stages, name = chain(repo_root() / a.module, a.run_id), Path(a.module).stem
        guard = [f"/opt/venv/bin/python -I -m objdel.infra.vast run_spec {shlex.quote(a.module)} "
                 f"/workspace/runs/{a.run_id}/run_spec.json"]  # the spec, code and versions, before any stage
        guard += [f"{{ test -f {r} || {{ echo 'missing G-vision reference {r}: run train_sft "
                 f"--vision-reference-only first' >&2; false; }}; }}"
                 for r in vision_refs(repo_root() / a.module, a.run_id)]  # fmt: skip
    elif RUNNABLE.fullmatch(a.module):
        stages, name = [["-m", a.module, *a.args]], "_".join(a.module.split(".")[-2:])
    else:
        raise InfraError(f"module must be {CLI_PACKAGE}.<group>.<command>, a tool script or a run spec "
                         "specs/<name>.yaml")
    cap = int(min(a.cap_hours * 3600, run["deadline"] - time.time() - MARGIN_S))
    if cap <= 0:
        raise InfraError("no time left before the instance cap")

    def job(s: list[str]) -> str:  # "~": a failure is logged and the chain goes on
        soft = s[0] == "~"
        cmd = shlex.join(s[soft:]) if s[soft] == "bash" else f"/opt/venv/bin/python -I {shlex.join(s[soft:])}"
        return f"{{ {cmd} || echo {shlex.quote('step failed, continuing: ' + cmd)}; }}" if soft else cmd

    rd = f"/workspace/runs/{a.run_id}"
    jobs = " && ".join([*guard, *map(job, stages)])
    # the session clock of train_rl --budget-usd (callbacks.BudgetStop): accrued(run, t) = fixed + rate x hours
    rate, fixed = run["dph"] + run.get("inet_up_cost", 0) * UP_GB_PER_H, run.get("inet_down_cost", 0) * DOWN_GB
    clock = f"OBJDEL_SESSION_START={run['started']!r} OBJDEL_RATE_USD_PER_H={rate!r} OBJDEL_FIXED_USD={fixed!r}"
    # EXIT is removed first: the watchdog's idle destroy needs it; `;` keeps cd/export in this shell
    inner = (f"cd /workspace/src && export OBJDEL_RUN_DIR={rd} {clock} && mkdir -p {rd} && rm -f {rd}/EXIT "
             "|| exit 1; "
             f"(while :; do date -u +%s > {rd}/HEARTBEAT; sleep 60; done) & "
             f"timeout -k 120 {cap}s bash -c {shlex.quote(jobs)} > {rd}/{name}.log 2>&1; "
             f"echo $? > {rd}/EXIT; kill %1")  # fmt: skip
    remote(run, f"tmux new-session -d -s {a.run_id}-{name} {shlex.quote('bash -c ' + shlex.quote(inner))}")
    ledger(event="run", run_id=a.run_id, instance_id=run["instance"], detail=f"{a.module} cap_s={cap}")


def cmd_run_spec(a) -> None:
    """On the box, the first job of a spec's chain: OUT = the spec, the code sha (objdel.core.io.code_sha:
    CODE_SHA, written by sync_up, else git; else "unknown"), the trl, torch and transformers versions and the reward
    fingerprint, as strict JSON. A different spec already in OUT stops the chain: a second run with the same id never
    overwrites it silently."""
    from importlib.metadata import version

    from objdel.verify.reward.spec import fingerprint

    rec = {"spec": load_spec(repo_root() / a.spec).model_dump(mode="json"),
           "git_sha": code_sha(repo_root(), required=False) or "unknown",
           **{p: version(p) for p in ("trl", "torch", "transformers")}, "reward_fingerprint": fingerprint()}
    out = Path(a.out)
    if out.is_file() and read_json(out)["spec"] != rec["spec"]:
        raise InfraError(f"{out} holds a different spec: use a new run id")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=1, allow_nan=False))


def cmd_pull(a) -> None:
    """rsync verifies every transferred file with a whole-file checksum; only then does PULLED_OK appear, carrying
    the pull's start time (PULLED_OK.next, touched before rsync), so a file written during the pull is newer."""
    for run_id in [a.run_id] if a.run_id else list(active()):
        with pull_lock(run_id):  # a teardown's pull or a second pull would race on PULLED_OK.next
            run, dst = get_run(run_id), Path(os.environ.get("VAST_RUNS", repo_root() / "runs")) / run_id
            dst.mkdir(parents=True, exist_ok=True)
            rd = f"/workspace/runs/{run_id}"
            remote(run, f"touch {rd}/PULLED_OK.next")
            rsync(run, f"root@{run['ip']}:{rd}/", f"{dst}/", opts=tuple(f"--exclude={x}" for x in PULL_EXCLUDE))
            remote(run, f"mv {rd}/PULLED_OK.next {rd}/PULLED_OK")


def hf_token() -> str:
    """The HF write token, refused unless it is fine-grained with write on the one repo only (it goes to a rented
    box over ssh stdin)."""
    token = secret("objdel-hf-write")
    if bad := hf.token_problems(hf.whoami(token), hf.REPO):
        raise InfraError(f"secret objdel-hf-write (.secrets or keychain) is too broad for a rented box: "
                         f"{'; '.join(bad)}. Recreate it at huggingface.co/settings/tokens as fine-grained with only "
                         f"{hf.REPO}: write (no other boxes ticked)")
    return token


def cmd_push_hf(a) -> None:
    """hf_push on the box for the run spec CONFIG's run (push_target: its last rl or sft stage); the token is checked
    here first."""
    run = get_run(a.run_id)
    path, out, rl = push_target(load_spec(repo_root() / a.config), a.run_id)
    kind = "--rl" if rl else "--run"
    token = hf_token()
    cmd = (f"cd /workspace/src && /opt/venv/bin/python -I -m {CLI_PACKAGE}.hub.push {kind} {path} "
           f"--repo {hf.REPO} --revision {out}")  # fmt: skip
    remote(run, cmd, token + "\n", "HF token from objdel.infra.secrets")


ROWS_CHECK = ("import hashlib, json, pathlib, sys; d, want = pathlib.Path(sys.argv[1]), sys.argv[2]; "
              "m = d / 'manifest.json'; ok = m.is_file() and hashlib.sha256(m.read_bytes()).hexdigest() == want "
              "and sum(p.is_file() for p in d.rglob('*')) - 1 == json.loads(m.read_text())['schedule']['files_count']; "
              "sys.exit(0 if ok else f'rows root {d}: no manifest.json, not SHA-256 {want}, or not files_count files')")


UNTAR = ("import os, sys, tarfile; t, d = sys.argv[1:]; f = tarfile.open(t); f.extractall(d, filter='data'); "
         "f.close(); os.remove(t)")  # hf.ROWS_TAR into the root, then gone (ROWS_CHECK counts the root's files)


def pinning(revision: str) -> list:
    """The rl stages of specs/ that pin REVISION as their rows revision."""
    return [s for p in sorted((repo_root() / "specs").glob("*.yaml")) for s in load_spec(p).stages
            if s.kind == "rl" and s.rows_revision == revision]


def rows_pin(revision: str) -> str:
    """The rows_manifest_sha256 that the rl stages of specs/ pin beside REVISION."""
    pins = {s.rows_manifest_sha256 for s in pinning(revision)}
    if len(pins) != 1:
        raise InfraError(f"--rows {revision}: the plans pin {len(pins)} manifest SHA-256 values for it "
                         "(need exactly one)")
    return pins.pop()


def cmd_fetch_hf(a) -> None:
    """hf_pull on the box: REVISION of the private HF repo to /workspace/runs/RUN_ID/NAME (--as; default ckpt, RL's
    start checkpoint); or --rows REVISION ROOT: the repo folder rows/ROOT/ to <box data>/rows/ROOT (the dataset of
    the stages that pin REVISION, chain.BOX_DATA), then
    the same ssh call unpacks hf.ROWS_TAR there (UNTAR) and checks the root (manifest.json with the SHA-256 the plans
    pin for REVISION, files_count files)."""
    if (a.revision is None) == (a.rows is None):
        raise InfraError("fetch_hf: give REVISION or --rows REVISION ROOT")
    rev, name = a.rows or (a.revision, a.name)
    if "/" in name or ".." in name or name in ("", "."):
        raise InfraError(f"fetch_hf: folder name {name!r} may not hold / or ..")
    pull = (f"cd /workspace/src && /opt/venv/bin/python -I -m {CLI_PACKAGE}.hub.pull --repo {hf.REPO} "
            f"--revision {shlex.quote(rev)}")
    if a.rows:
        pin = rows_pin(rev)  # refuses a revision no plan pins, or pins twice
        if len(found := {s.dataset for s in pinning(rev)}) != 1:
            raise InfraError(f"--rows {rev}: the plans that pin it name {len(found)} datasets (need exactly one)")
        box = BOX_DATA.format(dataset=found.pop())
        d = f"{box}/rows/{name}"
        cmd = (f"{pull} --out {box} --allow {shlex.quote(f'rows/{name}/**')} && /opt/venv/bin/python "
               f"-I -c {shlex.quote(UNTAR)} {shlex.quote(f'{d}/{hf.ROWS_TAR}')} {shlex.quote(d)} && "
               f"/opt/venv/bin/python -I -c {shlex.quote(ROWS_CHECK)} {shlex.quote(d)} {pin}")
    else:
        cmd = f"{pull} --out /workspace/runs/{a.run_id}/{name}"
    run, token = get_run(a.run_id), hf_token()
    remote(run, cmd, token + "\n", "HF token from objdel.infra.secrets")


def cmd_teardown(a) -> None:
    """teardown.sh: pull the run (three tries), then destroy the instance and record the session's cost."""
    run, lost = get_run(a.run_id), None
    try:
        for _ in range(3 if "port" in run else 0):  # no port: the box never got ssh, so nothing is on it to pull
            try:
                cmd_pull(a)
                lost = None
                break
            except InfraError as e:  # ssh refused, so the pull failed; the box must still go
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
        raise InfraError(f"final pull failed ({lost}); instance {run['instance']} destroyed")


def cmd_sweep_orphans(a) -> None:
    """sweep_orphans.sh: close runs whose instance is gone, destroy orphan or overdue instances, and destroy every
    instance at the hard stop."""
    listed, now = vast("show", "instances", "--raw", tries=3), time.time()
    for k, r in active().items():  # gone without our destroy (watchdog, host loss): close at its last state
        if r["instance"] not in {i["id"] for i in listed} and now - r["started"] > GRACE_S:
            alert(f"run {k}: instance {r['instance']} no longer listed; closing it")
            with locked():
                update_active(lambda runs, k=k: {j: x for j, x in runs.items() if j != k})
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
    "offers": {"gpu": {}, "--pick": {}, "--disk": {"type": int, "default": 80}},
    "launch": {"offer": {"type": int}, "run_id": {}, "--config": {"required": True, "help": "the run's config.yaml"},
               "--cap-hours": {"type": float, "default": 4},
               "--allow-parallel": FLAG, "--owner-approved": FLAG,
               "--idle-minutes": {"type": float, "default": 30, "help": "watchdog idle destroy after a pulled EXIT"}},
    "sync_up": {"run_id": {}, "--bundle": {"default": "lidc", "choices": sorted(datasets.DATASETS)},
                "--bundle-dir": {"type": Path, "default": repo_root() / "build" / "bundles"}, "--code-only": FLAG},
    "run": {"run_id": {}, "module": {}, "--cap-hours": {"type": float, "default": 24},
            "args": {"nargs": argparse.REMAINDER}},
    "run_spec": {"spec": {}, "out": {}},
    "pull": {"run_id": {"nargs": "?"}},
    "push_hf": {"run_id": {}, "config": {}},
    "fetch_hf": {"run_id": {}, "revision": {"nargs": "?"}, "--rows": {"nargs": 2, "metavar": ("REVISION", "ROOT")},
                 "--as": {"dest": "name", "default": "ckpt", "help": "folder under the run dir (default ckpt)"}},
    "teardown": {"run_id": {}},
    "sweep_orphans": {},
}  # fmt: skip


def main(argv: list[str] | None = None) -> None:
    """The vast CLI: one subcommand per shell shim (CLI table)."""
    ap = argparse.ArgumentParser(prog="vast", description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for cmd, args in CLI.items():
        p = sub.add_parser(cmd)
        for name, kw in args.items():
            p.add_argument(name, **kw)
    a = ap.parse_args(argv)
    try:
        globals()[f"cmd_{a.cmd}"](a)
    except InfraError as e:
        sys.exit(f"vast: {e}")


if __name__ == "__main__":
    main()
