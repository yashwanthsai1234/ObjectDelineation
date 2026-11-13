"""The box: vastai, ssh and rsync calls, destroy, the onstart watchdog, the code snapshot that sync_up sends with its
CODE_SHA, and the secrets those calls stream. Under DRY_RUN=1 ssh, rsync and every vastai call that rents or destroys
are printed, never run; vastai calls that only read still run."""

import base64
import gzip
import json
import os
import random
import re
import shlex
import subprocess
import sys
import time

from objdel.core.io import git_files
from objdel.core.provenance import repo_root
from objdel.core.spec import CLI
from objdel.infra.errors import InfraError, redact
from objdel.infra.secrets import require as secret  # noqa: F401  (re-exported: vast.__main__ and the tests)
from objdel.infra.vast.ledger import accrued, active, home, ledger, locked, update_active


def watchdog_script():
    """The onstart watchdog, from the code tree being launched (never from VAST_HOME, which holds state)."""
    return repo_root() / "ops/vast/watchdog.sh"

RUN_ID_RE = re.compile(r"[a-z0-9][a-z0-9_-]{0,39}")  # Vast label, shell, tmux name, local path; onstart size
BAD_STATUS = ("exited", "unknown", "offline")  # never reach running
# sync_up runs it on the box after the install: library modules and the box stack only, no tests.* or tools.*
IMPORT_CHECK = ("import objdel.geom.engine, objdel.verify, objdel.rl.trl_adapter, objdel.policy.sampler, pydantic, "
                "vllm")  # vllm and the sampler: the vllm_native image, one venv 
# The one list of live CLI modules; the CLI smoke and import tests and the module-string gate read it
LIVE_CLIS = tuple(f"{CLI}.{m}" for m in ("data.store", "data.targets", "data.fetch_lidc", "data.prepare", "data.views",
                                         "data.eval_views", "train.sft", "train.select", "train.rl", "train.probe",
                                         "eval.run", "box.preflight", "hub.push", "hub.pull"))
# MARGIN_S: run.sh's timeout ends this long before the watchdog cap, leaving time for the last pull. GRACE_S: sweep
# never closes a run younger than this as vanished (it may not be listed yet). PULL_S: added to --cap-hours, as the
# watchdog counts from launch; a measured 5.88 GB image pull took 5.7 min.
MARGIN_S = GRACE_S = PULL_S = 600
# Not pulled, and ignored by the watchdog's "nothing newer than PULLED_OK" (it gets this list in the onstart env):
# full checkpoints with optimizer state (29 GB; branches/ is model-only and is pulled) and the pull's own markers.
PULL_EXCLUDE = ("rolling/", "decay/", "stable/", "ckpt/", "PULLED_OK", "PULLED_OK.next")
# File names that never ship in the code snapshot, tracked or not
SECRET_NAME = re.compile(r"^(\.secrets|\.env|id_rsa|id_ed25519|\.netrc$)|\.(pem|key)$", re.I)
# Box paths and the box python, one place (the box image fixes them)
BOX_SRC = "/workspace/src"
BOX_PY = "/opt/venv/bin/python -I"
# The box sends nothing to vLLM, Hugging Face or any DO_NOT_TRACK-aware library
BOX_ENV = "VLLM_NO_USAGE_STATS=1 DO_NOT_TRACK=1 HF_HUB_DISABLE_TELEMETRY=1"
ONSTART_MAX = 4048  # Vast API spec
VAST_TIMEOUT_S = 120  # one `vastai` call: its API answers in seconds; a longer wait is a hang
# One ssh or rsync call. rsync also stops on 300 s without data (--timeout), so this only bounds a transfer that
# trickles on: 1 h moves a 3 GB model-only branch at 0.85 MB/s, and teardown makes up to 3 pulls.
SH_TIMEOUT_S = 3600
# ForwardAgent=no: a hostile box never gets our ssh-agent. No IdentitiesOnly: no key is pinned, the agent holds it.
SSH_OPTS = ("StrictHostKeyChecking=accept-new", "ForwardAgent=no", "ConnectTimeout=30", "ServerAliveInterval=30",
            "ServerAliveCountMax=4")
BACKOFF_CAP_S = 60
SPENDING = ("create", "destroy")  # the vastai verbs that rent or destroy an instance


def box_run_dir(run_id: str) -> str:
    """A run's folder on the box."""
    return f"/workspace/runs/{run_id}"


def dry_run() -> bool:
    """True under DRY_RUN=1."""
    return os.environ.get("DRY_RUN") == "1"


def backoff_s(i: int) -> float:
    """The wait after failed try i: VAST_POLL_S x 2**i, at most BACKOFF_CAP_S, times a random 0.5 to 1.5 (no two
    clients retry in step)."""
    return min(BACKOFF_CAP_S, float(os.environ.get("VAST_POLL_S", 10)) * 2**i) * random.uniform(0.5, 1.5)


def vast(*args: str, tries: int = 1):
    """Run `vastai ARGS` with the scoped API key, up to `tries` times (backoff_s between them); JSON for --raw, else
    the text; InfraError on failure. A call over VAST_TIMEOUT_S is a failed try. DRY_RUN=1 prints a spending call and
    returns {} (--raw) or "" without running it."""
    if dry_run() and args[0] in SPENDING:
        shown = ["<login>" if j and args[j - 1] == "--login" else a for j, a in enumerate(args)]  # the GHCR token
        print("DRY_RUN: vastai", shlex.join(shown))
        return {} if "--raw" in args else ""
    env = dict(os.environ, VAST_API_KEY=os.environ.get("VAST_API_KEY") or secret("objdel-vast-scoped"))
    for i in range(tries):
        try:
            p = subprocess.run(["vastai", *args], capture_output=True, text=True, env=env, stdin=subprocess.DEVNULL,
                               timeout=VAST_TIMEOUT_S)
        except subprocess.TimeoutExpired:
            why = f"no answer in {VAST_TIMEOUT_S} s"
        else:
            # an API error exits 0 with empty stdout
            if not p.returncode and (p.stdout.strip() or "--raw" not in args):
                return json.loads(p.stdout) if "--raw" in args else p.stdout
            why = redact(p.stderr.strip()[-300:])
        if i < tries - 1:
            time.sleep(backoff_s(i))
    raise InfraError(f"vastai {args[:2]} failed: {why}")


def sh(cmd: list[str], stdin: str | None = None, what: str = "") -> None:
    """ssh / rsync: printed under DRY_RUN=1 (stdin shown only as a label), otherwise run for at most SH_TIMEOUT_S; a
    failure or a timeout stops."""
    if dry_run():
        print("DRY_RUN:", shlex.join(cmd), f"<stdin: {what}>" if stdin is not None else "")
        return
    try:
        p = subprocess.run(cmd, input=stdin, text=True, stdin=None if stdin else subprocess.DEVNULL,
                           timeout=SH_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        raise InfraError(f"{cmd[0]} timed out after {SH_TIMEOUT_S} s: {redact(shlex.join(cmd[-2:]))[:160]}") from None
    if p.returncode:
        raise InfraError(f"{cmd[0]} failed: {redact(shlex.join(cmd[-2:]))[:160]}")


def ssh(run: dict) -> list[str]:
    """The ssh argv prefix for a run's box (its port, our known_hosts under home(), connect and keep-alive limits)."""
    opts = [*SSH_OPTS, f"UserKnownHostsFile={home() / 'known_hosts'}"]
    return ["ssh", "-p", str(run["port"]), *(x for o in opts for x in ("-o", o))]


def remote(run: dict, command: str, stdin: str | None = None, what: str = "") -> None:
    """Run `command` on a run's box over ssh (DRY_RUN=1 prints it)."""
    sh([*ssh(run), f"root@{run['ip']}", command], stdin, what)


def rsync(run: dict, *paths: str, opts: tuple = ()) -> None:
    """paths: sources, then the destination (two sources once made the .sha256 the destination). Stops after 300 s
    without data; a partly sent file is kept for the next try."""
    # --safe-links: a link that points outside the copied tree (planted by a host, or ours) is skipped, never followed
    sh(["rsync", "-az", "--timeout=300", "--partial", "--safe-links", "--exclude=data/raw/", *opts, "-e",
        shlex.join(ssh(run)), *paths])


def destroy(iid: int, run_id: str, event: str, inst: dict | None = None) -> None:
    """Destroy (never stop), verify the instance is no longer listed, record the session's cost (Vast's dph and
    start when it lists the instance, else the run record's). Six tries for the destroy itself: a box left running
    bills until the watchdog or the next sweep."""
    inst = inst or next((i for i in vast("show", "instances", "--raw", tries=3) if i["id"] == iid), {})
    run = next((r for r in active().values() if r["instance"] == iid), {})
    r = dict(run, **({"dph": inst.get("dph_total", 0), "started": inst["start_date"]} if "start_date" in inst else {}))
    cost = accrued(r, time.time()) if r else 0.0
    vast("destroy", "instance", str(iid), "-y", tries=6)
    if dry_run():
        return
    for _ in range(10):
        if all(i["id"] != iid for i in vast("show", "instances", "--raw", tries=3)):
            break
        time.sleep(float(os.environ.get("VAST_POLL_S", 10)))
    else:
        raise InfraError(f"instance {iid} still listed after destroy: destroy it in the console NOW")
    with locked():
        update_active(lambda runs: {k: x for k, x in runs.items() if x["instance"] != iid})
        ledger(event=event, run_id=run_id, instance_id=iid, cost_usd=round(cost, 4))


def onstart_cmd(run_id: str, cap_s: int, created: float | None = None, idle_s: int = 1800) -> str:
    """--onstart-cmd: unpack the watchdog and start it; carries no secret (CONTAINER_API_KEY is injected by Vast).
    Vast's API spec limits onstart to ONSTART_MAX characters."""
    # gzip: Debian essential
    b64 = base64.b64encode(gzip.compress(watchdog_script().read_bytes(), mtime=0)).decode()
    created = int(created or time.time())
    env = (f"RUN_ID={run_id} WATCHDOG_CAP_S={cap_s} WATCHDOG_IDLE_S={idle_s} WATCHDOG_CREATED={created} "
           f"WATCHDOG_PULLED_OK={box_run_dir(run_id)}/PULLED_OK WATCHDOG_SKIP={','.join(PULL_EXCLUDE)}")  # fmt: skip
    k = "/root/.ssh/authorized_keys"  # some hosts write it with a foreign owner and sshd refuses it (vast-cli #336)
    fix = (f"(for i in $(seq 60); do [ -f {k} ] && chown root:root {k} && chmod 600 {k} && chmod 700 /root/.ssh "
           f"&& chmod go-w /root; sleep 2; done) & ")
    tmux = "touch /root/.no_auto_tmux; "  # Vast's login tmux broke every `ssh host CMD` ("not a terminal")
    cmd = tmux + fix + (f"printf %s {b64} | base64 -d | gunzip > /root/watchdog.sh && {env} "
                        "nohup bash /root/watchdog.sh >> /root/watchdog.log 2>&1 &")
    if len(cmd) > ONSTART_MAX:
        raise InfraError(f"onstart is {len(cmd)} characters, over Vast's {ONSTART_MAX} "
                         "(shorten watchdog.sh or the run_id)")
    return cmd


def get_run(run_id: str) -> dict:
    """The active run `run_id`; InfraError when there is none."""
    if (run := active().get(run_id)) is None:
        raise InfraError(f"no active run {run_id}")
    return run


def code_files() -> list[str]:
    """The code snapshot (tracked and new files, but never a secret-shaped name, an untracked top-level .md or a link
    out of the tree), without the roots data/, docs/, archive/ and .claude/, the checker's
    raw results (results/checker/; its record results/checker.md ships, the checker eval writes into it), or
    RIGA-derived tables; only files on disk (rsync --files-from aborts on a missing one).
    tests/ ships: the chain runs tests/verify/test_gpu_scoring.py on the box. The records the box reads (final_eval's
    frozen view hashes, val_half's table; LIDC ids, counts and hashes, no pixels, no RIGA) are in reference/.
    docs/ is local only: nothing in it ships."""
    root = repo_root()
    files = git_files(root, "-co", "--exclude-standard")  # NUL split: paths with spaces stay whole
    tracked = set(git_files(root))
    riga = re.compile(r"riga[^/]*\.(csv|jsonl?)$", re.I)
    drop = re.compile(r"^(data|docs|archive|\.claude)/|^results/checker/|^ops/vast/\.lock$", re.I)
    # never a secret-shaped name, and never an untracked note at the top level (HANDOVER.md reached a box once)
    files = [f for f in dict.fromkeys(files)
             if not riga.search(f) and not drop.search(f) and not SECRET_NAME.search(f.rsplit("/", 1)[-1])
             and not (f not in tracked and "/" not in f and f.endswith(".md")) and (root / f).is_file()
             and (root / f).resolve().is_relative_to(root.resolve())]  # a link out of the tree never ships
    if "tests/verify/test_gpu_scoring.py" not in files:  # the chain runs it on the box before any scoring stage
        raise InfraError("tests/verify/test_gpu_scoring.py missing from the code snapshot")
    return files


def alert(msg: str) -> None:
    """Print an alert to stderr (and a macOS notification unless DRY_RUN=1)."""
    print(f"vast ALERT: {msg}", file=sys.stderr)
    if sys.platform == "darwin" and not dry_run():
        subprocess.run(["osascript", "-e", f'display notification {json.dumps(msg)} with title "objdel vast"'])
