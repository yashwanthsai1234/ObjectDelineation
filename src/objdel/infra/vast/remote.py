"""The box: vastai, ssh and rsync calls (ssh and rsync printed under DRY_RUN=1), destroy, the onstart watchdog, the
code snapshot that sync_up sends with its CODE_SHA, and the secrets those calls stream."""

import base64
import gzip
import json
import os
import re
import shlex
import subprocess
import sys
import time

from objdel.core.provenance import repo_root
from objdel.infra.vast.ledger import accrued, active, home, ledger


def watchdog_script():
    """The onstart watchdog, from the code tree being launched (never from VAST_HOME, which holds state)."""
    return repo_root() / "experiments/phase2/vast/watchdog.sh"

RUN_ID_RE = re.compile(r"[a-z0-9][a-z0-9_-]{0,39}")  # Vast label, shell, tmux name, local path (BK4-16); onstart size
BAD_STATUS = ("exited", "unknown", "offline")  # never reach running (V19)
# sync_up runs it on the box after the install: library modules and the box stack only, no tests.* or experiments.*
IMPORT_CHECK = ("import objdel.geom.engine, objdel.verify, objdel.rl.trl_adapter, objdel.policy.sampler, pydantic, "
                "vllm")  # vllm and the sampler: the vllm_native image, one venv (ADR 0023)
# The one list of live CLI modules (E8 step 2); the CLI smoke and import tests and the module-string gate read it
LIVE_CLIS = tuple(f"experiments.phase1.{m}" for m in ("build_data", "build_targets", "fetch_lidc_pixels",
                                                      "prepare_rl")) + tuple(
    f"experiments.phase2.{m}" for m in ("build_views", "build_eval_views", "train_sft", "select_ckpt", "train_rl",
                                        "g_probe", "final_eval", "preflight", "hf_push", "hf_pull"))
# MARGIN_S: run.sh's timeout ends this long before the watchdog cap, leaving time for the last pull. GRACE_S: sweep
# never closes a run younger than this as vanished (it may not be listed yet). PULL_S: added to --cap-hours, as the
# watchdog counts from launch; V0's 5.88 GB image pull took 5.7 min.
MARGIN_S = GRACE_S = PULL_S = 600
# Not pulled, and ignored by the watchdog's "nothing newer than PULLED_OK" (it gets this list in the onstart env):
# full checkpoints with optimizer state (R1: 29 GB; branches/ is model-only and is pulled) and the pull's own markers.
PULL_EXCLUDE = ("rolling/", "decay/", "stable/", "PULLED_OK", "PULLED_OK.next")
ONSTART_MAX = 4048  # Vast API spec


def die(msg: str):
    """Stop the CLI with a `vast: ` message (exit status 1)."""
    sys.exit(f"vast: {msg}")


def secret(service: str) -> str:
    """objdel.infra.secrets: <repo>/.secrets, else the macOS keychain; dies when neither has it."""
    from objdel.infra import secrets

    try:
        value = secrets.get(service)
    except PermissionError as e:
        die(str(e))
    return value or die(f"secret {service} missing: not in {secrets.path()}, not in the keychain")


def vast(*args: str, tries: int = 1):
    """Run `vastai ARGS` with the scoped API key, `tries` times; JSON for --raw, else the text; dies on failure."""
    env = dict(os.environ, VAST_API_KEY=os.environ.get("VAST_API_KEY") or secret("objdel-vast-scoped"))
    for i in range(tries):
        p = subprocess.run(["vastai", *args], capture_output=True, text=True, env=env, stdin=subprocess.DEVNULL)
        # an API error exits 0 with empty stdout (V0)
        if not p.returncode and (p.stdout.strip() or "--raw" not in args):
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
    """The ssh argv prefix for a run's box (its port, our known_hosts under home())."""
    known = f"UserKnownHostsFile={home() / 'known_hosts'}"
    return ["ssh", "-p", str(run["port"]), "-o", "StrictHostKeyChecking=accept-new", "-o", known]


def remote(run: dict, command: str, stdin: str | None = None, what: str = "") -> None:
    """Run `command` on a run's box over ssh (DRY_RUN=1 prints it)."""
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
    # gzip: Debian essential
    b64 = base64.b64encode(gzip.compress(watchdog_script().read_bytes(), mtime=0)).decode()
    created = int(created or time.time())
    env = (f"RUN_ID={run_id} WATCHDOG_CAP_S={cap_s} WATCHDOG_IDLE_S={idle_s} WATCHDOG_CREATED={created} "
           f"WATCHDOG_PULLED_OK=/workspace/runs/{run_id}/PULLED_OK WATCHDOG_SKIP={','.join(PULL_EXCLUDE)}")  # fmt: skip
    k = "/root/.ssh/authorized_keys"  # V1: some hosts write it with a foreign owner and sshd refuses it (vast-cli #336)
    fix = (f"(for i in $(seq 60); do [ -f {k} ] && chown root:root {k} && chmod 600 {k} && chmod 700 /root/.ssh "
           f"&& chmod go-w /root; sleep 2; done) & ")
    tmux = "touch /root/.no_auto_tmux; "  # V1: Vast's login tmux broke every `ssh host CMD` ("not a terminal")
    cmd = tmux + fix + (f"printf %s {b64} | base64 -d | gunzip > /root/watchdog.sh && {env} "
                        "nohup bash /root/watchdog.sh >> /root/watchdog.log 2>&1 &")
    if len(cmd) > ONSTART_MAX:
        die(f"onstart is {len(cmd)} characters, over Vast's {ONSTART_MAX} (shorten watchdog.sh or the run_id)")
    return cmd


def get_run(run_id: str) -> dict:
    """The active run `run_id`; dies when there is none."""
    return active().get(run_id) or die(f"no active run {run_id}")


def code_files() -> list[str]:
    """The code snapshot (tracked and new files), without the roots data/, docs/, archive/ and .claude/, lab or results
    dirs, or RIGA-derived tables (D16, BK4-13); only files on disk (rsync --files-from aborts on a missing one).
    tests/ ships: the chain runs tests/verify/test_gpu_scoring.py on the box. Kept from docs/: the two records the box
    reads
    (final_eval's frozen view hashes, val_half's table; LIDC ids, counts and hashes, no pixels, no RIGA; B1.6) and the
    kernel drafts K1-K3 without their build/ and results/ folders (the R-D1 runbooks run them on the box)."""
    git = ["git", "ls-files", "-z", "-co", "--exclude-standard"]  # -z: paths with spaces stay whole (V1)
    files = subprocess.run(git, cwd=repo_root(), capture_output=True, text=True, check=True).stdout.split("\0")[:-1]
    riga = re.compile(r"riga[^/]*\.(csv|jsonl?)$", re.I)
    drop = re.compile(r"^(data|docs|archive|\.claude)/|^experiments/[^/]+/(lab|results)/", re.I)
    keep = re.compile(r"^docs/reference/(frozen_views|val_halves)\.json$"
                      r"|^docs/implementation2\.0/kernels/K[123]/(?!build/|results/)")
    files = [f for f in dict.fromkeys(files)
             if not riga.search(f) and (keep.match(f) or not drop.search(f)) and (repo_root() / f).is_file()]
    if "tests/verify/test_gpu_scoring.py" not in files:  # the chain runs it on the box before any scoring stage
        die("tests/verify/test_gpu_scoring.py missing from the code snapshot")
    return files


def alert(msg: str) -> None:
    """Print an alert to stderr (and a macOS notification unless DRY_RUN=1)."""
    print(f"vast ALERT: {msg}", file=sys.stderr)
    if sys.platform == "darwin" and os.environ.get("DRY_RUN") != "1":
        subprocess.run(["osascript", "-e", f'display notification {json.dumps(msg)} with title "objdel vast"'])
