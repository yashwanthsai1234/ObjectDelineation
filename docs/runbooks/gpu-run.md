# Runbook: a GPU run on Vast AI

**When to use this runbook.** You want to run a chain config (SFT, G-probe, RL, evaluation, a speed test) on a rented GPU.
**Who.** The repository owner, or an engineer with the owner's approval. A launch near the cost alarm needs the owner's written approval.
**Time.** From launch to the first job: about 30 to 45 minutes. The bundle upload and the image pull dominate (`r5`: launch 15:35, sync done 16:06, run started 16:13; the `r2` header notes 45 minutes).
**Rollback.** `teardown.sh RUN_ID` destroys the box. Nothing on the box is needed after a final pull.

Commands live in `experiments/phase2/vast/`. The `*.sh` files are symlinks to one shim, `vast.sh`. The link name is the command. All logic is in `policy.py`.
Run every command from the repository root on the Mac.

## 0. Prerequisites

| Need | Detail |
|---|---|
| `vastai` CLI | In `PATH`. Called with a scoped API key |
| Secrets | Names in `.secrets` (mode 600) or the macOS keychain: `objdel-vast-scoped` (scoped Vast key), `objdel-ghcr-read` (read-only GHCR token), `objdel-hf-write` (fine-grained Hugging Face token, write on one repo only). `objdel-riga-age` only for RIGA. Values are never printed |
| Image | `image` in `experiments/phase2/vast/.vast-allowlist`: `ghcr.io/yashwanthsai1234/objdel-gpu:vllm_native-3ef6e9179051@sha256:...` |
| Bundle | `build/bundles/bundle-lidc.tar` and its `.sha256`, from `experiments.phase2.vast.bundle lidc` ([Data, section 7](../reference/data.md#7-bundle-experimentsphase2vastbundle)) |
| `rsync`, `ssh` | On the Mac |
| Sweeper | Optional: the launchd job `com.objdel.vast.sweep.plist` runs `sweep_orphans.sh` every hour. Install with `sed "s#__REPO__#$PWD#g" FILE > ~/Library/LaunchAgents/...` and `launchctl load` |
| Budget | Read the cost rules in section 6 first |

## 1. Prove everything on the Mac

Do this last, after your final code change. Any later edit makes the stamp stale.

```bash
make test
make lint
make preflight CONFIG=experiments/phase2/configs/r1_qwen3vl2b.yaml
```

Pass `CONFIG` as the training config that names the data view, not the chain config.
Preflight stops at the first failing stage: DATA, then CHECKER, then CODE. It writes `runs/preflight/<config>.json`.
`launch.sh` refuses when the stamp is missing, failed, or stale (code, data or config changed).

## 2. Find an offer

```bash
experiments/phase2/vast/offers.sh RTX_5090 --pick top
```

- The filter is in `policy.terms`: US or CA host, one GPU, verified, reliability above 0.98, at least 12 effective cores, at least 32 GB of CPU RAM, CUDA 13.0 or newer, price ceiling per card, direct ports.
- It prints the top five. `--pick top` saves the cheapest admissible offer. `--pick ID` saves a specific one.
- Cards: `RTX_3090`, `RTX_4090`, `RTX_5090`, `RTX_PRO_6000_WS`, `RTX_PRO_6000_S`, `A100_SXM4`. RTX 5090 ceiling: $1.50 per hour.
- `--disk GB` sets the disk (default 80). `--riga` is for RIGA only. It needs recorded permission and a datacenter offer.

## 3. Launch

```bash
experiments/phase2/vast/launch.sh OFFER_ID RUN_ID --config experiments/phase2/configs/r1_qwen3vl2b.yaml \
    --cap-hours 3.5 [--owner-approved] [--idle-minutes 30] [--allow-parallel]
```

`RUN_ID` must match `[a-z0-9][a-z0-9_-]{0,39}`. It becomes the Vast label, the tmux name and the local folder name.

What `launch.sh` checks, in order: the run id is free, the preflight stamp is fresh, the offer file exists, the offer still passes the filter, then the cost rule (section 6).
It creates the instance with the allowlisted image, an on-start watchdog and SSH. It waits up to 900 s for `running` with a public IP and an SSH port.
Any failure or interrupt after `create` destroys the instance, found by its label if the reply was lost.

`--cap-hours` (default 4) sets the watchdog's time limit. The watchdog adds 600 s for the image pull.
Choose the cap from the estimate in the config header. Cover the sync, which can be slow through a proxy.

## 4. Sync, fetch, run

```bash
experiments/phase2/vast/sync_up.sh RUN_ID                       # code + bundle
experiments/phase2/vast/fetch_hf.sh RUN_ID r1-qwen3vl2b         # a checkpoint onto the box, as runs/RUN_ID/ckpt
experiments/phase2/vast/run.sh RUN_ID experiments/phase2/configs/rl_full.yaml
```

`sync_up.sh`:

1. Refuses RIGA unless permission is recorded and the offer is a datacenter one. Refuses anything under `data/raw`.
2. Sends the files that `git ls-files -co --exclude-standard` lists (tracked, or untracked and not ignored), without `data/`, `docs/`, the `lab` and `results` folders, or RIGA-derived tables (`policy.code_files`). It installs the package on the box and runs an import check.
3. Runs the GPU check (`image_check.sh`) before the bundle, so a bad card shows before a 12-minute upload.
4. Sends the bundle and verifies its SHA-256 and `SHA256SUMS`.
5. `--code-only` sends only the code, for a fix in mid-session.

`fetch_hf.sh` streams the HF token over SSH stdin. It refuses a token that is too broad.
The default repo is `svenk029/objdel-sft` (private). The revision is the branch name.

`run.sh RUN_ID CONFIG` starts a tmux session `<run_id>-<name>` on the box. `--cap-hours` goes before the config. The run's timeout is `min(cap, deadline - now - 600 s)`.
The chain writes `/workspace/runs/RUN_ID/{<name>.log, HEARTBEAT, EXIT}`. `EXIT` holds the exit code.

## 5. Monitor, pull, tear down

Look at a run:

```bash
ssh -p PORT root@IP                       # PORT and IP are printed by launch.sh and kept in active_runs.json
tmux attach -t RUN_ID-rl_full             # or: tail -f /workspace/runs/RUN_ID/rl_full.log
nvidia-smi -l 5
```

When `EXIT` exists:

```bash
experiments/phase2/vast/pull.sh RUN_ID                         # rsync runs/RUN_ID/ to the Mac, checksummed
experiments/phase2/vast/push_hf.sh RUN_ID experiments/phase2/configs/r1_qwen3vl2b.yaml   # optional: selected branch to HF
experiments/phase2/vast/teardown.sh RUN_ID
```

- `pull.sh` excludes `rolling/`, `decay/`, `stable/` (full checkpoints with optimizer, 29 GB for R1), and its own markers. It writes `PULLED_OK` only after rsync verifies every file.
- `teardown.sh` tries the final pull up to three times. It destroys the instance in every case, also on Ctrl-C. It checks that the instance is gone, saves the invoice JSON and writes the ledger rows. Vast bills late, so the invoice row may say `pending`.
- `push_hf.sh` refuses a non-private repo, a run trained on anything but LIDC, and a run with no selected branch. It never uploads RIGA data.

The watchdog on the box destroys it in three cases. The cap passes. The GPU idles after `EXIT` and the pull. No job ever started, after twice the idle time. It is a backstop. Do not rely on it.

## 6. Cost rules

| Rule | Value | Where |
|---|---|---|
| Estimate | `dph_total x hours` plus bandwidth: 15 GB down once, 6 GB up per hour, at the offer's $/GB | `policy.accrued`, `DOWN_GB`, `UP_GB_PER_H` |
| Alarm | A launch whose committed money plus its worst case reaches $10.00 needs `--owner-approved` | `ALARM = 10.0` |
| Hard stop | A launch is refused when committed money plus its worst case passes $17.50 | `HARD_STOP = 17.5` |
| Committed money | Closed sessions in the ledger plus every active run to its deadline | `spend(committed=True)` |
| Sweeper | Destroys unknown or overdue instances. At `HARD_STOP` spent it destroys all. At `ALARM` it alerts | `cmd_sweep_orphans` |
| One box | A second instance is refused without `--allow-parallel` | `cmd_launch` |

The hard stop rose with owner decisions: $15, then $16.60, $17.00 and $17.50 (commits `cf2f507`, `74a0135`, `edf3129`). The constants decide.

Recorded costs. The ledger value is the estimate at teardown (`experiments/phase2/vast/ledger.csv`). Vast bills late, so the invoice differs. Billed values come from the result notes where they exist:

| Run | What | Ledger estimate (USD) | Billed (USD) |
|---|---|---|---|
| p1 | RL readiness on RTX 5090 | 1.46 | 1.253 (`perf/p1_results.md`) |
| r2 | G-probe | 0.63 | 0.550 (`issues_solved_2026-11-07.md`) |
| r3 | GPU test of the scorer | 0.50 | 0.440 (same) |
| rl | Full RL run, final session | 4.56 | not recorded in the files read |
| r4 | Speed test | 0.65 | 0.65 (`r4_speed_2026-11-08.md`) |
| r5 | Shared-prefix test | 1.07 | not recorded in the files read |

Several earlier sessions of `rl` ended after $0.00 to $0.23 each. The one-layer image timed out on those hosts (see [Troubleshooting](troubleshooting.md)).

How to estimate a run: take seconds per step from the latest result. It is 27.5 s mean over the full run, and 23.6 s with `--shared-prefix`. Multiply by steps. Add about 30 minutes for checks, loading and saves. Multiply by the hourly price.
Example from `rl_full.yaml`: 770 steps x 36.3 s is 7.8 h, plus 0.5 h, at $0.641 per hour is about $5.3. The run took 27.5 s per step and billed $4.56.

## 7. Failure recovery

| Symptom | Action |
|---|---|
| `launch.sh`: `preflight gate: ...` | Run `make preflight CONFIG=...` again. Do not edit code afterwards |
| `launch.sh`: `hard stop` or `alarm` | Check the ledger. Ask the owner. Pass `--owner-approved` only with written approval |
| `launch.sh`: instance `exited`, `unknown`, `offline` or timeout | The instance is destroyed for you. Pick another offer. Hosts that fail on a paid box go into `BAD_MACHINES` |
| SSH refused | Vast hosts sometimes write `authorized_keys` with a foreign owner. The on-start command repairs it for 120 s. Retry `teardown.sh` (it retries the pull three times) |
| Image pull never finishes | A big layer cannot resume. The image is built with small layers for this reason. Pick another host |
| Chain step failed | Read `/workspace/runs/RUN_ID/<name>.log`. A soft step (`"~"`) logs "step failed, continuing". A hard failure stops the chain and writes `EXIT` |
| Box lost during RL | Launch a new box. `sync_up.sh`, then sync the old `train_rl/` folder (the pulled checkpoints) up. Run the same chain. `--resume` continues from the last `checkpoint-N`. At most 50 steps are lost |
| Cap too small | The timeout kills the chain before it writes final files. Launch with a larger `--cap-hours`. Include the sync time |
| Forgotten instance | `sweep_orphans.sh`. It destroys unknown or overdue instances and closes vanished runs in the ledger |
| Pull incomplete | Run `pull.sh RUN_ID` again before `teardown.sh`. A pull is idempotent |

## 8. Dry run

`DRY_RUN=1` prints every `ssh` and `rsync` command and does not run it. `vastai` is still called. Tests put a fake `vastai` on `PATH` (`tests/test_vast.py`).
Environment: `VAST_HOME` (ledger, active runs, offers), `VAST_RUNS` (pull target, default `runs`), `VAST_POLL_S`, `VAST_LAUNCH_TIMEOUT_S`.
