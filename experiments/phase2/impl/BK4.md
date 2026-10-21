# BK4 break-it findings

Agent BK4, wave 2, 2026-10-20. Scope: Vast operations (`experiments/phase2/vast/`, `docker/`) and a CI simulation.
No Vast API call, no money spent, no image built or pushed.

Tests: `tests/test_p2_breakit_ops.py`, 341 lines. Results:

- Mac `.venv`: 30 passed, 22 strict xfails.
- `[dev]`-only CI venv: also 30 passed, 22 xfailed.
- `ruff check` and `ruff format --check` are clean on the file.

The ops tests reuse `tests/test_p2_vast.py`'s fakes and add two knobs:

- a transient `vastai` failure;
- a `running` reply without ports.

The CI tests block `torch`, `transformers`, `trl`, `datasets`, `peft`, `accelerate`, `pydicom` and `pylidc` through
`sys.modules`, so the `[dev]`-only job is reproduced on any machine.

## B. CI simulation

**Run 1.** `uv venv --python 3.12` + `-e ".[dev]"` in a scratch venv, then `OBJDEL_FIXTURES=/nonexistent
OBJDEL_STORE=/nonexistent pytest -q -m "not slow"` from the repo root.
- Result: 811 passed, 96 skipped, **3 failed, 4 errors** (134 s).

**Run 2.** A copy of the `git ls-files -co` snapshot without `data/` (the real CI has no `data/`), with a fresh `[dev]`
venv.
- `ruff check src tests experiments`: **clean**.
- pytest: 798 passed, 105 skipped, 7 failed, 4 errors.
- 5 of those failures are `test_p2_vast` tests. They fail only because the copy was not a git repo
  (`git ls-files` exit 128). After `git init`, all 47 `test_p2_vast` tests passed. CI's checkout is a git repo, so
  these are not CI breaks.

**CI breaks found:**

| Test | Cause | Finding |
|---|---|---|
| `tests/test_p2_bridge.py` (4 errors) | the fixture imports `datasets` without a skip | BK4-1 |
| `tests/test_p2_views.py::test_build_eval_views_help` | `build_eval_views` imports `datasets` at module level | BK4-2 |
| `tests/test_p2_eval.py::test_nsd_equals_the_deepmind_reference` | known issue S1-15 | — |

Root cause of S1-15, added here (not a new id):
- `nsd` returns -2.2e-16 at τ = 0 (Hypothesis falsifying example `rng=np.random.default_rng(703608005), side=39,
  tau=0.0`). This is floating-point rounding, not test order.
- Hypothesis draws new seeds on every run, so it breaks CI at random.
- Fix: clip the result to [0, 1], or compute the overlap as an integer ratio.

Not simulated within the time box: the `train` job (`[dev,data,train]`, CPU torch, `--cov-fail-under=90`).

## A. Findings

Severity "high" means a money, licence or security risk, or a CI break.

| id | sev | component | repro (test) | expected vs actual | suggested fix | owner |
|---|---|---|---|---|---|---|
| BK4-1 | **high** (CI) | tests / M2 | `test_ci_dev_only_job_skips_instead_of_failing[test_p2_bridge.py]`; CI `[dev]` job | expected: skip without `datasets`. Actual: 4 setup errors (`rows` fixture → `objdel.train.dataset` → `from datasets import …`, `dataset.py:64`) | `pytest.importorskip("datasets")` at the top of `test_p2_bridge.py` | T |
| BK4-2 | **high** (CI) | tests / M1 | `…[test_p2_views.py::test_build_eval_views_help]` | expected: skip. Actual: `ModuleNotFoundError: datasets` (`build_eval_views.py:34` imports `objdel.train.dataset` at module level) | `importorskip("datasets")` in the test, or a lazy import in `build_eval_views` (see S1-2) | V |
| BK4-3 | **high** (money) | M7 `policy.cmd_launch` | `test_launch_destroys_when_polling_fails` | expected: any failure after `create instance` destroys the instance (G8). Actual: one transient `vastai show instance` error ends in `die()` (`policy.py:177`) with the instance alive. It is in `active_runs`, so the sweep treats it as neither an orphan nor overdue, and it bills until the watchdog cap (default 4 h; PRO 6000 ≈ $10). Ctrl-C during the poll has the same effect | wrap everything after `create` (`policy.py:276-294`) in `try/except BaseException: destroy(iid, …); raise`; retry `show instance` a few times before giving up | O |
| BK4-4 | **high** (money) | M7 `policy.cmd_launch` | `test_launch_destroys_when_running_reply_lacks_ports` | expected: destroy, or keep polling. Actual: a `running` reply without `public_ipaddr` / `ports["22/tcp"]` raises a `KeyError` (`policy.py:287`); traceback, instance alive (same consequence as BK4-3). Vast often reports `running` before the port map is filled | treat a missing port map as "not running yet" (`inst.get("ports") or {}`), plus the fix for BK4-3 | O |
| BK4-5 | **high** (data loss / money) | M8 configs × M7 `run`/`pull` | `test_config_outputs_land_where_pull_reads` | expected: train / select / probe outputs end up under `/workspace/runs/<RUN_ID>/`, which `pull.sh` copies and the watchdog's `PULLED_OK` refers to. Actual: `run.sh` runs `cd /workspace/src` (`policy.py:363`) and every config has a relative `output_dir: runs/<name>` (`configs/*.yaml`). Checkpoints land in `/workspace/src/runs/…`, are **never pulled**, and are lost when teardown or the idle watchdog destroys the box. `pull.sh` reports success because `run.log` / `EXIT` / `HEARTBEAT` arrived | `cmd_run` passes the run dir explicitly (e.g. `cd /workspace/runs/<id>` and `--output_dir .`, or `OBJDEL_RUN_DIR`), or `pull` also syncs `/workspace/src/runs/`. Add a DRY_RUN test joining a config to `run.sh` | R / O |
| BK4-6 | medium | M8 configs × M7 `sync_up` | `test_config_data_root_is_the_box_bundle_dir` | expected: the data root on the box is `/workspace/data/lidc` (where `sync_up` extracts the bundle, `policy.py:347`). Actual: no config sets `data_root`. The default `data` resolves to `/workspace/src/data/…`, which does not exist. Same for the defaults of `k3_harness --root data/prepared` and of `select_ckpt` / `g_probe`. It fails fast, but each failure costs minutes of GPU time | set `data_root: /workspace/data/lidc` in the configs, or have `sync_up` symlink `/workspace/src/data` → `/workspace/data/lidc` (keep the `data/raw` rule) | R / O |
| BK4-7 | medium (money) | M7 `spend` / `sweep` | `test_sweep_reconciles_a_run_whose_instance_vanished` | expected: when an active run's instance is no longer listed (watchdog self-destroy G5, host loss), sweep records its cost and drops it. Actual: `active_runs.json` keeps it forever and `spend()` accrues `dph × (now - started)` without bound (`policy.py:140`). The real session cost never reaches the ledger. After a day the phantom accrual passes $15 and every launch is refused. `destroy()` with an unlisted instance also records cost 0 (`policy.py:207-208`) | in `cmd_sweep`, reconcile active runs whose instance is absent: cost = dph × (min(now, deadline) − started), event `vanished`, remove from `active_runs` | O |
| BK4-8 | medium (money) | M7 `cmd_launch` G10 | `test_parallel_launches_cannot_jointly_pass_the_hard_stop` | expected: the hard-stop projection counts the remaining worst case of every live run. Actual: it counts only their accrual so far (`policy.py:263`). Three `--allow-parallel` 12 h launches at $0.55/h are all admitted: $19.80 worst case > $15. The hourly sweep catches it ≤ 1 h after accrual passes $15, if the Mac is awake | `worst += Σ r.dph × max(0, r.deadline − now)` over `active()` | O |
| BK4-9 | medium (money / policy) | M7 `cmd_launch` G10 | `test_launch_projected_past_the_alarm_needs_approval` | PLAN §4.6: "passing either [$10 / $15] needs the owner's written approval". Actual: the alarm compares only money already spent (`policy.py:266`). With $9 spent, a 4 h launch (projected $11.20) needs no `--owner-approved`. **Coordinator decision**: if "alarm = spent only" is intended, record it in PLAN and delete this test | `if now_usd + worst >= ALARM and not a.owner_approved` | O (coordinator decides) |
| BK4-10 | medium | M7 `cmd_launch` | `test_launch_refuses_an_active_run_id` | expected: a run_id that is still active is refused. Actual: `active({**active(), run_id: run})` overwrites the first record (`policy.py:279`). The first instance's accrual disappears from `spend()`, and the next sweep destroys the first box as an "orphan" in the middle of its run | refuse when `a.run_id in active()` | O |
| BK4-11 | medium (data loss) | M7 `watchdog.sh` | `test_watchdog_keeps_box_when_outputs_are_newer_than_the_pull` | expected: the idle self-destroy happens only when everything the run wrote is on the Mac. Actual: `watchdog.sh:29` only requires `PULLED_OK` to be newer than the start of the GPU-idle period. CPU stages (k3_harness, scoring after vLLM sampling, `select_ckpt` writes, HF downloads) can keep writing for > 30 min of GPU idle. A launchd pull inside that window lets the watchdog destroy the box with outputs that were never pulled, or in the middle of a job | also require `/workspace/runs/<id>/EXIT` to exist and no file under the run dir newer than `PULLED_OK` (`find "$RUN" -newer "$PULLED" -print -quit` is empty). Together with BK4-5, the run dir must be where the outputs are | O |
| BK4-12 | low | M7 `reasons` | `test_offer_recheck_covers_the_whole_query[*]` (5) | the docstring says "re-checked on the returned JSON (defence in depth)". Actual: reliability, rentable, inet_down, direct_port_count and disk_space (also duration) are not re-checked (`policy.py:91-101`). Also, `reasons` converts with MB/1024 while the query uses Vast's units (UNVERIFIED: the CLI may multiply `gpu_ram`/`cpu_ram` by 1000). At 79 GB that is a 2.4 % difference | generate `query` and `reasons` from one table of (term, predicate) (see Simplicity) | O |
| BK4-13 | low (licence) | M7 `code_files` | `test_code_snapshot_has_no_riga_derived_results` | D16: no RIGA-derived data on community hosts. Actual: the code rsync ships these tracked files: `experiments/phase0/results/{riga,riga_outlines,medsam_riga.csv}`, `phase1/lab/t_data/riga/*.csv` (per-image areas, topology, vCDR) and `phase0/results/vlm_polygon/replies.jsonl` (VLM polygons predicted on RIGA images). They are statistics and model outputs, not crops or masks, but they serve no purpose on the box | `code_files()` keeps only `src/ experiments/*.py experiments/phase*/{configs,vast,docker} tests/ pyproject.toml CMakeLists.txt cpp/ …` (an allowlist), or excludes `*/results/`, `*/lab/`, `docs/` | O |
| BK4-14 | medium (licence) | M7 `bundle.walk` | `test_bundle_refuses_a_symlinked_view_directory` | docstring: "Symlinks … are refused". Actual: only entries *inside* a walked directory are checked (`bundle.py:34-38`). A symlinked top-level directory (`data/views/lidc-x → riga-train`, or `prepared/lidc` → elsewhere) is followed, and RIGA files enter `bundle-lidc.tar`, which may go to community hosts. `is_raw` catches only `data/raw` targets | in `walk`, refuse `top.is_symlink()`, or any `p.resolve()` outside `data.resolve()/<kind>/<name>` | O |
| BK4-15 | low | M7 `Dockerfile.phase2` | `test_dockerfile_extras_are_pinned` | D15 / §4.3: rebuild only when pins change, tag `p2-<pinhash>`. Actual: `flash-linear-attention`, `kernels` and `causal-conv1d` are unpinned (`Dockerfile.phase2:37-38`), so two builds with the same tag can differ (Triton kernel route, numerics). `/opt/venv.freeze.txt` records the versions only after the fact | pin them with `==` in an `extras.txt` that is part of the pin hash (`build_images.sh:40`) | O |
| BK4-16 | low (security hygiene) | M7 `run_id` | `test_launch_validates_run_id[*]` (3) | expected: run_id restricted to `[a-z0-9][a-z0-9_-]*`. Actual: anything is accepted. The run_id goes into Vast `--label` and, unquoted, into the `bash -c` text on the box (`mkdir -p /workspace/runs/{run_id}`, `policy.py:361-367`). It is also the tmux session name (tmux rejects `.` and `:`) and a local path (`VAST_RUNS/<id>`, `../` traversal) | `re.fullmatch` in `cmd_launch`; `shlex.quote` the run dir | O |
| BK4-17 | medium (money) | M7 `watchdog.sh` | no test (needs a state path the fix must choose) | expected: the cap counts from the first start. Actual: `start=$(date +%s)` (`watchdog.sh:21`). Vast re-runs the onstart command on every container restart (host reboot, Vast restart), so the cap resets each time and `/etc/environment` gets duplicate lines. The Mac sweep (hourly, Mac awake) is then the only cap. Also, the cap starts after the image pull (5–15 min for 7.5 GB) while the Mac deadline starts at `create` | persist the first start: `[ -f /root/watchdog.start ] \|\| date +%s > /root/watchdog.start; start=$(cat …)`; pass `WATCHDOG_CREATED=<epoch>` from `onstart_cmd` and use the minimum | O |
| BK4-18 | medium (money) | M7 filter / `spend` / teardown | no test (PLAN gap) | G1-G10 cover only `dph_total` × time. Bandwidth is charged separately per host (`inet_up_usd_tb_median` up to $22.48/TB in the A3 snapshot, and hosts set their own price). Pulls of 3–5 GB checkpoints every 20 min are not in the filter, not in `spend()`, and not reconciled: `invoices-v1` is saved verbatim (`policy.py:388-390`) and never parsed | add `inet_up_cost` / `inet_down_cost` ceilings to the filter and record them in the ledger; at teardown, parse the invoice and write the actual charge row (rows with `detail=invoice`) so `spend()` uses real money | O (after V0, see S1-5) |
| BK4-19 | info | M7 RIGA decrypt | no test (UNVERIFIED) | `/dev/shm` in Docker defaults to 64 MB; the RIGA plaintext is ≈ 280 MB (§4.4). If Vast does not enlarge shm, `age -d … \| tar -x -C /dev/shm/riga` fails. It fails safe: the `sha256sum -c` check still fails the step | check `df -h /dev/shm` in V0 / `image_check.sh` | O |
| BK4-20 | low (local CI only) | tests (phase 1) | `test_data_lidc.py::test_parity_with_pylidc_database_sample` | fails in a `[dev]`-only env when `data/raw/lidc_xml` exists (pylidc absent → `find_spec` returns None → `AttributeError`, `lidc.py:286`). It skips in real CI (no `data/`) | `pytest.importorskip("pylidc")` (find_spec only) in the test | any |

**Guards that pass today** (kept in the file):
- G10 projection for a single launch.
- Launch with a `None` status still destroys at the timeout.
- Ledger blank or partial cost rows.
- H-GOV permission × hosting matrix (18 cases).
- Per-card price / gpu_ram / cpu_ram / datacenter boundaries for all 6 cards: exact ceiling admitted, +$0.001 refused; `hi×1024` admitted, +1 MB refused; A100 `79×1024 − 1` refused.
- No `stop` verb anywhere: E2E + sweep hard stop + source; the watchdog uses `DELETE`.
- RIGA identity, across sync → run → pull → teardown → sweep with `DRY_RUN=0`: carried only on the one `age -d -i /dev/stdin` ssh stdin; never in argv, env, script output, ledger, offers or state files.
- `image_check.sh` peak lookup per card name (4090 → 165.2, A100-SXM4-80GB → 312, unknown → fail) and `tf >= 0.6 * peak`.
- The Mac freeze carries `scikit-build-core==1.1.1` and `nanobind==3.1.0` (needed by the box's `--no-build-isolation -e .`) and the D15 pins (torch 2.14.1, transformers 5.20.1, TRL 1.16.0, numpy 2.5.4, numba 0.68.0).

**Checked, no defect:**
- **Image versions vs CI.** The image is the Mac freeze by construction (`build_images.sh:28`). CI's `[dev]` resolve gave the same numpy 2.5.4, numba 0.68.0, scipy 1.18.1, opencv 5.0.0.93, pillow 13.0.0 and pyarrow 26.0.0, and differs only in hypothesis (6.170.5 vs 6.170.2). Caveat: `numpy>=2.0` / `scipy>=1.13` are unpinned in `pyproject.toml`, so CI drifts from the image whenever they release.
- **Box imports.** `cd / && .venv/bin/python -I -c "import experiments.phase2.{k3_harness,train_sft,select_ckpt,g_probe}"` succeeds through `_objdel_repo.pth`, the same mechanism `/workspace/src` uses on the box. `tests/` is a namespace package, and nothing in the freeze installs a top-level `tests`. `sample_vllm` is launched by script path under `/opt/vllm`, so it needs no `.pth`.
- **Orphan label collision.** Runs are matched on instance id (already tested by O).

## Simplicity audit (`policy.py`, 438 lines)

It cannot go below ~200 lines without dropping commands: there are 7 commands with real guards. About 70 lines (≈ 15 %) can go, and the first item also fixes BK4-12:

1. **Filter as data** (`policy.py:70-102`, 33 lines). Build one `FILTER = [(query_term, predicate), …]` list and derive both `query` and `reasons` from it. This also removes the datacenter test duplicated in `governance` (`:113`). About -12 lines, and the query and the re-check can no longer drift apart.
2. **Seven shims → one.** The 7 × 5-line `*.sh` files differ only in the subcommand. Use `vast.sh` with `exec … policy "$(basename "$0" .sh)" "$@"` plus symlinks, or call `python -I -m experiments.phase2.vast.policy <cmd>` directly from the plists and docs. About -30 lines across files.
3. **argparse help strings** (`:418-432`) repeat the module docstring's usage lines. Drop them (the docstring is the `--help` text). About -7 lines.
4. **`pull` + `cmd_pull`** (`:371-381`) and **`get_run`** (`:297-298`) can merge into `cmd_pull` and an inline `active()[id]`. About -5 lines.
5. **Four `data/raw` guards.** `under_raw` (`:301`), `code_files`' `data/` filter (`:309`), the rsync `--exclude=data/raw/` (`:201`) and `bundle.is_raw` (`bundle.py:27`) all guard the same thing. Move `is_raw` into `policy` and use it for both callers. About -3 lines.
6. **The `SECRET_RE` re-check of the offer JSON** (`:250-251`) is redundant: the same row then goes through `ledger()`, and Vast's offer JSON carries no secrets. -2 lines.
7. The fixes for BK4-3/4/7/8/10 add about 15 lines. Net is about 380 lines.
