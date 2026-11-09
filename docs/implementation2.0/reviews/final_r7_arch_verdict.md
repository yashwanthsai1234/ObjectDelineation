# Verdict on final review round 7, architecture (final_r7_arch.md)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item makes an NVIDIA product claim; no catalog lookup or install needed. |
| engineering:code-review | Opened every cited code line: `select_ckpt.py:359-379`, `g_probe.py:131-145`, `policy.py:425-490` (`chain`, `job`), `tests/test_p2_scripts.py:340-390`, `preflight.py:34`, `pyproject.toml`. |
| engineering:debug | Re-ran the reviewer probe `(planning probe, not kept: r7arch)/probe_select.py` in the r6 venv (pydantic-settings 2.16.2): `--procs 1` and `--no-vision` exit 2 on E6.1's `SelectStage`; `--timeout_s` exits 2, `--timeout-s` passes. Traced `job()` for an argv that starts with `pytest`. |
| engineering:system-design | Walked the V4 chain: a non-soft failing step 2 stops the `&&` chain, so step 3 (the only V4 eval) never runs. |
| engineering:architecture | Checked E2.1's "every E0 flag becomes a field" rule against E6.1's explicit field lists (arch-2) and A3's `--rescore`/`--hash-view` forms against `EvalStage` (arch-4). |
| engineering:documentation | Ran the check 14 regex on `docs/implementation2.0/kernels` (16 lines, 7 files) and read DOC1's kernel bullet and re-scan sentence. |
| huggingface-skills:huggingface-best | No model choice in this round; nothing changed. |
| WebSearch/WebFetch | Loaded; not needed: the parser facts were re-measured locally. |

## Verdicts

- arch-1: downgraded (minor). E2's code does import `CliSuppress` and the E2 Libraries line names only pydantic, but the very first E2 test fails with an immediate `ModuleNotFoundError` that names the fix, U42 already lists `pydantic-settings>=2.16,<3` as a base dependency mapped to E2, and the live box image (vllm_native, ADR 0018) has pydantic-settings 2.16.2 in its lock, so no outcome or GPU spend changes.
- arch-2: confirmed (major). E6.1's `SelectStage`/`ProbeStage` code has no `procs` (nor `k`, `seed`, `model`, `timeout_s`, `out`, `no_vision`) while `select_ckpt.py:369-376` and `g_probe.py:139` define them and the preflight dry run passes `--procs 1` to both [tests/test_p2_scripts.py:369,384; preflight.py:34]; re-run probe exits 2 on `--procs 1` and on `--no-vision`, so the E2.1 rule and the E6.1 code disagree and the `vision=None -> --no-vision` path cannot round-trip.
- arch-3: confirmed (major). E.8 V4 step 2 and E.9 write the box tests as `pytest -m gpu ...` and E6 step 8b says "CmdStages with the same argv"; `policy.py:486` turns any argv not starting with `bash` into `/opt/venv/bin/python -I pytest ...`, which fails, stops the `&&` chain before step 3, and U51 does not check the leading tokens (B.md I22a pins `-m, pytest` for V1, E does not).
- arch-4: confirmed (major). E6.1 puts `--hash-view` "out" of `EvalStage` (required `ckpt`, `extra="forbid"`) and E6.4 routes old flags through `cli_parse`, yet A3 needs `final_eval --hash-view VIEW --out FILE` (A.md A3, U54, B1.0/I36) and `--rescore DIR --view V --split S --out O` with no `--ckpt` (A.md:1705); both exit 2 under the stated parser and the doc gives no home for either mode.
- arch-5: confirmed (major). Check 14 scans `docs/implementation2.0/kernels/**`, but DOC1's kernel bullet replaces only absolute paths and the re-scan sentence omits the kernel folders; measured `grep -rnE` with the check 14 patterns finds `S[P]/impl` and `S[P]/[A-Z]` hits in `K1/GPU_RUN.md:20,25,31`, `K1/DESIGN.md:3,242,304`, `K2/DESIGN.md:3`, `K3/DESIGN.md:3,9`, so check 14 stays red at DOC4.
- arch-6: minor (as filed). Stale counts in DOC1; wording only.
- arch-7: minor (as filed). Teardown after a failed pull; the watchdog cap bounds cost.
- arch-8: minor (as filed). Round 6 minors carried over.
- arch-9: minor (as filed). E-G13/GP3 wording; follows arch-1.
