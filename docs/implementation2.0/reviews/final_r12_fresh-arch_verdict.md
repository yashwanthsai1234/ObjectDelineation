# Verdict: final_r12_fresh-arch

## Skills used
- nvidia-skills:nvidia-skill-finder: loaded; no item depends on an NVIDIA product, so no catalog lookup.
- engineering:code-review: read policy.py code_files() (drop regex line 379), cmd_run/cmd_pull/cmd_push_hf, watchdog.sh, pull plist.
- engineering:debug: traced item 1 (B1.6 lookup -> code sync drop -> E05 refusal) and item 2 (pull job -> PULLED_OK -> idle destroy).
- engineering:system-design: checked item 3 against the V1 time table and the session-clock BudgetStop design.
- engineering:architecture: items 5-6 are design trade-offs; they stay minor.
- engineering:documentation: every confirmed item quotes the source heading and the code line.
- huggingface-skills:huggingface-best: no model choice; read huggingface_hub 2.2.0 snapshot_download / hf_hub_download (Xet files skip HEAD when the tree listing is cached, file_download.py:1363-1366).

## Verdicts
- fresh-arch-1: CONFIRMED (blocking). B.md B1.6 "final_eval ... load frozen_views.json" at `docs/reference/frozen_views.json` and E05 "missing frozen_views.json gives a clear refusal", but code_files() drops `^(data|docs)/` [code: experiments/phase2/vast/policy.py:379] and the only planned keep rule (00_front section 5, D1 step 0) covers docs/implementation2.0/kernels/K1-K3 and lands after V1; every $0 check (6.1 step 4c, E05, I05, CT-30 imports only) runs where docs/ exists, so every V1 eval on a frozen view refuses on the box after spend.
- fresh-arch-2: DOWNGRADED (minor). 00_front 6.2 "Why the pull job" already states the exact constraint ("the watchdog destroys an idle box 30 min after a pulled EXIT ... a3_measure.yaml must start within minutes of the b1_rows EXIT") and the launch line says "start right after the b1_rows EXIT"; the deadline is given, only the polling tool is unstated and "as today" is a wording slip.
- fresh-arch-3: DOWNGRADED (minor). The download time is unmeasured, but the plan routes any overrun through the session-clock BudgetStop, which B1.5 justifies with "the rows-root download also comes before train_rl starts" and front 6.2 "If it trips, the adapter is saved and items 3-5 still run", so no spend is wasted and Gate B1 still gets its inputs (fewer steps only); the >15 min HEAD-per-file path assumes non-Xet files, while hf_hub_download skips HEAD for Xet files with a cached tree listing [code: .venv huggingface_hub/file_download.py:1363-1366].
- fresh-arch-4: minor as filed (no change).
- fresh-arch-5: minor as filed (no change).
- fresh-arch-6: minor as filed (no change).

Counts: confirmed 1, refuted 0, downgraded 2.
