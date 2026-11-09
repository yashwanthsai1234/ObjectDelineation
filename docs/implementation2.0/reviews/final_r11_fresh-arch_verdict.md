# Verdict on final review round 11, fresh-arch

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | None of the items depends on an NVIDIA product, so I did no catalog lookup. |
| engineering:code-review | I read `train_rl.py:1-70` (`final/` is written only after `trainer.train` returns), `vast/policy.py:466-494` (`cmd_run` joins steps with `&&`; `~` makes a step soft), `policy.py:59-61` (`UP_GB_PER_H = 6` for checkpoint pulls), `com.objdel.vast.pull.plist` (StartInterval 1200), `gpucheck_py.cpp:140-158` (`check_masks` returns rings) and `targets/build.py:80-97`. |
| engineering:debug | I traced the crash path. A non-zero `train_rl` exit skips step 4 (`&&`). The later `pull.sh` still brings `trainer/checkpoint-N` back, and `pull.sh` excludes only `rolling/decay/stable`. So a crash loses the B1 eval, and a host loss loses everything. |
| engineering:system-design | I checked the pairing of Gate B3. Before E, its arms differ only by host and by D7, and D7's byte checks cover D7, so that part is minor. After E, the sampler differs, and no $0 check covers it. |
| engineering:architecture | I weighed the polygonizer seam. D3 already plans `rings.py` over `_gpucheck.check_masks`. Today's `check_masks` is measured identical to the numpy path, so the seam can land before B1.4. |
| engineering:documentation | Each confirmed item below quotes the source heading and the code line or measurement. |
| huggingface-skills:huggingface-best | None of the items is about a model choice. No change. |

## Verdicts

- fresh-arch-1: confirmed (major). E4 step 3 says "`final_eval` calls `VllmSampler` + `evaluate`. It no longer calls `build_trainer`", and the E4 Gate makes the TRL `SamplingParams` comparison box-only. But 00_front "End state under (b')" lets V3 run after E0/E6, and V3 has no B1 re-eval. So on that path Gate B3 pairs a new-sampler B3 eval with V1's trainer-path B1 eval. E.8 says only that V4 checks "can ride" on a later B run, and front 6.3 does not price that. The pre-E part (host and D7) is minor on its own: D7 has byte-identity checks, and the Gate B1 fallback already accepts cross-host pairing.
- fresh-arch-2: confirmed (major). The V1 launch lines in 00_front 6.2 have one `pull.sh`, at the end. The pull plist exists (`StartInterval 1200`) but no plan step loads it (`launchctl list | grep objdel` is empty, measured 2026-11-09). r_s already pays for the pulls (`UP_GB_PER_H = 6`, "checkpoint pulls of 3-5 GB per 20 min", `policy.py:59-61`). The B1.5 "Resume" task (`--resume` from `trainer/checkpoint-N`) has nothing to read after a host loss. A `train_rl` crash skips the B1 eval through the `&&` join (`policy.py:490`). `BudgetStop` covers only the budget trip. FR12 covers only launch and image-pull failures. No earlier review raised this.
- fresh-arch-3: confirmed (major). B.3.3 says "Owner exception needed: the targets use the numpy polygonizer ... ends when Phase D K1 passes parity", and 99_back decision 9 sets the default to [accept]. That default breaks the binding owner rule "CUDA for the polygonizer". Yet `check_masks` already matches the numpy polygonizer on 37,093 masks with 0 mismatches (K1/results/cross_numpy_m3.txt: "check_masks IDENTICAL", "0 mismatches -> PASS"). It also matches on an RTX 5090 (issues_solved_2026-11-07.md:37). It takes 4.83 s for all masks (bench_cpu_m3.txt). D3 already wraps `check_masks` in `rings.py`. The round-1 handling (B log row 24, "Exception in ADR-0015 with end condition") rests on a premise that these measurements disprove. The end condition "K1 passes parity" is also undefined (Mac CPU parity, or a GPU run that (b') never funds). This is major, not blocking: the doc asks the owner for sign-off instead of silently breaking the rule.
- fresh-arch-4: minor (agreed). The tie of C1 adoption to the speed gate is stated and coherent. The suggested change is an owner option.
- fresh-arch-5: minor (agreed). The text is stale: 99_back row 12 still says "V2 item 4 (R-D1", the diagram says "V3: R-D2 first", and front 6.3 drops R-D1 and makes R-D2 optional. Front 6.3 is the binding table.
- fresh-arch-6: minor (agreed). D3 names phase0 and lab scripts that CLEAN5 wave 1 archives (CLEAN.md:495, 500). The glob then matches nothing and nothing fails.
- fresh-arch-7: minor (agreed). E.3.3 "after C has removed syncs" is false under (b'), but the gate records whatever count E0 has.

Counts: confirmed 3 (major 3), refuted 0, downgraded 0, minor as filed 4.
