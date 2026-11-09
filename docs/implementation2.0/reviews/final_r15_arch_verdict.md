# Verdict on final_r15_arch.md (fresh verifier, 2026-11-09)

## Skills used
- nvidia-skills:nvidia-skill-finder: no NVIDIA product claim in these items; no catalog lookup, no install.
- engineering:code-review: re-read `experiments/phase2/select_ckpt.py:88-96` and `final_eval.py:70-76` to check the caller claims in item 1.
- engineering:debug: walked the V4-failure branch (restore files -> who calls them -> U2/E-G14 -> which eval path is live).
- engineering:system-design / engineering:architecture: a rollback must restore a working path; checked E4 step 1 against the E.7 rollback row and E-G14.
- engineering:documentation: confirmed the stale E.7 mitigation "samplers archived only after V4" (E.md line 1837).
- huggingface-skills:huggingface-best: no model choice in scope; nothing changed.

## Verdicts
- arch-1: CONFIRMED (major). E.md "### E4" step 1 says "If V4 fails, the owner restores them from the tag", but E-G14 (E.md line 1817) and U2 (line 2031) say the two samplers exist only under `archive/phase2_samplers/` "no exception" / "with no V4 exception", step 1 itself switches their only caller (`select_ckpt.py:94-95`, `head = [vllm_python(), "-I", str(ROOT / "experiments/phase2/sample_vllm.py")]`) to `objdel.policy.sampler`, so the restore helps no path and leaves `final_eval` on the failed `VllmSampler`; the E.7 row still says "samplers archived only after V4" and "revert E4 only", and 00_front 6.3 repeats the restore, so the engineer must guess the real V4-failure action (the fresh-arch-1 handling in the round-14 log is itself the wrong text). Not blocking: only on a conditional, owner-approved V4 rental path.
- arch-2: minor as filed (not counted); a missed test parameter and commit split that R-3 `make test` catches at $0.
- arch-3: minor as filed (not counted); dead `vllm_python`/`VLLM_VENV` and HF carry-over detail, surfaced by CLEAN16 at $0.
- arch-4: minor as filed (not counted); `RlStage` lacks `merge_last` / `opsm_start_delta`, caught by the E6.0 spike, golden 5 and E-G9 at $0.
- arch-5: minor as filed (not counted); `LIVE_CLIS` list incomplete, caught by the E8 gate and `test_cli_smoke` at $0.

Counts: confirmed 1 (major), refuted 0, downgraded 0, minor as filed 4.
