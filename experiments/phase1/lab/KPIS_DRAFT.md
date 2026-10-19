# Phase 1 KPIs (draft v0, to be refined after research)

Each module is accepted only when its KPIs are met by a reproducible test in the lab.

| ID | Module | KPI | Target | Measured by |
|---|---|---|---|---|
| K1 | rasterizer (polygon → mask) | agreement with exact oracle, both fill rules, incl. self-intersecting / collinear / on-edge centres | 100 % of pixels on ≥ 10k random shapes | oracle.py |
| K2 | legality (graded) | label agreement with exact oracle | 100 % on ≥ 10k rings incl. grid-degenerate | oracle.ring_problems |
| K3 | reward throughput | wall time to score one GRPO batch (256 rollouts, realistic sizes) | ≤ 250 ms on 8 cores (≤ 5 % of a step) | benchmark |
| K4 | reward latency | p99 per rollout | ≤ 5 ms | benchmark |
| K5 | SFT targets | topology preserved after simplify → re-rasterize | 100 % of expert masks | checker |
| K6 | SFT targets | IoU of re-rasterized target vs source mask | min ≥ 0.98, median ≥ 0.99 | checker |
| K7 | SFT targets | tokens per object (chosen tokenizer) | median ≤ TBD after tokenizer test | tokenizer |
| K8 | scale-aware band | expert masks scored leave-one-rater-out: independence from object size | Spearman abs(rho) with area ≤ 0.1; median ≥ 0.9 | fixtures |
| K9 | reward discrimination | high-IoU wrong-topology perturbations ranked below topology-correct ones | ≥ 95 % of matched pairs | fixtures + perturbations |
| K10 | reward determinism | identical scores and hashes across runs, processes, thread counts | 100 % | test |
| K11 | data splits | patient / image leakage across splits; near-duplicates across splits | 0 / 0 | split tests |
| K12 | training stack | chosen framework runs 2 GRPO steps with our reward and the smallest model on CPU/MPS, logging each reward term | pass | smoke test |
| K13 | engineering | line coverage of src/objdel; CI; lint | ≥ 90 %; green; clean | CI |
