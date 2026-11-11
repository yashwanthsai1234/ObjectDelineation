# Experiments

Scripts that produced the published numbers. Run each as a module:

```bash
.venv/bin/python -I -m experiments.<folder>.<script>
```

`-I` is required for scripts that read `data/` (downloaded data is untrusted).

| Folder | What it holds | Results |
|---|---|---|
| `phase1/` | Data and target builders, reward gates, evaluation runner | [RESULTS.md](phase1/RESULTS.md) |
| `phase2/` | SFT and RL training, final evaluation, Vast rental tooling, run configs in `configs/`, Hub push and pull | [RESULTS.md](phase2/RESULTS.md) |
| `cuda_checker/` | The checker on a GPU box, compared with the CPU build | [RESULTS.md](cuda_checker/RESULTS.md) |

All three folders run today. What is frozen is the record, not the script:

- `RESULTS.md` files: frozen sections are never edited. A correction is a dated note below them.
- Plans and judge reports (`PLAN*.md`, `JUDGE_*.md`, `V0_V1.md`, `EXECUTION.md`, `RESEARCH.md`), every `impl/` folder and `cuda_checker/results/`: written once, kept as they were. The same holds for `docs/research/`.
- Code that no script runs any more moves to `archive/<name>/` (local only).

New scripts live in `experiments/<phase>/`. Plans for runs are typed run specs in `specs/`. The old chains and SFT configs of phase2/configs/ are archived, local only.

Headline numbers are in [reference/headline.json](../reference/headline.json). `make headline` rebuilds that file from a `final_eval` result.
Run specs are in [specs/](../specs/).
