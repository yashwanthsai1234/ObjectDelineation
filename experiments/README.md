# Experiments

Scripts that produced the published numbers. Run each as a module:

```bash
.venv/bin/python -I -m experiments.<folder>.<script>
```

`-I` is required for scripts that read `data/` (downloaded data is untrusted).

| Folder | What it holds | Results |
|---|---|---|
| `phase1/` | Data and target builders, evaluation runner, cross-arch hash tool | [RESULTS.md](phase1/RESULTS.md) |
| `phase2/` | SFT and RL training, final evaluation, Vast rental tooling, Hub push and pull, checker and data tools | [RESULTS.md](phase2/RESULTS.md) |
| `cuda_checker/` | The checker on a GPU box, compared with the CPU build | [RESULTS.md](cuda_checker/RESULTS.md) |

All three folders run today. What is frozen is the record, not the script:

- `RESULTS.md` files: frozen sections are never edited. A correction is a dated note below them.
- Plans and judge reports (`PLAN*.md`, `JUDGE_*.md`, `V0_V1.md`, `EXECUTION.md`, `RESEARCH.md`), every `impl/` folder and `cuda_checker/results/`: written once, kept as they were.
- Code that no script runs any more moves to `archive/<name>/` (local only). Finished one-time studies (the v4 reward gate tables, the B1 gate power simulation, the B1 view rebuild, the phase E file-move tool) and the run specs that `specs/b1_rows.yaml` replaced are there.

New scripts live in `experiments/<phase>/`. Plans for runs are typed run specs in `specs/`: `b1_rows.yaml` (the RL run), `r1_qwen3vl2b.yaml` (SFT) and the C speed and B2 fp16 arm specs.

Headline numbers are in [reference/headline.json](../reference/headline.json). `make headline` rebuilds that file from a `final_eval` result.
Run specs are in [specs/](../specs/).
