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

Headline numbers are in [docs/results/headline.json](../docs/results/headline.json). `make headline` rebuilds that file from a `final_eval` result.
Run configs are listed in [docs/reference/configs.md](../docs/reference/configs.md).
