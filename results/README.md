# Results

The frozen records of what produced the published numbers.

| File | What it records |
|---|---|
| [phase1.md](phase1.md) | Phase 1: the exact reward, the data and the training integration |
| [sft_r1.md](sft_r1.md) | Phase 2: SFT on Vast AI |
| [checker.md](checker.md), `checker/` | The checker on a GPU box, compared with the CPU build (`tools.checker.eval` writes both) |

- Frozen sections are never edited. A correction is a dated note below them.
- `checker/`: written once, kept as it was.
- Code that no script runs any more moves to `archive/<name>/` (local only).

The commands are `objdel <group> <command>` (src/objdel/cli), the tools are under [tools/](../tools/), the run specs
under [specs/](../specs/). Headline numbers are in [reference/headline.json](../reference/headline.json);
`make headline` rebuilds that file from a final eval result.
