# Custom kernels (drafts)

K1 (polygonizer), K2 (scoring) and K3 (distance/topology helpers) are design drafts with reference code.

- They are not built into the `objdel` package and are not imported by any library code.
- They are validated only as each `GPU_RUN.md` states: CPU builds and tests ran here; GPU runs are planned inside the Vast budget of `../Implementation.md` section 6.
- Treat any speed number as unmeasured until a GPU_RUN result file says otherwise.
- Each folder holds a `DESIGN.md` (read first) and a `GPU_RUN.md` (exact run steps).
