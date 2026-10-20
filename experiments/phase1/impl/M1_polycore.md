# M1: polycore, build system and CI

Gate: K1, K2, R9 (partial) (plan §7). Promotes `lab/t_polycore` (core, bindings, wrapper, tests, fuzz driver) and
`lab/common/oracle.py` (as test oracles).

**Success criteria, written first:**
- 0 wrong pixels against the exact oracle on ≥ 10k shapes, both fill rules (K1);
- 0 label mismatches on ≥ 10k rings (K2);
- input caps rejected with `ValueError`;
- sanitizer fuzz clean on macOS and Linux;
- wheels built and tested for manylinux x86-64 + aarch64 and macOS arm64;
- `make test` and `make lint` green after the build-backend switch;
- speed at least the lab's (batch of 256 in 8.5 ms; p99 0.38 ms).

## What was built

| Path | Contents |
|---|---|
| `cpp/polycore/polycore.hpp` | header-only C++17 core: `check_rings/check_roles/check_dims/check_grid`, `signed_area2`, `rasterize`, `validate`, `score_batch` (696 lines) |
| `cpp/polycore/tests/fuzz_main.cpp` | native fuzz/differential driver: spec in `__int128` per pixel, all-pairs edge scan, cap/malformed rejection, thread-count determinism |
| `src/objdel/_polycore.cpp` | nanobind module `objdel._polycore`; inputs (caps included) checked before the GIL is released |
| `src/objdel/checker/polycore.py` | Python API (contract in EXECUTION.md §3) |
| `CMakeLists.txt` | `nanobind_add_module(_polycore NB_STATIC NOMINSIZE STABLE_ABI …)`; in editable installs, also a `.pth` with the repo root |
| `pyproject.toml` | `[build-system]` scikit-build-core 1.1.1 + nanobind 3.1.0 (pinned); `[tool.scikit-build]`; `[tool.cibuildwheel]` |
| `Makefile` | new targets `build`, `test-all`, `cov`, `fuzz`, `wheels`; `test` now excludes `slow`; `setup` installs `.[models,train,data,dev]` (now `.[train,data,dev]`: the `models` extra and `src/objdel/models` were deleted, S4) |
| `.github/workflows/tests.yml` | jobs `python` (lint + pytest + coverage; ubuntu x86-64, ubuntu-24.04-arm, macOS arm64), `gates` (slow K1/K2), `sanitizers` (`make fuzz`, 3 OS), `wheels` (cibuildwheel, 3 OS); later + `train` (training stack, coverage K13; F5) |
| `tests/test_polycore.py`, `tests/test_polycore_oracle.py`, `tests/test_polycore_props.py` | unit, oracle (K1/K2), Hypothesis property tests |
| `tests/oracles/raster_oracle.py`, `tests/oracles/validity_oracle.py` | slow exact oracles and shape generators, derived from `lab/common/oracle.py` and the lab helpers |

## API (contract changes, EXECUTION.md §3 updated)

| Function | Notes |
|---|---|
| `DEFECTS`, `MAX_RINGS = 10`, `MAX_RING_VERTICES = 1000`, `MAX_DIM = 4096` | |
| `pack(rings) -> (xy int32 (N,2), off int32 (R+1,))` | **off is int32, not int64.** The caps keep it far below 2³¹, and one dtype means no conversion on the hot path. |
| `rasterize(shape, W, H, G=1000)` | `evenodd`, `nonzero`, `ring_wmin`/`ring_wmax`, `area_*`, `px_multi`, `rings_mixed_sign`. The lab's optional `winding` array is dropped: masks plus per-ring ranges are the spec, and the fuzz driver checks winding natively. |
| `validate(shape, roles=None, G=1000)` | `total`, **`per_ring`** (the lab called it `ring`), `area2`, `n_clean`, `pair_cross`, `pair_touch` |
| `signed_area2(ring) -> int` | exact (`__int128`, returned as hi/lo and joined in Python); not capped |
| `score_batch(shapes, Ws, Hs, roles=None, gts=None, G=1000, threads=0)` | takes shapes like `rasterize`, not the lab's parallel `xys, offs` lists |

Removed from the lab API:
- `score`: a single rollout is `score_batch` of 1;
- `ring_labels`: a test helper, now `validity_oracle.labels`;
- `hardware_threads`, `MAX_RINGS_VALIDATE`;
- `validate(method="brute")`: the brute pair scan lives on only as test code in the fuzz driver.

## Decided semantics (plan §4.1)

- **`wrong_winding` uses net area.** The net shoelace area of the cleaned ring decides, so a bow-tie is judged by its larger lobe, and a zero-net-area bow-tie is `zero_area` and not `wrong_winding`. Test: `test_wrong_winding_uses_net_area`.
- **A hole crossing its outer gets both flags**, `ring_ring_crossing` and `hole_outside_outer`. Test: `test_hole_crossing_its_outer_gets_both_flags`.
- **Caps.** At most 10 rings, at most 1,000 raw vertices per ring (before de-duplication), 1 ≤ W, H ≤ 4096, 1 ≤ G ≤ 2³⁰. Anything else raises `std::invalid_argument`, which nanobind maps to `ValueError`. Coordinates may be any int32: out-of-grid vertices are a defect (`out_of_range`), not a rejection. Tests: `test_caps_are_accepted_at_the_limit` and `test_beyond_caps_raise_value_error` (7 cases); the fuzz driver checks `check_rings` against an independent predicate.

## Research

**Build**
- scikit-build-core editable installs ([editable](https://scikit-build-core.readthedocs.io/en/latest/configuration/editable.html), [config](https://scikit-build-core.readthedocs.io/en/latest/configuration/index.html)):
  - the default redirect mode serves Python files from `src/` and the extension from site-packages;
  - `editable.rebuild` is experimental. It needs a persistent `build-dir` and the build tools in the runtime env (`--no-build-isolation`), and a rebuild only takes effect in new processes.
- nanobind ([CMake API](https://nanobind.readthedocs.io/en/latest/api_cmake.html), [GIL](https://nanobind.readthedocs.io/en/latest/functions.html#gil-management), [nanobind_example](https://github.com/wjakob/nanobind_example)):
  - `-Os` is the default and `NOMINSIZE` disables it;
  - `NB_STATIC` is the default;
  - `STABLE_ABI` needs Python ≥ 3.12, and with `wheel.py-api = "cp312"` it gives one abi3 wheel per platform;
  - `std::invalid_argument` is mapped to `ValueError` (checked in nanobind's `nb_internals.cpp`).
- cibuildwheel 4.3.1 ([options](https://cibuildwheel.pypa.io/en/stable/options/), [changelog](https://cibuildwheel.pypa.io/en/stable/changelog/)): the default image is manylinux_2_28 (since 4.0), the default frontend is `build`, abi3 wheels get `abi3audit`, and the `test-command` runs outside the source tree (hence `test-sources`).
- GitHub arm64 runners: `ubuntu-24.04-arm` is GA for public repos ([changelog 2025-08-07](https://github.blog/changelog/2025-08-07-arm64-hosted-runners-for-public-repositories-are-now-generally-available/)) and available to private repos since Jan 2026 ([report](https://dev.classmethod.jp/articles/github-actions-arm64-private-repos/)).
- Action versions: checkout v6, [setup-uv v10](https://github.com/astral-sh/setup-uv/releases), [upload-artifact v7](https://github.com/actions/upload-artifact/releases), cibuildwheel v4.3.1.

**Algorithms**
- **Rasterization.** The lab accumulates per-pixel crossing differences into H×(W+1) and H×W int32 buffers, the cell/cover approach of AGG and font rasterizers ([notes](https://rasterizer.readthedocs.io/en/stable/algorithm.html)). The classic alternative buckets edge crossings per scanline and walks them in x order as spans ([Wikipedia](https://en.wikipedia.org/wiki/Scanline_rendering), [Stanford CS248](https://graphics.stanford.edu/courses/cs248-01/scan/scan1.html)). Both apply the same exact integer crossing formula, so exactness is unchanged. The tie rule is the GPU top-left convention ([Direct3D rules](https://learn.microsoft.com/en-us/windows/win32/direct3d11/d3d10-graphics-programming-guide-rasterizer-stage-rules)).
- **Validity.**
  - Bentley–Ottmann ([Wikipedia](https://en.wikipedia.org/wiki/Bentley%E2%80%93Ottmann_algorithm)) costs O((n+k) log n). Every intersecting pair must be counted anyway (counts are reported per class), so it cannot beat the lab's sort-and-sweep ([sweep and prune](https://en.wikipedia.org/wiki/Sweep_and_prune)), which costs O(n log n + candidate pairs) and is already sub-quadratic on smooth rings (20 µs at 1,000 vertices against 232 µs brute). It would add robustness-critical event-queue code for no gain, so it was not adopted.
  - Shewchuk's adaptive predicates ([robust.html](https://www.cs.cmu.edu/~quake/robust.html)) exist to make floating-point predicates exact. Our inputs are integers, and int64/`__int128` orientation tests are exact by construction, so nothing needs adopting.

## Sandbox alternatives and numbers

All in the scratch dir. Native differential bench: lab header vs the new header, each in its own translation unit with the namespace renamed.

| Variant | Exactness | Rollouts (256, lab workload): p50 / p99 per rollout | Memory |
|---|---|---|---|
| lab accumulation buffers | (reference) | stats only 26.0 / 109.8 µs; with masks 39.5 / 169.2 µs | H×(W+1) + H×W int32 scratch per thread, kept at its maximum (134 MB at 4096²) |
| **span walk over sorted crossing keys (adopted)** | **0 mismatches vs lab on 10,256 cases × 2 modes (masks, stats, per-ring ranges)** | **17.1 / 50.1 µs; with masks 19.9 / 53.8 µs** | O(crossings), rows in bands of 512 (≤ 41 MB worst case) |

- **Pathological worst case at the caps** (10 rings × 1,000 vertices zigzagging the full height of a 4096² crop, 41 M crossings): lab 0.62 s and 162 MB RSS, new 1.03 s and 137 MB RSS. The sort costs more here. This is accepted: it is bounded by the caps and 8× smaller at realistic crop sizes.
- **Python, same process, same load** (lab `_polycore` vs new, `score_batch` of 256, median):

  | Threads | gt | lab | new |
  |---|---|---|---|
  | 1 | yes | 15.6 ms | 12.2 ms |
  | 1 | no | 11.7 ms | 5.8 ms |
  | 8 | yes | 5.4 ms | 4.0 ms |
  | 8 | no | 4.1 ms | 2.0 ms |

  Outputs are identical.
- **Not tried:** per-row incremental (DDA) crossing instead of one division per crossing. Divisions are a small share once the O(W·H) passes are gone, and the gates have a ~50× margin.
- **Editable rebuild on import** (`editable.rebuild=true`, `--no-build-isolation`): it works and rebuilds after a header edit (2.5 s). But it adds about 0.2 s to every process that imports `objdel._polycore`, races between parallel processes and agents, needs the build tools in the runtime venv, and is experimental. Not adopted. Instead, `make build`, and `make test` itself, reinstall incrementally when `cpp/polycore/polycore.hpp`, `src/objdel/_polycore.cpp` or `CMakeLists.txt` are newer than a stamp (`build/.polycore-stamp`). Isolated reinstall takes 4–6 s.
- **Stable ABI:** no measurable cost (a Python call to `rasterize` on 8×8 takes 3.3 µs), and the wheel count drops to one per platform. Adopted.

## Build-system switch (hatchling → scikit-build-core)

- **Editable install.** The default redirect mode. Verified in a sandbox mirror of the repo, then in the repo:
  - new modules and subpackages under `src/objdel` import without a reinstall;
  - `importlib.resources` sees both trees.
- **`experiments` importable.** Hatch's `dev-mode-dirs = ["src", "."]` is replaced by `_objdel_repo.pth`, which CMake writes and installs only when `SKBUILD_STATE == "editable"`. `python -I -m experiments.<phase>.<script>` works (`-I` still processes site-packages `.pth` files). A new `experiments/phase9/s.py` also ran without a reinstall.
- **Persistent build dir.** `build-dir = "build/{wheel_tag}"`, which is git-ignored and sits next to `make cuda`'s `build/` CMake tree without touching it.
- **The switch.** Done in one step: write `pyproject.toml`, then `uv pip install --no-deps -e .`. The wheel was built before the old one was removed, so `import objdel` was unavailable only for the swap, which took milliseconds. The full `make test` was green straight after (237 passed, 7 skipped, 6 deselected, 68 s).
- **sdist.** It excludes `data/`, `.venv` and `build/` (git-ignored); it is 3 MB. Note: the root `.gitignore` pattern `Makefile` (meant for CMake output) also drops the root Makefile from the sdist. This is harmless. *(The pattern is no longer in `.gitignore`.)*

## Tests

| File | Fast (default) | Slow (`-m slow`) |
|---|---|---|
| `test_polycore.py` | 31 tests:<br>- contract, keys and dtypes;<br>- top-left rule;<br>- shared-edge partition (no gaps or overlaps on cut lines);<br>- mixed-sign and multi-winding;<br>- 9 targeted validity cases;<br>- both decided semantics;<br>- exact `signed_area2` beyond int64;<br>- roles checks;<br>- caps at the limit and beyond (7 cases);<br>- malformed raw input;<br>- 3,000 garbage inputs (only `ValueError`/`OverflowError`/`TypeError`);<br>- batch = single calls (gt as bool, uint8 0/255 and int64);<br>- determinism across thread counts 0–64, Python threads and processes;<br>- GIL release (a spinner thread keeps ≥ 30 % of its solo rate during a kernel) | — |
| `test_polycore_oracle.py` | 11 tests:<br>- vectorised oracle = pure oracle (150 cases);<br>- K1 450 shapes;<br>- on-edge lattices (5 crop sizes);<br>- int32-extreme coordinates on the `__int128` path vs the pure oracle;<br>- K2 1,400 rings;<br>- fine classes on 300 rings;<br>- multi-ring on 300 shapes | 5 tests:<br>- K1 10,000 shapes;<br>- K1 400 shapes vs the pure-Python oracle;<br>- K2 10,000 rings;<br>- fine classes on 3,000 rings;<br>- multi-ring on 1,500 shapes |
| `test_polycore_props.py` | 8 Hypothesis properties:<br>- equals the pure oracle on lattices (incl. per-ring ranges);<br>- translation by whole pixels;<br>- mirror differs only at on-edge centres;<br>- reversal (same masks, negated ring ranges);<br>- start-vertex invariance;<br>- validity affine invariance;<br>- labels on lattices;<br>- `score_batch` determinism across threads | — |

Fast: 50 tests in about 7–11 s. Slow: 5 tests in 30 s (3 worker processes). Line coverage of `src/objdel/checker/polycore.py`: 98 %.

## Gate results

| Gate | Result |
|---|---|
| K1 | **0 / 849,802,293 pixels wrong**, both rules, on 10,000 shapes. Families: star 1,112; soup, grid, halfpixel, holes, multipart, tiny, degenerate and rect 1,111 each. The cases exercise 6,510 on-edge centres (10 % sample), 2,426 cases where even-odd ≠ nonzero, and 2,839 mixed-sign rings. Per-ring ranges and statistics: 0 wrong. Pure-Python oracle: 0 wrong on 400 shapes (34.5 M px). |
| K2 | **0 / 10,000 rings** mismatched. Labels: 4,480 self-crossing, 3,045 self-touching, 325 zero-area, 741 fewer-than-3, 3,672 valid. Fine classes vs the all-pairs reference: 0 / 3,000. Multi-ring (all 13 classes, pair matrices): 0 / 1,500. |
| Fuzz, macOS | `make fuzz FUZZ_N=200000` (Apple clang, ASan+UBSan, `-Werror`): clean. 156,995 shapes checked against the spec and the all-pairs scan, 43,005 malformed or over-cap inputs rejected; 343 s. |
| Fuzz, Linux | Docker `gcc:14` (g++ 14.4), ASan+UBSan, `-Werror`: arm64 50,000 iterations clean (39,248 checked, 10,752 rejected, 67 s); amd64 under emulation 10,000 iterations clean (7,814 checked). |
| Wheels | cibuildwheel 4.3.1 from the sdist in Docker: `objdel-0.1.0-cp312-abi3-manylinux_2_26_aarch64.manylinux_2_28_aarch64.whl` and `…manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl`. Both pass `abi3audit --strict` and the in-container test command (50 passed). macOS arm64: `uv build --wheel` with `MACOSX_DEPLOYMENT_TARGET=11.0` gives `objdel-0.1.0-cp312-abi3-macosx_11_0_arm64.whl` (minos 11.0), installed in a fresh venv: 50 passed. cibuildwheel itself was not run for macOS locally; it runs in the CI `wheels` job. |
| CI YAML | `rhysd/actionlint` (Docker): 0 errors. The workflow has not run on GitHub (no push from here). |
| Speed (K3/K4 workload, repo build, load average ≈ 8) | `score_batch` of 256 rollouts with gt: 19.4 / 11.0 / 7.2 / 4.9 ms median at 1 / 2 / 4 / 8 threads. Per rollout (`score_batch` of 1, incl. Python): p50 0.061 ms, p99 0.195 ms. `rasterize` with masks: p50 0.021 ms, p99 0.051 ms. |

## Open issues

1. **Counts are edge pairs.** One pinch counts 4, so the reward should use presence or a saturating transform (lab risk 2; for M6).
2. **The tie rule is not mirror-symmetric.** Flip augmentation must re-rasterize after flipping (lab risk 5; for M5).
3. **Pathological shapes at the caps** take about 1 s per call (they are legal input). If the reward ever sees such outputs, cap the crop side (≤ 512) earlier.
4. **No masks from `score_batch`.** M6 can call `rasterize` from Python threads (the GIL is released) or ask for a mask option.
5. **A `std::bad_alloc` inside a `score_batch` worker thread calls `std::terminate`.** Inputs are capped, so this needs about 41 MB per thread to fail.
6. **Not verified:** a real GitHub Actions run, and cibuildwheel on macOS.
