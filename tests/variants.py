"""CPU test variants of objdel._gpucheck (plan D2 step 7): the same source built with test-only compile definitions
(CMake cache string OBJDEL_TEST_HOOKS), each a non-editable install in build/variants/<name>/site.

A variant is never put on PYTHONPATH: `-I` ignores it, and without `-I` the editable install's finder (first on
sys.meta_path) still serves the installed extension. A variant runs in a `python -I -S` child whose sys.path is set
inside it (its site folder, then the repo root, then this interpreter's site-packages), and the child's load guard
checks that objdel._gpucheck comes from that site folder and carries exactly the variant's definitions (exit 3
otherwise), so a wrong load is red, never a silent pass.

`python -I -m tests.variants` (make test-order) builds every variant and applies the pass rules: order1, order2 and
poison pass the mask tests; alias (a planted scan alias under the poison) must fail them, which proves the poison
check can go red; mutfg4 is used by tests/geom/test_variants.py, run with OBJDEL_REQUIRE_VARIANTS=1. Not collected by
pytest (no test_ prefix)."""
import json
import os
import site
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = {"order1": ("K1_ORDER=1",), "order2": ("K1_ORDER=2",), "poison": ("GC_POISON_SCAN",),
            "alias": ("GC_POISON_SCAN", "GC_ALIAS_SCAN_FOR_TEST"), "mutfg4": ("K1_MUTATE_FG4",)}
MASK_TESTS = ("-q", "-p", "no:cacheprovider", "tests/geom/test_cuda_checker.py", "-k", "masks",
              "-m", "not slow and not data")
EXPECT = {"order1": 0, "order2": 0, "poison": 0, "alias": 1}  # pytest exit codes: 0 passed, 1 tests failed

# Run inside the child: set sys.path, check the load, then run pytest or a code string.
_BOOT = """
import json, sys
from pathlib import Path
cfg = json.loads(sys.argv[1])
sys.path[:0] = [cfg["site"], cfg["root"]]
sys.path += cfg["paths"]
import objdel._gpucheck as g
if not Path(g.__file__).resolve().is_relative_to(Path(cfg["site"]).resolve()) or \\
        set(g.test_hooks()) != set(cfg["hooks"]):
    print("variant load guard:", g.__file__, g.test_hooks(), file=sys.stderr)
    sys.exit(3)
if cfg["pytest"] is not None:
    import pytest
    sys.exit(pytest.main(cfg["pytest"]))
exec(cfg["code"])
"""


def site_dir(name: str) -> Path:
    return ROOT / "build" / "variants" / name / "site"


def built(name: str) -> bool:
    return any((site_dir(name) / "objdel").glob("_gpucheck*"))


def build(name: str, defs=None) -> None:
    """Non-editable install of the package with the variant's definitions into build/variants/<name>/site."""
    defs = VARIANTS[name] if defs is None else defs
    subprocess.run(["uv", "pip", "install", "--python", sys.executable, "--no-deps", "--reinstall-package", "objdel",
                    "--target", str(site_dir(name)), str(ROOT), "-C", f"build-dir={site_dir(name).parent / 'cmake'}",
                    "-C", "cmake.define.OBJDEL_TEST_HOOKS=" + ";".join(defs)], check=True, cwd=ROOT)


def _child(name: str, hooks, pytest_args=None, code: str = "") -> subprocess.CompletedProcess:
    cfg = dict(site=str(site_dir(name)), root=str(ROOT), paths=site.getsitepackages(), hooks=sorted(hooks),
               pytest=None if pytest_args is None else list(pytest_args), code=code)
    return subprocess.run([sys.executable, "-I", "-S", "-c", _BOOT, json.dumps(cfg)], cwd=ROOT, capture_output=True,
                          text=True)


def run(name: str, pytest_args=MASK_TESTS) -> int:
    """pytest in the variant's child: its exit code (3 = the load guard failed)."""
    p = _child(name, VARIANTS[name], pytest_args)
    sys.stdout.write(p.stdout[-3000:]), sys.stderr.write(p.stderr[-3000:])
    return p.returncode


def run_code(name: str, code: str, hooks) -> str:
    """A code string in the variant's child: its stdout. Raises on any failure, the load guard's included."""
    p = _child(name, hooks, code=code)
    if p.returncode:
        raise RuntimeError(f"variant {name} exited {p.returncode}: {p.stderr[-2000:]}")
    return p.stdout


def main() -> int:
    bad = []
    for name in VARIANTS:
        build(name)
    for name, want in EXPECT.items():
        got = run(name)
        print(f"variant {name}: pytest exit {got} (want {want})", flush=True)
        if got != want:
            bad.append(name)
    env = {**os.environ, "OBJDEL_REQUIRE_VARIANTS": "1"}
    got = subprocess.run([sys.executable, "-I", "-m", "pytest", "-q", "-p", "no:cacheprovider",
                          "tests/geom/test_variants.py"], cwd=ROOT, env=env).returncode
    print(f"variant mutfg4: tests/geom/test_variants.py exit {got} (want 0)")
    if got:
        bad.append("mutfg4")
    print("test-order:", "FAIL " + " ".join(bad) if bad else "OK")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
