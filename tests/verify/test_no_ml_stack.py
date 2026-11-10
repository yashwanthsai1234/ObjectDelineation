"""E3 gate (R-3): the verifier side (geom, verify, data) loads no ML stack and no training package."""
import subprocess
import sys

CODE = """
import importlib, pkgutil, sys
import objdel.data, objdel.geom, objdel.verify
for pkg in (objdel.geom, objdel.verify, objdel.data):
    for m in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
        if m.name != "objdel.verify.http":  # the FastAPI face (extra `serve`), not the library
            importlib.import_module(m.name)
ml = {"torch", "trl", "transformers", "vllm", "peft", "datasets", "accelerate"} & set(sys.modules)
up = {m for m in sys.modules if m.split(".")[:2] in (["objdel", "train"], ["objdel", "eval"])}
print(sorted(ml | up))
"""


def test_verifier_side_imports_no_ml_stack():
    out = subprocess.run([sys.executable, "-I", "-c", CODE], capture_output=True, text=True, check=True).stdout
    assert out.strip() == "[]"
