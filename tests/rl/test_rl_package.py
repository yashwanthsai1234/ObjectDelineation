"""U40 (E5.1): the objdel.rl package is lazy. Its __init__ holds a docstring only, so importing the rows module loads no
TRL and no vLLM (sampler-side and data-side callers import objdel.rl.rows without the trainer stack)."""
import ast
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_rl_init_holds_a_docstring_only():
    body = ast.parse((ROOT / "src/objdel/rl/__init__.py").read_text()).body
    assert len(body) == 1 and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant)
    assert isinstance(body[0].value.value, str)


def test_importing_rl_rows_loads_no_trl_and_no_vllm():
    code = "import sys, objdel.rl.rows; print(sorted({'trl', 'vllm'} & {m.split('.')[0] for m in sys.modules}))"
    out = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "[]"
