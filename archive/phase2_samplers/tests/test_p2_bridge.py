"""Archived with experiments/phase2/sample_vllm.py (E4 step 3): tests of the script itself."""


def test_sampler_imports_nothing_from_objdel():
    """It runs in /opt/vllm (vLLM's own torch / numpy pins, no objdel): top-level imports are the standard library;
    vLLM, torch and PIL are imported inside main."""
    tree = ast.parse(SCRIPT.read_text())
    mods = {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
    mods |= {n.module.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
    assert not any(isinstance(n, ast.ImportFrom) and n.level for n in ast.walk(tree))  # no relative imports
    assert mods <= set(sys.stdlib_module_names) | {"vllm", "torch", "PIL"}
    top = {a.name for n in tree.body if isinstance(n, ast.Import) for a in n.names}
    top |= {n.module for n in tree.body if isinstance(n, ast.ImportFrom)}
    assert top <= set(sys.stdlib_module_names)




def test_sampler_help_without_vllm():
    r = subprocess.run([sys.executable, "-I", str(SCRIPT), "--help"], capture_output=True, text=True, timeout=60)
    assert r.returncode == 0 and "--prompts" in r.stdout and "--gpu-memory-utilization" in r.stdout


