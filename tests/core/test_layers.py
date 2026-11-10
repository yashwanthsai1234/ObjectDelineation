"""E0: AST scans that import-linter and ruff cannot do (private TRL names, host syncs, module strings, LR)."""
import ast
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC, TESTS, EXPERIMENTS = ROOT / "src/objdel", ROOT / "tests", ROOT / "experiments"
GOLDENS = TESTS / "goldens/phase_e"

PRIVATE = {"_get_per_token_logps_and_entropies", "_generate_and_score_completions", "_prepare_inputs", "_compute_loss",
           "_generate", "_calculate_rewards", "shuffle_sequence_dict", "_iter_named_params", "_metrics", "_step",
           "_buffered_inputs"}
# Allowed to touch private TRL names. E5 leaves only rl/trl_adapter.py (Implementation.md E5 step 4).
TRL_ADAPTER_FILES = {"train/fast.py", "train/trl_adapter.py", "train/sample.py"}
TRL_BASES = {"GRPOTrainer", "VLLMGeneration"}
HOST_SYNC = re.compile(r"\.(item|cpu|tolist|numpy)\(\)")


def parse(path: Path) -> ast.AST:
    return ast.parse(path.read_text(), filename=str(path))


def py_files(root: Path):
    return sorted(p for p in root.rglob("*.py") if "__pycache__" not in p.parts)


def rel(path: Path, root: Path = SRC) -> str:
    return path.relative_to(root).as_posix()


def imported_roots(tree: ast.AST) -> set[str]:
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            out.add(node.module)
    return out


def trl_aliases(tree: ast.AST) -> set[str]:
    """Names bound to a module whose dotted path starts with `trl` (import trl.x as y; from trl... import module)."""
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names |= {a.asname or a.name.split(".")[0] for a in node.names if a.name.split(".")[0] == "trl"}
        elif isinstance(node, ast.ImportFrom) and (node.module or "").split(".")[0] == "trl":
            names |= {a.asname or a.name for a in node.names}
    return names


def private_trl_uses(tree: ast.AST) -> list[int]:
    """Lines that reach a private TRL name: self.<n> in a GRPOTrainer/VLLMGeneration subclass, trainer.<n>,
    super().<n>, or <trl module alias>.<n>."""
    aliases, lines = trl_aliases(tree), []

    def visit(node, in_trl_class):
        if isinstance(node, ast.ClassDef):
            bases = {b.attr if isinstance(b, ast.Attribute) else getattr(b, "id", "") for b in node.bases}
            in_trl_class = bool(bases & TRL_BASES)
        if isinstance(node, ast.Attribute) and node.attr in PRIVATE:
            v = node.value
            if isinstance(v, ast.Name) and (v.id == "trainer" or v.id in aliases or (v.id == "self" and in_trl_class)):
                lines.append(node.lineno)
            elif isinstance(v, ast.Call) and getattr(v.func, "id", "") == "super":
                lines.append(node.lineno)
        for child in ast.iter_child_nodes(node):
            visit(child, in_trl_class)

    visit(tree, False)
    return lines


def test_no_private_trl_names_outside_adapter():
    bad = {rel(p): private_trl_uses(parse(p)) for p in py_files(SRC) if rel(p) not in TRL_ADAPTER_FILES}
    assert {k: v for k, v in bad.items() if v} == {}


def test_private_trl_scan_flags_each_form():
    src = ("import trl.trainer.grpo_trainer as g\nfrom trl import GRPOTrainer\n"
           "class A(GRPOTrainer):\n    def f(self):\n        self._compute_loss()\n        super()._generate()\n"
           "class B:\n    def f(self):\n        self._compute_loss()\n"  # our own class: allowed
           "def h(trainer):\n    trainer._step\n    g.shuffle_sequence_dict\n")
    assert private_trl_uses(ast.parse(src)) == [5, 6, 11, 12]


def torch_modules() -> list[Path]:
    return [p for p in py_files(SRC) if any(m.split(".")[0] == "torch" for m in imported_roots(parse(p)))]


def host_sync_counts() -> dict[str, int]:
    counts = {rel(p): len(HOST_SYNC.findall(p.read_text())) for p in torch_modules()}
    return {k: v for k, v in counts.items() if v}


def test_host_sync_sites_not_raised():
    recorded = json.loads((GOLDENS / "host_syncs.json").read_text())
    now = host_sync_counts()
    assert {k: v for k, v in now.items() if v > recorded.get(k, 0)} == {}


def files_importing_experiments() -> list[str]:
    out = []
    for p in py_files(TESTS):
        if p.relative_to(TESTS).parts[0] != "experiments" and any(
            m.split(".")[0] == "experiments" for m in imported_roots(parse(p))
        ):
            out.append(rel(p, TESTS))
    return out


def test_tests_never_import_experiments():
    """The recorded list only shrinks (0 at E8)."""
    recorded = set(json.loads((GOLDENS / "tests_import_experiments.json").read_text()))
    assert set(files_importing_experiments()) <= recorded


def test_src_never_imports_experiments_or_tests():
    bad = {rel(p): sorted(m for m in imported_roots(parse(p)) if m.split(".")[0] in {"experiments", "tests", "oracles"})
           for p in py_files(SRC)}
    assert {k: v for k, v in bad.items() if v} == {}


def top_level_names(path: Path) -> set[str]:
    names = set()
    for node in parse(path).body:
        if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            names |= {t.id for t in node.targets if isinstance(t, ast.Name)}
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
        elif isinstance(node, ast.ImportFrom) and (path.name == "__init__.py" or path == SRC / "geom/engine.py"):
            # a package's re-exports are its names; geom/engine re-exports the CUDA extension (E3a)
            names |= {a.asname or a.name for a in node.names}
    return names


MODULE_STRING = re.compile(r"\bobjdel(?:\.\w+)+")


def module_imports(path: Path) -> set[str]:
    """Names a module binds with a top-level `import x` (a monkeypatch target may go through them)."""
    return {a.asname or a.name for n in parse(path).body if isinstance(n, ast.Import) for a in n.names}


def resolves(dotted: str) -> bool:
    """The longest prefix that is a module file exists, and the rest (if any) names a top-level definition in it."""
    parts = dotted.split(".")[1:]
    if parts[0] == "_gpucheck":  # the compiled extension
        return True
    for n in range(len(parts), 0, -1):
        base = SRC.joinpath(*parts[:n])
        for f in (base.with_suffix(".py"), base / "__init__.py"):
            if f.exists():
                rest = parts[n:]
                names = top_level_names(f) | module_imports(f)
                return not rest or rest[0] in names or (base / f"{rest[0]}.py").exists()
    return False


def string_tokens(path: Path):
    for node in ast.walk(parse(path)):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            yield from MODULE_STRING.findall(node.value)


def unresolved_module_strings() -> list[str]:
    files = [*py_files(SRC), *py_files(EXPERIMENTS), *py_files(TESTS)]
    live = [p for p in files if "lab" not in p.relative_to(ROOT).parts]
    return sorted({f"{p.relative_to(ROOT).as_posix()}: {t}" for p in live for t in string_tokens(p) if not resolves(t)})


def test_module_string_through_an_imported_module_resolves():
    """A monkeypatch target such as objdel.io.subprocess.run names a module that objdel.io imports."""
    assert resolves("objdel.io.subprocess.run") and not resolves(".".join(("objdel", "io", "no_such_name")))


def test_module_strings_resolve():
    """Every dotted objdel name in a string or docstring of src/, live experiments/ and tests/ names a module or a
    definition. The recorded stale tokens (docstrings of retired names) only shrink."""
    stale = set(json.loads((GOLDENS / "module_strings_stale.json").read_text()))
    assert set(unresolved_module_strings()) <= stale


def test_g_probe_lr_not_read():
    """The pre-registered record `prereg.G_PROBE` is read for no learning rate; `configs.rl_lr` is the one source."""
    bad = []
    for p in [*py_files(SRC), *py_files(EXPERIMENTS)]:
        if p.name == "prereg.py":
            continue
        for node in ast.walk(parse(p)):
            if (isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name) and node.value.id == "G_PROBE"
                    and isinstance(node.slice, ast.Constant) and node.slice.value in {"lr_full", "lr_lora"}):
                bad.append(f"{p.relative_to(ROOT).as_posix()}:{node.lineno}")
    assert bad == []


# E-G16 baseline (U10): the old modules whose top-level names must each find a home after the E moves. A listed
# module that does not exist yet (an A-D module of another lane) is recorded once it lands.
OLD_MODULES = [*(f"src/objdel/train/{m}.py" for m in ("dataset", "configs", "fast", "rewards", "callbacks", "sample",
                                                       "serve", "nvtx", "trl_adapter", "gpu_pixels", "vllm_sync",
                                                       "vllm_requests")),
               *(f"src/objdel/{m}.py" for m in ("eval/ceiling", "eval/select", "eval/device", "checker/rings",
                                                "geometry/box")),
               *(f"experiments/phase2/{m}.py" for m in ("select_ckpt", "train_sft", "g_probe", "train_rl",
                                                        "final_eval", "vast/policy", "vast/bundle", "preflight",
                                                        "hf_push", "hf_pull"))]


def old_names() -> dict[str, list[str]]:
    return {m: sorted(top_level_names(ROOT / m)) for m in OLD_MODULES if (ROOT / m).exists()}


def new_names() -> set[str]:
    """Top-level names of every library module (the homes of the moved names)."""
    return set().union(*(top_level_names(p) for p in py_files(SRC)))


def record_old_names():
    """Re-record the modules that still exist; a module the E moves removed keeps its recorded names."""
    path = GOLDENS / "old_names.json"
    path.write_text(json.dumps({**json.loads(path.read_text()), **old_names()}, indent=1) + "\n")


def test_old_names_recorded_on_this_tree():
    """Re-record after a merge that changes an old module: record_old_names() (goldens README)."""
    recorded, now = json.loads((GOLDENS / "old_names.json").read_text()), old_names()
    assert {m: v for m, v in recorded.items() if m in now} == now


def test_move_complete():
    """E-G16: each recorded name of an old module that the E moves removed has a home in a library module."""
    recorded, now, homes = json.loads((GOLDENS / "old_names.json").read_text()), old_names(), new_names()
    lost = {m: sorted(set(v) - homes) for m, v in recorded.items() if m not in now}
    assert {m: v for m, v in lost.items() if v} == {}
