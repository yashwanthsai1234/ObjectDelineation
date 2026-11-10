"""Hand edits that follow `e_moves` on the tree of commit a255988 (the cases it reports as MANUAL).

Each edit is an exact (old, new) pair tagged with its sub-move; an edit whose old text is missing is reported and
skipped, so a re-run on a later tree shows what still needs a hand edit. Not reported: an edit whose new text a later
due edit of the same file starts from (that edit carries it on, as E3c's audit import after E3b's). Run after e_moves
with the same stage:

    .venv/bin/python -I -m experiments.phase2.tools.e_moves_manual --stage E3a
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from experiments.phase2.tools.e_moves import STAGES

EDITS: list[tuple[str, str, str, str]] = [
    # fixtures that returned a split module
    ("E5", "tests/experiments/test_p2_breakit_data.py", 'pytest.importorskip("datasets")\n    return dataset',
     'pytest.importorskip("datasets")\n    from objdel.sft import rows\n\n    return rows'),
    ("E5", "tests/experiments/test_p2_breakit_scoring.py", 'pytest.importorskip("trl")\n    return sample',
     'pytest.importorskip("trl")\n    from objdel.policy import sampler\n\n    return sampler'),
    # module docstring of configs now lives in policy/models.py (docstring_to)
    ("E5", "tests/experiments/test_repo_hygiene.py", 'def test_one_stated_rl_learning_rate():\n\n    assert "prereg.G_PROBE" in configs.__doc__ '
     'and "1e-6 full, 1e-5 LoRA" not in configs.__doc__',
     'def test_one_stated_rl_learning_rate():\n    from objdel.policy import models\n\n    assert "prereg.G_PROBE" in '
     'models.__doc__ and "1e-6 full, 1e-5 LoRA" not in models.__doc__'),
    # names that are gone (negative checks) now look in the new homes
    ("E3b", "tests/experiments/test_train_rl_guards.py", 'not hasattr(serve, "VIEW_SEP")\n    assert not hasattr(serve, "_tagged")',
     'not hasattr(store, "VIEW_SEP")\n    assert not hasattr(store, "_tagged")'),
    ("E3b", "tests/experiments/test_train_rl_guards.py", '    from objdel.train import dataset\n\n    assert "--train-views"',
     '    from objdel.train import dataset\n    from objdel.verify import store\n\n    assert "--train-views"'),
    # the move drops a module import that is only partly due (the fixture still returns the old module until E5)
    ("E3b", "tests/experiments/test_p2_breakit_data.py", 'pytest.importorskip("datasets")\n    return dataset',
     'pytest.importorskip("datasets")\n    from objdel.train import dataset\n    return dataset'),
    ("E5", "tests/experiments/test_outline_lengths.py", 'assert not hasattr(dataset, "_tokens") and not hasattr(dataset, "_voting")',
     'assert not any(hasattr(m, n) for m in (policy_rows, verify_store) for n in ("_tokens", "_voting"))'),
    ("E3b", "tests/rl/test_train_trl.py",
     '    assert all(hasattr(objdel.eval, n) for n in ("hack_audit", "gates", "metrics", "summarize", "vision"))\n'
     '    assert not hasattr(objdel.eval, "alarms" + "_for")',
     '    from objdel.eval import audit\n\n    assert callable(audit.hack_audit)\n'
     '    assert not hasattr(audit, "alarms" + "_for")'),
    ("E3c", "tests/rl/test_train_trl.py", "    from objdel.eval import audit\n\n    assert callable(audit.hack_audit)\n",
     "    from objdel.data import audit\n\n    assert callable(audit.hack_audit)\n"),
    ("E3b", "tests/rl/test_train_trl.py", "from pathlib import Path\n\n    import objdel.eval\n    from objdel.train import callbacks\n",
     "from pathlib import Path\n\n    from objdel.train import callbacks\n"),
    # code inside strings (subprocess / planted sources)
    ("E3a", "tests/geom/test_variants.py", 'CERTIFY = f"""\nfrom objdel.checker import certify\n',
     'CERTIFY = f"""\nfrom objdel.geom.certify import certify\n'),
    ("E5", "tests/experiments/test_p2_train_rl.py", "from objdel.train import rewards\n", "from objdel.rl import reward as rewards\n"),
    ("E5", "tests/rl/test_trl_adapter.py", '"from objdel.train import trl_adapter as T\\n"', '"from objdel.rl import trl_adapter as T\\n"'),
    ("E5", "tests/experiments/test_p2_eval.py", '"from objdel.eval import vision; print(vision.decide(1, 1, 1)[\'action\'])"',
     '"from objdel.policy import vision; print(vision.decide(1, 1, 1)[\'action\'])"'),
    ("E5", "tests/experiments/test_nvtx.py", '"\\n\\ndef _planted():\\n    from objdel.train import nvtx\\n    return nvtx\\n"',
     '"\\n\\ndef _planted():\\n    from objdel.policy import marks\\n    return marks\\n"'),
    ("E5", "tests/experiments/test_nvtx.py", 'replace("from objdel.train import nvtx", "from ...train import nvtx")',
     'replace("from objdel.policy import marks", "from ...policy import marks")'),
    # the NVTX name scan walked train/; its code now sits in policy/ and rl/
    ("E5", "tests/experiments/test_nvtx.py", "    import objdel.train\n\n    names, bad = set(), []\n"
     "    for f in Path(objdel.train.__file__).parent.glob(\"*.py\"):",
     "    import objdel.policy\n    import objdel.rl\n\n    names, bad = set(), []\n"
     "    pkgs = (Path(objdel.policy.__file__).parent, Path(objdel.rl.__file__).parent)\n"
     "    for f in (f for d in pkgs for f in d.glob(\"*.py\")):"),
    # code in strings imported the package re-exports that E3b removes (R-1): name the defining modules
    ("E3b", "tests/verify/test_reward.py", 'from objdel.verify.reward import ExpertCache\\n"',
     'from objdel.verify.reward.cache import ExpertCache\\n"'),
    ("E3b", "tests/verify/test_reward.py", "from objdel.verify.reward import NUMERIC, ExpertCache, advantages_v4, score_batch\n",
     "from objdel.verify.reward.advantages import advantages_v4\nfrom objdel.verify.reward.cache import ExpertCache\n"
     "from objdel.verify.reward.score import NUMERIC, score_batch\n"),
    ("E3b", "tests/experiments/test_p2_eval.py", '"from objdel.eval import gates, metrics, vision; print(',
     '"from objdel.eval import vision; print('),
    # the seam test: the package API is geom (polygonize) plus the defining modules
    ("E3a", "tests/experiments/test_rings.py", "        from objdel.checker import bitquads  # noqa: F401\n",
     "        from objdel.geom import bitquads  # noqa: F401\n"),
    ("E3a", "tests/experiments/test_rings.py",
     "        assert name in objdel.checker.__all__ and hasattr(objdel.checker, name), name\n",
     "        assert any(hasattr(importlib.import_module(f\"objdel.geom.{m}\"), name)\n"
     "                   for m in (\"rings\", \"certify\", \"topology\")), name\n"
     "    assert objdel.geom.__all__ == [\"polygonize\"]\n"),
    ("E3a", "tests/experiments/test_rings.py", "    for owner in (rings, objdel.checker):", "    for owner in (rings, objdel.geom):"),
]

# README's example and a subprocess string imported package re-exports, which phase E removes (R-1): name the defining
# modules
EDITS.append(("E3c", "tests/data/test_targets.py", "from objdel.data.targets import build_target;",
              "from objdel.data.targets.build import build_target;"))
EDITS.append(("E3c", "README.md", "from objdel.data.targets import build_target, min_areas\n",
              "from objdel.data.targets.build import build_target, min_areas\n"))

# the one test without the fixture keeps the module alias
EDITS.append(("E5", "tests/experiments/test_p2_breakit_data.py", "t = D._no_holes(", "t = sft_rows._no_holes("))

# every occurrence: the tests take the fixture `D` but the alias rewrite named the module
REPLACE_ALL: list[tuple[str, str, str, str]] = [
    ("E5", "tests/experiments/test_p2_breakit_data.py", "sft_rows.view_rows(", "D.view_rows("),
    # the render file of a split module: window_record and prepare (the frames) now live in data/prepare.py; the
    # frozen records keep their recorded commits (freeze_views --check then lists the file as moved, which it is)
    ("E3c", "experiments/phase2/freeze_views.py", '"src/objdel/train/dataset.py"', '"src/objdel/data/prepare.py"'),
    ("E3c", "docs/reference/frozen_views.json", '"src/objdel/train/dataset.py"', '"src/objdel/data/prepare.py"'),
    # a skip on a split module would silently skip forever
    ("E3b", "tests/experiments/test_p2_probes.py", 'pytest.importorskip("objdel.train.serve")', 'pytest.importorskip("objdel.verify.store")'),
    # source scans that globbed the removed packages now glob their new homes
    ("E5", "tests/experiments/test_final_eval.py", '    """I28 (grep): the scoring code (experiments/phase2/*.py, objdel.eval, objdel.train) names the test view only in\n    final_eval (freeze_views hashes it, never scores it), and every config line that names it is a final_eval step."""\n    hits = [p.relative_to(ROOT).as_posix() for p in [*(ROOT / "experiments/phase2").glob("*.py"),\n                                                       *(ROOT / "src/objdel/eval").glob("*.py"),\n                                                       *(ROOT / "src/objdel/train").glob("*.py")]\n',
     '    """I28 (grep): the scoring code (experiments/phase2/*.py; objdel.verify, policy, rl, sft; data audit, probes and\n    prepare) names the test view only in final_eval (freeze_views hashes it, never scores it), and every config line\n    that names it is a final_eval step."""\n    lib = [*(f for d in ("verify", "policy", "rl", "sft") for f in (ROOT / "src/objdel" / d).rglob("*.py")),\n           *(ROOT / "src/objdel/data" / f for f in ("audit.py", "probes.py", "prepare.py"))]\n    hits = [p.relative_to(ROOT).as_posix() for p in [*(ROOT / "experiments/phase2").glob("*.py"), *lib]\n'),
    # S12: the verify layer never imports the ML layers (policy, rl, sft; import-linter checks it too)
    ("E5", "tests/experiments/test_nvtx.py", 'if any(m == "objdel.train" or m.startswith("objdel.train.") for m in mods):',
     'if any(m.split(".")[:2] in (["objdel", "policy"], ["objdel", "rl"], ["objdel", "sft"]) for m in mods):'),
    ("E5", "tests/experiments/test_nvtx.py", "# S12: objdel.verify.reward never imports objdel.train ----------",
     "# S12: objdel.verify.reward never imports objdel.policy/rl/sft"),
    ("E5", "tests/experiments/test_nvtx.py", '"""Imports of objdel.train anywhere', '"""Imports of objdel.policy, rl or sft anywhere'),
]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stage", choices=STAGES, default=STAGES[-1], help="apply the edits up to this sub-move")
    due = set(STAGES[: STAGES.index(ap.parse_args(argv).stage) + 1])
    root = Path.cwd()
    missed = 0
    for i, (stage, rel, old, new) in enumerate(EDITS):
        if stage not in due:
            continue
        p = root / rel
        text = p.read_text() if p.exists() else ""
        if new in text:  # done already (also when new contains old)
            continue
        if old in text:
            p.write_text(text.replace(old, new, 1))
        elif not any(s in due and r == rel and o in new for s, r, o, _ in EDITS[i + 1:]):  # a later edit took it on
            print(f"MISSING: {rel}: {old[:70]!r}")
            missed += 1
    for stage, rel, old, new in REPLACE_ALL:
        if stage not in due:
            continue
        p = root / rel
        text = p.read_text()
        if old not in text and new not in text:
            print(f"MISSING: {rel}: {old[:70]!r}")
            missed += 1
        p.write_text(text.replace(old, new))
    return 1 if missed else 0


if __name__ == "__main__":
    sys.exit(main())
