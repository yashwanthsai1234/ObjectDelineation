"""Gate B1 / Gate B3 on real final_eval outputs (Implementation.md B.6.2): a thin CLI over objdel.verify.rungates.

    .venv/bin/python -I -m objdel.cli.eval.gates b1 --b1 runs/rl-run/eval_val --b1-edge runs/rl-run/eval_edge \
        --r5 runs/r5/eval_val --view data/eval_views/lidc-val-shifted --edge-view data/eval_views/lidc-val-edge \
        [--a3 DIR --a3-edge DIR --identity holds|fails] [--pre runs/a3/val_pre_rebuild.json] [--partial-steps N] \
        --out runs/rl-run/gate_b1.json
    .venv/bin/python -I -m objdel.cli.eval.gates b3 --b3 DIR --b3-edge DIR --b1 DIR --b1-edge DIR \
        --view ... --edge-view ... --select val_select.json --out runs/b3/gate_b3.json

Each run folder is a final_eval output (items.jsonl + final_eval.json). --pre is the hash record of the view before its
rebuild: with no A3 baseline, Gate B1 uses R5's items only when its hash chain holds (rungates.hash_chain_ok); without --pre the chain
counts as broken. --select is a JSON list of val-select patient ids. Prints the result line, writes the full verdict.
Reads untrusted data: run with python -I.
"""

import argparse
import json
import sys
from pathlib import Path

from objdel.core.io import read_json, write_json
from objdel.policy.final_eval import FROZEN_VIEWS
from objdel.verify import rungates as G


def _run(path: Path | None, name: str) -> G.Run | None:
    if path is None:
        return None
    run = G.load_run(path)
    if run is None:
        sys.exit(f"run_gates: --{name} {path} holds no items.jsonl")
    return run


def main(argv=None) -> int:
    """The run_gates CLI (b1 or b3; see the module docstring)."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="gate", required=True)
    for gate in ("b1", "b3"):
        p = sub.add_parser(gate)
        p.add_argument("--view", type=Path, required=True, help="the frozen val view (lidc/val.jsonl inside)")
        p.add_argument("--edge-view", type=Path, required=True)
        p.add_argument("--frozen", type=Path, default=FROZEN_VIEWS)
        p.add_argument("--out", type=Path, required=True)
        p.add_argument("--b1", type=Path, required=True)
        p.add_argument("--b1-edge", type=Path, required=True)
    b1, b3 = sub.choices["b1"], sub.choices["b3"]
    b1.add_argument("--r5", type=Path, required=True)
    b1.add_argument("--a3", type=Path)
    b1.add_argument("--a3-edge", type=Path)
    b1.add_argument("--identity", choices=("holds", "fails"), help="whether A3 and R5 gave identical answers on the shared views")
    b1.add_argument("--pre", type=Path, help="pre-rebuild view hash record (the R5 hash chain)")
    b1.add_argument("--partial-steps", type=int)
    b3.add_argument("--b3", type=Path, required=True)
    b3.add_argument("--b3-edge", type=Path, required=True)
    b3.add_argument("--select", type=Path, required=True)
    a = ap.parse_args(argv)

    frozen = read_json(a.frozen)
    meta, edge_meta = G.view_meta(a.view), G.view_meta(a.edge_view)
    b1_run, b1_edge = _run(a.b1, "b1"), _run(a.b1_edge, "b1-edge")
    if a.gate == "b1":
        chain_ok = a.pre is not None and G.hash_chain_ok(read_json(a.pre), frozen.get(b1_run.view, {}))
        verdict = G.evaluate_b1(
            b1_run, b1_edge, a3=_run(a.a3, "a3"), a3_edge=_run(a.a3_edge, "a3-edge"), r5=_run(a.r5, "r5"),
            meta=meta, edge_meta=edge_meta, frozen=frozen, chain_ok=chain_ok,
            identity_ok=None if a.identity is None else a.identity == "holds", partial_steps=a.partial_steps)
    else:
        verdict = G.evaluate_b3(_run(a.b3, "b3"), _run(a.b3_edge, "b3-edge"), b1=b1_run, b1_edge=b1_edge, meta=meta,
                                edge_meta=edge_meta, frozen=frozen, select=set(read_json(a.select)))
    write_json(a.out, verdict)
    print(json.dumps({"gate": a.gate, "result": verdict["result"], "rules": verdict["rules"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
