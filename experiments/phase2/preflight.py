"""Preflight: everything that can be proven on the Mac, before any GPU is rented (the launch gate reads its stamp).

    python -I -m experiments.phase2.preflight --config experiments/phase2/configs/r1_lfm.yaml     (make preflight CONFIG=...)

Checks, in order; the stage (DATA / CHECKER / CODE) with a failure is finished, the later stages are not run:
  DATA     bundle integrity (rows.jsonl and context hashes, eval-view file hashes), patients shared between training and
           any exam (or probe) view, coverage of every exam view by the training view (objdel.eval.coverage), prompt
           text identical, CT window consistent, every training target parses and is legal for the reference checker
  CHECKER  tests/test_cuda_checker*.py when they exist, the scorers on the CUDA checker (tests/test_gpu_scoring.py:
           equal to the host scorer it replaced), then the checker's exactness tests (the slow exact oracles, the
           archived CPU checker's recorded results, topology, geometry)
  CODE     the private TRL names fast.py relies on (objdel.train.trl_adapter), a quick test subset, the CPU dry run
           of train -> select on a tiny model, the projected GPU hours and dollars
Writes runs/preflight/<config>.json: pass per check and the code, data and config hashes. vast launch (policy.py)
refuses unless that stamp passes and all three hashes still match: refusal(). Exit 0 = every check passed.
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# the exam views policy.chain passes to select_ckpt (gates: shifted, rescaled; report: rotated; probes)
EXAMS = ("lidc-val-shifted", "lidc-val-rescaled", "lidc-val-rotated", "probe-precision", "probe-topology")
CHECKER_TESTS = ["tests/test_checker.py", "tests/test_checker_oracle.py", "tests/test_checker_host.py",
                 "tests/test_topology.py", "tests/test_geometry.py", "tests/test_geometry_mask.py",
                 "tests/test_geometry_oracle.py"]  # fmt: skip
QUICK_TESTS = ["tests/test_eval_coverage.py", "tests/test_preflight.py", "tests/test_p2_scripts.py",
               "tests/test_train_configs.py", "tests/test_vast.py", "tests/test_p2_v1_select.py"]  # fmt: skip
DRY_RUN = "tests/test_p2_scripts.py::test_cpu_smoke_train_then_select"  # train -> select on a tiny model, CPU
HF_RATE = {"lfm2.5-vl-450m": "hf-rate-lfm", "qwen3-vl-2b": "hf-rate-qwen3vl"}  # runs/<run>/<dir>/hf_rate.json (R0)


def ceiling(gpu: str) -> float:
    """$/h price ceiling of the card (policy.CARDS): a bound, not a quote."""
    from experiments.phase2.vast import policy

    return policy.CARDS[gpu]["ceiling"]


def sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(cfg: Path):
    """(RunArgs, raw config dict) of a run config."""
    import yaml

    from experiments.phase2 import train_sft as TS

    raw = yaml.safe_load(Path(cfg).read_text())
    return TS.RunArgs(**raw), raw


def hashes(cfg: Path) -> dict:
    """code: as run manifests (train_sft.code_hashes); data: the training view hash, each exam view's manifest and the
    CUDA checker eval set's (it ships in the LIDC bundle);
    config: the file's bytes."""
    from experiments.phase2 import train_sft as TS

    a, _ = load(cfg)
    root = Path(a.data_root)
    view = json.loads((root / a.view / "manifest.json").read_text())["view_hash"]
    exams = [sha(root / "eval_views" / v / "manifest.json") for v in EXAMS] + [sha(root / "checker_eval" / "manifest.json")]
    return {"code": TS.code_hashes()["objdel_and_script_sha256"], "config": sha(cfg),
            "data": hashlib.sha256((view + "".join(exams)).encode()).hexdigest()}  # fmt: skip


def stamp_path(cfg: Path) -> Path:
    """runs/preflight/<config>.json under $OBJDEL_RUN_DIR (as train_sft)."""
    return ROOT / os.environ.get("OBJDEL_RUN_DIR", "runs") / "preflight" / f"{Path(cfg).stem}.json"


def refusal(cfg: Path) -> str | None:
    """Why a launch of this config must be refused (no stamp, a failed one, or code / data / config changed since), else None."""
    path = stamp_path(cfg)
    if not path.is_file():
        return f"no preflight stamp: run `make preflight CONFIG={cfg}` first"
    try:
        stamp, now = json.loads(path.read_text()), hashes(cfg)
    except (OSError, ValueError, KeyError) as e:
        return f"cannot check {path} against the current code, data and config: {e}"
    if not stamp["pass"]:
        return f"{path} records a FAILED preflight (failed: {[k for k, v in stamp['checks'].items() if not v['pass']]})"
    if stale := [k for k in now if now[k] != stamp["hashes"][k]]:
        return f"{path} is stale ({', '.join(stale)} changed since it passed): run `make preflight` again"
    return None


# --- DATA ---------------------------------------------------------------------------------------------------------


def lines(path):
    return [json.loads(x) for x in Path(path).read_text().splitlines()]


def view_integrity(a, cfg) -> list[str]:
    from objdel.train import dataset as D

    root, vdir = Path(a.data_root), Path(a.data_root) / a.view
    man, rows = json.loads((vdir / "manifest.json").read_text()), lines(vdir / "rows.jsonl")
    bad = []
    if sha(vdir / "rows.jsonl") != man["view_hash"]:
        bad.append(f"{vdir}/rows.jsonl does not match the manifest's view_hash")
    try:
        D._contexts(vdir, man, {r["item_id"] for r in rows})
    except ValueError as e:
        bad.append(str(e))
    for v in EXAMS:
        d = root / "eval_views" / v
        for name, f in json.loads((d / "manifest.json").read_text())["files"].items():
            if not (d / name).is_file() or sha(d / name) != f["sha256"]:
                bad.append(f"{d / name} does not match its manifest")
    return bad


def patients(a, cfg) -> list[str]:
    """No patient (group_id) and no CT series in both the training rows and an exam view; train and val prepared splits
    disjoint too."""
    from objdel.train import dataset as D

    root = Path(a.data_root)
    group = {r["item_id"]: r["group_id"] for s in ("train", "val") for r in D.records(root / "prepared", "lidc", s)}
    train_items = {r["item_id"] for r in lines(root / a.view / "rows.jsonl")}
    train_groups, series = {group.get(i) for i in train_items}, {i.split("/")[1] for i in train_items}
    bad = [f"{sum(i not in group for i in train_items)} training items have no patient"] if None in train_groups else []
    val = {r["group_id"] for r in D.records(root / "prepared", "lidc", "val")}
    if shared := train_groups & val:
        bad.append(f"{len(shared)} patients in training and prepared val, e.g. {sorted(shared)[:3]}")
    for v in EXAMS:
        rows = lines(root / "eval_views" / v / "lidc" / "val.jsonl")
        shared = train_groups & {r["group_id"] for r in rows}
        if shared or series & {r["item_id"].split("/")[1] for r in rows}:
            bad.append(f"{v}: {len(shared)} patients also in training, e.g. {sorted(shared)[:3]}")
    return bad


def exam_props(root: Path, v: str) -> list[dict]:
    from objdel.eval import coverage as C

    out = []
    for r in lines(root / "eval_views" / v / "lidc" / "val.jsonl"):
        out += [p for t in r["targets"] if (p := C.props(t, r["W"], r["view"].get("op", "id")))]
    return out


def coverage(a, cfg) -> list[str]:
    from objdel.eval import coverage as C

    root = Path(a.data_root)
    train = [p for r in lines(root / a.view / "rows.jsonl") if r["epoch"] == 0
             and (p := C.props(r["target"], r["side"], r["op"]))]  # fmt: skip
    td, bad = C.distribution(train), []
    for v in EXAMS:
        for g in C.gaps(td, C.distribution(exam_props(root, v))):
            bad.append(f"{v}: {g['key']}={g['category']} is {g['exam']:.1%} of the exam, {g['train']:.2%} of training "
                       f"(floor {C.TRAIN_MAX:.0%}, exam share >= {C.EVAL_MIN:.0%})")  # fmt: skip
    return bad


def prompts(a, cfg) -> list[str]:
    from objdel.train import dataset as D

    want = next(c["text"] for c in D._user("lidc")["content"] if c.get("type") == "text")
    return [f"{v}: exam prompt differs from the training prompt"
            for v in EXAMS if json.loads((Path(a.data_root) / "eval_views" / v / "manifest.json").read_text())
            .get("prompts", {"lidc": want})["lidc"] != want]  # fmt: skip  (probe manifests carry no prompt)


def ct_window(a, cfg) -> list[str]:
    """Every training window is a proper (lo < hi) window, and the exams' fixed lung window lies in the 1-99 percentile
    range the training jitter covers."""
    import numpy as np

    from objdel.targets.frame import LUNG_WINDOW

    w = np.array([r["window"] for r in lines(Path(a.data_root) / a.view / "rows.jsonl")])
    bad = ["training rows with lo >= hi"] if (w[:, 0] >= w[:, 1]).any() else []
    for i, name in enumerate(("lo", "hi")):
        lo, hi = np.percentile(w[:, i], [1, 99])
        if not lo <= LUNG_WINDOW[i] <= hi:
            bad.append(f"exam window {name} {LUNG_WINDOW[i]} outside training 1-99 % range [{lo:.0f}, {hi:.0f}]")
    return bad


def targets(a, cfg) -> list[str]:
    from objdel.geometry import parse, to_mask

    bad = []
    for r in lines(Path(a.data_root) / a.view / "rows.jsonl"):
        p = parse(r["target"], ["nodule"])
        if not p.ok or not to_mask(p.objs["nodule"], r["side"], r["side"]).legal:
            bad.append(f"{r['row_id']}: {p.error or 'illegal'}")
    return [f"{len(bad)} training targets rejected by the reference checker, e.g. {bad[:3]}"] if bad else []


# --- CHECKER, CODE ------------------------------------------------------------------------------------------------


def pytest(*args: str) -> list[str]:
    p = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", *args], cwd=ROOT,
                       capture_output=True, text=True)
    return [] if p.returncode == 0 else ["pytest " + " ".join(args)] + p.stdout.splitlines()[-15:]


def trl_problems() -> list[str]:
    """The private TRL names fast.py relies on (objdel.train.trl_adapter): [] or the contract error's message."""
    from objdel.train import trl_adapter

    try:
        trl_adapter.check()
    except trl_adapter.TRLContractError as e:
        return [str(e)]
    return []


def checker_tests(a, cfg) -> list[str]:
    cuda = sorted(p.relative_to(ROOT).as_posix() for p in (ROOT / "tests").glob("test_cuda_checker*.py"))
    return pytest("-m", "not slow", *cuda, "tests/test_gpu_scoring.py", *CHECKER_TESTS)  # CUDA equivalence first


def projection(a, cfg) -> list[str]:
    """project_r1's verdict on this config (R1 configs only), with the dollars it implies at the card's ceiling."""
    from experiments.phase2 import project_r1 as P

    if a.model not in P.R1 or P.R1[a.model][0] != Path(cfg).name:
        print(f"    no projection for {Path(cfg).name} (not an R1 config)")
        return []
    arms = [p for p in sorted(ROOT.glob("runs/*/r0-*/manifest.json"))
            if json.loads(p.read_text())["config"]["model"] == a.model]  # this model's R0 timing arm(s)
    rates = sorted(ROOT.glob(f"runs/*/{HF_RATE.get(a.model, '-')}/hf_rate.json"))
    if not arms or not rates:
        print(f"    no R0 measurement of {a.model} yet: R0 measures it on the box; project_r1 is the stop rule before R1")
        return []
    argv = ["--model", a.model, *map(str, arms), "--hf_rate", str(rates[-1])]
    import contextlib
    import io

    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = P.main(argv)
    doc = json.loads(out.getvalue().splitlines()[-1]) if code != 2 else {}
    if code == 2:
        return ["project_r1 could not project (see its message above)"]
    dph = ceiling(a.gpu)
    print(f"    projected {doc['projected_h']:.2f} h (stop {doc['stop_h']:.2f} h) = ${doc['projected_h'] * dph:.2f}"
          f" at ${dph:.2f}/h")  # fmt: skip
    return [] if doc["go"] else [out.getvalue().splitlines()[0]]


CHECKS = [
    ("DATA bundle integrity", view_integrity), ("DATA patient disjointness", patients),
    ("DATA prompt identity", prompts), ("DATA CT window", ct_window), ("DATA training targets legal", targets),
    ("DATA coverage of the exam views", coverage), ("CHECKER CPU exactness (+ CUDA equivalence)", checker_tests),
    ("CODE TRL private names", lambda a, cfg: trl_problems()),
    ("CODE quick tests", lambda a, cfg: pytest("-m", "not slow", *QUICK_TESTS)),
    ("CODE CPU dry run train -> select", lambda a, cfg: pytest(DRY_RUN)),
    ("CODE projected GPU hours and dollars", projection),
]  # fmt: skip


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--gpu", default="RTX_5090", help="card of the rental (policy.CARDS): sets the $/h ceiling")
    ns = ap.parse_args(argv)
    cfg = ns.config
    os.chdir(ROOT)
    a, _ = load(cfg)
    a.gpu = ns.gpu
    failed = ""
    stamp = {"config": cfg.as_posix(), "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "hashes": hashes(cfg),
             "checks": {}, "pass": True}  # fmt: skip
    for name, fn in CHECKS:
        if not stamp["pass"] and name.split()[0] != failed:  # fail fast by stage: all DATA problems show together
            break
        t0 = time.time()
        print(f"[{len(stamp['checks']) + 1}/{len(CHECKS)}] {name}")
        try:
            bad = fn(a, cfg)
        except Exception as e:  # a check that crashes is a failed check, never a pass
            bad = [f"check crashed: {type(e).__name__}: {e}"]
        stamp["checks"][name] = {"pass": not bad, "seconds": round(time.time() - t0, 1), "problems": bad[:40]}
        if bad:
            stamp["pass"], failed = False, name.split()[0]
            print(f"    FAIL ({len(bad)})\n" + "\n".join(f"    - {b}" for b in bad[:40]))
        else:
            print("    ok")
    path = stamp_path(cfg)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(stamp, indent=1))
    print(f"PREFLIGHT {'PASSED' if stamp['pass'] else 'FAILED'}: {cfg.name} -> {path}")
    return 0 if stamp["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
