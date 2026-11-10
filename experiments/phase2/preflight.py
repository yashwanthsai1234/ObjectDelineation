"""Preflight: everything that can be proven on the Mac, before any GPU is rented (the launch gate reads its stamp).

    python -I -m experiments.phase2.preflight --config experiments/phase2/configs/r1_qwen3vl2b.yaml     (make preflight CONFIG=...)

Checks, in order; the stage (DATA / CHECKER / CODE) with a failure is finished, the later stages are not run:
  DATA     bundle integrity (rows.jsonl and context hashes, eval-view file hashes), patients shared between training and
           any exam (or probe) view, coverage of every exam view by the training view (objdel.data.coverage), prompt
           text identical, CT window consistent, every training target parses and is legal for the reference checker;
           with `rows_plan` (B1.3): the plan's train_rl rows root built locally passes the B data gates
           (experiments.phase2.audit_data) with the manifest SHA-256 the plan pins, and the pinned HF revision holds
           exactly the local root's files
  CHECKER  tests/test_cuda_checker*.py when they exist, the scorers on the CUDA checker (tests/test_gpu_scoring.py:
           equal to the host scorer it replaced), then the checker's exactness tests (the slow exact oracles, the
           archived CPU checker's recorded results, topology, geometry)
  CODE     the private TRL names fast.py relies on (objdel.train.trl_adapter), a quick test subset, the CPU dry run
           of train -> select on a tiny model, the projected GPU hours and dollars
Writes runs/preflight/<config>.json: pass per check and the code, data and config hashes. vast launch (policy.py)
refuses unless that stamp passes and all three hashes still match: refusal(). Exit 0 = every check passed.
"""

import argparse
import contextlib
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

from objdel.hashing import file_sha256

ROOT = Path(__file__).resolve().parents[2]
# the exam views policy.chain passes to select_ckpt (gates: shifted, rescaled; report: rotated; probes)
EXAMS = ("lidc-val-shifted", "lidc-val-rescaled", "lidc-val-rotated", "probe-precision", "probe-topology")
CHECKER_TESTS = ["tests/test_checker.py", "tests/test_checker_oracle.py", "tests/test_checker_host.py",
                 "tests/test_topology.py", "tests/test_geometry.py", "tests/test_geometry_mask.py",
                 "tests/test_geometry_oracle.py"]  # fmt: skip
QUICK_TESTS = ["tests/test_data_coverage.py", "tests/test_preflight.py", "tests/test_p2_scripts.py",
               "tests/test_train_configs.py", "tests/test_vast.py", "tests/test_p2_v1_select.py"]  # fmt: skip
DRY_RUN = "tests/test_p2_scripts.py::test_cpu_smoke_train_then_select"  # train -> select on a tiny model, CPU
VASTAI_VERSION = "1.11.0"  # the CLI pin (CLEAN13): an isolated `uv tool`, never in pyproject.toml
HF_RATE = {"lfm2.5-vl-450m": "hf-rate-lfm", "qwen3-vl-2b": "hf-rate-qwen3vl"}  # runs/<run>/<dir>/hf_rate.json (R0)


def ceiling(gpu: str) -> float:
    """$/h price ceiling of the card (policy.CARDS): a bound, not a quote."""
    from experiments.phase2.vast import policy

    return policy.CARDS[gpu]["ceiling"]


def load(cfg: Path):
    """(RunArgs, raw config dict) of a run config."""
    import yaml

    from experiments.phase2 import train_sft as TS

    raw = yaml.safe_load(Path(cfg).read_text())
    return TS.RunArgs(**raw), raw


def hashes(cfg: Path) -> dict:
    """code: as run manifests (train_sft.code_hashes); data: the training view hash, each exam view's manifest, the
    CUDA checker eval set's (it ships in the LIDC bundle) and, with rows_plan, the rows root's manifest.json;
    config: the file's bytes."""
    from experiments.phase2 import train_sft as TS

    a, _ = load(cfg)
    root = Path(a.data_root)
    view = json.loads((root / a.view / "manifest.json").read_text())["view_hash"]
    exams = [file_sha256(root / "eval_views" / v / "manifest.json") for v in EXAMS] + [file_sha256(root / "checker_eval" / "manifest.json")]
    rows = ""
    if a.rows_plan:  # B1.3: an edited rows root refuses the launch
        with contextlib.suppress(OSError, ValueError, KeyError, StopIteration):
            rows = file_sha256(rows_plan(a)[0] / "manifest.json")
    return {"code": TS.code_hashes()["objdel_and_script_sha256"], "config": file_sha256(cfg),
            "data": hashlib.sha256((view + "".join(exams) + rows).encode()).hexdigest()}  # fmt: skip


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
    if not stamp.get("pass"):
        return f"{path} records a FAILED preflight (failed: {[k for k, v in stamp.get('checks', {}).items() if not v['pass']]})"
    if stale := [k for k in now if now[k] != stamp.get("hashes", {}).get(k)]:
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
    if file_sha256(vdir / "rows.jsonl") != man["view_hash"]:
        bad.append(f"{vdir}/rows.jsonl does not match the manifest's view_hash")
    try:
        D._contexts(vdir, man, {r["item_id"] for r in rows})
    except ValueError as e:
        bad.append(str(e))
    for v in EXAMS:
        d = root / "eval_views" / v
        for name, f in json.loads((d / "manifest.json").read_text())["files"].items():
            if not (d / name).is_file() or file_sha256(d / name) != f["sha256"]:
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
    from objdel.data import coverage as C

    out = []
    for r in lines(root / "eval_views" / v / "lidc" / "val.jsonl"):
        out += [p for t in r["targets"] if (p := C.props(t, r["W"], r["view"].get("op", "id")))]
    return out


def coverage(a, cfg) -> list[str]:
    from objdel.data import coverage as C

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


def local_rows_root(train_arg: str, data_root="data") -> Path:
    """The local folder B1.4 built for the box path `{data}/rows/<name>` (fetch_hf.sh --rows): <data_root>/views/<name>;
    ValueError for any other form."""
    m = re.fullmatch(r"\{data\}/rows/([^/]+)", train_arg)
    if not m or m[1] in (".", ".."):
        raise ValueError(f"--train {train_arg!r}: a rows plan names {{data}}/rows/<name>")
    return Path(data_root) / "views" / m[1]


def rows_plan(a) -> tuple[Path, str, str]:
    """(local rows root, --rows-manifest-sha256, --rows-revision) of the train_rl stage of a.rows_plan."""
    import yaml

    step = next(s for s in yaml.safe_load(Path(a.rows_plan).read_text())["steps"] if "experiments.phase2.train_rl" in s)
    train, manifest, revision = (str(step[step.index(f) + 1]) for f in ("--train", "--rows-manifest-sha256",
                                                                         "--rows-revision"))
    return local_rows_root(train, a.data_root), manifest, revision


def rows_root(a, cfg) -> list[str]:
    """B1.3: the local rows root has the manifest the plan pins and passes every hard B data gate."""
    if not a.rows_plan:
        return []
    from experiments.phase2 import audit_data
    from objdel.data import quality as Q

    root, want, _ = rows_plan(a)
    if not (root / "manifest.json").is_file():
        return [f"{root}: no manifest.json (build_rl_rows writes it only after the gates pass)"]
    if (got := file_sha256(root / "manifest.json")) != want:
        return [f"{root}/manifest.json SHA-256 {got} is not the plan's --rows-manifest-sha256 {want}"]
    data = Path(a.data_root)
    report = audit_data.audit(root, data / "store", data / "eval_views", model=a.model)
    return [f"{g.name}: {g.value:g} (limit {g.limit})" for g in Q.failures(Q.gates(report))]


def rows_revision(a, cfg) -> list[str]:
    """B1.3 2b: the pinned HF revision holds exactly the local root's files under rows/<name>/ and its manifest.json
    has the pinned SHA-256 ($0: two Hub API reads, no box)."""
    if not a.rows_plan:
        return []
    import huggingface_hub as hub

    from experiments.phase2.hf_push import REPO
    from objdel import secrets

    root, want, rev = rows_plan(a)
    prefix, token = f"rows/{root.name}/", secrets.get("objdel-hf-write")
    remote = sorted(f[len(prefix):] for f in hub.HfApi(token=token).list_repo_files(REPO, revision=rev)
                    if f.startswith(prefix))
    local = sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file())
    bad = [] if remote == local else [f"{REPO}@{rev} {prefix}: {len(remote)} files, local root {len(local)}; "
                                      f"differ in {sorted(set(local) ^ set(remote))[:3]}"]
    got = file_sha256(hub.hf_hub_download(REPO, prefix + "manifest.json", revision=rev, token=token))
    return bad + ([] if got == want else [f"{REPO}@{rev} {prefix}manifest.json SHA-256 {got} is not {want}"])


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


def vastai_problems() -> list[str]:
    """`vastai --version` must print the pinned token (exactly, not as a prefix): [] or one message. It is a local call,
    never the Vast API."""
    fix = f"install it with: uv tool install vastai=={VASTAI_VERSION}"
    try:
        p = subprocess.run(["vastai", "--version"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as e:
        return [f"vastai {VASTAI_VERSION} is required, but `vastai --version` did not run ({e}); {fix}"]
    if p.returncode == 0 and VASTAI_VERSION in p.stdout.split():
        return []
    return [f"vastai {VASTAI_VERSION} is required, found {p.stdout.strip() or 'no output'!r} "
            f"(exit {p.returncode}); {fix}"]  # fmt: skip


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


def runbooks(a, cfg, root: Path = ROOT) -> list[str]:
    """The kernel drafts run from the synced code (Implementation.md D1 step 0): no scratch path in any of their files,
    and no hand copy or stale box torch in their runbooks."""
    k = root / "docs/implementation2.0/kernels"
    files = [p for p in k.rglob("*") if p.is_file() and "build" not in p.relative_to(k).parts]
    bad = [f"{p.relative_to(root)}: scratch path" for p in files if re.search(r"/[p]rivate/tmp", p.read_text("latin-1"))]
    return bad + [f"{p.relative_to(root)}: {s.strip()!r}" for p in sorted(k.glob("*/GPU_RUN.md"))
                  for s in ("torch 2.14.1", "scp ", "rsync ") if s in p.read_text()]


CHECKS = [
    ("DATA bundle integrity", view_integrity), ("DATA patient disjointness", patients),
    ("DATA prompt identity", prompts), ("DATA CT window", ct_window), ("DATA training targets legal", targets),
    ("DATA coverage of the exam views", coverage), ("DATA ROWS root", rows_root), ("DATA ROWS revision", rows_revision),
    ("CHECKER CPU exactness (+ CUDA equivalence)", checker_tests),
    ("CODE TRL private names", lambda a, cfg: trl_problems()),
    ("CODE vastai version", lambda a, cfg: vastai_problems()),
    ("CODE kernel runbooks", runbooks), ("CODE quick tests", lambda a, cfg: pytest("-m", "not slow", *QUICK_TESTS)),
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
