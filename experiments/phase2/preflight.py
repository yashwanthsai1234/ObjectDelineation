"""Preflight: everything that can be proven on the Mac, before any GPU is rented (the launch gate reads its stamp).

    python -I -m experiments.phase2.preflight --config specs/b1_rows.yaml [--bundle-dir DIR]
    (make preflight CONFIG=specs/... [BUNDLE_DIR=DIR]; DIR as sync_up's --bundle-dir, default build/bundles)

The config is a run spec (objdel.core.spec.load_spec). On the Mac, {data} is the data root ($OBJDEL_DATA_ROOT, else
data; `local`). The training view T of each check is each sft stage's view, else the spec's bundle_view; a spec with
neither (no training) has no T. Checks, in order; the stage (DATA / CHECKER / CODE) with a failure is finished, the
later stages are not run:
  DATA     every {data}/... input of the rl, eval, select and probe stages exists; the upload bundle (DIR/bundle-lidc.tar)
           holds each of them but rows/<name> and a sample of its expert caches has the current reward fingerprint;
           bundle integrity (each T's rows.jsonl
           and context hashes, the eval-view file hashes), patients shared between T and any exam (or probe) view,
           coverage of every exam view by T (objdel.data.coverage), prompt text identical, T's CT window consistent,
           every T target parses and is legal for the reference checker; per rl stage trained on {data}/rows/<name>
           (B1.3; a train_rl run inside a cmd stage's `bash -c` string, the nsys-wrapped arms of the C specs, counts
           as the rl stage its words parse to, rl_stages): the root built locally passes the B data gates
           (experiments.phase2.audit_data) with the manifest
           SHA-256 the stage pins, and the pinned HF revision holds exactly the local root's files
  CHECKER  tests/test_cuda_checker*.py when they exist, the scorers on the CUDA checker (tests/test_gpu_scoring.py:
           equal to the host scorer it replaced), then the checker's exactness tests (the slow exact oracles, the
           archived CPU checker's recorded results, topology, geometry)
  CODE     the private TRL names fast.py relies on (objdel.train.trl_adapter), a quick test subset, the CPU dry run
           of train -> select on a tiny model
Writes runs/preflight/<spec>.json: pass per check and the code, data and config (spec bytes) hashes. vast launch
(policy.py) refuses unless that stamp passes and all three hashes still match: refusal(). Exit 0 = every check passed.
"""

import argparse
import contextlib
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

from objdel.hashing import file_sha256
from objdel.io import read_json

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
BUNDLE_DIR = ROOT / "build" / "bundles"  # sync_up's --bundle-dir default (policy.CLI); --bundle-dir here too
CACHE_SAMPLE = 16  # expert caches whose fingerprint the bundle check reads (of ~12.7k in the LIDC bundle)


@dataclass(frozen=True)
class Inputs:
    """What the checks read: the spec, the Mac data root, the training views T, the rl stages on a rows root and the
    LIDC bundle sync_up uploads."""

    spec: object
    root: Path
    views: tuple[Path, ...]
    rows: tuple
    bundle: Path


def data_root() -> Path:
    """The Mac value of {data}: $OBJDEL_DATA_ROOT, else data."""
    return Path(os.environ.get("OBJDEL_DATA_ROOT") or "data")


def local(template: str) -> Path | None:
    """The Mac folder of a box path: {data}/rows/<name> -> <data root>/views/<name> (what fetch_hf.sh --rows fetches),
    {data}/<rest> -> <data root>/<rest>; None for a run output {rd}/...; ValueError for anything else."""
    if template.startswith("{rd}/"):
        return None
    rest = template.removeprefix("{data}/")
    if rest == template or ".." in rest.split("/") or any(c in rest for c in "{}"):
        raise ValueError(f"{template!r}: not {{data}}/<path> or {{rd}}/<path>")
    if re.fullmatch(r"rows/[^/]+", rest):
        return local_rows_root(template, data_root())
    return data_root() / rest


def inputs(cfg: Path, bundle_dir: Path = BUNDLE_DIR) -> Inputs:
    """The spec's inputs on the Mac: the sft views, the rl rows stages and the data root."""
    from objdel.core.spec import load_spec

    spec, root = load_spec(cfg), data_root()
    views = [Path(s.data_root) / s.view if s.data_root else root / s.view for s in spec.stages if s.kind == "sft"]
    views = views or ([root / spec.bundle_view] if spec.bundle_view else [])
    rows = [s for s in rl_stages(spec) if s.train and re.fullmatch(r"\{data\}/rows/[^/]+", s.train)]
    return Inputs(spec, root, tuple(views), tuple(rows), Path(bundle_dir) / "bundle-lidc.tar")


def spec_inputs(spec) -> list[str]:
    """Every {data}/... string field (and tuple item) of the rl, eval, select and probe stages, in stage order."""
    out = []
    for s in spec.stages:
        if s.kind in ("rl", "eval", "select", "probe"):
            for v in s.model_dump(exclude={"kind", "name"}).values():
                out += [x for x in (v if isinstance(v, list | tuple) else [v])
                        if isinstance(x, str) and x.startswith("{data}")]
    return out


def hashes(cfg: Path) -> dict:
    """code: as run manifests (core.provenance.code_hashes); data: each T's view hash, each exam view's manifest, the CUDA
    checker eval set's (it ships in the LIDC bundle), each rows root's manifest.json and the manifest.json of each
    spec input that has one; config: the spec file's bytes."""
    from objdel.core.provenance import code_hashes

    s = inputs(cfg)
    views = [read_json(t / "manifest.json")["view_hash"] for t in s.views]
    exams = [file_sha256(s.root / "eval_views" / v / "manifest.json") for v in EXAMS]
    exams.append(file_sha256(s.root / "checker_eval" / "manifest.json"))
    manifests = []
    for t in [st.train for st in s.rows] + spec_inputs(s.spec):  # an edited rows root or input refuses the launch
        with contextlib.suppress(OSError, ValueError):
            if (p := local(t)) is not None and (p / "manifest.json").is_file():
                manifests.append(file_sha256(p / "manifest.json"))
    return {"code": code_hashes()["objdel_and_script_sha256"], "config": file_sha256(cfg),
            "data": hashlib.sha256("".join(views + exams + manifests).encode()).hexdigest()}  # fmt: skip


def stamp_path(cfg: Path) -> Path:
    """runs/preflight/<spec>.json under $OBJDEL_RUN_DIR (as train_sft); the stamp records the spec's path."""
    return ROOT / os.environ.get("OBJDEL_RUN_DIR", "runs") / "preflight" / f"{Path(cfg).stem}.json"


def refusal(cfg: Path) -> str | None:
    """Why a launch of this config must be refused (no stamp, a failed one, or code / data / config changed since), else
    None."""
    path = stamp_path(cfg)
    if not path.is_file():
        return f"no preflight stamp: run `make preflight CONFIG={cfg}` first"
    try:  # a stamp of the wrong shape (null, a list, hashes: null) is a refusal, never a traceback in cmd_launch
        stamp, now = read_json(path), hashes(cfg)
        if stamp.get("config") != str(Path(cfg).resolve()):
            return f"{path} was written for another spec path ({stamp.get('config')}): run `make preflight CONFIG={cfg}`"
        if not stamp.get("pass"):
            failed = [k for k, v in stamp.get("checks", {}).items() if not v["pass"]]
            return f"{path} records a FAILED preflight (failed: {failed})"
        if stale := [k for k in now if now[k] != stamp.get("hashes", {}).get(k)]:
            return f"{path} is stale ({', '.join(stale)} changed since it passed): run `make preflight` again"
    except (OSError, ValueError, KeyError, AttributeError, TypeError) as e:
        return f"cannot check {path} against the current code, data and config: {e}"
    return None


# --- DATA ---------------------------------------------------------------------------------------------------------


def lines(path):
    """The JSON records of a JSONL file."""
    return [json.loads(x) for x in Path(path).read_text().splitlines()]


def given_inputs(s: Inputs) -> list[str]:
    """Every {data}/... input of the rl, eval, select and probe stages is a folder on the Mac."""
    bad = []
    for t in spec_inputs(s.spec):
        try:
            if (p := local(t)) is not None and not p.is_dir():
                bad.append(f"{t}: no folder {p}")
        except ValueError as e:
            bad.append(str(e))
    return bad


def sampled(names: list[str]) -> list[str]:
    """About CACHE_SAMPLE of the names, evenly spread over their sorted order (the first always)."""
    return sorted(names)[::max(1, len(names) // CACHE_SAMPLE)]


def bundle(s: Inputs) -> list[str]:
    """The LIDC bundle sync_up uploads holds every {data} input of the spec (but rows/<name>, which fetch_hf.sh
    brings), and a sample of its expert caches (caches/**.npz) has the current reward fingerprint. Reads the tar
    headers and the sampled members only. Box C, 2026-11-10: a Nov 5 bundle held b1eb9f51 caches and no
    lidc-val-edge, and final_eval died with StaleCacheError after training."""
    import tarfile

    import numpy as np

    from objdel.verify.reward.spec import fingerprint

    if not s.bundle.is_file():
        return [f"no bundle {s.bundle}: build it (vast/bundle.py) or pass --bundle-dir"]
    with tarfile.open(s.bundle, "r:") as tar:
        members = {m.name: m for m in tar.getmembers()}
        bad = [f"{s.bundle}: no {rest}/ (a {{data}} input of the spec)" for t in dict.fromkeys(spec_inputs(s.spec))
               if not re.fullmatch(r"\{data\}/rows/[^/]+", t) and not any(n.startswith(
                   (rest := t.removeprefix("{data}/")) + "/") for n in members)]  # fmt: skip
        want = fingerprint()
        for name in sampled([n for n in members if "caches/" in n and n.endswith(".npz")]):
            with np.load(tar.extractfile(members[name]), allow_pickle=False) as z:
                if (fp := str(z["fingerprint"])) != want:
                    bad.append(f"{s.bundle}: {name} built under reward fingerprint {fp[:12]}, not {want[:12]}: "
                               "rebuild the bundle (vast/bundle.py)")
    return bad


def view_integrity(s: Inputs) -> list[str]:
    """Each training view's rows.jsonl matches its view_hash and its contexts load; each exam view's files match its
    manifest."""
    from objdel.train import dataset as D

    bad = []
    for vdir in s.views:
        man, rows = read_json(vdir / "manifest.json"), lines(vdir / "rows.jsonl")
        if file_sha256(vdir / "rows.jsonl") != man["view_hash"]:
            bad.append(f"{vdir}/rows.jsonl does not match the manifest's view_hash")
        try:
            D._contexts(vdir, man, {r["item_id"] for r in rows})
        except ValueError as e:
            bad.append(str(e))
    for v in EXAMS:
        d = s.root / "eval_views" / v
        for name, f in read_json(d / "manifest.json")["files"].items():
            if not (d / name).is_file() or file_sha256(d / name) != f["sha256"]:
                bad.append(f"{d / name} does not match its manifest")
    return bad


def patients(s: Inputs) -> list[str]:
    """No patient (group_id) and no CT series in both the training rows and an exam view; train and val prepared splits
    disjoint too."""
    from objdel.train import dataset as D

    if not s.views:
        return []
    group = {r["item_id"]: r["group_id"] for x in ("train", "val") for r in D.records(s.root / "prepared", "lidc", x)}
    val, bad = {r["group_id"] for r in D.records(s.root / "prepared", "lidc", "val")}, []
    for t in s.views:
        train_items = {r["item_id"] for r in lines(t / "rows.jsonl")}
        train_groups, series = {group.get(i) for i in train_items}, {i.split("/")[1] for i in train_items}
        if None in train_groups:
            bad.append(f"{sum(i not in group for i in train_items)} training items have no patient")
        if shared := train_groups & val:
            bad.append(f"{len(shared)} patients in training and prepared val, e.g. {sorted(shared)[:3]}")
        for v in EXAMS:
            rows = lines(s.root / "eval_views" / v / "lidc" / "val.jsonl")
            shared = train_groups & {r["group_id"] for r in rows}
            if shared or series & {r["item_id"].split("/")[1] for r in rows}:
                bad.append(f"{v}: {len(shared)} patients also in training, e.g. {sorted(shared)[:3]}")
    return bad


def exam_props(root: Path, v: str) -> list[dict]:
    """The coverage properties of each record of exam view `v`."""
    from objdel.data import coverage as C

    out = []
    for r in lines(root / "eval_views" / v / "lidc" / "val.jsonl"):
        out += [p for t in r["targets"] if (p := C.props(t, r["W"], r["view"].get("op", "id")))]
    return out


def coverage(s: Inputs) -> list[str]:
    """Each training view covers the exam views' properties (objdel.data.coverage)."""
    from objdel.data import coverage as C

    bad = []
    for t in s.views:
        train = [p for r in lines(t / "rows.jsonl") if r["epoch"] == 0
                 and (p := C.props(r["target"], r["side"], r["op"]))]  # fmt: skip
        td = C.distribution(train)
        for v in EXAMS:
            for g in C.gaps(td, C.distribution(exam_props(s.root, v))):
                bad.append(f"{v}: {g['key']}={g['category']} is {g['exam']:.1%} of the exam, {g['train']:.2%} of "
                           f"training (floor {C.TRAIN_MAX:.0%}, exam share >= {C.EVAL_MIN:.0%})")  # fmt: skip
    return bad


def prompts(s: Inputs) -> list[str]:
    """Each exam view's prompt text equals the training prompt."""
    from objdel.train import dataset as D

    want = next(c["text"] for c in D._user("lidc")["content"] if c.get("type") == "text")
    return [f"{v}: exam prompt differs from the training prompt"
            for v in EXAMS if read_json(s.root / "eval_views" / v / "manifest.json")
            .get("prompts", {"lidc": want})["lidc"] != want]  # fmt: skip  (probe manifests carry no prompt)


def ct_window(s: Inputs) -> list[str]:
    """Every training window is a proper (lo < hi) window, and the exams' fixed lung window lies in the 1-99 percentile
    range the training jitter covers."""
    import numpy as np

    from objdel.data.targets.frame import LUNG_WINDOW

    bad = []
    for t in s.views:
        w = np.array([r["window"] for r in lines(t / "rows.jsonl")])
        bad += [f"{t}: training rows with lo >= hi"] if (w[:, 0] >= w[:, 1]).any() else []
        for i, name in enumerate(("lo", "hi")):
            lo, hi = np.percentile(w[:, i], [1, 99])
            if not lo <= LUNG_WINDOW[i] <= hi:
                bad.append(f"{t}: exam window {name} {LUNG_WINDOW[i]} outside training 1-99 % range [{lo:.0f}, {hi:.0f}]")
    return bad


def targets(s: Inputs) -> list[str]:
    """Each training target parses and rasterizes to a legal mask (the reference checker)."""
    from objdel.geom.shape import to_mask
    from objdel.geom.text import parse

    bad = []
    for t in s.views:
        for r in lines(t / "rows.jsonl"):
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


def train_rl_argv(step) -> list[str] | None:
    """The words after `-m experiments.phase2.train_rl` of an argv: a stage of its own, or one inside a `bash -c`
    string (the nsys-wrapped arms of specs/c_speed.yaml); None for any other argv."""
    body = list(step[1:] if step[0] == "~" else step)
    words = shlex.split(body[2]) if body[:2] == ["bash", "-c"] and "experiments.phase2.train_rl" in body[2] else body
    i = next((i for i in range(len(words) - 1) if words[i:i + 2] == ["-m", "experiments.phase2.train_rl"]), None)
    return None if i is None else [str(w) for w in words[i + 2:]]


def rl_stages(spec) -> list:
    """Each rl stage, and each train_rl run inside a cmd stage (the nsys-wrapped arms of c_speed.yaml) as the RlStage
    its words parse to (named as the cmd stage), in stage order. Unparsable words raise ValueError."""
    from objdel.core import spec as S

    out = []
    for st in spec.stages:
        words = train_rl_argv([S.fill(w, "{rd}", "{data}") for w in st.argv]) if st.kind == "cmd" else None
        if st.kind == "rl":
            out.append(st)
        elif words is not None:
            try:
                out.append(S.cli_parse(S.RlStage, words).model_copy(update={"name": st.name}))
            except SystemExit as e:
                raise ValueError(f"stage {st.name}: its train_rl words do not parse as an rl stage") from e
    return out


def rows_root(s: Inputs) -> list[str]:
    """B1.3: each rl stage's local rows root has the manifest the stage pins and passes every hard B data gate."""
    from experiments.phase2 import audit_data
    from objdel.data import quality as Q

    bad = []
    for st in s.rows:
        root = local(st.train)
        if st.rows_manifest_sha256 is None:
            bad.append(f"stage {st.name}: a rows root needs rows_manifest_sha256")
        elif not (root / "manifest.json").is_file():
            bad.append(f"{root}: no manifest.json (build_rl_rows writes it only after the gates pass)")
        elif (got := file_sha256(root / "manifest.json")) != st.rows_manifest_sha256:
            bad.append(f"{root}/manifest.json SHA-256 {got} is not stage {st.name}'s {st.rows_manifest_sha256}")
        else:
            report = audit_data.audit(root, s.root / "store", s.root / "eval_views", model=st.model)
            bad += [f"{g.name}: {g.value:g} (limit {g.limit})" for g in Q.failures(Q.gates(report))]
    return bad


def rows_revision(s: Inputs) -> list[str]:
    """B1.3 2b: each pinned HF revision holds exactly the local root's files under rows/<name>/ and its manifest.json
    has the pinned SHA-256 ($0: two Hub API reads per root, no box)."""
    if not s.rows:
        return []
    import huggingface_hub as hub

    from objdel.infra import secrets
    from objdel.infra.hf import REPO

    bad, token = [], secrets.get("objdel-hf-write")
    for st in s.rows:
        root, want, rev = local(st.train), st.rows_manifest_sha256, st.rows_revision
        prefix = f"rows/{root.name}/"
        remote = sorted(f[len(prefix):] for f in hub.HfApi(token=token).list_repo_files(REPO, revision=rev)
                        if f.startswith(prefix))
        here = sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file())
        if remote != here:
            bad.append(f"{REPO}@{rev} {prefix}: {len(remote)} files, local root {len(here)}; "
                       f"differ in {sorted(set(here) ^ set(remote))[:3]}")
        got = file_sha256(hub.hf_hub_download(REPO, prefix + "manifest.json", revision=rev, token=token))
        bad += [] if got == want else [f"{REPO}@{rev} {prefix}manifest.json SHA-256 {got} is not {want}"]
    return bad


# --- CHECKER, CODE ------------------------------------------------------------------------------------------------


def pytest(*args: str) -> list[str]:
    """Run pytest with `args` in the repo; [] on a pass, else the command and the last 15 output lines."""
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


def checker_tests(s: Inputs) -> list[str]:
    """The CUDA checker tests, the GPU scoring test and CHECKER_TESTS (CUDA equivalence first)."""
    cuda = sorted(p.relative_to(ROOT).as_posix() for p in (ROOT / "tests").glob("test_cuda_checker*.py"))
    return pytest("-m", "not slow", *cuda, "tests/test_gpu_scoring.py", *CHECKER_TESTS)  # CUDA equivalence first


def runbooks(s: Inputs, root: Path = ROOT) -> list[str]:
    """The kernel drafts run from the synced code (Implementation.md D1 step 0): no scratch path in any of their files,
    and no hand copy or stale box torch in their runbooks."""
    k = root / "docs/implementation2.0/kernels"
    files = [p for p in k.rglob("*") if p.is_file() and "build" not in p.relative_to(k).parts]
    bad = [f"{p.relative_to(root)}: scratch path" for p in files if re.search(r"/[p]rivate/tmp",
        p.read_text("latin-1"))]
    return bad + [f"{p.relative_to(root)}: {s.strip()!r}" for p in sorted(k.glob("*/GPU_RUN.md"))
                  for s in ("torch 2.14.1", "scp ", "rsync ") if s in p.read_text()]


CHECKS = [
    ("DATA spec inputs", given_inputs), ("DATA upload bundle", bundle),
    ("DATA bundle integrity", view_integrity), ("DATA patient disjointness", patients),
    ("DATA prompt identity", prompts), ("DATA CT window", ct_window), ("DATA training targets legal", targets),
    ("DATA coverage of the exam views", coverage), ("DATA ROWS root", rows_root), ("DATA ROWS revision", rows_revision),
    ("CHECKER CPU exactness (+ CUDA equivalence)", checker_tests),
    ("CODE TRL private names", lambda s: trl_problems()),
    ("CODE vastai version", lambda s: vastai_problems()),
    ("CODE kernel runbooks", runbooks), ("CODE quick tests", lambda s: pytest("-m", "not slow", *QUICK_TESTS)),
    ("CODE CPU dry run train -> select", lambda s: pytest(DRY_RUN)),
]  # fmt: skip


def main(argv=None) -> int:
    """The preflight CLI: run every check on --config and write the stamp only when all pass."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", type=Path, required=True, help="a run spec (specs/<name>.yaml)")
    ap.add_argument("--bundle-dir", type=Path, default=BUNDLE_DIR,
                    help="the folder of bundle-lidc.tar, as sync_up's --bundle-dir (default build/bundles)")
    a = ap.parse_args(argv)
    cfg, bundle_dir = a.config.resolve(), a.bundle_dir.resolve()
    os.chdir(ROOT)
    s = inputs(cfg, bundle_dir)
    failed = ""
    stamp = {"config": str(cfg.resolve()), "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "hashes": hashes(cfg),
             "checks": {}, "pass": True}  # fmt: skip
    for name, fn in CHECKS:
        if not stamp["pass"] and name.split()[0] != failed:  # fail fast by stage: all DATA problems show together
            break
        t0 = time.time()
        print(f"[{len(stamp['checks']) + 1}/{len(CHECKS)}] {name}")
        try:
            bad = fn(s)
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
