"""Preflight: everything that can be proven on the Mac, before any GPU is rented (the launch gate reads its stamp).

    python -I -m objdel.cli.box.preflight --config specs/b1_rows.yaml [--bundle-dir DIR]
    (make preflight CONFIG=specs/... [BUNDLE_DIR=DIR]; DIR as sync_up's --bundle-dir, default build/bundles)

That script parses the flags and calls `run` here.

The config is a run spec (objdel.core.spec.load_spec). On the Mac, {data} is the data root ($OBJDEL_DATA_ROOT, else
data; `local`). The training view T of each check is each sft stage's view, else the spec's bundle_view; a spec with
neither (no training) has no T. Checks, in order; the stage (DATA / CHECKER / CODE) with a failure is finished, the
later stages are not run:
  DATA     every {data}/... input of the rl, eval, select and probe stages exists; the upload bundle
           (DIR/bundle-<dataset>.tar) holds each of them but rows/<name> and a sample of its expert caches has the
           current reward fingerprint; bundle integrity (each T's rows.jsonl and context hashes, the eval-view file
           hashes), patients shared between T and any exam (or probe) view,
           coverage of every exam view by T (objdel.data.coverage), prompt text identical, T's CT window consistent,
           every T target parses and is legal for the reference checker; per rl stage trained on {data}/rows/<name>
           (a train_rl run inside a cmd stage's `bash -c` string, the nsys-wrapped arms of the C specs, counts
           as the rl stage its words parse to, rl_stages): the root built locally passes the B data gates
           (objdel.data.quality.audit) with the manifest
           SHA-256 the stage pins, and the pinned HF revision holds exactly the local root's files
  CHECKER  tests/geom/test_cuda_checker*.py when they exist, the scorers on the CUDA checker
           (tests/verify/test_gpu_scoring.py:
           equal to the host scorer it replaced), then the checker's exactness tests (the slow exact oracles, the
           archived CPU checker's recorded results, topology, geometry)
  CODE     the private TRL names the adapter relies on (objdel.rl.trl_adapter), a quick test subset, the CPU dry run
           of train -> select on a tiny model
Writes runs/preflight/<spec>.json: pass per check and the code, data and config (spec bytes) hashes. vast launch
(objdel.infra.vast) refuses unless that stamp passes and all three hashes still match: refusal(). Exit 0 = every check
passed.
"""

import contextlib
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import time
import types
from dataclasses import dataclass
from pathlib import Path

from objdel.core.hashing import file_sha256
from objdel.core.io import read_json
from objdel.core.provenance import repo_root

CHECKER_TESTS = ["tests/geom/test_checker.py", "tests/geom/test_checker_oracle.py", "tests/geom/test_checker_host.py",
                 "tests/geom/test_topology.py", "tests/geom/test_geometry.py", "tests/geom/test_geometry_mask.py",
                 "tests/geom/test_geometry_oracle.py"]  # fmt: skip
QUICK_TESTS = ["tests/data/test_data_coverage.py", "tests/experiments/test_preflight.py",
               "tests/experiments/test_p2_scripts.py", "tests/policy/test_train_configs.py", "tests/infra/test_vast.py",
               "tests/experiments/test_p2_v1_select.py"]  # fmt: skip
DRY_RUN = "tests/experiments/test_p2_scripts.py::test_cpu_smoke_train_then_select"  # train -> select, tiny model, CPU
VASTAI_VERSION = "1.11.0"  # the CLI pin: an isolated `uv tool`, never in pyproject.toml
BUNDLE_DIR = repo_root() / "build" / "bundles"  # sync_up's --bundle-dir default (objdel.infra.vast CLI); --bundle-dir
CACHE_SAMPLE = 16  # expert caches whose fingerprint the bundle check reads (of ~12.7k in the LIDC bundle)


@dataclass(frozen=True)
class Inputs:
    """What the checks read: the spec, its dataset entry (objdel.core.datasets: the exam views, prompt, window and
    structure names), the Mac data root, the training views T, the rl stages on a rows root and the bundle sync_up
    uploads."""

    spec: object
    dataset: object
    root: Path
    views: tuple[Path, ...]
    rows: tuple
    bundle: Path


def rows_name(template: str) -> str | None:
    """<name> of a `{data}/rows/<name>` box path, else None."""
    m = re.fullmatch(r"\{data\}/rows/([^/]+)", template)
    return m[1] if m and m[1] not in (".", "..") else None


def data_root() -> Path:
    """The Mac value of {data}: $OBJDEL_DATA_ROOT, else data."""
    return Path(os.environ.get("OBJDEL_DATA_ROOT") or "data")


def local(template: str) -> Path | None:
    """The Mac folder of a box path: {data}/rows/<name> -> <data root>/views/<name> (what fetch_hf.sh --rows fetches),
    {data}/<rest> -> <data root>/<rest>; None for a run output {rd}/...; ValueError for anything else."""
    if template.startswith("{rd}/"):
        return None
    rest = template.removeprefix("{data}/")
    if rest == template or {".", ".."} & set(rest.split("/")) or any(c in rest for c in "{}"):
        raise ValueError(f"{template!r}: not {{data}}/<path> or {{rd}}/<path>")
    if name := rows_name(template):
        return data_root() / "views" / name
    return data_root() / rest


def inputs(cfg: Path, bundle_dir: Path = BUNDLE_DIR) -> Inputs:
    """The spec's inputs on the Mac: the sft views, the rl rows stages and the data root."""
    from objdel.core import datasets
    from objdel.core.spec import load_spec, spec_dataset

    spec, root = load_spec(cfg), data_root()
    ds = datasets.get(spec_dataset(spec))
    views = [Path(s.data_root) / s.view if s.data_root else root / s.view for s in spec.stages if s.kind == "sft"]
    views = views or ([root / spec.bundle_view] if spec.bundle_view else [])
    rows = [s for s in rl_stages(spec) if s.train and rows_name(s.train)]
    return Inputs(spec, ds, root, tuple(views), tuple(rows), Path(bundle_dir) / f"bundle-{ds.name}.tar")


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
    """code: as run manifests (core.provenance.code_hashes); data: each T's view hash, each exam view's manifest, the
    CUDA checker eval set's (it ships in the LIDC bundle), each rows root's manifest.json and the manifest.json of each
    spec input that has one; config: the spec file's bytes."""
    from objdel.core.provenance import code_hashes

    s = inputs(cfg)
    views = [read_json(t / "manifest.json")["view_hash"] for t in s.views]
    exams = [file_sha256(s.root / "eval_views" / v / "manifest.json") for v in s.dataset.exams]
    exams += [file_sha256(s.root / "checker_eval" / "manifest.json")] if s.dataset.checker_eval else []
    manifests = []
    for t in [st.train for st in s.rows] + spec_inputs(s.spec):  # an edited rows root or input refuses the launch
        with contextlib.suppress(OSError, ValueError):
            if (p := local(t)) is not None and (p / "manifest.json").is_file():
                manifests.append(file_sha256(p / "manifest.json"))
    return {"code": code_hashes()["objdel_and_script_sha256"], "config": file_sha256(cfg),
            "data": hashlib.sha256("".join(views + exams + manifests).encode()).hexdigest()}  # fmt: skip


def stamp_path(cfg: Path) -> Path:
    """runs/preflight/<spec>.json under $OBJDEL_RUN_DIR (as train_sft); the stamp records the spec's path."""
    return repo_root() / os.environ.get("OBJDEL_RUN_DIR", "runs") / "preflight" / f"{Path(cfg).stem}.json"


def refusal(cfg: Path) -> str | None:
    """Why a launch of this config must be refused (no stamp, a failed one, or code / data / config changed since), else
    None."""
    path = stamp_path(cfg)
    if not path.is_file():
        return f"no preflight stamp: run `make preflight CONFIG={cfg}` first"
    try:  # a stamp of the wrong shape (null, a list, hashes: null) is a refusal, never a traceback in cmd_launch
        stamp, now = read_json(path), hashes(cfg)
        if stamp.get("config") != str(Path(cfg).resolve()):
            return (f"{path} was written for another spec path ({stamp.get('config')}): "
                    f"run `make preflight CONFIG={cfg}`")
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
    """The bundle sync_up uploads holds every {data} input of the spec (but rows/<name>, which fetch_hf.sh
    brings), and a sample of its expert caches (caches/**.npz) has the current reward fingerprint. Reads the tar
    headers and the sampled members only. A stale bundle (old reward fingerprint, a missing eval view) fails here,
    not after training."""
    import tarfile

    import numpy as np

    from objdel.verify.reward.spec import fingerprint

    if not s.bundle.is_file():
        return [f"no bundle {s.bundle}: build it (vast/bundle.py) or pass --bundle-dir"]
    with tarfile.open(s.bundle, "r:") as tar:
        members = {m.name: m for m in tar.getmembers()}
        bad = [f"{s.bundle}: no {rest}/ (a {{data}} input of the spec)" for t in dict.fromkeys(spec_inputs(s.spec))
               if not rows_name(t) and not any(n.startswith(
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
    from objdel.sft import rows as sft_rows

    bad = []
    for vdir in s.views:
        man, rows = read_json(vdir / "manifest.json"), lines(vdir / "rows.jsonl")
        if file_sha256(vdir / "rows.jsonl") != man["view_hash"]:
            bad.append(f"{vdir}/rows.jsonl does not match the manifest's view_hash")
        try:
            sft_rows._contexts(vdir, man, {r["item_id"] for r in rows})
        except ValueError as e:
            bad.append(str(e))
    for v in s.dataset.exams:
        d = s.root / "eval_views" / v
        for name, f in read_json(d / "manifest.json")["files"].items():
            if not (d / name).is_file() or file_sha256(d / name) != f["sha256"]:
                bad.append(f"{d / name} does not match its manifest")
    return bad


def patients(s: Inputs) -> list[str]:
    """No patient (group_id) and no CT series in both the training rows and an exam view; train and val prepared splits
    disjoint too."""
    from objdel.policy import rows as policy_rows

    if not s.views:
        return []
    ds = s.dataset.name
    group = {r["item_id"]: r["group_id"] for x in ("train", "val") for r in policy_rows.records(s.root / "prepared",
        ds, x)}
    val, bad = {r["group_id"] for r in policy_rows.records(s.root / "prepared", ds, "val")}, []
    for t in s.views:
        train_items = {r["item_id"] for r in lines(t / "rows.jsonl")}
        train_groups, series = {group.get(i) for i in train_items}, {i.split("/")[1] for i in train_items}
        if None in train_groups:
            bad.append(f"{sum(i not in group for i in train_items)} training items have no patient")
        if shared := train_groups & val:
            bad.append(f"{len(shared)} patients in training and prepared val, e.g. {sorted(shared)[:3]}")
        for v in s.dataset.exams:
            rows = lines(s.root / "eval_views" / v / ds / "val.jsonl")
            shared = train_groups & {r["group_id"] for r in rows}
            if shared or series & {r["item_id"].split("/")[1] for r in rows}:
                bad.append(f"{v}: {len(shared)} patients also in training, e.g. {sorted(shared)[:3]}")
    return bad


def exam_props(root: Path, v: str, dataset: str) -> list[dict]:
    """The coverage properties of each record of exam view `v`."""
    from objdel.data import coverage as C

    out = []
    for r in lines(root / "eval_views" / v / dataset / "val.jsonl"):
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
        for v in s.dataset.exams:
            for g in C.gaps(td, C.distribution(exam_props(s.root, v, s.dataset.name))):
                bad.append(f"{v}: {g['key']}={g['category']} is {g['exam']:.1%} of the exam, {g['train']:.2%} of "
                           f"training (floor {C.TRAIN_MAX:.0%}, exam share >= {C.EVAL_MIN:.0%})")  # fmt: skip
    return bad


def prompts(s: Inputs) -> list[str]:
    """Each exam view's prompt text equals the training prompt (the dataset entry's)."""
    ds, want = s.dataset.name, s.dataset.prompt
    return [f"{v}: exam prompt differs from the training prompt"
            for v in s.dataset.exams if read_json(s.root / "eval_views" / v / "manifest.json")
            .get("prompts", {ds: want})[ds] != want]  # fmt: skip  (probe manifests carry no prompt)


def ct_window(s: Inputs) -> list[str]:
    """Every training window is a proper (lo < hi) window, and the exams' fixed window (the dataset entry's) lies in the
    1-99 percentile range the training jitter covers. A dataset without a window (uint8 images) has nothing to check."""
    import numpy as np

    if (window := s.dataset.window) is None:
        return []
    bad = []
    for t in s.views:
        w = np.array([r["window"] for r in lines(t / "rows.jsonl")])
        bad += [f"{t}: training rows with lo >= hi"] if (w[:, 0] >= w[:, 1]).any() else []
        for i, name in enumerate(("lo", "hi")):
            lo, hi = np.percentile(w[:, i], [1, 99])
            if not lo <= window[i] <= hi:
                bad.append(f"{t}: exam window {name} {window[i]} outside training 1-99 % range "
                           f"[{lo:.0f}, {hi:.0f}]")
    return bad


def targets(s: Inputs) -> list[str]:
    """Each training target parses and rasterizes to a legal mask (the reference checker)."""
    from objdel.geom.shape import to_mask
    from objdel.geom.text import parse

    bad, names = [], list(s.dataset.names)
    for t in s.views:
        for r in lines(t / "rows.jsonl"):
            p = parse(r["target"], names)
            if not p.ok or not all(to_mask(p.objs[n], r["side"], r["side"]).legal for n in names):
                bad.append(f"{r['row_id']}: {p.error or 'illegal'}")
    return [f"{len(bad)} training targets rejected by the reference checker, e.g. {bad[:3]}"] if bad else []


def train_rl_argv(step) -> list[str] | None:
    """The words after `-m <the rl module>` (objdel.core.spec.MODULE) of an argv: a stage of its own, or one inside a
    `bash -c` string (the nsys-wrapped arms of specs/c_speed.yaml); None for any other argv."""
    from objdel.core.spec import MODULE

    body = list(step[1:] if step[0] == "~" else step)
    words = shlex.split(body[2]) if body[:2] == ["bash", "-c"] and MODULE["rl"] in body[2] else body
    i = next((i for i in range(len(words) - 1) if words[i:i + 2] == ["-m", MODULE["rl"]]), None)
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
    """Each rl stage's local rows root has the manifest the stage pins and passes every hard B data gate."""
    from objdel.data import quality as Q
    from objdel.policy import models

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
            report = Q.audit(root, *models.tokenizer(st.model), s.root / "store", s.root / "eval_views")
            bad += [f"{g.name}: {g.value:g} (limit {g.limit})" for g in Q.failures(Q.gates(report))]
    return bad


def rows_revision(s: Inputs) -> list[str]:
    """Each pinned HF revision holds rows/<name>/ = hf.ROWS_TAR and manifest.json only, the tar's LFS SHA-256
    is that of hf.rows_tar over the local root (hashed as a stream, no copy on disk), and the manifest has the pinned
    SHA-256 ($0: three Hub API reads per root, no box)."""
    if not s.rows:
        return []
    import huggingface_hub as hub

    from objdel.infra import secrets
    from objdel.infra.hf import REPO, ROWS_TAR, rows_tar

    bad, token = [], secrets.get("objdel-hf-write")
    for st in s.rows:
        root, want, rev = local(st.train), st.rows_manifest_sha256, st.rows_revision
        prefix, api = f"rows/{root.name}/", hub.HfApi(token=token)
        remote = sorted(f[len(prefix):] for f in api.list_repo_files(REPO, revision=rev) if f.startswith(prefix))
        if remote != sorted([ROWS_TAR, "manifest.json"]):
            bad.append(f"{REPO}@{rev} {prefix}: holds {remote[:3]} ({len(remote)} files), not {ROWS_TAR} and "
                       "manifest.json")
            continue
        lfs = api.get_paths_info(REPO, [prefix + ROWS_TAR], revision=rev)[0].lfs
        h = hashlib.sha256()
        rows_tar(root, types.SimpleNamespace(write=lambda b, h=h: h.update(b) or len(b)))
        if lfs is None or lfs.sha256 != h.hexdigest():
            bad.append(f"{REPO}@{rev} {prefix}{ROWS_TAR}: SHA-256 {getattr(lfs, 'sha256', None)} is not the local "
                       f"root's {h.hexdigest()}")
        got = file_sha256(hub.hf_hub_download(REPO, prefix + "manifest.json", revision=rev, token=token))
        bad += [] if got == want else [f"{REPO}@{rev} {prefix}manifest.json SHA-256 {got} is not {want}"]
    return bad


# --- CHECKER, CODE ------------------------------------------------------------------------------------------------


def pytest(*args: str) -> list[str]:
    """Run pytest with `args` in the repo; [] on a pass, else the command and the last 15 output lines."""
    p = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", *args], cwd=repo_root(),
                       capture_output=True, text=True)
    return [] if p.returncode == 0 else ["pytest " + " ".join(args)] + p.stdout.splitlines()[-15:]


def trl_problems() -> list[str]:
    """The private TRL names the adapter relies on (objdel.rl.trl_adapter): [] or the contract error's message."""
    from objdel.rl import trl_adapter

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
    root = repo_root()
    cuda = sorted(p.relative_to(root).as_posix() for p in (root / "tests/geom").glob("test_cuda_checker*.py"))
    return pytest("-m", "not slow", *cuda, "tests/verify/test_gpu_scoring.py", *CHECKER_TESTS)  # CUDA equivalence first


CHECKS = [
    ("DATA spec inputs", given_inputs), ("DATA upload bundle", bundle),
    ("DATA bundle integrity", view_integrity), ("DATA patient disjointness", patients),
    ("DATA prompt identity", prompts), ("DATA CT window", ct_window), ("DATA training targets legal", targets),
    ("DATA coverage of the exam views", coverage), ("DATA ROWS root", rows_root), ("DATA ROWS revision", rows_revision),
    ("CHECKER CPU exactness (+ CUDA equivalence)", checker_tests),
    ("CODE TRL private names", lambda s: trl_problems()),
    ("CODE vastai version", lambda s: vastai_problems()),
    ("CODE quick tests", lambda s: pytest("-m", "not slow", *QUICK_TESTS)),
    ("CODE CPU dry run train -> select", lambda s: pytest(DRY_RUN)),
]  # fmt: skip


def run(cfg: Path, bundle_dir: Path = BUNDLE_DIR) -> int:
    """Run every check on the spec `cfg` (from the repo root) with the upload bundle in `bundle_dir` and write the
    stamp; 0 when all pass, else 1."""
    cfg, bundle_dir = cfg.resolve(), Path(bundle_dir).resolve()
    os.chdir(repo_root())
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
