"""Data bundles for the Vast box (PLAN §4.4, D16): one deterministic tar per dataset with a SHA256SUMS manifest.

    .venv/bin/python -I -m experiments.phase2.vast.bundle lidc|riga [--view NAME ...] [--data data] [--out build/bundles]

Members (sorted, mtime 0, uid/gid 0, mode 0644; SHA256SUMS first, so two builds of the same tree are byte-identical):
`prepared/{manifest.json, <ds>/, images/<ds>/, caches/<ds>/}`, `views/<view>/` (the --view names, the run's training
view; default every `<ds>-*`), and `eval_views/<ds>-*/` (LIDC also `probe-*`; LIDC requires lidc-val-shifted,
lidc-val-rescaled, probe-precision and probe-topology, and `checker_eval/`, the CUDA checker eval set). Symlinks (also a
symlinked view or prepared directory, BK4-14) and anything under data/raw are refused; `build_report.json` (wall clock)
is skipped (BK1-7); `prepared/manifest.json` is rewritten to list only the bundled dataset's files (BK1-2). LIDC writes `bundle-lidc.tar`; RIGA streams the tar straight into
`age -r <riga_recipient from .vast-allowlist>`, so only `bundle-riga.tar.age` touches the disk (age output is not
byte-identical: random file key; the tar inside is). Each bundle gets `<bundle>.sha256` (sha256sum -c format);
sync_up checks it and SHA256SUMS on the box.
"""

import argparse
import hashlib
import io
import json
import subprocess
import tarfile
from pathlib import Path

from experiments.phase2.vast.policy import ROOT, die, is_raw, load
from objdel.hashing import file_sha256

LIDC_EVAL = ("lidc-val-shifted", "lidc-val-rescaled", "probe-precision", "probe-topology")
SKIP = ("build_report.json",)  # wall-clock timings: two builds of one view would differ (BK1-7)


def walk(data: Path, top: Path, prefix: str) -> list[tuple[str, Path]]:
    """(archive name, path) of every file under `top`; refuses a symlink anywhere between data/ and the files."""
    if top.resolve() != data.resolve() / top.relative_to(data):  # a symlink between data/ and top (BK4-14)
        die(f"refusing symlinked path {top}")
    paths = list(top.rglob("*"))  # lists symlinks without following them
    if links := [p for p in paths if p.is_symlink()]:
        die(f"refusing symlink {links[0]}")
    return [(f"{prefix}/{p.relative_to(top).as_posix()}", p) for p in paths if not p.is_dir() and p.name not in SKIP]


def manifest(ds: str, prep: Path) -> bytes:
    """prepared/manifest.json with only `ds`'s files and sets: the box has no other dataset's records (BK1-2)."""
    m = json.loads((prep / "manifest.json").read_text())
    m["files"], m["sets"] = ({k: v for k, v in m[s].items() if k.startswith(f"{ds}/")} for s in ("files", "sets"))
    if missing := [k for k in m["files"] if not (prep / k).is_file()]:
        die(f"prepared/manifest.json lists missing files: {missing}")
    return (json.dumps(m, indent=1, sort_keys=True) + "\n").encode()


def members(ds: str, data: Path, views: list[str] | None = None) -> list[tuple[str, Path | bytes]]:
    """The bundle's (archive name, path or bytes) entries for dataset `ds`; refuses anything under data/raw."""
    if is_raw(data):
        die("refusing to bundle anything under data/raw")
    prep = data / "prepared"
    picked = [f for sub in (ds, f"images/{ds}", f"caches/{ds}") for f in walk(data, prep / sub, f"prepared/{sub}")]
    for kind in ("views", "eval_views"):
        names = set(views) if kind == "views" and views else {p.name for p in (data / kind).glob(f"{ds}-*")}
        if kind == "views" and not views:  # RL train windows (build_eval_views --split train): only when named
            names = {n for n in names if not n.startswith(f"{ds}-train-")}
        if missing := [n for n in names if not (data / kind / n).is_dir()]:
            die(f"{kind} missing: {missing}")
        if kind == "eval_views" and ds == "lidc":
            if missing := [n for n in LIDC_EVAL if not (data / kind / n).is_dir()]:
                die(f"eval views missing: {missing}")
            names |= set(LIDC_EVAL)
        if kind == "eval_views":  # held-out test views never ride in a train bundle (plan B1 rental 4c)
            names = {n for n in names if not n.startswith(f"{ds}-test-")}
        if not names:
            die(f"no {ds} {kind} under {data / kind}")
        picked += [f for n in sorted(names) for f in walk(data, data / kind / n, f"{kind}/{n}")]
    if ds == "lidc":  # the CUDA checker eval set (LIDC only; archive/cpu_checker/experiments/prepare_data.py)
        if not (data / "checker_eval" / "manifest.json").is_file():
            die("data/checker_eval missing: it holds the archived CPU checker's references and cannot be rebuilt; "
                "copy it from a checkout that has it")
        picked += walk(data, data / "checker_eval", "checker_eval")
    if bad := [p for _, p in picked if is_raw(p) or not p.is_file()]:
        die(f"refusing (data/raw or not a file): {bad[:3]}")
    return sorted([*picked, ("prepared/manifest.json", manifest(ds, prep))])


def sha(p: Path | bytes) -> str:
    """SHA-256 hex of a file or of bytes."""
    return hashlib.sha256(p).hexdigest() if isinstance(p, bytes) else file_sha256(p)


def write_tar(fileobj, files: list[tuple[str, Path | bytes]]) -> None:
    """Members are files on disk or bytes built here (the filtered manifest); TarInfo's defaults are mtime 0, uid/gid
    0, mode 0644 and no owner names."""
    sums = "".join(f"{sha(p)}  {n}\n" for n, p in files).encode()
    with tarfile.open(fileobj=fileobj, mode="w|", format=tarfile.PAX_FORMAT) as t:
        for name, src in [("SHA256SUMS", sums), *files]:
            info = tarfile.TarInfo(name)
            info.size = len(src) if isinstance(src, bytes) else src.stat().st_size
            with io.BytesIO(src) if isinstance(src, bytes) else src.open("rb") as f:
                t.addfile(info, f)


def build(ds: str, data: Path, out: Path, views: list[str] | None = None) -> Path:
    """Write the bundle tar (RIGA: age-encrypted) and its .sha256 under `out`; returns the tar path."""
    files = members(ds, data, views)
    out.mkdir(parents=True, exist_ok=True)
    target = out / ("bundle-riga.tar.age" if ds == "riga" else f"bundle-{ds}.tar")
    tmp = target.with_name(target.name + ".part")
    if ds == "riga":  # D16: plaintext never on disk
        if not (rcpt := load(".vast-allowlist").get("riga_recipient", "")).startswith("age1"):
            die("no age riga_recipient recorded in .vast-allowlist")
        with subprocess.Popen(["age", "-r", rcpt, "-o", str(tmp)], stdin=subprocess.PIPE) as age:
            write_tar(age.stdin, files)
            age.stdin.close()
        if age.returncode:
            tmp.unlink(missing_ok=True)
            die("age encryption failed")
    else:
        with tmp.open("wb") as f:
            write_tar(f, files)
    tmp.replace(target)
    target.with_name(target.name + ".sha256").write_text(f"{sha(target)}  {target.name}\n")
    print(f"{target}: {len(files)} files, {target.stat().st_size / 1e9:.2f} GB")
    return target


def main(argv: list[str] | None = None) -> None:
    """The bundle CLI (dataset, --view, --data, --out)."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dataset", choices=["lidc", "riga"])
    ap.add_argument("--view", action="append", help="training view to ship (repeatable; default every <ds>-*)")
    ap.add_argument("--data", type=Path, default=ROOT / "data")
    ap.add_argument("--out", type=Path, default=ROOT / "build" / "bundles")
    a = ap.parse_args(argv)
    build(a.dataset, a.data, a.out, a.view)


if __name__ == "__main__":
    main()
