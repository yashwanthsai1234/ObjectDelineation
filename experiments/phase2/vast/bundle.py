"""Data bundles for the Vast box (PLAN §4.4, D16): one deterministic tar per dataset with a SHA256SUMS manifest.

    .venv/bin/python -I -m experiments.phase2.vast.bundle lidc|riga [--data data] [--out build/bundles]

Members (sorted, mtime 0, uid/gid 0, mode 0644; SHA256SUMS first, so two builds of the same tree are byte-identical):
`prepared/{manifest.json, <ds>/, images/<ds>/, caches/<ds>/}`, `views/<ds>-*/`, and `eval_views/<ds>-*/` (LIDC also
`probe-*`; LIDC requires lidc-val-shifted, lidc-val-rescaled, probe-precision and probe-topology). Symlinks (also a
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
import os
import subprocess
import tarfile
from pathlib import Path

from experiments.phase2.vast import policy

LIDC_EVAL = ("lidc-val-shifted", "lidc-val-rescaled", "probe-precision", "probe-topology")


is_raw = policy.is_raw
SKIP = ("build_report.json",)  # wall-clock timings: two builds of one view would differ (BK1-7)


def walk(data: Path, top: Path, prefix: str) -> list[tuple[str, Path]]:
    if top.resolve() != data.resolve() / top.relative_to(data):  # a symlink between data/ and top (BK4-14)
        policy.die(f"refusing symlinked path {top}")
    out = []
    for dirpath, dirs, files in os.walk(top):
        for n in dirs + files:
            p = Path(dirpath) / n
            if p.is_symlink():
                policy.die(f"refusing symlink {p}")
        out += [(f"{prefix}/{(Path(dirpath) / f).relative_to(top).as_posix()}", Path(dirpath) / f)
                for f in files if f not in SKIP]  # fmt: skip
    return out


def manifest(ds: str, prep: Path) -> bytes:
    """prepared/manifest.json with only `ds`'s files and sets: the box has no other dataset's records (BK1-2)."""
    m = json.loads((prep / "manifest.json").read_text())
    m["files"], m["sets"] = ({k: v for k, v in m[s].items() if k.startswith(f"{ds}/")} for s in ("files", "sets"))
    if missing := [k for k in m["files"] if not (prep / k).is_file()]:
        policy.die(f"prepared/manifest.json lists missing files: {missing}")
    return (json.dumps(m, indent=1, sort_keys=True) + "\n").encode()


def members(ds: str, data: Path) -> list[tuple[str, Path]]:
    if is_raw(data):
        policy.die("refusing to bundle anything under data/raw")
    prep, picked = data / "prepared", []
    for sub in (ds, f"images/{ds}", f"caches/{ds}"):
        picked += walk(data, prep / sub, f"prepared/{sub}")
    for kind in ("views", "eval_views"):
        names = sorted(p.name for p in (data / kind).glob("*") if p.name.startswith(f"{ds}-"))
        if kind == "eval_views" and ds == "lidc":
            if missing := [n for n in LIDC_EVAL if not (data / kind / n).is_dir()]:
                policy.die(f"eval views missing: {missing}")
            names = sorted({*names, *LIDC_EVAL})
        if not names:
            policy.die(f"no {ds} {kind} under {data / kind}")
        for n in names:
            picked += walk(data, data / kind / n, f"{kind}/{n}")
    if bad := [p for _, p in picked if is_raw(p) or not p.is_file()]:
        policy.die(f"refusing (data/raw or not a file): {bad[:3]}")
    return sorted([*picked, ("prepared/manifest.json", manifest(ds, prep))])


def sha(p: Path | bytes) -> str:
    if isinstance(p, bytes):
        return hashlib.sha256(p).hexdigest()
    with p.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def write_tar(fileobj, files: list[tuple[str, Path | bytes]]) -> None:
    """Members are files on disk or bytes built here (the filtered manifest)."""
    sums = "".join(f"{sha(p)}  {n}\n" for n, p in files).encode()
    with tarfile.open(fileobj=fileobj, mode="w|", format=tarfile.PAX_FORMAT) as t:
        for name, src in [("SHA256SUMS", sums), *files]:
            info = tarfile.TarInfo(name)
            size = len(src) if isinstance(src, bytes) else src.stat().st_size
            info.size, info.mtime, info.mode, info.uname, info.gname = size, 0, 0o644, "", ""
            if isinstance(src, bytes):
                t.addfile(info, io.BytesIO(src))
            else:
                with src.open("rb") as f:
                    t.addfile(info, f)


def build(ds: str, data: Path, out: Path) -> Path:
    files = members(ds, data)
    out.mkdir(parents=True, exist_ok=True)
    target = out / ("bundle-riga.tar.age" if ds == "riga" else f"bundle-{ds}.tar")
    tmp = target.with_name(target.name + ".part")
    if ds == "riga":  # D16: plaintext never on disk
        rcpt = policy.allowlist().get("riga_recipient", "")
        if not rcpt.startswith("age1"):
            policy.die("no age riga_recipient recorded in .vast-allowlist")
        with subprocess.Popen(["age", "-r", rcpt, "-o", str(tmp)], stdin=subprocess.PIPE) as age:
            write_tar(age.stdin, files)
            age.stdin.close()
        if age.returncode:
            tmp.unlink(missing_ok=True)
            policy.die("age encryption failed")
    else:
        with tmp.open("wb") as f:
            write_tar(f, files)
    tmp.replace(target)
    target.with_name(target.name + ".sha256").write_text(f"{sha(target)}  {target.name}\n")
    print(f"{target}: {len(files)} files, {target.stat().st_size / 1e9:.2f} GB")
    return target


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dataset", choices=["lidc", "riga"])
    ap.add_argument("--data", type=Path, default=policy.ROOT / "data")
    ap.add_argument("--out", type=Path, default=policy.ROOT / "build" / "bundles")
    a = ap.parse_args(argv)
    build(a.dataset, a.data, a.out)


if __name__ == "__main__":
    main()
