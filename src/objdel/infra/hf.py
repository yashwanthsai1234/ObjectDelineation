"""The owner's private HF model repo: `push` uploads a run's selected branch (on the box; the Mac runs push_hf.sh),
`card_only` rewrites its card's evaluation sections, `pull` fetches one revision onto the box (fetch_hf.sh). The CLIs
are objdel.cli.hub.push and objdel.cli.hub.pull (argument parsing, the stdin token, one call here).

    python -I -m objdel.cli.hub.push --run RUN_DIR --repo svenk029/objdel-sft --revision r1-qwen3vl2b < token
    python -I -m objdel.cli.hub.push --rl TRAIN_RL_OUT --repo svenk029/objdel-sft --revision b1-rows < token
    python -I -m objdel.cli.hub.push --rows data/views/lidc-train-b1-780 --repo svenk029/objdel-sft \
        --revision main
    (on the Mac, after pull.sh: add --view data/views/<name>, the local copy of the box's view; its view_hash must
    be the run's)

The token arrives on stdin (push_hf.sh streams it from the Mac keychain over ssh: never on argv, in the environment or
on disk); push_hf.sh refuses it unless token_problems is empty: fine-grained, write access to this one repo only (the
owner's free HF account has no card; storage only, never inference, endpoints or jobs). Refused here: a repo that is
not private, a run trained on anything but LIDC (RIGA's licence forbids redistribution), a run with no selected
branch (select_ckpt: none passed its gates). One commit on branch `revision` (created if missing): the files of
<run>/branches/<selected>, the base model's LICENSE (Apache-2.0 for Qwen) and a model card
stating the change, the base model, the selection and the LIDC-IDRI attribution (CC BY 3.0).
--card-only --headline reference/headline.json (no --run or --rl): the card's Evaluation and Limitations sections
at <revision>, written by card_update; the weights stay. The full-push card() is unchanged.
--rl (a train_rl output, always LIDC): <out>/final at the revision's root (the plain model folder fetch_hf.sh and
final_eval --ckpt read) and <out>/final_adapter as final_adapter/ (the unmerged LoRA), with the licence and a card;
refused without either folder.

Reference: https://huggingface.co/docs/huggingface_hub/guides/upload#upload-files
"""

import json
import re
from importlib import resources
from pathlib import Path

from objdel.core import datasets
from objdel.core.io import read_json
from objdel.infra.errors import InfraError

REPO = "svenk029/objdel-sft"  # the owner's private model repo
TOKEN_OK = {"repo.content.read", "repo.write"}
ROWS_TAR = "root.tar"  # a rows root on the Hub: rows/<name>/root.tar + manifest.json (the repo holds <= 20,000 files)
LIDC = ("Trained on views derived from LIDC-IDRI (CC BY 3.0). Data citation: Armato III SG et al., Data From "
        "LIDC-IDRI, The Cancer Imaging Archive, DOI 10.7937/K9/TCIA.2015.LO9QL9SX. We acknowledge the National Cancer "
        "Institute and the Foundation for the National Institutes of Health, and their critical role in the creation "
        "of the free publicly available LIDC/IDRI Database used in this study.")  # fmt: skip


def token_problems(who: dict, repo: str) -> list[str]:
    """whoami() of a token -> why it may not go to a rented box ([] = fine-grained, write on `repo` only)."""
    t = (who.get("auth") or {}).get("accessToken") or {}
    if t.get("role") != "fineGrained":
        return [f"role {t.get('role')} (a fine-grained token is needed)"]
    fg, out, write = t.get("fineGrained") or {}, [], False
    out += [f"global permission {p}" for p in fg.get("global", [])]
    for s in fg.get("scoped", []):
        e, perms = s.get("entity", {}), set(s.get("permissions", []))
        if (e.get("type"), e.get("name")) != ("model", repo):
            out.append(f"scope on {e.get('type')} {e.get('name')} (only model {repo} is allowed)")
            continue
        out += [f"permission {p} on {repo}" for p in sorted(perms - TOKEN_OK)]
        write |= "repo.write" in perms
    return out + ([] if write else [f"no write access to {repo}"])


def whoami(token: str) -> dict:
    """huggingface_hub whoami() of a token."""
    from huggingface_hub import HfApi

    return HfApi(token=token).whoami()


def api(token: str):
    """An HfApi client for a token."""
    from huggingface_hub import HfApi

    return HfApi(token=token)




def require_private(hub, repo: str) -> None:
    """InfraError unless `repo` is a private model repo."""
    if not hub.repo_info(repo, repo_type="model").private:
        raise InfraError(f"{repo} is not private")


def base_card(base: str) -> tuple[bytes, dict]:
    """The base model's LICENSE and card metadata (public repos: no token). A repo with no LICENSE file whose card
    declares apache-2.0 (Qwen/Qwen3-VL-2B-Instruct) gets the Apache-2.0 text (Apache-2.0 §4(a): give a copy)."""
    from huggingface_hub import ModelCard, hf_hub_download
    from huggingface_hub.errors import EntryNotFoundError

    meta = ModelCard.load(base, token=False).data.to_dict()
    try:
        return Path(hf_hub_download(base, "LICENSE", token=False)).read_bytes(), meta
    except EntryNotFoundError:
        if meta.get("license") != "apache-2.0":
            raise
        # apache.org/licenses/LICENSE-2.0.txt, sha256 cfc7749b...
        return resources.files("objdel.infra").joinpath("licenses/Apache-2.0.txt").read_bytes(), meta


def model_card(base: str, meta: dict, title: str, how: str, facts: str) -> str:
    """A model card: the YAML header (the base model's licence fields), `title`, the modification notice (`how`), the
    bullet lines `facts` and the LIDC attribution."""
    import yaml

    head = {k: meta[k] for k in ("license", "license_name", "license_link") if meta.get(k)}
    head |= {"base_model": base, "tags": ["medical-imaging", "segmentation", "lidc-idri"]}
    return (f"---\n{yaml.safe_dump(head, sort_keys=False)}---\n\n# {title}\n\n"
            f"Modified from [{base}](https://huggingface.co/{base}) by the objdel project ({how}): every weight file "
            f"here differs from the base model's. The base licence is in LICENSE.\n\n{facts}\n\n"
            f"## Data\n\n{LIDC}\n\nResearch use only; not a medical device.\n")


def card(base: str, meta: dict, man: dict, sel: dict, run: Path) -> str:
    """The model card of a train_sft run's selected branch."""
    c = man["config"]
    failed = f" (failed: {', '.join(sel['gates_failed'])})" if sel["gates_failed"] else ""
    facts = (f"- Run {run.name}: seed {c['seed']}, WSD branches {c['branches']}; selected branch {sel['selected']} "
             f"(sampler {sel['sampler']}, K = {sel['k']}, {sel['items'] or 'all'} items); gates pass: "
             f"{sel['gates_pass']}{failed}\n"
             f"- Code SHA-256 {man['code']['objdel_and_script_sha256']}; "
             f"training view hash {man['view']['view_hash']}")
    return model_card(base, meta, f"objdel SFT {run.name}, branch {sel['selected']}",
                      "LoRA, merged into the base weights; every model is LoRA", facts)


def rl_card(base: str, meta: dict, doc: dict, out: Path) -> str:
    """The model card of a train_rl output."""
    rows = doc.get("rows") or {}
    facts = (f"- {doc['steps']} steps{' (partial: merged from the last checkpoint)' * bool(doc.get('partial'))} "
             f"from {doc['ckpt']}; rows revision {rows.get('revision')}, ids SHA-256 {rows.get('ids_sha256')}")
    return model_card(base, meta, f"objdel RL {out.name}",
                      "GRPO on a LoRA adapter, merged into the base weights; the unmerged adapter is in final_adapter/",
                      facts)


def rl_ops(out: Path, hub, repo: str) -> list:
    """--rl: the commit's files (final/ at the root, final_adapter/ under its name, LICENSE, README.md)."""
    from huggingface_hub import CommitOperationAdd

    from objdel.policy import models

    for d in ("final", "final_adapter"):
        if not (out / d).is_dir():
            raise InfraError(f"{out}/{d} missing (train_rl --merge-last writes both from the last checkpoint)")
    require_private(hub, repo)
    doc = read_json(out / "rl.json")
    base = models.get(doc["model"]).hf_id
    lic, meta = base_card(base)
    files = [(f.relative_to(out / "final"), f) for f in sorted((out / "final").rglob("*")) if f.is_file()]
    files += [(f.relative_to(out), f) for f in sorted((out / "final_adapter").rglob("*")) if f.is_file()]
    card = rl_card(base, meta, doc, out).encode()
    return [CommitOperationAdd(r.as_posix(), str(f)) for r, f in files] + [CommitOperationAdd("LICENSE", lic),
                                                                           CommitOperationAdd("README.md", card)]


METRIC_NAMES = {"iou_tc": "IoU (against the expert band)", "f5": "boundary F-score", "precision_vote": "precision",
                "recall_vote": "recall", "dice_vote": "Dice", "legal": "legal outlines"}  # fmt: skip


def card_sections(h: dict) -> str:
    """The '## Evaluation' and '## Limitations' sections of the model card, from reference/headline.json; its
    optional `limitations` list adds the run's own lines (checkpoint selection, gate results) to the general ones."""
    t, v = h["test"], h["val"]
    rows = "\n".join(f"| {name} | {t['metrics'][k]['mean']:.3f} ({t['metrics'][k]['ci95'][0]:.3f} to "
                     f"{t['metrics'][k]['ci95'][1]:.3f}) | {v['metrics'][k]['mean']:.3f} |"
                     for k, name in METRIC_NAMES.items())
    return (f"## Evaluation\n\nRun {h['run']}, seed {h['seed']}, K = {h['K']} answers per image, 95 % intervals "
            f"({h['ci_unit']}-level bootstrap). Test: {t['items']} items on `{t['view']}`. Validation: {v['items']} "
            f"items on `{v['view']}`.\n\n| Metric | Test (95 % CI) | Validation |\n|---|---|---|\n{rows}\n\n"
            f"Sources: `{t['source']}`, `{v['source']}`.\n\n"
            "## Limitations\n\n"
            "- The test uses shifted views from one synthetic shift law, and training matches that law.\n"
            "- An image counts as a miss when its mean IoU over the K answers is below 0.1.\n"
            "- The test split informed the diagnosis of the plateau.\n"
            + "".join(f"- {line}\n" for line in h.get("limitations", []))
            + "- Training-seed variance is not measured: a gain is a single-run result. The intervals hold only the "
            "K = 8 answer noise and the image sample.\n"
            "- LIDC-IDRI only: RIGA is not redistributed.\n")  # fmt: skip


def card_update(old: str, h: dict) -> str:
    """Replace the Evaluation and Limitations sections of `old` (or append them); every other line is kept."""
    kept = re.sub(r"^## (?:Evaluation|Limitations)\n.*?(?=^## |\Z)", "", old, flags=re.M | re.S).rstrip("\n")
    return f"{kept}\n\n{card_sections(h)}"


def card_only(repo: str, revision: str, headline: Path, token: str) -> int:
    """--card-only: replace the card's evaluation sections at one revision. No weights, no GPU."""
    from huggingface_hub import hf_hub_download

    if problems := token_problems(whoami(token), repo):
        raise InfraError(f"the token may not be used here: {'; '.join(problems)}")
    hub = api(token)
    require_private(hub, repo)
    old = Path(hf_hub_download(repo, "README.md", revision=revision, token=token)).read_text()
    new = card_update(old, read_json(headline))
    hub.upload_file(path_or_fileobj=new.encode(), path_in_repo="README.md", repo_id=repo, revision=revision,
                    repo_type="model", commit_message="card: headline from reference/headline.json")
    print(f"hf_push: card of {repo}@{revision} updated")
    return 0



def push(run: Path | None, rl: Path | None, repo: str, revision: str, view: Path | None, token: str) -> int:
    """One commit on branch `revision` of `repo`: a train_sft run's selected branch (`run`) or a train_rl output
    (`rl`), with the base model's licence and a card; refused as the module docstring says."""
    from objdel.policy import models

    if rl:
        hub = api(token)
        ops = rl_ops(rl, hub, repo)
        hub.create_branch(repo, branch=revision, repo_type="model", exist_ok=True)
        hub.create_commit(repo, operations=ops, commit_message=f"{rl.name}: train_rl final and final_adapter",
                          revision=revision, repo_type="model")  # fmt: skip
        print(f"hf_push: {rl} -> {repo}@{revision}")
        return 0
    sel = read_json(run / "select" / "selected.json")
    if not sel["selected"]:
        raise InfraError(f"{run}: no branch selected (none passed its gates); nothing to upload")
    man = read_json(run / "manifest.json")
    vm = read_json((view or Path(man["config"]["view"])) / "manifest.json")
    if vm["view_hash"] != man["view"]["view_hash"]:
        raise InfraError(f"{view}: view_hash {vm['view_hash']} is not the run's {man['view']['view_hash']}")
    if (ds := vm["dataset"]) not in datasets.DATASETS or datasets.get(ds).restricted:
        raise InfraError(f"{run} trained on {ds}: not uploaded (unregistered, or its licence forbids redistribution: "
                         "objdel.core.datasets restricted)")
    hub = api(token)
    require_private(hub, repo)
    from huggingface_hub import CommitOperationAdd

    base = models.get(man["config"]["model"]).hf_id
    lic, meta = base_card(base)
    src = run / "branches" / sel["selected"]
    ops = [CommitOperationAdd(f.relative_to(src).as_posix(), str(f)) for f in sorted(src.rglob("*")) if f.is_file()]
    readme = card(base, meta, man, sel, run).encode()
    ops += [CommitOperationAdd("LICENSE", lic), CommitOperationAdd("README.md", readme)]
    hub.create_branch(repo, branch=revision, repo_type="model", exist_ok=True)
    hub.create_commit(repo, operations=ops, commit_message=f"{run.name}: branch {sel['selected']}",
                      revision=revision, repo_type="model")  # fmt: skip
    print(f"hf_push: {src} -> {repo}@{revision}")
    return 0


def rows_tar(root: Path, out) -> None:
    """ROOT's files except manifest.json as an uncompressed tar streamed to the binary file OUT, the same bytes for the
    same files (sorted paths; mtime 0, owner 0, mode 0644, no names), so preflight can hash it without a copy on disk
    and compare it with the Hub's LFS SHA-256 of ROWS_TAR. Directories are implied by the paths."""
    import tarfile

    files = sorted((p.relative_to(root).as_posix(), p) for p in root.rglob("*") if p.is_file())
    with tarfile.open(fileobj=out, mode="w|", format=tarfile.GNU_FORMAT) as tf:
        for rel, p in files:
            if rel == "manifest.json":
                continue
            ti = tarfile.TarInfo(rel)  # TarInfo's defaults: mode 0644, mtime 0, uid/gid 0, no user or group name
            ti.size = p.stat().st_size
            with p.open("rb") as fh:
                tf.addfile(ti, fh)


def push_rows(root: Path, repo: str, revision: str, token: str) -> int:
    """--rows: a build_rl_rows root to rows/<root name>/ on branch `revision` of `repo`. Refused without
    its manifest.json (written only when every hard data gate passed), for a set that is not LIDC, or to a public
    repo. Two files, two commits: ROWS_TAR (rows_tar: every file but the manifest; a root has ~2 files per window and
    the free repo holds at most 20,000), then the manifest alone, last; so the printed rows_revision, the last
    commit, holds the whole root. It and rows_manifest_sha256 are the rl stage's pins (specs/b1_rows.yaml);
    fetch_hf.sh unpacks the tar on the box and checks the root there."""
    import tempfile

    from objdel.core.hashing import file_sha256

    man = root / "manifest.json"
    if not man.is_file():
        raise InfraError(f"{man} missing (build_rl_rows writes it only when every hard gate passes); nothing uploaded")
    held = [k for k in read_json(man)["sets"] if (ds := k.split("/")[0]) not in datasets.DATASETS
            or datasets.get(ds).restricted]
    if held:
        raise InfraError(f"{root} holds {held}: not uploaded (unregistered, or a licence that forbids redistribution: "
                         "objdel.core.datasets restricted)")
    hub = api(token)
    require_private(hub, repo)
    prefix = f"rows/{root.name}"
    hub.create_branch(repo, branch=revision, repo_type="model", exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        tar = Path(tmp) / ROWS_TAR
        with tar.open("wb") as fh:
            rows_tar(root, fh)
        hub.upload_file(path_or_fileobj=str(tar), path_in_repo=f"{prefix}/{ROWS_TAR}", repo_id=repo,
                        revision=revision, repo_type="model", commit_message=f"{prefix}: rows root ({ROWS_TAR})")
    info = hub.upload_file(path_or_fileobj=str(man), path_in_repo=f"{prefix}/manifest.json", repo_id=repo,
                           revision=revision, repo_type="model", commit_message=f"{prefix}: manifest")
    print(f"hf_push: {root} -> {repo}@{revision} {prefix}/\nrows_revision: {info.oid}\n"
          f"rows_manifest_sha256: {file_sha256(man)}")
    return 0


def pull(repo: str, revision: str, out: str, allow: str | None, token: str) -> int:
    """One revision of `repo` into `out` (only the repo paths matching `allow`, kept at their repo paths); the token is
    passed to huggingface_hub directly, so nothing is cached on the box. A full model pull (no `allow`) also writes
    <out>/hf_revision.json ({repo, revision}): the F check's comparator runs its output-level arms only on the --ckpt
    that names rl-full-3ep."""
    from huggingface_hub import snapshot_download

    print(snapshot_download(repo_id=repo, revision=revision, local_dir=out, token=token, allow_patterns=allow))
    if allow is None:
        Path(out).mkdir(parents=True, exist_ok=True)
        (Path(out) / "hf_revision.json").write_text(json.dumps({"repo": repo, "revision": revision}) + "\n")
    return 0
