"""Upload a run's selected branch to the owner's private HF model repo (on the box; the Mac runs push_hf.sh).

    python -I -m experiments.phase2.hf_push --run RUN_DIR --repo svenk029/objdel-sft --revision r1-lfm < token
    (on the Mac, after pull.sh: add --view data/views/<name>, the local copy of the box's view; its view_hash must
    be the run's)

The token arrives on stdin (push_hf.sh streams it from the Mac keychain over ssh: never on argv, in the environment or
on disk); push_hf.sh refuses it unless token_problems is empty: fine-grained, write access to this one repo only (the
owner's free HF account has no card; storage only, never inference, endpoints or jobs). Refused here: a repo that is
not private, a run trained on anything but LIDC (RIGA's licence forbids redistribution, D16), a run with no selected
branch (select_ckpt: none passed D12). One commit on branch `revision` (created if missing): the files of
<run>/branches/<selected>, the base model's LICENSE (LFM Open License v1.0 §4(a); Apache-2.0 for Qwen) and a model
card stating the change (§4(b)), the base model, the selection and the LIDC-IDRI attribution (CC BY 3.0).

Reference: https://huggingface.co/docs/huggingface_hub/guides/upload#upload-files
"""

import argparse
import json
import sys
from pathlib import Path

REPO = "svenk029/objdel-sft"  # the owner's private model repo (2026-11-02)
TOKEN_OK = {"repo.content.read", "repo.write"}
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
    from huggingface_hub import HfApi

    return HfApi(token=token).whoami()


def api(token: str):
    from huggingface_hub import HfApi

    return HfApi(token=token)


def base_card(base: str) -> tuple[bytes, dict]:
    """The base model's LICENSE and card metadata (public repos: no token)."""
    from huggingface_hub import ModelCard, hf_hub_download

    return Path(hf_hub_download(base, "LICENSE", token=False)).read_bytes(), ModelCard.load(base, token=False).data.to_dict()


def card(base: str, meta: dict, man: dict, sel: dict, run: Path) -> str:
    import yaml

    c = man["config"]
    head = {k: meta[k] for k in ("license", "license_name", "license_link") if meta.get(k)}
    head |= {"base_model": base, "tags": ["medical-imaging", "segmentation", "lidc-idri"]}
    how = "LoRA, merged into the base weights" if c.get("lora") else "full fine-tune"
    failed = f" (failed: {', '.join(sel['gates_failed'])})" if sel["gates_failed"] else ""
    return (f"---\n{yaml.safe_dump(head, sort_keys=False)}---\n\n# objdel SFT {run.name}, branch {sel['selected']}\n\n"
            f"Modified from [{base}](https://huggingface.co/{base}) by the objdel project ({how}): every weight file "
            f"here differs from the base model's. The base licence is in LICENSE.\n\n"
            f"- Run {run.name}: seed {c['seed']}, WSD branches {c['branches']}; selected branch {sel['selected']} "
            f"(sampler {sel['sampler']}, K = {sel['k']}, {sel['items'] or 'all'} items); gates pass: "
            f"{sel['gates_pass']}{failed}\n"
            f"- Code SHA-256 {man['code']['objdel_and_script_sha256']}; training view hash {man['view']['view_hash']}\n\n"
            f"## Data\n\n{LIDC}\n\nResearch use only; not a medical device.\n")  # fmt: skip


def main(argv=None) -> int:
    from objdel.train import configs

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", type=Path, required=True, help="train_sft output dir (select/ inside)")
    ap.add_argument("--repo", required=True)
    ap.add_argument("--revision", required=True, help="the repo branch, e.g. r1-lfm")
    ap.add_argument("--view", type=Path, help="local copy of the training view (default: the manifest's path)")
    a = ap.parse_args(argv)
    token = sys.stdin.readline().strip() or sys.exit("hf_push: no token on stdin")
    sel = json.loads((a.run / "select" / "selected.json").read_text())
    if not sel["selected"]:
        sys.exit(f"hf_push: {a.run}: no branch selected (none passed D12); nothing to upload")
    man = json.loads((a.run / "manifest.json").read_text())
    view = json.loads(((a.view or Path(man["config"]["view"])) / "manifest.json").read_text())
    if view["view_hash"] != man["view"]["view_hash"]:
        sys.exit(f"hf_push: {a.view}: view_hash {view['view_hash']} is not the run's {man['view']['view_hash']}")
    if (ds := view["dataset"]) != "lidc":
        sys.exit(f"hf_push: {a.run} trained on {ds}: only LIDC runs are uploaded (RIGA's licence forbids "
                 "redistribution, D16)")  # fmt: skip
    hub = api(token)
    if not hub.repo_info(a.repo, repo_type="model").private:
        sys.exit(f"hf_push: {a.repo} is not private")
    from huggingface_hub import CommitOperationAdd

    base = configs.MODELS[man["config"]["model"]]
    lic, meta = base_card(base)
    src = a.run / "branches" / sel["selected"]
    ops = [CommitOperationAdd(f.relative_to(src).as_posix(), str(f)) for f in sorted(src.rglob("*")) if f.is_file()]
    ops += [CommitOperationAdd("LICENSE", lic), CommitOperationAdd("README.md", card(base, meta, man, sel, a.run).encode())]
    hub.create_branch(a.repo, branch=a.revision, repo_type="model", exist_ok=True)
    hub.create_commit(a.repo, operations=ops, commit_message=f"{a.run.name}: branch {sel['selected']}",
                      revision=a.revision, repo_type="model")  # fmt: skip
    print(f"hf_push: {src} -> {a.repo}@{a.revision}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
