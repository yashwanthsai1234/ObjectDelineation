"""Prepared roots built from synthetic items."""


def prepare_items(items, base, *, splits=None, out="prepared", procs=1) -> dict:
    """Write `items` to the store `base/store`, run data.prepare into `base/<out>` (every split unless `splits`),
    and return its manifest. The imports are inside: a CPU run without `datasets` still collects the callers."""
    from objdel.data import prepare as data_prepare
    from objdel.data import store
    from objdel.policy import prompt as policy_prompt

    store.write(items, base / "store")
    extra = {} if splits is None else {"splits": splits}
    return data_prepare.prepare(base / "store", base / out, **extra, procs=procs, prompts=policy_prompt.PROMPTS)


def write_prepared_root(root, sets, *, manifest_extra=None):
    """Write the smallest prepared root the verifier reads: for each set key ("lidc/val") a `<key>.jsonl` of its
    records, `caches/<item_id>/<name>.npz` for each record's `caches` {name: rater stack (k, H, W)} (the record keeps
    the digest in their place), and a manifest.json with the SHA-256 of each jsonl plus `manifest_extra`."""
    import json

    from objdel.core.hashing import file_sha256
    from objdel.verify.reward.cache import ExpertCache

    files = {}
    for key, records in sets.items():
        lines = []
        for rec in records:
            caches = {}
            for name, stack in rec["caches"].items():
                (root / "caches" / rec["item_id"]).mkdir(parents=True, exist_ok=True)
                caches[name] = ExpertCache.build(stack).save(root / "caches" / rec["item_id"] / f"{name}.npz")
            lines.append(json.dumps(rec | {"caches": caches}))
        path = root / f"{key}.jsonl"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(line + "\n" for line in lines))
        files[f"{key}.jsonl"] = {"sha256": file_sha256(path)}
    (root / "manifest.json").write_text(json.dumps({"files": files} | (manifest_extra or {})))
    return root
