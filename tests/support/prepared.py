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
