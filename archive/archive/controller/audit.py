"""objdel.eval.audit.alarms_for as it was at 08e6455: hack_audit of a training evaluation, the Controller's
alarm input. hack_audit, metrics and ExpertCache stay live in objdel.eval and objdel.reward."""
from collections.abc import Callable, Iterable, Mapping, Sequence

from objdel.eval import hack_audit, metrics
from objdel.reward import ExpertCache


def alarms_for(records: Iterable[tuple[str, str, bool, int, int, Sequence[str]]],
               caches: Callable[[str], Mapping[str, ExpertCache]], baseline: dict | None = None) -> dict:
    """hack_audit of a training evaluation. records: one (key, text, truncated, W, H, names) per completion, in
    order (TRLReward's reward-call columns); key = item id starting with its dataset ("lidc/...", "riga/..."); a key's
    completions form one metrics item (split into GRPO groups of 8 there). caches(key) -> {object: ExpertCache}
    (objdel.train.dataset.cache_loader). baseline: an earlier alarms_for result on the same prompts (the SFT
    checkpoint's, e.g. the first evaluation). Returns the hack_audit dict; "fired" lists the alarms and flags."""
    by_key: dict[str, tuple] = {}
    for key, text, truncated, W, H, names in records:
        texts, flags, *_ = by_key.setdefault(key, ([], [], W, H, tuple(names)))
        texts.append(text)
        flags.append(bool(truncated))
    items = []
    for key, (texts, flags, W, H, names) in by_key.items():
        c = caches(key)
        items.append({"item_id": key, "dataset": key.split("/", 1)[0],
                      **metrics(texts, {n: c[n].raters for n in names}, W, H, names, flags, c)})
    return hack_audit(items, baseline)
