"""End-to-end demo of the evaluation harness on synthetic "models" over the validation store (M8).

    .venv/bin/python -I -m experiments.phase1.eval_demo [--store data/store] [--split val] [--out data/eval/demo] \
        [--lidc-every 1] [--procs 3]

Models, each written as a completions JSONL and evaluated by experiments.phase1.eval_run:
- raters: one sample per rater, the rater's own SFT target (objdel.targets.build_target with the reward's size rule,
  objdel.targets.min_areas, as build_targets and prepare build it);
- honest: one realistic honest rollout per rater, as in the R3 gate (the rater's mask with a smooth boundary
  perturbation of amplitude U(0, tau), correlation length max(3, tau), keeping the rater's size-filtered topology;
  tests/oracles/reward_harness.topology_jitter, else the exact mask), as its SFT target;
- <adversary>: the R3 adversary (tests/oracles/reward_harness.adversaries: dilated union, hull, union, core,
  in-band speck, holes stamped, minority hole) as a deterministic model, 64 samples per item; RIGA objects without
  the adversary keep the vote; items where no object has it are skipped;
- <adversary>+honest: the adversary as one rollout in a group with the honest rollouts (the R3 setting).
The baseline audit (component alarm, envelope rises) is the honest model's on the same items: the SFT stand-in.
Writes <out>/<model>/{completions.jsonl, report.json, report.md, items.jsonl} and <out>/demo.md (comparison).
The R3 harness comes from tests/oracles (shared with the pytest gates; ISSUES I19). Reads untrusted data: run with
python -I.

References:
- R3 adversaries and the realistic honest rollout: plan §7 (R3), experiments/phase1/lab/t_reward/scripts/
  iter4_adversaries.py, ported to tests/oracles/reward_harness.py (impl/M6_reward.md).
- Rater targets: objdel.targets.build_target / min_areas (EXECUTION.md §3, M5; impl/M5_targets.md).
- Evaluation: experiments.phase1.eval_run, objdel.eval (impl/M8_eval.md); seeds from SHA-256 of the item id:
  https://docs.python.org/3/library/hashlib.html ,
  https://numpy.org/doc/stable/reference/random/generator.html#numpy.random.default_rng
"""
import argparse
import hashlib
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from experiments.phase1 import eval_run
from objdel.eval import hack_audit
from objdel.geometry import SCHEMAS
from objdel.reward import ExpertCache
from objdel.targets import build_target, min_areas
from tests.oracles import reward_harness as H

ADVERSARIES = ("dilated_union", "hull", "union", "core", "band_speck", "holes_stamped", "minority_hole")
K_ADVERSARY = 64


def _perturbed(m, cache, sig, rng):
    """R3 'real' honest rollout of a rater mask (its size-filtered signature kept, else the exact mask)."""
    x = H.topology_jitter(m, sig, cache.min_area, float(rng.uniform(0, cache.tau)), max(3.0, cache.tau), rng)
    return m if x is None else x


def _texts(args) -> tuple[str, dict]:
    """(item_id, {model: samples}) of one framed item."""
    item_id, fr = args
    names, masks = SCHEMAS[fr["dataset"]], fr["masks"]
    caches = {n: ExpertCache.build(masks[n]) for n in names}
    ma = min_areas(masks)  # the size rule of build_targets and train.dataset.prepare (with the RIGA rim)
    rng = np.random.default_rng(int.from_bytes(hashlib.sha256(item_id.encode()).digest()[:8], "little"))
    k = caches[names[0]].k
    raters = [build_target({n: masks[n][r] for n in names}, names, min_area=ma).text for r in range(k)]
    honest = []
    for r in range(k):
        try:
            honest.append(build_target({n: _perturbed(masks[n][r], caches[n], caches[n].rater_sigs[r], rng)
                                        for n in names}, names, min_area=ma).text)
        except ValueError:  # sub-filter specks beyond the 10-ring text cap: the rater's own target instead
            honest.append(raters[r])
    adv = {n: H.adversaries(masks[n], caches[n], rng) for n in names}
    out = {"raters": raters, "honest": honest}
    for a in ADVERSARIES:
        if not any(a in adv[n] for n in names):
            continue
        try:
            text = build_target({n: adv[n].get(a, caches[n].vote) for n in names}, names, min_area=ma).text
        except ValueError:  # not certifiable as text (e.g. > 10 rings): skipped
            continue
        out[a] = [text] * K_ADVERSARY
        out[f"{a}+honest"] = [text] + honest
    return item_id, out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--store", type=Path, default=Path("data/store"))
    ap.add_argument("--split", default="val")
    ap.add_argument("--out", type=Path, default=Path("data/eval/demo"))
    ap.add_argument("--lidc-every", type=int, default=1, help="use every n-th LIDC item")
    ap.add_argument("--procs", type=int, default=3)
    a = ap.parse_args()
    frames = eval_run.load_frames(a.store, a.split, every={"lidc": a.lidc_every})
    with Pool(a.procs) as pool:
        texts = dict(pool.map(_texts, frames.items(), chunksize=8))
    models = {m: [] for m in ["honest", "raters", *ADVERSARIES, *(f"{x}+honest" for x in ADVERSARIES)]}
    for item_id, per in texts.items():
        for m, samples in per.items():
            models[m].append({"item_id": item_id, "samples": samples, "truncated": [False] * len(samples)})
    reports, honest = {}, None
    for m, rows in models.items():
        if not rows:
            continue
        out = a.out / m
        out.mkdir(parents=True, exist_ok=True)
        (out / "completions.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
        results = eval_run.evaluate(rows, frames, a.procs)
        honest = honest or results  # the SFT stand-in, audited on the same items as each model
        ids = {r["item_id"] for r in rows}
        reports[m] = eval_run.build_report(results, {"audit": hack_audit([r for r in honest if r["item_id"] in ids])})
        eval_run.write(reports[m], results, out, f"{m} ({a.split})")
        print(m, len(rows), "items, fired:", reports[m]["audit"]["fired"], flush=True)
    (a.out / "demo.md").write_text(comparison(reports))
    print(comparison(reports))


def comparison(reports: dict) -> str:
    """One row per model and dataset; "fired" lists the alarms and flags that fired on that dataset."""
    f = eval_run._f
    L = ["| model | dataset | items | format | legal | IoU_tc median | D_max | topo-any | pass | TV dis. | "
         "min-mode dis. | soft Dice | GED | RHR / NRFR | fired |", "|" + "---|" * 15]
    for m, r in reports.items():
        checks = r["audit"]["alarms"] | r["audit"]["flags"]
        for ds, d in r["summary"]["datasets"].items():
            o, v = d["objects"].values(), r["audit"]["datasets"][ds]
            per = lambda key: " / ".join(f(x[key]) for x in o)  # noqa: E731
            fired = [k for k, c in checks.items() if any(x.split("/")[0] == ds for x in c["on"])]
            L.append(f"| {m} | {ds} | {d['items']} | {f(d['format'], 1)} | {f(d['legal'], 1)} | {per('iou_tc_median')} | "
                     f"{per('d_max_mean')} | {f(d['topo_any'], 1)} | {f(d['pass_rate'], 1)} | "
                     f"{f(d['disagreement']['tv'])} | {f(d['disagreement']['min_mode'])} | {per('soft_dice')} | "
                     f"{per('ged')} | {f(v['rhr'], 1)} / {f(v['nrfr'], 1)} | {', '.join(fired) or '–'} |")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()
