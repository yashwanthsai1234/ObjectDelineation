"""B2.2 Mac check before the optional fp16 arm ($0, CPU): the whole model in float16, teacher-forced on saved answers.

    python -I -m tools.studies.fp16_check --ckpt runs/r1/r1-qwen3vl2b/branches/e2 --out RD/fp16_check.json

Records the max |activation| of every vision block and every text decoder layer and the count of inf/NaN in the
logits and the block outputs. Exit 1 (drop the fp16 arm) if any maximum exceeds 0.5 x 65,504 or anything is not
finite. Prompts and answers are b2_attribution's (`pick`, `batch`): the e2 completions of the R1 selection by default.
"""
import argparse
import json
import re
import sys
from pathlib import Path

import torch

LIMIT = 0.5 * 65504  # half the fp16 range
BLOCK = re.compile(r"(^|\.)(visual\.blocks|layers)\.\d+$")


def check(model, batches, forward) -> dict:
    """Runs forward(model, batch) on each batch under no_grad with a hook on every vision block and decoder layer.
    -> {max_abs: {block name: max |output|}, nonfinite_blocks, nonfinite_logits, over: [names above LIMIT], passed}"""
    top, bad = {}, [0]

    def hook(name):
        def f(_, __, out):
            x = (out[0] if isinstance(out, tuple) else out).detach().float()
            finite = torch.isfinite(x)
            bad[0] += int((~finite).sum())
            top[name] = max(top.get(name, 0.0), float(x.abs().where(finite, torch.zeros(())).max()))
        return f

    hooks = [m.register_forward_hook(hook(n)) for n, m in model.named_modules() if BLOCK.search(n)]
    logits = 0
    try:
        with torch.no_grad():
            for b in batches:
                out = forward(model, b)
                out = out.logits if hasattr(out, "logits") else out
                logits += int((~torch.isfinite(out)).sum())
    finally:
        for h in hooks:
            h.remove()
    over = [n for n, v in top.items() if v > LIMIT]
    return {"max_abs": top, "nonfinite_blocks": bad[0], "nonfinite_logits": logits, "over": over,
            "passed": not over and not bad[0] and not logits}  # fmt: skip


def exit_code(result: dict) -> int:
    return 0 if result["passed"] else 1


def main(argv=None) -> int:
    from tools.studies import b2_attribution as attr

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ckpt", type=Path, default=Path("runs/r1/r1-qwen3vl2b/branches/e2"))
    ap.add_argument("--answers", type=Path, default=attr.ANSWERS)
    ap.add_argument("--views", type=Path, default=Path("data/eval_views"))
    ap.add_argument("--prompts", type=int, default=32)
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    from transformers import AutoModelForImageTextToText, AutoProcessor

    from objdel.policy import models as configs

    proc = AutoProcessor.from_pretrained(a.ckpt, **configs.processor_kwargs("qwen3-vl-2b"))
    model = AutoModelForImageTextToText.from_pretrained(a.ckpt, dtype=torch.float16).eval()
    batches = [attr.batch(proc, p, ans, a.views, proc.tokenizer.pad_token_id)
               for p, ans in attr.pick(a.answers, a.prompts, a.k)]  # fmt: skip

    def forward(m, b):
        ids, mask, _, kw = b
        return m(input_ids=ids, attention_mask=mask, pixel_values=kw["pixel_values"].half(),
                 image_grid_thw=kw["image_grid_thw"], mm_token_type_ids=kw["mm_token_type_ids"])  # fmt: skip

    result = check(model, batches, forward) | {"ckpt": str(a.ckpt), "sequences": len(batches) * a.k}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "max_abs"}, indent=1))
    return exit_code(result)


if __name__ == "__main__":
    sys.exit(main())
