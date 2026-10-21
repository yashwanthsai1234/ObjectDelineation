"""B5 / H-A5-1, step 1 (training env, repo .venv): the prompts exactly as SFT and GRPO build them, for the vLLM check.

For N validation items per dataset (evenly spaced over the sorted item ids of data/prepared/<dataset>/val.jsonl) and
each model, with AutoProcessor.from_pretrained(repo, **objdel.train.configs.processor_kwargs(model, dataset)):
- gen: TRL 1.16 GRPOTrainer._tokenize_prompts: prepare_multimodal_messages([user turn], images=[image]) ->
  processor.apply_chat_template(add_generation_prompt=True, tokenize=True, return_dict=True) (chat_template None)
  -> input_ids; the same call with tokenize=False gives the text a standalone vLLM sampler sends with the image;
- sft: TRL 1.16 SFTTrainer's DataCollatorForVisionLanguageModeling(assistant_only_loss=True, chat_template=
  get_training_chat_template(processor)) on the SFT row (user turn + the item's first target): input_ids before the
  first loss token = the prompt the model was trained on. sft_prefix_equal says whether it equals gen;
  sft_starts_with_gen whether the SFT sequence starts with gen (same tokens, loss span may start earlier).
Writes JSON (ids, texts, image paths, processor kwargs) to argv[2]; the image files stay where they are (data/ never
leaves the Mac). Run with -I (reads data/):
    .venv/bin/python -I -m experiments.phase2.lab.B5.export_prompts data/prepared <out.json> [N]
"""
import json
import sys
from pathlib import Path

from PIL import Image
from transformers import AutoProcessor
from trl.chat_template_utils import get_training_chat_template
from trl.data_utils import prepare_multimodal_messages
from trl.trainer.sft_trainer import DataCollatorForVisionLanguageModeling

from objdel.train.configs import processor_kwargs
from objdel.train.dataset import _image, _user, records

MODELS = {"lfm2.5-vl-450m": "LiquidAI/LFM2.5-VL-450M", "qwen3.5-0.8b": "Qwen/Qwen3.5-0.8B",
          "smolvlm-256m": "HuggingFaceTB/SmolVLM-256M-Instruct"}


def pick(root: str, dataset: str, n: int) -> list[dict]:
    recs = sorted((r for r in records(root, dataset, "val") if any(r["targets"])), key=lambda r: r["item_id"])
    step = max(1, len(recs) // n)
    return recs[::step][:n]


def main(root: str, out: str, n: int) -> None:
    res = {"n_per_dataset": n, "models": {}}
    for model, repo in MODELS.items():
        res["models"][model] = {"repo": repo, "datasets": {}}
        for dataset in ("lidc", "riga"):
            kw = processor_kwargs(model, dataset)
            proc = AutoProcessor.from_pretrained(repo, **kw)
            try:
                tmpl = get_training_chat_template(proc)
                collator = DataCollatorForVisionLanguageModeling(proc, assistant_only_loss=True, chat_template=tmpl)
                sft_err = None
            except Exception as e:  # noqa: BLE001 - recorded: a model TRL cannot train with assistant-only loss
                collator, sft_err = None, f"{type(e).__name__}: {e}"[:400]
            items = []
            for r in pick(root, dataset, n):
                path = _image(root, r["item_id"])
                img = Image.open(path).convert("RGB")
                prompt = prepare_multimodal_messages([_user(dataset)], images=[img])
                tok = proc.apply_chat_template(conversation=[prompt], add_generation_prompt=True, tokenize=True,
                                               return_dict=True)
                gen_ids = [int(x) for x in tok["input_ids"][0]]
                text = proc.apply_chat_template(conversation=[prompt], add_generation_prompt=True, tokenize=False)
                text = text[0] if isinstance(text, list) else text
                rec = {"item_id": r["item_id"], "image": path, "image_size": list(img.size), "text": text,
                       "gen_ids": gen_ids, "target": next(t for t in r["targets"] if t)}
                if collator is not None:
                    answer = {"role": "assistant", "content": [{"type": "text", "text": rec["target"]}]}
                    b = collator([{"messages": [_user(dataset), answer], "image": img}])
                    ids, labels = b["input_ids"][0].tolist(), b["labels"][0].tolist()
                    first = next(i for i, v in enumerate(labels) if v != -100)
                    rec["sft_prefix_len"] = first
                    rec["sft_prefix_equal"] = ids[:first] == gen_ids
                    # same token sequence, loss span starting inside the generation prompt (Qwen3.5: the empty
                    # <think></think> block that the generation prompt already contains is trained as completion)
                    rec["sft_starts_with_gen"] = ids[:len(gen_ids)] == gen_ids
                    rec["gen_prompt_tokens_in_loss"] = max(0, len(gen_ids) - first)
                    if not rec["sft_prefix_equal"]:
                        rec["sft_prefix_tail"] = proc.tokenizer.decode(ids[max(0, first - 12):first])
                        rec["gen_tail"] = proc.tokenizer.decode(gen_ids[-12:])
                items.append(rec)
            res["models"][model]["datasets"][dataset] = {"processor_kwargs": kw, "sft_error": sft_err,
                                                          "items": items}
            eq = [x.get("sft_prefix_equal") for x in items]
            sw = [x.get("sft_starts_with_gen") for x in items]
            print(model, dataset, "prompt lens", sorted({len(x["gen_ids"]) for x in items}),
                  "sft prefix equal", sum(bool(e) for e in eq), "/", len(eq), "sft starts with gen prompt",
                  sum(bool(e) for e in sw), "/", len(sw), "gen tokens in loss",
                  sorted({x.get("gen_prompt_tokens_in_loss") for x in items}), sft_err or "", flush=True)
    Path(out).write_text(json.dumps(res))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 32)
