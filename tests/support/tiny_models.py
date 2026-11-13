"""Tiny offline language models for the TRL tests."""


def tiny_llama(n_words=40):
    """(model, tokenizer): a 16-wide 1-layer Llama (torch seed 0) over the word vocabulary
    ["<pad>", "<eos>", "prompt", "0" ... str(n_words - 1)]; <pad> = 0, <eos> = 1, left padding. The imports are inside:
    a run without transformers still collects the callers."""
    import torch
    from tokenizers import Tokenizer, models, pre_tokenizers
    from transformers import LlamaConfig, LlamaForCausalLM, PreTrainedTokenizerFast

    vocab = {w: i for i, w in enumerate(["<pad>", "<eos>", "prompt", *map(str, range(n_words))])}
    tk = Tokenizer(models.WordLevel(vocab, unk_token="<pad>"))
    tk.pre_tokenizer = pre_tokenizers.WhitespaceSplit()
    tok = PreTrainedTokenizerFast(tokenizer_object=tk, pad_token="<pad>", eos_token="<eos>", padding_side="left")
    torch.manual_seed(0)
    model = LlamaForCausalLM(LlamaConfig(vocab_size=len(vocab), hidden_size=16, intermediate_size=32,
                                         num_hidden_layers=1, num_attention_heads=2, num_key_value_heads=2,
                                         pad_token_id=0, eos_token_id=1, bos_token_id=1))
    return model, tok


def tiny_vlm(hub, text, vision, *, out=None, processor_kwargs=None):
    """A random small vision-language model from the cached config and processor of `hub` (no download; the test
    skips when they are not in the local HF cache). `text` and `vision` update the config's fields; torch seed 0.
    With `out`: save model and processor there and return `out`; else return (model, processor). HF_HUB_OFFLINE is
    set only while the cache is read."""
    import pytest

    torch = pytest.importorskip("torch")
    from transformers import AutoConfig, AutoModelForImageTextToText, AutoProcessor

    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("HF_HUB_OFFLINE", "1")
        try:
            cfg = AutoConfig.from_pretrained(hub, local_files_only=True)
            proc = AutoProcessor.from_pretrained(hub, local_files_only=True, **(processor_kwargs or {}))
        except OSError:
            pytest.skip(f"{hub} config / processor not in the local HF cache")
    cfg.text_config.update(text)
    cfg.vision_config.update(vision)
    torch.manual_seed(0)
    model = AutoModelForImageTextToText.from_config(cfg)
    if out is None:
        return model, proc
    model.save_pretrained(out)
    proc.save_pretrained(out)
    return out


QWEN_TINY = {  # model key -> (cached hub id, config fields of a random 2-layer model)
    "qwen3-vl-2b": ("Qwen/Qwen3-VL-2B-Instruct", dict(
        text=dict(hidden_size=64, intermediate_size=128, num_hidden_layers=2, num_attention_heads=4,
                  num_key_value_heads=2, head_dim=16,
                  rope_parameters={"rope_type": "default", "rope_theta": 5e6, "mrope_interleaved": True,
                                   "mrope_section": [2, 3, 3]}),
        vision=dict(hidden_size=32, intermediate_size=64, depth=2, num_heads=2, out_hidden_size=64,
                    deepstack_visual_indexes=[0, 1]))),
}


def tiny_qwen(out, model="qwen3-vl-2b"):
    """A random 2-layer `model` with its cached processor saved to `out` (see tiny_vlm); returns `out`."""
    hub, dims = QWEN_TINY[model]
    return tiny_vlm(hub, dims["text"], dims["vision"], out=out)
