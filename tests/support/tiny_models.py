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
