"""C6: vLLM requests with n = G (train_rl --group-requests, off by default): one request per block of G identical
prompts instead of G requests of n = 1. vLLM V1 fans `n` out to child requests after input processing, so the saving
is host side (input processing of 32 requests, not 256). Not bit-equal to n = 1 (other per-request RNG streams)."""


def group_requests(prompts: list[dict], G: int) -> list[dict]:
    """prompts[::G] after checking each positional block of G rows holds one prompt (same ids, same image object);
    two neighbouring blocks may be equal (a repeated sample)."""
    if G < 1:
        raise ValueError(f"G = {G}: needs at least 1 answer per prompt")
    if len(prompts) % G:
        raise ValueError(f"{len(prompts)} prompts are not a multiple of G = {G}")
    for k in range(0, len(prompts), G):
        head = prompts[k]
        img = head.get("multi_modal_data", {}).get("image")
        for row in prompts[k + 1:k + G]:
            if (row["prompt_token_ids"] != head["prompt_token_ids"]
                    or row.get("multi_modal_data", {}).get("image") is not img):  # fmt: skip
                raise ValueError(f"block {k // G} (rows {k}..{k + G - 1}) mixes prompts: n = {G} needs one per block")
    return prompts[::G]


def group_params(params: list, G: int) -> list:
    """One clone at n = G per block of a per-row SamplingParams list (C11 caps); a block with two max_tokens raises."""
    out = []
    for k in range(0, len(params), G):
        if any(p.max_tokens != params[k].max_tokens for p in params[k + 1:k + G]):
            raise ValueError(f"block {k // G} (rows {k}..{k + G - 1}) has two caps: n = {G} sends one max_tokens")
        out.append(params[k].clone())
        out[-1].n = G
    return out


def grouped(generate, G: int):
    """generate(prompts, sampling_params=...) sending group_requests(prompts, G) with a clone of the sampling params
    at n = G (a per-row list: one clone per block, group_params); the returned RequestOutputs (one per block, G
    completions each) flatten to TRL's row order k*G + j."""

    def call(prompts, *args, sampling_params, **kwargs):
        sent = group_requests(prompts, G)
        if isinstance(sampling_params, list):
            if len(sampling_params) != len(prompts):
                raise ValueError(f"{len(sampling_params)} sampling params for {len(prompts)} prompts")
            sp = group_params(sampling_params, G)
        else:
            sp = sampling_params.clone()
            sp.n = G
        return generate(sent, *args, sampling_params=sp, **kwargs)

    return call
