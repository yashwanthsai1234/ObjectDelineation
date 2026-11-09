"""C6: vLLM requests with n = G (train_rl --group-requests, off by default): one request per block of G identical
prompts instead of G requests of n = 1. vLLM V1 fans `n` out to child requests after input processing, so the saving
is host side (input processing of 32 requests, not 256). Not bit-equal to n = 1 (other per-request RNG streams)."""


def group_requests(prompts: list[dict], G: int) -> list[dict]:
    """prompts[::G] after checking each positional block of G rows holds one prompt (same ids, same image object);
    two neighbouring blocks may be equal (a repeated sample)."""
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


def grouped(generate, G: int):
    """generate(prompts, sampling_params=...) sending group_requests(prompts, G) with a clone of the sampling params
    at n = G; the returned RequestOutputs (one per block, G completions each) flatten to TRL's row order k*G + j."""

    def call(prompts, *args, sampling_params, **kwargs):
        sp = sampling_params.clone()
        sp.n = G
        return generate(group_requests(prompts, G), *args, sampling_params=sp, **kwargs)

    return call
