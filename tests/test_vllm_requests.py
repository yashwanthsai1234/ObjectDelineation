"""C6: vllm_requests.group_requests / grouped (one n = G request per block of G identical prompts) and the FastGRPO
wiring behind train_rl --group-requests (outside the one-pad check)."""

import copy

import pytest

from objdel.train.vllm_requests import group_requests, grouped


def rows(*blocks, G=2):
    """G rows per (ids, image) block, the same dict contents and the same image object in a block (as OpenOnce)."""
    return [{"prompt_token_ids": list(ids), "multi_modal_data": {"image": img}}
            for ids, img in blocks for _ in range(G)]


def test_one_request_per_block_and_equal_neighbours_stay_two_requests():
    a, b = object(), object()
    p = rows(([1, 2], a), ([1, 2], a), ([3], b))
    assert group_requests(p, 2) == [p[0], p[2], p[4]]


def test_a_mixed_block_or_a_ragged_length_raises():
    a, b = object(), object()
    with pytest.raises(ValueError, match="multiple of G"):
        group_requests(rows(([1], a))[:1] + rows(([1], a)), 2)
    mixed = rows(([1], a), ([2], b))
    mixed[1] = mixed[2]
    with pytest.raises(ValueError, match="block 0"):
        group_requests(mixed, 2)
    same_ids_other_image = rows(([1], a))
    same_ids_other_image[1] = {"prompt_token_ids": [1], "multi_modal_data": {"image": copy.copy(a)}}
    with pytest.raises(ValueError, match="block 0"):
        group_requests(same_ids_other_image, 2)


def test_edge_cases_g_zero_empty_g_one_and_a_run_of_eight_plus_fours():
    """U-GR1: G = 0 is a ValueError (not ZeroDivisionError), [] gives [], G = 1 returns the list, 16 identical rows with
    G = 8 give 2 requests, and block lengths 8 + 4 + 4 are not blocks of 8."""
    a, b, c = object(), object(), object()
    with pytest.raises(ValueError, match="G = 0"):
        group_requests(rows(([1], a)), 0)
    assert group_requests([], 8) == []
    p = rows(([1], a), ([2], b), G=1)
    assert group_requests(p, 1) == p
    same = rows(([1], a), ([1], a), G=8)
    assert group_requests(same, 8) == [same[0], same[8]]
    ragged = rows(([1], a), G=8) + rows(([2], b), G=4) + rows(([3], c), G=4)
    with pytest.raises(ValueError, match="block 1"):
        group_requests(ragged, 8)
    for n in (15, 25):
        with pytest.raises(ValueError, match="multiple of G"):
            group_requests(rows(([1], a), G=n), 8)


def test_equal_ids_with_different_image_objects_raise_but_the_same_image_with_other_ids_too():
    a = object()
    equal_content = rows(([1], a), G=2)
    equal_content[1] = {"prompt_token_ids": [1], "multi_modal_data": {"image": object()}}
    with pytest.raises(ValueError, match="block 0"):
        group_requests(equal_content, 2)
    other_ids = rows(([1], a), G=2)
    other_ids[1] = {"prompt_token_ids": [2], "multi_modal_data": {"image": a}}
    with pytest.raises(ValueError, match="block 0"):
        group_requests(other_ids, 2)


def test_multi_modal_data_reaches_vllm_untouched():
    """U-GR2: the grouped requests are the caller's own dicts (PIL-style images, never image_embeds)."""
    a, b = object(), object()
    p = rows(([1], a), ([2], b), G=2)
    before = [dict(r, multi_modal_data=dict(r["multi_modal_data"])) for r in p]
    sent = group_requests(p, 2)
    assert all(s is r for s, r in zip(sent, p[::2], strict=True)) and p == before
    assert all(set(s["multi_modal_data"]) == {"image"} and s["multi_modal_data"]["image"] is img
               for s, img in zip(sent, (a, b), strict=True))


class Params:
    def __init__(self, max_tokens=None):
        self.n = 1
        self.max_tokens = max_tokens

    def clone(self):
        return copy.copy(self)


def test_grouped_sends_n_g_on_a_clone_and_returns_the_outputs_unchanged():
    sent = []

    def generate(prompts, sampling_params):
        sent.append((prompts, sampling_params.n))
        return ["out"] * len(prompts)

    sp, a = Params(), object()
    assert grouped(generate, 4)(rows(([1], a), ([2], a), G=4), sampling_params=sp) == ["out", "out"]
    assert [len(p) for p, _ in sent] == [2] and sent[0][1] == 4 and sp.n == 1


def test_grouped_with_caps_sends_one_cap_per_block_and_refuses_two():
    """U-RC3: a per-row SamplingParams list is cut like the prompts; each block must hold one max_tokens."""
    sent, a = [], object()

    def generate(prompts, sampling_params):
        sent.append(sampling_params)
        return prompts

    p = rows(([1], a), ([2], a), G=2)
    caps = [Params(9), Params(9), Params(5), Params(5)]
    grouped(generate, 2)(p, sampling_params=caps)
    assert [(s.max_tokens, s.n) for s in sent[0]] == [(9, 2), (5, 2)] and all(c.n == 1 for c in caps)
    with pytest.raises(ValueError, match="block 1.*caps"):
        grouped(generate, 2)(p, sampling_params=[Params(9), Params(9), Params(5), Params(6)])
    with pytest.raises(ValueError, match="4 prompts"):
        grouped(generate, 2)(p, sampling_params=caps[:3])


def test_fast_grpo_groups_outside_the_one_pad_check():
    pytest.importorskip("torch")
    pytest.importorskip("trl")
    from test_train_fast import END, IMG, START, vllm_trainer

    from objdel.train import fast

    t, sent = vllm_trainer(expand=2)
    t.group_requests, t.num_generations, t.num_generations_eval = True, 2, 4
    t.model.training = True
    fast.FastGRPO.patch_vllm_prompts(t)
    fast.FastGRPO.patch_vllm_prompts(t)  # idempotent
    img = object()
    row = {"prompt_token_ids": [START, IMG, IMG, END], "multi_modal_data": {"image": img}}
    out = t.vllm_generation.llm.generate([row, dict(row)], sampling_params=Params())
    assert len(out) == 1 and sent == [[{"prompt_token_ids": [START, IMG, END], "multi_modal_data": {"image": img}}]]
    t.model.training = False  # eval: num_generations_eval rows per prompt
    assert len(t.vllm_generation.llm.generate([row] * 4, sampling_params=Params())) == 1


def test_fast_grpo_caps_then_groups_one_request_per_block_with_its_cap():
    """U-RC3 wiring: --runaway-cap and --group-requests together send one n = G request per prompt, max_tokens = cap."""
    pytest.importorskip("torch")
    pytest.importorskip("trl")
    from test_train_fast import END, IMG, START, vllm_trainer

    from objdel.train import fast

    t, sent = vllm_trainer(expand=2)
    t.group_requests, t.runaway_cap, t.num_generations = True, True, 2
    t.model.training = True
    t.rows = [{"cap": 7}, {"cap": 7}, {"cap": 9}, {"cap": 9}]
    llm, params, engine = t.vllm_generation.llm, [], t.vllm_generation.llm.generate
    llm.generate = lambda prompts, sampling_params=None, **kw: params.append(sampling_params) or engine(prompts, **kw)
    fast.FastGRPO.patch_vllm_prompts(t)
    a, b = object(), object()
    rows = [{"prompt_token_ids": [START, IMG, IMG, END], "multi_modal_data": {"image": img}} for img in (a, a, b, b)]
    assert len(llm.generate(rows, sampling_params=Params(1536))) == 2
    assert [(p.max_tokens, p.n) for p in params[0]] == [(7, 2), (9, 2)] and len(sent[0]) == 2
