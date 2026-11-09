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


class Params:
    def __init__(self):
        self.n = 1

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
