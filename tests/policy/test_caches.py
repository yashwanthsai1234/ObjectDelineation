"""objdel.policy.caches: the vision and rope caches give TRL's log-probabilities and compute each image or row once;
rope positions equal Qwen3-VL's get_rope_index (from the old test_train_fast.py, E5)."""

import functools

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from objdel.policy import caches as policy_caches  # noqa: E402
from tests.experiments.test_trl_adapter_equivalence import (  # noqa: E402
    END,
    IMG,
    PATCH,
    START,
    Fast,
    Plain,
    batch,
    tiny_model,
    trainer,
)


def counted_tower(m):
    """Counts the images the real vision tower encodes (grid rows per call)."""
    seen, orig = [], m.model.visual.forward

    def forward(pixel_values, grid_thw=None, **kw):
        seen.append(len(grid_thw))
        return orig(pixel_values, grid_thw=grid_thw, **kw)

    m.model.visual.forward = forward
    return seen


def test_vision_cache_gives_trls_logps_and_encodes_each_image_once():
    m = tiny_model()
    ids, mask, T, kw = batch(pads=[2, 0, 5, 1, 0], lengths=[3, 7, 1, 2, 5], images=[0, 0, 1, 1, 0])  # 2 images
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
    seen = counted_tower(m)
    cache = policy_caches.VisionCache.install(m)
    t = trainer(Fast, buckets=(4, 8, 12, 16))
    for _ in range(3):  # old-logprob, reference and update passes of one generation batch
        lp, _, _ = t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
        torch.testing.assert_close(lp, ref, atol=1e-5, rtol=1e-5)
    assert sum(seen) == 2  # 5 rows x 3 passes = 15 image uses, 2 distinct images encoded
    cache.clear()  # a new generation batch: the next pass encodes again
    t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
    assert sum(seen) == 4


def test_rope_cache_gives_trls_logps_and_computes_each_row_once():
    """C5: Qwen3-VL's 3-D rope positions (get_rope_index: per-row host reads) are computed once per row per generation
    batch at the full width and sliced per chunk, not once per chunk in each of the old, reference and update passes.
    The update's micro-batches are other row subsets (TRL shuffles): rows are found by content."""
    m = tiny_model()
    ids, mask, T, kw = batch(pads=[2, 0, 5, 1, 0], lengths=[3, 7, 1, 2, 5])
    ref, _, _ = trainer(Plain)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
    rows, orig = [], m.model.get_rope_index
    m.model.get_rope_index = lambda input_ids, **k: rows.append(len(input_ids)) or orig(input_ids, **k)
    t = trainer(Fast, buckets=(4, 8, 12, 16), rope=policy_caches.RopeCache.install(m))
    for _ in range(2):  # old-logprob and reference passes over the whole batch
        lp, _, _ = t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
        torch.testing.assert_close(lp, ref, atol=1e-5, rtol=1e-5)
    perm = [3, 0, 4]  # an update micro-batch: a shuffled subset
    sub = {"pixel_values": kw["pixel_values"].reshape(5, 16, PATCH)[perm].reshape(-1, PATCH), "num_images": [1] * 3,
           "image_grid_thw": kw["image_grid_thw"][perm], "mm_token_type_ids": kw["mm_token_type_ids"][perm]}
    lp, _, _ = t._get_per_token_logps_and_entropies(m, ids[perm], mask[perm], T, **sub)
    torch.testing.assert_close(lp, ref[perm], atol=1e-5, rtol=1e-5)
    assert sum(rows) == 5  # 5 rows x 3 passes, each row's positions computed once
    t.rope.clear()  # a new generation batch
    t._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=2, **kw)
    assert sum(rows) == 10 and t.rope.pending is None


def test_rope_cache_keys_rows_by_their_exact_tokens_and_image_grids():
    """The prompt text is fixed per dataset, so a row's ids carry only its pad count t*h/2*w/2: grids (1, 4, 8) and
    (1, 8, 4) give the same ids but other image positions. Rows are keyed exactly, ids, mask and grids."""

    class Base(torch.nn.Module):
        calls = []

        def get_rope_index(self, input_ids, image_grid_thw, **k):
            self.calls.append(image_grid_thw.tolist())
            pos = image_grid_thw[:, 1:2].float().expand(len(input_ids), input_ids.size(1))  # h as the position
            return pos.unsqueeze(0).expand(3, -1, -1), None

    cache = policy_caches.RopeCache(Base())
    ids, mask, mm = torch.ones(3, 4, dtype=torch.long), torch.ones(3, 4, dtype=torch.long), torch.zeros(3, 4)
    grid = torch.tensor([[1, 4, 8], [1, 8, 4], [1, 4, 8]])
    pos = cache.positions(ids, mask, mm, grid, [1, 1, 1])
    assert Base.calls == [[[1, 4, 8], [1, 8, 4]]]  # row 2 repeats row 0 exactly: computed once
    assert pos[0, :, 0].tolist() == [4.0, 8.0, 4.0]


def test_rope_positions_are_computed_on_host_copies_and_returned_on_the_rows_device():
    """get_rope_index reads the host several times per row (~9 syncs on CUDA tensors); on the host copies keys()
    already makes, the same function makes none, and one copy returns the positions to the rows' device."""
    seen = []

    class Base(torch.nn.Module):
        def get_rope_index(self, input_ids, image_grid_thw, attention_mask, mm_token_type_ids):
            seen.append({t.device.type for t in (input_ids, image_grid_thw, attention_mask, mm_token_type_ids)})
            return torch.zeros(3, *input_ids.shape), None

    ids = torch.ones(2, 4, dtype=torch.long)
    pos = policy_caches.RopeCache(Base()).positions(ids, torch.ones_like(ids), torch.zeros_like(ids),
                                           torch.tensor([[1, 4, 8], [1, 8, 4]]), [1, 1])  # fmt: skip
    assert seen == [{"cpu"}] and pos.device == ids.device and pos.shape == (3, 2, 4)


def rope_rows(prompt, answers, pads, W):
    """Rows laid out as TRL builds them: left pads, the shared prompt, an answer of text tokens, right pads."""
    ids, mask = torch.zeros(len(answers), W, dtype=torch.long), torch.zeros(len(answers), W, dtype=torch.long)
    for r, (a, p) in enumerate(zip(answers, pads, strict=True)):
        row = prompt + a
        ids[r, p:p + len(row)], mask[r, p:p + len(row)] = torch.tensor(row), 1
    return ids, mask, (ids == IMG).long() * mask


def test_rope_positions_equal_get_rope_index_and_flag_rows_off_the_template():
    """C2: one get_rope_index call on the bare prompt gives a template; every row's positions follow from its real
    token count by one cumsum, gather and where. Rows that leave the template set `bad` (no host read inside)."""
    m, g = tiny_model(), torch.Generator().manual_seed(3)
    prompt = [5, 6, START, *[IMG] * 4, END, 7, 8, 9]
    grid = torch.tensor([[1, 4, 4]])
    ans = lambda n: torch.randint(0, 99, (n,), generator=g).tolist()  # noqa: E731
    answers, pads = [[], ans(1), ans(13), ans(20), ans(9)], [3, 0, 2, 0, 1]
    ids, mask, mm = rope_rows(prompt, answers, pads, W=34)
    mask[4, 1 + len(prompt) + 4] = 0  # a hole in an answer
    mask[3, len(prompt):] = 0  # a truncated row whose answer is masked
    mm = mm * mask
    tmpl = policy_caches.rope_template(m.model, torch.tensor(prompt), (torch.tensor(prompt) == IMG).long(), grid)
    pos, bad = policy_caches.rope_positions(tmpl, ids, mask, mm)
    want, _ = m.model.get_rope_index(ids, mm_token_type_ids=mm, image_grid_thw=grid.repeat(5, 1), attention_mask=mask)
    assert torch.equal(pos, want) and pos.dtype == torch.int64 and not bad.any()
    odd = ids.clone()
    odd[0, 3 + 1] = 9  # a prompt token differs
    mm_odd = mm.clone()
    mm_odd[2, 2 + len(prompt) + 3] = 1  # an image token in an answer
    short = mask.clone()
    short[1, 4:] = 0  # fewer real tokens than the prompt
    assert policy_caches.rope_positions(tmpl, odd, mask, mm)[1].tolist() == [True, False, False, False, False]
    assert policy_caches.rope_positions(tmpl, ids, mask, mm_odd)[1].tolist() == [False, False, True, False, False]
    assert policy_caches.rope_positions(tmpl, ids, short, mm * short)[1].tolist() == [False, True, False, False, False]


@functools.cache
def rope_model():
    return tiny_model()


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(rows=st.lists(st.tuples(st.integers(0, 40), st.integers(0, 6), st.booleans()), min_size=1, max_size=6),
       seed=st.integers(0, 2**16))  # fmt: skip
def test_rope_positions_equal_get_rope_index_on_random_rows(rows, seed):
    """C2 property: any answer lengths, left pads, a mask hole or a masked (truncated) answer: rope_positions equals
    get_rope_index on every column and no row is bad."""
    m, g = rope_model(), torch.Generator().manual_seed(seed)
    prompt = [5, 6, START, *[IMG] * 4, END, 7, 8, 9]
    answers = [torch.randint(0, 99, (n,), generator=g).tolist() for n, _, _ in rows]
    ids, mask, mm = rope_rows(prompt, answers, [p for _, p, _ in rows], W=6 + len(prompt) + 40)
    for r, (n, p, odd) in enumerate(rows):
        if odd and n > 1:
            mask[r, p + len(prompt) + int(torch.randint(0, n, (1,), generator=g))] = 0  # a hole in the answer
        elif odd:
            mask[r, p + len(prompt):] = 0  # a truncated answer, masked
    mm = mm * mask
    grid = torch.tensor([[1, 4, 4]])
    tmpl = policy_caches.rope_template(m.model, torch.tensor(prompt), (torch.tensor(prompt) == IMG).long(), grid)
    pos, bad = policy_caches.rope_positions(tmpl, ids, mask, mm)
    want, _ = m.model.get_rope_index(ids, mm_token_type_ids=mm, image_grid_thw=grid.repeat(len(rows), 1),
                                     attention_mask=mask)  # fmt: skip
    assert torch.equal(pos, want) and not bad.any()


def test_rope_positions_equal_get_rope_index_on_the_real_lidc_prompt():
    """C2: the real processor's LIDC prompt (325 tokens, grid (1, 32, 32)) on a model with the real vision token ids:
    answers of 0, 1, 131 and 1,536 tokens, left pads, a mask hole and a truncated masked row equal get_rope_index;
    an image pad in an answer sets bad."""
    import PIL.Image
    from transformers import AutoConfig, AutoProcessor

    from objdel.policy import prompt as dataset

    try:
        name = "Qwen/Qwen3-VL-2B-Instruct"
        proc, cfg = AutoProcessor.from_pretrained(name, local_files_only=True), AutoConfig.from_pretrained(
            name, local_files_only=True)
    except OSError:
        pytest.skip("the Qwen3-VL processor is not in the local HF cache")
    text = proc.apply_chat_template([dataset._user("lidc")], add_generation_prompt=True, tokenize=False)
    one = proc(text=[text], images=[PIL.Image.new("RGB", (512, 512))], return_tensors="pt")
    prompt, grid = one["input_ids"][0].tolist(), one["image_grid_thw"]
    assert len(prompt) == 325 and grid.tolist() == [[1, 32, 32]]
    m = tiny_model(vocab_size=cfg.text_config.vocab_size, img=cfg.image_token_id, video=cfg.video_token_id,
                   start=cfg.vision_start_token_id, end=cfg.vision_end_token_id)  # fmt: skip
    g = torch.Generator().manual_seed(5)
    answers = [torch.randint(0, 151000, (n,), generator=g).tolist() for n in (0, 1, 131, 1536, 131)]
    ids, mask, _ = rope_rows(prompt, answers, [3, 0, 2, 0, 1], W=325 + 1536 + 3)
    mm = (ids == cfg.image_token_id).long() * mask
    mask[2, 2 + 325 + 7] = 0  # a hole
    mask[4, 1 + 325:] = 0  # truncated and masked
    mm = mm * mask
    tmpl = policy_caches.rope_template(m.model, torch.tensor(prompt), one["mm_token_type_ids"][0], grid)
    pos, bad = policy_caches.rope_positions(tmpl, ids, mask, mm)
    want, _ = m.model.get_rope_index(ids, mm_token_type_ids=mm, image_grid_thw=grid.repeat(5, 1), attention_mask=mask)
    assert torch.equal(pos, want) and not bad.any()
    mm[3, 325 + 10] = 1  # an image token in an answer
    assert policy_caches.rope_positions(tmpl, ids, mask, mm)[1].tolist() == [False, False, False, True, False]


def test_rope_bad_rows_routes_on_host_data_only():
    """U-RP7: a prompt off its template, an answer holding the image pad id or the vision start id; a clean row."""
    tmpl, vision = (1, 2, 3), (IMG, 101, START, END)
    bad = policy_caches.rope_bad_rows([(1, 2, 3), (1, 2, 4), (1, 2, 3), (1, 2, 3)], [[7, 8], [7], [7, IMG], [START]],
                             [tmpl] * 4, vision)  # fmt: skip
    assert bad.dtype == bool and bad.tolist() == [False, True, True, True]
