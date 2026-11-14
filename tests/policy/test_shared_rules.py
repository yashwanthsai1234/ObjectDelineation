"""Shared train/eval rules after the D10 merge of their copies: the teacher-forced batch extension (MT-04), the
outline length rule (MT-06) and the answers row and completions grouping (MT-08)."""
import pytest  # noqa: I001

from objdel.core.io import write_jsonl
from objdel.policy import rows as policy_rows
pytest.importorskip("pandas")  # the [dev]-only CI job has no ML stack
pytest.importorskip("trl")  # the [dev]-only CI job has no ML stack

from objdel.policy import sampler as S

torch = pytest.importorskip("torch")


def _old_vision_extend(b, tgt):  # policy.vision.next_token_probs before MT-04
    b = dict(b)
    b["input_ids"] = torch.cat([b["input_ids"], tgt], 1)
    b["attention_mask"] = torch.ones_like(b["input_ids"])
    if "mm_token_type_ids" in b:
        b["mm_token_type_ids"] = torch.cat([b["mm_token_type_ids"], torch.zeros_like(tgt)], 1)
    return b


def test_extend_equals_the_vision_copy_on_a_qwen_shaped_batch():
    from objdel.policy.vision import extend

    n = 6  # Qwen3-VL keys and dtypes (all int64 but pixel_values); image_grid_thw is (1, 3), not token-aligned
    b = {"input_ids": torch.arange(n)[None], "attention_mask": torch.ones(1, n, dtype=torch.long),
         "mm_token_type_ids": torch.tensor([[0, 1, 1, 1, 0, 0]]), "pixel_values": torch.rand(4, 3),
         "image_grid_thw": torch.tensor([[1, 2, 2]])}
    tgt = torch.tensor([[11, 12, 13]])
    got, want = extend(b, tgt), _old_vision_extend(b, tgt)
    assert got.keys() == want.keys()
    for k in got:
        assert got[k].dtype == want[k].dtype and torch.equal(got[k], want[k]), k
    assert got["mm_token_type_ids"].tolist() == [[0, 1, 1, 1, 0, 0, 0, 0, 0]]


def test_outline_tokens_is_target_tokens_plus_the_end_token():
    def tok(text, add_special_tokens=True):
        assert add_special_tokens is False
        return {"input_ids": list(text)}

    assert policy_rows.target_tokens(tok, "abcd") == 4 and policy_rows.outline_tokens(tok, "abcd") == 5


def _rec(i, k, ids=(1, 7)):
    return {"id": i, "k": k, "text": f"{i}{k}", "token_ids": list(ids)}


def test_samples_and_import_share_one_checked_grouping(tmp_path):
    lines = [{"id": "b", "n": 2}, {"id": "a", "n": 2}]
    recs = [_rec("a", 1), _rec("a", 0, (1, 2)), _rec("b", 0), _rec("b", 1)]
    got = S.samples(lines, recs, {7})
    assert [s.item_id for s in got] == ["b", "a"]  # line order
    assert got[1].row() == {"item_id": "a", "samples": ["a0", "a1"], "truncated": [True, False], "tokens": [2, 2]}
    path = tmp_path / "c.jsonl"
    write_jsonl(path, recs, sort_keys=False)
    assert S.import_completions(path, {7}, n=2) == [s.row() for s in got[::-1]]  # first-seen order
    with pytest.raises(ValueError, match="k = 0 twice"):
        S.samples(lines, [*recs, _rec("a", 0)], {7})
    with pytest.raises(ValueError, match="want 2 samples"):  # a missing sample: ValueError, no bare KeyError
        S.samples(lines, recs[1:], {7})


def test_chat_turns_keep_their_dicts():
    from objdel.core import datasets
    from objdel.policy.prompt import _user, assistant_turn, user_turn

    assert user_turn("p") == {"role": "user", "content": [{"type": "image"}, {"type": "text", "text": "p"}]}
    assert _user("lidc") == user_turn(datasets.get("lidc").prompt)  # RL / SFT / eval rows and G-vision: one turn
    assert assistant_turn("t") == {"role": "assistant", "content": [{"type": "text", "text": "t"}]}
