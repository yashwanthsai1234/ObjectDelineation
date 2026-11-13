"""objdel.core.io.write_json / read_json (X2-6): the bytes every manifest writer wrote before (indent 1, sorted keys,
one trailing newline), so manifest hashes and goldens do not move."""
import hashlib
import json

import pytest

from objdel.core.io import json_text, plain, read_json, read_jsonl, read_jsonl_glob, read_split, write_json, write_jsonl


def test_write_json_keeps_the_manifest_bytes(tmp_path):
    doc = {"b": [1, 2.5, None], "a": {"z": "é", "y": True}}
    want = json.dumps(doc, indent=1, sort_keys=True) + "\n"
    write_json(tmp_path / "new" / "m.json", doc)  # atomic: makes the folder, leaves no .tmp
    assert (tmp_path / "new" / "m.json").read_text() == want
    assert [p.name for p in (tmp_path / "new").iterdir()] == ["m.json"]
    assert read_json(str(tmp_path / "new" / "m.json")) == doc


def test_hashing_helpers_agree(tmp_path):
    import numpy as np

    from objdel.core.hashing import bytes_sha256, file_sha256, pixels_sha256

    (tmp_path / "f").write_bytes(b"abc")
    assert bytes_sha256(b"abc") == file_sha256(tmp_path / "f") == hashlib.sha256(b"abc").hexdigest()
    a = np.arange(6, dtype=np.uint8).reshape(2, 3)
    assert pixels_sha256(a) == bytes_sha256(repr(a.shape).encode() + a.tobytes())


ROWS = [{"b": 1, "a": "é"}, {"z": [1.5, None], "y": {"k": True}}]


def test_write_jsonl_keeps_the_bytes_of_the_old_copies(tmp_path):
    old = "".join(json.dumps(r, sort_keys=True) + "\n" for r in ROWS)
    sha = write_jsonl(tmp_path / "deep" / "a.jsonl", iter(ROWS))  # an iterator streams; the folder is made
    assert (tmp_path / "deep" / "a.jsonl").read_text() == old
    assert sha == hashlib.sha256(old.encode()).hexdigest()
    old = "".join(json.dumps(r) + "\n" for r in ROWS)  # the unsorted copies
    write_jsonl(tmp_path / "b.jsonl", ROWS, sort_keys=False, atomic=True)
    assert (tmp_path / "b.jsonl").read_text() == old and [p.name for p in tmp_path.glob("b*")] == ["b.jsonl"]
    with pytest.raises(ValueError):
        write_jsonl(tmp_path / "c.jsonl", [{"x": float("nan")}], allow_nan=False)


def test_json_text_options_give_each_old_byte_shape():
    doc = {"b": float(1), "a": [1, {"d": 2, "c": 3}]}
    assert json_text(doc) == json.dumps(doc, indent=1, sort_keys=True) + "\n"
    assert json_text(doc, sort_keys=False) == json.dumps(doc, indent=1) + "\n"
    assert json_text(doc, newline=False) == json.dumps(doc, indent=1, sort_keys=True)
    assert json_text({"s": {1}}, default=str) == json.dumps({"s": {1}}, indent=1, sort_keys=True, default=str) + "\n"
    with pytest.raises(ValueError):
        json_text({"x": float("inf")}, strict=True)


def test_write_json_options_and_plain(tmp_path):
    import numpy as np

    write_json(tmp_path / "m.json", {"b": np.float64(2), "a": np.arange(2)}, sort_keys=False, strict=True,
               default=plain)
    assert (tmp_path / "m.json").read_text() == '{\n "b": 2.0,\n "a": [\n  0,\n  1\n ]\n}\n'


def test_read_split_and_glob_read_the_rows_in_name_order(tmp_path):
    write_jsonl(tmp_path / "lidc" / "val.jsonl", ROWS)
    write_jsonl(tmp_path / "lidc" / "test.jsonl", ROWS[:1])
    write_jsonl(tmp_path / "riga" / "val.jsonl", ROWS[1:])
    assert read_split(tmp_path, "lidc", "val") == ROWS == read_jsonl(tmp_path / "lidc" / "val.jsonl")
    assert read_jsonl_glob(tmp_path, "lidc/*.jsonl") == ROWS[:1] + ROWS  # test.jsonl, then val.jsonl
    assert read_jsonl_glob(tmp_path, "*/val.jsonl") == ROWS + ROWS[1:]
    assert read_jsonl_glob(tmp_path, "none/*.jsonl") == []
