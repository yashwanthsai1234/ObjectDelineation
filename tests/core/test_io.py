"""objdel.core.io.write_json / read_json (X2-6): the bytes every manifest writer wrote before (indent 1, sorted keys,
one trailing newline), so manifest hashes and goldens do not move."""
import hashlib
import json

import pytest

from objdel.core.io import (
    files_under,
    json_text,
    plain,
    read_json,
    read_jsonl,
    read_jsonl_cut,
    read_jsonl_glob,
    read_split,
    write_json,
    write_jsonl,
)


def test_write_json_keeps_the_manifest_bytes(tmp_path):
    doc = {"b": [1, 2.5, None], "a": {"z": "é", "y": True}}
    want = json.dumps(doc, indent=1, sort_keys=True) + "\n"
    write_json(tmp_path / "new" / "m.json", doc)  # atomic: makes the folder, leaves no .tmp
    assert (tmp_path / "new" / "m.json").read_text() == want
    assert [p.name for p in (tmp_path / "new").iterdir()] == ["m.json"]
    assert read_json(str(tmp_path / "new" / "m.json")) == doc


def test_replace_dir_swaps_in_a_complete_folder(tmp_path):
    from objdel.core.io import replace_dir

    dest = tmp_path / "branches" / "e1"
    dest.mkdir(parents=True)
    (dest / "old.bin").write_text("old")
    (tmp_path / "branches" / ".e1.tmp").mkdir()  # a stale temp of a crashed write goes first
    (tmp_path / "branches" / ".e1.tmp" / "stale").write_text("x")

    def write(t):
        t.mkdir()
        (t / "new.bin").write_text("new")

    replace_dir(write, dest)
    assert sorted(p.name for p in dest.iterdir()) == ["new.bin"]
    assert sorted(p.name for p in dest.parent.iterdir()) == ["e1"]


def test_replace_dir_keeps_the_old_folder_when_the_write_fails(tmp_path):
    import pytest

    from objdel.core.io import replace_dir

    dest = tmp_path / "final"
    dest.mkdir()
    (dest / "w").write_text("kept")

    def write(t):
        t.mkdir()
        raise OSError("disk full")

    with pytest.raises(OSError):
        replace_dir(write, dest, tmp=tmp_path / "final.tmp")  # RL's save_final names its temp final.tmp
    assert (dest / "w").read_text() == "kept"


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


def test_read_json_default_is_for_a_missing_file_only(tmp_path):
    assert read_json(tmp_path / "none.json", {}) == {} and read_json(tmp_path / "none.json", None) is None
    assert read_json(tmp_path, "d") == "d"  # a folder is not a JSON file
    with pytest.raises(FileNotFoundError):
        read_json(tmp_path / "none.json")
    (tmp_path / "bad.json").write_text("{")
    with pytest.raises(ValueError):
        read_json(tmp_path / "bad.json", {})  # a broken file is an error, not the default


def test_read_jsonl_skip_blank(tmp_path):
    (tmp_path / "a.jsonl").write_text('{"a": 1}\n\n  \n{"a": 2}\n')
    assert read_jsonl(tmp_path / "a.jsonl", skip_blank=True) == [{"a": 1}, {"a": 2}]
    with pytest.raises(ValueError):
        read_jsonl(tmp_path / "a.jsonl")


def test_read_jsonl_cut_skips_a_cut_line_anywhere_and_keeps_the_raw_line(tmp_path):
    f = tmp_path / "c.jsonl"
    good1, good2 = '{"id": "a", "k": 0}', '{"id": "b", "k": 1, "x": [1]}'
    cut = ['{"id": "a", "k"', "[1, 2]", "7", '"s"', "", '{"id": "d", "k": 2']
    f.write_text("\n".join([good1, cut[0], cut[1], cut[2], cut[3], '{"id": "c"}', cut[4], good2, cut[5]]))
    assert read_jsonl_cut(f, ("id", "k")) == [(good1, {"id": "a", "k": 0}), (good2, {"id": "b", "k": 1, "x": [1]})]
    assert [c["id"] for _, c in read_jsonl_cut(f, ("id",))] == ["a", "c", "b"]  # fewer keys needed, one more line
    assert read_jsonl_cut(f) == [(x, json.loads(x)) for x in (good1, '{"id": "c"}', good2)]  # objects only
    assert read_jsonl_cut(tmp_path / "none.jsonl", ("id",)) == []
    with pytest.raises(ValueError):
        read_jsonl(f)  # the strict reader (select.complete) still raises on the same file


def test_files_under_sorts_by_the_posix_string_not_by_path(tmp_path):
    for rel in ("a/x", "a-b/x", "manifest.json", "a/sub/y", "z.txt"):
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel).write_text(rel)
    got = files_under(tmp_path)
    assert [r for r, _ in got] == sorted(r for r, _ in got) == ["a-b/x", "a/sub/y", "a/x", "manifest.json", "z.txt"]
    assert sorted(p for _, p in got) != [p for _, p in got]  # Path order differs: a/x before a-b/x
    assert [r for r, _ in files_under(tmp_path, skip=("manifest.json",))] == ["a-b/x", "a/sub/y", "a/x", "z.txt"]
    assert all(p == tmp_path / r for r, p in got)
