"""objdel.core.io.write_json / read_json (X2-6): the bytes every manifest writer wrote before (indent 1, sorted keys,
one trailing newline), so manifest hashes and goldens do not move."""
import json

from objdel.core.io import read_json, write_json


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
