"""objdel.io.write_json / read_json (X2-6): the bytes every manifest writer wrote before (indent 1, sorted keys, one
trailing newline), so manifest hashes and goldens do not move."""
import json

from objdel.io import read_json, write_json


def test_write_json_keeps_the_manifest_bytes(tmp_path):
    doc = {"b": [1, 2.5, None], "a": {"z": "é", "y": True}}
    want = json.dumps(doc, indent=1, sort_keys=True) + "\n"
    write_json(tmp_path / "new" / "m.json", doc)  # atomic: makes the folder, leaves no .tmp
    assert (tmp_path / "new" / "m.json").read_text() == want
    assert [p.name for p in (tmp_path / "new").iterdir()] == ["m.json"]
    write_json(tmp_path / "plain.json", doc, atomic=False)
    assert (tmp_path / "plain.json").read_text() == want
    assert read_json(tmp_path / "plain.json") == read_json(str(tmp_path / "new" / "m.json")) == doc
