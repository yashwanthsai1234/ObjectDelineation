"""A dataset the library does not know (registered in this test only: two raters, small uint8 RGB images, one
structure) goes store -> prepare -> RL prompt rows -> reward on CPU with no library change: the dataset entry is the
one place a dataset is described."""
import numpy as np
import pytest

from objdel.core import datasets
from objdel.data import prepare as data_prepare
from objdel.data import store
from objdel.data.store import Item
from objdel.policy import rows as policy_rows
from objdel.rl import run as rl_run
from objdel.rl.reward import TRLReward
from objdel.verify import store as serve
from objdel.verify.reward.spec import REWARD_V4

SYNTH = datasets.Dataset(
    name="synth",
    structures=(datasets.Structure("blob", nsd_tau=2.0, biou_d=2),),
    prompt='Outline the blob in this image. Answer with JSON only: {"blob":[{"outer":[[x,y],...],"holes":[]}]}, '
           'integers 0-1000 relative to the crop, outer rings clockwise.',
    pixels=("uint8", 3),
)
EOS = 1


def _items(disc_mask):
    """Four 96 x 96 images, two raters per item (a disc and a disc 1 px wider), one patient each."""
    rng = np.random.default_rng(0)
    out = []
    for i, (cy, cx, r) in enumerate([(40, 44, 12), (52, 48, 10), (48, 40, 14), (44, 52, 11)]):
        img = rng.integers(0, 256, (96, 96, 3), dtype=np.uint8)
        masks = {"blob": np.stack([disc_mask(96, 96, cy, cx, r + d) for d in (0, 1)])}
        out.append(Item(f"synth/img{i}", "synth", f"synth:p{i}", "train", img, masks, ["a", "b"], {}))
    return out


@pytest.fixture
def synth(monkeypatch):
    monkeypatch.setitem(datasets.DATASETS, SYNTH.name, SYNTH)
    return SYNTH


def test_a_new_dataset_goes_store_prepare_rows_reward_on_cpu(synth, disc_mask, tmp_path, clean_store):
    store.write(_items(disc_mask), tmp_path / "store")
    man = data_prepare.prepare(tmp_path / "store", tmp_path / "prep", splits=("train",), procs=1,
                               prompts={synth.name: synth.prompt})
    assert set(man["sets"]) == {"synth/train"}
    root = tmp_path / "prep"
    policy_rows.check_root(root, "synth", "train")
    assert rl_run.rl_refusal("synth") is None  # two raters: RL accepts it
    rows = policy_rows.prompts(root, "synth", "train")
    assert len(rows) == 4 and {r["key"] for r in rows} == {f"synth/img{i}" for i in range(4)}
    row = rows[0]
    assert row["names"] == ["blob"] and row["prompt"][0]["content"][1]["text"] == synth.prompt
    target = next(r for r in policy_rows.records(root, "synth", "train") if r["item_id"] == row["key"])["targets"][0]

    assert serve.preload(root, ["synth"]) == 4
    g = REWARD_V4.group
    texts = [target] * (g // 2) + ['{"blob":[]}'] * (g // 4) + ["not json"] * (g - g // 2 - g // 4)
    funcs, _ = TRLReward(serve.loader(root), {EOS}).funcs()
    kw = {"prompts": [row["prompt"]] * g, "completions": [[{"role": "assistant", "content": t}] for t in texts],
          "completion_ids": [[7, EOS]] * g, "key": [row["key"]] * g, "W": [row["W"]] * g, "H": [row["H"]] * g,
          "names": [row["names"]] * g}
    out = {f.__name__: f(**kw) for f in funcs}
    iou, adv = np.array(out["iou_tc"]), np.array(out[funcs[-1].__name__])
    assert (iou[: g // 2] > 0.8).all() and (iou[g // 2:] == 0).all()  # the rater's own outline scores
    assert np.isfinite(adv).all() and adv[: g // 2].min() > adv[g // 2:].max()
