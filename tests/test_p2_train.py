"""Phase 2 M2 / M9 (train): view_rows over a training view directory (C-ROWS), the Phase 2 configs, and the I45
dataset fixes (prepare checks cv2.imwrite, cache_loader reads the manifest, unbounded preload)."""
import hashlib
import json
import pickle
import re
import shutil
import sys
import types

import numpy as np
import pytest

pytest.importorskip("datasets")
pytest.importorskip("trl")

from objdel.data import store  # noqa: E402
from objdel.geom.shape import canonicalize, to_mask  # noqa: E402
from objdel.geom.text import SCHEMAS, parse, serialize  # noqa: E402
from objdel.targets import build_target, min_areas  # noqa: E402
from objdel.train import configs  # noqa: E402
from objdel.train import dataset as D  # noqa: E402
from objdel.verify import store as verify_store  # noqa: E402

# ------------------------------------------------------------------ a tiny training view directory (C-VIEW)


def _target(disc, dataset, *discs):
    names = SCHEMAS[dataset]
    stacks = {n: np.stack([disc(50, 50, *d) for d in ds]) for n, ds in zip(names, discs, strict=True)}
    ma = min_areas(stacks)
    return [build_target({n: stacks[n][r] for n in names}, names, None, ma).text for r in range(len(discs[0]))]


def _row(i, item, dataset, rater, target, op="id", rep=0, simple=True):
    return {"row_id": f"{item}#{i}", "epoch": 0, "item_id": item, "dataset": dataset, "rater": rater,
            "rater_index": int(rater) - 1, "op": op, "offset": [i, -i], "side": 50,
            "window": [-1000, 400] if dataset == "lidc" else None, "noise_sigma": 0.0, "noise_seed": i,
            "simple": simple, "rep": rep, "target": target}


@pytest.fixture(scope="module")
def view_dir(tmp_path_factory, disc_mask):
    """Three LIDC rows (a solid disc, a ring with a hole repeated twice) and one RIGA row."""
    root = tmp_path_factory.mktemp("view") / "lidc-train-e0"
    solid, ring = _target(disc_mask, "lidc", [(25, 25, 10), (25, 25, 10, 4)])
    (rim,) = _target(disc_mask, "riga", [(25, 25, 15)], [(25, 25, 6)])
    rows = [_row(0, "lidc/1.2.1/n000/z001", "lidc", "1", solid),
            _row(1, "lidc/1.2.1/n000/z001", "lidc", "2", ring, "r90", 0, False),
            _row(2, "lidc/1.2.1/n000/z001", "lidc", "2", ring, "fh", 1, False),
            _row(3, "riga/MESSIDOR/image001", "riga", "1", rim)]
    text = "".join(json.dumps(r, sort_keys=True) + "\n" for r in rows)
    (root / "context").mkdir(parents=True)
    (root / "rows.jsonl").write_text(text)
    rng = np.random.default_rng(0)
    contexts = {}
    for item, px in (("lidc/1.2.1/n000/z001", rng.integers(-1000, 400, (60, 60), dtype=np.int16)),
                     ("riga/MESSIDOR/image001", rng.integers(0, 256, (60, 60, 3), dtype=np.uint8))):
        path = root / "context" / f"{item}.npz"
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(path, pixels=px, x0=0, y0=0)
        contexts[item] = hashlib.sha256(path.read_bytes()).hexdigest()
    (root / "manifest.json").write_text(json.dumps({"view": root.name, "contexts": contexts, "view_hash":
                                                    hashlib.sha256(text.encode()).hexdigest()}))
    return root, rows


@pytest.fixture
def fake_augment(monkeypatch):
    """objdel.targets.augment per C-VIEW until agent V lands it: render(spec, context) -> uint8 (512, 512, 3)."""
    calls = []

    def render(spec, context):
        calls.append(spec)
        return np.full((512, 512, 3), (spec["offset"][0] * 40) % 256, np.uint8)

    def load_context(view_dir, item_id):
        with np.load(f"{view_dir}/context/{item_id}.npz") as z:
            return {k: z[k] for k in z}

    fake = types.SimpleNamespace(render=render, load_context=load_context, calls=calls)
    import objdel.targets
    monkeypatch.setitem(sys.modules, "objdel.targets.augment", fake)
    monkeypatch.setattr(objdel.targets, "augment", fake, raising=False)
    return fake


def chars(text, add_special_tokens=False):  # a one-token-per-character tokenizer
    return {"input_ids": list(text)}


def test_view_rows_schema_and_lazy_render(view_dir, fake_augment):
    root, rows = view_dir
    ds = D.view_rows(root, chars)
    assert len(ds) == 4 and not fake_augment.calls  # building renders nothing
    r = ds[1]
    assert {"messages", "image", "key", "rater", "row_id", "epoch", "length"} <= set(r) and "target" not in r
    assert ds["epoch"] == [0] * 4  # S1-7: the view epoch as a column
    user, answer = r["messages"]
    assert [c["type"] for c in user["content"]] == ["image", "text"]  # Arrow adds "text": None to the image part
    assert user["content"][1]["text"] == D.PROMPTS["lidc"]  # as sft_rows
    assert answer == {"role": "assistant", "content": [{"type": "text", "text": rows[1]["target"]}]}
    assert (r["key"], r["rater"], r["row_id"], r["length"]) == (rows[1]["item_id"], "2", rows[1]["row_id"],
                                                               len(rows[1]["target"]))
    assert r["image"].mode == "RGB" and r["image"].size == (512, 512) and r["image"].getpixel((0, 0))[0] == 40
    assert fake_augment.calls == [{k: v for k, v in rows[1].items() if k != "target"}]  # the row is the spec
    assert ds[3]["messages"][0]["content"][1]["text"] == D.PROMPTS["riga"]


def test_length_column_and_the_sampler_render_nothing(view_dir, fake_augment):
    """transformers' group_by_length sampler reads ds["length"] once per epoch; a plain with_transform would render
    every row there (impl note M2: 1,005 renders per 1,000 rows)."""
    from transformers.trainer_pt_utils import LengthGroupedSampler
    root, rows = view_dir
    ds = D.view_rows(root, chars)
    assert ds["length"] == [len(r["target"]) for r in rows] and "length" in ds.column_names
    assert sorted(LengthGroupedSampler(2, dataset=ds, lengths=ds["length"])) == [0, 1, 2, 3]
    assert not fake_augment.calls
    assert len(ds.__getitems__([0, 2])) == 2 and len(fake_augment.calls) == 2  # DataLoader's batched fetch
    assert "image" in next(iter(ds))  # what SFTTrainer checks to pick its vision collator


def test_view_rows_pickle_for_dataloader_workers(view_dir, fake_augment):
    root, _ = view_dir
    ds = pickle.loads(pickle.dumps(D.view_rows(root, chars)))
    assert type(ds).__name__ == "_ViewRows" and ds[0]["image"].size == (512, 512)


def test_view_rows_drop_answers_beyond_the_budget(view_dir, fake_augment):
    root, rows = view_dir
    n = sorted(len(r["target"]) for r in rows)
    kept = D.view_rows(root, chars, max_tokens=n[1] + 1)  # + 1: the end-of-turn token (I5)
    assert len(kept) == sum(k <= n[1] for k in n) < len(rows)
    assert D.view_rows(root)["length"] == [len(r["target"]) for r in rows]  # no tokenizer: characters
    with pytest.raises(ValueError, match="tokenizer"):
        D.view_rows(root, max_tokens=1024)  # BK1-5: a token budget needs the tokenizer


def test_a0_empties_holes_and_keeps_parts(view_dir, fake_augment):
    root, rows = view_dir
    a0 = D.view_rows(root, chars, a0=True)
    for r, row in zip(rows, a0, strict=True):
        names = SCHEMAS[r["dataset"]]
        text = row["messages"][1]["content"][0]["text"]
        before, after = parse(r["target"], names).objs, parse(text, names).objs
        assert serialize(canonicalize(after), names) == text  # canonical
        assert all(not hs for n in names for _, hs in after[n]) and '"holes":[[' not in text
        assert {n: [o.tolist() for o, _ in before[n]] for n in names} == {n: [o.tolist() for o, _ in after[n]]
                                                                           for n in names}
        assert row["length"] == len(text)
    assert '"holes":[[' in rows[1]["target"]  # the ring row had a hole to drop


def test_view_rows_refuse_a_changed_rows_file(view_dir, tmp_path, fake_augment):
    root, _ = view_dir
    shutil.copytree(root, tmp_path / "v")
    with open(tmp_path / "v" / "rows.jsonl", "a") as f:
        f.write("\n")
    with pytest.raises(ValueError, match="view_hash"):
        D.view_rows(tmp_path / "v")


def test_a0_drops_a_part_inside_another_parts_hole(disc_mask):
    """BK1-4: an island in a ring's hole would become a second outer ring nested in the ring's once holes are
    filled (illegal); A0 keeps only the outermost parts, so its mask is the filled union."""
    from scipy.ndimage import binary_fill_holes
    ring = disc_mask(50, 50, 25, 25, 14, 9) | disc_mask(50, 50, 25, 25, 3)
    blob = disc_mask(50, 50, 6, 44, 3)
    text = build_target({"nodule": ring | blob}, SCHEMAS["lidc"], None, {"nodule": 5}).text
    assert len(parse(text, SCHEMAS["lidc"]).objs["nodule"]) == 3
    parts = parse(D._no_holes(text, SCHEMAS["lidc"]), SCHEMAS["lidc"]).objs["nodule"]
    g = to_mask(parts, 50, 50)
    assert g.legal and len(parts) == 2 and np.array_equal(g.mask, binary_fill_holes(ring | blob))


def test_view_rows_check_every_context(view_dir, tmp_path, fake_augment):
    """BK1-3: each context file must hash to the manifest's contexts[item_id] (checked once at load: 0.3-1.2 s for
    the 7,341 LIDC train contexts); a missing, edited or unlisted context is refused before training starts."""
    root, _ = view_dir
    v = tmp_path / "v"
    for edit in (lambda m, p: p.write_bytes(p.read_bytes() + b"\0"), lambda m, p: p.unlink(),
                 lambda m, p: m["contexts"].pop("riga/MESSIDOR/image001"), lambda m, p: m.pop("contexts")):
        shutil.rmtree(v, ignore_errors=True)
        shutil.copytree(root, v)
        m = json.loads((v / "manifest.json").read_text())
        edit(m, v / "context" / "riga" / "MESSIDOR" / "image001.npz")
        (v / "manifest.json").write_text(json.dumps(m))
        with pytest.raises(ValueError, match="context"):
            D.view_rows(v)


def test_view_rows_on_the_real_view_builder(tmp_path, synthetic_items, disc_mask):
    """Agent V's builder (M1) on the conftest items plus an item with a sub-filter speck rater: rows render through
    the real augment.render, only voting raters are taught (I27 / B2-4), non-simple rows are repeated, A0 keeps
    every row."""
    from objdel.data import Item
    from objdel.targets import augment
    build = pytest.importorskip("experiments.phase2.build_views")
    rng = np.random.default_rng(1)
    speck = np.zeros((100, 100), bool)
    speck[50:52, 50:52] = True  # 4 px < min_area 5: abstains from the topology vote
    item = Item("lidc/1.2.9/n000/z001", "lidc", "LIDC-IDRI-0009", "train",
                rng.integers(-1000, 400, (100, 100), dtype=np.int16),
                {"nodule": np.stack([disc_mask(100, 100, 50, 50, 12), disc_mask(100, 100, 50, 50, 11), speck])},
                ["1", "2", "3"], {})
    store.write(sorted([*synthetic_items(), item], key=lambda i: (i.dataset, i.split, i.item_id)), tmp_path / "s")
    m = build.build(tmp_path / "s", tmp_path / "v", "lidc", "train", procs=1)
    ds = D.view_rows(tmp_path / "v")
    assert len(ds) == m["counts"]["rows"] > 0 and not any(r["key"] == item.item_id and r["rater"] == "3" for r in ds)
    by = {(r["key"], r["rater"]) for r in ds}
    assert {(item.item_id, "1"), (item.item_id, "2"), ("lidc/1.2.1/n000/z001", "2")} <= by
    spec = json.loads((tmp_path / "v" / "rows.jsonl").read_text().splitlines()[0])
    want = augment.render(spec, augment.load_context(tmp_path / "v", spec["item_id"]))
    assert np.array_equal(np.asarray(ds[0]["image"]), want) and want.shape == (512, 512, 3)
    assert len(D.view_rows(tmp_path / "v", a0=True)) == len(ds)


# ------------------------------------------------------------------ configs (plan D4-D7, D13, §6.0)

# Linear modules of each model (one layer per repeated block), recorded offline from
# AutoModelForImageTextToText.from_config on the meta device (transformers 5.20.1, impl note M2).
LINEAR = {
    "qwen3-vl-2b": [
        "lm_head", "model.language_model.layers.0.mlp.down_proj", "model.language_model.layers.0.mlp.gate_proj",
        "model.language_model.layers.0.mlp.up_proj", "model.language_model.layers.0.self_attn.k_proj",
        "model.language_model.layers.0.self_attn.o_proj", "model.language_model.layers.0.self_attn.q_proj",
        "model.language_model.layers.0.self_attn.v_proj", "model.visual.blocks.0.attn.proj",
        "model.visual.blocks.0.attn.qkv", "model.visual.blocks.0.mlp.linear_fc1",
        "model.visual.blocks.0.mlp.linear_fc2",
        "model.visual.deepstack_merger_list.0.linear_fc1", "model.visual.deepstack_merger_list.0.linear_fc2",
        "model.visual.merger.linear_fc1", "model.visual.merger.linear_fc2"],
}
VISION_PREFIX = ("model.visual", "model.vision_tower", "model.multi_modal_projector", "model.vision_model",
                 "model.connector")


def toy(names):
    """A module tree with an nn.Linear(2, 2) at every recorded name (and a non-linear vision parameter)."""
    import torch
    root = torch.nn.Module()
    for name in names:
        *path, leaf = name.split(".")
        m = root
        for p in path:
            if not hasattr(m, p):
                m.add_module(p, torch.nn.Module())
            m = getattr(m, p)
        m.add_module(leaf, torch.nn.Linear(2, 2))
    return root


def test_models_and_budgets():
    assert configs.MODELS["qwen3-vl-2b"] == "Qwen/Qwen3-VL-2B-Instruct"
    assert configs.MAX_COMPLETION["qwen3-vl-2b"] == 1536


@pytest.mark.parametrize("model", list(LINEAR))
def test_sampling_is_explicit_and_card_free(model):
    s = configs.sampling(model)
    assert s == {"temperature": 1.0, "top_p": 1.0, "top_k": 0, "min_p": 0.0, "repetition_penalty": 1.0,
                 "presence_penalty": 0.0, "frequency_penalty": 0.0, "max_tokens": configs.MAX_COMPLETION[model]}


def test_sampling_names_for_grpo_and_hf_generate(tmp_path):
    """S1-10: one source, three spellings. GRPOConfig and GenerationConfig have no presence / frequency penalty
    (both must be 0), the length field is renamed, HF samples with do_sample and stops at the exported ids."""
    from transformers import GenerationConfig
    m = "qwen3-vl-2b"
    common = {"temperature": 1.0, "top_p": 1.0, "top_k": 0, "min_p": 0.0, "repetition_penalty": 1.0}
    assert configs.sampling(m, "grpo") == common | {"max_completion_length": 1536}
    assert configs.sampling(m, "hf") == common | {"max_new_tokens": 1536, "do_sample": True}
    configs.grpo_config(m, str(tmp_path), use_cpu=True, **configs.sampling(m, "grpo"))  # accepted names
    exported = configs.sampling(m) | {"seed": 3, "stop_token_ids": [7]}  # a prompts JSONL "sampling" field
    hf = configs.as_names(exported, "hf")
    assert hf == common | {"max_new_tokens": 1536, "do_sample": True, "eos_token_id": [7]}
    GenerationConfig(**hf).validate(strict=True)
    with pytest.raises(ValueError, match="presence_penalty"):
        configs.as_names(exported | {"presence_penalty": 1.5}, "grpo")


@pytest.mark.parametrize("model", list(LINEAR))
def test_lora_targets_the_language_model_only(model):
    """PEFT's own injection on the recorded module names: no vision tower or projector module gets an adapter,
    every language-model linear layer does (B1 §2.5: the all-linear adapter collapsed vision)."""
    from peft import inject_adapter_in_model
    cfg = configs.lora("qwen3-vl-2b")  # the adapter config is the same for any model
    assert (cfg.r, cfg.lora_alpha) == (16, 32)
    names = [n for n in LINEAR[model] if n != "lm_head"]
    m = inject_adapter_in_model(cfg.__class__(**{**cfg.to_dict(), "target_modules": sorted(
        {n.rsplit(".", 1)[-1] for n in names})}), toy(LINEAR[model]))  # all-linear = every linear leaf name
    adapted = {n.removesuffix(".lora_A") for n, _ in m.named_modules() if n.endswith(".lora_A")}
    assert adapted == {n for n in names if not n.startswith(VISION_PREFIX)}
    assert all(re.fullmatch(cfg.exclude_modules, n) for n in names if n.startswith(VISION_PREFIX))


def test_every_model_is_lora():
    assert configs.lora("qwen3-vl-2b") is not None and configs.lora("qwen3-vl-2b").target_modules == "all-linear"
    assert configs._LORA == frozenset(configs.MODELS)  # ADR-0008: no full fine-tuning


@pytest.mark.parametrize("model", list(LINEAR))
def test_sft_config_bundle(model, tmp_path):
    c = configs.sft_config(model, str(tmp_path))
    assert c.lr_scheduler_type == "warmup_stable_decay" and c.lr_scheduler_kwargs["decay_type"] == "linear"
    assert c.train_sampling_strategy == "group_by_length" and c.length_column_name == "length"
    assert (c.dataloader_num_workers, c.dataloader_persistent_workers, c.dataloader_prefetch_factor) == (6, True, 2)
    assert c.gradient_checkpointing == (model == "qwen3-vl-2b")
    assert c.assistant_only_loss and c.max_length is None and not c.remove_unused_columns
    assert c.learning_rate == 5e-4  # LoRA
    assert configs.sft_config(model, str(tmp_path), learning_rate=1e-4).learning_rate == 1e-4  # override
    c = configs.sft_config(model, str(tmp_path), dataloader_num_workers=0)  # BK3-1: R0's workers-0 timing arms
    assert (c.dataloader_num_workers, c.dataloader_persistent_workers, c.dataloader_prefetch_factor) == (0, False,
                                                                                                         None)


def lrs(c, n):
    import torch
    from transformers import get_scheduler
    opt = torch.optim.SGD([torch.nn.Parameter(torch.zeros(1))], lr=c.learning_rate)
    s = get_scheduler(c.lr_scheduler_type, opt, num_warmup_steps=c.get_warmup_steps(c.max_steps),
                      num_training_steps=c.max_steps, scheduler_specific_kwargs=c.lr_scheduler_kwargs)
    out = []
    for _ in range(n):
        out.append(opt.param_groups[0]["lr"])
        opt.step()
        s.step()
    return out


def test_wsd_branches_follow_the_stable_run_then_decay_to_zero(tmp_path):
    """D5: 3 % warmup, constant LR, a branch at step s decays linearly to 0 over 15 % of s; up to s the branch's
    schedule is the stable run's, so a branch resumes the stable checkpoint without a jump."""
    stable = configs.sft_config("qwen3-vl-2b", str(tmp_path), steps=1000)
    branch = configs.sft_config("qwen3-vl-2b", str(tmp_path), steps=1000, branch_at=500)
    assert (stable.max_steps, branch.max_steps) == (1000, 575)
    a, b = lrs(stable, 1000), lrs(branch, 575)
    assert a[:500] == b[:500] and a[0] == 0 and a[30] == a[999] == 5e-4 and max(a) == 5e-4
    assert b[500] == 5e-4 and b[574] == pytest.approx(5e-4 / 75) and all(np.diff(b[500:]) < 0)
    with pytest.raises(ValueError):
        configs.sft_config("qwen3-vl-2b", str(tmp_path), branch_at=500)  # a branch needs the stable length


def test_qwen3_vl_trains_with_trls_own_template():
    """No template swap of ours: TRL maps the processor's Qwen3-VL template to qwen3_vl_training (assistant-only
    mask, stop token trained). Needs the Qwen3-VL-2B processor in the local HF cache."""
    from transformers import AutoProcessor
    from trl import chat_template_utils as ct

    try:
        proc = AutoProcessor.from_pretrained(configs.MODELS["qwen3-vl-2b"], local_files_only=True)
    except OSError:
        pytest.skip("Qwen3-VL-2B processor not in the local HF cache")
    train = ct.get_training_chat_template(proc)
    assert train == ct.qwen3_vl_training_chat_template and ct.has_generation_markers(train)
    assert ct.is_chat_template_stop_token_trained(proc, chat_template=train)


# ------------------------------------------------------------------ M9 (I45): prepare and cache_loader


@pytest.fixture(scope="module")
def prepared(tmp_path_factory, synthetic_items):
    root = tmp_path_factory.mktemp("prep9")
    store.write(synthetic_items(), root / "store")
    D.prepare(root / "store", root / "out", splits=("train",), procs=1)
    return root / "out"


def test_prepare_refuses_a_failed_image_write(tmp_path, synthetic_items, monkeypatch):
    """One prepare job in this process (prepare's Pool workers would not see the patch)."""
    from objdel.targets import frame
    it = next(i for i in synthetic_items() if i.dataset == "lidc" and i.image is not None)
    monkeypatch.setattr(D.cv2, "imwrite", lambda *a, **k: False)
    with pytest.raises(OSError, match="imwrite"):
        D._prepare((str(tmp_path), it.item_id, it.dataset, it.group_id, it.raters, frame(it)))


def test_cache_loader_reads_only_the_manifests_files(prepared, tmp_path):
    out = tmp_path / "out"
    shutil.copytree(prepared, out)
    (out / "lidc" / "stray.jsonl").write_text('{"item_id": "x", "rl": true, "caches": {"nodule": "0"}}\n')
    load = verify_store.cache_loader(out)
    with pytest.raises(KeyError):
        load("x")  # not in the manifest: never read
    assert load("lidc/1.2.1/n000/z001")["nodule"].k == 2
    with open(out / "lidc" / "train.jsonl", "a") as f:
        f.write("\n")
    with pytest.raises(ValueError, match="manifest"):
        verify_store.cache_loader(out)


def test_cache_loader_unbounded_preload(prepared, tmp_path):
    out = tmp_path / "out"
    shutil.copytree(prepared, out)
    load = verify_store.cache_loader(out, maxsize=None)
    shutil.rmtree(out / "caches")  # everything was loaded up front
    assert load("lidc/1.2.1/n000/z001")["nodule"].k == 2 and set(load("riga/MESSIDOR/image001")) == {"disc", "cup"}


def test_prompt_rows_with_and_without_the_rl_filter(prepared):
    """Public prompt rows from the records (select_ckpt, probes): rl_only keeps RL-eligible items (grpo_rows without
    oversampling); rl_only=False keeps every item with pixels (probe records are rl: false)."""
    recs = D.records(prepared, "lidc", "train")
    every, rl = D.prompt_rows(prepared, "lidc", "train", rl_only=False), D.prompt_rows(prepared, "lidc", "train")
    assert every["key"] == [r["item_id"] for r in recs] and rl["key"] == [r["item_id"] for r in recs if r["rl"]]
    assert len(rl) < len(every) and rl.to_list() == D.grpo_rows(prepared, "lidc", "train").to_list()
    assert every[0]["prompt"][0]["content"][1]["text"] == D.PROMPTS["lidc"] and every[0]["image"].mode == "RGB"
