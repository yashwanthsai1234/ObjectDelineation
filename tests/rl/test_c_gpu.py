"""Phase C box tests (plan C.9.4). G-RC0 and G-GR1 load the start checkpoint `$OBJDEL_RUN_DIR/ckpt` (HF revision
`b1-rows`, fetched by the plan's item 0) on the pinned vLLM and run by node id in their own plan steps; U-RC6 checks
on the Mac that the same model's config, tokenizer and processor load (marker `data`: network and an HF token).
`gpu` tests skip on the CPU build; a missing checkpoint fails them, never skips."""

import contextlib
import os
from pathlib import Path

import pytest

from objdel.geom import engine as _gpucheck
from objdel.infra.hf import REPO

CUDA = pytest.mark.skipif(_gpucheck.build() != "cuda", reason="CPU build: the CUDA build runs on the box")
REVISION = "b1-rows"
FILES = ("config.json", "generation_config.json", "tokenizer.json", "tokenizer_config.json", "processor_config.json",
         "chat_template.jinja")  # what U-RC6 downloads: no weights
PROMPT = "Describe a nodule in one sentence."


def check_model_files(path):
    """U-RC6: the folder's config, tokenizer and processor load, and the chat template renders one text prompt."""
    from transformers import AutoConfig, AutoProcessor, AutoTokenizer

    cfg = AutoConfig.from_pretrained(path)
    assert cfg.architectures == ["Qwen3VLForConditionalGeneration"]
    assert cfg.text_config.max_position_embeddings >= 256  # the engine's max_model_len
    tokenizer = AutoTokenizer.from_pretrained(path)
    assert tokenizer("Outline of object A:")["input_ids"] and tokenizer.pad_token_id is not None
    processor = AutoProcessor.from_pretrained(path)
    assert getattr(processor, "image_processor", None) is not None and processor.image_processor.patch_size == 16
    chat = [{"role": "user", "content": [{"type": "text", "text": PROMPT}]}]
    assert PROMPT in processor.apply_chat_template(chat, tokenize=False, add_generation_prompt=True)


@pytest.mark.data
def test_vllm_check_model_loads(tmp_path):
    """U-RC6: the model of G-RC0 and G-GR1 (HF revision b1-rows) loads on the Mac from its small files only."""
    from huggingface_hub import get_token, snapshot_download

    token = get_token()
    if token is None:
        pytest.skip("no HF token set: this check needs the owner's token for the private model repo")
    snapshot_download(REPO, revision=REVISION, allow_patterns=list(FILES), local_dir=tmp_path, token=token)
    check_model_files(tmp_path)


def engine():
    """The pinned vLLM on the start checkpoint, text-only prompts (plan C.9.4 model)."""
    from vllm import LLM

    ckpt = Path(os.environ["OBJDEL_RUN_DIR"]) / "ckpt"
    assert (ckpt / "config.json").is_file(), f"{ckpt}: no start checkpoint (fetch_hf.sh <run> {REVISION} puts it there)"
    return LLM(model=str(ckpt), max_model_len=256, gpu_memory_utilization=0.3, enforce_eager=True)


@pytest.mark.gpu
@CUDA
def test_g_rc0_per_row_caps():
    """G-RC0: a SamplingParams list gives each answer its own max_tokens (ignore_eos: every answer reaches its cap,
    so a shared cap would show)."""
    from vllm import SamplingParams

    caps = [4, 8, 16]
    params = [SamplingParams(max_tokens=c, ignore_eos=True, temperature=1.0) for c in caps]
    out = engine().generate([PROMPT] * len(caps), sampling_params=params)
    assert [len(o.outputs[0].token_ids) for o in out] == caps


@pytest.mark.gpu
@CUDA
def test_g_gr1_requests_independent():
    """G-GR1: one request with n = 8 per prompt gives independent samples: at temperature 1.0 a prompt's 8 answers
    average at least 2 distinct texts over 32 prompts."""
    from vllm import SamplingParams

    prompts = [f"{PROMPT} Case {i}." for i in range(32)]
    out = engine().generate(prompts, sampling_params=SamplingParams(n=8, temperature=1.0, max_tokens=32))
    assert all(len(o.outputs) == 8 for o in out)
    assert sum(len({c.text for c in o.outputs}) for o in out) / len(out) >= 2


# --- C1, C2, C4 box tests on synthetic inputs (plan C.9.4: V2 item 2, pytest) and the recorded-batch replay ---

H2D_BOUND = 32 * 512 * 512 * 3 + 2**20  # G-PX4: the 32 distinct uint8 images once, plus 1 MB


def processor():
    from transformers import Qwen2VLImageProcessor

    return Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)  # Qwen3-VL's patch and merge


def pil_images(n=32, size=512, seed=0):
    """n distinct seeded random RGB images (512 x 512: the training views' size)."""
    import numpy as np
    import PIL.Image

    rng = np.random.default_rng(seed)
    return [PIL.Image.fromarray(rng.integers(0, 256, (size, size, 3), dtype=np.uint8)) for _ in range(n)]


def gpu_memo(backend="torch", per_image=False):
    from objdel.policy.caches import ImageMemo

    return ImageMemo(processor(), gpu=(backend, "cuda"), per_image=per_image)


def rows_of(images, G=8):
    """A generation batch's pass-2 rows: each image object G times in a row, as OpenOnce hands them."""
    return [x for x in images for _ in range(G)]


def no_syncs(fn):
    """(fn's result, syncs counted around it); fn also runs under set_sync_debug_mode("error") (raises on a sync)."""
    import torch

    from objdel.policy.marks import count_syncs

    with count_syncs() as n:
        out = fn()
    prev = torch.cuda.get_sync_debug_mode()
    torch.cuda.set_sync_debug_mode("error")
    try:
        fn()
    finally:
        torch.cuda.set_sync_debug_mode(prev)
    return out, n.syncs


def median_ms(fn, n=50):
    """Median wall time of fn in ms over n warm runs, the device synchronized around each."""
    import statistics
    import time

    import torch

    fn()
    times = []
    for _ in range(n):
        torch.cuda.synchronize()
        t = time.perf_counter()
        fn()
        torch.cuda.synchronize()
        times.append((time.perf_counter() - t) * 1e3)
    return statistics.median(times)


@pytest.mark.parametrize("backend", ["hf", "torch"])
@pytest.mark.gpu
@CUDA
def test_g_px2_patchify_and_pass_two_make_no_host_sync(backend):
    """G-PX2: after warm-up, patchify of 32 uint8 images onto CUDA and one pass-2 ImageMemo call of 256 rows after
    prefill make no host sync (count_syncs 0, and set_sync_debug_mode("error") raises on none); the grid is a CPU
    tensor built from the shapes."""
    from objdel.policy import gpu_pixels

    ip, images = processor(), pil_images()
    gpu_pixels.patchify(images, ip, device="cuda", backend=backend)  # warm-up
    (pv, grid), syncs = no_syncs(lambda: gpu_pixels.patchify(images, ip, device="cuda", backend=backend))
    assert syncs == 0 and pv.is_cuda and grid.device.type == "cpu"
    memo, rows = gpu_memo(backend), rows_of(images)
    memo.prefill(images, return_tensors="pt")
    memo(rows, return_tensors="pt")  # warm-up of the pass-2 assembly
    out, syncs = no_syncs(lambda: memo(rows, return_tensors="pt"))
    assert syncs == 0 and out["image_grid_thw"].device.type == "cpu"


@pytest.mark.gpu
@CUDA
def test_g_px3_pixel_speed(record_property):
    """G-PX3 (soft gate): 32 images of 512 x 512: patchify (torch, cuda) under 5 ms median of 50; the pass-2 ImageMemo
    call under 10 ms for 32 rows and under 20 ms for 256 rows, warm; at least 5x faster than the CPU processor (its
    time is recorded as the baseline)."""
    if os.environ.get("CI"):
        pytest.skip("soft gate: timings on a shared CI host are not the box's")
    from objdel.policy import gpu_pixels

    ip, images = processor(), pil_images()
    gpu = median_ms(lambda: gpu_pixels.patchify(images, ip, device="cuda", backend="torch"))
    cpu = median_ms(lambda: ip(images, return_tensors="pt"), n=5)
    memo = gpu_memo()
    memo.prefill(images, return_tensors="pt")
    rows32, rows256 = median_ms(lambda: memo(images, return_tensors="pt")), median_ms(
        lambda: memo(rows_of(images), return_tensors="pt"))  # fmt: skip
    for k, v in {"patchify_ms": gpu, "cpu_processor_ms": cpu, "memo_32_ms": rows32, "memo_256_ms": rows256}.items():
        record_property(k, v)
    assert gpu < 5 and rows32 < 10 and rows256 < 20 and 5 * gpu <= cpu


@pytest.fixture
def h2d(monkeypatch):
    """Bytes copied host -> CUDA through Tensor.to, Tensor.cuda and Tensor.copy_ while the fixture is active (tensors
    built directly on the device, such as the processor's two (3,) mean/std vectors, are not counted)."""
    import torch

    seen = {"bytes": 0}
    real = {n: getattr(torch.Tensor, n) for n in ("to", "cuda", "copy_")}

    def moved(name):
        def spy(self, *a, **k):
            out = real[name](self, *a, **k)
            if name == "copy_":
                src = a[0] if a else k["src"]
                if isinstance(src, torch.Tensor) and src.device.type == "cpu" and self.is_cuda:
                    seen["bytes"] += src.nbytes
            elif self.device.type == "cpu" and isinstance(out, torch.Tensor) and out.is_cuda:
                seen["bytes"] += self.nbytes
            return out

        return spy

    for name in real:
        monkeypatch.setattr(torch.Tensor, name, moved(name))
    return seen


@pytest.mark.parametrize("backend", ["hf", "torch"])
@pytest.mark.gpu
@CUDA
def test_g_px4_pixel_pass_uploads_each_image_once(backend, h2d):
    """G-PX4: the pixel pass of a generation batch (prefill of its 32 images, then the 256-row pass-2 call) copies at
    most 32 x 512 x 512 x 3 + 1 MB = 25.2 MB host -> device: the uint8 images once, no 256-row upload."""
    images = pil_images()
    memo = gpu_memo(backend)
    memo.prefill(images, return_tensors="pt")
    out = memo(rows_of(images), return_tensors="pt")
    assert out["pixel_values"].is_cuda and 0 < h2d["bytes"] <= H2D_BOUND


@pytest.mark.gpu
@CUDA
def test_g_px5_out_of_memory_leaves_no_damage():
    """G-PX5: a patchify too large for the device (4,096 images of 512 x 512 or more, so the float32 rows exceed its
    memory) raises OutOfMemoryError; the next ImageMemo call equals the processor; after gc.collect() and
    empty_cache() the allocated memory is back at its level before the call (+ 1 MB)."""
    import gc
    import math

    import torch

    from objdel.policy import gpu_pixels

    ip, images = processor(), pil_images(8)
    n = max(4096, math.ceil(1.2 * torch.cuda.get_device_properties(0).total_memory / (1024 * 1536 * 4)))
    gc.collect()
    torch.cuda.empty_cache()
    before = torch.cuda.memory_allocated()
    raised = False
    try:
        gpu_pixels.patchify([images[0]] * n, ip, device="cuda", backend="torch")
    except torch.OutOfMemoryError:
        raised = True  # the traceback (and the frames' tensors) is dropped when this block ends
    assert raised, f"{n} images did not exhaust the device"
    memo = gpu_memo()
    out = memo(rows_of(images, 2), return_tensors="pt")
    want = ip(rows_of(images, 2), return_tensors="pt")
    assert torch.equal(out["pixel_values"].cpu(), want["pixel_values"])
    assert torch.equal(out["image_grid_thw"], want["image_grid_thw"])
    del memo, out
    gc.collect()
    torch.cuda.empty_cache()
    assert torch.cuda.memory_allocated() <= before + 2**20


@pytest.mark.parametrize("per_image", [False, True])
@pytest.mark.gpu
@CUDA
def test_g_im1_memo_pixels_on_cuda_grid_on_cpu(per_image):
    """G-IM1: the memo's pixel_values are on CUDA and its image_grid_thw on the CPU, in both per_image modes; the
    values equal the processor after .cpu() (per_image: each distinct image once, first-appearance order)."""
    import torch

    ip, images = processor(), pil_images(4)
    rows = [images[i] for i in (0, 0, 1, 2, 0, 1, 3, 3)]
    out = gpu_memo(per_image=per_image)(rows, return_tensors="pt")
    want = ip(images if per_image else rows, return_tensors="pt")
    assert out["pixel_values"].is_cuda and out["image_grid_thw"].device.type == "cpu"
    assert torch.equal(out["pixel_values"].cpu(), want["pixel_values"])
    assert torch.equal(out["image_grid_thw"], ip(rows, return_tensors="pt")["image_grid_thw"])


@pytest.mark.xfail(strict=True, raises=AssertionError,
                   reason="known gap: ImageMemo.prepare's fallback returns the CPU processor's pixel_values")
@pytest.mark.gpu
@CUDA
def test_g_im1_fallback_pixels_on_cuda():
    """G-IM1, fallback path: a call holding an image patchify cannot match (500 x 500) falls back to the processor
    for the whole call; its pixel_values should still be on CUDA and equal the processor."""
    import numpy as np
    import PIL.Image
    import torch

    ip = processor()
    a = pil_images(1)[0]
    bad = PIL.Image.fromarray(np.random.default_rng(1).integers(0, 256, (500, 500, 3), dtype=np.uint8))
    memo = gpu_memo()
    out = memo([a, bad, a], return_tensors="pt")
    assert memo.fallback_calls == 1
    assert torch.equal(out["pixel_values"].cpu(), ip([a, bad, a], return_tensors="pt")["pixel_values"])
    assert out["pixel_values"].is_cuda


def rope_case(device, B=256, answer=1536, seed=0):
    """The LIDC prompt layout on the tiny Qwen3-VL (Lp 325, grid (1, 32, 32): 256 image pads) and B rows of it with
    left pads and answers of up to `answer` text tokens: (template on `device`, ids, mask, mm types on `device`)."""
    import torch

    from objdel.policy import caches
    from tests.experiments.test_trl_adapter_equivalence import END, IMG, START, tiny_model

    g = torch.Generator().manual_seed(seed)
    prompt = torch.tensor([5, 6, START, *[IMG] * 256, END, *torch.randint(0, 99, (325 - 260,), generator=g).tolist()])
    W = 8 + len(prompt) + answer
    ids, mask = torch.zeros(B, W, dtype=torch.long), torch.zeros(B, W, dtype=torch.long)
    for r, n in enumerate(torch.randint(0, answer + 1, (B,), generator=g).tolist()):
        p = r % 8  # left pads 0..7
        ids[r, p:p + len(prompt)] = prompt
        ids[r, p + len(prompt):p + len(prompt) + n] = torch.randint(0, 99, (n,), generator=g)
        mask[r, p:p + len(prompt) + n] = 1
    mm = (ids == IMG).long() * mask
    tmpl = caches.rope_template(tiny_model().model, prompt.to(device), (prompt == IMG).long().to(device),
                                torch.tensor([[1, 32, 32]]))  # fmt: skip
    return tmpl, ids.to(device), mask.to(device), mm.to(device)


@pytest.mark.gpu
@CUDA
def test_g_rp1_rope_positions_on_the_device_without_a_sync():
    """G-RP1: rope_positions on CUDA rows (template on CUDA, as rope_template places it) makes no host sync
    (count_syncs 0, set_sync_debug_mode("error") raises nothing) and takes under 5 ms for 256 rows."""
    from objdel.policy import caches

    tmpl, ids, mask, mm = rope_case("cuda")
    (pos, bad), syncs = no_syncs(lambda: caches.rope_positions(tmpl, ids, mask, mm))
    assert syncs == 0 and pos.is_cuda and bad.is_cuda
    assert median_ms(lambda: caches.rope_positions(tmpl, ids, mask, mm)) < 5


@pytest.mark.xfail(strict=True, raises=RuntimeError,
                   reason="known gap: rope_positions indexes the template with the rows' index tensor, so the two "
                          "must share a device (U-RP5 asks for both mixes to work)")  # fmt: skip
@pytest.mark.parametrize("where", ["cpu template", "cuda template"])
@pytest.mark.gpu
@CUDA
def test_g_rp1_template_and_rows_on_different_devices(where):
    """G-RP1: a CPU template with CUDA rows, and the reverse, give the same positions without a sync."""
    import torch

    from objdel.policy import caches

    tmpl, ids, mask, mm = rope_case("cuda", B=8, answer=64)
    want = caches.rope_positions(tmpl, ids, mask, mm)[0].cpu()
    if where == "cpu template":
        tmpl = caches.RopeTemplate(tmpl.pos.cpu(), tmpl.end, tmpl.ids.cpu(), tmpl.types.cpu())
    else:
        ids, mask, mm = ids.cpu(), mask.cpu(), mm.cpu()
    (pos, _), syncs = no_syncs(lambda: caches.rope_positions(tmpl, ids, mask, mm))
    assert syncs == 0 and torch.equal(pos.cpu(), want)


def load_records(folder: Path) -> list[dict]:
    """The --record-batches files (step*.pt) or c_speed_check batches (b*.pt) of a folder, in name order."""
    import torch

    return [torch.load(p, weights_only=True) for p in sorted(folder.glob("*.pt"))]


def meta_of_record(rec):
    import numpy as np

    from objdel.policy.caches import RowMeta

    return RowMeta(**{k: np.asarray(v, dtype=np.int32) for k, v in rec["meta"].items()})


def host_calls(rec, meta, chunk=16) -> dict:
    """The C4 host functions on one recorded batch, each a no-argument call: row_meta (from the batch's id lists),
    window_from_meta and meta_runs of its first chunk, group_order and length_order."""
    from objdel.rl import forward

    P, T = rec["prompt_ids"].size(1), rec["completion_ids"].size(1)
    G = int((meta.group == meta.group[0]).sum())
    pm, cm = rec["prompt_mask"].bool(), rec["completion_mask"]
    prompts = [r[m].tolist() for r, m in zip(rec["prompt_ids"], pm, strict=True)]
    answers = [r[:int(n)].tolist() if n else r.tolist() for r, n in zip(rec["completion_ids"], cm.sum(1), strict=True)]
    ids = meta.image_id.tolist()
    return {"row_meta": lambda: forward.row_meta(prompts, answers, ids, G, False, (0,)),
            "window_from_meta": lambda: forward.window_from_meta(meta, slice(0, chunk), P, T),
            "meta_runs": lambda: forward.meta_runs(meta, slice(0, chunk)),
            "group_order": lambda: forward.group_order(meta, G),
            "length_order": lambda: forward.length_order(meta)}  # fmt: skip


@pytest.mark.gpu
@CUDA
def test_g_rm1_host_functions_make_no_sync_and_are_fast():
    """G-RM1 (sync and time part): on the 5 c_speed_check batches ($OBJDEL_RUN_DIR/c_check/batches) the C4 host
    functions make no host sync (count_syncs 0 across all of them) and each takes under 1 ms (row_meta under 3 ms)."""
    import statistics
    import time

    from objdel.policy.marks import count_syncs

    run = os.environ.get("OBJDEL_RUN_DIR")
    folder = Path(run or ".") / "c_check" / "batches"
    if not run or not any(folder.glob("b*.pt")):
        pytest.skip(f"no c_speed_check batches at {folder} (c_speed_check --pixels writes them)")
    for rec in load_records(folder):
        calls = host_calls(rec, meta_of_record(rec))
        with count_syncs() as n:
            for fn in calls.values():
                fn()
        assert n.syncs == 0
        for name, fn in calls.items():
            times = []
            for _ in range(20):
                t = time.perf_counter()
                fn()
                times.append((time.perf_counter() - t) * 1e3)
            assert statistics.median(times) < (3 if name == "row_meta" else 1), name


def v2_batches() -> list[dict]:
    """G-RP2 / G-RM1 equality inputs: the pulled V2 run's 5 c_speed_check batches and the adopted fast arm's check
    batches (arm_batches_b when item 3c ran, else arm_batches); skips without OBJDEL_V2_RUN, fails under 15."""
    name = os.environ.get("OBJDEL_V2_RUN")
    if not name:
        pytest.skip("no V2 run pulled (OBJDEL_V2_RUN unset)")
    check = Path(__file__).resolve().parents[2] / "runs" / name / "c_check"
    arm = check / "arm_batches_b" if (check / "arm_batches_b").is_dir() else check / "arm_batches"
    recs = load_records(check / "batches") + load_records(arm)
    assert len(recs) >= 15, f"{len(recs)} box batches under {check} (needs at least 15)"
    return recs


@pytest.mark.data
def test_g_rm1_host_meta_equals_the_device_derivation():
    """G-RM1 (equality part, G-EQ3): on every pulled box batch the host functions equal the device derivation they
    replace: each 16-row chunk's window, the length and group orders, the shared-prompt runs, and the lengths and
    pixel rows that packed_chunk's seg and idx are built from (check_meta)."""
    import numpy as np
    import torch

    from objdel.rl import forward

    for rec in v2_batches():
        meta = meta_of_record(rec)
        P, T = rec["prompt_ids"].size(1), rec["completion_ids"].size(1)
        G = int((meta.group == meta.group[0]).sum())
        am = torch.cat([rec["prompt_mask"], rec["completion_mask"]], 1)
        for a in range(0, len(am), 16):
            rows = slice(a, a + 16)
            assert forward.window(am[rows], T) == forward.window_from_meta(meta, rows, P, T)
            assert forward.meta_runs(meta, rows) == forward.prompt_groups(
                rec["prompt_ids"][rows], rec["prompt_mask"][rows], None, image_ids=meta.image_id[rows])
        n = am.sum(1)
        assert np.array_equal(forward.length_order(meta), torch.sort(n, descending=True, stable=True).indices.numpy())
        groups = torch.sort(n.view(-1, G).max(1).values, descending=True, stable=True).indices
        assert np.array_equal(forward.group_order(meta, G), (groups[:, None] * G + torch.arange(G)).flatten().numpy())
        forward.check_meta({"host_meta": forward.host_rows(meta), "prompt_mask": rec["prompt_mask"],
                            "completion_mask": rec["completion_mask"], "image_grid_thw": rec.get("image_grid_thw"),
                            "num_images": rec["num_images"].tolist()})  # fmt: skip


@pytest.mark.gpu
@CUDA
def test_g_rm2_packed_staging_survives_a_delayed_copy(monkeypatch):
    """G-RM2: two packed chunks back to back while each chunk's index copy waits behind a long kernel
    (torch.cuda._sleep queued just before it): the indices on the device equal the numpy plan each chunk staged."""
    import torch

    from tests.experiments.test_trl_adapter_equivalence import group_batch, shared, tiny_model

    m = tiny_model().cuda()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[2, 0])
    ids, mask = ids.cuda(), mask.cuda()
    kw = {k: v.cuda() if isinstance(v, torch.Tensor) else v for k, v in kw.items()}
    stages, seen = [], []
    real_empty, real_select = torch.empty, torch.index_select

    def empty(*a, pin_memory=False, **k):
        out = real_empty(*a, pin_memory=pin_memory, **k)
        if pin_memory and k.get("dtype") is torch.int64:
            stages.append(out)
            torch.cuda._sleep(2_000_000_000)  # holds the stream: this chunk's copy queues behind it
        return out

    def index_select(x, dim, index, *a, **k):
        if stages and (not seen or seen[-1][0] is not stages[-1]):  # the chunk's first gather: its src indices
            seen.append((stages[-1], stages[-1].numpy().copy(), index.clone()))
        return real_select(x, dim, index, *a, **k)

    monkeypatch.setattr(torch, "empty", empty)
    monkeypatch.setattr(torch, "index_select", index_select)
    shared(m)._get_per_token_logps_and_entropies(m, ids, mask, T, batch_size=4, **kw)
    torch.cuda.synchronize()
    assert len(seen) == 2
    for _, plan, src in seen:
        cw = src.numel()  # the stage holds seg (c x W), then src (c x W), then idx
        assert torch.equal(src.cpu(), torch.from_numpy(plan[cw:2 * cw]))


@pytest.mark.gpu
@CUDA
def test_g_dt1_repeats_streams_and_cudnn_benchmark_give_identical_bytes():
    """G-DT1: 10 repeats of patchify, rope_positions and window_from_meta give identical bytes; another stream and
    either cudnn.benchmark setting change nothing."""
    import numpy as np
    import torch

    from objdel.policy import caches, gpu_pixels
    from objdel.rl import forward

    ip, images = processor(), pil_images(8)
    tmpl, ids, mask, mm = rope_case("cuda", B=32, answer=256)
    lens = (mask.sum(1) - 325).tolist()  # the prompt is 325 real tokens; P = 8 pads + 325 columns
    meta = caches.RowMeta(*(np.asarray(x, dtype=np.int32) for x in (
        [325] * 32, lens, range(32), [i // 8 for i in range(32)], [0] * 32, [0] * 32, [1024] * 32)))  # fmt: skip

    def run():
        torch.cuda.synchronize()
        pv, grid = gpu_pixels.patchify(images, ip, device="cuda", backend="torch")
        pos, bad = caches.rope_positions(tmpl, ids, mask, mm)
        win = [forward.window_from_meta(meta, slice(a, a + 8), 333, ids.size(1) - 333) for a in range(0, 32, 8)]
        torch.cuda.synchronize()
        return [pv.cpu(), grid, pos.cpu(), bad.cpu()], win

    ref, win = run()
    bench = torch.backends.cudnn.benchmark
    try:
        for i in range(12):
            torch.backends.cudnn.benchmark = i == 11
            with torch.cuda.stream(torch.cuda.Stream()) if i >= 10 else contextlib.nullcontext():
                got, w = run()
            assert w == win and all(torch.equal(a, b) for a, b in zip(got, ref, strict=True)), f"run {i}"
    finally:
        torch.backends.cudnn.benchmark = bench


@pytest.mark.gpu
@CUDA
def test_g_ms1_count_syncs_counts_each_sync_exactly():
    """G-MS1 (count part): 1,000 .item() calls on one line count 1,000; .cpu() 1; .tolist() 1; an async kernel 0;
    torch.cuda.synchronize() 0 (pinned: torch's sync-debug mode does not report an explicit synchronize); the
    previous mode (default 0, error 2) comes back, also after an exception; the user's warning filters (here
    "ignore") neither suppress the count nor change."""
    import warnings

    import torch

    from objdel.policy.marks import count_syncs

    x = torch.ones(4, device="cuda")

    def counted(fn):
        with count_syncs() as n:
            fn()
        return n.syncs

    assert counted(lambda: [x[0].item() for _ in range(1000)]) == 1000
    assert (counted(lambda: x.cpu()), counted(lambda: x.tolist()), counted(lambda: x * 2)) == (1, 1, 0)
    assert counted(torch.cuda.synchronize) == 0
    prev = torch.cuda.get_sync_debug_mode()
    try:
        for mode in (0, 2):
            torch.cuda.set_sync_debug_mode(mode)
            counted(lambda: x * 2)
            assert torch.cuda.get_sync_debug_mode() == mode
            with pytest.raises(KeyError), count_syncs():
                raise KeyError
            assert torch.cuda.get_sync_debug_mode() == mode
    finally:
        torch.cuda.set_sync_debug_mode(prev)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        before = list(warnings.filters)
        assert counted(lambda: x[0].item()) == 1
        assert warnings.filters == before


@pytest.mark.xfail(strict=True, raises=AssertionError,
                   reason="known gap: a nested count_syncs records its block's warnings in its own list, so the "
                          "outer block misses them (fake-CUDA check on the Mac: outer 1, inner 2)")  # fmt: skip
@pytest.mark.gpu
@CUDA
def test_g_ms1_nested_blocks_count_independently():
    """G-MS1: nested blocks count independently, and an outer block also counts the inner block's syncs."""
    import torch

    from objdel.policy.marks import count_syncs

    x = torch.ones(4, device="cuda")
    with count_syncs() as outer:
        x.cpu()
        with count_syncs() as inner:
            x.cpu()
            x.cpu()
    assert inner.syncs == 2 and outer.syncs == 3
