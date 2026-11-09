"""The model frame (plan §4.4 crop policy): square crop around all raters' union with margin, side <= 500 px and a
divisor of 1000, integer block-majority downscale for large objects, LIDC model image upscaled to 512 px."""
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from PIL import Image

from objdel.data import Item
from objdel.targets import FRAMES, FrameConfig, build_target, frame

DIVISORS = {d for d in range(1, 501) if 1000 % d == 0}


def item(masks, image=None, dataset="riga"):
    k = next(iter(masks.values())).shape[0]
    return Item(item_id="x", dataset=dataset, group_id="g", split="train", image=image, masks=masks,
                raters=[str(i) for i in range(k)], meta={})


def blob(shape, y0, x0, y1, x1):
    m = np.zeros(shape, bool)
    m[y0:y1, x0:x1] = True
    return m


def union_box(stack):
    ys, xs = np.nonzero(stack.any(0))
    return xs.min(), ys.min(), xs.max() + 1, ys.max() + 1


def check_frame(it, cfg):
    f = frame(it, cfg)
    x0, y0, x1, y1 = f.box
    assert f.W == f.H and f.W in DIVISORS and cfg.min_side <= f.W <= cfg.max_side
    assert x1 - x0 == y1 - y0 == f.scale * f.W
    assert all(type(v) is int for v in (*f.box, f.scale, f.W, f.H))  # JSON-native
    for n, s in it.masks.items():
        assert f.masks[n].shape == (s.shape[0], f.H, f.W) and f.masks[n].dtype == bool
    # the crop holds the union of all raters and structures with the margin (unless the image border cuts it)
    u = np.stack([s.any(0) for s in it.masks.values()])
    ux0, uy0, ux1, uy1 = union_box(u)
    side = max(ux1 - ux0, uy1 - uy0)
    m = max(cfg.margin_px, -(-side * cfg.margin_pct // 100))
    Hs, Ws = u.shape[1:]
    assert x0 <= max(0, ux0 - m) and x1 >= min(Ws, ux1 + m) and y0 <= max(0, uy0 - m) and y1 >= min(Hs, uy1 + m)
    return f


def test_crop_is_exact_slice_at_scale_1():
    s = np.stack([blob((400, 600), 100, 200, 160, 290), blob((400, 600), 110, 190, 150, 300)])
    img = np.random.default_rng(0).integers(0, 256, (400, 600, 3), dtype=np.uint8)
    f = check_frame(item({"disc": s, "cup": s & blob((400, 600), 120, 220, 140, 260)}, img), FRAMES["riga"])
    x0, y0, x1, y1 = f.box
    assert f.scale == 1 and 0 <= x0 and x1 <= 600 and 0 <= y0 and y1 <= 400  # shifted inside the image, no padding
    assert np.array_equal(f.masks["disc"], s[:, y0:y1, x0:x1])
    assert np.array_equal(f.image, img[y0:y1, x0:x1])
    # side: union 110 px + 2 x 28 px margin = 166 -> 200, the next divisor of 1000
    assert f.W == 200


def test_grid_maps_to_integer_pixels():
    """A target vertex g on the frame grid is the item pixel corner box[0] + scale * g * W / 1000, exactly."""
    s = np.stack([blob((300, 300), 120, 130, 141, 171)])
    f = frame(item({"nodule": s}, dataset="lidc"), FRAMES["lidc"])
    assert 1000 % f.W == 0
    t = build_target({"nodule": f.masks["nodule"][0]}, ("nodule",), None, {"nodule": 5})
    xy = t.objs["nodule"][0][0]
    assert np.all(xy * f.W % 1000 == 0)
    src = xy * f.W // 1000 * f.scale + np.array(f.box[:2])
    assert sorted(map(tuple, src.tolist())) == sorted([(130, 120), (171, 120), (171, 141), (130, 141)])


def test_padding_when_the_object_touches_the_image_border():
    s = np.stack([blob((60, 70), 0, 0, 30, 25)])
    img = np.full((60, 70), 400, np.int16)  # Hounsfield units
    f = check_frame(item({"nodule": s}, img, "lidc"), FRAMES["lidc"])
    assert f.W == 50 and f.box[:2] == (0, 0)  # 30 + 2 x 8 = 46 -> 50; fits the 60 x 70 image after shifting
    s2 = np.stack([blob((40, 40), 0, 0, 30, 25)])
    f2 = frame(item({"nodule": s2}, np.full((40, 40), 400, np.int16), "lidc"), FRAMES["lidc"])
    x0, y0, x1, y1 = f2.box
    assert x1 - x0 == 50 and (x0 < 0 or x1 > 40)  # image smaller than the crop: padded
    assert f2.masks["nodule"].sum() == s2.sum()


def test_lidc_image_is_windowed_rgb_upscaled_to_512_masks_not():
    s = np.stack([blob((512, 512), 200, 200, 230, 240)] * 3)
    hu = np.full((512, 512), -1000, np.int16)
    hu[200:230, 200:240] = 400
    f = frame(item({"nodule": s}, hu, "lidc"), FRAMES["lidc"])
    assert f.image.shape == (512, 512, 3) and f.image.dtype == np.uint8
    assert f.image.min() == 0 and f.image.max() == 255  # lung window [-1000, 400]
    assert f.masks["nodule"].shape == (3, f.H, f.W) and f.W < 512


def test_lidc_upscale_is_pillow_bicubic():
    """F3 (M9 finding 3): the upscale is Pillow's fixed-point bicubic, bit-identical across architectures (OpenCV's
    INTER_CUBIC differed by +-1 on x86 with and without IPP); input is the windowed uint8 crop at the frame side."""
    rng = np.random.default_rng(3)
    s = np.stack([blob((300, 300), 120, 130, 160, 170)])
    hu = rng.integers(-1200, 600, (300, 300)).astype(np.int16)
    f = frame(item({"nodule": s}, hu, "lidc"), FRAMES["lidc"])
    small = frame(item({"nodule": s}, hu, "lidc"), FrameConfig(image_side=None)).image  # same crop, not upscaled
    assert small.shape == (f.H, f.W, 3) and f.image.shape == (512, 512, 3)
    expect = np.asarray(Image.fromarray(small).resize((512, 512), Image.Resampling.BICUBIC))
    assert np.array_equal(f.image, expect)


def test_missing_image_gives_none():
    s = np.stack([blob((512, 512), 200, 200, 230, 240)])
    assert frame(item({"nodule": s}, None, "lidc"), FRAMES["lidc"]).image is None


def test_large_object_is_downscaled_by_block_majority():
    H, W = 1600, 1700
    yy, xx = np.mgrid[:H, :W]
    d = (yy - 800) ** 2 + (xx - 850) ** 2 <= 430 ** 2
    c = (yy - 790) ** 2 + (xx - 860) ** 2 <= 200 ** 2
    img = np.random.default_rng(1).integers(0, 256, (H, W, 3), dtype=np.uint8)
    f = check_frame(item({"disc": np.stack([d, d]), "cup": np.stack([c, c])}, img), FRAMES["riga"])
    # 861 px + 2 x 216 px = 1293 -> factor 3 -> 431 -> 500 px
    assert f.scale == 3 and f.W == 500
    x0, y0, x1, y1 = f.box
    assert x0 >= 0 and y0 >= 0 and x1 <= W and y1 <= H
    votes = d[y0:y1, x0:x1].reshape(500, 3, 500, 3).sum((1, 3))
    assert np.array_equal(f.masks["disc"][0], votes >= 5)  # majority of 9
    assert not (f.masks["cup"] & ~f.masks["disc"]).any()  # nesting survives the vote
    block = img[y0:y0 + 3, x0:x0 + 3].reshape(-1, 3).astype(np.int64)
    assert f.image.shape == (500, 500, 3) and np.array_equal(f.image[0, 0], (block.sum(0) + 4) // 9)  # mean


def test_even_factor_ties_go_to_foreground():
    s = np.zeros((1, 1100, 1100), bool)
    s[0, :1000, :1000] = True
    s[0, 0, 0:2] = False  # block (0, 0): 2 of 4 pixels -> tie -> foreground
    s[0, 2, 2:4] = s[0, 3, 2] = False  # block (1, 1): 1 of 4 -> background
    f = frame(item({"nodule": s}), FrameConfig(margin_px=0, margin_pct=0))
    assert f.scale == 2 and f.W == 500 and f.box[:2] == (0, 0)
    assert f.masks["nodule"][0, 0, 0] and not f.masks["nodule"][0, 1, 1]


def test_reproducible():
    rng = np.random.default_rng(2)
    s = rng.random((4, 300, 300)) < 0.001
    s[:, 100:150, 120:160] = True
    it = item({"nodule": s}, None, "lidc")
    a, b = frame(it, FRAMES["lidc"]), frame(it, FRAMES["lidc"])
    assert a.box == b.box and a.scale == b.scale and np.array_equal(a.masks["nodule"], b.masks["nodule"])


def test_empty_union_raises():
    with pytest.raises(ValueError):
        frame(item({"nodule": np.zeros((2, 50, 50), bool)}), FRAMES["lidc"])


@settings(max_examples=60, deadline=None)
@given(h=st.integers(1, 700), w=st.integers(1, 700), y=st.integers(0, 1500), x=st.integers(0, 1500),
       ds=st.sampled_from(["riga", "lidc"]))
def test_policy_properties(h, w, y, x, ds):
    Hs, Ws = y + h + 20, x + w + 30
    s = np.stack([blob((Hs, Ws), y, x, y + h, x + w)])
    f = check_frame(item({"nodule": s}, None, ds), FRAMES[ds])
    if f.scale == 1:
        assert f.masks["nodule"].sum() == h * w


@pytest.mark.parametrize("shape", [(80, 80), (80, 80, 1), (80, 80, 3), (80, 80, 4)])
def test_gray_rgb_and_rgba_images_become_rgb(shape):
    """B1-6: gray (2-D or one channel) is repeated to RGB, RGB is kept, RGBA drops its alpha channel."""
    rng = np.random.default_rng(0)
    img = rng.integers(0, 256, shape, dtype=np.uint8)
    it = item({"disc": blob((80, 80), 30, 30, 50, 50)[None], "cup": blob((80, 80), 35, 35, 45, 45)[None]}, img)
    f = frame(it)
    x0, y0, x1, y1 = f.box
    rgb = np.atleast_3d(img)[..., :3]
    expect = np.zeros((f.H, f.W, rgb.shape[2]), np.uint8)
    expect[max(0, -y0):min(f.H, 80 - y0), max(0, -x0):min(f.W, 80 - x0)] = rgb[max(0, y0):y1, max(0, x0):x1]
    assert f.image.shape == (f.H, f.W, 3) and f.image.dtype == np.uint8
    assert np.array_equal(f.image, np.broadcast_to(expect, f.image.shape))


@pytest.mark.parametrize("shape", [(80, 80, 2), (80, 80, 5), (80, 80, 0), (3, 80, 80, 3), (80,)])
def test_other_image_layouts_are_refused(shape):
    it = item({"disc": blob((80, 80), 30, 30, 50, 50)[None], "cup": blob((80, 80), 35, 35, 45, 45)[None]},
              np.zeros(shape, np.uint8))
    with pytest.raises(ValueError, match="gray, RGB or RGBA"):
        frame(it)


def test_frame_masks_are_canonical_bools():
    """I1: the frame is the targets' and caches' mask boundary; True bytes other than 1 become 1."""
    m = np.stack([blob((80, 80), 30, 30, 50, 50), blob((80, 80), 32, 30, 52, 50)])
    nb = (m.view(np.uint8) * np.uint8(255)).view(bool)
    f, g = frame(item({"nodule": nb}, dataset="lidc")), frame(item({"nodule": m}, dataset="lidc"))
    assert f.masks["nodule"].view(np.uint8).max() == 1 and np.array_equal(f.masks["nodule"], g.masks["nodule"])
