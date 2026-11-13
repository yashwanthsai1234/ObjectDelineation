"""Outline text and rings for the verifier and reward tests."""
from objdel.geom.text import serialize


def square_ring(x0, y0, x1, y1, hole=False, scale=1) -> list:
    """The clockwise box (x0, y0)-(x1, y1) as four points, coordinates times `scale`; `hole` reverses it."""
    ring = [[x0 * scale, y0 * scale], [x1 * scale, y0 * scale], [x1 * scale, y1 * scale], [x0 * scale, y1 * scale]]
    return ring[::-1] if hole else ring


def outline_text(ring, holes=(), name="nodule") -> str:
    """The compact answer text of one object part: outer `ring` and the hole rings `holes`."""
    return serialize({name: [(ring, list(holes))]}, (name,))
