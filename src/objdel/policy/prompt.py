"""The prompt table: one user turn per dataset (PROMPTS), shared by the RL and SFT rows, G-vision and the
prepared-set manifest (prepare_rl passes it to objdel.data.prepare, which stores it)."""

PROMPTS = {
    "riga": 'Outline the optic disc and the optic cup in this fundus crop. Answer with JSON only: '
            '{"disc":[{"outer":[[x,y],...],"holes":[]}],"cup":[{"outer":[[x,y],...],"holes":[]}]}, '
            'integers 0-1000 relative to the crop, outer rings clockwise.',
    "lidc": 'Outline the lung nodule in this CT crop. Answer with JSON only: '
            '{"nodule":[{"outer":[[x,y],...],"holes":[[[x,y],...]]}]}, '
            'integers 0-1000 relative to the crop, outer rings clockwise, holes counter-clockwise.',
}


def _user(dataset: str) -> dict:
    return {"role": "user", "content": [{"type": "image"}, {"type": "text", "text": PROMPTS[dataset]}]}
