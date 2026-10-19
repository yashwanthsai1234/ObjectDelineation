"""MedSAM (SAM ViT-B fine-tuned on medical images), box-prompted, via transformers."""
import numpy as np
import torch
from transformers import SamModel, SamProcessor

MODEL = "flaviagiammarino/medsam-vit-base"


class MedSAM:
    def __init__(self, device=None):
        self.device = device or ("mps" if torch.backends.mps.is_available() else "cpu")
        self.model = SamModel.from_pretrained(MODEL).to(self.device).eval()
        self.processor = SamProcessor.from_pretrained(MODEL)

    @torch.no_grad()
    def segment(self, image_rgb, boxes, threshold=0.5):
        """One HxW bool mask per box [x0, y0, x1, y1] on the same image."""
        inputs = self.processor(image_rgb, input_boxes=[[list(map(float, b)) for b in boxes]], return_tensors="pt")
        inputs = {k: (v.to(self.device, torch.float32) if v.is_floating_point() else v.to(self.device))
                  for k, v in inputs.items()}
        out = self.model(**inputs, multimask_output=False)
        probs = self.processor.image_processor.post_process_masks(
            out.pred_masks.sigmoid().cpu(), inputs["original_sizes"].cpu(),
            inputs["reshaped_input_sizes"].cpu(), binarize=False)[0]
        return [p[0].numpy() > threshold for p in probs]


def box_of(mask, margin):
    """[x0, y0, x1, y1] of a mask, grown by `margin` pixels and clipped to the image."""
    ys, xs = np.nonzero(mask)
    h, w = mask.shape
    return [max(0, xs.min() - margin), max(0, ys.min() - margin),
            min(w - 1, xs.max() + margin), min(h - 1, ys.max() + margin)]
