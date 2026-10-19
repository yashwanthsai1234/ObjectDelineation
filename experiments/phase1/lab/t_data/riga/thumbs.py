"""Cache grey 576x384 thumbnails of the reference images: python -I thumbs.py <refs> <out.npz>"""
import sys, numpy as np
from pathlib import Path
from PIL import Image
files = sorted(Path(sys.argv[1]).glob("*.png"))
T = np.stack([np.asarray(Image.open(f).convert("L").resize((576, 384), Image.BILINEAR)) for f in files])
np.savez_compressed(sys.argv[2], thumbs=T, ids=np.array([f.name[:-4] for f in files]))
