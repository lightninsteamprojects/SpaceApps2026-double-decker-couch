"""Read-only diagnostic: checks how well-aligned the cached Pleiades frames are
after reprojection. Writes only img/pleiades_check.png."""
import json
import numpy as np
from PIL import Image


def phase_cross_correlation(ref, img):
    f1 = np.fft.fft2(ref)
    f2 = np.fft.fft2(img)
    cross_power = (f1 * np.conj(f2))
    cross_power /= np.abs(cross_power) + 1e-10
    corr = np.abs(np.fft.ifft2(cross_power))
    corr = np.fft.fftshift(corr)
    peak_y, peak_x = np.unravel_index(np.argmax(corr), corr.shape)
    dy = peak_y - ref.shape[0] // 2
    dx = peak_x - ref.shape[1] // 2
    return dx, dy


with open("targets.json") as f:
    targets = json.load(f)
entry = next(t for t in targets if t["id"] == "pleiades")
png_paths = [fr["url"] for fr in entry["frames"]]
print(f"{len(png_paths)} frames, in time order:")

arrays = [np.array(Image.open(p).convert("L"), dtype=float) for p in png_paths]
ref = arrays[0]
for path, arr in zip(png_paths, arrays):
    dx, dy = phase_cross_correlation(ref, arr)
    print(f"  {path}: dx={dx}, dy={dy}")

# contact sheet: each frame upscaled 4x nearest-neighbour, side by side
tiles = [Image.open(p).convert("L").resize(
    (Image.open(p).width * 4, Image.open(p).height * 4), Image.NEAREST) for p in png_paths]
w, h = tiles[0].size
sheet = Image.new("L", (w * len(tiles), h))
for i, tile in enumerate(tiles):
    sheet.paste(tile, (i * w, 0))
sheet.save("img/pleiades_check.png")
print("saved img/pleiades_check.png")
