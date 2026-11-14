"""Read-only diagnostic: no downloads, no writes. Inspects already-cached files."""
import os
import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from PIL import Image

target_dir_fits = "data/test"
target_dir_png = "img/test"

# pick 5 frames spread across what's already on disk: first few IRSA ones,
# plus 013 which we know came from the S3 fallback (no local .fits)
ids = ["000", "001", "002", "003", "013"]

for fid in ids:
    fits_path = f"{target_dir_fits}/{fid}.fits"
    png_path = f"{target_dir_png}/{fid}.png"

    if not os.path.exists(png_path):
        print(f"{fid}: no PNG on disk, skipping")
        continue

    with Image.open(png_path) as im:
        shape = im.size  # (width, height)

    if os.path.exists(fits_path):
        source = "IRSA cutout (local FITS present)"
        with fits.open(fits_path) as hdul:
            header = hdul["IMAGE"].header
        wcs = WCS(header)
        cd = wcs.pixel_scale_matrix  # deg/pixel, CD-equivalent (accounts for PC+CDELT)
        pa_north = np.degrees(np.arctan2(cd[0, 1], cd[1, 1]))
        det = np.linalg.det(cd)
        parity = "flipped (mirror)" if det > 0 else "standard (not mirrored)"
        print(f"{fid}: {source}, png size {shape}")
        print(f"      CD matrix: {cd.tolist()}")
        print(f"      position angle of north: {pa_north:.2f} deg, parity: {parity}")
    else:
        print(f"{fid}: S3 fallback (no local FITS/header was saved for this path -- "
              f"can't report WCS orientation without re-fetching), png size {shape}")
