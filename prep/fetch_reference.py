"""Fetch one reference image per field from NASA/IPAC IRSA's 2MASS image service.

(NASA SkyView was tried first but its server could not reach its 2MASS/WISE archives,
so images now come straight from IRSA. DSS is deliberately not used.)

For each field: ask IRSA which 2MASS tiles cover the position, download a K-band tile
that fully contains the 0.2 deg field, crop it, arcsinh-stretch it, and save
img/ref/<id>.jpg (grayscale, north up, 600 px wide). Also writes ref.json.
Run from the project root:  python prep/fetch_reference.py
"""
import io, json, os, re
import numpy as np
import requests
from PIL import Image
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from astropy.wcs import WCS

SIZE_DEG = 0.2   # same sky area as the SPHEREx field images (the "0p2" folders)
NPIX = 600
SIA = "https://irsa.ipac.caltech.edu/cgi-bin/2MASS/IM/nph-im_sia"
CREDIT = "2MASS images via NASA/IPAC IRSA"
# 2MASS tile file names start with the band letter (ki = K, hi = H, ji = J); K is closest to SPHEREx
BANDS = [("k", "2MASS K-band"), ("h", "2MASS H-band"), ("j", "2MASS J-band")]


def tile_urls(ra, dec, letter):
    r = requests.get(SIA, params=dict(POS=f"{ra},{dec}", SIZE=SIZE_DEG, FORMAT="image/fits"), timeout=60)
    return re.findall(rf"(https://[^<\]\s]*name={letter}i[^<\]\s]*fits)", r.text)


def mosaic(urls, ra, dec):
    """Same grid as the SPHEREx frames (build_target_wcs in test_one.py: TAN, north up, centred on the
    field, covering SIZE_DEG) but 600 px across instead of 128, so 2MASS detail is kept.
    The field usually straddles tile edges, so all overlapping tiles are mosaicked (with their
    background levels matched) and reprojected in one go."""
    w = WCS(naxis=2)
    w.wcs.crpix = [NPIX / 2 + 0.5, NPIX / 2 + 0.5]
    w.wcs.cdelt = [-SIZE_DEG / NPIX, SIZE_DEG / NPIX]
    w.wcs.crval = [ra, dec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    tiles = []
    for url in urls:
        hdu = fits.open(io.BytesIO(requests.get(url, timeout=120).content))[0]
        tiles.append((hdu.data, WCS(hdu.header)))
    out, _ = reproject_and_coadd(tiles, w, shape_out=(NPIX, NPIX), reproject_function=reproject_interp,
                                 match_background=True)
    if np.isnan(out).mean() > 0.02:
        raise ValueError("tiles do not cover the whole field")
    return out


def to_jpg(data, path):
    data = np.nan_to_num(data.astype(float))
    lo, hi = np.percentile(data, [1, 99.7])          # clip the dark sky and the brightest stars
    x = np.clip((data - lo) / (hi - lo), 0, 1)
    x = np.arcsinh(10 * x) / np.arcsinh(10)          # arcsinh stretch: faint stuff visible, stars not blown out
    x = np.flipud(x)                                 # array row 0 is south; image row 0 must be north (same flip as test_one.py)
    Image.fromarray((x * 255).astype(np.uint8)).save(path, quality=90)


os.makedirs("img/ref", exist_ok=True)
ref, failed = {}, []
for t in json.load(open("targets.json")):
    coord = SkyCoord(t["ra"], t["dec"], unit="deg", frame="icrs")
    done = False
    for letter, name in BANDS:
        try:
            data = mosaic(tile_urls(t["ra"], t["dec"], letter), t["ra"], t["dec"])
            to_jpg(data, f"img/ref/{t['id']}.jpg")
            ref[t["id"]] = {"file": f"img/ref/{t['id']}.jpg", "survey": name, "credit": CREDIT}
            print(t["id"], "->", name)
            done = True
        except Exception as e:
            print(t["id"], name, "failed:", e)
        if done:
            break
    if not done:
        failed.append(t["id"])

json.dump(ref, open("ref.json", "w"), indent=1)
print("skipped:", failed or "none")
