"""Read-only diagnostic: no downloads (offline=True), reuses cached D2 FITS for
pleiades. Measures the brightest star's centroid per frame after reprojection."""
import os
import numpy as np
from astropy.coordinates import SkyCoord
from test_one import build_target_wcs, fetch_reprojected, OUT_SHAPE, Irsa

ra, dec, size = 56.75, 24.12, 0.2
request_size = size * 1.5
coord = SkyCoord(ra=ra, dec=dec, unit="deg")
target_wcs = build_target_wcs(ra, dec, size, npix=OUT_SHAPE[0])

query = f"""
SELECT 'https://irsa.ipac.caltech.edu/' || a.uri || '?center={coord.ra.deg},{coord.dec.deg}d&size={request_size}' AS uri,
       a.uri AS raw_uri, p.time_bounds_lower, p.energy_bandpassname
FROM spherex.artifact a JOIN spherex.plane p ON a.planeid = p.planeid
WHERE 1 = CONTAINS(POINT('ICRS', {coord.ra.deg}, {coord.dec.deg}), p.poly)
ORDER BY p.time_bounds_lower
"""
results = Irsa.query_tap(query).to_table()
frames = results[results["energy_bandpassname"] == "SPHEREx-D2"]
idx = np.unique(np.linspace(0, len(frames) - 1, 40).round().astype(int))
frames = frames[idx]

s3_fallbacks_used = [0]
edge_drops = [0]
imgs = []
for row in frames:
    stable_id = os.path.splitext(os.path.basename(row["raw_uri"]))[0]
    fits_path = f"data/pleiades/0p2/D2/{stable_id}.fits"
    img = fetch_reprojected(row["uri"], row["raw_uri"], fits_path, coord, request_size, target_wcs,
                             s3_fallbacks_used, edge_drops, offline=True)
    if img is not None:
        imgs.append(img)

print(f"{len(imgs)} frames used")
stack = np.array(imgs)
median_img = np.nanmedian(stack, axis=0)
py, px = np.unravel_index(np.nanargmax(median_img), median_img.shape)
print(f"brightest pixel in median image: row={py}, col={px}")

box = 7
half = box // 2
centroids = []
for img in stack:
    sub = img[py - half:py + half + 1, px - half:px + half + 1]
    yy, xx = np.mgrid[py - half:py + half + 1, px - half:px + half + 1]
    w = np.nan_to_num(sub, nan=0.0)
    w = np.clip(w, 0, None)
    total = w.sum()
    if total <= 0:
        continue
    cy = (yy * w).sum() / total
    cx = (xx * w).sum() / total
    centroids.append((cx, cy))

centroids = np.array(centroids)
print(f"{len(centroids)} centroids measured")
print(f"x: std={centroids[:,0].std():.3f}px  range={centroids[:,0].max()-centroids[:,0].min():.3f}px")
print(f"y: std={centroids[:,1].std():.3f}px  range={centroids[:,1].max()-centroids[:,1].min():.3f}px")
