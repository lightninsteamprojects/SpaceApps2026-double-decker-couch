from astroquery.ipac.irsa import Irsa
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import astropy.units as u

coord = SkyCoord(ra=210.80227, dec=54.34895, unit="deg")
size = 0.1

query = f"""
SELECT a.uri AS raw_uri
FROM spherex.artifact a
JOIN spherex.plane p ON a.planeid = p.planeid
WHERE 1 = CONTAINS(POINT('ICRS', {coord.ra.deg}, {coord.dec.deg}), p.poly)
ORDER BY p.time_bounds_lower
"""
results = Irsa.query_tap(query).to_table()
raw_uri = results["raw_uri"][0]
key = raw_uri.split("ibe/data/spherex/", 1)[1]
s3_url = f"s3://nasa-irsa-spherex/{key}"
print("raw_uri:", raw_uri)
print("key:", key)
print("s3_url guess:", s3_url)

import time
import s3fs

total_bytes = 0
orig_fetch = s3fs.S3File._fetch_range

def counting_fetch(self, start, end):
    global total_bytes
    total_bytes += (end - start)
    print(f"  GET range {start}-{end} ({end-start} bytes)")
    return orig_fetch(self, start, end)

s3fs.S3File._fetch_range = counting_fetch

t0 = time.time()
fs = s3fs.S3FileSystem(anon=True)
s3_path = s3_url.replace("s3://", "")
with fs.open(s3_path, mode="rb", cache_type="readahead", block_size=200_000) as f:
    with fits.open(f, lazy_load_hdus=True) as hdul:
        hdul.info()
        hdu = hdul["IMAGE"]
        print("shape:", hdu.shape)
        wcs = WCS(hdu.header)
        cut = Cutout2D(hdu.section, coord, size * u.deg, wcs=wcs)
        print("cutout shape:", cut.data.shape)
print("elapsed:", time.time() - t0, "s")
print("TOTAL BYTES FETCHED FROM S3:", total_bytes)
