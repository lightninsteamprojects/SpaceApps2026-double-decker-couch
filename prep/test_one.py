from astroquery.ipac.irsa import Irsa
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
import numpy as np
import requests
import time
import os
import json
import shutil
import s3fs
from collections import Counter
from PIL import Image
from reproject import reproject_interp

S3_FS = s3fs.S3FileSystem(anon=True)  # reused across frames, no per-frame login cost

# S3 has no cutout service, so a failed frame falls back to a lazy, ranged read
# of just the region we need (block_size=200KB keeps it ~a few MB, not the full
# ~65MB exposure). If this budget is exhausted, remaining frames are skipped
# instead of silently downloading full files.
MAX_S3_FALLBACKS_PER_TARGET = 15

OUT_SHAPE = (128, 128)  # shared pixel grid every frame of a target gets reprojected onto

# one band for every target, so the whole site is consistent. Falls back to the
# target's most common band only if it has fewer than MIN_BAND_FRAMES in BAND.
BAND = "SPHEREx-D2"
MIN_BAND_FRAMES = 30


def download(url, out_path, tries=3):
    for i in range(tries):
        try:
            r = requests.get(url, timeout=120)
            if r.status_code == 200:
                with open(out_path, "wb") as f:
                    f.write(r.content)
                return True
            print(f"status {r.status_code}, retry {i+1}/{tries}")
        except requests.RequestException as e:
            print(f"error {e}, retry {i+1}/{tries}")
        time.sleep(5 * (i + 1))      # wait 5s, 10s, 15s
    return False


def raw_uri_to_s3_path(raw_uri):
    key = raw_uri.split("ibe/data/spherex/", 1)[1]
    return f"nasa-irsa-spherex/{key}"


def s3_lazy_cutout(raw_uri, coord, request_size):
    """Read only the cutout region from the full S3 frame via HTTP range requests."""
    s3_path = raw_uri_to_s3_path(raw_uri)
    with S3_FS.open(s3_path, mode="rb", cache_type="readahead", block_size=200_000) as f:
        with fits.open(f, lazy_load_hdus=True) as hdul:
            hdu = hdul["IMAGE"]
            wcs = WCS(hdu.header)
            cut = Cutout2D(hdu.section, coord, request_size * u.deg, wcs=wcs)
            return cut.data.astype(float), cut.wcs


def build_target_wcs(ra, dec, size, npix=128):
    """One shared TAN grid per target, north up, covering `size` degrees, so every
    frame can be reprojected onto the same pixels regardless of its own roll angle."""
    pixscale = size / npix
    w = WCS(naxis=2)
    w.wcs.crpix = [npix / 2 + 0.5, npix / 2 + 0.5]
    w.wcs.cdelt = [-pixscale, pixscale]
    w.wcs.crval = [ra, dec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return w


def make_png(url, raw_uri, fits_path, png_path, coord, request_size, target_wcs, s3_fallbacks_used, edge_drops):
    if os.path.exists(fits_path):
        with fits.open(fits_path) as hdul:
            img = hdul["IMAGE"].data.astype(float)
            wcs = WCS(hdul["IMAGE"].header)
    elif download(url, fits_path):
        with fits.open(fits_path) as hdul:
            img = hdul["IMAGE"].data.astype(float)
            wcs = WCS(hdul["IMAGE"].header)
    else:
        if s3_fallbacks_used[0] >= MAX_S3_FALLBACKS_PER_TARGET:
            print("S3 fallback budget exhausted, skipping frame")
            return False
        print("S3 fallback")
        s3_fallbacks_used[0] += 1
        try:
            img, wcs = s3_lazy_cutout(raw_uri, coord, request_size)
        except Exception as e:
            print(f"S3 lazy read failed ({e}), skipping frame (not downloading full file)")
            return False

    # put every frame on the same sky grid, so brightness/position changes
    # between frames are real, not just each exposure's own pointing/roll
    img, _ = reproject_interp((img, wcs), target_wcs, shape_out=OUT_SHAPE)

    nan_frac = np.isnan(img).mean()
    if nan_frac > 0.02:
        edge_drops[0] += 1
        print(f"frame hits detector edge (nan_frac={nan_frac:.1%}), skipping")
        return False

    img = np.nan_to_num(img, nan=np.nanmedian(img))
    lo, hi = np.percentile(img, [1, 99])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return False  # degenerate frame (e.g. all-NaN), skip it
    img = np.clip((img - lo) / (hi - lo), 0, 1)
    img = np.arcsinh(10 * img) / np.arcsinh(10)     # brings out faint stuff
    img = np.flipud(img)                            # FITS has y flipped
    Image.fromarray((img * 255).astype(np.uint8)).save(png_path)
    return True


def process_target(target_id, ra, dec, label, size=0.2, max_frames=40):
    coord = SkyCoord(ra=ra, dec=dec, unit="deg")
    request_size = size * 1.5  # request more than the display size so reprojection has margin
    target_wcs = build_target_wcs(ra, dec, size, npix=OUT_SHAPE[0])
    size_tag = f"{size:g}".replace(".", "p")  # e.g. 0.2 -> "0p2", keeps caches from different sizes apart

    query = f"""
    SELECT
      'https://irsa.ipac.caltech.edu/' || a.uri
        || '?center={coord.ra.deg},{coord.dec.deg}d&size={request_size}' AS uri,
      a.uri AS raw_uri,
      p.time_bounds_lower,
      p.energy_bandpassname
    FROM spherex.artifact a
    JOIN spherex.plane p ON a.planeid = p.planeid
    WHERE 1 = CONTAINS(POINT('ICRS', {coord.ra.deg}, {coord.dec.deg}), p.poly)
    ORDER BY p.time_bounds_lower
    """
    results = Irsa.query_tap(query).to_table()
    print(f"{target_id}: TAP found {len(results)} frames")

    # same bandpass only, so frame-to-frame brightness changes aren't just band effects
    band = BAND
    frames = results[results["energy_bandpassname"] == band]
    if len(frames) < MIN_BAND_FRAMES:
        band = Counter(results["energy_bandpassname"]).most_common(1)[0][0]
        frames = results[results["energy_bandpassname"] == band]
        print(f"{target_id}: only few {BAND} frames, falling back to {band}")
    print(f"{target_id}: using band {band}, {len(frames)} frames")

    if len(frames) > max_frames:
        # frames is already time-sorted; spread picks across the full span
        # instead of clustering near the start, so the blink viewer still
        # covers the whole observing window
        idx = np.unique(np.linspace(0, len(frames) - 1, max_frames).round().astype(int))
        frames = frames[idx]
        print(f"{target_id}: sampled down to {len(frames)} frames across full time range")

    os.makedirs(f"data/{target_id}/{size_tag}", exist_ok=True)
    os.makedirs(f"img/{target_id}/{size_tag}", exist_ok=True)

    s3_fallbacks_used = [0]  # mutable int, shared across the loop's make_png calls
    edge_drops = [0]
    out_frames = []
    for row in frames:
        # stable id from the filename IRSA/S3 already gave this exposure, so a
        # cached file always matches the frame it was downloaded for, even if
        # a later run samples a different subset of frames
        stable_id = os.path.splitext(os.path.basename(row["raw_uri"]))[0]
        fits_path = f"data/{target_id}/{size_tag}/{stable_id}.fits"
        png_path = f"img/{target_id}/{size_tag}/{stable_id}.png"
        ok = make_png(row["uri"], row["raw_uri"], fits_path, png_path, coord, request_size, target_wcs, s3_fallbacks_used, edge_drops)
        if ok:
            # time_bounds_lower is MJD (days); convert to ISO for the frontend
            time_iso = Time(row["time_bounds_lower"], format="mjd").isot
            out_frames.append({"url": png_path, "time": time_iso, "band": band})
        time.sleep(1)  # be polite to IRSA

    print(f"{target_id}: {len(out_frames)}/{len(frames)} frames succeeded, {edge_drops[0]} dropped for hitting detector edge")
    if not out_frames:
        return None
    return {"id": target_id, "label": label, "ra": ra, "dec": dec, "band": band, "frames": out_frames}


def count_frames(ra, dec):
    """TAP query only, no downloads -- frames per band, plus the time span of BAND frames."""
    coord = SkyCoord(ra=ra, dec=dec, unit="deg")
    query = f"""
    SELECT p.energy_bandpassname, p.time_bounds_lower
    FROM spherex.artifact a
    JOIN spherex.plane p ON a.planeid = p.planeid
    WHERE 1 = CONTAINS(POINT('ICRS', {coord.ra.deg}, {coord.dec.deg}), p.poly)
    """
    results = Irsa.query_tap(query).to_table()
    counts = Counter(results["energy_bandpassname"])
    in_band = results[results["energy_bandpassname"] == BAND]
    span = None
    if len(in_band):
        t = in_band["time_bounds_lower"]
        span = (Time(float(np.min(t)), format="mjd").isot[:10],
                Time(float(np.max(t)), format="mjd").isot[:10])
    return counts, span


def count_all_targets():
    for target_id, ra, dec, label in TARGETS:
        counts, span = count_frames(ra, dec)
        print(f"{target_id} ({label}):")
        for band, n in counts.most_common():
            print(f"  {band}: {n}")
        if span:
            print(f"  {BAND} observed {span[0]} to {span[1]}")


TARGETS_JSON = "targets.json"


def load_targets():
    if os.path.exists(TARGETS_JSON):
        with open(TARGETS_JSON) as f:
            return json.load(f)
    return []


def save_target_entry(entry):
    targets = [t for t in load_targets() if t["id"] != entry["id"]]
    targets.append(entry)
    with open(TARGETS_JSON, "w") as f:
        json.dump(targets, f, indent=2)
    return len(targets)


def clean():
    answer = input("Delete data/, img/ and targets.json? (y/N) ")
    if answer.strip().lower() != "y":
        print("Aborted.")
        return
    for path in ("data", "img", TARGETS_JSON):
        if os.path.isdir(path):
            shutil.rmtree(path)
            print(f"deleted {path}/")
        elif os.path.isfile(path):
            os.remove(path)
            print(f"deleted {path}")


TARGETS = [
    # (id, ra, dec, label)
    ("test", 210.80227, 54.34895, "Test field"),
    # ecliptic-plane fields, best odds of catching a passing asteroid
    ("pleiades", 56.75, 24.12, "Pleiades"),
    ("beehive", 130.10, 19.67, "Beehive Cluster (M44)"),
    ("ecliptic0", 0.0, 0.0, "Ecliptic point RA 0"),
    ("ecliptic180", 180.0, 0.0, "Ecliptic point RA 180"),
    ("galactic_center", 266.40, -28.94, "Galactic Center region"),
    # ecliptic poles, likely many more visits (check with --count)
    ("north_pole", 270.0, 66.56, "North ecliptic pole"),
    ("south_pole", 90.0, -66.56, "South ecliptic pole"),
]

# frames per target in a quick test run; pass --full to use the normal default (40)
TEST_MAX_FRAMES = 5

# Usage:
#   python prep/test_one.py --count              frames per band + time span, no downloads
#   python prep/test_one.py --clean              delete data/, img/, targets.json (asks first)
#   python prep/test_one.py --only pleiades      run just one target (5 frames)
#   python prep/test_one.py --only pleiades --full
#   python prep/test_one.py --full               all targets, up to 40 frames each
if __name__ == "__main__":
    import sys

    if "--clean" in sys.argv:
        clean()
        sys.exit()

    if "--count" in sys.argv:
        count_all_targets()
        sys.exit()

    only = None
    if "--only" in sys.argv:
        only = sys.argv[sys.argv.index("--only") + 1]

    kwargs = {} if "--full" in sys.argv else {"max_frames": TEST_MAX_FRAMES}
    for target_id, ra, dec, label in TARGETS:
        if only and target_id != only:
            continue
        entry = process_target(target_id, ra, dec, label, **kwargs)
        if entry:
            n = save_target_entry(entry)
            print(f"saved {target_id}, targets.json now has {n} targets")