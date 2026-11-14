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

BANDS = ["D1", "D2", "D3", "D4", "D5", "D6"]
MIN_BAND_FRAMES = 30       # D2 is the primary band, used by the current page
MIN_OTHER_BAND_FRAMES = 10  # the other five are extra, lower bar

# micron edges per IRSA's SPHEREx Data Explorer overview (rounded from the
# Explanatory Supplement's 0.744-5.002um full-mission coverage)
BAND_RANGES_UM = {
    "SPHEREx-D1": [0.75, 1.09],
    "SPHEREx-D2": [1.10, 1.62],
    "SPHEREx-D3": [1.63, 2.41],
    "SPHEREx-D4": [2.42, 3.82],
    "SPHEREx-D5": [3.83, 4.41],
    "SPHEREx-D6": [4.42, 5.00],
}


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


def migrate_flat_cache(target_id, size_tag, band_code="D2"):
    """Older runs (before per-band folders) cached files directly under
    data|img/<id>/<size_tag>/. Move them under .../<band_code>/ so they're
    reused instead of re-downloaded."""
    for root in ("data", "img"):
        flat_dir = f"{root}/{target_id}/{size_tag}"
        if not os.path.isdir(flat_dir):
            continue
        band_dir = f"{flat_dir}/{band_code}"
        loose_files = [f for f in os.listdir(flat_dir) if os.path.isfile(f"{flat_dir}/{f}")]
        if not loose_files:
            continue
        os.makedirs(band_dir, exist_ok=True)
        for fname in loose_files:
            shutil.move(f"{flat_dir}/{fname}", f"{band_dir}/{fname}")
        print(f"migrated {len(loose_files)} files: {flat_dir}/*  ->  {band_dir}/")


def fetch_reprojected(url, raw_uri, fits_path, coord, request_size, target_wcs, s3_fallbacks_used, edge_drops):
    """Download (or reuse the cached FITS) and reproject onto the target's shared
    grid. Returns the reprojected array, or None if the frame should be skipped."""
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
            return None
        print("S3 fallback")
        s3_fallbacks_used[0] += 1
        try:
            img, wcs = s3_lazy_cutout(raw_uri, coord, request_size)
        except Exception as e:
            print(f"S3 lazy read failed ({e}), skipping frame (not downloading full file)")
            return None

    # put every frame on the same sky grid, so brightness/position changes
    # between frames are real, not just each exposure's own pointing/roll
    img, _ = reproject_interp((img, wcs), target_wcs, shape_out=OUT_SHAPE)

    nan_frac = np.isnan(img).mean()
    if nan_frac > 0.02:
        edge_drops[0] += 1
        print(f"frame hits detector edge (nan_frac={nan_frac:.1%}), skipping")
        return None
    return img


def process_band(target_id, size_tag, band_code, coord, request_size, target_wcs, frames):
    """Fetch+reproject every frame of one band, then apply one shared
    background-subtract + percentile stretch across the whole stack, so
    brightness is consistent frame to frame instead of each frame being
    stretched to its own contrast."""
    band_full = f"SPHEREx-{band_code}"
    data_dir = f"data/{target_id}/{size_tag}/{band_code}"
    img_dir = f"img/{target_id}/{size_tag}/{band_code}"
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(img_dir, exist_ok=True)

    s3_fallbacks_used = [0]  # mutable int, shared across this band's fetch calls
    edge_drops = [0]
    collected = []  # (png_path, time_iso, img)
    for row in frames:
        # stable id from the filename IRSA/S3 already gave this exposure, so a
        # cached file always matches the frame it was downloaded for, even if
        # a later run samples a different subset of frames
        stable_id = os.path.splitext(os.path.basename(row["raw_uri"]))[0]
        fits_path = f"{data_dir}/{stable_id}.fits"
        png_path = f"{img_dir}/{stable_id}.png"
        img = fetch_reprojected(row["uri"], row["raw_uri"], fits_path, coord, request_size, target_wcs,
                                 s3_fallbacks_used, edge_drops)
        if img is not None:
            time_iso = Time(row["time_bounds_lower"], format="mjd").isot  # MJD days -> ISO for the frontend
            collected.append((png_path, time_iso, img))
        time.sleep(1)  # be polite to IRSA

    print(f"{target_id}/{band_code}: {len(collected)}/{len(frames)} frames succeeded, "
          f"{edge_drops[0]} dropped for hitting detector edge")
    if not collected:
        return []

    bg_subbed = [img - np.nanmedian(img) for _, _, img in collected]
    stack = np.array(bg_subbed)
    lo, hi = np.nanpercentile(stack, [1, 99.5])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        print(f"{target_id}/{band_code}: degenerate stack stats, skipping band")
        return []

    out_frames = []
    for (png_path, time_iso, _), img_sub in zip(collected, bg_subbed):
        scaled = np.clip((img_sub - lo) / (hi - lo), 0, 1)
        scaled = np.arcsinh(10 * scaled) / np.arcsinh(10)   # brings out faint stuff
        scaled = np.nan_to_num(scaled, nan=0.0)
        scaled = np.flipud(scaled)                          # FITS has y flipped
        Image.fromarray((scaled * 255).astype(np.uint8)).save(png_path)
        out_frames.append({"url": png_path, "time": time_iso, "band": band_full})
    return out_frames


def process_target(target_id, ra, dec, label, size=0.2, max_frames=40, bands=None):
    bands = bands if bands is not None else BANDS
    coord = SkyCoord(ra=ra, dec=dec, unit="deg")
    request_size = size * 1.5  # request more than the display size so reprojection has margin
    target_wcs = build_target_wcs(ra, dec, size, npix=OUT_SHAPE[0])
    size_tag = f"{size:g}".replace(".", "p")  # e.g. 0.2 -> "0p2", keeps caches from different sizes apart

    migrate_flat_cache(target_id, size_tag)  # old single-band cache -> .../D2/

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

    bands_out = {}
    for band_code in bands:
        band_full = f"SPHEREx-{band_code}"
        frames = results[results["energy_bandpassname"] == band_full]
        min_needed = MIN_BAND_FRAMES if band_code == "D2" else MIN_OTHER_BAND_FRAMES
        if len(frames) < min_needed:
            print(f"SKIP {target_id} {band_code}: only {len(frames)} frames (< {min_needed})")
            continue

        if len(frames) > max_frames:
            # frames is already time-sorted; spread picks across the full span
            # instead of clustering near the start, so the blink viewer still
            # covers the whole observing window
            idx = np.unique(np.linspace(0, len(frames) - 1, max_frames).round().astype(int))
            frames = frames[idx]
            print(f"{target_id}/{band_code}: sampled down to {len(frames)} frames across full time range")

        out_frames = process_band(target_id, size_tag, band_code, coord, request_size, target_wcs, frames)
        if out_frames:
            bands_out[band_full] = out_frames

    if not bands_out:
        return None

    return {
        "id": target_id,
        "label": label,
        "ra": ra,
        "dec": dec,
        "frames": bands_out.get("SPHEREx-D2", []),  # keeps the current page working off D2
        "bands": bands_out,
        "band_ranges": BAND_RANGES_UM,
    }


def count_frames(ra, dec):
    """TAP query only, no downloads -- frames per band, plus D2's time span."""
    coord = SkyCoord(ra=ra, dec=dec, unit="deg")
    query = f"""
    SELECT p.energy_bandpassname, p.time_bounds_lower
    FROM spherex.artifact a
    JOIN spherex.plane p ON a.planeid = p.planeid
    WHERE 1 = CONTAINS(POINT('ICRS', {coord.ra.deg}, {coord.dec.deg}), p.poly)
    """
    results = Irsa.query_tap(query).to_table()
    counts = Counter(results["energy_bandpassname"])
    in_d2 = results[results["energy_bandpassname"] == "SPHEREx-D2"]
    span = None
    if len(in_d2):
        t = in_d2["time_bounds_lower"]
        span = (Time(float(np.min(t)), format="mjd").isot[:10],
                Time(float(np.max(t)), format="mjd").isot[:10])
    return counts, span


def count_all_targets():
    for target_id, ra, dec, label in TARGETS:
        counts, span = count_frames(ra, dec)
        print(f"{target_id} ({label}):")
        for band_code in BANDS:
            n = counts.get(f"SPHEREx-{band_code}", 0)
            print(f"  {band_code}: {n}")
        if span:
            print(f"  D2 observed {span[0]} to {span[1]}")


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
    ("m101", 210.80227, 54.34895, "M101 (Pinwheel Galaxy)"),
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
#   python prep/test_one.py --count                   frames per band + D2 time span, no downloads
#   python prep/test_one.py --clean                    delete data/, img/, targets.json (asks first)
#   python prep/test_one.py --only pleiades             run just one target (5 frames/band)
#   python prep/test_one.py --only pleiades --full       all bands, up to 40 frames/band
#   python prep/test_one.py --only pleiades --bands D1,D2 --full   only those bands
#   python prep/test_one.py --full                      all targets, all bands
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

    bands = BANDS
    if "--bands" in sys.argv:
        bands = sys.argv[sys.argv.index("--bands") + 1].split(",")

    kwargs = {"bands": bands}
    if "--full" not in sys.argv:
        kwargs["max_frames"] = TEST_MAX_FRAMES

    for target_id, ra, dec, label in TARGETS:
        if only and target_id != only:
            continue
        entry = process_target(target_id, ra, dec, label, **kwargs)
        if entry:
            n = save_target_entry(entry)
            print(f"saved {target_id}, targets.json now has {n} targets")
