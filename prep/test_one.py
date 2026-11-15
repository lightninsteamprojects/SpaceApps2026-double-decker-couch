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
from scipy.ndimage import binary_dilation, median_filter, gaussian_filter, shift as ndi_shift, label as ndi_label
from astropy.stats import sigma_clip
from skimage.registration import phase_cross_correlation

S3_FS = s3fs.S3FileSystem(anon=True)  # reused across frames, no per-frame login cost

# S3 has no cutout service, so a failed frame falls back to a lazy, ranged read
# of just the region we need (block_size=200KB keeps it ~a few MB, not the full
# ~65MB exposure). If this budget is exhausted, remaining frames are skipped
# instead of silently downloading full files.
MAX_S3_FALLBACKS_PER_TARGET = 20

# pixels this bright in the static-sky median are masked to mid-grey in the diff
# images (a bright star's own PSF/alignment residual isn't a moving/changing source)
BRIGHT_MASK_PERCENTILE = 99.5  # dimmer tier: covers an ordinary star's halo
BRIGHT_MASK_DILATION_ITERATIONS = 6
BRIGHT_MASK_PERCENTILE_STRONG = 99.9  # brighter tier: wider dilation for a saturating star's larger residual
BRIGHT_MASK_DILATION_ITERATIONS_STRONG = 12

# how each frame's smooth background is estimated before building the median/diff:
# "poly2" fits a degree-2 2D polynomial with sigma-clipping (stars excluded as outliers);
# "median_filter" instead runs a large median filter over the frame.
BACKGROUND_METHOD = "poly2"
MEDIAN_FILTER_SIZE = 31  # pixels, only used when BACKGROUND_METHOD == "median_filter"
POLY_SIGMA_CLIP = 3      # only used when BACKGROUND_METHOD == "poly2"

# sub-pixel alignment of each frame against the median, before differencing
ALIGN_UPSAMPLE_FACTOR = 20
ALIGN_CLIP_PERCENTILE = 99.9  # clip the brightest 0.1% before cross-correlating, so a
                               # saturated core doesn't dominate the shift estimate

# Gaussian-smooths each frame before differencing, to soften the +/- dipoles that
# sub-pixel alignment residuals leave at point sources. 0 disables it.
GAUSSIAN_SIGMA_DIFF = 0.7

# registration guard: a frame's estimated shift is only trusted (and applied) if
# it's small and phase_cross_correlation itself is confident about it. error is
# skimage's normalized RMS registration error: 0 = perfect match, ~1 = no
# correlation at all, so 0.5 is a middling-permissive cutoff that still rejects
# clearly bad matches. Anything that fails either check gets shift=(0,0) instead.
REGISTRATION_MAX_SHIFT_PX = 2.0
REGISTRATION_MAX_ERROR = 0.5
REGISTRATION_MIN_SOURCES = 3     # below this many point sources in the median, skip alignment for the whole band
REGISTRATION_MIN_PASS_FRAC = 0.7  # below this fraction of confident frames, don't write diff images for the band
REGISTRATION_SOURCE_SIGMA = 5     # how many sigma above background counts as a "point source" for the count above

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


def save_cutout_fits(data, wcs, out_path):
    """Writes an S3 lazy cutout to disk as a normal FITS the rest of the code can
    read (same IMAGE extension/header shape as a full IRSA download), so later
    runs and --diff-only find it in the cache instead of re-fetching from S3."""
    hdu = fits.ImageHDU(data=data.astype(np.float32), header=wcs.to_header(), name="IMAGE")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(out_path, overwrite=True)


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


def fetch_reprojected(url, raw_uri, fits_path, coord, request_size, target_wcs, s3_fallbacks_used, edge_drops,
                       offline=False):
    """Download (or reuse the cached FITS) and reproject onto the target's shared
    grid. Returns the reprojected array, or None if the frame should be skipped.
    offline=True never downloads -- only cached FITS files are used (for --diff-only)."""
    if os.path.exists(fits_path):
        with fits.open(fits_path) as hdul:
            img = hdul["IMAGE"].data.astype(float)
            wcs = WCS(hdul["IMAGE"].header)
    elif offline:
        print("not cached, skipping (--diff-only)")
        return None
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
            save_cutout_fits(img, wcs, fits_path)  # cache it so future runs reuse this, not re-fetch
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


def fit_poly2_background(img):
    """Degree-2 2D polynomial fit, excluding stars via sigma-clipping."""
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    valid = np.isfinite(img)
    data, xs, ys = img[valid], xx[valid], yy[valid]
    keep = ~sigma_clip(data, sigma=POLY_SIGMA_CLIP, maxiters=5).mask
    xs, ys, data = xs[keep], ys[keep], data[keep]
    A = np.column_stack([np.ones_like(xs), xs, ys, xs**2, xs * ys, ys**2])
    coeffs, *_ = np.linalg.lstsq(A, data, rcond=None)
    return coeffs[0] + coeffs[1] * xx + coeffs[2] * yy + coeffs[3] * xx**2 + coeffs[4] * xx * yy + coeffs[5] * yy**2


def estimate_background(img):
    if BACKGROUND_METHOD == "median_filter":
        return median_filter(np.nan_to_num(img, nan=np.nanmedian(img)), size=MEDIAN_FILTER_SIZE)
    return fit_poly2_background(img)


def clip_bright(img, pct):
    img = np.nan_to_num(img, nan=0.0)
    return np.clip(img, None, np.nanpercentile(img, pct))


def count_point_sources(median_img, sigma_thresh=REGISTRATION_SOURCE_SIGMA):
    """Rough point-source count: connected components this many sigma above the
    median image's own background, used to decide whether there's enough to
    register against at all."""
    img = np.nan_to_num(median_img, nan=np.nanmedian(median_img))
    bg = np.nanmedian(img)
    sigma = 1.4826 * np.nanmedian(np.abs(img - bg))
    if not np.isfinite(sigma) or sigma <= 0:
        return 0
    _, n = ndi_label(img > bg + sigma_thresh * sigma)
    return n


def process_band(target_id, size_tag, band_code, coord, request_size, target_wcs, frames, offline=False):
    """Fetch+reproject every frame of one band, then apply one shared
    background-subtract + percentile stretch across the whole stack, so
    brightness is consistent frame to frame instead of each frame being
    stretched to its own contrast. Also builds a difference image per frame
    (frame minus the per-pixel time median) to make moving/changing sources
    pop out against the mostly-static star field."""
    band_full = f"SPHEREx-{band_code}"
    data_dir = f"data/{target_id}/{size_tag}/{band_code}"
    img_dir = f"img/{target_id}/{size_tag}/{band_code}"
    diff_dir = f"{img_dir}/diff"
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(img_dir, exist_ok=True)

    s3_fallbacks_used = [0]  # mutable int, shared across this band's fetch calls
    edge_drops = [0]
    collected = []  # (png_path, stable_id, time_iso, img)
    for row in frames:
        # stable id from the filename IRSA/S3 already gave this exposure, so a
        # cached file always matches the frame it was downloaded for, even if
        # a later run samples a different subset of frames
        stable_id = os.path.splitext(os.path.basename(row["raw_uri"]))[0]
        fits_path = f"{data_dir}/{stable_id}.fits"
        png_path = f"{img_dir}/{stable_id}.png"
        img = fetch_reprojected(row["uri"], row["raw_uri"], fits_path, coord, request_size, target_wcs,
                                 s3_fallbacks_used, edge_drops, offline=offline)
        if img is not None:
            time_iso = Time(row["time_bounds_lower"], format="mjd").isot  # MJD days -> ISO for the frontend
            collected.append((png_path, stable_id, time_iso, img))
        if not offline:
            time.sleep(1)  # be polite to IRSA

    print(f"{target_id}/{band_code}: {len(collected)}/{len(frames)} frames succeeded, "
          f"{edge_drops[0]} dropped for hitting detector edge")
    if not collected:
        return [], 0

    # subtract a smooth per-frame background (poly2 fit or large median filter)
    # instead of a single median value, so real structure -- not a sky gradient -- is
    # what's left to difference
    bg_subbed = [img - estimate_background(img) for _, _, _, img in collected]
    stack = np.array(bg_subbed)
    lo, hi = np.nanpercentile(stack, [1, 99.5])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        print(f"{target_id}/{band_code}: degenerate stack stats, skipping band")
        return [], 0

    # first-pass median (pre-alignment), used as the registration reference AND
    # as the fixed reference sigma for the before/after fairness report below
    median_img0 = np.nanmedian(stack, axis=0)
    diff_before = stack - median_img0
    mad_before = np.nanmedian(np.abs(diff_before - np.nanmedian(diff_before)))
    sigma_ref = 1.4826 * mad_before

    # skip alignment entirely for a band with too few point sources to register against
    n_sources = count_point_sources(median_img0)
    skip_registration = n_sources < REGISTRATION_MIN_SOURCES
    if skip_registration:
        print(f"{target_id}/{band_code}: only {n_sources} point source(s) in the median image "
              f"(< {REGISTRATION_MIN_SOURCES}), skipping registration for the whole band")

    # sub-pixel align every frame to that reference (brightest 0.1% clipped so a
    # saturated core doesn't dominate the cross-correlation). A shift is only
    # trusted if it's small AND phase_cross_correlation is confident about it;
    # otherwise shift=(0,0) is used for that frame instead.
    ref_for_align = clip_bright(median_img0, ALIGN_CLIP_PERCENTILE)
    shifts, confident, aligned, aligned_valid, raw_mags, errors = [], [], [], [], [], []
    for img in bg_subbed:
        valid = ~np.isnan(img)  # true per-pixel coverage, before NaN gets filled in for alignment
        filled = np.nan_to_num(img, nan=np.nanmedian(img))
        if skip_registration:
            shift_yx, ok, raw_mag, error = np.zeros(2), True, 0.0, 0.0
        else:
            raw_shift, error, _ = phase_cross_correlation(ref_for_align, clip_bright(filled, ALIGN_CLIP_PERCENTILE),
                                                            upsample_factor=ALIGN_UPSAMPLE_FACTOR)
            raw_mag = float(np.hypot(*raw_shift))
            ok = raw_mag <= REGISTRATION_MAX_SHIFT_PX and error <= REGISTRATION_MAX_ERROR
            shift_yx = raw_shift if ok else np.zeros(2)
        shifts.append(shift_yx)
        confident.append(ok)
        raw_mags.append(raw_mag)
        errors.append(error)
        aligned.append(ndi_shift(filled, shift_yx, order=3))
        # carry validity through the same shift (nearest-neighbor, so it stays a clean
        # boolean) instead of letting the filled-in value pass as real data downstream
        aligned_valid.append(ndi_shift(valid.astype(float), shift_yx, order=0) > 0.5)
    shifts = np.array(shifts)
    confident = np.array(confident)
    raw_mags = np.array(raw_mags)
    errors = np.array(errors)
    pass_frac = 1.0 if skip_registration else confident.mean()
    shift_mag = np.hypot(shifts[:, 0], shifts[:, 1])
    print(f"{target_id}/{band_code}: sub-pixel shift (px) mean={shift_mag.mean():.3f} "
          f"std={shift_mag.std():.3f} max={shift_mag.max():.3f}, "
          f"{int((~confident).sum())}/{len(confident)} frames failed the registration guard (shift=0 used)")

    # band-level fallback: if too many frames failed the guard, force shift=0 for
    # EVERY frame instead of leaving a mix of shifted/unshifted frames in the stack,
    # and keep the diffs (WCS alignment only) instead of dropping them
    registration_mode = "skipped" if skip_registration else "aligned"
    if not skip_registration and pass_frac < REGISTRATION_MIN_PASS_FRAC:
        failed = ~confident
        shift_fail = failed & (raw_mags > REGISTRATION_MAX_SHIFT_PX)
        error_fail = failed & (errors > REGISTRATION_MAX_ERROR)
        print(f"{target_id}/{band_code}: registration guard failures -- "
              f"{int(shift_fail.sum())} over {REGISTRATION_MAX_SHIFT_PX}px shift, "
              f"{int(error_fail.sum())} over {REGISTRATION_MAX_ERROR} error "
              f"(median shift={np.median(raw_mags):.3f}px, median error={np.median(errors):.3f})")
        print(f"{target_id}/{band_code}: only {pass_frac:.0%} of frames passed "
              f"(< {REGISTRATION_MIN_PASS_FRAC:.0%}) -- falling back to WCS-only alignment "
              f"(shift=0) for the whole band")
        shifts = np.zeros_like(shifts)
        aligned = [np.nan_to_num(img, nan=np.nanmedian(img)) for img in bg_subbed]
        aligned_valid = [~np.isnan(img) for img in bg_subbed]
        registration_mode = "wcs_only"
    aligned_stack = np.array(aligned)
    valid_stack = np.array(aligned_valid)

    # soften the +/- dipoles sub-pixel alignment residuals leave at point sources
    if GAUSSIAN_SIGMA_DIFF > 0:
        aligned_stack = np.array([gaussian_filter(f, GAUSSIAN_SIGMA_DIFF) for f in aligned_stack])

    # static sky = per-pixel median over time (post-alignment); diff = each frame minus that.
    # re-mask invalid pixels (lost when they were filled in for alignment) so a frame that
    # doesn't cover a spot never contributes a fake value to the median or the diff
    valid_aligned_stack = np.where(valid_stack, aligned_stack, np.nan)
    median_img = np.nanmedian(valid_aligned_stack, axis=0)
    diff_stack = valid_aligned_stack - median_img
    sigma_display = 1.4826 * np.nanmedian(np.abs(diff_stack - np.nanmedian(diff_stack)))  # for PNG scaling only

    # fairness: both counts use the SAME (pre-alignment) sigma_ref, not each version's own
    exceed_before = int((np.abs(diff_before) > 5 * sigma_ref).sum()) if sigma_ref > 0 else None
    exceed_after = int((np.abs(diff_stack) > 5 * sigma_ref).sum()) if sigma_ref > 0 else None
    print(f"{target_id}/{band_code}: pixels exceeding +/-5*sigma_ref ({sigma_ref:.4f}) out of "
          f"{diff_stack.size}: before={exceed_before} after={exceed_after}")

    # background-only std: excludes pixels within ~3px of the brightest 1% of the median,
    # so a bright star's own residual doesn't dominate the number
    bg_exclude = binary_dilation(median_img0 >= np.nanpercentile(median_img0, 99), iterations=3)
    bg_only = ~bg_exclude
    print(f"{target_id}/{band_code}: background-only std: "
          f"before={np.nanstd(diff_before[:, bg_only]):.4f} after={np.nanstd(diff_stack[:, bg_only]):.4f}")

    make_diffs = np.isfinite(sigma_display) and sigma_display > 0
    bright_mask = None
    masked_count = 0
    if make_diffs:
        os.makedirs(diff_dir, exist_ok=True)
        # bright stars leave their own PSF/alignment residual in the diff, not a
        # moving/changing source -- mask them out to mid-grey instead of hiding that fact.
        # two tiers: an ordinary bright star gets a modest dilation, a saturating one
        # (much larger residual) gets a wider dilation so its halo is fully covered
        bright_thresh = np.nanpercentile(median_img, BRIGHT_MASK_PERCENTILE)
        bright_mask = binary_dilation(median_img >= bright_thresh, iterations=BRIGHT_MASK_DILATION_ITERATIONS)
        strong_thresh = np.nanpercentile(median_img, BRIGHT_MASK_PERCENTILE_STRONG)
        strong_mask = binary_dilation(median_img >= strong_thresh, iterations=BRIGHT_MASK_DILATION_ITERATIONS_STRONG)
        bright_mask = bright_mask | strong_mask
        masked_count = int(bright_mask.sum())
        print(f"{target_id}/{band_code}: masking {masked_count}/{bright_mask.size} pixels per frame "
              f"({masked_count / bright_mask.size:.1%}, >= {BRIGHT_MASK_PERCENTILE}th/"
              f"{BRIGHT_MASK_PERCENTILE_STRONG}th pct of the static sky)")
    else:
        print(f"{target_id}/{band_code}: diff sigma is zero/non-finite, skipping diff images")

    out_frames = []
    for (png_path, stable_id, time_iso, _), img_sub, diff_img, valid_i in zip(collected, bg_subbed, diff_stack,
                                                                                valid_stack):
        scaled = np.clip((img_sub - lo) / (hi - lo), 0, 1)
        scaled = np.arcsinh(10 * scaled) / np.arcsinh(10)   # brings out faint stuff
        scaled = np.nan_to_num(scaled, nan=0.0)
        scaled = np.flipud(scaled)                          # FITS has y flipped
        Image.fromarray((scaled * 255).astype(np.uint8)).save(png_path)

        entry = {"url": png_path, "time": time_iso, "band": band_full, "registration": registration_mode}
        if make_diffs:
            diff_path = f"{diff_dir}/{stable_id}.png"
            # -5sigma..+5sigma -> 0..255, zero diff at mid-grey, one scale for the whole stack
            diff_norm = np.clip(diff_img / (5 * sigma_display), -1, 1)
            diff_norm = np.nan_to_num(diff_norm, nan=0.0)
            diff_gray = (diff_norm + 1) / 2 * 255
            diff_gray[bright_mask] = 127.5  # mid-grey over bright-star pixels
            invalid_mask = binary_dilation(~valid_i, iterations=1)  # +1px so a smoothed edge doesn't leak through
            diff_gray[invalid_mask] = 127.5
            diff_gray = np.flipud(diff_gray)
            Image.fromarray(diff_gray.astype(np.uint8)).save(diff_path)
            entry["diff_url"] = diff_path
        out_frames.append(entry)
    return out_frames, masked_count


def process_target(target_id, ra, dec, label, size=0.2, max_frames=40, max_frames_other=20, bands=None,
                    offline=False):
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
    total_masked = 0
    for band_code in bands:
        band_full = f"SPHEREx-{band_code}"
        frames = results[results["energy_bandpassname"] == band_full]
        min_needed = MIN_BAND_FRAMES if band_code == "D2" else MIN_OTHER_BAND_FRAMES
        if len(frames) < min_needed:
            print(f"SKIP {target_id} {band_code}: only {len(frames)} frames (< {min_needed})")
            continue

        # offline (--diff-only) rebuilds use every frame already cached for this band,
        # instead of resampling to the cap -- otherwise a rebuild can lose frames whose
        # cached file doesn't match this run's sampled subset
        band_max = max_frames if band_code == "D2" else max_frames_other
        if not offline and len(frames) > band_max:
            # frames is already time-sorted; spread picks across the full span
            # instead of clustering near the start, so the blink viewer still
            # covers the whole observing window
            idx = np.unique(np.linspace(0, len(frames) - 1, band_max).round().astype(int))
            frames = frames[idx]
            print(f"{target_id}/{band_code}: sampled down to {len(frames)} frames across full time range")

        t0 = time.time()
        out_frames, masked_count = process_band(target_id, size_tag, band_code, coord, request_size, target_wcs,
                                                  frames, offline=offline)
        print(f"{target_id}/{band_code}: took {time.time() - t0:.1f}s, masked {masked_count} pixels")
        total_masked += masked_count
        if out_frames:
            bands_out[band_full] = out_frames

    print(f"{target_id}: {total_masked} pixels masked total across {len(bands_out)} bands")
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


def merge_entry(existing, new):
    """Merge a run's entry into the saved one, per band: bands this run produced
    replace the saved ones, bands it didn't touch (e.g. D2 on a --bands D1,D3 run)
    are kept. A band that already has frames is never dropped."""
    if not existing:
        return new
    merged = {**existing, **{k: v for k, v in new.items() if k not in ("bands", "band_ranges", "frames")}}
    bands = dict(existing.get("bands", {}))
    for band, frames in new.get("bands", {}).items():
        if frames:
            bands[band] = frames
    merged["bands"] = bands
    merged["band_ranges"] = {**existing.get("band_ranges", {}), **new.get("band_ranges", {})}
    # top-level "frames" is the D2 view the page falls back on
    merged["frames"] = bands.get("SPHEREx-D2", existing.get("frames", []))
    return merged


def save_target_entry(entry):
    existing = next((t for t in load_targets() if t["id"] == entry["id"]), None)
    targets = [t for t in load_targets() if t["id"] != entry["id"]]
    targets.append(merge_entry(existing, entry))
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
#   python prep/test_one.py --only pleiades --full       D2 up to 40 frames, other bands up to 20
#   python prep/test_one.py --only pleiades --full --max-frames-other 10   fewer non-D2 frames
#   python prep/test_one.py --only pleiades --bands D1,D2 --full   only those bands
#   python prep/test_one.py --full                      all targets, all bands
#   python prep/test_one.py --only pleiades --bands D2 --full --diff-only   diffs from cache, no downloads
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

    max_frames_other = 20
    explicit_mfo = "--max-frames-other" in sys.argv
    if explicit_mfo:
        max_frames_other = int(sys.argv[sys.argv.index("--max-frames-other") + 1])

    kwargs = {"bands": bands, "offline": "--diff-only" in sys.argv}
    if "--full" in sys.argv:
        kwargs["max_frames_other"] = max_frames_other
    else:
        # quick test runs stay uniform across bands unless the user overrides it
        kwargs["max_frames"] = TEST_MAX_FRAMES
        kwargs["max_frames_other"] = max_frames_other if explicit_mfo else TEST_MAX_FRAMES

    for target_id, ra, dec, label in TARGETS:
        if only and target_id != only:
            continue
        entry = process_target(target_id, ra, dec, label, **kwargs)
        if entry:
            n = save_target_entry(entry)
            print(f"saved {target_id}, targets.json now has {n} targets")
