"""Build transparent red overlay PNGs from each cached cutout's FLAGS extension,
showing only instrument/detector problem bits (not source bits). Read-only against
the existing data/ cache -- never downloads anything. Re-runnable: a frame whose
flag PNG already exists, or whose FITS cutout has no FLAGS extension, is skipped.

Bit selection approved by the user; see the chat for the reasoning per bit.
Values are the FITS header's own MSKMnnnn bitmasks, not guessed.
"""
import glob
import json
import os
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from PIL import Image
from reproject import reproject_interp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from prep.test_one import TARGETS, OUT_SHAPE, build_target_wcs

RA_DEC = {t[0]: (t[1], t[2]) for t in TARGETS}
SIZE = 0.2  # every cached target uses this size (folder tag "0p2")

# instrument/detector problem bits only -- excludes DICHROIC (permanent low-efficiency
# region by design), FULLSAMPLE (informational, not a problem), PHANMISS/NONLINEAR
# (calibration-step applicability caveats, not declared defects), and SOURCE (marks a
# real astronomical source, not a defect)
INCLUDED_BITS = {
    "TRANSIENT": 1,
    "OVERFLOW": 2,
    "SUR_ERROR": 4,
    "NONFUNC": 64,
    "MISSING_DATA": 512,
    "HOT": 1024,
    "COLD": 2048,
    "PERSIST": 131072,
    "OUTLIER": 524288,
    "CROSSTALK": 1048576,
    "GHOST": 4194304,
    "GHOST_FPA": 8388608,
    "GHOST_EXT": 16777216,
    "STREAK": 33554432,
    "BLOOM": 67108864,
    "SNOWBALL": 134217728,
    "HALO": 268435456,
}
COMBINED_MASK = 0
for _bit in INCLUDED_BITS.values():
    COMBINED_MASK |= _bit

FLAGS_JSON = "flags.json"


def load_flags_json():
    if os.path.exists(FLAGS_JSON):
        with open(FLAGS_JSON) as f:
            return json.load(f)
    return {}


def main():
    flags_out = load_flags_json()
    no_flags_ext = []  # (target, band, stable_id) with no FLAGS extension at all

    # target -> band -> stats accumulated over frames analyzed in THIS run
    stats = {}

    for fits_path in sorted(glob.glob("data/*/0p2/*/*.fits")):
        parts = fits_path.replace("\\", "/").split("/")
        target_id, size_tag, band_code = parts[1], parts[2], parts[3]
        stable_id = os.path.splitext(parts[4])[0]

        flag_dir = f"img/{target_id}/{size_tag}/{band_code}/flags"
        flag_path = f"{flag_dir}/{stable_id}.png"

        band_stats = stats.setdefault(target_id, {}).setdefault(band_code, {
            "processed": 0, "skipped_exists": 0, "skipped_no_flags": 0,
            "skipped_empty": 0, "bit_counts": {name: 0 for name in INCLUDED_BITS},
            "combined_counts": 0,
        })

        if os.path.exists(flag_path):
            band_stats["skipped_exists"] += 1
            continue

        try:
            with fits.open(fits_path) as hdul:
                if "FLAGS" not in hdul:
                    band_stats["skipped_no_flags"] += 1
                    no_flags_ext.append((target_id, band_code, stable_id))
                    continue
                flags_data = hdul["FLAGS"].data.astype(np.int64)
                flags_wcs = WCS(hdul["FLAGS"].header)
        except OSError:
            # the all-bands run may still be writing this file -- treat as not
            # cached yet rather than crashing; a later re-run will pick it up
            continue

        ra, dec = RA_DEC[target_id]
        target_wcs = build_target_wcs(ra, dec, SIZE, npix=OUT_SHAPE[0])
        reproj, footprint = reproject_interp(
            (flags_data, flags_wcs), target_wcs, shape_out=OUT_SHAPE, order="nearest-neighbor"
        )
        reproj = np.nan_to_num(reproj, nan=0.0).astype(np.int64)
        covered = footprint > 0

        combined_flagged = (reproj & COMBINED_MASK) != 0
        combined_flagged &= covered
        n_flagged = int(combined_flagged.sum())

        for name, bit in INCLUDED_BITS.items():
            band_stats["bit_counts"][name] += int(((reproj & bit) != 0 & covered).sum())
        band_stats["combined_counts"] += n_flagged
        band_stats["processed"] += 1

        if n_flagged == 0:
            band_stats["skipped_empty"] += 1
            continue

        os.makedirs(flag_dir, exist_ok=True)
        rgba = np.zeros((*OUT_SHAPE, 4), dtype=np.uint8)
        rgba[combined_flagged] = [255, 0, 0, 255]
        rgba = np.flipud(rgba)  # match the frame PNGs: FITS has y flipped
        Image.fromarray(rgba, mode="RGBA").save(flag_path)

        flags_out.setdefault(target_id, {}).setdefault(band_code, {})[stable_id] = flag_path

    with open(FLAGS_JSON, "w") as f:
        json.dump(flags_out, f, indent=2)

    # ----- report -----
    total_pixels = OUT_SHAPE[0] * OUT_SHAPE[1]
    print(f"{'target':<16}{'band':<6}{'processed':<11}{'skip_exist':<11}{'skip_noflags':<13}"
          f"{'skip_empty':<11}{'avg_flagged%':<13}")
    for target_id in sorted(stats):
        for band_code in sorted(stats[target_id]):
            s = stats[target_id][band_code]
            n = s["processed"]
            avg_pct = (s["combined_counts"] / (n * total_pixels) * 100) if n else 0.0
            print(f"{target_id:<16}{band_code:<6}{n:<11}{s['skipped_exists']:<11}"
                  f"{s['skipped_no_flags']:<13}{s['skipped_empty']:<11}{avg_pct:<13.3f}")

    print("\nper-bit average flagged fraction (% of pixels), across frames processed this run:")
    bit_totals = {name: 0 for name in INCLUDED_BITS}
    n_total = 0
    for target_id in stats:
        for band_code in stats[target_id]:
            s = stats[target_id][band_code]
            n_total += s["processed"]
            for name in INCLUDED_BITS:
                bit_totals[name] += s["bit_counts"][name]
    if n_total:
        for name in INCLUDED_BITS:
            pct = bit_totals[name] / (n_total * total_pixels) * 100
            print(f"  {name:<14}{pct:.4f}%")
    else:
        print("  (no frames processed this run)")

    if no_flags_ext:
        print(f"\n{len(no_flags_ext)} frame(s) skipped for having no FLAGS extension:")
        for t, b, s in no_flags_ext:
            print(f"  {t}/{b}/{s}")


if __name__ == "__main__":
    main()
