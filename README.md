# SPHEREx Sky Blink

Made for the NASA Space Apps Challenge 2026 (SPHEREx challenge).
**Live site:** https://wonderful-douhua-5dfc7e.netlify.app/

SPHEREx photographs the whole sky again and again. This site shows the same patch of sky on different dates, and visitors flip through the photos and report whether a bright spot moved, changed brightness, or nothing changed.

## What it does

- Eight fields to explore: the Pleiades, M101, the Beehive Cluster (M44), two points on the ecliptic, the Galactic Center region, and the north and south ecliptic poles.
- A blink player: play, pause, step and scrub through the frames of a field, in order of date.
- Six wavelength bands (Band 1 to Band 6, which are SPHEREx detectors D1 to D6), with a short description of what each band shows.
- An experimental Difference view, which subtracts the static sky so changes stand out.
- A detector-flag overlay that tints pixels SPHEREx's own pipeline flagged (saturation, cosmic rays, bad pixels, ghosts), so visitors don't report artifacts.
- "Inspect a spot": click a star to get its sky coordinates, links to SIMBAD and WiseView, and a live check against SIMBAD.
- Hover coordinates, a scale bar and an N/E compass on the image.
- A reference image of the same patch of sky from 2MASS, for comparison.
- Report buttons, an optional comment, and a moderated Discussion tab for each field.
- Search by name or coordinates, Random field, and shareable links.

## How it works

- A Python pipeline (`prep/`) downloads SPHEREx image cutouts from NASA/IPAC IRSA, aligns the frames onto one sky grid, applies one shared brightness stretch per band, builds difference images, and makes the detector-flag overlays and 2MASS reference images. It writes PNG files and JSON describing each field.
- The site is a static page (HTML, CSS and JavaScript) hosted on Netlify. Because the images are prepared in advance, it loads quickly and does not depend on the archive being online.
- Visitors' answers and comments are stored in a Supabase database. Comments are reviewed by the team before they appear.

## Run it locally

See `docs/RUN_LOCALLY.md`. In short: `python serve.py`, then open http://localhost:8000. The database setup is in `docs/setup.sql`.

## Team

- Arjun Ninad Bendkhale with name: Python pipeline and deployment
- Susab Singh and Iqra Mirza with name: design, README and field descriptions

## Data and credits

- SPHEREx images: NASA/IPAC Infrared Science Archive (IRSA).
- Reference images: 2MASS (Two Micron All Sky Survey), via NASA/IPAC IRSA.
- Catalogue lookups: SIMBAD, CDS Strasbourg (with its Harvard mirror).
- Field descriptions were written by the team from these sources:
  - https://science.nasa.gov/mission/hubble/science/explore-the-night-sky/hubble-messier-catalog/messier-45/
  - https://science.nasa.gov/mission/hubble/science/explore-the-night-sky/hubble-messier-catalog/messier-101/
  - https://science.nasa.gov/asset/hubble/hubbles-largest-galaxy-portrait-offers-a-new-high-definition-view/
  - https://science.nasa.gov/mission/webb/science-overview/science-explainers/what-is-the-center-of-our-galaxy-like/
  - https://galacticcenter.astro.ucla.edu/learn-about-our-galaxy.html
  - https://esahubble.org/images/heic1606a/
  - https://www.sciencedirect.com/topics/physics-and-astronomy/galactic-center
  - http://hyperphysics.phy-astr.gsu.edu/hbase/eclip.html
