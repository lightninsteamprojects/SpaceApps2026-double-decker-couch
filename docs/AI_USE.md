# AI use in this project

NASA Space Apps Challenge 2026 asks teams to describe the AI tools they used. This file lists them and what each was used for. Edit any line that is not accurate for your team.

## Tools

- **Claude (Anthropic), in the Claude chat app.** Used for planning the project, comparing approaches (static site plus Supabase, which bands to use, how to align frames), reviewing example Space Apps repositories, drafting prompts, explaining errors, and reviewing results.
- **Claude, in the VS Code extension.** Used to write most of the code: the Python data pipeline (`prep/test_one.py`) and the web page (`index.html`), working from prompts written by the team. [CONFIRM: add anything else, e.g. other AI tools]

## What the team did

- Chose the challenge approach: a public web tool where visitors flip through SPHEREx images of the same patch of sky and report what moved or changed.
- Chose the target fields, the band choices and the overall design.
- Set up the data accounts and services: the Supabase database and security rules, GitHub, and Netlify hosting.
- Ran the pipeline on our own computer, and read and tested every change on the live site.
- Made the decisions about what to keep and what to drop (for example, removing the colour toggle and labelling the difference view as experimental). [CONFIRM]
- [Friend: design and styling of the page, README, testing on phones. CONFIRM]

## What AI helped write

- The Python pipeline: fetching SPHEREx cutouts from NASA/IPAC IRSA, retries, reprojecting frames onto a shared grid, brightness scaling, difference images and masks.
- The web page: layout, blink viewer, band selector, sharing buttons, comment box and saving answers.
- Documentation drafts.

## What we checked ourselves

- Alignment: cross-correlation of 26 frames of one field and 5 of another showed a shift of 0 pixels (whole-pixel precision only). Later, measuring the centroid of the brightest star in the Pleiades field showed about 0.3 pixels of scatter (up to about 1.3 pixels), which is why difference images show residuals around bright stars.
- Frame counts per field and band, and how many frames were dropped at the detector edge.
- That the same wavelength band is used within a field, so brightness comparisons make sense.
- The live site: loading each field, saving an answer, and checking the row in the database.
- Problems we found and fixed or documented: frames cut off at the detector edge, a black wedge in a Pleiades difference image caused by a very bright star, and registration that could not lock on in faint fields.

## Data and credits

- Images: SPHEREx spectral images from the NASA/IPAC Infrared Science Archive (IRSA).
- Reference images (comparison card in the Browse tab): 2MASS K-band images via NASA/IPAC IRSA, fetched by `prep/fetch_reference.py` (written with AI help), cropped to the same 0.2 degree field and arcsinh-stretched. NASA SkyView was tried first but its server could not reach the 2MASS archive, so it is not the source.
- Libraries: astropy, astroquery, reproject, numpy, scipy, Pillow, s3fs.
- Hosting: Netlify. Database: Supabase.

## Limits we know about

- Difference images are experimental. Bright stars leave residuals, and the telescope's pixels are large.
- Only some fields have all six wavelength bands. The others have one.
- A "user" in our consensus counts is a browser, so one person could fake agreement.
