# Use of AI

We used Claude (Anthropic) throughout this project, through the Claude chat app and the Claude extension for VS Code.

## What the team did

- Chose the idea, decided what to build, and decided what to cut.
- Ran the pipeline, deployed the site (GitHub, Netlify, Supabase), and tested it in the browser and on other people's screens.
- Wrote the description and fun fact for each field, with the sources listed in the About tab.

## What Claude helped with

- The Python pipeline: downloading SPHEREx cutouts from IRSA, aligning frames, the shared brightness stretch, difference images, detector-flag masks and the 2MASS reference images.
- The web page (HTML, CSS and JavaScript), its layout, and the Supabase setup.
- Drafts of the band descriptions, the per-band feature lists and the README structure, which the team reviewed and edited.
- Suggestions for features and for what to leave out.

## What we checked

- Frame alignment, centroid scatter and frame counts per field and band, with results reviewed by the team.
- That the 2MASS reference images match the orientation of the SPHEREx frames in all eight fields.
- That the page's hover and click coordinates match the pipeline's own sky coordinates (worst difference about 0.07 arcseconds).
- That the page works at several window sizes.

## Limits

- The Difference view is experimental.
- The per-band feature lists are indicative, and were not reviewed by an astronomer.
- AI can make mistakes, and we may have missed some.

## Data and credits

SPHEREx data via NASA/IPAC IRSA. Reference images from 2MASS via NASA/IPAC IRSA. Catalogue lookups use SIMBAD (CDS). Sources for the field descriptions are listed in the About tab.
