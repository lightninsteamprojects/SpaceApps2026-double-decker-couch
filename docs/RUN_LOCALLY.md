# Running SPHEREx Sky Blink locally

Commands are for a terminal opened in the project folder (the one containing `index.html`).
On Windows use `python`; on macOS/Linux use `python3` if `python` isn't found.

## 1. Run the site

```
python serve.py
```

Then open <http://localhost:8000>. Stop with Ctrl+C. `serve.py` only uses Python's standard library and tells the
browser not to cache anything, so after editing `index.html` a normal refresh shows the change.
(Use `python serve.py 8080` for another port.) Opening `index.html` by double-click does **not** work: the page loads
`targets.json` and other files with `fetch`, which browsers block on `file://` addresses.

## 2. Set up the data pipeline (only needed to make or update data)

```
python -m venv .venv
.venv\Scripts\activate          # Windows.  macOS/Linux:  source .venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` covers every package imported by `prep/test_one.py`, `prep/make_flags.py` and
`prep/fetch_reference.py` (astropy, astroquery, numpy, scipy, Pillow, requests, reproject, s3fs, scikit-image).
The pipeline needs internet access: it downloads SPHEREx cutouts from NASA/IPAC IRSA.

## 3. Add a new field

**a) Add it to the list.** In `prep/test_one.py`, add one line to `TARGETS` (id, RA in degrees, Dec in degrees, label):

```python
    ("orion", 83.82, -5.39, "Orion Nebula (M42)"),
```

The id is used for folder names and URLs, so keep it lowercase with no spaces.

**b) Download and process it** (all six bands, up to 40 frames for D2 and 20 for the others):

```
python prep/test_one.py --only orion --full
```

Other options (there is no `--help`; they are listed in the "Usage" comment just above `if __name__ == "__main__":`
in `prep/test_one.py`):

| option | what it does |
|---|---|
| `--count` | show frames available per band, no downloads |
| `--only ID` | process just that target |
| `--full` | D2 up to 40 frames, other bands up to 20 (without it: a quick 5-frame test) |
| `--bands D1,D2` | only those bands |
| `--max-frames-other N` | cap on frames for the non-D2 bands |
| `--diff-only` | rebuild difference images from cached data, no downloads |
| `--clean` | delete `data/`, `img/` and `targets.json` (asks first) |

**Careful:** running `prep/test_one.py` with no `--only` processes *every* target. This writes into `data/`,
`img/` and `targets.json`.

If IRSA answers with a 502 or a timeout, wait a few minutes and run the same command again; frames already
downloaded are cached.

**c) Detector-flag overlays** (reads the cached data, downloads nothing, skips frames that already have one):

```
python prep/make_flags.py
```

**d) Reference images** (2MASS images via NASA/IPAC IRSA, one per field in `targets.json`; this redoes all
fields and takes a few minutes):

```
python prep/fetch_reference.py
```

Optionally add a description for the new field in `fields.json` (`"blurb"` and `"fact"`).
Reload the page and the new field appears in the Browse list.

## 4. Use your own Supabase project

By default the site stores **answers** and **Discussion comments** in the team's Supabase project, so anything you
submit while testing goes there. To use your own:

1. Create a free project at <https://supabase.com>.
2. In the dashboard open **SQL Editor**, paste the contents of `docs/setup.sql` and run it. It creates the
   `comments` and `classifications` tables with Row Level Security. (Check the `classifications` part against the
   live table first, as noted at the top of that file.)
3. In **Project Settings -> API** copy the project URL and the publishable (anon) key.
4. In `index.html`, search for `supabase.createClient(`. The two values passed to it are the project URL and the
   publishable key. Replace both with yours. (Note the URL has no `/rest/v1/` on the end.)
5. Comments only show up on the site after you approve them: in **Table Editor -> comments**, tick `approved`.

The publishable key is meant to be visible in the page; what visitors can do is limited by the policies in
`docs/setup.sql`. Never put a `service_role` key in `index.html`.
