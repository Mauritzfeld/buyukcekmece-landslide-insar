# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A single-script scientific analysis pipeline that reproduces the data-derived figures for the
manuscript *"Satellite radar reveals rotational self-stabilization in earthquake-triggered
landslides"* (`manusript.tex`). It studies the Büyükçekmece landslide near Istanbul before and
after the 2019 Mw 5.7 Silivri earthquake using multi-geometry Sentinel-1 Persistent Scatterer
(PS) InSAR.

All code lives in `generate_figures.py`. See `README.md` for user-facing install/run docs.

## Running

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python generate_figures.py        # inputs resolved relative to the script (data/); runs from any CWD
```

The script has an `if __name__ == '__main__'` guard and runs headless (matplotlib `Agg`
backend). `main()` produces exactly four figures into the working directory:

- `landslide_analysis__figure1.pdf` / `.png` — Fig. 1, velocity time series
- `landslide_analysis__figure2.pdf` / `.png` — Fig. 2, phase-space trajectories
- `Acquisition_dates.pdf` — Fig. S3, acquisition timeline
- `synthetic_ps_locations.png` / `.pdf` — Fig. S4, synthetic PS locations

The reference copies as used in the manuscript live in `figures_manuscript/`. The remaining
manuscript figures (conceptual model, study-area map, geological cross-section, conceptual
phase-space diagram) are made outside this code and are not generated here.

## Architecture

`generate_figures.py` is organized top-to-bottom as: module docstring → config constants →
`load_data()` → `find_closest_points()` → `decompose()` + `average_components()` → the four
figure functions → `main()`.

The key design point: **all four figures derive from one computation.** `decompose(df_ascending,
df_descending)` pairs ascending/descending PS points by geodesic nearest neighbour, decomposes
line-of-sight velocity into horizontal/vertical components (fixed incidence angles — descending
43.8154°, ascending 39.2606° — and azimuth projection factors in the `_elosd/_hlosd/_elosa/_hlosa`
helpers), and returns the per-point velocity time-series lists plus the synthetic-point
(midpoint) frame. From there:

- Figs 1–2 use `average_components()` (cross-point mean/std per acquisition; the vertical mean is
  sign-flipped so positive = uplift).
- Fig. S4 attaches per-point medians to the synthetic-point locations and maps them.

## Data and important data quirk

Inputs live in `data/`. Excel files: `data/TSbulkdSA.xlsx`/`TSbulkaSA.xlsx` (post-earthquake
descending/ascending) and `data/newTS160dCSV.xlsx`/`newTS160aCSV.xlsx` (pre-earthquake
descending/ascending). Each row is a PS point with `LON`, `LAT`, `ID`. Shapefiles
`data/toe_boundary_QGIS.shp` (analysis area) and `data/shape_landslide_boundary.shp` (full extent)
define the regions; points are sliced with `gdf.within(...)`. Paths are resolved from a `DATA_DIR`
constant (`Path(__file__).parent / "data"`).

**Date-column quirk (do not "fix"):** each spreadsheet carries the acquisition dates *twice* — once
as `datetime`-typed column headers and once as string headers. `decompose()` deliberately uses the
`datetime` columns for the time axis (`isinstance(c, datetime)`) and the *string* columns for the
velocity values (`pd.to_datetime(c, errors='ignore') != c`). Both selections must stay; changing
either breaks the alignment and the numbers.

## Constraints / gotchas

- `pd.to_datetime(..., errors='ignore')` is used for date-column detection. It is **removed in
  pandas 3.0**, so `requirements.txt` caps `pandas<3`. A `FutureWarning` is silenced at import.
- Fig. S4 downloads an Esri World Imagery basemap via `contextily` at run time — needs internet,
  and the background tiles may differ over time (plotted points/boundaries are unaffected).
- When editing the figure code, note the velocity-axis sign conventions: horizontal is plotted as
  `-average_V_hor` and the vertical mean is pre-negated in `average_components()`.
