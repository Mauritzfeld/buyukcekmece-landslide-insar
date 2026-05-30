# Rotational self-stabilization in earthquake-triggered landslides — figure code

Code and data to reproduce the data-derived figures of the manuscript

> **Satellite radar reveals rotational self-stabilization in earthquake-triggered landslides**

The study analyses the Büyükçekmece landslide (western Istanbul, Turkey) before
and after the 2019 Mw 5.7 Silivri earthquake using multi-geometry Sentinel-1
Persistent Scatterer (PS) InSAR. Line-of-sight velocities from ascending and
descending geometries are paired by geodesic nearest neighbour and decomposed
into horizontal and vertical components, then analysed as velocity time series
and phase-space trajectories.

## Repository contents

| Path | Description |
|---|---|
| `generate_figures.py` | Single script that produces all data-derived figures. |
| `requirements.txt` | Python dependencies (tested with Python 3.9). |
| `CITATION.cff` | Citation metadata for this software/data deposit. |
| `LICENSE` | MIT license. |
| `data/TSbulkaSA.xlsx`, `data/TSbulkdSA.xlsx` | Post-earthquake PS time series (ascending / descending). |
| `data/newTS160aCSV.xlsx`, `data/newTS160dCSV.xlsx` | Pre-earthquake PS time series (ascending / descending). |
| `data/toe_boundary_QGIS.shp` (+ sidecars) | Analysis area (landslide toe) used to slice PS points. |
| `data/shape_landslide_boundary.shp` (+ sidecars) | Full landslide extent (for the location map). |
| `figures_manuscript/` | Reference copies of the figures as used in the manuscript. |

Each row of the time-series spreadsheets is a PS point with `LON`, `LAT`, `ID`
and the per-acquisition displacement values in dated columns. Shapefiles keep
all of their sidecar files (`.shp`/`.shx`/`.dbf`/`.prj`/`.cpg`/`.qmd`) together.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Running

```bash
python generate_figures.py
```

Inputs are resolved relative to the script's own location (`data/`), so it can
be run from any working directory. The figures are written to the current
directory.

## Outputs and manuscript figures

| Output file | Manuscript figure |
|---|---|
| `landslide_analysis__figure1.pdf` / `.png` | Fig. 1 — velocity time series |
| `landslide_analysis__figure2.pdf` / `.png` | Fig. 2 — phase-space trajectories |
| `Acquisition_dates.pdf` | Fig. S3 — Sentinel-1 acquisition timeline |
| `synthetic_ps_locations.png` / `.pdf` | Fig. S4 — synthetic PS locations |

The synthetic-PS map (Fig. S4) overlays the points on an Esri World Imagery
basemap downloaded at run time via `contextily`, so generating it **requires an
internet connection**. Because that basemap is a live tile service, the
satellite imagery in the background may differ slightly from the committed copy
over time; the plotted point locations and boundaries are unaffected.

### Figures not produced by this code

The following manuscript figures are created outside this pipeline (GIS /
illustration software, or reproduced from a cited publication) and are **not**
generated here: the conceptual model (Fig. 3), the study-area location map
(Fig. S2), the geological cross-section (Fig. S1, from Martino et al. 2018), and
the conceptual phase-space diagram (Fig. S5).

## Method summary

1. Load the four spreadsheets and build point GeoDataFrames; slice to the
   landslide-toe analysis area defined by `toe_boundary_QGIS.shp`.
2. Pair ascending and descending PS points by geodesic nearest neighbour
   (`find_closest_points`).
3. Decompose line-of-sight velocity into horizontal and vertical components
   using fixed incidence angles (descending 43.8154°, ascending 39.2606°) and
   azimuth-dependent projection factors (`decompose`).
4. Average across points per acquisition for the time-series and phase-space
   figures; use per-point medians at the synthetic-point locations for the map.

## Creating the archive deposit

The public code/data deposit contains only the source, inputs and reference
figures. To (re)build the upload archive from the repository root:

```bash
zip -r buyukcekmece_landslide_insar_v1.0.zip \
    README.md LICENSE CITATION.cff requirements.txt generate_figures.py \
    data figures_manuscript \
    -x '*/.DS_Store' '*/__pycache__/*'
```

This deliberately excludes editing/working files (e.g. the manuscript source
and editor configuration). After uploading to Zenodo, record the assigned DOI
in `CITATION.cff` and cite the version-independent (concept) DOI in the
manuscript.

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
