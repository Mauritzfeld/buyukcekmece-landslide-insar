#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure generation for the manuscript

    "Satellite radar reveals rotational self-stabilization in
     earthquake-triggered landslides"

This script reproduces the four figures in the manuscript that are produced
from the multi-geometry Sentinel-1 Persistent Scatterer (PS) InSAR data of the
Büyükçekmece landslide (Istanbul, Turkey) before and after the 2019 Mw 5.7
Silivri earthquake:

    landslide_analysis__figure1.pdf / .png   Fig. 1  velocity time series
    landslide_analysis__figure2.pdf / .png   Fig. 2  phase-space trajectories
    Acquisition_dates.pdf                     Fig. S3 acquisition timeline
    synthetic_ps_locations.png / .pdf         Fig. S4 synthetic PS locations

The remaining manuscript figures (conceptual model, study-area map, geological
cross-section and the conceptual phase-space diagram) are not produced here;
they are created outside this pipeline (GIS / illustration software / a cited
publication).

Method
------
Line-of-sight (LOS) velocities from ascending and descending geometries are
paired by geodesic nearest neighbour and decomposed into horizontal and
vertical components using fixed incidence angles (descending 43.8154°,
ascending 39.2606°) and azimuth-dependent projection factors. The per-point
time series are then averaged across points (Figs 1-2) or reduced to per-point
medians attached to the synthetic-point locations (Fig. S4).

Usage
-----
    pip install -r requirements.txt
    python generate_figures.py

Outputs are written to the current working directory. The synthetic-PS map
downloads an Esri World Imagery basemap via ``contextily`` and therefore
requires an internet connection.
"""

import math
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")  # headless / non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.lines import Line2D

import geopandas as gpd
from shapely.geometry import Point
from shapely.ops import unary_union
from geopy.distance import geodesic
from scipy.optimize import curve_fit
import contextily as ctx

# ``pd.to_datetime(col, errors='ignore')`` is used below to detect the (string)
# date columns; it is deprecated in pandas >= 2.2 but still functional. Silence
# the FutureWarning to keep the console output clean.
warnings.filterwarnings("ignore", category=FutureWarning)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Inputs live in data/ next to this script, so the script can be run from any
# working directory. Output figures are written to the current directory.
DATA_DIR = Path(__file__).resolve().parent / "data"

# Input time-series spreadsheets. Each row is a PS point with LON/LAT/ID and
# many date columns; the files carry the date columns twice (once as datetime
# headers used for the time axis, once as string headers carrying the values).
AFTER_DESCENDING = DATA_DIR / "TSbulkdSA.xlsx"    # post-earthquake, descending
AFTER_ASCENDING = DATA_DIR / "TSbulkaSA.xlsx"     # post-earthquake, ascending
BEFORE_ASCENDING = DATA_DIR / "newTS160aCSV.xlsx"  # pre-earthquake, ascending
BEFORE_DESCENDING = DATA_DIR / "newTS160dCSV.xlsx"  # pre-earthquake, descending

# Region-of-interest shapefiles.
TOE_BOUNDARY = DATA_DIR / "toe_boundary_QGIS.shp"          # analysis area (landslide toe)
LANDSLIDE_BOUNDARY = DATA_DIR / "shape_landslide_boundary.shp"  # full landslide extent

# LOS-decomposition geometry (radians).
INCIDENCE_D = math.radians(43.8154)
INCIDENCE_A = math.radians(39.2606)
AZIMUTH_D = math.radians(-9.3835 + 180)
AZIMUTH_A = math.radians(-169.99429)

# Acquisition index separating the roto-translational phase from the
# subsequent steady-translation phase (post-earthquake series).
ROTATION_INDEX = 135

# Output filename prefix for the two main-text figures.
SAVE_PREFIX = "landslide_analysis"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _to_gdf(df):
    """Build a GeoDataFrame of point geometries from a frame's LON/LAT."""
    geometry = [Point(lon, lat) for lon, lat in zip(df["LON"], df["LAT"])]
    return gpd.GeoDataFrame(df, geometry=geometry)


def load_data():
    """Load the four time-series spreadsheets and the two shapefiles, then
    slice the PS points to the landslide-toe analysis area.

    Returns a dict with the four sliced (ascending/descending, before/after)
    GeoDataFrames, the raw frames needed for the acquisition timeline, and the
    boundary shapefiles needed for the synthetic-PS map.
    """
    after_d = pd.read_excel(AFTER_DESCENDING)
    after_a = pd.read_excel(AFTER_ASCENDING)
    before_a = pd.read_excel(BEFORE_ASCENDING, header=0)
    before_d = pd.read_excel(BEFORE_DESCENDING, header=0)

    gdf_after_d = _to_gdf(after_d)
    gdf_after_a = _to_gdf(after_a)
    gdf_before_a = _to_gdf(before_a)
    gdf_before_d = _to_gdf(before_d)

    toe_shape = gpd.read_file(TOE_BOUNDARY)
    combined_toe = unary_union(list(toe_shape.geometry))
    landslide_shape = gpd.read_file(LANDSLIDE_BOUNDARY)

    return {
        # sliced points within the toe analysis area
        "after_descending": gdf_after_d[gdf_after_d.within(combined_toe)],
        "after_ascending": gdf_after_a[gdf_after_a.within(combined_toe)],
        "before_ascending": gdf_before_a[gdf_before_a.within(combined_toe)],
        "before_descending": gdf_before_d[gdf_before_d.within(combined_toe)],
        # raw frames for the acquisition timeline (full series, not sliced)
        "raw_after_descending": after_d,
        "raw_after_ascending": after_a,
        "raw_before_ascending": before_a,
        "raw_before_descending": before_d,
        # boundaries for the synthetic-PS map
        "toe_shape": toe_shape,
        "landslide_shape": landslide_shape,
    }


# ---------------------------------------------------------------------------
# Point pairing and LOS velocity decomposition
# ---------------------------------------------------------------------------

def find_closest_points(gdf1, gdf2):
    """Pair each point in ``gdf1`` with its geodesically nearest point in
    ``gdf2``. Returns a list of (id1, id2) pairs, 1-based to match the data."""
    closest_points = []

    for index1, row1 in gdf1.iterrows():
        point1 = (row1["LAT"], row1["LON"])
        closest_point = None
        min_distance = float("inf")

        for index2, row2 in gdf2.iterrows():
            point2 = (row2["LAT"], row2["LON"])
            distance = geodesic(point1, point2).meters

            if distance < min_distance:
                min_distance = distance
                closest_point = (index1, index2)

        if closest_point is not None:
            closest_points.append(closest_point)

    adjusted_closest_points = [(index1 + 1, index2 + 1) for index1, index2 in closest_points]
    return adjusted_closest_points


def _elosd(incidence_d, azi_d):
    return math.cos(1.571 - incidence_d) * math.cos(4.712 - azi_d)


def _hlosd(incidence_d):
    return math.cos(incidence_d)


def _elosa(incidence_a, azi_a):
    return math.cos(1.571 - incidence_a) * math.cos(4.712 - azi_a)


def _hlosa(incidence_a):
    return math.cos(incidence_a)


def decompose(df_ascending, df_descending):
    """Pair ascending/descending PS points and decompose their LOS velocities
    into horizontal and vertical components.

    Returns
    -------
    date_datetime : pandas.DatetimeIndex
        Acquisition dates (the per-increment time axis).
    v_hor_list, v_vert_list : list of arrays
        Per synthetic-point horizontal / vertical velocity time series.
    df_synthetic : pandas.DataFrame
        Synthetic points (midpoints of each ascending/descending pair) with
        LAT, LON and a Synthetic_ID column.
    """
    new_ascending = gpd.GeoDataFrame(columns=df_ascending.columns)
    new_descending = gpd.GeoDataFrame(columns=df_descending.columns)

    midpoints = []
    synthetic_ids = []

    closest_points = find_closest_points(
        df_ascending[["LAT", "LON", "ID"]],
        df_descending[["LAT", "LON", "ID"]],
    )

    for pair_index, (asc_id, desc_id) in enumerate(closest_points):
        asc_id = int(asc_id)
        desc_id = int(desc_id)

        match_asc = df_ascending[df_ascending["ID"].astype(int) == asc_id].copy()
        match_desc = df_descending[df_descending["ID"].astype(int) == desc_id].copy()

        mid_lat = np.mean([match_asc["LAT"].values[0], match_desc["LAT"].values[0]])
        mid_lon = np.mean([match_asc["LON"].values[0], match_desc["LON"].values[0]])
        midpoints.append((mid_lat, mid_lon))
        synthetic_ids.append(f"syn_{pair_index}")

        new_ascending = pd.concat([new_ascending, match_asc])
        new_descending = pd.concat([new_descending, match_desc])

    result_ascending = new_ascending.reset_index(drop=True)
    result_descending = new_descending.reset_index(drop=True)

    df_synthetic = pd.DataFrame(midpoints, columns=["LAT", "LON"])
    df_synthetic["Synthetic_ID"] = synthetic_ids

    # The date columns appear twice in the spreadsheets: as datetime-typed
    # headers (used for the time axis) and as string headers (carrying the
    # velocity values). ``date_index`` is taken from the former, the value
    # columns from the latter.
    date_index = pd.to_datetime([c for c in result_descending.columns if isinstance(c, datetime)])
    value_columns = [c for c in df_ascending.columns if pd.to_datetime(c, errors="ignore") != c]
    value_column_indices = [df_ascending.columns.get_loc(c) for c in value_columns]
    time_diffs_years = date_index.to_series().diff().dt.days / 365.25
    cumulative_time_years = time_diffs_years.cumsum()

    v_hor_list = []
    v_vert_list = []
    for i in range(len(result_descending)):
        v_d_row = result_descending.iloc[i, value_column_indices].values / cumulative_time_years
        v_a_row = result_ascending.iloc[i, value_column_indices].values / cumulative_time_years

        # Horizontal component.
        v_xtop = (v_d_row / _hlosd(INCIDENCE_D)) - (v_a_row / _hlosa(INCIDENCE_A))
        v_xbot = (_elosd(INCIDENCE_D, AZIMUTH_D) / _hlosd(INCIDENCE_D)) - (
            _elosa(INCIDENCE_A, AZIMUTH_A) / _hlosa(INCIDENCE_A)
        )
        v_hor_list.append(v_xtop / v_xbot)

        # Vertical component.
        v_ytop = (v_d_row / _elosd(INCIDENCE_D, AZIMUTH_D)) - (v_a_row / _elosa(INCIDENCE_A, AZIMUTH_A))
        v_ybot = (_hlosd(INCIDENCE_D) / _elosd(INCIDENCE_D, AZIMUTH_D)) - (
            _hlosa(INCIDENCE_A) / _elosa(INCIDENCE_A, AZIMUTH_A)
        )
        v_vert_list.append(v_ytop / v_ybot)

    date_datetime = pd.to_datetime(value_columns)
    return date_datetime, v_hor_list, v_vert_list, df_synthetic


def average_components(v_hor_list, v_vert_list):
    """Reduce per-point velocity time series to the cross-point mean and
    standard deviation per acquisition. The vertical mean is sign-flipped so
    that positive values denote uplift in the figures."""
    v_hor = np.array(v_hor_list).astype(np.float64)
    v_vert = np.array(v_vert_list).astype(np.float64)
    average_v_hor = np.mean(v_hor, axis=0)
    variance_x = np.std(v_hor, axis=0)
    average_v_vert = -np.mean(v_vert, axis=0)
    variance_y = np.std(v_vert, axis=0)
    return average_v_hor, average_v_vert, variance_x, variance_y


# ---------------------------------------------------------------------------
# Figure 1 -- velocity time series
# ---------------------------------------------------------------------------

def create_velocity_time_series_figure(
    date_datetime_before, average_V_hor_before, average_V_vert_before,
    variancex_per_increment_before, variancey_per_increment_before,
    date_datetime_after, average_V_hor_after, average_V_vert_after,
    variancex_per_increment_after, variancey_per_increment_after,
    rotation_index=ROTATION_INDEX, slice_name='', save_prefix=SAVE_PREFIX
):
    """Figure 1: horizontal and vertical velocity time series, before and
    after the earthquake, with a hyperbolic fit to the post-seismic transient.
    """
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']

    fig_width_mm = 174
    fig_height_mm = 120

    fig_width_in = fig_width_mm / 25.4
    fig_height_in = fig_height_mm / 25.4

    LABEL_SIZE = 9
    TICK_SIZE = 8

    plt.rc('font', family='sans-serif')
    plt.rc('axes', labelsize=LABEL_SIZE)
    plt.rc('xtick', labelsize=TICK_SIZE)
    plt.rc('ytick', labelsize=TICK_SIZE)
    plt.rc('legend', fontsize=TICK_SIZE)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(fig_width_in, fig_height_in), dpi=300)
    plt.subplots_adjust(hspace=0.4)

    # ============= HORIZONTAL VELOCITY TIME SERIES =============

    ax1.fill_between(date_datetime_before,
                    -average_V_hor_before - variancex_per_increment_before,
                    -average_V_hor_before + variancex_per_increment_before,
                    color='blue', alpha=0.2)

    ax1.fill_between(date_datetime_after[:rotation_index],
                    -average_V_hor_after[:rotation_index] - variancex_per_increment_after[:rotation_index],
                    -average_V_hor_after[:rotation_index] + variancex_per_increment_after[:rotation_index],
                    color='orange', alpha=0.2)

    ax1.fill_between(date_datetime_after[rotation_index:],
                    -average_V_hor_after[rotation_index:] - variancex_per_increment_after[rotation_index:],
                    -average_V_hor_after[rotation_index:] + variancex_per_increment_after[rotation_index:],
                    color='red', alpha=0.2)

    ax1.plot(date_datetime_before, -average_V_hor_before, 'b-', marker='.', markersize=4,
            label='Before Earthquake', linewidth=1.5, alpha=0.8)

    ax1.plot(date_datetime_after[:rotation_index], -average_V_hor_after[:rotation_index],
            'orange', marker='.', markersize=4, linewidth=1.5, alpha=0.8,
            label='Roto-Translational Phase')

    ax1.plot(date_datetime_after[rotation_index:], -average_V_hor_after[rotation_index:],
            'r-', marker='.', markersize=4, linewidth=1.5, alpha=0.8,
            label='Constant Translational Phase')

    # Hyperbolic fit
    peak_window = 10
    peak_idx = 2 + np.argmax(np.abs(-average_V_hor_after[2:2+peak_window]))
    post_eq_dates_numeric = mdates.date2num(date_datetime_after[peak_idx:rotation_index])

    def hyperbolic_decay(x, a, b, c):
        return a/(x-b) + c

    initial_guess = [1000, mdates.date2num(date_datetime_after[0])-10, -50]

    try:
        params, _ = curve_fit(
            hyperbolic_decay,
            post_eq_dates_numeric,
            -average_V_hor_after[peak_idx:rotation_index],
            p0=initial_guess,
            maxfev=10000
        )
        x_fit = np.linspace(
            mdates.date2num(date_datetime_after[peak_idx]),
            mdates.date2num(date_datetime_after[rotation_index-1]),
            100
        )
        y_fit = hyperbolic_decay(x_fit, *params)
        ax1.plot(
            mdates.num2date(x_fit), y_fit,
            color='black', linestyle=(0, (5, 2, 1, 2)),
            linewidth=1.2, alpha=0.7,
            label='Hyperbolic Fit'
        )
    except Exception as e:
        print(f"Curve fitting failed: {str(e)}")

    # Earthquake line
    ax1.axvline(pd.Timestamp('2019-09-26'), color='black', linestyle='--',
               linewidth=1.5, alpha=0.7, label='Earthquake')

    # Configure axes
    ax1.set_ylabel('Horizontal Velocity (mm/year)', fontsize=LABEL_SIZE)
    ax1.set_ylim(-500, 500)

    major_ticks = np.arange(-500, 501, 100)
    minor_ticks = np.arange(-500, 501, 20)

    ax1.set_yticks(major_ticks)
    ax1.set_yticks(minor_ticks, minor=True)

    ax1.grid(which='major', alpha=0.4, linestyle='-', linewidth=0.5)
    ax1.grid(which='minor', alpha=0.2, linestyle=':', linewidth=0.5)

    ax1.legend(loc='upper right', fontsize=TICK_SIZE-1, framealpha=0.7,
              ncol=2, handlelength=1.5, columnspacing=1.0,
              borderpad=0.3, labelspacing=0.4)

    ax1.xaxis.set_major_locator(mdates.YearLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax1.tick_params(axis='x', labelbottom=True)
    ax1.set_xlabel('Date', fontsize=LABEL_SIZE)

    ax1.text(0.02, 0.85, 'A', transform=ax1.transAxes,
            fontsize=LABEL_SIZE+2, fontweight='bold')

    # ============= VERTICAL VELOCITY TIME SERIES =============

    ax2.fill_between(date_datetime_before,
                    average_V_vert_before - variancey_per_increment_before,
                    average_V_vert_before + variancey_per_increment_before,
                    color='blue', alpha=0.2)

    ax2.fill_between(date_datetime_after[:rotation_index],
                    average_V_vert_after[:rotation_index] - variancey_per_increment_after[:rotation_index],
                    average_V_vert_after[:rotation_index] + variancey_per_increment_after[:rotation_index],
                    color='orange', alpha=0.2)

    ax2.fill_between(date_datetime_after[rotation_index:],
                    average_V_vert_after[rotation_index:] - variancey_per_increment_after[rotation_index:],
                    average_V_vert_after[rotation_index:] + variancey_per_increment_after[rotation_index:],
                    color='red', alpha=0.2)

    ax2.plot(date_datetime_before, average_V_vert_before, 'b-', marker='.', markersize=4,
            linewidth=1.5, alpha=0.8)

    ax2.plot(date_datetime_after[:rotation_index], average_V_vert_after[:rotation_index],
            'orange', marker='.', markersize=4, linewidth=1.5, alpha=0.8)

    ax2.plot(date_datetime_after[rotation_index:], average_V_vert_after[rotation_index:],
            'r-', marker='.', markersize=4, linewidth=1.5, alpha=0.8)

    ax2.axvline(pd.Timestamp('2019-09-26'), color='black', linestyle='--',
               linewidth=1.5, alpha=0.7)

    # Configure axes
    ax2.set_xlabel('Date', fontsize=LABEL_SIZE)
    ax2.set_ylabel('Vertical Velocity (mm/year)', fontsize=LABEL_SIZE)
    ax2.set_ylim(-100, 100)

    v_major_ticks = np.arange(-100, 101, 50)
    v_minor_ticks = np.arange(-100, 101, 10)

    ax2.set_yticks(v_major_ticks)
    ax2.set_yticks(v_minor_ticks, minor=True)

    ax2.grid(which='major', alpha=0.4, linestyle='-', linewidth=0.5)
    ax2.grid(which='minor', alpha=0.2, linestyle=':', linewidth=0.5)

    ax2.xaxis.set_major_locator(mdates.YearLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

    ax2.text(0.02, 0.85, 'B', transform=ax2.transAxes,
            fontsize=LABEL_SIZE+2, fontweight='bold')

    plt.tight_layout()
    plt.savefig(f'{save_prefix}_{slice_name}_figure1.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(f'{save_prefix}_{slice_name}_figure1.png', bbox_inches='tight', dpi=300)
    plt.close(fig)

    return fig


# ---------------------------------------------------------------------------
# Figure 2 -- phase-space trajectories
# ---------------------------------------------------------------------------

def create_phase_plots_figure(
    date_datetime_before, average_V_hor_before, average_V_vert_before,
    variancex_per_increment_before, variancey_per_increment_before,
    date_datetime_after, average_V_hor_after, average_V_vert_after,
    variancex_per_increment_after, variancey_per_increment_after,
    rotation_index=ROTATION_INDEX, slice_name='', save_prefix=SAVE_PREFIX
):
    """Figure 2: horizontal-vs-vertical velocity phase-space trajectories for
    the three kinematic phases (before, roto-translational, steady)."""
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']

    fig_width_mm = 174
    fig_height_mm = 100

    fig_width_in = fig_width_mm / 25.4
    fig_height_in = fig_height_mm / 25.4

    LABEL_SIZE = 9
    TICK_SIZE = 8

    plt.rc('font', family='sans-serif')
    plt.rc('axes', labelsize=LABEL_SIZE)
    plt.rc('xtick', labelsize=TICK_SIZE)
    plt.rc('ytick', labelsize=TICK_SIZE)
    plt.rc('legend', fontsize=TICK_SIZE)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(fig_width_in, fig_height_in), dpi=300)
    plt.subplots_adjust(wspace=0.3)

    def setup_phase_plot(ax, panel_label, show_y_axis_label=True):
        # Circular polar grid only — no square grid
        for r in [50, 100, 150, 200, 250, 300]:
            circle = plt.Circle((0, 0), r, fill=False, color='lightgray', linestyle='--', alpha=0.3)
            ax.add_artist(circle)

        # Cardinal direction lines
        for angle in [0, np.pi/2, np.pi, 3*np.pi/2]:
            ax.plot([0, 350*np.cos(angle)], [0, 350*np.sin(angle)], 'lightgray', linestyle='--', alpha=0.3)

        ax.scatter(0, 0, c='black', s=40, marker='+')

        ax.text(330, 0, 'E', fontsize=TICK_SIZE, ha='center', va='center')
        ax.text(-330, 0, 'W', fontsize=TICK_SIZE, ha='center', va='center')
        ax.text(0, 310, 'Up', fontsize=TICK_SIZE, ha='center', va='center')
        ax.text(0, -330, 'Down', fontsize=TICK_SIZE, ha='center', va='center')

        ax.set_aspect('equal')
        ax.set_xlim(-350, 350)
        ax.set_ylim(-350, 350)

        major_ticks = np.arange(-300, 301, 100)
        minor_ticks = np.arange(-300, 301, 50)

        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)

        ax.text(0.02, 0.92, panel_label, transform=ax.transAxes,
               fontsize=LABEL_SIZE+2, fontweight='bold')

        if show_y_axis_label:
            ax.set_ylabel('Vertical Velocity (mm/year)', fontsize=LABEL_SIZE)

        ax.set_xlabel('Horizontal Velocity (mm/year)', fontsize=LABEL_SIZE)

    # ============= BEFORE EARTHQUAKE (A) =============
    setup_phase_plot(ax1, 'A', show_y_axis_label=True)

    ax1.plot(-average_V_hor_before, average_V_vert_before, 'b-', alpha=0.4, linewidth=1.0)
    ax1.scatter(-average_V_hor_before, average_V_vert_before, c='blue', s=12, alpha=0.7)

    key_indices = [0, len(average_V_hor_before)//2, len(average_V_hor_before)-2]
    for i in key_indices:
        if i < len(average_V_hor_before)-1:
            x1 = -average_V_hor_before[i]
            y1 = average_V_vert_before[i]
            x2 = -average_V_hor_before[i+1]
            y2 = average_V_vert_before[i+1]
            dx, dy = x2 - x1, y2 - y1
            magnitude = np.sqrt(dx**2 + dy**2)
            if magnitude > 0:
                dx = dx / magnitude * 10
                dy = dy / magnitude * 10
            ax1.arrow(x1, y1, dx, dy,
                     head_width=3, head_length=3, fc='blue', ec='blue',
                     alpha=0.8, linewidth=0.5, length_includes_head=True)

    ax1.set_title('Before Earthquake', fontsize=LABEL_SIZE, pad=10)

    # ============= ROTO-TRANSLATIONAL PHASE (B) =============
    setup_phase_plot(ax2, 'B', show_y_axis_label=False)

    ax2.plot(-average_V_hor_after[1:11], average_V_vert_after[1:11],
             'darkred', alpha=0.5, linewidth=1.2)
    ax2.plot(-average_V_hor_after[10:rotation_index], average_V_vert_after[10:rotation_index],
             'orange', alpha=0.5, linewidth=1.2)

    ax2.scatter(-average_V_hor_after[1:11], average_V_vert_after[1:11],
               c='darkred', marker='^', s=15, alpha=0.9, label='Rotation')
    ax2.scatter(-average_V_hor_after[11:rotation_index], average_V_vert_after[11:rotation_index],
               c='orange', s=1, alpha=0.7, label='Translation')

    ax2.text(-average_V_hor_after[1]+50, average_V_vert_after[1]+25,
            'Start', color='darkred', fontsize=TICK_SIZE-1, ha='center')

    for i in [1, 3, 5, 7, 9]:
        if i < 10:
            x1 = -average_V_hor_after[i]
            y1 = average_V_vert_after[i]
            x2 = -average_V_hor_after[i+1]
            y2 = average_V_vert_after[i+1]
            dx, dy = x2 - x1, y2 - y1
            magnitude = np.sqrt(dx**2 + dy**2)
            if magnitude > 0:
                dx = dx / magnitude * 25
                dy = dy / magnitude * 25
            ax2.arrow(x1, y1, dx, dy,
                     head_width=5, head_length=5, fc='darkred', ec='darkred',
                     alpha=0.9, linewidth=0.7, length_includes_head=True)

    arrow_indices = [15, 30, 45, 60, 90, 120]
    for i in arrow_indices:
        if i < rotation_index-1:
            x1 = -average_V_hor_after[i]
            y1 = average_V_vert_after[i]
            x2 = -average_V_hor_after[i+5]
            y2 = average_V_vert_after[i+5]
            dx, dy = x2 - x1, y2 - y1
            magnitude = np.sqrt(dx**2 + dy**2)
            if magnitude > 0:
                dx = dx / magnitude * 15
                dy = dy / magnitude * 15
            ax2.arrow(x1, y1, dx, dy,
                     head_width=3, head_length=3, fc='orange', ec='orange',
                     alpha=0.9, linewidth=0.5, length_includes_head=True)

    ax2.legend(loc='upper right', fontsize=TICK_SIZE, framealpha=0.7,
              borderpad=0.3, handlelength=1.0, bbox_to_anchor=(1.035, 1.025))

    ax2.set_title('Roto-Translational Phase', fontsize=LABEL_SIZE, pad=10)

    # ============= CONSTANT TRANSLATIONAL PHASE (C) =============
    setup_phase_plot(ax3, 'C', show_y_axis_label=False)

    ax3.plot(-average_V_hor_after[rotation_index:], average_V_vert_after[rotation_index:],
             'r-', alpha=0.4, linewidth=1.0)
    ax3.scatter(-average_V_hor_after[rotation_index:], average_V_vert_after[rotation_index:],
               c='red', s=12, alpha=0.7)

    if len(average_V_hor_after) > rotation_index:
        i = rotation_index + (len(average_V_hor_after) - rotation_index) // 2
        if i < len(average_V_hor_after)-5:
            x1 = -average_V_hor_after[i]
            y1 = average_V_vert_after[i]
            x2 = -average_V_hor_after[i+5]
            y2 = average_V_vert_after[i+5]
            dx, dy = x2 - x1, y2 - y1
            magnitude = np.sqrt(dx**2 + dy**2)
            if magnitude > 0:
                dx = dx / magnitude * 10
                dy = dy / magnitude * 10
            ax3.arrow(x1, y1, dx, dy,
                     head_width=3, head_length=3, fc='red', ec='red',
                     alpha=0.8, linewidth=0.5, length_includes_head=True)

    ax3.set_title('Constant Translational Phase', fontsize=LABEL_SIZE, pad=10)

    plt.tight_layout()
    plt.savefig(f'{save_prefix}_{slice_name}_figure2.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(f'{save_prefix}_{slice_name}_figure2.png', bbox_inches='tight', dpi=300)
    plt.close(fig)

    return fig


# ---------------------------------------------------------------------------
# Figure S3 -- acquisition timeline
# ---------------------------------------------------------------------------

def plot_acquisition_dates(dataframes_list, TS):
    """Figure S3: Sentinel-1 acquisition timeline for the ascending and
    descending, pre- and post-earthquake time series.

    Parameters
    ----------
    dataframes_list : list of DataFrame
        Frames containing the acquisition-date columns.
    TS : list of str
        Per-frame labels; must contain 'Ascending'/'Descending' and
        'Pre'/'Post' to drive marker placement and styling.
    """
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']

    fig_width_mm = 174
    fig_height_mm = 50

    fig_width_in = fig_width_mm / 25.4
    fig_height_in = fig_height_mm / 25.4

    LABEL_SIZE = 10
    TICK_SIZE = 9

    fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in), dpi=300)

    color_pre = 'blue'  # Matches Figure 1 and 2
    color_post = '#FF8C00'  # darkorange hex — matches Figure 1's orange line rendering

    y_value_asc = 1   # Ascending on top
    y_value_desc = 0  # Descending on bottom

    def identify_date_columns(df):
        return [col for col in df.columns if pd.to_datetime(col, errors='coerce') is not pd.NaT]

    for df, ts_label in zip(dataframes_list, TS):
        date_columns = identify_date_columns(df)
        date_datetimes = pd.to_datetime(date_columns)

        y_values = np.full(len(date_datetimes), y_value_asc if 'Ascending' in ts_label else y_value_desc)
        color = color_pre if 'Pre' in ts_label else color_post
        is_pre = 'Pre' in ts_label

        if is_pre:
            # Filled circle for before earthquake
            ax.scatter(date_datetimes, y_values,
                      marker='o', color=color,
                      s=7, alpha=0.8, linewidths=0)
        else:
            # Empty square for after earthquake
            ax.scatter(date_datetimes, y_values,
                      marker='s', facecolors='none', edgecolors=color,
                      s=7, alpha=0.8, linewidths=0.9)

    # Configure x-axis
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.tick_params(axis='x', labelsize=TICK_SIZE)

    # Configure y-axis — tighter limits to reduce whitespace
    ax.set_yticks([y_value_desc, y_value_asc])
    ax.set_yticklabels(['Descending', 'Ascending'], fontsize=LABEL_SIZE)
    ax.set_ylim(-0.4, 1.4)

    # Earthquake line
    ax.axvline(pd.Timestamp('2019-09-26'), color='black',
               linestyle='--', linewidth=1.5, alpha=0.7)

    # Legend — filled circle / empty square to match markers
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue',
               markeredgecolor='blue', markersize=5,
               label='Before earthquake', linestyle='None'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='none',
               markeredgecolor=color_post, markersize=5,
               label='After earthquake', linestyle='None'),
        Line2D([0], [0], color='black', linestyle='--', linewidth=1.5,
               alpha=0.7, label='Earthquake')
    ]

    # Clean spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='y', length=0)

    ax.legend(handles=legend_elements, loc='upper center',
              bbox_to_anchor=(0.5, -0.3),
              ncol=3, fontsize=TICK_SIZE + 1,
              framealpha=0.9, handlelength=1.5)

    fig.tight_layout()
    plt.subplots_adjust(bottom=0.3)

    plt.savefig('Acquisition_dates.pdf', dpi=300, bbox_inches='tight')
    plt.close(fig)

    return fig, ax


# ---------------------------------------------------------------------------
# Figure S4 -- synthetic PS locations
# ---------------------------------------------------------------------------

def plot_synthetic_ps_location(combined_gdf, boundary_gdf, shape_file,
                               title="Selected area with location of 'Synthetic' PSs used for decomposition",
                               filename="analysis_location_map"):
    """Figure S4: publication-quality map of the synthetic PS locations within
    the landslide boundary, over an Esri World Imagery basemap (requires
    internet access)."""
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']

    fig_width_mm = 87
    fig_height_mm = 87

    fig_width_in = fig_width_mm / 25.4
    fig_height_in = fig_height_mm / 25.4

    LABEL_SIZE = 11
    LEGEND_SIZE = 10

    fig = plt.figure(figsize=(fig_width_in, fig_height_in), dpi=300, constrained_layout=False)

    # [left, bottom, width, height] - leave room for attribution
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])

    # Ensure proper CRS
    if combined_gdf.crs is None:
        combined_gdf.set_crs(epsg=4326, inplace=True)
    if boundary_gdf.crs is None:
        boundary_gdf.set_crs(epsg=4326, inplace=True)

    # Convert to web mercator for contextily
    combined_gdf = combined_gdf.to_crs(epsg=3857)
    boundary_gdf = boundary_gdf.to_crs(epsg=3857)
    shape_file = shape_file.to_crs(epsg=3857)

    # Plot boundary and analysis area
    boundary_gdf.boundary.plot(ax=ax, color="black", linewidth=1.8)
    shape_file.plot(ax=ax, color='skyblue', alpha=0.3, edgecolor='white', linewidth=1.2)

    # Plot synthetic PS points
    ax.scatter(combined_gdf.geometry.x, combined_gdf.geometry.y, color='yellow', s=15,
              edgecolor='black', linewidth=0.3, zorder=3, alpha=0.9)

    # Add basemap (downloads tiles — requires internet)
    ctx.add_basemap(ax, source=ctx.providers.Esri.WorldImagery, attribution_size=7)

    # Fix attribution position by moving it inside the plot
    attribution = ax.texts[-1]
    attribution.set_position((0.01, 0.01))
    attribution.set_ha('left')
    attribution.set_fontsize(5)

    # Set extent based on boundary
    total_bounds = boundary_gdf.total_bounds
    ax.set_xlim([total_bounds[0], total_bounds[2]])
    ax.set_ylim([total_bounds[1], total_bounds[3]])

    # Legend
    shape_legend = mpatches.Patch(color='skyblue', alpha=0.4, label='Analysis Area')
    boundary_handle = mlines.Line2D([], [], color='black', linewidth=1.5, label='Landslide Boundary')
    ps_handle = mlines.Line2D([], [], color='black', marker='o', markerfacecolor='yellow',
                            markersize=6, linewidth=0, label='"Synthetic" PSs')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    ax.legend(handles=[shape_legend, boundary_handle, ps_handle],
             loc='lower left', fontsize=LEGEND_SIZE, framealpha=0.8,
             handlelength=1.2, borderpad=0.6)

    # North arrow
    arrow_x = 0.96
    arrow_y = 0.86
    arrow_length = 0.06
    ax.add_patch(plt.Arrow(arrow_x, arrow_y, 0, arrow_length,
                          transform=ax.transAxes, width=0.03,
                          color='black'))
    ax.text(arrow_x, arrow_y + arrow_length + 0.01, 'N', transform=ax.transAxes,
           ha='center', fontsize=LABEL_SIZE, color='black')

    plt.subplots_adjust(bottom=0.05)

    plt.savefig(f"{filename}.pdf", bbox_inches='tight', dpi=300)
    plt.savefig(f"{filename}.png", bbox_inches='tight', dpi=300)
    plt.close(fig)

    return fig, ax


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def main():
    data = load_data()

    # --- Velocity decomposition (shared by Figs 1, 2 and S4) ---
    date_before, v_hor_before, v_vert_before, synthetic_before = decompose(
        data["before_ascending"], data["before_descending"]
    )
    date_after, v_hor_after, v_vert_after, _ = decompose(
        data["after_ascending"], data["after_descending"]
    )

    avg_hor_before, avg_vert_before, var_x_before, var_y_before = average_components(
        v_hor_before, v_vert_before
    )
    avg_hor_after, avg_vert_after, var_x_after, var_y_after = average_components(
        v_hor_after, v_vert_after
    )

    # --- Figure 1: velocity time series ---
    create_velocity_time_series_figure(
        date_before, avg_hor_before, avg_vert_before, var_x_before, var_y_before,
        date_after, avg_hor_after, avg_vert_after, var_x_after, var_y_after,
        rotation_index=ROTATION_INDEX, slice_name='', save_prefix=SAVE_PREFIX,
    )
    print("Wrote landslide_analysis__figure1.pdf / .png")

    # --- Figure 2: phase-space trajectories ---
    create_phase_plots_figure(
        date_before, avg_hor_before, avg_vert_before, var_x_before, var_y_before,
        date_after, avg_hor_after, avg_vert_after, var_x_after, var_y_after,
        rotation_index=ROTATION_INDEX, slice_name='', save_prefix=SAVE_PREFIX,
    )
    print("Wrote landslide_analysis__figure2.pdf / .png")

    # --- Figure S3: acquisition timeline ---
    acquisition_frames = [
        data["raw_before_ascending"],
        data["raw_before_descending"],
        data["raw_after_descending"],
        data["raw_after_ascending"],
    ]
    acquisition_labels = [
        'TS Pre- Ascending',
        'TS Pre- Descending',
        'TS Post- Descending',
        'TS Post- Ascending',
    ]
    plot_acquisition_dates(acquisition_frames, acquisition_labels)
    print("Wrote Acquisition_dates.pdf")

    # --- Figure S4: synthetic PS locations (before-earthquake points) ---
    synthetic_before = synthetic_before.copy()
    synthetic_before['Hor_Vel_Before'] = [np.median(v) for v in v_hor_before]
    synthetic_before['Vert_Vel_Before'] = [np.median(v) for v in v_vert_before]
    combined_before_gdf = gpd.GeoDataFrame(
        synthetic_before,
        geometry=gpd.points_from_xy(synthetic_before.LON, synthetic_before.LAT),
    )
    plot_synthetic_ps_location(
        combined_before_gdf, data["landslide_shape"], data["toe_shape"],
        filename="synthetic_ps_locations",
    )
    print("Wrote synthetic_ps_locations.pdf / .png")


if __name__ == '__main__':
    main()
