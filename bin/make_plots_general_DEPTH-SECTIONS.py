#!/usr/bin/env python

import sys
import shapefile
import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib.pyplot as plt
from tomo2plt import core as PYTMCR
from tomo2plt import io as PYTMIO

from obspy.geodetics.base import kilometers2degrees, degrees2kilometers

# from PROFILES_ALL import PROFILES_EK, PROFILES_TD, PROFILES_PAPER

if len(sys.argv) != 3:
    print("USAGE: %s  GRID_CSV  CONFIG_PATH" %
          Path(sys.argv[0]).name)
    sys.exit()

# ============================================================
# ==============================================  PREPARE FONTS

import matplotlib.font_manager as fm
fname_normal = "/System/Library/Fonts/Supplemental/Arial Narrow.ttf"
fname_italic = "/System/Library/Fonts/Supplemental/Arial Narrow Italic.ttf"
fname_bold = "/System/Library/Fonts/Supplemental/Arial Narrow Bold.ttf"
fm.fontManager.addfont(fname_normal)
fm.fontManager.addfont(fname_italic)
fm.fontManager.addfont(fname_bold)
custom_font = fm.FontProperties(fname=fname_normal)
custom_font_italic = fm.FontProperties(fname=fname_italic)
custom_font_bold = fm.FontProperties(fname=fname_bold)


# ============================================================
# ==============================================  METHODS

def get_parts(shapeObj):
    points = shapeObj.points
    parts = shapeObj.parts
    parts.append(len(points))  # ensure the last part is included
    return [points[parts[i]:parts[i + 1]] for i in range(len(parts) - 1)]


# ============================================================
# ==============================================  0. PREPARE

simulout_path = Path(sys.argv[1])   # After `extract_simulps_output`
assert simulout_path.exists()
config_path = Path(sys.argv[2])
assert config_path.exists()

configs = PYTMIO.read_configuration_file(config_path, check_version=True)
workpath = Path(configs["work_path"])

if not Path(configs["DEPTH_SECTIONS"]["store_dir"]).is_dir():
    Path(configs["DEPTH_SECTIONS"]["store_dir"]).mkdir()

# ============================================================
# ==============================================  1. READING

#  Init
grid = PYTMCR.TomoGrid(sys.argv[1])
min1d = np.array(configs["INPUT_MODEL_1D"])

if isinstance(configs["DEPTH_SECTIONS"]["profiles"], (tuple, list)):
    PROFILES = configs["DEPTH_SECTIONS"]["profiles"]
elif isinstance(configs["DEPTH_SECTIONS"]["profiles"], str):
    raise ValueError("String option for profile specification not yet implemented!")
else:
    raise ValueError("Must be either a string to a *py file oa list of list")

# ==============================================================
# ======================================  HELPER-MAP-SECTIONS


if configs["DEPTH_SECTIONS"]["helper_section_map"]:

    print("... Creating Helper MAP")
    start = (configs["DEPTH_SECTIONS"]["helper_section_map_bounds"][0][0],
             configs["DEPTH_SECTIONS"]["helper_section_map_bounds"][1][0])

    end = (configs["DEPTH_SECTIONS"]["helper_section_map_bounds"][0][1],
           configs["DEPTH_SECTIONS"]["helper_section_map_bounds"][1][1])
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1)

    # -----------------------------------------  Boundaries + Coastlines

    # 1. boundaries
    if configs["DEPTH_SLICES"]["plot_boundaries"]:
        shapefile_path = configs["DATASETS"]["shapefile_boundaries"]
        sf = shapefile.Reader(shapefile_path)
        shapes = sf.shapes()
        records = sf.records()
        for xx, (shape, record) in enumerate(zip(shapes, records)):
            if record[0] == 'International boundary (verify)' and record[10] == 'Land':
                # print(xx+1, record[0], record[10], shape.parts, len(shape.points))
                parts = get_parts(shape)
                for part in parts:
                    _x, _y = zip(*part)
                    ax.plot(_x, _y, lw=1.0, color="teal")

    # 2. coastline
    if configs["DEPTH_SLICES"]["plot_coastline"]:
        shapefile_path = configs["DATASETS"]["shapefile_coastlines"]
        sf = shapefile.Reader(shapefile_path)
        shapes = sf.shapes()
        for shape in shapes:
            _x, _y = zip(*shape.points)
            ax.plot(_x, _y, lw=1.0, color="teal")

    # # ---------------------- Events
    # if configs["DEPTH_SLICES"]["plot_events"]:
    #     events = np.genfromtxt(configs["tag"]+"_events.csv")
    #     events_xyz = events[:, [2, 3, 4]]
    #     #
    #     # mask = (events_xyz[:, -1] >= _dep_min) & (events_xyz[:, -1] < _dep_max)
    #     # events_xyz = events_xyz[mask]
    #     #
    #     ax.scatter(events_xyz[:, 0], events_xyz[:, 1], marker='o',
    #                linewidths=0.7, s=6, alpha=0.9,
    #                facecolor=(0/255, 250/255, 150/225), edgecolor="black")

    # # ---------------------- AlpArray BOUND
    # bound = np.genfromtxt(configs["DATASETS"]["alparray_bound"])
    # ax.plot(bound[:, 0], bound[:, 1], color='orange', lw=1.3)

    # ---------------------- SCALE
    scaleme = np.array([[2.5, 50.5],
                        [2.5+kilometers2degrees(100), 50.5]])
    ax.plot(scaleme[:, 0], scaleme[:, 1], color='black', lw=2)
    ax.text(3.25, 50.0, '100 km', ha='center', va='center',
            fontproperties=custom_font_italic, fontsize=10)

    # ---------------------- PROFILES
    for (_prof_name, _prof_start, _prof_end) in PROFILES:
        _coords = np.array([_prof_start, _prof_end])
        ax.plot(_coords[:, 0], _coords[:, 1], color="black", lw=2.1)

    # ----------------------Set AXES
    ax.set_xlim(start[0], end[0])
    ax.set_ylim(start[1], end[1])
    ax.set_aspect(1.5)

    ax.set_xlabel("longitude (dec.deg)", fontproperties=custom_font_bold, fontsize=12)
    ax.set_xticklabels(ax.get_xticklabels(), fontproperties=custom_font_italic, fontsize=11)

    ax.set_ylabel("latitude (dec.deg)", fontproperties=custom_font_bold, fontsize=12, fontweight="bold")
    ax.set_yticklabels(ax.get_yticklabels(), fontproperties=custom_font_italic, fontsize=11, fontweight="bold")

    # --------------- SAVE
    plt.savefig(("%s/HelperMap_profiles.png" % (
                    configs["DEPTH_SECTIONS"]["store_dir"])),
                format='png', bbox_inches='tight', dpi=310)
    plt.savefig(("%s/HelperMap_profiles.pdf" % (
                    configs["DEPTH_SECTIONS"]["store_dir"])),
                format='pdf', bbox_inches='tight', dpi=310)


# ========================================================
# ======================================  SECTIONS

for (profile_name, start, end) in PROFILES:

    print("SECTION:  %s --> Start: %5.3f / %5.3f  End: %5.3f / %5.3f" % (
           profile_name, *start[:2], *end[:2]))
    (fig, ax1, cbar) = grid.plot_depth_section(
                        start, end,
                        increment_x_km=configs["DEPTH_SECTIONS"]["x_km_sect"],
                        increment_y_km=configs["DEPTH_SECTIONS"]["y_km_sect"],
                        what=configs["DEPTH_SECTIONS"]["what_to_plot"],
                        smooth=configs["DEPTH_SECTIONS"]["gauss_smooth"],
                        mask_rde=configs["DEPTH_SECTIONS"]["mask_rde"],
                        interpolate="linear",
                        add_topography=(configs["DATASETS"]["dem_grid"]
                                        if configs["DEPTH_SECTIONS"]["plot_dem"]
                                        else False),
                        palettes=configs["PALETTES"]["absolute"],
                        isolines=configs["DEPTH_SECTIONS"]["isolines"])

    # ----------------------  EVENTS
    if configs["DEPTH_SLICES"]["plot_events"]:
        print("Plotting events ...")

        # ---- Read
        EVENTS_PATH = (Path(configs["work_path"]) / (configs["tag"] +
                                                     "_events.csv"))
        EVENTS_PD = pd.read_csv(str(EVENTS_PATH),
                                names=["OT_DATE", "OT_TIME",
                                       "LON", "LAT", "DEP",
                                       "MAG", "NOBS", "RMS", "X", "Y", "Z"],
                                sep="\s+")
        EVENTS_ARR = EVENTS_PD[["LON", "LAT", "DEP"]].values
        EVENTS = PYTMCR.Points_2D_Grid(EVENTS_ARR, tag="EVENTS")

        # ---- Project / plot
        (_prof, _depth) = EVENTS.project_points(
                        start, end,
                        project=configs["DEPTH_SECTIONS"]["project_events"])
        ax1.scatter(_prof, _depth, marker='o',
                    linewidths=0.7, s=6, alpha=0.9,
                    facecolor=(0/255, 250/255, 150/225), edgecolor="black")

    # ----------------------  INV.NODES depth
    for (_x, _y) in zip([80.0 for dd in range(min1d.shape[0])], min1d[:, 0]):
        if _y >= 0 and _y <= end[-1]:
            ax1.scatter(_x, _y, marker='o',
                        color="darkgray", edgecolors="black", s=30)

    # # ----------------------  Spada2013 MOHOs
    MOHO_EUR = PYTMCR.Plane_2D_Grid(
                configs["DATASETS"]["moho_Europe"], tag="European_Moho")
    MOHO_ADR = PYTMCR.Plane_2D_Grid(
                configs["DATASETS"]["moho_Adria"], tag="Adria_Moho")
    MOHO_TYR = PYTMCR.Plane_2D_Grid(
                configs["DATASETS"]["moho_Tirrenia"], tag="Tirrenia_Moho")
    MOHO_KOS = PYTMCR.Plane_2D_Grid(
                configs["DATASETS"]["moho_new"], tag="Kostantinos2023")

    for _xx, (_moho, _color, _interp) in enumerate(
                            #      0          1         2         3
                            zip([MOHO_EUR, MOHO_ADR, MOHO_TYR, MOHO_KOS],
                                ["white", "black", "red", "purple"],
                                [2.0, 2.0, 2.0, 20.0])):
        print("Plotting:  %s" % _moho.tag)
        try:
            (_prof, _depth) = _moho.get_profile(start, end, _interp,
                                                query_profile_only=True)
            ax1.plot(_prof, _depth, color=_color)
        except ValueError:
            print("    nothing to interp...skipping!")

    # --------------- FONT CHANGE
    # Title
    ax1.set_title(("Start: %.3f / %.3f - End: %.3f / %.3f" % (*start[:2], *end[:2])),
                  fontproperties=custom_font_bold, fontsize=14)

    # X-axis
    ax1.set_xlabel("distance along profile (km)", fontproperties=custom_font_bold, fontsize=12)
    ax1.set_xticklabels(ax1.get_xticklabels(), fontproperties=custom_font_italic, fontsize=11)

    # Y-axis
    ax1.set_ylabel("depth (km)", fontproperties=custom_font_bold, fontsize=12, fontweight="bold")
    ax1.set_yticklabels(ax1.get_yticklabels(), fontproperties=custom_font_italic, fontsize=11, fontweight="bold")

    # Colorbar
    cbar.ax.set_ylabel(
               "Vp (%)",
               fontproperties=custom_font_bold, fontsize=12, fontweight="bold")
    cbar.ax.set_yticklabels(
            cbar.ax.get_yticklabels(),
            fontproperties=custom_font_italic, fontsize=11, fontweight="bold")

    # --------------- SAVE
    plt.savefig(("%s/depth_SECTION__Start_%5.3f_%5.3f__End_%5.3f_%5.3f_Name_%s.png" % (
                 configs["DEPTH_SECTIONS"]["store_dir"], *start[:2], *end[:2],
                 profile_name)), format='png', bbox_inches='tight', dpi=310)
    plt.savefig(("%s/depth_SECTION__Start_%5.3f_%5.3f__End_%5.3f_%5.3f_Name_%s.pdf" % (
                 configs["DEPTH_SECTIONS"]["store_dir"], *start[:2], *end[:2],
                 profile_name)), format='pdf', bbox_inches='tight', dpi=310)

    if configs["DEPTH_SECTIONS"]["show"]:
        plt.show()

print("DONE SECTIONING ...")
