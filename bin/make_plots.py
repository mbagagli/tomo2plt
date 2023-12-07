#!/usr/bin/env python

import sys
import shapefile
from tqdm import tqdm
from pathlib import Path
from tomo2plt import core as PYTMCR
import matplotlib.pyplot as plt
import numpy as np
from PROFILES_ALL import PROFILES_EK, PROFILES_TD, PROFILES_PAPER

import matplotlib.font_manager as fm

from obspy.geodetics.base import kilometers2degrees, degrees2kilometers


fname_normal = "/System/Library/Fonts/Supplemental/Arial Narrow.ttf"
fname_italic = "/System/Library/Fonts/Supplemental/Arial Narrow Italic.ttf"
fname_bold = "/System/Library/Fonts/Supplemental/Arial Narrow Bold.ttf"
fm.fontManager.addfont(fname_normal)
fm.fontManager.addfont(fname_italic)
fm.fontManager.addfont(fname_bold)
custom_font = fm.FontProperties(fname=fname_normal)
custom_font_italic = fm.FontProperties(fname=fname_italic)
custom_font_bold = fm.FontProperties(fname=fname_bold)



# ================  FAST MIN1D
# -5.00      5.250
#  2.50      5.800
# 10.00      6.000
# 20.00      6.200
# 30.00      6.800
# 45.00      7.800
# 60.00      8.100
# 90.00      8.150
MIN1D_VELCT = [5.25, 5.80,  6.00,  6.20,  6.80,  7.80,  8.10,  8.15]
MIND1D_NODES = [-5.00, 2.50, 10.00, 20.00, 30.00, 45.00, 60.00, 90.00]

DO_SLICES = False
X_KM_SLICE = 15       # cell size in longitude (km)
Y_KM_SLICE = 15       # cell size in latitude (km)
SMOOTH_SLICE = 0.9    # GAUSSIAN KERNEL value or None/False
SLICE_PLOT_DIR = "./SLICES_%.01fx%.01f_%.01f" % (X_KM_SLICE, Y_KM_SLICE, SMOOTH_SLICE)

DO_HELPER_SECTION_MAP = True

DO_SECTIONS = True
X_KM_SECT = 9.0       # cell size in longitude (km)
Y_KM_SECT = 2.5       # cell size in latitude (km)
SMOOTH_SECT = 0.9    # GAUSSIAN KERNEL value or None/False
SECT_PLOT_DIR = "./SECTIONS_%.01fx%.01f_%.01f" % (X_KM_SECT, Y_KM_SECT, SMOOTH_SECT)
DEM_GRID = "/Users/matteo/ETH/Tomography_TD+EPOS+AARSC_FMTOMO/DEM_GRIDs_ALPS/aadem.1km.2.5-19-41.5-49.5.xyz"  # 500m / 1km

MASK_RDE = 0.2  # RDE value or None/False

#  Init
grid = PYTMCR.TomoGrid(sys.argv[1], parents=True, exist_ok=True)

# ===============================
def get_parts(shapeObj):
    points = shapeObj.points
    parts = shapeObj.parts
    parts.append(len(points))  # ensure the last part is included
    return [points[parts[i]:parts[i + 1]] for i in range(len(parts) - 1)]


# ========================================================
# ======================================  DEPTH

if DO_SLICES:
    if not Path(SLICE_PLOT_DIR).is_dir():
        Path(SLICE_PLOT_DIR).mkdir()

    start = (1, 41, 0.0)
    end = (21, 51, 70.0)

    for _dep, _vel in zip(MIND1D_NODES[1:-1], MIN1D_VELCT[1:-1]):
        (fig, ax) = grid.plot_depth_slice(
                            start, end,
                            increment_x_km=X_KM_SLICE, increment_y_km=Y_KM_SLICE,
                            depth=_dep,  what="delta",
                            smooth=SMOOTH_SLICE, mask_rde=MASK_RDE,
                            interpolate="linear")

        # ---------------------- COASTLINES + BOUNDARIES

        # 1. boundaries
        shapefile_path = "./datas/coastlines/ne_10m_admin_0_boundary_lines_land/ne_10m_admin_0_boundary_lines_land.shp"
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
        shapefile_path = "./datas/coastlines/ne_10m_coastline/ne_10m_coastline.shp"
        sf = shapefile.Reader(shapefile_path)
        shapes = sf.shapes()
        for shape in shapes:
            _x, _y = zip(*shape.points)
            ax.plot(_x, _y, lw=1.0, color="teal")

        # ---------------------- Inversion Grid
        nodes = np.genfromtxt("./datas/MOD4.grid.official")
        ax.scatter(nodes[:, 0], nodes[:, 1], marker='+',
                   linewidths=1, edgecolors='black',  # color='black',
                   s=3, alpha=0.25)

        # ---------------------- AlpArray BOUND
        bound = np.genfromtxt("./datas/AlpArray_BORDER_250km.dat")
        ax.plot(bound[:, 0], bound[:, 1], color='orange', lw=1.3)

        # ---------------------- SCALE
        scaleme = np.array([[2.5, 50.5],
                            [2.5+kilometers2degrees(100), 50.5]])
        ax.plot(scaleme[:, 0], scaleme[:, 1], color='black', lw=2)
        ax.text(3.25, 50.0, '100 km', ha='center', va='center',
                fontproperties=custom_font_italic, fontsize=10)

        # ----------------------Set AXES
        ax.set_xlim(start[0], end[0])
        ax.set_ylim(start[1], end[1])
        ax.set_aspect(1.5)

        # --------------- FONT CHANGE
        ax.set_title(('DEPTH SLICE @ %.2f km - %.2f km/s' % (_dep, _vel)),
                      fontproperties=custom_font_bold, fontsize=14)

        ax.set_xlabel("longitude (dec.deg)", fontproperties=custom_font_bold, fontsize=12)
        ax.set_xticklabels(ax.get_xticklabels(), fontproperties=custom_font_italic, fontsize=11)

        ax.set_ylabel("latitude (dec.deg)", fontproperties=custom_font_bold, fontsize=12, fontweight="bold")
        ax.set_yticklabels(ax.get_yticklabels(), fontproperties=custom_font_italic, fontsize=11, fontweight="bold")

        plt.savefig(('%s/depth_SLICE_plot_%05.2f_km_%.3f_kms.png' % (
                    SLICE_PLOT_DIR, _dep, _vel)),
                    format='png', bbox_inches='tight', dpi=310)

        plt.savefig(('%s/depth_SLICE_plot_%05.2f_km_%.3f_kms.pdf' % (
                    SLICE_PLOT_DIR, _dep, _vel)),
                    format='pdf', bbox_inches='tight', dpi=310)
        # plt.show()

    print("DONE SLICING ...")
    print("")


# ==============================================================
# ======================================  HELPER-MAP-SECTIONS

# PROFILES_TO_PLOT =  PROFILES_EK + PROFILES_TD
PROFILES_TO_PLOT = PROFILES_PAPER


if DO_HELPER_SECTION_MAP:

    if not Path(SECT_PLOT_DIR).is_dir():
        Path(SECT_PLOT_DIR).mkdir()

    print("... Creating Helper MAP")

    start = (1, 41, 0.0)
    end = (21, 51, 70.0)
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1)

    # ---------------------- COASTLINES + BOUNDARIES

    # 1. boundaries
    shapefile_path = "./datas/coastlines/ne_10m_admin_0_boundary_lines_land/ne_10m_admin_0_boundary_lines_land.shp"
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
    shapefile_path = "./datas/coastlines/ne_10m_coastline/ne_10m_coastline.shp"
    sf = shapefile.Reader(shapefile_path)
    shapes = sf.shapes()
    for shape in shapes:
        _x, _y = zip(*shape.points)
        ax.plot(_x, _y, lw=1.0, color="teal")

    # # ---------------------- Inversion Grid
    # nodes = np.genfromtxt("./datas/MOD4.grid.official")
    # ax.scatter(nodes[:, 0], nodes[:, 1], marker='+',
    #            linewidths=1, edgecolors='black',  # color='black',
    #            s=3, alpha=0.25)

    # ---------------------- AlpArray BOUND
    bound = np.genfromtxt("./datas/AlpArray_BORDER_250km.dat")
    ax.plot(bound[:, 0], bound[:, 1], color='orange', lw=1.3)

    # ---------------------- SCALE
    scaleme = np.array([[2.5, 50.5],
                        [2.5+kilometers2degrees(100), 50.5]])
    ax.plot(scaleme[:, 0], scaleme[:, 1], color='black', lw=2)
    ax.text(3.25, 50.0, '100 km', ha='center', va='center',
            fontproperties=custom_font_italic, fontsize=10)

    # ---------------------- PROFILES
    for (_prof_name, _prof_start, _prof_end) in PROFILES_PAPER:
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
    plt.savefig(("%s/HelperMap_profiles.png" % (SECT_PLOT_DIR)),
                 format='png', bbox_inches='tight', dpi=310)
    plt.savefig(("%s/HelperMap_profiles.pdf" % (SECT_PLOT_DIR)),
                 format='pdf', bbox_inches='tight', dpi=310)
    # plt.show()



# ========================================================
# ======================================  SECTIONS

if DO_SECTIONS:

    if not Path(SECT_PLOT_DIR).is_dir():
        Path(SECT_PLOT_DIR).mkdir()

    for (profile_name, start, end) in PROFILES_TO_PLOT:

        print("SECTION:  %s --> Start: %5.3f / %5.3f  End: %5.3f / %5.3f" % (
                profile_name, *start[:2], *end[:2]))
        (fig, ax1) = grid.plot_depth_section(
                            start, end,
                            increment_x_km=X_KM_SECT, increment_y_km=Y_KM_SECT,
                            what="pvel",
                            smooth=SMOOTH_SECT, mask_rde=MASK_RDE,
                            interpolate="linear",
                            add_topography=DEM_GRID)

        # ----------------------  INV.NODES depth
        ax1.scatter(np.full(len(MIND1D_NODES)-1, 80),
                    MIND1D_NODES[1:], marker='o',
                    color="darkgray",
                    edgecolors="black", s=30)

        # # ----------------------  Spada2013 MOHOs
        print("")
        MOHO_EUR = PYTMCR.Plane_2D_Grid("./datas/mohos/European.xyz",
                                        tag="European_Moho")
        MOHO_ADR = PYTMCR.Plane_2D_Grid("./datas/mohos/Adria.xyz",
                                        tag="Adria_Moho")
        MOHO_TYR = PYTMCR.Plane_2D_Grid("./datas/mohos/Tirrenia.xyz",
                                        tag="Tirrenia_Moho")
        MOHO_KOS = PYTMCR.Plane_2D_Grid("./datas/mohos/Konstantinos2023.csv",
                                        tag="Tirrenia_Kostantinos2023")

        for _moho, _color in zip([MOHO_EUR, MOHO_ADR, MOHO_TYR, MOHO_KOS],
                                 ["white", "black", "red", "purple"]):
            print("Plotting:  %s" % _moho.tag)
            try:
                (_prof, _depth) = _moho.get_profile(start, end, 2.0, query_profile_only=True)
                ax1.plot(_prof, _depth, color=_color)
            except ValueError:
                print("    nothing to interp...skipping!")

        # --------------- FONT CHANGE
        ax1.set_title(("Start: %.3f / %.3f - End: %.3f / %.3f" % (*start[:2], *end[:2])),
                      fontproperties=custom_font_bold, fontsize=14)

        ax1.set_xlabel("distance along profile (km)", fontproperties=custom_font_bold, fontsize=12)
        ax1.set_xticklabels(ax1.get_xticklabels(), fontproperties=custom_font_italic, fontsize=11)

        ax1.set_ylabel("depth (km)", fontproperties=custom_font_bold, fontsize=12, fontweight="bold")
        ax1.set_yticklabels(ax1.get_yticklabels(), fontproperties=custom_font_italic, fontsize=11, fontweight="bold")


        # --------------- SAVE
        plt.savefig(("%s/depth_SECTION__Start_%5.3f_%5.3f__End_%5.3f_%5.3f_Name_%s.png" % (
                     SECT_PLOT_DIR, *start[:2], *end[:2], profile_name)),
                     format='png', bbox_inches='tight', dpi=310)
        plt.savefig(("%s/depth_SECTION__Start_%5.3f_%5.3f__End_%5.3f_%5.3f_Name_%s.pdf" % (
                     SECT_PLOT_DIR, *start[:2], *end[:2], profile_name)),
                     format='pdf', bbox_inches='tight', dpi=310)
        # plt.show()

    print("DONE SECTIONING ...")
