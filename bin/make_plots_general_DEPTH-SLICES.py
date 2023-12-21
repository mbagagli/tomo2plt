#!/usr/bin/env python

import sys
import shapefile
import numpy as np
from tqdm import tqdm
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


def find_half_distances(arr, depth):
    """ Helper function to plot relevant events for layer definition """
    # Find the index of the row that matches the given depth
    index = np.where(arr[:, 0] == depth)[0][0]

    # Calculate half distance with the previous depth, if not the first element
    if index > 0:
        prev_half_distance = (arr[index, 0] + arr[index - 1, 0]) / 2
    else:
        prev_half_distance = None

    # Calculate half distance with the next depth, if not the last element
    if index < len(arr) - 1:
        next_half_distance = (arr[index, 0] + arr[index + 1, 0]) / 2
    else:
        next_half_distance = None

    return (prev_half_distance, next_half_distance)


# ============================================================
# ==============================================  0. PREPARE

simulout_path = Path(sys.argv[1])   # After `extract_simulps_output`
assert simulout_path.exists()
config_path = Path(sys.argv[2])
assert config_path.exists()

configs = PYTMIO.read_configuration_file(config_path, check_version=True)
workpath = Path(configs["work_path"])

if not Path(configs["DEPTH_SLICES"]["store_dir"]).is_dir():
    Path(configs["DEPTH_SLICES"]["store_dir"]).mkdir()

start = (configs["DEPTH_SLICES"]["lon_min"],
         configs["DEPTH_SLICES"]["lat_min"],
         configs["DEPTH_SLICES"]["dep_min"])
end = (configs["DEPTH_SLICES"]["lon_max"],
       configs["DEPTH_SLICES"]["lat_max"],
       configs["DEPTH_SLICES"]["dep_max"])

# ============================================================
# ==============================================  1. READING

#  Init
grid = PYTMCR.TomoGrid(sys.argv[1])
min1d = np.array(configs["INPUT_MODEL_1D"])

for _dep in np.unique(np.sort(grid.coordinates[:, -1])):
    mask = min1d[:, 0] == _dep
    _vel = min1d[mask]
    assert _vel.shape == (1, 2)
    _vel = _vel[0, 1]

    # -----------------------------------------  Core-Image
    (fig, ax, cbar) = grid.plot_depth_slice(
                        start, end,
                        increment_x_km=configs["DEPTH_SLICES"]["x_km_slice"],
                        increment_y_km=configs["DEPTH_SLICES"]["y_km_slice"],
                        depth=_dep,
                        what=configs["DEPTH_SLICES"]["what_to_plot"],
                        smooth=configs["DEPTH_SLICES"]["gauss_smooth"],
                        mask_rde=configs["DEPTH_SLICES"]["mask_rde"],
                        interpolate=configs["DEPTH_SLICES"]["interpolate"],
                        palettes=configs["PALETTES"]["relative"],
                        isolines=configs["DEPTH_SLICES"]["isolines"])

    # -----------------------------------------  Boundaries+Coastlines

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

    # ---------------------- Grid Nodes
    if configs["DEPTH_SLICES"]["plot_nodes"]:
        nodes = np.unique(grid.coordinates[:, :2], axis=0)
        ax.scatter(nodes[:, 0], nodes[:, 1], marker='+',
                   linewidths=0.5, facecolor='black',  # edgecolor="", color='black',
                   s=3, alpha=0.25)

    # ---------------------- AlpArray BOUND
    bound = np.genfromtxt(configs["DATASETS"]["alparray_bound"])
    ax.plot(bound[:, 0], bound[:, 1], color='orange', lw=1.3)

    # ---------------------- Events
    if configs["DEPTH_SLICES"]["plot_events"]:
        # Test the function with depth = 2.5
        (_dep_min, _dep_max) = find_half_distances(min1d, _dep)
        events = np.genfromtxt(configs["tag"]+"_events.csv")
        events_xyz = events[:, [2, 3, 4]]
        #
        mask = (events_xyz[:, -1] >= _dep_min) & (events_xyz[:, -1] < _dep_max)
        events_xyz = events_xyz[mask]
        #
        ax.scatter(events_xyz[:, 0], events_xyz[:, 1], marker='o',
                   linewidths=0.7, s=6, alpha=0.8,
                   facecolor=(0/255, 250/255, 150/225), edgecolor="black")

    # ---------------------- SCALE
    scaleme = np.array([[2.5, 50.5],
                        [2.5+kilometers2degrees(100), 50.5]])
    ax.plot(scaleme[:, 0], scaleme[:, 1], color='black', lw=2)
    ax.text(3.25, 50.0, '100 km', ha='center', va='center',
            fontproperties=custom_font_italic, fontsize=10)

    # ----------------------Set AXES
    ax.set_xlim(start[0], end[0])
    ax.set_ylim(start[1], end[1])
    ax.set_aspect(1.5)  # remember to change

    # --------------- FONT CHANGE
    # Title
    ax.set_title(('DEPTH SLICE @ %.2f km - %.2f km/s' % (_dep, _vel)),
                 fontproperties=custom_font_bold, fontsize=14)

    # X-axis
    ax.set_xlabel("longitude (dec.deg)",
                  fontproperties=custom_font_bold, fontsize=12)
    ax.set_xticklabels(ax.get_xticklabels(),
                       fontproperties=custom_font_italic, fontsize=11)

    # Y-axis
    ax.set_ylabel("latitude (dec.deg)",
                  fontproperties=custom_font_bold, fontsize=12, fontweight="bold")
    ax.set_yticklabels(ax.get_yticklabels(),
                       fontproperties=custom_font_italic, fontsize=11, fontweight="bold")

    # Colorbar
    cbar.ax.set_ylabel(
               "Vp (%)",
               fontproperties=custom_font_bold, fontsize=12, fontweight="bold")
    cbar.ax.set_yticklabels(
            cbar.ax.get_yticklabels(),
            fontproperties=custom_font_italic, fontsize=11, fontweight="bold")

    # --------------- SAVE
    plt.savefig(("%s/%s_DepthSlices_%05.1f_km.png" % (
                 configs["DEPTH_SLICES"]["store_dir"],
                 configs["tag"], _dep)),
                format='png', bbox_inches='tight', dpi=310)
    plt.savefig(("%s/%s_DepthSlices_%05.1f_km.pdf" % (
                 configs["DEPTH_SLICES"]["store_dir"],
                 configs["tag"], _dep)),
                format='pdf', bbox_inches='tight', dpi=310)

    if configs["DEPTH_SLICES"]["show"]:
        plt.show()


print("DONE SLICING ...")
print("")
