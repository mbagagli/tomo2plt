#!/usr/bin/env python

import os
import sys
import yaml
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import Path
from pathlib import Path as PathFile

from tomo2plt.io import MOD
matplotlib.use('TkAgg')


try:
    # Load data from the YAML file
    with open(sys.argv[1], 'r') as yaml_file:
        configs = yaml.load(yaml_file, Loader=yaml.FullLoader)
except IndexError:
    print("USAGE:  Make3DsimulPS_SYNTHETICS.py  CONFINGPATH")
    sys.exit()


sys.stdout.write(
    os.linesep.join([
        "# ==================================================================",
        "# Make 3D MOD synthetics",
        "#",
        "# Perturbation will be `layered`, by doing so, you can override",
        "# previously selected points.",
        "#",
        "# Please beware ALL NODES will be plotted..also the boundary nodes!",
        "# Follow the instructions on the figure for",
        "#",
        "# WARNING:  THE SOFTWARE TREAT MODELS as WEST-POSITIVE!",
        "#",
        "# AUTHOR: MatteoBagagli @ 10.2022",
        "# ==================================================================",
        "",
        ""
    ])+os.linesep)


# -------------  SIMUL
# short distance conversion factors.
# NB: It is origin dependent, copy them from your SIMUL output file)
#      one min lat   1.8525 km --> 1 degree (*60 min) 111.15 km
#      one min lon   1.2911 km --> 1 degree (*60 min) 77.466 km

# ORIGINLON = 10.5
# ORIGINLAT = 46.0
# DEGKMLAT = 111.15  # km  (short-distance conversion Y --> check SIMULPS output)
# DEGKMLON = 77.466  # km  (short-distance conversion X --> check SIMULPS output)


fig = plt.figure(figsize=(8, 7), dpi=150)


# ===========================================================
# ==============================================  DEF
def _in_or_out(vertices, pointspool):
    """ Vertices and Points mus tbe in the form of [(x1, y1), (x2, y2), ...] """
    hull_path = Path(np.array(vertices))
    return [(_x, _y) for _x, _y in pointspool if hull_path.contains_point((_x,_y))]


def get_parts(shapeObj):
    points = shapeObj.points
    parts = shapeObj.parts
    parts.append(len(points))  # ensure the last part is included
    return [points[parts[i]:parts[i + 1]] for i in range(len(parts) - 1)]


def plot_coastlines_and_boundaries(ax, boundaries, coastlines):
    # ---------------------- COASTLINES + BOUNDARIES
    # 1. boundaries

    sf = shapefile.Reader(boundaries)
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
    sf = shapefile.Reader(coastlines)
    shapes = sf.shapes()
    for shape in shapes:
        _x, _y = zip(*shape.points)
        ax.plot(_x, _y, lw=1.0, color="teal")

    return ax

# ===========================================================
# ==============================================  WORK

try:
    modobj = MOD(configs["model_path"])
    modobj.set_origin(configs["origin"][0], configs["origin"][1])
    modobj.calculate_lon_lat_grid()


except IndexError:
    print("USAGE: script.py  MODPATH ")
    sys.exit()

print("### Loaded model:  %s" % sys.argv[1])
print("")
print("  IDX   Depth  Velocity")
print("-------------------------")
for xx, (dd, vv) in enumerate(zip(modobj.zarr, modobj.velocities)):
    print("   %2d  %7.2f   %4.2f" % (xx, dd, vv))

# ------------------------------------------------------
# ----------------- Select regime (WEST - POSITIVE)
geoswitch = configs["use_geographic_coordinates"]
if geoswitch:
    print("... using geographic coordinates [dec.deg.]")
    XLIM = [min(modobj.lon_grid), max(modobj.lon_grid)]
    _xlim_tolerance = np.abs(XLIM[1]-XLIM[0])*0.05
    XLIM = [XLIM[0]-_xlim_tolerance, XLIM[1]+_xlim_tolerance]

    YLIM = [min(modobj.lat_grid), max(modobj.lat_grid)]
    _ylim_tolerance = np.abs(YLIM[1]-YLIM[0])*0.05
    YLIM = [YLIM[0]-_ylim_tolerance, YLIM[1]+_ylim_tolerance]

else:
    print("... using realative coordinates [km]")
    XLIM = [min(modobj.xarr*-1), max(modobj.xarr*-1)]
    _xlim_tolerance = np.abs(XLIM[1]-XLIM[0])*0.05
    XLIM = [XLIM[0]-_xlim_tolerance, XLIM[1]+_xlim_tolerance]

    YLIM = [min(modobj.yarr), max(modobj.yarr)]
    _ylim_tolerance = np.abs(YLIM[1]-YLIM[0])*0.05
    YLIM = [YLIM[0]-_ylim_tolerance, YLIM[1]+_ylim_tolerance]

#
if geoswitch:
    try:
        import shapefile

        geoswitch, geocoord = True, modobj.get_origin()
        x_axplot = [((_x*-1)/configs["lon_short_distance_conversion"])+geocoord[0] for _x in modobj.xarr]
        y_axplot = [(_y/configs["lat_short_distance_conversion"])+geocoord[1] for _y in modobj.yarr]

    except ImportError:
        print("    ... PySHP package Missing, automatically switched to cartesian!")
        geoswitch, geocoord = False, None
        x_axplot = modobj.xarr*-1
        y_axplot = modobj.yarr

else:
    print("    ... Switch to cartesian")
    geoswitch, geocoord = False, None
    x_axplot = modobj.xarr*-1
    y_axplot = modobj.yarr


xplt, yplt = [], []
for xx in x_axplot:
    for yy in y_axplot:
        xplt.append(xx)
        yplt.append(yy)


# ======================  Loop MODIFY

MODIFICATIONS = {str(xx): [] for xx in range(len(modobj.zarr))}

while True:
    laynum = input(" --> Which layer (IDX) you want to perturbate?  ([q] to quit)  ")
    if not laynum or laynum.lower() == "q":
        # End of loop
        break
    else:
        laynum = int(laynum)

    # =========================================================
    # ==================================================  GUI
    happy = False
    while not happy:

        if geoswitch:
            ax = plt.subplot(1, 1, 1)
            ax = plot_coastlines_and_boundaries(
                ax,
                "/Users/matteo/GEOPHYSICS/imaging_tomo/datasets/shapefiles/"
                "ne_10m_admin_0_boundary_lines_land/ne_10m_admin_0_boundary_lines_land.shp",
                "/Users/matteo/GEOPHYSICS/imaging_tomo/datasets/shapefiles/"
                "ne_10m_coastline/ne_10m_coastline.shp")

            ax.set_xlim(XLIM)
            ax.set_xlabel("longitude")

            ax.set_ylim(YLIM)
            ax.set_ylabel("latitude")

        else:
            ax = fig.add_subplot(1, 1, 1)
            ax.set_xlabel("x")
            ax.set_ylabel("y")

        ax.plot(xplt, yplt, marker=".", markersize=2,
                linestyle="None", color="black")  # "teal"

        fig.suptitle('Select area. ESC to stop', fontsize=15)
        ax.set_title(("Depth:  %d km - Velocity:  %.2f km/s" %
                      (int(modobj.zarr[laynum]),  modobj.velocities[laynum])))

        # -- Before going on with selection, check if modification are already present
        if MODIFICATIONS[str(laynum)]:
            # We have already polygons --> plot them
            for _pt in MODIFICATIONS[str(laynum)]:
                _col = "blue" if int(_pt[0]) > 0 else "red"
                _pt_plot, = ax.fill([cc[0] for cc in _pt[1]],
                                    [cc[1] for cc in _pt[1]],
                                    _col, lw=2, alpha=0.2)
                # --- Return selection among point
                _pt_plot_point = _in_or_out(np.array(_pt[1]), zip(xplt, yplt))

                ax.plot([xy[0] for xy in _pt_plot_point],
                        [xy[1] for xy in _pt_plot_point],
                        marker="x", markersize=5,
                        linestyle="None", color=_col)
                # -----------------
                plt.draw()
        #
        selection_edges = plt.ginput(n=-1, timeout=60, show_clicks=True,
                                     mouse_pop=None,
                                     mouse_stop='1')
        if not selection_edges:
            # I've closed the figure --> wanna go to layer selection
            fig = plt.figure(figsize=(7, 5), dpi=150)
            break
        #
        fillme, = ax.fill([cc[0] for cc in selection_edges],
                          [cc[1] for cc in selection_edges],
                          'gray', lw=2, alpha=0.5)
        plt.draw()
        fig.suptitle('Sure of the selection? YES: keyboard, NO: mouse click', fontsize=15)

        happy = plt.waitforbuttonpress()
        print(happy)
        # if happy == False:
        if not happy:
            # Clicked mouse --> not happy with drawing..redo
            fillme.remove()  # remove the artist
            plt.draw()       # force rendering
            continue
        else:
            # Keyboard press --> selection OK , go further
            plt.close(fig=fig)
            fig = plt.figure(figsize=(8, 7), dpi=150)
            # ==================================================  DO SOMETHING
            pert_perc = input(" --> --> How much you want to perturbate (%)? (+/- [0] to 100)   ")
            MODIFICATIONS[str(laynum)].append((pert_perc, selection_edges))

            # ==================================================
            happy = False
            break


# ========= SAVE OBJECT
MODIFICATIONS_LAYER = [(kk, vv) for kk, vv in MODIFICATIONS.items()]

print("")
if not any(MODIFICATIONS_LAYER):
    print(" ... No mods to be done ... exit!")
else:
    print(" ... Applying modifications")
    for idx, (_laynum, _pert) in enumerate(MODIFICATIONS_LAYER):
        if len(_pert) > 0:
            for (_amount, _vertices) in _pert:

                _pt_plot_point = _in_or_out(_vertices, zip(xplt, yplt))

                _amount = float(_amount)
                _velpert = modobj.velocities[int(_laynum)] + (
                           modobj.velocities[int(_laynum)] * (_amount/100.0))

                for (_x, _y) in _pt_plot_point:
                    modobj.matrix[int(_laynum),
                                  y_axplot.index(_y),
                                  x_axplot.index(_x)] = _velpert

    print("Exporting model")
    if configs["west_positive"]:
        modobj.write_mod(outfile="MOD.synth",
                         reverse_lon=False,
                         reverse_lat=False)
    else:
        print("WARNING: using East-Positive mode! "
              "Compatibility issue may occur with SIMUL")
        modobj.write_mod(outfile="MOD.synth",
                         reverse_lon=True,
                         reverse_lat=False)
