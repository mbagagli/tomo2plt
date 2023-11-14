#!/usr/bin/env python

import sys
import shapefile
from tqdm import tqdm
from tomo2plt import core as PYTMCR
import matplotlib.pyplot as plt
import numpy as np

# It's always  GOURAUD  SHADING !
X_KM_SLICE = 15       # cell size in longitude (km)
Y_KM_SLICE = 15       # cell size in latitude (km)
SMOOTH_SLICE = 0.9    # GAUSSIAN KERNEL value or None/False
SLICE_PLOT_DIR = "./SLICES_%.01fx%.01f_%.01f" % (X_KM_SLICE, Y_KM_SLICE, SMOOTH_SLICE)

MASK_RDE = 0.2  # RDE value or None/False

# ======================================  Init
grid = PYTMCR.TomoGrid(sys.argv[1])


# ======================================  DEPTH - SLICE
start_point = (1, 41, 0.0)
end_point = (21, 51, 70.0)

for xx, _dep in enumerate(tqdm(np.arange(0, 60, 1))):
    (fig, ax) = grid.plot_depth_slice(
                        start_point, end_point,
                        increment_x_km=X_KM_SLICE, increment_y_km=Y_KM_SLICE,
                        depth=_dep,  what="delta",
                        smooth=SMOOTH_SLICE, mask_rde=MASK_RDE,
                        interpolate="linear")
    plt.title('DEPTH SLICE @ %.2f km' % _dep)

    # ---------------------- COASTLINES
    shapefile_path = "./datas/coastlines/ne_10m_coastline.shp"
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

    # ----------------------Finisci e SAVE
    ax.set_xlim(start_point[0], end_point[0])
    ax.set_ylim(start_point[1], end_point[1])
    ax.set_aspect(1.5)

    plt.savefig(('%s/depth_SLICE_plot_%05.2f_km.png' % (SLICE_PLOT_DIR, _dep)),
                format='png', bbox_inches='tight', dpi=310)
