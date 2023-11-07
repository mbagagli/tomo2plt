import sys
import numpy as np
import pandas as pd
from pathlib import Path


ORIGINLON = 10.5
ORIGINLAT = 46.0
DEGKMLAT = 111.15  # km  (short-distance conversion Y --> check SIMULPS output)
DEGKMLON = 77.466  # km  (short-distance conversion X --> check SIMULPS output)


def __split_string__(string, chunk_size):
    return [string[i:i+chunk_size] for i in range(0, len(string), chunk_size)]


def __get_arr__(file_path, matchstr):
    idx_match, out_arr = [], []
    fp = Path(file_path)
    # 1. Find idx
    with open(str(fp), "r") as IN:
        for xx, line in enumerate(IN):
            line_fields = line.strip().split()
            if not line_fields: continue
            # --------------------------------
            if line_fields[0] == matchstr:
                idx_match.append(xx+1)
    #
    with open(str(fp), "r") as IN:
        for xx, line in enumerate(IN):
            if xx in idx_match:
                newline = line[3:-1]
                line_fields = __split_string__(newline, 7)
                out_arr.extend([float(de) for de in line_fields])

    return out_arr


def __extract_grid_nodes__(file_path, geographic=True):
    xarr = __get_arr__(file_path, "xgrid")
    yarr = __get_arr__(file_path, "ygrid")
    zarr = __get_arr__(file_path, "zgrid")
    if geographic:
        xarr = [((_x*-1)/DEGKMLON)+ORIGINLON for _x in xarr]
        yarr = [(_y / DEGKMLAT)+ORIGINLAT for _y in yarr]
    #
    return (xarr, yarr, zarr)


def __read_block__(file_path, start_fields, end_fields, elem_per_line):
    # start_fields = ["FINAL", "P-VELOCITY", "MODEL"]
    # end_fields = ["FINAL", "P-VELOCITY", "MODEL"]

    OUTMAT, zval_list = [], []
    layer_matrix, row = [], []
    fp = Path(file_path)
    #
    print("Start-Fields:  %r" % start_fields)
    print("End-Fields:  %r" % end_fields)

    # ===============================================================
    with open(str(fp), "r") as IN:

        record_on = False

        for xx, line in enumerate(IN):

            line_fields = line.strip().split()
            if not line_fields: continue

            # ========================================================

            if line_fields[:3] == start_fields:
                record_on = True

            if line_fields[:3] == end_fields:
                break

            if record_on:
                if line_fields[0] == "layer":
                    if layer_matrix:
                        OUTMAT.append(layer_matrix)

                    zval = float(line_fields[-1])
                    zval_list.append(zval)
                    print("Layer:  %5.1f" % zval)
                    layer_matrix, row = [], []
                    continue

                else:
                    try:
                        _row = [float(de) for de in line_fields]
                        row.extend(_row)
                        if len(_row) != elem_per_line:
                            layer_matrix.append(row)
                            row = []
                    except ValueError:
                        # Cannot convert float
                        continue
                #
            #
        #
    OUTMAT.append(layer_matrix)
    OUTMAT = np.array(OUTMAT)
    return (OUTMAT, zval_list)


def __read_block_res__(file_path, start_fields, end_fields, elem_per_line):
    # start_fields = ["FINAL", "P-VELOCITY", "MODEL"]
    # end_fields = ["FINAL", "P-VELOCITY", "MODEL"]

    OUTMAT, zval_list = [], []
    layer_matrix, row = [], []
    fp = Path(file_path)
    #
    print("Start-Fields RESOLUTION:  %r" % start_fields)
    print("End-Fields RESOLUTION:  %r" % end_fields)

    # ===============================================================
    with open(str(fp), "r") as IN:

        record_on = False

        for xx, line in enumerate(IN):

            line_fields = line.strip().split()
            if not line_fields: continue

            # ========================================================

            if line_fields[:3] == start_fields:
                record_on = True

            if line_fields[:3] == end_fields:
                break

            if record_on:
                if line_fields[0] == "layer":
                    if layer_matrix:
                        OUTMAT.append(layer_matrix)

                    zval = float(line_fields[-1])
                    zval_list.append(zval)
                    print("Layer:  %5.1f" % zval)
                    layer_matrix, row = [], []
                    continue

                else:
                    try:
                        if len(line) == 262:
                            newline = line[1:-1]
                        else:
                            newline = line[:-1]
                        _row = [float(de.split(":")[-1])for de in
                                __split_string__(newline, 13)]
                        row.extend(_row)
                        if len(_row) != elem_per_line:
                            layer_matrix.append(row)
                            row = []
                    except ValueError:
                        # Cannot convert float
                        continue
                #
            #
        #
    OUTMAT.append(layer_matrix)
    OUTMAT = np.array(OUTMAT)
    return (OUTMAT, zval_list)

# =================================================================
# =================================================================
# =================================================================
# =================================================================
# =================================================================
# =================================================================
# =================================================================
# =================================================================


def create_matrix(xarr, yarr, zarr, vel_mat, delta_mat, res_mat, out_file="outgrid.txt"):
    data = []
    with open(out_file, "w") as OUT:
        for _z, zz in enumerate(zarr):
            for _y, yy in enumerate(yarr):
                for _x, xx in enumerate(xarr):
                    OUT.write("%11.5f  %11.5f  %8.3f  %6.3f  %6.3f  %5.3f\n" % (
                                xx, yy, zz,
                                vel_mat[_z, _y, _x],
                                delta_mat[_z, _y, _x],
                                res_mat[_z, _y, _x]))
                    # Append data to the list
                    data.append([xx, yy, zz,
                                 vel_mat[_z, _y, _x],
                                 delta_mat[_z, _y, _x],
                                 res_mat[_z, _y, _x]])
    df = pd.DataFrame(data, columns=['X', 'Y', 'Z', 'Velocity', 'Delta', 'RDE'])
    return df


def read_simulps_output(file_path, store_path):

    (xarr, yarr, zarr) = __extract_grid_nodes__(file_path, geographic=True)

    print("... Reading VELOCITY")
    VEL, depths = __read_block__(
                        file_path,
                        start_fields=["FINAL", "P-VELOCITY", "MODEL"],
                        end_fields=["OBSERVATION", "MATRIX", "-"],
                        elem_per_line=20)

    # ==============================  TO BE IMPLEMENTED SOON
    # print("")
    # print("... Reading DELTA-MIN1D")
    # MIN1D, depths = __read_block__(
    #                     file_path,
    #                     start_fields=["velocity", "values", "on"],
    #                     end_fields=["iteration", "step", "0"],
    #                     elem_per_line=20)

    print("")
    print("... Reading RES")
    RES, depths = __read_block_res__(
                        file_path,
                        start_fields=["RESOLUTION", ":", "GRIDOINT"],
                        end_fields=["Computation", "finished", "at"],
                        elem_per_line=20)

    # ================  REMOVE BOUNDARY NODES
    xarr = xarr[1:-1]
    yarr = yarr[1:-1]
    VEL = VEL[:, 1:-1, 1:-1]

    # ================  FAST MIN1D
    # -5.00      5.250
    #  2.50      5.800
    # 10.00      6.000
    # 20.00      6.200
    # 30.00      6.800
    # 45.00      7.800
    # 60.00      8.100
    # 90.00      8.150
    MIN1D_PROFILE = [5.250, 5.800, 6.000, 6.200, 6.800, 7.800, 8.100, 8.150]
    MIN1D_VOLUME = []
    for _vel in MIN1D_PROFILE:
        MIN1D_VOLUME.append(np.full((VEL.shape[1], VEL.shape[2]), _vel))
    MIN1D_VOLUME = np.array(MIN1D_VOLUME)
    #
    DELTA = ((VEL - MIN1D_VOLUME) / MIN1D_VOLUME) * 100
    assert VEL.shape == DELTA.shape == RES.shape

    outpd = create_matrix(xarr, yarr, depths, VEL, DELTA, RES, out_file=store_path)
    return outpd


if __name__ == "__main__":
    read_simulps_output(sys.argv[1], sys.argv[2])


# ==============================================================
# ================================ REFERENCES
# ==============================================================
# ==============================================================
    # print("")
    # print("... Reading KHIT")
    # KHIT, depths = __read_block__(
    #                     file_path,
    #                     start_fields=["OBSERVATION", "MATRIX", "-"],
    #                     end_fields=["DERIVATIVE", "WEIGHT", "SUM"],
    #                     elem_per_line=19,
    #                     is_resolution=False)
    # print("")
    # print("... Reading DWS")
    # DWS, depths = __read_block__(
    #                     file_path,
    #                     start_fields=["DERIVATIVE", "WEIGHT", "SUM"],
    #                     end_fields=["RESOLUTION", ":", "GRIDOINT"],
    #                     elem_per_line=19,
    #                     is_resolution=False)
