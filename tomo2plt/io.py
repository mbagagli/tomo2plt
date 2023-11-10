import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

from tomo2plt import __version__


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



class MOD(object):
    def __init__(self, modfile, origin_lon=None, origin_lat=None):
        """ Class model for MOD file SIMULPS
        ### SIMUL
        # short distance conversion factors
        #      one min lat   1.8525 km --> 1 degree (*60 min) 111.15 km
        #      one min lon   1.2911 km --> 1 degree (*60 min) 77.466 km
        """

        # Constants
        self.GRDFMT = " %.1f"
        self.VELFMT = " %4.2f"
        self.DEGKMLAT = 111.15  # km for degree --> Check SIMULPS
        self.DEGKMLON = 77.466  # km for degree --> Check SIMULPS
        # Attributes
        self.filepath = None
        self.general = {}
        self.xarr = []
        self.yarr = []
        self.zarr = []
        self.matrix = []
        self.velocities = []
        self.origin = []
        self.lon_grid = []
        self.lat_grid = []
        if not origin_lon or not origin_lat:
            print("WARNING:  geographical coordinates will not be possible "
                  "without a decimal degree origin")
        else:
            self.origin = (origin_lon, origin_lat)
        #
        self.import_file(modfile)

    def import_file(self, filepath):
        """ Simply add """
        with open(filepath, "r") as IN:
            for idx, line in enumerate(IN):

                lines = line.strip().split()
                if idx == 0:
                    self.general["bldunits"] = float(lines[0])
                    self.general["nx"] = int(lines[1])
                    self.general["ny"] = int(lines[2])
                    self.general["nz"] = int(lines[3])

                elif idx == 1:
                    self.xarr = np.array([float(xx) for xx in lines])

                elif idx == 2:
                    self.yarr = np.array([float(yy) for yy in lines])

                elif idx == 3:
                    self.zarr = np.array([float(zz) for zz in lines])

                elif (idx-5)%self.general["ny"] == 0:
                    self.velocities.append(float(lines[0]))

                else:
                    pass

        # Close File, create matrix  --> Zx(YxX)
        _mat = []
        for _z in range(self.general["nz"]):
            _mat.append(np.full((self.general["ny"], self.general["nx"]),
                                self.velocities[_z], dtype="float32"))

        self.matrix = np.array(_mat)

    def calculate_lon_lat_grid(self):
        if self.origin:
            self.lon_grid = [(_x/self.DEGKMLON)+self.origin[0] for _x in self.xarr]
            self.lat_grid = [(_y/self.DEGKMLAT)+self.origin[1] for _y in self.yarr]
        else:
            logger.error("MISSING ORIGIN!")

    def set_origin(self, lon, lat):
        self.origin = (float(lon), float(lat))

    def get_origin(self):
        if self.origin:
            return self.origin
        else:
            return None

    def _reverse_east_west(self):
        """ If called it fill flip the X direction of the matrix """
        if len(self.matrix) == 0:
            raise AttributeError("Missing matrix! Please import file first")
        #
        self.matrix = self.matrix[:, :, ::-1]

    def _reverse_north_south(self):
        """ If called it fill flip the Y direction of the matrix """
        if len(self.matrix) == 0:
            raise AttributeError("Missing matrix! Please import file first")
        #
        self.matrix = self.matrix[:, ::-1, :]

    def write_mod(self, outfile="MOD", reverse_lon=False, reverse_lat=False):
        """ Writing out a MOD file for inversions.
            if reverse_lon specified the matrix will be flipped on the LON-axis
        """
        if reverse_lon:
            self._reverse_east_west()
        if reverse_lat:
            self._reverse_north_south()
        #
        print("Creating outputs")
        with open(outfile, "w") as OUT:
            # 1
            OUT.write((" %3.1f %2d %2d %2d      Created by TOMO2PLT - v_%s"+os.linesep) %
                      (self.general["bldunits"],
                       self.general["nx"],
                       self.general["ny"],
                       self.general["nz"], __version__))
            # 2
            for n, x in enumerate(self.xarr):
                OUT.write(self.GRDFMT % x)
                if n == self.general["nx"]-1:
                    OUT.write(os.linesep)
            # 3
            for n, y in enumerate(self.yarr):
                OUT.write(self.GRDFMT % y)
                if n == self.general["ny"]-1:
                    OUT.write(os.linesep)
            # 4
            for n, z in enumerate(self.zarr):
                OUT.write(self.GRDFMT % z)
                if n == self.general["nz"]-1:
                    OUT.write(os.linesep)
            # 5
            OUT.write("  0  0  0"+os.linesep)

            # VEL
            for zz in range(self.general["nz"]):
                for yy in range(self.general["ny"]):
                    for xx in range(self.general["nx"]):
                        OUT.write(self.VELFMT % self.matrix[zz, yy, xx])
                    #
                    OUT.write(os.linesep)

    def export_grid(self, outfile="MOG.grid", in_coordinates=False):
        """ Writing out an ascii file (x,y) for plot inversion grid"""

        if in_coordinates and self.origin:
            print("Creating MAP-GRIDPOINTS ... Geographical")
            xvals_grid = [(_x/self.DEGKMLON)+self.origin[0] for _x in self.xarr]
            yvals_grid = [(_y/self.DEGKMLAT)+self.origin[1] for _y in self.yarr]
        else:
            print("Creating MAP-GRIDPOINTS ... Cartesian")
            xvals_grid, yvals_grid = self.xarr, self.yarr

        with open(outfile, "w") as OUT:
            for xg in xvals_grid:
                for yg in yvals_grid:
                    OUT.write(("%9.5f  %9.5f"+os.linesep) % (xg, yg))


# if __name__ == "__main__":
#     read_simulps_output(sys.argv[1], sys.argv[2])


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
