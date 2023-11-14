import os
import sys
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

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


def __read_block__(file_path, start_fields, end_fields, elem_per_line,
                   skip=0, trim=0):
    # start_fields = ["FINAL", "P-VELOCITY", "MODEL"]
    # end_fields = ["FINAL", "P-VELOCITY", "MODEL"]

    OUTMAT, zval_list = [], []
    layer_matrix, row = [], []
    fp = Path(file_path)
    #
    print("Start-Fields:  %r" % start_fields)
    print("End-Fields:  %r" % end_fields)

    start_fields_idx, end_fields_idx = [], []
    with open(str(fp), "r") as IN:
        for _x, _line in enumerate(IN):
            line_fields = _line.strip().split()
            if line_fields[:3] == start_fields:
                start_fields_idx.append(_x)
            elif line_fields[:3] == end_fields:
                end_fields_idx.append(_x)
            else:
                continue
    #
    assert len(start_fields_idx) >= 1
    assert len(end_fields_idx) >= 1
    start_fields_idx = start_fields_idx[-1]
    end_fields_idx = end_fields_idx[-1]
    assert start_fields_idx < end_fields_idx

    # ===============================================================
    with open(str(fp), "r") as IN:

        record_on = False

        for xx, line in enumerate(IN):

            line_fields = line.strip().split()
            if not line_fields: continue

            # ========================================================

            if xx == start_fields_idx:
                record_on = True

            if xx == end_fields_idx:
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
                        # --- Skip FIRST n element
                        if int(skip) > 0:
                            _row = _row[skip:]

                        # --- Remove LAST n element
                        if int(trim) > 0:
                            _row = _row[:-trim]

                        # --- Append line to eixtingrow
                        row.extend(_row)
                        if len(_row) != (elem_per_line - skip - trim):
                            # --- Last line belonging to row
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


def __extract_events__(file_path):
    """
               FINAL LOCATIONS

     YRMODY HRMN  SEC  LATITUDE LONGITUDE  DEPTH    MAG  NO RMSRES     x      y      z
    Full line:
               |         |         |         |         |         |         |         |
     01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678
    '    1 960116 07 0 20.69 46N 51.87   7E 23.48  16.35   2.60  17  0.10 238.94  96.08  16.35'

    """

    def __degreedecimal2degreeminutes__(infloat):
        """ de """
        degr = np.int(str(infloat).split(".")[0])
        mint = np.float("0." + str(infloat).split(".")[1]) * 60.0
        return degr, mint

    def __degreeminute2decimaldegree__(instr):
        """
        Degrees Minutes.m to Decimal Degrees
        .d = M.m / 60
        Decimal Degrees = Degrees + .d
        """
        _ll = instr.strip()
        #
        _dl = float(_ll[-6:])
        _comp = _ll[-7]
        _deg = float(_ll[:-7])
        #
        decdeg = _deg + _dl/60.0
        if _comp.lower() in ("s", "w"):
            decdeg = -decdeg
        #
        return decdeg

    # ===============================================================

    __events_dictionary_keys__ = [
        "OT",
        "LONGITUDE",
        "LATITUDE",
        "DEPTH",
        "MAGNITUDE",
        "NPHASES",
        "RMS",
        "X", "Y", "Z"
    ]
    OUTDICT = {}
    fp = Path(file_path)

    # ===============================================================
    with open(str(fp), "r") as IN:

        record_on = False
        dict_idx = 0
        for xx, line in enumerate(IN):
            line_fields = line.strip().split()
            if not line_fields: continue

            # ========================================================

            if line_fields[:3] == ["YRMODY", "HRMN", "SEC"]:
                record_on = True
                continue

            if line_fields[0] == "RMSALL:":
                break

            if record_on:
                _indict = {kk: None for kk in __events_dictionary_keys__}
                OUTDICT[str(dict_idx)] = _indict

                # next line to avoid cross-century issues
                if int(line[6]) > 3:
                    _yr = 1900
                else:
                    _yr = 2000

                # --- Origin TIME
                year = _yr + int(line[6:8])
                month = int(line[ 8:10])
                day = int(line[10:12])

                hour = int(line[13:15])
                if hour == 24:
                    day += 1
                    hour = 0
                minute = int(line[15:17])
                if minute == 60:
                    hour += 1
                    minute = 0

                seconds = float(line[17:23])
                if seconds == 60:
                    minute += 1
                    seconds = 0

                ot = "%04d-%02d-%02d %02d:%02d:%05.2f" % (
                            year, month, day, hour, minute, seconds)

                # --------------- lon / lat / dep

                lat = __degreeminute2decimaldegree__(line[24:33])
                lon = __degreeminute2decimaldegree__(line[35:44])
                dep = float(line[45:51])
                mag = float(line[53:58])
                nop = int(line[58:62])
                rms = float(line[63:68])

                xxx = float(line[68:75])
                yyy = float(line[75:82])
                zzz = float(line[82:-1])

                # ===========================================  POPULATE
                OUTDICT[str(dict_idx)]["OT"] = datetime.strptime(
                                                   ot, '%Y-%m-%d %H:%M:%S.%f')
                OUTDICT[str(dict_idx)]["LONGITUDE"] = lon
                OUTDICT[str(dict_idx)]["LATITUDE"] = lat
                OUTDICT[str(dict_idx)]["DEPTH"] = dep
                OUTDICT[str(dict_idx)]["MAGNITUDE"] = mag
                OUTDICT[str(dict_idx)]["NPHASES"] = nop
                OUTDICT[str(dict_idx)]["RMS"] = rms
                OUTDICT[str(dict_idx)]["X"] = xxx
                OUTDICT[str(dict_idx)]["Y"] = yyy
                OUTDICT[str(dict_idx)]["Z"] = zzz
                #
                dict_idx += 1

    df = pd.DataFrame.from_dict(OUTDICT, orient="index")
    return (df, OUTDICT)

# =================================================================
# =================================================================
# =================================================================
# =================================================================
# =================================================================
# =================================================================
# =================================================================
# =================================================================


def read_configuration_file(config_path):
    # Can be either a `Path` or `str` instance
    with open(str(config_path), 'r') as yaml_file:
        configs = yaml.load(yaml_file, Loader=yaml.FullLoader)
    return configs


def create_matrix(xarr_geo, yarr_geo, zarr_geo, xarr, yarr,
                  vel_mat, delta_mat, res_mat, khit_mat, dws_mat,
                  out_file="outgrid.txt"):
    data = []
    with open(out_file, "w") as OUT:
        for _z, zz in enumerate(zarr_geo):
            for _y, yy in enumerate(yarr_geo):
                for _x, xx in enumerate(xarr_geo):

                    OUT.write("%11.5f  %11.5f  %8.1f  %8.1f  %8.3f  %6.3f  %6.3f  %5.3f  %5d  %5d\n" % (
                                xx, yy, xarr[_x], yarr[_y], zz,
                                vel_mat[_z, _y, _x],
                                delta_mat[_z, _y, _x],
                                res_mat[_z, _y, _x],
                                int(khit_mat[_z, _y, _x]),
                                int(dws_mat[_z, _y, _x])))
                    # Append data to the list
                    data.append([xx, yy, xarr[_x], yarr[_y], zz,
                                 vel_mat[_z, _y, _x],
                                 delta_mat[_z, _y, _x],
                                 res_mat[_z, _y, _x],
                                 khit_mat[_z, _y, _x],
                                 dws_mat[_z, _y, _x]])

    data_types = {
        'LONGITUDE': 'float32',
        'LATITUDE': 'float32',
        'X': 'float32',
        'Y': 'float32',
        'DEPTH': 'float32',
        'Velocity': 'float32',
        'Delta': 'float32',
        'RDE': 'float32',
        'KHIT': 'int32',
        'DWS': 'int32',
    }

    # Apply the type conversion
    df = pd.DataFrame(
                data,
                columns=data_types.keys())
    df = df.astype(data_types)

    return df


def read_simulps_output(file_path, store_path, config_path):
    store_path = Path(store_path)
    file_path = Path(file_path)
    configs = read_configuration_file(config_path)
    #
    (xarr, yarr, zarr) = __extract_grid_nodes__(
                                str(file_path), geographic=True)
    (xxx, yyy, zzz) = __extract_grid_nodes__(
                                str(file_path), geographic=False)

    print("... Reading VELOCITY")
    VEL, depths = __read_block__(
                        str(file_path),
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
    # Already without boundary nodes removed
    RES, depths = __read_block_res__(
                        str(file_path),
                        start_fields=["RESOLUTION", ":", "GRIDOINT"],
                        end_fields=["Computation", "finished", "at"],
                        elem_per_line=20)

    print("")
    print("... Reading RES")
    KHIT, depths = __read_block__(
                        str(file_path),
                        start_fields=["OBSERVATION", "MATRIX", "-"],
                        end_fields=["DERIVATIVE", "WEIGHT", "SUM"],
                        elem_per_line=21, skip=1, trim=0)

    print("")
    print("... Reading DWS")
    DWS, depths = __read_block__(
                        file_path,
                        start_fields=["DERIVATIVE", "WEIGHT", "SUM"],
                        end_fields=["RESOLUTION", ":", "GRIDOINT"],
                        elem_per_line=19, skip=1, trim=0)

    # ================  REMOVE BOUNDARY NODES and fix shapes
    xarr = xarr[1:-1]
    yarr = yarr[1:-1]
    VEL = VEL[:, 1:-1, 1:-1]
    KHIT = KHIT[:, :, :-1]  # removing last element
    DWS = DWS[:, :, :-1]  # removing last element

    MIN1D = configs["1d_input_model"]
    MIN1D_PROFILE = []
    for _dpt in depths:
        for (_depth, _pvalue) in MIN1D:
            if _dpt == _depth:
                MIN1D_PROFILE.append(_pvalue)

    MIN1D_VOLUME = []
    for _vel in MIN1D_PROFILE:
        MIN1D_VOLUME.append(np.full((VEL.shape[1], VEL.shape[2]), _vel))
    MIN1D_VOLUME = np.array(MIN1D_VOLUME)
    #
    DELTA = ((VEL - MIN1D_VOLUME) / MIN1D_VOLUME) * 100
    assert VEL.shape == DELTA.shape == RES.shape == KHIT.shape == DWS.shape

    matrixpd = create_matrix(xarr, yarr, depths, xxx, yyy,
                             VEL, DELTA, RES, KHIT, DWS,
                             out_file=str(store_path))
    matrixpd.to_csv(
            str(store_path.parent / (store_path.stem + '_pandas.csv')),
            float_format='%.5f',
            index=False)

    if configs["extract_events"]:
        print("")
        print("... Extracting EVENTS")
        eventpd, EVENTS = __extract_events__(file_path)
        eventpd.to_csv(
                str(store_path.parent / (store_path.stem + '_events.csv')),
                index=False)

    return (matrixpd, eventpd)


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
            OUT.write((" %3.1f %2d %2d %2d      Created by TOMO2PLT - v%s"+os.linesep) %
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

# # ========================================================
# # =======================  OLD-REFERENCE
# def create_matrix(xarr, yarr, zarr,
#                   vel_mat, delta_mat, res_mat,
#                   out_file="outgrid.txt"):
#     data = []
#     with open(out_file, "w") as OUT:
#         for _z, zz in enumerate(zarr):
#             for _y, yy in enumerate(yarr):
#                 for _x, xx in enumerate(xarr):
#                     OUT.write("%11.5f  %11.5f  %8.3f  %6.3f  %6.3f  %5.3f\n" % (
#                                 xx, yy, zz,
#                                 vel_mat[_z, _y, _x],
#                                 delta_mat[_z, _y, _x],
#                                 res_mat[_z, _y, _x]))
#                     # Append data to the list
#                     data.append([xx, yy, zz,
#                                  vel_mat[_z, _y, _x],
#                                  delta_mat[_z, _y, _x],
#                                  res_mat[_z, _y, _x]])

#     df = pd.DataFrame(
#                 data,
#                 columns=['X', 'Y', 'Z', 'Velocity', 'Delta', 'RDE'])
#     return df
# #
#     matrixpd = create_matrix(xarr, yarr, depths, VEL, DELTA, RES,
#                              out_file=str(store_path))
