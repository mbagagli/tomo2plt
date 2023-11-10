#!/usr/bin/env python

# Make 3D MOD
#
#

import os
import sys
import argparse
import numpy as np


# =================== Globals
GRDFMT = " %.1f"
VELFMT = " %4.2f"

### SIMUL
# short distance conversion factors
#      one min lat   1.8525 km --> 1 degree (*60 min) 111.15 km
#      one min lon   1.2911 km --> 1 degree (*60 min) 77.466 km

DEGKMLAT = 111.15  # km
DEGKMLON = 77.466  # km

# =================== Parsing
if len(sys.argv) == 1:
    sys.stderr.write(("USAGE: %s --help"+os.linesep) % (
                             sys.argv[0].split(os.sep)[-1]
                    ))
    sys.exit()

# =================== Parsing input
parser = argparse.ArgumentParser(
                description=("Script to create a 3D MOD file for SIMULPS, "
                             "starting from 1D velocity profile"))
parser.add_argument('-x', '--xspacing', type=str, nargs='+',
                    help='km spacing for inversion grid. X direction')
parser.add_argument('-y', '--yspacing', type=str, nargs='+',
                    help='km spacing for inversion grid. Y direction')
parser.add_argument('-z', '--zspacing', type=str, nargs='+',
                    help='km spacing for inversion grid. Z direction')
parser.add_argument('-v', '--vellayers', type=str, nargs='+',
                    help='velocity values')
parser.add_argument('-b', '--spacingunit', type=float,
                    help="minimal units for the inversion grid [1.0 or 0.1]")
parser.add_argument('-t', '--modeltag', type=str, default=None,
                    help="String to append in MOD's file firstline.")
parser.add_argument('-o', '--outputfile', type=str, default=None,
                    help="output filename (path) of the MOD file.")
parser.add_argument('-g', '--origingrid', type=float, nargs='+',
                    default=[10.5, 46.0],
                    help=("origin of the MAP.grid file associated with the spacing. For GMT plot."
                          "It must correspond to the decimal degree conversion of the origin listed in STNS file!!!"))
parser.add_argument('-m', '--gridoutputfile', type=str, default=None,
                    help="origin of the MAP.grid file associated with the spacing. For GMT plot")
parser.add_argument('-f', '--modprofile', type=str, default=None,
                    help="input filename (path) of containing 'depth vel' values.")

parser.add_argument(
                '--flip_horizontal', action='store_true', dest='flip',
                help=("by default SIMULPS is WEST-POSITIVE! Make sure to create the X "
                      "direction with the correct way! If specified, the X array in input "
                      "will be flipped, and so the resulting grid"))


parser.set_defaults(outputfile="MOD")
parser.set_defaults(gridoutputfile="MOD.grid")
parser.set_defaults(spacingunit=1.0)
parser.set_defaults(flip=False)
parser.set_defaults(modeltag="Created with %s" % sys.argv[0].split(os.sep)[-1])

args = parser.parse_args()

# Checks
if not args.xspacing or not args.xspacing:
    sys.stderr.write(" X-Y spacing needed!" +
                     os.linesep)
    sys.exit()
if (not args.zspacing and not args.vellayers) and not args.modprofile:
    sys.stderr.write(" Z and V values OR input file ('depth vel') needed!" +
                     os.linesep)
    sys.exit()

print("Creating arrays")
#  ------------------  Create X
xvals = []
for xx in args.xspacing:
    fields = xx.split(":")

    if len(fields) == 1:
        # singleValue
        xvals.append(np.float(xx))
    elif len(fields) == 3:
        # numpyarray
        _tmp = np.arange(np.float(fields[0]),   # start
                         np.float(fields[1]),   # stop
                         np.float(fields[2]))   # increment
        xvals.extend(_tmp)
    else:
        raise ValueError("Either VAL or START:STOP:DELTA syntax allowed!")

if args.flip:
    print("WARNING: Flipping horizontal grid!")
    xvals = sorted([_xx * -1 for _xx in xvals])


#  ------------------  Create Y
yvals = []
for yy in args.yspacing:
    fields = yy.split(":")

    if len(fields) == 1:
        # singleValue
        yvals.append(np.float(yy))
    elif len(fields) == 3:
        # numpyarray
        _tmp = np.arange(np.float(fields[0]),   # start
                         np.float(fields[1]),   # stop
                         np.float(fields[2]))   # increment
        yvals.extend(_tmp)
    else:
        raise ValueError("Either VAL or START:STOP:DELTA syntax allowed!")


if args.modprofile:
    print(" ---> Model 'depth vel' Profile:  %s " % args.modprofile)
    #
    zvals, vvals = [], []
    with open(args.modprofile, "r") as IN:
        for modline in IN:
            modlist = modline.strip().split()
            zvals.append(np.float(modlist[0]))
            vvals.append(np.float(modlist[1]))
else:
    #  ------------------  Create Z
    zvals = []
    for zz in args.zspacing:
        fields = zz.split(":")

        if len(fields) == 1:
            # singleValue
            zvals.append(np.float(zz))
        elif len(fields) == 3:
            # numpyarray
            _tmp = np.arange(np.float(fields[0]),   # start
                             np.float(fields[1]),   # stop
                             np.float(fields[2]))   # increment
            zvals.extend(_tmp)
        else:
            raise ValueError("Either VAL or START:STOP:DELTA syntax allowed!")

    #  ------------------  Create V
    vvals = []
    for vv in args.vellayers:
        fields = vv.split(":")

        if len(fields) == 1:
            # singleValue
            vvals.append(np.float(vv))
        elif len(fields) == 3:
            # numpyarray
            _tmp = np.arange(np.float(fields[0]),   # start
                             np.float(fields[1]),   # stop
                             np.float(fields[2]))   # increment
            vvals.extend(_tmp)
        else:
            raise ValueError("Either VAL or START:STOP:DELTA syntax allowed!")

nx, ny, nz, nv = len(xvals), len(yvals), len(zvals), len(vvals)
if nz != nv:
    raise ValueError("Depth and Velocity arrays must be EQUAL SIZE !!! "
                     "DEPTH: %d  -  VEL:  %d" % (nz, nv))

# ==========================  BUILDING

print("Creating outputs")
with open(args.outputfile, "w") as OUT:
    # 1
    OUT.write((" %3.1f %2d %2d %2d      %s"+os.linesep) %
              (args.spacingunit, nx, ny, nz, args.modeltag))
    # 2
    for n, x in enumerate(xvals):
        OUT.write(GRDFMT % x)
        if n == nx-1:
            OUT.write(os.linesep)
    # 3
    for n, y in enumerate(yvals):
        OUT.write(GRDFMT % y)
        if n == ny-1:
            OUT.write(os.linesep)
    # 4
    for n, z in enumerate(zvals):
        OUT.write(GRDFMT % z)
        if n == nz-1:
            OUT.write(os.linesep)
    # 5
    OUT.write("  0  0  0"+os.linesep)

    # VEL
    for vv in vvals:
        for yy in range(ny):
            for xx in range(nx):
                OUT.write(VELFMT % vv)
            #
            OUT.write(os.linesep)

# ----

print("Creating MAP-GRIDPOINTS --> %s" % args.gridoutputfile)

originLON, originLAT = args.origingrid
yvals_grid = [(_y/DEGKMLAT)+originLAT for _y in yvals]
if args.flip:
    xvals_grid = [((_x*-1)/DEGKMLON)+originLON for _x in xvals]
else:
    xvals_grid = [(_x/DEGKMLON)+originLON for _x in xvals]

with open(args.gridoutputfile, "w") as OUT:
    for xg in xvals_grid:
        for yg in yvals_grid:
            OUT.write(("%9.5f  %9.5f"+os.linesep) % (xg, yg))

print("DONE!")
