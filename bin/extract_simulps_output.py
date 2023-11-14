#!/usr/bin/env python

import sys
from pathlib import Path
from tomo2plt import io as PYTMIO

if len(sys.argv) != 4:
    print("USAGE: %s  SIMULPS_OUTPUT  TOMO2PLT_FILE  CONFIG_PATH" %
           Path(sys.argv[0]).name)
    sys.exit()

# ./SCRIPTNAME simulout_path store_path
(grid_df, events_df) = PYTMIO.read_simulps_output(*sys.argv[1:])
print("DONE")
