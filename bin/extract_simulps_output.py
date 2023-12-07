#!/usr/bin/env python

import sys
from pathlib import Path
from tomo2plt import io as PYTMIO

if len(sys.argv) != 3:
    print("USAGE: %s  SIMULPS_OUTPUT  CONFIG_PATH" %
          Path(sys.argv[0]).name)
    sys.exit()

simulout_path = Path(sys.argv[1])
assert simulout_path.exists()
config_path = Path(sys.argv[2])
assert config_path.exists()

# SPREAD function --> They call it “resolving width function” and the equation is shown in the appendix (eq A12).
configs = PYTMIO.read_configuration_file(config_path, check_version=True)
(grid_df, events_df, grid_matrix, event_matrix) = PYTMIO.read_simulps_output(
            simulout_path, configs)

print("DONE")
