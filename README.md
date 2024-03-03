# TOMO2PLT

This package contains all utilities to work with SIMULPS14/17 environments.
And is structured like this:

- `bin` folder: contains all executable and governor scripts for the package modules



### INSTALL

How to do it:

```
$ conda env create -f environment_file.yml
$ conda activate tomo2plt
$ pip install .
```

Ready to go.. now the executable in `bin` folder will be added to your path

### USAGE

```
$ extract_simulps_output.py  OUTPUT  TOMO2PLTCONFFILE
$ make_plots_general_DEPTH-SLICES.py  OUTPUT  TOMO2PLTCONFFILE
$ make_plots_general_DEPTH-SLICES_synthetics.py  OUTPUT  TOMO2PLTCONFFILE
$ make_plots_general_DEPTH-SECTIONS.py  OUTPUT  TOMO2PLTCONFFILE
```
