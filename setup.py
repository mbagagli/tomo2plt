from tomo2plt import __version__
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    required_list = f.read().splitlines()

setup(
    name="tomo2plt",
    version=__version__,
    author="Matteo Bagagli",
    author_email="matteo.bagagli@dst.unipi.it",
    description="Collection of modules and routines for handling SIMULPS outputs and plot them",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mbagagli/tomo2plt",
    python_requires='>=3.7',
    install_requires=required_list,
    setup_requires=['wheel'],
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        "Operating System :: MacOS",
        "Intended Audience :: Science/Research",
    ],
    include_package_data=True,
    zip_safe=False,
    scripts=[
        'bin/Make3DsimulPS_MOD.py',
        'bin/Make3DsimulPS_SYNTHETICS.py',
        'bin/extract_simulps_output.py',
        'bin/make_plots_general_DEPTH-SLICES.py',
        'bin/make_plots_general_DEPTH-SLICES_synthetics.py',
        'bin/make_plots_general_DEPTH-SECTIONS.py',]
)
