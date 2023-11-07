from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    required_list = f.read().splitlines()

setup(
    name="tomotools",
    version="0.1.2",
    author="Matteo Bagagli",
    author_email="matteo.bagagli@ingv.it",
    description="Collection of modules and bin for handling of SIMULPS and FMTOMO formats",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mbagagli/tomotools",
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
        'bin/velest2simulps_converter.py',
        'bin/Make3DsimulPS_MOD.py',
        'bin/Make3DsimulPS_SYNTHETICS.py',
        'bin/simulps2fmtomo.py',
        'bin/CreateAnalyze_EQKS.py',
        'bin/CreateError_df.py',
        'bin/CreateGMTfiles_tomotools.py']
)
