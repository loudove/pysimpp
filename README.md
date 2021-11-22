pysimpp
=======

pysimpp project is a collection of scripts for post processing molecular simulations trajectories.

The calculation intensive parts of the scipts are based on python binds of fortran and c/c++ imlementations. Python offers a welth of easyly accessible, high quality libraries reducing drasticly the imlementation time.

The infrastructure for accessing the trajectory data is mainly based on the MDAnalysis project. Numpy, networkx and pyvoro are also central modules used heavily in pysimpp scripts.

The reader is essentially a wrapper of MDAnalysis.Universe with an API that is rather handy from classical materials modeling perspective. In the case of LAMMPS molecular dynamics simulator, a more elaborate parser is available. Therefore, in principle, pysimpp's scripts could process trajectories from the simulation engines supported by MDAnalysis either directly or after minor changes.

Installation
------------

[//]: # (Recommended - installation via `pip`: pip3 install pysimpp)

The installation from source is the same as for any other python module. Clone pysimpp repository:

    git clone https://github.com/loudove/pysimpp.git

enter the pysimpp directory and use pip:
  
    pip3 install .

or `setup.py`:

    python3 setup.py install
    
to install pysimpp system-wide. Otherwise use:  

    pip3 install . --user

or:

    python setup.py install --user

to install it only for the current user.

Requirements
------------

Some of the functionality of pysimm requir the installation of the [pyvoro](https://github.com/joe-jordan/pyvoro) module. It is recommended to install a modified version of pyvoro where memory requirements for large systems have been significantly reduced . You should clone the feature/python3 branch for a fork of the project found in https://github.com/loudove/pyvoro

    git clone --branch feature/python3 https://github.com/loudove/pyvoro.git

and follow the installation instructions.

Use
---

After installing pysimm you can list the available pysimpp commands with:

    pysimpp -h

To access detailed help on a specific topic use:

    pysimpp command -h

A more detail description of the functionality for specific commands ([cluster](./doc/cluster.md)) can be found in [doc](./doc) folder.
