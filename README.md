pysimpp
=======

pysimpp project is a collection of scripts for post processing molecular simulations trajectories.

Python offers a wealth of easily accessible, high quality libraries reducing drastically the implementation time. Nevertheless, the calculation intensive parts of the scripts are based on python binds of Fortran and c/c++ implementations.

The infrastructure for accessing the trajectory data is mainly based on the [MDAnalysis](https://www.mdanalysis.org/) project. [Numpy](https://numpy.org/), [NetworkX](https://networkx.org/) and [pyvoro](https://github.com/joe-jordan/pyvoro) are also central modules used heavily in pysimpp scripts.

The basic trajectory reader is essentially a wrapper of MDAnalysis.Universe with an API that is rather handy from classical materials modeling perspective. In the case of LAMMPS molecular dynamics simulator, a more elaborate parser is available. Therefore, in principle, pysimpp's scripts could process trajectories from the simulation engines supported by MDAnalysis either directly or after minor changes.

Requirements
------------

Some of the functionality in pysimm, requires the installation of the [pyvoro](https://github.com/joe-jordan/pyvoro) module. It is recommended to install a modified version of pyvoro where the memory requirements for large systems are significantly reduced. You should clone the feature/python3 branch from a [fork](https://github.com/loudove/pyvoro) of the project:

    git clone --branch feature/python3 https://github.com/loudove/pyvoro.git

and follow the installation instructions.

To install pysimpp on windows hosts natively, a working version of `f2py` is needed. Otherwise, [WSL](https://docs.microsoft.com/en-us/windows/wsl/install) is an excalent alternative.

Installation
------------

[//]: # (Recommended - installation via `pip`: pip3 install pysimpp)

The installation from source is the same as for any other python module. Clone pysimpp repository:

    git clone https://github.com/loudove/pysimpp.git

enter the pysimpp directory and use pip:

    pip3 install .

to install pysimpp system-wide. Otherwise use:

    pip3 install . --user

to install it only for the current user.

Use
---

After installing pysimm you can list the available pysimpp commands with:

    pysimpp help

To access detailed help on a specific command use:

    pysimpp command -h

More details on [commands](./doc/README.md) with a more complex set of arguments can be found in the [doc](./doc) folder.

TODO
---

- Complete and improve the documentation.

- Add automated testing for python.
