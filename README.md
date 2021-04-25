pysimpp
=======

pysimpp project is a collection of scripts for post processing molecular simulations trajectories.

The calculation intensive parts of the scipts are based on python binds of fortran and c/c++ imlementations. Python offers a welth of easyly accessible, high quality libraries reducing drasticly the imlementation time.

The infrastructure for accessing the trajectory data is mainly based on the MDAnalysis project. Numpy, networkx and pyvoro are also central modules used heavily in pysimpp scripts.

The reader is essentially a wrapper of MDAnalysis.Universe with an API that is rather handy from classical materials modeling perspective. In the case of LAMMPS molecular dynamics simulator, a more elaborate parser is available. Therefore, in principle, pysimpp's scripts could process trajectories from the simulation engines supported by MDAnalysis either directly or after minor changes.

Installation
------------

Recommended - installation via `pip`:

    pip3 install pysimpp

Installation from source is the same as for any other python module. Issuing 
  
    pip3 install .

or:

    python3 setup.py install
    
will install pysimpp system-wide, while 

    pip3 install . --user

or

    python setup.py install --user

will install it only for the current user.

