Clusters
========

Trace the clusters formed in the simulation cell based on close contact analysis preformed with voronoi tesselation as implemented in [voro++](http://math.lbl.gov/voro++/) and its python bindings [pyvoro](https://github.com/joe-jordan/pyvoro).

Positional arguments
--------------------
`path`: the path to the simulation trajectory file. A topology file with the same based name should reside on the same directory, otherwise it should be provided with the `-topo` option.


Options
-------
- `-h, --help`:
show this help message and exit.

- `-start n` :
start processing form timestep n [inclusive] (step based for lammps and time based for gromacs).

-  `-end n` :
stop processing at timestep n [inclusive] (step based for lammps and time based for gromacs).

-  `-every n`: processing frequency (every n timesteps, step based for lammps and time based for gromacs).

-  `-vis n`:
visualization frequency (save gro files for the traced clusters every n timesteps iteration based for lammps and time based for gromacs).

-  `-molecules molid range`:
molecules to be used. A comma seperated list with the ranges of molecules ids, e.g., `1,2,3` or `1:10,20,30:100`

-  `-molnames molnames range`: 
types/names of clusters constituent molecules. A comma separated list with the names of the molecules, e.g., `Na,Cl`

-  `-ends list of end atoms`:
the pairs of atoms defining the end-to-end vectors of the molecular species participating in a cluster (see `-molnames`). For example, if the end-to-end vectors of clusters' constituent molecular types TIC and TOC are defined by atoms (1, 16) and (3,25), respectively, the arguments could be `TIC:1,16@TOC:3,25`.

-   `-excluded types range`:
the atom types to be exluded system wide from the closecontact analysis. A comma seperated list with the atoms' type name should be provided e.g.`HA,HW`.

-   `-specific atoms for each molname`:
the names of the atoms to be considered in the close contact analysis. A comma separated list for each molecular type in the `-molnames` argument should be provided. If a wildcard `"*"` is given for a residue then all the atoms of the molecular species will be considered. For example, if A, and B are the clusters' constituent molecular types, the argument could look like `*:C1,C2,C3` specifying  that all the atoms of species A and only the atoms C1,C2,C3 of species B should be considered in the analysis.

-   `-radii file with atoms' type/elemet radii`
 the file with the radii of the atoms. It can be element or type based. The first line of the file contains the keywords `(element|type) (r|d)`; the first, specifies the atom type identifier and the second if the radius (`r`) or the diameter (`d`) is given for each type. The rest of the lines contain the (`type`, `radius`) pairs. The `type` could be either a number (type id) or a string (type name).

-   `-hist list of neighbor group pairs`
the histograms of the number of neighbor pairs to be calculated. Each pair consists of a groups of atoms and a list of species. A list of pairs separated with "@" should be provided. Each pair contains the information for a histogram in the format `GROUPNAME:MOLNAME:ATOMSLIST:RESIDUESLIST`, where:
    - GROUPNAME:  the name of the group
    - MOLNAME:    the name of the molecule where the group belongs (all the atoms
                should belong to the same molecule)
    - ATOMSLIST:  a comma-separated list of atoms define the group. The atoms should belong to the same residue.
    - SPECIESLIST: a comma-separated list with molecule types. One histogramwill be calculated for each molecular type in the list

    The histogram is written in the file GROUPNAME_SPECIES_neighbors.dat in the simulation directory. For example, the argument: `C1:CTAC:C1,H31,H32,H33:CTAC,CL,SOL@C2:CTAC:C2,H1,H2:CTAC,CL,SOL` defines two groups and the ollowing pairs `(C1,H31,H32,H33)-CTAC`, `(C1 H31 H32 H33)-CL`, `(C1 H31 H32 H33)-SOL`, `(C2,H1,H2)-CTAC`, `(C2,H1,H2)-CL`, and `(C2,H1,H2)-SOL`; the files `C1_CTAC_neighbors.dat`, `C1_CL_neighbors.dat`, `C1_SOL_neighbors.dat`, `C2_CTAC_neighbors.dat`, `C2_CL_neighbors.dat`, and `C2_SOL_neighbors.dat`, will be wirtten in the simulation directory, respectively.

-  `--nbhist`:
for the molecular types participating in a cluster, calculate the histograms of their total number of neighbors with the other molecular types in the system. For example, consider a system consisting of molecular types A, B, and C. If he constituent species of the clusters are A and B, the script will calculate he histograms of the total number of neighbors pairs for `A-A`, `A-B`, `A-C`, `B-A`, `B-B` and `B-C` pairs. For each pair, the histogram will be written in `A_B_neighbors.at` file in the simulation directory.

-  `-hist2d list of pairs of residues`
for a cluster, calculates the conditional probability of having n neighbor species of type A-A given that m neighbor species of type B-B exist. The argument is a  column separated list of pairs, e.g., `A,B:C,D`. For each pair, the conditional probability is written in the file `A_B_TOTAL2D_neighbors.dat` in the simulation directory.

-  `-hist3d list of triplets of residues`
for a cluster, calculate the conditional probability of having n neighbor residues of type C given that m neighbor residues of type B exist for all the possible values of the number of neighbor residues of type A. The argument is a column separated list of triplet, e.g., `A,B,C:D,E,F`. For each triplet, the conditional probability is written in the file `A_B_C_TOTAL3D_neighbors.dat` in the simulation directory.

-  `-profiles list of atom names`
calculate the number density profile of the given groups of atoms. The distance considered for the profiles depends on the position of the atoms in the micelle and eventually from the shape of the cluster. If the atom is located in a spherical micelle or in the spherical caps of an elongated/rod-like micelle, the distance from the center of the sphere is used. If the atoms belong to a cylindrical column or in the body of an elongated micelle, the length of its projection on the column axis is taken. The argument is a set of lists of comma-separated atoms' names lists separated with "@". For example the argument `HEAD:CTAC:C17,N,C18,C19@WAT:SOL:OW` defines two named lists (groups); `HEAD` consists of atoms `C17`,`N`,`C18`, and `C19` that belong to `CTAC` molecules and `WAT` consists of `OW` atoms belong to `SOL` molecules. Atoms specified with the `-excluded` argument will be excluded also here. For each list, the density profile will be written in the file `profile_{listname}_{shape}.prof` in the simulation directory where list name is the name of the group of atoms and shape the shape of the cluster. Three types of clusters are considered:  spherical (`s`), cylindrical infinite periodic (`c`), and wedged/elongated (`e`). Therefore file `HEAD_s.prof` corresponds to the density profiles of atoms `C17`,`N`,`C18`, and `C19` belong to `CTAC` molecules, with respect to the center of mass for spherical clusters.

-  `-phist list of bin lengths`
the properties for which their histograms will be calculated. The name of the property together with the bin length for the histogram should be given. The following properties are supported:
    b : asphericity for both clusters and their constituent molecular species
    c : acylindricity for both clusters and their constituent molecular species
    sqk : anisotropy for both clusters and their constituent molecular species
    sqrg : square radius of gyration for both clusters and their constituent molecular species
    sqee : end-to-end distance for clusters' constituent molecular species
    bbox : bounding box
    qlong : global order
    qlocal : local order
The bin length of a poperty of a molecular species is named by apending `_{SPECIES SNAME}` at its name. For example, with argument `c:0.01,c_CTAC:0.05`, the bin length of the acylindricity distribution of the clusters is set to 0.01 and for the CTAC molecules to 0.05.
The special keywords `all`, `clusters`, and `species` can be given instead of a property's name. In this case, `all` the available properties, the properties of the traced `clusters`, or the properties of the molecular `species` will be calculated, respectively, and the corresponding histograms will be calculated.

-  `-phist2d list of properties`
provide the properties for which the conditional probability of having a specific value given the size of the cluster, will be calculated.

-  `-voropp voro++ executable`
provide the voro++ executable to be used instead of pyvoro ython module.
