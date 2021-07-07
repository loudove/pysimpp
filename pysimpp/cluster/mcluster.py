# -*- coding: utf-8 -*-

import os

filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
    with open(filename) as fobj:
        startup_file = fobj.read()
    exec(startup_file)
    import os

import sys
from collections import defaultdict, Counter
import subprocess
import numpy as np
from scipy import stats
import networkx as nx

import pysimpp.readers
# import pysimpp.readers.lammpsreader as lmp
# import pysimpp.readers.groreader as gmp

from pysimpp.utils.utils import chk_filename, parse_radii, isrange, isfile, IsList, IsListOfList, IsListOfNamedList
from pysimpp.utils.statisticsutils import Binning, Histogram, Histogram2D, Histogram3D
from pysimpp.fastpost import fastunwrapv, order_parameter # pylint: disable=no-name-in-module
from pysimpp.utils.clusterutils import Node, Connection
import pysimpp.utils.voropputils as voropputils
from .mclusterutils import MClusterTracker, MCluster # pylint: disable=import-error
from .mclusterlog import ClPropertiesLog, ClMolecularLog, ClOrderLog, ClDetailsLog # pylint: disable=import-error

# global flags
_critical_size = 2
_wholeused = True
_debug = False
_extend = True  # use the extended pyvoro version
_runext = True  # run voro++ cmd
_typeattribute = 'type'
_typebased = "types"
_probecritical = 0.0

def clusters(filename,
             finfo,
             start=-1,
             end=sys.maxsize,
             every=1,
             moleculesid=[],
             molnames=[],
             ends={},
             specific=[],
             excluded=[],
             voropp="",
             visfreq=sys.maxsize,
             histograms={},
             nbhist=False,
             histograms2d=[],
             histograms3d=[],
             profiles=[]):

    # check the input file
    reader = pysimpp.readers.create(filename)
    if not reader: sys.exit(0)

    dirname = reader.dir

    natoms = reader.natoms

    ## parse and keep info data
    print('>> reading info file ...')
    info = parse_radii(finfo)
    if not "types" in info.keys():
        print("type based radii where expected")
        sys.exit(0)

    # set voro++ stuff
    voropputils.set_voropp(voropp)
    run_voropp = voropputils.run_voropp
    extruct_close_neighbors = voropputils.extruct_close_neighbors

    # get molecular data and select molecules of interest

    atom_molecule = reader.get_atom_molecule() - 1  # atom molecule array (index @zero)
    nmolecules = atom_molecule.max() + 1  # number of molecules
    molecule_atoms = defaultdict(list)  # molecule atoms array
    for i, im in enumerate(atom_molecule):
        molecule_atoms[im].append(i)

    # if residues names are provided convert it to moleculesid (needed for lammps)
    molecule_name = None
    if len(molnames) > 0:
        molecule_name = reader.get_molecule_name()
        moleculesid = np.where(np.isin(molecule_name, molnames))[0] + 1

    # selected molecules
    selected = np.sort(moleculesid) - 1  # selected molecules (index @zero)
    nselected = selected.size  # number of selected molecules
    hasselected = nselected > 0  # selected flag

    # set what to read
    attributes = 'id x y z type'
    dmap = {0: ('x', 'ix'), 1: ('y', 'iy'), 2: ('z', 'iz')}
    dmapinv = {'x': 0, 'y': 1, 'z': 2}
    ndims = 3

    print('>> count configurations in dump file(s) ...')
    reader.set_attributes(attributes)

    # set the molecule selection flag
    molecule_selected = np.empty(shape=(nmolecules), dtype=np.bool)
    molecule_selected[:] = False
    molecule_selected[selected] = True
    # set the atom selection flag
    atom_selected = np.empty(shape=(natoms), dtype=np.bool)
    atom_selected[:] = False if len(selected) > 0 else True
    for im in selected:
        atom_selected[molecule_atoms[im]] = True

    # get atomic data
    atom_type = reader.get_atom_type()
    atom_name = reader.get_atom_name()
    atom_mass = reader.get_atom_mass()
    atom_element = reader.get_atom_element()
    atom_radius = np.array(
        tuple(
            map(lambda x: info[_typebased][str(x).upper()],
                atom_type)))  # this will not work for non arithmetic atom_type

    # LDP check
    # based on the excluded atoms fix the selected state for molecules and atoms
    if not len(excluded) == 0:
        atom_excluded = np.isin(atom_type, excluded)
        atom_selected = np.where(atom_excluded, False, atom_selected)
        for im in range(nmolecules):
            if molecule_selected[im]:
                _atoms = molecule_atoms[im]
                if np.all(atom_excluded[_atoms]):
                    molecule_selected[im] = False
    else:
        atom_excluded = np.array([False] * natoms)

    # restrain the cluster analysis using the specific atoms provided
    # essentially, all system atoms except the excluded will be used
    # will be used for the voronoi analysis. The selected atoms will
    # be used to retrieve the required statistics while the specific
    # atoms will be used foe the cluster analysis. note that:
    # 1 excluded is a subset of atoms
    # 2 selected is a subset of atoms
    # 3 selected can't be a subset of excluded
    # 4 specific intersect with selected is non emtpy and contains
    #     the atoms to be used fro the cluster analysis
    atom_molname = np.array( list(map( lambda im: molecule_name[im], atom_molecule)))
    if len(specific) > 0:
        if not len(specific) == len(molnames):
            print( "A specific list of atom names should be provided for each residue to be considered.")
            sys.exit(0)
        atom_specific = np.zeros((natoms), dtype=np.bool)
        for _resname, _specific in zip(molnames, specific):
            _atom_specific = atom_molname == _resname
            if not "*" in _specific:
                _atom_specific[:] = np.logical_and(_atom_specific, np.isin(atom_name, _specific))
            atom_specific[:] = np.logical_or(atom_specific, _atom_specific)
        # atom_specific = np.isin(atom_name, specific)
        # atom_specific[:] = np.logical_and(atom_selected, atom_specific)
        if not np.any(atom_specific):
            print( "There are no residues to consider. Check the -resname and 0specific arguments.")
            sys.exit(0)
    else:
        atom_specific = atom_selected

    # density profiles calculations
    doprofiles = len(profiles) > 0
    if doprofiles:
        bprofiles = defaultdict(
            lambda: defaultdict(lambda: Binning.free(0.25, 0.0, False)))
        profileid = np.zeros(natoms, dtype=int)
        profilemap = {}
        for ik, k in enumerate( sorted( profiles.keys())):
            v = profiles[k]
            for name in v[1]:
                profileid[ np.where((atom_name == name)&(np.isin(atom_molname,v[0]))) ] = ik+1
            profilemap[ ik+1] = k
        if excluded:
            profileid[ atom_excluded ] = 0

    # check ends (add missing residues)
    # TODO finish "-ends" implementation (remove or generalize)
    res_ends = {}
    if len(ends) > 0:
        for _res in molnames:
            _ends = ends.setdefault(_res,[])
            res_ends[_res] = tuple(np.array(_ends)-1)
    else:
        res_ends['CTAC'] = (0,15)

    # coordinates and conformation stuff
    rw = np.empty(shape=(natoms, ndims),
                  dtype=np.float32)  # wraped coordinates
    # if _wholeused:
    #     pass
    # else:
    #     r = np.empty( shape=(natoms,ndims), dtype=np.float32)  # unwraped coordinates
    #ip = np.empty( shape=(natoms,ndims), dtype=np.int32)   # periodic indexes
    periodic = [True, True, True]  # system periodicity
    steps = []  # loaded steps
    boxes = []  # loaded boxed

    # histograms and stuff
    # clusters = []     # clusters spotted per configurations
    clsize = []  # length of the clusters spotted per configurations
    clresnum = []  # number of residues in the clusters
    clnumber = []
    hclsize = Histogram.free(1, 0, addref=False)  # cluster size histogram
    hclnumber = Histogram.free(1, 0,
                               addref=False)  # number of clusters histogram
    # properties histograms
    hprop = defaultdict(lambda: Histogram.free(1, 0, addref=False))
    hprop['sqk'] = Histogram.free(0.1, 0, addref=False)
    hprop['sqk_mol'] = Histogram.free(0.1, 0, addref=False)
    hprop['qlong'] = Histogram.free(0.1, 0, addref=False)
    hprop['b_sqrg'] = Histogram.free(0.1, 0, addref=False)
    # properties,size histograms
    hprop2D = defaultdict(lambda: Histogram2D(1.0, (0, 0), addref=False))
    hprop2D['sqk'] = Histogram2D((0.1, 1.0), (0, 0), addref=False)
    hprop2D['sqk_mol'] = Histogram2D((0.1, 1.0), (0, 0), addref=False)
    hprop2D['qlong'] = Histogram2D((0.1, 1.0), (0, 0), addref=False)

    # tracer = MClusterTracker()

    # prepare on the fly log
    log = ClPropertiesLog(dirname)
    if True:
        mollog = ClMolecularLog(dirname)
        orderlog = ClOrderLog(dirname)
        detailslog = ClDetailsLog(dirname)

    # number of close contact neighbors histogram. one histogram is
    # calculated per residues pair encountered during the clusters
    # trace.

    # cluster residue - residue histograms
    ######################################
    hist2d = {}
    for _pair in histograms2d:
        hist2d["%s_%s" % tuple(_pair)] = _pair
    hist3d = {}
    for _triplet in histograms3d:
        hist3d["%s_%s_%s" % tuple(_triplet)] = _triplet
    hneighbors2d = defaultdict(lambda: Histogram2D(1.0, (0, 0), addref=False))
    hneighbors3d = defaultdict(
        lambda: Histogram3D(1.0, (0, 0, 0), addref=False))
    hneighbors = Histogram.free(1.0, 0, addref=False)

    has2dTotalHistograms = not len(histograms2d) == 0
    has3dTotalHistograms = not len(histograms3d) == 0
    hasTotalHistograms = nbhist or has2dTotalHistograms or has3dTotalHistograms

    hist_neighbors = defaultdict(
        lambda: defaultdict(lambda: Histogram.free(1, 0, addref=False)))

    # atoms group - residue histograms
    ##################################
    _chk_rnm = {}
    _chk_grps = {}
    _all_res = {}
    _s = set()
    for k, v in histograms.items():
        _chk_rnm[k] = v[0]
        _chk_grps[k] = v[1]
        _all_res[k] = v[2]
        _s.update(set(v[2]))
    __all_res = tuple(_s)

    _chk_grps_indxs = {}
    for k, v in _chk_grps.items():
        _atom_ok = atom_molname == _chk_rnm[k]
        _atom_ok[:] = np.logical_and(_atom_ok, np.isin(atom_name, v))
        _ndxs = np.where(_atom_ok)[0]
        _chk_grps_indxs[k] = _ndxs.reshape(int(_ndxs.size / len(v)), len(v))

    hasHistograms = not len(histograms) == 0
    ##################################

    vis = not visfreq == sys.maxsize

    print('>> reading dump file(s) ...')
    while (True):
        step, box, data = reader.read_next_frame()
        if step is None:
            break
        elif step < start:
            continue
        elif step > end:
            break

        if not step % every == 0:
            continue

        if box:
            if _debug and len(steps) == 0:
                print("{0:4s}  {1:4s}  {2:4s}  {3:4s}  {4:4s}  {5:4s}  {6:4s} \
 {7:4s}  {8:7s}  {9:7s}  {10:7s}  {11:5s}  {12:3s}"                                                   .format('b','b1','b2','c','c1','d','d1','sqk','bb1','bb2','bb3','ee','n'))

            steps.append(step)
            boxes.append(box)
            for k, v in dmap.items():
                np.copyto(rw[:, k], data[v[0]])
            #    np.copyto(ip[:, k] ,  data[v[1]])
            #r[:,:] = fastunwrapv( rw, ip, box.va, box.vb, box.vc)

            # get voronoi cells for each atom
            fields = run_voropp(rw,
                                box,
                                atom_radius,
                                atom_excluded,
                                periodic,
                                dirname=dirname)

            # get the volume and the neighbors of the voronoi cells
            # and create connections list of connected molecule pairs (im,jm)
            neighbors = [set() for i in range(nmolecules)
                         ]  # molecular neighbors per molecule
            atom_neighbors = [None for i in range(natoms)
                              ]  # molecular neighbors per atom
            _nodes = {}
            _connections = {}
            for i, field in enumerate(fields):
                if not atom_selected[i]:
                    continue
                im = atom_molecule[i]
                _neighbors = extruct_close_neighbors(field)
                if atom_specific[i]:
                    _cc = [j for j in _neighbors if atom_specific[j]]
                    _indexes = tuple(
                        map(lambda x: box.minimum_index(x), rw[_cc] - rw[i]))
                    # _connections = [ Connection(i,j,pbc=jp) for j, jp in zip(_cc,_indexes) ]
                    _nodes[i] = Node(i)
                    _connections[i] = zip(_cc, _indexes)
                # convert to molecular neighbors
                _neighbors = tuple(
                    set(jm for jm in [atom_molecule[j] for j in _neighbors]
                        if jm != im))
                atom_neighbors[i] = _neighbors
                neighbors[im].update(_neighbors)

            for i, n in _nodes.items():
                n.connections = [
                    Connection(n, _nodes[j], pbc=jp)
                    for j, jp in _connections[i]
                ]

            # free some memory
            del _connections
            del fields

            # update atom based neighbors histograms
            ########################################
            for k, v in _chk_grps_indxs.items():
                for _indxs in v:
                    nb = set()
                    im = atom_molecule[_indxs[0]]
                    for _i in _indxs:
                        if not atom_neighbors[_i] is None:
                            nb.update(atom_neighbors[_i])
                    _c = Counter(map(lambda x: molecule_name[x], nb))
                    for res in _all_res[k]:
                        _c[res] += 0
                    for jname in _c.keys():
                        hist_neighbors[k][jname].add(_c[jname])

            # update residue based neighbors histograms
            ###########################################
            if hasTotalHistograms:
                for im, nb in enumerate(neighbors):
                    iname = molecule_name[im]
                    if molecule_selected[im]:
                        _c = Counter(map(lambda x: molecule_name[x], nb))
                        for res in __all_res:  # all residues should be present
                            _c[res] += 0
                        if nbhist:
                            for jname in _c.keys():
                                hist_neighbors[iname][jname].add(_c[jname])
                        for k, v in hist3d.items():
                            hneighbors3d[k].add((_c[v[0]], _c[v[1]], _c[v[2]]))
                        for k, v in hist2d.items():
                            hneighbors2d[k].add((_c[v[0]], _c[v[1]]))
                        hneighbors.add(len(nb))

            # create the graph
            G = nx.Graph()
            G.add_nodes_from(_nodes.values())
            # add the non periodic connections
            _connections = []
            for n in _nodes.values():
                _connections += [
                    c for c in n.connections
                    if c[0] < c[1] and not c.is_periodic()
                ]
            G.add_edges_from(_connections)
            # identify the connected subgraphs (subclusters)
            _connected = sorted(nx.connected_components(G),
                                key=len,
                                reverse=True)
            _subclusters = [
                MCluster.create(i, x) for i, x in enumerate(_connected)
            ]
            for _c in _subclusters:
                _c.update_molecules(atom_molecule)
            _clusters = MCluster.defrag(_subclusters)

            # analyse
            # tracer.update(_clusters, step*reader.timestep)
            r = rw if _wholeused else reader.get_unwrapped()
            for i, _cl in enumerate(_clusters):
                MCluster.udpate_coordinates(_cl, box, rw, r, atom_molecule,
                                           molecule_atoms)
                if vis and step % visfreq == 0:
                    # if True:
                    # MCluster.write(_cl, r, fmt='xyz',
                    #     molecule_atoms=molecule_atoms, atom_element=atom_element,
                    #     fname="cluster%d_step%d.xyz"%(i,step), dirname=dirname)
                    MCluster.write(_cl,
                                  r,
                                  fmt='gro',
                                  molecule_atoms=molecule_atoms,
                                  molecule_name=molecule_name,
                                  atom_name=atom_name,
                                  box=box,
                                  fname="cluster%d_step%d.gro" % (i, step),
                                  dirname=dirname)
                # add calculated properties as cluster's attributes
                _cl._infinit = _cl.is_infinit()
                # _sqrg: square radious of gyration
                # _b: asphericity
                # _c: acylindricity
                # _sqk: relative shape anisotropy
                MCluster.order(_cl, r, atom_mass, molecule_atoms, molecule_name, neighbors, res_ends)
                MCluster.shape(_cl, r, atom_mass, molecule_atoms)
                if doprofiles:
                    MCluster.profile(_cl, r, box, profileid, profilemap, bprofiles)
                _size = len(_cl.molecules)
                # group q matrix eigenvalues based on three basic shapes:
                # spherical: all eigenvalues are close to zero (a threshold of 0.1 will be used)
                # cylindrical:  one negative (close to -0.5) and two positives (close to 0.25)
                # planar: two negative (close to -0.5) and one positive (close to 1.0)
                # other: anything else
                if not _cl._infinit:
                    hprop['b_sqrg'].add(_cl._b / _cl._sqrg)
                    hprop['sqrg'].add(_cl._sqrg)
                    hprop['b'].add(_cl._b)
                    hprop['c'].add(_cl._c)
                    hprop['sqk'].add(_cl._sqk)
                    hprop['bbox_x'].add(_cl._bbox[0])
                    hprop['bbox_y'].add(_cl._bbox[1])
                    hprop['bbox_z'].add(_cl._bbox[2])
                    hprop2D['sqrg'].add((_cl._sqrg, _size))
                    hprop2D['b'].add((_cl._b, _size))
                    hprop2D['c'].add((_cl._c, _size))
                    hprop2D['sqk'].add((_cl._sqk, _size))
                    hprop2D['bbox_x'].add((_cl._bbox[0], _size))
                    hprop2D['bbox_y'].add((_cl._bbox[1], _size))
                    hprop2D['bbox_z'].add((_cl._bbox[2], _size))

                    hprop['sqee_mol'].add(_cl._msqee)
                    hprop['sqrg_mol'].add(_cl._msqrg)
                    hprop['b_mol'].add(_cl._mb)
                    hprop['c_mol'].add(_cl._mc)
                    hprop['sqk_mol'].add(_cl._msqk)
                    hprop['qloc_mol'].add(_cl._qloc)

                    hprop2D['sqee_mol'].add((_cl._msqee, _size))
                    hprop2D['sqrg_mol'].add((_cl._msqrg, _size))
                    hprop2D['b_mol'].add((_cl._mb, _size))
                    hprop2D['c_mol'].add((_cl._mc, _size))
                    hprop2D['sqk_mol'].add((_cl._msqk, _size))
                    hprop2D['qloc_mol'].add((_cl._qloc, _size))

                    hprop['qlong'].add(_cl._qlong)

            # if any( [_cl._infinit for _cl in _clusters ]):
            if True:
                mollog.log(step, _clusters)
                orderlog.log(step, _clusters)
                detailslog.log(step, _clusters)
            log.log(step, _clusters)

            # check if more that two molecules with the first residue name given exist in the cluster
            connected = []
            _rescount = defaultdict(int)
            _first = molnames[0]
            _clusterslen = []
            for _cl in _clusters:
                # _res = [ molecule_name[ atom_molecule[i]] for i in _cl ]
                _res = [molecule_name[im] for im in _cl.molecules]
                _clusterslen.append(len(_res))
                if _res.count(_first) > 1:
                    connected.append(_cl)
                    for _r in _res:
                        _rescount[_r] += 1
            # LDP TODO: in order to reduce memory we should keep
            # directly the length of the clusters needed for the
            # final post process
            # clusters.append( _clusters)
            clsize.append(_clusterslen)
            clresnum.append(_rescount)
            # for cl in map( len, connected):
            for cl in map(lambda x: len(x.molecules), connected):
                hclsize.add(cl)
            _clnumber = len(connected)
            hclnumber.add(_clnumber)
            clnumber.append(_clnumber)

        else:
            break

    # complete the statistics
    # nconfs = len( clsize) # number of processed configurations
    # maxsize = hclsize.variable.max()
    # cnumbers = np.zeros( shape=(nconfs,maxsize+1), dtype=int)
    # for iconf, clusterlen in enumerate( clsize):
    #     a = cnumbers[iconf]
    #     # add to the size number distribution
    #     for ic in clusterlen:
    #         a[ ic] += 1

    if len(steps) == 0:
        print("No frames processed")
        return

    # write down the results
    for ik, iv in hist_neighbors.items():
        for jk, jv in iv.items():
            jv.write(dirname + os.sep + '%s_%s_neighbors.dat' % (ik, jk))

    for k, _ in hist3d.items():
        _h = hneighbors3d[k]
        _h.normalize()
        _range = _h.bin_range()
        for i in range(_range[0, 0], _range[0, 1]):
            __h = _h.get_h2d_x(i)
            __h.write(dirname + os.sep + "%s%d_TOTAL3D_neighbors.dat" % (k, i))

    for k, _ in hist2d.items():
        _h = hneighbors2d[k]
        _h.normalize()
        _h.write_conditional_(dirname + os.sep + "%s_TOTAL2D_neighbors.dat" % k)

    if hasTotalHistograms:
        hneighbors.write(dirname + os.sep + "TOTAL_neighbors.dat")

    hclsize.write(dirname + os.sep + 'hclsize.dat')
    hclnumber.write(dirname + os.sep + 'hclnumber.dat')
    # f=open(dirname+os.sep+'cnumbers.dat','w')
    # f.write( "# number of cluster per size (column) per frame (row)\n")
    # f.write( "# size:    %s  total\n" %  " ".join( map(str, range(1,maxsize+1))))
    # for i, j in enumerate(steps):
    #     _total = sum(cnumbers[i,1:])
    #     f.write(" %8d  %s  %s\n" % ( j, " ".join( map(str, cnumbers[i,1:])), str(_total) ))
    # f.close()
    f = open(dirname + os.sep + 'cmolecules.dat', 'w')
    f.write("# number of molecules in clusters (column) per frame (row)\n")
    f.write("# step    %s\n" % "   ".join(molnames))
    for _i, _res in zip(steps, clresnum):
        _s = "    ".join(map(str, [_res[_n] for _n in molnames]))
        f.write("%-7d    %s\n" % (_i, _s))
    f.close()
    # f.close()
    f = open(dirname + os.sep + 'cnumber.dat', 'w')
    f.write("# number of clusters per frame\n")
    f.write("# step    clusters\n")
    for _i, _num in zip(steps, clnumber):
        f.write("%-7d    %d\n" % (_i, _num))
    f.close()

    # write down properties histograms
    for k, v in hprop.items():
        _name = "h%s.dat" % k
        v.write(dirname + os.sep + _name)
    # write down 2d properties histograms
    for k, v in hprop2D.items():
        _name = "h%s2d.dat" % k
        v.normalize()
        v.write(dirname + os.sep + _name)

    # write down profiles
    if doprofiles:
        for k, v in bprofiles.items():
            for _k, _v in v.items():
                _name = dirname + os.sep + "%s_%s.prof"%(profilemap[_k],k)
                _v.write( _name, header="# Number density profile")

    if True:
        mollog.close()
        orderlog.close()
        detailslog.close()
    log.close()

def _short_description():
    return 'Trace the clusters formed in the system based on close contact analysis.'

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser( description=_short_description() )

    # add arguments
    string = '''
    the path to the simulation trajectory file. A topology file
    with the same based name should reside on the same directory,
    otherwise it will be provided with the -topo option. '''
    parser.add_argument('path', type=isfile, help=string)

    string = '''
    start processing form timestep n [inclusive] (step based for lammps
    and time based for gromacs). '''
    parser.add_argument('-start', nargs=1, type=int, metavar='n', \
        default=[-1], help=string)

    string = '''
    stop processing at timestep n [inclusive] (step based for lammps and
    time based for gromacs). '''
    parser.add_argument('-end', nargs=1, type=int, metavar='n', \
        default=[sys.maxsize], help=string)

    string = '''
    processing frequency (every dn timesteps, step based for lammps
    and time based for gromacs). '''
    parser.add_argument('-every', nargs=1, type=int, metavar='dn', \
        default=[1], help=string)

    string = '''
    visualization frequency (save gro files for the traced clusters every
    dn timesteps iteration based for lammps and time based for gromacs). '''
    parser.add_argument('-vis', nargs=1, type=int, metavar='dn', \
        default=[sys.maxsize], help=string)

    group = parser.add_mutually_exclusive_group()
    chktype = IsList("wrong molecules indexs range (check: %s)",itemtype=int,positive=True)
    string = '''
    molecules to be used. A comma seperated list with the ranges of molecules
    ids e.g. "1,2,3" or "1:10,20,30:100" '''
    group.add_argument('-molecules', nargs=1, type=chktype, default=[[]], \
        metavar='molid range', help=string)

    chktype = IsList("wrong residue names range (check: %s)",itemtype=str,positive=True)
    string = '''
    molecule types/names to be used. A comma seperated list with the names of
    the molecules e.g. "Na,Cl" '''
    group.add_argument('-molnames', nargs=1, type=chktype, default=[[]], \
        metavar='molnames range', help=string)

    chktype = IsListOfNamedList("wrong ends argument (check: %s)", itemtype=int,
        positive=True, llen=2)
    string = '''
    provide the two end atoms defining the end to end vectors for the molecular
    types participating in a cluster (see "-molnames" option). '''
    parser.add_argument('-ends', nargs=1, type=chktype, default=[{}], \
        metavar='list end atoms', help=string)

    chktype = IsList("wrong atom names range (check: %s)",itemtype=str,positive=True)
    string='''
    the atom types to be exluded system wide from the close contact analysis. A
    comma seperated list with the type names should be provided e.g. "HA,HW". '''
    parser.add_argument('-excluded', nargs=1, type=chktype, default=[[]], \
        metavar='types range', help=string)

    chktype = IsListOfList("wrong specific format (check: %s)")
    string = '''
    the names of the atoms to be considered in the close contact analysis. A comma
    separated list for each molecular type in the "-molnames" argument should be
    provided. If a wildcard "*" is given for a residue then all the atoms of the
    residue will be considered. For example, if the molecules A,B are given in
    "-molnames" argument then this argument should look like "C1,C2,C3:O1" meaning
    that atoms C1,C2,C3 should be considered for residue A and atom O1 for residue B. '''
    parser.add_argument('-specific', nargs=1, type=chktype, default=[[]], \
                    metavar='atoms for each molname', help=string)

    string = '''
    the file with the radii of the atoms. It can be element or type based.
    The first line of the file contains the keywords:
    [element|type] [r|d].
    The first, is the atom type identifier and the second if the radius (r)
    or the diameter (d) is given for each type. The rest of the lines contain
    pairs of values (type, rarious). The type could be either a number (type id)
    or a string (type name). '''
    parser.add_argument('-radii', nargs=1, type=argparse.FileType('r'),
        metavar='file with atoms\' type/elemet radii', help=string)

    string = '''
    the histograms of the number of neighbor pairs to be calculated. Each pair
    consists of a groups of atoms and a list of molecules. The distance is the
    distance between their centers of masses. A list of pairs separated with "@"
    should be provided. Each pair contains the information for a histogram
    in the following format:
        GROUPNAME:MOLNAME:ATOMSLIST:RESIDUESLIST
    with:
        GROUPNAME:  the name of the group
        MOLNAME:    the name of the molecule where the group belongs (all the atoms
                    should belong to the same molecule)
        ATOMSLIST:  a comma-separated list of atoms define the group. The
                    atoms should belong to the same residue.
        MOLECULSLIST: a comma-separated list with molecule types. One histogram
                      will be calculated for each molecular type in the list
    The histogram is written in the file GROUPNAME_RESNAME_neighbors.dat in the
    simulation directory.

    For example, the line:
        C1:CTAC:C1,H31,H32,H33:CTAC,CL,SOL@C2:CTAC:C2,H1,H2:CTAC,CL,SOL
    defines two groups and the following pairs: (C1,H31,H32,H33)-CTAC,
    (C1 H31 H32 H33)-CL, (C1 H31 H32 H33)-SOL, (C2,H1,H2)-CTAC,
     (C2,H1,H2)-CL, (C2,H1,H2)-SOL. The following files are wirtten in
     the simulation directory: C1_CTAC_neighbors.dat, C1_CL_neighbors.dat,
     C1_SOL_neighbors.dat, C2_CTAC_neighbors.dat, C2_CL_neighbors.dat,
     C2_SOL_neighbors.dat. '''
    def argshist(string):
        ''' check the -hist option argument. '''
        lines = string.split("@")
        ret = {}
        for line in lines:
            tk = line.split(":")
            if not len(tk) == 4:
                msg = "wrong pair group format (check: %s)" % line
                raise argparse.ArgumentTypeError(msg)
            rnm = tk[1]
            grp = list(map(lambda x: x.strip(), tk[2].split(",")))
            res = list(map(lambda x: x.strip(), tk[3].split(",")))
            if len(grp) == 0 or len(res) == 0:
                msg = "wrong pair group format (check: %s)" % line
                raise argparse.ArgumentTypeError(msg)
            if tk[0] in ret:
                msg = "group %s is previously defined (check: %s)" % (tk[0],
                                                                      line)
                raise argparse.ArgumentTypeError(msg)
            ret[tk[0]] = (rnm, grp, res)
        return ret

    parser.add_argument('-hist', nargs=1, type=argshist, default=[{}],
        metavar='list of neighbor group pairs', help=string)

    string = '''
    for the molecular types participating in a cluster, calculate the histograms of
    their total number of neighbors with the other molecular types in the system. For
    example, consider the system consisting of molecules with types A, B and C. If the
    clusters consist of types A and B, the script will calculate the histograms of the
    total number of neighbors for the pairs A-A, A-B, A-C, B-A, B-B and B-C. For each
    pair, the histogram is written in file A_B_neighbors.dat in the simulation directory. '''
    parser.add_argument('--nbhist', dest='nbhist', default=False, action='store_true', help=string)

    chktype = IsListOfList("wrong argument (check: %s)", itemtype=str, llen=2)
    string = '''
    for a cluster, calculates the conditional probability of having n neighbor
    residues of type A given that m neighbor residues of type B exist. The
    argument is a list of pairs separated with ":" e.g. A,B:C,D. For each pair,
    the conditional probability is written in the file A_B_TOTAL2D_neighbors.dat
    in the simulation directory.
    '''
    parser.add_argument('-hist2d', nargs=1, type=chktype, default=[[]],
        metavar='list of pairs of residues', help=string)

    chktype = IsListOfList("wrong argument (check: %s)", itemtype=str, llen=3)
    string = '''
    for a cluster, calculate the conditional probability of having n neighbor residues
    of type C given that m neighbor residues of type B exist for all the possible values
    of the number of neighbor residues of type A. The argument is a list of triplets
    separated with ":" e.g. A,B,C:D,E,F. For each triplet, the conditional probability
    is written in the file An_B_C_TOTAL3D_neighbors.dat in the simulation directory.
    '''
    parser.add_argument('-hist3d', nargs=1, type=chktype, default=[[]],
        metavar='list of triplets of residues', help=string)

    chktype = IsListOfNamedList("wrong profiles argument (check: %s)", klen=3, itemtype=str)
    string = '''
    calculate the number density profile of the given groups of atoms. The distance
    considered for the profiles depends on the position of the atoms in the micelle
    and eventually from the shape of the cluster. If the atom is located in a spherical
    micelle or in the spherical caps of an elongated/waged micelle, the distance from
    the center of the sphere is used. If the atoms belong to a cylindrical column or
    in the body of an elongated micelle, the length of its projection on the column
    axis is taken. The argument is a set of lists of comma-separated atoms name lists 
    separated with "@". For example the argument "HEAD:CTAC:C17,N,C18,C19@WAT:SOL:OW"
    defines two named lists (groups); HEAD cosist of atoms C17,N,C18, and C19 belong
    to CTAC molecules and WAT consists of OW atoms belong to SOL molecules. Atoms
    specified with the "-excluded" argument will be excluded also here. For each list,
    the density profile will be written in the file profile_{listname}_{shape}.prof in
    the simulation directory where list name is the name of the group of atom and shape
    the shape of the cluster. Three types of clusters are considered:  spherical (s),
    cylindrical infinite periodic (c), and wedged/elongated (e). Therefore file HEAD_s.prof
    corresponds to the density profiles of atoms C17,N,C18, and C19 belong to CTAC molecules,
    with the respect to the center of mass for spherical clusters.
    '''
    parser.add_argument('-profiles', nargs=1, type=chktype, default=[[]],
        metavar='list of atom names', help=string)

    parser.add_argument('-voropp', nargs=1, type=str, default=[""], metavar='voro++ executable', \
        help='provide a voro++ executable to be used instead of pyvoro python module.')

    # parse the arguments
    args = parser.parse_args()

    clusters(args.path,
             args.radii[0],
             start=args.start[0],
             end=args.end[0],
             every=args.every[0],
             moleculesid=args.molecules[0],
             molnames=args.molnames[0],
             ends=args.ends[0],
             specific=args.specific[0],
             excluded=args.excluded[0],
             voropp=args.voropp[0],
             visfreq=args.vis[0],
             histograms=args.hist[0],
             nbhist=args.nbhist,
             histograms2d=args.hist2d[0],
             histograms3d=args.hist3d[0],
             profiles=args.profiles[0])

if __name__ == '__main__':
    command()
