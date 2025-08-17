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

from pysimpp.utils.utils import chk_filename, parse_radii, isrange, isfile, IsList, IsListOfList, IsListOfNamedList
from pysimpp.utils.statisticsutils import Binning, Histogram, Histogram2D, Histogram3D, Variable
from pysimpp.fastpost import fastunwrapv, order_parameter # pylint: disable=no-name-in-module
from pysimpp.utils.clusterutils import Node, Connection, Assembly, EvolutionTracker, Reaction
import pysimpp.utils.voropputils as voropputils
from .mclusterutils import MCluster # pylint: disable=import-error
from .mclusterlog import ClPropertiesLog, ClMolecularLog, ClOrderLog, ClDetailsLog # pylint: disable=import-error

# global flags
_debug = False
_typebased = "types"
_profilebin = 0.25

def clusters(filename,
             fradii,
             start=-1,
             end=sys.maxsize,
             every=1,
             moleculesid=[],
             molnames=[],
             ends={},
             whole=False,
             phist={},
             phist2d=[],
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

    ## parse and keep radii data
    print('>> reading info file ...')
    radii = parse_radii(fradii)
    if not "types" in radii.keys():
        print("Problem parsiong the raddii file %s. A type based radii where expected" % fradii.name)
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
    reader.set_attributes(attributes)

    # set the molecule selection flag
    molecule_selected = np.empty(shape=(nmolecules), dtype=np.bool_)
    molecule_selected[:] = False
    molecule_selected[selected] = True
    # set the atom selection flag
    atom_selected = np.empty(shape=(natoms), dtype=np.bool_)
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
            map(lambda x: radii[_typebased][str(x).upper()],
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
        atom_specific = np.zeros((natoms), dtype=np.bool_)
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
            lambda: defaultdict(lambda: Binning.free(_profilebin, 0.0, False)))
        b2dprofiles = defaultdict(
            lambda: defaultdict(lambda: defaultdict( lambda: Binning.free(_profilebin, 0.0, False))))
        hprofiles = defaultdict(
            lambda: defaultdict(lambda: Histogram.free(_profilebin, 0.0, False)))
        h2dprofiles = defaultdict(
            lambda: defaultdict(lambda: Histogram2D((10.0, _profilebin), (0, 0), addref=False)))
        vprofiles = defaultdict( lambda: defaultdict( lambda: Variable() ))
        profileid = np.zeros(natoms, dtype=int)
        profilemap = {}
        for ik, k in enumerate( sorted( profiles.keys())):
            v = profiles[k]
            for name in v[1]:
                profileid[ np.where((atom_name == name)&(np.isin(atom_molname,v[0]))) ] = ik+1
            profilemap[ ik+1] = k
        if excluded:
            profileid[ atom_excluded ] = 0

    if True:
        tracer = EvolutionTracker()

    # check ends (add missing residues)
    mol_ends = {}
    hasends = False
    if len(ends) > 0:
        hasends = True
        for _res, _ends in ends.items():
            if _res in molnames:
                mol_ends[_res] = tuple(np.array(_ends)-1)

    # coordinates and conformation stuff
    rw = np.empty(shape=(natoms, 3),
                  dtype=np.float32)  # wraped coordinates
    reader.set_whole(whole)
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
    __properties = {
        'cluster':{ 'b':(),'c':(),'sqk':(),'sqrg':(),'bbox':(0,1,2)},
        'order':{ 'qlong':(),'qlocal':()},
        'species': { "%s_%s"%(_p,_ends):() for _ends in mol_ends for _p in ('b','c','sqk','sqee','sqrg')}
    }
    _properties_known = [ _k for _t in ('cluster','order','species') for _k in __properties[_t] ]
    __default_bins = {
        'cluster':{ 'b':0.05, 'c':0.05, 'sqk':0.05, 'sqee':1.0, 'sqrg':1.0, 'bbox':0.5 },
        'order':{ 'qlong':0.05, 'qlocal':0.05},
        'species':{  "%s_%s"%(_p[0],_ends):_p[1] for _ends in mol_ends for _p in (('b',0.05),('c',0.05),('sqk',0.05),('sqee',1.0),('sqrg',1.0))}
    }
    _cases = { 'all': ('cluster', 'order', 'species'), 'clusters': ('cluster', 'order'), 'species': ('order', 'species') }
    _given = [ _p.lower() for _p in phist]
    _special = [ _k for _k in ('all','clusters','species') if _k in _given]
    if len(_special) > 0:
        _bin = { _k:_v for _t in _cases[ _special[0]] for _k, _v in __default_bins[_t].items() }
        _properties = { _k:_v for _t in _cases[ _special[0]] for _k, _v in __properties[_t].items() }
    else:
        _bin = defaultdict( lambda: 0.1, { _k:_v[0] for _k, _v in phist.items() if _k in _properties_known })
    hprop = defaultdict(lambda: Histogram.free(1.0, 0, addref=False))
    for _k, _v in _bin.items():
        hprop[_k] = Histogram.free(_v, 0, addref=False)
    # hprop['qlocal_s'] = Histogram.free(0.05, 0, addref=False)
    # hprop['qlocal_ec'] = Histogram.free(0.05, 0, addref=False)
    # hprop['ld_s'] = Histogram.free(0.05, 0, addref=False)
    # hprop['ld_ec'] = Histogram.free(0.05, 0, addref=False)

    # size conditional properties histograms
    hprop2D = defaultdict(lambda: Histogram2D(1.0, (0, 0), addref=False))
    __default_bins_all = { _k:_v for _d in __default_bins.values() for _k, _v in _d.items() }
    for _k in phist2d:
        if _k in _properties_known:
            _v = _bin[_k] if _k in _bin else __default_bins_all[_k]
            hprop2D[_k] = Histogram2D((_v, 1.0), (0, 0), addref=False)

    # prepare on the fly log
    log = ClPropertiesLog(dirname)
    if True:
        if hasends: mollog = ClMolecularLog(dirname, mol_ends.keys())
        if hasends: orderlog = ClOrderLog(dirname)
        detailslog = ClDetailsLog(dirname, hasends)

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

    vis = not visfreq == sys.maxsize and not visfreq < 1
    if vis:
        _dir = dirname + os.sep + "clusters"
        if not os.path.exists(_dir): os.mkdir(_dir)    

    print('>> reading dump file(s) ...')
    if _debug:
        fndx = open(dirname + os.sep + "assemblies.ndx",'w')
    iframe = -1
    while (True):
        step, box, data = reader.read_next_frame()
        iframe+=1

        if step is None:
            break
        elif step < start:
            continue
        elif not iframe % every == 0:
            continue
        elif step > end:
            break

        steps.append(step)
        boxes.append(box)
        
        for k, v in {'x':0,'y':1,'z':2}.items():
            np.copyto(rw[:, v], data[k])

        # get voronoi cells for each atom
        fields = run_voropp(rw,
                            box,
                            atom_radius,
                            atom_excluded,
                            periodic,
                            dirname=dirname)

        # get the volume and the neighbors of the voronoi cells
        # and create connections list of connected molecule pairs (im,jm)
        neighbors = [set() for i in range(nmolecules) ]  # molecular neighbors per molecule
        atom_neighbors = [None for i in range(natoms) ]  # molecular neighbors per atom
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

        for i, _cl in enumerate(_clusters):
            MCluster.udpate_coordinates(_cl, box, rw, rw, atom_molecule,
                                        molecule_atoms)

            if _debug:
                if iframe == 0:
                    MCluster.writendx( fndx, i, _cl, molecule_atoms)

            # if vis and step % visfreq == 0:
            #     # if True:
            #     # MCluster.write(_cl, rw, fmt='xyz',
            #     #     molecule_atoms=molecule_atoms, atom_element=atom_element,
            #     #     fname="cluster%d_step%d.xyz"%(i,step), dirname=dirname)
            #     MCluster.write(_cl,
            #                   rw,
            #                   fmt='gro',
            #                   molecule_atoms=molecule_atoms,
            #                   molecule_name=molecule_name,
            #                   atom_name=atom_name,
            #                   box=box,
            #                   fname="cluster%d_step%d.gro" % (i, step),
            #                   dirname=dirname+os.sep+"clusters")
            # add calculated properties as cluster's attributes
            _cl._infinit = _cl.is_infinit()
            # _sqrg: square radious of gyration
            # _b: asphericity
            # _c: acylindricity
            # _sqk: relative shape anisotropy
            if hasends: MCluster.order( _cl, rw, atom_mass, molecule_atoms, molecule_name, neighbors, mol_ends)
            MCluster.shape( _cl, rw, atom_mass, molecule_atoms)
            if doprofiles:
                MCluster.profile(_cl, rw, box, profileid, profilemap, bprofiles, b2dprofiles, hprofiles, h2dprofiles, vprofiles)
            _size = len(_cl.molecules)
            # if not _cl._infinit:
            #     if _cl._shape == 's':
            #         hprop['sqrg_s'].add(_cl.sqrg)
            #         for _mname in mol_ends:
            #             hprop['sqee_%s_s' % _mname].add( _cl._msqee[ _mname])
            #         hprop['qlocal_s'].add(_cl.qlocal)
            #         hprop['ld_s'].add(_cl._ldensity[0])
            #     elif _cl._shape == 'ec':
            #         hprop['sqrg_ec'].add(_cl.sqrg)
            #         for _mname in mol_ends:
            #             hprop['sqee_%s_ec' % _mname].add( _cl._msqee[ _mname])
            #         hprop['qlocal_ec'].add(_cl.qlocal)
            #         hprop['ld_ec'].add(_cl._ldensity[0])
            #         hprop['ldn_ec'].add(_cl._ldensity[2])
            #         hprop['ldl_ec'].add(_cl._ldensity[3])

            for _p in _bin:
                hprop[_p].add( vars(_cl)[_p])
            for _p in phist2d:
                hprop2D[ _p].add( ( vars(_cl)[_p], _size))

        # analyse
        if True:
            _assemplies = [ _c.assembly_from_molecules(molecule_name, step) for _c in _clusters]
            tracer.update_current(_assemplies, step)
        #########
            for i, (_as, _cl) in enumerate(zip(_assemplies,_clusters)):
                if vis and step % visfreq == 0:
                    # if True:
                    # MCluster.write(_cl, rw, fmt='xyz',
                    #     molecule_atoms=molecule_atoms, atom_element=atom_element,
                    #     fname="cluster%d_step%d.xyz"%(i,step), dirname=dirname)
                    _i = _as.uid if hasattr(_as,'uid') else -i
                    MCluster.write(_cl,
                                rw,
                                fmt='gro',
                                molecule_atoms=molecule_atoms,
                                molecule_name=molecule_name,
                                atom_name=atom_name,
                                box=box,
                                fname="cluster%d_step%d.gro" % (_i, step),
                                dirname=dirname+os.sep+"clusters")

        # if any( [_cl._infinit for _cl in _clusters ]):
        if True:
            if hasends: mollog.log(step, _clusters, mol_ends)
            if hasends: orderlog.log(step, _clusters)
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

    # complete the statistics
    # nconfs = len( clsize) # number of processed configurations
    # maxsize = hclsize.variable.max()
    # cnumbers = np.zeros( shape=(nconfs,maxsize+1), dtype=int)
    # for iconf, clusterlen in enumerate( clsize):
    #     a = cnumbers[iconf]
    #     # add to the size number distribution
    #     for ic in clusterlen:
    #         a[ ic] += 1

    if _debug:
        fndx.close()

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
    f.write("# number of molecules in the clusters, per frame (row) and per species (column) \n")
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
        _dir = dirname + os.sep + "profiles"
        if not os.path.exists(_dir): os.mkdir(_dir)

        for k, v in bprofiles.items():
            _var = vprofiles[k]['total']
            for _k, _v in v.items():
                _name = _dir + os.sep + "%s_%s.bprof"%(profilemap[_k],k)
                _nsamples = int(_v.variable._n)
                _v.write( _name, header="# Number density profile of %s for clusters of molecular size %d(%d), from %d samples" %(k,_var.mean(),_var.std(),_nsamples))

        for k, v in b2dprofiles.items():
            for _s, _h in v.items():
                _var = vprofiles[k][_s]
                for _k, _v in _h.items():
                    _size="%d-%d_"%(_s*10,(_s+1)*10)
                    _name = _dir + os.sep + "%s_%s_%s.bprof"%(profilemap[_k],k,_size)
                    _nsamples = int(_v.variable._n)
                    _v.write( _name, header="# Number density profile of %s for clusters of molecular size %d(%d), from %d samples" %(k,_var.mean(),_var.std(),_nsamples))

        for k, v in hprofiles.items():
            _var = vprofiles[k]['total']
            for _k, _v in v.items():
                _name = _dir + os.sep + "%s_%s.hprof"%(profilemap[_k],k)
                _v.write( _name, header="# Number density profile of %s for clusters of molecular size %d(%d)" %(k,_var.mean(),_var.std()))
        for k, v in h2dprofiles.items():
            _var = vprofiles[k]['total']
            for _k, _v in v.items():
                _name = _dir + os.sep + "%s_%s.h2dprof"%(profilemap[_k],k)
                _v.normalize()
                _header = "# Number density profile of %s for clusters of molecular size %d(%d)" %(k,_var.mean(),_var.std())
                _v.write( _name, header=_header )

    if True:
        _dir = dirname + os.sep + "evolution"
        if not os.path.exists(_dir): os.mkdir(_dir)
        tracer.write(_dir)

    if True:
        if hasends: mollog.close()
        if hasends: orderlog.close()
        detailslog.close()
    log.close()

def _short_description():
    return 'Trace the clusters formed in the system based on close contact analysis.'

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser( description=_short_description())

    # add arguments
    string = '''
    the path to the simulation trajectory file. A topology file
    with the same based name should reside on the same directory,
    otherwise it will be provided with the -topo option. '''
    parser.add_argument('path', type=isfile, help=string)

    string = '''
    start and processing form time step STEP (inclusive). '''
    parser.add_argument('-start', nargs=1, type=int, metavar='STEP', \
        default=[-1], help=string)

    string = '''
    stop processing at time step STEP (inclusive). '''
    parser.add_argument('-end', nargs=1, type=int, metavar='STEP', \
        default=[sys.maxsize], help=string)

    string = '''
    process every EVERY frames (process frequency).'''
    parser.add_argument('-every', nargs=1, type=int, metavar='EVERY', \
        default=[1], help=string)

    string = '''
    visualization frequency (save gro files for the traced clusters every
    dn timesteps iteration based for lammps and time based for gromacs). '''
    parser.add_argument('-vis', nargs=1, type=int, metavar='n', \
        default=[sys.maxsize], help=string)

    string = '''
    reconstruct the molecules making them whole again, before spatial reconstruction 
    of the clusters. Use this option if the coordinates of the input trajectory are 
    wrapped in to the simulation cell and you want to correctly visualize the clusters. '''
    parser.add_argument('--whole', dest='whole', default=False, action='store_true', \
                       help=string)

    group = parser.add_mutually_exclusive_group(required=True)
    chktype = IsList("wrong molecules indexs range (check: %s)",itemtype=int,positive=True)
    string = '''
    indexes of clusters' constituent molecules. A comma separated list with the ranges of molecules
    ids e.g. "1,2,3" or "1:10,20,30:100" '''
    group.add_argument('-molecules', nargs=1, type=chktype, default=[[]], \
        metavar='<molecules\' index range>', help=string)

    chktype = IsList("wrong residue names range (check: %s)",itemtype=str,positive=True)
    string = '''
    types/names of clusters' constituent species. A comma seperated list with the names of
    the molecules e.g. "Na,Cl" '''
    group.add_argument('-molnames', nargs=1, type=chktype, default=[[]], \
        metavar='<species\' list>', help=string)

    chktype = IsListOfNamedList("wrong ends argument (check: %s)", itemtype=int,
        positive=True, llen=2)
    string = '''
    the pairs of atoms defining the end-to-end vectors of the molecular species
    participating in a cluster (see "-molnames"). For example, if the end-to-end
    vectors of clusters' constituent molecular types TIC and TOC are defined by
    atoms (1, 16) and (3,25), respectively, the arguments could be "TIC:1,16@TOC:3,25".
     '''
    parser.add_argument('-ends', nargs=1, type=chktype, default=[{}], \
        metavar='<list of end atoms per species>', help=string)

    chktype = IsList("wrong atom names range (check: %s)",itemtype=str,positive=True)
    string='''
    the atom types to be excluded system wide from the close contact analysis. A
    comma separated list with the atoms' type name should be provided e.g. "HA,HW". '''
    parser.add_argument('-excluded', nargs=1, type=chktype, default=[[]], \
        metavar='<types\' range>', help=string)

    chktype = IsListOfList("wrong specific format (check: %s)")
    string = '''
    the names of the atoms to be considered in the close contact analysis. A comma
    separated list for each molecular species in the `-molnames` argument should be 
    provided. If a wildcard "*" is given for a residue then all the atoms of the 
    molecular species will be considered. For example, if A, and B are the clusters' 
    constituent molecular types, the argument could look like "*:C1,C2,C3" specifying
    that all the atoms of species A and only the atoms C1,C2,C3 of species B should 
    be considered in the analysis. '''
    parser.add_argument('-specific', nargs=1, type=chktype, default=[[]], \
                    metavar='<list of atoms per species>', help=string)

    string = '''
    the file with the radii of the atoms. It can be element or type based. 
    The first line of the file contains the keywords "(elements|types) (r|d)"; 
    the former, specifies the atom type identifier and the latter if the 
    radius (r) or the diameter (d) is given for each type. The rest of the 
    lines contain the (type, radius) pairs. The type could be either a 
    number (type id) or a string (type name). '''
    parser.add_argument('-radii', nargs=1, type=argparse.FileType('r'), required=True,
        metavar='<file with atoms\' type/element radii>', help=string)

    string = '''
    the histograms of the number of neighbor pairs to be calculated. Each pair
    consists of a groups of atoms and a list of species. A list of pairs separated 
    with "@" should be provided. Each pair contains the information for a histogram
    in the following format:
        GROUPNAME:MOLNAME:ATOMSLIST:RESIDUESLIST
    with:
        GROUPNAME:  the name of the group
        MOLNAME:    the name of the molecule where the group belongs (all the atoms
                    should belong to the same molecule)
        ATOMSLIST:  a comma-separated list of atoms define the group. The
                    atoms should belong to the same residue.
        SPECIESLIST: a comma-separated list with molecule types. One histogram 
                      will be calculated for each molecular type in the list
    The histogram is written in the file GROUPNAME_SPECIES_neighbors.dat in the
    simulation directory.

    For example, the argument:
        C1:CTAC:C1,H31,H32,H33:CTAC,CL,SOL@C2:CTAC:C2,H1,H2:CTAC,CL,SOL
    defines two groups and the allowing pairs (C1,H31,H32,H33)-CTAC,
    (C1 H31 H32 H33)-CL, (C1 H31 H32 H33)-SOL, (C2,H1,H2)-CTAC,
    (C2,H1,H2)-CL, and (C2,H1,H2)-SOL; the files C1_CTAC_neighbors.dat, 
    C1_CL_neighbors.dat, C1_SOL_neighbors.dat, C2_CTAC_neighbors.dat, 
    C2_CL_neighbors.dat, and C2_SOL_neighbors.dat will be written in
    the simulation directory, respectively. '''
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
        metavar='<list of neighbor pairs\' groups>', help=string)

    string = '''
    for the molecular types participating in a cluster, calculate the histograms of
    their total number of neighbors with the other molecular types in the system. For
    example, consider a system consisting of molecular types A, B, and C. If the 
    constituent species of the clusters are A and B, the script will calculate the 
    histograms of the total number of neighbors pairs for A-A, A-B, A-C, B-A, B-B 
    and B-C pairs. For each pair, the histogram will be written in A_B_neighbors.dat
    file in the simulation directory. '''
    parser.add_argument('--nbhist', dest='nbhist', default=False, action='store_true', help=string)

    chktype = IsListOfList("wrong argument (check: %s)", itemtype=str, llen=2)
    string = '''
    for a cluster, calculates the conditional probability of having n neighbor 
    species of type A-A given that m neighbor species of type B-B exist. The 
    argument is a column separated list of pairs, e.g., A,B:C,D. For each pair, 
    the conditional probability is written in the file A_B_TOTAL2D_neighbors.dat 
    in the simulation directory.
    '''
    parser.add_argument('-hist2d', nargs=1, type=chktype, default=[[]],
        metavar='<list of species\' pairs>', help=string)

    chktype = IsListOfList("wrong argument (check: %s)", itemtype=str, llen=3)
    string = '''
    for a cluster, calculate the conditional probability of having n species residues
    of type C given that m species residues of type B exist for all the possible values
    of the number of neighbor species of type A. The argument is a column separated list
    of triplet, e.g., A,B,C:D,E,F. For each triplet, the conditional probability is 
    written in the file An_B_C_TOTAL3D_neighbors.dat in the simulation directory.
    '''
    parser.add_argument('-hist3d', nargs=1, type=chktype, default=[[]],
        metavar='<list of species\' triplets>', help=string)

    chktype = IsListOfNamedList("wrong profiles argument (check: %s)", klen=3, itemtype=str)
    string = '''
    calculate the number density profile of the given groups of atoms. The distance
    considered for the profiles depends on the position of the atoms in the micelle
    and eventually from the shape of the cluster. If the atom is located in a spherical
    micelle or in the spherical caps of an elongated/rod-like micelle, the distance from
    the center of the sphere is used. If the atoms belong to a cylindrical column or
    in the body of an elongated micelle, the length of its projection on the column
    axis is taken. The argument is a set of lists of comma-separated atoms name lists
    separated with "@". For example the argument "HEAD:CTAC:C17,N,C18,C19@WAT:SOL:OW"
    defines two named lists (groups); HEAD consists of atoms C17,N,C18, and C19 that belong
    to CTAC molecules and WAT consists of OW atoms belong to SOL molecules. Atoms
    specified with the "-excluded" argument will be excluded also here. For each list,
    the density profile will be written in the file profile_{listname}_{shape}.prof in
    the simulation directory where list name is the name of the group of atoms and shape
    the shape of the cluster. Three types of clusters are considered:  spherical (s),
    cylindrical infinite periodic (c), and wedged/elongated (e). Therefore file HEAD_s.prof
    corresponds to the density profiles of atoms C17,N,C18, and C19 belong to CTAC molecules,
    with  respect to the center of mass for spherical clusters.
    '''
    parser.add_argument('-profiles', nargs=1, type=chktype, default=[[]],
        metavar='<list of groups>', help=string)

    chktype = IsListOfNamedList("wrong phist argument (check: %s)", itemtype=float,
        positive=True, llen=1)
    string = '''
    provide the properties' histograms to be calculated together with the length of their bin.
    The following properties are supported:
      1) b : asphericity for both clusters and their constituent molecular species,
      2) c : acylindricity for both clusters and their constituent molecular species,
      3) sqk : anisotropy for both clusters and their constituent molecular species,
      4) sqrg : square radius of gyration for both clusters and their constituent molecular species,
      5) sqee : end-to-end distance for clusters' constituent molecular species,
      6) bbox : bounding box,
      7) qlong : global order, and
      8) qlocal : local order.
    The property of a molecular species is named by appending "_{SPECIES SNAME}" at the property's name.
    For example, with argument "c:0.01@c_CTAC:0.05", the bin length of the
    acylindricity distribution of the clusters is set to 0.01 and for the CTAC molecules to 0.05.
    The special keywords 'all', 'clusters', and 'species' can be given instead of a property's name.
    In this case, all the available properties, the properties of the traced clusters, or the 
    properties of the molecular species will be calculated, respectively, and the corresponding 
    histograms will be calculated.
    '''
    parser.add_argument('-phist', nargs=1, type=chktype, default=[{}], \
        metavar='<list of porperties histograms>', help=string)

    chktype = IsList("wrong phist2d argument (check: %s)",itemtype=str,positive=True)
    string = '''
    provide the properties for which the conditional probability of having
    a specific value, given the size of the cluster, will be calculated.
    '''
    parser.add_argument('-phist2d', nargs=1, type=chktype, default=[[]], \
        metavar='<list of properties>', help=string)

    parser.add_argument('-voropp', nargs=1, type=str, default=[""], metavar='<voro++ executable>', \
        help='provide the voro++ executable to be used instead of pyvoro python module.')

    # parse the arguments
    args = parser.parse_args()

    clusters(args.path,
             args.radii[0],
             start=args.start[0],
             end=args.end[0],
             every=args.every[0],
             whole=args.whole,
             moleculesid=args.molecules[0],
             molnames=args.molnames[0],
             ends=args.ends[0],
             phist=args.phist[0],
             phist2d=args.phist2d[0],
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
