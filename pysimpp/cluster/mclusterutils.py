#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from collections import defaultdict, Counter

import numpy as np
import networkx as nx

from pysimpp.utils.simulationbox import SimulationBox
from pysimpp.utils.statisticsutils import Variable
from pysimpp.utils.clusterutils import Node, Connection, Cluster
from pysimpp.fastpost import fastcom, gyration, inertia, order_parameter, order_parameter_local # pylint: disable=no-name-in-module

_debug = False

class MCluster(Cluster):
    ''' Implements a molecular cluster (i.e. an ensemble
        of connected molecular fragments.

        Attributes:
            molecules (set): a set with the global index of the  molecules
                in the cluster.
            molecules_first (dist): a dict with the first atom of each molecule
                in the cluster.
            molecules_shift: a dict with the shift (implemented as np.array) of
                each molecule in the cluster.
        '''

    def __init__( self, i, s):
        ''' Initialize a molecular cluster with the give set of atom indexes. '''
        super(MCluster, self).__init__( i, s)
        self.molecules = None
        self.molecules_shift = None
        self.molecules_first = None

    def update_molecules( self, atom_molecule):
        ''' Update the list with the global indexes of the molecules in the
            cluster. '''
        molecules_atoms = defaultdict(list)
        for i in self:
            molecules_atoms[ atom_molecule[i]].append( i)
        self.molecules = set( molecules_atoms.keys())
        self.molecules_first = { k:v[0] for k, v in molecules_atoms.items() }

    def update_molecules_shift( self, other):
        ''' Update molecules and their shifts when join with the parent cluster.
            This method should be called for a parent cluster otherwise nothing
            happens (use this api with care). '''
        if self.parent is self:
            if self.molecules_shift is None:
                # copy other to self
                self.molecules = set(other.molecules)
                self.molecules_shift = { im:np.zeros(3,dtype=np.float32) for im in other.molecules }
                self.molecules_first = dict(other.molecules_first)
            else:
                molecules_shift = self.molecules_shift
                molecules_first = self.molecules_first
                other_molecules_first = other.molecules_first
                shift = other.shift
                # update molecule shift and first belong only to other cluster
                for im in other.molecules-self.molecules:
                    molecules_shift[im] = shift
                    molecules_first[im] = other_molecules_first[im]
                # now update molecule
                self.molecules.update( other.molecules)

    @staticmethod
    def join( connected, clusterid):
        ''' Join the given list of connected clusters to obtain their union
            (parent cluster).
            Args:
                connected: list of connected clusters.
                clusterid: the parent cluster id.
            Returns:
                MCluster: the union (parent cluster) of the provided clusters
                    list. Returns an empty molecular cluster if the given
                    clusters are not connected.
        '''

        # prepare the clusters to make their union whole
        connected_sorted = Cluster.whole(connected)

        # join the given clusters to the larger one
        _parent = MCluster( clusterid, connected_sorted[0])
        _parent.set_parent( _parent)
        _parent.update_molecules_shift( connected_sorted[0])
        # join the nodes of each subclusters into the parent cluster
        # and also update the molecules list.
        for _c in connected_sorted[1:]:
            _parent.update( _c)
            _parent.update_molecules_shift( _c)
        # update the connections
        _parent.update_connect()

        return( _parent)

    @staticmethod
    def defrag( clusters):
        ''' Connect the clusters in the given set.
            Args:
                clusters: the list of clusters to be connected
            Returns:
                list: the list of the parent clusters.
            '''
                
        # Thε given set should be complete i.e. all the nodes of the system
        # should have been assigned in one of the clusters

        # Make sure that connections at each cluster is updated.
        for c in clusters:
            c.update_connect()
        # find clusters connectivity
        ### use your own code
        # _clusters = list(clusters)
        # _connected = []
        # current = []
        # buffer = [ _clusters[0] ]
        # while len(_clusters) > 0:
        #     buffer = [ _clusters[0] ]
        #     while len(buffer) > 0:
        #         c0 = buffer.pop(0)
        #         _clusters.remove(c0)
        #         if not c0 in current:
        #             current.append(c0)
        #         for c1 in c0.connectedto:
        #             if not c1 in buffer and not c1 in current:
        #                 buffer.append(c1)
        #     _connected.append(current)
        #     current = []
        ### or use nx
        _nodes = []
        _connections = set()
        for ic in clusters:
            _nodes.append(ic.i)
            _connections.update([ (ic.i,jc.i) for jc in ic.connectedto if ic.i<jc.i ])
        G = nx.Graph()
        G.add_nodes_from( _nodes)
        G.add_edges_from( _connections)
        # the pbc conditions have been removed. filter the groups of
        # clusters based on the number of the constituent molecules.
        # essentially at least two molecules should belong to a cluster
        # in order to be considered.
        _clusters = [ list( map( lambda y: clusters[y], x)) for x in sorted( nx.connected_components(G), key=len, reverse=True) ]
        _connected = []
        for _c in _clusters:
            _molecules = set()
            for _sc in _c:
                _molecules.update(_sc.molecules)
            if len(_molecules) > 4: _connected.append(_c)

        # and now join the connected clusters to obtain their union
        # (parent clusters)
        parents = []
        clusterid = len( clusters)+1
        for i, _c in enumerate(_connected):
            clusterid += i
            parents.append( MCluster.join( _c , clusterid))
        return sorted(parents,key=len, reverse=True)

    @staticmethod
    def udpate_coordinates( cluster, box, rw, ruw, atom_molecule, molecule_atoms):
        ''' update the given unwrapped molecules using the connectivity
            of the cluster specifed in Cluster.defrag method. '''
        shift_vector = box.shift_vector
        used=[]
        for i in cluster:
            im = atom_molecule[i]
            if im in used:
                continue
            imatoms = molecule_atoms[im]
            shift=cluster.molecules_shift[im]
            ifirst = cluster.molecules_first[im]
            dr = rw[ifirst]-ruw[ifirst]
            _dr = dr + shift_vector( shift)
            ruw[imatoms,:] += _dr
            used.append(im)

    @staticmethod
    def write( cluster, r, **kwargs):
        ''' write down the cluster. keyword arguments are:
           fmt: output format
           fname: output file name
           dirname: output directory
           molecule_atoms: atoms for each molecule
        '''
        fmt = kwargs['fmt'] if 'fmt' in kwargs else ""
        fname = kwargs['fname'] if 'fname' in kwargs else ""
        dirname = kwargs['dirname'] if 'dirname' in kwargs else "."
        if fmt == "xyz":
            if len(fname) == 0: fname="cluster%d.xyz" % cluster.i
            if not "molecule_atoms" in kwargs or \
               not "atom_element" in kwargs:
                return
            molecule_atoms = kwargs["molecule_atoms"]
            atom_element = kwargs["atom_element"]
            with open(dirname+os.sep+fname,'w') as f:
                _natoms = len(cluster.molecules)*20
                f.write("%d\n\n" % _natoms)
                for im in cluster.molecules:
                    for iat in molecule_atoms[im][:20]:
                        _r = r[iat]
                        f.write("%s %f %f %f\n" % (atom_element[iat],_r[0],_r[1],_r[2]))
        if fmt == "gro":
            if len(fname) == 0: fname="cluster%d.xyz" % cluster.i
            if not "molecule_atoms" in kwargs or \
               not "molecule_name" in kwargs or \
               not "atom_name" in kwargs or \
               not "box" in kwargs:
                return
            molecule_atoms = kwargs["molecule_atoms"]
            molecule_name = kwargs["molecule_name"]
            atom_name = kwargs["atom_name"]
            box = kwargs["box"]
            atomformat = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n"
            with open(dirname+os.sep+fname,'w') as f:
                _natoms = len(cluster.molecules)*62
                f.write("cluster %d\n" % cluster.i)
                f.write("%d\n" % _natoms)
                iatom = 0
                for i, im in enumerate(cluster.molecules):
                    for iat in molecule_atoms[im]:
                        x = r[iat] / 10.0
                        atomname = atom_name[iat]
                        iatom = iatom+1 if iatom < 99999 else 1
                        imol = i+1 if i < 99999 else 1
                        resname = molecule_name[im]
                        f.write( atomformat % ( imol, resname, atomname, iatom, x[0], x[1], x[2], 0.0, 0.0, 0.0))
                # f.write( "%f %f %f\n" % (box.a/10.0,box.b/10.0,box.c/10.0))
                a = box.va / 10.0
                b = box.vb / 10.0
                c = box.vc / 10.0
                f.write( "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n" %
                     ( a[0], b[1], c[2], a[1], a[2], b[0], b[2], c[0], c[1]))

    @staticmethod
    def order( cluster, r, atom_mass, molecule_atoms, neighbors, ends):
        ''' Calculate molecular shape and order characteristics. The calculated
            properties are added as cluster attributes. '''
        _catoms = []
        _molecules = []
        for i, im in enumerate(cluster.molecules):
            _atoms = molecule_atoms[im]
            _catoms += _atoms
            _molecules += [i]*len(_atoms)
        _r = r[ _catoms]
        _masses = atom_mass[ _catoms]
        _molecules = np.array( _molecules)
        _exclude = np.zeros( len(cluster.molecules), dtype=np.bool)

        ### LDP specific for CTAC and the connectivity used
        # (TODO remove or generalize)
        _matoms = [ molecule_atoms[im] for im in cluster.molecules]
        # _catoms = [ _at for _mat in _matoms for _at in _mat ]
        _start = [ _mat[0] for _mat in _matoms ]
        _end = [ _mat[15] for _mat in _matoms ]
        _ee = r[_end] - r[_start]
        # no pbc needed for coordinates since they are already whole
        _msqee = (_ee*_ee).sum(axis=1).mean()
        cluster._msqee = _msqee
        #######################

        _rp, _rg, _eigval, _eigvec, _ierr =  gyration(_r, _masses, _molecules, _exclude)
        _sqrg = _eigval.sum(axis=1)
        sqrg = _sqrg.mean() # square radious of gyration
        _b = _eigval[:,0] - 0.5*(_eigval[:,1]+_eigval[:,2])
        b = _b.mean()       # aspherisity
        _c = _eigval[:,1]-_eigval[:,2]
        c = _c.mean()       # acylindricity
        _sqk = (_b*_b+0.75*_c*_c)/(_sqrg*_sqrg)
        sqk = _sqk.mean()   # anisotropy
        cluster._msqrg = sqrg
        cluster._mb = b
        cluster._mc = c
        cluster._msqk = sqk
        # molecular axes array
        _v = _eigvec[:,0:3]
        # cluster order tensor
        _q, _eigval, _eigvec, _ierr = order_parameter( _v, _exclude)
        cluster._qval = _eigval  # order parameters
        cluster._qvec = _eigvec  # directors
        # cluster local order (map molecules in cluster local numbering)
        _map = { im:i for i,im in enumerate(cluster.molecules) }
        _nneighbors = []
        _neighbors = []
        for im in cluster.molecules:
            ngh = [ _map[x] for x in neighbors[im] if x in _map ]
            _nneighbors.append( len(ngh))
            _neighbors += ngh
        _nneighbors = np.array( _nneighbors)
        _neighbors = np.array( _neighbors)
        _bin = 0.1
        _nbins = int(1.5/_bin)+1
        # _h = np.zeros(_nbins, dtype=np.float32)
        _qloc, _qlochist, _ierr = order_parameter_local( _v, _nneighbors, _neighbors, _bin, _nbins )
        cluster._qloc = _qloc
        cluster._qlochist = _qlochist

    @staticmethod
    def shape( cluster, r, atom_mass, molecule_atoms):
        ''' Calculate cluster shape and order characteristics. The calculated
            properties are added as cluster attributes. '''
        _catoms = []
        for im in cluster.molecules:
            _catoms += molecule_atoms[im]

        # shape of the cluster
        _r = r[ list(_catoms)]
        _masses = atom_mass[ _catoms]
        _nomasses = np.ones( _masses.size, dtype=np.float32)
        _molecules = np.zeros(len(_catoms),dtype=np.int32)
        _exclude = np.array((False))
        _rp, _rg, _eigval, _eigvec, _ierr =  gyration(_r, _nomasses, _molecules, _exclude)
        _e = _eigval[0]
        sqrg = _e[0].sum()              # square radious of gyration
        b = _e[0]-0.5*(+_e[1]+_e[2])    # aspherisity
        c = _e[1]-_e[2]                 # acylindricity
        sqk = (b*b+0.75*c*c)/(sqrg*sqrg)# anisotropy
        # bounding ortho box: find the probability distribution along
        # each principal axis. check which range is higher that the
        # value of the corresponding uniform distribution scaled by f
        # and consider the length of the range as the characteristic
        # length for this principal direction.
        f = 0.1
        bbox = []
        on = []
        for i in (0,1,2):
            d = np.max(_rp[:,i])-np.min(_rp[:,i])
            uniform = 1.0 / d
            y, x = np.histogram(_rp[:,i], density=True)
            _on = x[ np.where(y > f*uniform)]
            bbox.append(_on[-1]-_on[0])
            on.append((_on[0],_on[-1] ))

        # estimate the total linear number molecular density
        # find the real molecule for each atom of the cluster
        if len(cluster.molecules) > 50:
            _rmolecules = []
            for i, im in enumerate(cluster.molecules):
                _atoms = molecule_atoms[im]
                _catoms += _atoms
                _rmolecules += [i]*len(_atoms)
            _cm = fastcom( _rp, _masses, _rmolecules, len(cluster.molecules))
            _x = _cm[:,0]
            _xL = 2.0
            _dL = 10.0
            _xmin = on[0][0]+_xL
            _xmax = on[0][1]-_xL-_dL
            _delta = _xmax - _xmin
            rand = np.random.rand
            _v = Variable()
            for i in range(100000):
                _x0 = _xmin + _delta*rand()
                _v.set( float(len(_x[ (_x>=_x0) & (_x<=_x0+_dL)]))/_dL)
            cluster._ldensity = (_v.mean(),_v.std())
        else:
            cluster._ldensity = (0.0,0.0)

        # gives the same results as in the case of gyration tensor
        # (ie. different eigenvalues but the same eigenvectors)
        # _inert, _eigval, _eigvec, _ierr =  inertia(_r, _masses, _molecules, _exclude)
        cluster._sqrg = sqrg
        cluster._b = b
        cluster._c = c
        cluster._sqk = sqk
        cluster._bbox = bbox
        cluster._rgval = _eigval[0]
        cluster._rgvec = _eigvec[0]

        # retrieve the q matrix for the cluster and check if any
        #  of the primitive axes coinsides with a director
        dirs = cluster._qvec.reshape(3,3)
        _crit = 10.0*np.pi/180.0
        # take the eigenvector corresponding to the larger eigenvalues
        # (longest molecular axis)
        _v = _eigvec[0][0:3]
        _dir = np.where(np.arccos( dirs.dot( _v)) < _crit)[0]
        cluster._qlong = cluster._qval[_dir[0]] if not _dir.size == 0 else 0.0
