#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from collections import defaultdict, Counter
from pysimpp.cluster import mcluster
import inspect

import numpy as np
import networkx as nx

from pysimpp.utils.simulationbox import SimulationBox
from pysimpp.utils.statisticsutils import Variable, Histogram
from pysimpp.utils.clusterutils import Node, Connection, Cluster, Assembly, EvolutionTracker
from pysimpp.utils.vectorutils import get_length, get_unit, get_projection, projection, get_vertical
from pysimpp.fastpost import fastcom, gyration, inertia, order_parameter, order_parameter_local # pylint: disable=no-name-in-module

_debug = False
_critical_size = 4
__iu = np.array((1.0,0.0,0.0))
__ju = np.array((0.0,1.0,0.0))
__ku = np.array((0.0,0.0,1.0))

class UnitHexagolal2D():
    ''' implements a 2d hexagonal grid. '''
    __sqrt23 = np.sqrt(2./3.)
    __onesf = np.ones(3, dtype=float)
    __iu = np.array((1.0,0.0,0.0))
    __ju = np.array((0.0,1.0,0.0))
    __ku = np.array((0.0,0.0,1.0))

    def __init__(self):
        ''' initialize unit cell. '''
        pass

    @classmethod
    def create(cls, axis, point, box, atol = 0.1):

        # access static members
        _ones = UnitHexagolal2D.__onesf
        _sqrt23 = UnitHexagolal2D.__sqrt23

        # calculate axis/box utility vectors
        sign = np.sign( axis)
        absaxis = np.abs( axis)
        v = np.array( ( box.a, box.b, box.c))
        u = get_unit( v)
        obj = None

        iscubic = np.isclose(box.a, box.b) and np.isclose(box.a, box.c)
        isortho = np.isclose( box.a*box.b*box.c, box.volume)

        if not isortho: return None

        # space diagonal
        if np.all( np.isclose(absaxis, u, atol = atol)) and iscubic:
            normal = u * sign # plane normal vector (space diagonal)

            # sign indicates the box apex where the space diagonal points
            # if the origin is located at the center of the box.
            # Get the box apex in box coordinates (i,j,k) with i,j,k := {0,1}
            # (with the origin being located at the start of the box)
            apex = np.array( (sign+_ones)/2, dtype=int)
            origin = v * apex
            _i, _j, _k = apex
            # find the coordinates of the neighboring/connected apexes
            iv = v * np.array( ((_i+1)%2, _j, _k))
            jv = v * np.array( (_i, (_j+1)%2, _k))
            # kv = v * np.array( (_i, _j, (_k+1)%2))
            # get their projections to the plance
            ip = get_vertical( normal, iv-origin)
            jp = get_vertical( normal, jv-origin)
            # force right handed (ip,jp,normal) frame
            if np.dot( np.cross(iv,jv), normal) < 0: ip, jp = jp, ip
            kp = ip + jp
            # box edges length
            a, b, c = v
            # the total length of the column equals box space diagonal
            length = np.sqrt(a*a+b*b+c*c)
            obj = cls()
            obj.pbc = obj._pbc_space_diagonal

        # face diagonal
        if obj is None:
            for mask in (np.array((1,1,0)),np.array((1,0,1)),np.array((0,1,1))):
                _v = v * mask
                cmask = (mask+1)%2 # complementary mask (for the normal to the plane)
                u = get_unit( _v)
                if np.all( np.isclose(absaxis, u, atol = atol)):
                    normal = u * sign # plane normal vector (face diagonal)

                    # find box apex where the diagonal points in box coordinates
                    # [ (i,j,k) with i,j,k := {0,1}]
                    sign[np.array(cmask,dtype=bool)] = -1 # force using the lower box planes
                    apex = np.array( (sign+_ones)/2, dtype=int)
                    origin = _v * apex
                    # find neighboring apexes
                    _i, _j, _k = apex
                    iv = v * np.array( ((_i+1)%2, _j, _k)) # 1st apex connected to plane point
                    jv = v * np.array( (_i, (_j+1)%2, _k)) # 2nd apex connected to plane point
                    kv = v * np.array( (_i, _j, (_k+1)%2)) # 3rd apex connected to plane point

                    _box = np.array( (iv, jv, kv))
                    iv, jv = _box[np.array( mask,dtype=bool)]
                    kv = v * cmask
                    ip = get_vertical( normal, iv - jv) * 0.5
                    # force (ip,kv,normal) to be right handed
                    if np.dot( np.cross(ip, kv), normal) < 0.0:
                        ip = get_vertical( normal, jv - iv) * 0.5
                    jp = kv

                    a, b = get_length(iv-origin), get_length(jv-origin)
                    c = get_length( kv)
                    diagonal = np.sqrt(a*a+b*b)
                    # not necessary but set a convinient length for the normal vector
                    normal *= diagonal
                    # the total length of the column equals box space diagonal
                    length = diagonal

                    kp = origin = np.zeros(3)
                    obj = cls()
                    obj.pbc = obj._pbc_face

        # face normal
        if obj is None:
            for mask in (np.array((1,0,0)),np.array((0,1,0)),np.array((0,0,1))):
                _v = v * mask
                cmask = (mask+1)%2 # complementary mask (for the normal to the plane)
                u = get_unit( _v)
                if np.all( np.isclose(absaxis, u, atol = atol)):
                    normal = u * sign # plane normal vector (face diagonal)

                    _box = np.array((box.va, box.vb, box.vc))
                    ip, jp = _box[ np.array( cmask, dtype=bool)]
                    if np.dot( np.cross( ip, jp), normal) < 0:
                        ip, jp = jp, ip
                    kp = _box[ np.array( mask, dtype=bool)]
                    a, b, c = get_length(ip), get_length(jp), get_length(kp)

                    # not necessary but set a convinient length for the normal vector
                    normal[:] *= c
                    # the total length of the column equals box face diagonal
                    length = c

                    kp = origin = np.zeros(3)
                    obj = cls()
                    obj.pbc = obj._pbc_face

        if not obj is None:
            # find transofmation matrix from fram to local coordinates (assuming row vectors)
            matrix = np.array( (ip, jp, normal))
            tolocal = np.linalg.inv( matrix)
            obj.axis = np.array(axis)
            obj.point = np.array(point)
            obj.normal = np.array(normal)
            obj.origin = np.array(origin)
            obj.u = np.array(ip)
            obj.v = np.array(jp)
            obj.w = np.array(kp)
            obj.tolocal = np.array( tolocal)
            obj.a = float(a)
            obj.b = float(b)
            obj.c = float(c)
            obj.length = float(length)
            obj.box = np.array( (get_length(ip),get_length(jp)))
            return obj

    def _pbc_space_diagonal(self, r):
        d = np.zeros(r.size//3)
        # wrap fractional coordinates into the cell
        r[:,0:2] -= np.floor(r[:,0:2])
        # convert the projections on the (v,u) plane to lab frame
        _r = np.outer(r[:,0],self.v)+np.outer(r[:,1],self.u)
        # get the projections on the v, u, and w vectors
        pv = np.dot(_r, self.v) / np.dot(self.v, self.v)
        pu = np.dot(_r, self.u) / np.dot(self.u, self.u)
        pw = np.dot(_r, self.w) / np.dot(self.w, self.w)
        # find the distances of the projections from the near apex
        # of the unit cell
        condv = pv <= 0.5
        condu = pu <= 0.5
        condw = pw <= 0.5
        where = condv & condu & condw
        remain = np.array(where)
        _x = _r[where]
        d[where] = np.sqrt( np.sum( _x*_x, axis=1))
        where = (~ condv) & (pu <= 0.0)
        remain |= where
        _x = _r[where]-self.v
        d[where] = np.sqrt( np.sum( _x*_x, axis=1))
        where = (~ condu) & (pv <= 0.0)
        remain |= where
        _x = _r[where]-self.u
        d[where] = np.sqrt( np.sum( _x*_x, axis=1))
        remain = ~ remain
        _x = _r[remain]-self.v-self.u
        d[remain] = np.sqrt( np.sum( _x*_x, axis=1))
        return d

    def _pbc_face(self, r):
        _box = self.box
        _r = ( r[:,:2] - np.rint( r[:,:2])) * _box
        d = np.sqrt( np.sum(_r*_r,axis=1))
        return d

    def distance(self, r):
        # center coordinates to the reference point and convert
        # them to fractional coordinates of the unit cell
        _r = np.dot(r - self.point, self.tolocal)
        return self.pbc( _r)


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
        self._shape = ''

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
            if len(_molecules) > _critical_size: _connected.append(_c)

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
    def writendx( fndx, icluster, cluster, molecule_atoms):
        _every = 20
        fndx.write("[ ASSMBL%d ]\n"%icluster)
        iatom = 0
        for im in cluster.molecules:
            for iat in molecule_atoms[im]:
                iatom+=1
                s = "\n" if iatom % _every == 0 else " "
                fndx.write("%d%s"%(iat,s))
        if not iatom % _every == 0: fndx.write("\n")
        fndx.write("\n")
        if _debug:
            _head=[16,17,18,19,53,54,55,56,57,58,59,60,61]
            # _head=[15,16,17,18,19,47,48,53,54,55,56,57,58,59,60,61]
            _tail=list(set(range(62))-set(_head))
            fndx.write("[ ASSMBL%d_HEAD ]\n"%icluster)
            iatom = 0
            for im in cluster.molecules:
                for iat in np.array(molecule_atoms[im])[_head]:
                    iatom+=1
                    s = "\n" if iatom % _every == 0 else " "
                    fndx.write("%d%s"%(iat,s))
            if not iatom % _every == 0: fndx.write("\n")
            fndx.write("\n")
            fndx.write("[ ASSMBL%d_TAIL ]\n"%icluster)
            iatom = 0
            for im in cluster.molecules:
                for iat in np.array(molecule_atoms[im])[_tail]:
                    iatom+=1
                    s = "\n" if iatom % _every == 0 else " "
                    fndx.write("%d%s"%(iat,s))
            if not iatom % _every == 0: fndx.write("\n")
            fndx.write("\n")

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
                _natoms = sum( [len(molecule_atoms[im]) for im in cluster.molecules])
                f.write("%d\n\n" % _natoms)
                for im in cluster.molecules:
                    for iat in molecule_atoms[im]:
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
                _natoms = sum( [len(molecule_atoms[im]) for im in cluster.molecules])
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
    def order( cluster, r, atom_mass, molecule_atoms, molecule_name, neighbors, ends):
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
        eigvec = np.zeros( (len(cluster.molecules),9), dtype=np.float64)
        exclude = np.zeros( len(cluster.molecules), dtype=np.bool_)
        _exclude = np.zeros( len(cluster.molecules), dtype=np.bool_)
        _molecule_name = molecule_name[ list(cluster.molecules)]

        # groups molecules per species
        _matoms = defaultdict(list)
        for im in cluster.molecules:
            _matoms[ molecule_name[im]].append( molecule_atoms[im])
        cluster.nspecies = { _k:len(_v) for _k, _v in _matoms.items() }

        exclude[:] = False
        _exclude[:] = True

        ### LDP specific for CTAC and the connectivity used
        if len(ends) >0 :
            cluster._msqee = {}
            cluster._msqrg = {}
            cluster._mb = {}
            cluster._mc = {}
            cluster._msqk = {}            
            for _mname, _ends in ends.items():
                _end0, _end1 = ends[ _mname]
                _start = [ _mat[_end0] for _mat in _matoms[_mname] ]
                _end = [ _mat[_end1] for _mat in _matoms[_mname] ]
                _ee = r[_end] - r[_start]
                # no pbc needed for coordinates since they are already whole
                _msqee = (_ee*_ee).sum(axis=1).mean()
                cluster._msqee[_mname] = _msqee
                setattr(cluster,'sqee_%s'%_mname,_msqee)
                #######################
                _mask = _molecule_name == _mname
                exclude[ _mask] = False
                _exclude[:] = True
                _exclude[ _mask] = False
                _rp, _rg, _eigval, _eigvec, _ierr =  gyration(_r, _masses, _molecules, _exclude)
                eigvec[ _mask,:] =  _eigvec[ _mask,:]               

                _sqrg = _eigval.sum(axis=1)
                sqrg = _sqrg.mean() # square radious of gyration
                _b = _eigval[:,0] - 0.5*(_eigval[:,1]+_eigval[:,2])
                b = _b.mean()       # aspherisity
                _c = _eigval[:,1]-_eigval[:,2]
                c = _c.mean()       # acylindricity
                _sqk = (_b*_b+0.75*_c*_c)/(_sqrg*_sqrg)
                sqk = _sqk.mean()   # anisotropy
                cluster._msqrg[_mname]  = sqrg
                cluster._mb[_mname]  = b
                cluster._mc[_mname]  = c
                cluster._msqk[_mname]  = sqk
                setattr(cluster,'b_%s'%_mname,b)
                setattr(cluster,'c_%s'%_mname,c)
                setattr(cluster,'sqk_%s'%_mname,sqk)
                setattr(cluster,'sqrg_%s'%_mname,sqrg)
        else:
            _rp, _rg, _eigval, eigvec, _ierr =  gyration(_r, _masses, _molecules, exclude)

        # molecular axes array
        _v = eigvec[:,0:3]
        # cluster order tensor
        _q, _eigval, _eigvec, _ierr = order_parameter( _v, exclude)
        cluster.qval = _eigval  # order parameters
        cluster.qvec = _eigvec  # directors
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
        _qlocal, _qlochist, _ierr = order_parameter_local( _v, _nneighbors, _neighbors, _bin, _nbins )
        cluster.qlocal = _qlocal
        cluster._qlochist = _qlochist

    @staticmethod
    def shape( cluster, r, atom_mass, molecule_atoms):
        ''' Calculate cluster shape and order characteristics. The calculated
            properties are added as cluster attributes. '''
        # cluster atoms based on molecucles
        _catoms = []
        for im in cluster.molecules:
            _catoms += molecule_atoms[im]
        cluster._catoms = _catoms
        # shape of the cluster
        _r = r[ _catoms]
        _masses = atom_mass[ _catoms]
        _nomasses = np.ones( _masses.size, dtype=np.float32)
        _molecules = np.zeros(len(_catoms),dtype=np.int32)
        _exclude = np.array((False))
        _rp, _rg, _eigval, _eigvec, _ierr =  gyration(_r, _nomasses, _molecules, _exclude)
        _com = fastcom( _r, _masses, _molecules, 1)
        _e = _eigval[0]
        sqrg = _e.sum()              # square radious of gyration
        b = _e[0]-0.5*(_e[1]+_e[2])    # aspherisity
        c = _e[1]-_e[2]                # acylindricity
        sqk = (b*b+0.75*c*c)/(sqrg*sqrg)# anisotropy
        # "normalize" asphericity dividing by Rg^2
        b = b / sqrg 
        # b = b / _e[0]  # dividing by max eigenvalue
        # or following gromacs implementation
        # b = (3./2.) * (np.sum((_e - sqrg/3.)**2))/sqrg**2
        # "normalize" acylindricity dividing by _e[1]
        c = c / sqrg
        # c = c / _e[1] #  dividing by max eigenvalue
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
        cluster._bounds = on

        # estimate the total linear number molecular density
        # find the real molecule for each atom of the cluster
        # _eps = np.sqrt( cluster._msqee['CTAC'])
        _eps = 20.0
        _epsxmin = on[0][0]+_eps
        _epsxmax = on[0][1]-_eps
        if len(cluster.molecules) > 50 and (_epsxmax-_epsxmin)>0.5*_eps:
            _rmolecules = []
            for i, im in enumerate(cluster.molecules):
                _atoms = molecule_atoms[im]
                _catoms += _atoms
                _rmolecules += [i]*len(_atoms)
            _cm = fastcom( _rp, _masses, _rmolecules, len(cluster.molecules))
            _x = _cm[:,0]
            _xL = _eps # 2.0
            _dL = 10.0
            _xmin = on[0][0]+_xL
            _xmax = on[0][1]-_xL-_dL
            _delta = _xmax - _xmin
            rand = np.random.rand
            _v = Variable()
            for i in range(100000):
                _x0 = _xmin + _delta*rand()
                _v.set( float(len(_x[ (_x>=_x0) & (_x<=_x0+_dL)]))/_dL)
            _natoms = len(_x[ (_x>=_epsxmin) & (_x<=_epsxmax)])
            cluster._ldensity = (_v.mean(),_v.std(),_natoms,_xmax-_xmin)
        else:
            cluster._ldensity = (0.0,0.0,0.0,0.0)

        # gives the same results as in the case of gyration tensor
        # (ie. different eigenvalues but the same eigenvectors)
        # _inert, _eigval, _eigvec, _ierr =  inertia(_r, _masses, _molecules, _exclude)
        cluster.sqrg = sqrg
        cluster.b = b
        cluster.c = c
        cluster.sqk = sqk
        cluster.bbox = bbox
        cluster.rgval = _eigval[0]
        cluster.rgvec = _eigvec[0]
        cluster.com = _com[0]

        # retrieve the q matrix for the cluster and check if any
        # of the primary axes coinsides with a director
        if 'qvec' in vars(cluster):
            dirs = cluster.qvec.reshape(3,3)
            _crit = 10.0*np.pi/180.0
            # take the eigenvector corresponding to the larger eigenvalues
            # (longest clusters extension axis)
            _v = _eigvec[0][0:3]
            _dir = np.where(np.arccos( dirs.dot( _v)) < _crit)[0]
            cluster.qlong = cluster.qval[_dir[0]] if not _dir.size == 0 else 0.0
        else:
            cluster.qlong = _eigvec[0][0:3]

    @staticmethod
    def chkaxis(cluster, box, atol=0.11):
        hbox = UnitHexagolal2D.create( cluster.rgvec[:3], cluster.com, box, atol=atol)
        return hbox

    @staticmethod
    def profile(cluster, r, box, profileid, profilemap, bprofiles, b2dprofiles, hprofiles, h2dprofiles, vprofiles):
        # TODO support triclinic boxes (work in fractional coordinates)
        # (currently on cubic boxes are supported)

        _profilebin = mcluster._profilebin
        # TODO add criteria/cases based on the molecular weight
        # cylindrical column
        if cluster._infinit:
            cluster._shape = 'c'
            _hbox = MCluster.chkaxis(cluster, box)
            _bs = bprofiles['c']
            _hs = defaultdict(lambda: Histogram.free(_profilebin, 0.0, False))

            _h1ds = hprofiles['c']
            # _h2ds = h2dprofiles['c']
            # _b2ds = b2dprofiles['c']
            _vs = vprofiles['c']['total']
            _size = len(cluster.molecules)
            # _sizebin = int(_size // 10.0)
            # _b2d = _b2ds[_sizebin]
            _vs.set(_size)

            _dist = _hbox.distance(r)
            for k in profilemap:
                if k == 0: continue
                _h = _hs[k]
                _h1d = _h1ds[k]
                # _h2d = _h2ds[k]
                _d = _dist[np.where(profileid == k)]
                if len(_d):
                    np.vectorize(_h.add)(_d) # TODO compile the universal function once
                    np.vectorize(_h1d.add)(_d)
                    # for _dx in _d: _h2d.add((_size,_dx))
                else:
                    print( 'WARNING: something is wrong with profiles calculation in %s(@%d)' % (inspect.currentframe().f_code.co_name, inspect.currentframe().f_lineno))
            for k, v in _hs.items():
                _bs[k].add_histogram(v, options={'type':'p','profile':'c', 'length':_hbox.length})
                # _b2d[k].add_histogram(v, options={'type':'p','profile':'c', 'length':_hbox.length})

        # spherical
        elif cluster.b < 0.2 and cluster.c < 0.3 and cluster.sqk < 0.2 and len(cluster.molecules) > 40:
            cluster._shape = 's'
            _bs = bprofiles['s']
            _hs = defaultdict(lambda: Histogram.free(_profilebin, 0.0, False))

            _h1ds = hprofiles['s']
            _h2ds = h2dprofiles['s']
            _b2ds = b2dprofiles['s']
            _size = len(cluster.molecules)
            _sizebin = int(_size // 10.0)
            _b2d = _b2ds[_sizebin]
            vprofiles['s']['total'].set(_size)
            vprofiles['s'][_sizebin].set(_size)

            _cm = cluster.com
            _dr = r - _cm
            box.set_to_minimum(_dr)
            _dist = np.sqrt( np.sum( _dr*_dr, axis=1))
            for k in profilemap:
                if k == 0: continue
                _h = _hs[k]
                _h1d = _h1ds[k]
                _h2d = _h2ds[k]
                _d = _dist[np.where(profileid == k)]
                if len(_d) > 0:
                    np.vectorize(_h.add)(_d) # TODO compile the universal function once
                    np.vectorize(_h1d.add)(_d)
                    for _dx in _d: _h2d.add((_size,_dx))
                else:
                    print( 'WARNING: something is wrong with profiles calculation in %s(@%d)' % (inspect.currentframe().f_code.co_name, inspect.currentframe().f_lineno))
            for k, v in _hs.items():
                _bs[k].add_histogram(v, options={'type':'p','profile':'s'})
                _b2d[k].add_histogram(v, options={'type':'p','profile':'s'})

        # rod-like
        elif 0.2 < cluster.b < 0.65 and cluster.c < 0.3 and 0.1 < cluster.sqk < 0.4:
            _rest = (cluster.bbox[1]+cluster.bbox[2])/4.0 # estimate end-to-end
            if cluster.bbox[0] > 2.0*_rest and len(cluster.molecules) > 100:
                cluster._shape = 'ec'
                _cm = cluster.com
                _axis = cluster.rgvec[:3]
                _r = r - _cm
                box.set_to_minimum(_r)
                # projection length on cylinder axis.
                _rp = np.dot(_r, _axis) / np.dot(_axis, _axis)
                _rp_cluster = _rp[cluster._catoms]
                _min, _max = np.min(_rp_cluster), np.max(_rp_cluster)
                _l0 = _min+_rest
                _cp0 = _l0 * _axis
                _l1 = _max-_rest
                _cp1 = _l1 * _axis
                _where = (_rp<_l0) | (_rp>_l1)
                _x = _r[ _where]
                _dr0 = _x - _cp0
                box.set_to_minimum(_dr0)
                _dist0 = np.sqrt( np.sum( _dr0*_dr0, axis=1))
                _dr1 = _x - _cp1
                box.set_to_minimum(_dr1)
                _dist1 = np.sqrt( np.sum( _dr1*_dr1, axis=1))
                _minimum = _dist0 < _dist1 # find the min(_dist0, _dist1)
                _dist = np.array( _dist1)
                _dist[_minimum] = _dist0[_minimum]

                _bs = bprofiles['es']
                _hs = defaultdict(lambda: Histogram.free(_profilebin, 0.0, False))

                _b2ds = b2dprofiles['es']
                _h1ds = hprofiles['es']
                _h2ds = h2dprofiles['es']
                _size = len(cluster.molecules)
                _sizebin = int(_size // 10.0)
                _b2d = _b2ds[_sizebin]
                vprofiles['es']['total'].set(_size)
                vprofiles['es'][_sizebin].set(_size)

                _profileid = profileid[ _where]
                for k in profilemap:
                    if k == 0: continue
                    _h = _hs[k]
                    _h1d = _h1ds[k]
                    _h2d = _h2ds[k]
                    _d = _dist[np.where(_profileid == k)]
                    if len(_d) > 0:
                        np.vectorize(_h.add)(_d) # TODO compile the universal function once
                        np.vectorize(_h1d.add)(_d)
                        for _dx in _d: _h2d.add((_size,_dx))
                    else:
                        print( 'WARNING: something is wrong with profiles calculation in %s(@%d)' % (inspect.currentframe().f_code.co_name, inspect.currentframe().f_lineno))

                for k, v in _hs.items():
                    _bs[k].add_histogram(v, options={'type':'p','profile':'s'})
                    _b2d[k].add_histogram(v, options={'type':'p','profile':'s'})

                _bs = bprofiles['ec']
                _hs = defaultdict(lambda: Histogram.free(_profilebin, 0.0, False))

                _b2ds = b2dprofiles['ec']
                _h1ds = hprofiles['ec']
                _h2ds = h2dprofiles['ec']
                _size = len(cluster.molecules)
                _sizebin = int(_size // 10.0)
                _b2d = _b2ds[_sizebin]
                vprofiles['ec']['total'].set(_size)
                vprofiles['ec'][_sizebin].set(_size)

                _where = ~ _where
                _profileid = profileid[ _where]
                _rpv = np.outer(_rp[_where], _axis)
                _dr = _r[_where] - _rpv
                box.set_to_minimum(_dr)
                _dist = np.sqrt( np.sum( _dr*_dr, axis=1))
                for k in profilemap:
                    if k == 0: continue
                    _h = _hs[k]
                    _h1d = _h1ds[k]
                    _h2d = _h2ds[k]
                    _d = _dist[ np.where(_profileid == k)]
                    if len(_d) > 0:
                        np.vectorize(_h.add)(_d) # TODO compile the universal function once
                        np.vectorize(_h1d.add)(_d)
                        for _dx in _d: _h2d.add((_size,_dx))
                    else:
                        print( 'WARNING: something is wrong with profiles calculation in %s(@%d)' % (inspect.currentframe().f_code.co_name, inspect.currentframe().f_lineno))

                for k, v in _hs.items():
                    _bs[k].add_histogram(v, options={'type':'p','profile':'c', 'length':_l1-_l0})
                    _b2d[k].add_histogram(v, options={'type':'p','profile':'c', 'length':_l1-_l0})

    def assembly_from_molecules(self, molecule_name, step):
        ''' Return an assembly consist of the molecules of the cluster. The kind
            for each molecules is defined from its name given in molecule_name
            argument. '''
        assembly = None
        if not self.molecules is None:
            _nodes = list( map( lambda i: Node(i,kind=molecule_name[i]), self.molecules))
            assembly = Assembly( _nodes, -1)
            _data = {}
            _data['n'] = len(_nodes)
            _data['sqrg'] = self.sqrg
            _data['b'] = self.b
            _data['c'] = self.c
            _data['sqk'] = self.sqk
            _data["diff_previous"] = 0.0
            _data["diff_cumulative"] = 0.0
            _data["diff_initial"] = 0.0
            assembly.data[step] = _data
        return assembly

