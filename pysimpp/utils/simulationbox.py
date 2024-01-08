# -*- coding: utf-8 -*-

from math import sqrt, cos, sin, radians, floor
import numpy as np
from collections import defaultdict
import networkx as nx

from .vectorutils import *


class SimulationBox():
    ''' Implements a periodic parallelepiped sumulation box. 

        Attributes:
            origin (np.array((3),float)): the origin of the box.
            ua, ub, uc (np.array((3),float)): box spaning unit
                vectors.
            a, b, c (float): box vectors length.
            alpha, beta, gamma: box angles betwen (vb,vc), (vc,va) 
                and (va,vb) vectors.
            volume (float): box volume.
            axb, bxc, cxa (np.array((3),float)): unit cross products
                of the box spaning vector as indicated by the variable
                name.
            a_rnorm, b_rnorm, c_rnorm (float): distances between the
                bc, ca and ab faces, also equal to the projections'
                lengths of va, vb and vc to bxc, cxa and axb respectively.
    '''

    def __init__(self):
        ''' Initialize a SimulationBox object. '''
        self.origin = np.zeros(3, dtype=np.float32)
        self.va = np.zeros(3, dtype=np.float32)
        self.vb = np.zeros(3, dtype=np.float32)
        self.vc = np.zeros(3, dtype=np.float32)
        self.ua = np.zeros(3, dtype=np.float32)
        self.ub = np.zeros(3, dtype=np.float32)
        self.uc = np.zeros(3, dtype=np.float32)
        self.a = 0.0
        self.b = 0.0
        self.c = 0.0
        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.volume = 0.0

        self.axb = np.zeros(3, dtype=np.float32)
        self.bxc = np.zeros(3, dtype=np.float32)
        self.cxa = np.zeros(3, dtype=np.float32)
        self.axb_ = 0.0
        self.bxc_ = 0.0
        self.cxa_ = 0.0
        self.a_hnorm = 0.0
        self.b_hnorm = 0.0
        self.c_hnorm = 0.0
        self.a_rnorm = 0.0
        self.b_rnorm = 0.0
        self.c_rnorm = 0.0
        self.f2c = np.matrix(np.zeros((3, 3), dtype=np.float32))
        self.c2f = np.matrix(np.zeros((3, 3), dtype=np.float32))

        # print "constructor I"

    def __str__(self):
        ''' Return str(self). lammps based printout of the simulation cell '''
        ret =  " %.6f %.6f  xlo xhi\n" % (self.origin[0], self.va[0])
        ret += " %.6f %.6f  ylo yhi\n" % (self.origin[1], self.vb[1])
        ret += " %.6f %.6f  zlo zhi" % (self.origin[2], self.vc[2])
        if not self.__is_ortho():
            ret += "\n %.6f %.6f %.6f  xy xz yz" % (
                self.vb[0], self.vc[0], self.vc[1])
        return ret

    def __is_ortho(self):
        ''' check if the simulation cell is orthogonal '''
        isclose = np.isclose
        hpi = 0.5 * np.pi
        return isclose(self.alpha, hpi) and isclose(self.beta, hpi) and isclose(self.gamma, hpi)

    def set_from_scalars(self, a, b, c, alpha=90.0, beta=90.0, gamma=90.0):
        ''' Set the dimensions of the unit cell using the
            lengths of the spaning vectors (a, b and c) and
            their angles (alpha, beta and gamma).
        '''
        av = [0.0, 0.0, 0.0]
        bv = [0.0, 0.0, 0.0]
        cv = [0.0, 0.0, 0.0]
        origin = [0.0, 0.0, 0.0]
        alpha_ = radians(alpha)
        beta_ = radians(beta)
        gamma_ = radians(gamma)

        # construct the vectors
        av[0] = a
        bv[0] = b * cos(gamma_)
        bv[1] = sqrt(b * b-bv[0] * bv[0])
        cv[0] = c * cos(beta_)
        cv[1] = (b * c * cos(alpha_)-bv[0] * cv[0]) / bv[1]
        cv[2] = sqrt(c * c-cv[0] * cv[0]-cv[1] * cv[1])

        self.set_from_vectors(origin, av, bv, cv)

    def set_from_vectors(self, origin, a, b, c):
        ''' Set the object data using the input arrays/lists. origin is
            the box origin and a, b and c the box spaning vectors.
        '''
        self.origin[:] = origin
        self.va[:] = a
        self.vb[:] = b
        self.vc[:] = c
        self._update()

    def set_from_lammps(self, xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz):
        ''' Set the box data using lammps data. in this case, the
            spaning vectors are:
                va = (xhi-xlo, 0.0,     0.0)
                vb = (xy,      yhi-ylo, 0.0)
                vc = (xz,      yz,      zhi-zlo)
            and the origin is assumed zero.
        '''
        (self.xlo, self.xhi, self.ylo, self.yhi, self.zlo, self.zhi, self.xy, self.xz, self.yz) = \
            (xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz)
        # print "set_from_vectors"
        zlo_ = zlo
        zhi_ = zhi
        ylo_ = max(ylo, ylo-yz)
        yhi_ = min(yhi, yhi-yz)
        xlo_ = max(xlo, max(max(xlo-xy, xlo-xz), xlo-xy-xz))
        xhi_ = min(xhi, min(min(xhi-xy, xhi-xz), xhi-xy-xz))
        o_ = [xlo_, ylo_, zlo_]
        a_ = [xhi_ - xlo_, 0.0, 0.0]
        b_ = [xy, yhi_ - ylo_, 0.0]
        c_ = [xz, yz, zhi_ - zlo_]
        self.set_from_vectors(o_, a_, b_, c_)

    def _update(self):
        ''' Calculate object's data for self.origin, self.va, self.vb, self.vc. '''
        self._find_units()
        self._find_lengths()
        self._find_angles()
        self._find_volume()
        self._set_transformation_matrix()
        self._find_normals()

    def _find_lengths(self):
        ''' Calculate box spaning vectors lengths. '''
        self.a = get_length(self.va)
        self.b = get_length(self.vb)
        self.c = get_length(self.vc)

    def _find_units(self):
        ''' Calculate box unit spaning vectors lengths. '''
        self.ua = get_unit(self.va)
        self.ub = get_unit(self.vb)
        self.uc = get_unit(self.vc)

    def _find_angles(self):
        ''' Calculate box angles. '''
        self.alpha = get_angle(self.vb, self.vc)
        self.beta = get_angle(self.vc, self.va)
        self.gamma = get_angle(self.va, self.vb)

    def _find_volume(self):
        ''' Calculate box volume. '''
        self.volume = np.abs(np.dot(self.va, np.cross(self.vb, self.vc)))

    def _find_normals(self):
        ''' Calculate box faces surface, normal vectors and relevant variables. '''
        self.axb[:] = np.cross(self.va, self.vb)
        self.bxc[:] = np.cross(self.vb, self.vc)
        self.cxa[:] = np.cross(self.vc, self.va)
        self.axb_ = get_length(self.axb)
        self.bxc_ = get_length(self.bxc)
        self.cxa_ = get_length(self.cxa)
        self.axb = self.axb / self.axb_
        self.bxc = self.bxc / self.bxc_
        self.cxa = self.cxa / self.cxa_
        self.a_rnorm = self.bxc_ / self.volume
        self.b_rnorm = self.cxa_ / self.volume
        self.c_rnorm = self.axb_ / self.volume
        self.a_hnorm = 0.5 * self.volume / self.bxc_
        self.b_hnorm = 0.5 * self.volume / self.cxa_
        self.c_hnorm = 0.5 * self.volume / self.axb_
        # for cubic box _rnorm equals (1.0/a,1.0/b,1.0/c)
        _rnorm = (self.a_rnorm, self.b_rnorm, self.c_rnorm)
        self._rnorm = np.array(_rnorm)
        self._norm = 1.0/np.array(_rnorm)
        self._normv = np.array((self.bxc, self.cxa, self.axb))

    def _set_transformation_matrix(self):
        ''' Calculate box transformation matrixes. '''
        a = self.a
        b = self.b
        c = self.c
        cosa = cos(self.alpha)
        sina = sin(self.alpha)
        cosb = cos(self.beta)
        sinb = sin(self.beta)
        cosg = cos(self.gamma)
        sing = sin(self.gamma)

        # fractional to cartesian transformation
        m0 = a
        m1 = b * cosg
        m2 = c * cosb
        m3 = sqrt(b * b - m1 * m1)
        m4 = (b * c * cosa - m1 * m2) / m3
        m5 = sqrt(c * c - m2 * m2 - m4 * m4)
        # f->r
        # f[0] = m0*f[0]+m1*f[1]+m2*f[2] + origin[0]
        # f[1] =         m3*f[1]+m4*f[2] + origin[1]
        # f[2] =                 m5*f[2] + origin[2]
        self.f2c[0, 0] = m0
        self.f2c[1, 0] = m1
        self.f2c[2, 0] = m2
        self.f2c[1, 1] = m3
        self.f2c[2, 1] = m4
        self.f2c[2, 2] = m5

        # cartesian to fractional transformation
        v = a*b*c*sqrt(1.0-cosa*cosa-cosb*cosb-cosg*cosg+2.0*cosa*cosb*cosg)
        m0 = 1.0/a
        m1 = -cosg/(a*sing)
        m2 = b*c*(cosg*(cosa-cosb*cosg)/sing-cosb*sing)/v
        m3 = 1.0/(b*sing)
        m4 = -a*c*(cosa-cosb*cosg)/(v*sing)
        m5 = a*b*sing/v
        # # r -> f
        # f = r - origin
        # r[0] = m0*r[0] + m1*r[1] + m2*r[2]
        # r[1] =           m3*r[1] + m4*r[2]
        # r[2] =                     m5*r[2]
        # and check for precission; if the values are close to 0.0 or 1.0
        # within epsilon force them to these values
        self.c2f[0, 0] = m0
        self.c2f[1, 0] = m1
        self.c2f[2, 0] = m2
        self.c2f[1, 1] = m3
        self.c2f[2, 1] = m4
        self.c2f[2, 2] = m5

    # def setToMinimumImage(self, v):
    #     self.set_to_minimum(v)

    def set_to_minimum(self, v):
        ''' Set the input vector(s) to the corresponding minimum 
            image(s) and return the periodic indexe(s).
        '''
        # project on the normals
        prj = np.inner(v, self._normv)
        # periodic indexes
        pi = np.rint(prj * self._rnorm).astype(dtype=int)
        # non zero periodic indexes
        xxx = (pi != 0)
        if len(v.shape) > 1:
            for i, w in enumerate((self.va, self.vb, self.vc)):
                xx = xxx[:, i]
                if np.any(xx):
                    v[xx] -= pi[xx, i, np.newaxis]*w
        else:
            if xxx[0]:
                v -= pi[0]*self.va
            if xxx[1]:
                v -= pi[1]*self.vb
            if xxx[2]:
                v -= pi[2]*self.vc
        return pi

    def minimum(self, v):
        ''' Return the minimum image of the given vector 
            and its periodic images.
        '''
        u = np.array(v)
        pindx = self.set_to_minimum(u)
        return u, pindx

    def minimum_index(self, v):
        '''  Returns a tuple of periodic indexes set the
             given vector to its minimum image.  '''
        # project on the normals
        prj = np.inner(v, self._normv)
        # periodic indexes
        pi = np.rint(prj * self._rnorm).astype(dtype=int)
        return tuple(pi)
        # only for cubic boxes
        # return (int( np.rint( v[0]/self.a)),
        #         int( np.rint( v[1]/self.b)),
        #         int( np.rint( v[2]/self.c)))

    def wrap1(self, v):
        ''' Wraps the given vector(s) into the box. '''
        r = v - self.origin
        f = np.inner(r, self.c2f)
        f -= np.floor(f)
        r = np.inner(f, self.f2c)
        return self.origin+r

    def wrap(self, v):
        ''' Wraps the given vector into the box. 
            Works only for orthogonal boxes. '''
        r = v - self.origin
        for j, (rd_, d_) in enumerate(zip(self._rnorm, self._norm)):
            r[j] -= d_ * np.floor(r[j]*rd_)
        return self.origin + r

    def shift_vector(self, ip):
        ''' Return shift vector using the given periodic indexes. '''
        return ip[0] * self.va + ip[1] * self.vb + ip[2] * self.vc

    def whole(self, rw, atom_molecule, molecule_atoms, bonds):
        ''' make molecules whole (CUBIC BOX).
            rw:             wrapped coordinates
            atom_molecule:  atom molecule
            molecule_atoms: atoms per molecule
            bonds:          bond pairs
        '''
        G = nx.Graph()
        box = np.array(self.va[0], self.vb[1], self.vc[2])
        # bond vectors
        bv = rw[bonds[1], :] - rw[bonds[0], :]
        # bond periodic index
        bpi = np.round(bv / box).astype(dtype=int)
        # bond boundary cross flag
        bf = np.any(bpi != 0, axis=1)
        # bonds per molecule
        bm = defaultdict(list)
        for ib, b in enumerate(bonds):
            bm[atom_molecule[b[0]]].append((b[0], b[1], bf[ib]))

        for im, atoms in enumerate(molecule_atoms):
            G.clear()
            # add atoms local local index as nodes
            G.add_nodes_from(list(range(len(atoms))))
            # map bonds from global to local index
            local = {at: i for i, at in enumerate(atoms)}
            _bonds = tuple(
                map(lambda b: (local[b[0]], local[b[1]], b[2]), bm[im]))
            # add the non crossing bonds as edges
            edges = [(b[0], b[1]) for b in _bonds if not b[2]]
            G.add_edges_from(edges)
            # sort the subgraphs (fragments) based on their size
            fragments = sorted(nx.connected_components(G),
                               key=len, reverse=True)
            # if one fragment, the molecule is already in the box (infinit periodic or not)
            if len(fragments) > 1:
                # set atom fragment
                atfrg = {}
                for i, f in enumerate(fragments):
                    atfrg[f] = i
                # find connections and order them
                _connections = []
                for b in _bonds:
                    if b[2]:
                        f0 = atfrg[b[0]]
                        f1 = atfrg[b[1]]
                        _connections.append(
                            (f0, f1, b[0], b[1]) if f0 > f1 else (f1, f0, b[1], b[0]))
                _connections = set(_connections)
                for i, f in enumerate(fragments):
                    _connectedto = filter(lambda x: x[0] == i, _connections)

            # _connections = set(map(lambda b: sorted(b), _mb[_crossing]))


class NeighborCellList():
    ''' Implements a basic neighbor cell list. '''

    def __int__(self):
        ''' Initialize the List object. '''
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.nxy = 0
        self.n = 0
        self.dx = 0.0
        self.dy = 0.0
        self.dz = 0.0
        self.cell = None
        self.atom = None

    def _index(self, i, j, k):
        return k * self.nxy + j * self.nx + i

    def _indexes(self, index):
        k = index / self.nxy
        t = index - k * self.nxy
        j = t / self.nx
        i = t - j * self.nx
        return i, j, k

    @staticmethod
    def check(box, dx, dy, dz):
        ok = dx * 4.0 < box.a or \
            dy * 4.0 < box.b or \
            dz * 4.0 < box.c
        return ok

    def set(self, box, coordinates, dx, dy, dz):
        ''' Split the box in cells where x >= dx, y >= dy, z >= dz.
            If this is not possble return False otherwise True.
            NOTE: works for cubic cells. '''
        self.box = box
        self.coordinates = coordinates
        self.nx = max(4, int(floor(box.a / dx)))
        self.ny = max(4, int(floor(box.b / dy)))
        self.nz = max(4, int(floor(box.c / dz)))
        self.dx = box.a / self.nx
        self.dy = box.b / self.ny
        self.dz = box.c / self.nz
        self.nxy = self.nx * self.ny
        self.nxy = self.nx * self.ny
        self.n = self.nz * self.nxy
        if self.dx >= dx and self.dy >= dy and self.dz >= dz:
            natoms = len(coordinates) / 3
            self.atom = np.zeros(natoms, dtype=int)
            self.cell = []
            for i in range(self.n):
                self.cell.append([])
            o = self.box.origin
            a = self.box.a
            b = self.box.b
            c = self.box.c
            for i in range(natoms):
                j = i * 3
                v = coordinates[j:j + 3] - o
                while v[0] < 0.0:
                    v[0] += a
                while v[0] > a:
                    v[0] -= a
                while v[1] < 0.0:
                    v[1] += b
                while v[1] > b:
                    v[1] -= b
                while v[2] < 0.0:
                    v[2] += c
                while v[2] > c:
                    v[2] -= c

                i_ = int(v[0] / self.dx)
                if i_ < 0:
                    i_ = 0
                elif i_ >= self.nx:
                    i_ = self.nx-1

                j_ = int(v[1] / self.dy)
                if j_ < 0:
                    j_ = 0
                elif j_ >= self.ny:
                    j_ = self.ny-1

                k_ = int(v[2] / self.dz)
                if k_ < 0:
                    k_ = 0
                elif k_ >= self.nz:
                    k_ = self.nz-1

                index = self._index(i_, j_, k_)

                self.atom[i] = self._index(i_, j_, k_)
                self.cell[index].append(i)
            return True
        else:
            return False

    def getNeighborCells(self, index):
        ''' Returns a tuple with the indexes of the index cell neighboring cells. '''
        (x, y, z) = self._indexes(index)
        _cells = []
        _list = [-1, 0, + 1]
        for i in _list:
            # apply pbc
            _i = x + i
            if _i == -1:
                _i = self.nx - 1
            elif _i == self.nx:
                _i = 0
            for j in _list:
                # apply pbc
                _j = y + j
                if _j == -1:
                    _j = self.ny - 1
                elif _j == self.ny:
                    _j = 0
                for k in _list:
                    _k = z + k
                    if _k == -1:
                        _k = self.nz - 1
                    elif _k == self.nz:
                        _k = 0
                    _cells.append(self._index(_i, _j, _k))
        return tuple(_cells)

    def getNeighborAtoms(self, index):
        ''' Returns a tuple with the indexes of the index atoms neighboring atoms
            i.e. atoms belongs the same or neighboring cells. '''
        cindex = self.atom[index]
        cells = self.getNeighborCells(cindex)
        _list = []
        for c in cells:
            _list.extend(self.cell[c])
        return tuple(_list)

# if __name__ == '__main__':
#    #cProfile.run('reader.read_dump()','dumpprof')
#    #box = Box()
#
#    box = SimulationBox()
#    #o=np.array([0,0,0])
#    #x=np.array([1,0,0])
#    #y=[0,1,0]
#    #z=[0,0,1]
#    #box.set_from_vectors(o, x, y, z)
#    box.set_from_lammps(0.0,1.0,0.0,1.0,0.0,1.0,0.0,0.0,0.0)
#
#    # check minimum image
#    vec = np.array([0.6, 0.6, 0.0])
#    print "vec : ", vec
#    box.set_to_minimum(vec)
#    print "min(vec) : ", vec
#
#    # check volume and fractional transformation matrix
#    print "volume : ", box.volume
#    print "transdormation matrix : \n", box.m
#
#    # set some atoms
#    n=8
#    c = np.zeros((n,3),dtype=np.float32)
#    c[0,:] = (0.1, 0.1, 0.1)
#    c[1,:] = (0.3, 0.3, 0.3)
#    c[2,:] = (0.5, 0.5, 0.5)
#    c[3,:] = (0.7, 0.7, 0.7)
#    c[4,:] = (0.9, 0.7, 0.7)
#    c[5,:] = (0.1, 0.1, 0.1)
#    c[6,:] = (0.1, 0.1, 0.9)
#    c[7,:] = (0.01, 0.01, 0.01)
#    print "c (%d*3) :\n" % (n), c
#
#    c_ = c.reshape(3*n)
#    print "c (%d) :\n" % (3*n), c_
#
#    nl = NeighborCellList()
#    nl.set(box,c_,0.2,0.2,0.2)
#
#    print "neighbor cells for cell 445: \n", nl.getNeighborCells(445)
#    print "neighbor cells for cell 0: \n", nl.getNeighborCells(0)
#
#    print "\n\n"
#    for i, c in enumerate(nl.cell):
#        if not len(c) == 0:
#            print "[", i, "] ", nl._indexes(i), " : ", c
#
#    print "\n\n"
#    print " neighbors list of atoms 0:\n", nl.getNeighborAtoms(0)
#    print " neighbors list of atoms 1:\n", nl.getNeighborAtoms(1)
#    print " neighbors list of atoms 2:\n", nl.getNeighborAtoms(2)
#    print " neighbors list of atoms 3:\n", nl.getNeighborAtoms(3)
#    print " neighbors list of atoms 4:\n", nl.getNeighborAtoms(4)
#


def test():
    v = (np.random.rand(1000, 3)*2.0-1.0)*20.0
    box = SimulationBox()
    box.set_from_scalars(2., 2., 2., 90., 90., 90.)
    for v_ in v:
        u = box.wrap(v_)
        # u = box.wrap_(v_)

        # u1 = np.array(v_)
        # u2 = np.array(v_)
        # u1 = box.wrap(u1)
        # u2 = box.wrap_(u2)
        # if not np.allclose( u1, u2, atol=1.e-12):
        #     print('probem')
