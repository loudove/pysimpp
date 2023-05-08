# -*- coding: utf-8 -*-

import os
from shutil import which
import subprocess
import pyvoro

__debug = False
__pyvoro_extended = True    # use the extended pyvoro version
__voropp_external = True    # run voro++ cmd
__voropp_exe = None         # run voro++ cmd

def set_voropp( voropp=""):
    ''' iinitialize voro++ method. if voropp is provided, then voro++
        will be executed externally and the results will be loaded from
        the resulted output file. if voropp is a valid executable will
        be used directly while in a different case system path will be
        checked for a valid voro++ executable. otherwise, the pyvoro module
        will be used. if the latest extension for improved memory handling
        is pressent it will be utilized.  '''
    global __debug
    global __pyvoro_extended
    global __voropp_external
    global __voropp_exe

    # by default use pyvoro module and if available the extented implementation
    __voropp_external = False
    __pyvoro_extended = hasattr(pyvoro,'compute_3d_voronoi')
    __voropp_exe = "none"

    # if a valid voro++ executable is provided the use it. otherwise check if
    # a valid executable can be found.
    if len( voropp) > 0: # check external executable
        if os.path.isfile(voropp) and os.access(voropp, os.X_OK):
            print("voropp: run using the given voro++ executalbe %s" % voropp)
            __voropp_external = True
            __voropp_exe = voropp
        else:
            _voropp = which("voro++")
            if os.path.isfile(_voropp) and os.access(_voropp, os.X_OK):
                print("voropp: run using voro++ executalbe %s" % _voropp)
                __voropp_external = True
                __voropp_exe = _voropp
    else:
        if __debug:
            print("DEBUG: use voro++ python bindings")

def run_voropp( rw, box, atom_radius, atom_excluded, periodic, dirname="."):
    ''' Run voro++ using the extended python binding or the
        external executable to reduce memory requirements. '''
    global __debug
    global __pyvoro_extended
    global __voropp_external
    global __voropp_exe

    if __voropp_external:

        ## write the input for voro++
        vfilename = dirname + os.sep + 'voro.atoms'
        f = open(vfilename, 'w')
        _map={}
        _id = 0
        for _r, _radii in zip( rw, atom_radius):
            _x, _y, _z = _r
            _id += 1
            f.write("%d %f %f %f %f\n" % (_id, _x, _y, _z, _radii))

        f.close()

        ## run voro++
        if os.name == "posix":
            prefix = "./"
            startup = None
        else:
            prefix = ""
            # Prevent shell window from showing up
            startup = subprocess.STARTUPINFO()
            startup.dwFlags |= subprocess.STARTF_USESHOWWINDOW

        #args = [_voroexe, '-p', '-r', '-y', '-yp', '-yv', '-c',
        args = [__voropp_exe, '-p', '-r', '-c',
                '\"%i %q %r %v %F %w %P %o %m %g %s %a %f %t %l %n\"']
        #args += shlex.split("'%i %w %P %o %m %g %s %a %f %t %n'")
        zero = "0.0"
        args += [zero, str(box.a), zero, str(box.b), zero, str(box.c), vfilename]
        strargs = " ".join( map(str, args))
        if __debug: print("DEBUG: run voro++ (%s) with: %s" % ( __voropp_exe, strargs))
        try:
            proc = subprocess.Popen(args, shell=False, startupinfo=startup)
            proc.wait()
            if __debug: print('DEBUG: voro++ run ... ok')
        except:
            print("problem calculating per atomm accessible/free volume surface")
            proc.kill()
            proc.wait()
            return

        f = open(vfilename + '.vol')
        lines = f.readlines()
        lines = map(lambda x: x[1:-2], lines)
        lines = sorted(lines, key=lambda x: int(x[:x.index(" ")]))
        f.close()
        fields = lines

    else:
        ## use voro++ from python.
        # it is faster but we need to parse the new data structures
        if __pyvoro_extended:
            fields = pyvoro.compute_3d_voronoi(rw,
                [ [ 0.0, box.a], [ 0.0, box.b],[ 0.0, box.c]],
                4.0,
                radii=atom_radius,
                periodic=periodic,
                excluded=atom_excluded,
                properties=['neighbors'])
        else:
            fields = pyvoro.compute_voronoi(rw,
                [ [ 0.0, box.a], [ 0.0, box.b],[ 0.0, box.c]],
                4.0,
                radii=atom_radius,
                periodic=periodic)

    return fields

def extruct_close_neighbors( field):
    ''' extruct the atom based neighbors list from the given field. '''
    # global __pyvoro_extended
    # global __voropp_external
    # global __voropp_exe

    if __voropp_external:
        tk = field.split()
        nv = int(tk[7])                 # the number of vertices in the Voronoi cell
        nf = int(tk[8 + 2 * nv + 2])    # The number of faces of the Voronoi cell
        neighbors = [ j for j in map(lambda x: int(x)-1, tk[8 + 2 * nv + 3 + 4 * nf:8 + 2 * nv + 3 + 5 * nf]) ]
    else:
        if __pyvoro_extended:
            neighbors = field['neighbors']
        else:
            neighbors = [ y for y in [ x['adjacent_cell'] for x in field['faces' ] ] ]

    return neighbors


def debug( ):
    ''' print debug info. '''
    # global __pyvoro_extended
    # global __voropp_external
    # global __voropp_exe
    print ("voropp: use extented implementation %r" % __pyvoro_extended)
    print ("voropp: use voro++ executable %r" % __voropp_external)
    print ("voropp: voro++ executable %s " % __voropp_exe)

import numpy as np
from collections import defaultdict
from enum import Enum
# from .vectrorutils import get_unit
from pysimpp.utils.vectorutils import get_unit

class ViewMode( Enum):
    ALWAYS = 0  # display allways
    ANY = 1      # display if vertex free
    ALL = 2      # display only if all vertexes free

class VoronoiCell():
    ''' Implements a voronoi cell.'''

    def __init__(self):
        ''' Default constructor '''
        self.i = -1
        self.shared = None

    @classmethod
    def fromrecord(cls, record, indx, network):
        ''' Ininitalize a voronoi cell from the record provided from voro++
        python bindings (hacked). '''
        obj = cls()
        obj.i = indx
        obj.x = record['original']
        obj.r = record['radius'] #+ network.probecritical
        obj.rsq = obj.r * obj.r
        obj.volume = record['volume']
        obj.surface = record['surface']
        obj.vrmax = 0.25 * record['rmaxsq'] # according to voro++ code

        _vstr = ["(%.8f,%.8f,%.8f)"%tuple(x) for x in record['vertices']]
        obj.vertexes = tuple( map(network.add_vertex, _vstr))
        obj.nv = len( obj.vertexes)

        _vfs = tuple( record['areas'])
        _faces = [x['vertices'] for x in record['faces']] # faces local indexes
        _vfl = tuple( [tuple([obj.vertexes[i] for i in x]) for x in _faces])
        obj.neighbors = [x['adjacent_cell'] for x in record['faces']]
        obj.normals = tuple( [tuple(x) for x in record['normals']])

        obj.faces = []
        obj.edges = set()
        for v, s, n in zip(_vfl, _vfs, obj.normals):
            face = network.add_face( v, s)
            obj.faces.append( face)
            obj.edges = obj.edges.union( face.edges)

        obj.avolume = obj.volume
        obj.asurface = obj.surface

        return obj

    @classmethod
    def fromline(cls, line, network):
        '''Initialize a voronoi cell from the output line of voro++ and update the data
           of the given network.'''
        obj = cls()
        tk = line.split()
        # %i ( The particle ID number)
        obj.i = int(tk[0]) - 1
        # %q ( The position vector of the particle, short for “%x %y %z”, wrapped)
        obj.x = list(map(float, tk[1:4]))
        # %r ( The radius of the particle (only printed if the polydisperse information is available))
        obj.r = float(tk[4]) # + network.probecritical
        obj.rsq = obj.r * obj.r
        # %v ( The volume of the Voronoi cell)
        obj.volume = float(tk[5])
        # %F ( The total surface area of the Voronoi cell)
        obj.surface = float(tk[6])
        # %w ( The number of vertices in the Voronoi cell)
        obj.nv = int(tk[7])
        nv = obj.nv  # keep the number of vertexes
        # keep only the indexes of each vertex
        # %P ( A list of the vertices of the Voronoi cell in the format (x,y,z), relative to the global coordinate system)
        obj.vertexes = tuple( map(network.add_vertex, tk[8:8 + nv]))
        # %o ( A list of the orders of each vertex) IS NOT NEEDED TODO: remove it
        #obj.vo = tuple( map( int, tk[8+nv:8+2*nv]))
        # %m ( The maximum radius squared of a vertex position, relative to the particle center)
        obj.vrmax = float(tk[8 + 2 * nv])
        # %g ( The number of edges of the Voronoi cell)
        obj.ne = int(tk[8 + 2 * nv + 1])
        # %s ( The number of faces of the Voronoi cell)
        obj.nf = int(tk[8 + 2 * nv + 2])
        nf = obj.nf  # keep the number of faces
        # %a ( A list of the orders of the faces, showing how many edges make up each face) IS NOT NEEDED TODO: remove it
        #vfo = tuple( map( int, tk[8+2*nv+3:8+2*nv+3+nf]))
        # %f ( A list of areas of each face)
        vfs = tuple(
            map(float, tk[8 + 2 * nv + 3 + nf:8 + 2 * nv + 3 + 2 * nf]))
        # %t ( A list of bracketed sequences of vertices that make up each face index @zero)
        vfl = tuple([tuple([obj.vertexes[int(i)] for i in x[
                    1:-1].split(",")]) for x in tk[8 + 2 * nv + 3 + 2 * nf:8 + 2 * nv + 3 + 3 * nf]])
        # %l ( A list of normal vectors for each face : can be calculated instead)
        obj.normals = tuple([tuple(map(
            float, x[1:-1].split(","))) for x in tk[8 + 2 * nv + 3 + 3 * nf:8 + 2 * nv + 3 + 4 * nf]])
        obj.faces = []
        obj.edges = set()
        for v, s, n in zip(vfl, vfs, obj.normals):
            face = network.add_face(v, s)
            obj.faces.append(face)
            obj.edges = obj.edges.union(face.edges)
        # %n ( A list of the neighboring particle or wall IDs corresponding to each face)
        obj.neighbors = tuple(
            [int(x)-1 for x in tk[8 + 2 * nv + 3 + 4 * nf:8 + 2 * nv + 3 + 5 * nf]])
        # accessible volume and surface
        obj.avolume = obj.volume
        obj.asurface = obj.surface

        return obj

    @staticmethod
    def get_neighbors( record):
        return extruct_close_neighbors( record)

    def face_point(self):
        '''Returns the points on the constituent faces.'''
        points = []
        for f in self.faces:
            points.append(list(f.vertexes)[0].x)
        return np.array(points)

    def resetshared(self):
        ''' Reset xplanes attribute. '''
        if not self.shared == None:
            del self.shared
            self.shared = None

    def updateshared(self):
        ''' Check if this cell is shared between two or more pores.
            Is should be called only after pore identification. '''
        pores = defaultdict( list)
        for v in [x for x in self.vertexes if not x.overlaped]:
            pores[ v.cluster].append(v)
        vects = {}
        for k, v in pores.items():
            vects[k] = np.sum( np.array([x.x for x in v]), axis=0)
        xplanes = defaultdict( list)
        keys = list(pores.keys())
        size = len( keys)
        for k1 in range( size-1):
            for k2 in range(k1+1, size):
                p1 = keys[k1]
                p2 = keys[k2]
                u = get_unit( vects[ p1] - vects[ p2])
                xplanes[p1].append( u)
                xplanes[p2].append( -u)
        self.shared = xplanes

    def tovmd(self, offset=0.1, viewmode = ViewMode.ALWAYS, particle=False,
              vertexes=False, vradius=0.2, edges=False, ewidth=2, VMDresolution = 17, color="red3", fcolor="blue2"):
        '''Retruns draw commands to display this cell to vmd.'''
        ret = "draw color %s\n" % color
        # draw cell center (i.e. the coresponding particle) using commnet if the partilce flag is False
        ret += "" if particle else ";"
        if particle:
            ret += "draw sphere {%f %f %f} radius %f resolution %d\n" % tuple( list(self.x) + [ self.r/2, VMDresolution])
        ret += "draw color %s\n" % fcolor
        # ret += "draw material Transparent\n"
        if viewmode == ViewMode.ALWAYS: # print without checking overlap flag
            for f, n in zip(self.faces, self.normals):
                ret += f.tovmd(self.x, -np.array(n), offset=offset)
            ret += "draw color %s\n" % color
            # ret += "draw material Opaque\n"
            if vertexes:
                for v in self.vertexes:
                    if not v.overlaped: ret += v.tovmd(radius=vradius)
            if edges:
                for e in self.edges:
                    if not e.overlaped: ret += e.tovmd(width=ewidth)
        else:
            if viewmode == ViewMode.ANY:   # print face only if any of the vertexes is not overlaped
                for f, n in zip(self.faces, self.normals):
                    if any( [not x.overlaped for x in f.vertexes]):
                        ret += f.tovmd(self.x, -np.array(n), offset=offset)
            elif viewmode == ViewMode.ALL: # print face only if all the vertexes are not overlaped
                for f, n in zip(self.faces, self.normals):
                    if all( [not x.overlaped for x in f.vertexes]):
                        ret += f.tovmd(self.x, -np.array(n), offset=offset)
            ret += "draw color %s\n" % color
            # ret += "draw material Opaque\n"
            if vertexes: # print vertex if is not overlaped
                for v in self.vertexes:
                    if not v.overlaped: ret += v.tovmd(radius=vradius)
            if edges:    # print edge if is not overlaped
                for e in self.edges:
                    if not e.overlaped: ret += e.tovmd(width=ewidth)
        return ret
