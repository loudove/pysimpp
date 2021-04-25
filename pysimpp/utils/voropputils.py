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
