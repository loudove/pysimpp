# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from collections import defaultdict
import itertools
from math import fabs, acos, copysign

import MDAnalysis

import pysimpp.readers

from pysimpp.utils.simulationbox import SimulationBox, NeighborCellList
from pysimpp.utils.statisticsutils import Histogram, Histogram2D
from pysimpp.utils.vectorutils import get_length, get_unit, get_angle, get_angle_unit, get_dihedral

def _is_command(): return True
def _short_description(): return 'Calculate orientation distribution (angle with z axis).'
def _command(): command()

# import scipy
# from scipy import interpolate
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt

__rad2deg = 180./np.pi

def dihedrals(filename, ndxfile, start, end, every):

    # check the input file
    reader = pysimpp.readers.create(filename)
    if not reader: sys.exit(0)

    dirname = reader.dir
    basename = reader.basename
    natoms = reader.natoms

    # set what to read
    attributes = 'id x y z type'
    reader.set_attributes(attributes)

    # read index file
    lines=ndxfile.readlines()
    ndxfile.close()
    # find types and keep angles quatrets per type
    _d=defaultdict(str)
    for line in lines:
        line=line.strip()
        if len(line) == 0:
            continue
        elif line.startswith("["):
            _type=line[1:-1].strip()
        else:
            _d[_type]+=" "+line
    d = {}
    h = {}
    for k in list(_d.keys()):
        a = np.array(list(map( int, _d[k].split()))) - 1
        _n = a.size
        if not _n % 4 == 0:
            print("Dihedral indexes should be given in quartets! Check %s type" % k)
            return
        d[k] = np.array(a).reshape(int(_n/4),4)
        h[k] = Histogram.free( 2.0, 0.0, addref=False)

    # buffer function
    dot = np.dot
    sqrt = np.sqrt
    nframes = 0
    r = np.empty(shape=(natoms, 3),
                  dtype=np.float32)  # coordinates
    print('>> reading dump file(s) ...')
    while (True):
        step, box, data = reader.read_next_frame()
        if step is None:
            break
        elif step < start:
            continue
        elif step > end:
            break

        nframes += 1

        if not step % every == 0:
            continue

        np.copyto(r[:, 0], data['x'])
        np.copyto(r[:, 1], data['y'])
        np.copyto(r[:, 2], data['z'])
        set_to_minimum = box.set_to_minimum

        for k in list(d.keys()):
            kd = d[k]
            kh = h[k]
            for _d in kd:
                _r = r[ _d,:]
                v1 = _r[1]-_r[0]
                v2 = _r[2]-_r[1]
                v3 = _r[3]-_r[2]
                set_to_minimum(v1)
                set_to_minimum(v2)
                set_to_minimum(v3)
                angle = get_dihedral( v1, v2, v3)*__rad2deg
                kh.add( angle)

    for k in list(h.keys()):
        h[k].write( dirname+os.sep+basename+"_%s.data" % k)

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(description=_short_description())

    # add arguments (self explaned)
    string = 'the path to the simulation trajectory file.A topology file' + \
             'should be present in the same directory (preferably a tpr file).'
    parser.add_argument('-ndx', nargs=1, type=argparse.FileType('r'))             
    parser.add_argument('path', default="."+os.sep,  \
                       help=string)
    parser.add_argument('-start', nargs=1, type=int, metavar='n', default=[-1], \
                       help='start processing form configuration n [inclusive]')
    parser.add_argument('-end', nargs=1, type=int, metavar='n', default=[sys.maxsize], \
                       help='stop processing at configuration n [inclusive]')
    parser.add_argument('-every', nargs=1, type=int, metavar='n', default=[1], \
                       help='processing frequency (every n configurations)')


    # parse the arguments
    args = parser.parse_args()

    print("INPUT")
    print("ndx : ", args.ndx[0].name)
    print("path : ", args.path)
    print("start : ", args.start[0])
    print("every : ", args.every[0])
    print("end : ", args.end[0])

    dihedrals( args.path, args.ndx[0], args.start[0], args.end[0], args.every[0])

if __name__ == '__main__':
    command()