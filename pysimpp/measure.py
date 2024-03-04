# -*- coding: utf-8 -*-
import os
import sys
import numpy as np

import pysimpp.readers
from pysimpp.fastpost import fastbonds, fastangles, fastdihedrals # pylint: disable=no-name-in-module
from pysimpp.utils.utils import IsList, read_ndx
from pysimpp.utils.statisticsutils import Histogram

def _is_command(): return True
def _short_description(): return 'Calculate orientation distribution (angle with z axis).'
def _command(): command()

# import scipy
# from scipy import interpolate
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt

__rad2deg = 180./np.pi
__deg2rad = np.pi/180.

def measure(filename, measurement, ndxfile, start, end, every, groups=[], doacf=False):

    # create the reader and set what to read
    reader = pysimpp.readers.create(filename)
    if not reader: sys.exit(0)
    attributes = 'id x y z type'
    reader.set_attributes(attributes)

    # access simulation directory and trajectory basename
    dirname = reader.dir
    basename = reader.basename

    # setup measurement type stuff
    _multi={'distanse':2, 'angle':3, 'dihedral':4}[measurement]
    _func={'distanse':fastbonds, 'angle':fastangles, 'dihedral':fastdihedrals}[measurement]
    _bin={'distanse':0.05, 'angle':2.0*__deg2rad, 'dihedral':2.0*__deg2rad}[measurement]

    # read the index file
    ndxs = read_ndx(ndxfile)
    ndxfile.close()
    # basic checks, convert serial to index, create fastpost methods
    # arguments and the histograms
    d = {}
    h = {}
    for k, v in ndxs.items():
        if not v.size % _multi == 0:
            print("ERROR: dihedral indexes should be given in quartets! Check %s type" % k)
            return
        if len(groups) == 0 or k in groups:
            d[k] = v.reshape((int(v.size/_multi),_multi),order='C')-1
            h[k] = Histogram.free( _bin, 0.0, addref=False)

    # buffer function
    nframes = 0
    r = np.empty(shape=(reader.natoms, 3), dtype=np.float32,order='C')
    while (True):
        step, box, data = reader.read_next_frame()
        nframes += 1
        if step is None:
            break
        elif nframes < start:
            continue
        elif nframes > end:
            break

        if not nframes % every == 0:
            continue

        np.copyto(r[:, 0], data['x'])
        np.copyto(r[:, 1], data['y'])
        np.copyto(r[:, 2], data['z'])

        for k in list(d.keys()):
            kh = h[k]
            for v in _func(r.T,box.va,box.vb,box.vc,d[k].T):
                kh.add( v)

    for k in list(h.keys()):
        h[k].write( dirname+os.sep+basename+"_%s.data" % k)

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(description=_short_description())

    # add arguments (self explaned)
    string = 'the path to the simulation trajectory file. A topology file ' + \
             'should be present in the same directory (preferably a tpr file).'
    parser.add_argument('path', default="."+os.sep, help=string)

    message='''
    type of the measurement. The number of indexes for each group in the index file should
    be multiple of 2, 3, and 4 for "distance", "angle" and "dihedral" measurement type
    respectively. '''
    parser.add_argument('-type',nargs=1, metavar='type', choices=['distanse', 'angle', 'dihedral'],
                        required=True, help=message)

    parser.add_argument('-start', nargs=1, type=int, metavar='n', default=[-1], \
                       help='start processing form frame n [inclusive].')

    parser.add_argument('-end', nargs=1, type=int, metavar='n', default=[sys.maxsize], \
                       help='stop processing at frame n [inclusive].')

    parser.add_argument('-every', nargs=1, type=int, metavar='n', default=[1], \
                       help='processing frequency (every n frames).')

    message='''
    a file with the gromacs style indexes for the bonded items to be considered.
    The calculated probability distribution is written in the trajectory file
    direcory in a file {trjbasename}_{group}.dat. '''
    parser.add_argument('-ndx', nargs=1, type=argparse.FileType('r'), metavar='file', required=True,
                        help=message)

    chktype = IsList("wrong group names (check: %s)",itemtype=str)
    message='''
    a comma separated list with the name of the groups in the index file to be considered.'''
    parser.add_argument('-group', nargs=1, type=chktype, metavar='group', default=[[]], \
                       help=message)

    message='''
    calculate the total autocorrelation function for each group considered. The calculated
    acf is written in the trajectory directory is the file trjbasename}_acf_{group}.dat. '''
    parser.add_argument('--acf', dest='doacf', default=False, action='store_true',
                        help=message)

    # parse the arguments
    args = parser.parse_args()

    print("INPUT")
    print("type : ", args.type[0])
    print("ndx : ", args.ndx[0].name)
    print("group(s) : ", ",".join(args.group[0]) if len(args.group[0])>0 else '-')
    print("acf :", 'yes' if args.doacf else 'no')
    print("path : ", args.path)
    print("start : ", args.start[0])
    print("end : ", "max" if args.end[0] == sys.maxsize else args.end[0])
    print("every : ", args.every[0])

    measure(args.path, args.type[0], args.ndx[0], args.start[0], args.end[0],
              args.every[0], groups=args.group[0], doacf=args.doacf)

if __name__ == '__main__':
    command()
