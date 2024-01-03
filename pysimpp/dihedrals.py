# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from collections import defaultdict

import pysimpp.readers

from pysimpp.fastpost import fastdihedrals # pylint: disable=no-name-in-module
from pysimpp.utils.utils import IsList, read_ndx
from pysimpp.utils.statisticsutils import Histogram, Histogram2D
from pysimpp.utils.vectorutils import get_dihedral

def _is_command(): return True
def _short_description(): return 'Calculate then dihedrals distribution (simple and joint).'
def _command(): command()

# import scipy
# from scipy import interpolate
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt

__rad2deg = 180./np.pi

def dihedrals(filename, ndxfile, jnt, start, end, every):

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
    ndxs = read_ndx(ndxfile)
    ndxfile.close()

    for k in jnt:
        if not k in ndxs:
            print("ERROR: group %s cannot be found in %s " % (k,ndxfile.name))
            return

    # identify the dihedrals to be calculated and assign them a unique id.
    # for each type keep a list with angles id and the the corresponding histogram.
    d = {}
    h = {}
    _phi_name = {} # dihedral {name:id} dictionary
    _phi_seq = [] # list of dihedral sequences; the index corresponds to the id
    for k, v in ndxs.items():
        # basic check 
        kinjnt = k in jnt 
        _n = v.size
        if not _n % 4 == 0:
            print("ERROR: indexes should be given in quartets! Check %s type" % k)
            return
        elif kinjnt and not _n % 4 == 0:
            print("ERROR: indexes should be given in pairs of quartets! Check %s type" % k)
            return
        # parse the quartets and assign them a id
        idl = [] # list of dihedrals id
        for _d in np.array(v-1).reshape(int(_n/4),4):
            # get dihedral unique name 
            _d0 = _d[0]; _d3 = _d[3]
            _name = "%d %d %d %d" % tuple(_d) if _d0 < _d3 else (_d3, _d[2], _d[1], _d0)
            # find the id
            if not _name in _phi_name:
                _id = len( _phi_name)
                _phi_name[ _name] = _id
                _phi_seq.append( _d)
            else:
                _id = _phi_name[ _name]
            idl.append( _id)
        # for each type keep the idexes of the dihedrals and the construct the histograms
        _bin = 2.0 * np.pi / 180.0
        if kinjnt:
            d[k] = np.array(idl).reshape(len(idl)//2,2)
            h[k] = Histogram2D((_bin, _bin), (0, 0), addref=False)
        else:
            d[k] = np.array(idl)
            h[k] = Histogram.free( _bin, _bin, addref=False)

    _phi_seq = np.array( _phi_seq)

    # buffer function
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

        _phi = fastdihedrals(r.T,box.va,box.vb,box.vc,_phi_seq.T)

        for k in list(d.keys()):
            if k in jnt:
                np.vectorize(h[k].add, signature='(n)->()')(_phi[d[k]])
            else:
                np.vectorize(h[k].add)(_phi[d[k]])

    for k in list(h.keys()):
        if k in jnt:
            h[k].normalize()
        h[k].write( dirname+os.sep+basename+"_%s.data" % k)

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(description=_short_description())

    # add arguments (self explaned)
    string = 'the path to the simulation trajectory file. A topology file' + \
             'should be present in the same directory (preferably a tpr file).'
    parser.add_argument('path', default="."+os.sep, help=string)

    message='''
    a file with the gromacs style indexes for the dihedrals to be considered.
    The calculated probability distribution is written in the trajectory file
    direcory in a file {trjbasename}_{group}.data. '''    
    parser.add_argument('-ndx', nargs=1, type=argparse.FileType('r'), metavar='file',
                        required=True, help=message)
    
    chktype = IsList("wrong group names (check: %s)",itemtype=str)
    message='''
    a comma separated list with the names of the groups in the index file is to
    be considered for a correlation check. The groups should contain pairs of
    dihedrals, the joint probability of which will be calculated. Also in this
    case the name of the output file will be {trjbasename}_{group}.data.'''
    parser.add_argument('-j', nargs=1, type=chktype, metavar='joint', default=[[]], \
                       help=message)
    
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
    print("joint ", ",".join( args.j[0]) if len(args.j[0]) > 0 else "-")
    print("path : ", args.path)
    print("start : ", args.start[0])
    print("every : ", args.every[0])
    print("end : ", "max" if args.end[0] == sys.maxsize else args.end[0])
    print()

    dihedrals( args.path, args.ndx[0], args.j[0], args.start[0], args.end[0], args.every[0])

if __name__ == '__main__':
    command()
