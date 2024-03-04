# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from collections import defaultdict

import pysimpp.readers
from pysimpp.fastpost import fastbonds # pylint: disable=no-name-in-module
from pysimpp.utils.statisticsutils import Histogram

def _is_command(): return True
def _short_description(): return 'Calculate bond lengths distribution and time evolution or the mean value.'
def _command(): command()

# import scipy
# from scipy import interpolate
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt

def bonds(filename, bin, start, end, every):

    # check the input file
    reader = pysimpp.readers.create( filename)
    if not reader: sys.exit(0)

    dirname = reader.dir
    basename = reader.basename
    natoms = reader.natoms

    # set what to read
    attributes = 'id type x y z'
    reader.set_attributes(attributes)

    # check if a lammps dump file exist for bonds
    bndhandler = None
    if reader.__class__.__name__ == "LammpsReader":
        bndhandler = pysimpp.readers.LammpsReader.create_bond_dump_handler( reader)
    # if not retrieve bonds that do not change during the simulation
    if bndhandler is None:
        bonds = reader.get_bonds() - 1
        bond_type = reader.get_bond_type()
        bond_types = np.unique( bond_type)

    # keep bonds time evolution and distribution per type
    trj = defaultdict( list)
    hst = defaultdict( lambda: Histogram.free( bin, 0.0, addref=False))

    frames = []
    r = np.empty(shape=(natoms, 3), dtype=np.float32)  # coordinates
    print('>> reading dump file(s) ...')
    while (True):
        step, box, data = reader.read_next_frame()

        if step is None:
            break

        # retrieve bonds per frame
        if not bndhandler is None:
            _, _, bonds_data = pysimpp.readers.LammpsReader.read_next_frame_handler( bndhandler, step, box)
            bonds = np.column_stack( (bonds_data['bead1']-1, bonds_data['bead2']-1))
            bond_type = bonds_data['type']
            bond_types = np.unique( bond_type)

        if step < start:
            continue
        elif step > end:
            break

        frames.append( step)

        if not step % every == 0:
            continue

        np.copyto(r[:, 0], data['x'])
        np.copyto(r[:, 1], data['y'])
        np.copyto(r[:, 2], data['z'])

        b = r[bonds[:,1]] - r[bonds[:,0]]
        bl = fastbonds( r.T, box.va, box.vb, box.vc, bonds.T)
        for _t in bond_types:
            _b = bl[bond_type == _t]
            np.vectorize( hst[_t].add)( _b)
            trj[_t].append( (np.mean(_b), np.std(_b)))

    print()

    # dump time evolution
    _types = sorted( trj.keys())
    header = "# %-12s" % 'frame' + " ".join( [ "mean_t%-6s std_t%-7s" % ((str(_t),)*2) for _t in _types])
    f = open( dirname+os.sep+basename+"_bonds.dat", 'w')
    lines = [ "%-12s" % str(x) for x in frames]
    for _t in _types:
        for i, _v in enumerate(trj[_t]):
            lines[i] += " %-12g %-12g" % _v
    f.write( header + '\n')
    for line in lines: f.write( line + '\n')
    f.close()

    # dump distributions
    for k in list( hst.keys()):
        hst[k].write( dirname+os.sep+basename+"_b%s.hst" % k)

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(description=_short_description())

    # add arguments (self explaned)
    string = 'the path to the simulation trajectory file. A topology file ' + \
             'should be present in the same directory (preferably a tpr file).'
    parser.add_argument('path', default="."+os.sep, help=string)

    parser.add_argument('-bin', nargs=1, type=float, metavar='bin', default=[0.2], \
                       help='bin length (Å) for the distribution')

    parser.add_argument('-start', nargs=1, type=int, metavar='n', default=[-1], \
                       help='start processing form configuration n [inclusive]')

    parser.add_argument('-end', nargs=1, type=int, metavar='n', default=[sys.maxsize], \
                       help='stop processing at configuration n [inclusive]')

    parser.add_argument('-every', nargs=1, type=int, metavar='n', default=[1], \
                       help='processing frequency (every n configurations)')

    # parse the arguments
    args = parser.parse_args()

    print("INPUT")
    print("path : ", args.path)
    print("bin : ", args.bin[0])
    print("start : ", args.start[0])
    print("every : ", args.every[0])
    print("end : ", "max" if args.end[0] == sys.maxsize else args.end[0])

    bonds( args.path, args.bin[0], args.start[0], args.end[0], args.every[0])

if __name__ == '__main__':
    command()