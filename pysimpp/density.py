#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from collections import defaultdict
from operator import add
import pickle

import numpy as np
from scipy import stats

import pysimpp.readers
from pysimpp.utils.utils import argparse_moleculestype
from pysimpp.utils.statisticsutils import Histogram
from pysimpp.fastpost import fastunwrapv, fastcom, fastcom_total

__debug = False

def _is_command(): return True

def _short_description(): return 'Calculate probability density profiles along x, y, and/or z axes.'

def _command(): command()

def _dimsbookkeep( dimensions, unwrap=True):
    ''' Given the dimensions i.e. a list of strings, creates the attributes
        string for the lammpsreader and the dimensions map for data storage. '''

    # set what to read and store bookkeeping
    attributes = 'id type'
    xstring = ""
    istring = ""
    dmap = {}
    dmapinv = {}
    ndims = 0
    for dim in [ 'x', 'y', 'z']:
        if dim in dimensions:
            xstring += " %s" % dim
            idim = "i%s" % dim
            istring += " %s" % idim
            dmap[ ndims] = (dim, idim)
            dmapinv[ dim] = ndims
            ndims += 1
    attributes = attributes + xstring
    if unwrap:
        attributes += istring
    return attributes, dmap, dmapinv

# @profile
def density( filename, bin, start=-1, end=sys.maxsize, every=1, dimensions=['z'], molids=(), dowrap=True, usecom=True):

    # create the reader
    reader = pysimpp.readers.create(filename)
    if not reader:
        print("ERROR: it was not possible to parse the given file.")
        return

    # remove center of mass flag
    rmcom = False
    # access wrapped coordinates
    reader.set_wrap( dowrap)
    _unwrap = not dowrap

    attributes = 'id mol x y z'
    reader.set_attributes( attributes)

    print('>> reading data file ...')
    reader.read_topology()

    # TODO create a utility method in pysimpp.readers for accessing molecular data
    natoms = reader.natoms
    # get molecular data and select molecules of interest
    # TODO refactor massesar
    types = reader.get_atom_type()              # atom type array (number or string based)
    masses = reader.get_type_mass()             # type mass
    massesar = np.array( [ masses[ types[iat]] for iat in range( natoms) ]) # atom mass array (accelerate)
    molecules = reader.get_atom_molecule() - 1  # atom molecule array (index @zero)
    nmolecules = molecules.max() + 1            # number of molecules
    # TODO add this functionality in reader
    mol_atoms = defaultdict( list)              # molecule atoms array
    for i, im in enumerate( molecules):
        mol_atoms[ im].append(i)

    # selected molecules
    selected = np.sort( molids) - 1             # selected molecules (index @zero)
    if not usecom:                              # convert to atoms
        satoms = []
        for im in selected:
            satoms += mol_atoms[im]
        selected = np.sort( satoms)

    nselected = selected.size                   # number of selected molecules
    hasselected = nselected > 0                 # selected flag

    # bookkeeping dimensions to read
    attributes, dmap, dmapinv = _dimsbookkeep( 'xyz', _unwrap)

    _dims = tuple( [ ['x','y','z'].index( _d) for _d in dimensions] )
    _dnames = { ['x','y','z'].index( _d):_d for _d in dimensions}
    reader.set_attributes( attributes)

    # allocate wrapped and unwrapped molecules coordinates
    n_ = nselected if hasselected else nmolecules if usecom else natoms
    cm = np.zeros( (n_, 3), dtype=np.float32)              # com coordinates
    r = np.empty( shape=(natoms,3), dtype=np.float32)      # unwrapped coordinates
    # initial system com
    cm0 = np.zeros((3),dtype=np.float32)
    # vector for com removal
    delta = np.zeros((3),dtype=np.float32)

    print('\n>> reading dump file(s) ...')
    steps = []
    boxes = []
    profiles = {}
    iconf = 0
    while( True):
        step, box, data = reader.read_next_frame()
        iconf += 1
        if step == None:
            break
        elif not iconf % every == 0:
            continue
        elif step < start:
            continue
        elif step > end:
            print(step)
            break
        if box:

            if len(profiles) == 0:
                for _d in _dims:
                    origin = box.origin[_d]
                    length = origin + (box.a,box.b,box.c)[_d]
                    # nbins=round(length/bin)
                    profiles[_d] = Histogram.fixed( origin, length, bin)

            steps.append( step)
            boxes.append( box)

            for k, v in dmap.items():
                np.copyto(r[:, k] ,  data[v[0]])

            if rmcom:
                cmtotal, masstotal = fastcom_total( r, massesar)
                if len(boxes) == 1:
                    cm0[:] = cmtotal
                else:
                    delta[:] = cmtotal-cm0
                    r -= delta

            # calculate center of masses and if regions
            # set also the wraped coordinates to ensure
            # the correct bining
            if hasselected:
                cm_ = fastcom( r, massesar, molecules, nmolecules) if usecom else r
                cm[:,:] = cm_[selected,:]
            else:
                cm[:,:] = fastcom( r, massesar, molecules, nmolecules) if usecom else r

            for _d in _dims:
                profile = profiles[_d]
                for _cm in cm:
                    profile.add(_cm[_d])

        else:
            break

    for _d, profile in profiles.items():
        _dname = "%s%s"%("com_" if usecom else "", _dnames[ _d])
        header="# %s_[A] probability_[1/A]" % ( _dname)
        profile.write(reader.dir+os.sep+"%s_prf.data" % _dname, header=header)

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(description='Calculate probability density profiles along x, y, and/or z axes.')

    # add arguments (self explaned)
    string = 'the path to the simulation file. In the case of gromacs' +\
             'simulation, a topology file should be present in the same' +\
             'directory (preferably a tpr file). In the case of lammps the' +\
             'data and dump files will be traced from the corresponding' +\
             'log records, otherwise a data and a dump file with the same' +\
             'base name as the log file should exist in the same directory.'

    parser.add_argument('path', default="."+os.sep,  \
                       help=string)

    parser.add_argument('-dim', nargs=1, choices=[ 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'], \
                       metavar='dim', default=['z'], \
                       help='produce the density for the given dimension(s)')

    parser.add_argument('-bin', nargs=1, type=float, metavar='bin', required=True, \
                       help='bin width in A')

    parser.add_argument('-start', nargs=1, type=int, metavar='n', default=[-1], \
                       help='start processing form configuration START [inclusive]')

    parser.add_argument('-end', nargs=1, type=int, metavar='n', default=[sys.maxsize], \
                       help='stop processing at configuration END [inclusive]')

    parser.add_argument('-every', nargs=1, type=int, metavar='n', default=[1], \
                       help='process every EVERY configuration')

    parser.add_argument('-molecules', nargs=1, type=argparse_moleculestype, default=[[]],  metavar='range', \
                       help='molecules to be used. A list with comma seperated id ranges should be provided e.g. "1,2,3" or "1:10,20,30:100"')

    parser.add_argument('--no-wrap', dest='dowrap', default=True, action='store_false', \
                       help="do not wrap molecules sine the provided trajectories provide wrapped coordinates.")

    parser.add_argument('--no-com', dest='usecom', default=True, action='store_false', \
                       help="calculate the distributions using atoms rather that molecules' centers of mass.")

    # parse the arguments
    args = parser.parse_args()

    print("INPUT")
    print("path  : ", args.path)
    print("start : ", args.start[0])
    print("end   : ", args.end[0])
    print("every : ", args.every[0])
    print("dim   : ", args.dim[0])
    print("wrap  : %s" % ("True" if args.dowrap else "False"))
    print("com   : %s" % ("True" if args.usecom else "False"))

    print("molecules : %d \n" % len(args.molecules[0]))
    if __debug:
        print(args.molecules[0])
    print()

    _dims = [_d for _d in args.dim[0]]
    density( args.path, args.bin[0], start=args.start[0], end=args.end[0], every=args.every[0], \
        dimensions=_dims, molids=args.molecules[0], dowrap=args.dowrap, usecom=args.usecom)

if __name__ == '__main__':
    command()
