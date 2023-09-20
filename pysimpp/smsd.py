# -*- coding: utf-8 -*-

import os
import sys
from collections import defaultdict
from operator import add
import pickle

import numpy as np
from scipy import stats

from pysimpp.msd import MSDException, MSDDataUtility, _dimsbookkeep, _msd_fft3, _msd_fft3x
import pysimpp.readers
from pysimpp.utils.utils import isrange, islist, read_ndx, ispositive, argparse_moleculestype
from pysimpp.utils.statisticsutils import Histogram

from pysimpp.fastpost import fasts1x, fastwrap, fastcom, fastcom_total

def _is_command(): return True
def _short_description(): return 'Calculate segmental mean square displacement.'
def _command(): command()

__debug = False
__rmcm = True

# @profile
def smsd( filename, dt=0.0, start=-1, end=sys.maxsize, every=1, dimensions=['xyz'],
         maxconfs=-1, segmentsfile=None, selectfile=None):

    # check the input file
    reader = pysimpp.readers.create(filename)
    if not reader: sys.exit(0)

    usecom = False

    # TODO check for gromacs: the provided trajectory is assumed allready unwrapped.
    _unwrap = True
    reader.set_unwrap( _unwrap)
    attributes = 'id x y z type'
    reader.set_attributes( attributes)

    natoms = reader.natoms
    masses = reader.get_atom_mass() # get atoms mass

    molecule = reader.get_atom_molecule() - 1   # atom molecule array (index @zero)
    nmolecules = molecule.max() + 1             # number of molecules
    # TODO add this functionality in reader
    mol_atoms = defaultdict( list)              # molecule atoms array
    for i, im in enumerate( molecule):
        mol_atoms[ im].append(i)

    # set selection
    selection = False
    if not selectfile is None:
        _basename = os.path.basename( selectfile.name)
        _selected = read_ndx( selectfile)
        for _name, _atoms in _selected.items(): # check valid indexes
            if np.any(_atoms>natoms) or np.any(_atoms<1):
                print("ERROR: group %s in %s contains invalid index(es)" % ( _name, _basename))
                return
        if not len( _selected) == 1: # check the number of groups included in the selection
            print('WARNING: more than one groups were found in %s' % _basename)
        _selected = list( set().union( *( _selected.values())))
        if len( _selected) == 0: # check for empty selection
            print('WARNING: selection was not found in %s' % _basename)
        else: # setup selection mask array and inform
            _selected = np.array( _selected) - 1 # zero based indes
            selected = np.full((natoms,), False)
            selected[ _selected] = True
            nselected = np.count_nonzero(selected)
            print("INFO: %d atoms were seleted" % nselected)
            selection = True

    # if not selection: # select everything
    #     selected = np.full((natoms,), True)
    #     nselected = natoms

    # set the segmens
    if not segmentsfile is None:
        _basename = os.path.basename( segmentsfile.name)
        _segments = read_ndx( segmentsfile)
        # check that all segments have the same size
        if not len( set( [len(x) for x in _segments.values()])) == 1:
            print('ERROR: all segments in %s should have the same size' % _basename)
            return
        # treat groups as individual molecules
        molecule = np.full( (natoms,), -1, dtype=np.int32)
        nmolecules = 0
        _chk = []

        for i, (k, v) in enumerate( _segments.items()):
            if np.any(v>natoms) or np.any(v<1):
                print("ERROR: segment %s in %s contains invalid index(es)" % ( k, _basename))
                return
            v = v - 1
            _n = np.count_nonzero( selected[ v]) if selection else v.size
            if _n == 0:
                print("WARNING: nothing is selected for segment %s in %s)" % ( k, _basename))
            else:
                molecule[ v] = nmolecules
                nmolecules += 1
                _chk.append( _n)

        if nmolecules == 0:
            print("ERROR: no valid segments were found in %s" % _basename)
            return
        if not np.all( np.array(_chk) == _chk[0]):
            print('ERROR: all segments in %s should have the same size' % _basename)
            return

        print("INFO: %d segments were found in %s" % ( nmolecules, _basename))
        if selection:
            print("INFO: slection is merged with segments")

        usecom = True

    # bookkeeping dimensions to read
    _dimensions = 'xyz'
    attributes, dmap, dmapinv = _dimsbookkeep( _dimensions, _unwrap)
    ndims = len( _dimensions)
    reader.set_attributes( attributes)

    print('>> count configurations in dump file(s) ...')
    nconfs = reader.count_frames( start, end) if  maxconfs == -1 else maxconfs

    # allocate unwrapped coordinates
    n_ = nmolecules if usecom else nselected if selection else natoms
    cm = np.zeros( (nconfs, n_, ndims), dtype=np.float64)
    r = np.empty( shape=(natoms,ndims), dtype=np.float64)  # unwrapped coordinates
    # initial system center of mass
    cm0 = np.zeros((3),dtype=np.float64)
    # vector for com move removal
    delta = np.zeros((3),dtype=np.float64)

    print('\n>> reading dump file(s) ...')
    steps = []
    boxes = []
    iconf = -1
    iconf_ = 0
    while( True):
        step, box, data = reader.read_next_frame()
        iconf_ += 1
        if step == None:
            break
        elif not iconf_%every == 0:
            continue
        elif step < start:
            continue
        elif step > end:
            print(step)
            break
        if box:
            iconf += 1
            steps.append( step)
            boxes.append( box)
            for k, v in dmap.items():
                np.copyto( r[:, k] ,  data[v[0]])

            if __rmcm:
                cmtotal, masstotal = fastcom_total( r, masses)
                if len(boxes) == 1:
                    cm0[:] = cmtotal
                else:
                    delta[:] = cmtotal-cm0
                    r -= delta

            # if needed calculate center of masses
            cm[iconf,:,:] = fastcom( r, masses, molecule, nmolecules) if usecom else r[selected] if selection else r

        else:
            break

    # number of configurations - resize if needed
    nconfs = len( steps)
    _s = cm.shape
    if not nconfs == _s[0]:
        print("to save some memory consider to use maxconfs = %d instead of %d" %( nconfs, _s[0]))
        cm = np.resize( cm, (nconfs, _s[1], _s[2]))

    # get the timestep in ps and handle some cases
    _dt = reader.timestep
    if _dt == 0.0:
        if dt == 0.0:
            dt = 1.0
            print("INFO: the timestep is set to %.2f fs. You shoud use '-dt' option to provide it." % dt)
    else:
        if dt > 0.0:
            print("INFO: the timestep is set to %.2f fs (provided with -dt option)." % dt)
            print("      Nevertheless, a timestep of %.2f fs is found in the trajectory files." % _dt)
        else:
            dt = _dt
    dt *= 1000.0 # TODO check for gromacs
    dt *= every

    print('\n>> calculate msds(s) ...')
    w = []

    msd, msds = _msd_fft3( cm)
    # msd, msds = _msd_fft3x( r, path=reader.dir+os.sep, w=w)

    print("msd was calculated")
    t = np.array( steps, dtype= np.float32) * dt
    du = MSDDataUtility( t, msd, msds, dimensions)
    f = open( reader.dir+os.sep+"msd.dat", 'w')
    du.write(f=f, var=0.10)
    f.close()

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(description='Calculate segmental mean square displacement.')

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
                       metavar='dim', default=['xyz'], \
                       help='produce the MSD for the given dimension(s)')

    parser.add_argument('-maxconfs', nargs=1, type=int, metavar='n', default=[-1], \
                       help='pre-allocate memmory for the given number of cofigurations.')

    parser.add_argument('-start', nargs=1, type=int, metavar='n', default=[-1], \
                       help='start processing form configuration START [inclusive]')

    parser.add_argument('-end', nargs=1, type=int, metavar='n', default=[sys.maxsize], \
                       help='stop processing at configuration END [inclusive]')

    parser.add_argument('-every', nargs=1, type=int, metavar='n', default=[1], \
                       help='process every EVERY configuration')

    def argdttype( string):
        val = ispositive( string, numbertype=float)
        if val is None:
            msg = "wrong integration timestep (check: %s)" % string
            raise argparse.ArgumentTypeError(msg)
        return val
    parser.add_argument('-dt', nargs=1, type=argdttype, default=[0.0], metavar='timestep', \
                       help='integration timeste in ps')

    message='a file with the gromacs style indexes for the segments to be considered. ' +\
            'The segments should be of the same kind since the msd of their center of ' +\
            'mass will be calculated. If no file is provided, the atom-based msd will ' +\
            'be calculated., '
    parser.add_argument('-segments', nargs=1, type=argparse.FileType('r'), metavar='file', required=False,
                        default=[None], help=message)

    message='a file with the gromacs style indexes for the atoms to be considered.'
    parser.add_argument('-select', nargs=1, type=argparse.FileType('r'), metavar='file', required=False,
                        default=[None], help=message)

    parser.add_argument('-gmxfit', dest='gmxfit', action='store_true',help="use gromacs approach for fitting msd(t) curve.")

    # parse the arguments
    args = parser.parse_args()

    print("INPUT")
    print("path : ", args.path)
    print("dt (ps): ", args.dt[0])
    print("start : ", args.start[0])
    print("end : ", "max" if args.end[0] == sys.maxsize else args.end[0])
    print("every : ", args.every[0])
    print("dim : ", args.dim[0])
    print("segments : ", "-" if args.segments[0] is None else args.segments[0].name)
    print("select : ", "-" if args.select[0] is None else args.select[0].name)
    print()

    smsd( args.path, dt=args.dt[0], start=args.start[0], end=args.end[0], every=args.every[0], dimensions=args.dim[0], \
         maxconfs=args.maxconfs[0], segmentsfile=args.segments[0], selectfile=args.select[0])

if __name__ == '__main__':
    command()
