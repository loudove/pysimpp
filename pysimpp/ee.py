#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
from scipy.signal import correlate

from pysimpp.utils.statisticsutils import Histogram
import pysimpp.utils.utils as utils
import pysimpp.readers
from pysimpp.fastpost import fastunwrapv, fastrg, pefastunwrap

__debug = False

def _is_command(): return True

def _short_description(): return 'Calculate end-to-end distances and stuff.'

def _command(): command()

def _FFT1D(x, normalize=False):
    ''' ACF using FFT. '''
    N=len(x)
    F = np.fft.fft(x, n=2*N)      # 2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res= (res[:N]).real           # now we have the autocorrelation in convention B
    if normalize:
        n=N*np.ones(N)-np.arange(0,N) # divide res(m) by (N-m)
        # n=np.arange(N, 0, -1)
        res = res/n
    return res                    # this is the autocorrelation in convention A

def __acf(t, data, fname):

    N = len(data)
    n = len(data[0])
    M = 3*n

    a=np.zeros( (M,N), dtype=float)
    for i, d in enumerate(data):
        ld=np.sqrt((d*d).sum(axis=1)) # convert to unit
        a_=d/ld[:,np.newaxis]
        a[:,i]=a_.reshape(M)
    acf=[]
    nrm_ = np.arange(N, 0, -1, dtype=float)
    for i, ai in enumerate(a):
        if i%3 == 0:
            if i: acf.append(s_/nrm_)
            s_=np.zeros(N, dtype=float)
        # tmp = np.correlate(ai, ai, mode='full')
        tmp = correlate(ai, ai, mode='full')
        # s_ += tmp[int(tmp.size/2):]/tmp[int(tmp.size/2)]
        s_ += tmp[int(tmp.size/2):]
        # s_ += _FFT1D(ai,normalize=False)

    acf=np.array(acf).mean(axis=0).transpose()

    f=open(fname,'w')
    for x1, x2 in np.vstack([t, acf]).T:
        f.write( "%.3f %.6f\n" % (x1, x2))
    f.close()

def endtoend(filename, 
    end1, 
    end2, 
    start=-1, 
    end=sys.maxsize, 
    molids=(), 
    unwrap=True, 
    camc=False):

    # check the ends
    if not camc:
        e1 = list(range(end1[0]-1, end1[1], end1[2]))
        e2 = list(range(end2[0]-1, end2[1], end2[2]))
        if not len(e1) == len(e2):
            print("ERROR: incompatible first/last atoms range.")
            return

    # create the reader. for the moment lammps, gromacs and camc (polybead) are supported.
    reader = pysimpp.readers.create(filename)

    if not reader: 
        print("ERROR: it was not possible to parse the given file.")
        return

    reader.set_unwrap( unwrap)
    attributes = 'id x y z'
    reader.set_attributes(attributes)


    print('>> reading data file ...')
    if not reader.read_topology():
        print("ERROR: there was no data file found.")
        return

    # get system  data
    natoms = reader.natoms
    types = reader.get_atom_type()  # atom type array (number or string based)
    # typemass = reader.get_type_mass()  # type mass
    # masses = np.array(
    #     [typemass[types[iat]] for iat in range(natoms)])  # atom mass array
    masses = reader.get_atom_mass()
    molecules = reader.get_atom_molecule() - 1  # atom molecule array (index @zero)
    nmolecules = molecules.max() + 1  # number of molecules

    # selected molecules
    selected = np.sort(molids) - 1  # selected molecules (index @zero)
    nselected = selected.size  # number of selected molecules
    hasselected = nselected > 0  # selected flag

    r = np.zeros((natoms, 3), dtype=np.float64)  # wrapped
    exclude = np.zeros((nmolecules), dtype=np.bool)  # all false
    cm = np.zeros((nmolecules, 3), dtype=np.float64)  # center of mass
    sqrg_ = np.zeros((nmolecules), dtype=np.float64)  # gyration tensors

    if hasselected:
        exclude[:] = True
        exclude[selected] = False
    else:
        selected = ()

    print('\n>> reading dump file(s) ...')
    steps = []
    boxes = []
    iconf = -1
    eevectors = []
    sqrg = []
    while (True):
        step, box, data = reader.read_next_frame()
        if step == None:
            break
        elif step < start:
            continue
        elif step > end:
            print(step)
            break
        if box:
            iconf += 1
            steps.append(step)
            boxes.append(box)
            # update the topology as needed in the case of camc
            if camc:
                natch, chat, e1_, e2_ = reader.get_topology()
                e1 = np.array(e1_)[selected]
                e2 = np.array(e2_)[selected]
                molecules = reader.get_atom_molecule() - 1
                masses = reader.get_atom_mass()
            r[:, 0] = data['x']
            r[:, 1] = data['y']
            r[:, 2] = data['z']
        else:
            break

        eev = r[e2] - r[e1]
        eevectors.append(eev)
        cm[:,:], sqrg__ = fastrg( r, masses, molecules, exclude)
        sqrg.append( sqrg__[selected])

    # number of configurations - resize if needed
    nconfs = len(steps)

    # create time array
    dt = reader.timestep / 1000.0 # form ps to ns
    time = np.array(steps, dtype=float) * dt

    # calculate the mean square end-to-end vector
    hsqree = Histogram.free(100.0, 0.0, addref=False)
    hree = Histogram.free(2.0, 0.0, addref=False)
    hmsqree = Histogram.free(100.0, 0.0, addref=False)
    hmree = Histogram.free(2.0, 0.0, addref=False)
    # and the mean square radius of gyration
    msqrg = []
    hsqrg = Histogram.free(20.0, 0.0, addref=False)
    hrg = Histogram.free(1.0, 0.0, addref=False)
    hmsqrg = Histogram.free(20.0, 0.0, addref=False)
    hmrg = Histogram.free(1.0, 0.0, addref=False)
    f1 = open(reader.dir+os.sep+"ree.data","w")
    f1.write("# time  msqree sqrt(msqree)\n")
    f2 = open(reader.dir+os.sep+"rg.data","w")
    f2.write("# time  msqree sqrt(msqrg)\n")
    for t, eev, sqrg_ in zip(time, eevectors, sqrg):
        # collect end-to-end vector
        eev_ = np.array(eev)
        sqeev_ = (eev_*eev_).sum(axis=1)
        leev_ = np.sqrt( sqeev_)
        for sqv, v in zip(sqeev_,leev_):
            hsqree.add( sqv)
            hree.add( v)
        msqv = np.mean(sqeev_)
        hmsqree.add( msqv)
        #mv = np.sqrt(msqv)
        mv = np.mean(leev_)
        hmree.add( mv)
        f1.write("%f %f %f\n" % (t, msqv, mv))
        # collect radious of gyration
        for val in sqrg_:
            hsqrg.add( val)
            hrg.add( np.sqrt(val))
        msqrg_ = sqrg_.mean()
        msqrg.append( msqrg_)
        hmsqrg.add( msqrg_)
        mrg_ = np.sqrt(msqrg_)
        hmrg.add( mrg_)
        f2.write("%f %f %f\n" % (t, msqrg_, mrg_))

    f1.close()
    f2.close()
    # output end-to-end data
    hsqree.write(reader.dir+os.sep+"sqree_hist.data")
    hree.write(reader.dir+os.sep+"ree_hist.data")
    hmsqree.write(reader.dir+os.sep+"msqree_hist.data")
    hmree.write(reader.dir+os.sep+"mree_hist.data")
    # output radious of gyration data
    hsqrg.write(reader.dir+os.sep+"sqrg_hist.data")
    hrg.write(reader.dir+os.sep+"rg_hist.data")
    hmsqrg.write(reader.dir+os.sep+"msqrg_hist.data")
    hmrg.write(reader.dir+os.sep+"mrg_hist.data")

    # mean square end-to-end vector autocorrelation
    __acf(time, eevectors, reader.dir+os.sep+"ree_acf.data")

    print()

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(
        description=_short_description(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=''' \
The output files are (located at the simulation directory):
    ree.data : mean square end-to-end and mean end-to-end length for each frame
    ree_acf.data : end-to-end vector autocorrelation function
    rg.data : mean square radius of gyration and mean radius of gyration for each frame
    sqree_hist.data : square end-to-end vector distribution
    msqree_hist.data : mean square end-to-end distribution
    ree_hist.data : end-to-end length distribution
    mree_hist.data : mean end-to-end length distribution
    sqrg_hist.data : square radius of gyration distribution
    msqrg_hist.data : mean square radius of gyration distribution
    rg_hist.data : radius of gyration distribution
    mrg_hist.data : mean radius of gyration distribution ''')

    # add arguments (self explaned)
    string = 'the path to the simulation log file. In the case of gromacs' +\
             'simulation, a topology file should be present in the same' +\
             'directory (preferably a tpr file). In the case of lammps the' +\
             'data and dump files will be traced from the corresponding' +\
             'log records, otherwise a data and a dump file with the same' +\
             'base name as the log file should exist in the same directory.'

    parser.add_argument('path', default="."+os.sep,  \
                       help=string)
    parser.add_argument('-start', nargs=1, type=int, metavar='n', default=[-1], \
                       help='start processing form configuration START [inclusive]')
    parser.add_argument('-end', nargs=1, type=int, metavar='n', default=[sys.maxsize], \
                       help='stop processing at configuration END [inclusive]')

    def argends(string):
        ''' check the "-end1" and "-end2" option arguments. '''
        numbers = utils.isrange(string, sep=":", positive=True)
        if len(numbers) != 3:
            msg = "wrong atoms range for (check: %s)" % string
            raise argparse.ArgumentTypeError(msg)
        return numbers
    parser.add_argument('-end1', nargs=1, type=argends, metavar='start atom', default=[()], \
                       help='range (start:end:step) of atoms to be used as molecules first atom')
    parser.add_argument('-end2', nargs=1, type=argends, metavar='end atom', default=[()], \
                       help='range (start:end:step) of atoms to be used as molecules last atom')

    def argmoleculestype(string):
        ''' check the "-molecules" option arguments. '''
        if len(string) == 0:
            return []
        numbers = utils.isrange(string, positive=True)
        if len(numbers) == 0:
            msg = "wrong molecules indexs range (check: %s)" % string
            raise argparse.ArgumentTypeError(msg)
        return numbers
    parser.add_argument('-molecules', nargs=1, type=argmoleculestype, default=[[]],  metavar='range', \
                       help='molecules to be used. A list with comma seperated id ranges should be provided e.g. "1,2,3" or "1:10,20,30:100"')

    parser.add_argument('--no-unwrap', dest='unwrap', default=True, action='store_false', \
                       help="do not unwrap molecules sine the provided trajectories provide unwrapped coordinates.")

    # parser.add_argument('--camc', dest='camc', default=False, action='store_true', \
    #                    help="process connectivity monte carlo output")

    # parse the arguments
    args = parser.parse_args()

    print("INPUT")
    print("path : ", args.path)
    print("start : ", args.start[0])
    print("end : ", args.end[0])
    print("molecules : %d" % len(args.molecules[0]))
    print("unwrap : %s \n" % ("True" if args.unwrap else "False"))
    # print("camc : %s \n" % ("True" if args.camc else "False"))
    if __debug:
        print(args.molecules[0])

    endtoend(
        args.path,
        args.end1[0],
        args.end2[0],
        start=args.start[0],
        end=args.end[0],
        molids=args.molecules[0],
        unwrap=args.unwrap)
        # camc=args.camc)

if __name__ == '__main__':
    command()
