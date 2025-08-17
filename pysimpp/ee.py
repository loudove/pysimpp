#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
from scipy.signal import correlate
from collections import defaultdict

from pysimpp.utils.statisticsutils import Histogram
import pysimpp.utils.utils as utils
import pysimpp.readers
# from pysimpp.msd import _FFT1D
from pysimpp.fastpost import fastrg

__debug = False

def _is_command(): return True

def _short_description(): return 'Calculate end-to-end distances and stuff.'

def _command(): command()

def _acf(t, data, fname):

    N = len(data)
    n = len(data[0])
    M = 3*n

    # conventional calculation for testing
    # # convert to unit
    # a=np.zeros((N,n,3),dtype=float)
    # for i, d in enumerate(data):
    #     l_=np.sqrt((d*d).sum(axis=1))
    #     a[i,:,:]=d/l_[:,np.newaxis]
    # acf=np.zeros(N,dtype=float)
    # count=np.zeros(N,dtype=float)
    # for i in range(N):
    #     ai=a[i]
    #     for j in range(i,N):
    #         aj=a[j]
    #         d=j-i
    #         for k in range(n):
    #             acf[d]+=np.sum(aj[k]*ai[k])
    #             count[d]+=1.0
    # acf /= count

    a=np.zeros( (M,N), dtype=float)
    for i, d in enumerate(data):
        ld=np.sqrt((d*d).sum(axis=1,keepdims=True)) # convert to unit
        a[:,i]=(d/ld).reshape(M)
    acf=[]
    nrm_ = np.arange(N, 0, -1, dtype=float)
    for i, ai in enumerate(a):
        if i%3 == 0:
            if i: acf.append(s_/nrm_)
            s_=np.zeros(N, dtype=float)
        # tmp = np.correlate(ai, ai, mode='full')
        tmp = correlate(ai, ai, mode='full')
        # s_ += tmp[int(tmp.size/2):]/tmp[int(tmp.size/2)]
        s_ += tmp[tmp.size//2:]
        # s_ += _FFT1D(ai,normalize=False)

    acf=np.array(acf).mean(axis=0).transpose()

    f=open(fname,'w')
    t0 = t[0]
    for x1, x2 in np.vstack([t, acf]).T:
        f.write( "%.3f %.6f\n" % (x1-t0, x2))
    f.close()

def endtoend(filename,
    end1,
    end2,
    dt=0.0,
    start=-1,
    every=1,
    end=sys.maxsize,
    molids=(),
    whole=True,
    persistence=False,
    camc=False,
    bins={}):

    # check the ends
    if not camc:
        e1 = list(range(end1[0]-1, end1[1]+1, end1[2]))
        e1p = np.array( range(end1[0], end1[1]+1, end1[2]))
        e2 = list(range(end2[0]-1, end2[1]+1, end2[2]))
        e2p = np.array( range(end2[0], end2[1]+1, end2[2]))
        if not len(e1) == len(e2):
            print("ERROR: incompatible first/last atoms range.")
            return

    # update the bins for the histograms
    _bins = {"pl":1.0,"bnd":0.2}
    for k, v in bins.items(): _bins[k] = v[0]

    # create the reader. for the moment lammps, gromacs and camc (polybead) are supported.
    reader = pysimpp.readers.create(filename)

    if not reader:
        print("ERROR: it was not possible to parse the given file.")
        return

    reader.set_whole( whole)
    attributes = 'id x y z'
    reader.set_attributes(attributes)

    # print('>> reading data file ...')
    # if not reader.read_topology():
    #     print("ERROR: there was no data file found.")
    #     return

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

    r = np.zeros((natoms, 3), dtype=np.float64)  # whole
    exclude = np.zeros((nmolecules), dtype=np.bool_)  # all false
    cm = np.zeros((nmolecules, 3), dtype=np.float64)  # center of mass
    # sqrg_ = np.zeros((nmolecules), dtype=np.float64)  # gyration tensors

    if hasselected:
        exclude[:] = True
        exclude[selected] = False
    else:
        selected = ()

    print('\n>> reading dump file(s) ...')
    steps = []
    boxes = []
    eevectors = []
    sqrg = []

    # persistence length stuff
    if persistence:
        bnd = np.zeros( (natoms,3), dtype=np.float64)
        natch_ = natoms//nmolecules
        nbndch_ = natch_ - 1
        bndacf = []
        hlp = Histogram.free(_bins["pl"], 0.0, addref=False)
        _hlp = defaultdict( lambda: Histogram.free(_bins["pl"], 0.0, addref=False))
        hbndl = Histogram.free(_bins["bnd"], 0.0, addref=False)
        nrm_ = np.arange( nbndch_, 0, -1, dtype=float)

    iframe = -1
    while (True):
        step, box, data = reader.read_next_frame()
        iframe += 1

        if step == None:
            break
        elif step < start:
            continue
        elif not iframe % every == 0:
            continue
        elif step > end:
            break

        steps.append(step)
        boxes.append(box)
    
        # update the topology as needed in the case of camc
        if camc:
            natch, chat, e1_, e2_ = reader.get_topology()
            e1 = np.array(e1_)[selected]
            e2 = np.array(e2_)[selected]
            e1p = e1+1
            e2p = e2+1
            molecules = reader.get_atom_molecule() - 1
            masses = reader.get_atom_mass()
        r[:, 0] = data['x']
        r[:, 1] = data['y']
        r[:, 2] = data['z']

        eev = r[e2] - r[e1]
        eevectors.append(eev)
        cm[:,:], sqrg__ = fastrg( r, masses, molecules, exclude)
        sqrg.append( sqrg__[selected])

        if persistence:
            bnd[1:,:] = r[1:,:] - r[:-1,:]
            bnd[e1,:] = eev
            bndl = np.sqrt( (bnd*bnd).sum(axis=1,keepdims=True))
            bnd[:] = bnd / bndl
            # flory definition
            np.vectorize(hlp.add)( ( eev*bnd[e1p,:]).sum(axis=1))
            for k in range(nbndch_):
                np.vectorize( _hlp[k].add)( ( eev*bnd[e1p+k,:]).sum(axis=1))
            # bond length distribution
            np.vectorize(hbndl.add)( bndl)
            # acf definition
            for i, j in zip(e1p,e2p):
                v = bnd[i:j].transpose()
                s_ = np.zeros(nbndch_, dtype=float)
                for i, vi in enumerate(v):
                    tmp = correlate(vi, vi, mode='full')
                    s_ += tmp[tmp.size//2:]
                s_/=nrm_
                bndacf.append(s_)

    # create time array
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
    dt /= 1000.0 # form ps to ns
    time = np.array(steps, dtype=float) * dt

    # precalculate quantities ot check the bins
    def setbin(a, name, n=20):
        if name in _bins:
            return _bins[name]
        else:
            bin = 2.5 * a.std() / n
            return np.round(bin,2) if bin > 1. else bin

    sqeev = np.array([ (eev_*eev_).sum(axis=1) for eev_ in eevectors ])
    msqeev = np.mean(sqeev, axis=1) # mean square end-to-end vector
    leev = np.sqrt(sqeev)   # end-to-end vector length
    mleev = np.mean(leev, axis=1) # mean end-to-end vector length
    sqrg = np.array(sqrg) # square radius of gyration
    msqrg = np.mean(sqrg, axis=1) # mean square radius of gyration
    rg = np.sqrt(sqrg) # radius of gyration
    mrg = np.sqrt(msqrg) # mean radius of gyration

    # calculate the mean square end-to-end vector
    hsqree = Histogram.free(setbin(sqeev,"sqree"), 0.0, addref=False)
    np.vectorize(hsqree.add)(sqeev.flat)
    hmsqree = Histogram.free(setbin(msqeev,"msqree"), 0.0, addref=False)
    np.vectorize(hmsqree.add)(msqeev)
    hree = Histogram.free(setbin(leev,"ree"), 0.0, addref=False)
    np.vectorize(hree.add)(leev.flat)
    hmree = Histogram.free(setbin(mleev,"mree"), 0.0, addref=False)
    np.vectorize(hmree.add)(mleev)

    # and the mean square radius of gyration
    hsqrg = Histogram.free(setbin(sqrg,"sqrg"), 0.0, addref=False)
    np.vectorize(hsqrg.add)(sqrg.flat)
    hmsqrg = Histogram.free(setbin(msqrg,"msqrg"), 0.0, addref=False)
    np.vectorize(hmsqrg.add)(msqrg)
    hrg = Histogram.free(setbin(rg,"rg"), 0.0, addref=False)
    np.vectorize(hrg.add)(rg.flat)
    hmrg = Histogram.free(setbin(mrg,"mrg"), 0.0, addref=False)    
    np.vectorize(hmrg.add)(mrg)

    # output time series data
    f1 = open(reader.dir+os.sep+"ree.data","w")
    f1.write("# time[ns]  msqree[Å²] sqrt(msqree[Å]\n")
    f2 = open(reader.dir+os.sep+"rg.data","w")
    f2.write("# time[ns]  msqree[Å²] sqrt(msqrg)[Å]\n")
    for t, msqeev_, mleev_, msqrg_, mrg_ in zip(time, msqeev, mleev, msqrg, mrg):
        f1.write("%f %f %f\n" % (t, msqeev_, mleev_))
        f2.write("%f %f %f\n" % (t, msqrg_, mrg_))
    f1.close()
    f2.close()

    # output end-to-end data
    hsqree.write(reader.dir+os.sep+"sqree_hist.data", fmt="%f %g", header="# %s"%str(hsqree.variable))
    hree.write(reader.dir+os.sep+"ree_hist.data", header="# %s"%str(hree.variable))
    hmsqree.write(reader.dir+os.sep+"msqree_hist.data", header="# %s"%str(hmsqree.variable))
    hmree.write(reader.dir+os.sep+"mree_hist.data", header="# %s"%str(hmree.variable))
    # output radious of gyration data
    hsqrg.write(reader.dir+os.sep+"sqrg_hist.data", fmt="%f %g", header="# %s"%str(hsqrg.variable))
    hrg.write(reader.dir+os.sep+"rg_hist.data", header="# %s"%str(hrg.variable))
    hmsqrg.write(reader.dir+os.sep+"msqrg_hist.data", header="# %s"%str(hmsqrg.variable))
    hmrg.write(reader.dir+os.sep+"mrg_hist.data", header="# %s"%str(hmrg.variable))

    # mean square end-to-end vector autocorrelation
    _acf(time, eevectors, reader.dir+os.sep+"ree_acf.data")

    if persistence:
        hlp.write(reader.dir+os.sep+"lp_hist.data", header="# %s"%str(hlp.variable))
        hbndl.write(reader.dir+os.sep+"bndl_hist.data", header="# %s"%str(hbndl.variable))

        # flory flory definition (bond dependent)
        f=open(reader.dir+os.sep+"lp_bnd.data",'w')
        f.write("# mean     std       min       max\n")
        for k in range(nbndch_):
            v = _hlp[k].variable
            f.write("%d %g %g %g %g\n"%(k, v.mean(), v.std(), v.min(), v.max()))
        f.close()

        # get the mean, normalize, and write down bndacf
        bndacf = np.array(bndacf).mean(axis=0).transpose()
        f=open(reader.dir+os.sep+"lp_acf.data",'w')
        for x, y in enumerate(bndacf): f.write("%f.1 %g\n"%(x*hbndl.variable.mean(), y))
        f.close()

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
    mrg_hist.data : mean radius of gyration distribution 
If the "--lp" option is enabled:
    bndl_hist.data : bond length distribution
    lp_hist.data : distribution of the projection length of the end-to-end vector in the direction
                   of the first bond (according to Flory's definition of the persistence length) 
    lp_bnd.data : presistance length calculated using Flory's definition for each of the backbone bonds
    lp_acf.data : backbone bond autocorrelation (--lp enabled) ''')

    # add arguments (self explaned)
    string = 'the path to the simulation file. In the case of gromacs' +\
             'simulation, a topology file should be present in the same' +\
             'directory (preferably a tpr file). In the case of lammps the' +\
             'data and dump files will be traced from the corresponding' +\
             'log records, otherwise a data and a dump file with the same' +\
             'base name as the log file should exist in the same directory.'

    parser.add_argument('path', default="."+os.sep,  \
                       help=string)

    parser.add_argument('-start', nargs=1, type=int, metavar='START', default=[-1], \
                       help='start processing form step START [inclusive]')

    parser.add_argument('-end', nargs=1, type=int, metavar='END', default=[sys.maxsize], \
                       help='stop processing at step END [inclusive]')

    parser.add_argument('-every', nargs=1, type=int, metavar='EVERY', default=[1], \
                       help='process every EVERY frames (process frequency)')
    
    argdttype = utils.chk_number("wrong integration time step",numbertype=float, positive=True)
    parser.add_argument('-dt', nargs=1, type=argdttype, default=[0.0], metavar='timestep', \
                       help='integration time step in ps')

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

    parser.add_argument('-molecules', nargs=1, type=utils.argparse_moleculestype, default=[[]],  metavar='range', \
                       help='molecules to be used. A list with comma seperated id ranges should be provided e.g. "1,2,3" or "1:10,20,30:100"')

    parser.add_argument('--no-whole', dest='whole', default=True, action='store_false', \
                       help="do not make molecules whole sine they are already in the trajectories.")

    chktype = utils.IsListOfNamedList("wrong -bins argument (check: %s)", itemtype=float,
        positive=True, llen=1, choices=("sqree","msqree","ree","mree","sqrg","msqrg","rg","mrg","pl","bnd"))
    string = '''
    provide the bins length for the histograms to be calculated. The following
    histograms are supported:
      - sqree : square end-to-end vector
      - msqree : mean square end-to-end
      - ree : end-to-end length
      - mree : mean end-to-end length
      - sqrg : square radius of gyration
      - msqrg : mean square radius of gyration
      - rg : radius of gyration
      - mrg : mean radius of gyration
    For example, with argument "ree:2.5@sqree:50", the bin length of the
    histogram for the end-to-end distance is set to 2.0 Å and for the square
    end-to-end distance to 50 Å².
    '''
    parser.add_argument('-bins', nargs=1, type=chktype, default=[{}], \
        metavar='<list of porperties histograms>', help=string)

    string = '''
    calculate the persistence length assuming a linear polybead topology with a single
    type of bead where end1 is the first atom and end2 is the last atom of the sequence.
    If this option is enalbed files bndl_hist.data, lp_hist.data, lp_bnd.data, and 
    lp_acf.data will be created in the simulation directory (see epilogue).
    '''
    parser.add_argument('--lp', dest='lp', default=False, action='store_true', \
                       help=string)

    parser.add_argument('--camc', dest='camc', default=False, action='store_true', \
                       help="process connectivity monte carlo output")

    # parse the arguments
    args = parser.parse_args()

    print("INPUT")
    print("path      : ", args.path)
    print("dt (ps)   : ", args.dt[0])
    print("start     : ", args.start[0])
    print("end       : ", "max" if args.end[0] == sys.maxsize else args.end[0])
    print("every : ", args.every[0])
    string = " ".join( map(str, args.molecules[0]))
    print("molecules : %s" % ('-' if len(string)==0 else string ) )
    print("whole     : %s" % ("True" if args.whole else "False"))
    print("persistence length : %s" % ("True" if args.camc else "False"))
    print("camc      : %s" % ("True" if args.camc else "False"))
    print()
    if __debug:
        print(args.molecules[0])

    endtoend(
        args.path,
        args.end1[0],
        args.end2[0],
        dt=args.dt[0],
        start=args.start[0],
        end=args.end[0],
        every=args.every[0],
        molids=args.molecules[0],
        whole=args.whole,
        bins=args.bins[0],
        persistence=args.lp,
        camc=args.camc)

if __name__ == '__main__':
    command()
