#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np

import pysimpp.utils.utils as utils
from pysimpp.utils.statisticsutils import Histogram, Histogram3D
import pysimpp.readers
from pysimpp.fastpost import fastunwrapv, gyration, pefastunwrap

__debug = False

def _is_command(): return True

def _short_description(): return 'Calculate gyration tensor and stuff.'

def _command(): command()

def shape(filename, 
    start=-1, 
    end=sys.maxsize, 
    molids=(), 
    unwrap=True,
    camc=False):

    # create the reader. for the moment lammps, gromacs and camc (polybead)
    reader = pysimpp.readers.create(filename)
    if not reader: 
        print("ERROR: it was not possible to parse the given file.")
        return

    reader.set_unwrap( unwrap)
    attributes = 'id mol x y z'
    reader.set_attributes( attributes)

    print('>> reading data file ...')
    reader.read_topology()

    # get system  data
    natoms = reader.natoms
    masses = reader.get_atom_mass()
    molecules = reader.get_atom_molecule() - 1  # atom molecule array (index @zero)
    nmolecules = molecules.max() + 1            # number of molecules

    # selected molecules
    selected = np.sort( molids) - 1             # selected molecules (index @zero)
    nselected = selected.size                   # number of selected molecules
    hasselected = nselected > 0                 # selected flag

    r = np.zeros( (natoms, 3), dtype=np.float64)    # wrapped
    rp = np.zeros( (natoms, 3), dtype=np.float64)   # principle frame coordinates
    exclude = np.zeros( (nmolecules), dtype=np.bool) # all false
    rg = np.zeros( (nmolecules,6), dtype=np.float64)        # gyration tensors
    eigval = np.zeros( (nmolecules,3), dtype=np.float64)    # gyration tensors eigenvalues
    eigvec = np.zeros( (nmolecules,9), dtype=np.float64)    # gyration tensors eigenvectors

    if hasselected:
        exclude[:] = True
        exclude[ selected ] = False
        atomselected = [ molecules[iat] in selected for iat in range(natoms)]
    else:
        selected = atomselected = ()

    print('\n>> reading dump file(s) ...')
    steps = []
    boxes = []
    iconf = -1

    # gyration tensor relevant histograms
    hsqrg = Histogram.free(0.25, 0.0, addref=False)
    hasph = Histogram.free(0.25, 0.0, addref=False)
    hacyl = Histogram.free(0.25, 0.0, addref=False)
    hanis = Histogram.free(0.25, 0.0, addref=False)
    hmsqrg = Histogram.free(0.25, 0.0, addref=False)
    hmasph = Histogram.free(0.25, 0.0, addref=False)
    hmacyl = Histogram.free(0.25, 0.0, addref=False)
    hmanis = Histogram.free(0.25, 0.0, addref=False)

    # 3D auto adjusted histogram
    hist = Histogram3D((1.,1.,1.),(0.,0.,0.),addref=False)
    # gyration tensor and eigen values output file
    f1 = open(reader.dir+os.sep+"gyration.data", 'w')
    f1.write( "#{:^14s} {:^14s} {:^14s} {:^14s} {:^12s} {:^12s} {:^12s} {:^12s} {:^12s}\n".format("rg_xx","rg_yy","rg_zz","rg_xy","rg_xz","rg_yz","ex","ey","ez") )
    f2 = open(reader.dir+os.sep+"descriptors.data", 'w')
    f2.write( "#{:^15s} {:^15s} {:^15s} {:^15s}\n".format("sqrg","asphericity","acylindricity","anisotropy") )

    while( True):
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
            steps.append( step)
            boxes.append( box)

            # update the topology as needed in the case of camc
            if camc:
                natch, chat, e1, e2 = reader.get_topology()
                molecules = reader.get_atom_molecule() - 1
                masses = reader.get_atom_mass()
                atomselected = [ molecules[iat] in selected for iat in range(natoms)] if hasselected else ()
            else:
                selected = ()
                atomselected = ()

            r[:,0] = data['x']
            r[:,1] = data['y']
            r[:,2] = data['z']
        else:
            break

        rp[:,:], rg[:,:], eigval[:,:], eigvec[:,:], ierr = gyration(r, masses, molecules, exclude)
        sqrg = []
        asph = []
        acyl = []
        anis = []
        for e_ in eigval[selected]: # LDP TODO: check this out
           # sqrg_ = (e_*e_).sum()
           sqrg_ = (e_).sum()
           sqrg.append( sqrg_)
           b_ = e_[0]-0.5*(e_[1]+e_[2])
           asph.append( b_)
           c_ = e_[1]-e_[2]
           acyl.append( c_)
           sqk_ = (b_*b_+0.75*c_*c_)/(sqrg_*sqrg_)
           anis.append( sqk_)
           hsqrg.add( sqrg_)
           hasph.add( b_)
           hacyl.add( c_)
           hanis.add( sqk_)
        msqrg = np.mean(sqrg)
        masph = np.mean(asph)
        macyl = np.mean(acyl)
        manis = np.mean(anis)
        hmsqrg.add( msqrg)
        hmasph.add( masph)
        hmacyl.add( macyl)
        hmanis.add( manis)

        for v, m in zip(rp[atomselected], masses[atomselected]):
            hist.add(v,m)
        f1.write("#step %d\n"%step)
        for v, e in zip(rg, eigval):
            f1.write( "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n" % (v[0],v[1],v[2],v[3],v[4],v[5],e[0],e[1],e[2]) )
        f2.write( "%15.6f %15.6f %15.6f %15.6f\n"%(msqrg,masph,macyl,manis))
    f1.close()
    f2.close()

    # number of configurations - resize if needed
    nconfs = len( steps)

    # write down the distributions
    hsqrg.write(reader.dir+os.sep+"gtsqrt_hist")
    hasph.write(reader.dir+os.sep+"gtasph_hist")
    hacyl.write(reader.dir+os.sep+"gtacyl_hist")
    hanis.write(reader.dir+os.sep+"gtanis_hist")

    hmsqrg.write(reader.dir+os.sep+"gtmsqrt_hist")
    hmasph.write(reader.dir+os.sep+"gtmasph_hist")
    hmacyl.write(reader.dir+os.sep+"gtmacyl_hist")
    hmanis.write(reader.dir+os.sep+"gtmanis_hist")

    hist.normalize()
    hist.write(reader.dir+os.sep+"dprof.cube",format='cube')

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(
        description='Calculate gyration tensor and shape.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=''' \
The output files are (located at the simulation directory):
    gtsqrt_hist : square rarious of gyration distrubtion
    gtasph_hist : asphericity distrubtion
    gtacyl_hist : acylindricity distrubtion
    gtanis_hist : anisotropy distrubtion
    gtmsqrt_hist : mean square rarious of gyration
    gtmasph_hist : mean asphericity distrubtion
    gtmacyl_hist : mean acylindricity distrubtion
    gtmanis_hist : mean anisotropy distrubtion
    dprof.cube : mass density spatian distribution(3d) ''')

    # add arguments (self explaned)
    string = 'the path to the simulation log file. In the case of gromacs simulation, a topology file' + \
             'should be present in the same directory (preferably a tpr file).'
    parser.add_argument('path', default="."+os.sep,  \
                       help=string)
    parser.add_argument('-start', nargs=1, type=int, metavar='n', default=[-1], \
                       help='start processing form configuration START [inclusive]')
    parser.add_argument('-end', nargs=1, type=int, metavar='n', default=[sys.maxsize], \
                       help='stop processing at configuration END [inclusive]')

    def argmoleculestype( string):
        ''' check the "-molecules" option arguments. '''
        if len( string) == 0:
            return []
        numbers = utils.isrange( string, positive=True)
        if len( numbers) == 0:
            msg = "wrong molecules indexs range (check: %s)" % string
            raise argparse.ArgumentTypeError(msg)
        return numbers
    parser.add_argument('-molecules', nargs=1, type=argmoleculestype, default=[[]],  metavar='range', \
                       help='molecules to be used. A list with comma seperated id ranges should be provided e.g. "1,2,3" or "1:10,20,30:100"')

    parser.add_argument('--no-unwrap', dest='unwrap', default=True, action='store_false', \
                       help="do not unwrap molecules sine the provided trajectories provide unwrapped coordinates.")

    # parser.add_argument('--camc', dest='camc', default=False, action='store_true', \
    #                    help="process connectivity monte carlo output.")

    # parse the arguments
    args = parser.parse_args()

    print("INPUT")
    print("path : ", args.path)
    print("start : ", args.start[0])
    print("end : ", args.end[0])
    print("molecules : %d \n" % len(args.molecules[0]))
    print("unwrap : %s \n" % ("True" if args.camc else "False"))
    # print("camc : %s \n" % ("True" if args.camc else "False"))

    if __debug:
        print(args.molecules[0])

    shape( args.path, 
        start=args.start[0], 
        end=args.end[0], 
        molids=args.molecules[0],
        unwrap=args.unwrap)
        # camc=args.camc)

if __name__ == '__main__':
    command()