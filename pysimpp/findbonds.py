# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
from collections import defaultdict

import MDAnalysis

import pysimpp.readers
from pysimpp.fastpost import fastfindbonds  # pylint: disable=no-name-in-module
from pysimpp.utils.statisticsutils import Histogram
from pysimpp.utils.utils import parse_radii, isnumber


def _is_command(): return True
def _short_description(
): return 'Find the bonds and perfom analysis (lengths distribution, time evolution or the mean value, etc.).'
def _command(): command()

# import scipy
# from scipy import interpolate
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt


_lmp_data_bond_local_header = '''\
ITEM: TIMESTEP
%d
ITEM: NUMBER OF ENTRIES
%d
ITEM: BOX BOUNDS pp pp ff
%s
ITEM: ENTRIES id type bead1 bead2
'''


_el_vdwradii = dict(MDAnalysis.topology.tables.vdwradii)


def findbonds(filename, fradii, scale, pairs, bin, start, end, every, dump=False):

    # check the input file
    reader = pysimpp.readers.create(filename)
    if not reader:
        sys.exit(0)

    dirname = reader.dir
    basename = reader.basename
    natoms = reader.natoms

    # set what to read
    attributes = 'id type x y z'
    reader.set_attributes(attributes)

    # find critical radii for the pairs of interest. if no radii file is
    # provided asume elements based analysis and retrieve the default vdw
    # radii from MDAnalysis.
    radii = parse_radii(fradii) if not fradii is None else {
        "elements": _el_vdwradii}

    # retrieve element or type based identifier
    types = reader.get_atom_element() if "elements" in radii else reader.get_atom_type()
    # map unique identifiers to id
    typesmap = {t: i for i, t in enumerate(set(types))}
    # map id to unique identifiers
    typesmaprev = {v: k for k, v in typesmap.items()}
    ntypes = len(typesmap)  # the number of the unique identifiers

    # number based atom id (element or type based)
    typesid = np.array([typesmap[t] for t in types])

    # setup the possible bonded pairs and their critical distances
    pairmap = {}
    pairid = np.zeros((ntypes, ntypes), dtype=np.int32)
    pairrc = np.zeros((ntypes, ntypes), dtype=np.float64)
    _radii = radii["elements"] if "elements" in radii else radii["types"]
    # Check: Inorganic Materials, Vol. 37, No. 9, 2001, pp. 871–885
    if not "CO" in _radii:
        _radii["CO"] = 2.30
    if not "TI" in _radii:
        _radii["TI"] = 2.16
    npairs = 0
    for i in range(ntypes):
        ti = typesmaprev[i]
        ri = _radii[ti]
        for j in range(i, ntypes):
            # check if the pair has explicitly requested/set assuming that
            # pairs are lexicographically ordered
            tj = typesmaprev[j]
            _pair = (ti, tj)
            pairrc[i, j] = pairs[_pair] if _pair in pairs else \
                scale * (ri + _radii[tj])
            if pairrc[i, j] > 0.0:
                pairid[i, j] = npairs
                pairmap[npairs] = "%s-%s" % _pair
                npairs += 1
            else:
                pairid[i, j] = -1

            if not i == j:
                pairid[j, i] = pairid[i, j]
                pairrc[j, i] = pairrc[i, j]

    # estimate the max number of bonds per frame
    maxnbd = natoms * 4

    # keep bond time evolution per type; the default zero tuples facilitate
    # the printout since not all the bond types are available at each frame.
    trj = defaultdict(lambda: defaultdict(lambda: (0, 0., 0.)))
    # keep bonds distribution per type
    hst = defaultdict(lambda: Histogram.free(bin, 0.0, addref=False))

    if dump:
        fd = open(dirname+os.sep+basename+".bonds.dump_local", 'w')

    frames = []
    r = np.empty(shape=(natoms, 3), dtype=np.float32)  # coordinates
    print('>> reading dump file(s) ...')
    while (True):
        step, box, data = reader.read_next_frame()

        if step is None:
            break

        if step < start:
            continue
        elif step > end:
            break

        frames.append(step)

        if not step % every == 0:
            continue

        np.copyto(r[:, 0], data['x'])
        np.copyto(r[:, 1], data['y'])
        np.copyto(r[:, 2], data['z'])

        nbd, bdid, bd, bdln, ierr = fastfindbonds(
            typesid, r.T, box.origin, box.va, box.vb, box.vc, pairid, pairrc, maxnbd)

        for _t in set(bdid):
            if not _t == -1:
                _b = bdln[bdid == _t]
                np.vectorize(hst[_t].add)(_b)
                trj[_t][step] = (len(_b), np.mean(_b), np.std(_b))

        if dump:
            fd.write(_lmp_data_bond_local_header % (step, nbd, str(box)))
            for _id, (_t, _at1, _at2) in enumerate(zip(bdid[:nbd], bd[0, :nbd], bd[1, :nbd])):
                fd.write("%d %d %d %d\n" % (_id+1, _t, _at1+1, _at2+1))

    print()

    if dump:  # close bonds local dump
        fd.close()

    # dump time evolution
    _types = sorted(trj.keys())
    header = "# %-12s" % 'frame' + \
        " ".join(["nbd_t%-7s mean_t%-6s std_t%-7s" %
                 ((pairmap[_t],)*3) for _t in _types])
    f = open(dirname+os.sep+basename+"_bonds.dat", 'w')
    lines = []
    for step in frames:
        line = "%-12s" % str(step)
        for _t in _types:
            line += " %-12d %-12g %-12g" % trj[_t][step]
        lines.append(line)
    f.write(header + '\n')
    for line in lines:
        f.write(line + '\n')
    f.close()

    # dump distributions
    for k in list(hst.keys()):
        hst[k].write(dirname+os.sep+basename+"_b%s.hst" % pairmap[k])


def command():

    # import sys
    # sys.argv += "--dump -scale 0.55 -bin 0.02 -rc A:B:2.0@A:C:2.1 /Users/loukas.peristeras/tmp/asma/sys1.dump".split()

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(description=_short_description())

    # add arguments (self explaned)
    string = 'the path to the simulation trajectory file. A topology file' + \
             'should be present in the same directory (preferably a tpr file).'
    parser.add_argument('path', default="."+os.sep, help=string)

    parser.add_argument('-start', nargs=1, type=int, metavar='n', default=[-1],
                        help='start processing form configuration n [inclusive]')

    parser.add_argument('-end', nargs=1, type=int, metavar='n', default=[sys.maxsize],
                        help='stop processing at configuration n [inclusive]')

    parser.add_argument('-every', nargs=1, type=int, metavar='n', default=[1],
                        help='processing frequency (every n configurations)')

    parser.add_argument('-bin', nargs=1, type=float, metavar='bin', default=[0.02],
                        help='bin length (Å) for the distributions')

    string = '''
    the file with the radii of the atoms. It can be element or type based. 
    The first line of the file contains the keywords "(elements|types) (r|d)"; 
    the former, specifies the atom type identifier and the latter if the 
    radius (r) or the diameter (d) is given for each type. The rest of the 
    lines contain the (type, radius) pairs. The type could be either a 
    number (type id) or a string (type name). If no radii file is provided,
    the analysis will be performed element-based with the vdw radii retrieved
    from the MDAnalysis tables. '''
    parser.add_argument('-radii', nargs=1, type=argparse.FileType('r'), required=False, default=[None],
                        metavar='<file with atoms\' type/element radii>', help=string)

    parser.add_argument('-scale', nargs=1, type=float, metavar='scale', default=[0.55],
                        help='scale the critical distance calculated based on the vdw radii.')

    string = '''
    set the pairs and their critical distance for bonds calculations. The argument is
    a set of column-separated records per pair, separated with "@". The first two record
    are the constituents elements/types of the pair and the third is the critical 
    distance (Å). For example the argument "AT1:AT1:rc1@AT1:AT3:rc2" provides the  
    critical distance rc1 and rc2 for the pairs AT1-AT1 and AT1-AT3, respectively. 
    For the definition of the pairs, either the types or the elements of the constituent
     particles should be provided; this should conform with the the format of the radii
    file its provided with "-radii" option. If 0.0 is provided as critical distance for
    a specific pair, the element/type based radii will be used for it specification, 
    while if the value is negative the pair will be ignored.
    '''

    def chktype(string):
        ''' check the -rc option argument. '''
        ret = {}
        lines = string.split("@")
        for line in lines:
            tk = line.split(":")
            if not len(tk) == 3 or isnumber(tk[2], float) is None:
                msg = "wrong pair's critical distance format (check: %s)" % line
                raise argparse.ArgumentTypeError(msg)
            # order the pairs
            ret[(tk[0], tk[1]) if tk[0] < tk[1]
                else (tk[1], tk[0])] = float(tk[2])
        return ret

    parser.add_argument('-rc', nargs=1, type=chktype, metavar='rc', default=[{}],
                        help=string)

    string = '''
    output a lammps local dump trajectory for bonds. The file will be located at
    the path of the input trajectory, name using the trajectory dump base name
    with extension ".bonds.dump_local".
    '''
    parser.add_argument('--dump', dest='dump', default=False, action='store_true',
                        help=string)

    # parse the arguments
    args = parser.parse_args()

    print("INPUT")
    # print("path : ", args.path)
    print("start : ", args.start[0])
    print("every : ", args.every[0])
    print("end : ", "max" if args.end[0] == sys.maxsize else args.end[0])
    print("bin : ", args.bin[0])
    print("radii : ", "-" if args.radii[0] is None else args.radii[0].name)
    print("scale : ", args.scale[0])
    print("dump : %s" % ("True" if args.dump else "False"))
    # print("rcrt : ", ",".join([ "%s:%f.4"%(k,v) for k v in args.rc[0].items()]))

    findbonds(args.path, args.radii[0], args.scale[0], args.rc[0],
              args.bin[0], args.start[0], args.end[0], args.every[0],
              dump=args.dump)


if __name__ == '__main__':
    command()
