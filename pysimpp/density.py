#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.constants import N_A as _NA  # in 1/mol
from scipy.constants import k as _kb  # in J/K
import os
import sys
from collections import defaultdict
from operator import add
import pickle

import numpy as np
from scipy import stats

import pysimpp.readers
from pysimpp.utils.utils import IsList, IsListOfList, chk_number, argparse_moleculestype
from pysimpp.utils.statisticsutils import Binning, Histogram
from pysimpp.fastpost import fastcom, fastcom_total, fast_localdensity


__debug = False


def _is_command(): return True


def _short_description(
): return 'Calculate probability density profiles along x, y, and/or z axes.'


def _command(): command()


def _dimsbookkeep(dimensions, unwrap=True):
    ''' Given the dimensions i.e. a list of strings, creates the attributes
        string for the lammpsreader and the dimensions map for data storage. '''

    # set what to read and store bookkeeping
    attributes = 'id type'
    xstring = ""
    istring = ""
    dmap = {}
    dmapinv = {'x': 0, 'y': 1, 'z': 2}
    dsurf = {}
    for dim in ['x', 'y', 'z']:
        if dim in dimensions:
            xstring += " %s" % dim
            idim = "i%s" % dim
            istring += " %s" % idim
            idim = dmapinv[dim]
            dmap[idim] = (dim, idim)
            dsurf[idim] = list(set((0, 1, 2))-set((idim,)))
    attributes = attributes + xstring
    if unwrap:
        attributes += istring
    return attributes, dmap, dsurf, dmapinv


class DensityBinning(Histogram):

    def __init__(self):
        self.h = None  # current histogram
        self.vH = None  # overall values histogram
        self.nH = None  # overall numbers histogram
        self.d = 1.0

    @classmethod
    def fixed(cls, min, max, d):
        obj = cls()
        obj.h = Binning.fixed(min, max, d)
        obj.vH = Binning.fixed(min, max, d)
        obj.nH = Binning.fixed(min, max, d)
        obj.d = d
        return obj

    def add(self, b, value=1.0, m=1):
        ''' Add a value in the histogram at the specific bin with multiplicity m. '''
        self.h.add(b, value, m)

    def write(self, filename, header=""):
        ''' Write histogram data to the file filename. '''
        x, n, nstd = self.nH.data()
        _, v, vstd = self.vH.data()
        # yn = n / np.sum(n) / self.d  # normalize
        yv = v / np.sum(v) / self.d  # normalize
        f = open(filename, 'w')
        if len(header) > 0:
            f.write("%s\n" % header)
        f.write("# %s\n" % " ".join(["%12s" % s for s in (
            "z[Å]", "n", "n_std", "rho[g/cm³]", "rho_std", "p(rho)[1/Å]")]))
        for x_, n_, nstd_, v_, vstd_, yv_ in zip(x, n, nstd, v, vstd, yv):
            f.write("%12g %12g %12g %12g %12g %12g\n" %
                    (x_, n_, nstd_, v_, vstd_, yv_))

    def append(self, surf=1.0):
        ''' Append the current values to the overall histogram. '''
        vol = surf * self.d * (_NA * 1.e-24)
        for n, v, yc in np.c_[self.nH.h, self.vH.h, self.h.h]:
            if yc._n > 0.0:
                v.set(yc._sum / vol)
                n.set(yc._n)
        self.h.reset()

# @profile
def density(filename, 
            bin, 
            temp, 
            start=-1, 
            end=sys.maxsize, 
            every=1, 
            dimensions=['z'], 
            molids=(), 
            types=(), 
            local=(), 
            dowrap=True, 
            usecom=True):

    # create the reader
    reader = pysimpp.readers.create(filename)
    if not reader:
        print("ERROR: it was not possible to parse the given file.")
        return

    # remove center of mass flag
    rmcom = False
    # access wrapped coordinates
    reader.set_wrap(dowrap)
    _unwrap = not dowrap

    attributes = 'id mol type x y z'
    reader.set_attributes(attributes)

    print('>> reading data file ...')
    reader.read_topology()

    # TODO create a utility method in pysimpp.readers for accessing molecular data
    natoms = reader.natoms
    # get molecular data and select molecules of interest
    masses = reader.get_atom_mass()
    molecule = reader.get_atom_molecule() - 1  # atom molecule array (index @zero)
    nmolecules = molecule.max() + 1            # number of molecules
    # TODO add this functionality in reader
    mol_atoms = defaultdict(list)             # molecule atoms array
    mol_mass = np.zeros(nmolecules, dtype=np.float32)
    for i, im in enumerate(molecule):
        mol_atoms[im].append(i)
        mol_mass[im] += masses[i]

    # selected molecules
    selected = np.sort(molids) - 1  # molecule based selection (index @zero)
    if not usecom:  # convert to atoms
        on = np.zeros(natoms, dtype=np.bool)
        for im in selected:
            on[mol_atoms[im]] = True
        selected = np.where(on)[0]
    nselected = selected.size  # number of selected molecules
    hasselected = nselected > 0  # selected flag

    # select types
    hastypes = False  # atom type bases select flag
    selected_types = {}  # type based selections (index @zero)
    # assuming that atom type does not change along the trajectory (constant
    # connectivity), type based selection arrays can be calculated now
    if not usecom:  # convert to atoms
        hastypes = len(types) > 0  # type bases select flag
        if hastypes:
            atom_type = reader.get_atom_type()
            if hasselected:
                on = np.zeros(natoms, dtype=np.bool)
                on[selected] = True
                for t in types:
                    selected_types[t] = np.where(on & (atom_type == t))[0]
            else:
                for t in types:
                    selected_types[t] = np.where(atom_type == t)[0]

    # bookkeeping dimensions to read
    attributes, dmap, _, _ = _dimsbookkeep('xyz', _unwrap)

    _dims = tuple([['x', 'y', 'z'].index(_d) for _d in dimensions])
    _dnames = {['x', 'y', 'z'].index(_d): _d for _d in dimensions}
    reader.set_attributes(attributes)

    # allocate wrapped and unwrapped molecules coordinates
    n_ = nselected if hasselected else nmolecules if usecom else natoms
    cm = np.zeros((n_, 3), dtype=np.float32)  # com coordinates
    mass = np.zeros(n_, dtype=np.float32)  # com coordinates
    r = np.empty(shape=(natoms, 3), dtype=np.float32)  # unwrapped coordinates
    # initial system com
    cm0 = np.zeros((3), dtype=np.float32)
    # vector for com removal
    delta = np.zeros((3), dtype=np.float32)

    print("\n>> reading dump file(s) ...")
    steps = []
    boxes = []
    profiles = {}
    profiles_type = {}

    dolocal = False
    if len(local) > 0:
        dolocal = True
        seed = np.array((43924342, 77928374, 46278346, 45329834, 6432847, 23984, 9384392, 234198), dtype=np.int32, )
        hld = defaultdict(lambda: Histogram.free(4.0, 0.0, addref=False))
        l, nl = zip(*local)

    iframe = -1
    while True:
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

        if len(profiles) == 0:
            for _d in _dims:
                origin = box.origin[_d]
                length = origin + (box.a, box.b, box.c)[_d]
                profiles[_d] = DensityBinning.fixed(origin, length, bin)
                profiles_type[_d] = {}
                for t in selected_types:
                    profiles_type[_d][t] = DensityBinning.fixed(origin, length, bin)

        steps.append(step)
        boxes.append(box)

        for k, v in dmap.items():
            np.copyto(r[:, k], data[v[0]])

        if rmcom:
            cmtotal, masstotal = fastcom_total(r, masses)
            if len(boxes) == 1:
                cm0[:] = cmtotal
            else:
                delta[:] = cmtotal - cm0
                r -= delta

        # calculate center of masses and if regions
        # set also the wrapped coordinates to ensure
        # the correct binning
        if hasselected:
            cm_, masses_ = (
                (fastcom(r, masses, molecule, nmolecules), mol_mass) if usecom else (r, masses))
            cm[:, :] = cm_[selected, :]
            mass[:] = masses_[selected]
        else:
            cm[:, :], mass = (
                (fastcom(r, masses, molecule, nmolecules), mol_mass) if usecom else (r, masses))

        for _d in _dims:
            profile = profiles[_d]
            np.vectorize(profile.add)(cm[:, _d], mass)
            surf = (box.bxc_, box.cxa_, box.axb_)[_d]
            profile.append(surf)

            profile_ = profiles_type[_d]
            for t in profile_.keys():
                profile_t = profile_[t]
                s = selected_types[t]
                np.vectorize(profile_t.add)(r[s, _d], masses[s])
                profile_t.append(surf)

        if dolocal:
            for l_, nl_ in zip(l, nl):
                ld = fast_localdensity( cm, box.origin, box.a, box.b, box.c, seed, l_, nl_)
                np.vectorize(hld[l_].add)(ld)

    for _d, profile in profiles.items():
        _dname = "%s%s" % ("com_" if usecom else "", _dnames[_d])
        header = "# %s_[Å] probability_[1/Å]" % (_dname,)
        profile.write(reader.dir + os.sep + "%s_prf.data" % _dname, header=header)

    for _d, profile in profiles_type.items():
        _dname = _dnames[_d]
        for t in profile.keys():
            header = "# number and density profiles along %s direction" % (_dname,)
            profile[t].write(reader.dir + os.sep + "%s_prf_%s.data" % (_dname, t), header=header)

    print()

    # print local density and stuff
    if dolocal:
        f = open(reader.dir + os.sep + "ld.data", "w")
        header = (
            "box[nm]         "
            + "mean[#beads]    "
            + "std[#beads]     "
            + "compressibility[1/bar]"
        )
        f.write("# %s\n" % header)
        for i, l_ in enumerate(l):
            h = hld[l_]
            h.write(reader.dir + os.sep + "ld_%d_hist.data" %
                    i, header="# %s l = %g vol = %g" % (str(h.variable), l_, l_ * l_ * l_), )
            mean = h.variable.mean()
            std = h.variable.std()
            vol = l_**3 * 1.0e-30
            k = std**2 * vol / _kb / temp / mean**2 * 1.0e5
            f.write(" %-16g %-16g %-16g %-16g\n" % (l_ / 10.0, mean, std, k))
        f.close()

    print()


def command():
    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(
        description="Calculate probability density profiles along x, y, and/or z axes.")

    # add arguments (self explaned)
    string = (
        "the path to the simulation file. In the case of gromacs"
        + "simulation, a topology file should be present in the same"
        + "directory (preferably a tpr file). In the case of lammps the"
        + "data and dump files will be traced from the corresponding"
        + "log records, otherwise a data and a dump file with the same"
        + "base name as the log file should exist in the same directory."
    )

    parser.add_argument("path", default="." + os.sep,
                        help=string)

    parser.add_argument("-dim", nargs=1, choices=["x", "y", "z", "xy", "xz", "yz", "xyz"],
                        metavar="dim", required=True,
                        help="produce the density for the given dimension(s)")

    parser.add_argument("-bin", nargs=1, type=float, metavar="bin", required=True,
                        help="bin width in Å")

    parser.add_argument("-start", nargs=1, type=int, metavar="START", default=[-1],
                        help="start processing form step START [inclusive]")

    parser.add_argument("-end", nargs=1, type=int, metavar="END", default=[sys.maxsize],
                        help="stop processing at step END [inclusive]")

    parser.add_argument("-every", nargs=1, type=int, metavar="EVERY", default=[1],
                        help="process every EVERY frames (process frequency)")

    parser.add_argument("-molecules", nargs=1, type=argparse_moleculestype, default=[[]], metavar="range",
                        help='molecules to be used. A list with comma separated id ranges should be provided e.g. "1,2,3" or "1:10,20,30:100"', )

    def arglocal(string):
        """check the "-local" option arguments."""
        chktype = IsListOfList("wrong local density argument (check: %s)")
        lst_ = chktype(string)
        lst = []
        for l_ in lst_:
            try:
                lst.append([float(l_[0]), int(l_[1])])
            except:
                msg = "wrong local density argument (check: %s)" % l_.join(":")
                raise argparse.ArgumentTypeError(msg)
        return lst

    string = """
    The length of the cubic box (in Å) and the number of trials to be used 
    for calculating of the local density. A list of pairs (length:trials)
    should be provided e.g. "25.0,2000:50.,1000".
    """
    parser.add_argument("-local", nargs=1, type=arglocal, metavar="<local density>", default=[()],
                        help=string)

    parser.add_argument("--no-wrap", dest="dowrap", default=True, action="store_false",
                        help="do not wrap molecules sine the provided trajectories provide wrapped coordinates")

    parser.add_argument("--no-com", dest="usecom", default=True, action="store_false",
                        help="calculate the distributions using atoms rather that molecules' centers of mass")

    chktype = IsList("wrong type names (check: %s)", itemtype=str)
    string = "a comma separated list with the atom type names to be considered."
    parser.add_argument("-types", nargs=1, type=chktype,
                        default=[()], help=string)

    argttype = chk_number("wrong temperature", numbertype=float, positive=True)
    parser.add_argument("-t", nargs=1, type=argttype, metavar="temperature",
                        help="the temperature in K")
    # parse the arguments
    args = parser.parse_args()

    # temperature is required only in local density calculations
    if len(args.local[0]) > 0:
        if args.t is None:
            parser.error(
                "The temperature is required (see the -t argument) to calculate the compressibility through the local density calculation (set with -local option)")
        else:
            temp = args.t[0]
    else:
        temp = 0.0

    print("INPUT")
    print("path  : ", args.path)
    print("start : ", args.start[0])
    print("end   : ", args.end[0])
    print("every : ", args.every[0])
    print("dim   : ", args.dim[0])
    print("temperature (K) : ", temp)
    print("local : ", args.local[0])
    print("wrap  : %s" % ("True" if args.dowrap else "False"))
    print("com   : %s" % ("True" if args.usecom else "False"))
    print("types : ", ",".join(args.types[0])
          if len(args.types[0]) > 0 else "-")

    print("molecules : %d \n" % len(args.molecules[0]))
    if __debug:
        print(args.molecules[0])
    print()

    _dims = [_d for _d in args.dim[0]]
    density(args.path, args.bin[0], 
            temp, start=args.start[0], 
            end=args.end[0], 
            every=args.every[0],
            dimensions=_dims, 
            molids=args.molecules[0], 
            types=args.types[0], 
            local=args.local[0], 
            dowrap=args.dowrap, 
            usecom=args.usecom, )


if __name__ == "__main__":
    command()
