# -*- coding: utf-8 -*-

import sys
#sys.setrecursionlimit(5500)
import os
import glob
import re
import inspect
from collections import defaultdict
from io import StringIO
import itertools

import numpy as np
import yaml

from pysimpp.utils.simulationbox import SimulationBox
import pysimpp.utils.statisticsutils as statutils
from pysimpp.fastpost import fastunwrapv, fastwhole # pylint: disable=no-name-in-module

from .reader import abcReader # pylint: disable=import-error

_inform = True
_excludeid = True
_fixcamc = False # Fix: the global indexing is screwed in old camc output

# if not _debug:
#     sys.tracebacklimit=0


class LammpsReaderException(Exception): ...

class LammpsDumpsHandler():
    ''' Implements a handler for the various dumps in a lammps simulation. '''

    def __init__(self, dumps={}):
        ''' Initiaze a LammpsDumpsHandler object.
            Args:
                dumps (dict): {dump name:dump object} pairs of dump names
                    with the corresponding LammpsDump objects.
        '''
        self.dumps = dumps
        # a dictioanary to connect the dump files with the attributes
        # to be retrieved {dump_file:attributes_list}
        self.toread = defaultdict(lambda: set())

    def set_attributes(self, attributes):
        ''' Set the attributes to read from the simulation dump files.
            Args:
                attributes (list): a list of atoms attributes (str)
                    to be retrieved from the dump files.
            Returns:
                bool: True if all all attributes are available in
                    the dump files.
        '''
        requested = set(
            attributes.split() if type(attributes) is str else attributes)
        # try to inject periodic indexes as extra attributes
        extra = False
        _poff = set(('ix','iy','iz'))
        if not _poff.issubset(requested):
            requested = requested.union(_poff)
            extra = True
        supported = self._update_canread( requested)
        # for v in self.dumps.values():
        #     supported = requested.intersection(set(v.attributes))
        #     if len(supported) > 0:
        #         self.toread[v].update(supported)
        #         requested = requested.difference(supported)
        if extra: # do not check for the extra attributes
            requested = requested.difference(_poff)
        # if coordinates are not pressent in the supported attributes
        # the try to see if unwrapped can be added
        self._make_link = {}
        if len( supported.intersection(set(('x','xu')))) == 0 :
            _tmp=(('x','y','z'),('xu','yu','zu'))
            _rm, _add = _tmp if ('x' in requested) else reversed(_tmp)
            requested.difference_update(_rm)
            requested.update(_add)
            supported.update( self._update_canread(requested))
            # for v in self.dumps.values():
            #     supported = requested.intersection(set(v.attributes))
            #     if len(supported) > 0:
            #         self.toread[v].update(supported)
            #         requested = requested.difference(supported)
            self._make_link = { _from:_to for _to, _from in zip(_rm, _add)}
            print("INFO: the requested attribute(s) %s are not present in the dump" % str(_rm))
            print("      file(s); instead they will be mapped into %s ttributes." % str(_add))

        self.canread = len(requested) == 0
        return self.canread

    def _update_canread(self, requested):
        ''' update self.toread based on the requested set of attributes
            and reture a set with the attribute that are supported by
            the registed dump objects. '''
        supported = set()
        for v in self.dumps.values():
            supported.update( requested.intersection(set(v.attributes)))
            if len(supported) > 0:
                self.toread[v].update(supported)
                requested.difference_update(supported)
        return supported

    def get_natoms(self):
        ''' Retrieve the number of atoms and check that all dumps in the
            handler have the same number of atoms.
            Returns:
                int: the number of atoms.
            Raises:
                LammpsReaderException: if the number of atom is not the same in
                    all dump files handles by the reader.
        '''
        natoms = 0
        if not len(self.dumps) == 0:
            for v in self.dumps.values():
                if natoms == 0:
                    natoms = v._natoms
                elif not v._natoms == natoms:
                    raise LammpsReaderException(
                        "different number of atoms in hundlers dumps")
        return natoms

    def count_frames(self, start=-1, end=sys.maxsize):
        ''' Count the number of configurations available in the dump
            files of the handler. Should be called after setting the
            attributes.
            Returns:
                int: the number of configurations.
        '''
        # just pick the first dump and count. '''
        nframes = 0
        if len(self.toread) > 0:
            _dump = list(self.toread.keys())[0]
            nframes = _dump.count_frames(start, end)
        return nframes

    def read_next(self, skip):
        ''' Read next frame available in the dump(s).
            Returns:
                int, SimulationBox, {}: return the frame step, the
                    frame SimulationBox object and a dict with the
                    data (value) for each attribute (key). If there
                    is no next frame to read the method returns the
                    (None, None, {}) tuple.
            Raises:
                LammpsReaderException: if wrong step sequence is traced.
        '''
        istep = None
        box = None
        data = {}
        if len(self.toread) > 0:
            # assume that steps and boxes are the same in all the dump files
            # and keep the ones from the first dump.
            for k, v in self.toread.items():
                _istep, _box, _data = k.read_next(skip)
                if _istep is None:
                    break
                if istep is None:
                    istep = _istep
                elif not istep == _istep:
                    msg = "wrong step sequence in %s (found step %d while expecting %d)" % (k.file, _istep, istep) # pylint: disable=bad-string-format-type
                    raise LammpsReaderException(msg)
                if box is None:
                    box = _box
                else: # free some memory
                    del _box
                for name in v:
                    if _excludeid and name == 'id': continue
                    data[name] = _data[name]
                del _data
                # rename attributes, i.e., create the links
                # decided in self.set_attributes
                for _from, _to in self._make_link.items():
                    data[_to] = data[_from]

        return (istep, box, data)

class LammpsDump():
    ''' Implements a lammps dump file handler. Operates on a list of
        dump files created by restarting the simulation. In this case,
        the dump from the ith restart should be named as declared in
        the lammps dump command, with the restart number appended. For
        example if the line:
            "dump           1 all atom 1000 c28_60.dump"
        in the lammps input script, will create a LammpsDump object
        handling, if available, the dump(s):
            [c28_60.dump, c28_60.dump.1, c28_60.dump.2, ... ]
    '''

    def __init__(self, line, path):
        ''' Initialize a dump from the corresponding line of
            lammps input script.
            Args:
                line (str): the dump line from the lammps input script.
                    or just the name of the dump file. In the later case
                    the basename of the file will be used also as the dump
                    name.
                path (str): the simulation path where the input script
                    is located.
        '''
        # initialize
        self.f = None               # current file unit
        self.nconfigurations = 0    # number of nconfigurations
        self.n_dump = []            # nconfigurations steps
        self.files = []             # dump files handled for this object
        self._stopnow = False       # force stop in thread excecution

        tk = line.split()            # passed striped
        if len(tk) == 1:
            self.name, _ = os.path.splitext( tk[0])
            self.filename = path + tk[0]
        else:
            self.name = tk[1]            # dump name
            self.filename = path + tk[5] # dump file name (absolute path)

        # check if the file exist
        if not os.path.isfile(self.filename):
            msg = '%s file not found.' % self.filename
            raise LammpsReaderException(msg)
        # if yes append it to the list of files to be handled
        self.files.append(self.filename)

        self._find_files()

        self._scan_first()

        self._next_file = self.__next_file()

    def has_attribute(self, attribute):
        ''' Check if this dump supports the given attribute. '''
        # return ' '+attribute.strip()+' ' in self.attributes
        return attribute.strip() in self.attributes

    def open(self):
        ''' Open the first dumpm file to start reading. '''
        self.file = next(self._next_file)
        if not os.path.isfile( self.file):
            msg = '%s file not found.'%self.file
            raise LammpsReaderException( msg)
        self.f = open( self.file, 'r')
        self.nconfigurations = 0

    def count_frames( self, start=-1, end=sys.maxsize):
        ''' Return the number of configuration in the dump file. '''
        cnt = 0
        nfiles = 0 # remove the first configuration from the second file and on
        regex = r'TIMESTEP[ ]*\n[ ]*([0-9]+)'
        for filename in self.files:
            if not filename is None:
                if _inform:
                    print('counting configs in', filename)
                with open(filename,'r') as f:
                    fsize = os.stat(filename).st_size
                    _cnt = cnt
                    batch = min(25000000, fsize) # read at most 25MB
                    nbytes = batch               # count bytes
                    while batch > 0:
                        txt=f.read(batch)
                        match = re.findall(regex, txt)
                        if match:
                            match = np.array(match,dtype=int)
                            # end frame has been exceeded, stop counting
                            if match[0] > end:
                                nbytes = fsize
                            # start frame has been exceeded, start counting
                            if not match[-1] < start:
                                cnt += np.count_nonzero( (start<match) & (match<end))
                        batch = min(batch,fsize-nbytes)
                        nbytes+=batch
                    # count the file only if contains frame of interest
                    nfiles += 1 if cnt > _cnt else 0
        return cnt - (nfiles-1)

    def read_next(self, skip=False):
        ''' Read the next cofiguration form the file of this
            dump and return the data. '''
        # open in first call or check later
        if not self.f:
            if  self.nconfigurations == 0: self.open()
            else: return ( None, None, None)
        # just skip
        if skip and self._skip_this():
            return ( None, None, None)

        # initialize
        natoms = 0
        istep = box = None
        ok = False
        for line in self.f:
            if self._stopnow: break
            line = line.strip()
            if len(line) == 0: continue

            if line.startswith("ITEM: TIMESTEP"):
                istep = int(next(self.f).strip())
                self.n_dump.append( istep)
                self.nconfigurations += 1
                if _inform:
                    print('>> ... reading %d step dump' % self.n_dump[-1], end="\r")
                    sys.stdout.flush()
                continue

            elif line.startswith("ITEM: NUMBER OF ATOMS"):
                natoms = int( next(self.f).strip())
                continue

            elif line.startswith("ITEM: BOX BOUNDS"):
                _list = next(self.f).strip().split()
                xlo = float(_list[0])
                xhi = float(_list[1])
                xy = float(_list[2]) if len(_list) == 3 else 0.0
                _list = next(self.f).strip().split()
                ylo = float(_list[0])
                yhi = float(_list[1])
                xz = float(_list[2]) if len(_list) == 3 else 0.0
                _list = next(self.f).strip().split()
                zlo = float(_list[0])
                zhi = float(_list[1])
                yz = float(_list[2]) if len(_list) == 3 else 0.0
                # support old lammps dumps. the next line will
                # used right after in any case
                line = next(self.f).strip()
                if not line.startswith("ITEM: ATOMS"):
                    _list = line.split()
                    xy = float(_list[0])
                    xz = float(_list[1])
                    yz = float(_list[2])
                box = SimulationBox()
                box.set_from_lammps(xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz)

            if line.startswith("ITEM: ATOMS"):
                # read the next natoms lines, convert them to ndarray and then
                # fill the dump data.
                _lines = list( itertools.islice(self.f, natoms))
                data = np.zeros( natoms, dtype={ 'names':self.attributes, 'formats':self.types})
                try:
                    if _fixcamc:
                        for i, _line in enumerate( _lines):
                            tk = list(map( float, _line.strip().split()))
                            tk[0] = i
                            data[ i] = tuple( tk)
                    else:
                        for _line in _lines:
                            tk = list(map( float, _line.strip().split()))
                            tk[0] = tk[0] - 1
                            data[ int(tk[0])] = tuple( tk)
                    ok = True
                except:
                    del(data)
                    ok = False
                break

        if not ok:
            # close the current and try to open the next file
            if self._open_next_file():
                # on success read the next configuration
                return self.read_next( skip)
            else:
                # we have finished
                return None, None, None

        return ( istep, box, data)

    def _open_next_file( self):
        ''' Close the current and opens the next dump file skiping
            the first configuration. Returns true if there is a next
            file.
        '''
        # if a unit is open close it
        if self.f: self.f.close()
        # get the next file name. If it is not None open it and read the
        # first saved configuration since the restart saves as first the
        # last configuration of the previous dump file
        self.file = next(self._next_file)
        hasnext = False
        if self.file:
            hasnext = True
            self.f = open( self.file, 'r')
            # skip the first configuration
            self._skip_this()
        return hasnext

    def _skip_this(self):
        ''' Skip the current configuration (self.f should be
            available and valid).
        '''
        skiped = False
        for line in self.f:
            line = line.strip()
            if line.startswith("ITEM: NUMBER OF ATOMS"):
                skiped = True
                natoms = int(next(self.f).strip())
                # read the box block and check if box tilt line exist
                lines = list( itertools.islice(self.f,5))  # check the fifth line for box tilts
                line = lines[-1].strip()
                n = natoms if line.startswith("ITEM") else natoms + 1
                itertools.islice(self.f, n)
                for _ in range(n):
                    line = next(self.f)
                break
        return skiped

    def __next_file(self):
        ''' Get the next dump file in self.files. '''
        for file in self.files:
            yield file

    def _find_files( self):
        ''' Find file with self.file+'n' for n = 1,2,... in the same path,
            adds them in self.files and terminates this list with None. '''
        i = 0
        while True:
            i += 1
            name = self.filename + ".%d" % i
            if not os.path.isfile( name) and i >=2 :
                break
            if not os.path.isfile( name) and i==1 : continue
            self.files.append( name)
        self.files.append( None)

    def _scan_first( self):
        f = open( self.files[0])
        _lines = list ( itertools.islice(f, 10))
        f.close()
        self._natoms = int( _lines[3].strip())
        for line in _lines[8:]:
            if line.startswith("ITEM: ATOMS"):
                tk = line.strip().split()
                self.attributes = tk[2:]
        self.types = tuple([LammpsReader._typesDump[x] for x in self.attributes])

class LammpsReader(abcReader):
    ''' Implements lammps reader functionality.
        Attributes:
            data ({}) : dict keeping the data extructed from lammps data files.
                The foolowing keys:values are supported:
                    'box':(list) [(Dx,xlo,xhi),(Dy,ylo,yhi),(Dz,zlo,zhi),(xy xz yz)]
                    'Masses':(dict) {type-ID:mass}
                    'Types':(dict) {type-ID:name}
                    'type':(np.array(str)) atoms type array
                    'molecule':(np.array(int)) atoms molecule array
                    'charge':(np.array(float)) atoms charge array
                    'localid':(np.array(int)) atoms local (molecular) index array
                    'r':(np.array[:,3](float)) atoms coordinates array
                    'image':(np.array[:,3](int)) atoms image) array
                    'Bonds':(dict) {bond-ID:(molecule, type, atom1, atom2)}
                    'Angles':(dict) {angle-ID:(molecule, type, atom1, atom2, atom3)}
                    'Dihedrals':(dict) {dihedral-ID:(molecule, type, atom1, atom2, atom3, atom4)}
                    'Impropers':(dict) {improper-ID:(molecule, type, atom1, atom2, atom3, atom4)}
                    'Pair Coeffs':(dict) {type-ID:list of tokenized parameter line}
                    'Bond Coeffs':(dict) {type-ID:list of tokenized parameter line}
                    'Angle Coeffs':(dict) {type-ID:list of tokenized parameter line}
                    'Dihedral Coeffs':(dict) {type-ID:list of tokenized parameter line}
                    'Improper Coeffs':(dict) {type-ID:list of tokenized parameter line}
                    'thermo_n':(list) list steps found in log
                    'thermo_t':(np.array(float)) array of times found in log

            options ({}) : dict keeping the values of the options found in the
                lammps input script. The following keys:values are supported:
                    'timestep':(float) integration timestep
                    'run':(int) simulation steps
                    'ismd':(bool) md simulation flag
                    'isopt':(bool) minimization simulation flag
                    'atom_style':(str) lammps atom style
                    'bond_style':(str) lammps bond style
                    'angle_style':(str) lammps angle style
                    'dihedral_style':(str) lammps dihedra style
                    'improper_style':(str) lammps improper style
                    'units':(str) lammps unit
                    'thermo':(int) lammps thermo frequency

            log ({}) : dict keeping the values of the thermo keywords
                found in the log.

            dir (str): simulation directory path
            infile (str): lammps script .in file
            datafile (str): lammps data .data file
            logfile (str): lammps log .log file
            dumpfile (str): lammps dump .dump file
            bondsfile (str): lammps bonds file (reaxff)
            vdumpfile (str): lammps velocities file
            fdumpfile (str): lammps forces file

            natoms: sumulation number of atoms
            timestep: simulation timestep
    '''

    # utility staff
    # known variables names in data file
    _var_names = ("atoms", "bonds", "angles", "dihedrals", "impropers",
                  "atom types", "bond types", "angle types", "dihedral types", "improper types")
    # known blocks nameσ in data file
    _blk_names = ("Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers",
                "Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs")
    # "BondBond Coeffs", "BondAngle Coeffs", "MiddleBondTorsion Coeffs", "EndBondTorsion Coeffs" ,
    # "AngleTorsion Coeffs", "AngleAngleTorsion Coeffs", "BondBond13 Coeffs", "AngleAngle Coeffs"

    # known dump attributes
    _keyDump = [' '+x+' ' for x in [
           'id', 'type', 'q', 'x', 'y', 'z', 'xs', 'ys', 'zs', 'xu', 'yu', 'zu',
           'ix', 'iy', 'iz', 'vx', 'vy', 'vz', 'fx', 'fy', 'fz', 'mol' ]]
    _ATTRIBUTES = [ x.strip() for x in _keyDump]
    # dump attribute types (default of float32 and then correct)
    _typesDump = defaultdict( lambda:'f4', { k:'f4' for k in _ATTRIBUTES })
    for k in ('id', 'type', 'ix', 'iy', 'iz', 'mol'):
        _typesDump[k] = 'i4'

    # supported atom styles
    _supportedAtomStyles = [ "full", "molecular", "charged", \
                             "atomic", "bond", "angle", "sphere" ]
    # TODO support more styles; check https://lammps.sandia.gov/doc/read_data.html
    # for the complete list
    _lammpsAtomStyles = [  \
              "angle", "atomic", "body", "bond", "charge", "dipole", "dpd", \
              "electron", "ellipsoid", "full", "line", "meso", "molecular", "peri", \
              "smd", "sphere", "template", "tri", "wavepacket" ]

    # key for stoping quich scan of log file
    _start_log_key = "Memory usage per processor"
    _yaml_required = "%s.%s requires a valid system information file (system.yaml) pressent in the simulation directory."

    def __init__(self):
        ''' Initialize a LammpsReader object.'''
        super(LammpsReader, self).__init__()
        self.data = {                  # simulation data from the data file
            'thermo_n':[], 'thermo_t': None }
        self.options = {               # simulation option from the script file
            'timestep':0.0, 'isopt':False, 'ismd':False,
            'units':"", 'atom_style':"", 'pair_style':"",
            'bond_style':"", 'angle_style':"",
            'dihedral_style':"", 'improper_style':"",
            'run':0, 'thermo':0 }

        self.options['isopt'] = False
        self.options['ismd'] = False
        self.options['run'] = 0
        self.log = {}                   # thermo data

        # simulation name and files
        self.basename = ""
        self.dir = ""       # simulation directory
        self.infile = ""    # lammps script .in file
        self.datafile = ""  # lammps data .data file
        self.logfile = ""   # lammps log .log file
        self.dumpfile = ""  # lammps dump .dump file
        self.bondsfile = "" # lammps bonds file
        self.vdumpfile = "" # lammps velocities file
        self.fdumpfile = "" # lammps forces file

        self.dumpshandler = None    # dumps handler
        self.fromMaps = True        # maps simulation setup flag
        self._stopnow = False       # interupt flag
        self.error = ""             # error flag
        self.natoms = 0             # number of atoms

        self.nmolecules = 0         # number of moleculs

        self.boxes = []             # dump(s) boxes
        self.coordinates = []       # dump(s) coordinates
        self.offsets = []           # dump(s) offsets

        self.nconfLog = 0    # number of configurations in log(s)
        self.nconfDump = 0
        self.nconfBond = 0

        self.logStdInfo = \
            {
                "TotEng":("Total energy", "Energy", "kcal/mol"), \
                "Enthalpy":("Enthalpy", "Energy", "kcal/mol"), \
                "KinEng":("Kinetic energy", "Energy", "kcal/mol"), \
                "PotEng":("Potential energy", "Energy", "kcal/mol"), \
                "E_bond":("Bond energy", "Energy", "kcal/mol"), \
                "E_angle":("Angle energy", "Energy", "kcal/mol"), \
                "E_dihed":("Dihedral energy", "Energy", "kcal/mol"), \
                "E_impro":("Out-of-plane energy", "Energy", "kcal/mol"), \
                "E_vdwl":("Van der Waals energy", "Energy", "kcal/mol"), \
                "E_coul":("Coulomb energy", "Energy", "kcal/mol"), \
                "E_long":("Long-range energy", "Energy", "kcal/mol"), \
                "E_tail":("Tail energy", "Energy", "kcal/mol"), \
                "E_mol":("Molecular energy", "Energy", "kcal/mol"), \
                "E_pair":("Pair energy", "Energy", "kcal/mol"), \
                "Temp":("Temperature", "Temperature", "K"), \
                "Press":("Pressure", "Pressure", "atm"), \
                "P_xx":("P xx", "Pressure", "atm"), \
                "P_yy":("P yy", "Pressure", "atm"), \
                "P_zz":("P zz", "Pressure", "atm"), \
                "P_xy":("P xy", "Pressure", "atm"), \
                "P_xz":("P xz", "Pressure", "atm"), \
                "P_yz":("P yz", "Pressure", "atm"), \
                "Lx":("x", "Length", "\u212B"), \
                "Ly":("y", "Length", "\u212B"), \
                "Lz":("z", "Length", "\u212B"), \
                "xy":("xy", "Length", "\u212B"), \
                "xz":("xz", "Length", "\u212B"), \
                "yz":("yz", "Length", "\u212B"), \
                "Volume":("Volume", "Volume", "\u212B\u00B3"), \
                "RotKEdip":("Rotational kinetic energy (dipolar atoms)", "Energy", "kcal/mol"), \
                "RotKEgrn":("Rotational kinetic energy (granular atoms)", "Energy", "kcal/mol"), \
                "T_ave":("Time-average temperature", "Temperature", "K"), \
                "P_ave":("Time-Average pressure", "Pressure", "atm"), \
                "E_ave":("Time-Average energy", "Energy", "kcal/mol"), \
                "PE_ave":("Time-Average potential energy", "Energy", "kcal/mol"), \
                "rffeb":("Bond energy", "Energy", "kcal/mol"), \
                "rffea":("Atom energy", "Energy", "kcal/mol"), \
                "rffelp":("Lone-pair energy", "Energy", "kcal/mol"), \
                "rffemol":("Molecule energy", "Energy", "kcal/mol"), \
                "rffev":("Valence angle energy", "Energy", "kcal/mol"), \
                "rffepen":("Double-bond valence angle penalty", "Energy", "kcal/mol"), \
                "rffecoa":("Valence angle conjugation energy", "Energy", "kcal/mol"), \
                "rffehb":("Hydrogen bond energy", "Energy", "kcal/mol"), \
                "rffet":("Torsion energy", "Energy", "kcal/mol"), \
                "rffeco":("Conjugation energy", "Energy", "kcal/mol"), \
                "rffew":("Van der Waals energy", "Energy", "kcal/mol"), \
                "rffep":("Coulomb energy", "Energy", "kcal/mol"), \
                "rffefi":("Electric field energy", "Energy", "kcal/mol"), \
                "rffeqeq":("Charge equilibration energy", "Energy", "kcal/mol") \
            }

        self.inter_groups = None
        self.ninter_groups = 0
        self.inter_etype = ["coulomb", "vdw", ""]

        self.timestep = 0.0   # simulation timestep

    def set_files(self, filename, topo=None):
        ''' Set the simulation files.
            Args:
                filename (str): the path to the simulation log, dump or data
                    file. If the log file is provided, it will be scaned and
                    the parameters together with the rest of the simulation
                    files will be traced. If the dump file is provided then,
                    either a data file with the same name should be located
                    in the same directory or it should be provided using the
                    topo argument.
                topo (str): the path to the simulation data file.
            Raises:
                LammpsReaderException: if the file with the given name is not valid.
                LammpsReaderException: if the simulation files can not be located.
        '''
        # check
        self._set_files_chk(filename, topo)

        print('INFO: set reader (%s)'%filename)
        kind, _ = LammpsReader.is_supported(filename)
        if len(kind) == 0:
            message = '%s file is not supported.' % filename
            raise LammpsReaderException(message)
        elif kind == 'log':
            self.logfile = filename
            # extract simulation info
            self._scan_log()
        elif kind == 'dump':
            self.dumpfile = filename
            dumps = { self.basename: LammpsDump( self.shortfilename, self.dir+os.sep)}
            self.dumpshandler = LammpsDumpsHandler( dumps)

        if kind == 'data':  # given directly
            self.dumpshandler = LammpsDumpsHandler()
            self.datafile = filename
        elif not topo is None:  # given with topo
            self.datafile = topo
            print('INFO: set reader topology (%s)' % self.datafile)
        elif len(self.datafile) == 0:   # guess
            _datafile = self.dir + os.sep + self.basename + ".data"
            if os.path.isfile(_datafile):
                self.datafile = _datafile
                print('INFO: set reader topology (%s)' % self.datafile)

        # get the number of atoms from the dumps
        self.natoms = self.dumpshandler.get_natoms() #if not self.dumpshandler is None else 0

        # if system.yaml exists read extra system info.
        self._read_yaml()

    def set_attributes(self, attributes):
        ''' Set the attributes to be readed from the dump files.
            Args:
                attributes (list(str)): a list of attributes. The
                    possible values and their meeing are detailed
                    in lammps dump command description.
        '''
        return( self.dumpshandler.set_attributes( attributes))

    def read_topology(self, atomstyle=None):
        ''' Read simulation data file.
            Args:
                atomstyle (str): the lammps atom_style for parsing the Atoms block.
            Returns:
                bool: success status
            Raises:
                LammpsReaderException: if the simulation data file is invalid.
                LammpsReaderException: if the length of the atoms or bonded terms
                    blocks do no much with the expexted one.
        '''

        if len(self.datafile) == 0: return False

        # parse the file
        lines, blocks = LammpsReader.parse_data( self.datafile)
        if lines is None:
            msg = "%s file in not accesible" % self.datafile
            raise LammpsReaderException(msg)
        elif not lines[0].startswith("LAMMPS data file"):
            msg = "%s file header is invalid." % self.datafile
            raise LammpsReaderException(msg)

        # check and, if neccessary, set atom style from lammps input script
        # currently only styles in _supportedAtomStyles are supported
        if atomstyle is None:
            atomstyle = self.options['atom_style']
        if len(atomstyle) == 0: atomstyle='full'
        it = 2 if atomstyle in ('angle', 'bond', 'full', 'line', 'molecular', 'tri') else 0 if atomstyle == 'charge' else 1
        im = 2 if atomstyle == 'template' else 1
        iq = 2 if atomstyle in ('dipole', 'electron') else 1 if atomstyle == 'charge' else 3
        read_charge = atomstyle in ('dipole', 'electron', 'full', 'charge')
        read_molecule = atomstyle in  ('angle', 'bond', 'full', 'line', 'molecular', 'template', 'tri')

        # parse box block:
        #  [ (lx, xlo, xhi), (ly, ylo, yhi), (lz, zlo, zhi), (xy, xz, yz)]
        _box = self.data['box'] = []
        for line in blocks[ 'box']:
            tk = line.split()
            if line.endswith("xlo xhi"):
                lst = list(map( float, tk[:2]))
                _box.append( ( lst[1]-lst[0], lst[0], lst[1]))
            elif line.endswith("ylo yhi"):
                lst = list(map( float, tk[:2]))
                _box.append( ( lst[1]-lst[0], lst[0], lst[1]))
            elif line.endswith("zlo zhi"):
                lst = list(map( float, tk[:2]))
                _box.append( ( lst[1]-lst[0], lst[0], lst[1]))
            elif line.endswith("xy xz yz"):
                _box.append( ( list(map(float,tk[:3]))))
            else:
                _box.append( [ [0.0, 0.0, 0.0]])

        # parse Masses block
        _masses = self.data['Masses'] = {}
        _types = {}
        for line in blocks[ 'Masses']:
            lst = line.split()
            # TODO mapping between lammps types and reaxff types.
            # check for inline comments contain the type names, assuming
            # the format: id mass #[ ]*type
            if len( lst) == 4:
                _types[ lst[0]] = lst[3]
            elif len(lst) == 3:
                _types[ lst[0]] = lst[2][1:]
            _masses[ lst[0]] = float( lst[1])
        if len(_types) > 0:
            self.data['Types'] = _types

        # parse Atoms block
        natoms = blocks[ 'atoms']
        if self.dumpshandler is None:
            self.natoms = natoms
        elif not natoms == self.natoms:
            msg = "The number of atoms (%d) in the data file  differs from the number of atoms (%d) in the dump file." % (
                natoms, self.natoms)
            raise LammpsReaderException(msg)
        nlines = len( blocks[  'Atoms'])
        if not natoms == nlines:
            msg = "Check Atoms block in %s file (%d lines found instead of %d)."
            msg = msg % ( self.datafile, nlines, natoms)
            raise LammpsReaderException(msg)
        _type = {}
        _molecule = self.data['molecule'] = np.zeros(natoms, dtype=int)
        _charge = self.data['charge'] = np.zeros(natoms, dtype=float)
        _localid = self.data['localid'] = np.zeros(natoms, dtype=int)
        _r = self.data['r'] = np.zeros((natoms,3), dtype=float)
        _image = self.data['image'] = np.zeros((natoms,3), dtype=int)
        # for all styles the last record will be either a coordinate (float)
        # or an image (int).
        _lines = blocks[  'Atoms']
        lastrecord_ = _lines[0].split()[-1] if len(_lines)>1 else ""
        has_image = not "." in lastrecord_
        previous_ = -1
        localid_ = 1
        for i, line in enumerate( _lines):
            tk = line.split()
            _id = int( tk[0])-1 # if it == 0 else i
            _type[_id] =  tk[it]
            molecule_ = int( tk[im]) if read_molecule else 0
            _molecule[ _id] = molecule_
            _charge[ _id] = float(tk[iq]) if read_charge else 0.0
            if not molecule_ == previous_:
                localid_ = 1
                previous_ = molecule_
            _localid [ _id] = localid_
            if has_image:
                r_ = tuple(map( float, tk[-6:-3]))
                image_ = tuple(map( int, tk[-3:]))
            else:
                r_ = tuple(map( float, tk[-3:]))
                image_ = (0,0,0)
            _image[ _id] = image_
            _r[ _id] = r_
            localid_ += 1
        self.data['type'] = np.array([_type[k] for k in sorted( _type.keys() )])

        # parse bonded terms blocks
        _atom_molecule = self.data['molecule']
        for k, v in {'Bonds':4, 'Angles':5, 'Dihedrals':6, 'Impropers':6}.items():
            _k = self.data[k] = {}
            for line in blocks[ k]:
                tk = line.split()
                # (id molecule type atom1 atom2 [atom3 [atom4]])
                # or if we want now to use indexes from
                #    tuple( map( lambda x: int(x)-1, tk(:3)))
                lst_ = list(map( int, tk[1:v]))
                lst_.insert(0, _atom_molecule[ lst_[1]-1]+1)
                _k[ int( tk[0])] = tuple( lst_)

        # parse coefficients
        for k in ('Pair Coeffs', 'Bond Coeffs', 'Angle Coeffs', 'Dihedral Coeffs', 'Improper Coeffs'):
            _k = self.data[k] = {}
            _s = k.split()[0].lower() + "_style"
            _style = self.options[_s] if _s in self.options else 'UNK'
            for line in blocks[k]:
                tk = line.split()
                lst_.insert(0, tk[1] if tk[1].isalpha() else _style)
                _k[ int(tk[0])] = tuple(tk[1:])

        if read_molecule:
            self.nmolecules = _atom_molecule.max()

        return True

    def read_next_frame(self, skip=False):
        ''' Read simulation's dump file.
            Args:
                skip (bool): if true the frame is skiped.
            Returns:
                int, Simulationbox, {}: the step, the simulation box and the
                    {attribute:data} dictionaly for the next frame.
        '''
        step, box, data = self.dumpshandler.read_next(skip)

        if self.do_whole:
            if 'x' in data and not 'xw' in data:
                molecule = self.get_atom_molecule() - 1
                bonds = self.get_bonds()
                r = np.zeros((data['x'].size,3),dtype=np.float32)
                for k, v in {'x':0,'y':1,'z':2}.items():
                    if k in data: np.copyto(r[:, v], data[k])
                rw = fastwhole(r, molecule, bonds, box.a, box.b, box.c)
                data['x'][:] = rw[:,0]
                data['y'][:] = rw[:,1]
                data['z'][:] = rw[:,2]
        elif self.do_unwrap:
            # check if unwrapped coordinates exist already and,  
            # if needed, create links to the normal coordinates
            if 'xu' in data and not 'x' in data:
                # just create a link
                data['x'] = data['xu']
                data['y'] = data['yu']
                data['z'] = data['zu']
                # or allocate memory and copy data
                # for k, v in {'x':'xu','y':'yu','z':'zu'}:
                #     data[v] = np.zeros(data[k].size, dtype=np.float32)
                #     np.copyto(data[v], data[k])
            elif 'x' in data and 'ix' in data:
                r = np.zeros((data['x'].size,3),dtype=np.float32)
                ip = np.zeros((data['ix'].size,3),dtype=np.int32)
                for k, v in {'x':0,'y':1,'z':2}.items():
                    if k in data: np.copyto(r[:, v], data[k])
                for k, v in {'ix':0,'iy':1,'iz':2}.items():
                    if k in data: np.copyto(ip[:, v], data[k])
                ru = fastunwrapv(r, ip, box.va, box.vb, box.vc)
                # a, b, c = box.va, box.vb, box.vc
                # for _r, ( i, j, k) in zip( r, ip):
                #     _r[:] += a * i + b * j + c * k
                data['x'][:] = ru[:,0]
                data['y'][:] = ru[:,1]
                data['z'][:] = ru[:,2]

        return (step, box, data)

    def get_type_mass(self):
        ''' Retrun the {type:mass} dictionary retrieved from data file.
            If the type names are provided from the comments in the data
            file or from the system.yaml file, type name  will be used
            as key. Otherwise, the type id (casting to str) will be used
            instead.
        '''
        typemass = self.data['Masses']
        if 'Types' in self.data and set(self.data['Types']) == set(typemass):
            typemass = { self.data['Types'][k]:typemass[k] for k in typemass.keys()}
        elif not self.yaml is None:
            names = list(
                itertools.chain.from_iterable([
                    s['atomtype'] * s['nmolecules'] for s in self.yaml['specie']
                ]))
            types = self.data['types']
            if len(names) == len(types):
                typemass = { n:t for n, t in zip(names,types)}
        return typemass

    def get_atom_molecule(self):
        ''' Retruns atom's molecule retrieved from data file. '''
        return self.data['molecule']

    def get_atom_mass(self):
        ''' Retruns atom's mass retrieved from data file. '''
        masses = self.data['Masses']
        types = self.data['type']
        return np.array( [ masses[ types[iat]] for iat in range( len(types)) ])

    def get_atom_type(self):
        ''' Retruns atom's types retrieved from data file. '''
        types = super(LammpsReader, self).get_atom_type()
        if len(types) == 0:
            types = self.data['type']
        return types

    def get_atom_name(self):
        ''' Retruns atom's names retrieved from yaml file. '''
        names = super(LammpsReader, self).get_atom_name()
        if len(names) == 0:
            msg = LammpsReader._yaml_required % (self.__class__.__name__,inspect.currentframe().f_code.co_name)
            raise LammpsReaderException( msg)
        return names

    def get_atom_element(self):
        ''' Retruns atom's element retrieved from data file. '''
        elements = super(LammpsReader, self).get_atom_element()
        if len(elements) == 0:
            # TODO add guess element based on mass functionality
            msg = LammpsReader._yaml_required % (self.__class__.__name__,inspect.currentframe().f_code.co_name)
            raise LammpsReaderException( msg)
        return elements

    def get_atom_charge(self):
        ''' Retruns atom's charge retrieved from data file. '''
        return self.data['charge']

    def get_molecule_name(self):
        ''' Return the molecues residue names '''
        names = super(LammpsReader, self).get_molecule_name()
        if len(names) == 0:
            msg = LammpsReader._yaml_required % (self.__class__.__name__,inspect.currentframe().f_code.co_name)
            raise LammpsReaderException( msg)
        return names

    def get_bonds(self):
        ''' Retruns bonds' atoms indexes. '''
        try:
            bonds = np.array( [ (v[2],v[3]) for k, v in sorted(self.data['Bonds'].items(), key=lambda item:item[0])])
        except:
            bonds = np.array((),dtype=np.int32)
        return bonds

    def get_topology(self):
        ''' Return the system topology (see McReader) '''
        message = '%s.%s is not yet supported.' % (self.__class__.__name__,inspect.currentframe().f_code.co_name)
        raise LammpsReaderException( message)

    def get_thermo_vars(self):
        ''' Return the thermodynamic variables names found in the log file.
            Returns:
                list(str): a list with thermo variables names
        '''
        if len(self.log) == 0: return []
        return list(self.log.keys())

    def get_thermo_info(self, var):
        ''' Return the information for the given thermodynamic variable.
            Returns:
                (description,name,units) : a tuple with the short descritpion of
                    description, the name and the units of the variable.
        '''
        try:
            (v, t, u) = self.logStdInfo.setdefault(var, ())
        except (ValueError):
            v = var
            t = ""
            u = ""
        return (v, t, u)

    def get_thermo_data(self, var=None):
        ''' Returns the data for the thermodynamic variables retrieved from the
            log file.
            Args:
                var (str) : the name of variable. If var is provided the
                    method returns only the data for the specific variable.
            Returns:
                [] : a list of values for the specific variable.
        '''
        if len(self.log) == 0: return []
        return self.log.setdefault(var, [])

    def print_simulation_info(self):
        ''' Print basic simulation info. '''

        _infile = self.infile if len(self.infile) > 0 else "-"
        _datafile = self.datafile if len(self.datafile) > 0 else  "-"
        _dumpfile = self.dumpfile if len(self.dumpfile) > 0 else  "-"
        _runtype = "opt" if self.options['isopt'] else "md" if self.options['ismd'] else "-"
        _nsteps = self.options['run']
        _nsteps = "-" if _nsteps == 0 else str(_nsteps)

        print("-"*20)
        print("engine        : %s" % "lammps")
        print("in script     : %s" % _infile)
        print("data file     : %s" % _datafile)
        print("dump file     : %s" % _dumpfile)
        print("# of atoms    : %d" % self.natoms)
        print("timestep (%s) : %f" % ("fs", self.options['timestep']))
        print("# of steps    : %s" % _nsteps)
        print("run type      : %s" % _runtype)
        print("-"*20)

    def thermo_tocsv(self, every=1, csvfilename=""):
        ''' Convert log to csv format and save the file in
            simulation directory with the given name.
            Args:
                every (int): dump every this many records.
                csvfilename (str): the cvs file name. If it
                is not provided, basename.csv name  will be
                used instead.
        '''
        if len(self.log) == 0: return

        csvfilename = self.dir + os.sep + self.basename + ".csv"
        f = open(csvfilename, 'w')
        keys = sorted( list( self.log.keys()))
        values = [ self.log[k] for k in keys]
        keys = [x.lower() for x in keys]

        # all the variables are of the same length
        n = len( values[0])

        # write header
        line = " ".ljust(2) + "step".ljust(9)
        for k in keys: line += k.ljust(15) + " "
        f.write(line[0:len(line)-1] + "\n")
        # write data
        for i in range(n):
            if not i%every == 0: continue
            line = " "+str(i).ljust(10)
            for v in values: line += str(v[i]).ljust(15) + " "
            f.write(line[0:len(line)-1] + "\n")

        f.close()

    def count_frames( self, start=-1, end=sys.maxsize):
        ''' Return the number of frames in lammps dump files. '''
        return self.dumpshandler.count_frames(start, end)

    def read_log(self):
        ''' Read simulation's log file.'''
        # check
        if len( self.logfile) == 0: return
        f = open( self.logfile, 'r')
        _thermo_n = self.data['thermo_n'] = []
        for line in f:
            if self._stopnow: break
            line = line.strip()
            if len(line) == 0: continue
            if line[0] == "#": continue
            if line.find(" Step   ") != -1:
                istep = int( line.split()[2])
                # skip the first step of a new log files if needed
                if len(_thermo_n) > 0 and istep == _thermo_n[-1]:
                    continue
                else:
                    _thermo_n.append( istep)
                    self._read_longlog(f)
            if line.find("Memory usage per processor") != -1 and self.shortlog:
                self._read_shortlog(f)
        f.close()
        # set the correspoding step time
        _thermo_n = self.data['thermo_n']
        _thermo_t = self.data['thermo_t'] = np.array( _thermo_n, dtype=np.float32) * self.options['timestep']
        # keep number of nconfigurations in log file
        self.nconfLog = len( _thermo_t)

    @staticmethod
    def is_supported(filename):
        ''' Check if the given file is a LAMMPS simulation file.
            The reader supports the types: 'log' for lammps log
            file, 'dump' for lammps trajectory files and 'data'
            for lammps data files.
            The data files (puropose 'top') are recongnized from
            their third line which should have the 'atoms' key as
            second record.
            The dump files (puropose 'trj) are recongnized form their
            first line which should start with the 'ITEM: ' key.
            The log file (puropose 'log') are recongnized from their
            first or second line which shoud start with 'LAMMPS ('
            key.
        '''
        s = None
        try:
            with open( filename, 'r') as f:
                line = f.readline().strip() # reade the first line
                if line.startswith("ITEM:"):
                    s =  ('dump','trj')
                elif line.startswith("LAMMPS ("):
                    s =  ('log','log')
                if s is None: # if nothing so far read the second
                    line = f.readline().strip()
                    if line.startswith("LAMMPS ("):
                        s =  ('log','log')
                    if s is None: # if nothing so far read the third
                        tk = f.readline().strip().split()
                        if len(tk)>1 and tk[1] == 'atoms':
                            s = ('data','top')
        except:
            s = None
        if s is None: s = ('','')
        return s

    def _scan_log( self):
        ''' Scan the log file to extract simulation info.
            Raises:
                LammpsReaderException: if the atom_style found in lammps input
                    script is not supported.
                LammpsReaderException: if no dump commands found in lammps input script.
        '''
        # check
        if len(self.logfile) == 0: return
        # initialize
        self.shortlog = True
        self.isopt = False
        self.ismd = False
        dumps = {}                  # dumps found in the log file
        _name = self.dir + os.sep   # use full path if possible
        # partial parsing of the log file to spot info of interest
        f = open( self.logfile, 'r')
        line = f.readline()
        if not line.startswith("LAMMPS"):
            msg = "%s file header is invalid." % self.logfile
            raise LammpsReaderException( msg)

        for line in f:
            # stop on start log block
            if line.startswith( LammpsReader._start_log_key): break
            # skip comments
            if line[0] == "#": continue
            # skip line with variables (lammps stepwise line parsing)
            if "${" in line: continue
            # strip and skip empty lines
            line = line.strip()
            if len(line) == 0: continue
            if line.startswith( "read_data"):
                self.datafile = _name + line.split()[1]
            elif line.startswith( "dump "):
                # if dump file dose not exist skip it
                try:
                    dump = LammpsDump( line, _name)
                    dumps[ dump.name] = dump
                except:
                    pass
            elif line.startswith( "thermo_modify"):
                tk = line.split()
                if "line" in tk:
                    i  = tk.index( "line")
                    if tk[i+1] == "multi": self.shortlog = False
            elif line.startswith( "timestep "): self.options['timestep'] = float( line.split()[1])
            elif line.startswith( "minimize "): self.options['isopt'] = True
            elif line.startswith( "units "): self.options['units'] = line.split()[1]
            elif line.startswith( "atom_style "): self.options['atom_style'] = line.split()[1]
            elif line.startswith( "pair_style"): self.options['pair_style'] = line.split()[1]
            elif line.startswith( "bond_style"): self.options['bond_style'] = line.split()[1]
            elif line.startswith( "angle_style"): self.options['angle_style'] = line.split()[1]
            elif line.startswith( "dihedral_style"): self.options['dihedral_style'] = line.split()[1]
            elif line.startswith( "improper_style"): self.options['improper_style'] = line.split()[1]
            elif line.startswith( "run "): self.options['run'] = int( line.split()[1])
            elif line.startswith( "thermo "): self.options['thermo'] = int( line.split()[1])

        self.timestep = self.options['timestep'] / 1000.0 # in ps (MDAnalysis compatiblity)

        # check atom style
        if not self.options['atom_style'] in LammpsReader._supportedAtomStyles:
            msg = " %s atom style is not supported." % self.options['atom_style']
            raise LammpsReaderException(msg)
        # setup dumps handler
        # if no dump files have been spotted rain an exception
        if len(dumps) == 0:
            msg = "no valid dump command found in %s." % self.logfile
            raise LammpsReaderException(msg)
        self.dumpshandler = LammpsDumpsHandler( dumps)


    def _read_longlog(self, f):
        ''' Read multiline thermo format. '''
        for line in f:
            if self._stopnow: break
            line = line.strip()
            if len(line) == 0: continue
            if line.find("vectors:") != -1: continue
            if line.find("Loop time") != -1: break
            if line.find("Step") != -1:
                istep = int( line.split()[2])
            else:
                s = line.split()
                for i in range(line.count("=")):
                    j = i * 3
                    # set the correct name for intermolecular energy variables
                    name = self._check_name(s[j])
                    self.log.setdefault(name, []).append(float(s[j + 2]))
                    #self.log.setdefault(s[j], []).append(float(s[j + 2]))

    def _read_shortlog(self, f):
        ''' Read single line thermo format. '''
        head = next(f)
        names = head.strip().split()
        variables = [[] for i in range(len(names))]
        icnt = 0
        _thremo_n = self.data['thermo_n']
        skip = len( _thremo_n) > 0
        for line in f:
            if self._stopnow: break
            if (len(line) == 0): continue
            if icnt == 0 and skip:
                icnt += 1
                continue
            skip = False
            line = line.strip()
            if line.find( "Loop time") != -1: break
            list_ = line.split()
            for v, d in zip( variables, list_):
                v.append( float( d))
            _thremo_n.append( int( list_[0])) # step is always the first row
        for n, v in zip( names, variables):
            # set the correct name for intermolecular energy variables
            name = self._check_name(n)
            self.log[name] = v

    def _check_name(self, name):
        ''' Create a proper name for intermolecular energy thermo variables. '''
        name_ = name
        i = name_.find("inter[")
        if i != -1:
            n = int(name[6:len(name)-1])
            (m, t) = divmod(n, 3)
            name_ = self.inter_etype[t] + " inter"
            if m < self.ninter_groups:
                name_ = self.inter_groups[m] + " " + name_ # pylint: disable=unsubscriptable-object
            name_.strip()
            if not name_ in self.logStdInfo:
                self.logStdInfo[name_] = (name_ + " energy", "Energy", "kcal/mol")
        return name_

    def _extract_coeffs(self, line):
        ''' Parse data file coefficient blocks. If a line comment
            exists, the name of the corresponding type is expected
            to be found (MAPS convention). '''
        h = line.find("#")
        # remove the part left to the first parenthesis
        p = line.find("(")
        comment = ''
        if p == -1:
            p = len(line)
        else:
            comment = line[p+1:-2]
        lst = line[:p].split()

        id_ = int( lst[0])
        types = []
        if h == -1:
            params = list(map( float, lst[1:]))
        else:
            indx = lst.index( "#")
            params = list(map( float, lst[1:indx]))
            types[:] = lst[indx+1:]
        # (id,  types, parameters)
        params.insert( 0, (types, comment))
        return id_, params

    @classmethod
    def create(cls, filename, topo=None):
        ''' Use the log or dump file to create a LammpsReader object. If the filename is
            a directory then the first .log file in this directory will be used instead.
            Args:
                filename (str): the name of the simulation file of of diretcory where the
                    simulation files reside. In the later case the filename suffix should
                    indicate their type i.e. '.log' for the log files, '.dump' for the dump
                    files, '.data' for the data files.
            Returns:
                LammpsReader
            Raises:
                LammpsReaderException: If the topology is not specified (i.e.,  neither a YAML
                    nor a data file was located).
                
        '''
        # check
        if filename is None or len(filename) == 0:
            filename = "."+os.sep # set the current directory

        # keep input
        _filename = os.path.expanduser( filename)
        _topo = topo
        # check if file
        if not os.path.isfile( _filename):
            # check if directory
            if os.path.isdir( _filename):
                _dirname = _filename
                filelist = glob.glob(_dirname + os.sep + "*.log")
                # if log files found use the first else check for dumps
                if len( filelist) == 0:
                    filelist = glob.glob(_dirname + os.sep + "*.dump")
                    # if dump files found use the first else fail
                    if len( filelist) == 0:
                        message = "no simulation files found in : %s" % filename
                        print(message)
                        return None
                    else:
                        _filename = filelist[0]     # set the dump file
                        if topo is None:
                            datalist = glob.glob(_dirname + os.sep + "*.data")
                            if not len(datalist) == 0:
                                _topo = datalist[0] # set the data file
                _filename = filelist[0]             # set the log file
        try:
            reader = cls()
            reader.set_files(_filename, topo=_topo)
            if not reader.read_topology():
                if reader.yaml is None:
                    raise LammpsReaderException("fail to read the topology file")
        except LammpsReaderException as err:
            reader = None
            print(err)

        

        return reader

    @staticmethod
    def parse_data(datafile):
        ''' Parse lammps data file.
            Args:
                datafile (str): data file full name.
            Returns:
                list, dict: lines list and blocks dictionary
        '''
        lines = []
        blocks = defaultdict(list)
        with open(datafile, 'r') as f:
            lines = list(map(lambda x: x.strip(), f.readlines()))
            f.close()
            blocks = defaultdict(list)
            for line in lines[1:]:
                if len(line) == 0:
                    continue
                found = False
                for v in LammpsReader._var_names:
                    if line.endswith(v):
                        blocks[v] = int(line.split()[0])
                        found = True
                        break
                if not found:
                    _line = line.split("#")[0].strip() # remove comments
                    if any(tuple(_line.endswith(s) for s in ("xhi", "yhi", "zhi", "yz") )):
                        blocks["box"].append(line)
                    elif _line in LammpsReader._blk_names:
                        blk = blocks[_line]
                    elif " Coeffs" in _line:
                        blk = blocks[_line]
                    else:
                        if not blk == None:
                            blk.append(line)
        return lines, blocks

    @staticmethod
    def write_data(blocks, datafile):
        ''' Parse lammps data file.
            Args:
                blocks (dict): dictionary with data blocks.
                datafile (str): data file full name.
        '''
        with open(datafile,'w') as f:
            # get the keywords
            keys = list( blocks.keys())

            lines = []
            # comment lines
            line = blocks['comment'] if 'comment' in blocks else "LAMMPS data file %s whole." % datafile
            lines.append( line)

            # integer value keywords (atoms, bonded terms and types)
            for key in ['',
                    'atoms', 'bonds', 'angles', 'dihedrals', 'impropers',
                    '',
                    'atom types', 'bond types', 'angle types', 'dihedral types', 'improper types',
                    '' ]:
                value = blocks[key] if key in blocks else 0
                lines.append( '' if key == '' else str(value) + ' '+ key)
                if key in keys: keys.remove( key)

            # box keyword
            lines += blocks['box']
            lines.append('')
            keys.remove( 'box')

            # topology and coefficients keywords
            common = ['Masses', 'Atoms', 'Bonds', 'Angles', 'Dihedrals', 'Pair Coeffs', 'Bond Coeffs', 'Angle Coeffs', 'Dihedral Coeffs', 'Improper Coeffs']
            extra = list( set(keys) - set(common))
            all = common + extra
            for key in all:
                if not key in blocks: continue
                lines.append(key)
                lines.append('')
                lines += blocks[key]
                keys.remove( key)
                lines.append('')

            # write down the lines
            for line in lines:
                f.write(line+"\n")

    def unwrap_frame(self, coords, images, box):
        ''' Unwrap the given frame.
            Args:
                coords (np.array((natoms,3), dtype=float)): atom coordinates array.
                images (np.array((natoms,3), dtype=int)): atom images array.
                box (SimulationBox) : the simulation box of the frame.
        '''
        a, b, c = box.va, box.vb, box.vc
        for r, (i,j,k) in zip(coords,images):
            r[:] += a * i + b * j + c * k
