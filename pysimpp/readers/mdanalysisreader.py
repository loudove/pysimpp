# -*- coding: utf-8 -*-

import os
import sys
import inspect
import numpy as np
from itertools import islice

import MDAnalysis

from .reader import abcReader # pylint: disable=import-error
from pysimpp.utils.simulationbox import SimulationBox
from pysimpp.fastpost import fastwhole # pylint: disable=no-name-in-module

_inform = True

class MDAnalysisReaderException(Exception): ...

class MDAnalysisReader(abcReader):
    ''' Implements reader functionality based on MDAnalysis.
        Attributes:
            u (MDAnalysis.Universe): MDAnalysis universe.
            topofile (str): the topology file (whatever MDAnanlysis supports for gromacs).
            trajfile (str): the trajectory file (whatever MDAnanlysis supports for gromacs).
            dir: the simulation directory
            frames: readed frames info
            natoms: sumulation number of atoms
            timestep: simulation timestep
    '''

    # set various dump keywords
    _KEYS = [ 'id', 'name', 'x', 'y', 'z', 'el', 'type']

    # set various dump keywords
    _TYPESKEYS = { k:'f4' for k in _KEYS }
    _TYPESKEYS[ 'id'] = 'i4'
    _TYPESKEYS[ 'name'] = 'U5' # set variable length for name
    _TYPESKEYS[ 'el'] = 'U2'
    _TYPESKEYS[ 'type'] = 'U5'

    def __init__(self):
        ''' Initialize a MDAnalysisReader reader object. '''
        super(MDAnalysisReader, self).__init__()

        self.u = None           # MDanalysis handler

        self.files = []
        self.topofile = None    # the topology file (whatever MDAnanlysis supports)
        self.trajfile = None    # the trajectory file (whatever MDAnanlysis supports)
        self.frames = []        # keep readed frames info
        self.data = {}          # gro file data
        self.natom = 0          # number of atoms
        self.nsteps = 0         # number of steprs read from trajectory file
        self.timeoffset = 0.0   # fix timestep (in case of time reset and multiple trajectory files)
        self.istepoffset = 0     # fix integration step

    def set_files(self, filename, topo=None):
        ''' Set the topology and the trajectory files. The convention is that
            both files have the same base name but different extension.
            Args:
                filename (str): the trajectory file name. A topology (tpr) file
                    with the same basename is expected to be found in the same
                    directory.
            Raises:
                MDAnalysisReaderException: if the file with the given name is not valid.
                MDAnalysisReaderException: if the topology file is not found.
        '''
        # check
        self._set_files_chk(filename, topo)
        # apply policy
        self._set_files_policy(filename, topo)
        # initialize
        self.filename = filename
        self.files.append( filename)
        self.__set_files_initialize()

    def _set_files_policy(self, filename, topo):
        ''' Apply readers policy for setting the trajectory (self.trajfile)
            and the topology files (self.topofile). '''
        # just set the files
        self.trajfile = filename
        self.topofile = topo

        print('INFO: set reader (%s)'%self.trajfile)
        print('INFO: set reader topology (%s)' % self.topofile)

    def __set_files_initialize(self):
        # initialize multiple trj files
        self._find_files()
        self._next_file = self.__next_file()

        # initialize the MDAnalysis univerce
        self.trajfile = next(self._next_file) # get the first file (always available)
        self._prepare()

    @staticmethod
    def is_supported(filename):
        ''' Check if MDAnalysis supports the topolog file.
            Args:
                filename (str) : topology file name.
            Returns:
                (type (str), puropose (str)): the type of the file and
                    its purpose  ('log', 'inout', 'top' or 'trj') if the
                    file is not suppoted both strings will be empty.
        '''
        parser = MDAnalysis.core._get_readers.get_parser_for(filename,format=None)
        if not parser is None:
            _, extension = os.path.splitext( filename)
            s = (extension[1:],'top')
        else:
            s = ('','')
        return s


    def set_attributes(self, attributes):
        ''' Set the attributes to be readed from the dump files. MDAnalysisReader
            supports the attributes:
                ['id', 'name', 'x', 'y', 'z', 'el', 'type']
            adn they should be set before calling the self.read_dump or
            self.read_next_frame methods.
        '''
        self.canread = all( [x in self._KEYS for x in attributes.split()])
        if self.canread:
            self.attributes = attributes
            self.names = attributes.split()
            self.types = [self._TYPESKEYS[x] for x in self.names]
        return self.canread

    def read_topology(self):
        ''' There is no data file in gromacs. Just set the number of atoms.
            Returns:
                bool: success status
        '''
        if not self.u is None:
            self.natoms = len( self.u.atoms)
        elif len(self.topofile):
            f = open(self.topofile,'r')
            _atoms = MDAnalysisReader.read_gro(f)
            _names = ['resid', 'resname', 'atomname', 'atomid', 'x', 'y', 'z', 'vx', 'vy', 'vz']
            _types = ['i4', '<U5', '<U5', 'i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4']
            data = np.array( _atoms, dtype={ 'names':_names, 'formats':_types})
            self.data = { k:data[k] for k in _names } # keep in LammpsReader self.data fashion
            self.natoms = len( _atoms)
        return (self.natom > 0)

    def read_next_frame(self, skip=False):
        ''' Read the next frame in the trajectory file.
            Args:
                skip (bool): if true the frame is skiped.
            Returns:
                int, Simulationbox, {}: the step, the simulation box and the
                    {attribute:data} dictionaly for the next frame.
        '''
        box = None
        data = None
        #istep = sys.maxsize
        istep = None
        self.nsteps += 1
        if self.nsteps == 1 and not self.do_whole and self.do_unwrap:
            print("WARNING: unwrapped trajectory is assumed.")

        if not self.u is None:
            ts = self._next_frame()

            atoms = self.u.atoms
            natoms = len(atoms)

            if not skip and natoms == self.natoms and not ts == None:
                # get the box
                boxdim = ts.dimensions
                box = SimulationBox()
                box.set_from_scalars( boxdim[0], boxdim[1], boxdim[2], boxdim[3], boxdim[4], boxdim[5])
                # get the positions
                _x = ts.positions
                # create the data (LDP TODO automate this based on the given attributes)
                data = np.zeros( natoms, dtype={ 'names':self.names, 'formats':self.types} )
                data['id'][:] = atoms.ids
                if 'type' in self.names: # LDP: check the attributes set
                    data['type'][:] = atoms.types
                # data['el'][:] = map( lambda x: MDAnalysis.topology.core.guess_atom_element(x), atoms.names)

                if self.do_whole:
                    molecule = self.get_atom_molecule() - 1
                    bonds = self.get_bonds()
                    rw = fastwhole(ts.positions, molecule, bonds, box.a, box.b, box.c)            
                    data['x'][:] = rw[:,0]
                    data['y'][:] = rw[:,1]
                    data['z'][:] = rw[:,2]
                # elif self.do_unwrap: # this is more like the "make whole"
                #     ru = self.u.atoms.unwrap(reference=None)
                #     data['x'][:] = ru[:,0]
                #     data['y'][:] = ru[:,1]
                #     data['z'][:] = ru[:,2]
                else:
                    data['x'][:] = ts.positions[:,0]
                    data['y'][:] = ts.positions[:,1]
                    data['z'][:] = ts.positions[:,2]

                # the step of the trajectory frame (inline with LammpsReader)
                _istepoffset = self.istepoffset
                istep = ts.data['step'] + _istepoffset
                self.frames.append( (ts.frame, ts.time+self.timeoffset, ts.data['step']+_istepoffset)) # keep frames info (debug info)

                if _inform:
                    print('>> ... reading %d step dump' % istep, end="\r")
                    sys.stdout.flush()

        return ( istep, box, data)

    def get_type_mass(self):
        ''' Retrun the {type:mass} dictionary retrieved from atoms
            type and masses.
        '''
        # find the types in the system
        typemass = {}
        atoms = self.u.atoms
        for atom in atoms: # pylint: disable=not-an-iterable
            atomtype = atom.type
            if not atomtype in typemass:
                typemass[ atomtype] = atom.mass
        return typemass

    def get_atom_molecule(self):
        ''' Retruns atom's molecule retrieved from resid. '''
        return self.u.atoms.molnums + 1 # fix the index to conform with lammpsreader

    def get_atom_mass(self):
        ''' Retruns atom's mass retrieved from the tpr file. '''
        return self.u.atoms.masses

    def get_atom_type(self):
        ''' Retruns atom's types retrieved from the tpr file. '''
        types = super(MDAnalysisReader, self).get_atom_type()
        if len(types) == 0:
            types = self.u.atoms.types
        return types

    def get_atom_name(self):
        ''' Retruns atom's names retrieved from the tpr file. '''
        names = super(MDAnalysisReader, self).get_atom_name()
        if len(names) == 0:
            names = self.u.atoms.names
        return self.u.atoms.names

    def get_atom_element(self):
        ''' Retruns atom's element guessed from atom names in the tpr file. '''
        elements = super(MDAnalysisReader, self).get_atom_element()
        if len(elements) == 0:
            elements = np.array( [ MDAnalysis.topology.guessers.guess_atom_element(name) for name in self.u.atoms.names])
        return elements

    def get_atom_charge(self):
        ''' Retruns atom's charge retrieved from data file. '''
        try:
            charges = self.u.atoms.charges
        except:
            charges = np.zeros( self.natoms, dtype=np.float64)
        return charges

    def get_molecule_name(self):
        ''' Return the molecues residue names. '''
        names = super(MDAnalysisReader, self).get_molecule_name()
        if len(names) == 0:
            _where = np.unique(self.u.atoms.molnums, return_index=True)
            names = self.u.atoms.moltypes[_where[1]]
        return names
        # return self.u.residues.resnames

    def get_bonds(self):
        ''' Retruns bonds' atoms indexes. '''
        try:
            bonds = self.u.bonds.indices
        except:
            bonds = np.array((),dtype=np.int32)
        return bonds
    
    def get_bond_type(self):
        ''' Retruns bonds' type. '''
        try:
            bpnd_type = self.u.bonds.type
        except:
            bpnd_type = np.array((),dtype=np.int32)
        return bpnd_type        

    def get_topology(self):
        ''' Return the system topology (see McReader) '''
        msg = '%s.%s not yet supported.' % (self.__class__.__name__,inspect.currentframe().f_code.co_name)
        raise MDAnalysisReaderException( msg)

    def count_frames( self, start=-1, end=sys.maxsize):
        ''' Return the number of configuration in the dump file. '''
        # implement a "hard count"
        icnt = 0
        if self.u:
            # count ...
            while True:
                ts = self._next_frame()
                if  ts == None:
                    break
                istep = ts.data['step']+self.istepoffset
                if istep >= start:
                    if istep <= end:
                        icnt+=1
                    else:
                        break
            # and reset
            self._next_file = self.__next_file()
            self.trajfile = next(self._next_file) # get the first file (always available)
            self._prepare()
            self.timeoffset = 0.0
            self.istepoffset = 0

        return  icnt

    def print_simulation_info(self):
        ''' Print basic simulation info. '''

        _tprfile = self.topofile if len(self.topofile) > 0 else '-'
        _trajfile = self.trajfile if len(self.trajfile) > 0 else "-"
        try:
            _dt = str(self.u.trajectory.dt)
        except:
            _dt = "-"

        print("-"*20)
        print("engine        : %s" % "gromacs")
        print("tpr file      : %s" % _tprfile)
        print("trj file      : %s" % _trajfile)
        print("# of atoms    : %d" % self.natoms)
        print("trj dt (%s)   : %s" % ("fs", _dt))
        print("-"*20)

    def _next_frame(self):
        ''' Return the next frame of the trajectorie(s) handled form the reader. '''
        ts = next(self._next_timestep)         # try to get the next frame
        if ts == None:                         # on fail
            filename = next(self._next_file)   # chek if there is a next trj file
            if not filename == None:           # if yes the
                self.trajfile = filename       # prepare it for reading
                # keep the last time-step
                lastframe = self.u.trajectory[-1]
                tlast = self.timeoffset + lastframe.time
                ilast = self.istepoffset + lastframe.data['step']
                self._prepare()
                # if the timestep is correct reset timeoffset
                firstframe = self.u.trajectory[0]
                tfirst = firstframe.time
                ifirst = firstframe.data['step']
                # assuming the same dump frequency with previous trajectory file
                self.timeoffset = tlast - tfirst
                self.istepoffset = ilast - ifirst

                ts = next(self._next_timestep) # skip the first frame
                ts = next(self._next_timestep) # and return the second
        return ts

    def __next_timestep(self):
        ''' Return the next trajectory timestep/frame from the current MDAnalysis univerce. '''
        for ts in self.u.trajectory:
            yield ts
        yield None

    def _find_files( self):
        ''' Find file with self.file+'n' for n = 1,2,... in the same path,
            adds them in self.files and terminates this list with None. '''
        i = 1
        while True:
            # name = self.filename + ".%d" % i
            name = self.dir + os.sep + self.basename + ".%d" % i + self.fileext
            if not os.path.isfile( name): break
            i += 1
            self.files.append( name)
        self.files.append( None)

    def __next_file(self):
        ''' Get the next trajectory file in self.files. '''
        for file in self.files:
            yield file

    def _prepare(self):
        ''' Scan the topology file to retrieve the number of atoms.
            Use the first configuration and assume that all
            configurations have the same number of atoms.
            Raises:
                MDAnalysisReaderException: if self.u fails to initialize. '''
        try:
            self.u = MDAnalysis.Universe( self.topofile, self.trajfile)
            delta = self.u.trajectory.dt # time between successive dumps
            # time of the first and secod dump
            t0, t1 = self.u.trajectory[0].data['step'], self.u.trajectory[1].data['step']
            self.timestep = delta/(t1-t0) # timestep in ps
        except:
            self.u = None
            raise MDAnalysisReaderException('Problem initializing MDAnalysis universe.')
        self._next_timestep = self.__next_timestep()

    def get_unwrapped(self):
        ''' Return the unwrap coordinates of the current frame. '''
        return self.u.atoms.unwrap(reference=None)

    @classmethod
    def create(cls, filename, topo):
        ''' Use the trajectory and the topology file names to create
            a MDAnalysisReader object.
            Args:
                filename (str): the trajectory (.trj or .xtc) file path.
                topo (str): the topology (.tpr) file path
        '''
        # keep input
        _filename = os.path.expanduser( filename)
        _topo = topo

        try:
            reader = cls()
            reader.set_files( _filename, _topo)
        except MDAnalysisReaderException as err:
            reader = None
            print(err)

        reader.read_topology()

        return reader

    @staticmethod
    def read_gro(f):
        ''' Parse gro file.
            Args:
                f: a valid file object to read
            Returns:
                list: a list of tuples with atoms info
                    [(resid, resname, atomname, atomid, x, y, z, vx, vy, vz)]
                    Gromacs units are assumed for coordinates and velocities.
        '''
        lines = f.readlines()
        atoms = []
        va = np.zeros(3, dtype=float)
        vb = np.zeros(3, dtype=float)
        vc = np.zeros(3, dtype=float)        
        for line in lines[2:2+int(lines[1].strip())]:
            resid = int(line[:5])
            resname = line[5:10]
            atomname = line[10:15]
            atomid = int(line[15:20])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            if len(line) > 45:
                vx = float(line[44:52])
                vy = float(line[52:60])
                vz = float(line[60:68])
            else:
                vx = vy = vz = 0.0
            atoms.append((resid, resname, atomname, atomid, x, y, z, vx, vy, vz))
        # get the box
        tk = lines[-1].strip().split()
        nbox = len(tk)
        if nbox == 9:
            va[0] = float(tk[0])
            va[1] = float(tk[3])
            va[2] = float(tk[4])
            vb[0] = float(tk[5])
            vb[1] = float(tk[1])
            vb[2] = float(tk[6])
            vc[0] = float(tk[7])
            vc[1] = float(tk[8])
            vc[2] = float(tk[2])
        elif nbox == 3:
            va[0] = float(tk[0])
            vb[1] = float(tk[1])
            vc[2] = float(tk[2])
        else:
            msg = "ERROR invalid box check the last line of the input gro file."
            raise MDAnalysisReaderException( msg)

        return atoms, (va,vb,vc)

def main():
    parser = MDAnalysis.core._get_readers.get_parser_for('/Users/loukas.peristeras/linux/files/kerogens/ucl_sructure/gromacs.methane/nvt1.tpr',format=None)

    reader = MDAnalysisReader()
    # reader.set_files('/Users/loukas.peristeras/work/demokritos/documents/paper.hydrates/loukas-P2000-T295-x05/pr2.trr')
    # reader.set_files('/Users/loukas.peristeras/linux/files/kerogens/ucl_sructure/gromacs.methane/nvt.trr')
    reader.set_files('/Users/loukas.peristeras/linux/files/kerogens/ucl_sructure/gromacs.methane/nvt1.trr')

    # print reader.count_frames()

    if reader.set_attributes('id x y z type el'):
        reader.read_topology()
        # i,b,d=reader.read_next_frame()
        # i,b,d=reader.read_next_frame()
        while True:
            i,b,d=reader.read_next_frame()
            if b == None:
                break
    else:
        print('Problem in setting GroReader attributes.')
