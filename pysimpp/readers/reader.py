# -*- coding: utf-8 -*-

import abc
import os
import yaml
import itertools

class ReaderException(Exception): ...


class abcReader(metaclass=abc.ABCMeta):
    ''' Prototype for reader class.
        Attributes:
            dir (str): the directory contains the simulation files. Usually, simulation
                files are divided in three categories:
                    a) input parameters, where the futures needed for the simulation
                    are provided,
                    b) system topology, where the topology of the system is provided
                    together with the force field parameters for the interactions
                    developed in the syste,
                    c) output, where thermodynamic data, system configurations (trajectory
                    frames) and simulation info are recorded during the simulations.

            yaml {}: if the file system.yaml exists in the simulation directory, it will
                be used to retrieve extra information for the system. This is a general
                mechanism available for all the derived classes, and it is useful in the
                the case where the simulation engine supported by the reader does not
                include this information in the simulation files. This attribute is the
                dict with the data of system.yaml file.

            timestep (int): simulation timestep.

            natoms (int): the number of particles in the system (atom is used in its literal
                sense i.e. a basic/unique entity)

            error (str): the error message from the last failed operation.
    '''

    @abc.abstractmethod
    def __init__(self):
        ''' Initialize Reader object. '''
        self.dir = ""
        self.natoms = 0
        self.timestep = 0.0
        self.yaml = None
        self.do_unwrap = False
        self.do_whole = False
        self.error = ""

    @abc.abstractmethod
    def set_files(self, filename, topo=None):
        ''' Set the simulation files.
            Args:
                filename (str): the name of a simulation file or of the
                    simulation directory. The reader should try to deduce
                    all the simulation files (input parameters, trajectory,
                    topology, log) and trace them. If no usable input/og
                    file is available then the simulation trajectory file
                    should be provided as argument.
                topo (str): the topology file. If provided, it will be used
                    and only on fail the reader will try to deduce it.

            Returns:
                bool: True if the needed files have been found and are valid.
                    and False otherwise.
            '''
        pass

    @staticmethod
    @abc.abstractmethod
    def is_supported(filename):
        ''' Check if the given file format is supported by the reader.
            Args:
                filename (str): the name of the file to check.
            Returns:
                (type (str), puropose (str)): the type of the file and
                    its purpose  ('log', 'inout', 'top' or 'trj') if the
                    file is not suppoted both strings will be empty.
        '''
        pass

    @abc.abstractmethod
    def set_attributes(self, attributes):
        ''' Set the attributes to be readed from the trajectory files. The
            attribute should be set prior to the call of self.read_next_frame
            method. Every reader should support a minimal set of attributes:
                ['id', 'type', 'x', 'y', 'z']
        '''
        pass

    @abc.abstractmethod
    def read_topology(self):
        ''' Read the topology of the system.
            Returns:
                bool: True on success.
        '''
        pass

    @abc.abstractmethod
    def read_next_frame(self, skip=False):
        ''' Read the next frame in the trajectory file.
            Returns:
                (int, SimulationBox, {attribute:data}): a tuple is the step
                    number,the simulation box and the dictionary with the
                    data of the next frame.
        '''
        pass

    @abc.abstractmethod
    def get_type_mass(self):
        ''' Retrun the {type (str):mass (float)} dictionary retrieved from the topology file. '''
        pass

    @abc.abstractmethod
    def get_atom_molecule(self):
        ''' Retruns atoms' molecules array retrieved from the topology file.
            Returns:
                np.array(..., dtype=int): atoms' molecules array.
        '''
        pass

    @abc.abstractmethod
    def get_atom_mass(self):
        ''' Retruns atoms' masses retrieved from the topology  file.
            Returns:
                np.array(..., dtype=float): atoms' masses array.
        '''
        pass

    @abc.abstractmethod
    def get_atom_type(self):
        ''' Retruns atoms' types retrieved from the topology  file.
            Returns:
                np.array(..., dtype=str): atoms' types array.
        '''
        return self._get_yaml_atom_property('atomtype')

    @abc.abstractmethod
    def get_atom_name(self):
        ''' Retruns atoms' names retrieved from the topology  file.
            Returns:
                np.array(..., dtype=str): atoms' names array.
        '''
        return self._get_yaml_atom_property('atomname')

    @abc.abstractmethod
    def get_atom_element(self):
        ''' Retruns atoms' elements retrieved from the topology  file.
            Returns:
                np.array(..., dtype=str): atoms' elements array.
        '''
        return self._get_yaml_atom_property('atomelement')

    @abc.abstractmethod
    def get_atom_charge(self):
        ''' Retruns atoms' charges retrieved from the topology  file.
            Returns:
                np.array(..., dtype=float): atoms' charges array.
        '''
        pass

    @abc.abstractmethod
    def get_bonds(self):
        ''' Retruns bonds' atoms indexes.
            Returns:
                np.array(..., dtype=int): bonds' atoms indexes (pair) array.
        '''
        pass

    @abc.abstractmethod
    def get_molecule_name(self):
        ''' Retruns molecule names retrieved from the topology  file.
            Returns:
                np.array(..., dtype=str): molecule names array.
        '''
        return self._get_yaml_atom_property('name')

    # @abc.abstractmethod
    # def get_topology(self):
    #     pass

    # @abc.abstractmethod
    # def get_forcefield(self):
    #     pass

    def set_unwrap(self, flag):
        ''' Set reader's unwrap flag.
            Args:
                flag(bool): if true the data dict returned form the read_next_frame
                    method reader will include the unwrapped coordinates of the frame
                    in the keys 'xu', 'yu', and 'zu'.
        '''
        self.do_unwrap = flag

    def set_whole(self, flag):
        ''' Set reader's whole flag.
            Args:
                flag(bool): if true the data dict returned form the read_next_frame
                    method reader will include the whole coordinates of the frame 
                    in the keys 'xw', 'yw', and 'zw'.
        '''
        self.do_whole = flag

    @abc.abstractmethod
    def count_frames(self):
        ''' Return the number of frames in the simulation trajectory.
            Returns:
                int: the number of frames.
        '''
        pass

    # @abc.abstractmethod
    # def print_simulation_info(self):
    #     ''' Print basic simulation info. '''
    #     pass

    def check_error(self, reset=False):
        ''' Check the reader error status.
            Args:
                reset (bool): if true the error will be cleared.
            Returns:
                bool : True if the error message has been set and
                    is non empty.
        '''
        ok = not len(self.error) == 0
        if reset:
            self.error = ""
        return ok

    def _read_yaml( self):
        ''' Check for the system.yaml file in the simulation directory.
            If exist read the extra data for the system.
            This file contains a dict with the folowing keys:
                ['name', 'comment', 'natoms', 'nmolecules', 'nspecie', 'specie']
                name: system name (str)
                comment: a comment (str)
                natoms: # of atoms (int)
                nmolecules: # ofmolecules (int)
                nspecie: # number of molecular types (int)
                specie: list of molecular types (list)
                atomtypes: list of the atom types name (list)
                bondtypes: list of the bond types name (list)
                angletypes: list of the bond types name (list)
                dihedraltypes: list of the bond types name (list)
                impropertypes: list of the bond types name (list)
            The 'specie' value is the list with the molecular types in the system.
            For each type the a dict with the following keys is provided:
                ['name', 'order', 'nmolecules', 'natoms', 'atomname', 'atomtype', 'atomres', ''atomelement']
                name: specie name
                order: the order of the listing in the data file
                nmolecules: # of molecules for the molecular type
                natom: # of atoms for the molecular type
                atomname: the name of the atoms
                atomtype: the type of the atoms
                atomres: the resname of the atoms
                atomelement: the element of the atoms
        '''
        fname = self.dir+os.sep+"system.yaml"
        if os.path.isfile(fname):
            with ( open(fname,'r')) as f:
                print("INFO: set reader info (%s)" % fname)
                self.yaml = yaml.load(f, Loader=yaml.FullLoader)

    def _get_yaml_atom_property( self, name):
        ''' Returns a per atom property retrieved from the yaml file.
            Args:
                name (str) : propery name.
            Returns:
                list (str): a list with propety for each atom. '''
        peratom = []
        try:
            if not self.yaml is None:
                _peratom = itertools.chain.from_iterable([
                        s[name] * s['nmolecules'] for s in self.yaml['specie'] ])
                peratom = list( _peratom)
        except:
            pass
        return peratom

    def _set_files_chk(self, filename, topo=None):
        ''' Check the set_files arguments and extract file(s) info. '''

        # check
        if not os.path.isfile( filename):
            message = '%s file not found.'%filename.strip()
            raise ReaderException(message)
        if not topo is None and not os.path.isfile(topo):
            message = '%s file not found.' % topo
            raise ReaderException(message)

        # extact path info
        ( self.dir, self.shortfilename ) = os.path.split( filename)
        ( self.basename, self.fileext ) = os.path.splitext( self.shortfilename)

        if len( self.dir) == 0:
            self.dir = "."

