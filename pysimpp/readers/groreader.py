# -*- coding: utf-8 -*-

import os
import sys
import re
import inspect
import numpy as np
from itertools import islice

import MDAnalysis

from .reader import abcReader, ReaderException # pylint: disable=import-error
from .mdanalysisreader import MDAnalysisReader, MDAnalysisReaderException
from pysimpp.utils.simulationbox import SimulationBox

_inform = True

class GromacsReaderException(Exception): ...


class GromacsReader(MDAnalysisReader):
    ''' Implements a gromacs reader from MDAnalysisReader by reimplementing
        the set_files and is_supported methods. '''

    def _set_files_policy(self, filename, topo):
        kind, puropose = GromacsReader.is_supported(filename)
        if puropose == 'trj':
            self.trajfile = filename
            if not topo is None:
                self.topofile = topo
            else:
                self.topofile = self.dir + os.sep + self.basename+".tpr"
                if not os.path.isfile( self.topofile):
                    print("WARNING: no type information will be abailable.")
                    self.topofile = self.dir + os.sep + self.basename+".pdb"
                    if not os.path.isfile( self.topofile):
                        self.topofile = self.dir + os.sep + self.basename+".gro"
                        if not os.path.isfile( self.topofile):
                            message = 'No valid topology file found for %s' % filename.strip()
                            raise GromacsReaderException( message)

            print('INFO: set reader (%s)'%self.trajfile)
            print('INFO: set reader topology (%s)' % self.topofile)
        else:
            message = 'The provided trajectory file is not valid %s' % filename.strip()
            raise ReaderException("message")

    @staticmethod
    def is_supported(filename):
        _, extension = os.path.splitext( filename)
        if extension in ('.trr', '.xtc'):
            s = (extension[1:],'trj')
        elif extension in ('.tpr'):
            s = ('tpr','top')
        elif extension in ('.gro'):
            s = ('gro','inout')
        else:
            s = ('','')
        return s
