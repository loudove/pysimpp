from .mdanalysisreader import MDAnalysisReader # pylint: disable=import-error
from .groreader import GromacsReader # pylint: disable=import-error
from .lammpsreader import LammpsReader # pylint: disable=import-error

def _is_command(): return False

# TODO automaticaly add the readers and update the file check
def create(filename, topo=None):
    ''' Create and return a reader for accessing the given
        situation trajectory file. T
        Args:
            filename (str): the full path to the trajectory file.
                For lammps, the full path to the simulation log file can
                also be provided. In this case, all the simulation details
                will be resolved from there.
            topo (str): the full path to the simulation topology file. For
                lammps and gromacs, if the topology file is not provided,
                the reader will check if a valid topology file resides in
                the same directory with the trajectory file having the same
                basename and if yes it will be used. 
    '''
    from .mdanalysisreader import MDAnalysisReaderException
    CHKREADER = {
        'lammps':LammpsReader.is_supported,
        'gromacs':GromacsReader.is_supported,
        'mdanalysis':MDAnalysisReader.is_supported
        }
    reader = None
    if any(CHKREADER['lammps'](filename)):
        reader = LammpsReader.create(filename,topo)
    elif any(CHKREADER['gromacs'](filename)):
        reader = GromacsReader.create(filename,topo)
    elif not topo is None:
        reader = MDAnalysisReader.create(filename,topo)
    return reader
