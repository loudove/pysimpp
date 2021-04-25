# -*- coding: utf-8 -*-

#from sys import *
from numpy import *
from collections import *

import pysimpp.utils.utils as utils

_simdbg = False

def _printdb_molecule(atoms, mol, conf):
    ''' Create the file conf_mol and print atoms info.
        *atoms* : type [ AtomInfo], in
        *mol* : type int, indent **in**
        *conf* : type int, indent **in** '''
    if not _simdbg: return
    f = open(str(conf) + "_" + str(mol), 'w')
    l = sorted(atoms, key=lambda x: x.i)
    for i in l:
        f.write(str(i.i) + " " + str(i.type) + " " + str(i.bonded) + "\n")
    f.close()

class AtomInfo():
    ''' Keep various atom info imported for reaxff simulation bond file. '''
    def __init__(self):
        ''' Initialize a AtomInfo object. '''
        self.i = 0       # atom index
        self.type = 0    # atom type
        # TODO save memory check if needed
        #self.q = 0.0     # charge
        #self.abo = 0.0   # abo
        #self.nlp = 0.0   # nlp
        self._initialize()

    @classmethod
    def create(cls, i, type):
        ''' Initialize a AtomInfo object. '''
        cls = cls()
        cls.i = i        # atom index
        cls.type = type  # atom type
        # TODO save memory check if needed
        #cls.q = 0.0     # charge
        #cls.abo = 0.0   # abo
        #cls.nlp = 0.0   # nlp
        cls._initialize()
        return cls

    def _initialize(self):
        self.bonded = [] # bonded partners list
        self.orders = [] # list with bonde orders
        self.used = False


    def set(self, line):
        ''' Parse the line (a record in lamms bonds file) string and
            set the atom info. '''
        list = line.split()
        self.i = int(list[0])-1
        self.type = int(list[1])
        nb = int(list[2])
        n_ = 3
        bonded = []
        for i in range(nb):
            bonded.append(int(list[n_ + i])-1)
        n_ = 4 + nb
        orders = []
        for i in range(nb):
            orders.append(float(list[n_ + i]))

        for j, o in zip(bonded, orders):
            if o > BondsInfo._bondOrderCutoff:
                self.bonded.append(j)
                self.orders.append(o)
        # TODO save memory check if needed
        #n_ = 4+2*nb
        #self.abo = float(list[n_])
        #self.nlp = float(list[n_+1])
        #self.q = float(list[n_+2])

    def set_from_hbonds(self, hbonds):
        ''' set the  bonded list and orders from the given hbonds list.
            the indexing should start from zero '''
        d = defaultdict(int)
        for i in hbonds: d[i] += 1
        for k, v in list(d.items()):
            self.bonded.append(k)
            self.orders.append(v)

class BondsInfo():
    ''' Keep the bonds info for a configuration. '''
    _bondOrderCutoff = 0.35

    def __init__(self):
        ''' Initialize a BondsInfo object. '''
        self.bonds = []
        self.orders = []
        self.q = None

    def set(self, line):
        ''' Parse the line (a record in lammps bonds file) string and update
            the configuration bonds. '''
        list = line.split()
        i = int(list[0])-1
        type = int(list[1])
        nb = int(list[2])
        n_ = 3
        bonded = []
        for i in range(nb):
            bonded.append(int(list[n_ + i])-1)
        n_ = 4 + nb
        orders = []
        for i in range(nb):
            orders.append(float(list[n_ + i]))
        for (j, o) in zip(bonded, orders):
            if i < j and o > BondsInfo._bondOrderCutoff:
                self.bonds.append(i)
                self.bonds.append(j)
                self.orders.append(o)
        self.q = float(list[6 + 2 * nb])

class Molecule():
    ''' Defines a molecule. '''
    _max_natoms = 200

    def __init__(self):
        ''' Initialize a Molecule object. '''
        self.i = -1        # molecule index (in configuration)
        self.atoms = []    # list of atoms
        self.sformula = "" # syntactic formula
        self.idKeys = None
        self.type = None   # molecular type
        self.life = 1      # life time in steps
        self.lifetime = 0.0    # life time in fs

        self.atoms_set = None

    def _find_formula(self, types):
        ''' Find the syntactic formula. types is the dict
            { id:element } e.g. {1:"C"}. '''
        elements = []
        list = ["C", "H", "O", "N", "S", "Si", "Al", "F", "Fe", "Au"]
        for atom in self.atoms:
            e = types[atom.type]
            elements.append(e)
            if not e in list: list.append(e)
        self.sformula = ""
        n = 0
        sum = 0
        for e in list:
            n = elements.count(e)
            if n != 0:
                sum += n
                if n == 1:
                    self.sformula += e
                else:
                    self.sformula += e + str(n)
        if sum < len(self.atoms):
            self.sformula += "X" + str(len(self.atoms)-n)

    def _find_string(self, types):
        ''' Return a string where molecule's connectivity is
            listed in detail. '''
        elements = []
        for atom in self.atoms:
            elements.append(types[atom.type])
        s = ""
        for i, atom in enumerate(self.atoms):
            if i != 0: s += "-"
            s += elements[i] + str(i) + "|" + "("
            n = len(atom.bonded)
            for j, b in enumerate(atom.bonded):
                s += elements[b] + str(b) + "|"
                if n < j: s += ","
            s += ")"
        return s

    # recursive
    #def trace(self, info, i, types):
    #    ''' Create the list of atoms for the molecule based on the connectivity info.
    #    '''
    #    atom = info[i]
    #    if atom.used: return
    #    atom.used = True
    #    self.atoms.append(atom)
    #    for b in atom.bonded:
    #        self.trace(info, b, types)

    # non recursive
    def trace(self, info, i, types):
        ''' Create the list of atoms for the molecule based on
            the connectivity info. '''
        atom = info[i]
        atom.used = True
        list = deque()
        list.append(atom)
        while len(list) != 0:
            atom = list.popleft()
            self.atoms.append(atom)
            for b in atom.bonded:
                atom_ = info[b]
                if not atom_.used:
                    atom_.used = True
                    list.append(atom_)
        self._find_formula(types)
        self.atoms_set = set( self.atoms)

    def normalize(self):
        ''' Normalize atoms numbering using Morgan's breadth-first scheme. The atoms list
            will be reordered and the bonded list of atoms will be  converted from global
            to local index list. '''
        self._convert_bonded_to_local()
        indxs = self._morgan()
        self._rearange(indxs)

    def get_bonds_list(self):
        ''' Return the bonds and bonds orders lists. Local indexs are
            uses for atoms. '''
        bonds = []
        orders = []
        for i, atom in enumerate(self.atoms):
            for (b, o) in zip(atom.bonded, atom.orders):
                if atom.i < self.atoms[b].i:
                    bonds.append((i, b))
                    orders.append(o)
        return bonds, orders

    def convert_atoms_list(self):
        ''' Convert the atoms list form AtomInfo object list to
            atoms index list. '''
        a = zeros(len(self.atoms), dtype=int32)
        for i, atom in enumerate(self.atoms): a[i] = atom.i
        del self.atoms
        self.atoms = a

    def _morgan(self):
        ''' Retunr a list with atoms normal numbering using Morgan's
            breadth-first scheme. Local indexes should be used in atoms
            bonded lists. '''
        labels = self._morgan_labeling()
        indxs = self._bfs_morgan_numbering(labels)
        self._rearange(indxs)
        return indxs

    def _rearange(self, indxs):
        ''' Rearange atoms list bases on the provided indxs list. Local
            indexes should be used in atoms bonded lists. '''
        # keep index maps in order to correct bonded lists
        d = {}
        for i in range(len(indxs)): d[i] = indxs[i]-1
        utils.sort_connected_lists(indxs, self.atoms, self.idKeys, False)
        # correct atoms bonded list
        for atom in self.atoms:
            for i in range(len(atom.bonded)): atom.bonded[i] = d[atom.bonded[i]]

    def _convert_bonded_to_local(self):
        ''' Use local atoms index in the AtomInfo bonded list. '''
        hash = {}
        for i, atom in enumerate(self.atoms):
            hash[atom.i] = i
        for atom in self.atoms:
            for i, j in enumerate(atom.bonded):
                atom.bonded[i] = hash[j]

    def _morgan_labeling(self):
        ''' Apply Morgan algorithm and find atoms labels. Return the list
            (array) with the labels. The atoms (AtomInfo) bonded list should
            use local numbering. '''
        a = None
        ap = None
        natoms = len(self.atoms)
        # avoid overflow: LDP check potential problems
        if natoms > self._max_natoms:
            a = zeros(natoms, dtype=float32)
            ap = zeros(natoms, dtype=float32)
            for i in range(natoms): ap[i] = len(self.atoms[i].bonded) / 1000.0
        else:
            a = zeros(natoms, dtype=int64)
            ap = zeros(natoms, dtype=int64)
            for i in range(natoms): ap[i] = len(self.atoms[i].bonded)

        np = -1
        n0 = -1
        N = 4     # ensure that after Nth iteration the number of different labels will not change
        cont = True
        while cont:
            a[:] = 0
            for i, atom in enumerate(self.atoms):
                for j in atom.bonded: a[i] += ap[j]
            d = defaultdict(int)
            for i in a: d[i] += 1
            n = len(d)
            if n == n0 or n == np:
                N -= 1
                if N == 0: cont = False
            np = n0
            n0 = n
            ap[:] = a[:]
        return a

    def _bfs_morgan_numbering(self, labels):
        ''' Use Morgan labels (output of the _MorganLabeling method) and
            perform bfs atoms numbering. Return a list (array) with atoms
            normal numbering. '''
        a = zeros(len(self.atoms), dtype=int32)
        __chk = -2
        a[:] = __chk

        self.idKeys = self._create_idlist(False)
        i_ = labels.argmax()     # locate the max value TODO take care the case of multiplicity
        num = 1
        q = deque()
        q.append(i_)
        a[i_] = -1
        while len(q) != 0:
            i_ = q.popleft()
            a[i_] = num
            num += 1
            ids = []
            lbls = []
            indx = []
            for i in self.atoms[i_].bonded:
                if a[i] == __chk:
                    ids.append(self.idKeys[i])
                    lbls.append(labels[i])
                    indx.append(i)
            if len(ids) != 0:
                utils.sort_connected_lists_check_multiplicity(lbls, ids, indx, True)
                for i in indx:
                    q.append(i)
                    a[i] = -1
        return a

    @classmethod
    def are_the_same(self, m1, m2):
        return m1.type == m2.type and m1.atoms_set == m2.atoms_set
#        return m1.type == m2.type and \
#          set( map( lambda x: x.i, m1.atoms)) == set( map( lambda x: x.i, m2.atoms))

    def _create_idlist(self, secondLevel):
        ''' Return a list with unique id's of molecule atoms used in
            the tie case of Morgan numbering algorithm. '''
        a = zeros(len(self.atoms), dtype=int32)
        idList = []
        # first level atom id
        for i, atom in enumerate(self.atoms):
            a[i]  = atom.type              # less than 10 atom id's should exist check
            a[i] += len(atom.bonded) * 10  # less that 10 bonds per atom should exist
            #a[i] += sum(atom.orders) * 100 # the sum of the bond orders should be less than 10
        # second level atom id
        if secondLevel:
            for i, atom in enumerate(self.atoms):
                list = []
                for b in atom.bonded: list.append(a[b])
                list.sort()
                s = ""
                for j in list: s += str(j)
                idList.append(str(a[i]) + s)
        else:
            for i in a: idList.append(i)

        return idList

class Configuration():
    ''' Define a configuration. '''
    def __init__(self, iconf):
        ''' Initialize a Configuration object. '''
        self.i = iconf              # configuration number
        self.step = 0               # simulation step (used in life time calculation
        #self.atoms = 0              # list of atoms info do not keep ! it check memory !
        self.bonds = None           # BondsInfo object
        self.coordinates = None     # array of coordinates
        self.periodicOffsets = None # array of periodic offsets
        self.box = None             # simulation box
        self.atomTypes = None       # list of atom types
        self.bondTypes = None       # list of bond types
        self.molecularTypes = []    # list of molecular types
        self.molecules = []         # list of molecules
        self.atomMolecule = None    # molecule where the atom belong
        self.moleculeType = []      # list with molecule type
        self.reactionsInfo = []     # list with reactions info [ (reaction, (:) reactunts molecules indexes, (:) products molecules indexes) ]

    def track_molecules(self, info, types):
        ''' Identify the molecules exist in the configuration using the
            connectivity information. '''
        self.molecules = []
        self.atomMolecule = zeros(len(info), dtype=int32)
        n = 0
        for i, atom in enumerate(info):
            if atom.used: continue
            m = Molecule()
            m.trace(info, i, types)
            self.molecules.append(m)
            m.i = n
            for at in m.atoms: self.atomMolecule[at.i] = n
            n += 1

    def set_step(self, step):
        ''' Set the step fo the configuration. '''
        self.step = step

    def remove_pbc(self, m):
        ''' Return a list with atoms corrdinates where the PBC have
            been removed. '''
        r = zeros(len(m.atoms) * 3, dtype=float32)
        if self.coordinates == None or len(self.coordinates) == 0 or self.box == None: return r

        #for (i_, i) in enumerate(m.atoms):
        #    r[3*i_:3*i_+3] = self.coordinates[3*i:3*i+3]
        used = zeros(len(m.atoms) * 3, dtype=int32)
        i = m.atoms[0].i
        used[0] = 1
        r[0:3] = self.coordinates[3 * i :3 * i + 3]

        ri = zeros(3, dtype=float32)
        dr = zeros(3, dtype=float32)
        q = deque()
        q.append(0)
        while len(q) != 0:
            i = q.popleft()
            ri[0:3] = r[3 * i :3 * i + 3]
            for b in m.atoms[i].bonded:
                if used[b] == 0:
                    used[b] = 1
                    q.append(b)
                    j = m.atoms[b].i
                    dr[0:3] = self.coordinates[3 * j :3 * j + 3] - ri[0:3]
                    self.box.set_to_minimum(dr)
                    r[3 * b :3 * b + 3] = ri + dr[0:3]
        return r

    def track_molecular_types(self, morgan=True):
        ''' Track the molecular types (track_molecules
             method should be called first). '''
        for m in self.molecules:
            exist = False
            if morgan:
                m.normalize()
            else:
                m._convert_bonded_to_local()
            for mt in self.molecularTypes:
                exist = mt.match(m, morgan)
                if exist: break
            if not exist:
                mt = MolecularType(m, 1)
                self.molecularTypes.append(mt)
                mt.name = "molecule%d" % len(self.molecularTypes)
            m.type = mt
            # IMPORTANT: saves memory from AtomInfo objects
            m.convert_atoms_list()
            self.moleculeType.append(mt)

class MolecularType():
    ''' Define a molecular type. '''
    def __init__(self, m, nconf):
        ''' Initialize a MolecularType object. '''
        self.name = ""                # the assigned name
        self.string = ""

        self.sformula = m.sformula    # molecule syntactic formula
        self.idKeys = m.idKeys        # molecule unique key
        self.types = []               # atoms type
        for atom in m.atoms: self.types.append(atom.type)
        self.bonds, self.orders = m.get_bonds_list() # bonds and their orders
        self.r = None                 # atoms coordinates

        # trajectory log. A dict/list with { step : # of molecules of the molecular type}
        self.trajectory = zeros(nconf, dtype=int32)

        # reactions log.
        self.creation = []
        self.destruction = []
        # list with the lifetime for the molecules of this molecular type
        self.lifetimes = []

    def match(self, m, morgan=True):
        ''' Return True is the given Molecule m is of the same molecular type. '''
        ret = False
        if m.sformula == self.sformula:
            if morgan and array_equal(m.idKeys, self.idKeys):
                ret = True
            else:
                ret = True
        return ret

    def add(self, m, conf, iconf):
        self.trajectory[iconf] += 1
        if self.r == None and conf.coordinates != None and len(conf.coordinates) != 0:
            self.r = conf.remove_pbc(m)  # atoms coordinates

    def first_appeared(self):
        n = 0
        for a in self.trajectory:
            if a != 0: break
            n += 1
        return n

    def life_time_statistics(self):
        n = len(self.lifetimes)
        if n == 0:
            return 0.0, 0.0
        elif n < 3:
            return mean(self.lifetimes), 0.0
        else:
            return mean(self.lifetimes), std(self.lifetimes)

class Reaction():
    ''' Defines a reaction. '''
    def __init__(self):
        self.reactants = None               # list with reactants
        self.products = None                # list with products
        self.step = defaultdict(lambda: 0)  # dict {step : number of reactions}

    def add_step(self, s):
        self.step[s] += 1

    def create_trajectory(self, nconfs):
        ''' convert self.step dict to array self.trajectory. '''
        self.trajectory = zeros(nconfs, dtype=int32)
        for k, v in list(self.step.items()):
            self.trajectory[k] = v

    def get_reactants_types(self):
        d = defaultdict(lambda: 0)
        for r in self.reactants: d[r] += 1
        return list(d.keys())

    def get_products_types(self):
        d = defaultdict(lambda: 0)
        for p in self.products: d[p] += 1
        return list(d.keys())

    def get_reactants_string(self):
        return self._get_string(self.reactants)

    def get_products_string(self):
        return self._get_string(self.products)

    def get_string(self):
        r = self._get_string(self.reactants)
        p = self._get_string(self.products)
        s = "%s \u2192 %s" % (r, p)
        return s

    def count(self):
        return sum(self.step.values())

    def compare(self, r1, r2):
        return r1.reactants == r2.reactants and r1.products == r2.products

    def _get_string(self, list):
        d = defaultdict(lambda: 0)
        for p in list: d[p] += 1
        s = ""
        i = 0
        n = len(d)
        for k, v in list(d.items()):
            if v != 1: s += str(v) + " "
            s += "%s [%s]" % (k.sformula, k.name)
            i += 1
            if i != n: s += " + "
        return s

    def contains(self, list):
        items = self.get_reactants_types()
        for item in items:
            if item in list: return True
        items = self.get_products_types()
        for item in items:
            if item in list: return True
        return False

    def is_reverse(self, r):
        return self.reactants == r.products and self.products == r.reactants

class System():
    ''' Define a simulation system. '''
    def __init__(self):
        ''' Initialize a System object. '''
        self.configurations = []    # list of configurations
        self.molecularTypes = []    # list of molecular types
        self.reactions = []         # list of reactions
        self._previous = None       # previous conformation
        self._current = None        # current conformation
        self._previousInfo = None   # previous conformation info
        self._currentInfo = None    # current conformation info
        self._iconf = 0             # current configuration index (TODO remove this if configurations are saved in the self. configurations list)

    def next_conformation_(self, nconfs, iconf, info, types, dt):
        ''' Given the next info Identify the molecular
            type and track the reactions took place. '''

        self._current = Configuration( self._iconf)
        self._currentInfo = info

        # If the size of the list is zero, we have process all the configurations.
        if len(self._currentInfo) == 0:
            # In the last configuration add the life times of the
            # existig molecules wich have not been processed.
            del self._currentInfo
            if not self._previousInfo == None:
                del self._previousInfo
            return False

        self._current.set_step( iconf)
        # Track the molecules
        self._current.track_molecules(self._currentInfo, types)
        # Track the molecular type
        self._track_molecular_types( types, nconfs, self._current, self._iconf)
        # Track reactions
        self._track_reactions( 1.0)
        # Store the configuration
        self.configurations.append(self._current)

        # Keep the last two info lists (current and previous configuration). Free some memory
        if not self._previous == None:
            del self._previous.molecules
            del self._previous.molecularTypes

        self._previous = self._current
        if not self._previousInfo == None: del self._previousInfo
        self._previousInfo = self._currentInfo
        self._current = None
        self._currentInfo = None

        ##print "number of reactions : ", len(self.reactions)

        # Increase the index for conformations numbering
        self._iconf += 1
        return True

    def next_configuration(self, reader, keep=False):
        ''' Read the next configuration from lammps bonds file.
            Identify the molecular type and track the reactions
            took place. '''
        ####if self._iconf == 10:
        ####    r = Reaction()
        ####    r.reactants = [self.molecularTypes[0]]
        ####    r.products = [self.molecularTypes[1]]
        ####    r.add_step(self._iconf-1)
        ####    self.reactions.append(r)
        ####    return False

        # When the simulation still running the number of configurations
        # in the bonds file could differ from the number of configurations
        # in the dump or log file.
        if self._iconf >= len(reader.t_bonds): return True

        # Create the current configuration
        self._current = Configuration(self._iconf)
        # Set the bond info
        if reader.bondsInfo != None and len(reader.bondsInfo) > self._iconf:
            self._current.bonds = reader.bondsInfo[self._iconf]
        # Set the coordinates
        iconf_ = reader._bonds_to_coordinates_configuration(self._iconf)
        ##print self._iconf, " -> ", iconf_
        if iconf_ != -1:
            self._current.coordinates = reader.coordinates[iconf_] # TODO check that the len of the lists are the same
            self._current.box = reader.boxes[iconf_]
            self._current.periodicOffsets = reader.offsets[iconf_]

        # Get the list with atoms info form the reader. This list (of AtomInfo objects) is memory intensive
        self._currentInfo = reader.read_bonds_next(keep)
        # If the size of the list is zero, we have process all the configurations in the bonds file
        if len(self._currentInfo) == 0 or len(reader.t_bonds) == self._iconf + 1:
            # In the last configuration add the life times of the existig molecules
            # have not been processed. LDP TODO check the behavior
            # Clean up the memory and return
            del self._currentInfo
            if not self._previousInfo == None: del self._previousInfo
            return False
        self._current.set_step(reader.t_bonds[self._iconf])
        # Track the molecules
        self._current.track_molecules(self._currentInfo, reader.types)
        # Track the molecular type
        self._track_molecular_types(reader.types, reader.configurationsInBondsFile, self._current, self._iconf)
        # Track reactions
        self._track_reactions( reader.timestep)
        # Store the configuration
        self.configurations.append(self._current)

        # Keep the last two info lists (current and previous configuration). Free some memory
        #if not self._previous == None: del self._previous
        if not self._previous == None:
            del self._previous.molecules
            del self._previous.molecularTypes
        self._previous = self._current
        if not self._previousInfo == None: del self._previousInfo
        self._previousInfo = self._currentInfo
        self._current = None
        self._currentInfo = None

        ##print "number of reactions : ", len(self.reactions)

        # Increase the index for conformations numbering
        self._iconf += 1
        return True

    def _track_reactions(self, timestep):
        ''' Track the reactions. '''
        if self._previous == None: return

        used0 = zeros(len(self._previous.molecules), dtype=int32)
        used1 = zeros(len(self._current.molecules), dtype=int32)
        molecules0 = self._previous.molecules
        molecules1 = self._current.molecules
        atomMolecule0 = self._previous.atomMolecule
        atomMolecule1 = self._current.atomMolecule

        dt = (self._current.step - self._previous.step) * timestep
        halfdt = 0.5 * dt

        for rm in molecules0:
            if used0[rm.i] != 0: continue

            i = rm.atoms[0]                     # first atom of the pm molecule
            pm = molecules1[atomMolecule1[i]]   # using this atom spot the rm molecule
            #if rm.type == pm.type:              # rm and pm have the same type, no reaction
            if Molecule.are_the_same(rm, pm): # rm and pm have the same type, no reaction
                pm.life += rm.life
                pm.lifetime = rm.lifetime + dt
            else:                               # reaction took place, using the atoms of rm and pm
                rlist = []                      # track reactants and products
                plist = []
                candidates = [pm.i]
                while len(candidates) != 0:
                    candidates = self._add_products_track_reactants(atomMolecule0, candidates, used1, molecules1, plist)
                    candidates = self._add_reactants_track_products(atomMolecule1, candidates, used0, molecules0, rlist)

                reactants = []
                for r in rlist: reactants.append(r.type)
                products = []
                for p in plist: products.append(p.type)

                for r, t in zip(rlist, reactants):
                    r.lifetime += halfdt
                    t.lifetimes.append(r.lifetime)

                for p in plist:
                    p.lifetime += halfdt

                reaction = self._much_reaction(reactants, products)

                rindex = []
                for r in rlist: rindex.append(r.i)
                pindex = []
                for p in plist: pindex.append(p.i)
                self._current.reactionsInfo.append((reaction, tuple(rindex), tuple(pindex)))

    def _much_reaction(self, reactants, products):
        reactants.sort()
        products.sort()
        ret = None
        for r in self.reactions:
            if r.reactants == reactants and r.products == products:
                ret = r
                break
        if ret == None:
            ret = Reaction()
            ret.reactants = reactants
            ret.products = products
            self.reactions.append(ret)
            self._update_molecular_type_creation_destruction(ret)

        ret.add_step(self._iconf)
        return ret

    def _update_molecular_type_creation_destruction(self, reaction):
        reactants = reaction.get_reactants_types()
        for r in reactants:
            if not reaction in r.destruction: r.destruction.append(reaction)
        products = reaction.get_products_types()
        for p in products:
            if not reaction in p.creation: p.creation.append(reaction)

    def _add_products_track_reactants(self, atomMolecule0, candidates, used1, molecules1, plist):
        d = defaultdict(lambda: 0)
        for i in candidates:                                 # loop over candidate products
            cm = molecules1[i]
            if used1[i] == 0:
                used1[i] = 1
                plist.append(cm)                             # add the products in the products list
                for j in cm.atoms: d[atomMolecule0[j]] += 1  # locate the candidate reactants
        list_ = list(d.keys())
        return list_

    def _add_reactants_track_products(self, atomMolecule1, candidates, used0, molecules0, rlist):
        d = defaultdict(lambda: 0)
        for i in candidates:                                 # loop over candidate reactants
            cm = molecules0[i]
            if used0[i] == 0:
                used0[i] = 1
                rlist.append(cm)                             # add the reactants in the reactants list
                for j in cm.atoms: d[atomMolecule1[j]] += 1  # locate the candidate products
        list_ = list(d.keys())
        return list_

    def _track_molecular_types(self, types, N, conf, iconf):
        ''' Track the molecular types. '''
        for m in conf.molecules:
            exist = False
            m.normalize()
            for mt in self.molecularTypes:
                exist = mt.match(m)
                if exist: break
            if not exist:
                mt = MolecularType(m, N)
                self.molecularTypes.append(mt)
                mt.name = "molecule%d" % len(self.molecularTypes)
            m.type = mt
            mt.add(m, conf, iconf)
            mt.string = m._find_string( types)
            # IMPORTANT: saves memory from AtomInfo objects
            m.convert_atoms_list()
            conf.molecularTypes.append(mt)
            conf.moleculeType.append(mt)

    def check_life_times_info(self):
        # Correct life times
        lifetime = defaultdict(list) # {molecular type:life times list}

        class BasicMoleculeLifeInfo():
            def __init__(self, type):
                self.first = None
                self.type = type
                self.life = 1

        previous_react = None
        previous_molecules = None
        previous = None
        current_react = None
        current_molecules = None

        for current in self.configurations:

            current_molecules = []
            for i, type in enumerate(current.moleculeType):
                m = BasicMoleculeLifeInfo(type)
                current_molecules.append(m)
                m.first = where(current.atomMolecule == i)[0][0]
            current_react = empty(len(current.moleculeType), dtype=bool)
            current_react[:] = False

            if not previous is None:

                for reaction in current.reactionsInfo:
                    print(reaction)
                    #for t, r in zip(reaction[0].reactants, reaction[1]):
                    for r in reaction[1]:
                        #lifetime[t].append(previous_molecules[r].life)
                        previous_react[r] = True # pylint: disable=unsupported-assignment-operation
                    for p in reaction[2]:
                        current_react[p] = True
                for i, m in enumerate(previous_molecules):
                    if previous_react[i]:   # pylint: disable=unsubscriptable-object
                        lifetime[m.type].append(m.life)
                    else:
                        m_ = current_molecules[current.atomMolecule[m.first]]
                        m_.life += m.life
                for i, m in enumerate(current_molecules):
                    if current_react[i]:
                        if not m.life == 1: print("problem ...A")
                    else:
                        m_ = previous_molecules[previous.atomMolecule[m.first]] # pylint: disable=unsubscriptable-object
                        if not m.life - m_.life == 1: print("problem ...B")

            previous = current
            previous_molecules = current_molecules
            previous_react = current_react
            previous_react[:] = False

        for t in self.molecularTypes:
            print(t.sformula, "   ", t.lifetimes, "   ", lifetime[t])

    def eliminate(self, reader):
        ''' Eliminate molecular type B given that:
        a) B was created form a simple reactions of the form A -> B
        b) B was destroied from the reverse reaction B -> A
        c) B appears only onse in the simulations. '''

        eliminate_reactions = []     # list with reaction to be eliminated
        eliminate_types = []         # list with types to be eliminated
        replace = {}                 # { eliminated : parent type}
        lifetime = defaultdict(list) # {molecular type:life times list}

        # Spot the reactions/types to be eliminated
        for mt in self.molecularTypes:
        #    print "type ", mt.sformula, ", sum :", sum(mt.trajectory)
        #    print "creation: ", len(mt.creation), ", desrtuction :", len(mt.destruction)
        #    if len(mt.creation) == 1 and len(mt.destruction) == 1:
        #        print "creation : ", mt.creation[0].get_string()
        #        print "destruction : ", mt.destruction[0].get_string()
        #        print "are reverse : ", mt.creation[0].is_reverse(mt.destruction[0])
            if sum(mt.trajectory) == 1 and \
                len(mt.creation) == 1 and len(mt.destruction) == 1 and \
                len(mt.creation[0].reactants) == 1 and len(mt.creation[0].products) == 1 and  \
                mt.creation[0].is_reverse(mt.destruction[0]):
                    eliminate_reactions.append(mt.creation[0])
                    eliminate_reactions.append(mt.destruction[0])
                    eliminate_types.append(mt)
                    replace[mt] = mt.creation[0].reactants[0]


        # Structure for keeping molecule basic info. Help the life time updates
        class BasicMoleculeLifeInfo():
            def __init__(self, type):
                self.first = None
                self.type = type
                self.life = 1
                self.lifetime = 0.0

        # Correct life times
        previous_react = None     # (prevous|current) array indicating if the molecule participate in a reaction
        previous_molecules = None # (prevous|current) list of molecule
        previous = None           # (prevous|current) configuration
        current_react = None
        current_molecules = None

        for current in self.configurations:

            # create the list of molecule for the current configuration
            current_molecules = []
            for i, type in enumerate(current.moleculeType):
                m = BasicMoleculeLifeInfo(type)
                current_molecules.append(m)
                m.first = where(current.atomMolecule == i)[0][0]
            # create and initialize to False the reaction flag array
            current_react = empty(len(current.moleculeType), dtype=bool)
            current_react[:] = False

            if previous:
                dt = (current.step - previous.step) * reader.timestep
                halfdt = 0.5 * dt
                #print previous.step, current.step, dt, halfdt

                # process the reactions info list of the current configuration
                for reaction in current.reactionsInfo:
                    # use the non eliminated reactions in order to set the reaction flag arrays for the
                    # previous and current configuration
                    if not reaction[0] in eliminate_reactions:
                        for r in reaction[1]:
                            previous_react[r] = True # pylint: disable=unsupported-assignment-operation
                        for p in reaction[2]:
                            current_react[p] = True
                # process the molecules in the previous configuration
                for i, m in enumerate(previous_molecules):
                    # if the molecule reacts then update the life time of the corresponding type otherwise
                    # increate the life index for the same molecule in the current configuration
                    if previous_react[i]:  # pylint: disable=unsubscriptable-object
                        m.lifetime += halfdt
                        lifetime[m.type].append(m.lifetime)
                        #lifetime[m.type].append(m.life)
                    else:
                        m_ = current_molecules[current.atomMolecule[m.first]]
                        m_.life += m.life
                # process the molecules in the current configuration and check if the update was correct.
                for i, m in enumerate(current_molecules):
                    if current_react[i]:
                        m.lifetime += halfdt
                        if not m.life == 1:
                            print("problem ...A")
                            print(m.type.sformula)
                    else:
                        m_ = previous_molecules[previous.atomMolecule[m.first]] # pylint: disable=unsubscriptable-object
                        m.lifetime = m_.lifetime + dt
                        if not m.life - m_.life == 1:
                            print("problem ...B")
                            print(m.type.sformula)
                            print(m_.type.sformula)
            # prepare the next iteration (update the previous with the current)
            previous = current
            previous_molecules = current_molecules
            previous_react = current_react
            previous_react[:] = False

        # correct molecular types: update trajectory info
        for mt in eliminate_types:
            #print "correct trajectory : %s -> %s" % (mt.sformula, mt.creation[0].reactants[0].sformula)
            mt.creation[0].reactants[0].trajectory += mt.trajectory

        # correct molecular type: udpate creation/destruction reactions
        for mt in self.molecularTypes:
            newlist = []
            for r in mt.creation:
                if not r in eliminate_reactions: newlist.append(r)
            mt.creation = newlist
            newlist = []
            for r in mt.destruction:
                if not r in eliminate_reactions: newlist.append(r)
            mt.destruction = newlist

        # correct reactions list: remove eliminated reactions
        for r in eliminate_reactions:
            #print "remove reaction:", r.get_string()
            self.reactions.remove(r)
        # correct types list: remove eliminated types
        for mt in eliminate_types:
            #print "remove type: ", mt.sformula
            self.molecularTypes.remove(mt)

        # correct confguration: update moleculeTypes and reactionsInfo
        for configuration in self.configurations:
            for i, t in enumerate(configuration.moleculeType):
                if t in eliminate_types: configuration.moleculeType[i] = replace[t]
            newinfo = []
            for info in configuration.reactionsInfo:
                if not info[0] in eliminate_reactions: newinfo.append(info)
            configuration.reactionsInfo = newinfo

        # correct molecule types: update names and life times
        for (i, mt) in enumerate(self.molecularTypes):
            mt.name = "molecule%d" % (i + 1)
            mt.lifetimes = lifetime[mt]
