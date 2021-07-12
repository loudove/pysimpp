# -*- coding: utf-8 -*-

import os
import threading
from collections import defaultdict, Counter

import numpy as np
import networkx as nx

_debug = True

class Node(int):
    ''' Implements a node as an integer equipped with extra properties
        i.e. the connections of the node and cluster where it belongs.
        These properties can be provided in the constructor of the object
        as "connections" and "cluster" keyword arguments.
    '''

    def __new__(cls, *args, **kwargs):
        ''' Create and return a new Node object.
        Args:
            *args: argument list for int construction.
            **kwargs: keyword arguments dict including
                connections: a list with the connections of the node
                cluster (Cluster): the cluster of the node
        Returns:
            Node: the constructed Node object'''
        if 'base' in kwargs:
            obj =  super(Node, cls).__new__(cls, *args, base=kwargs['base'])
        else:
            obj =  super(Node, cls).__new__(cls, *args)
        obj.connections = kwargs['connections'] if 'connections' in kwargs else []
        obj.cluster = kwargs['cluster'] if 'cluster' in kwargs else None
        obj.kind = kwargs['kind'] if 'kind' in kwargs else None
        return obj

class Connection(tuple):
    ''' Implements an undirected edge as a symmetric pair equipped with
        the periodic images of the corresponding spatial connection. '''

    def __new__(cls, i, j, pbc=(0,0,0)):
        ''' Create and return a new Connection object. The hash of
            the object will be the hash of the ordered (i,j) tuple.
        Args:
            i (Node): the first node
            j (Node): the second node
            pbc: a sequence of three int values corresponding to the
                periodic indexes of the connection due to the PBCs
                (periodic boundary conditions) imposed spatially.
        '''
        obj =  super(Connection, cls).__new__(cls, (i,j,))
        obj._hash = hash((i,j,)) if i < j else hash((j,i,))
        obj.ip = np.array(pbc)
        return obj

    def __eq__(self, other):
        ''' Return self==other. '''
        return self._hash == other._hash if type(other) is type(self) else False # pylint: disable=no-member

    def __hash__(self):
        ''' Return hash(self) '''
        return self._hash # pylint: disable=no-member

    def is_periodic(self):
        ''' Check if the connection is periodic i.e. crosses the boundaries.
        Returns:
            bool: true if any of the periodic indexes is non zero.
        '''
        return self.ip.any() # pylint: disable=no-member

class Cluster(set):
    ''' Implements a Cluster i.e. an ensemble of connected nodes. The
        connections of the nodes are inherited to the cluster. The non
        periodic nodes' connections define the cluster connections with
        other clusters.

        A cluster can represent a part of a molecule, a pore or an other
        3d object which, due to the periodic boundary conditions (PBCs),
        is fragmented (spatially interrupted).

        The class provides the functionality for join a set of clusters
        connected through PBCs, making their union (parent cluster) whole.

        Attributes:
            id: cluster global index used also as its hash
            connections: a list with the non-periodic connections of
                the constituent nodes
            periodic_connections: a list with the periodic connections
                of the constituent nodes
            connectedto: {Cluster:{(ix,iy,iz):n}} a dict with keys the
                connected clusters and  values the corresponing connection
                periodic info. The periodic info is itself a dictionaly
                whith keys the periodic index of the connections and values
                the number of the connections of the specific periodic index.
            shift: a squence with the translation indexes i.e. the number
                of boundary edges vectors to translate the cluster in order
                to recontruct its parent cluster.
            parent: the parent cluster
        '''

    def __init__( self, i, s):
        ''' Initialize a Cluster object.
        Args:
            i (int): cluster id. It will be used also as hash for the object.
            s : a container of Node objects.
        '''
        super(Cluster, self).__init__( s)
        self.i = i
        self.connections = []
        self.periodic_connections = []
        self.connectedto = defaultdict( lambda :defaultdict(int))
        self.shift = np.zeros(3,dtype=np.int32)
        self.parent = None

    def __str__( self):
        ''' Return str(self). '''
        return "%s %d\n nodes: %s" % (self.__class__.__name__, self.i, " ".join(map(str,self)) )

    def __repr__( self):
        ''' Return repr(self). '''
        return "%s %d" % (self.__class__.__name__, self.i)

    def __eq__( self, other):
        ''' Return self==other. '''
        return self.i == other.i if type(other) is type(self) else False

    def __hash__( self):
        ''' Return hash(self) '''
        return self.i

    @classmethod
    def create(cls, i, nodes):
        ''' Returns a Cluster object.
        Args:
            i: the cluster global id to be used also as hash
            node: a sequence with the nodes of the cluster
        '''
        obj = cls(i, nodes)
        for n in obj:
            n.cluster = obj
        return obj

    def set_parent( self, parent):
        ''' Set the parent cluster. '''
        self.parent = parent

    def is_infinit( self):
        ''' Check if the cluster is infinite i.e. is connected with it self through PBCs.
            Returns:
                bool: True if the cluster is infinit and False otherwise. '''
        return self in self.connectedto

    def update_connect( self):
        ''' Update the connections of the nodes and the
            connections with the other clusters. '''
        self.__update_nodes_connections()
        self.__update_cluster_connections()

    def __update_nodes_connections( self):
        ''' Update the connections of the nodes. '''
        connections=self.connections
        periodic_connections=self.periodic_connections
        connections.clear()
        periodic_connections.clear()

        # process the connections for each node
        for n in self:
            n.cluster = self
            _periodic = [ c.is_periodic() for c in n.connections ]
            # keep unique connections
            connections += [ c for pbc, c in zip(_periodic,n.connections) if not pbc and c[0] < c[1] ]
            # keep multiple connections (different in their periodic index)
            periodic_connections += [ c for pbc, c in zip(_periodic,n.connections) if pbc ]

    def __update_cluster_connections( self):
        ''' Update the connections with the other clusters. '''
        connectedto = self.connectedto
        connectedto.clear()
        for c in self.periodic_connections:
            connectedto[ c[1].cluster][tuple(c.ip)] += 1

    def __update_shift_self( self, parent, shift):
        ''' Update self cluster after joining the parent cluster.
            Essensially the nodes (connections and periodic_connections) and cluster
            (connectedto) connections are updated.
            Args:
                parent: the cluster with which self cluster has been joined.
                shift: the sequence with the translation indexes of self
                    cluster for joining with parent. '''

        self.parent = parent

        # update nodes connections:
        # turn periodic to non-periodic and update connection shift
        periodic_connections = self.periodic_connections
        connections = self.connections
        # update connection shift
        for c in periodic_connections:
            c.ip -= shift
        # check if are still periodic. if not move them to non-periodic
        # connections list
        toremove = [ c for c in periodic_connections if not c.is_periodic() ]
        while len(toremove) > 0:
            c = toremove.pop(0)
            periodic_connections.remove(c)
            connections.append(c)
        # update connectedto list:
        connectedto = self.connectedto
        for k, v in connectedto.items():
            # if the connection is periodic update the shift otherwise do nothing.
            _tmp = { tuple( np.array(_k)-shift) if any(_k) else _k : _v for _k, _v in v.items() }
            v.clear()
            v.update( _tmp)
        self.shift += shift

        # now update the shift of self cluster in the connected clusters
        # (for both for nodes and clusters connections).
        for _c in connectedto.keys():
            if not _c is self:
                _c.__update_shift_other(self,shift)

    def __update_shift_other( self, other, shift):
        ''' Update the shift of the other cluster in the nodes (connections and
            periodic_connections) and cluster (connectedto) connections lists of
            self cluster.
            Args:
                other: the clustes to update
                shift: the sequence with the translation indexes of self
                    cluster for joining with parent.
        '''

        # update nodes connections:
        # turn periodic to non-periodic and update connection shift
        periodic_connections = self.periodic_connections
        connections = self.connections
        # update connection shift, check if periodic and if not move it
        # to the non-periodic connections list
        toremove = []
        for c in periodic_connections:
            if c[1].cluster is other:
                c.ip += shift
                if not c.is_periodic():
                    toremove.append(c)
        while len(toremove) > 0:
            c = toremove.pop(0)
            periodic_connections.remove(c)
            connections.append(c)
        # update connectedto list:
        v = self.connectedto[other]
        # if the connection is periodic update the shift otherwise do nothing.
        _tmp = { tuple( np.array(_k)-shift) if any(_k) else _k : _v for _k, _v in v.items() }
        v.clear()
        v.update( _tmp)

    @staticmethod
    def whole( connected):
        ''' Update the given clusters to make their union whole again.
            After returning the connections lists (nodes and cluster)
            of the clusters and their shift have been updated.
            Args:
                connected: a list with connected clusters.
            Returns:
                list: the connected lidy sorted based on the size of
                    the clusters in descending order.
        '''
        # copy the list and in any case sort it
        connected_sorted = sorted( connected, key=len, reverse=True)
        _connected = list(connected_sorted)
        # and set the first clusters as the parent
        parent = _connected[0]
        parent.parent = parent
        while len(_connected) > 0:
            c0 = None
            # retrieve the first subcluster already joined in
            for i, _c in enumerate(_connected):
                if _c.parent is parent:
                    c0 = _connected.pop(i)
                    break
            if c0 is None:
                if _debug:
                    print("** the provided list of clusters to be joined are not connected!")
                return tuple()

            # join with the connected (by descending length order)
            shift = np.zeros(3,dtype=np.int32)
            for c1 in sorted(c0.connectedto, key=len, reverse=True):
                # if already joined continue
                if c1.parent is parent:
                    continue
                # select how c1 will be joined with c0 (the direction
                # with the larger number of periodic connections)
                directions = c0.connectedto[ c1]
                _directions = sorted( directions.items(), key=lambda x: x[1], reverse=True)
                ip = _directions[0][0]
                shift[:] = -np.array(ip)

                # update c1 and the clusters connected to it (including c0)
                c1.__update_shift_self( parent, shift)
        return connected_sorted

class Assembly(set):
    ''' Implements a collection of nodes suitable for tracing the evolution of
        a non connected cluster. To identify the specie of the assembly, the
        kind of the constituent node is used. '''

    def __init__( self, uid, s):
        ''' Initialize an Assembly object.
        Args:
            uid (int): cluster id. It will be used also as hash for the object.
            s : a container of Node like objects.
        '''
        super(Cluster, self).__init__( s)
        self.formula = Assembly.find_formula(s)
        self.uid = uid
        self.life = 0.0 # in number of configurations processed

    def __str__( self):
        ''' Return str(self). '''
        return "%s %d\n nodes: %s" % (self.__class__.__name__, self.uid, " ".join(map(str,self)) )

    def __repr__( self):
        ''' Return repr(self). '''
        return "%s %s %d" % (self.__class__.__name__, self.name, self.uid)

    def __eq__( self, other):
        ''' Return self==other. '''
        return self.uid == other.uid if type(other) is type(self) else False

    def __hash__( self):
        ''' Return hash(self) '''
        return self.uid

    @staticmethod
    def find_formula( assembly):
        try:
            _kinds = Counter( list( map( lambda x: x.kind, assembly)))
            _formula = ""
            for _k in sorted( _kinds):
                _formula += "%s_%s_" %(str(_k),str(_kinds[_k]))
        except AttributeError:
            _formula = "*_%s" % str(len(assembly))
        return _formula

    @staticmethod
    def find_node_assembly( assemblies):
        ''' Create the current state. The clusters list should be
            sorted based on cluster size (i.e. number of molecules). '''
        node_assembly = defaultdict(int)
        for _ic, _c  in enumerate(assemblies):
            for _m in _c:
                node_assembly[_m] = _ic+1
        return node_assembly

class EvolutionTracker():

    __uid = -1
    __lock = threading.Lock()

    def __get_uid(self):
        EvolutionTracker.__lock.acquire()
        EvolutionTracker.__uid += 1
        EvolutionTracker.__lock.release()
        return EvolutionTracker.__uid

    def __init__(self):
        # {assembly UID: [] }
        self.assemblies = {}
        # previous step configuration data
        self.previous_t = 0.0   # time
        self.previous_node_assembly = None  # {moleculeID:assemblyID}
        self.previous = None    # list of assemblies
        # current step configuration data
        self.current_t = 0.0
        self.current_node_assembly = None
        self.current = None

    def update_current(self, assemblies, t):

        # move current to previous
        self.previous = self.current
        self.previous_node_assembly = self.current_node_assembly
        self.previous_t = self.current_t

        # set current
        self.current = assemblies
        self.current_node_assembly = Assembly.find_node_assembly(assemblies)
        self.current_t = t

        if not self.previous is None:
            # track assemblies evolution from the previous to the current state. 
            # this evolution is modeled as a set if reactions where the assemblies
            # of the  previous state (reactants or left side assmblies) react to
            # give the assemblies of the current state (products or righ side assemblies).
            _assemblies = (self.previous, self.current)
            _node_assembly = (self.previous_node_assembly, self.current_node_assembly)
            _used = ( [False]*len(_assemblies[0]), [False]*len(_assemblies[1]))
            for _p in _assemblies[0]:
                _candidates = [_p]
                i, j = 0, 1
                _sides = ( [], [])
                while not len( _candidates) == 0:
                    _candidates = self._track(_candidates, _node_assembly[i], _assemblies[i], _node_assembly[j], _used[j], _sides[j])

        else:
            for _a in assemblies: # just set assemblies uid
                _a.uid = self.__get_uid()
                self.assemblies[_a.uid]._a


    def _track(self, ls_candidates, ls_node_assembly, rs_assemblies, rs_node_assembly, rs_used, rs_side):
        ''' Track the righ side (rs) assemblies containing nodes from the assemblies of the left side (ls). '''
        _candidates = []
        for _ls in ls_candidates:
            for pair in [ (k,v) for k, v in sorted( Counter([ rs_node_assembly[i] for i in _ls]).items(), key=lambda item: item[1]) if not k == 0 ]:
                _assembly = rs_assemblies[k]
                rs_side.append( _assembly)
                _candidates.append( _assembly)


class Reaction():
    ''' Implements a reaction. '''

    def __init__(self, reactants, products):
        ''' Initialize a reaction object. '''
        self.reactants = Counter(reactants) # stoichiometry dict for reactants {reactant:#}
        self.reactants_string = Reaction.__get_string(self.reactants)
        self.products = Counter(products)  # stoichiometry dict for products {products:#}
        self.products_string = Reaction.__get_string(self.reactants)
        self.steps = defaultdict( int)  # dict {step : number of reactions}

    @staticmethod
    def __get_string( side):
        ''' Return the string corespond to the given stoichiometry dict. '''
        return " + ".join( [ k if side[k] == 1 else "%d %s"% (side[k], k) for k in sorted( side.keys())])

    def __str__(self):
        ''' Return str(self). '''
        return "%s \u2192 %s" % (self.reactants_string, self.products_string)

    def __eq__(self, other):
        ''' Return self==other. '''
        return self.reactants_string == other.reactants_string and self.products_string == other.products_string
        # return self.reactants == other.reactants and self.products == other.products

    def is_reverse(self, other):
        ''' Return True if the other reaction is the reverse of self. '''
        # return self.reactants_string == other.products_string and self.products_string == other.reactants_string
        return self.reactants == other.products and self.products == other.reactants

    def add_step(self, step):
        ''' Update steps dict (record). '''
        self.steps[step] += 1
