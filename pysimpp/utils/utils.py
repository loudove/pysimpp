# -*- coding: utf-8 -*-

import os
import sys
import inspect
import time
import math
import numpy as np
from collections import defaultdict
from argparse import ArgumentTypeError
# from numpy import *

def inspect_classes( module):
    ''' returns a {name:class} dict with the classes exist in the given module. '''
    classes = {}
    for name, obj in inspect.getmembers( sys.modules[ module]):
        if inspect.isclass(obj):
            classes[name] = obj
    return classes

def sequences(a, w0):
    ''' returns a tuple of pairs with the ranges of the contignious
        sequences of the array a with values greater or equal with w0.
        >>> sequences( [4, 4, 0, 0, 0, 4, 4, 4, 4, 0, 0, 0, 0, 4, 4, 4, 4, 0, 4], 4)
        '''
    _a = np.where( np.array(a) >= w0)[0]
    d = np.diff( _a)
    w = np.where( d > 1)[0] + 1
    s = np.insert( w, 0, 0)
    s = np.insert(s, s.size, len(_a))
    # for i in range(len(s)-1): print a[_a[s[i]:s[i+1]]]
    return tuple( [(_a[s[i]],_a[s[i+1]-1]) for i in range(len(s)-1)])

def isfile(filename):
    ''' Check if the given file exists. '''
    try:
        with open(filename,'r') as f:
            pass
    except OSError as e:
        args = {'filename': filename, 'error': e}
        message = ("can't open '%(filename)s': %(error)s")
        raise ArgumentTypeError(message % args)
    return filename

def isnumber(string, numbertype=int):
    '''Check if the given string is of type numbertype and if yes return the number'''
    try:
        val = numbertype(string)
    except:
        val = None
    return val

def ispositive( string, numbertype=int):
    '''Check if the given string is possitive number of numbertype and if yes return the number'''
    val = isnumber( string, numbertype)
    if val and val > numbertype(0):
        return val
    else:
        return None

# this is not a general functionality so remove it from here
def checkstart(line, name, sep="="):
    '''Check if the given string line starts with name and if yes return a list tk contains
       the parsing of line "name tk[0] sep tk[0]" e.g. "#declare r=1.0;" if sep is "=". If
       this is not the case returns an empty pair [None, None]. '''
    if line.startswith( name):
        tk = line[len(name):-1].split( sep)
        if len( tk) == 2:
            return tk
    return [None, None]

class wildstr(str):
    ''' implements a wild string i.e. what ever is compared to it is equal. '''
    wildchar = ['*', 'X']
    def __init__(self,value):
        #str.__init__(self)
        self.iswild = value in wildstr.wildchar

    def __eq__(self, other):
        if self.iswild or other.iswild:
            return True
        return super(wildstr,self).__eq__(other)

class IsList(object):
    ''' Implements list check functionality with custom setting
        and error handling. Usefull for argparse type check. '''

    def __init__(self, message, itemtype=int, positive=False, sep=','):
        ''' Initialize the object and store settings.
            Args:
                message (str) : the error message.
                itemtype (type) : the type of the items in the list.
                    Currently int, float and str types are supported.
                positive (bool) : check if the numbers in the list are
                    positive.
                sep (str) : the separator.
        '''
        self.message = message
        self.itemtype = itemtype
        self.positive = positive
        self.sep = sep

    def __call__(self, string):
        ''' Check if the string is a list conforms with the specifications. '''
        if len(string) == 0:
            return []
        if self.itemtype is str:
            items = list(map(lambda x: x.strip(), string.split( self.sep)))
        else:
            items = islist(string, numbertype=self.itemtype,
                positive=self.positive, sep=self.sep)
        if len(items) == 0:
            message = self.message % ('"%s"'%string) if "%s" in self.message else self.message
            raise ArgumentTypeError(message)
        return items

class IsListOfList(object):
    ''' Implements list of lists check functionality with custom setting
        and error handling. Usefull for argparse type check.
        A list of list has the form:
            A,B,C:D,E,G:H,J,L,R
    '''
    def __init__(self, message, itemtype=str, positive=False,
                 llen=-1, sep1=':', sep2=','):
        ''' Initialize the object and store settings.
            Args:
                message (str) : the error message if one of the lists fails the
                    specifications (type or positiveness) or it is empty.
                itemtype (type) : the type of the items in the lists. Currently
                    int, float and str types are supported.
                positive (bool) : check if the numbers in the lists are positive.
                llen (int) : the length of the lists. If a possitive integer is
                    provided all the lists shoud have the same length and equal to
                    llen
                sep1 (str) : the separator for the list.
                sep2 (str) : the separator for the lists.
        '''
        self.message = message
        self.itemtype = itemtype
        self.positive = positive
        self.llen = llen
        self.sep1 = sep1
        self.sep2 = sep2

    def __call__(self, string):
        ''' Check if the string conforms with the specifications. '''
        ret = []
        if len(string) == 0:  # empty argument
            raise ArgumentTypeError("argument is missing ")
        lines = string.split(self.sep1)
        for line in lines:
            if self.itemtype is str:
                items = list(map(lambda x: x.strip(), line.split(self.sep2)))
            else:
                items = islist(line,
                               numbertype=self.itemtype,
                               positive=self.positive,
                               sep=self.sep2)
            if len(items) == 0 or self.llen>0 and not len(items) == self.llen:
                message = self.message % ('"%s"' % line) if "%s" in self.message else self.message
                raise ArgumentTypeError(message)
            ret.append(items)
        return ret

class IsListOfNamedList(IsListOfList):
    ''' Implements list of list check functionality with custom setting
        and error handling. Usefull for argparse type check.
        A list of list has the form:
            N:M:A,B,C@O:P:D,E@G:X:H,J,L,R
    '''
    def __init__(self, message, itemtype=str, positive=False,
                 klen=-1, llen=-1, sep='@', sep1=':', sep2=',',
                 choices=(wildstr('*'),)):
        ''' Initialize the object and store settings.
            Args:
                message (str) : the error message if one of the lists fails the
                    specifications (type or positiveness) or it is empty.
                itemtype (type) : the type of the items in the lists. Currently
                    int, float and str types are supported.
                positive (bool) : check if the numbers in the lists are positive.
                klen (int) : the length of the name-list
                llen (int) : the length of the named lists. If a possitive integer is
                    provided all the lists shoud have the same length and equal to
                    llen
                sep (str) : the separator for the list.
                sep1 (str) : the separator for the name-list.
                sep2 (str) : the separator for the named lists.
        '''
        super(IsListOfNamedList, self).__init__(message, itemtype, positive, llen, sep1, sep2)
        self.sep = sep
        self.klen = klen
        self.choices = choices

    def __call__(self, string):
        ''' Check if the string conforms with the specifications. '''
        ret = {}
        if len(string) == 0:  # empty argument
            raise ArgumentTypeError("argument is missing ")
        lines = string.split(self.sep)
        for line in lines:
            tk = tuple( map(lambda x: x.strip(), line.split(self.sep1)))
            if self.klen > 0 and not len(tk) == self.klen:
                message = self.message % str(line) if "%s" in self.message else self.message
                raise ArgumentTypeError(message)
            if not tk[0] in self.choices:
                message = self.message % ("%s - %s not supported" % (str(line), tk[0])) if "%s" in self.message else self.message
                raise ArgumentTypeError(message)
            items = [
                list(map(lambda x: x.strip(), _tk.split(self.sep2)))
                        if self.itemtype is str else
                islist(_tk, numbertype=self.itemtype, positive=self.positive, sep=self.sep2)
                        for _tk in tk[1:]
            ]
            lengths = list(map(len, items))
            if lengths.count(0) > 0 or self.llen>0 and not lengths.count(self.llen) == len(lengths):
                message = self.message % ('"%s"' % line) if "%s" in self.message else self.message
                raise ArgumentTypeError(message)
            ret[ tk[0]] = items[0] if len(tk) == 2 else items
            # if self.itemtype is str:
            #     items = list(map(lambda x: x.strip(), tk[1].split(self.sep2)))
            # else:
            #     items = islist(tk[1],numbertype=self.itemtype,
            #         positive=self.positive, sep=self.sep2)
            # if len(items) == 0 or self.llen>0 and not len(items) == self.llen:
            #     message = self.message % ('"%s"' % line) if "%s" in self.message else self.message
            #     raise ArgumentTypeError(message)
            # ret[ tk[0]] = items
        return ret

def islist(string, numbertype=int, positive=False, sep=','):
    '''Check if the given string is list of positive numbers (of numbertype) separated by sep.
       If yes returns the list of integers otherwise returns an empty list. If positive is True
       the number should also be positive. '''
    fnc =  ispositive if positive else isnumber
    numbers = [fnc(x, numbertype) for x in string.split( sep)]
    if None in numbers:
        numbers = []
    return numbers

def isrange(string, positive=True, sep=',', rangesep=':'):
    '''Check if the given string is list of positive numbers (of int type) ranges separated by sep.
       If yes returns the list of integers otherwise returns an empty list. If positive is True
       the number should also be positive. The range speperator can be specifide using the rangesep argument. '''
    fnc = ispositive if positive else isnumber
    tk = string.split( sep)
    numbers = []
    for t in tk:
        r = t.split( rangesep)
        r1 = fnc( r[0])
        r2 = r1 if len( r) == 1 else fnc( r[1])
        s = 1 if not len( r) == 3 else int(r[2])
        if r1 and r2:
            numbers += list(range( min(r1,r2), max(r1,r2)+1, s))
        else:
            numbers = []
            break
    return numbers

def enum(*sequential, **named):
    ''' Implements an enum in python (2). '''
    enums = dict(list(zip(sequential, list(range(len(sequential))))), **named)
    reverse = dict((value, key) for key, value in enums.items())
    enums['reverse'] = reverse
    return type('Enum', (), enums)

def call_with_timer(f, a):
    ''' Call function f with arguments a and time the execution. The results are
        printed in the console. '''
    print("Function ", f.__name__)
    t0 = os.times()
    print(list(map( type, a)))
    ret = f(*a)
    t1 = os.times()
    print("\t elapsed time: %f" % (t1[4]-t0[4]))
    print("\t user time: %f" % (t1[0]-t0[0]))
    print("\t system time: %f" % (t1[1]-t0[1]))
    print("\t cpu time: %f" % (t1[0] + t1[1]-t0[0]-t0[1]))
    print("\t cpu time and system call: %f" % (t1[2] + t1[3]-t0[2]-t0[3]))
    return ret

def sort_connected_lists(list1, list2, list3, reverse=False):
    ''' Sort ascending the given connected lists. If *reverce* is True the sort is descending.
        *list1* : type [], indent *in/out*. The sort will be prformed on this list
        *list2* : type [], indent *in/out*. Will be reordered according to the sort on *list1*
        *list3* : type [], indent *in/out*. Will be reordered according to the sort on *list1*
        *reverse* : type boolean '''
    zipped = list(zip(list1, list2, list3))
    zipped.sort()
    t1, t2, t3 = list(zip(*zipped))
    list1[:] = t1
    list2[:] = t2
    list3[:] = t3
    if reverse:
        list1.reverse()
        list2.reverse()
        list3.reverse()

def sort_connected_lists_check_multiplicity(list1, list2, list3, reverse=False):
    ''' Sort ascending the given connected lists. If *reverce* is True the sort is descending.
        The multiplicity of *list1* is handled using *list2* corresponding elements.
        *list1* : type [], indent *in/out*. The sort will be prformed on this list
        *list2* : type [], indent *in/out*. Will be reordered according to the sort on *list1*
        *list3* : type [], indent *in/out*. Will be reordered according to the sort on *list1*
        *reverse* : type boolean '''
    sort_connected_lists(list1, list2, list3, reverse)
    d = defaultdict(lambda: 0)
    for i in list1: d[i] += 1
    n = len(list1)
    i = 0
    while i < n:
        j = list1[i]
        m = d[j]
        if m > 1:
            i_ = i + m
            sl1 = list1[i:i_]
            sl2 = list2[i:i_]
            sl3 = list3[i:i_]
            sort_connected_lists(sl2, sl1, sl3, reverse)
            # or (is the same !)
            #sort_connected_lists(sl1, sl2, sl3, reverse)
            list1[i:i_] = sl1[0:]
            list2[i:i_] = sl2[0:]
            list3[i:i_] = sl3[0:]
            i += m
        else:
            i += 1

def chk_filename( filename):
    ''' Basic check of the given file. Returns the tuple:
        ( kind, DIR, BASENAME, EXT, TRJFILE, TOPFILE)
        kind: the kind of the given file ['lmp','gmx','pdb','unk']
        DIR: directory of the path
        BASENAME: file base name
        EXT: file extention
        TRJFILE: trajectory file name
        TOPFILE: topology file name
    '''


    # get files (clean-up this : partially contained in the readers)
    DIR=os.path.dirname(filename)

    ( DIR, TRJFILE ) = os.path.split( filename)
    ( BASENAME, EXT ) = os.path.splitext( TRJFILE)

    if EXT[1:] in ('trr','xtc'):
        TOPFILE=DIR+os.sep+BASENAME+".tpr"
        kind = 'gmx'
        if len(TOPFILE) == 0:
            TOPFILE=DIR+os.sep+BASENAME+".psf"
    elif EXT[1:] in ('dump'):
        TOPFILE=DIR+os.sep+BASENAME+".data"
        kind = 'lmp'
    elif EXT[1:] in ('pdb'):
        TOPFILE=""
        kind = 'pdb'
    else:
        kind = 'unk'

    # restore trj file using the full path given
    TRJFILE = filename

    return kind, DIR, BASENAME, EXT, TRJFILE, TOPFILE

def parse_radii( finfo):
    ''' Read and close the radii file (finfo). Returns a dict:
            { ["elements"|"types"] : [{[element|type]:[d/r]} }
        where the key specifies if the radii are element
        or type based and the value is a dictionary where
        the element or type (depending on the key) maps to
        the radius value.
    '''

    # open the file and read all the lines at once.
    lines = tuple( map( lambda x: x.strip(), finfo.readlines()))
    finfo.close()

    info = defaultdict( dict)
    # parse the first line:
    # [element|type] [d/r] c
    # the first token is the inderified for the atoms
    # the second indicate is the radius or the diameter is given
    # the third indicate the column of the radius
    tokens = lines[0].split()
    addto = info[tokens[0]]
    f = 1.0 if tokens[1] == 'r' else 0.5
    c = int(tokens[2])
    # now parse the data
    for line in lines[1:]:
        # skip comment and empty lines
        if line.startswith("!") or len(line) == 0: continue
        tokens = line.split()
        if len(tokens) > 1:
            addto[ tokens[0].upper()] = float(tokens[c]) * f
    return info

def read_ndx(f):
    ''' Read gromacs index file and return a dict where the keys are the group names and the
        values are the group index shaped as 1d np.array. '''
    lines = f.readlines()
    _d=defaultdict(str)
    _type = "default"
    for line in lines:
        line=line.strip()
        if len(line) == 0:
            continue
        elif line.startswith("["):
            _type=line[1:-1].strip()
        else:
            _d[_type]+=" "+line
    return { k:np.array(list(map( int, _d[k].split()))) for k in list(_d.keys()) }

def argparse_moleculestype( string):
    ''' check the "-molecules" option arguments. '''
    if len( string) == 0:
        return []
    numbers = isrange( string, positive=True)
    if len( numbers) == 0:
        msg = "wrong molecules indexs range (check: %s)" % string
        import argparse
        raise argparse.ArgumentTypeError(msg)
    return numbers