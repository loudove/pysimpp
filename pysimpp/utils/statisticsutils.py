# -*- coding: utf-8 -*-

import sys
from math import floor, sqrt
import numpy as np
from collections import defaultdict

_halfbin=0.0

def statistics( a, nblocks=1, start=0, end=-1 ):
    ''' Returns the (min, max, mean, std) of the variables
        a[start:end] using n blocks. '''
    if end == -1: end = len(a)
    N = end - start
    if N < 1: return (0.0, 0.0, 0.0, 0.0)
    atmp = np.array( a[start:end], dtype=np.float32)
    if nblocks == 1:
        blk = atmp
    else:
        # blocks should be of equal size
        n = int( N / nblocks)
        n0 =  end - n * nblocks
        blk = np.zeros( nblocks, dtype=np.float32)
        for i in range( nblocks):
            blk[i] = atmp[ n0 + i * n: n0 + (i+1) * n].mean()
    return (  blk.mean(), blk.std(), atmp.min(), atmp.max())

class Variable():
    ''' Implements a varialbe. '''

    def __init__(self):
        ''' Initialize the variable. '''
        self.reset()

    def reset(self):
        ''' Reset the variable. '''
        self._current = 0.0               # current value
        self._min = sys.float_info.min    # miminmum value
        self._max = -sys.float_info.max    # maximum value
        self._n = 0.0                     # number of values
        self._sum = 0.0                   # sum of values
        self._sumsq = 0.0                 # sum of values squares

    def set(self, value, m=1):
        ''' Set the current value with multiplicity m and udpate the info.'''
        self._current = value
        if (value < self._min):
            self._min = value
        elif (value > self._max):
            self._max = value
        self._n += 1.0 * m
        self._sum += value * m
        self._sumsq += value * value * m

    def value(self):
        ''' Get the current value. '''
        return self._current

    def min(self):
        ''' Get the minimum recorded value '''
        return self._min

    def max(self):
        ''' Get the maximum recorded value '''
        return self._max

    def mean(self):
        ''' Get the mean value. '''
        if self._n == 0.0: return 0.0
        return(self._sum / self._n)

    def std(self):
        ''' Get the standard deviation. '''
        if self._n < 2: return -1.0
        s = self._sumsq / (self._n-1.0)
        m = self._sum / self._n
        m_ = self._sum / (self._n-1.0)
        std_ = s - m * m_
        return sqrt(abs(std_))

class HistogramException(Exception): ...

class Binning():
    ''' Implements binning functionality. '''

    FIXED = 1
    FREE = 2

    def __init__(self):
        ''' Initialize the histogram. '''
        self.type = -1                # histogram type 1: fixed, 2:free
        self.n = 0.0                  # number of bin values in the histogram
        self.d = 0                    # bining
        self.min = sys.float_info.min # min value in fixed, reference value in free
        self.max = sys.float_info.max # max value in fixed
        self.range = 0.0              # range = max - min
        self.h = None                 # histogram data
        self.variable = Variable()    # global variable

    @classmethod
    def fixed(cls, min, max, d):
        ''' Initialize a fixed histogram. '''
        obj = cls()
        obj.type = obj.FIXED
        obj.min = min
        obj.max = max
        obj.d = d
        obj.range = max - min
        obj.n = int(obj.range / obj.d) + 1
        obj.h = np.array( [ Variable() for i in range(obj.n + 1)])
        return obj

    @classmethod
    def free(cls, d, ref, addref=True):
        ''' Initialize a free histogram. '''
        obj = cls()
        obj.type = obj.FREE
        obj.h = defaultdict(lambda: Variable())
        i = int(ref / d)
        obj.min = i * d
        obj.d = d
        if addref: obj.add(ref)  # there is no need to add it
        return obj

    def reset(self):
        if self.type == Histogram.FIXED:
            for v in self.h: v.reset()
        elif self.type == Histogram.FREE:
            self.h.clear()
        self.variable.reset()

    def add(self, b, value=1.0, m=1):
        ''' Add a value in the histogram with multiplicity m. '''
        if self.type == -1: return
        self.variable.set(value)
        n = int(floor((b-self.min) / self.d))
        if self.type == Histogram.FIXED:
            if (n < 0 or n > self.n): raise HistogramException("value out of range")
        self.h[n].set( value, m)

    def addat(self, n, value=1.0, m=1):
        ''' Add a value in the histogram at the specific bin with multiplicity m.. '''
        if self.type == -1: return
        self.variable.set(value)
        if self.type == Histogram.FIXED:
            if (n < 0 or n > self.n): raise HistogramException("value out of range")
        self.h[n].set( value, m)

    def add_histogram(self, other, options):
        ''' Add the given histogram . The options dict contains
            the keys defining the addition details. Currently
            works only for free type histograms.
        '''
        if not self.type == Histogram.FREE: return

        _h_other = other.h
        _h_self = self.h
        _add = self.addat

        # check for bins only in self or other histograms and
        # make sure that all bins have the same number of hits
        if len(_h_self):
            _self_set = set(_h_self)
            _other_set = set(_h_other)
            _nhits = next(iter(_h_self.values()))._n
            # for bins only in self histogram add 0.0
            for k in _self_set-_other_set: _add(k,value=0.0)
            # for bins only in other histogram (i.e. new bins) add
            # 0.0 a number of times since all bins should have the
            # same number of hits.
            for k in _other_set-_self_set: _add(k,value=0.0,m=_nhits)

        if options['type'] == 'r': # raw data, just add
            for k in _h_other.keys(): _add( k, value=_h_other[k])
        else:
            _profile = options['profile']
            _d = other.d
            _pi = np.pi
            if _profile == 's': # spherical
                for k in _h_other.keys():
                    r1 = (k*_d)
                    r2 = r1+_d
                    _v = 4.0 / 3.0 * _pi * (r2**3-r1**3)
                    _add( k, value=_h_other[k]/_v)
            elif _profile == 'c': # cylindrical
                _length = options['length']
                for k in _h_other.keys():
                    r1 = (k*_d)
                    r2 = r1+_d
                    _s = _pi * (r2**2-r1**2)
                    _add( k, value=_h_other[k]/(_s*_length))
            elif _profile == 'a': # axial
                _surface = options['surface']
                _v = _surface * _h_other.b
                for k in _h_other.keys(): _add( k, value=_h_other[k]/_v)

    def bin_range(self):
        ''' Get the bin range. '''
        if self.type == Histogram.FIXED:
            return (0, self.n)
        elif self.type == Histogram.FREE:
            keys = list(self.h.keys())
            return (min(keys), max(keys))
        return None

    def data(self):
        ''' Get histogram data. '''
        _h = self.h

        hd = _halfbin * self.d
        x = None
        y = None
        std = None
        if not _h == None:
            if self.type == Histogram.FIXED:
                x = np.zeros(self.n, dtype=np.float32)
                y = np.zeros(self.n, dtype=np.float32)
                std = np.zeros(self.n, dtype=np.float32)
                for i in range(self.n):
                    x[i] = self.min + i * self.d + hd
                    y[i] = _h[i].mean()
                    std[i] = _h[i].std()
            elif self.type == Histogram.FREE:
                (i0, iN) = self.bin_range()
                n = iN - i0 + 1
                x = np.zeros(n, dtype=np.float32)
                y = np.zeros(n, dtype=np.float32)
                std = np.zeros(n, dtype=np.float32)
                for i in range(i0, iN + 1):
                    i_ = i-i0
                    x[i_] = self.min + i * self.d + hd
                    if i in _h:
                        y[i_] = _h[i].mean()
                        std[i_] = _h[i].std()
        return (x, y, std)

    def write( self, filename, header=""):
        ''' Write histogram data to the file filename. '''
        if self.type == Histogram.FREE and len( self.h) == 0:
            return
        f = open( filename, 'w')
        _x, _y, _std = self.data()
        if len(header) > 0:
            f.write("%s\n" % header)
        for _xv, _yv, _stdv in zip( _x, _y, _std):
            f.write(   "%f  %f  %f\n"%(_xv, _yv, _stdv))
        f.close()

class Histogram():
    ''' Implements a histogram. '''

    FIXED = 1
    FREE = 2

    def __init__(self):
        ''' Initialize the histogram. '''
        self.type = -1                # histogram type 1: fixed, 2:free
        self.n = 0.0                  # number of bin values in the histogram
        self.d = 0                    # bining
        self.min = sys.float_info.min # min value in fixed, reference value in free
        self.max = sys.float_info.max # max value in fixed
        self.range = 0.0              # range = max - min
        self.h = None                 # histogram data
        self.variable = Variable()  # histogram variable
        self.ch = None              # cumulative histogram
        self.nch = 0                # number of instances added in the cumulative histogram

    @classmethod
    def fixed(cls, min, max, d):
        ''' Initialize a fixed histogram. '''
        obj = cls()
        obj.type = obj.FIXED
        obj.min = min
        obj.max = max
        obj.d = d
        obj.range = max - min
        obj.n = int(obj.range / obj.d) + 1
        obj.h = np.zeros(obj.n + 1, dtype=np.float32)
        return obj

    @classmethod
    def free(cls, d, ref, addref=True):
        ''' Initialize a free histogram. '''
        obj = cls()
        obj.type = obj.FREE
        obj.h = defaultdict(lambda: 0.0)
        i = int(ref / d)
        obj.min = i * d
        obj.d = d
        if addref: obj.add(ref)  # there is no need to add it
        return obj

    def reset(self):
        if self.type == Histogram.FIXED:
            self.h[:] = 0.0
        elif self.type == Histogram.FREE:
            self.h.clear()
        self.variable.reset()

    def cumulate(self):
        if self.type == Histogram.FIXED:
            if self.ch == None: self.ch = np.zeros(self.h.size, dtype=np.float32)
            self.ch[:] += self.h
        elif self.type == Histogram.FREE:
            if not self.ch: self.ch = defaultdict(lambda: 0.0)
            for k in list(self.h.keys()): self.ch[k] += self.h[k]
        self.nch += 1.0
        self.reset()

    def div_data(self, a):
        ''' Divide self histogram data with a array data. '''
        if type(a) in (int, float):
            if self.type == Histogram.FIXED:
                _h = self.h
                for k in _h.keys():
                    _h[k] /= a
            else:
                self.h /= a
        else:
            if not self.type == Histogram.FIXED or\
                not self.h.size == a.size: return
            for i in range(self.h.size):
                if not a[i] == 0.0: self.h[i] /= a[i]

    def add(self, value, multiplicity=1.0):
        ''' Add a value in the histogram. '''
        if self.type == -1: return
        self.variable.set(value)
        try:
            n = int(floor((value-self.min) / self.d))
        except:
            pass
        if self.type == Histogram.FIXED:
            if (n < 0 or n > self.n): raise HistogramException("value out of range")
        try:
            self.h[n] += multiplicity
        except:
            pass

    def bin_range(self):
        ''' Get the bin range. '''
        if self.type == Histogram.FIXED:
            return (0, self.n)
        elif self.type == Histogram.FREE:
            keys = list(self.h.keys())
            return (min(keys), max(keys))
        return None

    def data(self, cumulative=False):
        ''' Get histogram data. '''
        _h = None
        if cumulative and not self.ch == None:
            _h = self.ch
        else:
            _h = self.h

        hd = _halfbin * self.d
        x = None
        y = None
        if not _h == None:
            if self.type == Histogram.FIXED:
                x = np.zeros(self.n, dtype=np.float32)
                y = np.zeros(self.n, dtype=np.float32)
                for i in range(self.n):
                    x[i] = self.min + i * self.d + hd
                    y[i] = _h[i]
            elif self.type == Histogram.FREE:
                (i0, iN) = self.bin_range()
                n = iN - i0 + 1
                x = np.zeros(n, dtype=np.float32)
                y = np.zeros(n, dtype=np.float32)
                for i in range(i0, iN + 1):
                    i_ = i-i0
                    x[i_] = self.min + i * self.d + hd
                    if i in _h:
                        y[i_] = _h[i]
        return (x, y)

    def normalize(self, cumulative=False):
        ''' Normalize the histogram. '''
        if cumulative:
            if self.ch == None: return
            if self.type == Histogram.FIXED:
                self.ch[:] /= self.nch
            elif self.type == Histogram.FREE:
                for k in list(self.ch.keys()):  self.ch[k] /= self.nch
        else:
            if self.type == Histogram.FIXED:
                self.h[self.n-1] += self.h[self.n]
                self.h[self.n] = 0.0
                s = self.h.sum()
                if s == 0.0: return
                f = 1.0 / s / self.d
                self.h *= f
            elif self.type == Histogram.FREE:
                s = sum(self.h.values())
                if s == 0.0: return
                f = 1.0 / s / self.d
                for k in list(self.h.keys()):  self.h[k] *= f

    def write( self, filename, header="", normalize=True, fmt="%f  %f"):
        ''' Write histogram data to the file filename. If normalize is
            true the data will be normalized first. '''
        if self.type == Histogram.FREE and len( self.h) == 0:
            return
        f = open( filename, 'w')
        if normalize:
            self.normalize()
        _x, _y = self.data()
        if len(header) > 0:
            f.write("%s\n" % header)
        _fmt=fmt+"\n"
        for _xv, _yv in zip( _x, _y):
            f.write( _fmt % (_xv, _yv))
        f.close()

    def getf(self):
        ''' Return a dict with the real data. '''
        _x, _y = self.data()
        h = {}
        for k, v in zip( _x, _y):
            h[k]=v
        return h

class Histogram2D():
    ''' Implements a 2d Histogram. Data are kept in the
        structure:
            h{i0:h1{i1:v}}
    '''

    def __init__(self, d, ref, addref=True):
        ''' Initialize the histogram. '''
        self.n = np.array((0, 0))        # number of bin values
        self.d = np.array((0.0, 0.0))    # bining
        self.min = np.array((0.0, 0.0))  # reference values
        self.range = np.array((0.0, 0.0))# range = max - min
        self.h = defaultdict(lambda: defaultdict(float)) # histogram data
        self.variable = [Variable(), Variable()]        # histogram variable

        i = (np.array(ref) / d).astype(int)
        self.min[:] = i * d
        self.d[:] = d
        if addref: self.add(ref)

    def add(self, value, factor=1.0):
        ''' Add a value in the histogram. '''
        if not  len(value) == 2: return
        self.variable[0].set(value[0])
        self.variable[1].set(value[1])

        n = np.floor((value-self.min) / self.d).astype(int)
        self.h[n[0]] [n[1]] += factor

    def bin_range(self):
        ''' Get the bin range. '''
        range = np.zeros((2, 2), dtype=np.int32)
        keys = list(self.h.keys())
        range[0, 0] = min(keys)
        range[0, 1] = max(keys)
        min_ = sys.maxsize
        max_ = -sys.maxsize
        for v in list(self.h.values()):
            keys = list(v.keys())
            tmp = max(keys)
            if tmp > max_: max_ = tmp
            tmp = min(keys)
            if tmp < min_: min_ = tmp
        range[1, 0] = min_
        range[1, 1] = max_
        return range

    def write(self, fname, header="", include=True):
        ''' Write the histogram data. NOTE that the histogram should be normalized if
            the print out of the estimated probability distribution function is desired. '''
        if len(list(self.h.keys())) == 0: return
        hd = _halfbin * self.d

        i = self.bin_range()
        #n0=i[0,1]-i[0,0]+1
        #n1=i[1,1]-i[1,0]+1

        f = open(fname, 'w')
        f.write("%s\n" % header)
        if not include:
            for i0 in range(i[0, 0], i[0, 1] + 1):
                if i0 not in self.h: continue
                hash = self.h[i0]
                hi0 = self.min[0] + self.d[0] * i0 + hd[0]
                for i1 in range(i[1, 0], i[1, 1] + 1):
                    if i1 not in hash: continue
                    v = hash[i1]
                    hi1 = self.min[1] + self.d[1] * i1 + hd[1]
                    f.write( "%s %s %s\n" % ( str(hi0).ljust(15), str(hi1).ljust(15), str(v).ljust(15)) )
                f.write("\n")
        else:
            for i0 in range(i[0, 0], i[0, 1] + 1):
                hi0 = self.min[0] + self.d[0] * i0 + hd[0]
                if i0 not in self.h:
                    for i1 in range(i[1, 0], i[1, 1] + 1):
                        hi1 = self.min[1] + self.d[1] * i1 + hd[1]
                        f.write( "%s %s 0.0\n" % ( str(hi0).ljust(15), str(hi1).ljust(15)) )
                    f.write("\n")
                else:
                    hash = self.h[i0]
                    for i1 in range(i[1, 0], i[1, 1] + 1):
                        if i1 in hash:
                            v = hash[i1]
                        else:
                            v = 0.0
                        hi1 = self.min[1] + self.d[1] * i1 + hd[1]
                        f.write( "%s %s %s\n" % ( str(hi0).ljust(15), str(hi1).ljust(15), str(v).ljust(15)) )
                    f.write("\n")
        f.close()

    def write_conditional(self, fname, header="", include=True, condv=0, condp=None):
        ''' Write the histogram data. NOTE that the histogram should be normalized
            prior to the call of this method. '''
        if len(list(self.h.keys())) == 0: return
        hd = _halfbin * self.d

        i = self.bin_range()
        #n0=i[0,1]-i[0,0]+1
        #n1=i[1,1]-i[1,0]+1

        f = open(fname, 'w')
        f.write("%s\n" % header)
        if not include:
            for i0 in range(i[0, 0], i[0, 1] + 1):
                if i0 not in self.h: continue
                hash = self.h[i0]
                hi0 = self.min[0] + self.d[0] * i0 - hd[0]
                for i1 in range(i[1, 0], i[1, 1] + 1):
                    if i1 not in hash: continue

                    v = hash[i1]
                    if not condv == 0:
                        if condv==1: v/=condp[i0]
                        elif condv==2: v/=condp[i1]

                    hi1 = self.min[1] + self.d[1] * i1 - hd[1]
                    f.write( "%s %s %s\n" % ( str(hi0).ljust(15), str(hi1).ljust(15), str(v).ljust(15)) )
                f.write("\n")
        else:
            for i0 in range(i[0, 0], i[0, 1] + 1):
                hi0 = np.float32(self.min[0] + self.d[0] * i0 - hd[0])
                if i0 not in self.h:
                    for i1 in range(i[1, 0], i[1, 1] + 1):
                        hi1 = np.float32(self.min[1] + self.d[1] * i1 - hd[1])
                        f.write( "%s %s 0.0\n" % ( str(hi0).ljust(15), str(hi1).ljust(15)) )
                    f.write("\n")
                else:
                    hash = self.h[i0]
                    for i1 in range(i[1, 0], i[1, 1] + 1):
                        hi1 = np.float32( self.min[1] + self.d[1] * i1 - hd[1])
                        if i1 in hash:
                            v = hash[i1]
                            if not condv == 0:
                                if condv==1:
                                    if hi0 not in condp:
                                        raise HistogramException("error: conditional probability is not defined for x=%f" % hi0)
                                    _v = condp[hi0]
                                    if not _v == 0: v/=_v
                                    else: raise HistogramException("error: conditional probability is zero for x=%f" % hi0)
                                elif condv==2:
                                    if hi1 not in condp:
                                        raise HistogramException("error: conditional probability is not defined for x=%f" % hi1)
                                    _v = condp[hi1]
                                    if not _v == 0: v/=_v
                                    else: raise HistogramException("error: conditional probability is zero for x=%f" % hi1)
                        else:
                            v = 0.0
                        # line = str(hi0).ljust(15) + str(hi1).ljust(15) + str(v).ljust(15)
                        # f.write(line + "\n")
                        f.write( "%s %s %s\n" % ( str(hi0).ljust(15), str(hi1).ljust(15), str(v).ljust(15)) )
                    f.write("\n")
        f.close()

    def write_conditional_(self, fname, header="", include=True):
        ''' Write the histogram data with the conditional probabilities
            given the first dimension. '''

        if len(list(self.h.keys())) == 0: return
        hd = _halfbin * self.d
        # the probability of the events in the first dimension
        px = {}
        for x, vx in self.h.items():
            px[x] = sum( vx.values())

        i = self.bin_range()
        #n0=i[0,1]-i[0,0]+1
        #n1=i[1,1]-i[1,0]+1

        f = open(fname, 'w')
        f.write("%s\n" % header)
        if not include:
            for i0 in range(i[0, 0], i[0, 1] + 1):
                if i0 not in self.h: continue
                hash = self.h[i0]
                hi0 = self.min[0] + self.d[0] * i0 - hd[0]
                px0 = px[i0]
                for i1 in range(i[1, 0], i[1, 1] + 1):
                    if i1 not in hash: continue

                    v = hash[i1] / px0

                    hi1 = self.min[1] + self.d[1] * i1 - hd[1]
                    f.write( "%s %s %s\n" % ( str(hi0).ljust(15), str(hi1).ljust(15), str(v).ljust(15)) )
                f.write("\n")
        else:
            for i0 in range(i[0, 0], i[0, 1] + 1):
                hi0 = np.float32(self.min[0] + self.d[0] * i0 - hd[0])
                if i0 not in self.h:
                    for i1 in range(i[1, 0], i[1, 1] + 1):
                        hi1 = np.float32(self.min[1] + self.d[1] * i1 - hd[1])
                        f.write( "%s %s 0.0\n" % (str(hi0).ljust(15), str(hi1).ljust(15) ))
                    f.write("\n")
                else:
                    hash = self.h[i0]
                    px0 = px[i0]
                    for i1 in range(i[1, 0], i[1, 1] + 1):
                        hi1 = np.float32( self.min[1] + self.d[1] * i1 - hd[1])
                        if i1 in hash:
                            v = hash[i1] / px0
                        else:
                            v = 0.0
                        f.write( "%s %s %s\n" % ( str(hi0).ljust(15), str(hi1).ljust(15), str(v).ljust(15)) )
                    f.write("\n")
        f.close()


    def normalize(self):
        ''' Normalize the histogram. '''
        s = 0.0
        for hash in list(self.h.values()):
            s  += sum(hash.values())
        if s == 0.0: return
        f = 1.0 / s / (self.d[0] * self.d[1])
        for hash in list(self.h.values()):
            for key in hash:
                hash[key] *= f

class Histogram3D():
    ''' Implements a 3d Histogram. Data are kept in the
        structure:
            h{i0:h1{i1:h2{i2:v}}}
    '''

    def __init__(self, d, ref, addref=True):
        ''' Initialize the histogram. '''
        self.n = np.array((0, 0, 0))           # number of bin values
        self.d = np.array((0.0, 0.0, 0.0))     # bining
        self.min = np.array((0.0, 0.0, 0.0))   # reference values
        self.range = np.array((0.0, 0.0, 0.0)) # range = max - min
        self.h = defaultdict(lambda: defaultdict(lambda: defaultdict(float))) # histogram data
        self.variable = [Variable(), Variable(), Variable()]                  # histogram variable

        i = (np.array(ref) / np.array(d)).astype(int)
        self.min[:] = i * d
        self.d[:] = d
        if addref: self.add(ref)

    def add(self, value, factor=1.0):
        ''' Add a value in the histogram. '''
        if not  len(value) == 3: return
        self.variable[0].set(value[0])
        self.variable[1].set(value[1])
        self.variable[2].set(value[2])

        n = np.floor((value-self.min) / self.d).astype(int)
        self.h[n[0]] [n[1]] [n[2]] += factor

    def bin_range(self):
        ''' Get the bin range. '''
        range = np.zeros((3, 2), dtype=np.int32)
        keys = list(self.h.keys())
        range[0, 0] = min(keys)
        range[0, 1] = max(keys)

        min1 = min2 = sys.maxsize
        max1 = max2 = -sys.maxsize
        for v1 in list(self.h.values()):
            keys = list(v1.keys())
            max1 = max( max1, max(keys))
            min1 = min( min1, min(keys))
            for v2 in list(v1.values()):
                keys2 = list(v2.keys())
                max2 = max( max2, max(keys2))
                min2 = min( min2, min(keys2))
        range[1, 0] = min1
        range[1, 1] = max1
        range[2, 0] = min2
        range[2, 1] = max2
        return range

    def write(self, fname, format='', include=True):
        if format == 'cube':
            self._write_cube(fname, include)
        else:
            self._write(fname, include)

    def  _write_cube(self, fname, include=True):
        ''' Write histogram data in cube format. '''
        hd = 0.5 * self.d
        i = self.bin_range()
        f = open(fname, 'w')
        n = i[:,1]-i[:,0]+1

        f.write("HISTOGRAM3D CUBE FILE\n")
        f.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
        f.write("%5d %12.6f %12.6f %12.6f\n"%(1,0.,0.,0.))
        f.write("%5d %12.6f %12.6f %12.6f\n"%(n[0],self.d[0],0.,0.))
        f.write("%5d %12.6f %12.6f %12.6f\n"%(n[1],0.,self.d[1],0.))
        f.write("%5d %12.6f %12.6f %12.6f\n"%(n[2],0.,0.,self.d[2]))
        f.write("%5d %12.6f %12.6f %12.6f %12.6f\n"%(1,1.,0.,0.,0.))

        STRZERO = "%g "%0.0
        nn = 0
        for i0 in range(i[0, 0], i[0, 1] + 1):
            if i0 not in self.h:
                for i1 in range(i[1, 0], i[1, 1] + 1):
                    for i2 in range(i[2, 0], i[2, 1] + 1):
                        f.write(STRZERO)
                        nn += 1
                        if nn % 6 == 5: f.write("\n")
            else:
                h1 = self.h[i0]
                for i1 in range(i[1, 0], i[1, 1] + 1):
                    if i1 not in h1:
                        for i2 in range(i[2, 0], i[2, 1] + 1):
                            f.write(STRZERO)
                            nn += 1
                            if nn % 6 == 5: f.write("\n")
                    else:
                        h2 = h1[i1]
                        for i2 in range(i[2, 0], i[2, 1] + 1):
                            v = "%g "%h2[i2] if i2 in h2 else STRZERO
                            f.write(v)
                            nn += 1
                            if nn % 6 == 5: f.write("\n")
        f.close()

    def _write(self, fname, include=True):
        ''' Write histogram data in grid format. '''

        hd = 0.5 * self.d
        i = self.bin_range()
        f = open(fname, 'w')
        if not include:
            for i0 in range(i[0, 0], i[0, 1] + 1):
                if i0 not in self.h: continue
                hi0 = self.min[0] + self.d[0] * i0 + hd[0]
                h1 = self.h[i0]
                for i1 in range(i[1, 0], i[1, 1] + 1):
                    if i1 not in h1: continue
                    hi1 = self.min[1] + self.d[1] * i1 + hd[1]
                    h2 = h1[i1]
                    for i2 in range(i[2,0],i[2,1]+1):
                        if i2 not in h2: continue
                        hi2 = self.min[2] + self.d[2] * i2 + hd[2]
                        v = h2[i2]
                        f.write( "%s %s %s %s\n" % ( str(hi0).ljust(15), str(hi1).ljust(15), str(hi2).ljust(15) , str(v).ljust(15)) )
                    # f.write("\n")
                # f.write("\n")
        else:
            for i0 in range(i[0, 0], i[0, 1] + 1):
                hi0 = self.min[0] + self.d[0] * i0 + hd[0]
                if i0 not in self.h:
                    for i1 in range(i[1, 0], i[1, 1] + 1):
                        hi1 = self.min[1] + self.d[1] * i1 + hd[1]
                        for i2 in range(i[2, 0], i[2, 1] + 1):
                            hi2 = self.min[2] + self.d[2] * i1 + hd[2]
                            f.write( "%s %s %s 0.0\n" % ( str(hi0).ljust(15), str(hi1).ljust(15), str(hi2).ljust(15)) )
                        # f.write("\n")
                    # f.write("\n")
                else:
                    h1 = self.h[i0]
                    for i1 in range(i[1, 0], i[1, 1] + 1):
                        hi1 = self.min[1] + self.d[1] * i1 + hd[1]
                        if i1 not in h1:
                            for i2 in range(i[2, 0], i[2, 1] + 1):
                                hi2 = self.min[2] + self.d[2] * i1 + hd[2]
                                f.write( "%s %s %s 0.0\n" % ( str(hi0).ljust(15), str(hi1).ljust(15), str(hi2).ljust(15)) )
                            # f.write(line + "\n")
                        else:
                            h2 = h1[i1]
                            for i2 in range(i[2, 0], i[2, 1] + 1):
                                hi2 = self.min[2] + self.d[2] * i1 + hd[2]
                                v = h2[i2] if i2 in h2 else 0.0
                                f.write( "%s %s %s %s\n" % ( str(hi0).ljust(15), str(hi1).ljust(15), str(hi2).ljust(15) , str(v).ljust(15)) )
                            # f.write(line + "\n")
                    # f.write("\n")
        f.close()

    def _get_px(self):
        ''' Return the probability distribution function of the first dimension
            There is no need for the current histogram to be normalized. '''
        px = defaultdict(float)
        stot = 0.0
        for x, vx in self.h.items():
            s = 0.0
            for vy in vx.values():
                s += sum(vy.values())
            px[x]=s
            stot += s
        if not np.isclose(stot, 0.0):
            f = 1.0 / stot / self.d[0]
            for k in px.keys():
                px[k] *= f
        return px

    def get_h2d_x(self, x):
        ''' Returns the 2d histogram for a given value x of the first
            dimension. The current histogram should have been normalized.
            This result should be the conditional probabiliry:
              P[(y,z)|x]=P[(x,y,z)]/P[x]
            with P[(x,y,z)] been estimated by the current normalized
            histogram while P(x) is the estimated by the integral of
            P[(x,y,z)] with respect y and z variables:
              P[x] = Integral { dy dz P[(x,y,z)] } '''
        d = self.d[0]
        obj = Histogram2D(d, (0.0,0.0), addref=False)
        obj.min[:] = self.min[1:3]

        h2d = defaultdict(lambda: defaultdict(float))
        if x in self.h:
            _px = self._get_px()
            px = _px[x]
            obj.variable[:] = self.variable[1:3]
            vx = self.h[x]
            for y, vy in vx.items():
                for z, vz in vy.items():
                    h2d[y][z] = vz/px
        obj.h = h2d
        return obj

    def normalize(self):
        ''' Normalize the histogram. '''
        s = 0.0
        for hash1 in list(self.h.values()):
            for hash2 in list(hash1.values()):
                s  += sum(hash2.values())
        if s == 0.0: return
        f = 1.0 / s / (self.d[0] * self.d[1] * self.d[2])
        for hash1 in list(self.h.values()):
            for hash2 in list(hash1.values()):
                for key in hash2:
                    hash2[key] *= f
