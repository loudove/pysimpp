# -*- coding: utf-8 -*-

import os
import sys
from collections import defaultdict
from operator import add
import pickle

import numpy as np
from scipy import stats

import pysimpp.readers
from pysimpp.utils.utils import isrange, islist, read_ndx
from pysimpp.utils.statisticsutils import Histogram

from pysimpp.fastpost import fasts1, fasts1x, fastwrap, fastunwrapv, fastunwrap, fastcom, fastcom_total

def _is_command(): return True
def _short_description(): return 'Calculate mean square displacement.'
def _command(): command()

__debug = False
__rmcm = True

class MSDException(Exception):
    ''' Implements MSDException custom exception. '''

    def __init__(self, value):
        self.value = "MSDException exception: " + value

    def __str__(self):
        return repr(self.value)

class MSDDataUtility():
    ''' Implements MSD data fit and write utility class. '''

    def __init__(self, t, msd, msds, dimensions):
        ''' Initialize a MSDDataUtility object.'''
        # TODO check data validity
        # keep the data
        self.t = t
        self.msds =[]               # one msd per dimension
        dmap = {'x':0,'y':1,'z':2}
        for d in dimensions:
            self.msds.append( msds[ dmap[ d]])
        self.msd = sum( msds)    # total msd
        self.nconfs = len( self.msd)
        self.ndim = len( self.msds)
        self.dimensions = dimensions

        # work to the first set of non zero data (zero step excluded)
        # suppose that if msd is non zero all msds are also non zero
        self.upto = self.nconfs
        spotzero = np.where( msd[1:] == 0.0)[0]
        if not spotzero.size == 0:
            self.upto = spotzero[0]+1

        # print smds and their fits
        self.__printdims = True

    def __fit_gmx(self, weight=True):
        ''' atoms_msds[i,j,k] contains
        '''
        pass

    def fit(self, var=0.1):
        '''Fit the data to calculate self diffusion coefficient. Returns a list
           of lists with the fit results for the total msd (fist element) and for
           each of the dimensions (rest of the components). Each tuple contains
           the fit data.
             (slope, intercept, r_value, p_value, std_err, start, end)
        '''
        ret = []
        if self.ndim > 1:
            ret = [ self.__fit( self.msds[i]) for i in range( self.ndim)]
        ret.insert(0, self.__fit( msd))

    def __fit(self, msd, var=0.1):
        '''Fit the data to calculate self diffusion coefficient. Returns a
           tuple with fit results:
             (slope, intercept, r_value, p_value, std_err, start, end)
        '''

        upto = self.upto

        # log10(t) - exclude first time where msd is zero
        logt = np.log10( self.t[1:upto]-self.t[0]) # time zero is zero!
        # log10(msd)
        logmsd = np.log10( msd[ 1:upto])
        # gradient of log10(msd) with log10(t)
        grad = np.gradient( logmsd) / np.gradient( logt)

        # trace where the 1-var<grad<1+var
        _pone = 1.0 + var
        _mone = 1.0 - var
        _where = np.where( (grad<_mone)|(grad>_pone))
        ranges = np.arange( grad.size)+1
        ranges[ _where] = -1
        # remove the very first part
        try:
            _i0 = np.where( grad >=1.0 )[0][0]
            if _i0 > 0:
                ranges[:(_i0-1)] = -1
        except:
            pass

        ok = False
        istart = iend = 0
        slope = intercept = r_value = p_value = std_err = 0.0
        try:
            istart, iend = _max_consecutive( ranges)
            if iend - istart > 2: # more than two points are needed
                istart, iend = istart+1, iend+1  # increas by one to get the correct indexes
                slope, intercept, r_value, p_value, std_err = stats.linregress(self.t[istart:iend],msd[istart:iend])
                ok = True         # set success flag
        except:
            pass
        finally:
            if not ok:
                istart, iend = 0, 0

        return [slope, intercept, r_value, p_value, std_err, istart, iend, grad]

    def write( self, f=None, every=1, dofit=True, var=0.1):
        '''Writes the data to the give file unit. If dofit is true
           the fit is performed and the results are also written for
           the total msd and the msd for each dimension. Finally the
           volume of data can be reduced by increasing the value of
           argument every. '''

        # if no file unit is provided then write to the default
        closef = False
        if not f:
            f = open('msd.dat','w')
            closef = True

        # utility function for writing
        def __tostr( x):
            return " %s" % str(x)
        def __appendcolumn( lines, head, data):
            lines[0] += head
            for i in range(1,len(lines)): lines[i] += " %s" % str( data[i-1])
        def __printfit( slope, intercept, r_value, p_value, std_err, istart, var, dim, what):
            _m5u = 'x10\u207b\u2075 cm\u00b2/s'
            # _m10u = u'10\u207b\u00b9\u2070 m\u00b2/s'
            return  '# %s: points %d:%d (slope variance %f)\n' % ( what, istart, iend , var) + \
                    '# y = a x + b : a = %f, b = %f \n' % ( slope, intercept) + \
                    '#  r-value = %f , p-value = %f stderr = %f\n' % (r_value, p_value, std_err) + \
                    '# D = %f %s \n' % (slope/float(2*dim)*1.e4, _m5u)

        # do the print out
        lines = [""] * (self.upto+1)
        self.t[:] -= self.t[0]
        __appendcolumn( lines, "# time", self.t[:self.upto])
        __appendcolumn( lines, " msd", self.msd[:self.upto])
        if self.__printdims and self.ndim > 1:
            for i in range( self.ndim):
                __appendcolumn( lines, " msd[%d]"%(i+1), self.msds[i][:self.upto])

        fitlines = []
        if dofit:
            est = np.zeros( self.upto, dtype=np.float32)
            slope, intercept, r_value, p_value, std_err, istart, iend, grd = self.__fit( self.msd, var=var)
            fitlines.append( __printfit( slope, intercept, r_value, p_value, std_err, istart, var, self.ndim, "msd"))
            est[istart:iend] = slope*self.t[istart:iend]+intercept
            __appendcolumn( lines, " est", est)
            grd = np.insert( grd, 0, 0)
            __appendcolumn( lines, " grd", grd)

            if self.__printdims and self.ndim > 1:
                for i in range( self.ndim):
                    est[:] = 0.0
                    slope, intercept, r_value, p_value, std_err, istart, iend, grd = self.__fit( self.msds[i], var=var)
                    fitlines.append( __printfit( slope, intercept, r_value, p_value, std_err, istart, var, 1, "msd[%d]"%(i+1)))
                    est[istart:iend] = slope*self.t[istart:iend]+intercept
                    __appendcolumn( lines, " est[%d]"%(i+1), est)
                    grd = np.insert( grd, 0, 0)
                    __appendcolumn( lines, " grd[%d]"%(i+1), grd)

        f.write("# dimensions: %s\n" % self.dimensions)
        for line in fitlines: f.write( line)
        for line in lines: f.write( line+"\n")
        if closef: f.close()

def _dimsbookkeep( dimensions, unwrap=True):
    ''' Given the dimensions i.e. a list of strings creates the attributes
        string for the lammpsreader and the dimensions map for data storage.'''
    # set what to read and store bookkeeping
    attributes = 'id type'
    xstring = ""
    istring = ""
    dmap = {}
    dmapinv = {}
    ndims = 0
    for dim in [ 'x', 'y', 'z']:
        if dim in dimensions:
            xstring += " %s" % dim
            idim = "i%s" % dim
            istring += " %s" % idim
            dmap[ ndims] = (dim, idim)
            dmapinv[ dim] = ndims
            ndims += 1
    attributes = attributes + xstring
    if unwrap:
        attributes += istring
    return attributes, dmap, dmapinv

# @profile
def msd( filename, start=-1, end=sys.maxsize, every=1, dimensions=['xyz'], region=(), molids=(), maxconfs=-1, ndxfile=None):

    # check the input file
    reader = pysimpp.readers.create(filename)
    if not reader: sys.exit(0)

    usecom = True

    # TODO check for gromacs: the provided trajectory is assumed allready unwrapped.
    _unwrap = True 
    reader.set_unwrap( _unwrap)
    attributes = 'id x y z type'
    reader.set_attributes( attributes)

    dirname = reader.dir
    natoms = reader.natoms

    # get molecular data and select molecules of interest
    # TODO refactor massesar
    types = reader.get_atom_type()              # atom type array (number or string based)
    ntypes = len( set( types))
    masses = reader.get_type_mass()             # type mass
    massesar = np.array( [ masses[ types[iat]] for iat in range( natoms) ]) # atom mass array (accelerate)
    molecules = reader.get_atom_molecule() - 1  # atom molecule array (index @zero)
    nmolecules = molecules.max() + 1            # number of molecules
    # TODO add this functionality in reader
    mol_atoms = defaultdict( list)              # molecule atoms array
    for i, im in enumerate( molecules):
        mol_atoms[ im].append(i)

    # selected molecules
    if not ndxfile is None:
        _grps = read_ndx( ndxfile)
        _molids = set().union( molids, *(_grps.values()))
    else:
        _molids = molids
    
    selected = np.sort( _molids) - 1             # selected molecules (index @zero)
    if not usecom:                              # if needed convert to atoms
        satoms = []
        for im in selected:
            satoms += mol_atoms[im]
        selected = np.sort( satoms)

    nselected = selected.size                   # number of selected molecules
    hasselected = nselected > 0                 # selected flag

    # region based calculation
    hasregion = len( region) > 0

    # bookkeeping dimensions to read
    _dimensions = 'xyz'
    attributes, dmap, dmapinv = _dimsbookkeep( _dimensions, _unwrap)
    ndims = len( _dimensions)

    print('>> count configurations in dump file(s) ...')
    reader.set_attributes( attributes)
    nconfs = reader.count_frames( start, end) if  maxconfs == -1 else maxconfs

    # allocate wrapped and unwrapped molecules coordinates
    n_ = nselected if hasselected else nmolecules if usecom else natoms
    cmw = np.zeros( (nconfs, n_, ndims), dtype=np.float32) # wrapped
    cm = np.zeros( (nconfs, n_, ndims), dtype=np.float32)

    r = np.empty( shape=(natoms,ndims), dtype=np.float32)      # unwrapped coordinates
    if _unwrap:
        rw = np.empty( shape=(natoms,ndims), dtype=np.float32) # wrapped coordinates
        ip = np.empty( shape=(natoms,ndims), dtype=np.int32)   # periodic offsets
    steps = []
    boxes = []

    # initial system center of mass
    cm0 = np.zeros((3),dtype=np.float32)
    # vector for com move removal
    delta = np.zeros((3),dtype=np.float32)

    print('\n>> reading dump file(s) ...')
    iconf = -1
    iconf_ = 0
    while( True):
        step, box, data = reader.read_next_frame()
        iconf_ += 1
        if step == None:
            break
        elif not iconf_%every == 0:
            continue
        elif step < start:
            continue
        elif step > end:
            print(step)
            break
        if box:
            iconf += 1
            steps.append( step)
            boxes.append( box)
            for k, v in dmap.items():
                np.copyto(r[:, k] ,  data[v[0]])

            if __rmcm:
                cmtotal, masstotal = fastcom_total( r, massesar)
                if len(boxes) == 1:
                    cm0[:] = cmtotal                
                else:
                    delta[:] = cmtotal-cm0
                    r -= delta

            # calculate center of masses and if regions
            # set also the wraped coordinates to ensure
            # the correct bining
            if hasselected:
                cm_ = fastcom( r, massesar, molecules, nmolecules) if usecom else r
                cm[iconf,:,:] = cm_[selected,:]
                if hasregion:
                    # LDP TODO check the API of fastwarp
                    cmw_ = fastwrap( cm_.transpose(), box.origin, box.va, box.vb, box.vc)
                    cmw[iconf,:,:] = cmw_[:,selected].transpose()
            else:
                cm[iconf,:,:] = fastcom( r, massesar, molecules, nmolecules) if usecom else r
                if hasregion:
                    cmw_ = fastwrap( cm[iconf,:,:].transpose(), box.origin, box.va, box.vb, box.vc)
                    cmw[iconf,:,:] = cmw_.transpose()
        else:
            break

    # number of configurations - resize if needed
    nconfs = len( steps)
    _s = cm.shape
    if not nconfs == _s[0]:
        print("to save some memory consider to use maxconfs = %d instead of %d" %( nconfs, _s[0]))
        cm = np.resize( cm, (nconfs, _s[1], _s[2]))
        cmw = np.resize( cmw, (nconfs, _s[1], _s[2]))

    # get the timestep in fs
    dt = reader.timestep * 1000.0 # TODO check for gromacs
    dt *= every

    print('\n>> calculate msds(s) ...')
    w = []
    if hasregion == 0:
        print("no region")
        msd, msds = _msd_fft3( cm)
        # msd, msds = _msd_fft3x( r, path=reader.dir+os.sep, w=w)

        print("msd calculated")
        t = np.array( steps, dtype= np.float32) * dt
        du = MSDDataUtility( t, msd, msds, dimensions)
        f = open( reader.dir+os.sep+"msd.dat", 'w')
        du.write(f=f, var=0.10)
        f.close()
    else:
        REGIONDIM = dmapinv[ region[0]]
        LEVELS = region[1]
        NLEVELS = len( LEVELS)

        fr=open(reader.dir+os.sep+"reactions.dat",'w')

        # atom region flag (-1:not interested, 0:in region, 1...n (>0) subregion id)
        REGION = np.zeros( nconfs, dtype=np.int32) # region flag per atom
        REGIONDATA = {}
        for i in range( NLEVELS-1):
            # place it in a class later on:
            #                 msd, [msdx,msdy,msdz], normalize, timeshistogram]
            REGIONDATA[i] = [ np.zeros( nconfs, dtype=np.float32), \
                              np.zeros( (ndims, nconfs), dtype=np.float32), \
                              np.zeros( nconfs, dtype=np.int32), \
                              [], 0 ]

        f = np.zeros( nconfs, dtype=np.float32)
        fs = np.zeros( (ndims,nconfs), dtype=np.float32)
        zmax = 0
        threshold = 2
        for i in range( cm.shape[1]):
            x = cmw[:,i,:]                                          # x[nconfs:ndims] atom coordinataes for configurations
            zmax = max( zmax, x[:, REGIONDIM].max())
            REGION[:] = -1                                          # place in region
            ROLD = -1
            NR = 0
            for j in range( nconfs):
                for k in range( 1, NLEVELS):
                    if LEVELS[k-1] <= x[j, REGIONDIM] <= LEVELS[k]: # check REGIONDIM dimension
                        RNEW = k-1
                        REGION[j] = RNEW
                        if not RNEW == ROLD:
                            fr.write("%d#%d#%d,"%(j,ROLD,RNEW))
                            ROLD = RNEW
                            NR += 1
                            if NR%10==0:fr.write("\n")
                        break

            sequencies = _region_sequences( REGION, threshold=threshold)
            for seq in sequencies:
                if seq[2] == -1: continue # not in regions of interest
                length = seq[1]-seq[0]
                if length <= threshold: continue
                #f, fs = msd_fft( what[seq[0]:seq[1],i,:], atomic=True)
                f, fs = _msd_fft2a( cm[seq[0]:seq[1],i,:])
                if f.size != length:
                    print("ERROR")
                    sys.exit()
                _data = REGIONDATA[ seq[2]]
                _data[0][:length] += f
                _data[1][:,:length] += fs
                _data[2][:length] += 1
                _data[3].append( length) # cumulate residence time for this region
                _data[4] += 1

        for k, v in REGIONDATA.items():
            print("hits in %d : %f "% (k, v[4]))
            if ( v[4] == 0): continue
            msd = v[0]
            msds = v[1]
            norm = v[2]
            nz = np.where(norm!=0)[0]
            msd[ nz] /= norm[ nz]
            msds[ :, nz] /= norm[ nz]

            delta = dt*(steps[1]-steps[0])
            hist = Histogram.free( 1, v[3][0])
            for val in v[3][1:]: hist.add( val)
            f=open( reader.dir+os.sep+"thist%d.dat"%k, 'w')
            hist.normalize()
            _x, _y = hist.data()
            _x*=delta
            _y/=delta
            f.write("# mean residense time (timestep units): %f\n" % (np.sum(_x*_y)*delta))
            for _xv, _yv in zip( _x, _y):
                f.write(   "%f  %g\n"%(_xv, _yv))
            f.close()

            t = np.array( steps, dtype= np.float32) * dt
            du = MSDDataUtility( t, msd, msds, dimensions)
            f = open( reader.dir+os.sep+"msd%d.dat"%k, 'w')
            du.write(f=f, var=0.05)
            f.close()

# @profile
def _region_sequences( regions, threshold=1):
    ''' Given a time sequence of regions (in general id's), returns ranges of
        sub-sequences of the same region in the form of a list of tuples:
        [ (i0,i1,rid)]
        where i0 is the starting indexe, i1-1 is the end index and rid
        it the region of the sequence. A simple analysis is preformed so the
        sub-sequences of length smaller than a given therhold will be "absorbed"
        by the neighboring region. If the frequency of change between two regions
        is to high a warning will be given. If the change between regions is
        to drastic an error will be issued.'''

    # determine the sub-sequences (seq)
    # if we have n seq then we will need:
    # * n size array to store the starting position plus a last
    #   element to facilitate length calculation)
    # * n size array to store the corresponding region id
    _d = np.diff( regions)                  # regions[i+1]-regions[i] i=1,...,
    _w = np.where( _d != 0)[0]+1            # indexes of _d where seq breaks (_a1 non zero)
    s = np.insert( _w, 0, 0)                # starting inexes of each seq
                                            # prepend zero ( facilitate region id and length calculation)
    sid = regions[ s]                       # get the region id for each seq
    s  = np.insert(s, s.size, regions.size) # append regions size ( facilitate length calculation)
    sl = np.diff( s)                        # seq lengths
    ns = sl.size                            # number of seq
    nsm = ns - 1
    schk = np.where( sl <= threshold)[0]    # seq to be merged
    do = np.zeros(ns, dtype=np.int32)       # merge status:
                                            #   <0 merge with next
                                            #   =0 do nothing
                                            #   >0 merge with previous
    # just ignore them
    if True:
        for i in schk:
            if not do[i] == 0: continue
            if not sid[i] == -1:
                sid[i] = -1
                if i == 0:      # check only next
                    if sid[1] == -1:
                        sl[0] = sl[1] = sl[0] + sl[1]
                        do[0] = -1
                        do[1] = 1
                elif i == nsm:
                    im = i - 1  # check only previous
                    if sid[im] == -1:
                        sl[im] = sl[i] = sl[im] + sl[i]
                        do[im] = -1
                        do[i] = 1
                else:
                    im = i - 1
                    ip = i + 1
                    length = sl[i]
                    if sid[im] == -1:
                        length += sl[im]
                        do[im] = -1
                        do[i] = 1
                        sl[im] = sl[i] = length
                    if sid[ip] == -1:
                        length += sl[ip]
                        do[i] = -1
                        do[ip] = 1
                        sl[i] = sl[ip] = length
                        if sid[im] == -1:
                            sl[im] = length
                            sid[im] = -1
    else:

        # merge the sequences (serial). the out of region (-1) will be
        # be processes as a normal sequence id
        for i in schk:                  # loop over the seq need to be merged
            if not do[i] == 0: continue # this seq has been handled so skip it

            if i == 0:                  # first sequence
                do[0] = -1
                do[1] = 1
                # check if seq 0 will merged with seq 1 or vice versa
                _id = sid[1] if sl[1] > threshold or sl[1] > sl[0] else sid[0]
                sid[0] = sid[1] = _id
                sl[0] = sl[1] = sl[0]+sl[1]

            elif i == nsm:              # last sequence
                im = i - 1
                do[im] = -1
                do[i] = 1               # always merge with previous
                sid[i] = sid[im]
                sl[im] = sl[i] = sl[im]+sl[i]

            else:
                im = i - 1
                ip = i + 1
                if sid[im] == sid[ip]: # join adjusent seq with same region
                    do[im]= -1 # mark as merged with the next
                    do[i] = -1
                    do[ip]= 1  #(merged with the previous)
                    sid[i] = sid[im]
                    sl[im] = sl[i] = sl[ip] = sl[im]+sl[i]+sl[ip]
                else:
                    # join with the previous if the seq are of the same id or the
                    # previous seq length is greater that the next seq length
                    if sid[im] == sid[i] or  sl[im] >= sl[ip]:
                        do[im] = -1
                        do[i] = 1
                        sid[i] = sid[im]
                        sl[im] = sl[i] = sl[im]+sl[i]
                    elif sl[im] < sl[ip]: # join with the next if is of geater length
                        do[i] = -1
                        do[ip] = 1
                        sid[i] = sid[ip]
                        sl[ip] = sl[i] = sl[im]+sl[i]
                    else:
                        raise MSDException("error in region time sequnces analysis")

    # the positions where do is 0 or 1 indicate the end indexes of the sequences
    # after merging. A zero is prepended to facilitate the creation of seq reanges
    _w = np.where( (do==0)|(do==1))
    seq = np.insert( np.cumsum( sl[ _w]), 0, 0)
    _iw = np.insert( _w, 0, 0)
    res = []
    for i in range( 1, seq.size):
        im = i-1
        res.append( (seq[im], seq[i], sid[ _iw[i]]))
        #print "seq %d: [%d %d] is of id %d" % (i, _i[i-1], _i[i]-1, sid[_iw[i]])
    return res

def _consecutive(data, stepsize=1):
    ''' Retruns the groups of consecutive elements from the given array. '''
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

def _max_consecutive( data, stepsize=1):
    ''' Retruns the range with the maximum number of consecutive elements
        from the given array. '''
    _a1 = np.diff( data)            # data[i+1]-data[i] i=1,...,
    _a2 = np.where( _a1 != stepsize)[0]+1  # indexes of _a1 where sequences breaks (_a1 non zero)
    if _a2.size > 0:
        _a3 = np.insert( _a2, 0, 0)     # prepend zero and find ( facilitate length calculation)
        _a4 = np.diff( _a3)             # the squences lengths
        imaxsq = _a4.argmax()           # spot the sequence with the maximum length
        imax=_a3[ imaxsq+1]             # spot the end-index of the sequence in the data (_a3 contains sequences indxes)
        imin=_a3[ imaxsq]               # spot the start-index of the sequence in the data
        return (imin,imax)              # return the indexes
    else:
        return (1, data.size-1)
# LDP: this is kept only for archive 
def msd_straight( r):
    ''' Calculate the MDS using the conventional approach and a time
        based loop. The input array should contain atom's coordinates
        in different timesteps:
          r[nconfs,natoms,ndim]
        or:
          r[nconfs,natoms] if we are interested in one dimesion
        or:
          r[nconfs] if we are interested in one atom and one dimension
        where nconfs is the number of configurations, natoms the number
        of atoms and ndim the number of dimensions (dim={1,2,3}).
        Per atom call is this method is not efficient. '''

    n = len( r.shape)
    nconfs = len( r)
    natoms = ndim = 1.0
    if n > 1:
        natoms = np.float32( r.shape[1])
        if n == 3:
            ndim = np.float32( r.shape[2])

    msds = np.zeros( nconfs, dtype=np.float32)
    ncnt = np.zeros( nconfs, dtype=np.int32)
    for i in range( nconfs):
        for j in range( i+1,nconfs):
            d = j - i
            v = r[j] - r[i]
            msds[d] += np.sum( np.square(v))
            ncnt[d] += 1.0

    # normalize
    msds[1:] /= (natoms*ncnt[1:])
    return msds

# LDP: this is kept only for archive 
def msd_shift( r):
    ''' Calculate the MDS using the conventional approach and a shihft
        based loop. The input array
        should contain atom's coordinates in different timesteps:
          r[nconfs,natoms,ndim]
        or:
          r[nconfs,natoms] if we are interested in one dimesion
        or:
          r[nconfs] if we are interested in one atom and one dimension
        where nconfs is the number of configurations, natoms the number
        of atoms and ndim the number of dimensions (dim={1,2,3}).
        The lower is ndim the better is the performance of this
        implementations.
        When natoms is large the method is inefficient either will all
        atoms or per atom call. '''

    # identify input using the shape of the array
    n = len( r.shape)

    shifts = np.arange( len( r))
    msds = np.zeros( shifts.size)
    # Loop over the shifts in pythonic way
    for i, shift in enumerate(shifts):
        diffs = r[:-shift if shift else None] - r[shift:]
        sqdist = np.square(diffs) if n ==1 else np.square(diffs).sum(axis=1)
        msds[i] = sqdist.mean()

    # normalize
    natoms = ndim = 1.0
    if n > 1:
        natoms = np.float32( r.shape[1])
        if n == 3:
            ndim = np.float32( r.shape[2])
    msds *= ( ndim / natoms)
    return msds

def msd_fft( r, atomic=True):
    ''' Calculate the MDS using FFT formulation. The input array may contain
        atoms coordinates in different timesteps:
          r[nconfs,natoms,ndim]
        or:
          r[nconfs,atom] if we are interested in one dimesion (atomic=False)
          r[nconfs,ndim] if we are interested in one atom (atomic=True)
        or:
          r[nconfs] if we are interested in one atom and one dimension
        where nconfs is the number of configurations, natoms the number
        of atoms and ndim the number of dimensions (dim={1,2,3}).
        This method is far more efficient eitherm by all atoms or by per atom
        call.
        http://dx.doi.org/10.1051/sfn/201112010 '''

    # check the input
    n = len( r.shape)

    if n == 3:
        return _msd_fft3( r)
    elif n == 2:
        if atomic:
            return _msd_fft2a( r), None
        else:
            return _msd_fft2( r), None
    elif n == 1:
        return _msd_fft1( r), None

# LDP: this is kept only for archive 
def _msd_fft1( r):
    ''' One atom/One dimension. The argument is expected to be
        one dimension array (r[nconfs]).'''
    N=len(r)
    D = np.square(r)
    S2 = _FFT1D( r, normalize=True)
    Q=2*D.sum()
    S1=np.zeros(N)
    for m in range(N):
        _N = N-m
        Q = Q - D[m-1] - D[_N]
        S1[m] = Q / _N
    return (S1-2*S2)

# LDP: this is kept only for archive 
def _msd_fft2( r):
    ''' Many atoms/One dimension. The argument is expected to be
        two dimension array (r[nconfs,atoms]).'''
    N = len(r)
    natoms = np.float32( r.shape[1])
    D = np.append( np.square(r).sum(axis=1), 0)
    # calculate the ACF for each atom
    S2 = sum( [ _FFTND(r[:, i]) for i in range(r.shape[1])])
    Q=2*D.sum()
    S1=np.zeros(N)
    for m in range(N):
        _N = N-m
        Q = Q - D[m-1] - D[_N]
        S1[m] = Q / _N
    return (S1-2*S2)/natoms

# @profile
def _msd_fft2a( r):
    ''' One atom/Many dimension. The argument is expected to be
        two dimension array (r[nconfs,atoms,n]). The total msd and
        the msd per for each dimension are returned. '''

    N, n = r.shape
    Ds = np.square(r)
    # calculate the ACF for each dimension
    S2 = np.array( _FFTND( r))
    # normalize
    S1s = fasts1x( Ds)
    msds = S1s-2*S2

    return sum( msds), msds

# @profile
def _msd_fft3_chk( r):
    N = len(r)
    n = r.shape[1]
    msds=np.zeros((N,3),dtype=np.float32)
    cnt=np.zeros(N,dtype=np.float32)
    dx=np.zeros((n,r.shape[2]),dtype=np.float32)
    #for i, xi in enumerate(r[0:0]):
    for i, xi in enumerate(r[:-1]):
        for j, xj in enumerate(r[i+1:]):
            dx[:] = (xj-xi)**2
            jp = j+1
            msds[jp,:] += dx.sum(axis=0)
            cnt[jp] += n
    for k in range(r.shape[2]):
        msds[1:,k]/=cnt[1:]
    return  msds.sum(axis=1), (msds[:,0], msds[:,1], msds[:,2])

# @profile
def _msd_fft3x( r, path="", w=None):
    ''' Many atoms/Many dimension. The argument is expected to be
        two dimension array (r[nconfs,atoms,n]). The total msd and
        the msd per for each dimension are returned. In addition
        the file with the msds for each particle is writen in
        path+atomic.pkl pickle file and the weights w for each 
        timestep are calculated to be used in teh fitting of msd(t)
        performed ala gromacs. '''
    # check also:
    # https://stackoverflow.com/questions/34222272/computing-mean-square-displacement-using-python-and-fft/34222273

    N, natoms, n = r.shape

    # calculate position components squares and find their 
    # mean values (over particles)
    _Ds = np.square(r)
    Ds = _Ds.sum( axis=1)
    # calculate the ACF for each particle and find their 
    # mean value (over particles)
    _S2 = [ _FFTND(r[:, i]) for i in range(r.shape[1])]
    S2 = np.sum( _S2, axis=0)

    ###
    # normalize dimension components.-
    # the python code was to slow in per atom msd calculations (this method
    # is called for each atom and become bottleneck (wrap fortran and accelerate)

    # find the msd for each individual particle
    _msds = []
    for i in range(r.shape[1]):
        _S1s = fasts1x( _Ds[:,i,:])
        _msds.append( _S1s-2*np.array( _S2[i]))
    _msds = np.array(_msds)
    # save coordinates and msds for later use
    pickle.dump( ( r, _msds), open(path+"atomic.pkl", "w"))
    # find the fitting weights
    if not w == None:
        w.append( np.std( _msds, axis=0))
    # find the mean msds over particles
    S1s = fasts1x( Ds)
    msds = (S1s-2*S2) / np.float32(natoms)

    print("finished")
    return sum( msds), msds

# @profile
def _msd_fft3( r):
    ''' Many atoms/Many dimension. The argument is expected to be
        two dimension array (r[nconfs,atoms,n]). The total msd and
        the msd per for each dimension are returned. '''
    # check also:
    # https://stackoverflow.com/questions/34222272/computing-mean-square-displacement-using-python-and-fft/34222273
    N = len(r)
    natoms = r.shape[1]
    n = r.shape[2]

    #Ds = np.square(r).sum( axis=1)
    # or save some memory
    Ds = np.zeros((N,n), dtype=np.float32)
    for i in range(N):
        Ds[i,:] = np.square(r[i,:,:]).sum(axis=0)
    #D = Ds.sum(axis=1)

    # calculate the ACF for each particle
    #S2 = np.sum([  _FFTND(r[:, i]) for i in range(r.shape[1])], axis=0)
    # save memory in large systems
    S2 = np.zeros( (3,N), dtype=np.float32)
    for i in range(natoms):
        S2 += _FFTND(r[:, i])

    ###
    # normalize dimension components
    #Qs=2*np.sum( Ds, axis=0)
    #S1s=np.zeros( (n, N))
    ## total is not needed since it is calculated from components
    ##Q=2*D.sum()  # for the total
    ##S1=np.zeros( N)
    #for m in range(N): # check this out something is wrong
    #    _N = N-m
    #    _m = m-1
    #    #Q=Q-D[_m]-D[_N]
    #    #S1[m]=Q/_N
    #    Qs = Qs - Ds[_m,:]-Ds[_N,:]
    #    S1s[:,m] = Qs/_N
    ###
    # the python code was to slow in per atom msd calculations (this method
    # is called for each atom and become bottleneck
    # (wrap fortran and accelerate)
    S1s = fasts1x( Ds)

    msds = (S1s-2*S2) / natoms
    print("finished")
    return sum( msds), msds

# @profile
def _FFTND(x):
    ''' Calculate the ACF for each of the dimensions and add.
        (due to lineariry of fourier transformation). '''
    N=len(x)
    fi = [ _FFT1D(x[:, i]) for i in range(x.shape[1])]
    n = N*np.ones(N)-np.arange(0,N)
    # n = np.arange(N, 0, -1)
    for g in fi: g[:] = g / n # or fi/=n
    return fi

def _FFT1D(x, normalize=False):
    ''' ACF using FFT. '''
    N=len(x)
    F = np.fft.fft(x, n=2*N)      # 2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res= (res[:N]).real           # now we have the autocorrelation in convention B
    if normalize:
        n=N*np.ones(N)-np.arange(0,N) # divide res(m) by (N-m)
        # n=np.arange(N, 0, -1)
        res = res/n
    return res                    # this is the autocorrelation in convention A

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(description='Calculate mean square displacement.')

    # add arguments (self explaned)
    string = 'the path to the simulation log file. In the case of gromacs simulation, a topology file' + \
             'should be present in the same directory (preferably a tpr file).'
    parser.add_argument('path', default="."+os.sep,  \
                       help=string)
    parser.add_argument('-dim', nargs=1, choices=[ 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'], \
                       metavar='dim', default=['xyz'], \
                       help='produce the MSD for the given dimension(s)')
    parser.add_argument('-maxconfs', nargs=1, type=int, metavar='n', default=[-1], \
                       help='pre-allocate memmory for the given number of cofigurations.')
    parser.add_argument('-start', nargs=1, type=int, metavar='n', default=[-1], \
                       help='start processing form configuration START [inclusive]')
    parser.add_argument('-end', nargs=1, type=int, metavar='n', default=[sys.maxsize], \
                       help='stop processing at configuration END [inclusive]')
    parser.add_argument('-every', nargs=1, type=int, metavar='n', default=[1], \
                       help='process every EVERY configuration')

    def argmoleculestype( string):
        ''' check the "-molecules" option arguments. '''
        if len( string) == 0:
            return []
        numbers = isrange( string, positive=True)
        if len( numbers) == 0:
            msg = "wrong molecules indexs range (check: %s)" % string
            raise argparse.ArgumentTypeError(msg)
        return numbers
    parser.add_argument('-molecules', nargs=1, type=argmoleculestype, default=[[]],  metavar='range', \
                       help='molecules to be used. A list with comma seperated id ranges should be provided e.g. "1,2,3" or "1:10,20,30:100"')

    def argslabtype( string):
        ''' check the "-slab" option arguments. '''
        if len(string) == 0:
            return []
        tk = string.split(":")
        msg = ""
        if len( tk) == 2:
            dim = tk[0].lower()
            if dim == 'x' or dim == 'y' or dim == 'z':
                numbers = islist(tk[1], numbertype=float)
                if len( numbers) == 0:
                    msg = "wrong slab(s) position (check: %s)" % tk[1]
            else:
                msg = "wrong slab normal direction (check: %s)" % tk[0]
        else:
            msg = "wrong slab flag argument (check: %s)" % string

        if len( msg) > 0:
            raise argparse.ArgumentTypeError(msg)

        values = [ dim , numbers ]
        return values

    string = 'slabs to be used. The dimension perpedicular to the slabs and their positions should be ' + \
        'provided i.e. [x|y|z]:r0,r1,r2,... , where ri is the starting position of the ith slab.'
    parser.add_argument('-slab', nargs=1, type=argslabtype, default=[[]], metavar='slab', \
                       help=string)

    parser.add_argument('-gmxfit', dest='gmxfit', action='store_true',help="use gromacs approach for fitting msd(t) curve.")

    string = 'the file with the molecular indexes (starting from 1) to consideted in the calculation of the msd. ' + \
        'The file should conform with gromacs format. All the groups found in the file will be considered.'
    parser.add_argument('-ndx', nargs=1, type=argparse.FileType('r'), metavar='file', default=[None], required=False, help=string)             

    # parse the arguments
    args = parser.parse_args()

    print("INPUT")
    print("path : ", args.path)
    print("start : ", args.start[0])
    print("end : ", "max" if args.end[0] == sys.maxsize else args.end[0])
    print("every : ", args.every[0])
    print("dim : ", args.dim[0])
    print("slab : ", args.slab[0])
    string = " ".join( map(str, args.molecules[0]))
    print("molecules : %s" % ('-' if len(string)==0 else string ) )
    print("ndx : ", "-" if args.ndx[0] is None else args.ndx[0].name)

    if __debug:
        print(args.molecules[0])
    msd( args.path, start=args.start[0], end=args.end[0], every=args.every[0], dimensions=args.dim[0], \
         region=args.slab[0], molids=args.molecules[0], maxconfs=args.maxconfs[0], ndxfile=args.ndx[0])

if __name__ == '__main__':
    command()
