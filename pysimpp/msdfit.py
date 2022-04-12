# -*- coding: utf-8 -*-

import os
import sys
import re

import numpy as np
from scipy import stats

def _is_command(): return True
def _short_description(): return 'Calculate diffusion coefficient from msd.py or gromacs msd output.'
def _command(): command()

_m5u = 'x10\u207b\u2075 cm\u00b2/s'
_skip = False

# # utility functions for writing
# def __strfmt( n):
#     return " %"+ "%d"%n + "s "

# def __appendcolumn( lines, head, data):
#     if not len(lines) == len(data) + 1: return
#     _n = len( head)
#     _n = max( _n, np.max( data))
#     fmt = __strfmt( _n)
#     lines[0] += fmt % head
#     for i, line in enumerate( lines[1:]):
#         line += fmt % str( data[i])

#########################################################
# https://github.com/linsomniac/python-movingaverage/blob/master/movingaverage.py
def __moving_average(data, subset_size, data_is_list = None,
                    avoid_fp_drift = True):
    ''' Return the moving averages of the data, with a window size of
        `subset_size`.  `subset_size` must be an integer greater than 0 and
        less than the length of the input data, or a ValueError will be raised.
        `data_is_list` can be used to tune the algorithm for list or iteratable
        as an input.  The default value, `None` will auto-detect this.
        The algorithm used if `data` is a list is almost twice as fast as if
        it is an iteratable.
        `avoid_fp_drift`, if True (the default) sums every sub-set rather than
        keeping a "rolling sum" (which may be subject to floating-point drift).
        While more correct, it is also dramatically slower for subset sizes
        much larger than 20.
        NOTE: You really should consider setting `avoid_fp_drift = False` unless
        you are dealing with very small numbers (say, far smaller than 0.00001)
        or require extreme accuracy at the cost of execution time.  For
        `subset_size` < 20, the performance difference is very small.
    '''
    if subset_size < 1:
        raise ValueError('subset_size must be 1 or larger')

    if data_is_list is None:
        data_is_list = hasattr(data, '__getslice__')

    divisor = float(subset_size)
    if data_is_list:
        #  This only works if we can re-access old elements, but is much faster.
        #  In other words, it can't be just an iterable, it needs to be a list.

        if subset_size > len(data):
            raise ValueError('subset_size must be smaller than data set size')

        if avoid_fp_drift:
            for x in range(subset_size, len(data) + 1):
                yield sum(data[x - subset_size:x]) / divisor
        else:
            cur = sum(data[0:subset_size])
            yield cur / divisor
            for x in range(subset_size, len(data)):
                cur += data[x] - data[x - subset_size]
                yield cur / divisor
    else:
        #  Based on the recipe at:
        #     http://docs.python.org/library/collections.html#deque-recipes
        it = iter(data)
        d = deque(islice(it, subset_size))

        if subset_size > len(d):
            raise ValueError('subset_size must be smaller than data set size')

        if avoid_fp_drift:
            yield sum(d) / divisor
            for elem in it:
                d.popleft()
                d.append(elem)
                yield sum(d) / divisor
        else:
            s = sum(d)
            yield s / divisor
            for elem in it:
                s += elem - d.popleft()
                d.append(elem)
                yield s / divisor

def _moving_average(data_set, periods=3):
    weights = np.ones(periods) / periods
    return np.convolve(data_set, weights, mode='valid')

def _printfit( slope, intercept, r_value, p_value, std_err, istart, tstart, iend, tend, var, d, what):
    # _m10u = u'10\u207b\u00b9\u2070 m\u00b2/s'
    return  '# %s: points %d:%d (slope variance %.2f)\n' % ( what, istart, iend , var) + \
            '# y = a x + b : a = %f, b = %f\n' % ( slope, intercept) + \
            '# r-value = %f , p-value = %f stderr = %f\n' % (r_value, p_value, std_err) + \
            '# D = %f %s' % (slope/float(2*d)*1.e4, _m5u)

    #return  u'# %s: points %.0f:%.0f (slope variance %.2f)\n' % ( what, tstart, tend , var) + \
    #        u'# D = %f %s' % (slope/float(2*d)*1.e4, "-")

def _max_consecutive( data, stepsize=1):
    ''' Retruns the range with the maximum number of consecutive elements
        from the given array. '''
    _a1 = np.diff( data)            # data[i+1]-data[i] i=1,...,
    _a2 = np.where( _a1 != stepsize)[0]+1  # indexes of _a1 where sequences breaks (_a1 non zero)
    _a3 = np.insert( _a2, 0, 0)     # prepend zero and find ( facilitate length calculation)
    _a4 = np.diff( _a3)             # the squences lengths
    imaxsq = _a4.argmax()           # spot the sequence with the maximum length
    if _skip:
      _a4[imaxsq] = 1
      imaxsq = _a4.argmax()         # spot the second sequence 
    imax=_a3[ imaxsq+1]             # spot the end-index of the sequence in the data (_a3 contains sequences indxes)
    imin=_a3[ imaxsq]               # spot the start-index of the sequence in the data
    return (imin,imax)              # return the indexes

def _fit( t, msd, var=0.1, range=[], smooth=1):
    '''Fit the data to calculate self diffusion coefficient. Returns a
        tuple with fit results:
            (slope, intercept, r_value, p_value, std_err, start, end)
    '''
    # exclude first time where msd is zero
    _t = t - t[0]
    first = 0 if not t[0] == 0.0 else 1
    last = len(t)
    logt = np.log10( _t[first:last])
    logmsd = np.log10( msd[ first:last])
    _grad = np.gradient( logmsd, logt)

    # smooth the gradient setting the first values to zero 
    if smooth > _grad.size/2: smooth=1
    if smooth > 1:
        grad = np.zeros(_grad.size)    
        #grad[smooth-1:]=np.array(list(__moving_average(_grad,smooth)))
        grad[smooth-1:]=np.array(list(_moving_average(_grad,smooth)))
    else:
        grad=_grad

    # trace where the 1-var<grad<1+var
    _pone = 1.0 + var
    _mone = 1.0 - var
    _where = np.where( (grad<_mone)|(grad>_pone))
    ranges = np.arange( grad.size)+1
    ranges[ _where] = -1
    # remove the very first part
    if smooth == 1:
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
        if len(range) > 0:
            istart = np.where( t >= range[0])[0][0]
            iend = np.where( t >= range[1])[0][0]
            if iend > istart:
                slope, intercept, r_value, p_value, std_err = stats.linregress( t[istart:iend], msd[istart:iend])
                ok = True
        else:
            istart, iend = _max_consecutive( ranges)
            if iend - istart > 2: # more than two points are needed
                istart, iend = istart+1, iend+1  # increas by one to get the correct indexes
                slope, intercept, r_value, p_value, std_err = stats.linregress( t[istart:iend], msd[istart:iend])
                ok = True         # set success flag
    except:
        pass
    finally:
        if not ok:
            istart, iend = 0, 0

    return [slope, intercept, r_value, p_value, std_err, istart, iend, grad]

def msdfit( filename, dim, visualize, var, range=[], plrange=[], intype='gmx', every=1, column=1, smooth=1):

    # get file info
    DIR=os.path.dirname(filename)
    ( DIR, DATAFILE ) = os.path.split( filename)
    if len(DIR) == 0: DIR="."
    ( BASENAME, EXT ) = os.path.splitext( DATAFILE)


    # take care units stuff. in gmx depends on used input
    # TODO: add an option for the input time units
    if intype == 'gmx':
        _tcnv=1000.0  # set the default converter for time
        _dcnv=100.0
        _rcnv=1.0
    elif intype == 'msdout':
        _tcnv=1.0
        _dcnv=1.0
        _rcnv=1000.0
    if not len(range) == 0:
        range=tuple(np.array(range)*_rcnv)
    if not len(plrange) == 0:
        plrange=tuple(np.array(plrange)*_rcnv)        
    
    # check what is going to be calculated
    indxs=[]
    for d in dim:
        if d == 'x': indxs.append(2)
        elif d == 'y': indxs.append(3)
        elif d == 'z': indxs.append(4)
    dof = len(dim)
    if intype == 'gmx':
        indxs=[1]

    # read the msd data
    f = open(filename,'r',encoding="utf-8")
    lines_ = [x.strip() for x in f.readlines()]
    lines = [x for x in lines_ if not x.startswith(("#","@"))]

    if intype == "gmx":
        timere=re.compile(r'@[ ]*xaxis  label "Time \((ns|ps|fs)\)"')
        gmxtime=""
        for line in lines_:
            match = timere.match( line)
            if line.startswith(("#","@")):
                if match and match.groups(1):
                    gmxtime=match.groups(1)[0].lower()
                    if gmxtime=="ns":
                        _tcnv=1.0e6
                    elif gmxtime=="ps":
                        _ctnv=1.0e3
                    elif gmxtime=="fs":
                        _ctnv=1.0
                    else:
                        print("ERROR: unexpected time units in %s" % DATAFILE)
                        return
            else:
                break

    f.close()
    t = []
    msd = []
    if len(lines) < 5:
        print("to few data found in %s" % file)
        return
    if len(lines) < 10*every-smooth: every=1
    for i, line in enumerate(lines):
        if i%every == 0: 
            tk = np.array(list(map( float, line.split())))
            t.append( tk[0]*_tcnv)
            msd.append( np.sum(tk[indxs])*_dcnv)
    t = np.array( t)
    msd = np.array( msd)

    # fit the data
    slope, intercept, r_value, p_value, std_err, istart, iend, grad = _fit(t, msd, range=range, var=var, smooth=smooth)

    # print the results
    res = _printfit( slope, intercept, r_value, p_value, std_err, istart, t[istart], iend, t[iend], var, dof, "msd %s "%dim)
    print(res)

    # plot stuff
    if visualize:
        import matplotlib.pyplot as plt
        from matplotlib.ticker import NullFormatter  # useful for `logit` scale

        first=0
        last=t.size
        if not len(plrange) == 0:
            _a = np.where( t > plrange[0]) 
            if not len(_a) == 0: first = _a[0][0]
            _a = np.where( t < plrange[1])
            if not len(_a) == 0: last = _a[0][-1]+1
        length = last - first
    
        t_ = np.array(t)

        plt.figure(1,figsize=[11.8, 4.2],dpi=100)
        # msd(t)
        plt.subplot(141)
        y = np.zeros(t.size)
        y[istart:iend+1] = slope*t[istart:iend+1]+intercept
        t_-=t[0]
        t_/=1.e6
        plt.plot(t_[first:last], msd[first:last],'k',t_[first:last],y[first:last],'r')
        plt.yscale('linear')
        plt.xlabel('t [ns]')
        plt.ylabel('msd [A^2]')
        plt.title('msd(t)')
        plt.grid(True)

        # grad(t)
        plt.subplot(142)
        #x = np.log10(t_[1:])
        x = t_[1:]
        y = np.zeros(grad.size)
        y[istart:iend+1] = 1.0
        y1 = np.ones(grad.size)
        plt.plot(x[first:last], grad[first:last],'k',x[first:last],y[first:last],'r--',x[first:last],y1[first:last],'g')
        plt.xscale('linear')
        plt.yscale('linear')
        #plt.xlabel('log10(t) [t in ns]')
        plt.xlabel('t [ns]')
        # plt.xlabel('log10(t) [ps]')
        plt.ylabel('grad(Lmsd(Lt)) [ns^2/ps], L=log10')
        plt.title('grad(Lmsd(Lt))')
        plt.grid(True)

        # msd/6t
        plt.subplot(143)
        x = t_
        y0 = np.zeros(t.size)
        y0[1:] = msd[1:]/t[1:]/float(2*dof)*1.e4
        y0[0] = y0[1]
        y = np.ones(t_.size)
        y[:] = slope/float(2*dof)*1.e4
        plt.plot(x[first:last], y0[first:last],'k',x[first:last],y[first:last],'r--')
        plt.xscale('linear')
        plt.yscale('linear')
        plt.xlabel('t [t in ns]')
        plt.ylabel('msd(t)/(2*d*t) [%s]'%_m5u)
        #plt.title('msd(t)/(2*d*t)\nd=%.4f %s'%(y[0],_m5u))
        plt.title('d=%.4f %s'%(y[0],_m5u))
        plt.grid(True)

        # log(t)-log(msd)
        plt.subplot(144)
        x = np.log(t_[1:])
        y0 = np.log(msd[1:])
        y1 = np.log(slope*1.e6) + x # set the slope to the units of the diagram
        plt.plot(x[first:last], y0[first:last],'k',linewidth=2)
        plt.plot(x[first:last],y1[first:last],'r--',linewidth=1)
        #plt.plot(x[first:last], y0[first:last],'k', x[first:last],y1[first:last],'r--')
        plt.xscale('linear')
        plt.yscale('linear')
        plt.xlabel('ln(t) [t in ns]')
        plt.ylabel('ln(msd(t)) [msd in A^2]')
        plt.title('ln(msd)[ln(t)]')
        plt.grid(True)

        plt.gca().yaxis.set_minor_formatter(NullFormatter())

        plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25, wspace=0.35)

        plotname=DIR+os.sep+BASENAME+".png"
        print("(you can find the plot in: %s)" % plotname)
        plt.savefig(plotname, dpi=400, bbox_inches='tight')

        plt.show()

def command():

    import argparse

    # create an argument parser
    parser = argparse.ArgumentParser(description='Calculate diffusion coefficient from gromacs msd.')

    # add arguments (self explaned)
    parser.add_argument('path', default="."+os.sep, help='the path to the gromacs msd file.')
    def argdim( string):
        ret = ''
        for d in string:
            if not d in ('x','y','z'):
                msg = "wrong dim argument (%s)" % string
                raise argparse.ArgumentTypeError(msg)
            ret+=d
        if len(ret) == 0: ret='xyz'
        return ret               
    parser.add_argument('-dim', nargs=1, type=str, metavar='d', default=['xyz'],\
                        help='msd dimemsions.')
    parser.add_argument('-var', nargs=1, type=float, metavar='v', default=[0.1], \
                        help='slope divergence from one (defult value is 0.1).')
    parser.add_argument('-every', nargs=1, type=int, default=[1], \
                        help='use one every EVERY data points.')  
    parser.add_argument('-smooth', nargs=1, metavar='N', type=int, default=[1], \
                        help='use N number of point to comput the mooving average msd(t) derivative (smoothening).')                                       
    parser.add_argument('-vis', dest='visualize', action='store_true',help="visualize various msd plots.")
    parser.add_argument('-skip', dest='skip', action='store_true',help="skip the first region traced for msd(t) fitting.")
    parser.add_argument('-type', nargs=1, choices=['msdout','gmx'], default=['gmx'],\
                        help='the type of the input file to process.')
    def argrange( string):
        ''' check the "-range" option arguments. '''
        if len( string) == 0:
            return []
        ok = False
        try:
            numbers = list(map(float, string.split(':')))
            if len(numbers) == 2 and numbers[0] < numbers[1]: ok = True
        except:
            pass
        if not ok :
            msg = "wrong molecules indexs range (check: %s)" % string
            raise argparse.ArgumentTypeError(msg)
        return numbers
    parser.add_argument('-frange', nargs=1, type=argrange, default=[[]],  metavar='range', \
                       help='fit range (time scale in ps) to be used for the fit. A list with two column separated values should be provided e.g. 1:10.')
    parser.add_argument('-plrange', nargs=1, type=argrange, default=[[]],  metavar='range', \
                       help='fit range (time scale in ps) to be used for the plots. A list with two column separated values should be provided e.g. 1:10.')

    # parse the arguments
    args = parser.parse_args()

    _skip = args.skip

    msdfit( args.path, args.dim[0], args.visualize, args.var[0], range=args.frange[0], plrange=args.plrange[0], intype=args.type[0], every=args.every[0], smooth=args.smooth[0])

if __name__ == '__main__':
    command()
