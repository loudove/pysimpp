# -*- coding: utf-8 -*-

from math import fabs, acos, copysign, sqrt
import numpy as np

def get_square(a):
    ''' Vector square length. '''
    return np.dot(a, a)

def get_length(a):
    ''' Vector (iteratable) a length. '''
    #return sqrt(sum(a * a))
    return sqrt( np.dot(a, a))

def get_unit(a):
    ''' Return the unit vector .'''
    length = get_length(a)
    return np.array(a) / length

def get_angle(a, b):
    ''' Return the angle between vectors a and b. '''
    d = np.dot(a, b) / get_length(a) / get_length(b)
    abs_d = fabs(d)
    if abs_d > 1.0:
        if (abs_d - 1.0) < 0.000001:  d = copysign(1.0, d)
    return acos(d)

def get_angle_unit(a, b):
    ''' Return the angle between unit vectors a and b. '''
    d = np.dot(a, b)
    abs_d = fabs(d)
    if abs_d > 1.0:
        if (abs_d - 1.0) < 0.000001:  d = copysign(1.0, d)
    return acos(d)

def get_dihedral(a, b, c):
    ''' Return the dihedral between vectors a, b and c. '''
    # cis zero righ-handed
    axb = get_unit(np.cross(a, b))
    bxc = get_unit(np.cross(b, c))

    d = get_angle(axb, bxc)
    s = np.dot(axb, c)

    return copysign(d, s)

def get_projection(a, b):
    ''' Return the projection of vector b to a. '''
    p = np.dot(a, b) / np.dot(a, a)
    return np.array(a) * p

def get_vertical(a, b):
    ''' Return a vertor vertical to vector a. The new vector lies on
        plane defined form a and b toward. '''
    c = get_projection(a, b)
    return np.array(b) - c

def get_next(a, b, c, l, th, phi):
    pass
