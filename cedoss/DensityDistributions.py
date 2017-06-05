# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 13:03:37 2017

Various density distributions for use with creating thin plasma lens

@author: chris
"""

import numpy as np

#Creates a density distribution which is exponential axially and Gaussian
# radially.  The location of the gas jet is at z=-z_offset, and any z values
# below this threshold are automatically set to zero.
#  n0 - density at nozzle
#  Lz - scale length of exponential profile in z
#  Lr - scale length of Gaussian profile in r
#  z_offset - distance between nozzle and origin (center of beam, for example)
#  x,y,z - Cartesian coordinates for which to calculate density
def SimpleGasJet(n0,Lz,Lr,z_offset,x,y,z):
    if z<(-z_offset):
        return 0
    else:
        nr=np.exp(-(np.power(x,2)+np.power(y,2))/(2*np.power(Lr,2)))
        nz=np.exp(-(z+z_offset)/Lz)
        return n0*nr*nz