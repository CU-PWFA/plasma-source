#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 14:16:52 2022

I want to load a Matlab file and just print out all the bs

@author: chris
"""

import h5py as h5
import numpy as np
import os as os
import matplotlib.pyplot as plt
#from scipy import optimize
from scipy.optimize import curve_fit
import sys
import scipy.io as sio

superpath = '/media/chris/New Volume/SLACData/'
day = '20220812/'
dataset = 'E308_02493' #no gas
dataset = 'E308_02500' #24 psi
path = superpath + day + dataset + '/' + dataset +'.mat'


mat = sio.loadmat(path)

data = mat['data_struct']

"""
sys.stdout = open(superpath+day+dataset+'/'+dataset+'.txt','w')
for key,value in mat.items():
    sys.stdout.write(str(key)+ ':'+ str(value))
sys.stdout.close()
"""

backs = data['backgrounds'][0][0]
dtotr2_back = backs['DTOTR2']