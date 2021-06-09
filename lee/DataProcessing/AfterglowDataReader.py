#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 19:30:42 2020

@author: valentinalee
"""

#%%
import numpy as np;
import matplotlib.pyplot as plt
from PIL import Image

#%%
test= np.array(Image.open('2009210107/17583372_2009210107_0004.tiff'))
#%%
#testNew= test
testNew= test[50:250, 400:600]
pixel= 3.45
x= np.linspace(-testNew.shape[0]/2*pixel, testNew.shape[0]/2*pixel, testNew.shape[0])
y= np.linspace(-testNew.shape[1]/2*pixel, testNew.shape[1]/2*pixel, testNew.shape[1])
Nortest= testNew/(np.amax(testNew))
plt.pcolormesh(y, x, Nortest)
#plt.pcolormesh(y, x, testNew)
plt.axis('scaled')
plt.colorbar()
plt.xlabel('x ($\mu$m)')
plt.ylabel('y ($\mu$m)')
#%%
plt.pcolormesh(test)