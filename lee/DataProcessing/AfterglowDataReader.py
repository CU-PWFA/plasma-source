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
test1= np.array(Image.open('2109160180/17583372_2109160180_0000.tiff'))
BG1= np.array(Image.open('2109160241/17583372_2109160241_0000.tiff'))
test2= np.array(Image.open('2109160180/19423598_2109160180_0000.tiff'))
BG2= np.array(Image.open('2109160181/19423598_2109160181_0000.tiff'))
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
plt.figure(1)
plt.pcolormesh(test2)
plt.figure(2)
plt.pcolormesh(BG2)