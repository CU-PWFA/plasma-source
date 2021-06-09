#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 10:08:44 2019

@author: valentina_lee
"""

#%%
import numpy as np;
import matplotlib.pyplot as plt
from PIL import Image

#%%
im = Image.open('17583372_1910290007_0004.tiff')
PickOff= np.array(im)

#%%
pixel=3.45e-6
Energy=1.165e-3
NorE=PickOff*Energy/sum(sum(PickOff))
NorInt=(NorE/30e-15)/((pixel)**2)

xaxis=np.linspace(-1023, 1024, 2048)*pixel*10**3
yaxis=np.linspace(-767, 768, 1536)*pixel*10**3

NorInt=NorInt/10000/1e12
#%%
plt.figure(1)
plt.pcolormesh(xaxis, yaxis, NorInt)
plt.axis('scaled')
plt.colorbar(label='TW/cm^2')
plt.xlabel('mm')
plt.ylabel('mm')
#%%
np.save(PickOff, '2e17cm-3_170mJ')