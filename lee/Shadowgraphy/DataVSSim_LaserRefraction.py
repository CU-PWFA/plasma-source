#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 18:08:22 2020

@author: valentinalee
"""

#%%
import numpy as np;
import matplotlib.pyplot as plt
from PIL import Image

#%%
im064 = Image.open('064/19423601_2001230064.tiff')
im065 = Image.open('065/19423601_2001230065.tiff')
im080 = Image.open('080/19423601_2001230080.tiff')
On_8e16_500= np.array(im064) 
Off_8e16_500= np.array(im065)
On_2_5e16_500= np.array(im080) 

#%%
pixel=3.45e-6
xaxis=np.linspace(-1023, 1024, 2048)*pixel*10**3/0.4
yaxis=np.linspace(-767, 768, 1536)*pixel*10**3/0.4

#%%
lineout= On_8e16_500[:, 930]
plt.plot(lineout)
#%%
plt.figure(2)
#plt.pcolormesh(On_2_5e16_500)
plt.pcolormesh(xaxis, yaxis, On_8e16_500)
plt.axis('scaled')
