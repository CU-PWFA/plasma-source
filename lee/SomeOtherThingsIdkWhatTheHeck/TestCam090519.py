#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 15:33:16 2019

@author: valentina_lee
"""

#%%
import numpy as np;
import matplotlib.pyplot as plt
from PIL import Image

#%%
im029 = Image.open('029/19423601_2001230029.tiff')
















PickOff= np.array(im029)

#%%
plt.figure(1)
plt.pcolormesh(PickOff)
#plt.title('$2e17cm^{-3}$ 170mJ')
#plt.axis('scaled')

#%%
np.save(PickOff, '2e17cm-3_170mJ')