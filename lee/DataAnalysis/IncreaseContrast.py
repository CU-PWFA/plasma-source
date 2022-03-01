#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 18:57:52 2021

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from scipy.ndimage import gaussian_filter
#%%
Data= np.sqrt(np.array(Image.open('19423601_2009210001_0006.tiff')))
After_filter_Data= gaussian_filter(Data, sigma=3)

plt.figure(1)
plt.pcolormesh(Data)
plt.colorbar()

plt.figure(2)
plt.pcolormesh(After_filter_Data)
plt.colorbar()

#%%
DataInc= After_filter_Data**4
plt.figure(1)
plt.pcolormesh(Data)
plt.colorbar()

plt.figure(2)
plt.pcolormesh(DataInc)
plt.colorbar()
#%%
Data1D= Data[:, 312]/np.amax(Data[:, 312])
DataFT1D= After_filter_Data[:, 312]/np.amax(After_filter_Data[:, 312])
DataInc1D= DataInc[:, 312]/np.amax(DataInc[:, 312])

plt.figure(3)
plt.plot(Data1D, label= 'Raw')
plt.plot(DataInc1D, label= 'Processed Data')
plt.plot(DataFT1D, label= 'Gaussian Filtered Data')
plt.legend()