#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 16:54:11 2020

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift, fftfreq, fft2, ifft2;
from PIL import Image
from scipy.ndimage import gaussian_filter
from scipy import optimize
#%%
PeakIdxLArray= np.load('data.npy')
plt.plot(PeakIdxLArray[:, 0], PeakIdxLArray[:, 1], '.')
plt.xlabel('x (pixel#)')
plt.ylabel('y (pixel#)')

#%%
index= [i for i in range(int(PeakIdxLArray.shape[0]-1))]
your_data= PeakIdxLArray[:, 1]
count=0
NewLArray= np.zeros((1, 2))

while len(index)>0:
    print(len(index))
    index_processing = index.pop(0) 
    if your_data[index_processing+1]- your_data[index_processing]<12:
        NewLArray[count][0]= PeakIdxLArray[index_processing+1, 0]
        NewLArray[count][1]= PeakIdxLArray[index_processing+1, 1]
        NewLArray= np.append(NewLArray, [[0, 0]], axis=0)
        count= count+1
    else:
        mark = False
        while mark == False:
            if len(index)==0:
                break
            index_processing = index.pop(0)
            if your_data[index_processing+1]- your_data[index_processing]<12:
                NewLArray[count][0]= PeakIdxLArray[index_processing+1, 0]
                NewLArray[count][1]= PeakIdxLArray[index_processing+1, 1]
                NewLArray= np.append(NewLArray, [[0, 0]], axis=0)
                count= count+1
                mark = True
            
#%%
plt.plot(NewLArray[:,0], NewLArray[:,1], '.')
