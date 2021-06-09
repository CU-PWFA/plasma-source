#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 12:16:20 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
#%%
#Density= [3]
#Power= [5]
#PWidth= [50]
PixelSize= 10

Density= [2, 3, 5, 7, 8, 10, 20, 50, 80]
Power= [5]
PWidth= [30, 35, 40, 45, 50]

FPW_Results2D= np.zeros((len(Density), len(PWidth)))
PTValue_Results2D= np.zeros((len(Density), len(PWidth)))

denC= 0
PWC= 0
for den in Density:
    for PW in PWidth:
        for p in Power:
            try:
                dataE= np.load('PW'+str(PW)+'_GP'+str(p)+'_n'+str(den)+'.npy')
                MidLineOutI= abs(dataE[:, int(dataE.shape[1]/2)])**2

                x= MidLineOutI
                peaks, _ = find_peaks(x)
                peak1idx, peak2idx= peaks[np.argsort(x[peaks])[-1]], peaks[np.argsort(x[peaks])[-2]]
                peak1, peak2= x[peak1idx], x[peak2idx]
                
                if (np.argsort(x[peaks])[-1])-(np.argsort(x[peaks])[-2])==2:
                   
                    valleys= (np.diff(np.sign(np.diff(x[min(peak1idx, peak2idx): max(peak1idx, peak2idx)])))>0).nonzero()[0]+1\
                            +min(peak1idx, peak2idx)
                    valley1idx, valley2idx= valleys[np.argsort(x[valleys])[0:2]]
                    valley1, valley2= x[np.array([valley1idx, valley2idx])]
                else:
                    valleys= (np.diff(np.sign(np.diff(x[min(peak1idx, peak2idx): max(peak1idx, peak2idx)])))>0).nonzero()[0]+1\
                            +min(peak1idx, peak2idx)
                    valley1idx= valleys[np.argsort(x[valleys])[0]]
                    valley2idx= valley1idx
                    valley1, valley2= x[np.array([valley1idx, valley2idx])]
                    
                if (np.mean(np.array([peak1, peak2])))-(np.mean(np.array([valley1, valley2])))>0.5:
                    
                    plt.plot(x)
                    plt.plot(np.array([peak1idx, peak2idx]), x[np.array([peak1idx, peak2idx])], 'x')
                    plt.show()     
                    plt.plot(np.array([valley1idx, valley2idx]), x[np.array([valley1idx, valley2idx])], 'x')
    
                    FirstPeakWidth= abs(peak1idx-peak2idx)*PixelSize
                    PTValue= 100*(np.mean(np.array([peak1, peak2]))-np.mean(np.array([valley1, valley2])))/ \
                                ((np.mean(np.array([peak1, peak2]))+np.mean(np.array([valley1, valley2])))/2)
                    FPW_Results2D[denC, PWC]= FirstPeakWidth
                    PTValue_Results2D[denC, PWC]= PTValue
                    
            except:
                pass
        PWC= PWC+1
    denC= denC+1
    PWC= 0
    print(denC, PWC)

#%%
PTValue_Results2D_x= np.array(PWidth)
PTValue_Results2D_y= np.array(Density)

np.save('PTValue_Results2D', PTValue_Results2D)
np.save('PTValue_Results2D_x', PTValue_Results2D_x)
np.save('PTValue_Results2D_y', PTValue_Results2D_y)