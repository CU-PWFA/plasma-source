#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 13:58:20 2021

@author: valentinalee
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import chirp, find_peaks, peak_widths
#%%
def find(condition):
    res, = np.nonzero(np.ravel(condition))
    return res

def FWHM(X,Y):
    
    '''
    Args:
        X: 1D array, x variable
        Y: 1D array, y= f(x)
    
    Returns:
        value: the FWHM of signal f(x)
    
    '''

    half_max = max(Y) / 2
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    left_idx = find(d > 0)[0]
    right_idx = find(d < 0)[-1]
    print(half_max, X[right_idx], X[left_idx])
    return abs(X[right_idx] - X[left_idx])
#%%
c= 3e8
#%%
Start= -16.5
End= -16.75
Steps= 333
record= np.empty((int(Steps-1), 2))
record[:, 0]= [x for x in range(0, int(Steps-1))]
record[:, 1]= [Start+(End-Start)/Steps*x for x in range(0, int(Steps-1))]
#%%
Count= []
Delay= []
#path= '~/GitHub/plasma-source/lee/BobTheFrog'
path= ''
DataSetNumber= '2111050070'
for data_number in range(len(record)):
    data= np.load(path+ ''+ DataSetNumber+ '/michaelito_'+ DataSetNumber+ '_'+\
                  str(int(record[int(data_number), 0])).zfill(4)+ '.npy', allow_pickle=True)
    
    lam= data.item()['lambda']
    count= data.item()['I']
    Count.append(count)
    Delay.append(record[int(data_number)][1])
    
Count= np.array(Count)
Delay= np.array(Delay)
#%%
plt.pcolormesh(Count)

#%%
#plt.pcolormesh(lam, Delay, Count, vmax= 6000)
#plt.pcolormesh(lam[625: 875], Delay, Count[:, 625: 875], vmax= 6000)
#plt.figure(2)
plt.pcolormesh(lam[712:766], Delay*2*1e-3/c*1e15+110900, Count[:, 712:766])#, vmax= 12000)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Delay (fs)')
#plt.title('Grating Location= 2000')
plt.colorbar()
#%%
wavelength= lam[712:766]
time= Delay*2*1e-3/c*1e15+110900
Int= Count[:, 712:766]
#%%
for shot in range(len(time)-1):
    if np.sum(Int[shot, :])>np.sum(Int[shot+1, :])*1.2:
        Int[shot, :]= Int[shot+1, :]
#%%
plt.figure(2)
plt.pcolormesh(wavelength, time, Int)#, vmax= 12500)
plt.colorbar()
plt.xlabel('Wavelength (nm)')
plt.ylabel('Delay (fs)')
plt.title('Compressor at -15 steps')        
#%%
Auto= np.sum(Int, axis= 1)
x= Auto
peaks, _= find_peaks(x, height=26000, width=20)
results_half= peak_widths(x, peaks, rel_height=0.5)
print(results_half)
#%%
AutoResults=[]
#%%
AutoResults.append(Auto)

#%%
FWHMsave= []
#for w in range(115, 116):
for w in range(len(wavelength)):
    x= Int[:, w]
    peaks, _= find_peaks(x, height=3000, width=10)
    results_half= peak_widths(x, peaks, rel_height=0.5)
    if len(results_half[0])== 0:
        FWHMsave.append(1000)  # widths
    else:
        
        FWHMsave.append(float(results_half[0][0]))  # widths

#%%


#%%
x = Auto

peaks, properties = find_peaks(x, prominence=1, width=10, height= 250000)

properties["prominences"], properties["widths"]

plt.plot(x)

plt.plot(peaks, x[peaks], "x")

plt.vlines(x=peaks, ymin=x[peaks] - properties["prominences"],

           ymax = x[peaks], color = "C1")

plt.hlines(y=properties["width_heights"], xmin=properties["left_ips"],

           xmax=properties["right_ips"], color = "C1")
print(properties["widths"])
plt.show()
plt.title('-15')
#%%
print(FWHM(time, Int[:, 116]))
#%%
plt.plot(wavelength, FWHMsave)
plt.xlabel('nm')
plt.ylabel('FWHM')