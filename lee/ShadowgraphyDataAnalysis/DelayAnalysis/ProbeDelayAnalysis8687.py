#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 18 16:14:41 2021

@author: valentinalee
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from photutils.centroids import centroid_com
from PIL import Image
from scipy.ndimage.interpolation import shift
from scipy.ndimage import gaussian_filter
#%%
#load background data, find centroid 
BG= np.array(Image.open('2109170087/19423598_2109170087_0000.tiff'))
BG_cx, BG_cy = centroid_com(BG)

#move BG so that centroid goes to center
BG= shift(BG, ((300-BG_cy), (300-BG_cx)))

#load data, find centroid
Data= np.array(Image.open('2109170086/19423598_2109170086_0000.tiff'))
Data_cx, Data_cy = centroid_com(Data)
 #if it did work, use a low pass filter first 
#move data so that centroid goes to center
Data= shift(Data, ((300-Data_cy), (300-Data_cx)))

#sutraction 
Signal= gaussian_filter(Data-BG, sigma=5)
#%%
plt.figure(1)
plt.pcolormesh(BG)
plt.figure(2)
plt.pcolormesh(Data)
#%%
plt.figure(3)
plt.pcolormesh(Signal)

#%%
SignalRangeX= [190, 440]
SignalRangeY= [280, 360]
CentralIdx= (SignalRangeY[1]-SignalRangeY[0])/2
PeakIdx1= []
PeakIdx2= []
ValleyIdx= []
Peak1= []
Peak2= []
Valley= []
#%%
#for XIdx in range (300, 301):
for XIdx in range (SignalRangeX[0], SignalRangeX[1]):
    print(XIdx)
    x= Signal[SignalRangeY[0]: SignalRangeY[1], XIdx]
    peaks, _ = find_peaks(x, height=4000, distance=20)

    x2= -x
    valleys, _ = find_peaks(x2)

    try:
        valleyidx= valleys[np.argsort(abs(valleys-CentralIdx))[0]]
        valley= x[valleyidx]
    except:
        pass
    
    try:
        if peaks[np.argsort(abs(peaks-valleyidx))[0]] > peaks[np.argsort(abs(peaks-valleyidx))[1]]:
            peak1idx, peak2idx= peaks[np.argsort(abs(peaks-valleyidx))[0]], peaks[np.argsort(abs(peaks-valleyidx))[1]]
        else:
            peak2idx, peak1idx= peaks[np.argsort(abs(peaks-valleyidx))[0]], peaks[np.argsort(abs(peaks-valleyidx))[1]]
        peak1, peak2= x[peak1idx], x[peak2idx]
    except:
        pass
#    plt.plot(x)
#    plt.plot(np.array([peak1idx, peak2idx]), x[np.array([peak1idx, peak2idx])], 'x')
#    plt.plot(np.array([valleyidx]), x[np.array([valleyidx])], 'x')

    PeakIdx1.append([peak1idx, XIdx])
    PeakIdx2.append([peak2idx, XIdx])
    ValleyIdx.append([valleyidx, XIdx])
    Peak1.append([peak1, XIdx])
    Peak2.append([peak2, XIdx])
    Valley.append([valley, XIdx])
PeakIdx1= np.array(PeakIdx1)
PeakIdx2= np.array(PeakIdx2)
ValleyIdx= np.array(ValleyIdx)
Peak1= np.array(Peak1)
Peak2= np.array(Peak2)
Valley= np.array(Valley)
#%%
plt.pcolormesh(Signal[SignalRangeY[0]:SignalRangeY[1], SignalRangeX[0]:SignalRangeX[1]])
plt.plot(PeakIdx1[:, 1]-SignalRangeX[0], PeakIdx1[:, 0])
plt.plot(PeakIdx2[:, 1]-SignalRangeX[0], PeakIdx2[:, 0])
plt.plot(ValleyIdx[:, 1]-SignalRangeX[0], ValleyIdx[:, 0])

#%%
PeakWidth= np.mean(PeakIdx1[:, 0]-PeakIdx2[:, 0])
PeakToTrough= np.mean(((Peak1[:, 0]+Peak2[:, 0])/2-Valley[:, 0])/Valley[:, 0])
print('PeakWidth', PeakWidth)
print('PeakToTrough', PeakToTrough)